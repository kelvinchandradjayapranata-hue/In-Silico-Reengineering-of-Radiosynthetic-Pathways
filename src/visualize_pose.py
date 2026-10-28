"""Buat viewer 3D interaktif untuk pose ligan di kantung reseptor.

Output berupa HTML mandiri (kecuali library 3Dmol.js dari CDN) yang dapat
dibuka langsung di browser. Viewer menampilkan:

- reseptor sebagai cartoon,
- ligan pose terbaik sebagai stick,
- residu yang berinteraksi sebagai stick,
- garis putus-putus untuk H-bond, kontak hidrofobik, dan pi-pi.

Tidak membutuhkan py3Dmol. Koneksi internet diperlukan saat membuka HTML
untuk memuat 3Dmol.js dari CDN.
"""

from __future__ import annotations

import argparse
import html
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Sequence

from rdkit import Chem

import analyze_interactions as ai

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT_DIR = PROJECT_ROOT / "results" / "visualization"


def _slug(value: str) -> str:
    slug = re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_")
    return slug or "complex"


def _json_for_script(value: Any) -> str:
    return json.dumps(value, ensure_ascii=False).replace("</", "<\\/")


def _residue_groups(complex_data: dict[str, Any]) -> list[dict[str, Any]]:
    groups: dict[tuple[str, str], set[int]] = {}
    for key in ("hbonds", "hydrophobic_contacts", "pi_stacking"):
        for item in complex_data.get(key, []):
            chain = str(item.get("protein_chain", ""))
            resi = item.get("protein_resi")
            if resi is None:
                continue
            groups.setdefault((chain, str(item.get("protein_icode", ""))), set()).add(
                int(resi)
            )
    return [
        {
            "chain": chain,
            "icode": icode,
            "resi": sorted(residues),
        }
        for (chain, icode), residues in sorted(groups.items())
    ]


def _contact_lines(complex_data: dict[str, Any]) -> list[dict[str, Any]]:
    lines: list[dict[str, Any]] = []
    for item in complex_data.get("hbonds", []):
        lines.append(
            {
                "kind": "hbond",
                "start": item["ligand_position"],
                "end": item["protein_position"],
                "label": (
                    f"{item['protein_residue']} {item['protein_atom']}"
                    f"···{item['ligand_atom']}"
                ),
                "detail": (
                    f"{item['distance_angstrom']} Å, "
                    f"{item['angle_deg']}°"
                ),
                "color": "#ffd166",
            }
        )
    for item in complex_data.get("hydrophobic_contacts", []):
        lines.append(
            {
                "kind": "hydrophobic",
                "start": item["representative_ligand_position"],
                "end": item["representative_protein_position"],
                "label": f"{item['protein_residue']} (hidrofobik)",
                "detail": (
                    f"{item['contact_count']} kontak, min "
                    f"{item['min_distance_angstrom']} Å"
                ),
                "color": "#fb923c",
            }
        )
    for item in complex_data.get("pi_stacking", []):
        lines.append(
            {
                "kind": "pi_stacking",
                "start": item["ligand_centroid"],
                "end": item["protein_centroid"],
                "label": (
                    f"{item['protein_ring']} ··· {item['ligand_ring']}"
                ),
                "detail": (
                    f"{item['centroid_distance_angstrom']} Å centroid, "
                    f"{item['min_atom_distance_angstrom']} Å atom-min, "
                    f"{item['angle_deg']}°, {item['stacking_type']} "
                    f"({item['geometry']})"
                ),
                "color": "#e879f9",
            }
        )
    return lines


def _summary_rows(complex_data: dict[str, Any]) -> str:
    summary = complex_data.get("summary", {})
    return (
        f"<li><span>H-bond</span><strong>{summary.get('hbonds', 0)}</strong></li>"
        f"<li><span>Residu hidrofobik</span>"
        f"<strong>{summary.get('hydrophobic_residues', 0)}</strong></li>"
        f"<li><span>Pi-pi</span><strong>{summary.get('pi_stacking', 0)}</strong></li>"
    )


def _render_html(
    complex_data: dict[str, Any],
    *,
    pdb_text: str,
    sdf_text: str,
    generated_at: str,
) -> str:
    receptor = html.escape(str(complex_data.get("receptor", "-")))
    ligand = html.escape(str(complex_data.get("ligand", "-")))
    pdb_id = html.escape(str(complex_data.get("pdb_id") or "-"))
    target = html.escape(str(complex_data.get("target_accession") or "-"))
    chain = html.escape(str(complex_data.get("receptor_chain") or "-"))
    pose_file = html.escape(str(complex_data.get("pose_file") or "-"))
    receptor_file = html.escape(str(complex_data.get("receptor_file") or "-"))

    affinity = complex_data.get("affinity_kcal_mol")
    efficiency = complex_data.get("ligand_efficiency")
    affinity_text = f"{affinity:.2f} kcal/mol" if affinity is not None else "-"
    efficiency_text = (
        f"{efficiency:.3f} kcal/mol per heavy atom"
        if efficiency is not None
        else "-"
    )

    residue_groups = _residue_groups(complex_data)
    contact_lines = _contact_lines(complex_data)

    return f"""<!DOCTYPE html>
<html lang="id">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Pose 3D — {receptor} × {ligand}</title>
<style>
  :root {{
    color-scheme: dark;
    --bg: #070b14;
    --panel: #0d1524;
    --panel-2: #111c30;
    --border: #22324d;
    --text: #e5edf8;
    --muted: #94a3b8;
    --accent: #60a5fa;
    --hbond: #ffd166;
    --hydrophobic: #fb923c;
    --pipi: #e879f9;
  }}
  * {{ box-sizing: border-box; }}
  body {{
    margin: 0;
    background: radial-gradient(circle at 20% 0%, #13213a 0, var(--bg) 55%);
    color: var(--text);
    font-family: Inter, ui-sans-serif, system-ui, -apple-system, "Segoe UI", sans-serif;
    min-height: 100vh;
  }}
  .layout {{
    display: grid;
    grid-template-columns: minmax(300px, 360px) 1fr;
    gap: 18px;
    padding: 18px;
    min-height: 100vh;
  }}
  .panel {{
    background: linear-gradient(180deg, var(--panel), var(--panel-2));
    border: 1px solid var(--border);
    border-radius: 16px;
    padding: 18px;
    box-shadow: 0 18px 50px rgb(0 0 0 / 0.35);
  }}
  h1 {{ font-size: 1.15rem; line-height: 1.35; margin: 0 0 6px; }}
  h2 {{ font-size: .78rem; letter-spacing: .12em; text-transform: uppercase; color: var(--muted); margin: 20px 0 8px; }}
  p {{ color: var(--muted); line-height: 1.55; margin: 8px 0; }}
  code {{ color: #bfdbfe; font-size: .78rem; word-break: break-all; }}
  .badges {{ display: flex; flex-wrap: wrap; gap: 6px; margin-top: 10px; }}
  .badge {{
    border: 1px solid var(--border); border-radius: 999px;
    padding: 4px 8px; font-size: .72rem; color: #cbd5e1;
    background: rgb(15 23 42 / .65);
  }}
  ul.summary {{ list-style: none; padding: 0; margin: 0; display: grid; gap: 8px; }}
  ul.summary li {{ display: flex; justify-content: space-between; gap: 12px; padding: 9px 10px; background: rgb(2 6 23 / .45); border: 1px solid var(--border); border-radius: 10px; }}
  ul.summary span {{ color: var(--muted); font-size: .84rem; }}
  ul.summary strong {{ font-variant-numeric: tabular-nums; }}
  .legend {{ display: grid; gap: 7px; margin-top: 6px; }}
  .legend div {{ display: flex; align-items: center; gap: 8px; font-size: .82rem; color: var(--muted); }}
  .swatch {{ width: 22px; height: 0; border-top: 3px dashed currentColor; }}
  .controls {{ display: grid; grid-template-columns: 1fr 1fr; gap: 8px; margin-top: 14px; }}
  button {{
    border: 1px solid #2b4162; border-radius: 10px; padding: 9px 10px;
    background: #14233b; color: var(--text); font: inherit; font-size: .82rem;
    cursor: pointer;
  }}
  button:hover {{ background: #1b2f4e; border-color: #3b5d8a; }}
  button:focus-visible {{ outline: 2px solid var(--accent); outline-offset: 2px; }}
  .viewer-wrap {{
    position: relative; min-height: calc(100vh - 36px);
    border: 1px solid var(--border); border-radius: 16px; overflow: hidden;
    background: #070b14; box-shadow: 0 18px 50px rgb(0 0 0 / 0.35);
  }}
  #viewer {{ position: absolute; inset: 0; }}
  .hint {{
    position: absolute; left: 12px; bottom: 12px; z-index: 2;
    font-size: .72rem; color: #cbd5e1; background: rgb(2 6 23 / .72);
    border: 1px solid var(--border); border-radius: 8px; padding: 6px 8px;
    pointer-events: none;
  }}
  .note {{ font-size: .76rem; color: #94a3b8; margin-top: 12px; }}
  .note strong {{ color: #fca5a5; }}
  @media (max-width: 900px) {{
    .layout {{ grid-template-columns: 1fr; }}
    .viewer-wrap {{ min-height: 68vh; }}
  }}
</style>
</head>
<body>
<div class="layout">
  <aside class="panel">
    <h1>{receptor}<br><span style="color:var(--accent)">× {ligand}</span></h1>
    <p>Pose docking terbaik dan peta interaksi geometris. Warna residu dan garis mengikuti legenda.</p>
    <div class="badges">
      <span class="badge">PDB {pdb_id}</span>
      <span class="badge">target {target}</span>
      <span class="badge">chain {chain}</span>
      <span class="badge">{html.escape(affinity_text)}</span>
      <span class="badge">LE {html.escape(efficiency_text)}</span>
    </div>

    <h2>Ringkasan interaksi</h2>
    <ul class="summary">{_summary_rows(complex_data)}</ul>

    <h2>Legenda</h2>
    <div class="legend">
      <div style="color:var(--hbond)"><span class="swatch"></span>H-bond (kuning)</div>
      <div style="color:var(--hydrophobic)"><span class="swatch"></span>Kontak hidrofobik (oranye)</div>
      <div style="color:var(--pipi)"><span class="swatch"></span>Stacking pi-pi (magenta)</div>
    </div>

    <h2>Kontrol</h2>
    <div class="controls">
      <button id="resetBtn" type="button">Reset view</button>
      <button id="spinBtn" type="button">Putar</button>
      <button id="surfaceBtn" type="button">Surface pocket</button>
      <button id="contactsBtn" type="button">Sembunyikan kontak</button>
    </div>

    <h2>Berkas</h2>
    <p><code>{pose_file}</code></p>
    <p><code>{receptor_file}</code></p>
    <p class="note">Dibuat {html.escape(generated_at)}. Viewer memuat 3Dmol.js dari CDN; perlu koneksi internet. Interaksi adalah estimasi geometris, bukan energi ikat.</p>
    <p class="note"><strong>Catatan:</strong> reseptor adalah homolog surrogate, bukan struktur target C. cladosporioides.</p>
  </aside>

  <main class="viewer-wrap" aria-label="Viewer 3D pose ligan">
    <div id="viewer"></div>
    <div class="hint">Drag untuk rotasi · scroll untuk zoom · klik kanan untuk pan</div>
  </main>
</div>

<script src="https://3Dmol.org/build/3Dmol-min.js"></script>
<script type="text/plain" id="pdb-data">{pdb_text}</script>
<script type="text/plain" id="sdf-data">{sdf_text}</script>
<script>
const contactData = {_json_for_script(contact_lines)};
const residueGroups = {_json_for_script(residue_groups)};

function xyz(values) {{
  return {{x: values[0], y: values[1], z: values[2]}};
}}

const viewer = $3Dmol.createViewer("viewer", {{
  backgroundColor: "#070b14",
  antialias: true,
}});

const pdbText = document.getElementById("pdb-data").textContent;
const sdfText = document.getElementById("sdf-data").textContent;
viewer.addModel(pdbText, "pdb");
viewer.addModel(sdfText, "sdf");

viewer.setStyle({{model: 0}}, {{
  cartoon: {{color: "#8fa8c8", opacity: 0.88}},
}});

residueGroups.forEach((group) => {{
  const selection = {{model: 0, resi: group.resi}};
  if (group.chain) selection.chain = group.chain;
  viewer.addStyle(selection, {{
    stick: {{colorscheme: "grayCarbon", radius: 0.16}},
    cartoon: {{color: "#f2b84b", opacity: 0.95}},
  }});
}});

viewer.setStyle({{model: 1}}, {{
  stick: {{colorscheme: "greenCarbon", radius: 0.26}},
}});
viewer.addStyle({{model: 1, elem: "O"}}, {{stick: {{color: "#ff6b6b", radius: 0.30}}}});
viewer.addStyle({{model: 1, elem: "N"}}, {{stick: {{color: "#60a5fa", radius: 0.30}}}});

let contactShapes = [];
let contactLabels = [];

function drawContacts() {{
  contactData.forEach((contact) => {{
    const shape = viewer.addLine({{
      start: xyz(contact.start),
      end: xyz(contact.end),
      color: contact.color,
      dashed: true,
      linewidth: 2.5,
    }});
    contactShapes.push(shape);
  }});

  contactData.filter((item) => item.kind !== "hydrophobic").forEach((contact) => {{
    const mid = [
      (contact.start[0] + contact.end[0]) / 2,
      (contact.start[1] + contact.end[1]) / 2,
      (contact.start[2] + contact.end[2]) / 2,
    ];
    const label = viewer.addLabel(
      contact.label + "\\n" + contact.detail,
      {{
        position: xyz(mid),
        backgroundColor: "#0b1220",
        backgroundOpacity: 0.78,
        borderColor: contact.color,
        borderThickness: 1,
        fontColor: "#f8fafc",
        fontSize: 10,
        showBackground: true,
      }}
    );
    contactLabels.push(label);
  }});
}}

function clearContacts() {{
  contactShapes.forEach((shape) => viewer.removeShape(shape));
  contactLabels.forEach((label) => viewer.removeShape(label));
  contactShapes = [];
  contactLabels = [];
}}

drawContacts();

let surfaceShape = null;
let contactsVisible = true;
let spinning = false;

function resetView() {{
  viewer.zoomTo({{model: 1}});
  viewer.zoom(0.78);
  viewer.render();
}}

document.getElementById("resetBtn").addEventListener("click", resetView);

document.getElementById("spinBtn").addEventListener("click", () => {{
  spinning = !spinning;
  if (spinning) {{
    viewer.spin("y", 1);
  }} else {{
    viewer.spin(false);
  }}
  document.getElementById("spinBtn").textContent = spinning ? "Berhenti" : "Putar";
  viewer.render();
}});

document.getElementById("surfaceBtn").addEventListener("click", () => {{
  if (surfaceShape) {{
    viewer.removeSurface(surfaceShape);
    surfaceShape = null;
  }} else {{
    surfaceShape = viewer.addSurface(
      $3Dmol.SurfaceType.VDW,
      {{opacity: 0.18, color: "#60a5fa"}},
      {{model: 0, within: {{distance: 5.0, sel: {{model: 1}}}}}}
    );
  }}
  viewer.render();
}});

document.getElementById("contactsBtn").addEventListener("click", () => {{
  contactsVisible = !contactsVisible;
  if (contactsVisible) {{
    drawContacts();
  }} else {{
    clearContacts();
  }}
  document.getElementById("contactsBtn").textContent =
    contactsVisible ? "Sembunyikan kontak" : "Tampilkan kontak";
  viewer.render();
}});

resetView();
viewer.render();
</script>
</body>
</html>
"""


def _render_index(
    entries: Sequence[dict[str, Any]],
    *,
    generated_at: str,
    errors: Sequence[dict[str, Any]],
) -> str:
    rows = []
    for entry in entries:
        rows.append(
            "<tr>"
            f"<td><a href=\"{html.escape(entry['file'])}\">"
            f"{html.escape(entry['receptor'])} × {html.escape(entry['ligand'])}"
            "</a></td>"
            f"<td>{entry.get('hbonds', 0)}</td>"
            f"<td>{entry.get('hydrophobic_residues', 0)}</td>"
            f"<td>{entry.get('pi_stacking', 0)}</td>"
            "</tr>"
        )
    error_html = ""
    if errors:
        error_html = "<h2>Kompleks yang dilewati</h2><ul>" + "".join(
            f"<li>{html.escape(str(item['receptor']))} × "
            f"{html.escape(str(item['ligand']))}: "
            f"{html.escape(str(item['error']))}</li>"
            for item in errors
        ) + "</ul>"

    return f"""<!DOCTYPE html>
<html lang="id">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Visualisasi Pose Docking</title>
<style>
  body {{ margin: 0; padding: 32px; background: #070b14; color: #e5edf8; font-family: ui-sans-serif, system-ui, -apple-system, "Segoe UI", sans-serif; }}
  main {{ max-width: 980px; margin: 0 auto; }}
  h1 {{ font-size: 1.5rem; }}
  p {{ color: #94a3b8; line-height: 1.6; }}
  table {{ width: 100%; border-collapse: collapse; margin-top: 18px; background: #0d1524; border: 1px solid #22324d; border-radius: 12px; overflow: hidden; }}
  th, td {{ padding: 11px 12px; text-align: left; border-bottom: 1px solid #22324d; }}
  th {{ background: #111c30; color: #cbd5e1; font-size: .82rem; text-transform: uppercase; letter-spacing: .06em; }}
  tr:last-child td {{ border-bottom: 0; }}
  a {{ color: #7dd3fc; text-decoration: none; }}
  a:hover {{ text-decoration: underline; }}
  .note {{ font-size: .84rem; color: #94a3b8; }}
</style>
</head>
<body>
<main>
  <h1>Visualisasi 3D Pose Docking</h1>
  <p>Klik kompleks untuk membuka viewer interaktif. Viewer menampilkan pose ligan,
  residu yang berinteraksi, dan garis kontak. Dibuat {html.escape(generated_at)}.</p>
  <table>
    <thead><tr><th>Kompleks</th><th>H-bond</th><th>Residu hidrofobik</th><th>Pi-pi</th></tr></thead>
    <tbody>{''.join(rows)}</tbody>
  </table>
  {error_html}
  <p class="note">Analisis interaksi bersifat geometris. Reseptor adalah homolog
  surrogate; hasil bukan bukti interaksi target sebenarnya.</p>
</main>
</body>
</html>
"""


def _load_payload(args: argparse.Namespace) -> dict[str, Any]:
    if args.interactions:
        payload = ai._read_json(args.interactions)
        if not payload.get("complexes"):
            raise ai.AnalysisError(
                f"{args.interactions} tidak memiliki data complexes."
            )
        return payload
    return ai.build_payload(
        args.docking,
        receptors=args.receptor,
        ligands=args.ligand,
        chain=args.chain,
        cutoff_hbond=args.cutoff_hbond,
        min_hbond_angle=args.min_hbond_angle,
        cutoff_hydrophobic=args.cutoff_hydrophobic,
        cutoff_pipi=args.cutoff_pipi,
        max_pipi_centroid_offset=args.max_pipi_centroid,
        max_pipi_atom_distance=args.cutoff_pipi_atom,
    )


def generate_viewers(
    payload: dict[str, Any],
    out_dir: Path,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    out_dir.mkdir(parents=True, exist_ok=True)
    generated_at = datetime.now(timezone.utc).isoformat()
    entries: list[dict[str, Any]] = []
    errors: list[dict[str, Any]] = []

    for complex_data in payload.get("complexes", []):
        receptor = str(complex_data.get("receptor", "receptor"))
        ligand = str(complex_data.get("ligand", "ligand"))
        try:
            pose_path = ai.resolve_path(complex_data["pose_file"])
            receptor_path = ai.resolve_path(complex_data["receptor_file"])
            ligand_mol = ai.load_ligand_mol(pose_path)
            heavy_mol = Chem.RemoveHs(ligand_mol)
            sdf_text = Chem.MolToMolBlock(heavy_mol, confId=0)
            pdb_text = ai.receptor_pdb_text(
                receptor_path,
                chain_id=complex_data.get("receptor_chain"),
            )
            filename = f"{_slug(receptor)}__{_slug(ligand)}.html"
            (out_dir / filename).write_text(
                _render_html(
                    complex_data,
                    pdb_text=pdb_text,
                    sdf_text=sdf_text,
                    generated_at=generated_at,
                ),
                encoding="utf-8",
            )
            summary = complex_data.get("summary", {})
            entries.append(
                {
                    "receptor": receptor,
                    "ligand": ligand,
                    "file": filename,
                    "hbonds": summary.get("hbonds", 0),
                    "hydrophobic_residues": summary.get(
                        "hydrophobic_residues", 0
                    ),
                    "pi_stacking": summary.get("pi_stacking", 0),
                }
            )
        except (
            ai.AnalysisError,
            OSError,
            ValueError,
            KeyError,
            TypeError,
        ) as exc:
            errors.append(
                {
                    "receptor": receptor,
                    "ligand": ligand,
                    "error": str(exc),
                }
            )

    (out_dir / "index.html").write_text(
        _render_index(entries, generated_at=generated_at, errors=errors),
        encoding="utf-8",
    )
    return entries, errors


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Buat HTML 3D interaktif untuk pose ligan dan peta interaksinya."
        )
    )
    parser.add_argument(
        "--interactions",
        type=Path,
        default=None,
        help=(
            "Pakai results/interactions.json yang sudah ada; jika tidak, "
            "analisis dijalankan dari --docking."
        ),
    )
    parser.add_argument("--docking", type=Path, default=ai.DEFAULT_DOCKING)
    parser.add_argument("--receptor", nargs="+", default=None)
    parser.add_argument("--ligand", nargs="+", default=["DHN_monomer"])
    parser.add_argument("--chain", default=None)
    parser.add_argument("--cutoff-hbond", type=float, default=3.5)
    parser.add_argument("--min-hbond-angle", type=float, default=120.0)
    parser.add_argument("--cutoff-hydrophobic", type=float, default=4.0)
    parser.add_argument("--cutoff-pipi", type=float, default=5.5)
    parser.add_argument("--max-pipi-centroid", type=float, default=6.5)
    parser.add_argument("--cutoff-pipi-atom", type=float, default=4.0)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        payload = _load_payload(args)
        entries, errors = generate_viewers(payload, args.out_dir)
    except (
        ai.AnalysisError,
        OSError,
        ValueError,
        KeyError,
        TypeError,
    ) as exc:
        print(f"[ERROR] {exc}")
        return 1

    if not entries:
        print("[ERROR] Tidak ada viewer yang berhasil dibuat.")
        for error in errors:
            print(f"    {error['receptor']} x {error['ligand']}: {error['error']}")
        return 1

    print(f"[SUCCESS] {len(entries)} viewer ditulis ke: {args.out_dir}")
    print(f"    Index: {args.out_dir / 'index.html'}")
    for entry in entries:
        print(f"    - {entry['file']}")
    for error in errors:
        print(f"    [WARN] {error['receptor']} x {error['ligand']}: {error['error']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
