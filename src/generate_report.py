"""Kompilasi artefak hasil riset menjadi results/summary_report.md."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Sequence

PROJECT_ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = PROJECT_ROOT / "results"
RAW_DIR = PROJECT_ROOT / "data" / "raw"
SUMMARY_PATH = RESULTS_DIR / "summary_report.md"

EVIDENCE_LABELS = {
    "direct_sequence": "Urutan langsung",
    "ortholog_reference": "Ortholog referensi",
}


def _read_json(path: Path) -> dict[str, Any] | None:
    if not path.is_file():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return None


def _fmt(value: Any, digits: int = 3) -> str:
    if isinstance(value, bool):
        return "ya" if value else "tidak"
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        return f"{value:.{digits}f}"
    if value is None:
        return "-"
    return str(value)


def _table(headers: Sequence[str], rows: Sequence[Sequence[Any]]) -> str:
    lines = [
        "| " + " | ".join(headers) + " |",
        "|" + "|".join(["---"] * len(headers)) + "|",
    ]
    for row in rows:
        lines.append("| " + " | ".join(_fmt(cell) for cell in row) + " |")
    return "\n".join(lines)


def _section_fba(data: dict[str, Any] | None) -> list[str]:
    if not data:
        return []

    pathway = data.get("pathway_fluxes", {})
    balance = data.get("balance", {})
    lines = [
        "## 1. Simulasi FBA Radiosintesis",
        "",
        f"- Model: `{data.get('model', '-')}` (versi {data.get('model_version', '-')})",
        f"- Status solver: `{data.get('status', '-')}`",
        f"- Objective total ekspor melanin: **{_fmt(data.get('objective_value'))}**",
        f"- DHN-melanin export: **{_fmt(pathway.get('dhn_melanin_export'))}**",
        f"- DOPA-melanin export: **{_fmt(pathway.get('dopa_melanin_export'))}**",
        f"- Glucose terpakai: {_fmt(balance.get('glucose_used'))}",
        f"- Oksigen terpakai: {_fmt(balance.get('oxygen_used'))}",
        f"- Tirosin terpakai: {_fmt(balance.get('tyrosine_used'))}",
        f"- Elektron radiasi terpakai: {_fmt(balance.get('radiation_electron_used'))}",
        f"- Energi radiasi terpakai: {_fmt(balance.get('radiation_energy_used'))}",
        f"- Kontribusi ATP radiasi terhadap demand ATP: "
        f"**{_fmt(100.0 * balance.get('radiation_atp_share_of_demand', 0.0), 1)}%**",
        "",
    ]

    key_fluxes = data.get("fluxes") or {}
    if key_fluxes:
        lines.extend(
            [
                "### Fluks kunci",
                "",
                _table(
                    ["Reaksi", "Fluks"],
                    [[reaction_id, value] for reaction_id, value in key_fluxes.items()],
                ),
                "",
            ]
        )

    binding = data.get("binding_constraints") or []
    lines.extend(["### Batasan yang mengikat", ""])
    if binding:
        lines.extend([f"- {item}" for item in binding])
    else:
        lines.append("- Tidak ada batas uptake/input yang mengikat pada kondisi ini.")
    lines.append("")

    warnings = data.get("warnings") or []
    if warnings:
        lines.extend(["### Peringatan model", ""])
        lines.extend([f"- {item}" for item in warnings])
        lines.append("")

    assumptions = data.get("assumptions") or []
    if assumptions:
        lines.extend(["### Asumsi model", ""])
        lines.extend([f"- {item}" for item in assumptions])
        lines.append("")

    return lines


def _section_chem(data: dict[str, Any] | None) -> list[str]:
    if not data:
        return []

    monomers = data.get("monomers", {})
    if not monomers:
        return []

    rows = []
    for model_name, props in monomers.items():
        rows.append(
            [
                props.get("common_name", model_name),
                props.get("canonical_smiles", "-"),
                props.get("formula", "-"),
                props.get("molecular_weight_g_mol"),
                props.get("logP"),
                props.get("tpsa_A2"),
                props.get("h_bond_donors"),
                props.get("h_bond_acceptors"),
                props.get("gasteiger_charge_min"),
                props.get("gasteiger_charge_max"),
            ]
        )

    lines = [
        "## 2. Deskriptor Monomer Melanin",
        "",
        _table(
            [
                "Monomer",
                "SMILES kanonik",
                "Formula",
                "MW (g/mol)",
                "logP",
                "TPSA (Å²)",
                "HBD",
                "HBA",
                "Gasteiger min",
                "Gasteiger max",
            ],
            rows,
        ),
        "",
        f"Metode: {data.get('method', '-')}",
        "",
    ]

    assumptions = data.get("assumptions") or []
    if assumptions:
        lines.extend(["### Asumsi", ""])
        lines.extend([f"- {item}" for item in assumptions])
        lines.append("")
    return lines


def _section_provenance(manifest: dict[str, Any] | None) -> list[str]:
    if not manifest:
        return []

    targets = manifest.get("targets", [])
    rows = []
    for target in targets:
        rows.append(
            [
                target.get("target_id", "-"),
                target.get("gene", "-"),
                target.get("accession", "-"),
                target.get("expected_organism", "-"),
                EVIDENCE_LABELS.get(
                    target.get("evidence_level", ""), target.get("evidence_level", "-")
                ),
                target.get("length", "-"),
                target.get("status", "-"),
            ]
        )

    lines = [
        "## 3. Data dan Provenance",
        "",
        _table(
            [
                "Target",
                "Gen",
                "Aksesi",
                "Organisme",
                "Level bukti",
                "Panjang",
                "Status",
            ],
            rows,
        ),
        "",
    ]

    caveats = manifest.get("caveats") or []
    if caveats:
        lines.extend(["### Catatan ilmiah", ""])
        lines.extend([f"- {item}" for item in caveats])
        lines.append("")

    assemblies = manifest.get("reference_assemblies") or []
    if assemblies:
        lines.extend(["### Assembly referensi (tidak diunduh)", ""])
        for assembly in assemblies:
            lines.append(
                f"- `{assembly.get('accession', '-')}` — "
                f"{assembly.get('organism', '-')}: {assembly.get('description', '-')}"
            )
        lines.append("")

    return lines


def _section_sensitivity(data: dict[str, Any] | None) -> list[str]:
    if not data:
        return []

    rows = data.get("rows", [])
    if not rows:
        return []

    table_rows = []
    for row in rows:
        if row.get("status") == "error":
            table_rows.append(
                [
                    row.get("electron_limit"),
                    row.get("energy_limit"),
                    "error",
                    "-",
                    "-",
                    "-",
                    "-",
                    "-",
                ]
            )
        else:
            table_rows.append(
                [
                    row.get("electron_limit"),
                    row.get("energy_limit"),
                    row.get("status"),
                    row.get("objective_value"),
                    row.get("dhn_melanin_export"),
                    row.get("dopa_melanin_export"),
                    row.get("radiation_electron_used"),
                    row.get("radiation_energy_used"),
                ]
            )

    return [
        "## 4. Analisis Sensitivitas Input Radiasi",
        "",
        _table(
            [
                "e- limit",
                "Energi limit",
                "Status",
                "Objective",
                "DHN",
                "DOPA",
                "e- terpakai",
                "Energi terpakai",
            ],
            table_rows,
        ),
        "",
        "Catatan: status `error` (biasanya pada e- = 0 dan energi = 0) berarti "
        "kebutuhan ATP maintenance tidak dapat dipenuhi tanpa input radiasi.",
        "",
        f"{data.get('interpretation', '')}",
        "",
    ]


def _section_processed_stats(
    processed_files: Sequence[Path], root: Path = PROJECT_ROOT
) -> list[str]:
    if not processed_files:
        return []

    lines = ["## 5. Statistik Data Terproses", ""]
    for path in processed_files:
        data = _read_json(path)
        display_path = path.relative_to(root) if path.is_relative_to(root) else path
        lines.append(f"### `{display_path}`")
        lines.append("")
        if data is None:
            lines.append("Gagal membaca JSON.")
        else:
            lines.append("```json")
            lines.append(json.dumps(data, indent=2, ensure_ascii=False))
            lines.append("```")
        lines.append("")
    return lines


def generate_report(
    project_root: Path | None = None,
    output_path: Path | None = None,
) -> Path:
    """Bangun laporan markdown dari artefak yang tersedia."""
    root = Path(project_root) if project_root else PROJECT_ROOT
    results_dir = root / "results"
    raw_dir = root / "data" / "raw"
    processed_dir = root / "data" / "processed"
    destination = Path(output_path) if output_path else results_dir / "summary_report.md"

    fba = _read_json(results_dir / "metabolic_fba.json")
    chem = _read_json(results_dir / "melanin_chem.json")
    sensitivity = _read_json(results_dir / "fba_sensitivity.json")
    manifest = _read_json(raw_dir / "manifest.json")
    processed_files = sorted(processed_dir.glob("*_stats.json")) if processed_dir.is_dir() else []

    expected = [
        results_dir / "metabolic_fba.json",
        results_dir / "melanin_chem.json",
        raw_dir / "manifest.json",
    ]
    missing = [path for path in expected if not path.is_file()]

    lines: list[str] = [
        "# Ringkasan Hasil Riset Radiosintesis Melanin",
        "",
        f"> Dibuat: {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}",
        ">",
        "> Ringkasan ini mengompilasi artefak komputasi yang tersedia. Skor model",
        "> merupakan estimasi theoretical/in silico dan bukan pengganti validasi",
        "> eksperimental.",
        "",
    ]

    sections: list[list[str]] = [
        _section_fba(fba),
        _section_chem(chem),
        _section_provenance(manifest),
        _section_sensitivity(sensitivity),
        _section_processed_stats(processed_files, root),
    ]
    for section in sections:
        if section:
            lines.extend(section)

    lines.extend(["## Artefak yang belum tersedia / opsional", ""])
    if missing:
        lines.extend([f"- `{path.relative_to(root)}`" for path in missing])
    else:
        lines.append("- Semua artefak inti tersedia.")
    lines.append(
        "- `results/docking_results.json` — opsional; modul docking belum "
        "diimplementasikan."
    )
    lines.append(
        "- `data/processed/*_stats.json` — opsional; belum ada data terproses."
    )
    lines.append("")

    destination.parent.mkdir(parents=True, exist_ok=True)
    temp_path = destination.with_suffix(destination.suffix + ".tmp")
    temp_path.write_text("\n".join(lines), encoding="utf-8")
    temp_path.replace(destination)
    return destination


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Buat laporan ringkasan markdown.")
    parser.add_argument("--output", type=Path, default=SUMMARY_PATH)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    report_path = generate_report(output_path=args.output)
    print(f"[SUCCESS] Laporan ringkasan riset ditulis ke: {report_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
