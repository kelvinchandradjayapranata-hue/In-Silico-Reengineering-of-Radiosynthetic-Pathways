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
        lines.extend(["### Assembly referensi", ""])
        for assembly in assemblies:
            lines.append(
                f"- `{assembly.get('accession', '-')}` — "
                f"{assembly.get('organism', '-')}: {assembly.get('description', '-')}"
            )
        lines.extend(
            [
                "",
                "Genome assembly diunduh dan diverifikasi checksum-nya oleh "
                "`src/blast_validation.py` sesuai kebutuhan.",
                "",
            ]
        )

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


def _section_blast(data: dict[str, Any] | None, number: int) -> list[str]:
    if not data:
        return []

    summary = data.get("summary", {})
    counts = summary.get("classifications", {})
    assemblies = [
        assembly.get("accession", "-") for assembly in data.get("assemblies", [])
    ]
    best_hits = data.get("best_hits", [])
    parameters = data.get("parameters", {})

    lines = [
        f"## {number}. Validasi Homologi Ortholog (tblastn)",
        "",
        f"- Metode: {data.get('method', '-')}",
        f"- Versi BLAST: `{data.get('blast_version', '-')}`",
        f"- Parameter: e-value <= {parameters.get('evalue', '-')}, "
        f"max_target_seqs {parameters.get('max_target_seqs', '-')}, "
        f"max_hsps {parameters.get('max_hsps', '-')}, "
        f"SEG {parameters.get('seg', '-')}",
        f"- Total hit mentah: {data.get('hits_total', '-')}",
        "",
        "### Ringkasan klasifikasi",
        "",
        _table(
            ["Assembly", "Strong", "Moderate", "Weak", "Tidak terdeteksi"],
            [
                [
                    assembly,
                    summary.get("per_assembly", {}).get(assembly, {}).get("strong", 0),
                    summary.get("per_assembly", {}).get(assembly, {}).get("moderate", 0),
                    summary.get("per_assembly", {}).get(assembly, {}).get("weak", 0),
                    summary.get("per_assembly", {})
                    .get(assembly, {})
                    .get("not_detected", 0),
                ]
                for assembly in assemblies
            ],
        ),
        "",
        "### Hit terbaik per query",
        "",
    ]

    rows = []
    for query in data.get("queries", []):
        query_id = query.get("query_id", "-")
        row = [query_id.replace("|", "\\|"), query.get("gene", "-")]
        for assembly in assemblies:
            hit = next(
                (
                    item
                    for item in best_hits
                    if item.get("query_id") == query_id
                    and item.get("assembly") == assembly
                ),
                None,
            )
            if hit is None or hit.get("classification") == "not_detected":
                row.append("tidak terdeteksi")
            else:
                row.append(
                    "{classification} ({pident:.1f}%, cov {coverage:.2f})".format(
                        classification=hit.get("classification", "-"),
                        pident=hit.get("pident") or 0.0,
                        coverage=hit.get("query_coverage") or 0.0,
                    )
                )
        rows.append(row)

    lines.extend([_table(["Query", "Gen"] + assemblies, rows), ""])
    lines.extend(
        [
            "Detail subject, koordinat, dan e-value tersedia di "
            "`results/blast_validation.json` dan `results/blast_hits.tsv`.",
            "",
        ]
    )

    interpretation = data.get("interpretation")
    if interpretation:
        lines.extend([interpretation, ""])

    caveats = data.get("caveats") or []
    if caveats:
        lines.extend(["### Keterbatasan validasi BLAST", ""])
        lines.extend([f"- {item}" for item in caveats])
        lines.append("")
    return lines


def _section_docking(data: dict[str, Any] | None, number: int) -> list[str]:
    if not data:
        return []

    engine = data.get("engine", {})
    summary = data.get("summary", {})
    aggregation = summary.get("aggregation", {})
    multi_seed = aggregation.get("method") == "multi_seed"
    seeds = summary.get("seeds") or aggregation.get("seeds") or []
    title_kind = "Docking Multi-Seed" if multi_seed else "Docking Pilot"
    lines = [
        f"## {number}. {title_kind} DHN-melanin (AutoDock Vina)",
        "",
        f"- Status: `{data.get('status', '-')}`",
        f"- Engine: {engine.get('vina', '-')}; Meeko {engine.get('meeko', '-')}; "
        f"RDKit {engine.get('rdkit', '-')}",
    ]
    if seeds:
        lines.append(
            f"- Seed: {', '.join(str(seed) for seed in seeds)}; "
            f"{summary.get('runs_ok', 0)}/{summary.get('runs_total', 0)} "
            "kombinasi reseptor x ligan berhasil."
        )
    lines.append("")

    notice = data.get("surrogate_notice")
    if notice:
        lines.extend([f"> **Catatan:** {notice}", ""])

    matrix = summary.get("affinity_matrix", {})
    ligand_names = [
        ligand.get("name", "-") for ligand in data.get("ligands", [])
    ]
    affinity_label = (
        "Afinitas rata-rata antar-seed (kcal/mol)"
        if multi_seed
        else "Afinitas terbaik (kcal/mol)"
    )
    if matrix and ligand_names:
        rows = []
        for receptor, values in matrix.items():
            row = [receptor]
            for ligand_name in ligand_names:
                value = values.get(ligand_name)
                row.append(f"{value:.2f}" if value is not None else "-")
            rows.append(row)
        lines.extend(
            [
                f"### {affinity_label}",
                "",
                _table(["Reseptor (surrogate)"] + ligand_names, rows),
                "",
            ]
        )

    efficiency = summary.get("ligand_efficiency_matrix", {})
    if efficiency and ligand_names:
        efficiency_rows = []
        for receptor, values in efficiency.items():
            row = [receptor]
            for ligand_name in ligand_names:
                value = values.get(ligand_name)
                row.append(f"{value:.3f}" if value is not None else "-")
            efficiency_rows.append(row)
        efficiency_label = (
            "Ligand efficiency rata-rata antar-seed "
            "(kcal/mol per heavy atom)"
            if multi_seed
            else "Ligand efficiency (kcal/mol per heavy atom)"
        )
        basis = "rata-rata afinitas" if multi_seed else "afinitas terbaik"
        lines.extend(
            [
                f"### {efficiency_label}",
                "",
                _table(["Reseptor (surrogate)"] + ligand_names, efficiency_rows),
                "",
                f"Ligand efficiency = {basis} / jumlah atom berat. Metrik ini "
                "sebagian mengoreksi bias ukuran ligan pada skor Vina.",
                "",
            ]
        )
        if multi_seed:
            lines.extend(
                [
                    "Matriks nilai terbaik dan simpangan baku antar-seed "
                    "tersedia di `results/docking_results.json`.",
                    "",
                ]
            )

    receptor_rows = []
    for receptor in data.get("receptors", []):
        if receptor.get("status") == "error":
            continue
        receptor_rows.append(
            [
                receptor.get("name", "-"),
                receptor.get("pdb_id", "-"),
                receptor.get("target_accession", "-"),
                (
                    f"{receptor.get('identity_percent', '-')}% / "
                    f"{receptor.get('coverage_percent', '-')}%"
                ),
                receptor.get("box_center_source", "-"),
            ]
        )
    if receptor_rows:
        lines.extend(
            [
                "### Reseptor",
                "",
                _table(
                    [
                        "Reseptor",
                        "PDB",
                        "Target",
                        "Identitas / coverage",
                        "Sumber pusat kotak",
                    ],
                    receptor_rows,
                ),
                "",
            ]
        )

    interpretation = data.get("interpretation")
    if interpretation:
        lines.extend([interpretation, ""])

    caveats = data.get("caveats") or []
    if caveats:
        lines.extend(["### Keterbatasan docking", ""])
        lines.extend([f"- {item}" for item in caveats])
        lines.append("")
    return lines


def _section_processed_stats(
    processed_files: Sequence[Path],
    root: Path = PROJECT_ROOT,
    number: int = 5,
) -> list[str]:
    if not processed_files:
        return []

    lines = [f"## {number}. Statistik Data Terproses", ""]
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
    blast = _read_json(results_dir / "blast_validation.json")
    docking = _read_json(results_dir / "docking_results.json")
    manifest = _read_json(raw_dir / "manifest.json")
    processed_files = sorted(processed_dir.glob("*_stats.json")) if processed_dir.is_dir() else []

    expected = [
        results_dir / "metabolic_fba.json",
        results_dir / "melanin_chem.json",
        raw_dir / "manifest.json",
    ]
    missing = [path for path in expected if not path.is_file()]
    optional_missing = []
    if blast is None:
        optional_missing.append(
            "- `results/blast_validation.json` — jalankan "
            "`python src/blast_validation.py` untuk validasi homologi."
        )
    if docking is None:
        optional_missing.append(
            "- `results/docking_results.json` — jalankan "
            "`python src/melanin_docking.py` untuk docking pilot."
        )
    if not processed_files:
        optional_missing.append(
            "- `data/processed/*_stats.json` — belum ada data terproses."
        )

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
    ]
    next_number = sum(1 for section in sections if section) + 1
    if blast:
        sections.append(_section_blast(blast, next_number))
        next_number += 1
    if docking:
        sections.append(_section_docking(docking, next_number))
        next_number += 1
    sections.append(_section_processed_stats(processed_files, root, next_number))
    for section in sections:
        if section:
            lines.extend(section)

    lines.extend(["## Artefak yang belum tersedia / opsional", ""])
    if missing:
        lines.extend([f"- `{path.relative_to(root)}`" for path in missing])
    else:
        lines.append("- Semua artefak inti tersedia.")
    if optional_missing:
        lines.extend(optional_missing)
    else:
        lines.append("- Semua artefak opsional yang diharapkan tersedia.")
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
