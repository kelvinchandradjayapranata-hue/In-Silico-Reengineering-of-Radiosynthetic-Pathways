"""Generate README.md dari status artefak proyek."""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Sequence

PROJECT_ROOT = Path(__file__).resolve().parents[1]
README_PATH = PROJECT_ROOT / "README.md"
RESULTS_DIR = PROJECT_ROOT / "results"
RAW_DIR = PROJECT_ROOT / "data" / "raw"

EVIDENCE_LABELS = {
    "direct_sequence": "urutan langsung",
    "ortholog_reference": "ortholog referensi",
}


def _read_json(path: Path) -> dict[str, Any] | None:
    if not path.is_file():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return None


def _provenance_table(manifest: dict[str, Any] | None) -> str:
    if not manifest or not manifest.get("targets"):
        return (
            "Belum ada manifest data. Jalankan `python src/ncbi_fetcher.py --targets all`."
        )

    lines = [
        "| Target | Gen | Aksesi | Organisme | Level bukti |",
        "|---|---|---|---|---|",
    ]
    for target in manifest["targets"]:
        evidence = EVIDENCE_LABELS.get(
            target.get("evidence_level", ""), target.get("evidence_level", "-")
        )
        lines.append(
            "| {target_id} | {gene} | {accession} | {organism} | {evidence} |".format(
                target_id=target.get("target_id", "-"),
                gene=target.get("gene", "-"),
                accession=target.get("accession", "-"),
                organism=target.get("expected_organism", "-"),
                evidence=evidence,
            )
        )
    return "\n".join(lines)


def _reference_assembly_lines(manifest: dict[str, Any] | None) -> str:
    assemblies = manifest.get("reference_assemblies", []) if manifest else []
    if not assemblies:
        return "Belum ada metadata assembly referensi."
    lines = [
        f"- `{assembly.get('accession', '-')}` — "
        f"{assembly.get('organism', '-')}: {assembly.get('description', '-')}"
        for assembly in assemblies
    ]
    lines.extend(
        [
            "",
            "Genome assembly diunduh dan diverifikasi checksum-nya oleh "
            "`src/blast_validation.py` sesuai kebutuhan.",
        ]
    )
    return "\n".join(lines)


def _status_line(fba: dict[str, Any] | None) -> str:
    if not fba:
        return "FBA belum dijalankan."
    pathway = fba.get("pathway_fluxes", {})
    return (
        f"FBA `{fba.get('status', 'unknown')}`; DHN-melanin export "
        f"{pathway.get('dhn_melanin_export', 0):.3f}, DOPA-melanin export "
        f"{pathway.get('dopa_melanin_export', 0):.3f}."
    )


def _blast_status_line(blast: dict[str, Any] | None) -> str:
    if not blast:
        return "belum dijalankan (`python src/blast_validation.py`)."
    summary = blast.get("summary", {})
    counts = summary.get("classifications", {})
    total = summary.get("combinations_total", 0)
    assemblies = len(blast.get("assemblies", []))
    queries = len(blast.get("queries", []))
    return (
        f"tblastn pada {assemblies} assembly; {counts.get('strong', 0)} strong, "
        f"{counts.get('moderate', 0)} moderate, {counts.get('weak', 0)} weak, "
        f"{counts.get('not_detected', 0)} tidak terdeteksi dari {total} "
        f"kombinasi query x assembly ({queries} query ortholog)."
    )


def _docking_status_line(docking: dict[str, Any] | None) -> str:
    if not docking:
        return "belum dijalankan (`python src/melanin_docking.py`)."
    status = docking.get("status", "unknown")
    summary = docking.get("summary", {})
    runs_ok = summary.get("runs_ok", 0)
    runs_total = summary.get("runs_total", 0)
    seed_runs = summary.get("seed_runs_total")
    seed_suffix = f" ({seed_runs} seed-run)" if seed_runs else ""
    if status == "prepared":
        return (
            "persiapan ligan/reseptor selesai; Vina belum dijalankan "
            "(`--prepare-only`)."
        )
    ok_rows = [
        row for row in docking.get("results", []) if row.get("status") == "ok"
    ]
    if not ok_rows:
        return (
            f"{runs_ok}/{runs_total} kombinasi berhasil{seed_suffix}; "
            "periksa log docking."
        )
    best = min(ok_rows, key=lambda row: row["best_affinity_kcal_mol"])
    le_rows = [
        row for row in ok_rows if row.get("ligand_efficiency") is not None
    ]
    if not le_rows:
        return (
            f"{runs_ok}/{runs_total} kombinasi berhasil{seed_suffix}; afinitas "
            f"terbaik {best['best_affinity_kcal_mol']:.2f} kcal/mol untuk "
            f"{best['ligand']} pada {best['receptor']} (reseptor surrogate PDB)."
        )
    best_le = min(le_rows, key=lambda row: row["ligand_efficiency"])
    return (
        f"{runs_ok}/{runs_total} kombinasi berhasil{seed_suffix}; skor mentah "
        f"terbaik {best['best_affinity_kcal_mol']:.2f} kcal/mol untuk "
        f"{best['ligand']} pada {best['receptor']}, tetapi ligand efficiency "
        f"terbaik {best_le['ligand_efficiency']:.3f} kcal/mol per heavy atom "
        f"untuk {best_le['ligand']} pada {best_le['receptor']} "
        "(reseptor surrogate PDB)."
    )


def _colab_section() -> str:
    return r"""## Docking Presisi Tinggi di Google Colab (CPU)

AutoDock Vina hanya memakai CPU; GPU tidak mempercepat docking. Jalankan
runtime **CPU** Colab agar CPU lokal tidak terbebani. Script
`scripts/run_docking_colab.sh` meng-install dependency, mengunduh Vina 1.2.7
Linux, menjalankan beberapa seed, mengagregasi hasil, dan mengarsipkan pose/log.
Notebook siap pakai tersedia di `notebooks/docking_colab.ipynb`.

```bash
# Di sel Colab:
!git clone https://github.com/kelvinchandradjayapranata-hue/In-Silico-Reengineering-of-Radiosynthetic-Pathways.git
%cd In-Silico-Reengineering-of-Radiosynthetic-Pathways
!bash scripts/run_docking_colab.sh
```

Konfigurasi opsional (set sebagai environment variable sebelum perintah):

- `SEEDS="42 43 44"` — seed Vina; kurangi bila waktu terbatas.
- `EXHAUSTIVENESS=32` — naikkan ke 64 untuk pencarian lebih luas.
- `NUM_MODES=20`, `BOX_SIZE="24 24 24"`.
- `LIGANDS` / `RECEPTORS` — batasi kombinasi yang dijalankan.
- `CPU="$(nproc)"` — jumlah vCPU Colab.
- `OUTPUT_DIR="/content/drive/MyDrive/melanin_docking_runs"` — salin hasil ke
  Google Drive; set `MOUNT_DRIVE=1` untuk mount Drive dari script.

Contoh dua seed dengan penyimpanan Drive:

```bash
!SEEDS="42 43" EXHAUSTIVENESS=32 MOUNT_DRIVE=1 \
  OUTPUT_DIR="/content/drive/MyDrive/melanin_docking_runs" \
  bash scripts/run_docking_colab.sh
```

Hasil kanonik ditulis ke `results/docking_results.json` dan
`results/docking_summary.tsv`. Hasil per-seed, log, dan arsip pose disimpan di
`results/colab_runs/<run_id>/` (di-ignore Git). Script agregasi juga dapat
dijalankan terpisah:

```bash
python scripts/aggregate_docking_seeds.py \
  --input-dir results/colab_runs/<run_id> \
  --out results/docking_results.json \
  --summary results/docking_summary.tsv
```
"""


def build_readme() -> str:
    """Bangun isi README berdasarkan manifest dan hasil FBA jika tersedia."""
    manifest = _read_json(RAW_DIR / "manifest.json")
    fba = _read_json(RESULTS_DIR / "metabolic_fba.json")
    chem = _read_json(RESULTS_DIR / "melanin_chem.json")
    blast = _read_json(RESULTS_DIR / "blast_validation.json")
    docking = _read_json(RESULTS_DIR / "docking_results.json")
    generated = datetime.now(timezone.utc).strftime("%Y-%m-%d")

    caveats = manifest.get("caveats", []) if manifest else []
    caveat_lines = "\n".join(f"- {item}" for item in caveats) if caveats else (
        "- Belum ada manifest data NCBI."
    )

    chem_status = (
        "tersedia"
        if chem and chem.get("monomers")
        else "belum dijalankan"
    )

    return f"""# In Silico Reengineering of Radiosynthetic Pathways

**AI-Driven Bio-Computational Framework for Cosmic Radiation Energy Conversion Using Fungal Melanin Templates**

> README dibuat otomatis oleh `src/generate_readme.py` pada {generated}.
> Proyek ini masih tahap eksplorasi komputasi; semua hasil adalah estimasi in silico.

## Abstrak

Proyek ini mengeksplorasi kerangka kerja komputasi untuk mempelajari hipotesis
interaksi antara melanin jamur, protein, dan jalur metabolik
(*Cladosporium sphaerospermum*). Fokus awal ada pada jalur DHN-melanin dan
DOPA-melanin, konversi energi radiasi menjadi ATP/NADPH, dan simulator
stoikiometri tereduksi.

## Status Saat Ini

- Pipeline inti: berjalan end-to-end.
- Simulasi FBA: {_status_line(fba)}
- Deskriptor monomer: {chem_status}.
- Data NCBI: tervalidasi terhadap organisme dan kata kunci; lihat bagian
  **Data dan Provenance**.
- Validasi homologi: {_blast_status_line(blast)}
- Docking pilot: {_docking_status_line(docking)}

## Arsitektur Pipeline

1. **`src/ncbi_fetcher.py`** — mengambil dan memvalidasi FASTA dari NCBI.
2. **`src/melanin_chem.py`** — menghitung deskriptor RDKit monomer DHN/L-DOPA.
3. **`src/metabolic_fba.py`** — model FBA tereduksi jalur DHN/DOPA + input radiasi.
4. **`src/visualization.py`** — grafik fluks ke `results/plots/fba_fluxes.png`.
5. **`src/blast_validation.py`** — tblastn ortholog terhadap assembly
   *C. sphaerospermum*.
6. **`src/melanin_docking.py`** — docking pilot DHN ke struktur homolog
   PKS1/laccase dengan AutoDock Vina.
7. **`src/generate_report.py`** — laporan ke `results/summary_report.md`.
8. **`src/run_pipeline.py`** — orkestrasi tahap inti dan tahap opsional.

## Quickstart (Windows PowerShell)

```powershell
python -m venv venv
.\\venv\\Scripts\\Activate.ps1
pip install -r requirements.txt

# Unduh data NCBI (butuh jaringan; set NCBI_EMAIL agar sesuai kebijakan NCBI)
$env:NCBI_EMAIL = "nama@institusi.ac.id"
python src/ncbi_fetcher.py --targets all

# Jalankan pipeline tanpa mengunduh ulang
python src/run_pipeline.py --skip-fetch

# Pipeline inti + sensitivitas + tahap opsional (bila tool tersedia)
python src/run_pipeline.py --skip-fetch --sensitivity --with-blast --with-docking
```

Perintah manual per tahap:

```powershell
python src/melanin_chem.py
python src/metabolic_fba.py --sensitivity
python src/visualization.py
python src/generate_report.py
python src/generate_readme.py

# Opsional: validasi homologi (butuh BLAST+ dan genome; lihat catatan)
python src/blast_validation.py

# Opsional: docking pilot (butuh Meeko dan binary AutoDock Vina)
python src/melanin_docking.py
```

Catatan tool opsional:

- BLAST+ diletakkan di `tools/blast/ncbi-blast-*/bin`, atau atur `BLAST_HOME`,
  atau gunakan `--blast-bin`.
- Binary AutoDock Vina diletakkan di `tools/vina/` (mis.
  `vina_1.2.7_win.exe`), atau gunakan `--vina`.
- Genome assembly diunduh ke `data/genomes/`; struktur PDB ke
  `data/structures/`. Keduanya di-ignore Git karena berukuran besar.

{_colab_section()}

## Data dan Provenance

{_provenance_table(manifest)}

### Assembly referensi

{_reference_assembly_lines(manifest)}

### Catatan ilmiah

{caveat_lines}

Detail checksum dan metadata ada di `data/raw/manifest.json`. Manifest
membedakan **urutan langsung** dari *C. sphaerospermum* dan **ortholog
referensi** dari spesies *Cladosporium* lain.

## Asumsi dan Keterbatasan Model FBA

- Model adalah stoikiometri tereduksi, bukan genome-scale model.
- Jalur DHN: 1 acetyl-CoA + 4 malonyl-CoA -> THN, lalu 2 NADPH per DHN.
- Jalur DOPA: L-tyrosine -> L-DOPA -> dopaquinone -> DOPA-melanin.
- Radiasi dimodelkan sebagai input elektron (-> NADH) dan energi (-> ATP),
  bukan model fotokimia/kuantum.
- Koenzim seperti CoA dan H+ dihilangkan; neraca massa/charge tidak lengkap.
- Reaksi biomassa dan ATP maintenance adalah placeholder, bukan parameter
  terukur untuk *C. sphaerospermum*.
- Hasil FBA adalah estimasi in silico, bukan bukti eksperimental.

## Asumsi dan Keterbatasan Validasi dan Docking

- tblastn menunjukkan homologi sekuens, bukan bukti fungsi enzimatik,
  ekspresi gen, atau keterlibatan jalur DHN-melanin.
- Tidak adanya hit tidak membuktikan gen absent (gap assembly, divergensi,
  atau query parsial).
- Docking memakai struktur homolog PDB (1GW0 untuk laccase, 8CG4 untuk PKS1),
  bukan struktur *C. cladosporioides*; skor Vina bukan afinitas target
  sebenarnya.
- Ligan DHN dimer/trimer adalah model sintetis dengan linkage sederhana,
  bukan struktur oligomer melanin alami yang unik.
- PKS1 mendock DHN sebagai produk, bukan substrat; relevansi biologisnya
  spekulatif.
- Skor Vina tidak dinormalisasi terhadap ukuran ligan; bandingkan ligand
  efficiency, bukan hanya afinitas mentah.
- Tidak ada validasi eksperimental atau simulasi MD untuk pose docking.
- Semua hasil komputasi memerlukan validasi eksperimental.

## Output

- `results/metabolic_fba.json` — hasil FBA lengkap.
- `results/fba_sensitivity.json` — sweep batas elektron/energi radiasi.
- `results/melanin_chem.json` — deskriptor monomer.
- `results/plots/fba_fluxes.png` — grafik fluks.
- `results/blast_validation.json` — ringkasan validasi homologi tblastn.
- `results/blast_hits.tsv` — semua hit tblastn (gabungan).
- `results/blast/<assembly>.tsv` — hit tblastn per assembly.
- `results/docking_results.json` — hasil lengkap docking pilot.
- `results/docking_summary.tsv` — ringkasan afinitas docking.
- `results/summary_report.md` — laporan ringkasan.
- `models/melanin_radiosynthesis_fba_v2.xml` — model SBML.

## Struktur Direktori

```text
melanin-biotech/
├── data/
│   ├── raw/                 # FASTA NCBI + manifest + rejected/
│   ├── processed/           # data hasil pemrosesan (opsional)
│   ├── genomes/             # genome assembly + BLAST DB (di-ignore)
│   ├── structures/          # PDB + reseptor siap docking (di-ignore)
│   └── ligands/             # SDF/PDBQT ligan (di-ignore)
├── models/                  # model SBML
├── results/
│   ├── plots/
│   ├── blast/               # hit tblastn per assembly
│   ├── docking/             # pose dan log Vina (di-ignore)
│   ├── metabolic_fba.json
│   ├── fba_sensitivity.json
│   ├── melanin_chem.json
│   ├── blast_validation.json
│   ├── blast_hits.tsv
│   ├── docking_results.json
│   ├── docking_summary.tsv
│   └── summary_report.md
├── src/
│   ├── ncbi_fetcher.py
│   ├── melanin_chem.py
│   ├── metabolic_fba.py
│   ├── visualization.py
│   ├── blast_validation.py
│   ├── melanin_docking.py
│   ├── generate_report.py
│   ├── generate_readme.py
│   └── run_pipeline.py
├── scripts/
│   ├── run_docking_colab.sh  # runner Vina multi-seed di Colab CPU
│   └── aggregate_docking_seeds.py
├── tools/                   # BLAST+ dan Vina (di-ignore)
├── docs/                    # catatan proyek dan roadmap
├── notebooks/               # docking_colab.ipynb + rencana eksplorasi
├── tests/                   # placeholder test suite (belum ada test otomatis)
├── .gitignore
├── requirements.txt
└── README.md
```

## Disclaimer Ilmiah

Semua angka dalam proyek ini adalah estimasi komputasi. Validasi eksperimental
tetap diperlukan sebelum klaim biologis atau aplikasi rekayasa apa pun.
"""


def generate_readme(output_path: Path | None = None) -> Path:
    destination = Path(output_path) if output_path else README_PATH
    destination.parent.mkdir(parents=True, exist_ok=True)
    temp_path = destination.with_suffix(destination.suffix + ".tmp")
    temp_path.write_text(build_readme(), encoding="utf-8")
    temp_path.replace(destination)
    return destination


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Generate README.md proyek.")
    parser.add_argument("--output", type=Path, default=README_PATH)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    path = generate_readme(args.output)
    print(f"[SUCCESS] README ditulis ke: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
