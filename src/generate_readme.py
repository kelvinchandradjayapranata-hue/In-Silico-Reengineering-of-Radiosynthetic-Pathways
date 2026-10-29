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

Jika repositori masih **privat**, clone publik akan gagal. Pilihan yang aman:
jadikan repositori public, atau simpan token GitHub read-only di Colab Secrets
dan set `GITHUB_TOKEN` sebelum menjalankan script (script tidak menyimpan token
di `.git/config`).

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


def _interaction_section() -> str:
    return r"""## Analisis Interaksi dan Viewer 3D

Setelah pose PDBQT tersedia (lokal atau hasil unduhan dari Colab):

```bash
# Peta interaksi DHN monomer pada semua reseptor
python src/analyze_interactions.py --ligand DHN_monomer

# Viewer 3D interaktif (HTML mandiri; 3Dmol.js dari CDN)
python src/visualize_pose.py --interactions results/interactions.json
```

`analyze_interactions.py` mendukung JSON single-seed maupun agregasi
multi-seed. Script memilih pose terbaik yang dirujuk JSON, lalu menulis:

- `results/interactions.json` — detail semua kontak dan koordinat untuk viewer;
- `results/interactions.tsv` — tabel datar untuk analisis lanjutan;
- `results/interactions.md` — ringkasan per kompleks.

`visualize_pose.py` menulis satu HTML per kompleks dan `index.html` ke
`results/visualization/` (di-ignore Git). Viewer memuat 3Dmol.js dari CDN,
jadi perlu koneksi internet saat dibuka. Kontrol tersedia untuk rotasi,
spin, surface pocket, dan menyembunyikan garis kontak.

Kriteria default:

- H-bond: jarak donor-akseptor <=3.5 Å dan sudut D-H...A >=120°.
- Kontak hidrofobik: karbon side chain protein dengan karbon ligan <=4.0 Å.
- Pi-pi: centroid ring <=5.5 Å; kandidat parallel-displaced/T-shaped offset
  dengan centroid <=6.5 Å dan jarak atom-terdekat <=4.0 Å.

Reseptor laccase 1GW0 tidak memiliki atom H pada PDB bersih, sehingga sudut
H-bond tidak dapat dihitung; kontak polar yang cocok dilaporkan sebagai
kandidat `distance_only`. Untuk analisis produksi, tambahkan hidrogen dengan
PDBFixer/Reduce atau protonasi dengan PROPKA/H++ dan verifikasi dengan
PLIP/Arpeggio. Peta ini bersifat geometris, bukan energi ikat; validasi
memerlukan MD/MM-PBSA dan data eksperimental.
"""


def _clustering_section() -> str:
    return r"""## Klastering Pose (RMSD)

```bash
# Klaster DHN monomer pada semua reseptor (cutoff 2.0 Å)
python src/cluster_poses.py --ligand DHN_monomer

# Semua ligan, average linkage (lebih longgar; lihat caveat)
python src/cluster_poses.py --ligand all --linkage average
```

RMSD dihitung pada koordinat pose di **frame reseptor tanpa superposisi** dan
dikoreksi automorfisme graph RDKit. Superposisi ligan rigid seperti DHN akan
memberi RMSD ~0 untuk semua pose, sehingga mode `aligned` hanya disediakan
untuk diagnostik. Pengelompokan default memakai complete linkage agar setiap
pasangan dalam klaster berada <= cutoff.

Output:

- `results/pose_clusters.json` — klaster, populasi, afinitas, ligand
  efficiency, dan daftar anggota;
- `results/pose_clusters.tsv` — tabel ringkas per klaster;
- `results/pose_clusters.md` — ringkasan per kompleks.

Populasi klaster adalah frekuensi mode Vina yang dijalankan, bukan populasi
Boltzmann atau bukti kestabilan termodinamika.
"""


def _md_section() -> str:
    return r"""## Preparasi MD OpenMM (Colab/Linux)

Modul ini hanya menyiapkan sistem; simulasi produksi harus dijalankan di
Colab/Linux, bukan CPU lokal. Dependency berat dipisahkan di
`requirements-md.txt`:

```bash
python -m pip install -r requirements-md.txt
```

Preparasi sistem protein-ligan tersolvasi:

```bash
# Dry-run: validasi pose dan metadata tanpa membangun system
python src/prepare_openmm_md.py --ligand DHN_monomer --dry-run

# Bangun system untuk semua reseptor DHN monomer
python src/prepare_openmm_md.py --ligand DHN_monomer
```

Alur default: PDBFixer membuang heterogen non-air, menambahkan atom hilang dan
hidrogen pada pH 7.4, lalu `Modeller.addSolvent` menambahkan air TIP3P dan ion
Na+/Cl- pada 0.15 M. Protein memakai Amber ff14SB, ligan diparameterisasi
dengan OpenFF (`openff-2.2.1`, fallback ke versi sebelumnya) melalui
`openmmforcefields`; GAFF tersedia sebagai alternatif.

Setiap sistem ditulis ke `results/md_prep/<reseptor>__<ligan>/` (di-ignore Git)
berisi `system.xml`, `system.pdb`, `state.xml`, `md_config.json`, `run_md.py`,
`run_colab.sh`, dan `README.md`. Jalankan di Colab:

```bash
bash results/md_prep/PKS1_surrogate_8CG4__DHN_monomer/run_colab.sh
```

Runner menjalankan minimisasi, NVT 0.25 ns, NPT 0.25 ns, lalu produksi 10 ns
(timestep 2 fs; DCD, CSV, checkpoint). Override dengan `--production-ns`,
`--nvt-ns`, `--npt-ns`, `--platform CUDA|OpenCL|CPU`, atau `--minimize-only`.
Jalur OpenFF penuh perlu diverifikasi di Colab/Linux karena wheel OpenFF
tidak selalu tersedia untuk Windows atau versi Python lokal.
"""


def _md_trajectory_section() -> str:
    return r"""## Analisis Trajektori MD

Setelah `run_md.py` menyelesaikan tahap produksi, analisis trajektori dapat
dijalankan lokal dengan MDTraj (sudah ada di `requirements.txt`):

```bash
python src/analyze_md_trajectory.py \
  --system-dir results/md_prep/PKS1_surrogate_8CG4__DHN_monomer
```

Output default:

- `results/plots/md_trajectory_analysis.png` — empat panel: RMSD backbone,
  RMSD ligan di frame reseptor, RMSF per residu CA, dan occupancy H-bond;
- `results/md_trajectory_analysis.json` — nilai per frame, RMSF per residu,
  dan daftar H-bond;
- `results/md_trajectory_analysis.md` — ringkasan untuk audit.

Metrik dihitung setelah superposisi backbone ke frame pertama. RMSD ligan
memakai atom berat dan diukur di frame reseptor, sehingga mencerminkan drift
pose, bukan perubahan konformasi internal ligan. H-bond protein-ligan memakai
kriteria donor-akseptor <=3.5 Å dan sudut D-H...A >=120°; occupancy adalah
fraksi frame. Gunakan `--stride` untuk mengurangi beban memori, `--no-hbonds`
untuk melewati analisis H-bond, atau `--ligand-resname` bila nama residu ligan
berbeda. Analisis ini geometris dan bukan energi ikat.
"""


def _marcus_section() -> str:
    return r"""## Transfer Elektron Marcus (DFT → Laju)

`src/marcus_electron_transfer.py` mengubah parameter DFT/MD menjadi estimasi
laju transfer elektron non-adiabatik satu-elektron:

```text
k_ET = (2π/ħ) |H_AB|² / sqrt(4π λ k_B T)
       × exp[−(ΔG° + λ)² / (4 λ k_B T)]
```

Modul ini hanya menghitung persamaan laju (stdlib-only); ia tidak menjalankan
DFT. Parameter ΔG°, λ, dan H_AB harus berasal dari DFT/MD atau data
eksperimental.

```bash
# Self-test numerik
python src/marcus_electron_transfer.py --self-test

# Parameter langsung
python src/marcus_electron_transfer.py \
  --delta-g-ev -0.18 --lambda-inner-ev 0.45 --lambda-outer-ev 0.85 \
  --coupling 1.0 --coupling-unit meV

# Four-point λ dari DFT + outer-sphere kontinum
python src/marcus_electron_transfer.py \
  --four-point-json data/quantum/dhn_nad_four_point.json \
  --donor-radius-a 3.2 --acceptor-radius-a 3.4 --distance-a 6.5 \
  --eps-static 78.355 --eps-optical 1.776 \
  --delta-g-ev -0.18 --coupling 1.0 --coupling-unit meV
```

ΔG° juga dapat diturunkan dari potensial reduksi dengan `--donor-e0-v` dan
`--acceptor-e0-v` (referensi elektroda harus sama). Semua parameter dapat
diletakkan di satu JSON dan dijalankan dengan `--config`.

Output default: `results/marcus_electron_transfer.json` dan
`results/marcus_electron_transfer.md`. Protokol DFT lengkap (four-point,
SMD, GMH/CDFT, ensemble MD, validasi) ada di `docs/marcus_dft_protocol.md`.

**Caveat utama:** NAD+/NADP+ adalah akseptor hidrida dua-elektron dalam
metabolisme normal, bukan akseptor elektron luar-sphere satu-elektron.
Oksidasi DHN juga pH-dependent dan dapat berlangsung melalui PCET. Model
Marcus satu-elektron hanya sah untuk jalur yang didefinisikan eksplisit dan
tidak otomatis mewakili hidrida fisiologis. DFT produksi dijalankan di
Colab/HPC, bukan CPU lokal.
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
7. **`src/analyze_interactions.py`** — peta H-bond, kontak hidrofobik, dan
   pi-pi dari pose terbaik.
8. **`src/visualize_pose.py`** — viewer 3D interaktif pose ligan di kantung
   aktif (HTML + 3Dmol.js).
9. **`src/cluster_poses.py`** — klastering pose multi-seed berdasarkan RMSD
   di frame reseptor.
10. **`src/prepare_openmm_md.py`** — preparasi sistem protein-ligan tersolvasi
    untuk MD OpenMM eksplisit.
11. **`src/analyze_md_trajectory.py`** — RMSD, RMSF, dan occupancy H-bond
    dari trajektori MD.
12. **`src/marcus_electron_transfer.py`** — kalkulator laju transfer elektron
    non-adiabatik Marcus dari parameter DFT/MD.
13. **`src/generate_report.py`** — laporan ke `results/summary_report.md`.
14. **`src/run_pipeline.py`** — orkestrasi tahap inti dan tahap opsional.

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

# Peta interaksi dan viewer 3D dari pose yang tersedia
python src/analyze_interactions.py --ligand DHN_monomer
python src/visualize_pose.py --interactions results/interactions.json

# Klastering pose RMSD (tanpa superposisi, frame reseptor)
python src/cluster_poses.py --ligand DHN_monomer

# Preparasi MD OpenMM (butuh OpenMM/OpenFF; jalankan di Colab/Linux)
python src/prepare_openmm_md.py --ligand DHN_monomer

# Analisis trajektori MD setelah produksi selesai
python src/analyze_md_trajectory.py --system-dir results/md_prep/PKS1_surrogate_8CG4__DHN_monomer

# Estimasi laju transfer elektron Marcus dari parameter DFT/MD
python src/marcus_electron_transfer.py --self-test
python src/marcus_electron_transfer.py \\
  --delta-g-ev -0.18 --lambda-inner-ev 0.45 --lambda-outer-ev 0.85 \\
  --coupling 1.0 --coupling-unit meV
```

Catatan tool opsional:

- BLAST+ diletakkan di `tools/blast/ncbi-blast-*/bin`, atau atur `BLAST_HOME`,
  atau gunakan `--blast-bin`.
- Binary AutoDock Vina diletakkan di `tools/vina/` (mis.
  `vina_1.2.7_win.exe`), atau gunakan `--vina`.
- Genome assembly diunduh ke `data/genomes/`; struktur PDB ke
  `data/structures/`. Keduanya di-ignore Git karena berukuran besar.

{_colab_section()}

{_interaction_section()}

{_clustering_section()}

{_md_section()}

{_md_trajectory_section()}

{_marcus_section()}

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
- Peta interaksi memakai kriteria geometris (jarak/sudut), bukan energi
  ikat; hasilnya perlu diverifikasi dengan PLIP/Arpeggio dan MD/MM-PBSA.
- Laccase 1GW0 tidak memiliki hidrogen pada PDB bersih; sebagian H-bond hanya
  dapat dilaporkan sebagai kandidat distance-only.
- Belum ada simulasi MD produksi atau validasi eksperimental untuk pose
  docking; `prepare_openmm_md.py` baru menyiapkan sistem dan runner.
- Analisis trajektori MD memakai RMSD tanpa koreksi simetri dan occupancy
  H-bond geometris; keduanya bukan energi bebas ikat.
- Model Marcus satu-elektron tidak otomatis mewakili hidrida NAD+/NADP+
  fisiologis atau PCET oksidasi DHN; lihat `docs/marcus_dft_protocol.md`.
- Parameter DFT untuk four-point λ dan H_AB belum dijalankan; angka contoh
  di dokumentasi bukan hasil ilmiah.
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
- `results/interactions.json` — peta H-bond/hidrofobik/pi-pi.
- `results/interactions.tsv` — tabel interaksi untuk analisis lanjutan.
- `results/interactions.md` — ringkasan interaksi per kompleks.
- `results/pose_clusters.json` — klaster pose dan statistik per klaster.
- `results/pose_clusters.tsv` — tabel ringkas klaster pose.
- `results/pose_clusters.md` — ringkasan klaster per kompleks.
- `results/md_prep/` — sistem OpenMM tersolvasi dan runner MD (di-ignore Git).
- `results/plots/md_trajectory_analysis.png` — grafik RMSD/RMSF/H-bond MD.
- `results/md_trajectory_analysis.json` — data analisis trajektori.
- `results/md_trajectory_analysis.md` — ringkasan analisis trajektori.
- `results/marcus_electron_transfer.json` — parameter dan laju Marcus.
- `results/marcus_electron_transfer.md` — ringkasan audit estimasi Marcus.
- `results/visualization/index.html` — viewer 3D (di-ignore Git).
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
│   ├── interactions.json
│   ├── interactions.tsv
│   ├── interactions.md
│   ├── pose_clusters.json
│   ├── pose_clusters.tsv
│   ├── pose_clusters.md
│   ├── md_trajectory_analysis.json
│   ├── md_trajectory_analysis.md
│   ├── marcus_electron_transfer.json
│   ├── marcus_electron_transfer.md
│   ├── md_prep/              # sistem OpenMM + runner MD (di-ignore)
│   ├── visualization/       # HTML viewer 3D (di-ignore)
│   └── summary_report.md
├── src/
│   ├── ncbi_fetcher.py
│   ├── melanin_chem.py
│   ├── metabolic_fba.py
│   ├── visualization.py
│   ├── blast_validation.py
│   ├── melanin_docking.py
│   ├── analyze_interactions.py
│   ├── visualize_pose.py
│   ├── cluster_poses.py
│   ├── prepare_openmm_md.py
│   ├── analyze_md_trajectory.py
│   ├── marcus_electron_transfer.py
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
├── requirements-md.txt
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
