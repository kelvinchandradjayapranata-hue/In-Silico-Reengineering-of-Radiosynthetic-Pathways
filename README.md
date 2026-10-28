# In Silico Reengineering of Radiosynthetic Pathways

**AI-Driven Bio-Computational Framework for Cosmic Radiation Energy Conversion Using Fungal Melanin Templates**

> README dibuat otomatis oleh `src/generate_readme.py` pada 2026-10-28.
> Proyek ini masih tahap eksplorasi komputasi; semua hasil adalah estimasi in silico.

## Abstrak

Proyek ini mengeksplorasi kerangka kerja komputasi untuk mempelajari hipotesis
interaksi antara melanin jamur, protein, dan jalur metabolik
(*Cladosporium sphaerospermum*). Fokus awal ada pada jalur DHN-melanin dan
DOPA-melanin, konversi energi radiasi menjadi ATP/NADPH, dan simulator
stoikiometri tereduksi.

## Status Saat Ini

- Pipeline inti: berjalan end-to-end.
- Simulasi FBA: FBA `optimal`; DHN-melanin export 1.920, DOPA-melanin export 1.000.
- Deskriptor monomer: tersedia.
- Data NCBI: tervalidasi terhadap organisme dan kata kunci; lihat bagian
  **Data dan Provenance**.
- Validasi homologi: tblastn pada 3 assembly; 21 strong, 0 moderate, 0 weak, 0 tidak terdeteksi dari 21 kombinasi query x assembly (7 query ortholog).
- Docking pilot: 6/6 kombinasi berhasil; skor mentah terbaik -9.66 kcal/mol untuk DHN_trimer pada PKS1_surrogate_8CG4, tetapi ligand efficiency terbaik -0.493 kcal/mol per heavy atom untuk DHN_monomer pada PKS1_surrogate_8CG4 (reseptor surrogate PDB).

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
.\venv\Scripts\Activate.ps1
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

## Docking Presisi Tinggi di Google Colab (CPU)

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


## Data dan Provenance

| Target | Gen | Aksesi | Organisme | Level bukti |
|---|---|---|---|---|
| csph_pks_typeiii_gene | PKS tipe III (kandidat) | MK213995.1 | Cladosporium sphaerospermum | urutan langsung |
| csph_pks_typeiii_protein | PKS tipe III (kandidat) | QBZ73626.1 | Cladosporium sphaerospermum | urutan langsung |
| ccla_pks1 | PKS1 (non-reducing PKS) | KAN3049821.1 | Cladosporium cladosporioides | ortholog referensi |
| ccla_thn_reductase | THN reductase | KAN3049819.1 | Cladosporium cladosporioides | ortholog referensi |
| ccla_scytalone_dehydratase | Scytalone dehydratase | KAN3063391.1 | Cladosporium cladosporioides | ortholog referensi |
| ccla_laccase | Laccase | KAN3060494.1 | Cladosporium cladosporioides | ortholog referensi |
| ccla_tyrosinase | Tyrosinase | KAN3055251.1 | Cladosporium cladosporioides | ortholog referensi |
| chal_scytalone_dehydratase | Scytalone dehydratase | XP_069226016.1 | Cladosporium halotolerans | ortholog referensi |
| chal_thn_reductase | THN reductase | XP_069230222.1 | Cladosporium halotolerans | ortholog referensi |

### Assembly referensi

- `GCA_001931905.1` — Cladosporium sphaerospermum IMV 00045: Draft genome isolat dari ISS; file genomic.gbff tidak memiliki CDS teranotasi (hanya fitur source). Status NCBI: previous; digantikan oleh GCA_001931905.2.
- `GCA_001931905.2` — Cladosporium sphaerospermum IMV 00045: Assembly saat ini (current) untuk strain IMV 00045; menggantikan GCA_001931905.1.
- `GCA_000261425.2` — Cladosporium sphaerospermum UM 843: Assembly scaffold; direktori FTP NCBI tidak menyediakan protein FASTA maupun anotasi gen.

Genome assembly diunduh dan diverifikasi checksum-nya oleh `src/blast_validation.py` sesuai kebutuhan.

### Catatan ilmiah

- NCBI belum menyediakan set gen DHN-melanin teranotasi lengkap untuk Cladosporium sphaerospermum.
- Hanya MK213995.1 / QBZ73626.1 yang merupakan urutan langsung dari C. sphaerospermum; keduanya PKS tipe III parsial dan belum dikonfirmasi sebagai PKS jalur DHN-melanin.
- Target lain adalah ortholog referensi dari Cladosporium cladosporioides dan Cladosporium halotolerans, bukan urutan langsung C. sphaerospermum.
- Aksesi DQ000286 (Borrelia sp. 16S rRNA) yang sebelumnya tersimpan sebagai cladosporium_sample.fasta telah dipindahkan ke data/raw/rejected/.

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
