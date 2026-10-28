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

- Pipeline dasar: berjalan end-to-end.
- Simulasi FBA: FBA `optimal`; DHN-melanin export 1.920, DOPA-melanin export 1.000.
- Deskriptor monomer: tersedia.
- Data NCBI: tervalidasi terhadap organisme dan kata kunci; lihat bagian
  **Data dan Provenance**.

## Arsitektur Pipeline

1. **`src/ncbi_fetcher.py`** — mengambil dan memvalidasi FASTA dari NCBI.
2. **`src/melanin_chem.py`** — menghitung deskriptor RDKit monomer DHN/L-DOPA.
3. **`src/metabolic_fba.py`** — model FBA tereduksi jalur DHN/DOPA + input radiasi.
4. **`src/visualization.py`** — grafik fluks ke `results/plots/fba_fluxes.png`.
5. **`src/generate_report.py`** — laporan ke `results/summary_report.md`.
6. **`src/run_pipeline.py`** — orkestrasi semua tahap.

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
```

Perintah manual per tahap:

```powershell
python src/melanin_chem.py
python src/metabolic_fba.py --sensitivity
python src/visualization.py
python src/generate_report.py
python src/generate_readme.py
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
- Modul docking molekuler belum diimplementasikan; `results/docking_results.json`
  belum tersedia.
- Hasil FBA adalah estimasi in silico, bukan bukti eksperimental.

## Output

- `results/metabolic_fba.json` — hasil FBA lengkap.
- `results/fba_sensitivity.json` — sweep batas elektron/energi radiasi.
- `results/melanin_chem.json` — deskriptor monomer.
- `results/plots/fba_fluxes.png` — grafik fluks.
- `results/summary_report.md` — laporan ringkasan.
- `models/melanin_radiosynthesis_fba_v2.xml` — model SBML.

## Struktur Direktori

```text
melanin-biotech/
├── data/
│   ├── raw/                 # FASTA NCBI + manifest + rejected/
│   └── processed/           # data hasil pemrosesan (opsional)
├── models/                  # model SBML
├── results/
│   ├── plots/
│   ├── metabolic_fba.json
│   ├── fba_sensitivity.json
│   ├── melanin_chem.json
│   └── summary_report.md
├── src/
│   ├── ncbi_fetcher.py
│   ├── melanin_chem.py
│   ├── metabolic_fba.py
│   ├── visualization.py
│   ├── generate_report.py
│   ├── generate_readme.py
│   └── run_pipeline.py
├── docs/                    # catatan proyek dan roadmap
├── notebooks/               # placeholder notebook eksplorasi
├── tests/                   # placeholder test suite (belum ada test otomatis)
├── .gitignore
├── requirements.txt
└── README.md
```

## Disclaimer Ilmiah

Semua angka dalam proyek ini adalah estimasi komputasi. Validasi eksperimental
tetap diperlukan sebelum klaim biologis atau aplikasi rekayasa apa pun.
