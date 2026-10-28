# Dokumentasi Proyek

Dokumentasi utama saat ini tersedia di:

- `README.md` — gambaran proyek, quickstart, dan status pipeline.
- `results/summary_report.md` — laporan hasil terbaru (auto-generated).
- `data/raw/README.md` — provenance data NCBI dan catatan validasi.
- `data/raw/manifest.json` — metadata, checksum, dan level bukti setiap sekuens.
- `models/melanin_radiosynthesis_fba_v2.xml` — model SBML hasil FBA.
- `results/blast_validation.json` — hasil validasi homologi tblastn terhadap
  tiga assembly *C. sphaerospermum*.
- `results/docking_results.json` — hasil docking pilot DHN monomer/dimer/trimer
  ke struktur homolog PKS1 dan laccase.
- `scripts/run_docking_colab.sh` — runner AutoDock Vina multi-seed untuk
  Google Colab CPU (lihat README).
- `scripts/aggregate_docking_seeds.py` — agregasi statistik antar-seed.

## Menjalankan docking di Google Colab

Vina hanya memakai CPU. Buka runtime **CPU** Colab, clone repositori, lalu
jalankan:

```bash
!bash scripts/run_docking_colab.sh
```

Atur `SEEDS`, `EXHAUSTIVENESS`, `NUM_MODES`, `BOX_SIZE`, `LIGANDS`,
`RECEPTORS`, `CPU`, dan `OUTPUT_DIR` sebagai environment variable bila perlu.
Script mengunduh Vina 1.2.7 Linux, menjalankan beberapa seed, dan menulis
`results/docking_results.json` serta `results/docking_summary.tsv`.

## Pertanyaan riset awal

1. Bagaimana stoikiometri jalur DHN-melanin dan DOPA-melanin dapat
   direpresentasikan dalam model FBA tereduksi?
2. Sejauh mana input elektron/energi radiasi dapat menggantikan atau
   melengkapi ATP/NADPH dari metabolisme glukosa?
3. Parameter mana yang paling membatasi produksi melanin pada model ini?

## Roadmap

- [x] Validasi ortholog *C. sphaerospermum* dengan tblastn terhadap assembly
      GCA_001931905.1, GCA_001931905.2, dan GCA_000261425.2.
- [x] Menambahkan modul docking pilot (DHN monomer/dimer/trimer vs struktur
      homolog PKS1/laccase).
- [x] Menyediakan runner Colab CPU untuk docking presisi tinggi multi-seed.
- [ ] Kalibrasi parameter biomass dan ATP maintenance.
- [ ] Menambahkan analisis sensitivitas parameter (bukan hanya batas radiasi).
- [ ] Menambahkan test suite dan CI.
- [ ] Validasi eksperimental atau simulasi MD untuk pose docking (opsional).

## Catatan

Semua model masih bersifat in silico dan belum divalidasi eksperimental.
