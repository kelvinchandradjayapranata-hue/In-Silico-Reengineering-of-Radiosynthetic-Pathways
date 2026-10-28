# Dokumentasi Proyek

Dokumentasi utama saat ini tersedia di:

- `README.md` — gambaran proyek, quickstart, dan status pipeline.
- `results/summary_report.md` — laporan hasil terbaru (auto-generated).
- `data/raw/README.md` — provenance data NCBI dan catatan validasi.
- `data/raw/manifest.json` — metadata, checksum, dan level bukti setiap sekuens.
- `models/melanin_radiosynthesis_fba_v2.xml` — model SBML hasil FBA.

## Pertanyaan riset awal

1. Bagaimana stoikiometri jalur DHN-melanin dan DOPA-melanin dapat
   direpresentasikan dalam model FBA tereduksi?
2. Sejauh mana input elektron/energi radiasi dapat menggantikan atau
   melengkapi ATP/NADPH dari metabolisme glukosa?
3. Parameter mana yang paling membatasi produksi melanin pada model ini?

## Roadmap

- [ ] Validasi ortholog *C. sphaerospermum* dengan BLAST/tblastn terhadap
      assembly GCA_001931905.1 atau GCA_000261425.2.
- [ ] Kalibrasi parameter biomass dan ATP maintenance.
- [ ] Menambahkan analisis sensitivitas parameter (bukan hanya batas radiasi).
- [ ] Menambahkan modul docking bila diperlukan.
- [ ] Menambahkan test suite dan CI.

## Catatan

Semua model masih bersifat in silico dan belum divalidasi eksperimental.
