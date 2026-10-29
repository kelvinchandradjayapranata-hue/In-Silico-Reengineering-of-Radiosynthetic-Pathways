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
- `src/analyze_interactions.py` — peta H-bond, kontak hidrofobik, dan pi-pi
  dari pose terbaik DHN monomer.
- `src/visualize_pose.py` — viewer 3D interaktif pose ligan di kantung aktif.
- `src/cluster_poses.py` — klastering pose multi-seed berdasarkan RMSD di
  frame reseptor.
- `src/prepare_openmm_md.py` — preparasi sistem protein-ligan tersolvasi dan
  runner MD OpenMM untuk Colab.
- `src/analyze_md_trajectory.py` — RMSD backbone/ligan, RMSF per residu, dan
  occupancy H-bond dari trajektori MD.
- `src/marcus_electron_transfer.py` — kalkulator laju transfer elektron
  non-adiabatik Marcus dari parameter DFT/MD.
- `docs/marcus_dft_protocol.md` — protokol DFT, four-point reorganization,
  outer-sphere, coupling, dan caveat NAD+/NADP+.
- `requirements-md.txt` — dependency OpenMM/PDBFixer/OpenFF (Colab/Linux).

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

## Analisis interaksi dan visualisasi 3D

Setelah pose PDBQT tersedia (lokal atau hasil unduhan Colab), jalankan:

```bash
python src/analyze_interactions.py --ligand DHN_monomer
python src/visualize_pose.py --interactions results/interactions.json
```

`analyze_interactions.py` mendukung format JSON single-seed maupun agregasi
multi-seed. Script memilih pose terbaik yang dirujuk JSON, lalu menulis
`results/interactions.json`, `results/interactions.tsv`, dan
`results/interactions.md`. `visualize_pose.py` menulis HTML mandiri ke
`results/visualization/` (di-ignore Git) dengan `index.html` sebagai halaman
utama. Viewer memakai 3Dmol.js dari CDN, jadi perlu koneksi internet.

Kriteria default: H-bond jarak <=3.5 Å dan sudut D-H...A >=120°; kontak
hidrofobik karbon side chain <=4.0 Å; pi-pi centroid <=5.5 Å, ditambah kandidat
parallel-displaced/T-shaped offset dengan centroid <=6.5 Å dan jarak
atom-terdekat <=4.0 Å. Reseptor laccase 1GW0 tidak memiliki atom H, sehingga
sebagian kontak polar hanya dilaporkan sebagai kandidat `distance_only`.

## Klastering pose dan preparasi MD

Klastering pose memakai RMSD di frame reseptor tanpa superposisi (koreksi
automorfisme graph) dan complete linkage dengan cutoff default 2.0 Å:

```bash
python src/cluster_poses.py --ligand DHN_monomer
```

Output ditulis ke `results/pose_clusters.json`, `.tsv`, dan `.md`. Populasi
klaster adalah frekuensi mode Vina, bukan populasi Boltzmann.

Preparasi MD OpenMM (protein + ligan + TIP3P + ion) membutuhkan OpenMM,
PDBFixer, dan OpenFF, sehingga dijalankan di Colab/Linux:

```bash
python -m pip install -r requirements.txt -r requirements-md.txt
python src/prepare_openmm_md.py --ligand DHN_monomer --dry-run
python src/prepare_openmm_md.py --ligand DHN_monomer
bash results/md_prep/PKS1_surrogate_8CG4__DHN_monomer/run_colab.sh
```

Setiap sistem di `results/md_prep/<reseptor>__<ligan>/` berisi `system.xml`,
`system.pdb`, `state.xml`, `md_config.json`, `run_md.py`, dan `run_colab.sh`.
Runner menjalankan minimisasi, NVT/NPT singkat, lalu produksi 10 ns; semua
parameter dapat di-override dari CLI. `results/md_prep/` di-ignore Git.

Setelah produksi selesai, analisis trajektori dijalankan lokal:

```bash
python src/analyze_md_trajectory.py \
  --system-dir results/md_prep/PKS1_surrogate_8CG4__DHN_monomer
```

Output: `results/plots/md_trajectory_analysis.png`,
`results/md_trajectory_analysis.json`, dan
`results/md_trajectory_analysis.md`. RMSD dihitung setelah superposisi
backbone; RMSD ligan memakai frame reseptor; H-bond memakai kriteria
D...A <=3.5 Å dan sudut D-H...A >=120°. Gunakan `--stride` untuk mengurangi
memori pada trajektori panjang.

## Transfer elektron Marcus (DFT → laju)

`marcus_electron_transfer.py` hanya menghitung persamaan laju; ia tidak
menjalankan DFT. Parameter ΔG°, λ, dan H_AB harus berasal dari DFT/MD atau
data eksperimental. Modul ini stdlib-only dan dapat dijalankan lokal:

```bash
python src/marcus_electron_transfer.py --self-test
python src/marcus_electron_transfer.py \
  --delta-g-ev -0.18 --lambda-inner-ev 0.45 --lambda-outer-ev 0.85 \
  --coupling 1.0 --coupling-unit meV
```

Output default: `results/marcus_electron_transfer.json` dan
`results/marcus_electron_transfer.md`. Untuk four-point DFT, gunakan
`--four-point-json`; untuk parameter bola/kontinum, gunakan
`--donor-radius-a`, `--acceptor-radius-a`, `--distance-a`, `--eps-static`,
dan `--eps-optical`.

**Caveat utama:** NAD+/NADP+ adalah akseptor hidrida dua-elektron dalam
metabolisme normal, bukan akseptor elektron luar-sphere satu-elektron.
Oksidasi DHN juga pH-dependent dan dapat berlangsung melalui PCET. Model
Marcus satu-elektron hanya sah untuk jalur yang didefinisikan eksplisit dan
tidak otomatis mewakili hidrida fisiologis. Detail protokol dan validasi ada
di `docs/marcus_dft_protocol.md`.

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
- [x] Menambahkan peta interaksi dan viewer 3D pose ligan.
- [x] Menambahkan klastering pose RMSD dan preparasi sistem MD OpenMM.
- [x] Menambahkan analisis trajektori MD (RMSD, RMSF, occupancy H-bond).
- [x] Menambahkan kalkulator laju Marcus dan protokol DFT/Marcus.
- [ ] Menjalankan DFT produksi untuk four-point λ dan H_AB (Colab/HPC).
- [ ] Kalibrasi parameter biomass dan ATP maintenance.
- [ ] Menambahkan analisis sensitivitas parameter (bukan hanya batas radiasi).
- [ ] Menambahkan test suite dan CI.
- [ ] Menjalankan MD produksi dan validasi pose docking secara eksperimental
      (opsional; sistem dan runner sudah disiapkan).

## Catatan

Semua model masih bersifat in silico dan belum divalidasi eksperimental.
