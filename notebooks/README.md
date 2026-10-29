# Notebook

- `docking_colab.ipynb` — runner AutoDock Vina multi-seed untuk Google Colab
  CPU. Buka di Colab, pilih runtime CPU, lalu jalankan sel dari atas ke bawah.
- Setelah pose tersedia, buat peta interaksi dengan
  `python src/analyze_interactions.py --ligand DHN_monomer` dan viewer 3D
  dengan `python src/visualize_pose.py --interactions results/interactions.json`.
- Klastering pose RMSD: `python src/cluster_poses.py --ligand DHN_monomer`.
- Preparasi dan runner MD OpenMM di Colab:
  `python -m pip install -r requirements.txt -r requirements-md.txt`, lalu
  `python src/prepare_openmm_md.py --ligand DHN_monomer`. Jalankan sistem yang
  dihasilkan dengan `bash results/md_prep/<reseptor>__<ligan>/run_colab.sh`.
- Analisis trajektori MD setelah produksi:
  `python src/analyze_md_trajectory.py --system-dir results/md_prep/<reseptor>__<ligan>`
  (menulis grafik ke `results/plots/md_trajectory_analysis.png`).
- Kalkulator laju Marcus dari parameter DFT/MD:
  `python src/marcus_electron_transfer.py --self-test`, lalu lihat
  `docs/marcus_dft_protocol.md` untuk protokol four-point λ, outer-sphere,
  dan caveat NAD+/NADP+.

Rencana notebook eksplorasi tambahan:

1. `01_ncbi_validation.ipynb` — validasi aksesi dan perbandingan ortholog.
2. `02_fba_exploration.ipynb` — eksplorasi flux distribution dan sensitivity.
3. `03_melanin_descriptors.ipynb` — perbandingan deskriptor monomer/polimer.

Untuk notebook lokal, gunakan kernel dari `venv` proyek.
