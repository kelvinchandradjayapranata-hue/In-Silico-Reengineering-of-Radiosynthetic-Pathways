# Peta Interaksi Ligan-Reseptor

> Dibuat: 2026-10-28T15:30:57.078095+00:00
> Mode docking: single_seed
>
> Analisis geometris dari pose docking. Bukan energi bebas ikat dan
> bukan validasi eksperimental.

## Ringkasan

| Reseptor | Ligan | Afinitas | LE | H-bond | Residu hidrofobik | Pi-pi |
|---|---|---|---|---|---|---|
| Laccase_surrogate_1GW0 | DHN_monomer | -5.10 | -0.425 | 0 | 3 | 1 |
| PKS1_surrogate_8CG4 | DHN_monomer | -5.91 | -0.493 | 1 | 4 | 1 |

## Laccase_surrogate_1GW0 — DHN_monomer

- Pose: `results/docking/poses/Laccase_surrogate_1GW0__DHN_monomer.pdbqt`
- Reseptor: `data/structures/prepared/Laccase_surrogate_1GW0_clean.pdb` (chain A)
- PDB: 1GW0; target: KAN3060494.1
- Afinitas Vina: -5.10 kcal/mol
- Ligand efficiency: -0.425 kcal/mol per heavy atom

### Ikatan hidrogen

Tidak ada H-bond yang lolos kriteria jarak/sudut.

### Kandidat polar / distance-only

| Residu | Atom protein | Atom ligan | Jarak (Å) | Sudut (°) | Catatan |
|---|---|---|---|---|---|
| GLU235 | OE2 | O10 | 2.762 | 52.5 | sudut di bawah ambang H-bond |
| GLY299 | N | O0 | 3.133 | n/a | pasangan donor-akseptor cocok, tetapi atom H tidak tersedia untuk menghitung sudut |
| ALA297 | O | O10 | 3.154 | 117.1 | sudut di bawah ambang H-bond |
| GLU235 | OE2 | O0 | 3.155 | 111.2 | sudut di bawah ambang H-bond |
| GLY299 | N | O10 | 3.175 | n/a | pasangan donor-akseptor cocok, tetapi atom H tidak tersedia untuk menghitung sudut |

### Kontak hidrofobik

| Residu | Jumlah kontak | Jarak minimum (Å) | Atom protein |
|---|---|---|---|
| PHE194 | 3 | 3.357 | CD1, CE1 |
| LEU429 | 4 | 3.389 | CD1, CD2 |
| PRO192 | 7 | 3.567 | CD, CG |

### Stacking pi-pi

| Ring protein | Ring ligan | Jarak centroid (Å) | Jarak atom min (Å) | Sudut (°) | Tipe | Geometri |
|---|---|---|---|---|---|---|
| PHE194:ring | ligand_ring_1 | 5.962 | 3.357 | 80.2 | T-shaped | offset |

## PKS1_surrogate_8CG4 — DHN_monomer

- Pose: `results/docking/poses/PKS1_surrogate_8CG4__DHN_monomer.pdbqt`
- Reseptor: `data/structures/prepared/PKS1_surrogate_8CG4_clean.pdb` (chain A)
- PDB: 8CG4; target: KAN3049821.1
- Afinitas Vina: -5.91 kcal/mol
- Ligand efficiency: -0.493 kcal/mol per heavy atom

### Ikatan hidrogen

| Residu | Atom protein | Atom ligan | Jarak (Å) | Sudut (°) | Arah |
|---|---|---|---|---|---|
| SER791 | O | O10 | 3.005 | 136.3 | ligand_donor |

### Kandidat polar / distance-only

| Residu | Atom protein | Atom ligan | Jarak (Å) | Sudut (°) | Catatan |
|---|---|---|---|---|---|
| THR647 | O | O0 | 2.836 | 24.6 | sudut di bawah ambang H-bond |
| THR647 | O | O10 | 2.989 | 25.2 | sudut di bawah ambang H-bond |

### Kontak hidrofobik

| Residu | Jumlah kontak | Jarak minimum (Å) | Atom protein |
|---|---|---|---|
| THR680 | 3 | 3.482 | CG2 |
| PRO649 | 4 | 3.647 | CB, CG |
| PHE790 | 3 | 3.668 | CB, CD2 |
| THR682 | 1 | 3.931 | CB |

### Stacking pi-pi

| Ring protein | Ring ligan | Jarak centroid (Å) | Jarak atom min (Å) | Sudut (°) | Tipe | Geometri |
|---|---|---|---|---|---|---|
| PHE790:ring | ligand_ring_2 | 6.141 | 3.747 | 20.0 | parallel | offset |

## Keterbatasan

- H-bond ditentukan secara geometris dari pose rigid; sudut hanya dihitung bila atom hidrogen tersedia. Protonasi reseptor dan ligan tidak dioptimalkan pada pH fisiologis.
- Kontak hidrofobik adalah jarak atom karbon, bukan estimasi energi.
- Klasifikasi pi-pi memakai jarak centroid (<=5.5 A), sudut normal ring, dan kandidat offset (centroid <=6.5 A, atom-terdekat <=4.0 A); ambang ini heuristik dan perlu diverifikasi dengan PLIP/Arpeggio.
- Pose terbaik Vina belum tentu pose yang paling stabil secara dinamika; diperlukan klaster pose dan simulasi MD.
- Reseptor adalah struktur homolog surrogate, bukan struktur target C. cladosporioides; hasil tidak boleh dilaporkan sebagai interaksi target sebenarnya.
