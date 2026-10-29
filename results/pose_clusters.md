# Klastering Pose Docking (RMSD)

> Dibuat: 2026-10-29T12:52:51.287910+00:00
> Mode docking: single_seed
> Cutoff RMSD: 2.00 Å; linkage: complete
>
> RMSD dihitung di frame reseptor tanpa superposisi dan dikoreksi simetri. Populasi klaster bukan populasi Boltzmann.

## Ringkasan

| Reseptor | Ligan | Pose | Jumlah seed | Klaster | Afinitas terbaik | LE terbaik | Populasi terbesar |
|---|---|---|---|---|---|---|---|
| Laccase_surrogate_1GW0 | DHN_monomer | 9 | 1 | 5 | -5.102 | -0.425 | 3 |
| PKS1_surrogate_8CG4 | DHN_monomer | 9 | 1 | 6 | -5.911 | -0.493 | 2 |

## Laccase_surrogate_1GW0 — DHN_monomer

- Pose: 9 mode dari 1 seed
- Atom berat: 12
- Klaster: 5

| Klaster | Populasi | Fraksi | Afinitas terbaik | Rata-rata | SD | LE terbaik | Representatif | Max RMSD intra |
|---|---|---|---|---|---|---|---|---|
| 1 | 3 | 33.3% | -5.102 | -5.062 | 0.058 | -0.425 | 42/1 | 1.587 |
| 2 | 2 | 22.2% | -4.603 | -4.596 | 0.009 | -0.384 | 42/4 | 0.364 |
| 3 | 2 | 22.2% | -4.463 | -4.460 | 0.004 | -0.372 | 42/6 | 0.015 |
| 4 | 1 | 11.1% | -4.328 | -4.328 | 0.000 | -0.361 | 42/8 | 0.000 |
| 5 | 1 | 11.1% | -4.130 | -4.130 | 0.000 | -0.344 | 42/9 | 0.000 |

### Anggota klaster 1 (populasi terbesar)

| Pose | Seed | Model | Afinitas | LE | RMSD ke rep. |
|---|---|---|---|---|---|
| seed42/model1 | 42 | 1 | -5.102 | -0.425 | 0.000 |
| seed42/model2 | 42 | 2 | -5.089 | -0.424 | 0.051 |
| seed42/model3 | 42 | 3 | -4.995 | -0.416 | 1.587 |

## PKS1_surrogate_8CG4 — DHN_monomer

- Pose: 9 mode dari 1 seed
- Atom berat: 12
- Klaster: 6

| Klaster | Populasi | Fraksi | Afinitas terbaik | Rata-rata | SD | LE terbaik | Representatif | Max RMSD intra |
|---|---|---|---|---|---|---|---|---|
| 1 | 2 | 22.2% | -5.911 | -5.909 | 0.004 | -0.493 | 42/1 | 0.048 |
| 2 | 2 | 22.2% | -5.873 | -5.870 | 0.004 | -0.489 | 42/3 | 0.080 |
| 3 | 2 | 22.2% | -5.722 | -5.717 | 0.007 | -0.477 | 42/5 | 0.336 |
| 4 | 1 | 11.1% | -5.700 | -5.700 | 0.000 | -0.475 | 42/7 | 0.000 |
| 5 | 1 | 11.1% | -5.605 | -5.605 | 0.000 | -0.467 | 42/8 | 0.000 |
| 6 | 1 | 11.1% | -5.552 | -5.552 | 0.000 | -0.463 | 42/9 | 0.000 |

### Anggota klaster 1 (populasi terbesar)

| Pose | Seed | Model | Afinitas | LE | RMSD ke rep. |
|---|---|---|---|---|---|
| seed42/model1 | 42 | 1 | -5.911 | -0.493 | 0.000 |
| seed42/model2 | 42 | 2 | -5.906 | -0.492 | 0.048 |

## Keterbatasan

- RMSD default dihitung pada koordinat pose di frame reseptor tanpa superposisi, dengan koreksi automorfisme graph. Ini memisahkan orientasi/translasi pose di kantung; superposisi ligan rigid akan memberi RMSD ~0 untuk semua pose.
- Complete linkage menjamin setiap pasangan dalam klaster berada <= cutoff; average linkage tidak menjamin hal itu.
- Populasi klaster adalah frekuensi pose di antara mode Vina yang dijalankan, bukan populasi Boltzmann atau bukti kestabilan termodinamika.
- Nilai Vina bukan energi bebas ikat; klaster terbaik belum tentu pose yang paling stabil dalam MD.
- Reseptor adalah struktur homolog surrogate, bukan struktur target C. cladosporioides.
