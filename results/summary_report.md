# Ringkasan Hasil Riset Radiosintesis Melanin

> Dibuat: 2026-10-29 14:18 UTC
>
> Ringkasan ini mengompilasi artefak komputasi yang tersedia. Skor model
> merupakan estimasi theoretical/in silico dan bukan pengganti validasi
> eksperimental.

## 1. Simulasi FBA Radiosintesis

- Model: `melanin_radiosynthesis_fba_v2` (versi 2.0)
- Status solver: `optimal`
- Objective total ekspor melanin: **2.920**
- DHN-melanin export: **1.920**
- DOPA-melanin export: **1.000**
- Glucose terpakai: 5.000
- Oksigen terpakai: 3.950
- Tirosin terpakai: 1.000
- Elektron radiasi terpakai: 0.000
- Energi radiasi terpakai: 2.420
- Kontribusi ATP radiasi terhadap demand ATP: **13.3%**

### Fluks kunci

| Reaksi | Fluks |
|---|---|
| EX_GLC | 5.000 |
| GLYCOLYSIS | 4.900 |
| ACCARBOXYLASE | 7.680 |
| PKS_THN | 1.920 |
| MELANIN_DHN_EXPORT | 1.920 |
| TYROSINASE_CRESOLASE | 1.000 |
| MELANIN_DOPA_EXPORT | 1.000 |
| RESPIRATION | 2.980 |
| TRANS_HYDROGENASE | 3.840 |
| RAD_ELECTRON_IN | 0.000 |
| RAD_ENERGY_IN | 2.420 |
| ATPM | 10.000 |
| BIOMASS | 0.100 |

### Batasan yang mengikat

- Glucose uptake mencapai batas atas (5)
- Tyrosine uptake mencapai batas atas (1)

### Peringatan model

- Input elektron radiasi tidak terpakai; NADH dari glikolisis dan PPP cukup untuk kebutuhan redoks.

### Asumsi model

- Model adalah stoikiometri tereduksi, bukan genome-scale model.
- 1 acetyl-CoA + 4 malonyl-CoA menghasilkan 1 THN (PKS non-reducing, lumped).
- Jalur DHN: 2 NADPH per DHN (THN reductase dan 1,3,8-THN reductase).
- 1 ATP dan 1 CO2 dikonsumsi per malonyl-CoA (acetyl-CoA carboxylase).
- 1 glucose -> 2 acetyl-CoA + 2 NADH + 2 ATP + 2 CO2 (glycolysis lumped).
- 1 glucose -> 1 pentose phosphate + 2 NADPH + 1 CO2 (PPP oksidatif lumped).
- Respirasi lumped: 2 NADH + 0.5 O2 -> 2 NAD+ + 2 ATP (P/O = 1).
- Radiasi dimodelkan sebagai input elektron (-> NADH) dan energi (-> ATP), bukan model fotokimia/kuantum.
- Koenzim seperti CoA dan H+ dihilangkan; neraca massa/charge tidak lengkap.
- Reaksi biomassa adalah placeholder, bukan komposisi biomassa terukur.
- Jalur DOPA adalah reaksi tyrosinase + polimerisasi yang dilumpukan.
- Semua hasil adalah estimasi in silico, bukan validasi eksperimental.

## 2. Deskriptor Monomer Melanin

| Monomer | SMILES kanonik | Formula | MW (g/mol) | logP | TPSA (Å²) | HBD | HBA | Gasteiger min | Gasteiger max |
|---|---|---|---|---|---|---|---|---|---|
| 1,8-dihydroxynaphthalene (DHN) | Oc1cccc2cccc(O)c12 | C10H8O2 | 160.170 | 2.250 | 40.460 | 2 | 2 | -0.507 | 0.293 |
| L-DOPA | N[C@@H](Cc1ccc(O)c(O)c1)C(=O)O | C9H11NO4 | 197.190 | 0.050 | 103.780 | 4 | 4 | -0.504 | 0.320 |

Metode: RDKit 2D descriptors; Gasteiger partial charges after explicit hydrogen addition. SMILES diverifikasi terhadap PubChem.

### Asumsi

- DHN direpresentasikan sebagai monomer 1,8-dihydroxynaphthalene, bukan polimer.
- DOPA direpresentasikan sebagai L-DOPA, prekursor DOPA-melanin.
- Muatan Gasteiger dihitung setelah penambahan hidrogen eksplisit.
- Deskriptor 2D tidak menangkap polimerisasi, agregasi, atau konformasi 3D.
- Semua nilai adalah estimasi in silico, bukan hasil pengukuran eksperimental.

## 3. Data dan Provenance

| Target | Gen | Aksesi | Organisme | Level bukti | Panjang | Status |
|---|---|---|---|---|---|---|
| csph_pks_typeiii_gene | PKS tipe III (kandidat) | MK213995.1 | Cladosporium sphaerospermum | Urutan langsung | 705 | cached |
| csph_pks_typeiii_protein | PKS tipe III (kandidat) | QBZ73626.1 | Cladosporium sphaerospermum | Urutan langsung | 235 | cached |
| ccla_pks1 | PKS1 (non-reducing PKS) | KAN3049821.1 | Cladosporium cladosporioides | Ortholog referensi | 2180 | cached |
| ccla_thn_reductase | THN reductase | KAN3049819.1 | Cladosporium cladosporioides | Ortholog referensi | 264 | cached |
| ccla_scytalone_dehydratase | Scytalone dehydratase | KAN3063391.1 | Cladosporium cladosporioides | Ortholog referensi | 184 | cached |
| ccla_laccase | Laccase | KAN3060494.1 | Cladosporium cladosporioides | Ortholog referensi | 588 | cached |
| ccla_tyrosinase | Tyrosinase | KAN3055251.1 | Cladosporium cladosporioides | Ortholog referensi | 593 | cached |
| chal_scytalone_dehydratase | Scytalone dehydratase | XP_069226016.1 | Cladosporium halotolerans | Ortholog referensi | 183 | cached |
| chal_thn_reductase | THN reductase | XP_069230222.1 | Cladosporium halotolerans | Ortholog referensi | 265 | cached |

### Catatan ilmiah

- NCBI belum menyediakan set gen DHN-melanin teranotasi lengkap untuk Cladosporium sphaerospermum.
- Hanya MK213995.1 / QBZ73626.1 yang merupakan urutan langsung dari C. sphaerospermum; keduanya PKS tipe III parsial dan belum dikonfirmasi sebagai PKS jalur DHN-melanin.
- Target lain adalah ortholog referensi dari Cladosporium cladosporioides dan Cladosporium halotolerans, bukan urutan langsung C. sphaerospermum.
- Aksesi DQ000286 (Borrelia sp. 16S rRNA) yang sebelumnya tersimpan sebagai cladosporium_sample.fasta telah dipindahkan ke data/raw/rejected/.

### Assembly referensi

- `GCA_001931905.1` — Cladosporium sphaerospermum IMV 00045: Draft genome isolat dari ISS; file genomic.gbff tidak memiliki CDS teranotasi (hanya fitur source). Status NCBI: previous; digantikan oleh GCA_001931905.2.
- `GCA_001931905.2` — Cladosporium sphaerospermum IMV 00045: Assembly saat ini (current) untuk strain IMV 00045; menggantikan GCA_001931905.1.
- `GCA_000261425.2` — Cladosporium sphaerospermum UM 843: Assembly scaffold; direktori FTP NCBI tidak menyediakan protein FASTA maupun anotasi gen.

Genome assembly diunduh dan diverifikasi checksum-nya oleh `src/blast_validation.py` sesuai kebutuhan.

## 4. Analisis Sensitivitas Input Radiasi

| e- limit | Energi limit | Status | Objective | DHN | DOPA | e- terpakai | Energi terpakai |
|---|---|---|---|---|---|---|---|
| 0.000 | 0.000 | error | - | - | - | - | - |
| 0.000 | 2.500 | optimal | 2.920 | 1.920 | 1.000 | 0.000 | 2.420 |
| 0.000 | 5.000 | optimal | 2.920 | 1.920 | 1.000 | 0.000 | 2.420 |
| 0.000 | 7.500 | optimal | 2.920 | 1.920 | 1.000 | 0.000 | 2.420 |
| 0.000 | 10.000 | optimal | 2.920 | 1.920 | 1.000 | 0.000 | 2.420 |
| 2.500 | 0.000 | optimal | 2.920 | 1.920 | 1.000 | 2.420 | 0.000 |
| 2.500 | 2.500 | optimal | 2.920 | 1.920 | 1.000 | 0.000 | 2.420 |
| 2.500 | 5.000 | optimal | 2.920 | 1.920 | 1.000 | 0.000 | 2.420 |
| 2.500 | 7.500 | optimal | 2.920 | 1.920 | 1.000 | 0.000 | 2.420 |
| 2.500 | 10.000 | optimal | 2.920 | 1.920 | 1.000 | 0.000 | 2.420 |
| 5.000 | 0.000 | optimal | 2.920 | 1.920 | 1.000 | 2.420 | 0.000 |
| 5.000 | 2.500 | optimal | 2.920 | 1.920 | 1.000 | 0.000 | 2.420 |
| 5.000 | 5.000 | optimal | 2.920 | 1.920 | 1.000 | 0.000 | 2.420 |
| 5.000 | 7.500 | optimal | 2.920 | 1.920 | 1.000 | 0.000 | 2.420 |
| 5.000 | 10.000 | optimal | 2.920 | 1.920 | 1.000 | 0.000 | 2.420 |
| 7.500 | 0.000 | optimal | 2.920 | 1.920 | 1.000 | 2.420 | 0.000 |
| 7.500 | 2.500 | optimal | 2.920 | 1.920 | 1.000 | 0.000 | 2.420 |
| 7.500 | 5.000 | optimal | 2.920 | 1.920 | 1.000 | 0.000 | 2.420 |
| 7.500 | 7.500 | optimal | 2.920 | 1.920 | 1.000 | 0.000 | 2.420 |
| 7.500 | 10.000 | optimal | 2.920 | 1.920 | 1.000 | 0.000 | 2.420 |
| 10.000 | 0.000 | optimal | 2.920 | 1.920 | 1.000 | 2.420 | 0.000 |
| 10.000 | 2.500 | optimal | 2.920 | 1.920 | 1.000 | 0.000 | 2.420 |
| 10.000 | 5.000 | optimal | 2.920 | 1.920 | 1.000 | 0.000 | 2.420 |
| 10.000 | 7.500 | optimal | 2.920 | 1.920 | 1.000 | 0.000 | 2.420 |
| 10.000 | 10.000 | optimal | 2.920 | 1.920 | 1.000 | 0.000 | 2.420 |

Catatan: status `error` (biasanya pada e- = 0 dan energi = 0) berarti kebutuhan ATP maintenance tidak dapat dipenuhi tanpa input radiasi.

Sweep ini menunjukkan kondisi ketika input elektron/energi radiasi menjadi faktor pembatas. Nilai adalah estimasi in silico.

## 5. Validasi Homologi Ortholog (tblastn)

- Metode: tblastn protein query vs six-frame translation of genomic assembly nucleotide BLAST database
- Versi BLAST: `tblastn: 2.17.0+`
- Parameter: e-value <= 1e-05, max_target_seqs 25, max_hsps 10, SEG yes
- Total hit mentah: 543

### Ringkasan klasifikasi

| Assembly | Strong | Moderate | Weak | Tidak terdeteksi |
|---|---|---|---|---|
| GCA_001931905.1 | 7 | 0 | 0 | 0 |
| GCA_001931905.2 | 7 | 0 | 0 | 0 |
| GCA_000261425.2 | 7 | 0 | 0 | 0 |

### Hit terbaik per query

| Query | Gen | GCA_001931905.1 | GCA_001931905.2 | GCA_000261425.2 |
|---|---|---|---|---|
| ccla_pks1\|KAN3049821.1 | PKS1 (non-reducing PKS) | strong (92.1%, cov 0.74) | strong (92.1%, cov 0.74) | strong (91.6%, cov 0.74) |
| ccla_thn_reductase\|KAN3049819.1 | THN reductase | strong (93.9%, cov 0.99) | strong (93.9%, cov 0.99) | strong (92.0%, cov 0.99) |
| ccla_scytalone_dehydratase\|KAN3063391.1 | Scytalone dehydratase | strong (89.2%, cov 0.91) | strong (89.2%, cov 0.91) | strong (86.9%, cov 0.91) |
| ccla_laccase\|KAN3060494.1 | Laccase | strong (74.6%, cov 1.00) | strong (74.6%, cov 1.00) | strong (83.9%, cov 0.80) |
| ccla_tyrosinase\|KAN3055251.1 | Tyrosinase | strong (75.4%, cov 0.99) | strong (75.4%, cov 0.99) | strong (72.4%, cov 0.99) |
| chal_scytalone_dehydratase\|XP_069226016.1 | Scytalone dehydratase | strong (95.2%, cov 0.91) | strong (95.2%, cov 0.91) | strong (100.0%, cov 0.91) |
| chal_thn_reductase\|XP_069230222.1 | THN reductase | strong (94.3%, cov 1.00) | strong (94.3%, cov 1.00) | strong (100.0%, cov 1.00) |

Detail subject, koordinat, dan e-value tersedia di `results/blast_validation.json` dan `results/blast_hits.tsv`.

Dari 21 kombinasi query x assembly, 21 menunjukkan hit strong, 0 moderate, 0 weak, dan 0 tidak terdeteksi. Query dengan bukti strong/moderate di setidaknya satu assembly: ccla_laccase|KAN3060494.1, ccla_pks1|KAN3049821.1, ccla_scytalone_dehydratase|KAN3063391.1, ccla_thn_reductase|KAN3049819.1, ccla_tyrosinase|KAN3055251.1, chal_scytalone_dehydratase|XP_069226016.1, chal_thn_reductase|XP_069230222.1. Status ini adalah bukti homologi sekuens, bukan konfirmasi fungsi enzimatik atau ekspresi gen.

### Keterbatasan validasi BLAST

- Tidak adanya hit dapat disebabkan gap assembly, divergensi sekuens, atau query parsial; bukan bukti gen absent.
- Koordinat subject bersifat spesifik assembly dan tidak dinormalisasi antar assembly.
- Hasil ini tidak membuktikan ekspresi, aktivitas enzim, atau keterlibatan jalur DHN-melanin.
- GCA_001931905.1 berstatus previous dan digantikan oleh GCA_001931905.2; keduanya disertakan untuk transparansi.

## 6. Docking Pilot DHN-melanin (AutoDock Vina)

- Status: `completed`
- Engine: AutoDock Vina v1.2.7; Meeko 0.8.0; RDKit 2026.9.1

> **Catatan:** Reseptor default adalah struktur homolog PDB (1GW0 untuk laccase, 8CG4 untuk PKS1), bukan struktur C. cladosporioides. Skor Vina tidak boleh ditafsirkan sebagai energi ikat target sebenarnya.

### Afinitas terbaik (kcal/mol)

| Reseptor (surrogate) | DHN_monomer | DHN_dimer | DHN_trimer |
|---|---|---|---|
| Laccase_surrogate_1GW0 | -5.10 | -3.16 | 2.27 |
| PKS1_surrogate_8CG4 | -5.91 | -8.89 | -9.66 |

### Ligand efficiency (kcal/mol per heavy atom)

| Reseptor (surrogate) | DHN_monomer | DHN_dimer | DHN_trimer |
|---|---|---|---|
| Laccase_surrogate_1GW0 | -0.425 | -0.132 | 0.063 |
| PKS1_surrogate_8CG4 | -0.493 | -0.370 | -0.268 |

Ligand efficiency = afinitas terbaik / jumlah atom berat. Metrik ini sebagian mengoreksi bias ukuran ligan pada skor Vina.

### Reseptor

| Reseptor | PDB | Target | Identitas / coverage | Sumber pusat kotak |
|---|---|---|---|---|
| Laccase_surrogate_1GW0 | 1GW0 | KAN3060494.1 | 44.7% / 91.7% | T1 copper serial 9143 (Cys503 SG, 2.2 A) |
| PKS1_surrogate_8CG4 | 8CG4 | KAN3049821.1 | 36.4% / 94.7% | SG atom Cys543 (katalitik KS motif DTACSSSL) |

6 dari 6 run docking berhasil. Skor mentah terbaik -9.663 kcal/mol untuk DHN_trimer pada PKS1_surrogate_8CG4; ligand efficiency terbaik -0.493 kcal/mol per heavy atom untuk DHN_monomer pada PKS1_surrogate_8CG4. Skor Vina tidak dinormalisasi terhadap ukuran ligan, sehingga monomer, dimer, dan trimer harus dibandingkan dengan ligand efficiency, bukan skor mentah saja. Semua reseptor adalah surrogate homolog, sehingga angka ini bukan afinitas target C. cladosporioides dan hanya untuk perbandingan internal.

### Keterbatasan docking

- Skor Vina bukan energi bebas ikat eksperimental (bukan ΔG).
- Tidak ada validasi eksperimental atau simulasi MD untuk pose terbaik.
- PKS1 mendock DHN sebagai produk, bukan substrat; relevansi biologisnya spekulatif.
- Model dimer/trimer tidak mewakili heterogenitas melanin alami.
- Skor Vina tidak dinormalisasi terhadap ukuran ligan; ligan lebih besar cenderung memperoleh skor lebih favorable. Bandingkan ligand efficiency, bukan hanya afinitas mentah.
- Perbedaan kecil skor (<1 kcal/mol) tidak dapat dianggap signifikan.

## Artefak yang belum tersedia / opsional

- Semua artefak inti tersedia.
- `data/processed/*_stats.json` — belum ada data terproses.
