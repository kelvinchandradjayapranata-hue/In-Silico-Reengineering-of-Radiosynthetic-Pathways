# Ringkasan Hasil Riset Radiosintesis Melanin

> Dibuat: 2026-10-28 11:55 UTC
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

### Assembly referensi (tidak diunduh)

- `GCA_001931905.1` — Cladosporium sphaerospermum IMV 00045: Draft genome isolat dari ISS; file genomic.gbff tidak memiliki CDS teranotasi (hanya fitur source).
- `GCA_000261425.2` — Cladosporium sphaerospermum UM 843: Assembly scaffold; direktori FTP NCBI tidak menyediakan protein FASTA maupun anotasi gen.

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

## Artefak yang belum tersedia / opsional

- Semua artefak inti tersedia.
- `results/docking_results.json` — opsional; modul docking belum diimplementasikan.
- `data/processed/*_stats.json` — opsional; belum ada data terproses.
