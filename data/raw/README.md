# Data mentah NCBI

Direktori ini berisi sekuens referensi yang diunduh dan divalidasi oleh
`src/ncbi_fetcher.py`. Lihat `manifest.json` untuk checksum dan metadata lengkap.

## Target

| Target | Gen | Aksesi | Organisme | Level bukti | Status | Panjang |
|---|---|---|---|---|---|---|
| csph_pks_typeiii_gene | PKS tipe III (kandidat) | MK213995.1 | Cladosporium sphaerospermum | urutan langsung | cached | 705 |
| csph_pks_typeiii_protein | PKS tipe III (kandidat) | QBZ73626.1 | Cladosporium sphaerospermum | urutan langsung | cached | 235 |
| ccla_pks1 | PKS1 (non-reducing PKS) | KAN3049821.1 | Cladosporium cladosporioides | ortholog referensi | cached | 2,180 |
| ccla_thn_reductase | THN reductase | KAN3049819.1 | Cladosporium cladosporioides | ortholog referensi | cached | 264 |
| ccla_scytalone_dehydratase | Scytalone dehydratase | KAN3063391.1 | Cladosporium cladosporioides | ortholog referensi | cached | 184 |
| ccla_laccase | Laccase | KAN3060494.1 | Cladosporium cladosporioides | ortholog referensi | cached | 588 |
| ccla_tyrosinase | Tyrosinase | KAN3055251.1 | Cladosporium cladosporioides | ortholog referensi | cached | 593 |
| chal_scytalone_dehydratase | Scytalone dehydratase | XP_069226016.1 | Cladosporium halotolerans | ortholog referensi | cached | 183 |
| chal_thn_reductase | THN reductase | XP_069230222.1 | Cladosporium halotolerans | ortholog referensi | cached | 265 |

## Catatan ilmiah

- NCBI belum menyediakan set gen DHN-melanin teranotasi lengkap untuk Cladosporium sphaerospermum.
- Hanya MK213995.1 / QBZ73626.1 yang merupakan urutan langsung dari C. sphaerospermum; keduanya PKS tipe III parsial dan belum dikonfirmasi sebagai PKS jalur DHN-melanin.
- Target lain adalah ortholog referensi dari Cladosporium cladosporioides dan Cladosporium halotolerans, bukan urutan langsung C. sphaerospermum.
- Aksesi DQ000286 (Borrelia sp. 16S rRNA) yang sebelumnya tersimpan sebagai cladosporium_sample.fasta telah dipindahkan ke data/raw/rejected/.

## Assembly referensi (tidak diunduh)

- `GCA_001931905.1` — Cladosporium sphaerospermum IMV 00045: Draft genome isolat dari ISS; file genomic.gbff tidak memiliki CDS teranotasi (hanya fitur source).
- `GCA_000261425.2` — Cladosporium sphaerospermum UM 843: Assembly scaffold; direktori FTP NCBI tidak menyediakan protein FASTA maupun anotasi gen.

## Rekaman yang ditolak

- `rejected/DQ000286.1_borrelia_16S.fasta` — aksesi DQ000286 adalah
  *Borrelia* sp. IK/23 16S rRNA, bukan *Cladosporium sphaerospermum*.

