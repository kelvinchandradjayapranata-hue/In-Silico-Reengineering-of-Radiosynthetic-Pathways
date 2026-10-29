# Estimasi Laju Transfer Elektron Marcus

> Dibuat: 2026-10-29T13:22:14.413984+00:00
> Model non-adiabatik Marcus untuk transfer elektron satu-elektron.

## Parameter

| Parameter | Nilai |
|---|---|
| ΔG° | -0.180000 eV |
| λ inner-sphere | 0.450000 eV |
| λ outer-sphere | 0.850000 eV |
| λ total | 1.300000 eV |
| H_AB | 0.100000 eV (806.554 cm⁻¹) |
| Temperatur | 310.15 K |
| k_B T | 0.026727 eV |

## Hasil

| Besaran | Nilai |
|---|---|
| Energi aktivasi ΔG‡ | 0.241231 eV |
| Prefaktor | 1.444658e+14 s⁻¹ |
| Faktor eksponensial | 1.202608e-04 |
| Laju k_ET | 1.737357e+10 s⁻¹ |
| log10(k_ET / s⁻¹) | 10.24 |

Persamaan yang dipakai:

```text
k_ET = (2π/ħ) |H_AB|² / sqrt(4π λ k_B T)
       × exp[-(ΔG° + λ)² / (4 λ k_B T)]
```

## Asumsi dan Keterbatasan

- Model ini adalah Marcus non-adiabatik satu-elektron untuk pasangan donor-akseptor pada jarak dan orientasi tetap.
- NAD+/NADP+ secara fisiologis menerima **hidrida** (2e⁻/H⁺) pada cincin nikotinamida; transfer elektron satu-elektron ke NAD+ menghasilkan radikal NAD• yang reaktif. Model ini tidak otomatis mewakili jalur hidrida.
- Oksidasi DHN bersifat pH-dependent dan dapat berlangsung melalui PCET (deprotonasi + oksidasi), bukan ET luar-sphere sederhana.
- λ inner-sphere dari four-point DFT bergantung pada fungsi fungsional, basis set, koreksi ZPE, dan model solvasi.
- H_AB dari DFT sangat sensitif terhadap jarak, orientasi, dan tingkat teori; nilai kecil (≤10 meV) berada di batas akurasi.
- Persamaan klasik mengabaikan nuclear tunneling, efek kuantum ikatan hidrogen, dan dinamika pelarut non-kontinum.
- k_ET adalah laju unimolekuler untuk pasangan pada jarak tertentu, bukan konstanta bimolekuler; untuk larutan encer perlu koreksi difusi/encounter.
- Angka ini adalah estimasi komputasi, bukan bukti eksperimental.
