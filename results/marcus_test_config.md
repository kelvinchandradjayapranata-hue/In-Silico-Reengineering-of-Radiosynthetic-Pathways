# Estimasi Laju Transfer Elektron Marcus

> Dibuat: 2026-10-29T13:22:26.730718+00:00
> Model non-adiabatik Marcus untuk transfer elektron satu-elektron.

## Parameter

| Parameter | Nilai |
|---|---|
| ΔG° | -0.180000 eV |
| λ inner-sphere | 0.650000 eV |
| λ outer-sphere | 1.981034 eV |
| λ total | 2.631034 eV |
| H_AB | 0.001000 eV (8.066 cm⁻¹) |
| Temperatur | 298.15 K |
| k_B T | 0.025693 eV |

## Hasil

| Besaran | Nilai |
|---|---|
| Energi aktivasi ΔG‡ | 0.570837 eV |
| Prefaktor | 1.035719e+10 s⁻¹ |
| Faktor eksponensial | 2.243144e-10 |
| Laju k_ET | 2.323268e+00 s⁻¹ |
| log10(k_ET / s⁻¹) | 0.37 |

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
