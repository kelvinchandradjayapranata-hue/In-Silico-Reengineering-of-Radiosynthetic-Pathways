# Protokol DFT/Marcus untuk Transfer Elektron DHN → NAD(P)+

Dokumen ini menjelaskan cara memperoleh parameter teori transfer elektron
(ΔG°, λ, H_AB) dan menghitung laju Marcus dengan
`src/marcus_electron_transfer.py`. Dokumen ini adalah protokol komputasi,
bukan klaim bahwa jalur transfer elektron satu-elektron DHN → NAD+ terjadi
secara biologis.

## 1. Kesimpulan utama sebelum mulai

1. **NAD+/NADP+ bukan akseptor elektron luar-sphere yang sederhana.** Dalam
   metabolisme normal, NAD+/NADP+ menerima **hidrida** (2e⁻ + H⁺) pada posisi
   C4 cincin nikotinamida. Transfer elektron satu-elektron menghasilkan radikal
   NAD• yang reaktif dan dapat dimerisasi atau tereduksi lebih lanjut.
2. **Oksidasi DHN bukan ET luar-sphere murni.** DHN (1,8-dihidroksinaftalena)
   teroksidasi menjadi naftokuinon melalui 2e⁻/2H⁺; pada pH fisiologis sebagian
   gugus fenol dapat terdeprotonasi. Jalur yang relevan adalah **PCET**
   (proton-coupled electron transfer), bukan hanya transfer elektron.
3. **Marcus tetap berguna sebagai model pembanding.** Model ini memberi skala
   orde laju dan menunjukkan parameter mana yang mengontrol (ΔG°, λ, H_AB,
   jarak). Gunakan untuk:
   - akseptor satu-elektron yang benar (mis. ferisianida, viologen, sitokrom c,
     permukaan elektrode);
   - skenario NAD+/NAD• yang dilabeli eksplisit sebagai **non-fisiologis**;
   - studi sensitivitas, bukan prediksi jalur biologis.
4. **Untuk klaim hidrida/PCET**, gunakan model DFT hidrida atau PCET (EVB,
   QM/MM, CDFT terkopel proton), bukan persamaan Marcus satu-elektron.
5. **Semua angka DFT harus divalidasi** terhadap elektrokimia/kinetika
   eksperimental sebelum dipakai untuk klaim energi radiasi.

## 2. Tiga jalur yang harus dibedakan

| Jalur | Reaksi | Model yang tepat | Catatan |
|---|---|---|---|
| Hidrida fisiologis | DHN–H + NAD+ → DHN(ox) + NADH | DFT hidrida affinity, EVB, QM/MM | 2e⁻/H⁺; bukan Marcus 1e⁻ |
| PCET oksidasi DHN | DHN ⇌ DHN• + H⁺ + e⁻ | PCET/Marcus terkopel proton, CDFT | pH-dependent; semikuinon |
| ET luar-sphere 1e⁻ | D + A → D⁺ + A⁻ | Marcus non-adiabatik | hanya untuk akseptor 1e⁻ yang sah |

Jika tetap ingin memakai NAD+ dalam model Marcus, tulis di laporan:
"model satu-elektron NAD+/NAD•, bukan jalur hidrida fisiologis". Jangan
menyamakan ΔG° model ini dengan ΔG° NAD+/NADH (−0,320 V vs NHE, 2e⁻/H⁺).

## 3. Termodinamika: ΔG°

Untuk D + A → D⁺ + A⁻ pada satu elektron:

```text
ΔG° = F (E°_D − E°_A)  [J/mol]
ΔG° = n (E°_D − E°_A)  [eV per molekul]
```

- `E°_D` = potensial oksidasi donor, yaitu E°(D⁺/D).
- `E°_A` = potensial reduksi akseptor, yaitu E°(A/A⁻).
- Keduanya harus terhadap referensi elektroda yang sama (NHE/SHE).
- Konversi referensi: E(vs NHE) = E(vs SCE) + 0,241 V;
  E(vs NHE) = E(vs Ag/AgCl, KCl jenuh) + 0,197 V.
- Jika yang diukur adalah puncak oksidasi CV, itu bukan otomatis E°;
  gunakan potensial formal dari analisis reversibel atau ekstrapolasi.
- Untuk proses PCET dengan m proton dan n elektron pada 298 K:
  `E = E° − (0,05916 V) (m/n) pH`. Jangan memakai koreksi ini untuk proses
  yang bukan PCET.

**Nilai literatur yang sering disalahgunakan.** E°′ NAD+/NADH = −0,320 V vs
NHE pada pH 7 adalah proses **dua elektron + satu proton**. Pasangan
satu-elektron NAD+/NAD• kira-kira berada sekitar −1,1 V vs NHE (nilai
bervariasi antar kondisi; verifikasi dengan data primer). Konsekuensinya,
transfer 1e⁻ dari DHN (E° oksidasi kira-kira 0,0 sampai +0,5 V vs NHE pada
pH 7, tergantung mikrostate) ke NAD+ biasanya **sangat uphill** dan tidak
kompetitif dengan jalur hidrida. Ini alasan fisik, bukan sekadar formalitas,
untuk memisahkan model.

**Rekomendasi:** tentukan E° DHN dan E° pasangan akseptor yang benar dengan
voltametri siklik pada pH dan buffer yang sama, atau hitung dengan DFT +
SMD dan kalibrasi terhadap sistem referensi. Laporkan ketidakpastian ±50 mV
(≈0,05 eV) atau lebih.

## 4. Reorganisasi inner-sphere: metode empat titik

Untuk donor yang teroksidasi (D → D⁺) dan akseptor yang tereduksi (A → A⁻):

```text
λ_D = [E_D+(g_D)  − E_D+(g_D+)] + [E_D(g_D+)  − E_D(g_D)]
λ_A = [E_A−(g_A)  − E_A−(g_A−)] + [E_A(g_A−)  − E_A(g_A)]
λ_i = λ_D + λ_A
```

`g_D` = geometri teroptimasi state netral; `g_D+` = geometri teroptimasi
state kationik; analog untuk akseptor anionik. Semua energi harus dari level
teori, basis, solvasi, dan koreksi termal yang sama.

Praktik yang disarankan:

- Optimasi + frekuensi di fase larutan (SMD air) untuk setiap state.
- Pastikan tidak ada frekuensi imajiner (minimum sejati).
- Gunakan energi elektronik + ZPE, atau Gibbs, secara **konsisten** untuk
  keempat titik. Jangan mencampur E elektronik dengan G.
- Untuk anion, gunakan basis dengan fungsi difus (def2-TZVPD, ma-def2-TZVP,
  atau aug-cc-pVTZ). Tanpa fungsi difus, energi anion dan λ_A bias.
- Periksa kontaminasi spin (`<S²>` mendekati 0,75 untuk doublet), terutama
  untuk radikal NAD• dan semikuinon DHN.
- Koreksi BSSE dengan counterpoise hanya untuk kompleks dimer; untuk monomer
  four-point biasanya tidak diperlukan, tetapi basis harus konsisten.
- Standar state: DFT gas-phase memakai 1 atm; larutan 1 M. Jika menghitung
  ΔG absolut, tambahkan koreksi 1,89 kcal/mol per spesies terlarut pada 298 K.
  Untuk potensial redoks, koreksi ini sering saling menghilangkan, tetapi
  dokumentasikan.

Contoh JSON yang diterima `src/marcus_electron_transfer.py`:

```json
{
  "donor": {
    "neutral_at_neutral_eV": 0.0,
    "cation_at_neutral_eV": 7.00,
    "cation_at_cation_eV": 6.80,
    "neutral_at_cation_eV": 0.10
  },
  "acceptor": {
    "neutral_at_neutral_eV": 0.0,
    "anion_at_neutral_eV": -1.00,
    "anion_at_anion_eV": -1.25,
    "neutral_at_anion_eV": 0.10
  }
}
```

Nilai di atas hanya contoh sintetis; jangan dipakai sebagai hasil ilmiah.

## 5. Reorganisasi outer-sphere: kontinum Marcus

```text
λ_s = (e² / 8π ε₀) (1/ε_op − 1/ε_s) (1/r_D + 1/r_A − 1/R_DA)
```

- `e²/(4π ε₀) = 14,399645 eV·Å`, sehingga faktor depan = 7,1998225 eV·Å.
- `ε_op = n²` (kuadrat indeks bias optik). Air: n = 1,333 → ε_op ≈ 1,776;
  ε_s ≈ 78,355 pada 25 °C.
- `r_D`, `r_A` = radius cavity efektif (Å), biasanya 3,0–3,5 Å untuk molekul
  kecil; lebih baik diturunkan dari volume cavity SMD/PCM DFT:
  `r = (3V/(4π))^(1/3)`.
- `R_DA` = jarak pusat donor–akseptor (Å), harus ≥ r_D + r_A.
- Untuk interior protein, ε_s efektif kira-kira 4–10 dan ε_op kira-kira
  2,0–2,5; jangan memakai ε air 78 untuk lingkungan protein tanpa pembenaran.
- Model ini menganggap pelarut kontinum linear dan bola keras. Ia tidak
  menangkap reorganisasi spesifik protein, water wires, atau gerak gugus
  polar. Gunakan MD/energetika gap untuk koreksi (lihat §8).

## 6. Kopling elektronik H_AB

Pilihan metode, dari yang paling praktis:

1. **GMH dua-state (TD-DFT).** Hitung state adiabatik 1 dan 2 untuk dimer
   donor–akseptor. Dengan `ΔE_12` (eV), `μ_12` (transition dipole, Debye), dan
   `Δμ_12` (selisih dipole adiabatik):
   `H_AB = ΔE_12 μ_12 / sqrt(Δμ_12² + 4 μ_12²)`.
   Untuk sistem simetris, `H_AB = ΔE_12/2`.
2. **CDFT (constrained DFT).** Q-Chem menyediakan CDFT; definisikan fragmen
   donor dan akseptor, hitung state charge-localized, lalu ekstrak kopling.
   Laporkan definisi constraint dan functional.
3. **FCD (fragment charge difference).** Implementasi di PySCF/Q-Chem;
   berguna untuk dimer yang lebih besar.
4. **Model eksponensial jarak.** Jika DFT dimer tidak memungkinkan:
   `H_AB(R) = H_0 exp[−β(R − R_0)/2]`, dengan β kira-kira 1,0–1,4 Å⁻¹ untuk
   protein. Ini hanya skala, bukan prediksi.

Praktik:

- Optimasi dimer dengan dispersion correction (D3(BJ)); tanpa koreksi
  dispersi, π-stacking DHN oligomer akan terlalu longgar.
- Hitung pada beberapa jarak/orientasi dari MD; jangan hanya satu pose docking.
- Periksa apakah state tereksitasi benar-benar CT (charge-transfer), bukan
  eksitasi lokal. Lihat perubahan densitas dan dipole.
- Untuk dimer dengan kopling kuat (H_AB sebanding λ), model non-adiabatik
  mungkin tidak berlaku; pertimbangkan hopping adiabatik/Sumi–Marcus.
- Laporkan rentang H_AB dari beberapa metode; perbedaan faktor 2–3 lazim.

## 7. Level teori yang disarankan

| Tahap | Metode | Catatan |
|---|---|---|
| Screening konformer | GFN2-xTB / CREST | cepat; pilih <3 kcal/mol |
| Optimasi + frekuensi | ωB97X-D3/def2-TZVP, SMD(air) | produksi; RIJCOSX |
| Single-point anion | def2-TZVPD atau ma-def2-TZVP | fungsi difus wajib |
| Kopling | TD-DFT/GMH atau CDFT | functional sama bila mungkin |
| Benchmark | M06-2X/def2-TZVP, ωB97X-D/def2-TZVP | bandingkan λ dan E° |
| Solvasi | SMD (air); ε protein dari MD | jangan campur model sembarangan |

Catatan functional:

- ωB97X-D3 (range-separated + D3) adalah pilihan awal yang baik untuk CT,
  tetapi tidak ada functional yang universal.
- M06-2X sering baik untuk termokimia/konformer organik, tetapi dapat
  melebih-lebihkan CT jarak jauh.
- Untuk anion dan radikal, periksa stabilitas SCF dan kontaminasi spin.
- Jangan memakai B3LYP tanpa koreksi dispersi untuk π-stacking DHN.
- Untuk sistem ini (C, H, O, N), efek relativistik dapat diabaikan.

Contoh input ORCA (optimasi + frekuensi, state netral):

```text
! wB97X-D3 def2-TZVP def2/J RIJCOSX SMD(water) Opt Freq TightSCF
%pal nprocs 8 end
%maxcore 4000
* xyzfile 0 1 dhn_neutral.xyz
```

Untuk kation donor: `* xyzfile 1 2 dhn_cation.xyz`.
Untuk anion akseptor: `* xyzfile -1 2 nad_anion.xyz`.
Untuk single-point TD-DFT dimer, tambahkan `%tddft nroots 10 end` dan
pastikan geometri dimer sudah teroptimasi.

Contoh input Gaussian (alternatif):

```text
#p wb97xd/def2tzvp scrf=(smd,solvent=water) opt freq
```

Periksa manual versi software untuk nama basis dan keyword SMD yang tepat.

## 8. Dari MD ke parameter: ensemble dan reorganisasi

Model Marcus klasik menganggap satu geometri. Protein dan oligomer DHN
memiliki distribusi konformasi, sehingga laju harus dirata-ratakan:

1. Ekstrak 100–500 snapshot dari trajektori produksi (setelah ekuilibrasi).
2. Hitung `R_DA`, orientasi, dan kontak π-stacking.
3. Hitung `H_AB(R)` dengan model eksponensial atau DFT pada subset snapshot.
4. Hitung λ_s per snapshot dari fluktuasi energy gap (pendekatan linear
   response):

   ```text
   λ_s = Var(ΔE) / (2 k_B T)
   ΔG  = <ΔE>_reactant − λ_s
   ```

   dengan `ΔE = E(product) − E(reactant)` dievaluasi pada ensemble reactant.
   Ini membutuhkan parameterisasi state redoks (QM/MM atau EVB); jangan
   memakai muatan force field ground state untuk state teroksidasi tanpa
   validasi.
5. Hitung `k_i` per snapshot, lalu laporkan `<k>`, median, dan interval
   5–95%. **Jangan** menghitung `k(<ΔG>, <λ>)`; rata-rata parameter lalu
   memasukkan ke eksponensial memberi bias (Jensen).
6. Jika distribusi `R_DA` lebar, laju efektif didominasi konformasi jarak
   pendek. Laporkan distribusi, bukan hanya satu angka.

Untuk oligomer DHN:

- Modelkan dimer/trimer eksplisit; jangan menganggap monomer identik.
- Periksa lokalisasi hole/radikal (spin density). Jika hole terdelokalisasi
  kuat, model hopping situs tunggal gagal.
- Coupling antarmonomer sangat sensitif terhadap jarak π–π dan slip; rata-rata
  dari MD, bukan satu struktur docking.
- Untuk film melanin terhidrasi, transport dapat melibatkan hopping yang
  dibantu air dan reaksi komproportionasi kuinon/hidrokuinon; model Marcus
  sederhana hanya bagian kecil dari gambaran.

## 9. Menghitung laju dengan skrip

Persamaan yang dipakai:

```text
k_ET = (2π/ħ) |H_AB|² / sqrt(4π λ k_B T)
       × exp[−(ΔG° + λ)² / (4 λ k_B T)]
```

Jika `−ΔG° > λ`, sistem berada di **inverted region**: laju menurun lagi
ketika driving force bertambah. Ini adalah perilaku fisik Marcus, bukan
kesalahan numerik; pastikan tanda ΔG° dan definisi donor/akseptor benar.

Jika ada mode vibrasi frekuensi tinggi (mis. C=O atau C–H) yang tidak
termodulasi pelarut, gunakan Marcus–Levich–Jortner:

```text
k_MLJ = (2π/ħ) |H_AB|² / sqrt(4π λ_s k_B T)
        × Σ_n (e^−S S^n / n!)
        × exp[−(ΔG° + λ_s + n ħω)² / (4 λ_s k_B T)]
S = λ_v / (ħω)
```

`λ_s` adalah reorganisasi outer-sphere/pelarut dan `λ_v` adalah reorganisasi
mode kuantum. Skrip `marcus_electron_transfer.py` memakai bentuk klasik;
untuk MLJ, hitung suku-suku `n` secara terpisah atau gunakan skrip analisis
sendiri.

Contoh parameter langsung:

```bash
python src/marcus_electron_transfer.py \
  --delta-g-ev -0.18 \
  --lambda-inner-ev 0.45 \
  --lambda-outer-ev 0.85 \
  --coupling 1.0 --coupling-unit meV
```

Contoh dari potensial redoks + outer-sphere:

```bash
python src/marcus_electron_transfer.py \
  --donor-e0-v -0.50 --acceptor-e0-v -0.32 \
  --lambda-inner-ev 0.45 \
  --donor-radius-a 3.2 --acceptor-radius-a 3.4 --distance-a 6.5 \
  --eps-static 78.355 --eps-optical 1.776 \
  --coupling 1.0 --coupling-unit meV
```

Contoh dari four-point DFT:

```bash
python src/marcus_electron_transfer.py \
  --four-point-json data/quantum/dhn_nad_four_point.json \
  --lambda-outer-ev 0.85 \
  --delta-g-ev -0.18 \
  --coupling 1.0 --coupling-unit meV
```

Semua parameter juga bisa diletakkan di satu file config:

```json
{
  "delta_g_ev": -0.18,
  "four_point": {
    "donor": {
      "neutral_at_neutral_eV": 0.0,
      "cation_at_neutral_eV": 7.00,
      "cation_at_cation_eV": 6.80,
      "neutral_at_cation_eV": 0.10
    },
    "acceptor": {
      "neutral_at_neutral_eV": 0.0,
      "anion_at_neutral_eV": -1.00,
      "anion_at_anion_eV": -1.25,
      "neutral_at_anion_eV": 0.10
    }
  },
  "donor_radius_angstrom": 3.2,
  "acceptor_radius_angstrom": 3.4,
  "distance_angstrom": 6.5,
  "dielectric_static": 78.355,
  "dielectric_optical": 1.776,
  "coupling_ev": 0.001
}
```

```bash
python src/marcus_electron_transfer.py --config data/quantum/marcus_config.json
```

Dalam config, `coupling_ev` selalu dalam eV; gunakan `coupling` bersama
`coupling_unit` (`eV`, `meV`, `cm-1`) untuk satuan lain.

Output default: `results/marcus_electron_transfer.json` dan
`results/marcus_electron_transfer.md`. Gunakan `--json`/`--md` untuk path lain.

## 10. Sensitivitas dan pelaporan

Karena `k` bergantung eksponensial pada ΔG° dan λ, kesalahan kecil
menghasilkan orde besaran:

```text
∂ln k/∂ΔG° = −(ΔG° + λ) / (2 λ k_B T)
∂ln k/∂λ   = −(ΔG° + λ)(λ − ΔG°) / (4 λ² k_B T)
```

Pada ΔG° = −0,18 eV, λ = 1,30 eV, T = 298 K:
±0,05 eV pada ΔG° ≈ faktor 2,3 pada k; ±0,10 eV pada λ ≈ faktor 2,6.
H_AB tidak pasti faktor 2 → k berubah faktor 4. Laporkan rentang, bukan satu
angka.

Checklist pelaporan:

- [ ] Reaksi dan state redoks yang dimodelkan (1e⁻ vs 2e⁻/H⁺).
- [ ] Level teori, basis, functional, dispersion, solvasi.
- [ ] Charge/multiplicity, `<S²>`, stabilitas SCF.
- [ ] Geometri sumber (DFT, docking, MD snapshot) dan protonasi.
- [ ] Four-point: empat energi per state, ZPE/G, BSSE.
- [ ] λ_i, λ_s, λ total, dan asal setiap komponen.
- [ ] ΔG° dan referensi elektroda/pH.
- [ ] H_AB: metode, jarak, orientasi, rentang.
- [ ] Ensemble: jumlah snapshot, rata-rata, median, interval.
- [ ] Sensitivitas dan ketidakpastian.
- [ ] Pernyataan bahwa hasil bukan bukti eksperimental.

## 11. Validasi eksperimental yang disarankan

1. **Elektrokimia:** CV/DPV DHN dan NAD(P)+ pada pH, buffer, dan elektrode
   yang sama; tentukan E° formal dan reversibilitas. Bandingkan dengan DFT
   (target ±50 mV).
2. **Kinetika:** stopped-flow atau flash photolysis dengan akseptor 1e⁻ yang
   diketahui; ukur k_ET pada beberapa suhu dan jarak.
3. **Efek isotop:** H/D pada DHN untuk membedakan PCET dari ET murni.
4. **pH dependence:** kemiringan E vs pH menunjukkan stoikiometri proton.
5. **Ionic strength:** uji apakah laju mengikuti model elektrostatik atau
   kontrol konformasi.
6. **Benchmark internal:** hitung sistem referensi (mis. ferrosen/ferrisianida,
   sitokrom c) dengan protokol yang sama; jika meleset >1–2 orde, jangan
   percaya angka DHN.
7. **Kontrol negatif:** NADH vs NAD+; donor tanpa gugus kuinon; oligomer tanpa
   stacking.

## 12. Batas komputasi dan tempat eksekusi

- DFT produksi (opt+freq, TD-DFT, CDFT) berjalan di **CPU**, bukan GPU.
- Molekul kecil (DHN, NAD) masih mungkin di Colab CPU untuk screening, tetapi
  batas waktu/CPU/memori Colab tidak cocok untuk produksi dimer/oligomer.
- Gunakan Colab/HPC untuk DFT; jangan menjalankan DFT berat di CPU lokal.
- ORCA gratis untuk akademik; Gaussian/Q-Chem komersial; PySCF/PSI4
  open-source tetapi butuh setup SMD/CDFT tambahan.
- Simpan input, output, dan versi software di `data/quantum/` (di-ignore Git
  bila besar) dan catat checksum untuk audit.

## 13. Referensi kunci

- R. A. Marcus, *J. Chem. Phys.* **24**, 966 (1956); **67**, 853 (1964).
- R. A. Marcus & N. Sutin, *Biochim. Biophys. Acta* **811**, 265 (1985).
- J. Jortner, *J. Chem. Phys.* **64**, 4860 (1976) — mode kuantum.
- R. J. Cave & M. D. Newton, *Chem. Phys. Lett.* **249**, 15 (1996);
  *J. Chem. Phys.* **106**, 9213 (1997) — GMH.
- A. V. Marenich, C. J. Cramer & D. G. Truhlar, *J. Phys. Chem. B* **113**,
  6378 (2009) — SMD.
- J.-D. Chai & M. Head-Gordon, *J. Chem. Phys.* **128**, 084106 (2008);
  **132**, 154104 (2010) — ωB97X/ωB97X-D.
- S. Grimme et al., *J. Chem. Phys.* **132**, 154104 (2010);
  *J. Comput. Chem.* **32**, 1456 (2011) — D3/D3(BJ).
- Y. Zhao & D. G. Truhlar, *Theor. Chem. Acc.* **120**, 215 (2008) — M06-2X.
- F. Neese, *WIREs Comput. Mol. Sci.* **2**, 73 (2012); ORCA manual.
- Y. Shao et al., *Mol. Phys.* **113**, 184 (2015) — Q-Chem.
- D. Blumberger, *Chem. Rev.* **115**, 11191 (2015) — review ET biologis.
- M. J. Frisch et al., Gaussian 16 (2016).

Untuk nilai E° NAD+/NAD• dan DHN, gunakan data primer/elektrokimia terbaru;
nilai dalam dokumen ini hanya jangkar orde besaran dan harus diverifikasi.
