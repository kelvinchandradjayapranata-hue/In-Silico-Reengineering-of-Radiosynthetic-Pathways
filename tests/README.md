# Tests

Belum ada test suite otomatis.

Rencana pengujian:

- Validasi manifest NCBI: jumlah target, organisme, dan tidak ada aksesi
  *Borrelia*/non-target.
- Validasi descriptor RDKit: formula DHN = C10H8O2 dan DOPA = C9H11NO4.
- Smoke test FBA: status optimal, DHN/DOPA export > 0, dan SBML dapat dibaca
  ulang.
- Regresi: objective value berubah jika batas input radiasi/glukosa diubah.

Jalankan secara manual setelah ditambahkan, misalnya:

```powershell
python -m unittest discover -s tests
```
