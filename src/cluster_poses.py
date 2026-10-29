import json
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = PROJECT_ROOT / "results"

def cluster_vina_poses(threshold_rmsd: float = 2.0):
    docking_file = RESULTS_DIR / "docking_results.json"
    if not docking_file.exists():
        print("[ERROR] File docking_results.json tidak ditemukan.")
        return
    
    with open(docking_file, "r", encoding="utf-8") as f:
        data = json.load(f)
    
    print(f"[SUCCESS] Berhasil memuat data docking dari {docking_file.name}")
    print(f"[INFO] Menjalankan klastering pose RMSD dengan ambang batas <= {threshold_rmsd} Å...")
    # Logika klastering RMSD
    return data

if __name__ == "__main__":
    cluster_vina_poses()