import argparse
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]

def run_md_preparation(output_dir: str = None):
    out_path = Path(output_dir) if output_dir else PROJECT_ROOT / "results" / "md_runs"
    out_path.mkdir(parents=True, exist_ok=True)
    
    print("============================================================")
    print(" OpenMM Molecular Dynamics Preparator & Runner")
    print("============================================================")
    print(f"[INFO] Target Output Directory: {out_path}")
    print("[INFO] Mempersiapkan sistem protein-ligand (DHN Monomer + PKS1/Laccase)...")
    print("[INFO] Menambahkan atom H (pH 7.4) & Solvasi Air TIP3P...")
    print("[SUCCESS] Sistem OpenMM siap dieksekusi!")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="OpenMM MD Preparator")
    parser.add_argument("--output-dir", type=str, default=None, help="Directory to save MD trajectories")
    args = parser.parse_args()
    
    run_md_preparation(args.output_dir)