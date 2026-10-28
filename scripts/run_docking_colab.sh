#!/usr/bin/env bash
# =============================================================================
# Docking presisi tinggi AutoDock Vina di Google Colab (CPU runtime)
# =============================================================================
# Vina tidak memakai GPU. Jalankan notebook ini dengan Runtime > Change runtime
# type > CPU agar CPU lokal tidak terbebani.
#
# Pemakaian cepat di Colab:
#   !git clone https://github.com/kelvinchandradjayapranata-hue/In-Silico-Reengineering-of-Radiosynthetic-Pathways.git
#   %cd In-Silico-Reengineering-of-Radiosynthetic-Pathways
#   !bash scripts/run_docking_colab.sh
#
# Konfigurasi opsional (set sebelum memanggil script):
#   SEEDS="42 43 44"          # seed Vina; kurangi bila waktu terbatas
#   EXHAUSTIVENESS=32         # naikkan ke 64 untuk pencarian lebih luas
#   NUM_MODES=20
#   BOX_SIZE="24 24 24"       # ukuran kotak (Angstrom)
#   LIGANDS="DHN_monomer DHN_dimer DHN_trimer"
#   RECEPTORS=""              # kosong = semua reseptor di config
#   CPU="$(nproc)"            # jangan melebihi jumlah vCPU Colab
#   FORCE=0                   # 1 = regenerasi ligan/reseptor
#   OUTPUT_DIR="/content/drive/MyDrive/melanin_docking_runs"
#   MOUNT_DRIVE=1             # mount Google Drive sebelum menyalin hasil
#   SKIP_PULL=1               # lewati git pull
#   GITHUB_TOKEN="..."        # token read-only bila repo masih privat
#
# Contoh:
#   !SEEDS="42 43" EXHAUSTIVENESS=32 OUTPUT_DIR="/content/drive/MyDrive/melanin_docking_runs" bash scripts/run_docking_colab.sh
# =============================================================================

set -euo pipefail

REPO_URL="${REPO_URL:-https://github.com/kelvinchandradjayapranata-hue/In-Silico-Reengineering-of-Radiosynthetic-Pathways.git}"
GITHUB_TOKEN="${GITHUB_TOKEN:-}"
WORKDIR="${WORKDIR:-/content/In-Silico-Reengineering-of-Radiosynthetic-Pathways}"
SEEDS="${SEEDS:-42 43 44}"
EXHAUSTIVENESS="${EXHAUSTIVENESS:-32}"
NUM_MODES="${NUM_MODES:-20}"
BOX_SIZE="${BOX_SIZE:-24 24 24}"
LIGANDS="${LIGANDS:-DHN_monomer DHN_dimer DHN_trimer}"
RECEPTORS="${RECEPTORS:-}"
CPU="${CPU:-$(nproc)}"
FORCE="${FORCE:-0}"
OUTPUT_DIR="${OUTPUT_DIR:-}"
MOUNT_DRIVE="${MOUNT_DRIVE:-0}"
SKIP_PULL="${SKIP_PULL:-0}"
VINA_URL="${VINA_URL:-https://github.com/forlilab/AutoDock-Vina/releases/download/v1.2.7/vina_1.2.7_linux_x86_64}"
VINA_BIN="${VINA_BIN:-}"

echo "============================================================"
echo " AutoDock Vina di Colab (CPU)"
echo "============================================================"
echo " Working dir : $WORKDIR"
echo " Seeds       : $SEEDS"
echo " Exhaustive  : $EXHAUSTIVENESS"
echo " Num modes   : $NUM_MODES"
echo " Box size    : $BOX_SIZE"
echo " CPU         : $CPU"
echo "============================================================"

if command -v nvidia-smi >/dev/null 2>&1; then
    echo "[INFO] GPU terdeteksi, tetapi Vina tetap memakai CPU. Ini normal."
fi

# ---------------------------------------------------------------------------
# 1. Clone / update repository
# ---------------------------------------------------------------------------
if [ ! -d "$WORKDIR/.git" ]; then
    mkdir -p "$(dirname "$WORKDIR")"
    echo "[1/5] Meng-clone repository..."
    CLONE_URL="$REPO_URL"
    if [ -n "$GITHUB_TOKEN" ]; then
        # Gunakan token read-only dari Colab Secrets untuk repo privat.
        CLONE_URL="${REPO_URL/https:\/\/github.com\//https://x-access-token:${GITHUB_TOKEN}@github.com/}"
    fi
    git clone "$CLONE_URL" "$WORKDIR"
    # Jangan simpan token di .git/config.
    git -C "$WORKDIR" remote set-url origin "$REPO_URL"
fi
cd "$WORKDIR"

if [ "$SKIP_PULL" != "1" ] && [ -z "$GITHUB_TOKEN" ]; then
    echo "[1/5] Memperbarui repository (fast-forward)..."
    git pull --ff-only || echo "[WARN] git pull dilewati (tree lokal tidak bersih)."
elif [ -n "$GITHUB_TOKEN" ]; then
    echo "[1/5] Melewati git pull (clone privat memakai token)."
else
    echo "[1/5] Melewati git pull."
fi

# ---------------------------------------------------------------------------
# 2. Dependency Python untuk docking
# ---------------------------------------------------------------------------
echo "[2/5] Memeriksa dependency docking..."
python - <<'PY'
import importlib
import importlib.metadata
import subprocess
import sys

# Hanya pasang paket yang belum tersedia agar environment Colab tidak
# diubah lebih dari yang diperlukan.
required = {
    "Bio": "biopython==1.88",
    "rdkit": "rdkit==2026.9.1",
    "meeko": "meeko==0.8.0",
    "gemmi": "gemmi==0.7.5",
    "scipy": "scipy==1.18.1",
}
missing = []
for module, spec in required.items():
    try:
        importlib.import_module(module)
    except Exception:
        missing.append(spec)

if missing:
    print("  Meng-install:", ", ".join(missing))
    subprocess.check_call([sys.executable, "-m", "pip", "install", "-q", *missing])

for package in ("rdkit", "meeko", "gemmi", "scipy", "biopython"):
    try:
        print(f"  {package}: {importlib.metadata.version(package)}")
    except importlib.metadata.PackageNotFoundError:
        raise SystemExit(f"[ERROR] Package {package} tidak terpasang.")
PY

# ---------------------------------------------------------------------------
# 3. Binary AutoDock Vina (Linux x86_64, static)
# ---------------------------------------------------------------------------
if [ -z "$VINA_BIN" ]; then
    VINA_BIN="$WORKDIR/tools/vina/vina"
fi
if [ ! -x "$VINA_BIN" ]; then
    mkdir -p "$(dirname "$VINA_BIN")"
    echo "[3/5] Mengunduh AutoDock Vina 1.2.7..."
    if command -v curl >/dev/null 2>&1; then
        curl -fL --retry 3 --retry-delay 2 -o "$VINA_BIN" "$VINA_URL"
    else
        wget -O "$VINA_BIN" "$VINA_URL"
    fi
    chmod +x "$VINA_BIN"
else
    echo "[3/5] Binary Vina sudah tersedia."
fi
"$VINA_BIN" --version

# ---------------------------------------------------------------------------
# 4. Jalankan beberapa seed
# ---------------------------------------------------------------------------
RUN_ID="${RUN_ID:-$(date -u +%Y%m%dT%H%M%SZ)}"
RUN_DIR="$WORKDIR/results/colab_runs/$RUN_ID"
mkdir -p "$RUN_DIR"

read -r -a SEED_ARRAY <<< "$SEEDS"
if [ "${#SEED_ARRAY[@]}" -eq 0 ]; then
    echo "[ERROR] SEEDS kosong."
    exit 1
fi

read -r -a BOX_ARRAY <<< "$BOX_SIZE"
if [ "${#BOX_ARRAY[@]}" -ne 3 ]; then
    echo "[ERROR] BOX_SIZE harus berisi 3 angka, contoh: BOX_SIZE='24 24 24'."
    exit 1
fi

read -r -a LIGAND_ARRAY <<< "$LIGANDS"
read -r -a RECEPTOR_ARRAY <<< "$RECEPTORS"

FIRST_SEED="${SEED_ARRAY[0]}"
for seed in "${SEED_ARRAY[@]}"; do
    echo
    echo "------------------------------------------------------------"
    echo "[4/5] Menjalankan seed $seed"
    echo "------------------------------------------------------------"
    cmd=(
        python src/melanin_docking.py
        --exhaustiveness "$EXHAUSTIVENESS"
        --num-modes "$NUM_MODES"
        --cpu "$CPU"
        --seed "$seed"
        --box-size "${BOX_ARRAY[@]}"
        --vina "$VINA_BIN"
        --out "$RUN_DIR/docking_seed${seed}.json"
        --summary "$RUN_DIR/docking_seed${seed}.tsv"
    )
    if [ "${#LIGAND_ARRAY[@]}" -gt 0 ]; then
        cmd+=(--ligands "${LIGAND_ARRAY[@]}")
    fi
    if [ "${#RECEPTOR_ARRAY[@]}" -gt 0 ]; then
        cmd+=(--receptors "${RECEPTOR_ARRAY[@]}")
    fi
    # Regenerasi ligan/reseptor cukup sekali; seed berikutnya memakai file yang sama.
    if [ "$FORCE" = "1" ] && [ "$seed" = "$FIRST_SEED" ]; then
        cmd+=(--force)
    fi

    "${cmd[@]}" 2>&1 | tee "$RUN_DIR/seed${seed}.log"
done

# ---------------------------------------------------------------------------
# 5. Agregasi dan arsip
# ---------------------------------------------------------------------------
echo
echo "[5/5] Mengagregasi hasil multi-seed..."
python scripts/aggregate_docking_seeds.py \
    --input-dir "$RUN_DIR" \
    --out "$WORKDIR/results/docking_results.json" \
    --summary "$WORKDIR/results/docking_summary.tsv"

cp "$WORKDIR/results/docking_results.json" "$RUN_DIR/"
cp "$WORKDIR/results/docking_summary.tsv" "$RUN_DIR/"

echo "Mengarsipkan pose dan log..."
tar -czf "$RUN_DIR/docking_poses_logs.tar.gz" \
    -C "$WORKDIR/results" docking/poses docking/logs

if [ -n "$OUTPUT_DIR" ]; then
    if [ "$MOUNT_DRIVE" = "1" ]; then
        python -c "from google.colab import drive; drive.mount('/content/drive')"
    fi
    mkdir -p "$OUTPUT_DIR/$RUN_ID"
    cp -r "$RUN_DIR"/. "$OUTPUT_DIR/$RUN_ID"/
    echo "[INFO] Hasil disalin ke: $OUTPUT_DIR/$RUN_ID"
fi

echo
echo "============================================================"
echo " Selesai"
echo " Run dir : $RUN_DIR"
echo " JSON    : $WORKDIR/results/docking_results.json"
echo " TSV     : $WORKDIR/results/docking_summary.tsv"
echo " Arsip   : $RUN_DIR/docking_poses_logs.tar.gz"
echo "============================================================"
