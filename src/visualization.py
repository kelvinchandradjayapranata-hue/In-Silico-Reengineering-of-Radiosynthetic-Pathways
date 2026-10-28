"""Visualisasi fluks FBA radiosintesis melanin."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Sequence

import matplotlib.pyplot as plt

PROJECT_ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = PROJECT_ROOT / "results"
PLOTS_DIR = RESULTS_DIR / "plots"
FBA_RESULTS_PATH = RESULTS_DIR / "metabolic_fba.json"
DEFAULT_OUTPUT_PATH = PLOTS_DIR / "fba_fluxes.png"


def _color_for_reaction(reaction_id: str) -> str:
    """Warna berdasarkan kelompok jalur."""
    if reaction_id.startswith("RAD_"):
        return "#e76f51"  # orange: radiasi
    if reaction_id.startswith(("TYROSINASE", "DOPA", "MELANIN_DOPA")):
        return "#7b2cbf"  # ungu: DOPA-melanin
    if reaction_id.startswith(
        ("PKS", "THN", "SCYTALONE", "T3HN", "VERMELONE", "LACCASE", "MELANIN_DHN")
    ):
        return "#2a9d8f"  # hijau: DHN-melanin
    if reaction_id in {"EX_GLC", "GLYCOLYSIS", "PPP", "ACCARBOXYLASE"}:
        return "#457b9d"  # biru: karbon pusat
    return "#6c757d"  # abu-abu: maintenance/respirasi/lainnya


def plot_fba_fluxes(
    fba_json: Path = FBA_RESULTS_PATH,
    output_path: Path = DEFAULT_OUTPUT_PATH,
) -> Path | None:
    """Buat grafik batang horizontal untuk fluks kunci FBA."""
    fba_json = Path(fba_json)
    if not fba_json.is_file():
        print(f"[ERROR] File {fba_json} tidak ditemukan.")
        return None

    with fba_json.open("r", encoding="utf-8") as handle:
        data = json.load(handle)

    fluxes = data.get("fluxes") or data.get("all_fluxes") or {}
    if not fluxes:
        print("[ERROR] Tidak ada data fluks pada JSON.")
        return None

    # Sembunyikan fluks nol agar grafik tetap terbaca.
    nonzero = {key: value for key, value in fluxes.items() if abs(float(value)) > 1e-9}
    if not nonzero:
        print("[ERROR] Semua fluks bernilai nol.")
        return None

    items = sorted(nonzero.items(), key=lambda item: float(item[1]))
    reaction_ids = [item[0] for item in items]
    values = [float(item[1]) for item in items]
    colors = [_color_for_reaction(reaction_id) for reaction_id in reaction_ids]

    height = max(4.5, 0.42 * len(reaction_ids) + 1.8)
    fig, ax = plt.subplots(figsize=(10, height))
    bars = ax.barh(reaction_ids, values, color=colors)

    ax.set_xlabel("Fluks (mmol gDW$^{-1}$ h$^{-1}$)")
    objective_value = data.get("objective_value")
    status = data.get("status", "unknown")
    title = f"Fluks Radiosintesis Melanin (status: {status}"
    if objective_value is not None:
        title += f", objective: {float(objective_value):.3g}"
    title += ")"
    ax.set_title(title)
    ax.grid(axis="x", linestyle="--", alpha=0.45)
    ax.axvline(0.0, color="black", linewidth=0.8)

    for bar, value in zip(bars, values):
        ax.text(
            value + max(values) * 0.01,
            bar.get_y() + bar.get_height() / 2,
            f"{value:.2f}",
            va="center",
            fontsize=8,
        )

    x_padding = max(values) * 0.12 if max(values) > 0 else 0.1
    ax.set_xlim(left=min(0.0, min(values) * 1.1), right=max(values) + x_padding)

    fig.tight_layout()
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    print(f"[SUCCESS] Grafik fluks FBA disimpan ke: {output_path}")
    return output_path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Plot fluks FBA radiosintesis melanin.")
    parser.add_argument("--input", type=Path, default=FBA_RESULTS_PATH)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT_PATH)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    result = plot_fba_fluxes(args.input, args.output)
    return 0 if result is not None else 1


if __name__ == "__main__":
    raise SystemExit(main())
