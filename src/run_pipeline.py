"""Orkestrasi pipeline radiosintesis melanin end-to-end."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Sequence

from generate_readme import generate_readme
from generate_report import generate_report
from melanin_chem import save_results as save_chem_results
from metabolic_fba import run_fba, run_sensitivity
from ncbi_fetcher import fetch_all
from visualization import plot_fba_fluxes

PROJECT_ROOT = Path(__file__).resolve().parents[1]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Jalankan pipeline radiosintesis melanin secara berurutan."
    )
    parser.add_argument(
        "--skip-fetch",
        action="store_true",
        help="Lewati unduhan NCBI dan gunakan data/raw yang sudah ada.",
    )
    parser.add_argument(
        "--force-fetch",
        action="store_true",
        help="Unduh ulang semua target NCBI meskipun cache ada.",
    )
    parser.add_argument("--electron-limit", type=float, default=10.0)
    parser.add_argument("--energy-limit", type=float, default=10.0)
    parser.add_argument("--glucose-limit", type=float, default=5.0)
    parser.add_argument("--oxygen-limit", type=float, default=20.0)
    parser.add_argument("--tyrosine-limit", type=float, default=1.0)
    parser.add_argument("--atp-maintenance", type=float, default=10.0)
    parser.add_argument("--biomass", type=float, default=0.1)
    parser.add_argument("--sensitivity", action="store_true")
    parser.add_argument("--no-sbml", action="store_true")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    if not args.skip_fetch:
        print("\n[1/6] Mengambil data NCBI...")
        results = fetch_all(force=args.force_fetch)
        errors = [row for row in results if row.get("status") == "error"]
        if errors:
            print(
                f"[WARN] {len(errors)} target gagal diunduh; pipeline tetap "
                "dilanjutkan dengan data yang tersedia."
            )
    else:
        print("\n[1/6] Fetch dilewati (--skip-fetch).")

    print("\n[2/6] Menghitung deskriptor monomer melanin...")
    chem = save_chem_results()
    print(f"       {len(chem['monomers'])} monomer diproses.")

    print("\n[3/6] Menjalankan FBA radiosintesis...")
    fba_kwargs = {
        "electron_input_limit": args.electron_limit,
        "energy_input_limit": args.energy_limit,
        "glucose_limit": args.glucose_limit,
        "oxygen_limit": args.oxygen_limit,
        "tyrosine_limit": args.tyrosine_limit,
        "atp_maintenance": args.atp_maintenance,
        "biomass_requirement": args.biomass,
    }
    if args.no_sbml:
        fba_kwargs["sbml_path"] = None
    fba = run_fba(**fba_kwargs)
    print(
        "       Status: {status}; DHN={dhn}, DOPA={dopa}".format(
            status=fba["status"],
            dhn=fba["pathway_fluxes"]["dhn_melanin_export"],
            dopa=fba["pathway_fluxes"]["dopa_melanin_export"],
        )
    )

    if args.sensitivity:
        print("\n[4/6] Analisis sensitivitas input radiasi...")
        payload = run_sensitivity(
            glucose_limit=args.glucose_limit,
            oxygen_limit=args.oxygen_limit,
            tyrosine_limit=args.tyrosine_limit,
            atp_maintenance=args.atp_maintenance,
            biomass_requirement=args.biomass,
        )
        print(f"       {len(payload['rows'])} kombinasi diuji.")
    else:
        print("\n[4/6] Analisis sensitivitas dilewati (--sensitivity untuk mengaktifkan).")

    print("\n[5/6] Membuat grafik fluks...")
    plot_fba_fluxes()

    print("\n[6/6] Menulis laporan dan README...")
    report_path = generate_report()
    readme_path = generate_readme()

    print("\n=== Pipeline selesai ===")
    print(f"Laporan : {report_path}")
    print(f"README  : {readme_path}")
    print(f"Hasil   : {PROJECT_ROOT / 'results'}")
    print(f"Model   : {PROJECT_ROOT / 'models'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
