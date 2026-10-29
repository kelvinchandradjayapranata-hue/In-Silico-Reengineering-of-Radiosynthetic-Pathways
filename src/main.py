import argparse
from pathlib import Path
from ncbi_fetcher import fetch_target
from melanin_chem import characterize_monomer_models
from metabolic_fba import run_fba
from visualization import plot_fba_fluxes
from generate_report import generate_report
from generate_readme import generate_readme

def main():
    parser = argparse.ArgumentParser(description='Melanin Radiosynthesis Bio-Computational Engine')
    parser.add_argument('--fetch', action='store_true', help='Fetch sample NCBI FASTA')
    parser.add_argument('--chem', action='store_true', help='Run RDKit melanin chemical descriptors')
    parser.add_argument('--fba', action='store_true', help='Run COBRApy metabolic FBA simulation')
    parser.add_argument('--visualize', action='store_true', help='Generate Matplotlib flux plots')
    parser.add_argument('--report', action='store_true', help='Generate research summary report')
    parser.add_argument('--doc', action='store_true', help='Generate README documentation')
    parser.add_argument('--all', action='store_true', help='Run full pipeline end-to-end')

    args = parser.parse_args()

    if not any(vars(args).values()):
        args.all = True

    print('==================================================')
    print('  MELANIN RADIOSYNTHESIS BIO-COMPUTATIONAL ENGINE ')
    print('==================================================')

    if args.all or args.fetch:
        print('\n[1/5] Fetching NCBI genomic data...')
        try:
            fetch_target('pks1')
        except Exception as e:
            print(f'  - Fetch skipped/note: {e}')

    if args.all or args.chem:
        print('\n[2/5] Calculating melanin chemical descriptors...')
        chem_res = characterize_monomer_models()
        print('  - DHN & DOPA models processed successfully.')

    if args.all or args.fba:
        print('\n[3/5] Running metabolic Flux Balance Analysis...')
        fba_res = run_fba()
        status_val = fba_res['status']
        print(f'  - FBA Status: {status_val}')

    if args.all or args.visualize:
        print('\n[4/5] Generating visualization plots...')
        plot_fba_fluxes()

    if args.all or args.report or args.doc:
        print('\n[5/5] Compiling report and documentation...')
        rep_path = generate_report()
        doc_path = generate_readme()
        print(f'  - Report: {rep_path.name}')
        print(f'  - Documentation: {doc_path.name}')

    print('\n[SUCCESS] Pipeline execution finished successfully!')

if __name__ == '__main__':
    main()