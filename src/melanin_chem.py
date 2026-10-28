"""Deskriptor kimia untuk monomer/model melanin (DHN dan L-DOPA).

Model ini adalah representasi monomer, bukan polimer melanin utuh. Nilai
deskriptor dihitung dengan RDKit 2D dan muatan parsial Gasteiger.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors, rdPartialCharges

PROJECT_ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = PROJECT_ROOT / "results"
CHEM_RESULTS_PATH = RESULTS_DIR / "melanin_chem.json"

# SMILES diverifikasi terhadap PubChem:
# - DHN: 1,8-dihydroxynaphthalene, CID 68438, C10H8O2.
# - DOPA: L-DOPA, CID 6047, C9H11NO4.
# Catatan: versi lama memakai SMILES 1,5-dihydroxynaphthalene, bukan 1,8-DHN.
MONOMER_MODELS: dict[str, dict[str, Any]] = {
    "DHN_monomer_model": {
        "common_name": "1,8-dihydroxynaphthalene (DHN)",
        "pubchem_cid": 68438,
        "smiles": "Oc1cccc2cccc(O)c12",
        "notes": (
            "Monomer model untuk DHN-melanin. Bukan representasi polimer atau "
            "interaksi supramolekul."
        ),
    },
    "DOPA_monomer_model": {
        "common_name": "L-DOPA",
        "pubchem_cid": 6047,
        "smiles": "N[C@@H](Cc1ccc(O)c(O)c1)C(=O)O",
        "notes": (
            "Prekursor DOPA-melanin dalam bentuk L. Deskriptor tidak mencakup "
            "dopaquinone atau polimer melanin."
        ),
    },
}

ASSUMPTIONS: tuple[str, ...] = (
    "DHN direpresentasikan sebagai monomer 1,8-dihydroxynaphthalene, bukan polimer.",
    "DOPA direpresentasikan sebagai L-DOPA, prekursor DOPA-melanin.",
    "Muatan Gasteiger dihitung setelah penambahan hidrogen eksplisit.",
    "Deskriptor 2D tidak menangkap polimerisasi, agregasi, atau konformasi 3D.",
    "Semua nilai adalah estimasi in silico, bukan hasil pengukuran eksperimental.",
)


def characterize_smiles(smiles: str) -> dict[str, Any]:
    """Hitung deskriptor molekuler untuk satu SMILES."""
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError(f"SMILES tidak valid: {smiles}")

    molecule_with_hydrogens = Chem.AddHs(molecule)
    rdPartialCharges.ComputeGasteigerCharges(molecule_with_hydrogens)

    charges = []
    for atom in molecule_with_hydrogens.GetAtoms():
        charge = float(atom.GetProp("_GasteigerCharge"))
        if not math.isnan(charge) and not math.isinf(charge):
            charges.append(charge)

    return {
        "canonical_smiles": Chem.MolToSmiles(molecule),
        "formula": rdMolDescriptors.CalcMolFormula(molecule),
        "molecular_weight_g_mol": round(Descriptors.MolWt(molecule), 2),
        "logP": round(Crippen.MolLogP(molecule), 2),
        "tpsa_A2": round(rdMolDescriptors.CalcTPSA(molecule), 2),
        "h_bond_donors": int(rdMolDescriptors.CalcNumHBD(molecule)),
        "h_bond_acceptors": int(rdMolDescriptors.CalcNumHBA(molecule)),
        "rotatable_bonds": int(Descriptors.NumRotatableBonds(molecule)),
        "heavy_atoms": int(molecule.GetNumHeavyAtoms()),
        "gasteiger_charge_max": round(max(charges), 4) if charges else None,
        "gasteiger_charge_min": round(min(charges), 4) if charges else None,
    }


def characterize_monomer_models() -> dict[str, dict[str, Any]]:
    """Hitung deskriptor untuk semua monomer yang terdaftar."""
    results: dict[str, dict[str, Any]] = {}
    for name, metadata in MONOMER_MODELS.items():
        descriptors = characterize_smiles(metadata["smiles"])
        results[name] = {
            "common_name": metadata["common_name"],
            "pubchem_cid": metadata["pubchem_cid"],
            "source_smiles": metadata["smiles"],
            "notes": metadata["notes"],
            **descriptors,
        }
    return results


def save_results(output_path: Path | None = None) -> dict[str, Any]:
    """Hitung deskriptor dan simpan ke results/melanin_chem.json."""
    payload: dict[str, Any] = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "method": (
            "RDKit 2D descriptors; Gasteiger partial charges after explicit "
            "hydrogen addition. SMILES diverifikasi terhadap PubChem."
        ),
        "monomers": characterize_monomer_models(),
        "assumptions": list(ASSUMPTIONS),
    }

    destination = Path(output_path) if output_path else CHEM_RESULTS_PATH
    destination.parent.mkdir(parents=True, exist_ok=True)
    temp_path = destination.with_suffix(destination.suffix + ".tmp")
    temp_path.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    temp_path.replace(destination)
    return payload


def main() -> int:
    payload = save_results()
    print("=== Model Monomer Melanin ===")
    for model_name, properties in payload["monomers"].items():
        print(f"\n[{model_name}] {properties['common_name']}")
        for key in (
            "canonical_smiles",
            "formula",
            "molecular_weight_g_mol",
            "logP",
            "tpsa_A2",
            "h_bond_donors",
            "h_bond_acceptors",
            "rotatable_bonds",
            "gasteiger_charge_max",
            "gasteiger_charge_min",
        ):
            print(f"  - {key}: {properties[key]}")
    print(f"\n[SUCCESS] Hasil kimia disimpan ke: {CHEM_RESULTS_PATH}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
