"""Peta interaksi ligan-reseptor dari pose AutoDock Vina.

Modul ini membaca pose terbaik dari ``results/docking_results.json`` (format
single-seed maupun agregasi multi-seed), mengonversi pose PDBQT ke RDKit
melalui Meeko, lalu memetakan:

- ikatan hidrogen (H-bond) dengan kriteria jarak dan sudut,
- kontak hidrofobik atom karbon ligan dengan karbon side chain protein,
- stacking pi-pi antara ring aromatik ligan dan residu Phe/Tyr/Trp/His.

Analisis bersifat geometris dan hanya untuk triase. Hasilnya bukan energi
ikat, bukan bukti eksperimental, dan tidak menggantikan analisis MD/MM-PBSA
atau alat seperti PLIP/Arpeggio. Reseptor default adalah struktur homolog
surrogate, bukan struktur target *C. cladosporioides*.
"""

from __future__ import annotations

import argparse
import json
import math
import re
from dataclasses import dataclass
from datetime import datetime, timezone
from io import StringIO
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np
from Bio.PDB import NeighborSearch, PDBIO, PDBParser
from meeko import PDBQTMolecule, RDKitMolCreate
from rdkit import Chem

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DOCKING = PROJECT_ROOT / "results" / "docking_results.json"
DEFAULT_JSON = PROJECT_ROOT / "results" / "interactions.json"
DEFAULT_TSV = PROJECT_ROOT / "results" / "interactions.tsv"
DEFAULT_MD = PROJECT_ROOT / "results" / "interactions.md"
POSE_DIR = PROJECT_ROOT / "results" / "docking" / "poses"

POLAR_ELEMENTS = {"N", "O", "S"}
HYDROPHOBIC_ELEMENTS = {"C"}

# Ring aromatik protein yang dipakai untuk klasifikasi pi-pi.
PROTEIN_RINGS: dict[str, list[tuple[str, tuple[str, ...]]]] = {
    "PHE": [("ring", ("CG", "CD1", "CD2", "CE1", "CE2", "CZ"))],
    "TYR": [("ring", ("CG", "CD1", "CD2", "CE1", "CE2", "CZ"))],
    "HIS": [("ring", ("CG", "ND1", "CD2", "CE1", "NE2"))],
    "TRP": [
        ("pyrrole", ("CG", "CD1", "NE1", "CE2", "CD2")),
        ("benzene", ("CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2")),
    ],
}


class AnalysisError(RuntimeError):
    """Kesalahan analisis yang dapat ditampilkan tanpa traceback."""


@dataclass
class ReceptorContext:
    name: str
    path: Path
    chain_id: str
    residues: list[Any]
    heavy_atoms: list[Any]
    all_atoms: list[Any]
    neighbor_search: NeighborSearch
    hydrogen_search: NeighborSearch | None


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _relative(path: Path) -> str:
    try:
        return path.resolve().relative_to(PROJECT_ROOT.resolve()).as_posix()
    except ValueError:
        return str(path)


def resolve_path(value: str | Path | None, *, must_exist: bool = True) -> Path:
    """Resolve path relatif terhadap root proyek dan cari fallback struktur."""

    if value is None or str(value).strip() == "":
        raise AnalysisError("Path tidak diberikan.")
    text = str(value).replace("<PROJECT_ROOT>", str(PROJECT_ROOT))
    raw = Path(text)
    candidates: list[Path] = []
    if raw.is_absolute():
        candidates.append(raw)
    else:
        candidates.append(PROJECT_ROOT / raw)
    if raw.name:
        candidates.extend(
            [
                PROJECT_ROOT / "data" / "structures" / raw.name,
                PROJECT_ROOT / "data" / "structures" / "prepared" / raw.name,
            ]
        )
    for candidate in candidates:
        if candidate.exists():
            return candidate
    if must_exist:
        raise AnalysisError(f"File tidak ditemukan: {value}")
    return candidates[0]


def _read_json(path: Path) -> dict[str, Any]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise AnalysisError(f"Gagal membaca {path}: {exc}") from exc
    if not isinstance(payload, dict):
        raise AnalysisError(f"Isi {path} bukan objek JSON.")
    return payload


def _write_json_atomic(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp_path = path.with_suffix(path.suffix + ".tmp")
    temp_path.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    temp_path.replace(path)


def _pdbqt_first_affinity(path: Path) -> float | None:
    try:
        text = path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return None
    match = re.search(r"REMARK VINA RESULT:\s+(-?\d+(?:\.\d+)?)", text)
    return float(match.group(1)) if match else None


def _candidate_pose_paths(
    row: dict[str, Any], receptor: str, ligand: str
) -> list[Path]:
    candidates: list[Path] = []
    per_seed = row.get("per_seed") or []
    for value in [row.get("pose_file")] + [
        item.get("pose_file")
        for item in per_seed
        if isinstance(item, dict)
    ]:
        if not value:
            continue
        try:
            candidates.append(resolve_path(value))
        except AnalysisError:
            continue

    patterns = [
        POSE_DIR / f"{receptor}__{ligand}__seed*.pdbqt",
        POSE_DIR / f"{receptor}__{ligand}.pdbqt",
    ]
    for pattern in patterns:
        candidates.extend(sorted(POSE_DIR.glob(pattern.name)))

    unique: list[Path] = []
    seen: set[str] = set()
    for candidate in candidates:
        key = str(candidate.resolve())
        if key not in seen and candidate.is_file():
            seen.add(key)
            unique.append(candidate)
    return unique


def find_pose_path(row: dict[str, Any]) -> Path:
    receptor = str(row.get("receptor", ""))
    ligand = str(row.get("ligand", ""))

    # Prioritaskan pose yang dirujuk langsung oleh JSON. Pada agregasi
    # multi-seed, field ini adalah pose terbaik yang sudah dipilih agregator;
    # jangan biarkan pose pilot lama di folder yang sama mengalahkannya.
    row_pose = row.get("pose_file")
    if row_pose:
        try:
            path = resolve_path(row_pose)
            if path.is_file():
                return path
        except AnalysisError:
            pass

    per_seed_paths: list[Path] = []
    for item in row.get("per_seed") or []:
        if not isinstance(item, dict) or not item.get("pose_file"):
            continue
        try:
            path = resolve_path(item["pose_file"])
        except AnalysisError:
            continue
        if path.is_file():
            per_seed_paths.append(path)
    if per_seed_paths:
        ranked = [
            (path, _pdbqt_first_affinity(path))
            for path in per_seed_paths
        ]
        ranked = [(path, value) for path, value in ranked if value is not None]
        if ranked:
            return min(ranked, key=lambda item: item[1])[0]
        return per_seed_paths[0]

    candidates = _candidate_pose_paths(row, receptor, ligand)
    if not candidates:
        raise AnalysisError(
            f"Pose PDBQT tidak ditemukan untuk {receptor} x {ligand}. "
            "Jalankan docking atau unduh results/docking/poses/."
        )
    if len(candidates) == 1:
        return candidates[0]
    ranked = [
        (path, _pdbqt_first_affinity(path))
        for path in candidates
    ]
    ranked = [(path, value) for path, value in ranked if value is not None]
    if ranked:
        return min(ranked, key=lambda item: item[1])[0]
    return candidates[0]


def _pose_seed(row: dict[str, Any], pose_path: Path) -> Any:
    """Seed yang menghasilkan pose terpilih (bila metadata tersedia)."""

    target = pose_path.resolve()
    matches: list[dict[str, Any]] = []
    for item in row.get("per_seed") or []:
        if not isinstance(item, dict) or not item.get("pose_file"):
            continue
        try:
            candidate = resolve_path(item["pose_file"])
        except AnalysisError:
            continue
        if candidate.resolve() == target:
            matches.append(item)
    if matches:
        def _affinity(item: dict[str, Any]) -> float:
            value = item.get(
                "best_affinity_kcal_mol",
                item.get("affinity_kcal_mol", item.get("mean_affinity_kcal_mol")),
            )
            return float(value) if value is not None else float("inf")

        best = min(matches, key=_affinity)
        return best.get("seed")

    target_affinity = row.get("best_affinity_kcal_mol")
    if target_affinity is not None:
        for item in row.get("per_seed") or []:
            if not isinstance(item, dict):
                continue
            for key in (
                "best_affinity_kcal_mol",
                "affinity_kcal_mol",
                "mean_affinity_kcal_mol",
            ):
                value = item.get(key)
                if value is None:
                    continue
                try:
                    if abs(float(value) - float(target_affinity)) < 1e-6:
                        return item.get("seed")
                except (TypeError, ValueError):
                    continue

    if row.get("best_seed") is not None:
        return row.get("best_seed")
    seeds = row.get("seeds") or []
    return seeds[0] if len(seeds) == 1 else row.get("seed")


def load_ligand_mol(pose_path: Path) -> Chem.Mol:
    """Baca pose PDBQT dan kembalikan RDKit mol dengan semua pose sebagai conformer."""

    try:
        pdbqt = PDBQTMolecule.from_file(str(pose_path), skip_typing=True)
        mols = RDKitMolCreate.from_pdbqt_mol(pdbqt)
    except Exception as exc:
        raise AnalysisError(
            f"Meeko gagal mengonversi {pose_path.name}: {exc}"
        ) from exc
    if not mols:
        raise AnalysisError(f"Tidak ada molekul ligan dari {pose_path.name}.")
    mol = mols[0]
    if mol.GetNumConformers() == 0:
        raise AnalysisError(f"Pose {pose_path.name} tidak memiliki conformer.")
    return mol


def _protein_hydrogens(ctx: ReceptorContext, atom: Any) -> list[Any]:
    if ctx.hydrogen_search is None:
        return []
    return list(ctx.hydrogen_search.search(atom.coord, 1.25))


def _protein_roles(ctx: ReceptorContext, atom: Any) -> set[str]:
    roles: set[str] = set()
    element = atom.element.upper()
    name = atom.get_name()
    resname = atom.get_parent().get_resname()

    if element == "O":
        roles.add("acceptor")
        if name in {"OG", "OG1", "OH"}:
            roles.add("donor")
    elif element == "N":
        if name == "N" and resname == "PRO":
            pass
        else:
            roles.add("donor")
        if resname == "HIS" and name in {"ND1", "NE2"}:
            roles.add("acceptor")
    elif element == "S":
        roles.update({"donor", "acceptor"})

    if _protein_hydrogens(ctx, atom):
        roles.add("donor")
    return roles


def _ligand_roles(atom: Any) -> set[str]:
    roles: set[str] = set()
    symbol = atom.GetSymbol()
    if symbol not in POLAR_ELEMENTS:
        return roles

    has_h = any(neighbor.GetSymbol() == "H" for neighbor in atom.GetNeighbors())
    if has_h:
        roles.add("donor")
    if symbol == "O":
        roles.add("acceptor")
    elif symbol == "N":
        # Amida/amina protonated tetap dapat menjadi donor; akseptor hanya
        # untuk N netral yang bukan amida sederhana.
        if atom.GetFormalCharge() <= 0 and not _is_amide_nitrogen(atom):
            roles.add("acceptor")
    elif symbol == "S":
        roles.update({"donor", "acceptor"})
    return roles


def _is_amide_nitrogen(atom: Any) -> bool:
    for neighbor in atom.GetNeighbors():
        if neighbor.GetSymbol() != "C":
            continue
        for bond in neighbor.GetBonds():
            other = bond.GetOtherAtom(neighbor)
            if bond.GetBondType() == Chem.BondType.DOUBLE and other.GetSymbol() == "O":
                return True
    return False


def _angle_degrees(
    donor: Sequence[float], hydrogen: Sequence[float], acceptor: Sequence[float]
) -> float:
    vec_dh = np.asarray(donor, dtype=float) - np.asarray(hydrogen, dtype=float)
    vec_ah = np.asarray(acceptor, dtype=float) - np.asarray(hydrogen, dtype=float)
    denominator = float(np.linalg.norm(vec_dh) * np.linalg.norm(vec_ah))
    if denominator == 0.0:
        return 0.0
    cosine = float(np.dot(vec_dh, vec_ah) / denominator)
    return float(math.degrees(math.acos(max(-1.0, min(1.0, cosine)))))


def _ligand_atom_label(mol: Chem.Mol, index: int) -> str:
    atom = mol.GetAtomWithIdx(index)
    try:
        tripos = atom.GetProp("_TriposAtomName")
    except KeyError:
        tripos = ""
    return f"{atom.GetSymbol()}{index}" + (f"({tripos})" if tripos else "")


def _residue_label(atom: Any) -> str:
    residue = atom.get_parent()
    insertion = str(residue.id[2]).strip()
    return f"{residue.get_resname()}{residue.id[1]}{insertion}"


def detect_hbonds(
    mol: Chem.Mol,
    receptor: ReceptorContext,
    *,
    cutoff: float = 3.5,
    min_angle: float = 120.0,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Kembalikan (H-bonds, polar contacts non-H-bond)."""

    conformer = mol.GetConformer(0)
    coordinates = conformer.GetPositions()
    hbonds: dict[tuple[Any, ...], dict[str, Any]] = {}
    polar: dict[tuple[Any, ...], dict[str, Any]] = {}

    for ligand_index, atom in enumerate(mol.GetAtoms()):
        if atom.GetSymbol() not in POLAR_ELEMENTS:
            continue
        roles = _ligand_roles(atom)
        ligand_position = coordinates[ligand_index].tolist()
        for protein_atom in receptor.neighbor_search.search(
            coordinates[ligand_index], cutoff
        ):
            if protein_atom.element not in POLAR_ELEMENTS:
                continue
            distance = float(
                np.linalg.norm(
                    coordinates[ligand_index] - np.asarray(protein_atom.coord)
                )
            )
            protein_roles = _protein_roles(receptor, protein_atom)
            candidates: list[tuple[float, str, str]] = []

            if "donor" in roles and "acceptor" in protein_roles:
                for hydrogen in [
                    neighbor
                    for neighbor in atom.GetNeighbors()
                    if neighbor.GetSymbol() == "H"
                ]:
                    angle = _angle_degrees(
                        coordinates[ligand_index],
                        coordinates[hydrogen.GetIdx()],
                        protein_atom.coord,
                    )
                    candidates.append(
                        (
                            angle,
                            "ligand_donor",
                            _ligand_atom_label(mol, hydrogen.GetIdx()),
                        )
                    )

            if "acceptor" in roles and "donor" in protein_roles:
                for hydrogen in _protein_hydrogens(receptor, protein_atom):
                    angle = _angle_degrees(
                        protein_atom.coord,
                        hydrogen.coord,
                        coordinates[ligand_index],
                    )
                    candidates.append(
                        (
                            angle,
                            "protein_donor",
                            f"{protein_atom.get_name()}-H",
                        )
                    )

            role_pair = (
                "donor" in roles and "acceptor" in protein_roles
            ) or (
                "acceptor" in roles and "donor" in protein_roles
            )
            key = (
                ligand_index,
                protein_atom.get_parent().get_resname(),
                protein_atom.get_parent().id[1],
                str(protein_atom.get_parent().id[2]).strip(),
                protein_atom.get_name(),
            )
            base = {
                "ligand_atom": _ligand_atom_label(mol, ligand_index),
                "ligand_element": atom.GetSymbol(),
                "ligand_position": [round(value, 3) for value in ligand_position],
                "protein_residue": _residue_label(protein_atom),
                "protein_chain": receptor.chain_id,
                "protein_resname": protein_atom.get_parent().get_resname(),
                "protein_resi": int(protein_atom.get_parent().id[1]),
                "protein_icode": str(protein_atom.get_parent().id[2]).strip(),
                "protein_atom": protein_atom.get_name(),
                "protein_element": protein_atom.element,
                "protein_position": [
                    round(float(value), 3) for value in protein_atom.coord
                ],
                "distance_angstrom": round(distance, 3),
            }

            if not candidates:
                if role_pair:
                    direction = (
                        "protein_donor"
                        if "acceptor" in roles and "donor" in protein_roles
                        else "ligand_donor"
                    )
                    polar[key] = {
                        **base,
                        "direction": direction,
                        "hydrogen": None,
                        "angle_deg": None,
                        "geometry_ok": None,
                        "confidence": "distance_only",
                        "reason": (
                            "pasangan donor-akseptor cocok, tetapi atom H tidak "
                            "tersedia untuk menghitung sudut"
                        ),
                    }
                else:
                    polar[key] = {
                        **base,
                        "direction": None,
                        "hydrogen": None,
                        "angle_deg": None,
                        "geometry_ok": None,
                        "confidence": "none",
                        "reason": "tidak ada pasangan donor-akseptor yang cocok",
                    }
                continue

            angle, direction, hydrogen_label = max(
                candidates, key=lambda item: item[0]
            )
            record = {
                **base,
                "direction": direction,
                "hydrogen": hydrogen_label,
                "angle_deg": round(angle, 1),
                "geometry_ok": angle >= min_angle,
                "confidence": "strong" if angle >= min_angle else "weak",
            }
            if angle >= min_angle:
                hbonds[key] = record
            else:
                polar[key] = {
                    **record,
                    "reason": "sudut di bawah ambang H-bond",
                }

    return (
        sorted(
            hbonds.values(),
            key=lambda item: (item["distance_angstrom"], item["protein_residue"]),
        ),
        sorted(
            polar.values(),
            key=lambda item: (item["distance_angstrom"], item["protein_residue"]),
        ),
    )


def detect_hydrophobic_contacts(
    mol: Chem.Mol,
    receptor: ReceptorContext,
    *,
    cutoff: float = 4.0,
) -> list[dict[str, Any]]:
    conformer = mol.GetConformer(0)
    coordinates = conformer.GetPositions()
    groups: dict[tuple[Any, ...], dict[str, Any]] = {}

    for ligand_index, atom in enumerate(mol.GetAtoms()):
        if atom.GetSymbol() not in HYDROPHOBIC_ELEMENTS:
            continue
        for protein_atom in receptor.neighbor_search.search(
            coordinates[ligand_index], cutoff
        ):
            if protein_atom.element not in HYDROPHOBIC_ELEMENTS:
                continue
            # Hanya karbon side chain; kerangka N/CA/C/O bukan kontak
            # hidrofobik yang interpretable.
            if protein_atom.get_name() in {"N", "CA", "C", "O", "OXT"}:
                continue
            distance = float(
                np.linalg.norm(
                    coordinates[ligand_index] - np.asarray(protein_atom.coord)
                )
            )
            residue = protein_atom.get_parent()
            key = (
                residue.get_resname(),
                residue.id[1],
                str(residue.id[2]).strip(),
            )
            group = groups.setdefault(
                key,
                {
                    "protein_residue": _residue_label(protein_atom),
                    "protein_chain": receptor.chain_id,
                    "protein_resname": residue.get_resname(),
                    "protein_resi": int(residue.id[1]),
                    "protein_icode": str(residue.id[2]).strip(),
                    "contact_count": 0,
                    "ligand_atoms": set(),
                    "protein_atoms": set(),
                    "min_distance_angstrom": distance,
                    "representative_ligand_position": [
                        round(float(value), 3)
                        for value in coordinates[ligand_index]
                    ],
                    "representative_protein_position": [
                        round(float(value), 3) for value in protein_atom.coord
                    ],
                },
            )
            group["contact_count"] += 1
            group["ligand_atoms"].add(_ligand_atom_label(mol, ligand_index))
            group["protein_atoms"].add(protein_atom.get_name())
            if distance < group["min_distance_angstrom"]:
                group["min_distance_angstrom"] = distance
                group["representative_ligand_position"] = [
                    round(float(value), 3) for value in coordinates[ligand_index]
                ]
                group["representative_protein_position"] = [
                    round(float(value), 3) for value in protein_atom.coord
                ]

    result = []
    for group in groups.values():
        group["min_distance_angstrom"] = round(
            group["min_distance_angstrom"], 3
        )
        group["ligand_atoms"] = sorted(group["ligand_atoms"])
        group["protein_atoms"] = sorted(group["protein_atoms"])
        result.append(group)
    return sorted(
        result,
        key=lambda item: (
            item["min_distance_angstrom"],
            item["protein_residue"],
        ),
    )


def _plane_geometry(points: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    centroid = points.mean(axis=0)
    centered = points - centroid
    _, _, vh = np.linalg.svd(centered)
    normal = vh[-1]
    norm = np.linalg.norm(normal)
    if norm == 0:
        raise ValueError("Normal ring tidak valid.")
    return centroid, normal / norm


def _ring_angle(normal_a: np.ndarray, normal_b: np.ndarray) -> float:
    cosine = abs(float(np.dot(normal_a, normal_b)))
    return float(math.degrees(math.acos(max(-1.0, min(1.0, cosine)))))


def detect_pi_stacking(
    mol: Chem.Mol,
    receptor: ReceptorContext,
    *,
    cutoff: float = 5.5,
    max_centroid_offset: float = 6.5,
    max_atom_distance: float = 4.0,
) -> list[dict[str, Any]]:
    conformer = mol.GetConformer(0)
    coordinates = conformer.GetPositions()
    ligand_rings: list[dict[str, Any]] = []
    for ring_index, ring in enumerate(mol.GetRingInfo().AtomRings(), start=1):
        if len(ring) not in {5, 6}:
            continue
        if not all(mol.GetAtomWithIdx(index).GetIsAromatic() for index in ring):
            continue
        ring_coordinates = np.asarray([coordinates[index] for index in ring])
        try:
            centroid, normal = _plane_geometry(ring_coordinates)
        except ValueError:
            continue
        ligand_rings.append(
            {
                "label": f"ligand_ring_{ring_index}",
                "atoms": [_ligand_atom_label(mol, index) for index in ring],
                "coordinates": ring_coordinates,
                "centroid": centroid,
                "normal": normal,
            }
        )
    if not ligand_rings:
        return []

    protein_rings: list[dict[str, Any]] = []
    for residue in receptor.residues:
        definitions = PROTEIN_RINGS.get(residue.get_resname())
        if not definitions:
            continue
        for ring_name, atom_names in definitions:
            atoms = []
            for atom_name in atom_names:
                if not residue.has_id(atom_name):
                    atoms = []
                    break
                atoms.append(residue[atom_name])
            if len(atoms) != len(atom_names):
                continue
            ring_coordinates = np.asarray([atom.coord for atom in atoms])
            try:
                centroid, normal = _plane_geometry(ring_coordinates)
            except ValueError:
                continue
            protein_rings.append(
                {
                    "label": f"{_residue_label(atoms[0])}:{ring_name}",
                    "residue": _residue_label(atoms[0]),
                    "resname": atoms[0].get_parent().get_resname(),
                    "resi": int(atoms[0].get_parent().id[1]),
                    "icode": str(atoms[0].get_parent().id[2]).strip(),
                    "atoms": [atom.get_name() for atom in atoms],
                    "coordinates": ring_coordinates,
                    "centroid": centroid,
                    "normal": normal,
                }
            )

    results = []
    for ligand_ring in ligand_rings:
        for protein_ring in protein_rings:
            distance = float(
                np.linalg.norm(
                    ligand_ring["centroid"] - protein_ring["centroid"]
                )
            )
            if distance > max_centroid_offset:
                continue
            angle = _ring_angle(
                ligand_ring["normal"], protein_ring["normal"]
            )
            if angle <= 30.0:
                stacking_type = "parallel"
            elif angle >= 60.0:
                stacking_type = "T-shaped"
            else:
                # Sudut 30-60 derajat tidak cukup jelas untuk diklaim sebagai
                # stacking; biarkan analisis lain (PLIP) yang memutuskan.
                continue

            # Jarak centroid saja melewatkan parallel-displaced stacking.
            # Tambahkan kriteria jarak atom-terdekat untuk kasus offset.
            min_atom_distance = min(
                float(np.linalg.norm(ligand_atom - protein_atom))
                for ligand_atom in ligand_ring["coordinates"]
                for protein_atom in protein_ring["coordinates"]
            )
            is_direct = distance <= cutoff
            is_offset = (
                distance <= max_centroid_offset
                and min_atom_distance <= max_atom_distance
            )
            if not (is_direct or is_offset):
                continue

            results.append(
                {
                    "ligand_ring": ligand_ring["label"],
                    "ligand_ring_atoms": ligand_ring["atoms"],
                    "protein_ring": protein_ring["label"],
                    "protein_residue": protein_ring["residue"],
                    "protein_resname": protein_ring["resname"],
                    "protein_resi": protein_ring["resi"],
                    "protein_icode": protein_ring["icode"],
                    "protein_ring_atoms": protein_ring["atoms"],
                    "centroid_distance_angstrom": round(distance, 3),
                    "min_atom_distance_angstrom": round(min_atom_distance, 3),
                    "angle_deg": round(angle, 1),
                    "stacking_type": stacking_type,
                    "geometry": "direct" if is_direct else "offset",
                    "ligand_centroid": [
                        round(float(value), 3)
                        for value in ligand_ring["centroid"]
                    ],
                    "protein_centroid": [
                        round(float(value), 3)
                        for value in protein_ring["centroid"]
                    ],
                }
            )
    return sorted(
        results,
        key=lambda item: (
            item["centroid_distance_angstrom"],
            item["protein_ring"],
        ),
    )


def load_receptor_context(
    path: Path, *, chain_id: str | None = None, name: str = "receptor"
) -> ReceptorContext:
    try:
        structure = PDBParser(QUIET=True).get_structure(name, str(path))
        model = next(structure.get_models())
    except Exception as exc:
        raise AnalysisError(f"Gagal membaca reseptor {path}: {exc}") from exc

    chains = list(model.get_chains())
    if not chains:
        raise AnalysisError(f"Reseptor {path} tidak memiliki chain.")
    chain = None
    if chain_id:
        chain = next((item for item in chains if item.id == chain_id), None)
    if chain is None:
        chain = next(
            (
                item
                for item in chains
                if any(residue.id[0] == " " for residue in item)
            ),
            chains[0],
        )

    residues = [residue for residue in chain if residue.id[0] == " "]
    all_atoms = [
        atom
        for residue in residues
        for atom in residue
        if atom.get_altloc() in {" ", "A"}
    ]
    heavy_atoms = [atom for atom in all_atoms if atom.element != "H"]
    if not heavy_atoms:
        raise AnalysisError(f"Reseptor {path} tidak memiliki atom berat.")
    hydrogens = [atom for atom in all_atoms if atom.element == "H"]
    return ReceptorContext(
        name=name,
        path=path,
        chain_id=chain.id,
        residues=residues,
        heavy_atoms=heavy_atoms,
        all_atoms=all_atoms,
        neighbor_search=NeighborSearch(heavy_atoms),
        hydrogen_search=NeighborSearch(hydrogens) if hydrogens else None,
    )


def receptor_pdb_text(
    path: Path, *, chain_id: str | None = None
) -> str:
    """Ambil teks PDB reseptor tanpa hidrogen untuk viewer 3D."""

    structure = PDBParser(QUIET=True).get_structure("receptor", str(path))
    model = next(structure.get_models())
    if chain_id and chain_id not in {chain.id for chain in model}:
        chain_id = None

    class _Select:
        def accept_model(self, model: Any) -> bool:
            return True

        def accept_chain(self, chain: Any) -> bool:
            return chain_id is None or chain.id == chain_id

        def accept_residue(self, residue: Any) -> bool:
            return residue.id[0] == " "

        def accept_atom(self, atom: Any) -> bool:
            return (
                atom.get_altloc() in {" ", "A"}
                and atom.element != "H"
            )

    output = StringIO()
    io = PDBIO()
    io.set_structure(model)
    io.save(output, _Select())
    return output.getvalue()


def _receptor_path(metadata: dict[str, Any]) -> Path:
    for key in ("clean_pdb", "pdb_path"):
        value = metadata.get(key)
        if not value:
            continue
        try:
            return resolve_path(value)
        except AnalysisError:
            continue

    name = str(metadata.get("name", "")).strip()
    pdb_id = str(metadata.get("pdb_id", "")).strip()
    fallbacks = []
    if name:
        fallbacks.extend(
            [
                PROJECT_ROOT / "data" / "structures" / "prepared" / f"{name}_clean.pdb",
                PROJECT_ROOT / "data" / "structures" / "prepared" / f"{name}.pdb",
                PROJECT_ROOT / "data" / "structures" / f"{name}.pdb",
            ]
        )
    if pdb_id:
        fallbacks.append(PROJECT_ROOT / "data" / "structures" / f"{pdb_id}.pdb")
    for candidate in fallbacks:
        if candidate.is_file():
            return candidate
    raise AnalysisError(
        f"File reseptor untuk {metadata.get('name', '?')} tidak ditemukan."
    )


def _select_rows(
    payload: dict[str, Any],
    receptor_names: Sequence[str] | None,
    ligand_names: Sequence[str] | None,
) -> list[dict[str, Any]]:
    rows = [
        row
        for row in (payload.get("results") or [])
        if row.get("status") in {"ok", "completed"}
        or (
            row.get("status") is None
            and row.get("best_affinity_kcal_mol") is not None
        )
    ]
    if receptor_names and "all" not in {name.lower() for name in receptor_names}:
        selected = set(receptor_names)
        rows = [row for row in rows if row.get("receptor") in selected]
    if ligand_names and "all" not in {name.lower() for name in ligand_names}:
        selected = set(ligand_names)
        rows = [row for row in rows if row.get("ligand") in selected]
    if not rows:
        raise AnalysisError(
            "Tidak ada kombinasi reseptor x ligan yang berhasil dan cocok "
            "dengan filter."
        )
    return rows


def _receptor_metadata(
    payload: dict[str, Any], name: str
) -> dict[str, Any]:
    for receptor in payload.get("receptors") or []:
        if receptor.get("name") == name:
            return dict(receptor)
    return {"name": name}


def analyze_complex(
    row: dict[str, Any],
    receptor_metadata: dict[str, Any],
    *,
    chain_override: str | None = None,
    cutoff_hbond: float = 3.5,
    min_hbond_angle: float = 120.0,
    cutoff_hydrophobic: float = 4.0,
    cutoff_pipi: float = 5.5,
    max_pipi_centroid_offset: float = 6.5,
    max_pipi_atom_distance: float = 4.0,
) -> dict[str, Any]:
    receptor_name = str(row.get("receptor", ""))
    ligand_name = str(row.get("ligand", ""))
    pose_path = find_pose_path(row)
    ligand_mol = load_ligand_mol(pose_path)
    receptor_path = _receptor_path(receptor_metadata)
    chain_id = chain_override or receptor_metadata.get("chain") or "A"
    receptor = load_receptor_context(
        receptor_path, chain_id=chain_id, name=receptor_name
    )

    hbonds, polar_contacts = detect_hbonds(
        ligand_mol,
        receptor,
        cutoff=cutoff_hbond,
        min_angle=min_hbond_angle,
    )
    hydrophobic = detect_hydrophobic_contacts(
        ligand_mol, receptor, cutoff=cutoff_hydrophobic
    )
    pi_stacking = detect_pi_stacking(
        ligand_mol,
        receptor,
        cutoff=cutoff_pipi,
        max_centroid_offset=max_pipi_centroid_offset,
        max_atom_distance=max_pipi_atom_distance,
    )

    heavy_mol = Chem.RemoveHs(ligand_mol)
    affinity = row.get("best_affinity_kcal_mol")
    if affinity is None:
        affinity = _pdbqt_first_affinity(pose_path)
    heavy_atoms = heavy_mol.GetNumAtoms()
    ligand_efficiency = row.get("ligand_efficiency")
    if ligand_efficiency is None and affinity is not None and heavy_atoms:
        ligand_efficiency = round(float(affinity) / heavy_atoms, 3)
    seeds = row.get("seeds") or []
    seed = _pose_seed(row, pose_path)
    if seed is None:
        seed = seeds[0] if seeds else row.get("seed")

    return {
        "status": "ok",
        "receptor": receptor_name,
        "ligand": ligand_name,
        "pdb_id": receptor_metadata.get("pdb_id"),
        "target_accession": receptor_metadata.get("target_accession"),
        "identity_percent": receptor_metadata.get("identity_percent"),
        "coverage_percent": receptor_metadata.get("coverage_percent"),
        "surrogate": bool(receptor_metadata.get("surrogate", True)),
        "seed": seed,
        "seeds": seeds,
        "affinity_kcal_mol": affinity,
        "ligand_efficiency": ligand_efficiency,
        "heavy_atoms": heavy_atoms,
        "ligand_smiles": Chem.MolToSmiles(heavy_mol),
        "pose_file": _relative(pose_path),
        "receptor_file": _relative(receptor_path),
        "receptor_chain": receptor.chain_id,
        "hbonds": hbonds,
        "polar_contacts": polar_contacts,
        "hydrophobic_contacts": hydrophobic,
        "pi_stacking": pi_stacking,
        "summary": {
            "hbonds": len(hbonds),
            "polar_contacts": len(polar_contacts),
            "hydrophobic_residues": len(hydrophobic),
            "hydrophobic_atom_pairs": sum(
                item["contact_count"] for item in hydrophobic
            ),
            "pi_stacking": len(pi_stacking),
        },
    }


def build_payload(
    docking_path: Path = DEFAULT_DOCKING,
    *,
    receptors: Sequence[str] | None = None,
    ligands: Sequence[str] | None = None,
    chain: str | None = None,
    cutoff_hbond: float = 3.5,
    min_hbond_angle: float = 120.0,
    cutoff_hydrophobic: float = 4.0,
    cutoff_pipi: float = 5.5,
    max_pipi_centroid_offset: float = 6.5,
    max_pipi_atom_distance: float = 4.0,
) -> dict[str, Any]:
    payload = _read_json(docking_path)
    rows = _select_rows(payload, receptors, ligands)
    complexes: list[dict[str, Any]] = []
    errors: list[dict[str, Any]] = []

    fallback_seed = (payload.get("parameters") or {}).get("seed")
    for row in rows:
        row = dict(row)
        if row.get("seed") is None and not row.get("seeds"):
            row["seed"] = fallback_seed
        metadata = _receptor_metadata(payload, str(row.get("receptor", "")))
        try:
            complexes.append(
                analyze_complex(
                    row,
                    metadata,
                    chain_override=chain,
                    cutoff_hbond=cutoff_hbond,
                    min_hbond_angle=min_hbond_angle,
                    cutoff_hydrophobic=cutoff_hydrophobic,
                    cutoff_pipi=cutoff_pipi,
                    max_pipi_centroid_offset=max_pipi_centroid_offset,
                    max_pipi_atom_distance=max_pipi_atom_distance,
                )
            )
        except (AnalysisError, OSError, ValueError, KeyError, TypeError) as exc:
            errors.append(
                {
                    "receptor": row.get("receptor"),
                    "ligand": row.get("ligand"),
                    "error": str(exc),
                }
            )

    if not complexes:
        status = "failed"
    elif errors:
        status = "partial"
    else:
        status = "completed"

    aggregation = (payload.get("summary") or {}).get("aggregation")
    return {
        "generated_at": _now(),
        "tool": "melanin_biotech_interaction_analysis",
        "status": status,
        "docking_file": _relative(docking_path),
        "docking_status": payload.get("status"),
        "docking_aggregation": aggregation,
        "docking_mode": (
            "multi_seed"
            if isinstance(aggregation, dict)
            and aggregation.get("method") == "multi_seed"
            else "single_seed"
        ),
        "parameters": {
            "cutoff_hbond_angstrom": cutoff_hbond,
            "min_hbond_angle_deg": min_hbond_angle,
            "cutoff_hydrophobic_angstrom": cutoff_hydrophobic,
            "cutoff_pipi_angstrom": cutoff_pipi,
            "max_pipi_centroid_offset_angstrom": max_pipi_centroid_offset,
            "max_pipi_atom_distance_angstrom": max_pipi_atom_distance,
            "receptor_filter": list(receptors or []),
            "ligand_filter": list(ligands or []),
            "chain_override": chain,
        },
        "complexes": complexes,
        "errors": errors,
        "summary": {
            "complexes_total": len(complexes),
            "complexes_failed": len(errors),
            "hbonds_total": sum(
                item["summary"]["hbonds"] for item in complexes
            ),
            "hydrophobic_residues_total": sum(
                item["summary"]["hydrophobic_residues"] for item in complexes
            ),
            "pi_stacking_total": sum(
                item["summary"]["pi_stacking"] for item in complexes
            ),
        },
        "caveats": [
            "H-bond ditentukan secara geometris dari pose rigid; sudut hanya "
            "dihitung bila atom hidrogen tersedia. Protonasi reseptor dan ligan "
            "tidak dioptimalkan pada pH fisiologis.",
            "Kontak hidrofobik adalah jarak atom karbon, bukan estimasi energi.",
            "Klasifikasi pi-pi memakai jarak centroid (<=5.5 A), sudut normal "
            "ring, dan kandidat offset (centroid <=6.5 A, atom-terdekat "
            "<=4.0 A); ambang ini heuristik dan perlu diverifikasi dengan "
            "PLIP/Arpeggio.",
            "Pose terbaik Vina belum tentu pose yang paling stabil secara "
            "dinamika; diperlukan klaster pose dan simulasi MD.",
            "Reseptor adalah struktur homolog surrogate, bukan struktur target "
            "C. cladosporioides; hasil tidak boleh dilaporkan sebagai "
            "interaksi target sebenarnya.",
        ],
    }


def _tsv_rows(payload: dict[str, Any]) -> Iterable[list[str]]:
    for complex_data in payload.get("complexes", []):
        complex_name = f"{complex_data['receptor']}__{complex_data['ligand']}"
        for item in complex_data.get("hbonds", []):
            yield [
                complex_name,
                complex_data["receptor"],
                complex_data["ligand"],
                "hbond",
                item["protein_residue"],
                item["protein_chain"],
                item["protein_atom"],
                item["ligand_atom"],
                str(item["distance_angstrom"]),
                str(item["angle_deg"]),
                f"{item['direction']}; confidence={item['confidence']}",
            ]
        for item in complex_data.get("polar_contacts", []):
            yield [
                complex_name,
                complex_data["receptor"],
                complex_data["ligand"],
                "polar_contact",
                item["protein_residue"],
                item["protein_chain"],
                item["protein_atom"],
                item["ligand_atom"],
                str(item["distance_angstrom"]),
                ""
                if item.get("angle_deg") is None
                else str(item["angle_deg"]),
                f"{item.get('reason', '')}; confidence="
                f"{item.get('confidence', 'none')}",
            ]
        for item in complex_data.get("hydrophobic_contacts", []):
            yield [
                complex_name,
                complex_data["receptor"],
                complex_data["ligand"],
                "hydrophobic",
                item["protein_residue"],
                item["protein_chain"],
                ",".join(item["protein_atoms"]),
                ",".join(item["ligand_atoms"]),
                str(item["min_distance_angstrom"]),
                "",
                f"contact_count={item['contact_count']}",
            ]
        for item in complex_data.get("pi_stacking", []):
            yield [
                complex_name,
                complex_data["receptor"],
                complex_data["ligand"],
                "pi_stacking",
                item["protein_ring"],
                complex_data["receptor_chain"],
                ",".join(item["protein_ring_atoms"]),
                ",".join(item["ligand_ring_atoms"]),
                str(item["centroid_distance_angstrom"]),
                str(item["angle_deg"]),
                f"{item['stacking_type']}; geometry={item['geometry']}; "
                f"min_atom_distance={item['min_atom_distance_angstrom']}",
            ]


def write_tsv(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    header = [
        "complex",
        "receptor",
        "ligand",
        "interaction_type",
        "protein_residue",
        "chain",
        "protein_atom",
        "ligand_atom",
        "distance_angstrom",
        "angle_deg",
        "detail",
    ]
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("\t".join(header) + "\n")
        for row in _tsv_rows(payload):
            handle.write(
                "\t".join(
                    value.replace("\t", " ").replace("\n", " ")
                    for value in row
                )
                + "\n"
            )


def _md_table(headers: Sequence[str], rows: Sequence[Sequence[Any]]) -> str:
    lines = [
        "| " + " | ".join(headers) + " |",
        "|" + "|".join(["---"] * len(headers)) + "|",
    ]
    for row in rows:
        lines.append(
            "| " + " | ".join(str(value).replace("|", "\\|") for value in row) + " |"
        )
    return "\n".join(lines)


def write_markdown(path: Path, payload: dict[str, Any]) -> None:
    lines = [
        "# Peta Interaksi Ligan-Reseptor",
        "",
        f"> Dibuat: {payload['generated_at']}",
        f"> Mode docking: {payload.get('docking_mode', 'unknown')}",
        ">",
        "> Analisis geometris dari pose docking. Bukan energi bebas ikat dan",
        "> bukan validasi eksperimental.",
        "",
        "## Ringkasan",
        "",
    ]
    rows = []
    for item in payload.get("complexes", []):
        rows.append(
            [
                item["receptor"],
                item["ligand"],
                f"{item['affinity_kcal_mol']:.2f}"
                if item["affinity_kcal_mol"] is not None
                else "-",
                f"{item['ligand_efficiency']:.3f}"
                if item["ligand_efficiency"] is not None
                else "-",
                item["summary"]["hbonds"],
                item["summary"]["hydrophobic_residues"],
                item["summary"]["pi_stacking"],
            ]
        )
    if rows:
        lines.extend(
            [
                _md_table(
                    [
                        "Reseptor",
                        "Ligan",
                        "Afinitas",
                        "LE",
                        "H-bond",
                        "Residu hidrofobik",
                        "Pi-pi",
                    ],
                    rows,
                ),
                "",
            ]
        )
    else:
        lines.extend(["Tidak ada kompleks yang berhasil dianalisis.", ""])

    for item in payload.get("complexes", []):
        lines.extend(
            [
                f"## {item['receptor']} — {item['ligand']}",
                "",
                f"- Pose: `{item['pose_file']}`",
                f"- Reseptor: `{item['receptor_file']}` (chain {item['receptor_chain']})",
                f"- PDB: {item.get('pdb_id') or '-'}; target: "
                f"{item.get('target_accession') or '-'}",
                f"- Afinitas Vina: "
                f"{item['affinity_kcal_mol']:.2f} kcal/mol"
                if item["affinity_kcal_mol"] is not None
                else "- Afinitas Vina: -",
                f"- Ligand efficiency: "
                f"{item['ligand_efficiency']:.3f} kcal/mol per heavy atom"
                if item["ligand_efficiency"] is not None
                else "- Ligand efficiency: -",
                "",
            ]
        )

        lines.extend(["### Ikatan hidrogen", ""])
        if item["hbonds"]:
            lines.extend(
                [
                    _md_table(
                        [
                            "Residu",
                            "Atom protein",
                            "Atom ligan",
                            "Jarak (Å)",
                            "Sudut (°)",
                            "Arah",
                        ],
                        [
                            [
                                hbond["protein_residue"],
                                hbond["protein_atom"],
                                hbond["ligand_atom"],
                                hbond["distance_angstrom"],
                                hbond["angle_deg"],
                                hbond["direction"],
                            ]
                            for hbond in item["hbonds"]
                        ],
                    ),
                    "",
                ]
            )
        else:
            lines.extend(
                [
                    "Tidak ada H-bond yang lolos kriteria jarak/sudut.",
                    "",
                ]
            )

        lines.extend(["### Kandidat polar / distance-only", ""])
        if item["polar_contacts"]:
            lines.extend(
                [
                    _md_table(
                        [
                            "Residu",
                            "Atom protein",
                            "Atom ligan",
                            "Jarak (Å)",
                            "Sudut (°)",
                            "Catatan",
                        ],
                        [
                            [
                                contact["protein_residue"],
                                contact["protein_atom"],
                                contact["ligand_atom"],
                                contact["distance_angstrom"],
                                contact["angle_deg"]
                                if contact.get("angle_deg") is not None
                                else "n/a",
                                contact.get("reason", "-"),
                            ]
                            for contact in item["polar_contacts"]
                        ],
                    ),
                    "",
                ]
            )
        else:
            lines.extend(
                [
                    "Tidak ada kontak polar tambahan pada cutoff.",
                    "",
                ]
            )

        lines.extend(["### Kontak hidrofobik", ""])
        if item["hydrophobic_contacts"]:
            lines.extend(
                [
                    _md_table(
                        [
                            "Residu",
                            "Jumlah kontak",
                            "Jarak minimum (Å)",
                            "Atom protein",
                        ],
                        [
                            [
                                contact["protein_residue"],
                                contact["contact_count"],
                                contact["min_distance_angstrom"],
                                ", ".join(contact["protein_atoms"]),
                            ]
                            for contact in item["hydrophobic_contacts"]
                        ],
                    ),
                    "",
                ]
            )
        else:
            lines.extend(["Tidak ada kontak hidrofobik pada cutoff.", ""])

        lines.extend(["### Stacking pi-pi", ""])
        if item["pi_stacking"]:
            lines.extend(
                [
                    _md_table(
                        [
                            "Ring protein",
                            "Ring ligan",
                            "Jarak centroid (Å)",
                            "Jarak atom min (Å)",
                            "Sudut (°)",
                            "Tipe",
                            "Geometri",
                        ],
                        [
                            [
                                stack["protein_ring"],
                                stack["ligand_ring"],
                                stack["centroid_distance_angstrom"],
                                stack["min_atom_distance_angstrom"],
                                stack["angle_deg"],
                                stack["stacking_type"],
                                stack["geometry"],
                            ]
                            for stack in item["pi_stacking"]
                        ],
                    ),
                    "",
                ]
            )
        else:
            lines.extend(["Tidak ada stacking pi-pi pada cutoff.", ""])

    lines.extend(["## Keterbatasan", ""])
    lines.extend([f"- {caveat}" for caveat in payload.get("caveats", [])])
    lines.append("")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Petakan H-bond, kontak hidrofobik, dan pi-pi dari pose docking "
            "AutoDock Vina."
        )
    )
    parser.add_argument(
        "--docking",
        type=Path,
        default=DEFAULT_DOCKING,
        help="File results/docking_results.json.",
    )
    parser.add_argument(
        "--receptor",
        nargs="+",
        default=None,
        help="Nama reseptor (default: semua). Gunakan 'all' untuk semua.",
    )
    parser.add_argument(
        "--ligand",
        nargs="+",
        default=["DHN_monomer"],
        help="Nama ligan (default: DHN_monomer; 'all' untuk semua).",
    )
    parser.add_argument("--chain", default=None, help="Override chain reseptor.")
    parser.add_argument("--cutoff-hbond", type=float, default=3.5)
    parser.add_argument("--min-hbond-angle", type=float, default=120.0)
    parser.add_argument("--cutoff-hydrophobic", type=float, default=4.0)
    parser.add_argument("--cutoff-pipi", type=float, default=5.5)
    parser.add_argument(
        "--max-pipi-centroid",
        type=float,
        default=6.5,
        help="Jarak centroid maksimum untuk kandidat pi-pi offset.",
    )
    parser.add_argument(
        "--cutoff-pipi-atom",
        type=float,
        default=4.0,
        help="Jarak atom-terdekat maksimum untuk kandidat pi-pi offset.",
    )
    parser.add_argument("--json", type=Path, default=DEFAULT_JSON)
    parser.add_argument("--tsv", type=Path, default=DEFAULT_TSV)
    parser.add_argument("--md", type=Path, default=DEFAULT_MD)
    parser.add_argument("--list", action="store_true")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        if args.list:
            payload = _read_json(args.docking)
            for row in (payload.get("results") or []):
                if row.get("status") in {"ok", "completed"} or (
                    row.get("status") is None
                    and row.get("best_affinity_kcal_mol") is not None
                ):
                    print(
                        f"{row.get('receptor')} x {row.get('ligand')} "
                        f"({row.get('best_affinity_kcal_mol')} kcal/mol)"
                    )
            return 0

        payload = build_payload(
            args.docking,
            receptors=args.receptor,
            ligands=args.ligand,
            chain=args.chain,
            cutoff_hbond=args.cutoff_hbond,
            min_hbond_angle=args.min_hbond_angle,
            cutoff_hydrophobic=args.cutoff_hydrophobic,
            cutoff_pipi=args.cutoff_pipi,
            max_pipi_centroid_offset=args.max_pipi_centroid,
            max_pipi_atom_distance=args.cutoff_pipi_atom,
        )
        _write_json_atomic(args.json, payload)
        write_tsv(args.tsv, payload)
        write_markdown(args.md, payload)
        print(
            f"[SUCCESS] {payload['summary']['complexes_total']} kompleks; "
            f"{payload['summary']['hbonds_total']} H-bond, "
            f"{payload['summary']['hydrophobic_residues_total']} residu "
            f"hidrofobik, {payload['summary']['pi_stacking_total']} pi-pi."
        )
        print(f"    JSON : {args.json}")
        print(f"    TSV  : {args.tsv}")
        print(f"    MD   : {args.md}")
        for error in payload.get("errors", []):
            print(f"    [WARN] {error['receptor']} x {error['ligand']}: {error['error']}")
        return 0 if payload["status"] != "failed" else 1
    except (AnalysisError, OSError, ValueError, KeyError, TypeError) as exc:
        print(f"[ERROR] {exc}")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
