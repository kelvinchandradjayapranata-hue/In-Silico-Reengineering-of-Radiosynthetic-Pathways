"""Docking pilot DHN-melanin terhadap reseptor PKS1 dan laccase.

Pipeline:
1. Bangun model ligan DHN monomer, dimer, dan trimer dengan RDKit.
2. Siapkan ligan PDBQT dengan Meeko.
3. Siapkan reseptor PDBQT dengan Meeko (reseptor harus PDB; default memakai
   struktur homolog PDB 1GW0 untuk laccase dan 8CG4 untuk PKS1).
4. Jalankan AutoDock Vina pada kotak aktif-site yang telah ditentukan.
5. Simpan skor ke results/docking_results.json dan results/docking_summary.tsv.

PENTING: struktur target *C. cladosporioides* (KAN3049821.1 dan KAN3060494.1)
belum tersedia di PDB/AlphaFold DB. Modul ini karena itu memakai struktur
homolog sebagai *surrogate*:
- Laccase: PDB 1GW0 (Melanocarpus albomyces laccase), identitas ~44.7%
  terhadap query, coverage ~91.7%.
- PKS1: PDB 8CG4 (Aspergillus parasiticus norsolorinic acid synthase),
  identitas ~36.4%, coverage ~94.7% (konstruk KS-AT, residu 3-1254).

Skor Vina dari reseptor surrogate BUKAN afinitas target sebenarnya. Angka ini
hanya pilot komputasi dan tidak boleh dilaporkan sebagai energi ikat
C. cladosporioides.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import os
import re
import shutil
import subprocess
import sys
import time
import urllib.error
import urllib.request
from datetime import datetime, timezone
from math import dist
from pathlib import Path
from typing import Any, Sequence

from Bio.PDB import PDBParser
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = PROJECT_ROOT / "data"
LIGAND_DIR = DATA_DIR / "ligands"
STRUCTURE_DIR = DATA_DIR / "structures"
PREPARED_DIR = STRUCTURE_DIR / "prepared"
RESULTS_DIR = PROJECT_ROOT / "results"
DOCKING_DIR = RESULTS_DIR / "docking"
DEFAULT_JSON = RESULTS_DIR / "docking_results.json"
DEFAULT_TSV = RESULTS_DIR / "docking_summary.tsv"
TOOLS_DIR = PROJECT_ROOT / "tools"

USER_AGENT = "melanin-biotech-docking/1.0"
DHN_SMILES = "Oc1cccc2cccc(O)c12"

LIGAND_SPECS: tuple[tuple[str, int, str], ...] = (
    (
        "DHN_monomer",
        1,
        "1,8-dihydroxynaphthalene monomer (PubChem CID 68438)",
    ),
    (
        "DHN_dimer",
        2,
        "Model dimer 2,2'-linked; bukan struktur oligomer alami yang unik",
    ),
    (
        "DHN_trimer",
        3,
        "Model trimer linear: unit tengah terhubung di posisi 2 dan 7 "
        "(2,2'/2,7'-linked; posisi 2/7 ekuivalen secara simetri); bukan "
        "struktur oligomer alami yang unik",
    ),
)

DEFAULT_RECEPTORS: tuple[dict[str, Any], ...] = (
    {
        "name": "Laccase_surrogate_1GW0",
        "target_id": "ccla_laccase",
        "target_accession": "KAN3060494.1",
        "pdb_id": "1GW0",
        "chain": "A",
        "box_center": None,
        "box_center_source": "T1 copper (Cys SG within 2.8 A)",
        "box_size": [22.0, 22.0, 22.0],
        "surrogate": True,
        "identity_percent": 44.7,
        "coverage_percent": 91.7,
        "note": (
            "Laccase dari Melanocarpus albomyces; homolog terdekat yang "
            "tersedia di PDB. Tembaga T1 dipakai untuk pusat kotak, tetapi "
            "logam dihapus dari reseptor PDBQT (reseptor apo)."
        ),
    },
    {
        "name": "PKS1_surrogate_8CG4",
        "target_id": "ccla_pks1",
        "target_accession": "KAN3049821.1",
        "pdb_id": "8CG4",
        "chain": "A",
        "box_center": [169.803, 164.0, 134.493],
        "box_center_source": "SG atom Cys543 (katalitik KS motif DTACSSSL)",
        "box_size": [22.0, 22.0, 22.0],
        "surrogate": True,
        "identity_percent": 36.4,
        "coverage_percent": 94.7,
        "note": (
            "Norsolorinic acid synthase dari Aspergillus parasiticus; "
            "konstruk KS-AT (residu 3-1254). DHN adalah produk, bukan substrat "
            "PKS, sehingga docking ini hanya eksplorasi kualitatif."
        ),
    },
)


class DockingError(RuntimeError):
    """Kesalahan yang dapat ditampilkan ke pengguna tanpa traceback."""


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _download(url: str, dest: Path, timeout: int = 120, retries: int = 3) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    last_error: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(request, timeout=timeout) as response:
                with tmp.open("wb") as handle:
                    shutil.copyfileobj(response, handle, length=1024 * 1024)
            if tmp.stat().st_size < 1024:
                raise DockingError(f"Unduhan terlalu kecil: {url}")
            tmp.replace(dest)
            return
        except (OSError, urllib.error.URLError) as exc:
            last_error = exc
            tmp.unlink(missing_ok=True)
            if attempt < retries:
                time.sleep(2 * attempt)
    raise DockingError(f"Gagal mengunduh {url}: {last_error}")


def _portable_paths(value: Any) -> Any:
    """Ubah path absolut di dalam PROJECT_ROOT menjadi path relatif."""

    if isinstance(value, dict):
        return {key: _portable_paths(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_portable_paths(item) for item in value]
    if isinstance(value, str):
        root = PROJECT_ROOT.resolve()
        try:
            path = Path(value)
            if path.is_absolute():
                resolved = path.resolve()
                if resolved.is_relative_to(root):
                    return resolved.relative_to(root).as_posix()
        except (OSError, ValueError):
            pass
        for root_text in {str(root), root.as_posix()}:
            value = value.replace(f"{root_text}{os.sep}", "").replace(
                f"{root_text}/", ""
            ).replace(root_text, "<PROJECT_ROOT>")
    return value


def _write_json_atomic(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    tmp.replace(path)


def _run(
    command: Sequence[str | Path],
    *,
    log_path: Path | None = None,
    timeout: int | None = None,
) -> subprocess.CompletedProcess[str]:
    """Jalankan proses; jika log_path diberikan, stdout/stderr dialihkan ke file."""

    if log_path is not None:
        log_path.parent.mkdir(parents=True, exist_ok=True)
        with log_path.open("w", encoding="utf-8") as log_handle:
            return subprocess.run(
                [str(part) for part in command],
                stdout=log_handle,
                stderr=subprocess.STDOUT,
                text=True,
                encoding="utf-8",
                errors="replace",
                timeout=timeout,
                check=False,
            )
    return subprocess.run(
        [str(part) for part in command],
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        timeout=timeout,
        check=False,
    )


def _couple(
    mol_a: Chem.Mol, mol_b: Chem.Mol, index_a: int, index_b: int
) -> Chem.Mol:
    """Gabungkan dua unit DHN dengan ikatan C-C tunggal (model oligomer)."""

    combined = Chem.CombineMols(mol_a, mol_b)
    editable = Chem.RWMol(combined)
    editable.AddBond(index_a, mol_a.GetNumAtoms() + index_b, Chem.BondType.SINGLE)
    mol = editable.GetMol()
    Chem.SanitizeMol(mol)
    return mol


def build_ligand_molecules() -> list[tuple[str, int, str, Chem.Mol]]:
    """Bangun monomer, dimer, dan trimer DHN secara deterministik."""

    monomer = Chem.MolFromSmiles(DHN_SMILES)
    if monomer is None:
        raise DockingError("SMILES DHN tidak dapat diparsing.")
    dimer = _couple(monomer, monomer, 2, 2)
    trimer = _couple(dimer, monomer, 8, 2)
    molecules = {
        "DHN_monomer": monomer,
        "DHN_dimer": dimer,
        "DHN_trimer": trimer,
    }
    return [
        (name, units, description, molecules[name])
        for name, units, description in LIGAND_SPECS
    ]


def _find_meeko_cli(name: str, explicit: Path | None = None) -> Path:
    candidates: list[Path] = []
    if explicit is not None:
        candidates.append(explicit)
    found = shutil.which(name)
    if found:
        candidates.append(Path(found))
    scripts_dir = Path(sys.executable).resolve().parent
    candidates.append(scripts_dir / f"{name}.exe")
    candidates.append(scripts_dir / name)
    for candidate in candidates:
        if candidate.exists():
            return candidate
    raise DockingError(
        f"{name} tidak ditemukan. Install Meeko di venv atau gunakan opsi path."
    )


def _find_vina(explicit: Path | None = None) -> Path:
    candidates: list[Path] = []
    if explicit is not None:
        candidates.append(explicit)
    found = shutil.which("vina")
    if found:
        candidates.append(Path(found))
    if (TOOLS_DIR / "vina").is_dir():
        candidates.extend(sorted((TOOLS_DIR / "vina").glob("vina_*_win.exe")))
        candidates.extend(sorted((TOOLS_DIR / "vina").glob("vina*")))
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    raise DockingError(
        "AutoDock Vina tidak ditemukan. Letakkan binary di tools/vina/ atau "
        "gunakan --vina."
    )


def _version(command: Sequence[str | Path]) -> str:
    process = _run(command)
    output = (process.stdout or process.stderr).strip()
    return output.splitlines()[0] if output else "unknown"


def _ligand_descriptors(mol: Chem.Mol) -> dict[str, Any]:
    return {
        "canonical_smiles": Chem.MolToSmiles(mol),
        "formula": rdMolDescriptors.CalcMolFormula(mol),
        "molecular_weight": round(Descriptors.MolWt(mol), 2),
        "logp": round(Descriptors.MolLogP(mol), 3),
        "tpsa": round(Descriptors.TPSA(mol), 2),
        "h_bond_donors": rdMolDescriptors.CalcNumHBD(mol),
        "h_bond_acceptors": rdMolDescriptors.CalcNumHBA(mol),
        "rotatable_bonds": rdMolDescriptors.CalcNumRotatableBonds(mol),
        "heavy_atoms": mol.GetNumHeavyAtoms(),
    }


def _torsdof(pdbqt_path: Path) -> int | None:
    try:
        text = pdbqt_path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return None
    for line in text.splitlines():
        if line.startswith("TORSDOF"):
            match = re.search(r"(\d+)", line)
            if match:
                return int(match.group(1))
    return text.count("BRANCH")


def prepare_ligands(
    names: Sequence[str],
    meeko_ligand: Path,
    *,
    force: bool = False,
) -> list[dict[str, Any]]:
    """Siapkan SDF dan PDBQT ligan; kembalikan metadata."""

    LIGAND_DIR.mkdir(parents=True, exist_ok=True)
    prepared: list[dict[str, Any]] = []
    available = {name for name, _, _ in LIGAND_SPECS}
    unknown = sorted(set(names) - available)
    if unknown:
        raise DockingError(f"Ligan tidak dikenal: {', '.join(unknown)}")

    for name, units, description, mol in build_ligand_molecules():
        if name not in names:
            continue

        sdf_path = LIGAND_DIR / f"{name}.sdf"
        pdbqt_path = LIGAND_DIR / f"{name}.pdbqt"
        if force or not sdf_path.exists():
            mol_with_h = Chem.AddHs(mol)
            params = AllChem.ETKDGv3()
            params.randomSeed = 42
            if AllChem.EmbedMolecule(mol_with_h, params) != 0:
                raise DockingError(f"Gagal membangun konformer 3D untuk {name}.")
            AllChem.MMFFOptimizeMolecule(mol_with_h, maxIters=500)
            mol_with_h.SetProp("_Name", name)
            writer = Chem.SDWriter(str(sdf_path))
            writer.write(mol_with_h)
            writer.close()

        if force or not pdbqt_path.exists():
            process = _run([meeko_ligand, "-i", sdf_path, "-o", pdbqt_path])
            if process.returncode != 0 or not pdbqt_path.exists():
                raise DockingError(
                    f"Meeko gagal menyiapkan ligan {name}: "
                    f"{(process.stderr or process.stdout).strip()}"
                )

        record = {
            "name": name,
            "units": units,
            "description": description,
            "sdf": str(sdf_path),
            "pdbqt": str(pdbqt_path),
            "torsdof": _torsdof(pdbqt_path),
            "model_note": (
                "Model oligomer sintetis untuk eksplorasi; linkage bukan "
                "struktur DHN-melanin alami yang tunggal."
                if units > 1
                else "Monomer DHN."
            ),
        }
        record.update(_ligand_descriptors(mol))
        prepared.append(record)
    return prepared


def load_receptors(config_path: Path | None) -> list[dict[str, Any]]:
    if config_path is None:
        return [dict(entry) for entry in DEFAULT_RECEPTORS]
    if not config_path.is_file():
        raise DockingError(f"Config reseptor tidak ditemukan: {config_path}")
    payload = json.loads(config_path.read_text(encoding="utf-8"))
    receptors = payload.get("receptors")
    if not isinstance(receptors, list) or not receptors:
        raise DockingError("Config harus memiliki list 'receptors' yang tidak kosong.")
    return receptors


def ensure_receptor_pdb(
    receptor: dict[str, Any], *, force: bool = False
) -> dict[str, Any]:
    """Gunakan PDB lokal atau unduh dari RCSB berdasarkan pdb_id."""

    name = str(receptor.get("name", "")).strip()
    if not name:
        raise DockingError("Setiap reseptor harus memiliki 'name'.")

    local_path = receptor.get("pdb_path")
    if local_path:
        pdb_path = Path(str(local_path))
        if not pdb_path.is_absolute():
            pdb_path = PROJECT_ROOT / pdb_path
        if not pdb_path.is_file():
            raise DockingError(f"PDB lokal tidak ditemukan: {pdb_path}")
        source_url = None
        pdb_id = receptor.get("pdb_id")
    else:
        pdb_id = str(receptor.get("pdb_id", "")).strip().upper()
        if not pdb_id:
            raise DockingError(f"Reseptor {name} tidak punya pdb_id/pdb_path.")
        pdb_path = STRUCTURE_DIR / f"{pdb_id}.pdb"
        source_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        if force or not pdb_path.exists():
            print(f"    mengunduh struktur {pdb_id} ...")
            _download(source_url, pdb_path)

    if not pdb_path.read_text(encoding="utf-8", errors="replace").lstrip().startswith(
        ("HEADER", "ATOM", "REMARK", "TITLE", "CRYST1", "MODEL")
    ):
        raise DockingError(f"File bukan PDB yang valid: {pdb_path}")

    metadata = {
        "name": name,
        "target_id": receptor.get("target_id"),
        "target_accession": receptor.get("target_accession"),
        "pdb_id": pdb_id,
        "pdb_path": str(pdb_path),
        "source_url": source_url,
        "sha256": _sha256(pdb_path),
        "size_bytes": pdb_path.stat().st_size,
        "surrogate": bool(receptor.get("surrogate", False)),
        "identity_percent": receptor.get("identity_percent"),
        "coverage_percent": receptor.get("coverage_percent"),
        "note": receptor.get("note", ""),
    }
    return metadata


def _clean_receptor_pdb(pdb_path: Path, chain: str, out_path: Path) -> int:
    kept: list[str] = []
    for line in pdb_path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith("ATOM"):
            continue
        if chain and line[21:22] != chain:
            continue
        altloc = line[16:17]
        if altloc not in (" ", "A"):
            continue
        kept.append(line)
    if not kept:
        raise DockingError(
            f"Tidak ada atom chain {chain!r} pada struktur {pdb_path.name}."
        )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(kept) + "\nTER\nEND\n", encoding="utf-8")
    return len(kept)


def _find_t1_copper(pdb_path: Path, chain: str) -> dict[str, Any] | None:
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("receptor", str(pdb_path))
    model = next(structure.get_models())
    if chain not in model:
        return None
    chain_obj = model[chain]
    cys_sg = [
        atom
        for residue in chain_obj
        if residue.get_resname() == "CYS"
        for atom in residue
        if atom.name == "SG"
    ]
    if not cys_sg:
        return None

    candidates: list[tuple[float, Any, Any]] = []
    for atom in chain_obj.get_atoms():
        if atom.element != "CU":
            continue
        for sulfur in cys_sg:
            distance = dist(atom.coord, sulfur.coord)
            if distance <= 2.8:
                candidates.append((distance, atom, sulfur))
    if not candidates:
        return None

    distance, copper, sulfur = min(candidates, key=lambda item: item[0])
    return {
        "center": [round(float(value), 3) for value in copper.coord],
        "distance_to_cys_sg": round(float(distance), 3),
        "copper_serial": copper.get_serial_number(),
        "cys_residue": sulfur.get_parent().id[1],
    }


def prepare_receptor(
    receptor: dict[str, Any],
    metadata: dict[str, Any],
    meeko_receptor: Path,
    *,
    force: bool = False,
) -> dict[str, Any]:
    """Bersihkan PDB, tentukan kotak, dan siapkan PDBQT reseptor."""

    name = str(receptor["name"])
    chain = str(receptor.get("chain", "A"))
    pdb_path = Path(metadata["pdb_path"])
    clean_path = PREPARED_DIR / f"{name}_clean.pdb"
    prefix = PREPARED_DIR / name
    pdbqt_path = PREPARED_DIR / f"{name}.pdbqt"
    box_path = PREPARED_DIR / f"{name}_box.txt"

    atoms_kept = _clean_receptor_pdb(pdb_path, chain, clean_path)
    box_center = receptor.get("box_center")
    box_center_source = receptor.get("box_center_source", "config")
    t1_info = None
    if box_center is None:
        if receptor.get("box_center_source", "").startswith("T1 copper"):
            t1_info = _find_t1_copper(pdb_path, chain)
            if t1_info is None:
                raise DockingError(
                    f"Tembaga T1 tidak ditemukan pada {pdb_path.name}; "
                    "isi box_center secara eksplisit di config."
                )
            box_center = t1_info["center"]
            box_center_source = (
                f"T1 copper serial {t1_info['copper_serial']} "
                f"(Cys{t1_info['cys_residue']} SG, "
                f"{t1_info['distance_to_cys_sg']} A)"
            )
        else:
            raise DockingError(
                f"Reseptor {name} tidak punya box_center dan bukan mode T1 copper."
            )

    box_size = receptor.get("box_size") or [22.0, 22.0, 22.0]
    if len(box_center) != 3 or len(box_size) != 3:
        raise DockingError(f"box_center/box_size {name} harus berisi 3 angka.")
    box_center = [float(value) for value in box_center]
    box_size = [float(value) for value in box_size]

    if force or not pdbqt_path.exists():
        command = [
            meeko_receptor,
            "--read_pdb",
            clean_path,
            "-o",
            prefix,
            "-p",
            pdbqt_path,
            "--box_center",
            *[f"{value:.3f}" for value in box_center],
            "--box_size",
            *[f"{value:.3f}" for value in box_size],
            "--write_vina_box",
            box_path,
        ]
        process = _run(command)
        if process.returncode != 0 or not pdbqt_path.exists():
            raise DockingError(
                f"Meeko gagal menyiapkan reseptor {name}: "
                f"{(process.stderr or process.stdout).strip()}"
            )

    record = dict(metadata)
    record.update(
        {
            "chain": chain,
            "atoms_kept": atoms_kept,
            "clean_pdb": str(clean_path),
            "pdbqt": str(pdbqt_path),
            "box_file": str(box_path),
            "box_center": box_center,
            "box_size": box_size,
            "box_center_source": box_center_source,
            "t1_copper": t1_info,
        }
    )
    return record


def parse_vina_scores(log_text: str) -> list[dict[str, Any]]:
    """Parse tabel affinity dari output Vina 1.2.x."""

    poses: list[dict[str, Any]] = []
    in_table = False
    for line in log_text.splitlines():
        if "mode |" in line and "affinity" in line:
            in_table = True
            continue
        if not in_table:
            continue
        match = re.match(
            r"^\s*(\d+)\s+(-?\d+(?:\.\d+)?)\s+"
            r"(-?\d+(?:\.\d+)?)\s+(-?\d+(?:\.\d+)?)\s*$",
            line,
        )
        if match:
            poses.append(
                {
                    "mode": int(match.group(1)),
                    "affinity_kcal_mol": float(match.group(2)),
                    "rmsd_lb": float(match.group(3)),
                    "rmsd_ub": float(match.group(4)),
                }
            )
    return poses


def run_docking(
    *,
    receptors: list[dict[str, Any]],
    ligand_names: Sequence[str],
    meeko_ligand: Path,
    meeko_receptor: Path,
    vina: Path,
    force: bool = False,
    prepare_only: bool = False,
    exhaustiveness: int = 8,
    num_modes: int = 9,
    cpu: int = 4,
    seed: int = 42,
) -> dict[str, Any]:
    """Jalankan seluruh kombinasi reseptor x ligan."""

    DOCKING_DIR.mkdir(parents=True, exist_ok=True)
    (DOCKING_DIR / "poses").mkdir(parents=True, exist_ok=True)
    (DOCKING_DIR / "logs").mkdir(parents=True, exist_ok=True)

    ligands = prepare_ligands(ligand_names, meeko_ligand, force=force)
    results: list[dict[str, Any]] = []
    receptors_out: list[dict[str, Any]] = []

    for receptor in receptors:
        try:
            metadata = ensure_receptor_pdb(receptor, force=force)
            prepared = prepare_receptor(
                receptor, metadata, meeko_receptor, force=force
            )
            receptors_out.append(prepared)
        except DockingError as exc:
            receptors_out.append(
                {
                    "name": receptor.get("name"),
                    "status": "error",
                    "error": str(exc),
                }
            )
            continue

        if prepare_only:
            continue

        for ligand in ligands:
            stem = f"{receptor['name']}__{ligand['name']}__seed{seed}"
            pose_path = DOCKING_DIR / "poses" / f"{stem}.pdbqt"
            log_path = DOCKING_DIR / "logs" / f"{stem}.log"
            command = [
                vina,
                "--receptor",
                prepared["pdbqt"],
                "--ligand",
                ligand["pdbqt"],
                "--center_x",
                f"{prepared['box_center'][0]:.3f}",
                "--center_y",
                f"{prepared['box_center'][1]:.3f}",
                "--center_z",
                f"{prepared['box_center'][2]:.3f}",
                "--size_x",
                f"{prepared['box_size'][0]:.3f}",
                "--size_y",
                f"{prepared['box_size'][1]:.3f}",
                "--size_z",
                f"{prepared['box_size'][2]:.3f}",
                "--exhaustiveness",
                str(exhaustiveness),
                "--num_modes",
                str(num_modes),
                "--cpu",
                str(cpu),
                "--seed",
                str(seed),
                "--out",
                pose_path,
                "--verbosity",
                "1",
            ]
            try:
                process = _run(command, log_path=log_path, timeout=3600)
                log_text = log_path.read_text(encoding="utf-8", errors="replace")
                poses = parse_vina_scores(log_text)
                if process.returncode != 0 or not poses:
                    raise DockingError(
                        (process.stderr or log_text).strip().splitlines()[-1]
                        if (process.stderr or log_text).strip()
                        else "Vina tidak menghasilkan pose."
                    )
                results.append(
                    {
                        "receptor": receptor["name"],
                        "ligand": ligand["name"],
                        "status": "ok",
                        "best_affinity_kcal_mol": poses[0]["affinity_kcal_mol"],
                        "poses": poses,
                        "pose_file": str(pose_path),
                        "log": str(log_path),
                    }
                )
            except (DockingError, OSError, subprocess.TimeoutExpired) as exc:
                results.append(
                    {
                        "receptor": receptor["name"],
                        "ligand": ligand["name"],
                        "status": "error",
                        "error": str(exc),
                        "pose_file": str(pose_path),
                        "log": str(log_path),
                    }
                )

    heavy_atoms = {ligand["name"]: ligand["heavy_atoms"] for ligand in ligands}
    for row in results:
        if row.get("status") == "ok":
            heavy = heavy_atoms.get(row["ligand"])
            if heavy:
                row["ligand_efficiency"] = round(
                    row["best_affinity_kcal_mol"] / heavy, 3
                )

    ok = [row for row in results if row.get("status") == "ok"]
    failed = [row for row in results if row.get("status") != "ok"]
    matrix: dict[str, dict[str, float | None]] = {}
    efficiency_matrix: dict[str, dict[str, float | None]] = {}
    for receptor in receptors_out:
        receptor_name = receptor.get("name", "?")
        matrix[receptor_name] = {}
        efficiency_matrix[receptor_name] = {}
        for ligand in ligands:
            match = next(
                (
                    row
                    for row in results
                    if row["receptor"] == receptor_name
                    and row["ligand"] == ligand["name"]
                    and row["status"] == "ok"
                ),
                None,
            )
            matrix[receptor_name][ligand["name"]] = (
                match["best_affinity_kcal_mol"] if match else None
            )
            efficiency_matrix[receptor_name][ligand["name"]] = (
                match.get("ligand_efficiency") if match else None
            )

    if prepare_only:
        status = "prepared"
    elif failed and ok:
        status = "partial"
    elif failed and not ok:
        status = "failed"
    else:
        status = "completed"

    if ok:
        best = min(ok, key=lambda row: row["best_affinity_kcal_mol"])
        best_le = min(ok, key=lambda row: row["ligand_efficiency"])
        interpretation = (
            f"{len(ok)} dari {len(results)} run docking berhasil. Skor mentah "
            f"terbaik {best['best_affinity_kcal_mol']:.3f} kcal/mol untuk "
            f"{best['ligand']} pada {best['receptor']}; ligand efficiency "
            f"terbaik {best_le['ligand_efficiency']:.3f} kcal/mol per heavy atom "
            f"untuk {best_le['ligand']} pada {best_le['receptor']}. Skor Vina "
            "tidak dinormalisasi terhadap ukuran ligan, sehingga monomer, dimer, "
            "dan trimer harus dibandingkan dengan ligand efficiency, bukan skor "
            "mentah saja. Semua reseptor adalah surrogate homolog, sehingga "
            "angka ini bukan afinitas target C. cladosporioides dan hanya untuk "
            "perbandingan internal."
        )
    elif prepare_only:
        interpretation = (
            "Persiapan ligan dan reseptor selesai; docking Vina dilewati "
            "(--prepare-only)."
        )
    else:
        interpretation = (
            "Tidak ada run docking yang berhasil. Periksa log dan dependency."
        )

    return {
        "generated_at": _now(),
        "tool": "melanin_biotech_docking",
        "status": status,
        "engine": {
            "vina": _version([vina, "--version"]),
            "meeko": importlib.metadata.version("meeko"),
            "rdkit": importlib.metadata.version("rdkit"),
            "python": sys.version.split()[0],
        },
        "surrogate_notice": (
            "Reseptor default adalah struktur homolog PDB (1GW0 untuk laccase, "
            "8CG4 untuk PKS1), bukan struktur C. cladosporioides. Skor Vina "
            "tidak boleh ditafsirkan sebagai energi ikat target sebenarnya."
        ),
        "parameters": {
            "exhaustiveness": exhaustiveness,
            "num_modes": num_modes,
            "cpu": cpu,
            "seed": seed,
            "box_size_default": [22.0, 22.0, 22.0],
            "prepare_only": prepare_only,
        },
        "ligands": ligands,
        "receptors": receptors_out,
        "results": results,
        "summary": {
            "runs_total": len(results),
            "runs_ok": len(ok),
            "runs_failed": len(failed),
            "affinity_matrix": matrix,
            "ligand_efficiency_matrix": efficiency_matrix,
        },
        "assumptions": [
            "Ligan DHN dimer/trimer adalah model sintetis dengan linkage C-C sederhana, bukan struktur oligomer alami yang unik.",
            "Reseptor default adalah homolog PDB, bukan struktur target C. cladosporioides.",
            "Logam tembaga dan glikan dihapus dari reseptor laccase; reseptor bersifat apo.",
            "Kotak docking berpusat pada situs aktif T1 laccase atau Cys katalitik KS PKS1.",
            "Vina mensimulasikan docking rigid-receptor; fleksibilitas side chain tidak dimodelkan.",
        ],
        "caveats": [
            "Skor Vina bukan energi bebas ikat eksperimental (bukan ΔG).",
            "Tidak ada validasi eksperimental atau simulasi MD untuk pose terbaik.",
            "PKS1 mendock DHN sebagai produk, bukan substrat; relevansi biologisnya spekulatif.",
            "Model dimer/trimer tidak mewakili heterogenitas melanin alami.",
            "Skor Vina tidak dinormalisasi terhadap ukuran ligan; ligan lebih besar cenderung memperoleh skor lebih favorable. Bandingkan ligand efficiency, bukan hanya afinitas mentah.",
            "Perbedaan kecil skor (<1 kcal/mol) tidak dapat dianggap signifikan.",
        ],
        "interpretation": interpretation,
    }


def write_summary_tsv(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(
            "receptor\tligand\tstatus\tbest_affinity_kcal_mol\t"
            "ligand_efficiency_kcal_mol_per_heavy_atom\tmode1_rmsd_lb\terror\n"
        )
        for row in payload.get("results", []):
            poses = row.get("poses") or []
            rmsd = poses[0].get("rmsd_lb") if poses else ""
            handle.write(
                "\t".join(
                    [
                        str(row.get("receptor", "")),
                        str(row.get("ligand", "")),
                        str(row.get("status", "")),
                        str(row.get("best_affinity_kcal_mol", "")),
                        str(row.get("ligand_efficiency", "")),
                        str(rmsd),
                        str(row.get("error", "")).replace("\t", " "),
                    ]
                )
                + "\n"
            )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Docking pilot DHN monomer/dimer/trimer ke reseptor PKS1/laccase "
            "(default: struktur homolog PDB)."
        )
    )
    parser.add_argument(
        "--ligands",
        nargs="+",
        default=[name for name, _, _ in LIGAND_SPECS],
        help="Ligan yang dipakai.",
    )
    parser.add_argument(
        "--receptors",
        nargs="+",
        default=None,
        help="Nama reseptor dari config (default: semua).",
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=None,
        help="JSON config reseptor; default memakai registry bawaan.",
    )
    parser.add_argument("--force", action="store_true", help="Regenerasi file.")
    parser.add_argument(
        "--prepare-only",
        action="store_true",
        help="Siapkan ligan/reseptor tanpa menjalankan Vina.",
    )
    parser.add_argument("--exhaustiveness", type=int, default=8)
    parser.add_argument("--num-modes", type=int, default=9)
    parser.add_argument(
        "--box-size",
        nargs=3,
        type=float,
        metavar=("X", "Y", "Z"),
        default=None,
        help=(
            "Override ukuran kotak docking (Angstrom) untuk semua reseptor. "
            "Default mengikuti config per reseptor (22 A)."
        ),
    )
    parser.add_argument(
        "--cpu", type=int, default=min(4, os.cpu_count() or 1)
    )
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--vina", type=Path, default=None)
    parser.add_argument("--meeko-ligand", type=Path, default=None)
    parser.add_argument("--meeko-receptor", type=Path, default=None)
    parser.add_argument("--out", type=Path, default=DEFAULT_JSON)
    parser.add_argument("--summary", type=Path, default=DEFAULT_TSV)
    parser.add_argument(
        "--list",
        action="store_true",
        help="Tampilkan ligan dan reseptor yang tersedia, lalu keluar.",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    try:
        receptors = load_receptors(args.config)
        if args.receptors:
            selected = set(args.receptors)
            receptors = [
                receptor
                for receptor in receptors
                if receptor.get("name") in selected
            ]
            if not receptors:
                raise DockingError(
                    "Tidak ada reseptor yang cocok dengan --receptors."
                )

        if args.box_size:
            if any(value <= 0 for value in args.box_size):
                raise DockingError("--box-size harus berisi angka positif.")
            for receptor in receptors:
                receptor["box_size"] = [float(value) for value in args.box_size]

        if args.list:
            print("Ligan:")
            for name, units, description in LIGAND_SPECS:
                print(f"  - {name} ({units} unit): {description}")
            print("\nReseptor:")
            for receptor in receptors:
                print(
                    f"  - {receptor.get('name')}  target={receptor.get('target_id')} "
                    f"pdb={receptor.get('pdb_id', receptor.get('pdb_path'))}"
                )
            return 0

        meeko_ligand = _find_meeko_cli("mk_prepare_ligand", args.meeko_ligand)
        meeko_receptor = _find_meeko_cli(
            "mk_prepare_receptor", args.meeko_receptor
        )
        vina = _find_vina(args.vina)

        print("[1/2] Menyiapkan ligan dan reseptor...")
        payload = run_docking(
            receptors=receptors,
            ligand_names=args.ligands,
            meeko_ligand=meeko_ligand,
            meeko_receptor=meeko_receptor,
            vina=vina,
            force=args.force,
            prepare_only=args.prepare_only,
            exhaustiveness=args.exhaustiveness,
            num_modes=args.num_modes,
            cpu=args.cpu,
            seed=args.seed,
        )
        _write_json_atomic(args.out, _portable_paths(payload))
        write_summary_tsv(args.summary, payload)

        print("[2/2] Selesai.")
        print(f"    JSON    : {args.out}")
        print(f"    Summary : {args.summary}")
        print(f"    Status  : {payload['status']}")
        print(f"    {payload['interpretation']}")
        return 0

    except (DockingError, OSError, ValueError, json.JSONDecodeError) as exc:
        print(f"[ERROR] {exc}")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
