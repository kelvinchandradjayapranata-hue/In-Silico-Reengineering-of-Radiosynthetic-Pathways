"""Preparasi sistem protein-ligan untuk MD eksplisit dengan OpenMM.

Modul ini menyiapkan pose docking Vina terbaik menjadi sistem MD eksplisit:

1. membaca pose PDBQT dan reseptor clean PDB dari ``results/docking_results.json``;
2. menambahkan hidrogen/protonasi protein pada pH 7.4 dengan PDBFixer;
3. menambahkan hidrogen ligan dan memparameterisasinya dengan OpenFF
   (default ``openff-2.2.1``) melalui ``openmmforcefields``;
4. menyolvasi dengan TIP3P dan NaCl 0.15 M;
5. menulis ``system.xml``, ``system.pdb``, ``state.xml``, ``md_config.json``,
   ``run_md.py``, dan ``run_colab.sh``.

Pekerjaan berat (minimisasi, ekuilibrasi, produksi MD) dijalankan oleh
``run_md.py`` di Google Colab GPU. Modul ini sendiri hanya melakukan preparasi
struktur dan system setup; ia tidak menjalankan simulasi panjang.

Catatan ilmiah penting:

- PKS1 8CG4 adalah konstruk KS-AT dari norsolorinic acid synthase
  *Aspergillus parasiticus*; DHN adalah produk, bukan substrat PKS.
- Laccase 1GW0 adalah struktur homolog *Melanocarpus albomyces*. Clean PDB
  proyek tidak mengandung atom tembaga, sehingga sistem MD laccase adalah
  **apo** dan tidak merepresentasikan pusat katalitik T1/T2/T3.
- Parameter ligan berasal dari OpenFF/GAFF dan bukan parameter kuantum khusus
  untuk DHN-melanin. Hasil MD adalah eksplorasi kestabilan pose, bukan bukti
  afinitas atau reaktivitas enzim.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
import os
import re
import shutil
import statistics
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Sequence

from rdkit import Chem
from rdkit.Chem import AllChem

try:  # jalankan sebagai skrip (python src/prepare_openmm_md.py)
    import analyze_interactions as ai
except ImportError:  # impor sebagai paket (python -m src.prepare_openmm_md)
    from . import analyze_interactions as ai

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DOCKING = PROJECT_ROOT / "results" / "docking_results.json"
DEFAULT_OUT_DIR = PROJECT_ROOT / "results" / "md_prep"
DEFAULT_RECEPTORS = ("PKS1_surrogate_8CG4", "Laccase_surrogate_1GW0")
DEFAULT_LIGAND = "DHN_monomer"
LIGAND_RESIDUE_NAME = "LIG"

_RUNNER_TEMPLATE = r'''#!/usr/bin/env python3
"""Runner MD OpenMM yang dihasilkan oleh src/prepare_openmm_md.py.

Jalankan di Google Colab GPU:

    python run_md.py --config md_config.json

atau di CPU:

    python run_md.py --config md_config.json --platform CPU

Runner ini membaca ``system.xml`` dan ``state.xml``, lalu menjalankan
minimisasi, ekuilibrasi NVT, ekuilibrasi NPT, dan produksi. Trajektori DCD,
log CSV, checkpoint, dan ringkasan JSON ditulis ke subdirektori ``trajectory``
dan ``logs`` di samping file konfigurasi.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

DEFAULT_CONFIG = Path(__file__).with_name("md_config.json")


def _read_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict):
        raise ValueError(f"Isi {path} bukan objek JSON.")
    return payload


def _path(base: Path, value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else base / path


def _steps(ns: float, timestep_fs: float) -> int:
    if ns <= 0:
        return 0
    return max(1, int(round(float(ns) * 1_000_000.0 / float(timestep_fs))))


def _interval_steps(ps: float, timestep_fs: float) -> int:
    if ps <= 0:
        return max(1, _steps(1.0, timestep_fs))
    return max(1, int(round(float(ps) * 1000.0 / float(timestep_fs))))


def _select_platform(requested: str):
    from openmm import Platform

    if requested and requested.lower() != "auto":
        return Platform.getPlatformByName(requested)
    last_error: Exception | None = None
    for name in ("CUDA", "OpenCL", "CPU"):
        try:
            return Platform.getPlatformByName(name)
        except Exception as exc:  # platform tidak tersedia
            last_error = exc
    raise RuntimeError(f"Tidak ada platform OpenMM yang tersedia: {last_error}")


def _platform_properties(platform) -> dict[str, str]:
    name = platform.getName()
    if name == "CUDA":
        return {"CudaPrecision": "mixed", "CudaDeviceIndex": "0"}
    if name == "OpenCL":
        return {"OpenCLPrecision": "mixed", "OpenCLDeviceIndex": "0"}
    return {}


def _write_state(path: Path, context, get_velocities: bool = True) -> None:
    from openmm import XmlSerializer

    state = context.getState(
        getPositions=True,
        getVelocities=get_velocities,
        getEnergy=True,
        enforcePeriodicBox=True,
    )
    path.write_text(XmlSerializer.serialize(state), encoding="utf-8")


def _write_pdb(path: Path, topology, positions) -> None:
    from openmm.app import PDBFile

    with path.open("w", encoding="utf-8") as handle:
        PDBFile.writeFile(topology, positions, handle, keepIds=True)


def _reporters(
    simulation,
    *,
    prefix: str,
    out_dir: Path,
    report_interval: int,
    checkpoint_interval: int,
    total_steps: int,
) -> None:
    from openmm.app import (
        CheckpointReporter,
        DCDReporter,
        StateDataReporter,
    )

    trajectory_dir = out_dir / "trajectory"
    log_dir = out_dir / "logs"
    trajectory_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)

    simulation.reporters.append(
        DCDReporter(
            str(trajectory_dir / f"{prefix}.dcd"),
            report_interval,
            enforcePeriodicBox=True,
        )
    )
    simulation.reporters.append(
        StateDataReporter(
            str(log_dir / f"{prefix}.csv"),
            report_interval,
            step=True,
            time=True,
            potentialEnergy=True,
            kineticEnergy=True,
            totalEnergy=True,
            temperature=True,
            volume=True,
            density=True,
            speed=True,
            progress=True,
            remainingTime=True,
            totalSteps=total_steps,
            separator=",",
        )
    )
    simulation.reporters.append(
        CheckpointReporter(
            str(out_dir / f"{prefix}.chk"),
            checkpoint_interval,
        )
    )


def _state_snapshot(simulation):
    return simulation.context.getState(
        getPositions=True,
        getVelocities=True,
        getEnergy=True,
        enforcePeriodicBox=True,
    )


def _set_state(simulation, state) -> None:
    simulation.context.setPositions(state.getPositions())
    try:
        velocities = state.getVelocities()
    except Exception:
        velocities = None
    if velocities is not None:
        simulation.context.setVelocities(velocities)
    try:
        box = state.getPeriodicBoxVectors()
    except Exception:
        box = None
    if box is not None:
        simulation.context.setPeriodicBoxVectors(*box)


def _make_simulation(system, topology, integrator, platform, properties):
    from openmm.app import Simulation

    return Simulation(topology, system, integrator, platform, properties)


def _minimize(simulation, max_iterations: int) -> None:
    from openmm import LocalEnergyMinimizer

    print(f"[MD] Minimisasi energi (max {max_iterations} iterasi)...", flush=True)
    LocalEnergyMinimizer.minimize(
        simulation.context, maxIterations=max_iterations
    )


def _run_phase(
    simulation,
    *,
    prefix: str,
    out_dir: Path,
    steps: int,
    report_interval: int,
    checkpoint_interval: int,
) -> None:
    if steps <= 0:
        return
    simulation.reporters.clear()
    _reporters(
        simulation,
        prefix=prefix,
        out_dir=out_dir,
        report_interval=report_interval,
        checkpoint_interval=checkpoint_interval,
        total_steps=steps,
    )
    print(f"[MD] {prefix}: {steps:,} langkah...", flush=True)
    simulation.step(steps)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Runner MD OpenMM eksplisit.")
    parser.add_argument(
        "--config",
        type=Path,
        default=DEFAULT_CONFIG,
        help="File md_config.json (default: di samping run_md.py).",
    )
    parser.add_argument("--platform", default=None, help="CUDA/OpenCL/CPU/auto.")
    parser.add_argument("--production-ns", type=float, default=None)
    parser.add_argument("--nvt-ns", type=float, default=None)
    parser.add_argument("--npt-ns", type=float, default=None)
    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument(
        "--minimize-only",
        action="store_true",
        help="Hanya minimisasi, lalu tulis state/PDB minimized.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Tampilkan rencana langkah tanpa memuat OpenMM.",
    )
    args = parser.parse_args(argv)

    config_path = args.config.resolve()
    if not config_path.is_file():
        print(f"[ERROR] Konfigurasi tidak ditemukan: {config_path}")
        return 2
    config = _read_json(config_path)
    base = config_path.parent
    md = config.get("md") or {}
    files = config.get("files") or {}

    timestep_fs = float(md.get("timestep_fs", 2.0))
    nvt_ns = float(args.nvt_ns if args.nvt_ns is not None else md.get("nvt_ns", 0.25))
    npt_ns = float(args.npt_ns if args.npt_ns is not None else md.get("npt_ns", 0.25))
    production_ns = float(
        args.production_ns
        if args.production_ns is not None
        else md.get("production_ns", 10.0)
    )
    report_interval = _interval_steps(
        float(md.get("report_interval_ps", 10.0)), timestep_fs
    )
    checkpoint_interval = _interval_steps(
        float(md.get("checkpoint_interval_ps", 100.0)), timestep_fs
    )
    nvt_steps = _steps(nvt_ns, timestep_fs)
    npt_steps = _steps(npt_ns, timestep_fs)
    production_steps = _steps(production_ns, timestep_fs)
    seed = int(args.seed if args.seed is not None else md.get("seed", 20261029))
    platform_requested = args.platform or str(md.get("platform", "auto"))

    system_path = _path(base, files.get("system_xml", "system.xml"))
    topology_path = _path(base, files.get("topology_pdb", "system.pdb"))
    state_path = _path(base, files.get("state_xml", "state.xml"))

    plan = {
        "config": str(config_path),
        "system_xml": str(system_path),
        "topology_pdb": str(topology_path),
        "state_xml": str(state_path),
        "platform_requested": platform_requested,
        "temperature_k": md.get("temperature_k", 300.0),
        "timestep_fs": timestep_fs,
        "nvt_steps": nvt_steps,
        "npt_steps": npt_steps,
        "production_steps": production_steps,
        "production_ns": round(production_steps * timestep_fs / 1_000_000.0, 4),
        "report_interval_steps": report_interval,
        "checkpoint_interval_steps": checkpoint_interval,
        "barostat": bool((md.get("barostat") or {}).get("enabled", True)),
    }
    if args.dry_run:
        print(json.dumps(plan, indent=2, ensure_ascii=False))
        return 0

    try:
        import openmm as mm
        from openmm import unit
        from openmm.app import PDBFile
    except ImportError as exc:
        print(
            "[ERROR] OpenMM belum terpasang. Di Colab jalankan:\n"
            '  python -m pip install "openmm>=8.1" "pdbfixer>=1.9" '
            '"openmmforcefields>=0.14" "openff-toolkit>=0.16" '
            '"openff-forcefields>=2024.1"\n'
            f"Detail: {exc}"
        )
        return 2

    for label, path in (
        ("system.xml", system_path),
        ("system.pdb", topology_path),
        ("state.xml", state_path),
    ):
        if not path.is_file():
            print(f"[ERROR] File {label} tidak ditemukan: {path}")
            return 2

    out_dir = base
    started = time.perf_counter()
    result: dict[str, Any] = {
        "status": "running",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "config": str(config_path),
        "plan": plan,
    }
    try:
        platform = _select_platform(platform_requested)
        properties = _platform_properties(platform)
        print(
            f"[MD] Platform: {platform.getName()} {properties}", flush=True
        )

        system = mm.XmlSerializer.deserialize(system_path.read_text(encoding="utf-8"))
        state = mm.XmlSerializer.deserialize(state_path.read_text(encoding="utf-8"))
        topology = PDBFile(str(topology_path)).topology
        if topology.getNumAtoms() != system.getNumParticles():
            raise RuntimeError(
                "Jumlah atom topology tidak cocok dengan system: "
                f"{topology.getNumAtoms()} vs {system.getNumParticles()}."
            )

        temperature = float(md.get("temperature_k", 300.0)) * unit.kelvin
        friction = float(md.get("friction_per_ps", 1.0)) / unit.picosecond
        timestep = timestep_fs * unit.femtosecond
        integrator = mm.LangevinMiddleIntegrator(temperature, friction, timestep)
        simulation = _make_simulation(
            system, topology, integrator, platform, properties
        )
        _set_state(simulation, state)
        simulation.context.setVelocitiesToTemperature(temperature, seed)

        _minimize(simulation, int(md.get("minimization_steps", 5000)))
        minimized_state = _state_snapshot(simulation)
        _write_pdb(out_dir / "minimized.pdb", topology, minimized_state.getPositions())
        _write_state(out_dir / "minimized.xml", simulation.context)

        if args.minimize_only:
            result.update(
                {
                    "status": "completed",
                    "mode": "minimize_only",
                    "elapsed_seconds": round(time.perf_counter() - started, 2),
                }
            )
            (out_dir / "md_run.json").write_text(
                json.dumps(result, indent=2, ensure_ascii=False) + "\n",
                encoding="utf-8",
            )
            print("[MD] Minimisasi selesai.")
            return 0

        if nvt_steps:
            _run_phase(
                simulation,
                prefix="nvt",
                out_dir=out_dir,
                steps=nvt_steps,
                report_interval=report_interval,
                checkpoint_interval=checkpoint_interval,
            )
        nvt_state = _state_snapshot(simulation)

        barostat_config = md.get("barostat") or {}
        barostat_enabled = bool(barostat_config.get("enabled", True))
        if barostat_enabled:
            # System dasar tidak mengandung barostat. Muat ulang lalu tambahkan
            # MonteCarloBarostat agar NVT di atas benar-benar NVT.
            npt_system = mm.XmlSerializer.deserialize(
                system_path.read_text(encoding="utf-8")
            )
            pressure = float(barostat_config.get("pressure_bar", 1.0)) * unit.bar
            frequency = int(barostat_config.get("frequency_steps", 25))
            npt_system.addForce(
                mm.MonteCarloBarostat(pressure, temperature, frequency)
            )
            npt_integrator = mm.LangevinMiddleIntegrator(
                temperature, friction, timestep
            )
            npt_simulation = _make_simulation(
                npt_system, topology, npt_integrator, platform, properties
            )
            _set_state(npt_simulation, nvt_state)
            _run_phase(
                npt_simulation,
                prefix="npt",
                out_dir=out_dir,
                steps=npt_steps,
                report_interval=report_interval,
                checkpoint_interval=checkpoint_interval,
            )
            production_simulation = npt_simulation
        else:
            production_simulation = simulation

        _run_phase(
            production_simulation,
            prefix="production",
            out_dir=out_dir,
            steps=production_steps,
            report_interval=report_interval,
            checkpoint_interval=checkpoint_interval,
        )
        final_state = _state_snapshot(production_simulation)
        _write_pdb(
            out_dir / "final.pdb",
            topology,
            final_state.getPositions(),
        )
        _write_state(out_dir / "final.xml", production_simulation.context)

        elapsed = time.perf_counter() - started
        result.update(
            {
                "status": "completed",
                "platform": platform.getName(),
                "platform_properties": properties,
                "elapsed_seconds": round(elapsed, 2),
                "nvt_steps": nvt_steps,
                "npt_steps": npt_steps,
                "production_steps": production_steps,
                "production_ns": round(
                    production_steps * timestep_fs / 1_000_000.0, 4
                ),
                "final_energy": {
                    "potential_kj_mol": float(
                        final_state.getPotentialEnergy().value_in_unit(
                            unit.kilojoule_per_mole
                        )
                    ),
                    "kinetic_kj_mol": float(
                        final_state.getKineticEnergy().value_in_unit(
                            unit.kilojoule_per_mole
                        )
                    ),
                },
            }
        )
        (out_dir / "md_run.json").write_text(
            json.dumps(result, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
        print(
            f"[MD] Selesai dalam {elapsed / 3600.0:.2f} jam; "
            f"produksi {result['production_ns']} ns."
        )
        return 0
    except Exception as exc:
        result.update(
            {
                "status": "failed",
                "error": f"{type(exc).__name__}: {exc}",
                "elapsed_seconds": round(time.perf_counter() - started, 2),
            }
        )
        (out_dir / "md_run.json").write_text(
            json.dumps(result, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
        print(f"[ERROR] MD gagal: {type(exc).__name__}: {exc}")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
'''


class PrepError(RuntimeError):
    """Kesalahan preparasi yang dapat ditampilkan tanpa traceback."""


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _relative(path: Path) -> str:
    try:
        return path.resolve().relative_to(PROJECT_ROOT.resolve()).as_posix()
    except ValueError:
        return str(path)


def _read_json(path: Path) -> dict[str, Any]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise PrepError(f"Gagal membaca {path}: {exc}") from exc
    if not isinstance(payload, dict):
        raise PrepError(f"Isi {path} bukan objek JSON.")
    return payload


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    temp.replace(path)


def _resolve(value: str | Path | None, *, must_exist: bool = True) -> Path:
    try:
        return ai.resolve_path(value, must_exist=must_exist)
    except ai.AnalysisError as exc:
        raise PrepError(str(exc)) from exc


def _first_conformer(mol: Chem.Mol) -> Chem.Mol:
    conformers = list(mol.GetConformers())
    if not conformers:
        raise PrepError("Ligan tidak memiliki conformer 3D.")
    copy = Chem.Mol(mol)
    keep = conformers[0].GetId()
    for conformer in list(copy.GetConformers()):
        if conformer.GetId() != keep:
            copy.RemoveConformer(conformer.GetId())
    return copy


def _ligand_template_path(
    payload: dict[str, Any], ligand_name: str
) -> Path | None:
    for ligand in payload.get("ligands") or []:
        if ligand.get("name") != ligand_name:
            continue
        value = ligand.get("sdf")
        if not value:
            break
        try:
            path = _resolve(value)
        except PrepError:
            break
        if path.is_file():
            return path
    fallback = PROJECT_ROOT / "data" / "ligands" / f"{ligand_name}.sdf"
    return fallback if fallback.is_file() else None


def prepare_ligand(
    pose_path: Path,
    ligand_name: str,
    *,
    template_path: Path | None = None,
) -> tuple[Chem.Mol, dict[str, Any], list[str]]:
    """Baca pose, pilih mode terbaik, tambahkan H, dan siapkan OpenFF input."""

    warnings: list[str] = []
    try:
        mol = ai.load_ligand_mol(pose_path)
    except ai.AnalysisError as exc:
        raise PrepError(str(exc)) from exc

    mol = _first_conformer(mol)
    mol = Chem.RemoveHs(mol)
    if template_path and template_path.is_file():
        template = Chem.MolFromMolFile(
            str(template_path), removeHs=True, sanitize=True
        )
        if template is None:
            warnings.append(
                f"Template SDF {template_path.name} tidak dapat dibaca; "
                "bond order dari Meeko dipakai."
            )
        else:
            try:
                mol = AllChem.AssignBondOrdersFromTemplate(template, mol)
            except Exception as exc:
                warnings.append(
                    "AssignBondOrdersFromTemplate gagal "
                    f"({type(exc).__name__}: {exc}); bond order Meeko dipakai."
                )
    else:
        warnings.append(
            "Template SDF ligan tidak ditemukan; bond order dari Meeko dipakai."
        )

    mol = Chem.AddHs(mol, addCoords=True)
    if mol.GetNumConformers() == 0:
        raise PrepError("Ligan kehilangan conformer setelah penambahan H.")
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    mol.SetProp("_Name", LIGAND_RESIDUE_NAME)
    heavy_atoms = Chem.RemoveHs(mol).GetNumAtoms()
    if heavy_atoms == 0:
        raise PrepError("Ligan tidak memiliki atom berat.")
    info = {
        "smiles": Chem.MolToSmiles(Chem.RemoveHs(mol)),
        "heavy_atoms": heavy_atoms,
        "formal_charge": Chem.GetFormalCharge(mol),
        "total_atoms_with_h": mol.GetNumAtoms(),
        "template_sdf": _relative(template_path) if template_path else None,
    }
    return mol, info, warnings


def _make_openff_ligand(mol: Chem.Mol):
    try:
        from openff.toolkit import Molecule as OffMolecule
        from openff.toolkit import Topology as OffTopology
        from openff.units.openmm import to_openmm as offquantity_to_openmm
    except ImportError as exc:
        raise PrepError(
            "OpenFF Toolkit belum terpasang. Di Colab jalankan "
            '`python -m pip install "openff-toolkit>=0.16" '
            '"openff-forcefields>=2024.1"`. '
            f"Detail: {exc}"
        ) from exc

    try:
        off_mol = OffMolecule.from_rdkit(mol, allow_undefined_stereo=True)
    except Exception as exc:
        raise PrepError(
            f"OpenFF Molecule.from_rdkit gagal: {type(exc).__name__}: {exc}"
        ) from exc
    off_mol.name = LIGAND_RESIDUE_NAME
    if off_mol.n_conformers == 0:
        raise PrepError("OpenFF Molecule tidak memiliki conformer 3D.")
    off_topology = OffTopology.from_molecules(molecules=[off_mol])
    omm_topology = off_topology.to_openmm()
    for residue in omm_topology.residues():
        residue.name = LIGAND_RESIDUE_NAME
    positions = offquantity_to_openmm(off_mol.conformers[0])
    return off_mol, omm_topology, positions


def _make_template_generator(off_mol, requested: str):
    try:
        from openmmforcefields.generators import (
            GAFFTemplateGenerator,
            SMIRNOFFTemplateGenerator,
        )
    except ImportError as exc:
        raise PrepError(
            "openmmforcefields belum terpasang. Di Colab jalankan "
            '`python -m pip install "openmmforcefields>=0.14"`. '
            f"Detail: {exc}"
        ) from exc

    if requested.lower().startswith("gaff"):
        try:
            generator = GAFFTemplateGenerator(
                molecules=[off_mol], forcefield=requested
            )
        except Exception as exc:
            raise PrepError(
                f"GAFFTemplateGenerator gagal untuk {requested}: "
                f"{type(exc).__name__}: {exc}"
            ) from exc
        return generator, requested, "gaff"

    candidates = [requested]
    if requested == "openff-2.2.1":
        candidates.extend(["openff-2.2.0", "openff-2.1.0", "openff-2.0.0"])
    errors: list[str] = []
    for candidate in candidates:
        try:
            generator = SMIRNOFFTemplateGenerator(
                molecules=[off_mol], forcefield=candidate
            )
            return generator, candidate, "openff"
        except Exception as exc:
            errors.append(f"{candidate}: {type(exc).__name__}: {exc}")
    raise PrepError(
        "SMIRNOFFTemplateGenerator gagal untuk semua kandidat force field. "
        "Pastikan `openff-forcefields` terpasang. " + " | ".join(errors)
    )


def _add_disulfide_bonds(topology, positions, *, max_distance: float = 2.3) -> int:
    """Tambahkan ikatan S-S yang tidak dibuat otomatis oleh PDBFixer."""

    import numpy as np
    from openmm import unit

    def _position(index: int) -> np.ndarray:
        value = positions[index]
        if hasattr(value, "value_in_unit"):
            value = value.value_in_unit(unit.nanometer)
        return np.asarray(value, dtype=float)

    cysteine_sg = [
        atom
        for chain in topology.chains()
        for residue in chain.residues()
        if residue.name == "CYS"
        for atom in residue.atoms()
        if atom.name == "SG"
    ]
    existing = {
        frozenset((bond[0].index, bond[1].index))
        for bond in topology.bonds()
    }
    added = 0
    for index, first in enumerate(cysteine_sg):
        for second in cysteine_sg[index + 1 :]:
            key = frozenset((first.index, second.index))
            if key in existing:
                continue
            distance = float(
                np.linalg.norm(_position(first.index) - _position(second.index))
            )
            if 1.8 <= distance <= max_distance:
                topology.addBond(first, second)
                existing.add(key)
                added += 1
    return added


def prepare_protein(
    receptor_path: Path,
    *,
    ph: float,
    fill_missing_residues: bool,
) -> tuple[Any, list[str], dict[str, Any]]:
    """PDBFixer: nonstandard residues, missing atoms, dan H pada pH 7.4."""

    try:
        from pdbfixer import PDBFixer
    except ImportError as exc:
        raise PrepError(
            "PDBFixer belum terpasang. Di Colab jalankan "
            '`python -m pip install "pdbfixer>=1.9"`. '
            f"Detail: {exc}"
        ) from exc

    warnings: list[str] = []
    try:
        fixer = PDBFixer(filename=str(receptor_path))
    except Exception as exc:
        raise PrepError(
            f"PDBFixer gagal membaca {receptor_path.name}: "
            f"{type(exc).__name__}: {exc}"
        ) from exc

    raw_text = receptor_path.read_text(encoding="utf-8", errors="replace")
    heterogens = sum(1 for line in raw_text.splitlines() if line.startswith("HETATM"))
    if heterogens:
        warnings.append(
            f"{heterogens} baris HETATM ditemukan; PDBFixer.removeHeterogens "
            "akan membuangnya. Logam/cofactor perlu ditangani terpisah."
        )

    try:
        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()
        fixer.removeHeterogens(keepWater=False)
        if fill_missing_residues:
            fixer.findMissingResidues()
            if fixer.missingResidues:
                warnings.append(
                    f"PDBFixer akan membangun {len(fixer.missingResidues)} "
                    "segmen residu yang hilang; periksa hasil loop."
                )
        else:
            fixer.missingResidues = {}
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(ph)
    except Exception as exc:
        raise PrepError(
            f"PDBFixer gagal memproses {receptor_path.name}: "
            f"{type(exc).__name__}: {exc}"
        ) from exc

    return fixer, warnings, {
        "heterogens_removed": heterogens,
        "missing_residues_filled": bool(fill_missing_residues),
        "ph": ph,
    }


def _system_stats(modeller) -> dict[str, Any]:
    topology = modeller.topology
    residue_names = [residue.name for residue in topology.residues()]
    waters = sum(
        name in {"HOH", "WAT", "H2O", "TIP3", "SOL"} for name in residue_names
    )
    ions = sum(
        name.upper() in {"NA", "CL", "NA+", "CL-", "SOD", "CLA"}
        for name in residue_names
    )
    box = topology.getPeriodicBoxVectors()
    box_vectors = None
    if box is not None:
        box_vectors = [
            [round(float(vector.x), 4), round(float(vector.y), 4), round(float(vector.z), 4)]
            for vector in box
        ]
    return {
        "atoms": topology.getNumAtoms(),
        "residues": topology.getNumResidues(),
        "waters": waters,
        "ions": ions,
        "box_vectors_nm": box_vectors,
    }


def _make_system(
    fixer,
    ligand_off,
    ligand_topology,
    ligand_positions,
    *,
    protein_forcefield: str,
    water_forcefield: str,
    ligand_forcefield: str,
    water_model: str,
    padding_nm: float,
    ionic_strength_m: float,
    nonbonded_cutoff_nm: float,
) -> tuple[Any, Any, dict[str, Any], list[str]]:
    try:
        from openmm import app, unit
        from openmm.app import ForceField, Modeller
    except ImportError as exc:
        raise PrepError(
            "OpenMM belum terpasang. Di Colab jalankan "
            '`python -m pip install "openmm>=8.1"`. '
            f"Detail: {exc}"
        ) from exc

    warnings: list[str] = []
    generator, actual_ligand_ff, ligand_family = _make_template_generator(
        ligand_off, ligand_forcefield
    )
    if actual_ligand_ff != ligand_forcefield:
        warnings.append(
            f"Force field ligan diminta {ligand_forcefield}, dipakai "
            f"{actual_ligand_ff}."
        )

    disulfides = _add_disulfide_bonds(fixer.topology, fixer.positions)
    if disulfides:
        warnings.append(f"{disulfides} ikatan disulfida S-S ditambahkan.")

    try:
        forcefield = ForceField(protein_forcefield, water_forcefield)
        forcefield.registerTemplateGenerator(generator.generator)
        modeller = Modeller(fixer.topology, fixer.positions)
        modeller.add(ligand_topology, ligand_positions)
        modeller.addSolvent(
            forcefield,
            model=water_model,
            padding=padding_nm * unit.nanometer,
            ionicStrength=ionic_strength_m * unit.molar,
            positiveIon="Na+",
            negativeIon="Cl-",
            neutralize=True,
        )
        system = forcefield.createSystem(
            modeller.topology,
            nonbondedMethod=app.PME,
            nonbondedCutoff=nonbonded_cutoff_nm * unit.nanometer,
            constraints=app.HBonds,
            rigidWater=True,
            removeCMMotion=True,
        )
    except Exception as exc:
        raise PrepError(
            "Pembentukan system OpenMM gagal: "
            f"{type(exc).__name__}: {exc}"
        ) from exc

    metadata = {
        "protein_forcefield": protein_forcefield,
        "water_forcefield": water_forcefield,
        "ligand_forcefield_requested": ligand_forcefield,
        "ligand_forcefield_used": actual_ligand_ff,
        "ligand_family": ligand_family,
        "water_model": water_model,
        "padding_nm": padding_nm,
        "ionic_strength_m": ionic_strength_m,
        "nonbonded_cutoff_nm": nonbonded_cutoff_nm,
        "disulfide_bonds_added": disulfides,
    }
    return system, modeller, metadata, warnings


def _serialize_system_and_state(
    system,
    modeller,
    *,
    system_path: Path,
    topology_path: Path,
    state_path: Path,
) -> None:
    try:
        import openmm as mm
        from openmm.app import PDBFile
    except ImportError as exc:
        raise PrepError(f"OpenMM tidak tersedia saat menulis system: {exc}") from exc

    system_path.write_text(
        mm.XmlSerializer.serialize(system), encoding="utf-8"
    )
    with topology_path.open("w", encoding="utf-8") as handle:
        PDBFile.writeFile(
            modeller.topology, modeller.positions, handle, keepIds=True
        )

    try:
        platform = mm.Platform.getPlatformByName("Reference")
    except Exception:
        platform = mm.Platform.getPlatformByName("CPU")
    integrator = mm.VerletIntegrator(0.001)
    context = mm.Context(system, integrator, platform)
    context.setPositions(modeller.positions)
    box = modeller.topology.getPeriodicBoxVectors()
    if box is not None:
        context.setPeriodicBoxVectors(*box)
    state = context.getState(getPositions=True, enforcePeriodicBox=True)
    state_path.write_text(
        mm.XmlSerializer.serialize(state), encoding="utf-8"
    )


def _md_config(
    *,
    system_name: str,
    receptor_name: str,
    ligand_name: str,
    receptor_path: Path,
    pose_path: Path,
    row: dict[str, Any],
    ligand_info: dict[str, Any],
    system_stats: dict[str, Any],
    forcefield_meta: dict[str, Any],
    prep_meta: dict[str, Any],
    args: argparse.Namespace,
) -> dict[str, Any]:
    caveats = [
        "Pose awal berasal dari docking Vina; MD menguji kestabilan pose, "
        "bukan membuktikan afinitas atau mekanisme enzim.",
        "Protein diprotonasi dengan PDBFixer pada pH 7.4; protonasi His/Glu/Asp "
        "belum dioptimalkan dengan PROPKA/H++ dan tidak ada sampling pH.",
        "Ligan DHN diperlakukan netral; pKa kedua OH naftalena sekitar 9-10, "
        "sehingga bentuk netral wajar pada pH 7.4 tetapi tetap perlu verifikasi.",
        "Force field protein Amber ff14SB, air TIP3P, dan ligan OpenFF/GAFF "
        "adalah kombinasi standar untuk eksplorasi, bukan parameter kuantum "
        "spesifik untuk DHN-melanin.",
        "Sistem tidak mengandung membran, metaloenzym cofactor, atau kondisi "
        "sitoplasma penuh; hasil hanya berlaku untuk model eksplisit sederhana.",
        "Durasi produksi default 10 ns terlalu pendek untuk klaim "
        "termodinamika; jalankan beberapa replika dan analisis RMSD/RMSF/"
        "H-bond sebelum menyimpulkan.",
        "Reseptor adalah struktur homolog surrogate, bukan struktur target "
        "C. cladosporioides.",
    ]
    if "Laccase" in receptor_name:
        caveats.insert(
            0,
            "Laccase 1GW0 clean PDB tidak mengandung tembaga T1/T2/T3. Sistem "
            "ini adalah laccase apo, sehingga pusat katalitik dan jalur transfer "
            "elektron tidak direpresentasikan.",
        )
    if "PKS1" in receptor_name:
        caveats.insert(
            0,
            "PKS1 8CG4 adalah konstruk KS-AT norsolorinic acid synthase; DHN "
            "adalah produk, bukan substrat PKS. Docking/MD ini hanya eksplorasi "
            "kualitatif kantung, bukan validasi katalisis.",
        )

    return {
        "generated_at": _now(),
        "tool": "melanin_biotech_openmm_md_prep",
        "system_name": system_name,
        "status": "prepared",
        "source": {
            "docking_file": _relative(Path(args.docking)),
            "receptor": receptor_name,
            "ligand": ligand_name,
            "receptor_pdb": _relative(receptor_path),
            "pose_pdbqt": _relative(pose_path),
            "pose_affinity_kcal_mol": row.get("best_affinity_kcal_mol"),
            "pose_ligand_efficiency": row.get("ligand_efficiency"),
            "pose_seed": row.get("best_seed") or row.get("seed"),
            "pose_seeds": row.get("seeds") or [],
        },
        "ligand": ligand_info,
        "system": system_stats,
        "forcefields": forcefield_meta,
        "preparation": prep_meta,
        "files": {
            "system_xml": "system.xml",
            "topology_pdb": "system.pdb",
            "state_xml": "state.xml",
            "runner": "run_md.py",
            "colab_runner": "run_colab.sh",
            "config": "md_config.json",
        },
        "md": {
            "temperature_k": args.temperature,
            "friction_per_ps": 1.0,
            "timestep_fs": args.timestep_fs,
            "nonbonded_cutoff_nm": args.nonbonded_cutoff_nm,
            "constraints": "HBonds",
            "rigid_water": True,
            "barostat": {
                "enabled": True,
                "pressure_bar": args.pressure_bar,
                "frequency_steps": 25,
            },
            "minimization_steps": args.minimization_steps,
            "nvt_ns": args.nvt_ns,
            "npt_ns": args.npt_ns,
            "production_ns": args.production_ns,
            "report_interval_ps": args.report_interval_ps,
            "checkpoint_interval_ps": args.checkpoint_interval_ps,
            "seed": args.seed,
            "platform": args.platform,
        },
        "caveats": caveats,
    }


def _write_runner_files(
    system_dir: Path,
    *,
    system_name: str,
    config: dict[str, Any],
) -> None:
    runner_path = system_dir / "run_md.py"
    runner_path.write_text(_RUNNER_TEMPLATE, encoding="utf-8", newline="\n")

    shell = f"""#!/usr/bin/env bash
# Runner Colab untuk {system_name}
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${{BASH_SOURCE[0]}}")" && pwd)"
ROOT="${{PROJECT_ROOT:-$(cd -- "$SCRIPT_DIR/../../.." && pwd)}}"
cd "$ROOT"

echo "============================================================"
echo " OpenMM MD: {system_name}"
echo "============================================================"
python -m pip install -q \\
  "openmm>=8.1" \\
  "pdbfixer>=1.9" \\
  "openmmforcefields>=0.14" \\
  "openff-toolkit>=0.16" \\
  "openff-forcefields>=2024.1"

if command -v nvidia-smi >/dev/null 2>&1; then
  nvidia-smi --query-gpu=name,memory.total --format=csv,noheader || true
else
  echo "[INFO] GPU tidak terdeteksi; gunakan --platform CPU atau runtime GPU."
fi

python "$SCRIPT_DIR/run_md.py" --config "$SCRIPT_DIR/md_config.json" "$@"
"""
    shell_path = system_dir / "run_colab.sh"
    shell_path.write_text(shell, encoding="utf-8", newline="\n")
    try:
        shell_path.chmod(shell_path.stat().st_mode | 0o111)
    except OSError:
        pass

    readme = f"""# MD preparation: {system_name}

Sistem ini disiapkan oleh `src/prepare_openmm_md.py` untuk MD eksplisit
OpenMM. Semua simulasi panjang harus dijalankan di Google Colab GPU, bukan
di CPU lokal.

## File

- `system.xml` — OpenMM `System` (tanpa barostat; barostat ditambahkan runner).
- `system.pdb` — topology + koordinat sistem tersolvasi.
- `state.xml` — posisi dan box awal.
- `md_config.json` — parameter MD, force field, provenance pose.
- `run_md.py` — runner minimisasi/NVT/NPT/produksi.
- `run_colab.sh` — instalasi dependency dan pemanggilan runner di Colab.

## Menjalankan di Colab

Dari root repositori:

```bash
bash results/md_prep/{system_name}/run_colab.sh
```

GPU runtime disarankan. Untuk CPU:

```bash
bash results/md_prep/{system_name}/run_colab.sh --platform CPU
```

Hasil ditulis ke `trajectory/`, `logs/`, `minimized.pdb`, `final.pdb`,
`final.xml`, dan `md_run.json` di direktori ini.

## Peringatan ilmiah

{chr(10).join("- " + item for item in config.get("caveats", []))}
"""
    (system_dir / "README.md").write_text(readme, encoding="utf-8", newline="\n")


def _system_name(receptor: str, ligand: str) -> str:
    safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", f"{receptor}__{ligand}")
    return safe


def prepare_one(
    row: dict[str, Any],
    receptor_metadata: dict[str, Any],
    payload: dict[str, Any],
    *,
    args: argparse.Namespace,
    out_dir: Path,
) -> dict[str, Any]:
    receptor_name = str(row.get("receptor"))
    ligand_name = str(row.get("ligand"))
    try:
        receptor_path = ai._receptor_path(receptor_metadata)
        pose_path = ai.find_pose_path(row)
    except (ai.AnalysisError, OSError) as exc:
        raise PrepError(str(exc)) from exc

    ligand_template = _ligand_template_path(payload, ligand_name)
    ligand_mol, ligand_info, ligand_warnings = prepare_ligand(
        pose_path,
        ligand_name,
        template_path=ligand_template,
    )
    ligand_off, ligand_topology, ligand_positions = _make_openff_ligand(
        ligand_mol
    )
    fixer, protein_warnings, prep_meta = prepare_protein(
        receptor_path,
        ph=args.ph,
        fill_missing_residues=args.fill_missing_residues,
    )
    system, modeller, forcefield_meta, system_warnings = _make_system(
        fixer,
        ligand_off,
        ligand_topology,
        ligand_positions,
        protein_forcefield=args.protein_forcefield,
        water_forcefield=args.water_forcefield,
        ligand_forcefield=args.ligand_forcefield,
        water_model=args.water_model,
        padding_nm=args.padding_nm,
        ionic_strength_m=args.ionic_strength_m,
        nonbonded_cutoff_nm=args.nonbonded_cutoff_nm,
    )
    stats = _system_stats(modeller)
    system_name = _system_name(receptor_name, ligand_name)
    system_dir = out_dir / system_name
    system_dir.mkdir(parents=True, exist_ok=True)

    system_path = system_dir / "system.xml"
    topology_path = system_dir / "system.pdb"
    state_path = system_dir / "state.xml"
    _serialize_system_and_state(
        system,
        modeller,
        system_path=system_path,
        topology_path=topology_path,
        state_path=state_path,
    )

    prep_meta = {
        **prep_meta,
        "ligand_template_sdf": ligand_info.get("template_sdf"),
        "ligand_warnings": ligand_warnings,
        "protein_warnings": protein_warnings,
        "system_warnings": system_warnings,
        "warnings": [
            *ligand_warnings,
            *protein_warnings,
            *system_warnings,
        ],
    }
    config = _md_config(
        system_name=system_name,
        receptor_name=receptor_name,
        ligand_name=ligand_name,
        receptor_path=receptor_path,
        pose_path=pose_path,
        row=row,
        ligand_info=ligand_info,
        system_stats=stats,
        forcefield_meta=forcefield_meta,
        prep_meta=prep_meta,
        args=args,
    )
    _write_json(system_dir / "md_config.json", config)
    _write_runner_files(
        system_dir, system_name=system_name, config=config
    )

    return {
        "status": "prepared",
        "system_name": system_name,
        "system_dir": _relative(system_dir),
        "receptor": receptor_name,
        "ligand": ligand_name,
        "pose_file": _relative(pose_path),
        "receptor_file": _relative(receptor_path),
        "system": stats,
        "forcefields": forcefield_meta,
        "ligand": ligand_info,
        "files": config["files"],
        "warnings": prep_meta["warnings"],
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Siapkan sistem protein-ligan untuk MD eksplisit OpenMM: "
            "PDBFixer pH 7.4, TIP3P, dan parameter ligan OpenFF/GAFF."
        )
    )
    parser.add_argument("--docking", type=Path, default=DEFAULT_DOCKING)
    parser.add_argument(
        "--receptor",
        nargs="+",
        default=list(DEFAULT_RECEPTORS),
        help="Nama reseptor (default: PKS1 dan laccase; 'all' untuk semua).",
    )
    parser.add_argument(
        "--ligand",
        nargs="+",
        default=[DEFAULT_LIGAND],
        help="Nama ligan (default: DHN_monomer; 'all' untuk semua).",
    )
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--ph", type=float, default=7.4)
    parser.add_argument("--padding-nm", type=float, default=1.0)
    parser.add_argument("--ionic-strength-m", type=float, default=0.15)
    parser.add_argument("--nonbonded-cutoff-nm", type=float, default=1.0)
    parser.add_argument("--water-model", default="tip3p")
    parser.add_argument(
        "--protein-forcefield", default="amber/protein.ff14SB.xml"
    )
    parser.add_argument(
        "--water-forcefield", default="amber/tip3p_standard.xml"
    )
    parser.add_argument(
        "--ligand-forcefield",
        default="openff-2.2.1",
        help="Nama OpenFF (mis. openff-2.2.1) atau GAFF (gaff-2.11).",
    )
    parser.add_argument("--temperature", type=float, default=300.0)
    parser.add_argument("--pressure-bar", type=float, default=1.0)
    parser.add_argument("--timestep-fs", type=float, default=2.0)
    parser.add_argument("--minimization-steps", type=int, default=5000)
    parser.add_argument("--nvt-ns", type=float, default=0.25)
    parser.add_argument("--npt-ns", type=float, default=0.25)
    parser.add_argument("--production-ns", type=float, default=10.0)
    parser.add_argument("--report-interval-ps", type=float, default=10.0)
    parser.add_argument("--checkpoint-interval-ps", type=float, default=100.0)
    parser.add_argument("--seed", type=int, default=20261029)
    parser.add_argument(
        "--platform",
        default="auto",
        help="Platform runner MD: auto/CUDA/OpenCL/CPU.",
    )
    parser.add_argument(
        "--fill-missing-residues",
        action="store_true",
        help="Biarkan PDBFixer membangun loop/residu yang hilang (risiko loop arbitrer).",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Timpa system yang sudah ada (trajectory lama tidak dihapus).",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Validasi input dan tampilkan rencana tanpa OpenMM/PDBFixer.",
    )
    parser.add_argument("--list", action="store_true")
    return parser


def _missing_runtime() -> list[str]:
    modules = {
        "openmm": "openmm>=8.1",
        "pdbfixer": "pdbfixer>=1.9",
        "openmmforcefields": "openmmforcefields>=0.14",
        "openff.toolkit": "openff-toolkit>=0.16",
        "openff.units": "openff-units",
    }
    missing = []
    for module, package in modules.items():
        try:
            found = importlib.util.find_spec(module) is not None
        except (ImportError, ValueError):
            found = False
        if not found:
            missing.append(package)
    return missing


def _dry_run_one(
    row: dict[str, Any],
    receptor_metadata: dict[str, Any],
    payload: dict[str, Any],
    *,
    args: argparse.Namespace,
    out_dir: Path,
) -> dict[str, Any]:
    receptor_name = str(row.get("receptor"))
    ligand_name = str(row.get("ligand"))
    try:
        receptor_path = ai._receptor_path(receptor_metadata)
        pose_path = ai.find_pose_path(row)
    except (ai.AnalysisError, OSError) as exc:
        raise PrepError(str(exc)) from exc
    ligand_mol, ligand_info, warnings = prepare_ligand(
        pose_path,
        ligand_name,
        template_path=_ligand_template_path(payload, ligand_name),
    )
    return {
        "status": "dry_run",
        "system_name": _system_name(receptor_name, ligand_name),
        "system_dir": _relative(out_dir / _system_name(receptor_name, ligand_name)),
        "receptor": receptor_name,
        "ligand": ligand_name,
        "receptor_file": _relative(receptor_path),
        "pose_file": _relative(pose_path),
        "receptor_size_bytes": receptor_path.stat().st_size,
        "ligand": ligand_info,
        "warnings": warnings,
        "planned_forcefields": {
            "protein": args.protein_forcefield,
            "water": args.water_forcefield,
            "ligand": args.ligand_forcefield,
        },
        "planned_box": {
            "padding_nm": args.padding_nm,
            "ionic_strength_m": args.ionic_strength_m,
            "water_model": args.water_model,
        },
    }


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        payload = _read_json(args.docking)
        rows = ai._select_rows(payload, args.receptor, args.ligand)
    except (PrepError, ai.AnalysisError) as exc:
        print(f"[ERROR] {exc}")
        return 1

    if args.list:
        for row in rows:
            try:
                metadata = ai._receptor_metadata(
                    payload, str(row.get("receptor", ""))
                )
                receptor_path = ai._receptor_path(metadata)
                pose_path = ai.find_pose_path(row)
                print(
                    f"{row.get('receptor')} x {row.get('ligand')}: "
                    f"pose={_relative(pose_path)}; receptor={_relative(receptor_path)}"
                )
            except (PrepError, ai.AnalysisError, OSError) as exc:
                print(
                    f"[WARN] {row.get('receptor')} x {row.get('ligand')}: {exc}"
                )
        return 0

    out_dir = args.out_dir.resolve()
    missing = _missing_runtime()
    if args.dry_run:
        print(
            "[INFO] Dry run: tidak menulis system dan tidak menjalankan "
            "PDBFixer/OpenMM."
        )
        if missing:
            print(
                "[INFO] Dependency runtime belum ada: "
                + ", ".join(missing)
                + " (pasang di Colab)."
            )
    elif missing:
        print(
            "[ERROR] Dependency runtime belum terpasang: "
            + ", ".join(missing)
            + ".\n        Di Colab jalankan:\n"
            '        python -m pip install "openmm>=8.1" "pdbfixer>=1.9" '
            '"openmmforcefields>=0.14" "openff-toolkit>=0.16" '
            '"openff-forcefields>=2024.1"'
        )
        return 2

    prepared: list[dict[str, Any]] = []
    errors: list[dict[str, Any]] = []
    skipped: list[dict[str, Any]] = []
    for row in rows:
        receptor_name = str(row.get("receptor"))
        ligand_name = str(row.get("ligand"))
        system_name = _system_name(receptor_name, ligand_name)
        system_dir = out_dir / system_name
        if (
            not args.force
            and not args.dry_run
            and (system_dir / "system.xml").is_file()
        ):
            skipped.append(
                {
                    "system_name": system_name,
                    "system_dir": _relative(system_dir),
                    "reason": "system.xml sudah ada; gunakan --force untuk menimpa.",
                }
            )
            print(f"[SKIP] {system_name}: system.xml sudah ada.")
            continue

        metadata = ai._receptor_metadata(payload, receptor_name)
        try:
            if args.dry_run:
                result = _dry_run_one(
                    row,
                    metadata,
                    payload,
                    args=args,
                    out_dir=out_dir,
                )
            else:
                print(f"\n[PREP] {system_name}")
                result = prepare_one(
                    row,
                    metadata,
                    payload,
                    args=args,
                    out_dir=out_dir,
                )
            prepared.append(result)
            print(f"[OK] {system_name}: {result['system_dir']}")
            for warning in result.get("warnings", []):
                print(f"    [WARN] {warning}")
        except (PrepError, OSError, ValueError, KeyError, TypeError) as exc:
            errors.append(
                {
                    "receptor": receptor_name,
                    "ligand": ligand_name,
                    "system_name": system_name,
                    "error": str(exc),
                }
            )
            print(f"[ERROR] {system_name}: {exc}")

    status = "completed"
    if errors and (prepared or skipped):
        status = "partial"
    elif errors and not prepared and not skipped:
        status = "failed"

    summary = {
        "generated_at": _now(),
        "tool": "melanin_biotech_openmm_md_prep",
        "status": status,
        "dry_run": bool(args.dry_run),
        "docking_file": _relative(args.docking),
        "out_dir": _relative(out_dir),
        "parameters": {
            "ph": args.ph,
            "padding_nm": args.padding_nm,
            "ionic_strength_m": args.ionic_strength_m,
            "water_model": args.water_model,
            "ligand_forcefield": args.ligand_forcefield,
            "protein_forcefield": args.protein_forcefield,
            "water_forcefield": args.water_forcefield,
            "temperature_k": args.temperature,
            "production_ns": args.production_ns,
        },
        "systems": prepared,
        "skipped": skipped,
        "errors": errors,
    }
    if not args.dry_run:
        _write_json(out_dir / "md_prep_summary.json", summary)

    print(
        f"\n[SUMMARY] prepared={len(prepared)}, skipped={len(skipped)}, "
        f"errors={len(errors)}, status={status}"
    )
    if not args.dry_run:
        print(f"    Summary: {out_dir / 'md_prep_summary.json'}")
    if args.dry_run:
        print(json.dumps(summary, indent=2, ensure_ascii=False))
    return 0 if status != "failed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
