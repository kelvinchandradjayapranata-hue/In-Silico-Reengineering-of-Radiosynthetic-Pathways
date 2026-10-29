"""Klastering pose AutoDock Vina berdasarkan RMSD ligan.

Modul ini membaca semua pose dari ``results/docking_results.json`` (single-seed
maupun agregasi multi-seed), lalu mengelompokkan pose dengan ambang RMSD
default 2.0 Å. Untuk pose docking, RMSD dihitung **di frame reseptor tanpa
superposisi** dan dikoreksi terhadap simetri molekul. Ini penting: ligan rigid
seperti DHN monomer akan selalu memberi RMSD ~0 bila pose di-superposiskan
satu sama lain, sehingga superposisi tidak dapat dipakai untuk memisahkan
pose yang berbeda di kantung ikat.

Keluaran:

- ``results/pose_clusters.json`` — detail klaster dan anggota,
- ``results/pose_clusters.tsv`` — satu baris per klaster,
- ``results/pose_clusters.md`` — ringkasan yang dapat dibaca manusia.

Catatan: Vina adalah alat docking, bukan penentu energi bebas ikat. Populasi
klaster menunjukkan konsistensi pencarian pose, bukan populasi Boltzmann
eksperimental.
"""

from __future__ import annotations

import argparse
import json
import math
import re
import statistics
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np
from meeko import PDBQTMolecule, RDKitMolCreate
from rdkit import Chem
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DOCKING = PROJECT_ROOT / "results" / "docking_results.json"
DEFAULT_JSON = PROJECT_ROOT / "results" / "pose_clusters.json"
DEFAULT_TSV = PROJECT_ROOT / "results" / "pose_clusters.tsv"
DEFAULT_MD = PROJECT_ROOT / "results" / "pose_clusters.md"
POSE_DIR = PROJECT_ROOT / "results" / "docking" / "poses"

_MODEL_RE = re.compile(r"^MODEL\s+(\d+)\s*$", re.MULTILINE)
_VINA_RESULT_RE = re.compile(
    r"^REMARK VINA RESULT:\s+"
    r"(-?\d+(?:\.\d+)?)\s+"
    r"(-?\d+(?:\.\d+)?)\s+"
    r"(-?\d+(?:\.\d+)?)",
    re.MULTILINE,
)
_SEED_RE = re.compile(r"__seed(-?\d+)(?:\.pdbqt)?$", re.IGNORECASE)


class ClusterError(RuntimeError):
    """Kesalahan yang dapat ditampilkan tanpa traceback."""


@dataclass
class PoseRecord:
    receptor: str
    ligand: str
    seed: Any
    model: int
    conformer_index: int
    affinity: float | None
    heavy_atoms: int
    ligand_efficiency: float | None
    pose_file: Path
    source: str
    mol: Chem.Mol = field(repr=False, compare=False)
    symmetry_maps: list[tuple[int, ...]] = field(
        default_factory=list, repr=False, compare=False
    )

    @property
    def label(self) -> str:
        seed = "?" if self.seed is None else str(self.seed)
        return f"seed{seed}/model{self.model}"


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _relative(path: Path) -> str:
    try:
        return path.resolve().relative_to(PROJECT_ROOT.resolve()).as_posix()
    except ValueError:
        return str(path)


def resolve_path(value: str | Path | None, *, must_exist: bool = True) -> Path:
    if value is None or str(value).strip() == "":
        raise ClusterError("Path tidak diberikan.")
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
        raise ClusterError(f"File tidak ditemukan: {value}")
    return candidates[0]


def _read_json(path: Path) -> dict[str, Any]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ClusterError(f"Gagal membaca {path}: {exc}") from exc
    if not isinstance(payload, dict):
        raise ClusterError(f"Isi {path} bukan objek JSON.")
    return payload


def _write_json_atomic(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp_path = path.with_suffix(path.suffix + ".tmp")
    temp_path.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    temp_path.replace(path)


def _seed_from_path(path: Path) -> int | None:
    match = _SEED_RE.search(path.name)
    return int(match.group(1)) if match else None


def parse_vina_models(path: Path) -> list[dict[str, Any]]:
    """Ambil nomor model dan afinitas dari setiap blok MODEL PDBQT.

    Jumlah model di PDBQT adalah sumber kebenaran untuk pose yang benar-benar
    tersedia. Field ``poses`` di JSON dapat berasal dari log yang berbeda atau
    terpotong, sehingga hanya dipakai sebagai fallback bila afinitas PDBQT
    tidak lengkap.
    """

    try:
        text = path.read_text(encoding="utf-8", errors="replace")
    except OSError as exc:
        raise ClusterError(f"Gagal membaca pose {path}: {exc}") from exc

    starts = list(_MODEL_RE.finditer(text))
    result: list[dict[str, Any]] = []
    if starts:
        for index, match in enumerate(starts):
            block_start = match.end()
            block_end = (
                starts[index + 1].start() if index + 1 < len(starts) else len(text)
            )
            block = text[block_start:block_end]
            affinity_match = _VINA_RESULT_RE.search(block)
            result.append(
                {
                    "model": int(match.group(1)),
                    "affinity": float(affinity_match.group(1))
                    if affinity_match
                    else None,
                }
            )
        return result

    affinity_match = _VINA_RESULT_RE.search(text)
    return [
        {
            "model": 1,
            "affinity": float(affinity_match.group(1)) if affinity_match else None,
        }
    ]


def load_ligand_mol(path: Path) -> Chem.Mol:
    try:
        pdbqt = PDBQTMolecule.from_file(str(path), skip_typing=True)
        mols = RDKitMolCreate.from_pdbqt_mol(pdbqt)
    except Exception as exc:
        raise ClusterError(f"Meeko gagal mengonversi {path.name}: {exc}") from exc
    if not mols:
        raise ClusterError(f"Tidak ada molekul ligan dari {path.name}.")
    mol = mols[0]
    if mol.GetNumConformers() == 0:
        raise ClusterError(f"Pose {path.name} tidak memiliki conformer.")
    return mol


def _single_conformer_heavy(mol: Chem.Mol, conformer_index: int) -> Chem.Mol:
    """Salin satu conformer, buang H, dan pertahankan koordinat 3D."""

    conformers = list(mol.GetConformers())
    if not 0 <= conformer_index < len(conformers):
        raise ClusterError(
            f"Conformer index {conformer_index} di luar rentang "
            f"0..{len(conformers) - 1}."
        )
    copy = Chem.Mol(mol)
    keep_id = conformers[conformer_index].GetId()
    for conformer in list(copy.GetConformers()):
        if conformer.GetId() != keep_id:
            copy.RemoveConformer(conformer.GetId())
    return Chem.RemoveHs(copy)


def _symmetry_maps(
    reference: Chem.Mol, probe: Chem.Mol
) -> list[tuple[int, ...]]:
    """Peta atom probe -> atom reference untuk semua automorfisme graph."""

    if probe.GetNumAtoms() != reference.GetNumAtoms():
        return []
    try:
        matches = reference.GetSubstructMatches(
            probe,
            uniquify=False,
            useChirality=False,
            maxMatches=100000,
        )
    except Exception:
        return []
    return [tuple(int(index) for index in match) for match in matches]


def _direct_rmsd(
    probe: Chem.Mol,
    reference: Chem.Mol,
    maps: Sequence[Sequence[int]],
) -> float:
    """RMSD koordinat tanpa superposisi, memakai peta simetri terbaik."""

    if not maps:
        return float("inf")
    probe_positions = np.asarray(probe.GetConformer().GetPositions(), dtype=float)
    reference_positions = np.asarray(
        reference.GetConformer().GetPositions(), dtype=float
    )
    best = float("inf")
    for atom_map in maps:
        if len(atom_map) != len(probe_positions):
            continue
        diff = probe_positions - reference_positions[list(atom_map)]
        value = math.sqrt(float(np.mean(np.sum(diff * diff, axis=1))))
        if value < best:
            best = value
    return best


def _aligned_rmsd(probe: Chem.Mol, reference: Chem.Mol) -> float:
    """RMSD setelah superposisi; hanya untuk mode pembanding internal."""

    from rdkit.Chem import rdMolAlign

    try:
        return float(rdMolAlign.GetBestRMS(Chem.Mol(probe), Chem.Mol(reference)))
    except Exception:
        return float("inf")


def _collect_pose_files(
    row: dict[str, Any],
) -> list[tuple[Path, Any, str]]:
    """Kumpulkan file pose unik beserta seed dan asal metadata."""

    receptor = str(row.get("receptor", ""))
    ligand = str(row.get("ligand", ""))
    entries: list[tuple[Path, Any, str]] = []
    seen: set[str] = set()

    def add(value: Any, seed: Any, source: str) -> None:
        if not value:
            return
        try:
            path = resolve_path(value)
        except ClusterError:
            return
        if not path.is_file():
            return
        key = str(path.resolve())
        if key in seen:
            return
        seen.add(key)
        entries.append((path, seed, source))

    per_seed = [item for item in (row.get("per_seed") or []) if isinstance(item, dict)]
    for item in per_seed:
        add(item.get("pose_file"), item.get("seed"), "per_seed")

    add(row.get("pose_file"), row.get("seed") or row.get("best_seed"), "row")

    # Fallback untuk agregasi lama yang tidak menulis per_seed lengkap.
    patterns = [
        POSE_DIR / f"{receptor}__{ligand}__seed*.pdbqt",
        POSE_DIR / f"{receptor}__{ligand}.pdbqt",
    ]
    for pattern in patterns:
        for path in sorted(POSE_DIR.glob(pattern.name)):
            add(path, _seed_from_path(path), "glob")

    return entries


def _json_pose_affinities(row: dict[str, Any]) -> dict[int, float]:
    result: dict[int, float] = {}
    for item in row.get("poses") or []:
        if not isinstance(item, dict):
            continue
        try:
            mode = int(item.get("mode"))
            affinity = float(item.get("affinity_kcal_mol"))
        except (TypeError, ValueError):
            continue
        result[mode] = affinity
    return result


def build_pose_records(
    row: dict[str, Any],
    *,
    fallback_seed: Any = None,
) -> tuple[list[PoseRecord], list[str]]:
    receptor = str(row.get("receptor", ""))
    ligand = str(row.get("ligand", ""))
    warnings: list[str] = []
    sources = _collect_pose_files(row)
    if not sources:
        raise ClusterError(
            f"Tidak ada file pose PDBQT untuk {receptor} x {ligand}."
        )

    json_affinities = _json_pose_affinities(row)
    records: list[PoseRecord] = []
    reference_mol: Chem.Mol | None = None
    json_count_checked = False

    for path, seed, source in sources:
        if seed is None:
            seed = _seed_from_path(path)
        if seed is None:
            seed = fallback_seed
        models = parse_vina_models(path)
        if (
            not json_count_checked
            and source == "row"
            and json_affinities
            and len(json_affinities) != len(models)
        ):
            warnings.append(
                f"{path.name}: JSON mencantumkan {len(json_affinities)} mode "
                f"tetapi PDBQT berisi {len(models)} model; klaster memakai "
                "model PDBQT aktual."
            )
        if source == "row":
            json_count_checked = True
        mol = load_ligand_mol(path)
        n_conformers = mol.GetNumConformers()
        if len(models) != n_conformers:
            warnings.append(
                f"{path.name}: {len(models)} model PDBQT tetapi "
                f"{n_conformers} conformer Meeko; dipakai "
                f"{min(len(models), n_conformers)} pose."
            )
        heavy = Chem.RemoveHs(mol)
        heavy_atoms = heavy.GetNumAtoms()
        if heavy_atoms == 0:
            warnings.append(f"{path.name}: tidak ada atom berat; dilewati.")
            continue

        if reference_mol is None:
            reference_mol = heavy

        for index, model_info in enumerate(models):
            if index >= n_conformers:
                break
            model = int(model_info.get("model", index + 1))
            affinity = model_info.get("affinity")
            if affinity is None:
                affinity = json_affinities.get(model)
            if affinity is None:
                warnings.append(
                    f"{path.name} model {model}: afinitas tidak tersedia."
                )
            single = _single_conformer_heavy(mol, index)
            maps = _symmetry_maps(reference_mol, single)
            if not maps:
                warnings.append(
                    f"{path.name} model {model}: graph/atom order tidak cocok "
                    "dengan pose referensi; RMSD pasangan ini akan diabaikan."
                )
            ligand_efficiency = (
                round(float(affinity) / heavy_atoms, 3)
                if affinity is not None and heavy_atoms
                else None
            )
            records.append(
                PoseRecord(
                    receptor=receptor,
                    ligand=ligand,
                    seed=seed,
                    model=model,
                    conformer_index=index,
                    affinity=float(affinity) if affinity is not None else None,
                    heavy_atoms=heavy_atoms,
                    ligand_efficiency=ligand_efficiency,
                    pose_file=path,
                    source=source,
                    mol=single,
                    symmetry_maps=maps,
                )
            )

    if not records:
        raise ClusterError(
            f"Tidak ada pose valid untuk {receptor} x {ligand}."
        )
    return records, warnings


def _pairwise_rmsd(
    records: Sequence[PoseRecord],
    *,
    mode: str,
) -> np.ndarray:
    n = len(records)
    matrix = np.zeros((n, n), dtype=float)
    valid = [bool(record.symmetry_maps) for record in records]
    for i in range(n):
        for j in range(i + 1, n):
            if mode == "aligned":
                value = _aligned_rmsd(records[i].mol, records[j].mol)
            elif valid[i] and valid[j]:
                value = _direct_rmsd(
                    records[i].mol,
                    records[j].mol,
                    records[i].symmetry_maps,
                )
            else:
                value = float("inf")
            matrix[i, j] = matrix[j, i] = value
    return matrix


def _cluster_labels(
    matrix: np.ndarray, *, cutoff: float, method: str
) -> np.ndarray:
    n = matrix.shape[0]
    if n == 0:
        return np.array([], dtype=int)
    if n == 1:
        return np.ones(1, dtype=int)
    finite = np.where(np.isfinite(matrix), matrix, 1e6)
    condensed = squareform(finite, checks=False)
    if condensed.size == 0:
        return np.ones(n, dtype=int)
    tree = linkage(condensed, method=method)
    return fcluster(tree, t=cutoff, criterion="distance")


def _mean(values: Sequence[float]) -> float | None:
    return round(float(statistics.mean(values)), 3) if values else None


def _sd(values: Sequence[float]) -> float | None:
    if not values:
        return None
    if len(values) == 1:
        return 0.0
    return round(float(statistics.stdev(values)), 3)


def _member_payload(
    record: PoseRecord,
    *,
    rmsd_to_representative: float,
) -> dict[str, Any]:
    return {
        "label": record.label,
        "seed": record.seed,
        "model": record.model,
        "conformer_index": record.conformer_index,
        "affinity_kcal_mol": record.affinity,
        "ligand_efficiency": record.ligand_efficiency,
        "heavy_atoms": record.heavy_atoms,
        "pose_file": _relative(record.pose_file),
        "source": record.source,
        "rmsd_to_representative_angstrom": (
            round(rmsd_to_representative, 3)
            if math.isfinite(rmsd_to_representative)
            else None
        ),
    }


def _build_clusters(
    records: Sequence[PoseRecord],
    matrix: np.ndarray,
    labels: np.ndarray,
    *,
    cutoff: float,
    method: str,
    rmsd_mode: str,
) -> list[dict[str, Any]]:
    groups: dict[int, list[int]] = {}
    for index, label in enumerate(labels.tolist()):
        groups.setdefault(int(label), []).append(index)

    clusters: list[dict[str, Any]] = []
    for raw_label, indices in groups.items():
        affinities = [
            records[index].affinity
            for index in indices
            if records[index].affinity is not None
        ]
        efficiencies = [
            records[index].ligand_efficiency
            for index in indices
            if records[index].ligand_efficiency is not None
        ]
        best_index = min(
            indices,
            key=lambda index: (
                records[index].affinity
                if records[index].affinity is not None
                else float("inf"),
                records[index].seed if records[index].seed is not None else 10**9,
                records[index].model,
            ),
        )
        representative = records[best_index]

        # Medoid = anggota dengan total RMSD terkecil ke anggota lain.
        medoid_index = min(
            indices,
            key=lambda index: (
                sum(
                    matrix[index, other]
                    for other in indices
                    if other != index and math.isfinite(matrix[index, other])
                ),
                records[index].affinity
                if records[index].affinity is not None
                else float("inf"),
            ),
        )

        intra: list[float] = []
        for pos, i in enumerate(indices):
            for j in indices[pos + 1 :]:
                value = float(matrix[i, j])
                if math.isfinite(value):
                    intra.append(value)
        max_intra = max(intra) if intra else 0.0
        mean_intra = _mean(intra) if intra else 0.0
        members = []
        for index in sorted(
            indices,
            key=lambda item: (
                records[item].affinity
                if records[item].affinity is not None
                else float("inf"),
                records[item].seed if records[item].seed is not None else 10**9,
                records[item].model,
            ),
        ):
            members.append(
                _member_payload(
                    records[index],
                    rmsd_to_representative=matrix[best_index, index],
                )
            )

        seeds = sorted(
            {
                records[index].seed
                for index in indices
                if records[index].seed is not None
            }
        )
        clusters.append(
            {
                "raw_label": raw_label,
                "population": len(indices),
                "fraction": round(len(indices) / len(records), 4),
                "seeds": seeds,
                "representative": _member_payload(
                    representative, rmsd_to_representative=0.0
                ),
                "medoid": _member_payload(
                    records[medoid_index],
                    rmsd_to_representative=matrix[best_index, medoid_index],
                ),
                "best_affinity_kcal_mol": (
                    round(min(affinities), 3) if affinities else None
                ),
                "mean_affinity_kcal_mol": _mean(affinities),
                "sd_affinity_kcal_mol": _sd(affinities),
                "best_ligand_efficiency": (
                    min(efficiencies) if efficiencies else None
                ),
                "mean_ligand_efficiency": _mean(efficiencies),
                "sd_ligand_efficiency": _sd(efficiencies),
                "max_intra_rmsd_angstrom": round(max_intra, 3),
                "mean_intra_rmsd_angstrom": mean_intra,
                "members": members,
            }
        )

    # Urutkan klaster: afinitas terbaik (paling negatif), lalu populasi.
    clusters.sort(
        key=lambda item: (
            item["best_affinity_kcal_mol"]
            if item["best_affinity_kcal_mol"] is not None
            else float("inf"),
            -item["population"],
            item["representative"]["seed"]
            if item["representative"]["seed"] is not None
            else 10**9,
        )
    )
    for cluster_id, cluster in enumerate(clusters, start=1):
        cluster["cluster_id"] = cluster_id
        cluster["selection"] = (
            "best_affinity_member"
            if cluster["population"] > 1
            else "singleton"
        )

    if method == "average":
        for cluster in clusters:
            cluster["linkage_note"] = (
                "Average linkage tidak menjamin semua pasangan anggota "
                "berada <= cutoff; lihat max_intra_rmsd."
            )
    return clusters


def cluster_complex(
    row: dict[str, Any],
    *,
    cutoff: float,
    method: str,
    rmsd_mode: str,
    fallback_seed: Any = None,
) -> dict[str, Any]:
    records, warnings = build_pose_records(row, fallback_seed=fallback_seed)
    matrix = _pairwise_rmsd(records, mode=rmsd_mode)
    if rmsd_mode == "receptor":
        finite = matrix[np.isfinite(matrix) & (matrix > 0)]
        if finite.size == 0:
            warnings.append(
                "Semua RMSD pasangan ~0. Untuk ligan rigid, ini bisa berarti "
                "pose benar-benar identik; untuk ligan fleksibel, periksa "
                "kembali konversi PDBQT."
            )
    labels = _cluster_labels(matrix, cutoff=cutoff, method=method)
    clusters = _build_clusters(
        records,
        matrix,
        labels,
        cutoff=cutoff,
        method=method,
        rmsd_mode=rmsd_mode,
    )
    heavy_atoms = sorted({record.heavy_atoms for record in records})
    seeds = sorted(
        {record.seed for record in records if record.seed is not None}
    )
    best_cluster = clusters[0] if clusters else None
    return {
        "status": "ok",
        "receptor": str(row.get("receptor", "")),
        "ligand": str(row.get("ligand", "")),
        "pose_files": sorted({_relative(record.pose_file) for record in records}),
        "n_poses": len(records),
        "n_seeds": len(seeds),
        "seeds": seeds,
        "heavy_atoms": heavy_atoms,
        "n_clusters": len(clusters),
        "clusters": clusters,
        "best_cluster_id": best_cluster.get("cluster_id") if best_cluster else None,
        "best_affinity_kcal_mol": (
            best_cluster.get("best_affinity_kcal_mol") if best_cluster else None
        ),
        "best_ligand_efficiency": (
            best_cluster.get("best_ligand_efficiency") if best_cluster else None
        ),
        "warnings": warnings,
    }


def _select_rows(
    payload: dict[str, Any],
    receptors: Sequence[str] | None,
    ligands: Sequence[str] | None,
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
    if receptors and "all" not in {name.lower() for name in receptors}:
        selected = set(receptors)
        rows = [row for row in rows if row.get("receptor") in selected]
    if ligands and "all" not in {name.lower() for name in ligands}:
        selected = set(ligands)
        rows = [row for row in rows if row.get("ligand") in selected]
    if not rows:
        raise ClusterError(
            "Tidak ada kombinasi reseptor x ligan yang berhasil dan cocok "
            "dengan filter."
        )
    return rows


def build_payload(
    docking_path: Path = DEFAULT_DOCKING,
    *,
    receptors: Sequence[str] | None = None,
    ligands: Sequence[str] | None = None,
    cutoff: float = 2.0,
    method: str = "complete",
    rmsd_mode: str = "receptor",
) -> dict[str, Any]:
    payload = _read_json(docking_path)
    rows = _select_rows(payload, receptors, ligands)
    fallback_seed = (payload.get("parameters") or {}).get("seed")
    complexes: list[dict[str, Any]] = []
    errors: list[dict[str, Any]] = []

    for row in rows:
        try:
            complexes.append(
                cluster_complex(
                    row,
                    cutoff=cutoff,
                    method=method,
                    rmsd_mode=rmsd_mode,
                    fallback_seed=fallback_seed,
                )
            )
        except (ClusterError, OSError, ValueError, KeyError, TypeError) as exc:
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

    aggregation = (payload.get("summary") or {}).get("aggregation") or {}
    mode = "multi_seed" if aggregation.get("method") == "multi_seed" else "single_seed"
    # Mode juga dapat disimpulkan dari jumlah file per kompleks, tetapi
    # metadata agregator tetap menjadi sumber utama.
    return {
        "generated_at": _now(),
        "tool": "melanin_biotech_pose_clustering",
        "status": status,
        "docking_file": _relative(docking_path),
        "docking_status": payload.get("status"),
        "docking_mode": mode,
        "docking_aggregation": aggregation or None,
        "parameters": {
            "rmsd_cutoff_angstrom": cutoff,
            "linkage": method,
            "rmsd_mode": rmsd_mode,
            "receptor_filter": list(receptors or []),
            "ligand_filter": list(ligands or []),
            "heavy_atom_rmsd": True,
            "rmsd_note": (
                "RMSD dihitung di frame reseptor tanpa superposisi dan "
                "dikoreksi simetri; bukan RMSD internal setelah alignment."
            ),
        },
        "complexes": complexes,
        "errors": errors,
        "summary": {
            "complexes_total": len(complexes),
            "complexes_failed": len(errors),
            "poses_total": sum(item["n_poses"] for item in complexes),
            "clusters_total": sum(item["n_clusters"] for item in complexes),
            "largest_cluster_population": max(
                (
                    cluster["population"]
                    for item in complexes
                    for cluster in item["clusters"]
                ),
                default=0,
            ),
        },
        "caveats": [
            "RMSD default dihitung pada koordinat pose di frame reseptor "
            "tanpa superposisi, dengan koreksi automorfisme graph. Ini "
            "memisahkan orientasi/translasi pose di kantung; superposisi "
            "ligan rigid akan memberi RMSD ~0 untuk semua pose.",
            "Complete linkage menjamin setiap pasangan dalam klaster "
            "berada <= cutoff; average linkage tidak menjamin hal itu.",
            "Populasi klaster adalah frekuensi pose di antara mode Vina "
            "yang dijalankan, bukan populasi Boltzmann atau bukti "
            "kestabilan termodinamika.",
            "Nilai Vina bukan energi bebas ikat; klaster terbaik belum "
            "tentu pose yang paling stabil dalam MD.",
            "Reseptor adalah struktur homolog surrogate, bukan struktur "
            "target C. cladosporioides.",
        ],
    }


def _tsv_rows(payload: dict[str, Any]) -> Iterable[list[str]]:
    for complex_data in payload.get("complexes", []):
        for cluster in complex_data.get("clusters", []):
            representative = cluster.get("representative") or {}
            yield [
                f"{complex_data['receptor']}__{complex_data['ligand']}",
                complex_data["receptor"],
                complex_data["ligand"],
                str(cluster["cluster_id"]),
                str(cluster["population"]),
                f"{cluster['fraction']:.4f}",
                str(representative.get("seed", "")),
                str(representative.get("model", "")),
                _format_float(cluster.get("best_affinity_kcal_mol")),
                _format_float(cluster.get("mean_affinity_kcal_mol")),
                _format_float(cluster.get("sd_affinity_kcal_mol")),
                _format_float(cluster.get("best_ligand_efficiency")),
                _format_float(cluster.get("mean_ligand_efficiency")),
                _format_float(cluster.get("max_intra_rmsd_angstrom")),
                _format_float(cluster.get("mean_intra_rmsd_angstrom")),
                ",".join(str(seed) for seed in cluster.get("seeds", [])),
                representative.get("pose_file", ""),
            ]


def _format_float(value: Any) -> str:
    return f"{float(value):.3f}" if value is not None else ""


def write_tsv(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    header = [
        "complex",
        "receptor",
        "ligand",
        "cluster_id",
        "population",
        "fraction",
        "representative_seed",
        "representative_model",
        "best_affinity_kcal_mol",
        "mean_affinity_kcal_mol",
        "sd_affinity_kcal_mol",
        "best_ligand_efficiency",
        "mean_ligand_efficiency",
        "max_intra_rmsd_angstrom",
        "mean_intra_rmsd_angstrom",
        "seeds",
        "representative_pose_file",
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


def _table(headers: Sequence[str], rows: Sequence[Sequence[Any]]) -> str:
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
        "# Klastering Pose Docking (RMSD)",
        "",
        f"> Dibuat: {payload['generated_at']}",
        f"> Mode docking: {payload.get('docking_mode', 'unknown')}",
        f"> Cutoff RMSD: {payload['parameters']['rmsd_cutoff_angstrom']:.2f} Å; "
        f"linkage: {payload['parameters']['linkage']}",
        ">",
        "> RMSD dihitung di frame reseptor tanpa superposisi dan dikoreksi "
        "simetri. Populasi klaster bukan populasi Boltzmann.",
        "",
        "## Ringkasan",
        "",
    ]
    summary_rows = []
    for item in payload.get("complexes", []):
        summary_rows.append(
            [
                item["receptor"],
                item["ligand"],
                item["n_poses"],
                item["n_seeds"],
                item["n_clusters"],
                _format_float(item.get("best_affinity_kcal_mol")),
                _format_float(item.get("best_ligand_efficiency")),
                max(
                    (cluster["population"] for cluster in item["clusters"]),
                    default=0,
                ),
            ]
        )
    if summary_rows:
        lines.extend(
            [
                _table(
                    [
                        "Reseptor",
                        "Ligan",
                        "Pose",
                        "Jumlah seed",
                        "Klaster",
                        "Afinitas terbaik",
                        "LE terbaik",
                        "Populasi terbesar",
                    ],
                    summary_rows,
                ),
                "",
            ]
        )
    else:
        lines.extend(["Tidak ada kompleks yang berhasil diklaster.", ""])

    for item in payload.get("complexes", []):
        lines.extend(
            [
                f"## {item['receptor']} — {item['ligand']}",
                "",
                f"- Pose: {item['n_poses']} mode dari {item['n_seeds']} seed",
                f"- Atom berat: {', '.join(str(value) for value in item['heavy_atoms']) or '-'}",
                f"- Klaster: {item['n_clusters']}",
                "",
            ]
        )
        cluster_rows = []
        for cluster in item["clusters"]:
            representative = cluster["representative"]
            cluster_rows.append(
                [
                    cluster["cluster_id"],
                    cluster["population"],
                    f"{cluster['fraction'] * 100:.1f}%",
                    _format_float(cluster.get("best_affinity_kcal_mol")),
                    _format_float(cluster.get("mean_affinity_kcal_mol")),
                    _format_float(cluster.get("sd_affinity_kcal_mol")),
                    _format_float(cluster.get("best_ligand_efficiency")),
                    f"{representative.get('seed', '?')}/"
                    f"{representative.get('model', '?')}",
                    _format_float(cluster.get("max_intra_rmsd_angstrom")),
                ]
            )
        if cluster_rows:
            lines.extend(
                [
                    _table(
                        [
                            "Klaster",
                            "Populasi",
                            "Fraksi",
                            "Afinitas terbaik",
                            "Rata-rata",
                            "SD",
                            "LE terbaik",
                            "Representatif",
                            "Max RMSD intra",
                        ],
                        cluster_rows,
                    ),
                    "",
                ]
            )

        # Tabel anggota untuk klaster terbesar/terbaik (maksimal 12 baris).
        largest = max(
            item["clusters"], key=lambda cluster: cluster["population"], default=None
        )
        if largest and largest["population"] > 1:
            lines.extend(
                [
                    f"### Anggota klaster {largest['cluster_id']} "
                    "(populasi terbesar)",
                    "",
                    _table(
                        ["Pose", "Seed", "Model", "Afinitas", "LE", "RMSD ke rep."],
                        [
                            [
                                member["label"],
                                member["seed"],
                                member["model"],
                                _format_float(member["affinity_kcal_mol"]),
                                _format_float(member["ligand_efficiency"]),
                                _format_float(
                                    member["rmsd_to_representative_angstrom"]
                                ),
                            ]
                            for member in largest["members"][:12]
                        ],
                    ),
                    "",
                ]
            )
            if largest["population"] > 12:
                lines.extend(
                    [
                        f"*(Menampilkan 12 dari {largest['population']} anggota; "
                        "detail lengkap ada di JSON.)*",
                        "",
                    ]
                )

        if item.get("warnings"):
            lines.extend(["### Peringatan", ""])
            lines.extend([f"- {warning}" for warning in item["warnings"]])
            lines.append("")

    lines.extend(["## Keterbatasan", ""])
    lines.extend([f"- {caveat}" for caveat in payload.get("caveats", [])])
    lines.append("")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Klasterkan pose Vina berdasarkan RMSD ligan di frame reseptor "
            "(default cutoff 2.0 Å)."
        )
    )
    parser.add_argument("--docking", type=Path, default=DEFAULT_DOCKING)
    parser.add_argument(
        "--receptor",
        nargs="+",
        default=None,
        help="Nama reseptor (default: semua; 'all' untuk semua).",
    )
    parser.add_argument(
        "--ligand",
        nargs="+",
        default=["DHN_monomer"],
        help="Nama ligan (default: DHN_monomer; 'all' untuk semua).",
    )
    parser.add_argument(
        "--rmsd-cutoff",
        type=float,
        default=2.0,
        help="Ambang klaster RMSD dalam Å (default: 2.0).",
    )
    parser.add_argument(
        "--linkage",
        choices=("complete", "average"),
        default="complete",
        help=(
            "Complete linkage (default) menjamin semua pasangan anggota "
            "<= cutoff; average hanya rata-rata."
        ),
    )
    parser.add_argument(
        "--rmsd-mode",
        choices=("receptor", "aligned"),
        default="receptor",
        help=(
            "receptor: RMSD tanpa superposisi (default, sesuai pose docking); "
            "aligned: RMSD setelah superposisi (ligan rigid bisa ~0)."
        ),
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
            rows = _select_rows(payload, args.receptor, args.ligand)
            fallback_seed = (payload.get("parameters") or {}).get("seed")
            for row in rows:
                sources = _collect_pose_files(row)
                records, _ = build_pose_records(row, fallback_seed=fallback_seed)
                print(
                    f"{row.get('receptor')} x {row.get('ligand')}: "
                    f"{len(records)} pose dari {len(sources)} file; "
                    f"afinitas terbaik "
                    f"{min((r.affinity for r in records if r.affinity is not None), default=None)}"
                )
            return 0

        payload = build_payload(
            args.docking,
            receptors=args.receptor,
            ligands=args.ligand,
            cutoff=args.rmsd_cutoff,
            method=args.linkage,
            rmsd_mode=args.rmsd_mode,
        )
        _write_json_atomic(args.json, payload)
        write_tsv(args.tsv, payload)
        write_markdown(args.md, payload)
        print(
            f"[SUCCESS] {payload['summary']['complexes_total']} kompleks; "
            f"{payload['summary']['poses_total']} pose; "
            f"{payload['summary']['clusters_total']} klaster."
        )
        print(f"    JSON : {args.json}")
        print(f"    TSV  : {args.tsv}")
        print(f"    MD   : {args.md}")
        for error in payload.get("errors", []):
            print(
                f"    [WARN] {error['receptor']} x {error['ligand']}: "
                f"{error['error']}"
            )
        return 0 if payload["status"] != "failed" else 1
    except (ClusterError, OSError, ValueError, KeyError, TypeError) as exc:
        print(f"[ERROR] {exc}")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
