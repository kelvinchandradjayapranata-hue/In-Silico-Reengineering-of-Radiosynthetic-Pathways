"""Analisis trajektori MD OpenMM protein-ligan dengan MDTraj.

Modul ini membaca trajektori produksi OpenMM (``production.dcd``) dan topology
``system.pdb`` dari direktori ``results/md_prep/<sistem>/``, lalu menghitung:

- RMSD backbone protein setelah superposisi ke frame pertama,
- RMSD ligan di frame reseptor (setelah superposisi backbone),
- RMSF per residu berdasarkan atom CA,
- occupany ikatan hidrogen protein-ligan sepanjang trajektori.

Hasil utama adalah grafik empat panel di
``results/plots/md_trajectory_analysis.png``. Ringkasan numerik ditulis ke
JSON dan Markdown agar dapat diaudit.

Analisis ini bersifat geometris. Occupancy H-bond bukan energi ikat, dan RMSD
ligan tidak dikoreksi simetri graph. Hasil MD tetap memerlukan validasi
eksperimental.
"""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Sequence

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MD_PREP = PROJECT_ROOT / "results" / "md_prep"
DEFAULT_SYSTEM_DIR = (
    DEFAULT_MD_PREP / "PKS1_surrogate_8CG4__DHN_monomer"
)
DEFAULT_OUT = PROJECT_ROOT / "results" / "plots" / "md_trajectory_analysis.png"
DEFAULT_JSON = PROJECT_ROOT / "results" / "md_trajectory_analysis.json"
DEFAULT_MD = PROJECT_ROOT / "results" / "md_trajectory_analysis.md"

LIGAND_RESNAME = "LIG"
POLAR_ELEMENTS = {"N", "O"}
BACKBONE_NAMES = {"N", "CA", "C", "O"}
PROTEIN_RESIDUES = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "HID", "HIE", "HIP", "ASH", "GLH", "LYN", "CYM", "CYX",
}

# Palet kategorikal referensi (slot 1-3) dan warna chrome grafik.
COLOR_BLUE = "#2a78d6"
COLOR_ORANGE = "#eb6834"
COLOR_AQUA = "#1baf7a"
COLOR_INK = "#0b0b0b"
COLOR_SECONDARY = "#52514e"
COLOR_MUTED = "#898781"
COLOR_GRID = "#e1e0d9"
COLOR_BASELINE = "#c3c2b7"
COLOR_SURFACE = "#fcfcfb"


class TrajectoryError(RuntimeError):
    """Kesalahan analisis yang dapat ditampilkan tanpa traceback."""


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _relative(path: Path) -> str:
    try:
        return path.resolve().relative_to(PROJECT_ROOT.resolve()).as_posix()
    except ValueError:
        return str(path)


def _read_json(path: Path) -> dict[str, Any]:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError:
        return {}
    except json.JSONDecodeError as exc:
        raise TrajectoryError(f"JSON tidak valid: {path} ({exc})") from exc


def _element(topology: Any, index: int) -> str:
    element = topology.atom(int(index)).element
    symbol = getattr(element, "symbol", "") if element is not None else ""
    return str(symbol).upper()


def _atom_label(topology: Any, index: int) -> str:
    atom = topology.atom(int(index))
    residue = atom.residue
    chain = residue.chain
    chain_id = getattr(chain, "chain_id", None)
    prefix = f"{residue.name}{residue.resSeq}"
    if chain_id:
        prefix = f"{chain_id}:{prefix}"
    return f"{prefix}:{atom.name}"


def _residue_label(residue: Any) -> str:
    chain_id = getattr(residue.chain, "chain_id", None)
    prefix = f"{residue.name}{residue.resSeq}"
    return f"{chain_id}:{prefix}" if chain_id else prefix


def _list_systems(base: Path = DEFAULT_MD_PREP) -> list[Path]:
    if not base.is_dir():
        return []
    systems = []
    for path in sorted(base.iterdir()):
        if (path / "system.pdb").is_file() and (
            path / "trajectory" / "production.dcd"
        ).is_file():
            systems.append(path)
    return systems


def _resolve_inputs(args: argparse.Namespace) -> tuple[Path, Path, Path, Path]:
    system_dir = (
        Path(args.system_dir).expanduser()
        if args.system_dir
        else DEFAULT_SYSTEM_DIR
    )
    system_dir = system_dir.resolve()
    topology = (
        Path(args.topology).expanduser().resolve()
        if args.topology
        else system_dir / "system.pdb"
    )
    trajectory = (
        Path(args.trajectory).expanduser().resolve()
        if args.trajectory
        else system_dir / "trajectory" / "production.dcd"
    )
    config = (
        Path(args.config).expanduser().resolve()
        if args.config
        else system_dir / "md_config.json"
    )
    if not topology.is_file():
        raise TrajectoryError(f"Topology tidak ditemukan: {topology}")
    if not trajectory.is_file():
        raise TrajectoryError(
            f"Trajektori tidak ditemukan: {trajectory}\n"
            "Jalankan run_md.py di Colab sampai tahap produksi selesai."
        )
    return system_dir, topology, trajectory, config


def _load_trajectory(
    trajectory: Path, topology: Path, stride: int
) -> Any:
    try:
        import mdtraj as md
    except ImportError as exc:
        raise TrajectoryError(
            "MDTraj belum terpasang. Jalankan "
            "`python -m pip install -r requirements.txt`."
        ) from exc

    try:
        traj = md.load(str(trajectory), top=str(topology), stride=stride)
    except Exception as exc:
        raise TrajectoryError(
            f"Gagal membaca trajektori {trajectory.name}: "
            f"{type(exc).__name__}: {exc}"
        ) from exc

    if traj.n_frames < 1:
        raise TrajectoryError("Trajektori tidak memiliki frame.")
    if traj.n_atoms != traj.topology.n_atoms:
        raise TrajectoryError(
            f"Jumlah atom DCD ({traj.n_atoms}) tidak cocok dengan topology "
            f"({traj.topology.n_atoms})."
        )
    return traj


def _frame_times(
    n_frames: int,
    config: dict[str, Any],
    stride: int,
) -> tuple[np.ndarray, float | None, str]:
    """Kembalikan waktu relatif (ns), interval per frame, dan satuan label."""
    md_config = config.get("md", {}) if isinstance(config, dict) else {}
    report_ps = md_config.get("report_interval_ps")
    if isinstance(report_ps, (int, float)) and report_ps > 0:
        interval = float(report_ps) * max(1, int(stride))
        return np.arange(n_frames, dtype=float) * interval / 1000.0, interval, "ns"

    production_ns = md_config.get("production_ns")
    if isinstance(production_ns, (int, float)) and production_ns > 0:
        duration = float(production_ns)
        return np.linspace(0.0, duration, n_frames, endpoint=False), None, "ns"

    return np.arange(n_frames, dtype=float), None, "frame"


def _protein_and_ligand_indices(
    topology: Any, ligand_resname: str
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    protein: list[int] = []
    backbone: list[int] = []
    ca: list[int] = []
    ligand: list[int] = []
    for atom in topology.atoms:
        residue = atom.residue
        is_protein = residue.is_protein or residue.name.upper() in PROTEIN_RESIDUES
        if is_protein:
            protein.append(atom.index)
            if atom.name in BACKBONE_NAMES:
                backbone.append(atom.index)
            if atom.name == "CA":
                ca.append(atom.index)
        elif residue.name == ligand_resname:
            ligand.append(atom.index)

    if not protein:
        raise TrajectoryError("Tidak ada atom protein pada topology.")
    if not backbone:
        raise TrajectoryError("Tidak ada atom backbone (N/CA/C/O) pada protein.")
    if not ca:
        raise TrajectoryError("Tidak ada atom CA pada protein.")
    if not ligand:
        raise TrajectoryError(
            f"Residu ligan '{ligand_resname}' tidak ditemukan. "
            "Gunakan --ligand-resname sesuai system.pdb."
        )

    if not any(_element(topology, index) != "H" for index in ligand):
        raise TrajectoryError("Ligan tidak memiliki atom berat.")

    return (
        np.asarray(protein, dtype=int),
        np.asarray(backbone, dtype=int),
        np.asarray(ca, dtype=int),
        np.asarray(ligand, dtype=int),
    )


def _donor_hydrogens(
    topology: Any,
    atom_indices: np.ndarray,
    coords0: np.ndarray,
) -> list[tuple[int, int]]:
    """Pasangan (donor N/O, H) untuk atom dalam seleksi.

    Ikatan topology dipakai lebih dulu. Bila PDB tidak menulis CONECT dan
    MDTraj tidak membangun bond, H tanpa pasangan kovalen dipetakan ke N/O
    terdekat (<=1.35 Å).
    """
    index_set = set(int(index) for index in atom_indices)
    polar = np.asarray(
        [i for i in atom_indices if _element(topology, int(i)) in POLAR_ELEMENTS],
        dtype=int,
    )
    hydrogens = np.asarray(
        [i for i in atom_indices if _element(topology, int(i)) == "H"],
        dtype=int,
    )
    if len(polar) == 0 or len(hydrogens) == 0:
        return []

    adjacency: dict[int, list[int]] = {}
    for bond in topology.bonds:
        first = bond.atom1.index
        second = bond.atom2.index
        adjacency.setdefault(first, []).append(second)
        adjacency.setdefault(second, []).append(first)

    pairs: list[tuple[int, int]] = []
    for hydrogen in hydrogens:
        neighbors = adjacency.get(int(hydrogen), [])
        polar_neighbors = [
            neighbor
            for neighbor in neighbors
            if neighbor in index_set
            and _element(topology, neighbor) in POLAR_ELEMENTS
        ]
        if polar_neighbors:
            for donor in polar_neighbors:
                pairs.append((int(donor), int(hydrogen)))
            continue

        # Jika H sudah terikat pada atom berat non-polar (mis. C-H), jangan
        # dipetakan ke N/O terdekat; itu akan menghasilkan donor palsu.
        if any(
            _element(topology, neighbor) != "H" for neighbor in neighbors
        ):
            continue

        distances = np.linalg.norm(
            coords0[polar] - coords0[int(hydrogen)], axis=1
        )
        nearest = int(np.argmin(distances))
        if distances[nearest] <= 0.135:  # nm; 1.35 Å
            pairs.append((int(polar[nearest]), int(hydrogen)))
    return pairs


def _hbond_occupancy(
    traj: Any,
    topology: Any,
    protein_indices: np.ndarray,
    ligand_indices: np.ndarray,
    *,
    cutoff_angstrom: float,
    min_angle_deg: float,
    min_occupancy: float,
) -> tuple[list[dict[str, Any]], dict[str, int]]:
    """Hitung occupancy H-bond protein-ligan.

    Hanya ikatan antarmolekul yang dihitung (donor protein -> akseptor ligan
    dan donor ligan -> akseptor protein). Air dan ion tidak disertakan.
    """
    import mdtraj as md

    coords0 = traj.xyz[0]
    protein_donors = _donor_hydrogens(topology, protein_indices, coords0)
    ligand_donors = _donor_hydrogens(topology, ligand_indices, coords0)

    protein_acceptors = [
        int(index)
        for index in protein_indices
        if _element(topology, int(index)) in POLAR_ELEMENTS
    ]
    ligand_acceptors = [
        int(index)
        for index in ligand_indices
        if _element(topology, int(index)) in POLAR_ELEMENTS
    ]

    triplets: list[tuple[int, int, int]] = []
    kinds: list[str] = []
    for donor, hydrogen in protein_donors:
        for acceptor in ligand_acceptors:
            triplets.append((donor, hydrogen, acceptor))
            kinds.append("protein->ligand")
    for donor, hydrogen in ligand_donors:
        for acceptor in protein_acceptors:
            triplets.append((donor, hydrogen, acceptor))
            kinds.append("ligand->protein")

    stats = {
        "protein_donors": len(protein_donors),
        "ligand_donors": len(ligand_donors),
        "protein_acceptors": len(protein_acceptors),
        "ligand_acceptors": len(ligand_acceptors),
        "candidate_triplets": len(triplets),
    }
    if not triplets:
        return [], stats

    angle_indices = np.asarray(triplets, dtype=int)
    distance_indices = angle_indices[:, [0, 2]]
    periodic = traj.unitcell_lengths is not None
    distances = md.compute_distances(
        traj, distance_indices, periodic=periodic
    )
    angles = md.compute_angles(traj, angle_indices, periodic=periodic)

    cutoff_nm = cutoff_angstrom / 10.0
    min_angle = math.radians(min_angle_deg)
    present = (distances <= cutoff_nm) & (angles >= min_angle)
    occupancy = present.mean(axis=0)

    bonds: list[dict[str, Any]] = []
    for position, (donor, hydrogen, acceptor) in enumerate(triplets):
        fraction = float(occupancy[position])
        if fraction <= 0 or fraction < min_occupancy:
            continue
        present_frames = present[:, position]
        mean_distance = float(distances[present_frames, position].mean())
        mean_angle = float(angles[present_frames, position].mean())
        bonds.append(
            {
                "donor": _atom_label(topology, donor),
                "acceptor": _atom_label(topology, acceptor),
                "donor_residue": _residue_label(
                    topology.atom(donor).residue
                ),
                "acceptor_residue": _residue_label(
                    topology.atom(acceptor).residue
                ),
                "direction": kinds[position],
                "occupancy": round(fraction, 4),
                "mean_distance_angstrom": round(mean_distance * 10.0, 3),
                "mean_angle_deg": round(math.degrees(mean_angle), 2),
                "frames_present": int(present_frames.sum()),
                "frames_total": int(len(traj)),
            }
        )

    bonds.sort(key=lambda item: (-item["occupancy"], item["donor"]))
    return bonds, stats


def _window_means(values: np.ndarray) -> tuple[float, float]:
    if len(values) == 0:
        return 0.0, 0.0
    window = max(1, int(round(len(values) * 0.1)))
    return float(np.mean(values[:window])), float(np.mean(values[-window:]))


def _describe(values: np.ndarray, *, temporal: bool = True) -> dict[str, Any]:
    summary: dict[str, Any] = {
        "mean_angstrom": round(float(np.mean(values)), 3),
        "sd_angstrom": round(float(np.std(values)), 3),
        "min_angstrom": round(float(np.min(values)), 3),
        "max_angstrom": round(float(np.max(values)), 3),
    }
    if temporal:
        first, last = _window_means(values)
        summary.update(
            {
                "first_10pct_mean_angstrom": round(first, 3),
                "last_10pct_mean_angstrom": round(last, 3),
                "drift_angstrom": round(last - first, 3),
            }
        )
    return summary


def _style_axis(axis: Any) -> None:
    axis.set_facecolor(COLOR_SURFACE)
    axis.grid(axis="y", color=COLOR_GRID, linewidth=0.8)
    axis.set_axisbelow(True)
    for spine in ("top", "right"):
        axis.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        axis.spines[spine].set_color(COLOR_BASELINE)
        axis.spines[spine].set_linewidth(0.8)
    axis.tick_params(colors=COLOR_MUTED, labelsize=9)
    axis.xaxis.label.set_color(COLOR_SECONDARY)
    axis.yaxis.label.set_color(COLOR_SECONDARY)
    axis.title.set_color(COLOR_INK)


def _plot_analysis(
    *,
    system_name: str,
    times: np.ndarray,
    time_unit: str,
    backbone_rmsd: np.ndarray,
    ligand_rmsd: np.ndarray,
    rmsf: np.ndarray,
    rmsf_residues: list[Any],
    hbonds: list[dict[str, Any]],
    hbond_stats: dict[str, int],
    hbond_threshold: float,
    out_path: Path,
    dpi: int,
    top_hbonds: int,
) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise TrajectoryError(
            "Matplotlib belum terpasang. Jalankan "
            "`python -m pip install -r requirements.txt`."
        ) from exc

    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Segoe UI", "DejaVu Sans"],
            "axes.titlesize": 11,
            "axes.labelsize": 10,
            "figure.facecolor": COLOR_SURFACE,
            "savefig.facecolor": COLOR_SURFACE,
        }
    )

    fig, axes = plt.subplots(
        2,
        2,
        figsize=(12.5, 9.0),
        facecolor=COLOR_SURFACE,
        constrained_layout=False,
    )
    fig.subplots_adjust(
        top=0.86,
        bottom=0.08,
        left=0.075,
        right=0.98,
        hspace=0.42,
        wspace=0.28,
    )
    axis_bb, axis_lig, axis_rmsf, axis_hb = axes.ravel()
    x_label = "Waktu (ns)" if time_unit == "ns" else "Frame"

    # Panel 1: RMSD backbone.
    axis_bb.plot(times, backbone_rmsd, color=COLOR_BLUE, linewidth=1.6)
    mean_bb = float(np.mean(backbone_rmsd))
    axis_bb.axhline(
        mean_bb, color=COLOR_MUTED, linewidth=0.9, linestyle=(0, (4, 3))
    )
    axis_bb.set_title("RMSD backbone protein")
    axis_bb.set_xlabel(x_label)
    axis_bb.set_ylabel("RMSD (Å)")
    axis_bb.annotate(
        f"rata-rata {mean_bb:.2f} Å",
        xy=(0.98, 0.92),
        xycoords="axes fraction",
        ha="right",
        va="top",
        fontsize=9,
        color=COLOR_SECONDARY,
    )

    # Panel 2: RMSD ligan di frame reseptor.
    axis_lig.plot(times, ligand_rmsd, color=COLOR_ORANGE, linewidth=1.6)
    mean_lig = float(np.mean(ligand_rmsd))
    axis_lig.axhline(
        mean_lig, color=COLOR_MUTED, linewidth=0.9, linestyle=(0, (4, 3))
    )
    axis_lig.set_title("RMSD ligan DHN (frame reseptor)")
    axis_lig.set_xlabel(x_label)
    axis_lig.set_ylabel("RMSD (Å)")
    axis_lig.annotate(
        f"rata-rata {mean_lig:.2f} Å",
        xy=(0.98, 0.92),
        xycoords="axes fraction",
        ha="right",
        va="top",
        fontsize=9,
        color=COLOR_SECONDARY,
    )

    # Panel 3: RMSF per residu (CA).
    residue_x = np.arange(len(rmsf), dtype=float)
    axis_rmsf.plot(residue_x, rmsf, color=COLOR_AQUA, linewidth=1.4)
    peak = int(np.argmax(rmsf))
    axis_rmsf.scatter(
        [residue_x[peak]],
        [rmsf[peak]],
        s=28,
        color=COLOR_AQUA,
        edgecolor=COLOR_SURFACE,
        linewidth=1.0,
        zorder=3,
    )
    axis_rmsf.set_title(
        f"RMSF per residu (CA) · maks {rmsf[peak]:.2f} Å"
    )
    axis_rmsf.set_xlabel("Residu")
    axis_rmsf.set_ylabel("RMSF (Å)")
    tick_count = min(8, len(rmsf))
    if tick_count > 0:
        tick_positions = np.unique(
            np.linspace(0, len(rmsf) - 1, tick_count).astype(int)
        )
        axis_rmsf.set_xticks(tick_positions)
        axis_rmsf.set_xticklabels(
            [_residue_label(rmsf_residues[i]) for i in tick_positions],
            rotation=40,
            ha="right",
            fontsize=8,
        )

    # Panel 4: occupancy H-bond protein-ligan.
    if hbonds:
        shown = hbonds[: max(1, top_hbonds)]
        labels = [
            f"{item['donor']} → {item['acceptor']}" for item in shown
        ]
        values = [item["occupancy"] * 100.0 for item in shown]
        order = np.argsort(values)
        axis_hb.barh(
            np.arange(len(shown)),
            [values[i] for i in order],
            color=COLOR_BLUE,
            height=0.62,
        )
        axis_hb.set_yticks(np.arange(len(shown)))
        axis_hb.set_yticklabels([labels[i] for i in order], fontsize=8)
        for position, index in enumerate(order):
            axis_hb.text(
                values[index] + 1.0,
                position,
                f"{values[index]:.0f}%",
                va="center",
                fontsize=8,
                color=COLOR_SECONDARY,
            )
        axis_hb.set_xlim(0, max(105.0, max(values) * 1.18))
        axis_hb.set_title(
            f"H-bond protein–ligan (≥{hbond_threshold * 100:.0f}% frame)"
        )
        axis_hb.set_xlabel("Occupancy (% frame)")
    else:
        if hbond_stats.get("skipped"):
            message = "Analisis H-bond dilewati\n(--no-hbonds)"
        else:
            message = (
                "Tidak ada H-bond protein–ligan\n"
                "yang memenuhi ambang\n"
                f"({hbond_stats.get('candidate_triplets', 0)} kandidat)"
            )
        axis_hb.text(
            0.5,
            0.5,
            message,
            ha="center",
            va="center",
            transform=axis_hb.transAxes,
            fontsize=10,
            color=COLOR_SECONDARY,
        )
        axis_hb.set_title("H-bond protein–ligan")
        axis_hb.set_xticks([])
        axis_hb.set_yticks([])

    for axis in axes.ravel():
        _style_axis(axis)

    duration = times[-1] if len(times) else 0.0
    fig.suptitle(
        f"Analisis Trajektori MD — {system_name}",
        fontsize=14,
        color=COLOR_INK,
        x=0.01,
        y=0.995,
        ha="left",
        va="top",
    )
    subtitle = (
        f"{len(times):,} frame dianalisis"
        if time_unit == "frame"
        else f"{len(times):,} frame · 0–{duration:.2f} ns"
    )
    fig.text(
        0.01,
        0.950,
        subtitle,
        fontsize=9,
        color=COLOR_SECONDARY,
        ha="left",
        va="top",
    )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def _write_markdown(
    path: Path,
    *,
    system_name: str,
    inputs: dict[str, str],
    trajectory_info: dict[str, Any],
    backbone: dict[str, Any],
    ligand: dict[str, Any],
    rmsf_summary: dict[str, Any],
    hbonds: list[dict[str, Any]],
    hbond_stats: dict[str, int],
    warnings: list[str],
) -> None:
    lines = [
        f"# Analisis Trajektori MD — {system_name}",
        "",
        f"> Dibuat: {_now()}",
        f"> Trajektori: `{inputs['trajectory']}`",
        f"> Topology: `{inputs['topology']}`",
        f"> Frame dianalisis: {trajectory_info['n_frames_analyzed']:,} "
        f"(stride {trajectory_info['stride']})",
        "",
        "## Ringkasan",
        "",
        "| Metrik | Rata-rata | SD | Minimum | Maksimum |",
        "|---|---|---|---|---|",
        (
            f"| RMSD backbone (Å) | {backbone['mean_angstrom']:.3f} | "
            f"{backbone['sd_angstrom']:.3f} | {backbone['min_angstrom']:.3f} | "
            f"{backbone['max_angstrom']:.3f} |"
        ),
        (
            f"| RMSD ligan (Å) | {ligand['mean_angstrom']:.3f} | "
            f"{ligand['sd_angstrom']:.3f} | {ligand['min_angstrom']:.3f} | "
            f"{ligand['max_angstrom']:.3f} |"
        ),
        (
            f"| RMSF CA (Å) | {rmsf_summary['mean_angstrom']:.3f} | "
            f"{rmsf_summary['sd_angstrom']:.3f} | "
            f"{rmsf_summary['min_angstrom']:.3f} | "
            f"{rmsf_summary['max_angstrom']:.3f} |"
        ),
        "",
        "Perubahan RMSD dari 10% frame awal ke 10% frame akhir:",
        "",
        f"- Backbone: {backbone['drift_angstrom']:+.3f} Å",
        f"- Ligan: {ligand['drift_angstrom']:+.3f} Å",
        "",
        "## H-bond protein–ligan",
        "",
        f"Kandidat triplet: {hbond_stats.get('candidate_triplets', 0)}; "
        f"H-bond dilaporkan: {len(hbonds)}.",
        "",
    ]
    if hbonds:
        lines.extend(
            [
                "| Donor | Akseptor | Arah | Occupancy | Jarak (Å) | Sudut (°) |",
                "|---|---|---|---|---|---|",
            ]
        )
        for item in hbonds[:20]:
            direction = (
                "protein→ligan"
                if item["direction"] == "protein->ligand"
                else "ligan→protein"
            )
            lines.append(
                f"| {item['donor']} | {item['acceptor']} | {direction} | "
                f"{item['occupancy'] * 100:.1f}% | "
                f"{item['mean_distance_angstrom']:.2f} | "
                f"{item['mean_angle_deg']:.1f} |"
            )
        lines.append("")
    else:
        lines.extend(
            [
                "Tidak ada H-bond protein–ligan yang memenuhi ambang.",
                "",
            ]
        )

    lines.extend(
        [
            "## Keterbatasan",
            "",
            "- RMSD ligan dihitung di frame reseptor setelah superposisi "
            "backbone, tanpa koreksi simetri graph.",
            "- Occupancy adalah fraksi frame yang memenuhi kriteria geometris "
            "(jarak dan sudut), bukan energi ikat atau laju asosiasi.",
            "- H-bond yang dihitung hanya antarmolekul protein–ligan; "
            "H-bond intramolekul protein tidak dilaporkan.",
            "- Waktu frame dihitung dari `report_interval_ps` di "
            "`md_config.json`; DCD tidak menyimpan timestamp.",
            "- Hasil MD bergantung pada force field, protonasi, dan sampling. "
            "Interpretasi biologis memerlukan validasi eksperimental.",
            "",
        ]
    )
    if warnings:
        lines.extend(["## Peringatan", ""])
        lines.extend(f"- {warning}" for warning in warnings)
        lines.append("")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Analisis RMSD, RMSF, dan occupancy H-bond dari trajektori MD "
            "OpenMM."
        )
    )
    parser.add_argument(
        "--system-dir",
        type=Path,
        default=None,
        help=(
            "Direktori sistem MD (default: results/md_prep/"
            "PKS1_surrogate_8CG4__DHN_monomer)."
        ),
    )
    parser.add_argument("--trajectory", type=Path, default=None)
    parser.add_argument("--topology", type=Path, default=None)
    parser.add_argument("--config", type=Path, default=None)
    parser.add_argument("--ligand-resname", default=LIGAND_RESNAME)
    parser.add_argument("--stride", type=int, default=1)
    parser.add_argument("--hbond-cutoff", type=float, default=3.5)
    parser.add_argument("--hbond-angle", type=float, default=120.0)
    parser.add_argument("--hbond-occupancy", type=float, default=0.05)
    parser.add_argument("--top-hbonds", type=int, default=10)
    parser.add_argument("--no-hbonds", action="store_true")
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--json", type=Path, default=DEFAULT_JSON)
    parser.add_argument("--md", type=Path, default=DEFAULT_MD)
    parser.add_argument("--dpi", type=int, default=200)
    parser.add_argument("--list", action="store_true")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.list:
        systems = _list_systems()
        if not systems:
            print("[INFO] Belum ada sistem MD dengan production.dcd.")
            return 0
        print("Sistem MD yang tersedia:")
        for path in systems:
            print(f"- {_relative(path)}")
        return 0

    if args.stride < 1:
        print("[ERROR] --stride harus >= 1.")
        return 2
    if args.hbond_cutoff <= 0 or not 0 < args.hbond_angle < 180:
        print("[ERROR] Parameter H-bond tidak valid.")
        return 2
    if not 0 <= args.hbond_occupancy <= 1:
        print("[ERROR] --hbond-occupancy harus antara 0 dan 1.")
        return 2

    try:
        system_dir, topology_path, trajectory_path, config_path = (
            _resolve_inputs(args)
        )
        config = _read_json(config_path) if config_path.is_file() else {}
        traj = _load_trajectory(
            trajectory_path, topology_path, max(1, args.stride)
        )
        topology = traj.topology
        protein_idx, backbone_idx, ca_idx, ligand_all_idx = (
            _protein_and_ligand_indices(topology, args.ligand_resname)
        )
        ligand_idx = np.asarray(
            [
                index
                for index in ligand_all_idx
                if _element(topology, int(index)) != "H"
            ],
            dtype=int,
        )

        times, frame_interval_ps, time_unit = _frame_times(
            traj.n_frames, config, args.stride
        )

        # RMSD backbone: superposisi setiap frame ke frame pertama.
        import mdtraj as md

        reference = traj[0]
        backbone_rmsd_nm = np.asarray(
            md.rmsd(
                traj,
                reference,
                atom_indices=backbone_idx,
                parallel=False,
            )
        )
        aligned = traj.superpose(
            reference, atom_indices=backbone_idx, parallel=False
        )

        # RMSD ligan di frame reseptor: tanpa superposisi ulang pada ligan.
        ligand_delta = (
            aligned.xyz[:, ligand_idx, :] - reference.xyz[0, ligand_idx, :]
        )
        ligand_rmsd_nm = np.sqrt(
            np.mean(np.sum(ligand_delta * ligand_delta, axis=2), axis=1)
        )

        # RMSF per residu memakai CA setelah superposisi backbone. Dihitung
        # manual terhadap rata-rata koordinat CA agar tidak bergantung pada
        # API RMSF yang berbeda antar versi MDTraj.
        ca_xyz = aligned.xyz[:, ca_idx, :]
        ca_mean = ca_xyz.mean(axis=0)
        rmsf_nm = np.sqrt(
            np.mean(np.sum((ca_xyz - ca_mean) ** 2, axis=2), axis=0)
        )
        rmsf_residues = [topology.atom(index).residue for index in ca_idx]
        rmsf_angstrom = rmsf_nm * 10.0

        backbone_angstrom = backbone_rmsd_nm * 10.0
        ligand_angstrom = ligand_rmsd_nm * 10.0

        warnings: list[str] = []
        if args.no_hbonds:
            hbonds: list[dict[str, Any]] = []
            hbond_stats = {"candidate_triplets": 0, "skipped": 1}
            warnings.append("Analisis H-bond dilewati (--no-hbonds).")
        else:
            try:
                hbonds, hbond_stats = _hbond_occupancy(
                    traj,
                    topology,
                    protein_idx,
                    ligand_all_idx,
                    cutoff_angstrom=args.hbond_cutoff,
                    min_angle_deg=args.hbond_angle,
                    min_occupancy=args.hbond_occupancy,
                )
            except TrajectoryError as exc:
                warnings.append(str(exc))
                hbonds, hbond_stats = [], {"candidate_triplets": 0}
            except Exception as exc:
                warnings.append(
                    "Perhitungan H-bond gagal: "
                    f"{type(exc).__name__}: {exc}"
                )
                hbonds, hbond_stats = [], {"candidate_triplets": 0}
            if hbond_stats.get("candidate_triplets", 0) == 0:
                warnings.append(
                    "Tidak ada kandidat H-bond; periksa atom H pada topology."
                )

        backbone_summary = _describe(backbone_angstrom)
        ligand_summary = _describe(ligand_angstrom)
        rmsf_summary = _describe(rmsf_angstrom, temporal=False)

        out_path = Path(args.out).expanduser().resolve()
        _plot_analysis(
            system_name=system_dir.name,
            times=times,
            time_unit=time_unit,
            backbone_rmsd=backbone_angstrom,
            ligand_rmsd=ligand_angstrom,
            rmsf=rmsf_angstrom,
            rmsf_residues=rmsf_residues,
            hbonds=hbonds,
            hbond_stats=hbond_stats,
            hbond_threshold=args.hbond_occupancy,
            out_path=out_path,
            dpi=args.dpi,
            top_hbonds=args.top_hbonds,
        )

        payload = {
            "tool": "melanin_biotech_md_trajectory_analysis",
            "generated_at": _now(),
            "system": system_dir.name,
            "inputs": {
                "trajectory": _relative(trajectory_path),
                "topology": _relative(topology_path),
                "config": _relative(config_path)
                if config_path.is_file()
                else None,
            },
            "parameters": {
                "stride": args.stride,
                "ligand_resname": args.ligand_resname,
                "hbond_cutoff_angstrom": args.hbond_cutoff,
                "hbond_angle_deg": args.hbond_angle,
                "hbond_occupancy_threshold": args.hbond_occupancy,
                "hbond_analysis": not args.no_hbonds,
            },
            "trajectory": {
                "stride": args.stride,
                "n_frames_analyzed": int(traj.n_frames),
                "n_atoms": int(traj.n_atoms),
                "protein_atoms": int(len(protein_idx)),
                "backbone_atoms": int(len(backbone_idx)),
                "ca_atoms": int(len(ca_idx)),
                "ligand_heavy_atoms": int(len(ligand_idx)),
                "frame_interval_ps": round(frame_interval_ps, 3)
                if frame_interval_ps
                else None,
                "analyzed_duration_ns": round(float(times[-1]), 4)
                if time_unit == "ns"
                else None,
                "time_basis": time_unit,
            },
            "rmsd": {
                "backbone": {
                    **backbone_summary,
                    "values_angstrom": [
                        round(float(value), 4) for value in backbone_angstrom
                    ],
                },
                "ligand": {
                    **ligand_summary,
                    "values_angstrom": [
                        round(float(value), 4) for value in ligand_angstrom
                    ],
                },
            },
            "rmsf": {
                **rmsf_summary,
                "residues": [
                    {
                        "residue": _residue_label(residue),
                        "resname": residue.name,
                        "resseq": int(residue.resSeq),
                        "chain": getattr(residue.chain, "chain_id", None),
                        "rmsf_angstrom": round(float(value), 4),
                    }
                    for residue, value in zip(rmsf_residues, rmsf_angstrom)
                ],
            },
            "hbonds": {
                **hbond_stats,
                "reported": len(hbonds),
                "bonds": hbonds,
            },
            "warnings": warnings,
        }

        json_path = Path(args.json).expanduser().resolve()
        json_path.parent.mkdir(parents=True, exist_ok=True)
        json_path.write_text(
            json.dumps(payload, indent=2, ensure_ascii=False),
            encoding="utf-8",
        )

        md_path = Path(args.md).expanduser().resolve()
        _write_markdown(
            md_path,
            system_name=system_dir.name,
            inputs=payload["inputs"],
            trajectory_info=payload["trajectory"],
            backbone=backbone_summary,
            ligand=ligand_summary,
            rmsf_summary=rmsf_summary,
            hbonds=hbonds,
            hbond_stats=hbond_stats,
            warnings=warnings,
        )

        print(f"[SUCCESS] Sistem: {system_dir.name}")
        print(
            f"    Frame: {traj.n_frames:,}; RMSD backbone "
            f"{backbone_summary['mean_angstrom']:.2f} ± "
            f"{backbone_summary['sd_angstrom']:.2f} Å; RMSD ligan "
            f"{ligand_summary['mean_angstrom']:.2f} ± "
            f"{ligand_summary['sd_angstrom']:.2f} Å."
        )
        print(
            f"    RMSF maks: {rmsf_summary['max_angstrom']:.2f} Å; "
            f"H-bond dilaporkan: {len(hbonds)}."
        )
        print(f"    PNG  : {out_path}")
        print(f"    JSON : {json_path}")
        print(f"    MD   : {md_path}")
        return 0
    except TrajectoryError as exc:
        print(f"[ERROR] {exc}")
        return 1
    except Exception as exc:
        print(f"[ERROR] {type(exc).__name__}: {exc}")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
