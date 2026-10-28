"""Agregasi hasil docking multi-seed dari ``src/melanin_docking.py``.

Script ini membaca beberapa file JSON hasil docking (satu file per seed),
menghitung statistik antar-seed, lalu menulis ulang dua artefak kanonik:

- ``results/docking_results.json`` — skema yang sama dengan output modul,
  dengan matriks rata-rata dan nilai terbaik antar-seed.
- ``results/docking_summary.tsv`` — ringkasan panjang per kombinasi
  reseptor x ligan.

Contoh pemakaian setelah menjalankan beberapa seed di Colab:

    python scripts/aggregate_docking_seeds.py \
        --input-dir results/colab_runs/20261028T120000Z \
        --out results/docking_results.json \
        --summary results/docking_summary.tsv
"""

from __future__ import annotations

import argparse
import json
import statistics
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Sequence

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_JSON = PROJECT_ROOT / "results" / "docking_results.json"
DEFAULT_TSV = PROJECT_ROOT / "results" / "docking_summary.tsv"


def _read_json(path: Path) -> dict[str, Any]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise SystemExit(f"[ERROR] Gagal membaca {path}: {exc}") from exc
    if not isinstance(payload, dict):
        raise SystemExit(f"[ERROR] Isi {path} bukan objek JSON.")
    return payload


def _write_json_atomic(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp_path = path.with_suffix(path.suffix + ".tmp")
    temp_path.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    temp_path.replace(path)


def _round(value: float | None, digits: int = 3) -> float | None:
    return round(float(value), digits) if value is not None else None


def _mean(values: Sequence[float]) -> float | None:
    return _round(statistics.mean(values)) if values else None


def _sd(values: Sequence[float]) -> float | None:
    if len(values) < 2:
        return 0.0 if values else None
    return _round(statistics.stdev(values))


def _unique(items: Iterable[dict[str, Any]], key: str) -> list[dict[str, Any]]:
    seen: set[str] = set()
    result: list[dict[str, Any]] = []
    for item in items:
        name = str(item.get(key, "")).strip()
        if not name or name in seen:
            continue
        seen.add(name)
        result.append(item)
    return result


def _collect_payloads(input_dir: Path, pattern: str) -> list[tuple[Path, dict[str, Any]]]:
    paths = sorted(input_dir.glob(pattern))
    if not paths:
        raise SystemExit(
            f"[ERROR] Tidak ada file '{pattern}' di {input_dir}. "
            "Jalankan minimal satu seed terlebih dahulu."
        )
    payloads = [(path, _read_json(path)) for path in paths]
    payloads.sort(
        key=lambda item: (
            item[1].get("parameters", {}).get("seed", 0),
            item[0].name,
        )
    )
    return payloads


def _build_interpretation(
    rows: Sequence[dict[str, Any]], n_seeds: int
) -> str:
    ok_rows = [
        row
        for row in rows
        if row.get("status") == "ok"
        and row.get("best_affinity_kcal_mol") is not None
    ]
    if not ok_rows:
        return (
            "Tidak ada run docking yang berhasil. Periksa log, dependency "
            "Meeko/Vina, dan struktur reseptor."
        )

    best = min(ok_rows, key=lambda row: row["best_affinity_kcal_mol"])
    le_rows = [
        row for row in ok_rows if row.get("best_ligand_efficiency") is not None
    ]
    if not le_rows:
        return (
            f"Dari {n_seeds} seed, afinitas mentah terbaik "
            f"{best['best_affinity_kcal_mol']:.2f} kcal/mol untuk "
            f"{best['ligand']} pada {best['receptor']} (reseptor surrogate PDB)."
        )

    best_le = min(le_rows, key=lambda row: row["best_ligand_efficiency"])
    return (
        f"Agregasi {n_seeds} seed: skor mentah terbaik "
        f"{best['best_affinity_kcal_mol']:.2f} kcal/mol untuk {best['ligand']} "
        f"pada {best['receptor']} (rata-rata {best['mean_affinity_kcal_mol']:.2f} "
        f"± {best['sd_affinity_kcal_mol']:.2f} kcal/mol). Ligand efficiency "
        f"terbaik {best_le['best_ligand_efficiency']:.3f} kcal/mol per heavy "
        f"atom untuk {best_le['ligand']} pada {best_le['receptor']} (rata-rata "
        f"{best_le['mean_ligand_efficiency']:.3f} ± "
        f"{best_le['sd_ligand_efficiency']:.3f}). Karena skor Vina tidak "
        "dinormalisasi terhadap ukuran ligan, bandingkan ligand efficiency, "
        "bukan hanya afinitas mentah. Semua reseptor adalah surrogate homolog; "
        "angka ini bukan afinitas target C. cladosporioides."
    )


def aggregate(
    payloads: Sequence[tuple[Path, dict[str, Any]]],
) -> dict[str, Any]:
    if not payloads:
        raise SystemExit("[ERROR] Tidak ada payload untuk diagregasi.")

    first = payloads[0][1]
    ligs = _unique(
        [ligand for _, payload in payloads for ligand in payload.get("ligands", [])],
        "name",
    )
    recs = _unique(
        [
            receptor
            for _, payload in payloads
            for receptor in payload.get("receptors", [])
        ],
        "name",
    )
    ligand_names = [str(ligand["name"]) for ligand in ligs]
    receptor_names = [str(receptor["name"]) for receptor in recs]
    seeds = [
        payload.get("parameters", {}).get("seed")
        for _, payload in payloads
        if payload.get("parameters", {}).get("seed") is not None
    ]

    runs_by_key: dict[tuple[str, str], list[dict[str, Any]]] = {}
    for path, payload in payloads:
        seed = payload.get("parameters", {}).get("seed")
        for row in payload.get("results", []):
            receptor = str(row.get("receptor", "")).strip()
            ligand = str(row.get("ligand", "")).strip()
            if not receptor or not ligand:
                continue
            enriched = dict(row)
            enriched["_seed"] = seed
            enriched["_source_file"] = path.name
            runs_by_key.setdefault((receptor, ligand), []).append(enriched)

    aggregated: list[dict[str, Any]] = []
    for receptor in receptor_names:
        for ligand in ligand_names:
            runs = runs_by_key.get((receptor, ligand), [])
            ok_runs = [
                run
                for run in runs
                if run.get("status") == "ok"
                and isinstance(run.get("best_affinity_kcal_mol"), (int, float))
            ]
            affinities = [
                float(run["best_affinity_kcal_mol"]) for run in ok_runs
            ]
            efficiencies = [
                float(run["ligand_efficiency"])
                for run in ok_runs
                if isinstance(run.get("ligand_efficiency"), (int, float))
            ]

            per_seed = [
                {
                    "seed": run.get("_seed"),
                    "status": run.get("status", "error"),
                    "best_affinity_kcal_mol": run.get(
                        "best_affinity_kcal_mol"
                    ),
                    "ligand_efficiency": run.get("ligand_efficiency"),
                    "pose_file": run.get("pose_file"),
                    "log": run.get("log"),
                    "error": run.get("error"),
                }
                for run in runs
            ]

            if not ok_runs:
                errors = [
                    str(run.get("error", "")).strip()
                    for run in runs
                    if run.get("error")
                ]
                aggregated.append(
                    {
                        "receptor": receptor,
                        "ligand": ligand,
                        "status": "error",
                        "n_seeds": 0,
                        "seeds": [],
                        "per_seed": per_seed,
                        "error": "; ".join(dict.fromkeys(errors))
                        or "Tidak ada hasil docking yang berhasil.",
                    }
                )
                continue

            best_run = min(
                ok_runs, key=lambda run: float(run["best_affinity_kcal_mol"])
            )
            mean_affinity = _mean(affinities)
            mean_efficiency = _mean(efficiencies)
            aggregated.append(
                {
                    "receptor": receptor,
                    "ligand": ligand,
                    "status": "ok",
                    "n_seeds": len(ok_runs),
                    "seeds": [
                        run.get("_seed") for run in ok_runs
                    ],
                    "best_affinity_kcal_mol": _round(min(affinities)),
                    "mean_affinity_kcal_mol": mean_affinity,
                    "sd_affinity_kcal_mol": _sd(affinities),
                    "best_ligand_efficiency": _round(min(efficiencies))
                    if efficiencies
                    else None,
                    "mean_ligand_efficiency": mean_efficiency,
                    "sd_ligand_efficiency": _sd(efficiencies),
                    # Field ini dipakai generator README/report sebagai
                    # ringkasan satu nilai per kombinasi.
                    "ligand_efficiency": mean_efficiency,
                    "poses": best_run.get("poses") or [],
                    "pose_file": best_run.get("pose_file"),
                    "log": best_run.get("log"),
                    "per_seed": per_seed,
                }
            )

    def matrix(field: str) -> dict[str, dict[str, Any]]:
        result: dict[str, dict[str, Any]] = {}
        for receptor in receptor_names:
            result[receptor] = {}
            for ligand in ligand_names:
                row = next(
                    (
                        item
                        for item in aggregated
                        if item["receptor"] == receptor
                        and item["ligand"] == ligand
                    ),
                    None,
                )
                result[receptor][ligand] = row.get(field) if row else None
        return result

    ok_rows = [row for row in aggregated if row.get("status") == "ok"]
    failed_rows = [row for row in aggregated if row.get("status") != "ok"]
    if not failed_rows:
        status = "completed"
    elif ok_rows:
        status = "partial"
    else:
        status = "failed"

    seed_runs_total = sum(len(runs) for runs in runs_by_key.values())
    seed_runs_ok = sum(
        1
        for runs in runs_by_key.values()
        for run in runs
        if run.get("status") == "ok"
    )
    method = "multi_seed" if len(payloads) > 1 else "single_seed"

    caveats = list(first.get("caveats") or [])
    caveats.extend(
        [
            "Nilai matriks adalah rata-rata antar-seed; field best_* adalah "
            "nilai minimum (paling favorable) antar-seed.",
            "Sebaran antar-seed mencerminkan stochasticitas pencarian Vina, "
            "bukan interval kepercayaan populasi.",
        ]
    )
    if len(payloads) == 1:
        caveats.append(
            "Hanya satu seed yang diagregasi; variabilitas antar-seed belum "
            "diestimasi."
        )
    caveats = list(dict.fromkeys(caveats))

    return {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "tool": "melanin_biotech_docking_aggregate",
        "status": status,
        "engine": {
            **first.get("engine", {}),
            "seeds": seeds,
        },
        "surrogate_notice": first.get("surrogate_notice", ""),
        "parameters": {
            **first.get("parameters", {}),
            "seeds": seeds,
            "aggregation": method,
        },
        "ligands": ligs,
        "receptors": recs,
        "results": aggregated,
        "summary": {
            "runs_total": len(aggregated),
            "runs_ok": len(ok_rows),
            "runs_failed": len(failed_rows),
            "seed_runs_total": seed_runs_total,
            "seed_runs_ok": seed_runs_ok,
            "n_seeds": len(payloads),
            "seeds": seeds,
            "aggregation": {
                "method": method,
                "statistic": "mean",
                "best": "minimum",
                "n_seeds": len(payloads),
            },
            "affinity_matrix": matrix("mean_affinity_kcal_mol"),
            "best_affinity_matrix": matrix("best_affinity_kcal_mol"),
            "sd_affinity_matrix": matrix("sd_affinity_kcal_mol"),
            "ligand_efficiency_matrix": matrix("mean_ligand_efficiency"),
            "best_ligand_efficiency_matrix": matrix(
                "best_ligand_efficiency"
            ),
            "sd_ligand_efficiency_matrix": matrix("sd_ligand_efficiency"),
        },
        "assumptions": list(first.get("assumptions") or []),
        "caveats": caveats,
        "interpretation": _build_interpretation(aggregated, len(payloads)),
        "seed_files": [path.name for path, _ in payloads],
    }


def write_summary_tsv(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    header = (
        "receptor\tligand\tstatus\tn_seeds\tbest_affinity_kcal_mol\t"
        "mean_affinity_kcal_mol\tsd_affinity_kcal_mol\t"
        "best_ligand_efficiency\tmean_ligand_efficiency\t"
        "sd_ligand_efficiency\tseeds\terror\n"
    )
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(header)
        for row in payload.get("results", []):
            values = [
                row.get("receptor", ""),
                row.get("ligand", ""),
                row.get("status", ""),
                row.get("n_seeds", 0),
                row.get("best_affinity_kcal_mol", ""),
                row.get("mean_affinity_kcal_mol", ""),
                row.get("sd_affinity_kcal_mol", ""),
                row.get("best_ligand_efficiency", ""),
                row.get("mean_ligand_efficiency", ""),
                row.get("sd_ligand_efficiency", ""),
                ",".join(str(seed) for seed in row.get("seeds", [])),
                str(row.get("error", "")).replace("\t", " ").replace(
                    "\n", " "
                ),
            ]
            handle.write("\t".join(str(value) for value in values) + "\n")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Agregasi hasil docking multi-seed menjadi artefak kanonik."
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        required=True,
        help="Direktori berisi file JSON per seed.",
    )
    parser.add_argument(
        "--pattern",
        default="docking_seed*.json",
        help="Pola nama file input (default: docking_seed*.json).",
    )
    parser.add_argument("--out", type=Path, default=DEFAULT_JSON)
    parser.add_argument("--summary", type=Path, default=DEFAULT_TSV)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    payloads = _collect_payloads(args.input_dir, args.pattern)
    payload = aggregate(payloads)
    _write_json_atomic(args.out, payload)
    write_summary_tsv(args.summary, payload)
    print(f"[SUCCESS] Agregasi {len(payloads)} seed selesai.")
    print(f"    JSON    : {args.out}")
    print(f"    Summary : {args.summary}")
    print(f"    Status  : {payload['status']}")
    print(f"    {payload['interpretation']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
