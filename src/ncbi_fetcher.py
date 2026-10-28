"""Pengambilan sekuens referensi jalur melanin dari NCBI.

Prinsip penting:
- Setiap aksesi divalidasi terhadap nama organisme dan kata kunci deskripsi
  sebelum ditulis ke disk. Ini mencegah kejadian seperti DQ000286 (Borrelia
  sp. 16S rRNA) yang sebelumnya tersimpan sebagai "cladosporium_sample.fasta".
- NCBI belum memiliki anotasi gen DHN-melanin yang lengkap untuk
  Cladosporium sphaerospermum. Hanya MK213995.1 / QBZ73626.1 yang merupakan
  urutan langsung dari spesies tersebut (PKS tipe III, belum tentu PKS DHN).
  Target lain ditandai sebagai ortholog referensi dari Cladosporium
  cladosporioides dan Cladosporium halotolerans.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Sequence

from Bio import Entrez, SeqIO

PROJECT_ROOT = Path(__file__).resolve().parents[1]
RAW_DIR = PROJECT_ROOT / "data" / "raw"
MANIFEST_PATH = RAW_DIR / "manifest.json"
RAW_README_PATH = RAW_DIR / "README.md"

TOOL_NAME = "melanin_biotech_fetcher"
DEFAULT_EMAIL = "researcher@example.com"

EVIDENCE_LABELS = {
    "direct_sequence": "urutan langsung",
    "ortholog_reference": "ortholog referensi",
}


@dataclass(frozen=True)
class FetchTarget:
    """Satu target sekuens yang akan diunduh dan divalidasi."""

    target_id: str
    gene: str
    pathway: str
    evidence_level: str  # "direct_sequence" atau "ortholog_reference"
    db: str  # "nucleotide" atau "protein"
    accession: str
    expected_organism: str
    expected_keyword: str
    output_name: str
    notes: str = ""


# Registry target yang sudah diverifikasi manual terhadap NCBI (Oktober 2026).
# Jangan menambahkan aksesi tanpa memeriksa organism + deskripsi terlebih dahulu.
TARGETS: tuple[FetchTarget, ...] = (
    FetchTarget(
        target_id="csph_pks_typeiii_gene",
        gene="PKS tipe III (kandidat)",
        pathway="PKS / melanin",
        evidence_level="direct_sequence",
        db="nucleotide",
        accession="MK213995.1",
        expected_organism="Cladosporium sphaerospermum",
        expected_keyword="polyketide synthase",
        output_name="cladosporium_sphaerospermum_typeIII_pks_MK213995.1.fasta",
        notes=(
            "Satu-satunya urutan nukleotida PKS langsung dari C. sphaerospermum "
            "di NCBI. Bersifat parsial dan belum dikonfirmasi sebagai PKS DHN-melanin."
        ),
    ),
    FetchTarget(
        target_id="csph_pks_typeiii_protein",
        gene="PKS tipe III (kandidat)",
        pathway="PKS / melanin",
        evidence_level="direct_sequence",
        db="protein",
        accession="QBZ73626.1",
        expected_organism="Cladosporium sphaerospermum",
        expected_keyword="polyketide synthase",
        output_name="cladosporium_sphaerospermum_typeIII_pks_QBZ73626.1.fasta",
        notes="Terjemahan protein dari MK213995.1; parsial (235 aa).",
    ),
    FetchTarget(
        target_id="ccla_pks1",
        gene="PKS1 (non-reducing PKS)",
        pathway="DHN-melanin",
        evidence_level="ortholog_reference",
        db="protein",
        accession="KAN3049821.1",
        expected_organism="Cladosporium cladosporioides",
        expected_keyword="polyketide synthase",
        output_name="cladosporium_cladosporioides_PKS1_KAN3049821.1.fasta",
        notes=(
            "Ortholog referensi dari C. cladosporioides. Gen ini berdekatan dengan "
            "THN reductase (KAN3049819.1), konsisten dengan klaster DHN-melanin."
        ),
    ),
    FetchTarget(
        target_id="ccla_thn_reductase",
        gene="THN reductase",
        pathway="DHN-melanin",
        evidence_level="ortholog_reference",
        db="protein",
        accession="KAN3049819.1",
        expected_organism="Cladosporium cladosporioides",
        expected_keyword="tetrahydroxynaphthalene reductase",
        output_name="cladosporium_cladosporioides_THN_reductase_KAN3049819.1.fasta",
        notes="Ortholog referensi; langkah reduksi THN -> scytalone.",
    ),
    FetchTarget(
        target_id="ccla_scytalone_dehydratase",
        gene="Scytalone dehydratase",
        pathway="DHN-melanin",
        evidence_level="ortholog_reference",
        db="protein",
        accession="KAN3063391.1",
        expected_organism="Cladosporium cladosporioides",
        expected_keyword="scytalone dehydratase",
        output_name="cladosporium_cladosporioides_SCD1_KAN3063391.1.fasta",
        notes="Ortholog referensi; langkah dehidrasi scytalone -> 1,3,8-THN.",
    ),
    FetchTarget(
        target_id="ccla_laccase",
        gene="Laccase",
        pathway="DHN-melanin",
        evidence_level="ortholog_reference",
        db="protein",
        accession="KAN3060494.1",
        expected_organism="Cladosporium cladosporioides",
        expected_keyword="laccase",
        output_name="cladosporium_cladosporioides_laccase_KAN3060494.1.fasta",
        notes="Ortholog referensi; oksidasi/polimerisasi DHN menjadi melanin.",
    ),
    FetchTarget(
        target_id="ccla_tyrosinase",
        gene="Tyrosinase",
        pathway="DOPA-melanin",
        evidence_level="ortholog_reference",
        db="protein",
        accession="KAN3055251.1",
        expected_organism="Cladosporium cladosporioides",
        expected_keyword="tyrosinase",
        output_name="cladosporium_cladosporioides_tyrosinase_KAN3055251.1.fasta",
        notes="Ortholog referensi; konversi L-tyrosine -> L-DOPA -> dopaquinone.",
    ),
    FetchTarget(
        target_id="chal_scytalone_dehydratase",
        gene="Scytalone dehydratase",
        pathway="DHN-melanin",
        evidence_level="ortholog_reference",
        db="protein",
        accession="XP_069226016.1",
        expected_organism="Cladosporium halotolerans",
        expected_keyword="scytalone dehydratase",
        output_name="cladosporium_halotolerans_SCD1_XP_069226016.1.fasta",
        notes="Ortholog referensi tambahan dari C. halotolerans.",
    ),
    FetchTarget(
        target_id="chal_thn_reductase",
        gene="THN reductase",
        pathway="DHN-melanin",
        evidence_level="ortholog_reference",
        db="protein",
        accession="XP_069230222.1",
        expected_organism="Cladosporium halotolerans",
        expected_keyword="tetrahydroxynaphthalene reductase",
        output_name="cladosporium_halotolerans_THN_reductase_XP_069230222.1.fasta",
        notes="Ortholog referensi tambahan dari C. halotolerans.",
    ),
)

# Metadata assembly yang relevan. Tidak diunduh otomatis karena tidak memiliki
# anotasi gen/protein lengkap di NCBI dan ukurannya besar.
REFERENCE_ASSEMBLIES: tuple[dict[str, str], ...] = (
    {
        "accession": "GCA_001931905.1",
        "organism": "Cladosporium sphaerospermum IMV 00045",
        "description": (
            "Draft genome isolat dari ISS; file genomic.gbff tidak memiliki CDS "
            "teranotasi (hanya fitur source)."
        ),
        "url": "https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_001931905.1/",
    },
    {
        "accession": "GCA_000261425.2",
        "organism": "Cladosporium sphaerospermum UM 843",
        "description": (
            "Assembly scaffold; direktori FTP NCBI tidak menyediakan protein FASTA "
            "maupun anotasi gen."
        ),
        "url": "https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_000261425.2/",
    },
)

CAVEATS: tuple[str, ...] = (
    "NCBI belum menyediakan set gen DHN-melanin teranotasi lengkap untuk "
    "Cladosporium sphaerospermum.",
    "Hanya MK213995.1 / QBZ73626.1 yang merupakan urutan langsung dari "
    "C. sphaerospermum; keduanya PKS tipe III parsial dan belum dikonfirmasi "
    "sebagai PKS jalur DHN-melanin.",
    "Target lain adalah ortholog referensi dari Cladosporium cladosporioides "
    "dan Cladosporium halotolerans, bukan urutan langsung C. sphaerospermum.",
    "Aksesi DQ000286 (Borrelia sp. 16S rRNA) yang sebelumnya tersimpan sebagai "
    "cladosporium_sample.fasta telah dipindahkan ke data/raw/rejected/.",
)


class FetchValidationError(ValueError):
    """Record NCBI tidak cocok dengan target yang diminta."""


def configure_entrez(email: str | None = None, api_key: str | None = None) -> str:
    """Set identitas Entrez. Kembalikan email yang dipakai."""
    resolved_email = email or os.getenv("NCBI_EMAIL") or DEFAULT_EMAIL
    if resolved_email == DEFAULT_EMAIL:
        print(
            "[WARN] NCBI_EMAIL belum diatur. Gunakan --email atau variabel "
            "lingkungan NCBI_EMAIL agar sesuai kebijakan NCBI."
        )
    Entrez.email = resolved_email
    Entrez.tool = TOOL_NAME

    resolved_key = api_key or os.getenv("NCBI_API_KEY")
    if resolved_key:
        Entrez.api_key = resolved_key
    return resolved_email


def _source_url(target: FetchTarget) -> str:
    if target.db == "nucleotide":
        return f"https://www.ncbi.nlm.nih.gov/nuccore/{target.accession}"
    return f"https://www.ncbi.nlm.nih.gov/protein/{target.accession}"


def _validate_record(record: SeqIO.SeqRecord, target: FetchTarget) -> None:
    """Pastikan record benar-benar organisme dan gen yang diminta."""
    description = str(record.description)
    haystack = description.lower()

    if target.expected_organism.lower() not in haystack:
        raise FetchValidationError(
            f"{target.accession}: organism tidak cocok. Diharapkan "
            f"'{target.expected_organism}', deskripsi NCBI: '{description}'"
        )
    if target.expected_keyword.lower() not in haystack:
        raise FetchValidationError(
            f"{target.accession}: kata kunci '{target.expected_keyword}' tidak "
            f"ditemukan pada deskripsi: '{description}'"
        )

    accession_root = target.accession.split(".")[0].upper()
    if accession_root not in str(record.id).upper():
        raise FetchValidationError(
            f"{target.accession}: record ID '{record.id}' tidak cocok dengan "
            "aksesi yang diminta."
        )
    if len(record.seq) < 30:
        raise FetchValidationError(
            f"{target.accession}: sekuens terlalu pendek ({len(record.seq)} bp/aa)."
        )


def _efetch_with_retry(target: FetchTarget, retries: int = 3) -> SeqIO.SeqRecord:
    last_error: Exception | None = None
    for attempt in range(1, retries + 1):
        handle = None
        try:
            handle = Entrez.efetch(
                db=target.db,
                id=target.accession,
                rettype="fasta",
                retmode="text",
            )
            record = SeqIO.read(handle, "fasta")
            _validate_record(record, target)
            return record
        except Exception as exc:  # network, parser, atau validasi
            last_error = exc
            if attempt < retries:
                time.sleep(float(attempt))
        finally:
            if handle is not None:
                handle.close()

    raise RuntimeError(
        f"Gagal mengambil {target.accession} setelah {retries} percobaan: {last_error}"
    ) from last_error


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _read_cached(path: Path, target: FetchTarget) -> SeqIO.SeqRecord | None:
    if not path.is_file():
        return None
    try:
        record = SeqIO.read(path, "fasta")
        _validate_record(record, target)
        return record
    except Exception:
        return None


def fetch_target(
    target: FetchTarget,
    raw_dir: Path = RAW_DIR,
    force: bool = False,
    dry_run: bool = False,
    retries: int = 3,
) -> dict[str, Any]:
    """Unduh satu target (atau pakai cache yang sudah tervalidasi)."""
    raw_dir.mkdir(parents=True, exist_ok=True)
    output_path = raw_dir / target.output_name

    base: dict[str, Any] = {
        "target_id": target.target_id,
        "gene": target.gene,
        "pathway": target.pathway,
        "evidence_level": target.evidence_level,
        "db": target.db,
        "accession": target.accession,
        "expected_organism": target.expected_organism,
        "notes": target.notes,
        "source_url": _source_url(target),
        "path": str(output_path.relative_to(PROJECT_ROOT)),
    }

    if dry_run:
        return {
            **base,
            "status": "planned",
            "action": "would_refetch" if output_path.exists() and force else "would_fetch",
        }

    cached = None if force else _read_cached(output_path, target)
    if cached is not None:
        return {
            **base,
            "status": "cached",
            "record_id": cached.id,
            "description": cached.description,
            "length": len(cached.seq),
            "sha256": _sha256(output_path),
            "fetched_at": datetime.fromtimestamp(
                output_path.stat().st_mtime, tz=timezone.utc
            ).isoformat(),
        }

    record = _efetch_with_retry(target, retries=retries)

    temp_path = output_path.with_suffix(output_path.suffix + ".tmp")
    try:
        with temp_path.open("w", encoding="utf-8", newline="\n") as handle:
            SeqIO.write(record, handle, "fasta")
        temp_path.replace(output_path)
    finally:
        if temp_path.exists():
            temp_path.unlink()

    return {
        **base,
        "status": "fetched",
        "record_id": record.id,
        "description": record.description,
        "length": len(record.seq),
        "sha256": _sha256(output_path),
        "fetched_at": datetime.now(timezone.utc).isoformat(),
    }


def _select_targets(target_ids: Sequence[str] | None) -> tuple[FetchTarget, ...]:
    if not target_ids:
        return TARGETS

    by_id = {target.target_id: target for target in TARGETS}
    selected: list[FetchTarget] = []
    unknown: list[str] = []
    for target_id in target_ids:
        target = by_id.get(target_id)
        if target is None:
            unknown.append(target_id)
        else:
            selected.append(target)

    if unknown:
        known = ", ".join(sorted(by_id))
        raise SystemExit(f"Target tidak dikenal: {', '.join(unknown)}\nTarget valid: {known}")
    return tuple(selected)


def _write_json_atomic(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp_path = path.with_suffix(path.suffix + ".tmp")
    temp_path.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    temp_path.replace(path)


def _write_raw_readme(results: Iterable[dict[str, Any]], path: Path = RAW_README_PATH) -> None:
    rows = [result for result in results]
    lines = [
        "# Data mentah NCBI",
        "",
        "Direktori ini berisi sekuens referensi yang diunduh dan divalidasi oleh",
        "`src/ncbi_fetcher.py`. Lihat `manifest.json` untuk checksum dan metadata lengkap.",
        "",
        "## Target",
        "",
        "| Target | Gen | Aksesi | Organisme | Level bukti | Status | Panjang |",
        "|---|---|---|---|---|---|---|",
    ]
    for row in rows:
        length = row.get("length")
        lines.append(
            "| {target_id} | {gene} | {accession} | {organism} | {evidence} | "
            "{status} | {length} |".format(
                target_id=row.get("target_id", "-"),
                gene=row.get("gene", "-"),
                accession=row.get("accession", "-"),
                organism=row.get("expected_organism", "-"),
                evidence=EVIDENCE_LABELS.get(
                    row.get("evidence_level", ""), row.get("evidence_level", "-")
                ),
                status=row.get("status", "-"),
                length=f"{length:,}" if isinstance(length, int) else "-",
            )
        )

    lines.extend(
        [
            "",
            "## Catatan ilmiah",
            "",
            *[f"- {caveat}" for caveat in CAVEATS],
            "",
            "## Assembly referensi (tidak diunduh)",
            "",
            *[
                f"- `{assembly['accession']}` — {assembly['organism']}: "
                f"{assembly['description']}"
                for assembly in REFERENCE_ASSEMBLIES
            ],
            "",
            "## Rekaman yang ditolak",
            "",
            "- `rejected/DQ000286.1_borrelia_16S.fasta` — aksesi DQ000286 adalah",
            "  *Borrelia* sp. IK/23 16S rRNA, bukan *Cladosporium sphaerospermum*.",
            "",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    temp_path = path.with_suffix(path.suffix + ".tmp")
    temp_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    temp_path.replace(path)


def fetch_all(
    target_ids: Sequence[str] | None = None,
    raw_dir: Path = RAW_DIR,
    force: bool = False,
    dry_run: bool = False,
    email: str | None = None,
    api_key: str | None = None,
    retries: int = 3,
) -> list[dict[str, Any]]:
    """Unduh semua target terpilih dan tulis manifest + README data."""
    selected = _select_targets(target_ids)
    resolved_email = configure_entrez(email=email, api_key=api_key)

    results: list[dict[str, Any]] = []
    delay = 0.1 if (api_key or os.getenv("NCBI_API_KEY")) else 0.4

    for index, target in enumerate(selected, start=1):
        print(
            f"[{index}/{len(selected)}] {target.target_id}: "
            f"{target.db} {target.accession} ({target.expected_organism})"
        )
        try:
            result = fetch_target(
                target,
                raw_dir=raw_dir,
                force=force,
                dry_run=dry_run,
                retries=retries,
            )
            print(f"    -> {result['status']}")
        except Exception as exc:
            result = {
                "target_id": target.target_id,
                "gene": target.gene,
                "pathway": target.pathway,
                "evidence_level": target.evidence_level,
                "db": target.db,
                "accession": target.accession,
                "expected_organism": target.expected_organism,
                "notes": target.notes,
                "source_url": _source_url(target),
                "path": str((raw_dir / target.output_name).relative_to(PROJECT_ROOT)),
                "status": "error",
                "error": str(exc),
            }
            print(f"    -> ERROR: {exc}")
        results.append(result)

        if not dry_run and index < len(selected):
            time.sleep(delay)

    if not dry_run:
        manifest = {
            "generated_at": datetime.now(timezone.utc).isoformat(),
            "tool": TOOL_NAME,
            "entrez_email": resolved_email,
            "targets": results,
            "reference_assemblies": list(REFERENCE_ASSEMBLIES),
            "caveats": list(CAVEATS),
        }
        _write_json_atomic(raw_dir / MANIFEST_PATH.name, manifest)
        _write_raw_readme(results, raw_dir / RAW_README_PATH.name)

    return results


def _print_target_table() -> None:
    print(f"{'TARGET ID':34s} {'DB':11s} {'ACCESSION':15s} {'EVIDENCE':20s} ORGANISM")
    print("-" * 110)
    for target in TARGETS:
        print(
            f"{target.target_id:34s} {target.db:11s} {target.accession:15s} "
            f"{target.evidence_level:20s} {target.expected_organism}"
        )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Unduh sekuens referensi melanin dari NCBI dengan validasi organisme "
            "dan kata kunci. Cache yang sudah valid tidak diunduh ulang."
        )
    )
    parser.add_argument(
        "--targets",
        default="all",
        help="Daftar target dipisahkan koma, atau 'all' (default).",
    )
    parser.add_argument("--list", action="store_true", help="Tampilkan target lalu keluar.")
    parser.add_argument("--force", action="store_true", help="Unduh ulang meski cache ada.")
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Tampilkan rencana tanpa mengakses NCBI atau menulis file.",
    )
    parser.add_argument("--email", default=None, help="Email kontak NCBI (atau NCBI_EMAIL).")
    parser.add_argument("--api-key", default=None, help="API key NCBI (opsional).")
    parser.add_argument("--retries", type=int, default=3, help="Jumlah percobaan ulang.")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    if args.list:
        _print_target_table()
        return 0

    if args.targets.strip().lower() == "all":
        target_ids = None
    else:
        target_ids = [item.strip() for item in args.targets.split(",") if item.strip()]

    results = fetch_all(
        target_ids=target_ids,
        force=args.force,
        dry_run=args.dry_run,
        email=args.email,
        api_key=args.api_key,
        retries=max(1, args.retries),
    )

    fetched = sum(1 for row in results if row["status"] == "fetched")
    cached = sum(1 for row in results if row["status"] == "cached")
    errors = sum(1 for row in results if row["status"] == "error")
    planned = sum(1 for row in results if row["status"] == "planned")

    if args.dry_run:
        print(f"\nDry run: {planned} target direncanakan.")
        return 0

    print(f"\nSelesai: {fetched} diunduh, {cached} dari cache, {errors} gagal.")
    print(f"Manifest: {MANIFEST_PATH}")
    print(f"README data: {RAW_README_PATH}")
    return 0 if errors == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
