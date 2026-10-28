"""Validasi homology ortholog melanin terhadap assembly Cladosporium.

Modul ini menjalankan tblastn (protein query vs enam-frame translation genom)
untuk memetakan tujuh ortholog referensi dari *C. cladosporioides* dan
*C. halotolerans* ke assembly *C. sphaerospermum*:

- GCA_001931905.1 (assembly lama, digantikan oleh .2)
- GCA_001931905.2 (assembly saat ini)
- GCA_000261425.2 (strain UM 843)

Hasil BLAST adalah bukti homologi, bukan bukti fungsi. Hit strong/moderate
menunjukkan bahwa sekuens mirip hadir di assembly; tidak adanya hit tidak
membuktikan gen tersebut absent (gap assembly, divergensi, atau query parsial).

Catatan operasional:
- Genome dan BLAST+ diunduh ke folder yang di-ignore Git (`data/genomes/`,
  `tools/`). File hasil ringkas disimpan di `results/`.
- Checksum md5 dari NCBI diverifikasi bila tersedia.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import shutil
import subprocess
import time
import urllib.error
import urllib.request
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Sequence

from Bio import SeqIO

PROJECT_ROOT = Path(__file__).resolve().parents[1]
RAW_DIR = PROJECT_ROOT / "data" / "raw"
MANIFEST_PATH = RAW_DIR / "manifest.json"
GENOME_DIR = PROJECT_ROOT / "data" / "genomes"
DB_DIR = GENOME_DIR / "blast_db"
RESULTS_DIR = PROJECT_ROOT / "results"
DEFAULT_JSON = RESULTS_DIR / "blast_validation.json"
DEFAULT_TSV = RESULTS_DIR / "blast_hits.tsv"
DEFAULT_QUERIES = RESULTS_DIR / "blast_queries.faa"
BLAST_TOOLS_DIR = PROJECT_ROOT / "tools" / "blast"

USER_AGENT = "melanin-biotech-blast-validation/1.0"

TBLASTN_FIELDS = (
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "qlen",
    "slen",
)

TBLASTN_OUTFMT = "6 " + " ".join(TBLASTN_FIELDS)


class BlastValidationError(RuntimeError):
    """Kesalahan yang dapat ditampilkan ke pengguna tanpa traceback."""


@dataclass(frozen=True)
class Assembly:
    """Satu assembly genom yang akan dipakai sebagai target BLAST."""

    accession: str
    assembly_name: str
    organism: str
    strain: str
    url: str
    note: str = ""
    superseded_by: str | None = None

    @property
    def filename(self) -> str:
        return self.url.rsplit("/", 1)[-1]


ASSEMBLIES: tuple[Assembly, ...] = (
    Assembly(
        accession="GCA_001931905.1",
        assembly_name="ASM193190v1",
        organism="Cladosporium sphaerospermum",
        strain="IMV 00045",
        url=(
            "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/931/905/"
            "GCA_001931905.1_ASM193190v1/"
            "GCA_001931905.1_ASM193190v1_genomic.fna.gz"
        ),
        note=(
            "Assembly yang disebut pada rencana awal; status NCBI 'previous' "
            "dan telah digantikan oleh GCA_001931905.2."
        ),
        superseded_by="GCA_001931905.2",
    ),
    Assembly(
        accession="GCA_001931905.2",
        assembly_name="ASM193190v2",
        organism="Cladosporium sphaerospermum",
        strain="IMV 00045",
        url=(
            "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/931/905/"
            "GCA_001931905.2_ASM193190v2/"
            "GCA_001931905.2_ASM193190v2_genomic.fna.gz"
        ),
        note="Assembly saat ini (current) untuk strain IMV 00045.",
    ),
    Assembly(
        accession="GCA_000261425.2",
        assembly_name="UM843Vel2.0",
        organism="Cladosporium sphaerospermum",
        strain="UM 843",
        url=(
            "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/261/425/"
            "GCA_000261425.2_UM843Vel2.0/"
            "GCA_000261425.2_UM843Vel2.0_genomic.fna.gz"
        ),
        note="Assembly alternatif strain UM 843.",
    ),
)


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


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


def _md5_file(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _download(url: str, dest: Path, timeout: int = 180, retries: int = 3) -> None:
    """Unduh file dengan retry dan penulisan atomik via file .part."""

    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    last_error: Exception | None = None

    for attempt in range(1, retries + 1):
        try:
            request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(request, timeout=timeout) as response:
                with tmp.open("wb") as handle:
                    shutil.copyfileobj(response, handle, length=1024 * 1024)
            if tmp.stat().st_size == 0:
                raise BlastValidationError(f"Unduhan kosong: {url}")
            tmp.replace(dest)
            return
        except (OSError, urllib.error.URLError) as exc:
            last_error = exc
            tmp.unlink(missing_ok=True)
            if attempt < retries:
                time.sleep(2 * attempt)

    raise BlastValidationError(f"Gagal mengunduh {url}: {last_error}")


def _fetch_expected_md5(assembly: Assembly) -> str | None:
    """Ambil md5 dari md5checksums.txt NCBI untuk file genom tertentu."""

    checksum_url = assembly.url.rsplit("/", 1)[0] + "/md5checksums.txt"
    try:
        request = urllib.request.Request(
            checksum_url, headers={"User-Agent": USER_AGENT}
        )
        with urllib.request.urlopen(request, timeout=60) as response:
            text = response.read().decode("utf-8", errors="replace")
    except (OSError, urllib.error.URLError):
        return None

    for line in text.splitlines():
        parts = line.split()
        if len(parts) >= 2 and parts[-1].lstrip("./") == assembly.filename:
            return parts[0].lower()
    return None


def _decompress_gzip(src: Path, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    with gzip.open(src, "rb") as source, tmp.open("wb") as target:
        shutil.copyfileobj(source, target, length=1024 * 1024)
    tmp.replace(dest)


def _fasta_stats(path: Path) -> tuple[int, int]:
    sequences = 0
    total_bp = 0
    for record in SeqIO.parse(str(path), "fasta"):
        sequences += 1
        total_bp += len(record.seq)
    if sequences == 0:
        raise BlastValidationError(f"FASTA kosong atau tidak valid: {path}")
    return sequences, total_bp


def ensure_genome(
    assembly: Assembly,
    *,
    skip_download: bool = False,
    force_download: bool = False,
) -> dict[str, Any]:
    """Pastikan file genom dan provenance tersedia; kembalikan metadata."""

    out_dir = GENOME_DIR / assembly.accession
    archive = out_dir / assembly.filename
    genome = out_dir / "genome.fna"
    provenance_path = out_dir / "provenance.json"

    if provenance_path.exists() and genome.exists() and not force_download:
        try:
            cached = json.loads(provenance_path.read_text(encoding="utf-8"))
            if cached.get("sequences") and cached.get("length_bp"):
                cached["genome_path"] = str(genome)
                return cached
        except (OSError, json.JSONDecodeError):
            pass

    if force_download or not archive.exists():
        if skip_download:
            raise BlastValidationError(
                f"--skip-download diberikan tetapi {archive} tidak ada."
            )
        print(f"    mengunduh {assembly.accession} ...")
        _download(assembly.url, archive)

    expected_md5 = _fetch_expected_md5(assembly)
    observed_md5 = _md5_file(archive)
    if expected_md5 and expected_md5 != observed_md5:
        archive.unlink(missing_ok=True)
        raise BlastValidationError(
            f"md5 tidak cocok untuk {assembly.accession}: "
            f"diharapkan {expected_md5}, diperoleh {observed_md5}"
        )

    if force_download or not genome.exists():
        print(f"    mengekstrak {assembly.accession} ...")
        _decompress_gzip(archive, genome)

    sequences, total_bp = _fasta_stats(genome)
    provenance = {
        "accession": assembly.accession,
        "assembly_name": assembly.assembly_name,
        "organism": assembly.organism,
        "strain": assembly.strain,
        "source_url": assembly.url,
        "note": assembly.note,
        "superseded_by": assembly.superseded_by,
        "archive_path": str(archive),
        "genome_path": str(genome),
        "archive_size_bytes": archive.stat().st_size,
        "md5_expected": expected_md5,
        "md5_observed": observed_md5,
        "md5_verified": bool(expected_md5 and expected_md5 == observed_md5),
        "sequences": sequences,
        "length_bp": total_bp,
        "prepared_at": _now(),
    }
    _write_json_atomic(provenance_path, provenance)
    return provenance


def _exe(directory: Path, name: str) -> Path:
    suffix = ".exe" if os.name == "nt" else ""
    return directory / f"{name}{suffix}"


def find_blast_tools(explicit: Path | None = None) -> dict[str, Path]:
    """Cari tblastn dan makeblastdb dari argumen, PATH, atau tools/."""

    candidates: list[Path] = []
    if explicit is not None:
        candidates.append(explicit)
    env_home = os.environ.get("BLAST_HOME")
    if env_home:
        candidates.append(Path(env_home) / "bin")
    if BLAST_TOOLS_DIR.exists():
        candidates.extend(sorted(BLAST_TOOLS_DIR.glob("ncbi-blast-*/bin")))
    found = shutil.which("tblastn")
    if found:
        candidates.append(Path(found).parent)

    checked: list[str] = []
    for directory in candidates:
        tblastn = _exe(directory, "tblastn")
        makeblastdb = _exe(directory, "makeblastdb")
        checked.append(str(directory))
        if tblastn.exists() and makeblastdb.exists():
            return {"tblastn": tblastn, "makeblastdb": makeblastdb}

    raise BlastValidationError(
        "BLAST+ tidak ditemukan. Letakkan di tools/blast/ncbi-blast-*/bin, "
        "atur BLAST_HOME, atau gunakan --blast-bin. Lokasi diperiksa: "
        + ", ".join(checked)
    )


def _blast_version(tblastn: Path) -> str:
    process = subprocess.run(
        [str(tblastn), "-version"],
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        check=False,
    )
    output = (process.stdout or process.stderr).strip()
    return output.splitlines()[0] if output else "unknown"


def build_blast_db(
    assembly: Assembly,
    genome_path: Path,
    makeblastdb: Path,
    *,
    force: bool = False,
) -> dict[str, Any]:
    """Buat database nukleotida BLAST untuk satu assembly."""

    db_dir = DB_DIR / assembly.accession
    prefix = db_dir / assembly.accession
    marker = Path(str(prefix) + ".nsq")

    db_dir.mkdir(parents=True, exist_ok=True)
    if marker.exists() and not force:
        return {
            "accession": assembly.accession,
            "db_prefix": str(prefix),
            "built": False,
        }

    command = [
        str(makeblastdb),
        "-in",
        str(genome_path),
        "-dbtype",
        "nucl",
        "-out",
        str(prefix),
        "-title",
        f"{assembly.organism} {assembly.strain} {assembly.accession}",
    ]
    process = subprocess.run(
        command,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        check=False,
    )
    if process.returncode != 0 or not marker.exists():
        raise BlastValidationError(
            f"makeblastdb gagal untuk {assembly.accession}: "
            f"{(process.stderr or process.stdout).strip()}"
        )

    return {
        "accession": assembly.accession,
        "db_prefix": str(prefix),
        "built": True,
        "makeblastdb_stdout": process.stdout.strip(),
    }


def _resolve_manifest_path(value: str) -> Path:
    path = Path(value.replace("\\", "/"))
    if not path.is_absolute():
        path = PROJECT_ROOT / path
    return path


def select_queries(
    *,
    include_direct: bool = False,
    output_faa: Path = DEFAULT_QUERIES,
) -> list[dict[str, Any]]:
    """Pilih query protein ortholog (opsional: termasuk urutan langsung)."""

    if not MANIFEST_PATH.exists():
        raise BlastValidationError(f"Manifest tidak ditemukan: {MANIFEST_PATH}")

    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    selected: list[dict[str, Any]] = []

    for row in manifest.get("targets", []):
        evidence = row.get("evidence_level")
        if row.get("db") != "protein":
            continue
        if evidence != "ortholog_reference" and not (
            include_direct and evidence == "direct_sequence"
        ):
            continue

        path = _resolve_manifest_path(str(row.get("path", "")))
        if not path.exists():
            raise BlastValidationError(f"File query tidak ditemukan: {path}")

        records = list(SeqIO.parse(str(path), "fasta"))
        if len(records) != 1:
            raise BlastValidationError(
                f"Query harus berisi tepat satu record: {path} "
                f"(ditemukan {len(records)})"
            )
        record = records[0]
        query_id = f"{row['target_id']}|{row['accession']}"
        selected.append(
            {
                "target_id": row["target_id"],
                "query_id": query_id,
                "accession": row["accession"],
                "gene": row.get("gene", ""),
                "pathway": row.get("pathway", ""),
                "evidence_level": evidence,
                "length": len(record.seq),
                "path": str(path),
                "sequence": str(record.seq),
            }
        )

    if not selected:
        raise BlastValidationError("Tidak ada query protein yang dipilih.")

    output_faa.parent.mkdir(parents=True, exist_ok=True)
    with output_faa.open("w", encoding="utf-8", newline="\n") as handle:
        for query in selected:
            handle.write(f">{query['query_id']}\n")
            sequence = query["sequence"]
            for start in range(0, len(sequence), 60):
                handle.write(sequence[start : start + 60] + "\n")

    return selected


def run_tblastn(
    query_faa: Path,
    db_prefix: Path,
    tblastn: Path,
    out_tsv: Path,
    *,
    evalue: float,
    max_target_seqs: int,
    max_hsps: int,
    threads: int,
) -> str:
    """Jalankan tblastn dan kembalikan stderr (untuk provenance)."""

    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    command = [
        str(tblastn),
        "-query",
        str(query_faa),
        "-db",
        str(db_prefix),
        "-out",
        str(out_tsv),
        "-outfmt",
        TBLASTN_OUTFMT,
        "-evalue",
        f"{evalue:g}",
        "-max_target_seqs",
        str(max_target_seqs),
        "-max_hsps",
        str(max_hsps),
        "-num_threads",
        str(threads),
        "-seg",
        "yes",
    ]
    process = subprocess.run(
        command,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        check=False,
    )
    stderr = (process.stderr or "").strip()
    if process.returncode != 0:
        raise BlastValidationError(
            f"tblastn gagal untuk {db_prefix.name}: "
            f"{stderr or process.stdout.strip()}"
        )
    return stderr


def parse_hits(path: Path, assembly: str) -> list[dict[str, Any]]:
    """Parse output tabular tblastn menjadi list dict."""

    hits: list[dict[str, Any]] = []
    if not path.exists():
        return hits

    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) != len(TBLASTN_FIELDS):
            continue
        row = dict(zip(TBLASTN_FIELDS, parts))
        try:
            pident = float(row["pident"])
            alignment_length = int(row["length"])
            qlen = int(row["qlen"])
            qstart = int(row["qstart"])
            qend = int(row["qend"])
            sstart = int(row["sstart"])
            send = int(row["send"])
            evalue = float(row["evalue"])
            bitscore = float(row["bitscore"])
        except ValueError:
            continue

        query_coverage = (abs(qend - qstart) + 1) / qlen if qlen else 0.0
        hits.append(
            {
                "assembly": assembly,
                "query_id": row["qseqid"],
                "subject_id": row["sseqid"],
                "pident": pident,
                "alignment_length": alignment_length,
                "mismatch": int(row["mismatch"]),
                "gapopen": int(row["gapopen"]),
                "query_start": qstart,
                "query_end": qend,
                "subject_start": min(sstart, send),
                "subject_end": max(sstart, send),
                "strand": "plus" if sstart <= send else "minus",
                "evalue": evalue,
                "bitscore": bitscore,
                "query_length": qlen,
                "subject_length": int(row["slen"]),
                "query_coverage": round(query_coverage, 4),
            }
        )
    return hits


def classify_hit(hit: dict[str, Any] | None) -> str:
    """Klasifikasi heuristik kekuatan bukti homologi."""

    if hit is None:
        return "not_detected"
    if (
        hit["evalue"] <= 1e-20
        and hit["query_coverage"] >= 0.70
        and hit["pident"] >= 30.0
    ):
        return "strong"
    if (
        hit["evalue"] <= 1e-5
        and hit["query_coverage"] >= 0.50
        and hit["pident"] >= 25.0
    ):
        return "moderate"
    return "weak"


def _best_hit(rows: Iterable[dict[str, Any]]) -> dict[str, Any] | None:
    rows = list(rows)
    if not rows:
        return None
    return min(
        rows,
        key=lambda hit: (
            hit["evalue"],
            -hit["bitscore"],
            -hit["pident"],
            -hit["query_coverage"],
        ),
    )


def build_best_hits(
    queries: list[dict[str, Any]],
    assemblies: list[dict[str, Any]],
    hits: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    """Ambil hit terbaik per kombinasi query x assembly."""

    grouped: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for hit in hits:
        grouped[(hit["query_id"], hit["assembly"])].append(hit)

    best: list[dict[str, Any]] = []
    for query in queries:
        for assembly in assemblies:
            key = (query["query_id"], assembly["accession"])
            hit = _best_hit(grouped.get(key, []))
            classification = classify_hit(hit)
            if hit is None:
                best.append(
                    {
                        "query_id": query["query_id"],
                        "query_accession": query["accession"],
                        "gene": query["gene"],
                        "query_length": query["length"],
                        "assembly": assembly["accession"],
                        "classification": classification,
                        "subject_id": None,
                        "pident": None,
                        "alignment_length": None,
                        "query_coverage": None,
                        "evalue": None,
                        "bitscore": None,
                        "subject_start": None,
                        "subject_end": None,
                        "strand": None,
                    }
                )
            else:
                best.append(
                    {
                        "query_id": query["query_id"],
                        "query_accession": query["accession"],
                        "gene": query["gene"],
                        "query_length": query["length"],
                        "assembly": assembly["accession"],
                        "classification": classification,
                        "subject_id": hit["subject_id"],
                        "pident": hit["pident"],
                        "alignment_length": hit["alignment_length"],
                        "query_coverage": hit["query_coverage"],
                        "evalue": hit["evalue"],
                        "bitscore": hit["bitscore"],
                        "subject_start": hit["subject_start"],
                        "subject_end": hit["subject_end"],
                        "strand": hit["strand"],
                    }
                )
    return best


def write_hits_tsv(path: Path, hits: list[dict[str, Any]]) -> None:
    """Tulis semua hit mentah ke TSV."""

    columns = (
        "assembly",
        "query_id",
        "subject_id",
        "pident",
        "alignment_length",
        "mismatch",
        "gapopen",
        "query_start",
        "query_end",
        "subject_start",
        "subject_end",
        "strand",
        "evalue",
        "bitscore",
        "query_length",
        "subject_length",
        "query_coverage",
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("\t".join(columns) + "\n")
        for hit in sorted(
            hits,
            key=lambda row: (
                row["assembly"],
                row["query_id"],
                row["evalue"],
                -row["bitscore"],
            ),
        ):
            handle.write(
                "\t".join(str(hit.get(column, "")) for column in columns) + "\n"
            )


def build_summary(
    queries: list[dict[str, Any]],
    assemblies: list[dict[str, Any]],
    best_hits: list[dict[str, Any]],
) -> dict[str, Any]:
    """Ringkas klasifikasi per query dan per assembly."""

    classification_counts = Counter(hit["classification"] for hit in best_hits)
    per_assembly: dict[str, dict[str, int]] = {}
    for assembly in assemblies:
        accession = assembly["accession"]
        rows = [hit for hit in best_hits if hit["assembly"] == accession]
        per_assembly[accession] = dict(Counter(hit["classification"] for hit in rows))

    per_query: dict[str, dict[str, str]] = {}
    for query in queries:
        rows = [hit for hit in best_hits if hit["query_id"] == query["query_id"]]
        per_query[query["query_id"]] = {
            hit["assembly"]: hit["classification"] for hit in rows
        }

    queries_with_evidence = sorted(
        {
            hit["query_id"]
            for hit in best_hits
            if hit["classification"] in {"strong", "moderate"}
        }
    )
    queries_without_evidence = sorted(
        query["query_id"]
        for query in queries
        if query["query_id"] not in queries_with_evidence
    )

    return {
        "classifications": dict(classification_counts),
        "per_assembly": per_assembly,
        "per_query": per_query,
        "queries_with_evidence": queries_with_evidence,
        "queries_without_evidence": queries_without_evidence,
        "combinations_total": len(best_hits),
    }


def build_interpretation(summary: dict[str, Any]) -> str:
    total = summary["combinations_total"]
    counts = summary["classifications"]
    strong = counts.get("strong", 0)
    moderate = counts.get("moderate", 0)
    weak = counts.get("weak", 0)
    absent = counts.get("not_detected", 0)
    with_evidence = summary["queries_with_evidence"]
    without_evidence = summary["queries_without_evidence"]

    parts = [
        (
            f"Dari {total} kombinasi query x assembly, {strong} menunjukkan hit "
            f"strong, {moderate} moderate, {weak} weak, dan {absent} tidak "
            "terdeteksi."
        ),
        (
            "Query dengan bukti strong/moderate di setidaknya satu assembly: "
            + (", ".join(with_evidence) if with_evidence else "tidak ada")
            + "."
        ),
    ]
    if without_evidence:
        parts.append(
            "Query tanpa bukti strong/moderate: "
            + ", ".join(without_evidence)
            + ". Ini tidak membuktikan gen absent; dapat disebabkan gap assembly, "
            "divergensi sekuens, atau query yang parsial."
        )
    parts.append(
        "Status ini adalah bukti homologi sekuens, bukan konfirmasi fungsi "
        "enzimatik atau ekspresi gen."
    )
    return " ".join(parts)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Validasi homology ortholog melanin terhadap assembly "
            "C. sphaerospermum dengan tblastn."
        )
    )
    parser.add_argument(
        "--assemblies",
        nargs="+",
        choices=[assembly.accession for assembly in ASSEMBLIES],
        default=None,
        help="Assembly yang dipakai (default: ketiga assembly).",
    )
    parser.add_argument(
        "--include-direct",
        action="store_true",
        help="Sertakan juga query protein langsung C. sphaerospermum.",
    )
    parser.add_argument(
        "--skip-download",
        action="store_true",
        help="Gunakan genome yang sudah ada; jangan unduh.",
    )
    parser.add_argument(
        "--force-download",
        action="store_true",
        help="Unduh ulang genome meskipun sudah ada.",
    )
    parser.add_argument(
        "--force-db",
        action="store_true",
        help="Bangun ulang database BLAST.",
    )
    parser.add_argument("--evalue", type=float, default=1e-5)
    parser.add_argument("--max-target-seqs", type=int, default=25)
    parser.add_argument("--max-hsps", type=int, default=10)
    parser.add_argument(
        "--threads",
        type=int,
        default=min(8, os.cpu_count() or 1),
        help="Jumlah thread tblastn.",
    )
    parser.add_argument(
        "--blast-bin",
        type=Path,
        default=None,
        help="Direktori bin BLAST+ (berisi tblastn dan makeblastdb).",
    )
    parser.add_argument("--out", type=Path, default=DEFAULT_JSON)
    parser.add_argument("--tsv", type=Path, default=DEFAULT_TSV)
    parser.add_argument("--queries", type=Path, default=DEFAULT_QUERIES)
    parser.add_argument(
        "--list",
        action="store_true",
        help="Tampilkan assembly dan query yang akan dipakai, lalu keluar.",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    selected_accessions = args.assemblies or [
        assembly.accession for assembly in ASSEMBLIES
    ]
    selected_assemblies = [
        assembly
        for assembly in ASSEMBLIES
        if assembly.accession in selected_accessions
    ]

    if args.list:
        print("Assembly:")
        for assembly in selected_assemblies:
            print(
                f"  - {assembly.accession}  {assembly.organism} "
                f"{assembly.strain}  ({assembly.assembly_name})"
            )
        print("\nQuery:")
        for query in select_queries(
            include_direct=args.include_direct, output_faa=args.queries
        ):
            print(
                f"  - {query['query_id']}  {query['gene']}  "
                f"({query['length']} aa)"
            )
        return 0

    try:
        blast_tools = find_blast_tools(args.blast_bin)
        queries = select_queries(
            include_direct=args.include_direct, output_faa=args.queries
        )
        print(f"[1/4] {len(queries)} query protein disiapkan: {args.queries}")

        prepared_assemblies: list[dict[str, Any]] = []
        print("[2/4] Menyiapkan genome assembly...")
        for assembly in selected_assemblies:
            provenance = ensure_genome(
                assembly,
                skip_download=args.skip_download,
                force_download=args.force_download,
            )
            prepared_assemblies.append(provenance)
            print(
                f"    {assembly.accession}: {provenance['sequences']:,} "
                f"sekuens, {provenance['length_bp']:,} bp, "
                f"md5_verified={provenance['md5_verified']}"
            )

        print("[3/4] Membangun database BLAST dan menjalankan tblastn...")
        all_hits: list[dict[str, Any]] = []
        db_details: list[dict[str, Any]] = []
        for assembly, provenance in zip(selected_assemblies, prepared_assemblies):
            db_info = build_blast_db(
                assembly,
                Path(provenance["genome_path"]),
                blast_tools["makeblastdb"],
                force=args.force_db,
            )
            db_details.append(db_info)
            per_assembly_tsv = args.tsv.parent / "blast" / (
                f"{assembly.accession}.tsv"
            )
            run_tblastn(
                args.queries,
                Path(db_info["db_prefix"]),
                blast_tools["tblastn"],
                per_assembly_tsv,
                evalue=args.evalue,
                max_target_seqs=args.max_target_seqs,
                max_hsps=args.max_hsps,
                threads=args.threads,
            )
            assembly_hits = parse_hits(per_assembly_tsv, assembly.accession)
            all_hits.extend(assembly_hits)
            print(
                f"    {assembly.accession}: {len(assembly_hits)} hit "
                f"(evalue <= {args.evalue:g})"
            )

        best_hits = build_best_hits(queries, prepared_assemblies, all_hits)
        summary = build_summary(queries, prepared_assemblies, best_hits)
        write_hits_tsv(args.tsv, all_hits)

        payload = {
            "generated_at": _now(),
            "tool": "melanin_biotech_blast_validation",
            "method": (
                "tblastn protein query vs six-frame translation of genomic "
                "assembly nucleotide BLAST database"
            ),
            "blast_version": _blast_version(blast_tools["tblastn"]),
            "parameters": {
                "evalue": args.evalue,
                "max_target_seqs": args.max_target_seqs,
                "max_hsps": args.max_hsps,
                "threads": args.threads,
                "seg": "yes",
                "outfmt": TBLASTN_OUTFMT,
            },
            "assemblies": prepared_assemblies,
            "database_builds": db_details,
            "queries": [
                {key: value for key, value in query.items() if key != "sequence"}
                for query in queries
            ],
            "hits_tsv": str(args.tsv),
            "hits_total": len(all_hits),
            "best_hits": best_hits,
            "summary": summary,
            "classification_rules": {
                "strong": "evalue <= 1e-20 AND query_coverage >= 0.70 AND pident >= 30",
                "moderate": "evalue <= 1e-5 AND query_coverage >= 0.50 AND pident >= 25",
                "weak": "hit di bawah ambang strong/moderate",
                "not_detected": "tidak ada hit pada evalue yang dipakai",
            },
            "assumptions": [
                "Query adalah ortholog referensi dari C. cladosporioides/C. halotolerans, bukan sekuens C. sphaerospermum.",
                "tblastn membandingkan protein query dengan enam frame translation genom; hit adalah bukti homologi, bukan bukti fungsi.",
                "Ambang strong/moderate adalah heuristik untuk triase, bukan standar filogenetik.",
                "Query parsial (mis. PKS tipe III) memiliki coverage yang secara intrinsik terbatas.",
                "Database BLAST dibangun ulang dari FASTA NCBI yang checksum-nya diverifikasi bila tersedia.",
            ],
            "caveats": [
                "Tidak adanya hit dapat disebabkan gap assembly, divergensi sekuens, atau query parsial; bukan bukti gen absent.",
                "Koordinat subject bersifat spesifik assembly dan tidak dinormalisasi antar assembly.",
                "Hasil ini tidak membuktikan ekspresi, aktivitas enzim, atau keterlibatan jalur DHN-melanin.",
                "GCA_001931905.1 berstatus previous dan digantikan oleh GCA_001931905.2; keduanya disertakan untuk transparansi.",
            ],
            "interpretation": build_interpretation(summary),
        }
        _write_json_atomic(args.out, _portable_paths(payload))

        print("[4/4] Selesai.")
        print(f"    JSON : {args.out}")
        print(f"    TSV  : {args.tsv}")
        print(f"    {payload['interpretation']}")
        return 0

    except BlastValidationError as exc:
        print(f"[ERROR] {exc}")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
