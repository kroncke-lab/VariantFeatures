"""REVEL local-file batch handler for normalized variants.

REVEL v1.3 is distributed as a large CSV inside a zip archive:
``revel_with_transcript_ids``. The file is keyed by genomic coordinate and
Ensembl transcript ID, so the handler claims pending REVEL jobs, indexes their
canonical coordinates/transcripts, then scans the REVEL file once.
"""

from __future__ import annotations

import csv
import gzip
import io
import os
import zipfile
from contextlib import contextmanager
from pathlib import Path
from typing import Iterable, Optional, TextIO

from .clingen_ar import normalize_chromosome


SOURCE = "revel"
PREDICTOR_VERSION = "1.3"
SOURCE_LABEL = "revel_zenodo_7072866"
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
DEFAULT_CANDIDATES = (
    PROJECT_ROOT / "data" / "revel" / "revel_with_transcript_ids",
    PROJECT_ROOT / "data" / "revel" / "revel-v1.3_all_chromosomes.zip",
)


class HandlerError(Exception):
    pass


def resolve_file_path(explicit: Optional[str] = None) -> Path:
    if explicit:
        return Path(explicit)
    env_path = os.environ.get("REVEL_FILE")
    if env_path:
        return Path(env_path)
    for candidate in DEFAULT_CANDIDATES:
        if candidate.exists():
            return candidate
    return DEFAULT_CANDIDATES[0]


def file_present(file_path: Optional[str] = None) -> bool:
    return resolve_file_path(file_path).exists()


def run_batch(db, *, file_path: Optional[str] = None, limit: Optional[int] = None) -> dict:
    """Drain pending REVEL jobs by scanning the local REVEL file once."""
    path = resolve_file_path(file_path)
    if not path.exists():
        raise HandlerError(
            f"REVEL file not found at {path}. Set REVEL_FILE or place the "
            "revel_with_transcript_ids file/zip under data/revel/."
        )

    pending = _load_pending_index(db, limit=limit)
    claimed = len(pending["jobs"])
    if claimed == 0:
        return {"claimed": 0, "matched": 0, "failed": 0, "lines_scanned": 0}

    by_key = pending["by_key"]
    matched_job_ids: set[int] = set()
    matched = 0
    lines_scanned = 0

    with _open_revel(path) as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            lines_scanned += 1
            key = _record_key(row)
            if key is None:
                continue
            jobs = by_key.get(key)
            if not jobs:
                continue
            transcript = _transcript_base(row.get("Ensembl_transcriptid"))
            try:
                score = float(row["REVEL"])
            except (KeyError, TypeError, ValueError):
                continue
            for job in jobs:
                wanted_transcripts = job.get("transcript_bases") or set()
                if wanted_transcripts and transcript not in wanted_transcripts:
                    continue
                db.upsert_pathogenicity(
                    job["variant_id"],
                    "revel",
                    predictor_version=PREDICTOR_VERSION,
                    score=score,
                    source=SOURCE_LABEL,
                )
                db.mark_job_done(job["job_id"])
                matched_job_ids.add(job["job_id"])
                matched += 1

            if len(matched_job_ids) == claimed:
                break

    failed = 0
    for job in pending["jobs"]:
        if job["job_id"] in matched_job_ids:
            continue
        db.mark_job_failed(job["job_id"], "no REVEL entry for this coordinate/transcript")
        failed += 1

    return {"claimed": claimed, "matched": matched, "failed": failed, "lines_scanned": lines_scanned}


def parse_rows(path: str | Path) -> Iterable[dict]:
    """Yield raw REVEL CSV rows from a plain, gzipped, or zipped file."""
    with _open_revel(Path(path)) as fh:
        yield from csv.DictReader(fh)


@contextmanager
def _open_revel(path: Path) -> Iterable[TextIO]:
    if path.suffix == ".zip":
        with zipfile.ZipFile(path) as zf:
            names = [n for n in zf.namelist() if not n.endswith("/")]
            if not names:
                raise HandlerError(f"REVEL zip archive has no file entries: {path}")
            with zf.open(names[0], "r") as raw:
                yield io.TextIOWrapper(raw, encoding="utf-8", newline="")
        return
    if path.suffix == ".gz":
        with gzip.open(path, "rt", newline="") as fh:
            yield fh
        return
    with open(path, "r", newline="") as fh:
        yield fh


def _load_pending_index(db, *, limit: Optional[int]) -> dict:
    jobs = db.claim_pending_jobs(source=SOURCE, limit=limit if limit is not None else 1_000_000)
    if not jobs:
        return {"jobs": [], "by_key": {}}

    ids = [j["id"] for j in jobs]
    placeholders = ",".join("?" * len(ids))
    cur = db.conn.execute(
        f"""
        SELECT
            j.id AS job_id,
            j.variant_id AS variant_id,
            v.chromosome AS chromosome,
            v.position AS position,
            v.ref AS ref,
            v.alt AS alt,
            c.transcript_id AS transcript_id,
            c.consequence AS consequence
        FROM annotation_jobs j
        JOIN variants v ON v.id = j.variant_id
        LEFT JOIN variant_consequences c
            ON c.variant_id = j.variant_id
           AND c.source = 'enumerated'
        WHERE j.id IN ({placeholders})
        """,
        ids,
    )

    by_job: dict[int, dict] = {}
    by_key: dict[tuple, list[dict]] = {}
    for row in cur.fetchall():
        d = dict(row)
        job = by_job.setdefault(
            d["job_id"],
            {
                "job_id": d["job_id"],
                "variant_id": d["variant_id"],
                "transcript_bases": set(),
            },
        )
        if d.get("transcript_id"):
            job["transcript_bases"].add(_transcript_base(d["transcript_id"]))
        key = (
            normalize_chromosome(d["chromosome"]),
            int(d["position"]),
            d["ref"],
            d["alt"],
        )
        by_key.setdefault(key, []).append(job)

    return {"jobs": list(by_job.values()), "by_key": by_key}


def _record_key(row: dict) -> Optional[tuple]:
    chrom = row.get("chr")
    pos = row.get("grch38_pos")
    ref = row.get("ref")
    alt = row.get("alt")
    if not (chrom and pos and ref and alt):
        return None
    try:
        return (normalize_chromosome(chrom), int(pos), ref, alt)
    except (TypeError, ValueError):
        return None


def _transcript_base(transcript_id: Optional[str]) -> str:
    return str(transcript_id or "").split(".", 1)[0]
