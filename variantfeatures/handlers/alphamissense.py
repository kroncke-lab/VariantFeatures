"""AlphaMissense local-file handler.

The AlphaMissense data file (~1.1GB compressed, ~4GB uncompressed) lists every
possible single amino-acid substitution for ~20,000 human proteins. Per-variant
lookups via per-row scans are wasteful — the right shape is a single pass over
the file matching against the set of pending jobs at once.

This module exposes `run_batch(db, ...)` that:
1. Reads all pending `alphamissense` jobs into a (chrom, pos, ref, alt) index;
2. Streams through `AlphaMissense_aa_substitutions.tsv.gz` once;
3. For each match, writes annotations_pathogenicity + a hgvs_p alias and marks
   the job done.
4. Variants present in the queue but absent from the file are marked failed
   (e.g. the variant isn't a missense — AM only covers missense).

File format (header-prefixed by `#`):
    #CHROM  POS  REF  ALT  genome  uniprot_id  transcript_id  protein_variant  am_pathogenicity  am_class

Configure via `ALPHAMISSENSE_FILE` env var, or pass `file_path=` to run_batch.
Default search path: `data/alphamissense/AlphaMissense_aa_substitutions.tsv.gz`.

Download: https://storage.googleapis.com/dm_alphamissense/AlphaMissense_aa_substitutions.tsv.gz
"""

from __future__ import annotations

import gzip
import os
from pathlib import Path
from typing import Optional

from .clingen_ar import normalize_chromosome


SOURCE = "alphamissense"

DEFAULT_REL_PATH = Path("data") / "alphamissense" / "AlphaMissense_aa_substitutions.tsv.gz"
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent


class HandlerError(Exception):
    pass


def resolve_file_path(explicit: Optional[str] = None) -> Path:
    """Pick the AlphaMissense TSV path: explicit arg > env var > project default."""
    if explicit:
        return Path(explicit)
    env_path = os.environ.get("ALPHAMISSENSE_FILE")
    if env_path:
        return Path(env_path)
    return PROJECT_ROOT / DEFAULT_REL_PATH


def file_present(file_path: Optional[str] = None) -> bool:
    return resolve_file_path(file_path).exists()


# ---------------------------------------------------------------------------
# Batch runner
# ---------------------------------------------------------------------------

def run_batch(
    db,
    *,
    file_path: Optional[str] = None,
    limit: Optional[int] = None,
    progress_every: int = 1_000_000,
    progress_callback=None,
) -> dict:
    """Drain pending alphamissense jobs by single-pass over the AM TSV.

    Returns: {"claimed": int, "matched": int, "failed": int, "lines_scanned": int}
    """
    path = resolve_file_path(file_path)
    if not path.exists():
        raise HandlerError(
            f"AlphaMissense file not found at {path}. "
            f"Set ALPHAMISSENSE_FILE or download from "
            f"https://storage.googleapis.com/dm_alphamissense/AlphaMissense_aa_substitutions.tsv.gz"
        )

    # Pull pending jobs into memory, keyed by (chrom, pos, ref, alt).
    pending = _load_pending_index(db, limit=limit)
    if not pending:
        return {"claimed": 0, "matched": 0, "failed": 0, "lines_scanned": 0}

    matched = 0
    lines_scanned = 0

    with _open_maybe_gzip(path) as fh:
        header = None
        for line in fh:
            lines_scanned += 1
            if progress_callback and lines_scanned % progress_every == 0:
                progress_callback({"lines": lines_scanned, "matched": matched, "remaining": len(pending)})

            # AlphaMissense files start with copyright comments (#) then a header line that's
            # also `#`-prefixed (#CHROM  POS  ...). Treat the *first* commented line that splits
            # into the right number of fields as the header.
            if line.startswith("#"):
                if header is None:
                    fields = line.lstrip("#").strip().split("\t")
                    if len(fields) >= 9:
                        header = fields
                continue
            if header is None:
                # File without a header at all — assume the canonical column order.
                header = ["CHROM", "POS", "REF", "ALT", "genome", "uniprot_id", "transcript_id", "protein_variant", "am_pathogenicity", "am_class"]

            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            row = dict(zip(header, cols))
            try:
                key = (
                    normalize_chromosome(row["CHROM"]),
                    int(row["POS"]),
                    row["REF"],
                    row["ALT"],
                )
            except (KeyError, ValueError):
                continue

            jobs = pending.get(key)
            if not jobs:
                continue

            try:
                score = float(row["am_pathogenicity"])
            except (KeyError, ValueError):
                continue
            am_class = row.get("am_class") or None
            protein_variant = row.get("protein_variant") or ""
            transcript_id = row.get("transcript_id") or None

            # Map AlphaMissense class label -> short category code.
            category = _AM_CLASS_TO_CODE.get(am_class, am_class)

            for job in jobs:
                _persist_match(
                    db,
                    job=job,
                    score=score,
                    category=category,
                    protein_variant=protein_variant,
                    transcript_id=transcript_id,
                )
                matched += 1
            del pending[key]
            if not pending:
                break

    # Anything left in `pending` was not in the AM file (likely not a missense, or off-canonical
    # transcript). Mark those jobs failed with a clear message — the worker view stays accurate.
    failed = 0
    for jobs in pending.values():
        for job in jobs:
            db.mark_job_failed(job["job_id"], "no AlphaMissense entry for this (chrom,pos,ref,alt)")
            failed += 1

    return {"claimed": matched + failed, "matched": matched, "failed": failed, "lines_scanned": lines_scanned}


# ---------------------------------------------------------------------------
# Internals
# ---------------------------------------------------------------------------

_AM_CLASS_TO_CODE = {
    "likely_benign": "B",
    "ambiguous": "A",
    "likely_pathogenic": "P",
}


def _open_maybe_gzip(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def _load_pending_index(db, *, limit: Optional[int]) -> dict[tuple, list[dict]]:
    """Claim pending alphamissense jobs and index them by GRCh38 (chrom, pos, ref, alt)."""
    jobs = db.claim_pending_jobs(source=SOURCE, limit=limit if limit is not None else 1_000_000)
    if not jobs:
        return {}
    ids = [j["id"] for j in jobs]
    placeholders = ",".join("?" * len(ids))
    cur = db.conn.execute(
        f"""
        SELECT j.id AS job_id, j.variant_id AS variant_id,
               v.chromosome AS chromosome, v.position AS position, v.ref AS ref, v.alt AS alt
        FROM annotation_jobs j JOIN variants v ON v.id = j.variant_id
        WHERE j.id IN ({placeholders})
        """,
        ids,
    )
    out: dict[tuple, list[dict]] = {}
    for r in cur.fetchall():
        key = (normalize_chromosome(r["chromosome"]), int(r["position"]), r["ref"], r["alt"])
        out.setdefault(key, []).append(dict(r))
    return out


def _persist_match(db, *, job: dict, score: float, category: Optional[str], protein_variant: str, transcript_id: Optional[str]) -> None:
    db.upsert_pathogenicity(
        job["variant_id"],
        "alphamissense",
        score=score,
        category=category,
        source=SOURCE,
    )
    # AM's protein_variant looks like "KCNH2_A614V". Build a hgvs_p alias if we can.
    if protein_variant and "_" in protein_variant:
        aa_part = protein_variant.split("_", 1)[1]
        if transcript_id:
            db.add_aliases(
                job["variant_id"],
                [{"alias_type": "hgvs_p", "alias_value": f"{transcript_id}:p.{aa_part}", "source": SOURCE}],
            )
    db.mark_job_done(job["job_id"])
