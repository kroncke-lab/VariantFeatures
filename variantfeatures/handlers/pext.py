"""gnomAD pext local-file importer.

gnomAD's current pext data is distributed through browser/Hail-table downloads,
so this module intentionally consumes a local, region-filtered CSV/TSV export
instead of trying to hide that large-data dependency.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Optional

from .tabular_utils import (
    COORD_COLUMNS,
    coerce_float,
    first_present,
    load_variant_index,
    match_variant_ids,
    parse_variant_key,
    read_dicts,
)


SOURCE = "pext"
DEFAULT_DATASET = "gnomad_pext_v10"
UCSC_BIGWIG_DATASET = "ucsc_gnomad_pext_hg38"
UCSC_BIGWIG_SOURCE = "ucsc_hg38_gnomad_pext_bigwig"


class HandlerError(Exception):
    pass


def import_file(
    db,
    path: str | Path,
    *,
    gene_filter: Optional[str] = None,
    dataset: str = DEFAULT_DATASET,
    variant_ids: Optional[list[int]] = None,
) -> dict:
    """Import pext rows matched by coordinate.

    Supports either one-value-per-row (`tissue`, `pext`) or one-tissue-per-column
    exports. Rows without ref/alt are matched to all variants at the base.
    """
    p = Path(path)
    if not p.exists():
        raise HandlerError(f"pext file not found: {p}")

    index = load_variant_index(db, gene_filter=gene_filter, variant_ids=variant_ids)
    counts = {"rows": 0, "matched_rows": 0, "annotations": 0}

    for row in read_dicts(p):
        counts["rows"] += 1
        key = parse_variant_key(row)
        if not key:
            continue
        ids = match_variant_ids(index, key)
        if not ids:
            continue
        counts["matched_rows"] += 1
        values = _extract_pext_values(row)
        transcript_id = first_present(row, "transcript", "transcript_id", "Transcript", "transcript_id_version")
        for variant_id in ids:
            for tissue, score in values:
                db.upsert_expression(
                    variant_id,
                    "pext",
                    dataset=dataset,
                    tissue=tissue,
                    transcript_id=str(transcript_id or ""),
                    score=score,
                    source=SOURCE,
                )
                counts["annotations"] += 1
    return counts


def run_batch(
    db,
    *,
    file_path: str,
    limit: Optional[int] = None,
    gene_filter: Optional[str] = None,
    dataset: str = DEFAULT_DATASET,
) -> dict:
    """Drain queued pext jobs using a local pext TSV/CSV."""
    jobs = db.claim_pending_jobs(source=SOURCE, limit=limit if limit is not None else 1_000_000)
    if not jobs:
        return {"claimed": 0, "annotated": 0, "failed": 0, "rows": 0}
    job_ids_by_variant = {j["variant_id"]: j["id"] for j in jobs}
    counts = import_file(
        db,
        file_path,
        gene_filter=gene_filter,
        dataset=dataset,
        variant_ids=list(job_ids_by_variant),
    )

    annotated_variants = {
        row["variant_id"]
        for row in db.conn.execute(
            f"""
            SELECT DISTINCT variant_id
            FROM annotations_expression
            WHERE metric = 'pext' AND variant_id IN ({','.join('?' * len(job_ids_by_variant))})
            """,
            list(job_ids_by_variant),
        )
    }
    failed = 0
    for variant_id, job_id in job_ids_by_variant.items():
        if variant_id in annotated_variants:
            db.mark_job_done(job_id)
        else:
            db.mark_job_failed(job_id, "no matching pext row found")
            failed += 1

    return {
        "claimed": len(jobs),
        "annotated": len(annotated_variants),
        "failed": failed,
        "rows": counts["rows"],
    }


def import_bigwig_dir(
    db,
    path: str | Path,
    *,
    gene_filter: Optional[str] = None,
    dataset: str = UCSC_BIGWIG_DATASET,
    source_label: str = UCSC_BIGWIG_SOURCE,
    bigwig_average_bin: Optional[str] = None,
) -> dict:
    """Import UCSC/gnomAD pext bigWig tracks at each variant position.

    The repo's `data/pext/ucsc_hg38` cache is one bigWig per tissue. For
    single-base SNVs, `bigWigAverageOverBed`'s mean0 column is the pext value,
    with uncovered positions represented as 0.0.
    """
    p = Path(path)
    if not p.exists() or not p.is_dir():
        raise HandlerError(f"pext bigWig directory not found: {p}")
    bw_files = sorted(p.glob("*.bw"))
    if not bw_files:
        raise HandlerError(f"No .bw files found in pext directory: {p}")

    binary = _resolve_bigwig_average(bigwig_average_bin)
    variants = _load_variants_for_bigwig(db, gene_filter=gene_filter)
    if not variants:
        return {"variants": 0, "tissues": len(bw_files), "annotations": 0}

    annotations = 0
    with tempfile.TemporaryDirectory(prefix="pext_bigwig_") as tmp:
        bed_path = Path(tmp) / "variants.bed"
        _write_variant_bed(variants, bed_path)

        for bw in bw_files:
            tissue = bw.stem
            out_path = Path(tmp) / f"{tissue}.tab"
            subprocess.run(
                [binary, str(bw), str(bed_path), str(out_path)],
                check=True,
                capture_output=True,
                text=True,
            )
            rows = []
            with open(out_path, "r") as fh:
                for line in fh:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 5:
                        continue
                    try:
                        variant_id = int(parts[0])
                        score = float(parts[4])
                    except ValueError:
                        continue
                    rows.append((variant_id, "pext", dataset, tissue, "", score, source_label))
            annotations += _bulk_upsert_expression(db, rows)

    return {"variants": len(variants), "tissues": len(bw_files), "annotations": annotations}


def _extract_pext_values(row: dict) -> list[tuple[str, float]]:
    # Long form: one row per tissue, score in a pext-ish column.
    tissue = first_present(row, "tissue", "Tissue")
    pext_value = first_present(row, "pext", "mean_pext", "mean", "score", "value")
    score = coerce_float(pext_value)
    if tissue and score is not None:
        return [(str(tissue), score)]

    # Wide form: all numeric non-coordinate columns are treated as tissue tracks.
    out: list[tuple[str, float]] = []
    for key, value in row.items():
        if key is None or key.lower() in COORD_COLUMNS:
            continue
        score = coerce_float(value)
        if score is not None:
            out.append((key, score))
    return out


def _resolve_bigwig_average(explicit: Optional[str]) -> str:
    candidates = [
        explicit,
        os.environ.get("BIGWIG_AVERAGE_OVER_BED"),
        shutil.which("bigWigAverageOverBed"),
        str(Path.home() / "tools" / "bin" / "bigWigAverageOverBed"),
    ]
    for cand in candidates:
        if cand and Path(cand).exists():
            return str(cand)
    raise HandlerError(
        "bigWigAverageOverBed not found. Set BIGWIG_AVERAGE_OVER_BED or pass --bigwig-average-bin."
    )


def _load_variants_for_bigwig(db, *, gene_filter: Optional[str]) -> list[dict]:
    params: list = []
    sql = """
        SELECT DISTINCT v.id, v.chromosome, v.position
        FROM variants v
    """
    if gene_filter:
        sql += " JOIN variant_consequences c ON c.variant_id = v.id WHERE c.gene_symbol = ?"
        params.append(gene_filter.upper())
    cur = db.conn.execute(sql, params)
    return [dict(r) for r in cur.fetchall()]


def _write_variant_bed(variants: list[dict], path: Path) -> None:
    with open(path, "w") as f:
        for v in variants:
            chrom = str(v["chromosome"])
            if not chrom.startswith("chr"):
                chrom = f"chr{chrom}"
            pos = int(v["position"])
            f.write(f"{chrom}\t{pos - 1}\t{pos}\t{int(v['id'])}\n")


def _bulk_upsert_expression(db, rows: list[tuple]) -> int:
    if not rows:
        return 0
    sql = """
        INSERT INTO annotations_expression
            (variant_id, metric, dataset, tissue, transcript_id, score, source)
        VALUES (?, ?, ?, ?, ?, ?, ?)
        ON CONFLICT(variant_id, metric, dataset, tissue, transcript_id) DO UPDATE SET
            score = excluded.score,
            source = excluded.source,
            fetched_at = CURRENT_TIMESTAMP
    """
    db.conn.executemany(sql, rows)
    db.conn.commit()
    return len(rows)
