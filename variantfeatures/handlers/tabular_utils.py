"""Shared helpers for local tabular annotation imports."""

from __future__ import annotations

import csv
import gzip
import re
from pathlib import Path
from typing import Iterable, Optional

from .clingen_ar import normalize_chromosome


COORD_COLUMNS = {
    "variant", "chrom", "chr", "#chrom", "chromosome", "contig",
    "pos", "position", "start", "ref", "reference", "alt", "alternate",
    "gene", "gene_id", "gene_symbol", "symbol", "tissue", "transcript",
    "transcript_id", "dataset",
}


def open_text(path: str | Path):
    p = Path(path)
    if str(p).endswith(".gz"):
        return gzip.open(p, "rt", newline="")
    return open(p, "r", newline="")


def read_dicts(path: str | Path) -> Iterable[dict]:
    """Read CSV/TSV rows with delimiter inferred from the header line."""
    with open_text(path) as fh:
        first = fh.readline()
        if not first:
            return
        delimiter = "\t" if "\t" in first else ","
        fh.seek(0)
        yield from csv.DictReader(fh, delimiter=delimiter)


def parse_variant_key(row: dict) -> Optional[tuple[str, int, Optional[str], Optional[str]]]:
    """Parse a row into (chrom, pos, ref_or_None, alt_or_None)."""
    variant = first_present(row, "variant", "Variant", "variant_id", "var")
    if variant:
        parsed = parse_variant_string(str(variant))
        if parsed:
            return parsed

    chrom = first_present(row, "chrom", "chr", "#chrom", "#CHROM", "chromosome", "contig")
    pos = first_present(row, "pos", "position", "start", "POS", "Position")
    if not chrom or pos is None:
        return None
    try:
        pos_i = int(float(str(pos)))
    except (TypeError, ValueError):
        return None
    ref = first_present(row, "ref", "reference", "REF")
    alt = first_present(row, "alt", "alternate", "ALT")
    return normalize_chromosome(str(chrom)), pos_i, _clean_allele(ref), _clean_allele(alt)


def parse_variant_string(value: str) -> Optional[tuple[str, int, Optional[str], Optional[str]]]:
    """Parse common variant strings, e.g. `chr7:140753336:A>T`."""
    text = value.strip()
    patterns = [
        r"^(?:chr)?([A-Za-z0-9]+):(\d+):([ACGTN]+)>([ACGTN]+)$",
        r"^(?:chr)?([A-Za-z0-9]+):(\d+):([ACGTN]+):([ACGTN]+)$",
        r"^(?:chr)?([A-Za-z0-9]+)-(\d+)-([ACGTN]+)-([ACGTN]+)$",
    ]
    for pat in patterns:
        m = re.match(pat, text, flags=re.IGNORECASE)
        if not m:
            continue
        chrom, pos, ref, alt = m.groups()
        return normalize_chromosome(chrom), int(pos), ref.upper(), alt.upper()
    return None


def load_variant_index(
    db,
    *,
    gene_filter: Optional[str] = None,
    variant_ids: Optional[Iterable[int]] = None,
) -> dict:
    """Index variants by exact allele and by position."""
    params: list = []
    where = []
    join = ""
    if gene_filter:
        join = "JOIN variant_consequences c ON c.variant_id = v.id"
        where.append("c.gene_symbol = ?")
        params.append(gene_filter.upper())
    ids = list(variant_ids or [])
    if ids:
        where.append(f"v.id IN ({','.join('?' * len(ids))})")
        params.extend(ids)
    sql = f"""
        SELECT DISTINCT v.id, v.chromosome, v.position, v.ref, v.alt
        FROM variants v
        {join}
    """
    if where:
        sql += " WHERE " + " AND ".join(where)
    cur = db.conn.execute(sql, params)

    by_exact: dict[tuple, list[int]] = {}
    by_pos: dict[tuple, list[int]] = {}
    for row in cur.fetchall():
        chrom = normalize_chromosome(row["chromosome"])
        pos = int(row["position"])
        ref = row["ref"]
        alt = row["alt"]
        by_exact.setdefault((chrom, pos, ref, alt), []).append(row["id"])
        by_pos.setdefault((chrom, pos), []).append(row["id"])
    return {"by_exact": by_exact, "by_pos": by_pos}


def match_variant_ids(index: dict, key: tuple[str, int, Optional[str], Optional[str]]) -> list[int]:
    chrom, pos, ref, alt = key
    if ref and alt:
        return list(index["by_exact"].get((chrom, pos, ref, alt), []))
    return list(index["by_pos"].get((chrom, pos), []))


def first_present(row: dict, *keys: str):
    lower_map = {str(k).lower(): v for k, v in row.items()}
    for key in keys:
        if key in row and row[key] not in (None, "", ".", "NA", "nan"):
            return row[key]
        val = lower_map.get(key.lower())
        if val not in (None, "", ".", "NA", "nan"):
            return val
    return None


def coerce_float(value) -> Optional[float]:
    if value in (None, "", ".", "NA", "nan", "NaN", "-"):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _clean_allele(value) -> Optional[str]:
    if value in (None, "", ".", "NA", "-"):
        return None
    return str(value).upper()
