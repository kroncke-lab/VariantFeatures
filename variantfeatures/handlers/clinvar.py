"""Normalized ClinVar importer backed by local variant_summary.txt.gz."""

from __future__ import annotations

from pathlib import Path
from typing import Optional

from ..fetchers.clinvar import DEFAULT_DATA_DIR, fetch_clinvar
from .tabular_utils import load_variant_index, match_variant_ids, normalize_chromosome


SOURCE = "clinvar"
DEFAULT_DATA_FILE = DEFAULT_DATA_DIR / "variant_summary.txt.gz"


class HandlerError(Exception):
    pass


def file_present(file_path: Optional[str] = None) -> bool:
    return _resolve_data_file(file_path).exists()


def import_gene(
    db,
    gene_symbol: str,
    *,
    data_file: Optional[str | Path] = None,
    assembly: str = "GRCh38",
    include_all_types: bool = False,
) -> dict:
    """Import ClinVar clinical assertions for one built gene.

    Rows are matched to canonical variants by GRCh38 chromosome, position, ref,
    and alt, then written to ``annotations_clinical`` with ``source='clinvar'``.
    """
    path = _resolve_data_file(data_file)
    if not path.exists():
        raise HandlerError(
            f"ClinVar variant_summary file not found: {path}. "
            "Download NCBI ClinVar variant_summary.txt.gz into data/."
        )

    gene = gene_symbol.upper()
    index = load_variant_index(db, gene_filter=gene)
    counts = {"rows": 0, "matched_rows": 0, "annotations": 0, "aliases": 0}

    for record in fetch_clinvar(
        gene,
        data_file=path,
        assembly=assembly,
        include_all_types=include_all_types,
    ):
        counts["rows"] += 1
        key = _record_key(record)
        if key is None:
            continue
        variant_ids = match_variant_ids(index, key)
        if not variant_ids:
            continue
        counts["matched_rows"] += 1
        for variant_id in variant_ids:
            db.upsert_clinical(
                variant_id,
                SOURCE,
                record_id=_vcv_id(record.get("clinvar_id")),
                classification=record.get("clinvar_significance"),
                review_status=record.get("clinvar_review_status"),
                stars=record.get("clinvar_stars"),
                last_evaluated=record.get("clinvar_last_evaluated"),
                conditions=None,
            )
            counts["annotations"] += 1
            counts["aliases"] += db.add_aliases(
                variant_id,
                _aliases(record),
                source=SOURCE,
            )
    return counts


def _resolve_data_file(data_file: Optional[str | Path]) -> Path:
    return Path(data_file) if data_file else DEFAULT_DATA_FILE


def _record_key(record: dict) -> Optional[tuple[str, int, str, str]]:
    chrom = record.get("chromosome")
    pos = record.get("position")
    ref = record.get("ref")
    alt = record.get("alt")
    if not (chrom and pos and ref and alt):
        return None
    return normalize_chromosome(str(chrom)), int(pos), str(ref).upper(), str(alt).upper()


def _vcv_id(value) -> str:
    if value in (None, ""):
        return ""
    try:
        return f"VCV{int(value):09d}"
    except (TypeError, ValueError):
        return str(value)


def _aliases(record: dict) -> list[dict]:
    aliases = []
    vcv = _vcv_id(record.get("clinvar_id"))
    if vcv:
        aliases.append({"alias_type": "clinvar_vcv", "alias_value": vcv})
    allele_id = record.get("clinvar_allele_id")
    if allele_id not in (None, ""):
        aliases.append({"alias_type": "clinvar_allele", "alias_value": str(allele_id)})
    return aliases
