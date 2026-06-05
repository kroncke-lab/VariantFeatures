"""AbSplice / AbExp local-output importer."""

from __future__ import annotations

from pathlib import Path
from typing import Optional

from .tabular_utils import (
    coerce_float,
    first_present,
    load_variant_index,
    match_variant_ids,
    parse_variant_key,
    read_dicts,
)


SOURCE = "absplice"
DEFAULT_DATASET = "absplice"


class HandlerError(Exception):
    pass


ABSPLICE_COLUMNS = {
    "absplice_dna": "absplice_dna",
    "absplice_rna": "absplice_rna",
    "absplice": "absplice",
}


def import_file(
    db,
    path: str | Path,
    *,
    gene_filter: Optional[str] = None,
    dataset: str = DEFAULT_DATASET,
    source_label: str = SOURCE,
) -> dict:
    """Import AbSplice tabular output.

    The canonical AbSplice unique key is `(variant, gene_id, tissue)`, with
    scores such as `AbSplice_DNA`. AbExp-like numeric columns are stored in
    `annotations_expression`.
    """
    p = Path(path)
    if not p.exists():
        raise HandlerError(f"AbSplice file not found: {p}")

    index = load_variant_index(db, gene_filter=gene_filter)
    counts = {"rows": 0, "matched_rows": 0, "splice": 0, "expression": 0}

    for row in read_dicts(p):
        counts["rows"] += 1
        key = parse_variant_key(row)
        if not key:
            continue
        ids = match_variant_ids(index, key)
        if not ids:
            continue
        counts["matched_rows"] += 1
        tissue = str(first_present(row, "tissue", "Tissue") or "overall")
        transcript_id = str(
            first_present(row, "transcript", "transcript_id", "Transcript", "gene_id") or ""
        )

        for variant_id in ids:
            for predictor, score_type, score in _extract_absplice_scores(row, default_score_type=tissue):
                db.upsert_splice(
                    variant_id,
                    predictor,
                    score_type=score_type,
                    score=score,
                    source=source_label,
                )
                counts["splice"] += 1
            for metric, metric_tissue, score in _extract_abexp_scores(row, default_tissue=tissue):
                db.upsert_expression(
                    variant_id,
                    metric,
                    dataset=dataset,
                    tissue=metric_tissue,
                    transcript_id=transcript_id,
                    score=score,
                    source=source_label,
                )
                counts["expression"] += 1

    return counts


def _extract_absplice_scores(row: dict, *, default_score_type: str) -> list[tuple[str, str, float]]:
    out: list[tuple[str, str, float]] = []
    lower_to_key = {str(k).lower(): k for k in row}
    for lower_name, predictor in ABSPLICE_COLUMNS.items():
        key = lower_to_key.get(lower_name)
        if key:
            score = coerce_float(row.get(key))
            if score is not None:
                out.append((predictor, default_score_type, score))

    for key, value in row.items():
        if key is None:
            continue
        lower = str(key).lower()
        if lower in ABSPLICE_COLUMNS:
            continue
        for prefix, predictor in (
            ("absplice_dna_", "absplice_dna"),
            ("absplice_rna_", "absplice_rna"),
            ("absplice_", "absplice"),
        ):
            if lower.startswith(prefix) and lower != prefix.rstrip("_"):
                score = coerce_float(value)
                if score is not None:
                    out.append((predictor, str(key)[len(prefix):], score))
                break
        else:
            if lower == "delta_score":
                score = coerce_float(value)
                if score is not None:
                    out.append(("spliceai_delta_score", "overall", score))
    return out


def _extract_abexp_scores(row: dict, *, default_tissue: str) -> list[tuple[str, str, float]]:
    out: list[tuple[str, str, float]] = []
    for key, value in row.items():
        if key is None:
            continue
        lower = str(key).lower()
        metric = None
        tissue = default_tissue
        if lower.startswith("delta_logit_psi_"):
            metric = "absplice_delta_logit_psi"
            tissue = str(key)[len("delta_logit_psi_"):]
        elif lower.startswith("delta_psi_"):
            metric = "absplice_delta_psi"
            tissue = str(key)[len("delta_psi_"):]
        elif lower.startswith("splice_site_is_expressed_"):
            metric = "splice_site_is_expressed"
            tissue = str(key)[len("splice_site_is_expressed_"):]
        elif lower.startswith("abexp"):
            metric = lower

        if metric:
            score = coerce_float(value)
            if score is not None:
                out.append((metric, tissue, score))
    return out
