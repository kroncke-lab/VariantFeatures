"""Importer for external NMD predictors such as NMDEP / NMDetective."""

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


SOURCE = "nmd_external"


class HandlerError(Exception):
    pass


DEFAULT_SCORE_COLUMNS = (
    "nmdep_score",
    "nmd_efficiency",
    "nmd_probability",
    "NMDetectiveA",
    "NMDetectiveB",
    "score",
    "probability",
)


def import_file(
    db,
    path: str | Path,
    *,
    predictor: str,
    gene_filter: Optional[str] = None,
    score_column: Optional[str] = None,
    category_column: Optional[str] = None,
) -> dict:
    """Import coordinate-keyed NMD predictor output into annotations_pathogenicity."""
    p = Path(path)
    if not p.exists():
        raise HandlerError(f"NMD predictor file not found: {p}")
    if not predictor:
        raise HandlerError("predictor is required, e.g. 'nmdep' or 'nmdetective'")

    index = load_variant_index(db, gene_filter=gene_filter)
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
        scores = _extract_scores(row, predictor=predictor, score_column=score_column)
        category = first_present(row, category_column) if category_column else first_present(
            row, "prediction", "class", "category", "nmd_class", "label"
        )
        for variant_id in ids:
            for pred_name, score in scores:
                db.upsert_pathogenicity(
                    variant_id,
                    pred_name,
                    score=score,
                    category=str(category) if category is not None else None,
                    source=SOURCE,
                )
                counts["annotations"] += 1

    return counts


def _extract_scores(row: dict, *, predictor: str, score_column: Optional[str]) -> list[tuple[str, float]]:
    if score_column:
        score = coerce_float(first_present(row, score_column))
        return [(predictor, score)] if score is not None else []

    out: list[tuple[str, float]] = []
    lower_to_key = {str(k).lower(): k for k in row}
    for candidate in DEFAULT_SCORE_COLUMNS:
        key = lower_to_key.get(candidate.lower())
        if not key:
            continue
        score = coerce_float(row.get(key))
        if score is None:
            continue
        name = predictor
        if candidate.lower() not in {"score", "probability", f"{predictor.lower()}_score"}:
            name = f"{predictor}_{candidate.lower()}"
        out.append((name, score))
    return out
