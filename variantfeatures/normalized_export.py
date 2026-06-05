"""Normalized exports from the canonical variant + annotation tables.

The normalized schema deliberately keeps source families in separate tables:
pathogenicity, population, clinical, splice, structure, expression, and
conservation. This module turns that normalized storage into user-facing export
formats without hiding where a feature came from.
"""

from __future__ import annotations

import csv
import json
import re
from collections import defaultdict
from pathlib import Path
from typing import Iterable, Optional


BASE_COLUMNS = [
    "gene",
    "variant_id",
    "chromosome",
    "position",
    "ref",
    "alt",
    "variant_type",
    "ca_id",
    "vrs_id",
    "hgvs_g",
    "transcript_id",
    "consequence",
    "hgvs_c",
    "hgvs_p",
    "aa_pos",
    "aa_ref",
    "aa_alt",
    "codon_pos",
    "codon_ref",
    "codon_alt",
    "consequence_source",
    "is_canonical",
    "is_mane_select",
    "is_mane_plus_clinical",
    "aliases.rsid",
    "aliases.clinvar_vcv",
    "aliases.gnomad_id",
]

FEATURE_GROUPS = (
    "pathogenicity",
    "population",
    "clinical",
    "splice",
    "expression",
    "structure",
    "conservation",
    "gene_constraint",
)

LONG_COLUMNS = BASE_COLUMNS[:-3] + [
    "feature_group",
    "feature_name",
    "version",
    "dataset",
    "population",
    "tissue",
    "feature_transcript_id",
    "score_type",
    "field",
    "value",
    "source",
]

_SAFE_COMPONENT_RE = re.compile(r"[^A-Za-z0-9_-]+")


class ExportError(Exception):
    pass


def parse_groups(groups: Optional[str | Iterable[str]]) -> set[str]:
    """Parse a comma-separated group list, accepting ``all`` as a shortcut."""
    if groups is None:
        return set(FEATURE_GROUPS)
    if isinstance(groups, str):
        raw = [g.strip() for g in groups.split(",") if g.strip()]
    else:
        raw = [str(g).strip() for g in groups if str(g).strip()]
    if not raw or "all" in {g.lower() for g in raw}:
        return set(FEATURE_GROUPS)
    parsed = {g.lower() for g in raw}
    unknown = parsed - set(FEATURE_GROUPS)
    if unknown:
        raise ExportError(
            f"Unknown feature group(s): {', '.join(sorted(unknown))}. "
            f"Valid groups: {', '.join(FEATURE_GROUPS)}"
        )
    return parsed


def export_gene(
    db,
    gene: str,
    output_path: str | Path,
    *,
    layout: str = "wide",
    groups: Optional[str | Iterable[str]] = None,
    include_provenance: bool = False,
) -> dict:
    """Write a normalized gene export.

    ``layout='wide'`` writes one row per canonical variant with namespaced feature
    columns. ``layout='long'`` writes one row per feature field and is useful for
    auditing where feature families live.
    """
    selected_groups = parse_groups(groups)
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    if layout == "wide":
        rows, fieldnames = build_wide_rows(
            db,
            gene,
            groups=selected_groups,
            include_provenance=include_provenance,
        )
        _write_csv(output, rows, fieldnames)
        return {"layout": "wide", "variants": len(rows), "columns": len(fieldnames), "path": output}

    if layout == "long":
        rows = list(iter_long_rows(db, gene, groups=selected_groups))
        _write_csv(output, rows, LONG_COLUMNS)
        return {"layout": "long", "rows": len(rows), "columns": len(LONG_COLUMNS), "path": output}

    raise ExportError(f"Unknown export layout: {layout!r}")


def build_wide_rows(
    db,
    gene: str,
    *,
    groups: Optional[set[str]] = None,
    include_provenance: bool = False,
) -> tuple[list[dict], list[str]]:
    """Build one-row-per-variant export rows and field names."""
    selected = groups or set(FEATURE_GROUPS)
    rows = _base_variant_rows(db, gene)
    if not rows:
        return [], list(BASE_COLUMNS)

    fieldnames = list(BASE_COLUMNS)
    by_variant = {row["variant_id"]: row for row in rows}
    _add_alias_columns(db, gene, by_variant)

    if "gene_constraint" in selected:
        _add_gene_constraint(db, gene, rows, fieldnames, include_provenance=include_provenance)
    if "pathogenicity" in selected:
        _add_pathogenicity(db, gene, by_variant, fieldnames, include_provenance=include_provenance)
    if "population" in selected:
        _add_population(db, gene, by_variant, fieldnames, include_provenance=include_provenance)
    if "clinical" in selected:
        _add_clinical(db, gene, by_variant, fieldnames)
    if "splice" in selected:
        _add_splice(db, gene, by_variant, fieldnames, include_provenance=include_provenance)
    if "expression" in selected:
        _add_expression(db, gene, by_variant, fieldnames, include_provenance=include_provenance)
    if "structure" in selected:
        _add_structure(db, gene, by_variant, fieldnames, include_provenance=include_provenance)
    if "conservation" in selected:
        _add_conservation(db, gene, by_variant, fieldnames, include_provenance=include_provenance)

    return rows, fieldnames


def iter_long_rows(db, gene: str, *, groups: Optional[set[str]] = None) -> Iterable[dict]:
    """Yield normalized feature rows in long form."""
    selected = groups or set(FEATURE_GROUPS)
    base_rows = _base_variant_rows(db, gene)
    if not base_rows:
        return
    base_by_variant = {row["variant_id"]: {k: row.get(k) for k in BASE_COLUMNS[:-3]} for row in base_rows}

    if "gene_constraint" in selected:
        for row in _query_gene_constraint(db, gene):
            for field in (
                "pli",
                "lof_z",
                "mis_z",
                "syn_z",
                "oe_lof",
                "oe_lof_lower",
                "oe_lof_upper",
                "oe_mis",
                "oe_mis_upper",
                "exp_lof",
                "obs_lof",
            ):
                value = row.get(field)
                if value is None:
                    continue
                for base in base_by_variant.values():
                    yield _long_row(
                        base,
                        feature_group="gene_constraint",
                        feature_name=field,
                        dataset=row["dataset"],
                        field="value",
                        value=value,
                        source=row.get("source"),
                    )

    if "pathogenicity" in selected:
        for row in _query_gene_annotation(db, gene, "annotations_pathogenicity p", "p.*"):
            base = base_by_variant[row["variant_id"]]
            feature = row["predictor"]
            version = row["predictor_version"]
            for field in ("score", "rank_score", "category"):
                if row[field] is not None:
                    yield _long_row(
                        base,
                        feature_group="pathogenicity",
                        feature_name=feature,
                        version=version,
                        field=field,
                        value=row[field],
                        source=row.get("source"),
                    )

    if "population" in selected:
        for row in _query_gene_annotation(db, gene, "annotations_population p", "p.*"):
            base = base_by_variant[row["variant_id"]]
            for field in ("af", "ac", "an", "n_homozygotes", "filter_status"):
                if row[field] is not None:
                    yield _long_row(
                        base,
                        feature_group="population",
                        feature_name=field,
                        dataset=row["dataset"],
                        population=row["pop"],
                        field=field,
                        value=row[field],
                        source=row.get("source"),
                    )

    if "clinical" in selected:
        for row in _query_gene_annotation(db, gene, "annotations_clinical c", "c.*"):
            base = base_by_variant[row["variant_id"]]
            for field in ("classification", "review_status", "stars", "last_evaluated", "conditions"):
                if row[field] is not None:
                    yield _long_row(
                        base,
                        feature_group="clinical",
                        feature_name=row["source"],
                        dataset=row["record_id"],
                        field=field,
                        value=row[field],
                        source=row["source"],
                    )

    if "splice" in selected:
        for row in _query_gene_annotation(db, gene, "annotations_splice s", "s.*"):
            base = base_by_variant[row["variant_id"]]
            for field in ("score", "distance"):
                if row[field] is not None:
                    yield _long_row(
                        base,
                        feature_group="splice",
                        feature_name=row["predictor"],
                        version=row["predictor_version"],
                        score_type=row["score_type"],
                        field=field,
                        value=row[field],
                        source=row.get("source"),
                    )

    if "expression" in selected:
        for row in _query_gene_annotation(db, gene, "annotations_expression e", "e.*"):
            base = base_by_variant[row["variant_id"]]
            if row["score"] is not None:
                yield _long_row(
                    base,
                    feature_group="expression",
                    feature_name=row["metric"],
                    dataset=row["dataset"],
                    tissue=row["tissue"],
                    feature_transcript_id=row["transcript_id"],
                    field="score",
                    value=row["score"],
                    source=row.get("source"),
                )

    if "structure" in selected:
        for row in _query_gene_annotation(db, gene, "annotations_structure s", "s.*"):
            base = base_by_variant[row["variant_id"]]
            for field in ("score", "category", "residue_number", "protein_accession"):
                if row[field] is not None:
                    yield _long_row(
                        base,
                        feature_group="structure",
                        feature_name=row["feature"],
                        version=row["feature_version"],
                        field=field,
                        value=row[field],
                        source=row.get("source"),
                    )

    if "conservation" in selected:
        for row in _query_gene_annotation(db, gene, "annotations_conservation c", "c.*"):
            base = base_by_variant[row["variant_id"]]
            for field in ("score", "rank_score"):
                if row[field] is not None:
                    yield _long_row(
                        base,
                        feature_group="conservation",
                        feature_name=row["metric"],
                        field=field,
                        value=row[field],
                        source=row.get("source"),
                    )


def _base_variant_rows(db, gene: str) -> list[dict]:
    sql = """
    WITH ranked AS (
        SELECT
            ? AS gene,
            v.id AS variant_id,
            v.chromosome,
            v.position,
            v.ref,
            v.alt,
            v.variant_type,
            v.ca_id,
            v.vrs_id,
            v.hgvs_g,
            c.transcript_id,
            c.consequence,
            c.hgvs_c,
            c.hgvs_p,
            c.aa_pos,
            c.aa_ref,
            c.aa_alt,
            c.codon_pos,
            c.codon_ref,
            c.codon_alt,
            c.source AS consequence_source,
            c.is_canonical,
            c.is_mane_select,
            c.is_mane_plus_clinical,
            ROW_NUMBER() OVER (
                PARTITION BY v.id
                ORDER BY
                    CASE c.source
                        WHEN 'enumerated' THEN 0
                        WHEN 'vep' THEN 1
                        WHEN 'annovar' THEN 2
                        ELSE 3
                    END,
                    c.is_mane_select DESC,
                    c.is_canonical DESC,
                    c.transcript_id
            ) AS rn
        FROM variants v
        JOIN variant_consequences c ON c.variant_id = v.id
        WHERE UPPER(c.gene_symbol) = ?
    )
    SELECT
        gene, variant_id, chromosome, position, ref, alt, variant_type,
        ca_id, vrs_id, hgvs_g, transcript_id, consequence, hgvs_c, hgvs_p,
        aa_pos, aa_ref, aa_alt, codon_pos, codon_ref, codon_alt,
        consequence_source, is_canonical, is_mane_select, is_mane_plus_clinical,
        '' AS "aliases.rsid",
        '' AS "aliases.clinvar_vcv",
        '' AS "aliases.gnomad_id"
    FROM ranked
    WHERE rn = 1
    ORDER BY
        CASE
            WHEN chromosome IN ('X', 'x') THEN 23
            WHEN chromosome IN ('Y', 'y') THEN 24
            WHEN chromosome IN ('M', 'MT', 'm', 'mt') THEN 25
            ELSE CAST(chromosome AS INTEGER)
        END,
        position,
        ref,
        alt
    """
    cur = db.conn.execute(sql, [gene.upper(), gene.upper()])
    return [dict(row) for row in cur.fetchall()]


def _add_alias_columns(db, gene: str, by_variant: dict[int, dict]) -> None:
    sql = """
    WITH gene_vars AS (
        SELECT DISTINCT variant_id
        FROM variant_consequences
        WHERE UPPER(gene_symbol) = ?
    )
    SELECT a.variant_id, a.alias_type, a.alias_value
    FROM variant_aliases a
    JOIN gene_vars g ON g.variant_id = a.variant_id
    WHERE a.alias_type IN ('rsid', 'clinvar_vcv', 'gnomad_id')
    ORDER BY a.alias_type, a.alias_value
    """
    grouped: dict[tuple[int, str], list[str]] = defaultdict(list)
    for row in db.conn.execute(sql, [gene.upper()]):
        grouped[(row["variant_id"], row["alias_type"])].append(row["alias_value"])
    for (variant_id, alias_type), values in grouped.items():
        col = f"aliases.{alias_type}"
        row = by_variant.get(variant_id)
        if row is not None:
            row[col] = "|".join(_unique(values))


def _add_gene_constraint(
    db,
    gene: str,
    rows: list[dict],
    fieldnames: list[str],
    *,
    include_provenance: bool,
) -> None:
    for constraint in _query_gene_constraint(db, gene):
        dataset = _safe_component(constraint["dataset"])
        prefix = f"gene_constraint.{dataset}"
        for field in (
            "pli",
            "lof_z",
            "mis_z",
            "syn_z",
            "oe_lof",
            "oe_lof_lower",
            "oe_lof_upper",
            "oe_mis",
            "oe_mis_upper",
            "exp_lof",
            "obs_lof",
        ):
            _set_all(rows, fieldnames, f"{prefix}.{field}", constraint.get(field))
        if include_provenance:
            _set_all(rows, fieldnames, f"{prefix}.source", constraint.get("source"))


def _add_pathogenicity(
    db,
    gene: str,
    by_variant: dict[int, dict],
    fieldnames: list[str],
    *,
    include_provenance: bool,
) -> None:
    for row in _query_gene_annotation(db, gene, "annotations_pathogenicity p", "p.*"):
        prefix = _versioned_prefix("pathogenicity", row["predictor"], row["predictor_version"])
        dest = by_variant.get(row["variant_id"])
        if dest is None:
            continue
        _set(dest, fieldnames, f"{prefix}.score", row["score"])
        _set(dest, fieldnames, f"{prefix}.rank_score", row["rank_score"])
        _set(dest, fieldnames, f"{prefix}.category", row["category"])
        if include_provenance:
            _set(dest, fieldnames, f"{prefix}.source", row["source"])


def _add_population(
    db,
    gene: str,
    by_variant: dict[int, dict],
    fieldnames: list[str],
    *,
    include_provenance: bool,
) -> None:
    for row in _query_gene_annotation(db, gene, "annotations_population p", "p.*"):
        prefix = f"population.{_safe_component(row['dataset'])}.{_safe_component(row['pop'])}"
        dest = by_variant.get(row["variant_id"])
        if dest is None:
            continue
        _set(dest, fieldnames, f"{prefix}.af", row["af"])
        _set(dest, fieldnames, f"{prefix}.ac", row["ac"])
        _set(dest, fieldnames, f"{prefix}.an", row["an"])
        _set(dest, fieldnames, f"{prefix}.homozygotes", row["n_homozygotes"])
        _set(dest, fieldnames, f"{prefix}.filter_status", row["filter_status"])
        if include_provenance:
            _set(dest, fieldnames, f"{prefix}.source", row["source"])


def _add_clinical(db, gene: str, by_variant: dict[int, dict], fieldnames: list[str]) -> None:
    grouped: dict[tuple[int, str], list[dict]] = defaultdict(list)
    for row in _query_gene_annotation(db, gene, "annotations_clinical c", "c.*"):
        grouped[(row["variant_id"], row["source"])].append(dict(row))

    for (variant_id, source), rows in grouped.items():
        dest = by_variant.get(variant_id)
        if dest is None:
            continue
        rows = sorted(
            rows,
            key=lambda r: (
                -(r.get("stars") or 0),
                -_classification_severity(r.get("classification")),
                r.get("record_id") or "",
            ),
        )
        prefix = f"clinical.{_safe_component(source)}"
        _set(dest, fieldnames, f"{prefix}.classification", "|".join(_unique(r.get("classification") for r in rows)))
        _set(dest, fieldnames, f"{prefix}.review_status", "|".join(_unique(r.get("review_status") for r in rows)))
        _set(dest, fieldnames, f"{prefix}.stars_max", max((r.get("stars") or 0 for r in rows), default=None))
        _set(dest, fieldnames, f"{prefix}.last_evaluated_max", max((r.get("last_evaluated") or "" for r in rows), default=""))
        _set(dest, fieldnames, f"{prefix}.record_ids", "|".join(_unique(r.get("record_id") for r in rows)))
        _set(dest, fieldnames, f"{prefix}.conditions", "|".join(_unique(r.get("conditions") for r in rows)))


def _add_splice(
    db,
    gene: str,
    by_variant: dict[int, dict],
    fieldnames: list[str],
    *,
    include_provenance: bool,
) -> None:
    for row in _query_gene_annotation(db, gene, "annotations_splice s", "s.*"):
        prefix = _versioned_prefix("splice", row["predictor"], row["predictor_version"])
        prefix = f"{prefix}.{_safe_component(row['score_type'])}"
        dest = by_variant.get(row["variant_id"])
        if dest is None:
            continue
        _set(dest, fieldnames, f"{prefix}.score", row["score"])
        _set(dest, fieldnames, f"{prefix}.distance", row["distance"])
        if include_provenance:
            _set(dest, fieldnames, f"{prefix}.source", row["source"])


def _add_expression(
    db,
    gene: str,
    by_variant: dict[int, dict],
    fieldnames: list[str],
    *,
    include_provenance: bool,
) -> None:
    for row in _query_gene_annotation(db, gene, "annotations_expression e", "e.*"):
        components = [
            "expression",
            row["metric"],
            row["dataset"] or "dataset_unspecified",
            row["tissue"] or "tissue_unspecified",
        ]
        if row["transcript_id"]:
            components.append(row["transcript_id"])
        prefix = ".".join(_safe_component(c) for c in components)
        dest = by_variant.get(row["variant_id"])
        if dest is None:
            continue
        _set(dest, fieldnames, f"{prefix}.score", row["score"])
        if include_provenance:
            _set(dest, fieldnames, f"{prefix}.source", row["source"])


def _add_structure(
    db,
    gene: str,
    by_variant: dict[int, dict],
    fieldnames: list[str],
    *,
    include_provenance: bool,
) -> None:
    for row in _query_gene_annotation(db, gene, "annotations_structure s", "s.*"):
        prefix = _versioned_prefix("structure", row["feature"], row["feature_version"])
        dest = by_variant.get(row["variant_id"])
        if dest is None:
            continue
        _set(dest, fieldnames, f"{prefix}.score", row["score"])
        _set(dest, fieldnames, f"{prefix}.category", row["category"])
        _set(dest, fieldnames, f"{prefix}.residue_number", row["residue_number"])
        _set(dest, fieldnames, f"{prefix}.protein_accession", row["protein_accession"])
        if include_provenance:
            _set(dest, fieldnames, f"{prefix}.source", row["source"])


def _add_conservation(
    db,
    gene: str,
    by_variant: dict[int, dict],
    fieldnames: list[str],
    *,
    include_provenance: bool,
) -> None:
    for row in _query_gene_annotation(db, gene, "annotations_conservation c", "c.*"):
        prefix = f"conservation.{_safe_component(row['metric'])}"
        dest = by_variant.get(row["variant_id"])
        if dest is None:
            continue
        _set(dest, fieldnames, f"{prefix}.score", row["score"])
        _set(dest, fieldnames, f"{prefix}.rank_score", row["rank_score"])
        if include_provenance:
            _set(dest, fieldnames, f"{prefix}.source", row["source"])


def _query_gene_annotation(db, gene: str, table_expr: str, select_expr: str):
    sql = f"""
    WITH gene_vars AS (
        SELECT DISTINCT variant_id
        FROM variant_consequences
        WHERE UPPER(gene_symbol) = ?
    )
    SELECT {select_expr}
    FROM {table_expr}
    JOIN gene_vars g ON g.variant_id = {table_expr.split()[1]}.variant_id
    ORDER BY {table_expr.split()[1]}.variant_id
    """
    cur = db.conn.execute(sql, [gene.upper()])
    for row in cur:
        yield dict(row)


def _query_gene_constraint(db, gene: str) -> list[dict]:
    cur = db.conn.execute(
        """
        SELECT *
        FROM gene_constraint
        WHERE UPPER(gene_symbol) = ?
        ORDER BY dataset
        """,
        [gene.upper()],
    )
    return [dict(row) for row in cur.fetchall()]


def _long_row(base: dict, **fields) -> dict:
    row = {col: "" for col in LONG_COLUMNS}
    row.update(base)
    row.update(fields)
    return row


def _set(row: dict, fieldnames: list[str], key: str, value) -> None:
    if value is None:
        return
    if key not in fieldnames:
        fieldnames.append(key)
    existing = row.get(key)
    if existing in (None, ""):
        row[key] = value
    elif existing != value:
        row[key] = "|".join(_unique([existing, value]))


def _set_all(rows: list[dict], fieldnames: list[str], key: str, value) -> None:
    if value is None:
        return
    if key not in fieldnames:
        fieldnames.append(key)
    for row in rows:
        row[key] = value


def _safe_component(value) -> str:
    text = str(value or "").strip()
    if not text:
        return "unspecified"
    text = _SAFE_COMPONENT_RE.sub("_", text)
    text = re.sub(r"_+", "_", text).strip("_")
    return text or "unspecified"


def _versioned_prefix(group: str, name: str, version: Optional[str]) -> str:
    prefix = f"{group}.{_safe_component(name)}"
    if version:
        prefix += f".version_{_safe_component(version)}"
    return prefix


def _unique(values: Iterable) -> list[str]:
    seen = set()
    out = []
    for value in values:
        if value in (None, ""):
            continue
        text = str(value)
        if text in seen:
            continue
        seen.add(text)
        out.append(text)
    return out


def _classification_severity(value: Optional[str]) -> int:
    if not value:
        return 0
    text = value.lower()
    if "pathogenic" in text and "conflict" not in text:
        return 5 if "likely" not in text else 4
    if "uncertain" in text or "vus" in text:
        return 3
    if "benign" in text:
        return 2 if "likely" in text else 1
    return 0


def _write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def feature_schema() -> dict[str, dict]:
    """Return the public map of feature groups to storage tables and column names."""
    return {
        "pathogenicity": {
            "table": "annotations_pathogenicity",
            "examples": ["alphamissense", "revel", "cadd_phred", "dbnsfp predictors", "nmd_rule"],
            "wide_prefix": "pathogenicity.<predictor>[.version_<version>].<score|rank_score|category>",
        },
        "population": {
            "table": "annotations_population",
            "examples": ["gnomAD exome/genome AF, AC, AN, homozygotes"],
            "wide_prefix": "population.<dataset>.<population>.<af|ac|an|homozygotes|filter_status>",
        },
        "clinical": {
            "table": "annotations_clinical",
            "examples": ["ClinVar classification, review status, stars"],
            "wide_prefix": "clinical.<source>.<classification|review_status|stars_max|record_ids>",
        },
        "splice": {
            "table": "annotations_splice",
            "examples": ["SpliceAI", "dbscSNV", "MaxEntScan", "AbSplice"],
            "wide_prefix": "splice.<predictor>[.version_<version>].<score_type>.<score|distance>",
        },
        "expression": {
            "table": "annotations_expression",
            "examples": ["gnomAD pext", "AbExp", "splice-site expression"],
            "wide_prefix": "expression.<metric>.<dataset>.<tissue>[.<transcript>].score",
        },
        "structure": {
            "table": "annotations_structure",
            "examples": ["AlphaFold pLDDT", "InterPro domains"],
            "wide_prefix": "structure.<feature>[.version_<version>].<score|category|residue_number>",
        },
        "conservation": {
            "table": "annotations_conservation",
            "examples": ["phyloP", "phastCons", "GERP++", "SiPhy"],
            "wide_prefix": "conservation.<metric>.<score|rank_score>",
        },
        "gene_constraint": {
            "table": "gene_constraint",
            "examples": ["gnomAD pLI, LOEUF, missense Z"],
            "wide_prefix": "gene_constraint.<dataset>.<metric>",
        },
    }


def feature_schema_json() -> str:
    return json.dumps(feature_schema(), indent=2, sort_keys=True)
