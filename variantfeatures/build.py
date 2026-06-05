"""Normalized gene build orchestration."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import requests

from .database import VariantDB
from .enumerate import populate_for_transcript
from .transcripts import TranscriptError, fetch_canonical_transcript


JOB_SOURCES = {
    "clingen_ar",
    "myvariant",
    "gnomad",
    "alphamissense",
    "revel",
    "alphafold",
    "annovar",
    "vep",
    "cadd",
}

DIRECT_SOURCES = {
    "gene_constraint",
    "nmd_rule",
    "pext_bigwig",
}

SOURCE_ALIASES = {
    "core": {
        "myvariant",
        "gnomad",
        "alphamissense",
        "revel",
        "alphafold",
        "gene_constraint",
        "nmd_rule",
        "pext_bigwig",
    },
    "all": JOB_SOURCES | DIRECT_SOURCES,
    "identity": {"clingen_ar"},
    "aliases": {"clingen_ar"},
    "pathogenicity": {"myvariant", "alphamissense", "revel"},
    "dbnsfp": {"myvariant"},
    "revel": {"revel"},
    "alphamissense": {"alphamissense"},
    "cadd": {"cadd"},
    "population": {"gnomad"},
    "gnomad": {"gnomad"},
    "clinical": {"myvariant"},
    "clinvar": {"myvariant"},
    "structure": {"alphafold"},
    "alphafold": {"alphafold"},
    "splice": {"vep", "annovar"},
    "expression": {"pext_bigwig"},
    "pext": {"pext_bigwig"},
    "pext_bigwig": {"pext_bigwig"},
    "constraint": {"gene_constraint"},
    "gene_constraint": {"gene_constraint"},
    "nmd": {"nmd_rule", "vep"},
    "nmd_rule": {"nmd_rule"},
    "lof": {"vep", "nmd_rule"},
    "myvariant": {"myvariant"},
    "clingen_ar": {"clingen_ar"},
    "annovar": {"annovar"},
    "vep": {"vep"},
}

MISSENSE_ONLY = {"alphamissense", "revel"}
PROTEIN_POSITION_SOURCES = {"alphafold"}
ALL_CODING_SNV_SOURCES = {"clingen_ar", "myvariant", "gnomad", "annovar", "vep", "cadd"}


class BuildError(Exception):
    pass


@dataclass
class SourceRun:
    source: str
    status: str
    summary: dict = field(default_factory=dict)
    error: Optional[str] = None


@dataclass
class GeneBuildResult:
    gene: str
    transcript_label: str
    db_path: Path
    enumeration: dict
    sources: list[SourceRun] = field(default_factory=list)


def parse_sources(sources: str | None) -> set[str]:
    """Expand source aliases into canonical normalized build source names."""
    raw = [s.strip().lower() for s in (sources or "core").split(",") if s.strip()]
    if not raw:
        raw = ["core"]
    expanded: set[str] = set()
    unknown = []
    for item in raw:
        aliases = SOURCE_ALIASES.get(item)
        if aliases is None:
            unknown.append(item)
            continue
        expanded |= set(aliases)
    if unknown:
        valid = ", ".join(sorted(SOURCE_ALIASES))
        raise BuildError(f"Unknown build source(s): {', '.join(unknown)}. Valid sources: {valid}")
    return expanded


def sources_by_consequence(selected_sources: set[str]) -> dict[str, tuple[str, ...]]:
    """Return the annotation queue map for the selected sources."""
    job_sources = selected_sources & JOB_SOURCES
    by_consequence: dict[str, list[str]] = {}

    for consequence in (
        "missense_variant",
        "stop_gained",
        "stop_lost",
        "start_lost",
        "synonymous_variant",
    ):
        sources: list[str] = []
        for source in sorted(job_sources):
            if source in MISSENSE_ONLY and consequence != "missense_variant":
                continue
            if source in PROTEIN_POSITION_SOURCES and consequence == "synonymous_variant":
                continue
            if source in ALL_CODING_SNV_SOURCES or source in MISSENSE_ONLY or source in PROTEIN_POSITION_SOURCES:
                sources.append(source)
        by_consequence[consequence] = tuple(sources)

    return by_consequence


def build_gene(
    gene: str,
    *,
    db_path: Optional[Path] = None,
    sources: str | None = None,
    types: str = "missense,nonsense,stop_lost,start_lost,synonymous",
    enumerate_limit: Optional[int] = None,
    run_annotations: bool = True,
    annotation_limit: Optional[int] = None,
    strict: bool = False,
    pext_bigwig_dir: str | Path = Path("data") / "pext" / "ucsc_hg38",
    pext_dataset: str = "ucsc_gnomad_pext_hg38",
) -> GeneBuildResult:
    """Build/update a normalized database for one gene."""
    selected = parse_sources(sources)
    symbol = gene.upper()

    try:
        transcript = fetch_canonical_transcript(symbol)
    except TranscriptError:
        raise
    except requests.RequestException as e:
        raise BuildError(f"Network error talking to Ensembl REST: {e}") from e

    db = VariantDB(db_path)
    transcript_label = transcript.refseq_match or transcript.transcript_id_versioned
    db.upsert_gene(
        symbol,
        ensembl_id=transcript.gene_ensembl,
        canonical_transcript=transcript_label,
    )

    enumeration = populate_for_transcript(
        transcript,
        db,
        types=types,
        sources_by_consequence=sources_by_consequence(selected),
        enqueue=True,
        limit=enumerate_limit,
    )

    result = GeneBuildResult(
        gene=symbol,
        transcript_label=transcript_label,
        db_path=db.db_path,
        enumeration=enumeration,
    )

    if run_annotations:
        result.sources.extend(
            run_selected_sources(
                db,
                symbol,
                selected,
                annotation_limit=annotation_limit,
                strict=strict,
                pext_bigwig_dir=pext_bigwig_dir,
                pext_dataset=pext_dataset,
            )
        )

    return result


def run_selected_sources(
    db,
    gene: str,
    selected: set[str],
    *,
    annotation_limit: Optional[int],
    strict: bool,
    pext_bigwig_dir: str | Path,
    pext_dataset: str,
) -> list[SourceRun]:
    """Drain selected annotation sources in dependency-aware order."""
    runs: list[SourceRun] = []

    def record(source: str, func):
        try:
            summary = func()
            runs.append(SourceRun(source=source, status="done", summary=summary))
        except Exception as e:  # noqa: BLE001 - build should report all unavailable sources.
            if strict:
                raise
            runs.append(SourceRun(source=source, status="skipped", error=str(e)))

    if "gene_constraint" in selected:
        def _constraint():
            from .handlers import gnomad_constraint

            ok = gnomad_constraint.annotate_gene(db, gene)
            return {"annotated": int(bool(ok))}

        record("gene_constraint", _constraint)

    # MyVariant has an efficient batch endpoint and covers CADD/dbNSFP/ClinVar/
    # many population fallbacks; run it before narrower sources.
    if "myvariant" in selected:
        def _myvariant():
            from .handlers import myvariant

            return myvariant.run_batch(db, limit=annotation_limit)

        record("myvariant", _myvariant)

    for source in ("clingen_ar", "gnomad", "cadd"):
        if source not in selected:
            continue

        def _worker(src=source):
            from .worker import run_pending

            return run_pending(db, source=src, max_jobs=annotation_limit)

        record(source, _worker)

    if "alphamissense" in selected:
        def _alphamissense():
            from .handlers import alphamissense

            return alphamissense.run_batch(db, limit=annotation_limit)

        record("alphamissense", _alphamissense)

    if "revel" in selected:
        def _revel():
            from .handlers import revel

            return revel.run_batch(db, limit=annotation_limit)

        record("revel", _revel)

    if "alphafold" in selected:
        def _alphafold():
            from .handlers import alphafold

            return alphafold.run_batch(db, limit=annotation_limit)

        record("alphafold", _alphafold)

    if "annovar" in selected:
        def _annovar():
            from .handlers import annovar

            return annovar.run_batch(db, limit=annotation_limit)

        record("annovar", _annovar)

    if "vep" in selected:
        def _vep():
            from .handlers import vep

            return vep.run_batch(db, limit=annotation_limit)

        record("vep", _vep)

    if "pext_bigwig" in selected:
        def _pext_bigwig():
            from .handlers import pext

            return pext.import_bigwig_dir(
                db,
                pext_bigwig_dir,
                gene_filter=gene,
                dataset=pext_dataset,
            )

        record("pext_bigwig", _pext_bigwig)

    if "nmd_rule" in selected:
        def _nmd_rule():
            from .handlers import nmd_rules

            return nmd_rules.annotate_gene(db, gene)

        record("nmd_rule", _nmd_rule)

    return runs
