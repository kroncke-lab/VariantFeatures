"""Normalized gene build orchestration."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import requests

from .database import VariantDB
from .enumerate import populate_for_transcript
from .transcripts import TranscriptError, fetch_gene_transcripts


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
    "clinvar",
    "gene_constraint",
    "nmd_rule",
    "pext_bigwig",
    "frameshift_proxy",
}

SOURCE_ALIASES = {
    "core": {
        "myvariant",
        "clinvar",
        "gnomad",
        "alphamissense",
        "revel",
        "alphafold",
        "gene_constraint",
        "nmd_rule",
        "pext_bigwig",
        "frameshift_proxy",
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
    "clinical": {"myvariant", "clinvar"},
    "clinvar": {"myvariant", "clinvar"},
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
    "lof": {"vep", "nmd_rule", "frameshift_proxy"},
    "myvariant": {"myvariant"},
    "clingen_ar": {"clingen_ar"},
    "annovar": {"annovar"},
    "vep": {"vep"},
    "frameshift": {"frameshift_proxy"},
    "frameshift_proxy": {"frameshift_proxy"},
}

MISSENSE_ONLY = {"alphamissense", "revel"}
PROTEIN_POSITION_SOURCES = {"alphafold"}
ALL_CODING_SNV_SOURCES = {"clingen_ar", "myvariant", "gnomad", "annovar", "vep", "cadd"}
NONSENSE_REQUIRED_SOURCES = {"nmd_rule", "frameshift_proxy"}


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
    transcript_labels: list[str]
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
    frameshift_max_steps: int = 20,
    isoforms: str = "canonical",
    max_isoforms: Optional[int] = None,
) -> GeneBuildResult:
    """Build/update a normalized database for one gene."""
    selected = parse_sources(sources)
    symbol = gene.upper()

    try:
        transcripts = fetch_gene_transcripts(
            symbol,
            include=isoforms,
            max_transcripts=max_isoforms,
        )
    except TranscriptError:
        raise
    except requests.RequestException as e:
        raise BuildError(f"Network error talking to Ensembl REST: {e}") from e

    db = VariantDB(db_path)
    primary_transcript = transcripts[0]
    transcript_labels = [
        transcript.refseq_match or transcript.transcript_id_versioned
        for transcript in transcripts
    ]
    transcript_label = transcript_labels[0]
    db.upsert_gene(
        symbol,
        ensembl_id=primary_transcript.gene_ensembl,
        canonical_transcript=transcript_label,
    )

    enumeration = _empty_enumeration_summary()
    for transcript in transcripts:
        db.upsert_transcript(transcript)
        summary = populate_for_transcript(
            transcript,
            db,
            types=_types_with_required_nonsense(types, selected),
            sources_by_consequence=sources_by_consequence(selected),
            enqueue=True,
            limit=enumerate_limit,
        )
        _merge_enumeration_summary(
            enumeration,
            summary,
            transcript_label=transcript.refseq_match or transcript.transcript_id_versioned,
        )

    result = GeneBuildResult(
        gene=symbol,
        transcript_label=transcript_label,
        transcript_labels=transcript_labels,
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
                frameshift_max_steps=frameshift_max_steps,
                transcripts=transcripts,
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
    frameshift_max_steps: int,
    transcripts: Optional[list] = None,
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

    if "clinvar" in selected:
        def _clinvar():
            from .handlers import clinvar

            return clinvar.import_gene(db, gene)

        record("clinvar", _clinvar)

    # MyVariant has an efficient batch endpoint and covers CADD/dbNSFP/ClinVar/
    # many population fallbacks; run it before narrower sources.
    if "myvariant" in selected:
        def _myvariant():
            from .handlers import myvariant

            return myvariant.run_batch(db, limit=annotation_limit, gene_filter=gene)

        record("myvariant", _myvariant)

    for source in ("clingen_ar", "gnomad", "cadd"):
        if source not in selected:
            continue

        def _worker(src=source):
            from .worker import run_pending

            return run_pending(db, source=src, max_jobs=annotation_limit, gene_filter=gene)

        record(source, _worker)

    if "alphamissense" in selected:
        def _alphamissense():
            from .handlers import alphamissense

            return alphamissense.run_batch(db, limit=annotation_limit, gene_filter=gene)

        record("alphamissense", _alphamissense)

    if "revel" in selected:
        def _revel():
            from .handlers import revel

            return revel.run_batch(db, limit=annotation_limit, gene_filter=gene)

        record("revel", _revel)

    if "alphafold" in selected:
        def _alphafold():
            from .handlers import alphafold

            return alphafold.run_batch(db, limit=annotation_limit, gene_filter=gene)

        record("alphafold", _alphafold)

    if "annovar" in selected:
        def _annovar():
            from .handlers import annovar

            return annovar.run_batch(db, limit=annotation_limit, gene_filter=gene)

        record("annovar", _annovar)

    if "vep" in selected:
        def _vep():
            from .handlers import vep

            return vep.run_batch(db, limit=annotation_limit, gene_filter=gene)

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

    if "frameshift_proxy" in selected:
        def _frameshift_proxy():
            from . import frameshift

            selected_transcripts = transcripts or [None]
            total = {
                "positions": 0,
                "mapped_positions": 0,
                "unmapped_positions": 0,
                "mappings": 0,
                "nonsense_variants": 0,
                "jobs_queued": 0,
                "transcripts": len(selected_transcripts),
            }
            for transcript in selected_transcripts:
                summary = frameshift.annotate_gene(
                    db,
                    gene,
                    transcript=transcript,
                    max_steps=frameshift_max_steps,
                )
                for key in (
                    "positions",
                    "mapped_positions",
                    "unmapped_positions",
                    "mappings",
                    "nonsense_variants",
                    "jobs_queued",
                ):
                    total[key] += summary.get(key, 0)
            return total

        record("frameshift_proxy", _frameshift_proxy)

    return runs


def _types_with_required_nonsense(types: str | list[str] | tuple[str, ...], selected: set[str]):
    if not (selected & NONSENSE_REQUIRED_SOURCES):
        return types
    if isinstance(types, str):
        parsed = [t.strip() for t in types.split(",") if t.strip()]
    else:
        parsed = list(types)
    lowered = {t.strip().lower() for t in parsed}
    if "nonsense" not in lowered and "stop_gained" not in lowered:
        parsed.append("nonsense")
    return parsed


def _empty_enumeration_summary() -> dict:
    return {
        "variants": 0,
        "consequences": 0,
        "jobs_queued": 0,
        "by_consequence": {},
        "by_transcript": {},
    }


def _merge_enumeration_summary(total: dict, summary: dict, *, transcript_label: str) -> None:
    total["variants"] += summary.get("variants", 0)
    total["consequences"] += summary.get("consequences", 0)
    total["jobs_queued"] += summary.get("jobs_queued", 0)
    for consequence, count in summary.get("by_consequence", {}).items():
        total["by_consequence"][consequence] = total["by_consequence"].get(consequence, 0) + count
    total["by_transcript"][transcript_label] = summary
