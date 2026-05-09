"""Annotation handlers — one per data source.

Each handler exports a callable `handle(db, variant_id, payload)` that performs
exactly one annotation job. The worker (`variantfeatures.worker`) dispatches
queue rows to the correct handler by `source` name.
"""

from . import clingen_ar, myvariant, gnomad, alphamissense, annovar, vep  # noqa: F401

# Per-job handlers: one variant per call. Used by the generic worker loop.
# Each module exposes SOURCE, DEFAULT_RATE_LIMIT_SEC, handle(db, variant_id, payload).
HANDLERS: dict = {
    clingen_ar.SOURCE: clingen_ar,
    myvariant.SOURCE: myvariant,
    gnomad.SOURCE: gnomad,
}

# Batch handlers: one tool invocation processes many variants. The worker loop
# can't dispatch these one row at a time, so they're surfaced as separate CLI
# subcommands and skipped by `annotate-pending`. Each module exposes
# SOURCE and run_batch(db, **kwargs) -> dict.
BATCH_HANDLERS: dict = {
    alphamissense.SOURCE: alphamissense,
    annovar.SOURCE: annovar,
    vep.SOURCE: vep,
}
