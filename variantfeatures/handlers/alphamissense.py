"""AlphaMissense local-file handler.

Reads the published `AlphaMissense_aa_substitutions.tsv.gz` (~1.2GB) which is
keyed by `(uniprot_id, protein_variant)` rather than genomic coordinates. We
match pending jobs by joining `annotation_jobs` -> `variants` ->
`variant_consequences` to recover `(gene_symbol, aa_ref, aa_pos, aa_alt)`,
then convert to the AlphaMissense short form (e.g. `KCNH2 + A614V`).

For per-gene UniProt resolution we use a small hardcoded map for the cardiac
channelopathy Phase 1 genes (sufficient for the R01 work), with the option
to extend via the `gene_uniprot` argument or `ALPHAMISSENSE_UNIPROT` env var
(comma-separated `GENE:UNIPROT` pairs).

File format:
    uniprot_id  protein_variant  am_pathogenicity  am_class

Configure file path via `ALPHAMISSENSE_FILE` env var or `file_path=` kwarg.
Default: `data/alphamissense/AlphaMissense_aa_substitutions.tsv.gz`.

Download: https://storage.googleapis.com/dm_alphamissense/AlphaMissense_aa_substitutions.tsv.gz
"""

from __future__ import annotations

import gzip
import os
from pathlib import Path
from typing import Optional


SOURCE = "alphamissense"

DEFAULT_REL_PATH = Path("data") / "alphamissense" / "AlphaMissense_aa_substitutions.tsv.gz"
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent


# Cardiac channelopathy genes covered out of the box. Extend as needed.
_BUILTIN_GENE_TO_UNIPROT: dict[str, str] = {
    "KCNH2": "Q12809",
    "KCNQ1": "P51787",
    "SCN5A": "Q14524",
    "RYR2": "Q92736",
    "CACNA1C": "Q13936",
    "KCNJ2": "P63252",
    "KCNE1": "P15382",
    "KCNE2": "Q9Y6J6",
    # Common ACMG SF v3.2 genes (incomplete; expand as the project moves to ACMG-81)
    "BRCA1": "P38398",
    "BRCA2": "P51587",
    "TP53": "P04637",
    "MYH7": "P12883",
    "MYBPC3": "Q14896",
    "TTN": "Q8WZ42",
    "LDLR": "P01130",
    "APOB": "P04114",
    "PCSK9": "Q8NBP7",
    "BRAF": "P15056",  # for the test fixture
    # ACMG SF v3.2 + amyloidosis / metabolic
    "TTR": "P02766",     # transthyretin (familial amyloid polyneuropathy)
    "ALPL": "P05186",    # tissue-nonspecific alkaline phosphatase (hypophosphatasia)
    "RB1": "P06400",
    "VHL": "P40337",
    "ATP7B": "P35670",
    "RYR1": "P21817",
    "CACNA1S": "Q13698",
    "SDHB": "P21912",
    "SDHC": "Q99643",
    "SDHD": "O14521",
    "MEN1": "O00255",
    "RET": "P07949",
    "PTEN": "P60484",
    "STK11": "Q15831",
    "MUTYH": "Q9UIF7",
    "MLH1": "P40692",
    "MSH2": "P43246",
    "MSH6": "P52701",
    "PMS2": "P54278",
    "APC": "P25054",
}


class HandlerError(Exception):
    pass


# ---------------------------------------------------------------------------
# Path / config helpers
# ---------------------------------------------------------------------------

def resolve_file_path(explicit: Optional[str] = None) -> Path:
    if explicit:
        return Path(explicit)
    env_path = os.environ.get("ALPHAMISSENSE_FILE")
    if env_path:
        return Path(env_path)
    return PROJECT_ROOT / DEFAULT_REL_PATH


def file_present(file_path: Optional[str] = None) -> bool:
    return resolve_file_path(file_path).exists()


def gene_uniprot_map(extra: Optional[dict[str, str]] = None) -> dict[str, str]:
    """Return GENE -> UniProt accession map. Pulls in env-configured overrides."""
    out = dict(_BUILTIN_GENE_TO_UNIPROT)
    env_pairs = os.environ.get("ALPHAMISSENSE_UNIPROT")
    if env_pairs:
        for pair in env_pairs.split(","):
            if ":" in pair:
                g, u = pair.split(":", 1)
                out[g.strip().upper()] = u.strip()
    if extra:
        for g, u in extra.items():
            out[g.upper()] = u
    return out


# ---------------------------------------------------------------------------
# Batch runner
# ---------------------------------------------------------------------------

def run_batch(
    db,
    *,
    file_path: Optional[str] = None,
    limit: Optional[int] = None,
    progress_every: int = 5_000_000,
    progress_callback=None,
    gene_uniprot: Optional[dict[str, str]] = None,
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

    pending = _load_pending_index(db, limit=limit, gene_uniprot=gene_uniprot_map(gene_uniprot))
    if not pending["by_key"]:
        # Either no pending jobs, or none of them matched a known gene -> mark them failed.
        failed = _mark_unmappable(db, pending["unmappable"])
        return {"claimed": failed, "matched": 0, "failed": failed, "lines_scanned": 0}

    matched = 0
    lines_scanned = 0

    by_key = pending["by_key"]

    with _open_maybe_gzip(path) as fh:
        header = None
        for line in fh:
            lines_scanned += 1
            if progress_callback and lines_scanned % progress_every == 0:
                progress_callback({"lines": lines_scanned, "matched": matched, "remaining": len(by_key)})

            # Skip license/copyright comment lines.
            if line.startswith("#"):
                continue
            if header is None:
                header = line.rstrip("\n").split("\t")
                # Sanity check: expect (uniprot_id, protein_variant, am_pathogenicity, am_class)
                continue

            cols = line.rstrip("\n").split("\t")
            if len(cols) < 4:
                continue
            uniprot_id, protein_variant, score_str, am_class = cols[0], cols[1], cols[2], cols[3]
            key = (uniprot_id, protein_variant)
            jobs = by_key.get(key)
            if not jobs:
                continue
            try:
                score = float(score_str)
            except ValueError:
                continue
            category = _AM_CLASS_TO_CODE.get(am_class, am_class)
            for job in jobs:
                _persist_match(db, job=job, score=score, category=category, am_class=am_class)
                matched += 1
            del by_key[key]
            if not by_key:
                break

    # Anything still in by_key was not in the AM file (extremely rare). Mark failed.
    failed = 0
    for jobs in by_key.values():
        for job in jobs:
            db.mark_job_failed(job["job_id"], "no AlphaMissense entry for this (uniprot, AA-substitution)")
            failed += 1

    # Plus jobs we couldn't even map to a uniprot to begin with (unknown gene).
    failed += _mark_unmappable(db, pending["unmappable"])

    return {"claimed": matched + failed, "matched": matched, "failed": failed, "lines_scanned": lines_scanned}


# ---------------------------------------------------------------------------
# Internals
# ---------------------------------------------------------------------------

_AM_CLASS_TO_CODE = {
    # Newer AlphaMissense releases use bare labels; older releases prefixed with `likely_`.
    "benign": "B",
    "likely_benign": "B",
    "ambiguous": "A",
    "pathogenic": "P",
    "likely_pathogenic": "P",
}


def _open_maybe_gzip(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def _load_pending_index(db, *, limit: Optional[int], gene_uniprot: dict[str, str]) -> dict:
    """Claim pending alphamissense jobs and index them by (uniprot, AM-style variant).

    Jobs whose gene_symbol isn't in the uniprot map are returned in the
    'unmappable' list so the caller can mark them failed without scanning the
    big TSV.
    """
    jobs = db.claim_pending_jobs(source=SOURCE, limit=limit if limit is not None else 1_000_000)
    if not jobs:
        return {"by_key": {}, "unmappable": []}

    ids = [j["id"] for j in jobs]
    placeholders = ",".join("?" * len(ids))
    cur = db.conn.execute(
        f"""
        SELECT j.id AS job_id, j.variant_id AS variant_id,
               c.gene_symbol AS gene_symbol,
               c.aa_ref AS aa_ref, c.aa_pos AS aa_pos, c.aa_alt AS aa_alt,
               c.consequence AS consequence
        FROM annotation_jobs j
        JOIN variant_consequences c ON c.variant_id = j.variant_id
                                       AND c.source = 'enumerated'
        WHERE j.id IN ({placeholders})
        """,
        ids,
    )

    by_key: dict[tuple, list[dict]] = {}
    unmappable: list[dict] = []
    for r in cur.fetchall():
        d = dict(r)
        consequence = d.get("consequence")
        # AM only covers missense.
        if consequence != "missense_variant":
            d["_reason"] = f"alphamissense covers missense only (consequence={consequence!r})"
            unmappable.append(d)
            continue
        gene = (d.get("gene_symbol") or "").upper()
        uniprot = gene_uniprot.get(gene)
        if not uniprot:
            d["_reason"] = f"no UniProt mapping for gene {gene!r} (extend ALPHAMISSENSE_UNIPROT)"
            unmappable.append(d)
            continue
        if not (d.get("aa_ref") and d.get("aa_pos") and d.get("aa_alt")):
            d["_reason"] = "consequence missing aa_ref / aa_pos / aa_alt"
            unmappable.append(d)
            continue
        key = (uniprot, f"{d['aa_ref']}{d['aa_pos']}{d['aa_alt']}")
        by_key.setdefault(key, []).append(d)

    return {"by_key": by_key, "unmappable": unmappable}


def _persist_match(db, *, job: dict, score: float, category: Optional[str], am_class: str) -> None:
    db.upsert_pathogenicity(
        job["variant_id"],
        "alphamissense",
        score=score,
        category=category,
        source=SOURCE,
    )
    db.mark_job_done(job["job_id"])


def _mark_unmappable(db, unmappable: list[dict]) -> int:
    for j in unmappable:
        db.mark_job_failed(j["job_id"], j.get("_reason", "unmappable"))
    return len(unmappable)
