"""ClinGen Allele Registry annotation handler.

Takes one queued job, looks up the variant by GRCh38 (chrom, pos, ref, alt),
builds an NC-style g.HGVS, and resolves it through the ClinGen Allele Registry
to populate the canonical CA-ID and the full alias inventory (rsID, ClinVar
VCV / Allele ID, gnomAD ID, MyVariant.info ID, c./p. HGVS, etc.).

Idempotent: re-running over a variant whose aliases are already present is a
no-op (alias inserts use INSERT OR IGNORE).
"""

from __future__ import annotations

from typing import Optional

from .. import identity


SOURCE = "clingen_ar"
DEFAULT_RATE_LIMIT_SEC = 0.20  # ClinGen AR is fine at ~5 req/s


# NCBI RefSeq accessions for the GRCh38 primary assembly.
GRCH38_NC_ACCESSIONS: dict[str, str] = {
    "1": "NC_000001.11",
    "2": "NC_000002.12",
    "3": "NC_000003.12",
    "4": "NC_000004.12",
    "5": "NC_000005.10",
    "6": "NC_000006.12",
    "7": "NC_000007.14",
    "8": "NC_000008.11",
    "9": "NC_000009.12",
    "10": "NC_000010.11",
    "11": "NC_000011.10",
    "12": "NC_000012.12",
    "13": "NC_000013.11",
    "14": "NC_000014.9",
    "15": "NC_000015.10",
    "16": "NC_000016.10",
    "17": "NC_000017.11",
    "18": "NC_000018.10",
    "19": "NC_000019.10",
    "20": "NC_000020.11",
    "21": "NC_000021.9",
    "22": "NC_000022.11",
    "X": "NC_000023.11",
    "Y": "NC_000024.10",
    "MT": "NC_012920.1",
    "M": "NC_012920.1",
}


class HandlerError(Exception):
    """Raised when this handler can't run for a variant (caller marks job failed)."""


def normalize_chromosome(chrom: str) -> str:
    """Strip 'chr' prefix and normalize 'M' to 'MT' for the accession lookup."""
    c = chrom.upper().lstrip("CHR") if chrom else ""
    return c


def hgvs_g_for_snv(chromosome: str, position: int, ref: str, alt: str) -> str:
    """Build NC_xxxxx.x:g.PPPR>A HGVS for a GRCh38 SNV."""
    chrom = normalize_chromosome(chromosome)
    accession = GRCH38_NC_ACCESSIONS.get(chrom)
    if accession is None:
        raise HandlerError(f"No GRCh38 NC accession for chromosome {chromosome!r}")
    if len(ref) != 1 or len(alt) != 1:
        raise HandlerError(
            f"clingen_ar handler currently only handles SNVs (got ref={ref!r}, alt={alt!r})"
        )
    return f"{accession}:g.{position}{ref}>{alt}"


def handle(db, variant_id: int, payload: Optional[str] = None, *, timeout: int = 30) -> None:
    """Resolve one variant via ClinGen Allele Registry and persist aliases."""
    cur = db.conn.execute(
        "SELECT chromosome, position, ref, alt FROM variants WHERE id = ?",
        [variant_id],
    )
    row = cur.fetchone()
    if row is None:
        raise HandlerError(f"variant_id {variant_id} not found in variants table")

    hgvs_g = hgvs_g_for_snv(row["chromosome"], row["position"], row["ref"], row["alt"])

    rv = identity.resolve(hgvs_g, kind="hgvs_g", timeout=timeout)
    identity.persist(rv, db)
