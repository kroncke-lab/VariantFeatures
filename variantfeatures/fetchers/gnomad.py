"""Fetch gnomAD allele frequencies.

gnomAD v4 provides population allele frequencies for genetic variants.
Uses the GraphQL API at gnomad.broadinstitute.org.

API docs: https://gnomad.broadinstitute.org/api
"""

import time
from typing import Iterator, Optional, Dict
import requests

from ..handlers import gnomad as gnomad_handler

# gnomAD GraphQL API endpoint
GNOMAD_API = "https://gnomad.broadinstitute.org/api"

# Rate limiting
RATE_LIMIT_DELAY = 0.5  # seconds between requests


def _make_graphql_request(query: str) -> Dict:
    """Make a GraphQL request to gnomAD API with rate limiting."""
    response = requests.post(
        GNOMAD_API,
        json={"query": query},
        headers={"Content-Type": "application/json"},
        timeout=120
    )
    response.raise_for_status()
    
    data = response.json()
    if "errors" in data:
        raise RuntimeError(f"gnomAD API error: {data['errors']}")
    
    time.sleep(RATE_LIMIT_DELAY)
    return data


def fetch_gnomad(
    gene: str,
    dataset: str = "gnomad_r4",
    reference_genome: str = "GRCh38"
) -> Iterator[dict]:
    """
    Fetch gnomAD frequencies for variants in a gene.
    
    Args:
        gene: Gene symbol (e.g., "KCNH2")
        dataset: gnomAD dataset version ("gnomad_r4", "gnomad_r3", etc.)
        reference_genome: "GRCh38" or "GRCh37"
        
    Yields:
        Dicts with:
            - variant_id: gnomAD variant ID (chr-pos-ref-alt)
            - hgvs_p: protein change (e.g., "p.Arg528His") or None
            - hgvs_c: coding change or None
            - gnomad_af: global allele frequency
            - gnomad_af_popmax: max population AF
            - gnomad_homozygotes: count of homozygotes
            - gnomad_an: allele number (coverage)
    """
    # Build query with inline values (gnomAD API prefers this)
    query = f'''
    {{
      gene(gene_symbol: "{gene.upper()}", reference_genome: {reference_genome}) {{
        gene_id
        symbol
        variants(dataset: {dataset}) {{
          variant_id
          hgvsc
          hgvsp
          consequence
          exome {{
            ac
            an
            af
            homozygote_count
          }}
          genome {{
            ac
            an
            af
            homozygote_count
          }}
        }}
      }}
    }}
    '''
    
    print(f"Fetching gnomAD variants for {gene}...")
    
    try:
        data = _make_graphql_request(query)
    except Exception as e:
        print(f"  Error: {e}")
        return
    
    gene_data = data.get("data", {}).get("gene")
    if not gene_data:
        print(f"  No data found for gene: {gene}")
        return
    
    variants = gene_data.get("variants", [])
    print(f"  Found {len(variants)} variants")
    
    for var in variants:
        # Combine exome and genome data (prefer exome if available)
        exome = var.get("exome") or {}
        genome = var.get("genome") or {}
        
        # Get AF (prefer exome)
        af = exome.get("af") or genome.get("af")
        an = exome.get("an") or genome.get("an")
        homozygotes = (exome.get("homozygote_count") or 0) + (genome.get("homozygote_count") or 0)
        
        yield {
            "variant_id": var.get("variant_id"),
            "hgvs_p": var.get("hgvsp"),
            "hgvs_c": var.get("hgvsc"),
            "consequence": var.get("consequence"),
            "gnomad_af": af,
            "gnomad_homozygotes": homozygotes,
            "gnomad_an": an,
        }


def fetch_single_variant(
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    dataset: str = "gnomad_r4",
    reference_genome: str = "GRCh38"
) -> Optional[dict]:
    """
    Fetch gnomAD data for a single variant by coordinates.

    This compatibility wrapper delegates to the normalized gnomAD handler, whose
    GraphQL query matches the current public schema:
    `variant(variantId: ..., dataset: ...)`. The single-variant endpoint no
    longer accepts a `referenceGenome` argument and does not expose `hgvsc` /
    `hgvsp` directly on `VariantDetails`; transcript consequences live under
    `transcript_consequences`.
    
    Args:
        chrom: Chromosome (e.g., "7" or "chr7")
        pos: Position (1-based)
        ref: Reference allele
        alt: Alternate allele
        dataset: gnomAD dataset version
        reference_genome: "GRCh38" or "GRCh37"
        
    Returns:
        Dict with variant data, or None if not found
    """
    if reference_genome.upper() not in {"GRCH38", "GRCH37"}:
        raise ValueError("reference_genome must be 'GRCh38' or 'GRCh37'")

    try:
        var = gnomad_handler.fetch_gnomad(chrom, pos, ref, alt, dataset=dataset)
    except Exception as e:
        variant_id = f"{chrom.replace('chr', '')}-{pos}-{ref}-{alt}"
        print(f"Error fetching {variant_id}: {e}")
        return None

    if not var:
        return None

    exome = var.get("exome") or {}
    genome = var.get("genome") or {}
    exome_af = exome.get("af")
    genome_af = genome.get("af")

    return {
        "variant_id": var.get("variant_id"),
        "rsids": var.get("rsids") or [],
        "hgvs_p": None,
        "hgvs_c": None,
        "gnomad_af": exome_af if exome_af is not None else genome_af,
        "gnomad_homozygotes": (exome.get("homozygote_count") or 0) + (genome.get("homozygote_count") or 0),
        "gnomad_an": exome.get("an") if exome.get("an") is not None else genome.get("an"),
    }


if __name__ == "__main__":
    # Test with KCNH2
    import sys
    import itertools
    gene = sys.argv[1] if len(sys.argv) > 1 else "KCNH2"
    
    print(f"Testing gnomAD fetcher for {gene}")
    
    variants = list(itertools.islice(fetch_gnomad(gene), 10))
    print(f"\nSample variants ({len(variants)} shown):")
    for v in variants:
        af_str = f"{v['gnomad_af']:.2e}" if v['gnomad_af'] else "N/A"
        print(f"  {v['variant_id']}: AF={af_str}, hgvs_p={v['hgvs_p']}")
