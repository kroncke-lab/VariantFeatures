#!/usr/bin/env python3
"""
Fetch CADD scores for genetic variants.

CADD (Combined Annotation Dependent Depletion) scores predict deleteriousness
of variants based on ~60 genomic features.

Data source: https://cadd.gs.washington.edu/
API: https://cadd.gs.washington.edu/api

CADD PHRED scores:
- 10 = top 10% most deleterious
- 20 = top 1% most deleterious  
- 30 = top 0.1% most deleterious
- Suggested pathogenic threshold: 15-20
"""

import time
import requests
from typing import Iterator, Optional, Dict, List
from pathlib import Path

# CADD API endpoints
# Note: API can be unreliable; consider downloading pre-computed files for production
# Download: https://cadd.gs.washington.edu/download (per-chromosome TSV files)
CADD_API_BASE = "https://cadd.gs.washington.edu/api/v1.0"
CADD_VERSION = "GRCh38-v1.6"  # Use v1.6 for stability

# Rate limiting
RATE_LIMIT_DELAY = 0.2  # seconds between requests (be gentle with public API)


def fetch_cadd_single(
    chromosome: str,
    position: int,
    ref: str = None,
    alt: str = None,
    version: str = CADD_VERSION
) -> Optional[dict]:
    """
    Fetch CADD score for a single variant by genomic coordinates.
    
    Args:
        chromosome: Chromosome (e.g., "7" or "chr7")
        position: Genomic position (1-based)
        ref: Reference allele (optional, filters results)
        alt: Alternate allele (optional, filters results)
        version: CADD version/build (default: GRCh38-v1.7)
    
    Returns:
        Dict with cadd_raw and cadd_phred, or None if not found
    """
    # Normalize chromosome
    chrom = chromosome.replace("chr", "")
    
    # CADD API expects position range
    url = f"{CADD_API_BASE}/{version}/{chrom}:{position}-{position}"
    
    try:
        response = requests.get(url, timeout=30)
        response.raise_for_status()
        data = response.json()
        
        time.sleep(RATE_LIMIT_DELAY)
        
        if not data or len(data) < 2:
            return None
        
        # First row is header, rest are variants
        header = data[0]
        
        for row in data[1:]:
            variant = dict(zip(header, row))
            
            # Filter by ref/alt if provided
            if ref and variant.get('Ref') != ref:
                continue
            if alt and variant.get('Alt') != alt:
                continue
            
            return {
                'chromosome': variant.get('Chrom'),
                'position': int(variant.get('Pos')),
                'ref': variant.get('Ref'),
                'alt': variant.get('Alt'),
                'cadd_raw': float(variant.get('RawScore', 0)),
                'cadd_phred': float(variant.get('PHRED', 0))
            }
        
        return None
        
    except Exception as e:
        print(f"Error fetching CADD for {chrom}:{position}: {e}")
        return None


def fetch_cadd_batch(
    variants: List[Dict],
    version: str = CADD_VERSION
) -> Iterator[dict]:
    """
    Fetch CADD scores for a batch of variants.
    
    Args:
        variants: List of dicts with 'chromosome', 'position', 'ref', 'alt'
        version: CADD version/build
    
    Yields:
        Dicts with original variant info plus cadd_raw and cadd_phred
    """
    print(f"Fetching CADD scores for {len(variants)} variants...")
    
    for i, var in enumerate(variants):
        result = fetch_cadd_single(
            chromosome=var['chromosome'],
            position=var['position'],
            ref=var.get('ref'),
            alt=var.get('alt'),
            version=version
        )
        
        if result:
            # Merge with original variant info
            var.update(result)
            yield var
        else:
            # Return original with None scores
            var['cadd_raw'] = None
            var['cadd_phred'] = None
            yield var
        
        if (i + 1) % 100 == 0:
            print(f"  Processed {i + 1}/{len(variants)} variants...")


def fetch_cadd_for_gene(
    gene: str,
    variants_db_path: Optional[Path] = None
) -> Iterator[dict]:
    """
    Fetch CADD scores for all variants of a gene in the VariantFeatures database.
    
    Requires genomic coordinates to be populated in the database.
    
    Args:
        gene: Gene symbol (e.g., "KCNH2")
        variants_db_path: Path to variants.db (uses default if not provided)
    
    Yields:
        Dicts with variant info and CADD scores
    """
    import sqlite3
    
    if variants_db_path is None:
        variants_db_path = Path(__file__).parent.parent.parent / "data" / "variants.db"
    
    if not variants_db_path.exists():
        raise FileNotFoundError(f"Database not found: {variants_db_path}")
    
    conn = sqlite3.connect(variants_db_path)
    conn.row_factory = sqlite3.Row
    
    # Get variants with genomic coordinates
    cursor = conn.execute("""
        SELECT hgvs_p, chromosome, position, ref, alt
        FROM variants_missense
        WHERE gene = ? AND chromosome IS NOT NULL AND position IS NOT NULL
    """, (gene.upper(),))
    
    variants = [dict(row) for row in cursor.fetchall()]
    conn.close()
    
    if not variants:
        print(f"No variants with genomic coordinates found for {gene}")
        return
    
    print(f"Found {len(variants)} variants for {gene} with genomic coordinates")
    
    yield from fetch_cadd_batch(variants)


if __name__ == "__main__":
    import sys
    
    # Test single lookup
    print("Testing CADD API...")
    
    # Test with a known KCNH2 variant position
    # KCNH2 is on chr7, example position
    result = fetch_cadd_single("7", 150945368, "A", "C")
    
    if result:
        print(f"Success! CADD PHRED = {result['cadd_phred']}")
        print(f"  Raw: {result['cadd_raw']}")
    else:
        print("No result (variant may not be in CADD database)")
    
    # Test batch if gene provided
    if len(sys.argv) > 1:
        gene = sys.argv[1]
        print(f"\nFetching CADD for {gene}...")
        
        for i, var in enumerate(fetch_cadd_for_gene(gene)):
            if var.get('cadd_phred'):
                print(f"  {var['hgvs_p']}: PHRED={var['cadd_phred']:.1f}")
            if i >= 9:  # Limit to 10 for testing
                print("  ... (limited to 10 for testing)")
                break
