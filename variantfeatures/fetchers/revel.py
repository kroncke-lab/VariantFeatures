"""Fetch REVEL scores for missense variants.

REVEL (Rare Exome Variant Ensemble Learner) is an ensemble method for predicting
pathogenicity of missense variants using 13 different tools.

Data source: https://sites.google.com/site/revelgenomics/
File: revel_with_transcript_ids (TSV, ~6GB uncompressed)

REVEL scores range from 0 to 1:
- Higher scores = more likely pathogenic
- Suggested threshold: 0.5 (can be adjusted per application)
"""

import csv
from pathlib import Path
from typing import Iterator, Optional, Dict

# Default REVEL data file location
DEFAULT_REVEL_PATH = Path(__file__).parent.parent.parent / "data" / "revel" / "revel_with_transcript_ids"

# Gene to Ensembl transcript mappings (canonical transcripts)
# These are the primary transcripts used for variant annotation
GENE_TRANSCRIPTS = {
    "KCNH2": ["ENST00000262186", "ENST00000392968", "ENST00000330883"],
    "KCNQ1": ["ENST00000155840", "ENST00000621062"],
    "SCN5A": ["ENST00000333535", "ENST00000413689"],
    "RYR2": ["ENST00000366574"],
    # Add more genes as needed
}


def fetch_revel(
    gene: str,
    revel_path: Optional[Path] = None,
    transcript_ids: Optional[list] = None,
) -> Iterator[dict]:
    """
    Fetch REVEL scores for all missense variants in a gene.
    
    Args:
        gene: Gene symbol (e.g., "KCNH2")
        revel_path: Path to REVEL data file (uses default if not provided)
        transcript_ids: List of Ensembl transcript IDs to filter by
                       (uses GENE_TRANSCRIPTS mapping if not provided)
    
    Yields:
        Dicts with:
            - chromosome: Chromosome number
            - hg19_pos: hg19/GRCh37 position
            - grch38_pos: GRCh38 position
            - ref: Reference allele
            - alt: Alternate allele
            - aaref: Reference amino acid
            - aaalt: Alternate amino acid
            - revel_score: REVEL score (0-1)
            - transcript_ids: List of matching transcript IDs
            - hgvs_p: Protein change string (e.g., "p.Arg528His")
    """
    revel_file = revel_path or DEFAULT_REVEL_PATH
    
    if not revel_file.exists():
        raise FileNotFoundError(
            f"REVEL data file not found at {revel_file}. "
            "Download from: https://sites.google.com/site/revelgenomics/downloads"
        )
    
    # Get transcript IDs to search for
    if transcript_ids is None:
        transcript_ids = GENE_TRANSCRIPTS.get(gene.upper(), [])
        if not transcript_ids:
            raise ValueError(
                f"No transcript IDs configured for gene {gene}. "
                "Provide transcript_ids parameter or add to GENE_TRANSCRIPTS."
            )
    
    print(f"Fetching REVEL scores for {gene} (transcripts: {transcript_ids})...")
    
    # Amino acid 3-letter to 1-letter mapping
    AA_MAP = {
        'A': 'Ala', 'C': 'Cys', 'D': 'Asp', 'E': 'Glu', 'F': 'Phe',
        'G': 'Gly', 'H': 'His', 'I': 'Ile', 'K': 'Lys', 'L': 'Leu',
        'M': 'Met', 'N': 'Asn', 'P': 'Pro', 'Q': 'Gln', 'R': 'Arg',
        'S': 'Ser', 'T': 'Thr', 'V': 'Val', 'W': 'Trp', 'Y': 'Tyr'
    }
    
    count = 0
    with open(revel_file, 'r') as f:
        reader = csv.DictReader(f)
        
        for row in reader:
            # Check if this row matches any of our transcripts
            row_transcripts = row.get('Ensembl_transcriptid', '').split(';')
            matching = [t for t in transcript_ids if t in row_transcripts]
            
            if not matching:
                continue
            
            # Parse score
            try:
                score = float(row['REVEL'])
            except (ValueError, KeyError):
                continue
            
            # Build HGVS protein notation
            aaref = row.get('aaref', '')
            aaalt = row.get('aaalt', '')
            
            # We don't have position in the amino acid sequence from this file,
            # so we'll need to match by ref/alt later or get from another source
            # For now, just return what we have
            
            yield {
                'chromosome': row.get('chr', ''),
                'hg19_pos': row.get('hg19_pos', ''),
                'grch38_pos': row.get('grch38_pos', ''),
                'ref': row.get('ref', ''),
                'alt': row.get('alt', ''),
                'aaref': aaref,
                'aaalt': aaalt,
                'revel_score': score,
                'transcript_ids': matching,
            }
            count += 1
    
    print(f"  Found {count} REVEL scores for {gene}")


def lookup_revel_by_position(
    chromosome: str,
    position: int,
    ref: str,
    alt: str,
    revel_path: Optional[Path] = None,
    genome_build: str = "GRCh38"
) -> Optional[float]:
    """
    Look up REVEL score for a specific variant by genomic coordinates.
    
    Args:
        chromosome: Chromosome (e.g., "7" or "chr7")
        position: Genomic position
        ref: Reference allele
        alt: Alternate allele
        revel_path: Path to REVEL file
        genome_build: "GRCh38" or "hg19"
    
    Returns:
        REVEL score or None if not found
    """
    revel_file = revel_path or DEFAULT_REVEL_PATH
    
    if not revel_file.exists():
        return None
    
    # Normalize chromosome
    chrom = chromosome.replace("chr", "")
    pos_col = "grch38_pos" if genome_build == "GRCh38" else "hg19_pos"
    
    with open(revel_file, 'r') as f:
        reader = csv.DictReader(f)
        
        for row in reader:
            if (row.get('chr') == chrom and 
                row.get(pos_col) == str(position) and
                row.get('ref') == ref and 
                row.get('alt') == alt):
                try:
                    return float(row['REVEL'])
                except (ValueError, KeyError):
                    return None
    
    return None


if __name__ == "__main__":
    import sys
    import itertools
    
    gene = sys.argv[1] if len(sys.argv) > 1 else "KCNH2"
    
    print(f"Testing REVEL fetcher for {gene}")
    print(f"REVEL file: {DEFAULT_REVEL_PATH}")
    print(f"File exists: {DEFAULT_REVEL_PATH.exists()}")
    
    if DEFAULT_REVEL_PATH.exists():
        variants = list(itertools.islice(fetch_revel(gene), 10))
        print(f"\nSample variants ({len(variants)} shown):")
        for v in variants:
            print(f"  chr{v['chromosome']}:{v['grch38_pos']} {v['ref']}>{v['alt']} "
                  f"({v['aaref']}>{v['aaalt']}) REVEL={v['revel_score']:.3f}")
