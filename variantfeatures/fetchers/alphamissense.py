"""Fetch AlphaMissense predictions.

AlphaMissense provides pathogenicity scores for all possible single amino acid
substitutions in human proteins. Scores range from 0 (benign) to 1 (pathogenic).

Data source: https://storage.googleapis.com/dm_alphamissense/
Paper: Cheng et al., Science 2023
"""

import gzip
import os
from pathlib import Path
from typing import Iterator, Optional
import urllib.request

# Default cache location
DEFAULT_CACHE_DIR = Path(__file__).parent.parent.parent / "data" / "alphamissense"

# AlphaMissense data URL (amino acid substitutions, ~4GB compressed)
ALPHAMISSENSE_URL = "https://storage.googleapis.com/dm_alphamissense/AlphaMissense_aa_substitutions.tsv.gz"

# Gene name to UniProt mapping (common cardiac genes)
# Can be extended or loaded from file
GENE_TO_UNIPROT = {
    "KCNH2": "Q12809",
    "KCNQ1": "P51787",
    "SCN5A": "Q14524",
    "RYR2": "Q92736",
    "CACNA1C": "Q13936",
    "KCNJ2": "P63252",
    "KCNE1": "P15382",
    "KCNE2": "Q9Y6J6",
}


def get_cache_path() -> Path:
    """Get the cache directory, creating if needed."""
    cache_dir = Path(os.environ.get("ALPHAMISSENSE_CACHE", DEFAULT_CACHE_DIR))
    cache_dir.mkdir(parents=True, exist_ok=True)
    return cache_dir


def download_alphamissense(force: bool = False) -> Path:
    """
    Download AlphaMissense data file if not cached.
    
    Args:
        force: Re-download even if file exists
        
    Returns:
        Path to the downloaded file
    """
    cache_dir = get_cache_path()
    local_path = cache_dir / "AlphaMissense_aa_substitutions.tsv.gz"
    
    if local_path.exists() and not force:
        return local_path
    
    print(f"Downloading AlphaMissense data (~4GB)...")
    print(f"  From: {ALPHAMISSENSE_URL}")
    print(f"  To: {local_path}")
    
    urllib.request.urlretrieve(ALPHAMISSENSE_URL, local_path)
    print("  Done.")
    
    return local_path


def fetch_alphamissense(
    gene: str,
    uniprot_id: Optional[str] = None,
    cache_dir: Optional[Path] = None
) -> Iterator[dict]:
    """
    Fetch AlphaMissense scores for a gene.
    
    Args:
        gene: Gene symbol (e.g., "KCNH2")
        uniprot_id: Optional UniProt ID. If not provided, looks up from GENE_TO_UNIPROT
        cache_dir: Optional custom cache directory
        
    Yields:
        Dicts with:
            - hgvs_p: e.g., "p.Ala123Val"
            - position: int, amino acid position
            - ref_aa: reference amino acid (1-letter)
            - alt_aa: alternate amino acid (1-letter)
            - alphamissense_score: float 0-1
            - alphamissense_class: "likely_benign" | "ambiguous" | "likely_pathogenic"
    """
    # Resolve UniProt ID
    if uniprot_id is None:
        uniprot_id = GENE_TO_UNIPROT.get(gene.upper())
        if uniprot_id is None:
            raise ValueError(
                f"Unknown gene: {gene}. Provide uniprot_id or add to GENE_TO_UNIPROT mapping."
            )
    
    # Download if needed
    data_path = download_alphamissense()
    
    # Parse and filter
    # File format (TSV):
    # #CHROM  POS     REF     ALT     genome  uniprot_id  transcript_id   protein_variant am_pathogenicity        am_class
    # 1       69094   G       A       hg38    Q8NH21      ENST00000335137 OR4F5_A61T      0.1066  likely_benign
    
    print(f"Filtering AlphaMissense for {gene} ({uniprot_id})...")
    count = 0
    
    with gzip.open(data_path, 'rt') as f:
        header = None
        for line in f:
            # Skip comment lines (copyright notice)
            if line.startswith('#'):
                continue
            
            # First non-comment line is the header
            if header is None:
                header = line.strip().split('\t')
                continue
                
            fields = line.strip().split('\t')
            row = dict(zip(header, fields))
            
            # Filter by UniProt ID
            if row.get('uniprot_id') != uniprot_id:
                continue
            
            # Parse protein variant (e.g., "KCNH2_A123V" or "OR4F5_A61T")
            prot_var = row.get('protein_variant', '')
            if '_' in prot_var:
                var_part = prot_var.split('_', 1)[1]  # Get "A123V" part
            else:
                var_part = prot_var
            
            # Parse ref, position, alt from variant notation
            if len(var_part) >= 3:
                ref_aa = var_part[0]
                alt_aa = var_part[-1]
                try:
                    position = int(var_part[1:-1])
                except ValueError:
                    continue
            else:
                continue
            
            # Build HGVS notation
            # Convert 1-letter to 3-letter codes
            AA_MAP = {
                'A': 'Ala', 'C': 'Cys', 'D': 'Asp', 'E': 'Glu', 'F': 'Phe',
                'G': 'Gly', 'H': 'His', 'I': 'Ile', 'K': 'Lys', 'L': 'Leu',
                'M': 'Met', 'N': 'Asn', 'P': 'Pro', 'Q': 'Gln', 'R': 'Arg',
                'S': 'Ser', 'T': 'Thr', 'V': 'Val', 'W': 'Trp', 'Y': 'Tyr',
            }
            
            ref_3 = AA_MAP.get(ref_aa, ref_aa)
            alt_3 = AA_MAP.get(alt_aa, alt_aa)
            hgvs_p = f"p.{ref_3}{position}{alt_3}"
            
            try:
                score = float(row.get('am_pathogenicity', 0))
            except ValueError:
                score = None
            
            am_class = row.get('am_class', 'unknown')
            
            count += 1
            yield {
                'hgvs_p': hgvs_p,
                'position': position,
                'ref_aa': ref_aa,
                'alt_aa': alt_aa,
                'alphamissense_score': score,
                'alphamissense_class': am_class,
            }
    
    print(f"  Found {count} variants for {gene}")


def create_gene_cache(gene: str, output_dir: Optional[Path] = None) -> Path:
    """
    Pre-filter AlphaMissense data for a specific gene and cache as smaller file.
    
    This is useful for repeated queries on the same gene.
    
    Args:
        gene: Gene symbol
        output_dir: Where to save the filtered file
        
    Returns:
        Path to the gene-specific cache file
    """
    import csv
    
    if output_dir is None:
        output_dir = get_cache_path()
    
    output_path = output_dir / f"alphamissense_{gene.lower()}.tsv"
    
    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=[
            'hgvs_p', 'position', 'ref_aa', 'alt_aa', 
            'alphamissense_score', 'alphamissense_class'
        ], delimiter='\t')
        writer.writeheader()
        
        for variant in fetch_alphamissense(gene):
            writer.writerow(variant)
    
    print(f"Created gene cache: {output_path}")
    return output_path


if __name__ == '__main__':
    # Test with KCNH2
    import sys
    gene = sys.argv[1] if len(sys.argv) > 1 else "KCNH2"
    
    print(f"Testing AlphaMissense fetcher for {gene}")
    variants = list(fetch_alphamissense(gene))
    print(f"Found {len(variants)} variants")
    
    if variants:
        print("\nSample variants:")
        for v in variants[:5]:
            print(f"  {v['hgvs_p']}: {v['alphamissense_score']:.3f} ({v['alphamissense_class']})")
