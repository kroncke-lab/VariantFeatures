"""Fetch ClinVar annotations from variant_summary.txt.gz."""

import gzip
import re
from pathlib import Path
from typing import Iterator, Optional
from datetime import datetime

# ClinVar FTP: https://ftp.ncbi.nlm.nih.gov/pub/clinvar/
# variant_summary.txt.gz - all variants with clinical significance

CLINVAR_SUMMARY_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
DEFAULT_DATA_DIR = Path(__file__).parent.parent.parent / "data"

# Columns (0-indexed) from variant_summary.txt
COL_ALLELE_ID = 0
COL_TYPE = 1
COL_NAME = 2
COL_GENE_SYMBOL = 4
COL_CLINICAL_SIG = 6
COL_LAST_EVALUATED = 8
COL_ASSEMBLY = 16
COL_CHROMOSOME = 18
COL_START = 19
COL_REVIEW_STATUS = 24
COL_VARIATION_ID = 30
COL_POS_VCF = 31
COL_REF_VCF = 32
COL_ALT_VCF = 33

# Review status to stars mapping (ClinVar convention)
REVIEW_STATUS_STARS = {
    "practice guideline": 4,
    "reviewed by expert panel": 3,
    "criteria provided, multiple submitters, no conflicts": 2,
    "criteria provided, conflicting classifications": 1,
    "criteria provided, single submitter": 1,
    "no assertion for the individual variant": 0,
    "no assertion criteria provided": 0,
    "no classification for the single variant": 0,
    "no classifications from unflagged records": 0,
    "no classification provided": 0,
}


def parse_protein_change(name: str) -> Optional[str]:
    """
    Extract protein change (p.XXX) from ClinVar Name field.
    
    Examples:
        "NM_000238.4(KCNH2):c.1682C>T (p.Ala561Val)" -> "p.Ala561Val"
        "NM_000238.4(KCNH2):c.2893G>A (p.Gly965Arg)" -> "p.Gly965Arg"
    """
    # Match protein annotation in parentheses
    match = re.search(r'\(p\.([A-Za-z]{3}\d+[A-Za-z]{3})\)', name)
    if match:
        return f"p.{match.group(1)}"
    
    # Also try to match nonsense/stop variants like p.Arg534Ter
    match = re.search(r'\(p\.([A-Za-z]{3}\d+(?:Ter|\*))\)', name)
    if match:
        aa_change = match.group(1).replace('*', 'Ter')
        return f"p.{aa_change}"
    
    return None


def parse_hgvs_c(name: str) -> Optional[str]:
    """
    Extract coding change (c.XXX) from ClinVar Name field.
    
    Examples:
        "NM_000238.4(KCNH2):c.1682C>T (p.Ala561Val)" -> "c.1682C>T"
    """
    match = re.search(r':c\.([^\s(]+)', name)
    if match:
        return f"c.{match.group(1)}"
    return None


def parse_date(date_str: str) -> Optional[str]:
    """Parse ClinVar date to ISO format."""
    if not date_str or date_str == '-':
        return None
    try:
        dt = datetime.strptime(date_str, "%b %d, %Y")
        return dt.strftime("%Y-%m-%d")
    except ValueError:
        return None


def get_review_stars(review_status: str) -> int:
    """Convert review status to star rating (0-4)."""
    return REVIEW_STATUS_STARS.get(review_status.lower(), 0)


def fetch_clinvar(
    gene: str,
    data_file: Optional[Path] = None,
    assembly: str = "GRCh38",
    include_all_types: bool = False,
) -> Iterator[dict]:
    """
    Fetch ClinVar annotations for a gene from local variant_summary.txt.gz.
    
    Args:
        gene: Gene symbol (e.g., "KCNH2")
        data_file: Path to variant_summary.txt.gz (default: data/variant_summary.txt.gz)
        assembly: Genome assembly to filter on ("GRCh38" or "GRCh37")
        include_all_types: Include all variant types, not just SNVs
    
    Yields dicts with:
        - hgvs_p: e.g., "p.Ala561Val"
        - hgvs_c: e.g., "c.1682C>T"
        - clinvar_id: int (VariationID)
        - clinvar_significance: e.g., "Pathogenic", "Likely benign"
        - clinvar_review_status: e.g., "criteria provided, single submitter"
        - clinvar_stars: int (0-4)
        - clinvar_last_evaluated: ISO date string
        - chromosome: str
        - position: int (VCF position)
        - ref: str
        - alt: str
    """
    if data_file is None:
        data_file = DEFAULT_DATA_DIR / "variant_summary.txt.gz"
    
    if not data_file.exists():
        raise FileNotFoundError(
            f"ClinVar data file not found: {data_file}\n"
            f"Download with: wget {CLINVAR_SUMMARY_URL} -O {data_file}"
        )
    
    seen_variants = set()  # Track by (gene, hgvs_p) to avoid duplicates
    
    with gzip.open(data_file, 'rt', encoding='utf-8') as f:
        # Skip header
        next(f)
        
        for line in f:
            fields = line.strip().split('\t')
            
            # Filter by gene
            if fields[COL_GENE_SYMBOL] != gene:
                continue
            
            # Filter by assembly
            if fields[COL_ASSEMBLY] != assembly:
                continue
            
            # Filter by variant type (default: SNVs only for missense)
            if not include_all_types and fields[COL_TYPE] != "single nucleotide variant":
                continue
            
            # Parse protein change
            hgvs_p = parse_protein_change(fields[COL_NAME])
            if hgvs_p is None:
                continue
            
            # Skip duplicates (same variant reported multiple times)
            variant_key = (gene, hgvs_p)
            if variant_key in seen_variants:
                continue
            seen_variants.add(variant_key)
            
            # Parse other fields
            hgvs_c = parse_hgvs_c(fields[COL_NAME])
            review_status = fields[COL_REVIEW_STATUS]
            
            # Get position/ref/alt if available
            pos_str = fields[COL_POS_VCF]
            ref = fields[COL_REF_VCF]
            alt = fields[COL_ALT_VCF]
            
            try:
                position = int(pos_str) if pos_str and pos_str != 'na' else None
            except ValueError:
                position = None
            
            yield {
                "hgvs_p": hgvs_p,
                "hgvs_c": hgvs_c,
                "clinvar_id": int(fields[COL_VARIATION_ID]) if fields[COL_VARIATION_ID] else None,
                "clinvar_significance": fields[COL_CLINICAL_SIG],
                "clinvar_review_status": review_status,
                "clinvar_stars": get_review_stars(review_status),
                "clinvar_last_evaluated": parse_date(fields[COL_LAST_EVALUATED]),
                "chromosome": fields[COL_CHROMOSOME] if fields[COL_CHROMOSOME] != 'na' else None,
                "position": position,
                "ref": ref if ref and ref != 'na' else None,
                "alt": alt if alt and alt != 'na' else None,
            }


def load_clinvar_to_db(
    db,
    genes: list[str],
    data_file: Optional[Path] = None,
    verbose: bool = False,
) -> dict[str, int]:
    """
    Load ClinVar data for multiple genes into the database.
    
    Args:
        db: VariantDB instance
        genes: List of gene symbols
        data_file: Path to variant_summary.txt.gz
        verbose: Print progress
    
    Returns:
        Dict of gene -> count of variants loaded
    """
    counts = {}
    
    for gene in genes:
        if verbose:
            print(f"Loading ClinVar data for {gene}...")
        
        count = 0
        for variant in fetch_clinvar(gene, data_file):
            db.upsert_missense(
                gene=gene,
                hgvs_p=variant["hgvs_p"],
                hgvs_c=variant["hgvs_c"],
                clinvar_id=variant["clinvar_id"],
                clinvar_significance=variant["clinvar_significance"],
                clinvar_review_status=variant["clinvar_review_status"],
                clinvar_stars=variant["clinvar_stars"],
                clinvar_last_evaluated=variant["clinvar_last_evaluated"],
                chromosome=variant["chromosome"],
                position=variant["position"],
                ref=variant["ref"],
                alt=variant["alt"],
            )
            count += 1
        
        counts[gene] = count
        if verbose:
            print(f"  Loaded {count} variants")
    
    return counts


if __name__ == "__main__":
    # Quick test
    import sys
    gene = sys.argv[1] if len(sys.argv) > 1 else "KCNH2"
    
    print(f"Fetching ClinVar data for {gene}...")
    count = 0
    sample = None
    for v in fetch_clinvar(gene):
        count += 1
        if count == 1:
            sample = v
    
    print(f"Total variants: {count}")
    if sample:
        print(f"Sample: {sample}")
