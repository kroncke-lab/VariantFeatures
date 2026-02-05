#!/usr/bin/env python3
"""Load ClinVar data for cardiac genes - single-pass optimized version."""

import sys
import gzip
import re
from pathlib import Path
from datetime import datetime

sys.path.insert(0, str(Path(__file__).parent.parent))

from variantfeatures.database import VariantDB

CARDIAC_GENES = {'KCNH2', 'KCNQ1', 'SCN5A', 'RYR2'}
DATA_FILE = Path(__file__).parent.parent / "data" / "variant_summary.txt.gz"

# Column indices
COL_TYPE = 1
COL_NAME = 2
COL_GENE_SYMBOL = 4
COL_CLINICAL_SIG = 6
COL_LAST_EVALUATED = 8
COL_ASSEMBLY = 16
COL_CHROMOSOME = 18
COL_REVIEW_STATUS = 24
COL_VARIATION_ID = 30
COL_POS_VCF = 31
COL_REF_VCF = 32
COL_ALT_VCF = 33

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

def parse_protein_change(name):
    match = re.search(r'\(p\.([A-Za-z]{3}\d+[A-Za-z]{3})\)', name)
    if match:
        return f"p.{match.group(1)}"
    match = re.search(r'\(p\.([A-Za-z]{3}\d+(?:Ter|\*))\)', name)
    if match:
        return f"p.{match.group(1).replace('*', 'Ter')}"
    return None

def parse_hgvs_c(name):
    match = re.search(r':c\.([^\s(]+)', name)
    return f"c.{match.group(1)}" if match else None

def parse_date(date_str):
    if not date_str or date_str == '-':
        return None
    try:
        return datetime.strptime(date_str, "%b %d, %Y").strftime("%Y-%m-%d")
    except ValueError:
        return None

def main():
    print(f"Parsing {DATA_FILE}...")
    print(f"Looking for genes: {', '.join(sorted(CARDIAC_GENES))}")
    print()
    
    db = VariantDB()
    counts = {gene: 0 for gene in CARDIAC_GENES}
    seen = set()  # Track (gene, hgvs_p) to avoid duplicates
    total_lines = 0
    
    with gzip.open(DATA_FILE, 'rt', encoding='utf-8') as f:
        next(f)  # Skip header
        
        for line in f:
            total_lines += 1
            if total_lines % 500000 == 0:
                print(f"  Processed {total_lines:,} lines...")
            
            fields = line.strip().split('\t')
            gene = fields[COL_GENE_SYMBOL]
            
            # Filter by gene
            if gene not in CARDIAC_GENES:
                continue
            
            # Filter by assembly (GRCh38)
            if fields[COL_ASSEMBLY] != "GRCh38":
                continue
            
            # Filter by variant type (SNVs only)
            if fields[COL_TYPE] != "single nucleotide variant":
                continue
            
            # Parse protein change
            hgvs_p = parse_protein_change(fields[COL_NAME])
            if hgvs_p is None:
                continue
            
            # Skip duplicates
            key = (gene, hgvs_p)
            if key in seen:
                continue
            seen.add(key)
            
            # Parse fields
            hgvs_c = parse_hgvs_c(fields[COL_NAME])
            review_status = fields[COL_REVIEW_STATUS]
            pos_str = fields[COL_POS_VCF]
            ref = fields[COL_REF_VCF]
            alt = fields[COL_ALT_VCF]
            
            try:
                position = int(pos_str) if pos_str and pos_str != 'na' else None
            except ValueError:
                position = None
            
            # Insert into database
            db.upsert_missense(
                gene=gene,
                hgvs_p=hgvs_p,
                hgvs_c=hgvs_c,
                clinvar_id=int(fields[COL_VARIATION_ID]) if fields[COL_VARIATION_ID] else None,
                clinvar_significance=fields[COL_CLINICAL_SIG],
                clinvar_review_status=review_status,
                clinvar_stars=REVIEW_STATUS_STARS.get(review_status.lower(), 0),
                clinvar_last_evaluated=parse_date(fields[COL_LAST_EVALUATED]),
                chromosome=fields[COL_CHROMOSOME] if fields[COL_CHROMOSOME] != 'na' else None,
                position=position,
                ref=ref if ref and ref != 'na' else None,
                alt=alt if alt and alt != 'na' else None,
            )
            counts[gene] += 1
    
    print(f"\nProcessed {total_lines:,} total lines")
    print("\n=== Results ===")
    for gene in sorted(CARDIAC_GENES):
        print(f"  {gene}: {counts[gene]} variants")
    print(f"  Total: {sum(counts.values())}")
    
    # Verify
    print("\n=== Verification from DB ===")
    for gene in sorted(CARDIAC_GENES):
        variants = db.get_gene_missense(gene)
        with_clinvar = sum(1 for v in variants if v.get('clinvar_id'))
        print(f"  {gene}: {len(variants)} total, {with_clinvar} with ClinVar IDs")
    
    db.close()
    print("\nDone!")

if __name__ == "__main__":
    main()
