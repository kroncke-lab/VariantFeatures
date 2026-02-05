#!/usr/bin/env python3
"""Load ClinVar data for cardiac genes into the database."""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

from variantfeatures.fetchers.clinvar import fetch_clinvar, load_clinvar_to_db
from variantfeatures.database import VariantDB

CARDIAC_GENES = ['KCNH2', 'KCNQ1', 'SCN5A', 'RYR2']

def main():
    # Count variants per gene first
    print('Counting missense variants with protein annotations...')
    for gene in CARDIAC_GENES:
        count = sum(1 for _ in fetch_clinvar(gene))
        print(f'  {gene}: {count}')
    
    # Load into database
    print()
    print('Loading into database...')
    db = VariantDB()
    counts = load_clinvar_to_db(db, CARDIAC_GENES, verbose=True)
    
    print()
    print('Summary:')
    total = 0
    for gene, count in counts.items():
        print(f'  {gene}: {count} variants loaded')
        total += count
    print(f'  Total: {total}')
    
    # Verify
    print()
    print('Verification from DB:')
    for gene in CARDIAC_GENES:
        variants = db.get_gene_missense(gene)
        with_clinvar = sum(1 for v in variants if v.get('clinvar_id'))
        print(f'  {gene}: {len(variants)} total, {with_clinvar} with ClinVar IDs')
    
    db.close()

if __name__ == "__main__":
    main()
