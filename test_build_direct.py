#!/usr/bin/env python3
"""Test building database directly without CLI runner."""
import sys
sys.path.insert(0, '/mnt/temp2/kronckbm/gitrepos/VariantFeatures')

from variantfeatures.database import VariantDB
from variantfeatures.fetchers.gnomad import fetch_gnomad

print("Testing gnomAD build for KCNH2...", flush=True)

vdb = VariantDB()
print(f"Database: {vdb.db_path}", flush=True)

count = 0
skipped = 0

print("Fetching from gnomAD API...", flush=True)

for variant in fetch_gnomad('KCNH2'):
    if variant.get('hgvs_p'):
        vdb.upsert_missense(
            gene='KCNH2',
            hgvs_p=variant['hgvs_p'],
            hgvs_c=variant.get('hgvs_c'),
            gnomad_af=variant.get('gnomad_af'),
            gnomad_homozygotes=variant.get('gnomad_homozygotes'),
        )
        count += 1
    else:
        skipped += 1

print(f"Loaded {count} gnomAD variants (skipped {skipped} without hgvs_p)", flush=True)

# Check results
cur = vdb.conn.execute("""
    SELECT COUNT(*) as total,
           SUM(CASE WHEN gnomad_af IS NOT NULL THEN 1 ELSE 0 END) as with_af
    FROM variants_missense
    WHERE gene='KCNH2'
""")
row = cur.fetchone()
print(f"\nKCNH2 in DB: {row[0]} total, {row[1]} with gnomAD AF", flush=True)
