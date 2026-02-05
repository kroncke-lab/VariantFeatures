#!/usr/bin/env python3
"""Test gnomAD fetcher - check for variants with hgvs_p."""
import sys
sys.path.insert(0, '/mnt/temp2/kronckbm/gitrepos/VariantFeatures')

from variantfeatures.fetchers.gnomad import fetch_gnomad

print("Testing gnomAD fetcher for KCNH2...")
print("Looking for variants WITH hgvs_p annotations...")
print()

try:
    missense_count = 0
    total_count = 0
    sample_missense = []
    
    for v in fetch_gnomad('KCNH2'):
        total_count += 1
        if v.get('hgvs_p'):
            missense_count += 1
            if len(sample_missense) < 15:
                sample_missense.append(v)
    
    print(f"\nTotal variants: {total_count}")
    print(f"Variants with hgvs_p: {missense_count}")
    print(f"\nSample missense variants:")
    for v in sample_missense:
        af_str = f"{v['gnomad_af']:.2e}" if v['gnomad_af'] else "N/A"
        print(f"  {v['hgvs_p']}: AF={af_str}, consequence={v.get('consequence')}")

except Exception as e:
    print(f"Error: {type(e).__name__}: {e}")
    import traceback
    traceback.print_exc()
