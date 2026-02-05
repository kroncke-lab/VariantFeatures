#!/usr/bin/env python3
"""Test gnomAD fetcher."""
import sys
import itertools
sys.path.insert(0, '/mnt/temp2/kronckbm/gitrepos/VariantFeatures')

from variantfeatures.fetchers.gnomad import fetch_gnomad

print("Testing gnomAD fetcher for KCNH2...")
print("This will make a network request to gnomad.broadinstitute.org")
print()

try:
    variants = list(itertools.islice(fetch_gnomad('KCNH2'), 10))
    print(f"\nSuccess! Found variants")
    print(f"Sample ({len(variants)} shown):")
    for v in variants:
        af_str = f"{v['gnomad_af']:.2e}" if v['gnomad_af'] else "N/A"
        print(f"  {v['variant_id']}: AF={af_str}, hgvs_p={v['hgvs_p']}")
except Exception as e:
    print(f"Error: {type(e).__name__}: {e}")
    import traceback
    traceback.print_exc()
