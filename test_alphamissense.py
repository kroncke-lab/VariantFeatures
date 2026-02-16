#!/usr/bin/env python3

import sys
import os

# Add VariantFeatures to Python path
sys.path.insert(0, '/mnt/temp2/kronckbm/gitrepos/VariantFeatures')

from variantfeatures.fetchers.alphamissense import fetch_alphamissense

def test_alphamissense(variant_str):
    """Test AlphaMissense for a specific variant."""
    
    # Parse variant string (KCNH2 G628S)
    parts = variant_str.split(' ')
    if len(parts) != 2:
        print("Usage: python test_alphamissense.py <gene> <variant>")
        return
        
    gene = parts[0]
    
    # Extract position and amino acid change
    # For G628S: ref='G', pos=628, alt='S'
    var = parts[1]
    ref_aa = var[0]
    pos = int(var[1:-1])
    alt_aa = var[-1]
    
    print(f"Testing AlphaMissense for {gene} {var}")
    print(f"Gene: {gene}, Position: {pos}, Ref: {ref_aa}, Alt: {alt_aa}")
    print("-" * 50)
    
    # Get all variants for this gene
    variants = []
    print(f"Fetching AlphaMissense data for {gene}...")
    for v in fetch_alphamissense(gene):
        variants.append(v)
    
    print(f"Found {len(variants)} total variants for {gene}")
    
    # Find the specific variant
    specific_variants = [v for v in variants if v['position'] == pos and v['alt_aa'] == alt_aa]
    
    if specific_variants:
        for var in specific_variants:
            print(f"\nFound variant:")
            print(f"  HGVS: {var['hgvs_p']}")
            print(f"  Position: {var['position']}")
            print(f"  Change: {var['ref_aa']}{var['position']}{var['alt_aa']}")
            print(f"  AlphaMissense Score: {var['alphamissense_score']}")
            print(f"  AlphaMissense Class: {var['alphamissense_class']}")
            
            # Check if this is the G628S we're looking for
            if var['position'] == pos and var['alt_aa'] == alt_aa:
                print("\n✅ SPECIFIC VARIANT FOUND!")
                break
    else:
        print(f"\n❌ Variant {var} not found in AlphaMissense data")
        print(f"Available variants near position {pos}:")
        near_by = [v for v in variants if abs(v['position'] - pos) <= 5]
        for v in near_by[:10]:
            print(f"  {v['hgvs_p']}: {v['alphamissense_score']:.3f}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        test_alphamissense("KCNH2 G628S")
    else:
        test_alphamissense(f"{sys.argv[1]} {sys.argv[2]}")