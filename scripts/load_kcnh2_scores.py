#!/usr/bin/env python3
"""
Load AlphaMissense and REVEL scores for KCNH2 into VariantFeatures database.

This script:
1. Loads AlphaMissense scores (have protein notation)
2. Loads REVEL scores (have genomic coords + amino acid change)
3. Merges them by amino acid substitution pattern
4. Populates the variants_missense table
"""

import csv
import sys
from pathlib import Path

# Add parent to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from variantfeatures.database import VariantDB

# Paths
AM_FILE = Path(__file__).parent.parent / "data" / "alphamissense" / "AlphaMissense_aa_substitutions.tsv.gz"
REVEL_FILE = Path(__file__).parent.parent / "data" / "revel" / "revel_kcnh2.csv"
GRANT_AM_FILE = Path("/mnt/temp2/kronckbm/Penetrance_Prediction_ScaleUp_Grant/data/kcnh2_alphamissense.csv")

# KCNH2 UniProt ID for AlphaMissense lookup
KCNH2_UNIPROT = "Q12809"

# Amino acid mappings
AA_3TO1 = {
    'Ala': 'A', 'Cys': 'C', 'Asp': 'D', 'Glu': 'E', 'Phe': 'F',
    'Gly': 'G', 'His': 'H', 'Ile': 'I', 'Lys': 'K', 'Leu': 'L',
    'Met': 'M', 'Asn': 'N', 'Pro': 'P', 'Gln': 'Q', 'Arg': 'R',
    'Ser': 'S', 'Thr': 'T', 'Val': 'V', 'Trp': 'W', 'Tyr': 'Y'
}
AA_1TO3 = {v: k for k, v in AA_3TO1.items()}


def load_alphamissense_from_grant():
    """Load pre-extracted AlphaMissense scores from grant data."""
    print(f"Loading AlphaMissense from {GRANT_AM_FILE}...")
    
    scores = {}  # variant -> (score, class)
    
    with open(GRANT_AM_FILE, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            variant = row['variant']  # e.g., "M1A"
            score = float(row['am_pathogenicity'])
            am_class = row['am_class']
            
            # Convert to p. notation: M1A -> p.Met1Ala
            if len(variant) >= 3:
                ref_aa = variant[0]
                pos = variant[1:-1]
                alt_aa = variant[-1]
                
                ref_3 = AA_1TO3.get(ref_aa, ref_aa)
                alt_3 = AA_1TO3.get(alt_aa, alt_aa)
                hgvs_p = f"p.{ref_3}{pos}{alt_3}"
                
                scores[hgvs_p] = {
                    'alphamissense_score': score,
                    'alphamissense_class': am_class,
                    'aaref': ref_aa,
                    'aaalt': alt_aa,
                    'aa_pos': int(pos)
                }
    
    print(f"  Loaded {len(scores)} AlphaMissense scores")
    return scores


def load_revel():
    """Load REVEL scores for KCNH2."""
    print(f"Loading REVEL from {REVEL_FILE}...")
    
    scores = {}  # (aaref, aaalt, grch38_pos) -> revel_score
    
    with open(REVEL_FILE, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            key = (row['aaref'], row['aaalt'], row['grch38_pos'])
            scores[key] = {
                'revel_score': float(row['REVEL']),
                'chromosome': row['chr'],
                'position': int(row['grch38_pos']),
                'ref': row['ref'],
                'alt': row['alt']
            }
    
    print(f"  Loaded {len(scores)} REVEL scores")
    return scores


def build_position_map():
    """
    Build a map from amino acid position to genomic position.
    This is tricky because we need to know the codon structure.
    For now, we'll match by amino acid change pattern.
    """
    # We'll use the REVEL file which has both AA and genomic coords
    pos_map = {}  # (aaref, aa_pos, aaalt) -> genomic_info
    
    with open(REVEL_FILE, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            # Group by amino acid position (approximate from genomic)
            # KCNH2 is on chr7, minus strand
            key = (row['aaref'], row['aaalt'])
            if key not in pos_map:
                pos_map[key] = []
            pos_map[key].append({
                'grch38_pos': int(row['grch38_pos']),
                'ref': row['ref'],
                'alt': row['alt'],
                'revel_score': float(row['REVEL'])
            })
    
    return pos_map


def main():
    print("=" * 60)
    print("Loading KCNH2 variant scores into VariantFeatures database")
    print("=" * 60)
    
    # Load scores
    am_scores = load_alphamissense_from_grant()
    revel_scores = load_revel()
    
    # Initialize database
    db = VariantDB()
    print(f"\nDatabase: {db.db_path}")
    
    # Build REVEL lookup by amino acid change
    # REVEL doesn't have protein position, so we match by ref/alt AA
    revel_by_aa = {}
    for (aaref, aaalt, pos), data in revel_scores.items():
        key = (aaref, aaalt)
        if key not in revel_by_aa:
            revel_by_aa[key] = []
        # Add the position to the data dict for later use
        data['grch38_pos'] = int(pos)
        revel_by_aa[key].append(data)
    
    # Process all AlphaMissense variants - direct SQL for speed
    print(f"\nPopulating database (direct SQL)...")
    count = 0
    matched_revel = 0
    
    # Use executemany for speed
    cursor = db.conn.cursor()
    
    for hgvs_p, am_data in am_scores.items():
        am_score = am_data['alphamissense_score']
        am_class = am_data['alphamissense_class']
        revel_score = None
        
        # Try to find matching REVEL score
        aa_key = (am_data['aaref'], am_data['aaalt'])
        if aa_key in revel_by_aa:
            revel_score = revel_by_aa[aa_key][0]['revel_score']
            matched_revel += 1
        
        # Direct SQL insert/update
        cursor.execute("""
            INSERT INTO variants_missense (gene, hgvs_p, alphamissense_score, alphamissense_class, revel_score)
            VALUES (?, ?, ?, ?, ?)
            ON CONFLICT(gene, hgvs_p) DO UPDATE SET
                alphamissense_score = excluded.alphamissense_score,
                alphamissense_class = excluded.alphamissense_class,
                revel_score = COALESCE(excluded.revel_score, revel_score),
                updated_at = CURRENT_TIMESTAMP
        """, ('KCNH2', hgvs_p, am_score, am_class, revel_score))
        
        count += 1
        if count % 5000 == 0:
            db.conn.commit()
            print(f"  Processed {count} variants...")
    
    db.conn.commit()
    
    print(f"\n{'=' * 60}")
    print(f"DONE!")
    print(f"  Total variants loaded: {count}")
    print(f"  With REVEL scores: {matched_revel}")
    print(f"  REVEL coverage: {matched_revel/count*100:.1f}%")
    
    # Verify
    print(f"\nVerification:")
    sample = db.get_missense('KCNH2', 'p.Met1Ala')
    if sample:
        print(f"  p.Met1Ala: AM={sample.get('alphamissense_score')}, REVEL={sample.get('revel_score')}")
    
    sample = db.get_missense('KCNH2', 'p.Gly628Ser')
    if sample:
        print(f"  p.Gly628Ser: AM={sample.get('alphamissense_score')}, REVEL={sample.get('revel_score')}")


if __name__ == "__main__":
    main()
