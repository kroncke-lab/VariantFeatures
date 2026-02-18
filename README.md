# VariantFeatures

**Aggregate pathogenicity scores and variant annotations from multiple sources into a unified SQLite database.**

The goal: give it a gene name, get a complete annotation database ready for downstream penetrance modeling and clinical interpretation.

## Quick Start

```bash
# Query existing data
python -m variantfeatures query --gene KCNH2
python -m variantfeatures query --gene KCNH2 --format csv > kcnh2.csv

# Check database coverage
python -m variantfeatures stats

# Export for downstream analysis
python -m variantfeatures export --gene KCNH2 --output kcnh2_features.csv
```

## What's Included

### Pathogenicity Scores
| Source | Description | Coverage |
|--------|-------------|----------|
| **AlphaMissense** | Deep learning missense pathogenicity | All possible missense |
| **REVEL** | Ensemble of 13 pathogenicity tools | Observed variants |
| **CADD** | Deleteriousness based on 60+ features | All possible SNVs |

### Clinical Annotations
| Source | Description |
|--------|-------------|
| **ClinVar** | Clinical classifications + review status |
| **gnomAD** | Population allele frequencies |

### Planned
- **AlphaFold** — 3D structure, per-residue confidence (pLDDT)
- **Protein domains** — Functional region annotations

## Current Coverage (KCNH2)

| Source | Variants |
|--------|----------|
| AlphaMissense | 22,021 |
| REVEL | 9,522 |
| CADD | On-demand API |
| ClinVar | 699 |

## Target Architecture

```bash
# Future: single command builds everything for any gene
variantfeatures build --gene BRCA1

# This will:
# 1. Look up canonical transcript, genomic coordinates
# 2. Extract AlphaMissense scores for the gene
# 3. Extract REVEL scores for the gene
# 4. Fetch CADD scores via API
# 5. Fetch gnomAD frequencies via GraphQL
# 6. Download AlphaFold structure
# 7. Output unified SQLite database + CSV export
```

## Database Schema

```sql
-- Main table: missense variants with all annotations
variants_missense (
    gene, hgvs_p, hgvs_c,
    chromosome, position, ref, alt, genome_build,
    
    -- Pathogenicity scores
    alphamissense_score, alphamissense_class,
    revel_score,
    cadd_phred, cadd_raw,
    
    -- Clinical
    clinvar_significance, clinvar_stars,
    gnomad_af, gnomad_homozygotes,
    
    -- Structural (planned)
    alphafold_plddt, domain
)

-- Loss-of-function variants
variants_lof (gene, hgvs_c, lof_type, loftee_confidence, nmd_escape, ...)

-- Gene-level constraint
genes (symbol, pli, loeuf, canonical_transcript)

-- Penetrance estimates (from downstream modeling)
penetrance_estimates (variant_id, penetrance_mean, ci_lower, ci_upper, ...)
```

## Project Structure

```
VariantFeatures/
├── variantfeatures/
│   ├── cli.py              # Command-line interface
│   ├── database.py         # SQLite operations
│   └── fetchers/
│       ├── alphamissense.py # ~1.1GB TSV
│       ├── revel.py         # ~6.1GB TSV
│       ├── cadd.py          # REST API
│       ├── clinvar.py       # Bulk XML
│       ├── gnomad.py        # GraphQL API
│       └── lof.py           # LOFTEE annotations
├── scripts/
│   └── load_kcnh2_scores.py # Example loader
├── data/
│   ├── variants.db          # Main database
│   ├── alphamissense/       # Cached TSV
│   └── revel/               # Cached TSV
└── PIPELINE.md              # Full workflow docs
```

## Data Sources & Sizes

| Source | Method | Size | Download |
|--------|--------|------|----------|
| AlphaMissense | TSV | 1.1GB | [Google Cloud](https://console.cloud.google.com/storage/browser/dm_alphamissense) |
| REVEL | TSV | 6.1GB | [REVEL Site](https://sites.google.com/site/revelgenomics/) |
| CADD | API | On-demand | [CADD API](https://cadd.gs.washington.edu/api) |
| ClinVar | XML | ~1GB | [NCBI FTP](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/) |
| gnomAD | API | On-demand | [gnomAD API](https://gnomad.broadinstitute.org/api) |

## Related Projects

This is part of the **Kroncke Lab Variant Interpretation Pipeline**:

```
Literature → GeneVariantFetcher → VariantFeatures → BayesianPenetranceEstimator → Variant_Browser
                                       ↑
                                   (you are here)
```

| Repo | Purpose |
|------|---------|
| [GeneVariantFetcher](https://github.com/kroncke-lab/GeneVariantFetcher) | Mine literature for carriers/phenotypes |
| **VariantFeatures** | Aggregate pathogenicity scores |
| [BayesianPenetranceEstimator](https://github.com/kroncke-lab/BayesianPenetranceEstimator) | Model penetrance from features |
| [Variant_Browser](https://github.com/kroncke-lab/Variant_Browser) | Clinical-facing Django app |

## Requirements

- Python 3.11+
- ~10GB disk space for cached data files
- Internet for API-based fetchers (gnomAD, CADD)

## Installation

```bash
git clone https://github.com/kroncke-lab/VariantFeatures.git
cd VariantFeatures
pip install -e .
```

## License

MIT
