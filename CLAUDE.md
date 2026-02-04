# VariantFeatures

## Objective
Aggregate predictive features for genetic variants from multiple sources into a unified SQLite database.

## Current Status
- **Phase 1 genes:** KCNH2, KCNQ1, SCN5A, RYR2 (cardiac channelopathies)
- **Phase 2 target:** ACMG SF v3.2 — 81 genes for secondary findings
- **Deadline:** June 2026 (R01 grant)

## Architecture
```
VariantFeatures/
├── variantfeatures/
│   ├── cli.py              # Command-line interface
│   ├── database.py         # SQLite operations
│   └── fetchers/
│       ├── alphamissense.py # AlphaMissense pathogenicity scores
│       ├── clinvar.py       # ClinVar classifications
│       ├── gnomad.py        # Population frequencies
│       └── lof.py           # Loss-of-function annotations (LOFTEE, NMD)
├── scripts/
│   └── load_clinvar.py     # Bulk ClinVar loading
├── data/                   # Local data cache
└── PIPELINE.md             # End-to-end workflow documentation
```

## Key Commands
```bash
# Build database for a gene
python -m variantfeatures build --gene KCNH2

# Query variants
python -m variantfeatures query --gene KCNH2 --format csv > kcnh2.csv
python -m variantfeatures query --gene KCNH2 --format json > kcnh2.json
```

## Data Sources
| Source | Method | Features |
|--------|--------|----------|
| AlphaMissense | TSV (~4GB) | Pathogenicity score, class |
| ClinVar | Bulk XML | Classification, review status, stars |
| gnomAD | GraphQL API | AF, AF_popmax, homozygotes |
| REVEL | TSV | Ensemble missense score |
| CADD | TSV | PHRED + raw scores |

## Schema Overview
- `variants_missense` - Missense variants with scores
- `variants_lof` - Loss-of-function variants (nonsense, frameshift, splice)
- `genes` - Gene-level annotations (pLI, LOEUF)
- `penetrance_estimates` - From BayesianPenetranceEstimator

## Related Repos
- **GeneVariantFetcher** - Upstream: literature-derived carriers/phenotypes
- **BayesianPenetranceEstimator** - Downstream: uses features as model inputs
- **Variant_Browser** - Downstream: displays aggregated features

## Next Steps
See TASKS.md for detailed task list. Priority:
1. AlphaMissense integration
2. gnomAD v4 integration
3. Fix CLI build/query commands
4. REVEL score integration
5. Phase 1 gene validation
