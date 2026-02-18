# VariantFeatures

## Mission
**Single command to build a complete variant annotation database for any gene.**

The goal is: `variantfeatures build --gene KCNH2` → SQLite database with all pathogenicity scores, population frequencies, structural features, and functional annotations ready for downstream analysis.

## Current Status (Feb 2026)
- **Phase 1 genes:** KCNH2, KCNQ1, SCN5A, RYR2 (cardiac channelopathies)
- **Phase 2 target:** ACMG SF v3.2 — 81 genes for secondary findings
- **Deadline:** June 2026 (R01 grant submission)

### Data Sources Implemented
| Source | Status | Data |
|--------|--------|------|
| AlphaMissense | ✅ Done | 1.1GB TSV, pathogenicity scores for all missense |
| REVEL | ✅ Done | 6.1GB TSV, ensemble missense scores |
| CADD | ✅ Done | REST API, deleteriousness scores |
| ClinVar | ✅ Done | Classifications, review status |
| gnomAD | ⚠️ Fetcher exists | Population frequencies via GraphQL |
| AlphaFold | ❌ Planned | 3D structure, pLDDT confidence |

### Current KCNH2 Coverage
- 22,021 variants with AlphaMissense scores
- 9,522 variants with REVEL scores
- CADD: API-based lookup (on demand)

## Architecture
```
VariantFeatures/
├── variantfeatures/
│   ├── cli.py              # Command-line interface
│   ├── database.py         # SQLite operations (upsert, query, export)
│   └── fetchers/
│       ├── alphamissense.py # AlphaMissense pathogenicity (TSV)
│       ├── revel.py         # REVEL ensemble scores (TSV)
│       ├── cadd.py          # CADD deleteriousness (REST API)
│       ├── clinvar.py       # ClinVar classifications (XML)
│       ├── gnomad.py        # Population frequencies (GraphQL)
│       └── lof.py           # Loss-of-function annotations
├── scripts/
│   └── load_kcnh2_scores.py # Example: populate KCNH2 scores
├── data/
│   ├── variants.db          # Main SQLite database
│   ├── alphamissense/       # AlphaMissense TSV (~1.1GB)
│   └── revel/               # REVEL data (~6.1GB uncompressed)
└── PIPELINE.md              # End-to-end workflow
```

## Key Commands
```bash
# Query variants
python -m variantfeatures query --gene KCNH2
python -m variantfeatures query --gene KCNH2 --format csv > kcnh2.csv

# Database stats
python -m variantfeatures stats

# Export for downstream analysis
python -m variantfeatures export --gene KCNH2 --output kcnh2_features.csv
```

## Target Architecture (TODO)
```bash
# Future: single command builds everything
variantfeatures build --gene BRCA1
# → Looks up transcripts, coords
# → Downloads/extracts AlphaMissense, REVEL scores
# → Fetches CADD scores via API
# → Fetches gnomAD frequencies
# → Downloads AlphaFold structure
# → Outputs unified SQLite + CSV
```

## Database Schema
```sql
variants_missense:
  - gene, hgvs_p, hgvs_c
  - chromosome, position, ref, alt, genome_build
  - alphamissense_score, alphamissense_class
  - revel_score
  - cadd_phred, cadd_raw
  - clinvar_significance, clinvar_stars
  - gnomad_af, gnomad_homozygotes
  - alphafold_plddt, domain

genes:
  - symbol, ensembl_id, canonical_transcript
  - pli, loeuf (constraint metrics)

penetrance_estimates:
  - variant linkage
  - penetrance_mean, ci_lower, ci_upper
  - model_version, n_carriers
```

## Related Repos
| Repo | Relationship |
|------|-------------|
| **GeneVariantFetcher** | Upstream — literature-mined carriers/phenotypes |
| **BayesianPenetranceEstimator** | Downstream — uses features as model inputs |
| **Variant_Browser** | Downstream — Django app displaying features |

## Contribution Notes
- Python 3.11+ required
- Large data files (AlphaMissense, REVEL) not in git — download on first use
- Rate limit external APIs (gnomAD, CADD)
- See TASKS.md for detailed work items
