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
│   └── full_gene_pipeline.sh # End-to-end per-gene annotation pipeline
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
python -m variantfeatures export --gene KCNH2 --layout transcript-wide --output kcnh2_isoforms.csv

# Publish gene-scoped SQLite artifacts for Variant_Browser import
python -m variantfeatures publish --gene KCNH2 --dry-run
python -m variantfeatures publish --gene KCNH2
```

## Target Architecture (TODO)
```bash
# Future: single command builds everything
variantfeatures build --gene BRCA1
# → Looks up selected transcripts/isoforms, coords
# → Downloads/extracts AlphaMissense, REVEL scores
# → Fetches CADD scores via API
# → Fetches gnomAD frequencies
# → Uses pext/expression as an isoform relevance adjudicator
# → Downloads AlphaFold structure
# → Outputs unified SQLite + CSV
```

## Database Schema (normalized)
The old flat `variants_missense` / `variants_lof` tables were removed; data now
lives in a normalized, identity-centric schema:
```sql
variants:                 -- canonical identity, one row per GRCh38 SNV
  - chromosome, position, ref, alt, variant_type, hgvs_g, ca_id, vrs_id

variant_consequences:     -- one row per transcript per source
  - variant_id, gene_symbol, gene_ensembl, transcript_id, consequence
  - hgvs_c, hgvs_p, aa_pos, is_mane_select, is_canonical, source

transcripts:              -- isoforms selected during build/enumeration
  - transcript_id, refseq_match, protein_id, biotype, cds_length
  - is_mane_select, is_mane_plus_clinical, is_canonical, appris

annotations_*:            -- per feature family, keyed by variant_id + version
  - pathogenicity, population, clinical, conservation, splice, expression, structure

gene_constraint:          -- per-gene gnomAD metrics (pli, oe_lof; LOEUF = oe_lof_upper)
genes:                    -- gene identity (symbol, ensembl_id, canonical_transcript)
variant_aliases:          -- rsID / ClinVar / gnomAD / CA-ID / HGVS aliases

penetrance_estimates:     -- downstream BayesianPenetranceEstimator output
  - variant_id, penetrance_mean, ci_lower, ci_upper, model_version, n_carriers
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
