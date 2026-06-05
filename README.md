# VariantFeatures

**Aggregate pathogenicity scores and variant annotations from multiple sources into a unified SQLite database.**

The goal: give it a gene name, get a complete annotation database ready for downstream penetrance modeling and clinical interpretation.

## Quick Start

```bash
# Build normalized canonical variants + queued annotations
.venv/bin/python -m variantfeatures build --gene KCNH2 --no-run

# Query existing data
.venv/bin/python -m variantfeatures query --gene KCNH2
.venv/bin/python -m variantfeatures query --gene KCNH2 --format csv > kcnh2.csv

# Check database coverage
.venv/bin/python -m variantfeatures stats

# Export normalized feature matrix for downstream analysis
.venv/bin/python -m variantfeatures export --gene KCNH2 --output kcnh2_features.csv

# Show where normalized feature families live
.venv/bin/python -m variantfeatures feature-schema
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

### LoF / Splice / Structure
| Source | Description |
|--------|-------------|
| **AlphaFold DB** | Per-residue pLDDT at missense or truncation position |
| **SpliceAI** | Parsed from VEP plugin output (DS/DP + overall max) |
| **dbscSNV / MaxEntScan** | Parsed from VEP plugin output |
| **LOFTEE / VEP NMD / nmd-rule** | LoF confidence and NMD escape/trigger annotations |
| **gnomAD pext / AbSplice / AbExp / NMDEP / NMDetective** | Local tabular importers for large or external model outputs |

## Current Coverage (KCNH2)

| Source | Variants |
|--------|----------|
| AlphaMissense | 22,021 |
| REVEL | 9,522 |
| CADD | On-demand API |
| ClinVar | 699 |

## Target Architecture

```bash
# Normalized canonical build for any gene
variantfeatures build --gene BRCA1

# This:
# 1. Looks up the canonical/MANE transcript and enumerates coding SNVs
# 2. Queues/runs normalized sources such as MyVariant/dbNSFP, gnomAD,
#    AlphaMissense, REVEL, AlphaFold, pext, and NMD rules
# 3. Stores features in family-specific normalized tables
# 4. Exports one-row-per-variant wide CSVs or long audit tables
```

## Database Schema

The normalized path uses `variants` for canonical GRCh38 alleles,
`variant_consequences` for transcript effects, and feature-family tables:

| Feature family | Table | Export prefix |
|---|---|---|
| Pathogenicity / functional scores | `annotations_pathogenicity` | `pathogenicity.<predictor>...` |
| Population frequencies | `annotations_population` | `population.<dataset>.<pop>...` |
| Clinical assertions | `annotations_clinical` | `clinical.<source>...` |
| Splice predictors | `annotations_splice` | `splice.<predictor>...` |
| pext / expression features | `annotations_expression` | `expression.<metric>.<dataset>.<tissue>...` |
| Structure / domains | `annotations_structure` | `structure.<feature>...` |
| Conservation | `annotations_conservation` | `conservation.<metric>...` |
| Gene constraint | `gene_constraint` | `gene_constraint.<dataset>...` |

Legacy `variants_missense` and `variants_lof` tables are still present for old
loaders. Use `variantfeatures query --legacy`, `variantfeatures stats --legacy`,
or `variantfeatures export --legacy` when that old shape is needed.

## On-Disk Storage

VariantFeatures stores the built annotation database in a local SQLite file:

```text
data/variants.db
```

The database file is ignored by git (`data/` and `*.db` are in `.gitignore`)
because it is a generated data artifact. The current local snapshot is
8,404,180,992 bytes, reported by `ls -lh` as 7.8G.

The normalized schema separates variant identity, aliases, transcript effects,
and annotation features:

| Table | What it stores | Key columns |
|---|---|---|
| `variants` | One canonical row per allele, using GRCh38 normalized VCF-style coordinates | `id`, `chromosome`, `position`, `ref`, `alt`, `variant_type`, `hgvs_g`, `ca_id`, `vrs_id` |
| `variant_aliases` | All alternate names that resolve to the same allele | `variant_id`, `alias_type`, `alias_value`, `source`, `fetched_at` |
| `variant_consequences` | Per-transcript gene/cDNA/protein effects | `variant_id`, `transcript_id`, `gene_symbol`, `gene_ensembl`, `consequence`, `hgvs_c`, `hgvs_p`, `aa_pos`, `aa_ref`, `aa_alt`, `is_canonical`, `is_mane_select`, `source` |
| `genes` | Gene-level metadata used by builds and exports | `symbol`, `ensembl_id`, `ncbi_id`, `canonical_transcript`, `pli`, `loeuf` |
| `annotation_jobs` | Work queue for deferred annotation handlers | `variant_id`, `source`, `status`, `attempts`, `last_error` |

Alias rows support lookup by multiple external or human-facing names without
duplicating the feature rows. Common `alias_type` values include `hgvs_g`,
`hgvs_c`, `hgvs_p`, `rsid`, `ca_id`, `clinvar_vcv`, `clinvar_rcv`,
`clinvar_allele`, `gnomad_id`, and `myvariant_id`.

All variant-level features are stored in family-specific annotation tables and
linked back to `variants.id` through `variant_id`:

| Table | Feature data | Representative columns |
|---|---|---|
| `annotations_pathogenicity` | AlphaMissense, REVEL, CADD, dbNSFP predictors, and related functional scores | `variant_id`, `predictor`, `predictor_version`, `score`, `rank_score`, `category`, `source` |
| `annotations_population` | gnomAD and other allele frequency datasets | `variant_id`, `dataset`, `pop`, `af`, `ac`, `an`, `n_homozygotes`, `filter_status`, `source` |
| `annotations_clinical` | ClinVar and other clinical assertions | `variant_id`, `source`, `record_id`, `classification`, `review_status`, `stars`, `last_evaluated`, `conditions` |
| `annotations_splice` | SpliceAI, MaxEntScan, AbSplice, dbscSNV, and related splice predictors | `variant_id`, `predictor`, `predictor_version`, `score_type`, `score`, `distance`, `source` |
| `annotations_expression` | pext, AbSplice expression, ABExp, and tissue/transcript expression features | `variant_id`, `metric`, `dataset`, `tissue`, `transcript_id`, `score`, `source` |
| `annotations_structure` | AlphaFold pLDDT, InterPro/domain, and protein-position features | `variant_id`, `feature`, `feature_version`, `protein_accession`, `residue_number`, `score`, `category`, `source` |
| `annotations_conservation` | GERP, PhastCons, PhyloP, SiPhy, and other conservation scores | `variant_id`, `metric`, `score`, `rank_score`, `source` |
| `gene_constraint` | Gene-level pLI, LOEUF, observed/expected, and z-score metrics | `gene_symbol`, `dataset`, `pli`, `lof_z`, `mis_z`, `syn_z`, `oe_lof`, `oe_mis`, `source` |

Raw source datasets are kept outside the normalized SQLite tables under `data/`
when needed, for example `data/alphamissense/` and `data/revel/`. The CLI reads
from those raw files or APIs, then persists queryable results into the SQLite
tables above.

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
│   ├── variants.db          # Local generated SQLite database, ignored by git
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

## Deferred Predictor Commands

```bash
# AlphaFold pLDDT for queued variants
.venv/bin/python -m variantfeatures queue --source alphafold --gene KCNH2 --db data/variants.db
.venv/bin/python -m variantfeatures alphafold-run --db data/variants.db

# VEP plugin config string for SpliceAI / dbscSNV / MaxEntScan / LOFTEE / NMD
.venv/bin/python -m variantfeatures vep-plugin-config \
  --spliceai-snv data/vep_plugins/spliceai_scores.masked.snv.ensembl_mane.grch38.vcf.gz \
  --spliceai-indel data/vep_plugins/spliceai_scores.masked.indel.hg38.vcf.gz \
  --dbscsnv data/vep_plugins/dbscSNV1.1_GRCh38.txt.gz \
  --maxentscan-dir "$HOME/tools/maxentscan/fordownload" \
  --loftee-dir "$HOME/tools/loftee" \
  --loftee-data-dir data/loftee_data \
  --include-nmd

# Large/external predictors exported to CSV/TSV
.venv/bin/python -m variantfeatures pext-import pext_region.tsv --gene KCNH2 --db data/variants.db
.venv/bin/python -m variantfeatures absplice-import absplice.tsv --gene KCNH2 --db data/variants.db
.venv/bin/python -m variantfeatures nmd-import nmdep.tsv --predictor nmdep --gene KCNH2 --db data/variants.db
```

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
.venv/bin/python -m pip install -e .
```

## License

MIT
