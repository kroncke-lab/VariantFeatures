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

# Export one row per variant/transcript consequence for isoform-aware models
.venv/bin/python -m variantfeatures export --gene KCNH2 --layout transcript-wide --output kcnh2_isoforms.csv

# Publish a gene-scoped SQLite slice and manifest for Variant_Browser import
.venv/bin/python -m variantfeatures publish --gene KCNH2 --dry-run

# Map an observed frameshift to finite stop-gained SNV proxies
.venv/bin/python -m variantfeatures frameshift-map --gene KCNH2 --aa-position 1155 --queue-source cadd

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
| **ClinVar** | Clinical classifications + review status from local ClinVar summary and MyVariant fallback |
| **gnomAD** | Population allele frequencies |

### LoF / Splice / Structure
| Source | Description |
|--------|-------------|
| **AlphaFold DB** | Per-residue pLDDT at missense or truncation position |
| **Frameshift proxy map** | Maps frameshift amino-acid positions to nearby stop-gained SNV proxies |
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

## Required Fully Populated Genes

The current must-have target set is:

```text
APOE, MYBPC3, BRCA1, BRCA2, KCNH2, LDLR
```

`MYBPC3` is the canonical HGNC symbol for the cardiomyopathy gene sometimes
mistyped as `MYPBC3`. These targets are additive: keep the existing populated
genes in `data/variants.db` and fill any missing coverage for this set.

```bash
.venv/bin/python -m variantfeatures build \
  --genes APOE,MYBPC3,BRCA1,BRCA2,KCNH2,LDLR \
  --db data/variants.db
```

## Target Architecture

```bash
# Normalized canonical build for any gene
variantfeatures build --gene BRCA1

# Include alternate coding isoforms when isoform-specific effects matter
variantfeatures build --gene BRCA1 --isoforms all

# This:
# 1. Looks up the selected transcript set (canonical by default; MANE/all optional)
#    and enumerates coding SNVs per isoform
# 2. Queues/runs normalized sources such as MyVariant/dbNSFP, gnomAD,
#    AlphaMissense, REVEL, AlphaFold, pext, and NMD rules
# 3. Stores features in family-specific normalized tables
# 4. Exports one-row-per-variant, one-row-per-transcript, or long audit CSVs
```

Isoform-aware builds use one canonical GRCh38 allele row in `variants` and one
row per affected transcript in `variant_consequences`. Transcript metadata is
stored in `transcripts`. Transcript-scoped pext/expression rows stay keyed by
`transcript_id`, which lets downstream models prefer biologically relevant
isoforms and down-weight effects seen only in weakly expressed isoforms.

## Database Schema

The normalized path uses `variants` for canonical GRCh38 alleles,
`variant_consequences` for transcript effects, and feature-family tables:

| Feature family | Table | Export prefix |
|---|---|---|
| Pathogenicity / functional scores | `annotations_pathogenicity` | `pathogenicity.<predictor>...` |
| Population frequencies | `annotations_population` | `population.<dataset>.<pop>...` |
| Clinical assertions | `annotations_clinical` | `clinical.<source>...` |
| Splice predictors | `annotations_splice` | `splice.<predictor>...` |
| pext / expression features | `annotations_expression` | `expression.<metric>.<dataset>.<tissue>[.<transcript>]...` |
| Structure / domains | `annotations_structure` | `structure.<feature>...` |
| Conservation | `annotations_conservation` | `conservation.<metric>...` |
| Gene constraint | `gene_constraint` | `gene_constraint.<dataset>...` |

The legacy `variants_missense` / `variants_lof` tables have been removed — all
data now lives in the normalized `variants` + `variant_consequences` +
`annotations_*` schema above. There is no longer a `--legacy` mode on `query`,
`stats`, or `export`.

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
| `transcripts` | Gene isoforms selected during build/enumeration | `transcript_id`, `gene_symbol`, `refseq_match`, `protein_id`, `biotype`, `cds_length`, `is_canonical`, `is_mane_select`, `is_mane_plus_clinical`, `transcript_support_level`, `appris` |
| `variant_consequences` | Per-transcript gene/cDNA/protein effects | `variant_id`, `transcript_id`, `gene_symbol`, `gene_ensembl`, `consequence`, `hgvs_c`, `hgvs_p`, `aa_pos`, `aa_ref`, `aa_alt`, `is_canonical`, `is_mane_select`, `source` |
| `frameshift_nonsense_mappings` | Frameshift amino-acid positions mapped to finite stop-gained SNV proxy variants | `gene_symbol`, `transcript_id`, `frameshift_aa_pos`, `direction`, `n_steps_wrt_frameshift`, `proxy_variant_id` |
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

## Publishing to Azure Blob

Publishing is separate from `build` and is always opt-in. It creates per-gene
SQLite slices under `dist/publish/`, writes a provenance manifest, and then
uploads to Azure Blob when configured:

```bash
.venv/bin/python -m pip install -e '.[azure]'

export VF_BLOB_ACCOUNT_URL="https://<account>.blob.core.windows.net"
export VF_BLOB_CONTAINER="variantfeatures"

.venv/bin/python -m variantfeatures publish \
  --gene KCNH2 \
  --prefix pipeline/variantfeatures
```

Blob layout:

```text
pipeline/variantfeatures/{YYYYMMDD-HHMM}__{gitsha7}/genes/KCNH2.db
pipeline/variantfeatures/{YYYYMMDD-HHMM}__{gitsha7}/manifest.json
pipeline/variantfeatures/latest.json
```

`latest.json` is a per-gene pointer used by Variant_Browser. Each entry includes
the current slice path, SHA-256, build timestamp, and schema version. The pointer
is updated with an Azure Blob ETag write after versioned artifacts upload.

Authentication uses `DefaultAzureCredential`; no account keys or secrets are
stored in this repo. Run `az login` for local publishing, or use managed
identity/workload identity in cloud jobs. If `VF_BLOB_ACCOUNT_URL` or
`VF_BLOB_CONTAINER` is unset, `publish` automatically behaves like a dry run and
prints the blob paths, hashes, and `latest.json` diff without contacting Azure.

## Project Structure

```
VariantFeatures/
├── variantfeatures/
│   ├── cli.py              # Command-line interface
│   ├── database.py         # SQLite operations
│   ├── frameshift.py       # Frameshift-to-nonsense proxy mapping
│   └── fetchers/
│       ├── alphamissense.py # ~1.1GB TSV
│       ├── revel.py         # ~6.1GB TSV
│       ├── cadd.py          # REST API
│       ├── clinvar.py       # Bulk XML
│       ├── gnomad.py        # GraphQL API
│       └── lof.py           # LOFTEE annotations
├── scripts/
│   └── full_gene_pipeline.sh # End-to-end per-gene annotation pipeline
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
# Frameshift-to-nonsense proxy mapping. Proxy variants are ordinary stop_gained
# SNVs, so queued sources such as CADD attach through annotations_pathogenicity.
.venv/bin/python -m variantfeatures frameshift-map \
  --gene KCNH2 \
  --aa-position 1155 \
  --n-steps 20 \
  --queue-source cadd \
  --db data/variants.db

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
.venv/bin/python -m variantfeatures clinvar-import variant_summary.txt.gz --gene KCNH2 --db data/variants.db
.venv/bin/python -m variantfeatures pext-import pext_region.tsv --gene KCNH2 --db data/variants.db
.venv/bin/python -m variantfeatures absplice-import absplice.tsv --gene KCNH2 --db data/variants.db
.venv/bin/python -m variantfeatures dbscsnv-import hg38_dbscsnv11.txt --gene KCNH2 --db data/variants.db
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
