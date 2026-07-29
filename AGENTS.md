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
| gnomAD | ✅ Done | Population frequencies via GraphQL, wired through `build` |
| AlphaFold | ✅ Done | Per-residue pLDDT via AlphaFold DB |

Additional handlers/importers live in `variantfeatures/handlers/`: MyVariant/dbNSFP,
gnomAD constraint (pLI/LOEUF), VEP plugins (SpliceAI, dbscSNV, MaxEntScan, LOFTEE,
NMD), ANNOVAR, gnomAD pext, AbSplice/AbExp, and ClinGen.

### Current KCNH2 Coverage
- 22,021 variants with AlphaMissense scores
- 9,522 variants with REVEL scores
- CADD: API-based lookup (on demand)

## Architecture
```
VariantFeatures/
├── variantfeatures/
│   ├── cli.py               # Command-line interface (all subcommands)
│   ├── database.py          # Normalized SQLite schema + operations
│   ├── build.py             # Orchestrates enumerate → queue → run per gene
│   ├── enumerate.py         # Enumerates coding SNVs per transcript
│   ├── identity.py          # Canonical GRCh38 variant identity / normalization
│   ├── transcripts.py       # Transcript/isoform selection (MANE/canonical)
│   ├── uniprot.py           # Gene → UniProt resolution
│   ├── frameshift.py        # Frameshift → stop-gained proxy mapping
│   ├── normalized_export.py # Wide / transcript-wide / long feature export
│   ├── publish.py           # Per-gene SQLite slices + Azure Blob publish
│   ├── worker.py            # Annotation job-queue worker
│   ├── fetchers/            # Low-level source access (used by handlers)
│   │   ├── alphamissense.py revel.py cadd.py clinvar.py gnomad.py lof.py
│   └── handlers/            # Normalized, queue-driven annotation handlers
│       ├── alphamissense.py revel.py cadd.py clinvar.py gnomad.py
│       ├── gnomad_constraint.py myvariant.py alphafold.py vep.py annovar.py
│       ├── absplice.py pext.py clingen_ar.py nmd_rules.py nmd_external.py
│       └── tabular_utils.py
├── scripts/
│   └── full_gene_pipeline.sh # End-to-end per-gene annotation pipeline
├── tests/                    # pytest suite (pytest.ini: testpaths = tests)
├── data/                     # Generated SQLite + cached sources (gitignored)
│   ├── variants.db          # Main SQLite database
│   ├── alphamissense/       # AlphaMissense TSV (~1.1GB)
│   └── revel/               # REVEL data (~6.1GB uncompressed)
└── PIPELINE.md              # End-to-end workflow
```

## Local Data Storage

On Brett's current workstation, the large ignored datasets are backed by the
external APFS volume named `Ezekers`:

| Stable repo path | Physical storage |
|---|---|
| `data/` | `/Volumes/Ezekers/ResearchData/variantFeatures/data` |
| `annovar/humandb/` | `/Volumes/Ezekers/ResearchData/variantFeatures/annovar/humandb` |

Both repo paths are absolute symlinks. Keep code and commands on the stable
repo-relative paths (`data/variants.db`, `data/<source>/`, and
`annovar/humandb/`). Before any build, import, export, or database job, run:

```bash
test -L data && test -d data
test -L annovar/humandb && test -d annovar/humandb
```

If a check fails, mount `Ezekers` at `/Volumes/Ezekers`; do not rename the
volume, replace either broken symlink, or create a fallback local directory.
The symlinks are intentionally local-only and untracked, so a fresh checkout
needs them recreated after both external targets have been verified.

`.gitignore` ignores `data` with an anchored, slashless `/data` rule. That form
matters: a trailing-slash `data/` pattern matches directories only, so a
recreated symlink would show up as untracked in a fresh clone.

### The fail-closed guard

`variantfeatures/local_storage.py` enforces the contract above in code, because
a *missing* `data` link does not fail on its own — `mkdir(parents=True)` would
happily create a real local `data/`, and the job would then build a second
`variants.db` (or re-download a multi-GB TSV) onto the internal disk that no
later job reads. `require_external_storage()` refuses that write instead, and is
called from `VariantDB.__init__` and the AlphaMissense cache path.

The guard only covers the repo-root `data/` tree; a `tmp_path` database or an
explicit `--db` elsewhere passes straight through, so tests are unaffected. On a
machine that genuinely has no external volume (CI, a collaborator's laptop), set
`VARIANTFEATURES_ALLOW_LOCAL_DATA=1` to opt into plain local storage.

### Direct ANNOVAR commands

`scripts/full_gene_pipeline.sh` sets `ANNOVAR_HOME` and `ANNOVAR_DB` together so
the binaries and `humandb` stay in step. Running ANNOVAR outside that script
means setting both explicitly — `ANNOVAR_DB` defaults to `$ANNOVAR_HOME/humandb`,
which is the symlinked path, so leaving it unset while pointing `ANNOVAR_HOME`
somewhere else silently annotates against the wrong database:

```bash
ANNOVAR_HOME="$PWD/annovar" ANNOVAR_DB="$PWD/annovar/humandb" \
  python -m variantfeatures annovar-run --db data/variants.db \
    --build hg38 --protocols refGeneWithVer,clinvar_20240611 --operations g,f
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
