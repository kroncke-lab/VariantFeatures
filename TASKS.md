# Active Tasks

## 🎯 TARGET ARCHITECTURE (per Brett, 2026-02-18)

**Goal:** Single command to build complete gene database
```bash
variantfeatures build --gene KCNH2
```

**Should automatically:**
1. Look up gene → get canonical transcript(s), UniProt ID, genomic coords
2. Download/extract AlphaMissense scores for that gene
3. Download/extract REVEL scores for that gene  
4. Fetch gnomAD allele frequencies via API
5. Download AlphaFold structure (AF2 or newer public version)
6. Build SQLite database with all features linked

**Data sources to integrate:**
- [x] AlphaMissense (have 1.1GB file)
- [x] REVEL (local normalized batch handler)
- [x] gnomAD v4 API (normalized handler wired through build)
- [x] AlphaFold DB (https://alphafold.ebi.ac.uk/) — pLDDT fetcher/handler implemented
- [x] Gene metadata lookup (Ensembl canonical/MANE transcript enumeration)

**Current state:** Normalized build/export layer is implemented. `build` now
enumerates canonical variants into `variants` + `variant_consequences`, queues
family-specific normalized annotations, and can run available annotators. `export`
now defaults to normalized wide/long output with explicit feature-family prefixes.

---

## ✅ Task 1: AlphaMissense Integration — FETCHER IMPLEMENTED

**Status:** Fetcher implemented, needs data file download (~4GB)

**Objective** → Download AlphaMissense data and implement fetcher to populate `alphamissense_score` and `alphamissense_class` for missense variants.

**What's done:**
- `fetch_alphamissense(gene)` fully implemented in `variantfeatures/fetchers/alphamissense.py`
- Parses TSV file and yields variant dicts with scores
- Gene-to-UniProt mapping for Phase 1 genes (KCNH2, KCNQ1, SCN5A, RYR2)
- CLI `build` command wired up to call fetcher

**What's needed:**
- Download AlphaMissense data file (~4GB): `AlphaMissense_aa_substitutions.tsv.gz`
- Cache location: `data/alphamissense/`
- Run: `/usr/bin/python3 run_cli.py build -g KCNH2 --sources alphamissense`

**Verification**
```bash
# 1. Data file exists (~4GB)
ls -lh data/alphamissense/AlphaMissense_aa_substitutions.tsv.gz

# 2. Fetcher works
/usr/bin/python3 -c "from variantfeatures.fetchers.alphamissense import fetch_alphamissense; print(list(fetch_alphamissense('KCNH2'))[:3])"

# 3. DB has scores
sqlite3 data/variants.db "SELECT COUNT(DISTINCT p.variant_id) FROM annotations_pathogenicity p JOIN variant_consequences c ON c.variant_id=p.variant_id WHERE c.gene_symbol='KCNH2' AND p.predictor='alphamissense'"
```

---

## ✅ Task 2: gnomAD v4 Integration — FETCHER IMPLEMENTED

**Status:** Fetcher implemented, needs network access to gnomAD API

**Objective** → Implement gnomAD fetcher to populate allele frequencies (`gnomad_af`, `gnomad_af_popmax`, `gnomad_homozygotes`).

**What's done:**
- `fetch_gnomad(gene)` fully implemented in `variantfeatures/fetchers/gnomad.py`
- GraphQL query for gene variants
- Rate limiting (0.5s delay between requests)
- Single-variant lookup function `fetch_single_variant()`
- CLI `build` command wired up to call fetcher

**What's needed:**
- Network access to gnomAD API (`https://gnomad.broadinstitute.org/api`)
- Run: `/usr/bin/python3 run_cli.py build -g KCNH2 --sources gnomad`

**Verification**
```bash
# 1. Fetcher returns data
/usr/bin/python3 -c "from variantfeatures.fetchers.gnomad import fetch_gnomad; import itertools; print(list(itertools.islice(fetch_gnomad('KCNH2'), 5)))"

# 2. DB has AF values
sqlite3 data/variants.db "SELECT c.hgvs_p, pop.af FROM annotations_population pop JOIN variant_consequences c ON c.variant_id=pop.variant_id WHERE c.gene_symbol='KCNH2' AND pop.pop='all' LIMIT 5"
```

---

## ✅ Task 4: REVEL Score Integration

**Status:** ✅ DONE

**Objective** → Create REVEL fetcher and populate `revel_score` for missense variants.

**Success Measures** →
- REVEL data downloaded/cached (or use per-variant lookup) ✅
- `revel-run` drains queued normalized jobs from local REVEL file/zip ✅
- Scores stored in `annotations_pathogenicity` as `predictor='revel'` ✅

**Notes:**
- REVEL scores available from: https://sites.google.com/site/revelgenomics/
- ~3GB compressed TSV file
- Similar implementation pattern to AlphaMissense fetcher

---

## Task 5: Phase 1 Gene Validation

**Status:** Partially complete — ClinVar loaded, awaiting AlphaMissense/gnomAD data

**Current Coverage:**
```
Gene       Total    ClinVar    AlphaMissense   gnomAD     REVEL     
----------------------------------------------------------------------
KCNH2      699      699        0               0          0         
SCN5A      654      654        0               0          0         
KCNQ1      394      394        0               0          0         
RYR2       376      376        0               0          0         
```

**What's done:**
- ClinVar data loaded for all Phase 1 genes ✅
- Export pipeline working ✅ (e.g., `data/exports/kcnh2_features.csv`)

**What's needed:**
- Download and load AlphaMissense data
- Run gnomAD fetcher (needs network)
- Implement and run REVEL fetcher

---

## Backlog

- CADD score integration hardening (direct normalized API handler exists; MyVariant/dbNSFP remains preferred for bulk)
- LOF variant pipeline hardening (LOFTEE install/data checks; NMD-rule implemented, VEP NMD parsed)
- Structural features beyond AlphaFold pLDDT (protein domains, disorder, solvent accessibility)
- Expand to ACMG81 genes
- Nextflow pipeline for batch processing

## Normalized Build / Export Added

**Status:** Implemented.

- `variantfeatures build --gene KCNH2` is now the normalized canonical build path.
- `--sources` accepts source groups such as `core`, `all`, `population`,
  `splice`, `expression`, `structure`, `pathogenicity`, `pext`, `revel`, and `cadd`.
- `--no-run` supports deterministic enumerate/queue-only builds.
- `variantfeatures export` now defaults to normalized output.
- `--layout wide` emits one row per variant with namespaced columns.
- `--layout long` emits one row per feature field for auditing feature provenance.
- `variantfeatures feature-schema` documents where pext, splice, population,
  clinical, structure, conservation, pathogenicity, and gene-constraint features live.

## Deferred Predictor Wiring Added

**Status:** Implemented as first-class handlers/importers where external data remains large or model-specific.

- `alphafold-run`: fetches AlphaFold DB confidence JSON and stores per-variant `alphafold_plddt` in `annotations_structure`.
- VEP parser now persists SpliceAI overall/directional scores, dbscSNV ADA/RF, MaxEntScan scores, and VEP NMD plugin output.
- `pext-import` / `pext-run`: imports local gnomAD pext CSV/TSV exports into `annotations_expression`.
- `absplice-import`: imports AbSplice_DNA/RNA into `annotations_splice` and AbExp columns into `annotations_expression`.
- `nmd-import`: imports coordinate-keyed NMDEP/NMDetective-style outputs into `annotations_pathogenicity`.
- `vep-plugin-config`: prints a semicolon-separated `VEP_PLUGINS` string for SpliceAI, dbscSNV, MaxEntScan, LOFTEE, and NMD.

---

## Completed Tasks ✅

### Task 3: Fix CLI Build & Query Commands — COMPLETE

**Status:** ✅ DONE

**Objective** → Make CLI functional: `build` should call fetchers, `query` should work without errors.

**Completed:**
- Fixed CLI `query` command (was calling nonexistent `get_gene_variants`, now uses `get_gene_missense`)
- Implemented CLI `build` command to call AlphaMissense and gnomAD fetchers
- Added `stats` command for database coverage reporting
- Added `export` command for CSV export
- Added `__main__.py` for `python -m variantfeatures` support
- Created `run_cli.py` wrapper for easier invocation

**CLI Commands:**
```bash
# Stats
/usr/bin/python3 run_cli.py stats

# Query (table, csv, json)
/usr/bin/python3 run_cli.py query -g KCNH2
/usr/bin/python3 run_cli.py query -g KCNH2 --format csv
/usr/bin/python3 run_cli.py query -g KCNH2 --format json

# Export to file
/usr/bin/python3 run_cli.py export -g KCNH2 -o data/exports/kcnh2.csv

# Build (when data available)
/usr/bin/python3 run_cli.py build -g KCNH2 --sources alphamissense,gnomad
```
