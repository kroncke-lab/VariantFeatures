# Active Tasks

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
sqlite3 data/variants.db "SELECT COUNT(*) FROM variants_missense WHERE gene='KCNH2' AND alphamissense_score IS NOT NULL"
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
sqlite3 data/variants.db "SELECT hgvs_p, gnomad_af FROM variants_missense WHERE gene='KCNH2' AND gnomad_af IS NOT NULL LIMIT 5"
```

---

## ✅ Task 3: Fix CLI Build & Query Commands — COMPLETE

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

---

## Task 4: REVEL Score Integration

**Status:** Not started

**Objective** → Create REVEL fetcher and populate `revel_score` for missense variants.

**Success Measures** →
- REVEL data downloaded/cached (or use per-variant lookup)
- `fetch_revel(gene)` returns scores
- Phase 1 genes have REVEL scores in DB

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

- CADD score integration
- LOF variant pipeline (LOFTEE, NMD prediction)
- Structural features (protein domains, AlphaFold pLDDT)
- Expand to ACMG81 genes
- Nextflow pipeline for batch processing
