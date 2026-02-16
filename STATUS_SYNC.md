# VariantFeatures Sync Status

**Last checked:** 2026-02-05 16:30 CST

## Git Status

| Check | Status |
|-------|--------|
| Uncommitted changes | ✅ None |
| Local ahead of remote | ✅ Pushed |
| Pushed commit | `bf4dbb8` - AlphaMissense parser header detection fix |

## Pending Tasks (from TASKS.md)

### Highest Priority — Data Loading

1. **AlphaMissense** (Task 1) — Fetcher DONE, needs ~4GB data file
   - Download: `AlphaMissense_aa_substitutions.tsv.gz`
   - Location: `data/alphamissense/`

2. **gnomAD v4** (Task 2) — Fetcher DONE, needs network access
   - API: `https://gnomad.broadinstitute.org/api`

3. **REVEL Integration** (Task 4) — Not started
   - ~3GB data file needed

### Current Database Coverage

| Gene | Total | ClinVar | AlphaMissense | gnomAD | REVEL |
|------|-------|---------|---------------|--------|-------|
| KCNH2 | 699 | 699 | 0 | 0 | 0 |
| SCN5A | 654 | 654 | 0 | 0 | 0 |
| KCNQ1 | 394 | 394 | 0 | 0 | 0 |
| RYR2 | 376 | 376 | 0 | 0 | 0 |

## What's Done ✅

- CLI commands working (stats, query, export, build)
- ClinVar fetcher implemented and loaded
- AlphaMissense fetcher implemented
- gnomAD fetcher implemented
- Export pipeline working

## Next Steps

1. Download AlphaMissense data file (~4GB)
2. Run AlphaMissense build: `/usr/bin/python3 run_cli.py build -g KCNH2 --sources alphamissense`
3. Run gnomAD build (when network available)
4. Implement REVEL fetcher

## Backlog

- CADD score integration
- LOF variant pipeline
- Structural features (AlphaFold pLDDT)
- Expand to ACMG81 genes
