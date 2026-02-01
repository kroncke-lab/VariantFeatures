# Active Tasks

## Task 1: AlphaMissense Integration

**Objective** → Download AlphaMissense data and implement fetcher to populate `alphamissense_score` and `alphamissense_class` for missense variants.

**Success Measures** →
- AlphaMissense TSV cached locally in `data/`
- `fetch_alphamissense(gene)` returns scored variants (no more NotImplementedError)
- All 4 Phase 1 genes (KCNH2, KCNQ1, SCN5A, RYR2) have AlphaMissense annotations in DB

**Verification**
```bash
# 1. Data file exists (~4GB)
ls -lh data/AlphaMissense*.tsv.gz

# 2. Fetcher works
python -c "from variantfeatures.fetchers.alphamissense import fetch_alphamissense; print(list(fetch_alphamissense('KCNH2'))[:3])"

# 3. DB has scores
sqlite3 data/variants.db "SELECT COUNT(*) FROM variants_missense WHERE gene='KCNH2' AND alphamissense_score IS NOT NULL"
# Expected: >100 variants with scores
```

---

## Task 2: gnomAD v4 Integration

**Objective** → Implement gnomAD fetcher to populate allele frequencies (`gnomad_af`, `gnomad_af_popmax`, `gnomad_homozygotes`).

**Success Measures** →
- `fetch_gnomad(gene)` queries gnomAD GraphQL API successfully
- Rate limiting / pagination handled
- All Phase 1 genes have gnomAD AF annotations for existing variants

**Verification**
```bash
# 1. Fetcher returns data
python -c "from variantfeatures.fetchers.gnomad import fetch_gnomad; v = list(fetch_gnomad('KCNH2'))[:5]; print(len(v), v[0] if v else 'empty')"

# 2. DB has AF values
sqlite3 data/variants.db "SELECT hgvs_p, gnomad_af FROM variants_missense WHERE gene='KCNH2' AND gnomad_af IS NOT NULL LIMIT 5"

# 3. Coverage check
sqlite3 data/variants.db "SELECT COUNT(*) AS with_af, (SELECT COUNT(*) FROM variants_missense WHERE gene='KCNH2') AS total FROM variants_missense WHERE gene='KCNH2' AND gnomad_af IS NOT NULL"
```

---

## Task 3: Fix CLI Build & Query Commands

**Objective** → Make CLI functional: `build` should call fetchers, `query` should work without errors.

**Success Measures** →
- `python -m variantfeatures build --genes KCNH2` loads data from all available sources
- `python -m variantfeatures query --gene KCNH2` outputs variant table (fix missing `get_gene_variants` method)
- CSV/JSON export formats work

**Verification**
```bash
cd /mnt/temp2/kronckbm/gitrepos/VariantFeatures

# 1. Build command runs without error
python -m variantfeatures build --genes KCNH2 2>&1 | head -20

# 2. Query command works
python -m variantfeatures query --gene KCNH2 --format table | head -10

# 3. Export formats
python -m variantfeatures query --gene KCNH2 --format csv > /tmp/kcnh2.csv && wc -l /tmp/kcnh2.csv
python -m variantfeatures query --gene KCNH2 --format json > /tmp/kcnh2.json && python -c "import json; print(len(json.load(open('/tmp/kcnh2.json'))))"
```

---

## Task 4: REVEL Score Integration

**Objective** → Create REVEL fetcher and populate `revel_score` for missense variants.

**Success Measures** →
- REVEL data downloaded/cached (or use per-variant lookup)
- `fetch_revel(gene)` returns scores
- Phase 1 genes have REVEL scores in DB

**Verification**
```bash
# 1. Fetcher exists and works
python -c "from variantfeatures.fetchers.revel import fetch_revel; print(list(fetch_revel('KCNH2'))[:3])"

# 2. DB coverage
sqlite3 data/variants.db "SELECT COUNT(*) AS with_revel FROM variants_missense WHERE gene='KCNH2' AND revel_score IS NOT NULL"
# Expected: >50% of missense variants have REVEL scores
```

---

## Task 5: Phase 1 Gene Validation

**Objective** → Complete feature population for all 4 Phase 1 genes (KCNH2, KCNQ1, SCN5A, RYR2) and validate data quality.

**Success Measures** →
- Each gene has ≥100 missense variants with ClinVar annotations
- ≥50% of variants have AlphaMissense scores
- ≥30% of variants have gnomAD AF data
- Export to CSV works for downstream pipeline

**Verification**
```bash
# Coverage report
sqlite3 data/variants.db "
SELECT 
    gene,
    COUNT(*) as total,
    SUM(CASE WHEN clinvar_id IS NOT NULL THEN 1 ELSE 0 END) as clinvar,
    SUM(CASE WHEN alphamissense_score IS NOT NULL THEN 1 ELSE 0 END) as alphamissense,
    SUM(CASE WHEN gnomad_af IS NOT NULL THEN 1 ELSE 0 END) as gnomad,
    SUM(CASE WHEN revel_score IS NOT NULL THEN 1 ELSE 0 END) as revel
FROM variants_missense 
WHERE gene IN ('KCNH2', 'KCNQ1', 'SCN5A', 'RYR2')
GROUP BY gene
"

# Export for penetrance pipeline
python -m variantfeatures query --gene KCNH2 --format csv > data/exports/kcnh2_features.csv
ls -l data/exports/
```

---

## Backlog

- CADD score integration
- LOF variant pipeline (LOFTEE, NMD prediction)
- Structural features (protein domains, AlphaFold pLDDT)
- Expand to ACMG81 genes
- Nextflow pipeline for batch processing
