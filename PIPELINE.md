# Variant Interpretation Pipeline

End-to-end workflow for determining what genetic variants mean to the people who carry them.

## Target Scale

**ACMG SF v3.2: 81 genes** for secondary findings reporting.

Initial focus: Cardiac channelopathies (KCNH2, KCNQ1, SCN5A, RYR2)

## Pipeline Stages

```
┌─────────────────────────────────────────────────────────────────────┐
│                        STAGE 1: LITERATURE                          │
│                      (GeneVariantFetcher)                           │
├─────────────────────────────────────────────────────────────────────┤
│  • PubMed search by gene symbol                                     │
│  • Full-text PDF/HTML retrieval                                     │
│  • Supplement extraction (Nature, Elsevier, Springer, Oxford)       │
│  • Variant mention extraction (HGVS parsing)                        │
│  • Individual-level data extraction (phenotype + genotype)          │
│                                                                     │
│  OUTPUT: variants_literature.tsv                                    │
│          individuals_literature.tsv                                 │
└─────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────┐
│                      STAGE 2: ANNOTATION                            │
│                      (VariantFeatures)                              │
├─────────────────────────────────────────────────────────────────────┤
│  MISSENSE VARIANTS:                                                 │
│  • AlphaMissense pathogenicity scores                               │
│  • REVEL ensemble scores                                            │
│  • CADD (PHRED + raw)                                               │
│  • ClinVar classification + review status                           │
│  • gnomAD allele frequencies                                        │
│  • Structural features (domain, AlphaFold proximity)                │
│                                                                     │
│  LOF VARIANTS (nonsense, frameshift, splice):                       │
│  • gnomAD pLI / LOEUF scores (gene-level)                           │
│  • LOFTEE confidence (HC/LC)                                        │
│  • NMD escape prediction                                            │
│  • Truncation position (% protein remaining)                        │
│  • ClinVar classification                                           │
│  • gnomAD allele frequencies                                        │
│                                                                     │
│  OUTPUT: variants_missense.tsv                                      │
│          variants_lof.tsv                                           │
└─────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────┐
│                      STAGE 3: PENETRANCE                            │
│                   (BayesianPenetranceEstimator)                     │
├─────────────────────────────────────────────────────────────────────┤
│  INPUTS:                                                            │
│  • Variant features (from Stage 2)                                  │
│  • Phenotype counts (cases/carriers per variant)                    │
│  • Literature-derived individual data (from Stage 1)                │
│                                                                     │
│  MODELS:                                                            │
│  • Missense model: Feature-weighted Beta-Bernoulli                  │
│  • LOF model: Gene-level hierarchical pooling                       │
│                                                                     │
│  OUTPUT: penetrance_estimates.tsv                                   │
│          - penetrance_mean, penetrance_median                       │
│          - credible intervals (ci_lower, ci_upper)                  │
│          - posterior samples (optional)                             │
└─────────────────────────────────────────────────────────────────────┘
                                    │
                                    ▼
┌─────────────────────────────────────────────────────────────────────┐
│                      STAGE 4: CLINICAL                              │
│                      (VariantBrowser)                               │
├─────────────────────────────────────────────────────────────────────┤
│  • Web interface for clinicians                                     │
│  • Penetrance visualization with uncertainty                        │
│  • Evidence summary (literature, ClinVar, population)               │
│  • Decision support (actionability thresholds)                      │
│                                                                     │
│  URL: variantbrowser.org                                            │
└─────────────────────────────────────────────────────────────────────┘
```

## Data Flow

```
PubMed/Literature ──► GeneVariantFetcher ──┐
                                           │
ClinVar (bulk XML) ────────────────────────┼──► VariantFeatures ──┐
gnomAD v4 (VCF/Hail) ──────────────────────┤                      │
AlphaMissense (TSV) ───────────────────────┤                      │
REVEL (TSV) ───────────────────────────────┤                      │
CADD (TSV) ────────────────────────────────┘                      │
                                                                   │
                                                                   ▼
BioVU / Cohort data ──────────────────────────► BayesianPenetrance
International collaborations (LQT2, etc.) ────►    Estimator
                                                       │
                                                       ▼
                                               VariantBrowser
                                               (Clinical UI)
```

## Variant Type Handling

### Missense Variants
- Single amino acid substitutions
- Use: AlphaMissense, REVEL, CADD, structural features
- Penetrance model: Per-variant feature weighting

### LOF Variants (Loss of Function)
- **Nonsense**: Premature stop codons
- **Frameshift**: Insertions/deletions disrupting reading frame
- **Splice**: Canonical splice site disruptions (+/- 1,2)
- Use: pLI, LOEUF, LOFTEE, NMD prediction
- Penetrance model: Gene-level hierarchical pooling

### Other (future)
- In-frame indels
- Synonymous (splice effects)
- Non-coding regulatory

## Gene Targets

### Phase 1: Cardiac Channelopathies (4 genes)
| Gene | Syndrome | Priority |
|------|----------|----------|
| KCNH2 | LQT2 | Active (breakthrough-events collab) |
| KCNQ1 | LQT1 | High |
| SCN5A | LQT3/Brugada | High |
| RYR2 | CPVT | Medium |

### Phase 2: ACMG SF Expansion (81 genes)
Full list at: https://www.ncbi.nlm.nih.gov/clinvar/docs/acmg/

Categories:
- Cardiomyopathies (TTN, MYH7, MYBPC3, etc.)
- Cancer susceptibility (BRCA1, BRCA2, Lynch genes)
- Familial hypercholesterolemia (LDLR, APOB, PCSK9)
- Metabolic conditions
- Other actionable conditions

## Execution

### Per-Gene Workflow
```bash
# 1. Literature mining
cd GeneVariantFetcher
python -m gvf run --gene KCNH2 --output data/kcnh2/

# 2. Feature annotation
cd VariantFeatures
python -m variantfeatures build --gene KCNH2 --variants data/kcnh2/variants.tsv

# 3. Penetrance estimation
cd BayesianPenetranceEstimator
python -m bpe estimate \
    --features data/kcnh2/features_missense.tsv \
    --phenotypes data/kcnh2/phenotype_counts.tsv \
    --output results/kcnh2/

# 4. Load to browser
cd VariantBrowser
python manage.py import_penetrance results/kcnh2/penetrance.tsv
```

### Batch (All 81 genes)
```bash
nextflow run pipeline.nf --genes acmg_sf_v3.2.txt --output results/
```

## Repository Links

| Repo | Purpose | Status |
|------|---------|--------|
| [GeneVariantFetcher](https://github.com/kroncke-lab/GeneVariantFetcher) | Literature mining | Active |
| [VariantFeatures](https://github.com/kroncke-lab/VariantFeatures) | Feature aggregation | Scaffolded |
| [BayesianPenetranceEstimator](https://github.com/kroncke-lab/BayesianPenetranceEstimator) | Penetrance modeling | Active |
| [VariantBrowser](https://github.com/kroncke-lab/Variant_Browser) | Clinical UI | Production |

## Validation Strategy

1. **Calibration**: Predicted penetrance vs observed event rates
2. **Discrimination**: Separation of pathogenic vs benign
3. **Clinical utility**: Decision curve analysis vs ACMG guidelines
4. **External validation**: Hold-out genes or international cohorts

## R01 Narrative

**Problem**: Variants of uncertain significance (VUS) create clinical uncertainty
**Solution**: Bayesian penetrance estimation with calibrated uncertainty
**Impact**: Reclassify 30%+ of VUS to actionable categories
**Generalizability**: Pipeline scales to all 81 ACMG SF genes
