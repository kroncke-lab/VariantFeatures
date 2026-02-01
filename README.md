# VariantFeatures

Aggregate predictive features for genetic variants from multiple sources into a unified SQLite database.

## Data Sources

| Source | Features | URL |
|--------|----------|-----|
| **AlphaMissense** | Pathogenicity scores from DeepMind | [alphafold.ebi.ac.uk](https://alphafold.ebi.ac.uk/download) |
| **REVEL** | Ensemble method for missense variants | [sites.google.com/site/revelgenomics](https://sites.google.com/site/revelgenomics/) |
| **CADD** | Combined Annotation Dependent Depletion | [cadd.gs.washington.edu](https://cadd.gs.washington.edu/) |
| **ClinVar** | Clinical significance annotations | [ncbi.nlm.nih.gov/clinvar](https://www.ncbi.nlm.nih.gov/clinvar/) |
| **gnomAD** | Population allele frequencies | [gnomad.broadinstitute.org](https://gnomad.broadinstitute.org/) |

## Schema

```sql
CREATE TABLE variants (
    id INTEGER PRIMARY KEY,
    gene TEXT NOT NULL,
    hgvs_p TEXT,           -- e.g., p.Arg528His
    hgvs_c TEXT,           -- e.g., c.1583G>A
    chromosome TEXT,
    position INTEGER,
    ref TEXT,
    alt TEXT,
    
    -- Scores
    alphamissense_score REAL,
    alphamissense_class TEXT,  -- likely_benign, ambiguous, likely_pathogenic
    revel_score REAL,
    cadd_phred REAL,
    cadd_raw REAL,
    
    -- ClinVar
    clinvar_id INTEGER,
    clinvar_significance TEXT,
    clinvar_review_status TEXT,
    clinvar_last_evaluated DATE,
    
    -- gnomAD
    gnomad_af REAL,        -- global allele frequency
    gnomad_af_popmax REAL, -- max population AF
    gnomad_homozygotes INTEGER,
    
    UNIQUE(gene, hgvs_p)
);

CREATE INDEX idx_gene ON variants(gene);
CREATE INDEX idx_clinvar ON variants(clinvar_significance);
CREATE INDEX idx_gnomad_af ON variants(gnomad_af);
```

## Usage

```bash
# Fetch and build database
python -m variantfeatures build --genes KCNH2,KCNQ1,SCN5A

# Query
python -m variantfeatures query --gene KCNH2 --format csv > kcnh2_features.csv
```

## Project Structure

```
variantfeatures/
├── __init__.py
├── cli.py              # Command-line interface
├── database.py         # SQLite operations
├── fetchers/
│   ├── __init__.py
│   ├── alphamissense.py
│   ├── revel.py
│   ├── cadd.py
│   ├── clinvar.py
│   └── gnomad.py
└── utils.py
```

## Related Projects

- [GeneVariantFetcher](https://github.com/kroncke-lab/GeneVariantFetcher) — Literature mining for variants
- [BayesianPenetranceEstimator](https://github.com/kroncke-lab/BayesianPenetranceEstimator) — Penetrance modeling

## License

MIT
