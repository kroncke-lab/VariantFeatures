# VariantFeatures

Aggregate predictive features for genetic variants from multiple sources into a unified SQLite database.

**Part of the Variant Interpretation Pipeline** — see [PIPELINE.md](PIPELINE.md) for end-to-end workflow.

## Target Scale

- **Phase 1**: Cardiac channelopathies (KCNH2, KCNQ1, SCN5A, RYR2)
- **Phase 2**: ACMG SF v3.2 — 81 genes for secondary findings

## Variant Types

### Missense Variants
| Source | Features |
|--------|----------|
| **AlphaMissense** | Pathogenicity score, class |
| **REVEL** | Ensemble missense score |
| **CADD** | PHRED + raw scores |
| **ClinVar** | Classification, review status, stars |
| **gnomAD** | AF, AF_popmax, homozygotes |

### LOF Variants (Nonsense, Frameshift, Splice)
| Source | Features |
|--------|----------|
| **gnomAD** | pLI, LOEUF (gene-level constraint) |
| **LOFTEE** | Confidence (HC/LC), flags |
| **NMD** | Escape prediction |
| **ClinVar** | Classification, review status |
| **Position** | Truncation position (% protein remaining) |

## Schema

```sql
-- Missense variants
CREATE TABLE variants_missense (
    gene TEXT NOT NULL,
    hgvs_p TEXT,                 -- p.Arg528His
    hgvs_c TEXT,                 -- c.1583G>A
    chromosome TEXT,
    position INTEGER,
    ref TEXT, alt TEXT,
    genome_build TEXT DEFAULT 'GRCh38',
    transcript_id TEXT,          -- NM_ accession
    
    -- Scores
    alphamissense_score REAL,
    alphamissense_class TEXT,
    revel_score REAL,
    cadd_phred REAL,
    
    -- ClinVar
    clinvar_significance TEXT,
    clinvar_stars INTEGER,
    
    -- gnomAD
    gnomad_af REAL,
    gnomad_af_popmax REAL,
    
    UNIQUE(gene, hgvs_p),
    UNIQUE(chromosome, position, ref, alt, genome_build)
);

-- LOF variants
CREATE TABLE variants_lof (
    gene TEXT NOT NULL,
    hgvs_p TEXT,                 -- p.Arg528Ter
    hgvs_c TEXT,
    lof_type TEXT NOT NULL,      -- nonsense, frameshift, splice_donor, splice_acceptor
    
    -- LOF-specific
    loftee_confidence TEXT,      -- HC, LC
    nmd_escape INTEGER,
    truncation_position REAL,    -- Fraction of protein remaining
    gene_pli REAL,
    gene_loeuf REAL,
    
    -- ClinVar + gnomAD (same as missense)
    ...
    
    UNIQUE(gene, hgvs_c)
);

-- Gene-level annotations
CREATE TABLE genes (
    symbol TEXT NOT NULL UNIQUE,
    pli REAL,
    loeuf REAL,
    canonical_transcript TEXT
);

-- Penetrance estimates (from BayesianPenetranceEstimator)
CREATE TABLE penetrance_estimates (
    variant_type TEXT,           -- missense, lof
    variant_id INTEGER,
    gene TEXT,
    penetrance_mean REAL,
    penetrance_median REAL,
    ci_lower REAL,
    ci_upper REAL
);
```

## Usage

```bash
# Build database for a gene
python -m variantfeatures build --gene KCNH2

# Query missense variants
python -m variantfeatures query --gene KCNH2 --type missense --format csv

# Query LOF variants
python -m variantfeatures query --gene KCNH2 --type lof --format csv
```

## Data Sources

| Source | Method | Size |
|--------|--------|------|
| ClinVar | Bulk XML (monthly) | ~1GB |
| gnomAD v4 | Hail MatrixTable or VCF | ~500GB |
| AlphaMissense | TSV from GCS | ~4GB |
| REVEL | TSV | ~10GB |
| CADD | TSV (per-chromosome) | ~80GB |

## Project Structure

```
variantfeatures/
├── __init__.py
├── cli.py              # Command-line interface
├── database.py         # SQLite operations (missense + LOF)
├── fetchers/
│   ├── alphamissense.py
│   ├── clinvar.py
│   ├── gnomad.py
│   └── lof.py          # LOFTEE, NMD prediction
└── utils.py
```

## Related Projects

| Repo | Purpose |
|------|---------|
| [GeneVariantFetcher](https://github.com/kroncke-lab/GeneVariantFetcher) | Literature mining |
| [BayesianPenetranceEstimator](https://github.com/kroncke-lab/BayesianPenetranceEstimator) | Penetrance modeling |
| [Variant_Browser](https://github.com/kroncke-lab/Variant_Browser) | Clinical UI |

## License

MIT
