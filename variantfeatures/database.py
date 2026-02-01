"""SQLite database operations for variant features."""

import sqlite3
from pathlib import Path
from typing import Optional, Literal

DEFAULT_DB = Path(__file__).parent.parent / "data" / "variants.db"

SCHEMA = """
-- Gene-level annotations (pLI, LOEUF for LOF interpretation)
CREATE TABLE IF NOT EXISTS genes (
    id INTEGER PRIMARY KEY,
    symbol TEXT NOT NULL UNIQUE,
    
    -- gnomAD constraint metrics
    pli REAL,                    -- Probability of LoF intolerance
    loeuf REAL,                  -- LoF observed/expected upper bound
    loeuf_lower REAL,
    loeuf_upper REAL,
    
    -- Gene info
    ensembl_id TEXT,
    ncbi_id TEXT,
    canonical_transcript TEXT,   -- NM_ accession
    
    -- Metadata
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

-- Missense variants
CREATE TABLE IF NOT EXISTS variants_missense (
    id INTEGER PRIMARY KEY,
    
    -- Identifiers
    gene TEXT NOT NULL,
    hgvs_p TEXT,                 -- e.g., p.Arg528His
    hgvs_c TEXT,                 -- e.g., c.1583G>A
    
    -- Genomic coordinates (VCF-style, normalized)
    chromosome TEXT,
    position INTEGER,
    ref TEXT,
    alt TEXT,
    genome_build TEXT DEFAULT 'GRCh38',
    
    -- Transcript
    transcript_id TEXT,          -- NM_ accession
    
    -- Pathogenicity scores
    alphamissense_score REAL,
    alphamissense_class TEXT,    -- likely_benign, ambiguous, likely_pathogenic
    revel_score REAL,
    cadd_phred REAL,
    cadd_raw REAL,
    
    -- Structural features
    domain TEXT,                 -- Protein domain
    alphafold_plddt REAL,        -- AlphaFold confidence at position
    
    -- ClinVar
    clinvar_id INTEGER,
    clinvar_significance TEXT,
    clinvar_review_status TEXT,
    clinvar_stars INTEGER,
    clinvar_last_evaluated DATE,
    
    -- gnomAD
    gnomad_af REAL,              -- Global allele frequency
    gnomad_af_popmax REAL,       -- Max population AF
    gnomad_homozygotes INTEGER,
    gnomad_version TEXT,
    
    -- Provenance
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    
    UNIQUE(gene, hgvs_p),
    UNIQUE(chromosome, position, ref, alt, genome_build)
);

-- Loss-of-function variants (nonsense, frameshift, splice)
CREATE TABLE IF NOT EXISTS variants_lof (
    id INTEGER PRIMARY KEY,
    
    -- Identifiers
    gene TEXT NOT NULL,
    hgvs_p TEXT,                 -- e.g., p.Arg528Ter, p.Gly123fs
    hgvs_c TEXT,
    
    -- Variant type
    lof_type TEXT NOT NULL,      -- nonsense, frameshift, splice_donor, splice_acceptor
    
    -- Genomic coordinates
    chromosome TEXT,
    position INTEGER,
    ref TEXT,
    alt TEXT,
    genome_build TEXT DEFAULT 'GRCh38',
    
    -- Transcript
    transcript_id TEXT,
    
    -- LOF-specific annotations
    loftee_confidence TEXT,      -- HC (high confidence), LC (low confidence)
    loftee_flags TEXT,           -- Comma-separated flags
    nmd_escape INTEGER,          -- 1 if likely to escape NMD, 0 otherwise
    truncation_position REAL,    -- Fraction of protein remaining (0-1)
    last_exon INTEGER,           -- 1 if in last exon (NMD escape)
    
    -- Gene-level constraint (denormalized for convenience)
    gene_pli REAL,
    gene_loeuf REAL,
    
    -- ClinVar
    clinvar_id INTEGER,
    clinvar_significance TEXT,
    clinvar_review_status TEXT,
    clinvar_stars INTEGER,
    clinvar_last_evaluated DATE,
    
    -- gnomAD
    gnomad_af REAL,
    gnomad_af_popmax REAL,
    gnomad_homozygotes INTEGER,
    gnomad_version TEXT,
    
    -- Provenance
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    
    UNIQUE(gene, hgvs_c),
    UNIQUE(chromosome, position, ref, alt, genome_build)
);

-- Penetrance estimates (populated by BayesianPenetranceEstimator)
CREATE TABLE IF NOT EXISTS penetrance_estimates (
    id INTEGER PRIMARY KEY,
    
    -- Link to variant
    variant_type TEXT NOT NULL,  -- missense, lof
    variant_id INTEGER NOT NULL, -- FK to variants_missense or variants_lof
    gene TEXT NOT NULL,
    hgvs_p TEXT,
    
    -- Estimates
    penetrance_mean REAL,
    penetrance_median REAL,
    ci_lower REAL,               -- 95% credible interval
    ci_upper REAL,
    
    -- Model metadata
    model_version TEXT,
    n_cases INTEGER,
    n_carriers INTEGER,
    
    -- Provenance
    estimated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    
    UNIQUE(variant_type, variant_id)
);

-- Indexes
CREATE INDEX IF NOT EXISTS idx_missense_gene ON variants_missense(gene);
CREATE INDEX IF NOT EXISTS idx_missense_clinvar ON variants_missense(clinvar_significance);
CREATE INDEX IF NOT EXISTS idx_missense_coords ON variants_missense(chromosome, position);
CREATE INDEX IF NOT EXISTS idx_lof_gene ON variants_lof(gene);
CREATE INDEX IF NOT EXISTS idx_lof_type ON variants_lof(lof_type);
CREATE INDEX IF NOT EXISTS idx_lof_clinvar ON variants_lof(clinvar_significance);
CREATE INDEX IF NOT EXISTS idx_penetrance_gene ON penetrance_estimates(gene);
"""


class VariantDB:
    """SQLite database for variant features."""
    
    def __init__(self, db_path: Optional[Path] = None):
        self.db_path = db_path or DEFAULT_DB
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        self.conn = sqlite3.connect(self.db_path)
        self.conn.row_factory = sqlite3.Row
        self._init_schema()
    
    def _init_schema(self):
        self.conn.executescript(SCHEMA)
        self.conn.commit()
    
    def upsert_missense(self, gene: str, hgvs_p: str, **features):
        """Insert or update a missense variant."""
        return self._upsert("variants_missense", gene, hgvs_p, **features)
    
    def upsert_lof(self, gene: str, hgvs_c: str, lof_type: str, **features):
        """Insert or update a LOF variant."""
        features["lof_type"] = lof_type
        columns = ["gene", "hgvs_c"] + list(features.keys())
        placeholders = ", ".join(["?"] * len(columns))
        updates = ", ".join(f"{k} = excluded.{k}" for k in features.keys())
        updates += ", updated_at = CURRENT_TIMESTAMP"
        
        sql = f"""
        INSERT INTO variants_lof ({", ".join(columns)})
        VALUES ({placeholders})
        ON CONFLICT(gene, hgvs_c) DO UPDATE SET {updates}
        """
        self.conn.execute(sql, [gene, hgvs_c] + list(features.values()))
        self.conn.commit()
    
    def upsert_gene(self, symbol: str, **features):
        """Insert or update gene-level annotations."""
        columns = ["symbol"] + list(features.keys())
        placeholders = ", ".join(["?"] * len(columns))
        updates = ", ".join(f"{k} = excluded.{k}" for k in features.keys())
        updates += ", updated_at = CURRENT_TIMESTAMP"
        
        sql = f"""
        INSERT INTO genes ({", ".join(columns)})
        VALUES ({placeholders})
        ON CONFLICT(symbol) DO UPDATE SET {updates}
        """
        self.conn.execute(sql, [symbol] + list(features.values()))
        self.conn.commit()
    
    def _upsert(self, table: str, gene: str, hgvs_p: str, **features):
        """Generic upsert for variant tables."""
        columns = ["gene", "hgvs_p"] + list(features.keys())
        placeholders = ", ".join(["?"] * len(columns))
        updates = ", ".join(f"{k} = excluded.{k}" for k in features.keys())
        updates += ", updated_at = CURRENT_TIMESTAMP"
        
        sql = f"""
        INSERT INTO {table} ({", ".join(columns)})
        VALUES ({placeholders})
        ON CONFLICT(gene, hgvs_p) DO UPDATE SET {updates}
        """
        self.conn.execute(sql, [gene, hgvs_p] + list(features.values()))
        self.conn.commit()
    
    def get_missense(self, gene: str, hgvs_p: str) -> Optional[dict]:
        """Get a missense variant."""
        cur = self.conn.execute(
            "SELECT * FROM variants_missense WHERE gene = ? AND hgvs_p = ?",
            [gene, hgvs_p]
        )
        row = cur.fetchone()
        return dict(row) if row else None
    
    def get_lof(self, gene: str, hgvs_c: str) -> Optional[dict]:
        """Get a LOF variant."""
        cur = self.conn.execute(
            "SELECT * FROM variants_lof WHERE gene = ? AND hgvs_c = ?",
            [gene, hgvs_c]
        )
        row = cur.fetchone()
        return dict(row) if row else None
    
    def get_gene_missense(self, gene: str) -> list[dict]:
        """Get all missense variants for a gene."""
        cur = self.conn.execute(
            "SELECT * FROM variants_missense WHERE gene = ? ORDER BY hgvs_p",
            [gene]
        )
        return [dict(row) for row in cur.fetchall()]
    
    def get_gene_lof(self, gene: str) -> list[dict]:
        """Get all LOF variants for a gene."""
        cur = self.conn.execute(
            "SELECT * FROM variants_lof WHERE gene = ? ORDER BY position",
            [gene]
        )
        return [dict(row) for row in cur.fetchall()]
    
    def get_gene_all(self, gene: str) -> dict:
        """Get all variants for a gene, separated by type."""
        return {
            "missense": self.get_gene_missense(gene),
            "lof": self.get_gene_lof(gene),
        }
    
    def store_penetrance(self, variant_type: Literal["missense", "lof"], 
                         variant_id: int, gene: str, hgvs_p: str,
                         penetrance_mean: float, penetrance_median: float,
                         ci_lower: float, ci_upper: float,
                         model_version: str = None,
                         n_cases: int = None, n_carriers: int = None):
        """Store penetrance estimate."""
        sql = """
        INSERT INTO penetrance_estimates 
            (variant_type, variant_id, gene, hgvs_p, 
             penetrance_mean, penetrance_median, ci_lower, ci_upper,
             model_version, n_cases, n_carriers)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        ON CONFLICT(variant_type, variant_id) DO UPDATE SET
            penetrance_mean = excluded.penetrance_mean,
            penetrance_median = excluded.penetrance_median,
            ci_lower = excluded.ci_lower,
            ci_upper = excluded.ci_upper,
            model_version = excluded.model_version,
            n_cases = excluded.n_cases,
            n_carriers = excluded.n_carriers,
            estimated_at = CURRENT_TIMESTAMP
        """
        self.conn.execute(sql, [
            variant_type, variant_id, gene, hgvs_p,
            penetrance_mean, penetrance_median, ci_lower, ci_upper,
            model_version, n_cases, n_carriers
        ])
        self.conn.commit()
    
    def close(self):
        self.conn.close()
