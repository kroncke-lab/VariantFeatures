"""SQLite database operations for variant features."""

import sqlite3
from pathlib import Path
from typing import Optional

DEFAULT_DB = Path(__file__).parent.parent / "data" / "variants.db"

SCHEMA = """
CREATE TABLE IF NOT EXISTS variants (
    id INTEGER PRIMARY KEY,
    gene TEXT NOT NULL,
    hgvs_p TEXT,
    hgvs_c TEXT,
    chromosome TEXT,
    position INTEGER,
    ref TEXT,
    alt TEXT,
    
    -- Scores
    alphamissense_score REAL,
    alphamissense_class TEXT,
    revel_score REAL,
    cadd_phred REAL,
    cadd_raw REAL,
    
    -- ClinVar
    clinvar_id INTEGER,
    clinvar_significance TEXT,
    clinvar_review_status TEXT,
    clinvar_last_evaluated DATE,
    
    -- gnomAD
    gnomad_af REAL,
    gnomad_af_popmax REAL,
    gnomad_homozygotes INTEGER,
    
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    
    UNIQUE(gene, hgvs_p)
);

CREATE INDEX IF NOT EXISTS idx_gene ON variants(gene);
CREATE INDEX IF NOT EXISTS idx_clinvar ON variants(clinvar_significance);
CREATE INDEX IF NOT EXISTS idx_gnomad_af ON variants(gnomad_af);
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
    
    def upsert_variant(self, gene: str, hgvs_p: str, **features):
        """Insert or update a variant with features."""
        columns = ["gene", "hgvs_p"] + list(features.keys())
        placeholders = ", ".join(["?"] * len(columns))
        updates = ", ".join(f"{k} = excluded.{k}" for k in features.keys())
        updates += ", updated_at = CURRENT_TIMESTAMP"
        
        sql = f"""
        INSERT INTO variants ({", ".join(columns)})
        VALUES ({placeholders})
        ON CONFLICT(gene, hgvs_p) DO UPDATE SET {updates}
        """
        self.conn.execute(sql, [gene, hgvs_p] + list(features.values()))
        self.conn.commit()
    
    def get_variant(self, gene: str, hgvs_p: str) -> Optional[dict]:
        """Get a variant by gene and HGVS."""
        cur = self.conn.execute(
            "SELECT * FROM variants WHERE gene = ? AND hgvs_p = ?",
            [gene, hgvs_p]
        )
        row = cur.fetchone()
        return dict(row) if row else None
    
    def get_gene_variants(self, gene: str) -> list[dict]:
        """Get all variants for a gene."""
        cur = self.conn.execute(
            "SELECT * FROM variants WHERE gene = ? ORDER BY hgvs_p",
            [gene]
        )
        return [dict(row) for row in cur.fetchall()]
    
    def close(self):
        self.conn.close()
