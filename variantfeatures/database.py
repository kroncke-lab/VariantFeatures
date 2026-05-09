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

-- Canonical variant identity (build-agnostic, single row per allele)
-- Coordinates stored as GRCh38, normalized (left-aligned) per VCF spec.
CREATE TABLE IF NOT EXISTS variants (
    id INTEGER PRIMARY KEY,
    ca_id TEXT UNIQUE,                     -- ClinGen Allele Registry ID, e.g., CA12345678
    vrs_id TEXT UNIQUE,                    -- GA4GH VRS computed identifier (optional)

    chromosome TEXT,
    position INTEGER,
    ref TEXT,
    alt TEXT,

    variant_type TEXT,                     -- SNV, MNV, INS, DEL, DELINS, COMPLEX
    hgvs_g TEXT,                           -- canonical g.HGVS on GRCh38

    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,

    UNIQUE(chromosome, position, ref, alt)
);

CREATE INDEX IF NOT EXISTS idx_variants_coords ON variants(chromosome, position);

-- All known names for a variant. One variant -> many aliases.
CREATE TABLE IF NOT EXISTS variant_aliases (
    variant_id INTEGER NOT NULL REFERENCES variants(id) ON DELETE CASCADE,
    alias_type TEXT NOT NULL,              -- rsid, hgvs_g, hgvs_c, hgvs_p, clinvar_vcv, clinvar_rcv, clinvar_allele, gnomad_id, spdi, myvariant_id, ca_id, vrs
    alias_value TEXT NOT NULL,
    source TEXT,                           -- clingen_ar, clinvar, user_input, inferred
    fetched_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    PRIMARY KEY (variant_id, alias_type, alias_value)
);

CREATE INDEX IF NOT EXISTS idx_aliases_lookup ON variant_aliases(alias_type, alias_value);

-- Per-transcript consequence (one variant -> many transcripts)
-- Source 'enumerated' = predicted by saturation mutagenesis; 'vep'/'annovar'/'clingen_ar' = external annotator.
CREATE TABLE IF NOT EXISTS variant_consequences (
    variant_id INTEGER NOT NULL REFERENCES variants(id) ON DELETE CASCADE,
    transcript_id TEXT NOT NULL,           -- e.g. NM_000238.4 or ENST00000262186.10
    gene_symbol TEXT,
    gene_ensembl TEXT,
    consequence TEXT,                       -- SO term: missense_variant, stop_gained, synonymous_variant, intron_variant, ...
    hgvs_c TEXT,
    hgvs_p TEXT,
    aa_pos INTEGER,
    aa_ref TEXT,                            -- 1-letter
    aa_alt TEXT,                            -- 1-letter (or '*' for stop)
    codon_pos INTEGER,                      -- codon number (1-based)
    codon_ref TEXT,
    codon_alt TEXT,
    is_canonical INTEGER DEFAULT 0,
    is_mane_select INTEGER DEFAULT 0,
    is_mane_plus_clinical INTEGER DEFAULT 0,
    source TEXT,                            -- enumerated, vep, annovar, clingen_ar
    fetched_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    PRIMARY KEY (variant_id, transcript_id, source)
);

CREATE INDEX IF NOT EXISTS idx_consequences_gene ON variant_consequences(gene_symbol, consequence);
CREATE INDEX IF NOT EXISTS idx_consequences_transcript ON variant_consequences(transcript_id);

-- Async annotation queue: one row per (variant, source) pair.
CREATE TABLE IF NOT EXISTS annotation_jobs (
    id INTEGER PRIMARY KEY,
    variant_id INTEGER NOT NULL REFERENCES variants(id) ON DELETE CASCADE,
    source TEXT NOT NULL,                   -- alphamissense, revel, cadd, clinvar, gnomad, clingen_ar, myvariant, vep, annovar, ...
    status TEXT NOT NULL DEFAULT 'pending', -- pending, running, done, failed, skipped
    priority INTEGER DEFAULT 100,           -- lower = higher priority
    attempts INTEGER NOT NULL DEFAULT 0,
    last_attempted_at TIMESTAMP,
    error TEXT,
    payload TEXT,                           -- JSON for source-specific parameters
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    UNIQUE(variant_id, source)
);

CREATE INDEX IF NOT EXISTS idx_jobs_status_priority ON annotation_jobs(status, priority, id);
CREATE INDEX IF NOT EXISTS idx_jobs_source_status ON annotation_jobs(source, status);

-- ----------------------------------------------------------------------------
-- Generic per-predictor annotation tables.
-- Adding a new predictor is a data load, not a schema migration.
-- predictor_version is part of the PK with DEFAULT '' so multiple versions can
-- coexist without NULL-comparison surprises in ON CONFLICT.
-- ----------------------------------------------------------------------------

-- In silico pathogenicity / functional predictors (CADD, REVEL, AlphaMissense,
-- dbNSFP suite, PrimateAI, etc.)
CREATE TABLE IF NOT EXISTS annotations_pathogenicity (
    variant_id INTEGER NOT NULL REFERENCES variants(id) ON DELETE CASCADE,
    predictor TEXT NOT NULL,                -- alphamissense, revel, cadd_phred, cadd_raw, metasvm, primateai, ...
    predictor_version TEXT NOT NULL DEFAULT '',
    score REAL,
    rank_score REAL,                        -- 0-1 normalized within predictor (when source provides one)
    category TEXT,                          -- predictor-specific class label, e.g. 'P'/'B'/'A', 'D'/'T'
    source TEXT,                            -- where the value came from: myvariant, dbnsfp, alphamissense_file, vep_plugin
    fetched_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    PRIMARY KEY (variant_id, predictor, predictor_version)
);

CREATE INDEX IF NOT EXISTS idx_path_predictor ON annotations_pathogenicity(predictor);

-- Population allele frequencies (gnomAD exomes/genomes, TOPMed, ExAC, etc.)
CREATE TABLE IF NOT EXISTS annotations_population (
    variant_id INTEGER NOT NULL REFERENCES variants(id) ON DELETE CASCADE,
    dataset TEXT NOT NULL,                  -- gnomad_exomes_v4, gnomad_genomes_v4, exac, topmed_bravo, ...
    pop TEXT NOT NULL DEFAULT 'all',        -- 'all', 'afr', 'amr', 'asj', 'eas', 'fin', 'nfe', 'sas', 'oth', 'popmax'
    af REAL,
    ac INTEGER,
    an INTEGER,
    n_homozygotes INTEGER,
    filter_status TEXT,                     -- 'PASS' or filter codes
    source TEXT,
    fetched_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    PRIMARY KEY (variant_id, dataset, pop)
);

CREATE INDEX IF NOT EXISTS idx_pop_dataset ON annotations_population(dataset, pop);

-- Clinical assertions (ClinVar, HGMD, LOVD, etc.)
CREATE TABLE IF NOT EXISTS annotations_clinical (
    variant_id INTEGER NOT NULL REFERENCES variants(id) ON DELETE CASCADE,
    source TEXT NOT NULL,                   -- clinvar, hgmd, lovd
    record_id TEXT NOT NULL DEFAULT '',     -- VCV/RCV/HGMD ID; '' when source has just one record per variant
    classification TEXT,                    -- 'Pathogenic', 'Likely benign', 'VUS', etc.
    review_status TEXT,
    stars INTEGER,
    last_evaluated DATE,
    conditions TEXT,                        -- JSON or comma-separated
    fetched_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    PRIMARY KEY (variant_id, source, record_id)
);

CREATE INDEX IF NOT EXISTS idx_clin_source ON annotations_clinical(source, classification);

-- Cross-species conservation (phyloP, phastCons, GERP++, SiPhy)
CREATE TABLE IF NOT EXISTS annotations_conservation (
    variant_id INTEGER NOT NULL REFERENCES variants(id) ON DELETE CASCADE,
    metric TEXT NOT NULL,                   -- phylop100way_vertebrate, phastcons100way_vertebrate, gerp_pp_rs, siphy_29way_logodds
    score REAL,
    rank_score REAL,
    source TEXT,
    fetched_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    PRIMARY KEY (variant_id, metric)
);

-- Splice-effect predictors (SpliceAI, dbscSNV, MaxEntScan)
-- score_type is part of the PK because SpliceAI emits 4 directional scores per variant.
CREATE TABLE IF NOT EXISTS annotations_splice (
    variant_id INTEGER NOT NULL REFERENCES variants(id) ON DELETE CASCADE,
    predictor TEXT NOT NULL,                -- spliceai, dbscsnv_ada, dbscsnv_rf, maxentscan
    predictor_version TEXT NOT NULL DEFAULT '',
    score_type TEXT NOT NULL DEFAULT 'overall', -- 'donor_gain', 'donor_loss', 'acceptor_gain', 'acceptor_loss', 'overall'
    score REAL,
    distance INTEGER,                       -- e.g. SpliceAI distance from splice site
    source TEXT,
    fetched_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    PRIMARY KEY (variant_id, predictor, predictor_version, score_type)
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
    
    # ------------------------------------------------------------------
    # Canonical variant identity
    # ------------------------------------------------------------------

    def upsert_variant(
        self,
        chromosome: str,
        position: int,
        ref: str,
        alt: str,
        ca_id: Optional[str] = None,
        vrs_id: Optional[str] = None,
        variant_type: Optional[str] = None,
        hgvs_g: Optional[str] = None,
    ) -> int:
        """Insert or update a canonical variant row. Returns variant_id.

        Coordinates must already be GRCh38, left-aligned, VCF-normalized.
        """
        sql = """
        INSERT INTO variants (chromosome, position, ref, alt, ca_id, vrs_id, variant_type, hgvs_g)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?)
        ON CONFLICT(chromosome, position, ref, alt) DO UPDATE SET
            ca_id = COALESCE(excluded.ca_id, variants.ca_id),
            vrs_id = COALESCE(excluded.vrs_id, variants.vrs_id),
            variant_type = COALESCE(excluded.variant_type, variants.variant_type),
            hgvs_g = COALESCE(excluded.hgvs_g, variants.hgvs_g),
            updated_at = CURRENT_TIMESTAMP
        """
        self.conn.execute(sql, [chromosome, position, ref, alt, ca_id, vrs_id, variant_type, hgvs_g])
        self.conn.commit()
        cur = self.conn.execute(
            "SELECT id FROM variants WHERE chromosome = ? AND position = ? AND ref = ? AND alt = ?",
            [chromosome, position, ref, alt],
        )
        return cur.fetchone()["id"]

    def add_aliases(self, variant_id: int, aliases: list[dict], source: Optional[str] = None) -> int:
        """Insert aliases for a variant. Each alias dict has alias_type and alias_value.

        Returns count of new rows inserted (existing duplicates are ignored).
        """
        if not aliases:
            return 0
        rows = [
            (variant_id, a["alias_type"], a["alias_value"], a.get("source") or source)
            for a in aliases
        ]
        cur = self.conn.executemany(
            "INSERT OR IGNORE INTO variant_aliases (variant_id, alias_type, alias_value, source) VALUES (?, ?, ?, ?)",
            rows,
        )
        self.conn.commit()
        return cur.rowcount

    def get_variant_by_alias(self, alias_type: str, alias_value: str) -> Optional[dict]:
        """Find a canonical variant by any alias. Returns variant row + aliases list, or None."""
        cur = self.conn.execute(
            """
            SELECT v.* FROM variants v
            JOIN variant_aliases a ON a.variant_id = v.id
            WHERE a.alias_type = ? AND a.alias_value = ?
            LIMIT 1
            """,
            [alias_type, alias_value],
        )
        row = cur.fetchone()
        if not row:
            return None
        result = dict(row)
        result["aliases"] = self.get_aliases(result["id"])
        return result

    def get_aliases(self, variant_id: int) -> list[dict]:
        """Return all aliases for a variant_id."""
        cur = self.conn.execute(
            "SELECT alias_type, alias_value, source, fetched_at FROM variant_aliases WHERE variant_id = ? ORDER BY alias_type, alias_value",
            [variant_id],
        )
        return [dict(row) for row in cur.fetchall()]

    # ------------------------------------------------------------------
    # Per-transcript consequences
    # ------------------------------------------------------------------

    def upsert_consequence(self, variant_id: int, transcript_id: str, source: str, **fields) -> None:
        """Insert or update a consequence row for (variant, transcript, source)."""
        all_fields = {"variant_id": variant_id, "transcript_id": transcript_id, "source": source, **fields}
        columns = list(all_fields.keys())
        placeholders = ", ".join(["?"] * len(columns))
        update_cols = [c for c in columns if c not in ("variant_id", "transcript_id", "source")]
        updates = ", ".join(f"{c} = excluded.{c}" for c in update_cols)
        if updates:
            updates += ", fetched_at = CURRENT_TIMESTAMP"
        else:
            updates = "fetched_at = CURRENT_TIMESTAMP"
        sql = f"""
        INSERT INTO variant_consequences ({", ".join(columns)})
        VALUES ({placeholders})
        ON CONFLICT(variant_id, transcript_id, source) DO UPDATE SET {updates}
        """
        self.conn.execute(sql, list(all_fields.values()))
        self.conn.commit()

    def get_consequences(self, variant_id: int) -> list[dict]:
        cur = self.conn.execute(
            "SELECT * FROM variant_consequences WHERE variant_id = ? ORDER BY transcript_id, source",
            [variant_id],
        )
        return [dict(row) for row in cur.fetchall()]

    # ------------------------------------------------------------------
    # Annotation job queue
    # ------------------------------------------------------------------

    def enqueue_job(self, variant_id: int, source: str, *, priority: int = 100, payload: Optional[str] = None) -> None:
        """Add a pending annotation job. Idempotent: if (variant, source) already exists, leave its current status alone."""
        sql = """
        INSERT INTO annotation_jobs (variant_id, source, priority, payload)
        VALUES (?, ?, ?, ?)
        ON CONFLICT(variant_id, source) DO NOTHING
        """
        self.conn.execute(sql, [variant_id, source, priority, payload])
        self.conn.commit()

    def enqueue_source_for_all_variants(self, source: str, *, priority: int = 100, gene_filter: Optional[str] = None) -> int:
        """Enqueue a job for every variant in the DB (optionally filtered to one gene).

        Useful when you add a new annotator after the canonical variants have already
        been populated. Returns the number of new pending jobs created.
        """
        if gene_filter:
            sql = """
            INSERT OR IGNORE INTO annotation_jobs (variant_id, source, priority)
            SELECT DISTINCT v.id, ?, ?
            FROM variants v JOIN variant_consequences c ON c.variant_id = v.id
            WHERE c.gene_symbol = ?
            """
            cur = self.conn.execute(sql, [source, priority, gene_filter])
        else:
            sql = """
            INSERT OR IGNORE INTO annotation_jobs (variant_id, source, priority)
            SELECT id, ?, ? FROM variants
            """
            cur = self.conn.execute(sql, [source, priority])
        self.conn.commit()
        return cur.rowcount

    def claim_pending_jobs(self, source: Optional[str] = None, limit: int = 100) -> list[dict]:
        """Mark up to `limit` pending jobs as running and return them.

        If `source` is given, only that source's jobs are claimed.
        """
        if source is not None:
            cur = self.conn.execute(
                """
                SELECT id, variant_id, source, priority, attempts, payload
                FROM annotation_jobs
                WHERE status = 'pending' AND source = ?
                ORDER BY priority, id
                LIMIT ?
                """,
                [source, limit],
            )
        else:
            cur = self.conn.execute(
                """
                SELECT id, variant_id, source, priority, attempts, payload
                FROM annotation_jobs
                WHERE status = 'pending'
                ORDER BY priority, id
                LIMIT ?
                """,
                [limit],
            )
        rows = [dict(r) for r in cur.fetchall()]
        if rows:
            ids = [r["id"] for r in rows]
            placeholders = ",".join("?" * len(ids))
            self.conn.execute(
                f"UPDATE annotation_jobs SET status = 'running', attempts = attempts + 1, last_attempted_at = CURRENT_TIMESTAMP, updated_at = CURRENT_TIMESTAMP WHERE id IN ({placeholders})",
                ids,
            )
            self.conn.commit()
        return rows

    def mark_job_done(self, job_id: int) -> None:
        self.conn.execute(
            "UPDATE annotation_jobs SET status = 'done', error = NULL, updated_at = CURRENT_TIMESTAMP WHERE id = ?",
            [job_id],
        )
        self.conn.commit()

    def mark_job_failed(self, job_id: int, error: str) -> None:
        self.conn.execute(
            "UPDATE annotation_jobs SET status = 'failed', error = ?, updated_at = CURRENT_TIMESTAMP WHERE id = ?",
            [error[:1000], job_id],
        )
        self.conn.commit()

    def job_status_counts(self) -> list[dict]:
        """Return rows of (status, source, n) grouped counts for the job queue."""
        cur = self.conn.execute(
            "SELECT status, source, COUNT(*) AS n FROM annotation_jobs GROUP BY status, source ORDER BY status, source"
        )
        return [dict(r) for r in cur.fetchall()]

    # ------------------------------------------------------------------
    # Generic annotation upserts
    # ------------------------------------------------------------------

    def upsert_pathogenicity(
        self,
        variant_id: int,
        predictor: str,
        *,
        score: Optional[float] = None,
        predictor_version: str = "",
        rank_score: Optional[float] = None,
        category: Optional[str] = None,
        source: Optional[str] = None,
    ) -> None:
        sql = """
        INSERT INTO annotations_pathogenicity
            (variant_id, predictor, predictor_version, score, rank_score, category, source)
        VALUES (?, ?, ?, ?, ?, ?, ?)
        ON CONFLICT(variant_id, predictor, predictor_version) DO UPDATE SET
            score = excluded.score,
            rank_score = excluded.rank_score,
            category = excluded.category,
            source = excluded.source,
            fetched_at = CURRENT_TIMESTAMP
        """
        self.conn.execute(sql, [variant_id, predictor, predictor_version, score, rank_score, category, source])
        self.conn.commit()

    def upsert_population(
        self,
        variant_id: int,
        dataset: str,
        pop: str = "all",
        *,
        af: Optional[float] = None,
        ac: Optional[int] = None,
        an: Optional[int] = None,
        n_homozygotes: Optional[int] = None,
        filter_status: Optional[str] = None,
        source: Optional[str] = None,
    ) -> None:
        sql = """
        INSERT INTO annotations_population
            (variant_id, dataset, pop, af, ac, an, n_homozygotes, filter_status, source)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
        ON CONFLICT(variant_id, dataset, pop) DO UPDATE SET
            af = excluded.af,
            ac = excluded.ac,
            an = excluded.an,
            n_homozygotes = excluded.n_homozygotes,
            filter_status = excluded.filter_status,
            source = excluded.source,
            fetched_at = CURRENT_TIMESTAMP
        """
        self.conn.execute(sql, [variant_id, dataset, pop, af, ac, an, n_homozygotes, filter_status, source])
        self.conn.commit()

    def upsert_clinical(
        self,
        variant_id: int,
        source: str,
        record_id: str = "",
        *,
        classification: Optional[str] = None,
        review_status: Optional[str] = None,
        stars: Optional[int] = None,
        last_evaluated: Optional[str] = None,
        conditions: Optional[str] = None,
    ) -> None:
        sql = """
        INSERT INTO annotations_clinical
            (variant_id, source, record_id, classification, review_status, stars, last_evaluated, conditions)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?)
        ON CONFLICT(variant_id, source, record_id) DO UPDATE SET
            classification = excluded.classification,
            review_status = excluded.review_status,
            stars = excluded.stars,
            last_evaluated = excluded.last_evaluated,
            conditions = excluded.conditions,
            fetched_at = CURRENT_TIMESTAMP
        """
        self.conn.execute(sql, [variant_id, source, record_id, classification, review_status, stars, last_evaluated, conditions])
        self.conn.commit()

    def upsert_conservation(
        self,
        variant_id: int,
        metric: str,
        *,
        score: Optional[float] = None,
        rank_score: Optional[float] = None,
        source: Optional[str] = None,
    ) -> None:
        sql = """
        INSERT INTO annotations_conservation
            (variant_id, metric, score, rank_score, source)
        VALUES (?, ?, ?, ?, ?)
        ON CONFLICT(variant_id, metric) DO UPDATE SET
            score = excluded.score,
            rank_score = excluded.rank_score,
            source = excluded.source,
            fetched_at = CURRENT_TIMESTAMP
        """
        self.conn.execute(sql, [variant_id, metric, score, rank_score, source])
        self.conn.commit()

    def upsert_splice(
        self,
        variant_id: int,
        predictor: str,
        score_type: str = "overall",
        *,
        score: Optional[float] = None,
        predictor_version: str = "",
        distance: Optional[int] = None,
        source: Optional[str] = None,
    ) -> None:
        sql = """
        INSERT INTO annotations_splice
            (variant_id, predictor, predictor_version, score_type, score, distance, source)
        VALUES (?, ?, ?, ?, ?, ?, ?)
        ON CONFLICT(variant_id, predictor, predictor_version, score_type) DO UPDATE SET
            score = excluded.score,
            distance = excluded.distance,
            source = excluded.source,
            fetched_at = CURRENT_TIMESTAMP
        """
        self.conn.execute(sql, [variant_id, predictor, predictor_version, score_type, score, distance, source])
        self.conn.commit()

    def close(self):
        self.conn.close()
