"""Tests for the AlphaMissense protein-keyed batch handler."""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from variantfeatures.database import VariantDB
from variantfeatures.handlers import alphamissense as am


# Synthetic AlphaMissense file content matching the actual published format:
# 4 columns: uniprot_id, protein_variant, am_pathogenicity, am_class.
SAMPLE_TSV = (
    "# Copyright DeepMind\n"
    "uniprot_id\tprotein_variant\tam_pathogenicity\tam_class\n"
    # KCNH2 = Q12809 (resolved dynamically to this accession in tests)
    "Q12809\tA614V\t0.92\tlikely_pathogenic\n"
    "Q12809\tP2T\t0.30\tambiguous\n"
    "Q12809\tR5W\t0.95\tlikely_pathogenic\n"
    # An entry that's not requested by any pending job:
    "Q12809\tX9Y\t0.50\tambiguous\n"
    # LMNA = P02545 (never hardcoded; proves gene-agnostic resolution):
    "P02545\tR100C\t0.88\tlikely_pathogenic\n"
)


@pytest.fixture
def db(tmp_path: Path) -> VariantDB:
    return VariantDB(tmp_path / "test.db")


@pytest.fixture
def am_file(tmp_path: Path) -> Path:
    p = tmp_path / "AlphaMissense_aa_substitutions.tsv.gz"
    with gzip.open(p, "wt") as f:
        f.write(SAMPLE_TSV)
    return p


@pytest.fixture(autouse=True)
def offline_uniprot(monkeypatch):
    """Keep gene->UniProt resolution offline and deterministic.

    KCNH2 resolves to its real accession; any other unlisted gene is treated as
    a UniProt REST miss so the gene-agnostic path is exercised without network.
    """
    known = {"KCNH2": "Q12809"}

    def fake_resolve(gene, *, timeout=30, extra=None):
        g = gene.upper()
        if g in known:
            return known[g]
        raise am.UniProtError(f"No reviewed human UniProt entry found for gene {g!r}")

    monkeypatch.setattr(am, "resolve_uniprot_accession", fake_resolve)


def _add_kcnh2_variant(db, *, position, ref, alt, aa_ref, aa_pos, aa_alt, consequence="missense_variant"):
    """Insert a canonical variant + a KCNH2 enumerated consequence."""
    vid = db.upsert_variant(chromosome="7", position=position, ref=ref, alt=alt, variant_type="SNV")
    db.upsert_consequence(
        variant_id=vid,
        transcript_id="ENST00000262186.10",
        source="enumerated",
        gene_symbol="KCNH2",
        consequence=consequence,
        aa_ref=aa_ref, aa_pos=aa_pos, aa_alt=aa_alt,
    )
    return vid


# ---------------------------------------------------------------------------
# Path resolution
# ---------------------------------------------------------------------------

def test_resolve_file_path_explicit_wins(tmp_path):
    p = tmp_path / "x.tsv.gz"
    assert am.resolve_file_path(str(p)) == p


def test_resolve_file_path_env(monkeypatch, tmp_path):
    p = tmp_path / "from_env.tsv.gz"
    monkeypatch.setenv("ALPHAMISSENSE_FILE", str(p))
    assert am.resolve_file_path() == p


def test_gene_uniprot_map_has_no_hardcoded_genes():
    # The handler is gene-agnostic: the override map is empty unless configured.
    assert am.gene_uniprot_map() == {}


def test_gene_uniprot_map_env_override(monkeypatch):
    monkeypatch.setenv("ALPHAMISSENSE_UNIPROT", "FOO:Q11111,BAR:Q22222")
    m = am.gene_uniprot_map()
    assert m["FOO"] == "Q11111"
    assert m["BAR"] == "Q22222"


# ---------------------------------------------------------------------------
# Batch runner
# ---------------------------------------------------------------------------

def test_run_batch_matches_pending_missense_jobs(db, am_file):
    v_match1 = _add_kcnh2_variant(db, position=1, ref="C", alt="T", aa_ref="A", aa_pos=614, aa_alt="V")
    v_match2 = _add_kcnh2_variant(db, position=2, ref="C", alt="A", aa_ref="P", aa_pos=2, aa_alt="T")
    # In file but not enqueued -> ignored
    # Enqueued but not in file -> failed
    v_no_match = _add_kcnh2_variant(db, position=3, ref="A", alt="C", aa_ref="L", aa_pos=999, aa_alt="P")
    for v in (v_match1, v_match2, v_no_match):
        db.enqueue_job(v, "alphamissense")

    result = am.run_batch(db, file_path=str(am_file))
    assert result["claimed"] == 3
    assert result["matched"] == 2
    assert result["failed"] == 1

    cur = db.conn.execute("SELECT score, category FROM annotations_pathogenicity WHERE variant_id = ?", [v_match1])
    row = dict(cur.fetchone())
    assert row["score"] == 0.92
    assert row["category"] == "P"

    cur = db.conn.execute("SELECT score, category FROM annotations_pathogenicity WHERE variant_id = ?", [v_match2])
    row = dict(cur.fetchone())
    assert row["score"] == 0.30
    assert row["category"] == "A"

    cur = db.conn.execute("SELECT variant_id, status, error FROM annotation_jobs ORDER BY variant_id")
    statuses = {r["variant_id"]: (r["status"], r["error"]) for r in cur.fetchall()}
    assert statuses[v_match1][0] == "done"
    assert statuses[v_match2][0] == "done"
    assert statuses[v_no_match][0] == "failed"


def test_run_batch_resolves_arbitrary_gene_dynamically(db, am_file):
    """Any gene works without being hardcoded -- here LMNA via an explicit override."""
    vid = db.upsert_variant(chromosome="1", position=500, ref="C", alt="T", variant_type="SNV")
    db.upsert_consequence(
        variant_id=vid, transcript_id="NM_LMNA", source="enumerated",
        gene_symbol="LMNA", consequence="missense_variant",
        aa_ref="R", aa_pos=100, aa_alt="C",
    )
    db.enqueue_job(vid, "alphamissense")

    result = am.run_batch(db, file_path=str(am_file), gene_uniprot={"LMNA": "P02545"})
    assert result["matched"] == 1
    row = dict(db.conn.execute(
        "SELECT score, category FROM annotations_pathogenicity WHERE variant_id = ?", [vid]
    ).fetchone())
    assert row["score"] == 0.88
    assert row["category"] == "P"


def test_run_batch_marks_non_missense_jobs_failed(db, am_file):
    """AM only covers missense; stop_gained jobs should be failed without scanning the file."""
    vid = _add_kcnh2_variant(db, position=10, ref="A", alt="T",
                              aa_ref="Q", aa_pos=11, aa_alt="*",
                              consequence="stop_gained")
    db.enqueue_job(vid, "alphamissense")

    result = am.run_batch(db, file_path=str(am_file))
    assert result["matched"] == 0
    assert result["failed"] == 1
    cur = db.conn.execute("SELECT status, error FROM annotation_jobs WHERE variant_id = ?", [vid])
    row = dict(cur.fetchone())
    assert row["status"] == "failed"
    assert "missense" in (row["error"] or "")


def test_run_batch_unknown_gene_fails_cleanly(db, am_file):
    """A gene without a UniProt mapping fails without scanning the file."""
    vid = db.upsert_variant(chromosome="X", position=1, ref="A", alt="G", variant_type="SNV")
    db.upsert_consequence(
        variant_id=vid, transcript_id="NM_X", source="enumerated",
        gene_symbol="MADEUPGENE", consequence="missense_variant",
        aa_ref="A", aa_pos=1, aa_alt="V",
    )
    db.enqueue_job(vid, "alphamissense")

    result = am.run_batch(db, file_path=str(am_file))
    assert result["matched"] == 0
    assert result["failed"] == 1
    cur = db.conn.execute("SELECT error FROM annotation_jobs WHERE variant_id = ?", [vid])
    err = cur.fetchone()["error"]
    assert "MADEUPGENE" in err and "UniProt" in err


def test_run_batch_blank_gene_symbol_fails_without_resolution(db, am_file):
    """A blank/whitespace gene_symbol fails cleanly, without attempting a UniProt lookup."""
    vid = db.upsert_variant(chromosome="2", position=7, ref="A", alt="G", variant_type="SNV")
    db.upsert_consequence(
        variant_id=vid, transcript_id="NM_BLANK", source="enumerated",
        gene_symbol="   ", consequence="missense_variant",
        aa_ref="A", aa_pos=1, aa_alt="V",
    )
    db.enqueue_job(vid, "alphamissense")

    result = am.run_batch(db, file_path=str(am_file))
    assert result["matched"] == 0
    assert result["failed"] == 1
    err = db.conn.execute(
        "SELECT error FROM annotation_jobs WHERE variant_id = ?", [vid]
    ).fetchone()["error"]
    assert "gene_symbol" in err


def test_run_batch_no_pending_jobs_is_noop(db, am_file):
    result = am.run_batch(db, file_path=str(am_file))
    assert result == {"claimed": 0, "matched": 0, "failed": 0, "lines_scanned": 0}


def test_run_batch_missing_file_raises(db, tmp_path):
    vid = _add_kcnh2_variant(db, position=1, ref="A", alt="G", aa_ref="A", aa_pos=1, aa_alt="V")
    db.enqueue_job(vid, "alphamissense")
    with pytest.raises(am.HandlerError):
        am.run_batch(db, file_path=str(tmp_path / "no_such_file.tsv.gz"))


def test_run_batch_respects_limit(db, am_file):
    v1 = _add_kcnh2_variant(db, position=1, ref="C", alt="T", aa_ref="A", aa_pos=614, aa_alt="V")
    v2 = _add_kcnh2_variant(db, position=2, ref="C", alt="A", aa_ref="P", aa_pos=2, aa_alt="T")
    db.enqueue_job(v1, "alphamissense")
    db.enqueue_job(v2, "alphamissense")
    result = am.run_batch(db, file_path=str(am_file), limit=1)
    assert result["claimed"] == 1
