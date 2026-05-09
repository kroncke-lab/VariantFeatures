"""Tests for the generic annotation tables (pathogenicity / population / clinical / conservation / splice)."""

from __future__ import annotations

from pathlib import Path

import pytest

from variantfeatures.database import VariantDB


@pytest.fixture
def db(tmp_path: Path) -> VariantDB:
    return VariantDB(tmp_path / "test.db")


@pytest.fixture
def vid(db) -> int:
    return db.upsert_variant(chromosome="7", position=1, ref="A", alt="G")


# ---------------------------------------------------------------------------
# pathogenicity
# ---------------------------------------------------------------------------

def test_upsert_pathogenicity_inserts(db, vid):
    db.upsert_pathogenicity(vid, "alphamissense", score=0.95, category="P", source="myvariant")
    cur = db.conn.execute("SELECT predictor, score, category, source FROM annotations_pathogenicity WHERE variant_id = ?", [vid])
    row = dict(cur.fetchone())
    assert row == {"predictor": "alphamissense", "score": 0.95, "category": "P", "source": "myvariant"}


def test_upsert_pathogenicity_updates_on_rerun(db, vid):
    db.upsert_pathogenicity(vid, "revel", score=0.5, source="myvariant")
    db.upsert_pathogenicity(vid, "revel", score=0.7, rank_score=0.9, source="myvariant")
    cur = db.conn.execute("SELECT score, rank_score FROM annotations_pathogenicity WHERE variant_id = ? AND predictor = 'revel'", [vid])
    row = dict(cur.fetchone())
    assert row["score"] == 0.7
    assert row["rank_score"] == 0.9


def test_upsert_pathogenicity_versions_are_distinct(db, vid):
    db.upsert_pathogenicity(vid, "alphamissense", predictor_version="1.0", score=0.5)
    db.upsert_pathogenicity(vid, "alphamissense", predictor_version="2.0", score=0.7)
    cur = db.conn.execute(
        "SELECT predictor_version, score FROM annotations_pathogenicity WHERE variant_id = ? AND predictor = 'alphamissense' ORDER BY predictor_version",
        [vid],
    )
    rows = [dict(r) for r in cur.fetchall()]
    assert rows == [{"predictor_version": "1.0", "score": 0.5}, {"predictor_version": "2.0", "score": 0.7}]


# ---------------------------------------------------------------------------
# population
# ---------------------------------------------------------------------------

def test_upsert_population_inserts_and_updates(db, vid):
    db.upsert_population(vid, "gnomad_exomes_v4", "all", af=1e-5, ac=1, an=100000, n_homozygotes=0, filter_status="PASS", source="myvariant")
    db.upsert_population(vid, "gnomad_exomes_v4", "afr", af=2e-5, ac=1, an=50000, source="myvariant")
    db.upsert_population(vid, "gnomad_exomes_v4", "all", af=2e-5, ac=2, an=200000)  # update
    cur = db.conn.execute(
        "SELECT pop, af, ac, an FROM annotations_population WHERE variant_id = ? ORDER BY pop",
        [vid],
    )
    rows = [dict(r) for r in cur.fetchall()]
    assert rows == [
        {"pop": "afr", "af": 2e-5, "ac": 1, "an": 50000},
        {"pop": "all", "af": 2e-5, "ac": 2, "an": 200000},
    ]


# ---------------------------------------------------------------------------
# clinical
# ---------------------------------------------------------------------------

def test_upsert_clinical_distinct_records_coexist(db, vid):
    db.upsert_clinical(vid, "clinvar", "VCV1", classification="Pathogenic", stars=2)
    db.upsert_clinical(vid, "clinvar", "VCV2", classification="Likely benign", stars=1)
    cur = db.conn.execute("SELECT record_id, classification FROM annotations_clinical WHERE variant_id = ? ORDER BY record_id", [vid])
    rows = [dict(r) for r in cur.fetchall()]
    assert rows == [
        {"record_id": "VCV1", "classification": "Pathogenic"},
        {"record_id": "VCV2", "classification": "Likely benign"},
    ]


def test_upsert_clinical_updates_existing(db, vid):
    db.upsert_clinical(vid, "clinvar", "VCV1", classification="VUS", stars=1)
    db.upsert_clinical(vid, "clinvar", "VCV1", classification="Likely pathogenic", stars=2)
    cur = db.conn.execute("SELECT classification, stars FROM annotations_clinical WHERE variant_id = ?", [vid])
    row = dict(cur.fetchone())
    assert row == {"classification": "Likely pathogenic", "stars": 2}


# ---------------------------------------------------------------------------
# conservation
# ---------------------------------------------------------------------------

def test_upsert_conservation(db, vid):
    db.upsert_conservation(vid, "phylop100way_vertebrate", score=7.5, rank_score=0.99)
    db.upsert_conservation(vid, "gerp_pp_rs", score=5.2)
    cur = db.conn.execute("SELECT metric, score FROM annotations_conservation WHERE variant_id = ? ORDER BY metric", [vid])
    rows = [dict(r) for r in cur.fetchall()]
    assert rows == [
        {"metric": "gerp_pp_rs", "score": 5.2},
        {"metric": "phylop100way_vertebrate", "score": 7.5},
    ]


# ---------------------------------------------------------------------------
# splice
# ---------------------------------------------------------------------------

def test_upsert_splice_distinct_score_types(db, vid):
    db.upsert_splice(vid, "spliceai", "donor_gain", score=0.01, distance=-15)
    db.upsert_splice(vid, "spliceai", "donor_loss", score=0.02, distance=-15)
    db.upsert_splice(vid, "spliceai", "acceptor_gain", score=0.03)
    db.upsert_splice(vid, "spliceai", "acceptor_loss", score=0.99, distance=-2)
    cur = db.conn.execute(
        "SELECT score_type, score FROM annotations_splice WHERE variant_id = ? ORDER BY score_type",
        [vid],
    )
    rows = [dict(r) for r in cur.fetchall()]
    assert {r["score_type"] for r in rows} == {"donor_gain", "donor_loss", "acceptor_gain", "acceptor_loss"}


def test_annotation_cascade_on_variant_delete(db, vid):
    db.upsert_pathogenicity(vid, "revel", score=0.5)
    db.upsert_population(vid, "gnomad_exomes_v4", "all", af=1e-5)
    db.upsert_clinical(vid, "clinvar", "VCV1", classification="Pathogenic")
    db.upsert_conservation(vid, "phylop100way_vertebrate", score=7.5)
    db.upsert_splice(vid, "spliceai", "overall", score=0.0)

    db.conn.execute("PRAGMA foreign_keys = ON")
    db.conn.execute("DELETE FROM variants WHERE id = ?", [vid])
    db.conn.commit()

    for table in ("annotations_pathogenicity", "annotations_population", "annotations_clinical", "annotations_conservation", "annotations_splice"):
        cur = db.conn.execute(f"SELECT COUNT(*) AS n FROM {table} WHERE variant_id = ?", [vid])
        assert cur.fetchone()["n"] == 0, f"{table} did not cascade delete"
