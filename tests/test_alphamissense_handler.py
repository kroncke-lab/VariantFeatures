"""Tests for the AlphaMissense local-file batch handler."""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from variantfeatures.database import VariantDB
from variantfeatures.handlers import alphamissense as am


# Synthetic AlphaMissense file content. Matches the actual schema.
SAMPLE_TSV = (
    "# Copyright DeepMind\n"
    "#CHROM\tPOS\tREF\tALT\tgenome\tuniprot_id\ttranscript_id\tprotein_variant\tam_pathogenicity\tam_class\n"
    # Two rows that match canonical KCNH2 saturation variants:
    "7\t150977910\tG\tT\thg38\tQ12809\tENST00000262186.10\tKCNH2_P2T\t0.85\tlikely_pathogenic\n"
    "7\t150977899\tC\tA\thg38\tQ12809\tENST00000262186.10\tKCNH2_R5S\t0.20\tlikely_benign\n"
    # A row that does NOT correspond to any pending job (different position):
    "7\t999999\tA\tG\thg38\tQ12809\tENST00000262186.10\tKCNH2_X9Z\t0.50\tambiguous\n"
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


# ---------------------------------------------------------------------------
# resolve_file_path
# ---------------------------------------------------------------------------

def test_resolve_file_path_explicit_wins(tmp_path):
    p = tmp_path / "x.tsv.gz"
    assert am.resolve_file_path(str(p)) == p


def test_resolve_file_path_env(monkeypatch, tmp_path):
    p = tmp_path / "from_env.tsv.gz"
    monkeypatch.setenv("ALPHAMISSENSE_FILE", str(p))
    assert am.resolve_file_path() == p


# ---------------------------------------------------------------------------
# run_batch
# ---------------------------------------------------------------------------

def test_run_batch_matches_pending_jobs_and_marks_done(db, am_file):
    v_match1 = db.upsert_variant(chromosome="7", position=150977910, ref="G", alt="T", variant_type="SNV")
    v_match2 = db.upsert_variant(chromosome="7", position=150977899, ref="C", alt="A", variant_type="SNV")
    v_no_match = db.upsert_variant(chromosome="7", position=12345, ref="A", alt="C", variant_type="SNV")
    for v in (v_match1, v_match2, v_no_match):
        db.enqueue_job(v, "alphamissense")

    result = am.run_batch(db, file_path=str(am_file))
    assert result["claimed"] == 3
    assert result["matched"] == 2
    assert result["failed"] == 1

    # Pathogenicity rows
    cur = db.conn.execute("SELECT score, category FROM annotations_pathogenicity WHERE variant_id = ? AND predictor = 'alphamissense'", [v_match1])
    row = dict(cur.fetchone())
    assert row["score"] == 0.85
    assert row["category"] == "P"

    cur = db.conn.execute("SELECT score, category FROM annotations_pathogenicity WHERE variant_id = ? AND predictor = 'alphamissense'", [v_match2])
    row = dict(cur.fetchone())
    assert row["score"] == 0.20
    assert row["category"] == "B"

    # Job statuses
    cur = db.conn.execute("SELECT variant_id, status, error FROM annotation_jobs ORDER BY variant_id")
    statuses = {r["variant_id"]: (r["status"], r["error"]) for r in cur.fetchall()}
    assert statuses[v_match1][0] == "done"
    assert statuses[v_match2][0] == "done"
    assert statuses[v_no_match][0] == "failed"
    assert "no AlphaMissense entry" in (statuses[v_no_match][1] or "")


def test_run_batch_writes_hgvs_p_alias(db, am_file):
    vid = db.upsert_variant(chromosome="7", position=150977910, ref="G", alt="T", variant_type="SNV")
    db.enqueue_job(vid, "alphamissense")
    am.run_batch(db, file_path=str(am_file))
    aliases = {(a["alias_type"], a["alias_value"]) for a in db.get_aliases(vid)}
    assert ("hgvs_p", "ENST00000262186.10:p.P2T") in aliases


def test_run_batch_no_pending_jobs_is_noop(db, am_file):
    result = am.run_batch(db, file_path=str(am_file))
    assert result == {"claimed": 0, "matched": 0, "failed": 0, "lines_scanned": 0}


def test_run_batch_missing_file_raises(db, tmp_path):
    vid = db.upsert_variant(chromosome="7", position=1, ref="A", alt="G")
    db.enqueue_job(vid, "alphamissense")
    with pytest.raises(am.HandlerError):
        am.run_batch(db, file_path=str(tmp_path / "no_such_file.tsv.gz"))


def test_run_batch_respects_limit(db, am_file):
    v1 = db.upsert_variant(chromosome="7", position=150977910, ref="G", alt="T", variant_type="SNV")
    v2 = db.upsert_variant(chromosome="7", position=150977899, ref="C", alt="A", variant_type="SNV")
    db.enqueue_job(v1, "alphamissense")
    db.enqueue_job(v2, "alphamissense")
    result = am.run_batch(db, file_path=str(am_file), limit=1)
    assert result["claimed"] == 1
