"""Tests for the VEP wrapper.

Live runs require VEP installed and a cache; these are gated behind RUN_INTEGRATION
and are intentionally NOT included here. Once VEP is on the system, add a live test
following the pattern in test_gnomad_handler.py / test_myvariant_handler.py.
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import pytest

from variantfeatures.database import VariantDB
from variantfeatures.handlers import vep


SAMPLE_VEP_RECORD = {
    "input": "7 140753336 . A T . . .",
    "seq_region_name": "7",
    "start": 140753336,
    "end": 140753336,
    "allele_string": "A/T",
    "transcript_consequences": [
        {
            "transcript_id": "ENST00000288602",
            "gene_id": "ENSG00000157764",
            "gene_symbol": "BRAF",
            "consequence_terms": ["missense_variant"],
            "hgvsc": "ENST00000288602.11:c.1799T>A",
            "hgvsp": "ENSP00000288602.7:p.Val600Glu",
            "protein_start": 600,
            "canonical": 1,
            "mane_select": "NM_004333.6",
            # Plugin output:
            "am_pathogenicity": 0.998,
            "am_class": "likely_pathogenic",
            "revel_score": "0.927",
            "cadd_phred": 32.0,
            "cadd_raw": 5.5,
            "primateai_score": 0.85,
            "spliceai_pred_ds_ag": "0.01",
            "spliceai_pred_ds_dl": 0.99,
            "spliceai_pred_dp_dl": -2,
            "ada_score": "0.91",
            "rf_score": 0.82,
            "MaxEntScan_ref": 8.2,
            "MaxEntScan_alt": 2.2,
            "MaxEntScan_diff": 6.0,
            "NMD": "NMD_escaping_variant",
        },
        {
            "transcript_id": "ENST00000644969",
            "gene_id": "ENSG00000157764",
            "gene_symbol": "BRAF",
            "consequence_terms": ["missense_variant"],
            "hgvsc": "ENST00000644969.2:c.1919T>A",
            "hgvsp": "ENSP00000493543.1:p.Val640Glu",
            "protein_start": 640,
            "am_pathogenicity": 0.994,
            "revel_score": 0.92,
        },
    ],
    "colocated_variants": [
        {"id": "rs113488022"},
    ],
}


@pytest.fixture
def db(tmp_path: Path) -> VariantDB:
    return VariantDB(tmp_path / "test.db")


@pytest.fixture
def vid(db) -> int:
    return db.upsert_variant(chromosome="7", position=140753336, ref="A", alt="T", variant_type="SNV")


# ---------------------------------------------------------------------------
# Config helpers
# ---------------------------------------------------------------------------

def test_is_installed_when_no_binary(monkeypatch, tmp_path):
    monkeypatch.delenv("VEP_BIN", raising=False)
    monkeypatch.setattr(vep.Path, "home", lambda: tmp_path)  # no ~/tools/ensembl-vep/vep
    with patch("variantfeatures.handlers.vep.shutil.which", return_value=None):
        assert vep.is_installed() is False


def test_vep_cache_dir_default(monkeypatch):
    monkeypatch.delenv("VEP_CACHE_DIR", raising=False)
    assert vep.vep_cache_dir() == Path.home() / ".vep"


def test_configured_plugins(monkeypatch):
    # Plugin args within one plugin are comma-separated, so we split on `;` between plugins.
    monkeypatch.setenv("VEP_PLUGINS", "AlphaMissense,file=/x.tsv.gz; REVEL,file=/y.tsv.gz")
    out = vep.configured_plugins()
    assert len(out) == 2
    assert out[0] == "AlphaMissense,file=/x.tsv.gz"
    assert out[1] == "REVEL,file=/y.tsv.gz"


# ---------------------------------------------------------------------------
# _persist_record (the parser)
# ---------------------------------------------------------------------------

def test_persist_record_writes_consequence_pathogenicity_splice_alias(db, vid):
    vep._persist_record(db, vid, SAMPLE_VEP_RECORD)

    # Per-transcript consequence (MANE Select picked)
    cur = db.conn.execute("SELECT transcript_id, consequence, hgvs_p, is_mane_select FROM variant_consequences WHERE source = 'vep' AND variant_id = ?", [vid])
    row = dict(cur.fetchone())
    assert row["transcript_id"] == "ENST00000288602"
    assert row["consequence"] == "missense_variant"
    assert row["hgvs_p"] == "ENSP00000288602.7:p.Val600Glu"
    assert row["is_mane_select"] == 1

    # Pathogenicity scores: AlphaMissense gets the max across transcripts
    cur = db.conn.execute("SELECT score, category FROM annotations_pathogenicity WHERE variant_id = ? AND predictor = 'alphamissense'", [vid])
    row = dict(cur.fetchone())
    assert row["score"] == 0.998
    assert row["category"] == "P"

    cur = db.conn.execute("SELECT score FROM annotations_pathogenicity WHERE variant_id = ? AND predictor = 'cadd_phred'", [vid])
    assert dict(cur.fetchone())["score"] == 32.0

    cur = db.conn.execute("SELECT score FROM annotations_pathogenicity WHERE variant_id = ? AND predictor = 'revel'", [vid])
    assert dict(cur.fetchone())["score"] == 0.927  # string form gets coerced

    # SpliceAI
    cur = db.conn.execute("SELECT score_type, score, distance FROM annotations_splice WHERE variant_id = ? AND predictor = 'spliceai'", [vid])
    rows = [dict(r) for r in cur.fetchall()]
    by_type = {r["score_type"]: r for r in rows}
    assert by_type["acceptor_gain"]["score"] == 0.01
    assert by_type["donor_loss"]["score"] == 0.99
    assert by_type["donor_loss"]["distance"] == -2
    assert by_type["overall"]["score"] == 0.99

    cur = db.conn.execute(
        "SELECT predictor, score_type, score FROM annotations_splice WHERE variant_id = ? AND predictor IN ('dbscsnv_ada', 'dbscsnv_rf', 'maxentscan')",
        [vid],
    )
    splice_rows = {(r["predictor"], r["score_type"]): r["score"] for r in cur.fetchall()}
    assert splice_rows[("dbscsnv_ada", "overall")] == 0.91
    assert splice_rows[("dbscsnv_rf", "overall")] == 0.82
    assert splice_rows[("maxentscan", "diff")] == 6.0

    cur = db.conn.execute("SELECT score, category FROM annotations_pathogenicity WHERE variant_id = ? AND predictor = 'nmd_vep'", [vid])
    row = dict(cur.fetchone())
    assert row == {"score": 0.0, "category": "NMD_escaping_variant"}

    # rsID alias from colocated_variants
    aliases = {(a["alias_type"], a["alias_value"]) for a in db.get_aliases(vid)}
    assert ("rsid", "rs113488022") in aliases


def test_pick_primary_csq_prefers_mane_select():
    csqs = [
        {"transcript_id": "T1", "canonical": 1},
        {"transcript_id": "T2", "mane_select": "NM_xxx"},
    ]
    assert vep._pick_primary_csq(csqs)["transcript_id"] == "T2"


def test_pick_primary_csq_falls_back_to_canonical():
    csqs = [
        {"transcript_id": "T1"},
        {"transcript_id": "T2", "canonical": 1},
    ]
    assert vep._pick_primary_csq(csqs)["transcript_id"] == "T2"


# ---------------------------------------------------------------------------
# run_batch refuses to run when VEP isn't installed
# ---------------------------------------------------------------------------

def test_run_batch_no_binary_raises(db, monkeypatch, tmp_path):
    monkeypatch.delenv("VEP_BIN", raising=False)
    monkeypatch.setattr(vep.Path, "home", lambda: tmp_path)
    with patch("variantfeatures.handlers.vep.shutil.which", return_value=None):
        with pytest.raises(vep.HandlerError, match="VEP binary not found"):
            vep.run_batch(db)


def test_run_batch_no_cache_raises(db, monkeypatch, tmp_path):
    monkeypatch.setenv("VEP_BIN", "/usr/bin/true")  # pretend a binary exists
    monkeypatch.setenv("VEP_CACHE_DIR", str(tmp_path / "no_such_cache"))
    with pytest.raises(vep.HandlerError, match="VEP cache directory"):
        vep.run_batch(db)
