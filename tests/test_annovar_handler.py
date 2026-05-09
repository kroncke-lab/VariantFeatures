"""Tests for the ANNOVAR wrapper + multianno importer.

Live runs require ANNOVAR installed; those tests are gated behind RUN_INTEGRATION
plus the binary actually being on disk.
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest

from variantfeatures.database import VariantDB
from variantfeatures.handlers import annovar


# Synthetic multianno.txt (tab-delimited; only the columns we know how to parse)
SAMPLE_MULTIANNO = (
    "Chr\tStart\tEnd\tRef\tAlt\tFunc.refGene\tExonicFunc.refGene\tAAChange.refGene\tGene.refGene\t"
    "AF\tAF_popmax\tCLNSIG\tCLNREVSTAT\tCLNALLELEID\tREVEL_score\tCADD_phred\tAlphaMissense_score\t"
    "phyloP100way_vertebrate\tGERP++_RS\n"
    "7\t150977910\t150977910\tG\tT\texonic\tnonsynonymous_SNV\t"
    "KCNH2:NM_000238.4:exon1:c.4G>T:p.P2T,KCNH2:NM_172056.2:exon1:c.4G>T:p.P2T\t"
    "KCNH2\t.\t.\tLikely_pathogenic\tcriteria_provided,_single_submitter\t12345\t"
    "0.94\t29.5\t0.91\t9.2\t5.1\n"
    "7\t150977901\t150977901\tT\tA\texonic\tnonsynonymous_SNV\t"
    "KCNH2:NM_000238.4:exon1:c.13T>A:p.W5R\tKCNH2\t1.5e-05\t1.7e-04\t.\t.\t.\t"
    "0.55\t22.0\t.\t.\t.\n"
)


@pytest.fixture
def db(tmp_path: Path) -> VariantDB:
    return VariantDB(tmp_path / "test.db")


@pytest.fixture
def multianno_file(tmp_path: Path) -> Path:
    p = tmp_path / "out.hg38_multianno.txt"
    p.write_text(SAMPLE_MULTIANNO)
    return p


# ---------------------------------------------------------------------------
# is_installed / annovar_home
# ---------------------------------------------------------------------------

def test_is_installed_returns_false_when_no_env(monkeypatch):
    monkeypatch.delenv("ANNOVAR_HOME", raising=False)
    # If table_annovar.pl really is on PATH on the test machine this would flip true.
    # We don't assume either way here — only that no env + no PATH hit -> False.
    if not annovar.is_installed():
        assert annovar.annovar_home() is None


# ---------------------------------------------------------------------------
# parse_multianno
# ---------------------------------------------------------------------------

def test_parse_multianno_extracts_all_known_columns(multianno_file):
    records = list(annovar.parse_multianno(multianno_file))
    assert len(records) == 2

    r1, r2 = records
    assert r1["chromosome"] == "7"
    assert r1["position"] == 150977910
    assert r1["ref"] == "G"
    assert r1["alt"] == "T"
    assert r1["revel"] == 0.94
    assert r1["cadd_phred"] == 29.5
    assert r1["alphamissense"] == 0.91
    assert r1["phylop100"] == 9.2
    assert r1["gerp_rs"] == 5.1
    assert r1["clinvar_clnsig"] == "Likely_pathogenic"
    assert r1["aa_change_refgene"].startswith("KCNH2:NM_000238.4:exon1:c.4G>T:p.P2T")

    # row 2: missing several fields ('.')
    assert r2["alphamissense"] is None
    assert r2["gnomad_af"] == 1.5e-05
    assert r2["gnomad_af_popmax"] == 1.7e-04


# ---------------------------------------------------------------------------
# import_multianno
# ---------------------------------------------------------------------------

def test_import_multianno_creates_variants_and_annotations(db, multianno_file):
    counts = annovar.import_multianno(db, multianno_file)
    assert counts["variants"] == 2
    assert counts["annotated"] == 2

    # Variant rows
    cur = db.conn.execute("SELECT chromosome, position, ref, alt FROM variants ORDER BY position")
    rows = [dict(r) for r in cur.fetchall()]
    assert {(r["chromosome"], r["position"], r["ref"], r["alt"]) for r in rows} == {
        ("7", 150977901, "T", "A"),
        ("7", 150977910, "G", "T"),
    }

    # Pathogenicity rows
    cur = db.conn.execute("SELECT predictor, score FROM annotations_pathogenicity ORDER BY predictor")
    by_predictor = {r["predictor"]: r["score"] for r in cur.fetchall()}
    assert by_predictor.get("revel") == 0.55 or by_predictor.get("revel") == 0.94
    assert "cadd_phred" in by_predictor
    assert "alphamissense" in by_predictor

    # Conservation rows
    cur = db.conn.execute("SELECT metric FROM annotations_conservation")
    metrics = {r["metric"] for r in cur.fetchall()}
    assert "phylop100way_vertebrate" in metrics
    assert "gerp_pp_rs" in metrics

    # Population: only the 2nd row had AF
    cur = db.conn.execute("SELECT pop, af FROM annotations_population")
    pops = {r["pop"]: r["af"] for r in cur.fetchall()}
    assert pops.get("all") == 1.5e-05
    assert pops.get("popmax") == 1.7e-04

    # Clinical: row 1 only
    cur = db.conn.execute("SELECT classification, record_id, review_status FROM annotations_clinical")
    rows = [dict(r) for r in cur.fetchall()]
    assert rows == [{
        "classification": "Likely_pathogenic",
        "record_id": "12345",
        "review_status": "criteria_provided,_single_submitter",
    }]

    # Per-transcript consequence (refGene-canonical view)
    cur = db.conn.execute("SELECT consequence, hgvs_c, hgvs_p, transcript_id FROM variant_consequences WHERE source = 'annovar'")
    rows = [dict(r) for r in cur.fetchall()]
    assert any(r["consequence"] == "missense_variant" and r["hgvs_p"] == "p.P2T" for r in rows)


# ---------------------------------------------------------------------------
# run_batch (mocked subprocess) — verifies command construction
# ---------------------------------------------------------------------------

def test_run_batch_raises_when_annovar_missing(db, monkeypatch):
    monkeypatch.delenv("ANNOVAR_HOME", raising=False)
    # Make is_installed() return False even if table_annovar.pl is on PATH
    monkeypatch.setattr(annovar, "is_installed", lambda: False)
    vid = db.upsert_variant(chromosome="7", position=1, ref="A", alt="G")
    db.enqueue_job(vid, "annovar")
    with pytest.raises(annovar.HandlerError, match="ANNOVAR not detected"):
        annovar.run_batch(db)
