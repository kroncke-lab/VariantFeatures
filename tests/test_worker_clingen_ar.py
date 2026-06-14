"""Tests for the clingen_ar handler + the generic worker loop."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import pytest

from variantfeatures.database import VariantDB
from variantfeatures.handlers import clingen_ar
from variantfeatures.worker import run_pending


SAMPLE_PAYLOAD = {
    "@id": "https://reg.genome.network/allele/CA999",
    "genomicAlleles": [
        {
            "referenceGenome": "GRCh38",
            "chromosome": "7",
            "hgvs": ["NC_000007.14:g.150951925A>G"],
            "coordinates": [
                {"start": 150951924, "end": 150951925, "referenceAllele": "A", "allele": "G"}
            ],
        }
    ],
    "transcriptAlleles": [
        {"hgvs": ["NM_000238.4:c.1841T>C"], "proteinEffect": {"hgvs": "NP_000229.1:p.Met614Thr"}}
    ],
    "externalRecords": {
        "dbSNP": [{"rs": 199472953}],
        "ClinVarVariations": [{"variationId": 67244}],
        "gnomAD": [{"id": "7-150951925-A-G"}],
    },
}


@pytest.fixture
def db(tmp_path: Path) -> VariantDB:
    return VariantDB(tmp_path / "test.db")


# ---------------------------------------------------------------------------
# hgvs_g building
# ---------------------------------------------------------------------------

def test_hgvs_g_for_snv_chr7():
    assert clingen_ar.hgvs_g_for_snv("7", 150951925, "A", "G") == "NC_000007.14:g.150951925A>G"


def test_hgvs_g_for_snv_strips_chr_prefix():
    assert clingen_ar.hgvs_g_for_snv("chr7", 1, "A", "G") == "NC_000007.14:g.1A>G"


def test_hgvs_g_for_snv_x_chromosome():
    assert clingen_ar.hgvs_g_for_snv("X", 1, "A", "G") == "NC_000023.11:g.1A>G"


def test_hgvs_g_for_snv_unknown_chromosome_raises():
    with pytest.raises(clingen_ar.HandlerError):
        clingen_ar.hgvs_g_for_snv("99", 1, "A", "G")


def test_hgvs_g_for_snv_rejects_indel():
    with pytest.raises(clingen_ar.HandlerError):
        clingen_ar.hgvs_g_for_snv("7", 1, "A", "AT")  # insertion
    with pytest.raises(clingen_ar.HandlerError):
        clingen_ar.hgvs_g_for_snv("7", 1, "AT", "A")  # deletion


# ---------------------------------------------------------------------------
# handle() with mocked HTTP
# ---------------------------------------------------------------------------

class _FakeResp:
    def __init__(self, status_code=200, payload=None):
        self.status_code = status_code
        self._payload = payload

    def json(self):
        return self._payload

    def raise_for_status(self):
        if self.status_code >= 400:
            from requests import HTTPError
            raise HTTPError(f"{self.status_code}")


def test_handle_persists_ca_id_and_aliases(db):
    vid = db.upsert_variant(chromosome="7", position=150951925, ref="A", alt="G", variant_type="SNV")

    with patch("variantfeatures.identity.requests.get", return_value=_FakeResp(200, SAMPLE_PAYLOAD)):
        clingen_ar.handle(db, vid)

    cur = db.conn.execute("SELECT ca_id, hgvs_g FROM variants WHERE id = ?", [vid])
    row = cur.fetchone()
    assert row["ca_id"] == "CA999"
    assert row["hgvs_g"] == "NC_000007.14:g.150951925A>G"

    aliases = {(a["alias_type"], a["alias_value"]) for a in db.get_aliases(vid)}
    assert ("ca_id", "CA999") in aliases
    assert ("rsid", "rs199472953") in aliases
    assert ("clinvar_vcv", "VCV000067244") in aliases
    assert ("gnomad_id", "7-150951925-A-G") in aliases


def test_handle_unknown_variant_raises_handler_error(db):
    with pytest.raises(clingen_ar.HandlerError):
        clingen_ar.handle(db, variant_id=9999)  # not in table


# ---------------------------------------------------------------------------
# run_pending dispatch loop
# ---------------------------------------------------------------------------

def test_run_pending_drains_clingen_ar(db):
    v1 = db.upsert_variant(chromosome="7", position=1, ref="A", alt="G")
    v2 = db.upsert_variant(chromosome="7", position=2, ref="C", alt="T")
    db.enqueue_job(v1, "clingen_ar")
    db.enqueue_job(v2, "clingen_ar")

    with patch("variantfeatures.identity.requests.get", return_value=_FakeResp(200, SAMPLE_PAYLOAD)):
        summary = run_pending(db, source="clingen_ar", rate_limit_sec=0)

    assert summary["claimed"] == 2
    assert summary["done"] == 2
    assert summary["failed"] == 0
    cur = db.conn.execute("SELECT status, COUNT(*) AS n FROM annotation_jobs GROUP BY status")
    assert {r["status"]: r["n"] for r in cur} == {"done": 2}


def test_run_pending_marks_failed_on_handler_error(db):
    """If ClinGen AR returns 400 (e.g. reference mismatch), job is marked failed but worker keeps going."""
    v1 = db.upsert_variant(chromosome="7", position=1, ref="A", alt="G")
    db.enqueue_job(v1, "clingen_ar")

    error_resp = _FakeResp(400, {"errorType": "IncorrectReferenceAllele",
                                  "description": "Reference allele does not match"})
    with patch("variantfeatures.identity.requests.get", return_value=error_resp):
        summary = run_pending(db, source="clingen_ar", rate_limit_sec=0)

    assert summary["done"] == 0
    assert summary["failed"] == 1
    cur = db.conn.execute("SELECT status, error FROM annotation_jobs WHERE variant_id = ?", [v1])
    row = cur.fetchone()
    assert row["status"] == "failed"
    assert "Reference allele" in (row["error"] or "")


def test_run_pending_respects_max_jobs(db):
    v1 = db.upsert_variant(chromosome="7", position=1, ref="A", alt="G")
    v2 = db.upsert_variant(chromosome="7", position=2, ref="C", alt="T")
    v3 = db.upsert_variant(chromosome="7", position=3, ref="G", alt="A")
    for v in (v1, v2, v3):
        db.enqueue_job(v, "clingen_ar")

    with patch("variantfeatures.identity.requests.get", return_value=_FakeResp(200, SAMPLE_PAYLOAD)):
        summary = run_pending(db, source="clingen_ar", max_jobs=2, rate_limit_sec=0)

    assert summary["claimed"] == 2
    cur = db.conn.execute("SELECT status, COUNT(*) AS n FROM annotation_jobs GROUP BY status")
    counts = {r["status"]: r["n"] for r in cur}
    assert counts.get("pending", 0) == 1
    assert counts.get("done", 0) == 2


def test_run_pending_filters_by_gene(db):
    kcnh2_variant = db.upsert_variant(chromosome="7", position=1, ref="A", alt="G")
    brca1_variant = db.upsert_variant(chromosome="17", position=2, ref="C", alt="T")
    db.upsert_consequence(
        kcnh2_variant,
        "NM_KCNH2.1",
        "enumerated",
        gene_symbol="KCNH2",
        consequence="missense_variant",
    )
    db.upsert_consequence(
        brca1_variant,
        "NM_BRCA1.1",
        "enumerated",
        gene_symbol="BRCA1",
        consequence="missense_variant",
    )
    db.enqueue_job(kcnh2_variant, "clingen_ar")
    db.enqueue_job(brca1_variant, "clingen_ar")

    with patch("variantfeatures.identity.requests.get", return_value=_FakeResp(200, SAMPLE_PAYLOAD)):
        summary = run_pending(db, source="clingen_ar", gene_filter="BRCA1", rate_limit_sec=0)

    assert summary["claimed"] == 1
    rows = {
        row["variant_id"]: row["status"]
        for row in db.conn.execute("SELECT variant_id, status FROM annotation_jobs")
    }
    assert rows == {
        kcnh2_variant: "pending",
        brca1_variant: "done",
    }


def test_run_pending_skips_unknown_source(db):
    v1 = db.upsert_variant(chromosome="7", position=1, ref="A", alt="G")
    db.enqueue_job(v1, "imaginary_source")

    summary = run_pending(db, source="imaginary_source", rate_limit_sec=0)
    assert summary["claimed"] == 1
    assert summary["skipped"] == 1
    assert summary["done"] == 0


def test_run_pending_idempotent_when_rerun(db):
    """A done job is not re-claimed on a second run."""
    v1 = db.upsert_variant(chromosome="7", position=1, ref="A", alt="G")
    db.enqueue_job(v1, "clingen_ar")

    with patch("variantfeatures.identity.requests.get", return_value=_FakeResp(200, SAMPLE_PAYLOAD)):
        run_pending(db, source="clingen_ar", rate_limit_sec=0)
        second = run_pending(db, source="clingen_ar", rate_limit_sec=0)

    assert second["claimed"] == 0
    assert second["done"] == 0
