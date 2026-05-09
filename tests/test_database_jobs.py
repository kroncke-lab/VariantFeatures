"""Tests for variant_consequences and annotation_jobs in VariantDB."""

from __future__ import annotations

from pathlib import Path

import pytest

from variantfeatures.database import VariantDB


@pytest.fixture
def db(tmp_path: Path) -> VariantDB:
    return VariantDB(tmp_path / "test.db")


@pytest.fixture
def variant_id(db) -> int:
    return db.upsert_variant(
        chromosome="7", position=1, ref="A", alt="G",
        ca_id="CA1", variant_type="SNV",
    )


# ---------------------------------------------------------------------------
# variant_consequences
# ---------------------------------------------------------------------------

def test_upsert_consequence_insert(db, variant_id):
    db.upsert_consequence(
        variant_id=variant_id,
        transcript_id="NM_000238.4",
        source="enumerated",
        gene_symbol="KCNH2",
        consequence="missense_variant",
        hgvs_c="NM_000238.4:c.1841T>C",
        hgvs_p="NP_000229.1:p.Met614Thr",
        aa_pos=614,
        aa_ref="M",
        aa_alt="T",
        codon_pos=614,
        codon_ref="ATG",
        codon_alt="ACG",
        is_canonical=1,
        is_mane_select=1,
    )
    rows = db.get_consequences(variant_id)
    assert len(rows) == 1
    assert rows[0]["consequence"] == "missense_variant"
    assert rows[0]["aa_pos"] == 614


def test_upsert_consequence_idempotent_on_pk(db, variant_id):
    db.upsert_consequence(
        variant_id=variant_id, transcript_id="NM_x", source="enumerated",
        consequence="missense_variant", aa_pos=1,
    )
    db.upsert_consequence(
        variant_id=variant_id, transcript_id="NM_x", source="enumerated",
        consequence="missense_variant", aa_pos=2,  # value updated
    )
    rows = db.get_consequences(variant_id)
    assert len(rows) == 1
    assert rows[0]["aa_pos"] == 2


def test_upsert_consequence_distinct_sources(db, variant_id):
    """Same (variant, transcript) but different sources stay as separate rows."""
    db.upsert_consequence(
        variant_id=variant_id, transcript_id="NM_x", source="enumerated",
        consequence="missense_variant",
    )
    db.upsert_consequence(
        variant_id=variant_id, transcript_id="NM_x", source="vep",
        consequence="missense_variant",
    )
    rows = db.get_consequences(variant_id)
    assert {r["source"] for r in rows} == {"enumerated", "vep"}


# ---------------------------------------------------------------------------
# annotation_jobs
# ---------------------------------------------------------------------------

def test_enqueue_job_inserts_pending(db, variant_id):
    db.enqueue_job(variant_id, source="alphamissense", priority=10)
    cur = db.conn.execute("SELECT status, priority, attempts FROM annotation_jobs WHERE variant_id = ?", [variant_id])
    row = cur.fetchone()
    assert row["status"] == "pending"
    assert row["priority"] == 10
    assert row["attempts"] == 0


def test_enqueue_job_idempotent(db, variant_id):
    db.enqueue_job(variant_id, source="alphamissense")
    db.enqueue_job(variant_id, source="alphamissense")  # second call no-ops
    cur = db.conn.execute("SELECT COUNT(*) AS n FROM annotation_jobs WHERE variant_id = ?", [variant_id])
    assert cur.fetchone()["n"] == 1


def test_claim_pending_jobs_marks_running(db):
    v1 = db.upsert_variant(chromosome="7", position=1, ref="A", alt="G")
    v2 = db.upsert_variant(chromosome="7", position=2, ref="A", alt="G")
    db.enqueue_job(v1, "alphamissense", priority=10)
    db.enqueue_job(v2, "alphamissense", priority=20)

    claimed = db.claim_pending_jobs(source="alphamissense", limit=10)
    assert len(claimed) == 2
    assert claimed[0]["variant_id"] == v1  # priority 10 first
    assert claimed[0]["attempts"] == 0  # the row returned reflects pre-claim attempts

    # now nothing pending
    again = db.claim_pending_jobs(source="alphamissense", limit=10)
    assert again == []

    # rows in DB should be 'running' with attempts=1
    cur = db.conn.execute("SELECT status, attempts FROM annotation_jobs ORDER BY id")
    statuses = [(r["status"], r["attempts"]) for r in cur.fetchall()]
    assert statuses == [("running", 1), ("running", 1)]


def test_claim_pending_jobs_filters_by_source(db, variant_id):
    db.enqueue_job(variant_id, "alphamissense")
    db.enqueue_job(variant_id, "revel")
    am_jobs = db.claim_pending_jobs(source="alphamissense", limit=10)
    assert len(am_jobs) == 1
    assert am_jobs[0]["source"] == "alphamissense"


def test_mark_job_done(db, variant_id):
    db.enqueue_job(variant_id, "alphamissense")
    [job] = db.claim_pending_jobs(source="alphamissense", limit=10)
    db.mark_job_done(job["id"])
    cur = db.conn.execute("SELECT status, error FROM annotation_jobs WHERE id = ?", [job["id"]])
    row = cur.fetchone()
    assert row["status"] == "done"
    assert row["error"] is None


def test_mark_job_failed_records_error(db, variant_id):
    db.enqueue_job(variant_id, "alphamissense")
    [job] = db.claim_pending_jobs(source="alphamissense", limit=10)
    db.mark_job_failed(job["id"], "boom")
    cur = db.conn.execute("SELECT status, error FROM annotation_jobs WHERE id = ?", [job["id"]])
    row = cur.fetchone()
    assert row["status"] == "failed"
    assert row["error"] == "boom"


def test_job_status_counts(db):
    v1 = db.upsert_variant(chromosome="7", position=1, ref="A", alt="G")
    v2 = db.upsert_variant(chromosome="7", position=2, ref="A", alt="G")
    db.enqueue_job(v1, "alphamissense")
    db.enqueue_job(v2, "alphamissense")
    db.enqueue_job(v1, "revel")

    [job] = db.claim_pending_jobs(source="revel", limit=10)
    db.mark_job_done(job["id"])

    counts = db.job_status_counts()
    by_key = {(c["status"], c["source"]): c["n"] for c in counts}
    assert by_key[("pending", "alphamissense")] == 2
    assert by_key[("done", "revel")] == 1
