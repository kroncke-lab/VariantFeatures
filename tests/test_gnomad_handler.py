"""Tests for the gnomAD GraphQL handler."""

from __future__ import annotations

import os
from pathlib import Path
from unittest.mock import patch

import pytest

from variantfeatures.database import VariantDB
from variantfeatures.handlers import gnomad as gn


SAMPLE_VARIANT = {
    "variant_id": "7-140753336-A-T",
    "rsids": ["rs113488022"],
    "exome": {
        "ac": 1,
        "an": 251260,
        "af": 4e-06,
        "homozygote_count": 0,
        "filters": [],
        "populations": [
            {"id": "afr", "ac": 0, "an": 16252, "homozygote_count": 0},
            {"id": "amr", "ac": 0, "an": 34528, "homozygote_count": 0},
            {"id": "sas", "ac": 1, "an": 30000, "homozygote_count": 0},
            # Sex-stratified entries that we should skip:
            {"id": "afr_XX", "ac": 0, "an": 10070, "homozygote_count": 0},
            {"id": "afr_XY", "ac": 0, "an": 6182, "homozygote_count": 0},
            # Bare aggregate sex entries (no pop prefix) — also skip:
            {"id": "XX", "ac": 0, "an": 130000, "homozygote_count": 0},
            {"id": "XY", "ac": 1, "an": 121260, "homozygote_count": 0},
        ],
        "faf95": {"popmax": 1.5e-05, "popmax_population": "sas"},
    },
    "genome": None,
    "joint": {
        "ac": 1,
        "an": 251260,
        "homozygote_count": 0,
        "filters": [],
        "populations": [
            {"id": "afr", "ac": 0, "an": 16252, "homozygote_count": 0},
            {"id": "sas", "ac": 1, "an": 30000, "homozygote_count": 0},
        ],
    },
}


@pytest.fixture
def db(tmp_path: Path) -> VariantDB:
    return VariantDB(tmp_path / "test.db")


@pytest.fixture
def vid(db) -> int:
    return db.upsert_variant(chromosome="7", position=140753336, ref="A", alt="T", variant_type="SNV")


# ---------------------------------------------------------------------------
# parse_gnomad_variant
# ---------------------------------------------------------------------------

def test_parse_gnomad_variant_emits_per_pop_rows():
    rows = list(gn.parse_gnomad_variant(SAMPLE_VARIANT))
    by_key = {(r["dataset"], r["pop"]): r for r in rows}

    assert ("gnomad_exomes_v4", "all") in by_key
    assert by_key[("gnomad_exomes_v4", "all")]["af"] == 4e-06
    assert by_key[("gnomad_exomes_v4", "all")]["ac"] == 1
    assert by_key[("gnomad_exomes_v4", "all")]["an"] == 251260
    assert by_key[("gnomad_exomes_v4", "all")]["filter_status"] == "PASS"

    assert ("gnomad_exomes_v4", "popmax") in by_key
    assert by_key[("gnomad_exomes_v4", "popmax")]["af"] == 1.5e-05
    assert by_key[("gnomad_exomes_v4", "popmax")]["filter_status"] == "sas"

    assert ("gnomad_exomes_v4", "afr") in by_key
    assert by_key[("gnomad_exomes_v4", "afr")]["af"] == 0.0  # 0 / 16252
    assert ("gnomad_exomes_v4", "sas") in by_key
    assert by_key[("gnomad_exomes_v4", "sas")]["ac"] == 1
    assert by_key[("gnomad_exomes_v4", "sas")]["an"] == 30000

    # Sex breakdowns get skipped (both per-pop and aggregate forms)
    assert ("gnomad_exomes_v4", "afr_xx") not in by_key
    assert ("gnomad_exomes_v4", "afr_xy") not in by_key
    assert ("gnomad_exomes_v4", "xx") not in by_key
    assert ("gnomad_exomes_v4", "xy") not in by_key

    # No genome dataset present -> no rows
    assert not any(r["dataset"] == "gnomad_genomes_v4" for r in rows)

    assert ("gnomad_r4_joint", "all") in by_key
    assert by_key[("gnomad_r4_joint", "all")]["ac"] == 1
    assert by_key[("gnomad_r4_joint", "all")]["an"] == 251260
    assert by_key[("gnomad_r4_joint", "all")]["af"] == 1 / 251260
    assert ("gnomad_r4_joint", "sas") in by_key


def test_parse_gnomad_filters_are_serialized():
    payload = dict(SAMPLE_VARIANT)
    payload["exome"] = dict(SAMPLE_VARIANT["exome"], filters=["AC0", "RF"])
    payload["joint"] = dict(SAMPLE_VARIANT["joint"], filters=["AC0", "RF"])
    rows = list(gn.parse_gnomad_variant(payload))
    assert all(r["filter_status"] == "AC0,RF" for r in rows if r["pop"] != "popmax")


def test_parse_aliases_includes_rsid_and_gnomad_id():
    pairs = {(a["alias_type"], a["alias_value"]) for a in gn.parse_aliases(SAMPLE_VARIANT)}
    assert ("gnomad_id", "7-140753336-A-T") in pairs
    assert ("rsid", "rs113488022") in pairs


# ---------------------------------------------------------------------------
# persist + handle
# ---------------------------------------------------------------------------

def test_persist_writes_population_rows(db, vid):
    counts = gn.persist(db, vid, SAMPLE_VARIANT)
    assert counts["population"] >= 4
    assert counts["aliases"] >= 1
    cur = db.conn.execute("SELECT pop, af FROM annotations_population WHERE variant_id = ? AND pop = 'sas'", [vid])
    assert dict(cur.fetchone()) == {"pop": "sas", "af": 1 / 30000}
    cur = db.conn.execute(
        """
        SELECT ac, an, af
        FROM annotations_population
        WHERE variant_id = ? AND dataset = 'gnomad_r4_joint' AND pop = 'all'
        """,
        [vid],
    )
    row = dict(cur.fetchone())
    assert row["ac"] == 1
    assert row["an"] == 251260
    assert row["af"] == 1 / 251260


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


def test_handle_persists_rows(db, vid):
    body = {"data": {"variant": SAMPLE_VARIANT}}
    with patch("variantfeatures.handlers.gnomad.requests.post", return_value=_FakeResp(200, body)):
        gn.handle(db, vid)
    cur = db.conn.execute("SELECT COUNT(*) FROM annotations_population WHERE variant_id = ?", [vid])
    assert cur.fetchone()[0] >= 4


def test_handle_404_in_gnomad_is_silent(db, vid):
    body = {"data": {"variant": None}}
    with patch("variantfeatures.handlers.gnomad.requests.post", return_value=_FakeResp(200, body)):
        gn.handle(db, vid)
    cur = db.conn.execute("SELECT COUNT(*) FROM annotations_population WHERE variant_id = ?", [vid])
    assert cur.fetchone()[0] == 0


def test_handle_unknown_variant_raises(db):
    with pytest.raises(gn.HandlerError):
        gn.handle(db, variant_id=9999)


# ---------------------------------------------------------------------------
# Live integration
# ---------------------------------------------------------------------------
@pytest.mark.skipif(
    os.environ.get("RUN_INTEGRATION") != "1",
    reason="Set RUN_INTEGRATION=1 to hit the live gnomAD API.",
)
def test_handle_braf_v600e_live(tmp_path):
    """Live gnomAD lookup for BRAF V600E."""
    d = VariantDB(tmp_path / "live.db")
    vid = d.upsert_variant(chromosome="7", position=140753336, ref="A", alt="T", variant_type="SNV")
    gn.handle(d, vid, timeout=60)
    cur = d.conn.execute(
        "SELECT pop, af FROM annotations_population WHERE variant_id = ? AND pop = 'all'",
        [vid],
    )
    row = cur.fetchone()
    assert row is not None
    # BRAF V600E is very rare; expect AF less than 1e-3 (or genome-only entry)
    assert row["af"] is not None
