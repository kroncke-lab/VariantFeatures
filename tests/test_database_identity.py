"""Tests for the canonical-variant tables in VariantDB."""

from __future__ import annotations

from pathlib import Path

import pytest

from variantfeatures.database import VariantDB


@pytest.fixture
def db(tmp_path: Path) -> VariantDB:
    return VariantDB(tmp_path / "test.db")


def test_upsert_variant_inserts_and_returns_id(db):
    vid = db.upsert_variant(
        chromosome="7", position=151251182, ref="A", alt="G",
        ca_id="CA12345678", variant_type="SNV",
        hgvs_g="NC_000007.14:g.151251182A>G",
    )
    assert vid > 0


def test_upsert_variant_is_idempotent_on_coords(db):
    vid1 = db.upsert_variant(chromosome="7", position=1, ref="A", alt="G", ca_id="CA1")
    vid2 = db.upsert_variant(chromosome="7", position=1, ref="A", alt="G", ca_id="CA1")
    assert vid1 == vid2


def test_upsert_variant_does_not_overwrite_with_null(db):
    """COALESCE preserves prior values when newer call passes None."""
    vid = db.upsert_variant(
        chromosome="7", position=1, ref="A", alt="G",
        ca_id="CA1", variant_type="SNV", hgvs_g="NC_000007.14:g.1A>G",
    )
    db.upsert_variant(chromosome="7", position=1, ref="A", alt="G")  # no metadata

    cur = db.conn.execute("SELECT ca_id, variant_type, hgvs_g FROM variants WHERE id = ?", [vid])
    row = cur.fetchone()
    assert row["ca_id"] == "CA1"
    assert row["variant_type"] == "SNV"
    assert row["hgvs_g"] == "NC_000007.14:g.1A>G"


def test_add_aliases_inserts(db):
    vid = db.upsert_variant(chromosome="7", position=1, ref="A", alt="G")
    inserted = db.add_aliases(vid, [
        {"alias_type": "rsid", "alias_value": "rs1"},
        {"alias_type": "hgvs_c", "alias_value": "NM_x:c.1A>G"},
    ], source="clingen_ar")
    assert inserted == 2

    aliases = db.get_aliases(vid)
    assert {(a["alias_type"], a["alias_value"]) for a in aliases} == {
        ("rsid", "rs1"),
        ("hgvs_c", "NM_x:c.1A>G"),
    }


def test_add_aliases_dedups(db):
    vid = db.upsert_variant(chromosome="7", position=1, ref="A", alt="G")
    db.add_aliases(vid, [{"alias_type": "rsid", "alias_value": "rs1"}])
    again = db.add_aliases(vid, [{"alias_type": "rsid", "alias_value": "rs1"}])
    assert again == 0
    assert len(db.get_aliases(vid)) == 1


def test_get_variant_by_alias_returns_canonical_row_and_alias_list(db):
    vid = db.upsert_variant(
        chromosome="7", position=1, ref="A", alt="G",
        ca_id="CA1", variant_type="SNV",
    )
    db.add_aliases(vid, [
        {"alias_type": "rsid", "alias_value": "rs1"},
        {"alias_type": "hgvs_c", "alias_value": "NM_x:c.1A>G"},
    ])

    found = db.get_variant_by_alias("rsid", "rs1")
    assert found is not None
    assert found["id"] == vid
    assert found["ca_id"] == "CA1"
    assert {(a["alias_type"], a["alias_value"]) for a in found["aliases"]} == {
        ("rsid", "rs1"),
        ("hgvs_c", "NM_x:c.1A>G"),
    }


def test_get_variant_by_alias_misses(db):
    assert db.get_variant_by_alias("rsid", "rs_nonexistent") is None


def test_aliases_cascade_on_variant_delete(db):
    vid = db.upsert_variant(chromosome="7", position=1, ref="A", alt="G")
    db.add_aliases(vid, [{"alias_type": "rsid", "alias_value": "rs1"}])
    db.conn.execute("PRAGMA foreign_keys = ON")
    db.conn.execute("DELETE FROM variants WHERE id = ?", [vid])
    db.conn.commit()
    assert db.get_aliases(vid) == []
