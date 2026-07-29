"""A database must carry the schema shape it was written with.

Before stamping existed there was no way to tell a migrated database from a stale
one. `variants.db` is ~35 GiB and the old flat `variants_missense`/`variants_lof`
tables were already removed once, so reading a pre-migration file as if it were
current would return confidently wrong answers rather than failing.

These use tmp_path databases, which sit outside the external-storage guard, so
the two concerns stay independently testable.
"""

from __future__ import annotations

import sqlite3

import pytest

from variantfeatures import database
from variantfeatures.database import (
    SCHEMA_VERSION,
    SchemaVersionError,
    VariantDB,
)


def _stamp(path, version: int) -> None:
    conn = sqlite3.connect(path)
    conn.execute(f"PRAGMA user_version = {version}")
    conn.commit()
    conn.close()


# ---------------------------------------------------------------------------
# Stamping
# ---------------------------------------------------------------------------


def test_a_new_database_is_stamped_with_the_current_version(tmp_path):
    db = VariantDB(tmp_path / "new.db")
    try:
        assert db.schema_version == SCHEMA_VERSION
    finally:
        db.conn.close()


def test_reopening_a_stamped_database_is_accepted(tmp_path):
    path = tmp_path / "new.db"
    VariantDB(path).conn.close()

    db = VariantDB(path, initialize=False, read_only=True)
    try:
        assert db.schema_version == SCHEMA_VERSION
    finally:
        db.conn.close()


# ---------------------------------------------------------------------------
# Legacy adoption — the real variants.db is unstamped (version 0)
# ---------------------------------------------------------------------------


def test_unstamped_but_current_database_is_adopted(tmp_path):
    """Brett's existing database has user_version 0; it must keep opening."""
    path = tmp_path / "legacy.db"
    VariantDB(path).conn.close()
    _stamp(path, 0)

    db = VariantDB(path, initialize=False, read_only=True)
    try:
        assert db.schema_version == 0, "a read-only open must not rewrite the file"
    finally:
        db.conn.close()


def test_adopted_database_gets_stamped_on_the_next_write_open(tmp_path):
    path = tmp_path / "legacy.db"
    VariantDB(path).conn.close()
    _stamp(path, 0)

    db = VariantDB(path)  # write path -> _init_schema stamps it
    try:
        assert db.schema_version == SCHEMA_VERSION
    finally:
        db.conn.close()


def test_pre_normalized_database_is_refused(tmp_path):
    """An unstamped database still carrying the removed flat tables is stale."""
    path = tmp_path / "old.db"
    VariantDB(path).conn.close()
    conn = sqlite3.connect(path)
    conn.execute("CREATE TABLE variants_missense (id INTEGER PRIMARY KEY)")
    conn.execute("PRAGMA user_version = 0")
    conn.commit()
    conn.close()

    with pytest.raises(SchemaVersionError, match="predates the normalized schema"):
        VariantDB(path, initialize=False, read_only=True)


def test_pre_normalized_database_is_refused_before_a_build_can_write_to_it(tmp_path):
    """The dangerous direction: refusing must happen before _init_schema runs."""
    path = tmp_path / "old.db"
    VariantDB(path).conn.close()
    conn = sqlite3.connect(path)
    conn.execute("CREATE TABLE variants_lof (id INTEGER PRIMARY KEY)")
    conn.execute("PRAGMA user_version = 0")
    conn.commit()
    conn.close()

    with pytest.raises(SchemaVersionError):
        VariantDB(path)

    # Still unstamped: the refusal must not have relabelled a stale file.
    conn = sqlite3.connect(path)
    assert conn.execute("PRAGMA user_version").fetchone()[0] == 0
    conn.close()


# ---------------------------------------------------------------------------
# Version mismatches
# ---------------------------------------------------------------------------


def test_newer_database_is_refused(tmp_path):
    """Older code must not read a file written by a newer schema."""
    path = tmp_path / "future.db"
    VariantDB(path).conn.close()
    _stamp(path, SCHEMA_VERSION + 5)

    with pytest.raises(SchemaVersionError, match="newer VariantFeatures"):
        VariantDB(path, initialize=False, read_only=True)


def test_older_stamped_database_is_refused_pending_migration(tmp_path, monkeypatch):
    """A stamped-but-older database needs a migration, not a silent write.

    At SCHEMA_VERSION 1 nothing can be both stamped and older, so raise the code's
    version to reach the branch — this is exactly the state the first real bump
    will produce, and it must refuse rather than mix two schemas in one file.
    """
    path = tmp_path / "old.db"
    VariantDB(path).conn.close()  # stamped at the current version
    monkeypatch.setattr(database, "SCHEMA_VERSION", SCHEMA_VERSION + 1)

    with pytest.raises(SchemaVersionError, match="migration is required"):
        VariantDB(path, initialize=False, read_only=True)

    with pytest.raises(SchemaVersionError, match="migration is required"):
        VariantDB(path)  # the write path must refuse too, before _init_schema

    conn = sqlite3.connect(path)
    assert conn.execute("PRAGMA user_version").fetchone()[0] == SCHEMA_VERSION, (
        "a refused open must not restamp the file"
    )
    conn.close()


def test_version_check_ignores_an_empty_file(tmp_path):
    """A brand-new file has no tables yet and must not be judged as stale."""
    db = VariantDB(tmp_path / "fresh.db")
    try:
        assert db.schema_version == SCHEMA_VERSION
    finally:
        db.conn.close()
