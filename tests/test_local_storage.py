"""The `data/` external-storage guard must fail closed, not create local storage.

`data/` is an absolute symlink to an external volume on Brett's workstation. The
failure this file pins is the *missing* link (a fresh checkout, or a rename): the
old code called `mkdir(parents=True, exist_ok=True)`, which would silently build
a second `variants.db` on the internal disk that no later job reads.

The guard is anchored on the real repo root, so these tests relocate that anchor
to `tmp_path` rather than touching the developer's actual `data/` link.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from variantfeatures import local_storage
from variantfeatures.database import VariantDB
from variantfeatures.local_storage import (
    ALLOW_DB_CREATE_ENV,
    ALLOW_LOCAL_ENV,
    LocalStorageError,
    MissingDatabaseError,
    external_path_state,
    require_existing_database,
    require_external_storage,
)


@pytest.fixture
def fake_repo(monkeypatch, tmp_path):
    """Repoint the guard's repo-root anchor at an empty tmp_path checkout."""
    monkeypatch.setattr(local_storage, "_ANCHORS", (tmp_path,))
    monkeypatch.delenv(ALLOW_LOCAL_ENV, raising=False)
    monkeypatch.delenv(ALLOW_DB_CREATE_ENV, raising=False)
    return tmp_path


@pytest.fixture
def mounted_repo(fake_repo):
    """A fake checkout whose `data` link resolves — the volume is mounted."""
    target = fake_repo / "external"
    target.mkdir()
    (fake_repo / "data").symlink_to(target)
    return fake_repo


# ---------------------------------------------------------------------------
# external_path_state
# ---------------------------------------------------------------------------


def test_state_is_linked_for_a_mounted_symlink(fake_repo):
    target = fake_repo / "external"
    target.mkdir()
    (fake_repo / "data").symlink_to(target)

    assert external_path_state("data") == "linked"


def test_state_is_linked_for_a_real_local_directory(fake_repo):
    (fake_repo / "data").mkdir()

    assert external_path_state("data") == "linked"


def test_state_is_dangling_when_the_volume_is_unmounted(fake_repo):
    (fake_repo / "data").symlink_to(fake_repo / "never_mounted")

    assert external_path_state("data") == "dangling"


def test_state_is_absent_on_a_fresh_checkout(fake_repo):
    assert external_path_state("data") == "absent"


# ---------------------------------------------------------------------------
# require_external_storage
# ---------------------------------------------------------------------------


def test_absent_link_is_refused(fake_repo):
    """The regression: a fresh checkout must not grow a local data/ tree."""
    with pytest.raises(LocalStorageError, match="no 'data' entry"):
        require_external_storage(fake_repo / "data")

    assert not (fake_repo / "data").exists(), "the guard must not create anything"


def test_dangling_link_is_refused_with_the_mount_reason(fake_repo):
    (fake_repo / "data").symlink_to(fake_repo / "never_mounted")

    with pytest.raises(LocalStorageError, match="target is unreachable"):
        require_external_storage(fake_repo / "data" / "alphamissense")


def test_error_names_the_volume_and_the_escape_hatch(fake_repo):
    with pytest.raises(LocalStorageError) as excinfo:
        require_external_storage(fake_repo / "data")

    message = str(excinfo.value)
    assert "Ezekers" in message
    assert ALLOW_LOCAL_ENV in message


def test_nested_paths_under_a_mounted_link_pass(fake_repo):
    target = fake_repo / "external"
    target.mkdir()
    (fake_repo / "data").symlink_to(target)

    require_external_storage(fake_repo / "data" / "revel")


def test_paths_outside_the_guarded_tree_pass(fake_repo, tmp_path):
    """A tmp_path database or an explicit --db elsewhere is not the guard's business."""
    require_external_storage(tmp_path / "somewhere_else")
    require_external_storage(fake_repo / "dist" / "publish")
    require_external_storage(Path("/Volumes/Ezekers/ResearchData/variantFeatures/data"))


def test_data_prefixed_sibling_is_not_treated_as_data(fake_repo):
    """`is_relative_to` must not match a sibling that merely starts with 'data'."""
    require_external_storage(fake_repo / "database_backups")


def test_escape_hatch_allows_local_storage(fake_repo, monkeypatch):
    monkeypatch.setenv(ALLOW_LOCAL_ENV, "1")

    require_external_storage(fake_repo / "data")


@pytest.mark.parametrize("raw", ["0", "", "no", "maybe"])
def test_escape_hatch_stays_closed_for_non_truthy_values(fake_repo, monkeypatch, raw):
    monkeypatch.setenv(ALLOW_LOCAL_ENV, raw)

    with pytest.raises(LocalStorageError):
        require_external_storage(fake_repo / "data")


# ---------------------------------------------------------------------------
# VariantDB wiring
# ---------------------------------------------------------------------------


def test_variantdb_refuses_to_create_a_second_internal_database(fake_repo):
    with pytest.raises(LocalStorageError):
        VariantDB(fake_repo / "data" / "variants.db")

    assert not (fake_repo / "data").exists()


def test_variantdb_still_builds_a_database_outside_the_guarded_tree(tmp_path):
    """The offline suite builds throwaway databases; the guard must stay out of it."""
    db = VariantDB(tmp_path / "scratch" / "variants.db")

    assert db.db_path.exists()
    db.conn.close()


# ---------------------------------------------------------------------------
# require_existing_database — the mounted-but-empty case
# ---------------------------------------------------------------------------


def test_missing_database_on_a_mounted_volume_is_refused(mounted_repo):
    """The expensive case: `data/` resolves fine, `variants.db` simply is not there.

    Nothing above catches this — the link is healthy — so without this guard
    sqlite3 creates the file and the run re-annotates from zero.
    """
    db_path = mounted_repo / "data" / "variants.db"

    with pytest.raises(MissingDatabaseError, match="does not exist"):
        require_existing_database(db_path)

    assert not db_path.exists(), "the guard must not create the database"


def test_zero_byte_database_is_refused(mounted_repo):
    """A truncated remnant would otherwise get a fresh schema written into it."""
    db_path = mounted_repo / "data" / "variants.db"
    db_path.touch()

    with pytest.raises(MissingDatabaseError, match="empty"):
        require_existing_database(db_path)


def test_existing_database_passes(mounted_repo):
    db_path = mounted_repo / "data" / "variants.db"
    db_path.write_bytes(b"SQLite format 3\x00")

    require_existing_database(db_path)


def test_missing_database_error_explains_the_rebuild_and_the_opt_in(mounted_repo):
    with pytest.raises(MissingDatabaseError) as excinfo:
        require_existing_database(mounted_repo / "data" / "variants.db")

    message = str(excinfo.value)
    assert "from scratch" in message
    assert ALLOW_DB_CREATE_ENV in message
    assert "Ezekers" in message


def test_databases_outside_the_guarded_tree_are_still_created(tmp_path):
    require_existing_database(tmp_path / "scratch.db")


def test_opt_in_allows_a_genuine_first_build(mounted_repo, monkeypatch):
    """Brett or a collaborator rebuilding from scratch on purpose must still work."""
    monkeypatch.setenv(ALLOW_DB_CREATE_ENV, "1")

    require_existing_database(mounted_repo / "data" / "variants.db")


@pytest.mark.parametrize("raw", ["0", "", "no", "maybe"])
def test_opt_in_stays_closed_for_non_truthy_values(mounted_repo, monkeypatch, raw):
    monkeypatch.setenv(ALLOW_DB_CREATE_ENV, raw)

    with pytest.raises(MissingDatabaseError):
        require_existing_database(mounted_repo / "data" / "variants.db")


def test_variantdb_refuses_to_rebuild_a_missing_database_from_scratch(mounted_repo):
    """The end-to-end shape of an accidental run: mounted volume, absent database."""
    db_path = mounted_repo / "data" / "variants.db"

    with pytest.raises(MissingDatabaseError):
        VariantDB(db_path)

    assert not db_path.exists()


def test_variantdb_read_only_also_refuses_with_the_clear_error(mounted_repo):
    """A read of an absent database should name the cause, not raise OperationalError."""
    with pytest.raises(MissingDatabaseError):
        VariantDB(mounted_repo / "data" / "variants.db", initialize=False, read_only=True)


def test_variantdb_opens_an_existing_database_in_the_guarded_tree(mounted_repo):
    """The normal mounted case must stay completely unaffected."""
    # Build a real database outside the guarded tree, then move it into place, so
    # the guarded path holds a genuine SQLite file rather than a stub.
    seed_path = mounted_repo / "seed.db"
    seed = VariantDB(seed_path)
    seed.conn.close()
    db_path = mounted_repo / "data" / "variants.db"
    seed_path.rename(db_path)

    db = VariantDB(db_path, initialize=False, read_only=True)

    assert db.db_path == db_path
    assert db.conn.execute("SELECT COUNT(*) FROM genes").fetchone()[0] == 0
    db.conn.close()
