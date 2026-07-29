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
    ALLOW_LOCAL_ENV,
    LocalStorageError,
    external_path_state,
    require_external_storage,
)


@pytest.fixture
def fake_repo(monkeypatch, tmp_path):
    """Repoint the guard's repo-root anchor at an empty tmp_path checkout."""
    monkeypatch.setattr(local_storage, "_ANCHORS", (tmp_path,))
    monkeypatch.delenv(ALLOW_LOCAL_ENV, raising=False)
    return tmp_path


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
