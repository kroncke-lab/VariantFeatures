"""`doctor` must tell you the state of things, and tell the right person the right thing.

A collaborator cloning this repo has no external volume and no database. The
failure mode being pinned here is advice that only makes sense on Brett's
machine: telling someone to mount a drive they have never heard of is worse than
saying nothing, because it reads like a broken setup rather than a normal fresh
clone.

So the guidance is derived from the symlink actually present — a dangling link
names its own volume, and a checkout with no link gets the setup path instead.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from click.testing import CliRunner

from variantfeatures import local_storage
from variantfeatures.cli import main
from variantfeatures.database import SCHEMA_VERSION, VariantDB
from variantfeatures.local_storage import (
    ALLOW_DB_CREATE_ENV,
    ALLOW_LOCAL_ENV,
    LocalStorageError,
    external_link_target,
    external_volume_name,
    require_external_storage,
)


@pytest.fixture
def fake_repo(monkeypatch, tmp_path):
    monkeypatch.setattr(local_storage, "_ANCHORS", (tmp_path,))
    monkeypatch.delenv(ALLOW_LOCAL_ENV, raising=False)
    monkeypatch.delenv(ALLOW_DB_CREATE_ENV, raising=False)
    return tmp_path


# ---------------------------------------------------------------------------
# Volume detection — the basis for who gets told what
# ---------------------------------------------------------------------------


def test_volume_name_comes_from_the_link_not_a_hardcoded_string(fake_repo):
    """A drive named something else must be named correctly."""
    (fake_repo / "data").symlink_to("/Volumes/SomeOtherDisk/research/data")

    assert external_volume_name("data") == "SomeOtherDisk"


def test_link_target_is_readable_while_the_volume_is_unmounted(fake_repo):
    """readlink, not resolve: an unmounted target still reports where it points."""
    (fake_repo / "data").symlink_to("/Volumes/NotMounted/research/data")

    assert external_link_target("data") == Path("/Volumes/NotMounted/research/data")


def test_no_link_means_no_volume_to_name(fake_repo):
    assert external_volume_name("data") is None
    assert external_link_target("data") is None


def test_a_non_volumes_link_names_no_volume(fake_repo):
    """A link to plain local storage should not be described as a mounted drive."""
    target = fake_repo / "elsewhere"
    target.mkdir()
    (fake_repo / "data").symlink_to(target)

    assert external_volume_name("data") is None


# ---------------------------------------------------------------------------
# Guidance targeting
# ---------------------------------------------------------------------------


def test_detached_drive_names_the_volume_to_attach(fake_repo):
    (fake_repo / "data").symlink_to("/Volumes/Ezekers2/ResearchData/data")

    with pytest.raises(LocalStorageError) as excinfo:
        require_external_storage(fake_repo / "data")

    message = str(excinfo.value)
    assert "Ezekers2" in message
    assert "/Volumes/Ezekers2" in message


def test_fresh_clone_is_not_told_to_mount_someone_elses_drive(fake_repo):
    """The regression: no link to read means no volume may be named."""
    with pytest.raises(LocalStorageError) as excinfo:
        require_external_storage(fake_repo / "data")

    message = str(excinfo.value)
    assert "Volumes" not in message
    assert "mount" not in message.lower()
    # ...and it must say what to actually do instead.
    assert "fresh clone" in message
    assert ALLOW_LOCAL_ENV in message
    assert "ln -s" in message


def test_every_refusal_points_at_doctor(fake_repo):
    with pytest.raises(LocalStorageError) as excinfo:
        require_external_storage(fake_repo / "data")

    assert "variantfeatures doctor" in str(excinfo.value)


# ---------------------------------------------------------------------------
# The doctor command
# ---------------------------------------------------------------------------


def test_doctor_reports_a_fresh_clone_as_not_ready(fake_repo):
    result = CliRunner().invoke(main, ["doctor", "--db", str(fake_repo / "data" / "variants.db")])

    assert result.exit_code == 1, result.output
    assert "Not ready" in result.output
    assert "absent   data" in result.output
    assert "no database file" in result.output


def test_doctor_reports_a_working_setup_as_ready(fake_repo):
    storage = fake_repo / "external"
    storage.mkdir()
    (fake_repo / "data").symlink_to(storage)
    db_path = fake_repo / "data" / "variants.db"
    seed = fake_repo / "seed.db"
    VariantDB(seed).conn.close()
    seed.rename(db_path)

    result = CliRunner().invoke(main, ["doctor", "--db", str(db_path)])

    assert result.exit_code == 0, result.output
    assert "Ready." in result.output
    assert f"version {SCHEMA_VERSION}" in result.output


def test_doctor_flags_a_detached_volume_by_name(fake_repo):
    (fake_repo / "data").symlink_to("/Volumes/Ezekers2/ResearchData/data")

    result = CliRunner().invoke(main, ["doctor", "--db", str(fake_repo / "data" / "x.db")])

    assert result.exit_code == 1, result.output
    assert "MISSING  data" in result.output
    assert "Ezekers2" in result.output


def test_doctor_reports_active_overrides(fake_repo, monkeypatch):
    """An override silently in effect is exactly what you want surfaced."""
    monkeypatch.setenv(ALLOW_DB_CREATE_ENV, "1")

    result = CliRunner().invoke(main, ["doctor", "--db", str(fake_repo / "data" / "x.db")])

    assert ALLOW_DB_CREATE_ENV in result.output
    assert "from scratch is permitted" in result.output


def test_doctor_never_raises_on_a_zero_byte_database(fake_repo):
    """doctor is the thing you run when broken; it must not itself blow up."""
    storage = fake_repo / "external"
    storage.mkdir()
    (fake_repo / "data").symlink_to(storage)
    (storage / "variants.db").touch()

    result = CliRunner().invoke(main, ["doctor", "--db", str(fake_repo / "data" / "variants.db")])

    assert result.exit_code == 1, result.output
    assert "empty" in result.output
    assert result.exception is None or isinstance(result.exception, SystemExit)
