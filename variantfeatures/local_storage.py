"""Fail-closed guard for the repo paths that are backed by external storage.

On Brett's workstation the repo path ``data/`` is an absolute symlink to an
external APFS volume (see AGENTS.md > "Local Data Storage"); ``variants.db`` and
the multi-GB source caches live there, not on the internal disk. Any write that
would otherwise ``mkdir(parents=True)`` its way into ``data/`` has to answer one
question first: is the link actually present?

A *dangling* link already fails on its own — ``mkdir`` refuses to follow one. A
*missing* link does not: ``mkdir`` cheerfully creates a real local ``data/``, and
the job then builds a second database (or re-downloads a 6GB TSV) onto the
internal disk that no later job reads. That split is silent, slow to notice, and
expensive to undo, so the guard refuses the write instead.
"""

from __future__ import annotations

import os
from pathlib import Path

ALLOW_LOCAL_ENV = "VARIANTFEATURES_ALLOW_LOCAL_DATA"
ALLOW_DB_CREATE_ENV = "VARIANTFEATURES_ALLOW_DB_CREATE"

# Repo-relative paths that belong on external storage and must never be created.
EXTERNAL_PATHS = ("data",)

#: The checkout this code belongs to, for reporting.
REPO_ROOT = Path(__file__).resolve().parent.parent

# Both anchors for the same tree: the unresolved one matches how callers build
# their paths (``DEFAULT_DB`` is anchored on a bare ``__file__``), the resolved
# one matches the same location reached through a symlinked checkout.
_ANCHORS = tuple(
    dict.fromkeys(
        (
            Path(__file__).parent.parent,
            REPO_ROOT,
        )
    )
)

_TRUTHY = frozenset({"1", "true", "yes", "on"})


class LocalStorageError(Exception):
    """Raised when a write would create local storage that belongs on an external volume."""


class MissingDatabaseError(LocalStorageError):
    """Raised when the expected database is absent and creating it would rebuild from scratch."""


def local_storage_allowed() -> bool:
    """True when the operator has opted into plain local storage for ``data/``.

    The escape hatch exists for a machine that legitimately has no external
    volume (CI, a collaborator's laptop, a throwaway checkout).
    """
    return os.environ.get(ALLOW_LOCAL_ENV, "").strip().lower() in _TRUTHY


def external_path_state(name: str) -> str:
    """Classify repo-relative path ``name`` as ``linked``, ``dangling``, or ``absent``.

    ``linked`` covers both a symlink whose target is reachable and a real local
    directory — either way the storage exists and a caller may write into it.
    """
    for anchor in _ANCHORS:
        path = anchor / name
        if path.is_dir():
            return "linked"
        if path.is_symlink():
            return "dangling"
    return "absent"


def external_link_target(name: str) -> Path | None:
    """The symlink target of repo-relative ``name``, whether or not it resolves.

    Read with ``os.readlink`` rather than ``resolve()`` so an unmounted volume
    still reports where it was pointed — that string is what tells an operator
    which drive to attach.
    """
    for anchor in _ANCHORS:
        path = anchor / name
        if path.is_symlink():
            try:
                return Path(os.readlink(path))
            except OSError:
                return None
    return None


def external_volume_name(name: str) -> str | None:
    """The macOS volume a repo-relative link points into, e.g. ``Ezekers``.

    Derived from the link itself rather than hardcoded, so this says the right
    thing on a machine whose drive is named something else — and says nothing at
    all on a collaborator's checkout that has no link to read.
    """
    target = external_link_target(name)
    if target is None:
        return None
    parts = target.parts
    try:
        return parts[parts.index("Volumes") + 1]
    except (ValueError, IndexError):
        return None


def _remedy_lines(name: str) -> list:
    """Guidance for restoring ``name``, matched to what the checkout looks like.

    A dangling link means a known drive is detached — say which one. Nothing to
    read means a fresh clone, where telling someone to mount a volume they have
    never heard of is worse than saying nothing: they need the setup path.
    """
    volume = external_volume_name(name)
    if volume:
        return [
            f"Attach the volume {volume!r} and mount it at /Volumes/{volume}, then retry.",
            f"The {name}/ symlink points at {external_link_target(name)}.",
        ]
    return [
        f"This checkout has no {name}/ at all, which is the normal state of a "
        f"fresh clone.",
        "Set one up as either an absolute symlink to your own storage:",
        f"    ln -s /path/to/your/storage {name}",
        f"or a plain local directory, with {ALLOW_LOCAL_ENV}=1 set:",
        f"    mkdir {name}",
    ]


def _guarded_root(target: Path) -> str | None:
    """Return the guarded repo-relative name containing ``target``, if any.

    ``abspath`` normalizes without resolving symlinks, so a path under the repo's
    ``data/`` link keeps its repo-relative shape and still matches.
    """
    try:
        absolute = Path(os.path.abspath(target))
    except (OSError, ValueError):
        return None
    for name in EXTERNAL_PATHS:
        for anchor in _ANCHORS:
            if absolute.is_relative_to(anchor / name):
                return name
    return None


def require_external_storage(target: Path) -> None:
    """Refuse to create ``target`` when its external-storage link is not in place.

    Paths outside the guarded trees — a pytest ``tmp_path``, an explicit
    ``--db`` somewhere else, an already-resolved path on the volume itself — are
    not this guard's business and pass straight through.
    """
    name = _guarded_root(target)
    if name is None:
        return

    state = external_path_state(name)
    if state == "linked" or local_storage_allowed():
        return

    found = {
        "absent": f"there is no {name}/ in the repo root",
        "dangling": f"{name}/ is a symlink whose target is unreachable",
    }[state]
    lines = [
        f"cannot use {name}/ — {found}.",
        "",
        f"Refusing to create it, because {target} on the internal disk would be a "
        f"second copy that no later job reads.",
        "",
        *_remedy_lines(name),
        "",
        "For the full picture of what is and is not set up:",
        "    python -m variantfeatures doctor",
    ]
    raise LocalStorageError("\n".join(lines))


def db_creation_allowed() -> bool:
    """True when the operator has opted into building a database from scratch."""
    return os.environ.get(ALLOW_DB_CREATE_ENV, "").strip().lower() in _TRUTHY


def require_existing_database(db_path: Path) -> None:
    """Refuse to conjure a new database where an existing one was expected.

    The external-storage guard catches an unmounted volume, but not the case that
    actually costs the most: the volume is mounted, `data/` resolves, and
    `variants.db` simply is not there. `sqlite3.connect` creates the file,
    `_init_schema` writes the schema, and the run starts re-enumerating and
    re-annotating from zero — days of API calls to rebuild what already existed.

    So an absent (or truncated, zero-byte) database inside the guarded tree is a
    hard stop. Rebuilding from scratch stays possible, but has to be asked for.
    Databases outside the guarded tree — a pytest ``tmp_path``, an explicit
    ``--db`` elsewhere — are created as before.
    """
    if _guarded_root(db_path) is None or db_creation_allowed():
        return

    try:
        exists = db_path.is_file()
        size = db_path.stat().st_size if exists else 0
    except OSError:
        exists, size = False, 0

    if exists and size > 0:
        return

    problem = (
        f"{db_path} exists but is empty (0 bytes)"
        if exists
        else f"{db_path} does not exist"
    )
    volume = external_volume_name("data")
    where = (
        f"The storage is reachable, so the database was moved, removed, or never "
        f"copied onto the volume {volume!r}."
        if volume
        else "The storage is reachable — it just does not contain a database yet."
    )
    lines = [
        f"no database to read — {problem}.",
        "",
        f"{where} Refusing to create it here, because that would start a full "
        f"re-enumeration and re-annotation from zero rather than reading the "
        f"database it should have found.",
        "",
        "If you expect an existing database, put it in place first:",
        "    ls -l data/variants.db",
        "",
        "If you are starting fresh and *want* to build one (this takes hours to "
        "days per gene and makes many API calls), opt in explicitly:",
        f"    {ALLOW_DB_CREATE_ENV}=1 python -m variantfeatures build --gene KCNH2",
        "",
        "For the full picture of what is and is not set up:",
        "    python -m variantfeatures doctor",
    ]
    raise MissingDatabaseError("\n".join(lines))
