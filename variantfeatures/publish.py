"""Publish gene-scoped SQLite artifacts to Azure Blob Storage."""

from __future__ import annotations

import hashlib
import json
import os
import sqlite3
import subprocess
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping

from .database import SCHEMA


class PublishError(Exception):
    """Raised when a publish artifact cannot be built or uploaded."""


@dataclass(frozen=True)
class PublishArtifact:
    """Local artifact to include in a publish release."""

    gene: str | None
    path: Path
    kind: str = "gene_slice"


BlobServiceClient = None
DefaultAzureCredential = None
ContentSettings = None
MatchConditions = None
ResourceExistsError = None
ResourceModifiedError = None
ResourceNotFoundError = None


DATA_TABLES = (
    "genes",
    "gene_constraint",
    "transcripts",
    "variants",
    "variant_consequences",
    "variant_aliases",
    "penetrance_estimates",
    "frameshift_nonsense_mappings",
)

GNOMAD_R4_EXOME_DATASET = "gnomad_r4_exome_gene"
GNOMAD_R4_GENOME_DATASET = "gnomad_r4_genome_gene"
GNOMAD_R4_JOINT_DATASET = "gnomad_r4_joint"
GNOMAD_R4_JOINT_SOURCE = "gnomad_gene_api_joint_v4_1"
GNOMAD_API = "https://gnomad.broadinstitute.org/api"


def export_gene_slice(full_db: Path, gene: str, out_path: Path) -> dict[str, Any]:
    """Write a standalone SQLite database containing only rows for ``gene``.

    Variant membership is resolved from ``variant_consequences``. Variant-level
    annotations and aliases are copied for those variant IDs; gene-level rows are
    copied by symbol. Queue state is intentionally not copied.
    """

    full_db = Path(full_db)
    out_path = Path(out_path)
    gene = gene.upper()

    if not full_db.exists():
        raise PublishError(f"Source database does not exist: {full_db}")
    if full_db.resolve() == out_path.resolve():
        raise PublishError("Refusing to overwrite the source database with a slice")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    if out_path.exists():
        out_path.unlink()

    conn = sqlite3.connect(out_path)
    conn.row_factory = sqlite3.Row
    try:
        conn.execute("PRAGMA foreign_keys = OFF")
        conn.executescript(SCHEMA)
        conn.execute("ATTACH DATABASE ? AS src", [str(full_db)])
        conn.execute("CREATE TEMP TABLE gene_variant_ids (id INTEGER PRIMARY KEY)")
        conn.execute(
            """
            INSERT OR IGNORE INTO gene_variant_ids (id)
            SELECT DISTINCT variant_id
            FROM src.variant_consequences
            WHERE UPPER(gene_symbol) = ?
            """,
            [gene],
        )

        annotation_tables = _annotation_tables(conn)
        copied_tables = [
            "genes",
            "gene_constraint",
            "transcripts",
            "variants",
            "variant_consequences",
            "variant_aliases",
            *annotation_tables,
            "penetrance_estimates",
            "frameshift_nonsense_mappings",
        ]

        for table in copied_tables:
            if not _table_exists(conn, "src", table) or not _table_exists(conn, "main", table):
                continue
            _copy_gene_table(conn, table, gene)

        counts = _assert_slice_counts(conn, gene, copied_tables)
        added_joint_rows = _add_gnomad_r4_joint_population_rows(conn, gene)
        if "annotations_population" in counts:
            counts["annotations_population"] = _count(
                conn, "annotations_population", "1 = 1", []
            )
        conn.commit()
        conn.execute("DETACH DATABASE src")
        conn.commit()
        conn.execute("VACUUM")
        summary = {"gene": gene, "path": out_path, "row_counts": counts}
        if added_joint_rows:
            summary["added_rows"] = {
                f"annotations_population.{GNOMAD_R4_JOINT_DATASET}": added_joint_rows
            }
        return summary
    finally:
        conn.close()


def build_manifest(
    db: Path | sqlite3.Connection | Any,
    genes: list[str] | tuple[str, ...],
    built_at: str,
    slices: Mapping[str, Path] | list[PublishArtifact] | tuple[PublishArtifact, ...],
) -> dict[str, Any]:
    """Build a deterministic provenance manifest for a publish release."""

    conn, should_close = _open_read_conn(db)
    try:
        normalized_genes = [gene.upper() for gene in genes]
        artifacts = _normalize_artifacts(slices)
        manifest = {
            "producer": "variantfeatures",
            "git_sha": get_git_sha(),
            "built_at": built_at,
            "schema_version": _read_schema_version(conn),
            "source_versions": _source_versions(conn),
            "genes": {
                gene: {
                    "row_counts": _gene_row_counts(conn, gene),
                    "slice": f"genes/{gene}.db",
                }
                for gene in normalized_genes
            },
            "artifacts": {},
        }

        for artifact in artifacts:
            key = _artifact_relative_path(artifact)
            manifest["artifacts"][key] = {
                "path": key,
                "local_path": str(artifact.path),
                "sha256": sha256_file(artifact.path),
                "bytes": artifact.path.stat().st_size,
                "kind": artifact.kind,
                "gene": artifact.gene,
            }

        return manifest
    finally:
        if should_close:
            conn.close()


def upload(
    account_url: str | None,
    container: str | None,
    prefix: str,
    built_at: str,
    git_sha: str | None,
    slices: Mapping[str, Path] | list[PublishArtifact] | tuple[PublishArtifact, ...],
    manifest: dict[str, Any],
    *,
    dry_run: bool,
) -> dict[str, Any]:
    """Upload artifacts and atomically merge per-gene entries into latest.json."""

    artifacts = _normalize_artifacts(slices)
    release = release_name(built_at, git_sha or manifest.get("git_sha"))
    prefix = _clean_prefix(prefix)
    planned = [_planned_artifact(prefix, release, artifact) for artifact in artifacts]
    manifest_blob_path = f"{prefix}/{release}/manifest.json"
    latest_blob_path = f"{prefix}/latest.json"
    manifest_bytes = _manifest_bytes(manifest)
    latest_entries = _latest_entries(prefix, release, built_at, manifest, artifacts)
    empty_latest: dict[str, Any] = {}
    latest_after = {**empty_latest, **latest_entries}
    latest_diff = _latest_diff(empty_latest, latest_after, latest_entries)

    if dry_run:
        return {
            "dry_run": True,
            "release": release,
            "artifacts": planned,
            "manifest": {
                "blob_path": manifest_blob_path,
                "sha256": sha256_bytes(manifest_bytes),
                "bytes": len(manifest_bytes),
            },
            "latest_blob_path": latest_blob_path,
            "latest_before": empty_latest,
            "latest_after": latest_after,
            "latest_diff": latest_diff,
        }

    if not account_url or not container:
        raise PublishError("Azure account URL and container are required unless dry_run=True")

    deps = _load_azure()
    credential = deps["DefaultAzureCredential"]()
    service = deps["BlobServiceClient"](account_url=account_url, credential=credential)
    container_client = service.get_container_client(container)

    uploaded = []
    for plan in planned:
        blob_client = _blob_client(service, container_client, container, plan["blob_path"])
        if _blob_metadata_sha256(blob_client, deps) == plan["sha256"]:
            plan = {**plan, "status": "skipped"}
        else:
            with Path(plan["local_path"]).open("rb") as f:
                blob_client.upload_blob(
                    f,
                    overwrite=True,
                    metadata={"sha256": plan["sha256"], "producer": "variantfeatures"},
                    content_settings=deps["ContentSettings"](
                        content_type="application/x-sqlite3"
                    ),
                )
            plan = {**plan, "status": "uploaded"}
        uploaded.append(plan)

    manifest_client = _blob_client(service, container_client, container, manifest_blob_path)
    manifest_sha = sha256_bytes(manifest_bytes)
    if _blob_metadata_sha256(manifest_client, deps) == manifest_sha:
        manifest_status = "skipped"
    else:
        manifest_client.upload_blob(
            manifest_bytes,
            overwrite=True,
            metadata={"sha256": manifest_sha, "producer": "variantfeatures"},
            content_settings=deps["ContentSettings"](content_type="application/json"),
        )
        manifest_status = "uploaded"

    latest_result = _merge_latest_json(
        service,
        container_client,
        container,
        latest_blob_path,
        latest_entries,
        deps,
    )

    return {
        "dry_run": False,
        "release": release,
        "artifacts": uploaded,
        "manifest": {
            "blob_path": manifest_blob_path,
            "sha256": manifest_sha,
            "bytes": len(manifest_bytes),
            "status": manifest_status,
        },
        "latest_blob_path": latest_blob_path,
        **latest_result,
    }


def get_git_sha() -> str | None:
    """Return the current short git SHA, falling back to GIT_SHA."""

    try:
        result = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
        )
        sha = result.stdout.strip()
        if sha:
            return sha
    except (OSError, subprocess.CalledProcessError):
        pass
    return os.environ.get("GIT_SHA")


def release_name(built_at: str, git_sha: str | None) -> str:
    """Return the version directory name: YYYYMMDD-HHMM__gitsha7."""

    dt = parse_built_at(built_at)
    short_sha = (git_sha or "unknown")[:7]
    return f"{dt.strftime('%Y%m%d-%H%M')}__{short_sha}"


def parse_built_at(value: str) -> datetime:
    """Parse an ISO-8601 timestamp and normalize it to UTC."""

    raw = value.strip()
    if raw.endswith("Z"):
        raw = f"{raw[:-1]}+00:00"
    try:
        dt = datetime.fromisoformat(raw)
    except ValueError as e:
        raise PublishError(f"Invalid --built-at timestamp: {value!r}") from e
    if dt.tzinfo is None:
        dt = dt.replace(tzinfo=timezone.utc)
    return dt.astimezone(timezone.utc)


def format_built_at(dt: datetime) -> str:
    """Format a UTC ISO-8601 timestamp for manifests and release paths."""

    if dt.tzinfo is None:
        dt = dt.replace(tzinfo=timezone.utc)
    return dt.astimezone(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _copy_gene_table(conn: sqlite3.Connection, table: str, gene: str) -> None:
    where, params = _copy_where(table, gene)
    _insert_select(conn, table, where, params)


def _copy_where(table: str, gene: str) -> tuple[str, list[Any]]:
    if table == "genes":
        return "UPPER(symbol) = ?", [gene]
    if table == "gene_constraint":
        return "UPPER(gene_symbol) = ?", [gene]
    if table == "transcripts":
        return (
            """
            UPPER(gene_symbol) = ?
            OR transcript_id IN (
                SELECT transcript_id
                FROM src.variant_consequences
                WHERE UPPER(gene_symbol) = ?
            )
            """,
            [gene, gene],
        )
    if table == "variants":
        return "id IN (SELECT id FROM gene_variant_ids)", []
    if table == "variant_consequences":
        return "UPPER(gene_symbol) = ?", [gene]
    if table == "variant_aliases":
        return "variant_id IN (SELECT id FROM gene_variant_ids)", []
    if table == "penetrance_estimates":
        return "variant_id IN (SELECT id FROM gene_variant_ids)", []
    if table == "frameshift_nonsense_mappings":
        return (
            "UPPER(gene_symbol) = ? AND proxy_variant_id IN (SELECT id FROM gene_variant_ids)",
            [gene],
        )
    if table.startswith("annotations_"):
        return "variant_id IN (SELECT id FROM gene_variant_ids)", []
    raise PublishError(f"No gene-scope copy rule for table {table!r}")


def _insert_select(
    conn: sqlite3.Connection,
    table: str,
    where: str,
    params: list[Any],
) -> None:
    src_cols = set(_table_columns(conn, "src", table))
    dest_cols = [col for col in _table_columns(conn, "main", table) if col in src_cols]
    if not dest_cols:
        return
    columns = ", ".join(_quote_ident(col) for col in dest_cols)
    conn.execute(
        f"""
        INSERT OR IGNORE INTO {_quote_ident(table)} ({columns})
        SELECT {columns}
        FROM src.{_quote_ident(table)}
        WHERE {where}
        """,
        params,
    )


def _add_gnomad_r4_joint_population_rows(conn: sqlite3.Connection, gene: str) -> int:
    """Add v4.1 all-sample joint AC/AN rows to a gene slice.

    The source warehouse may already carry exome/genome rows, but joint AN is not
    always recoverable by summing those rows when a variant is absent from genomes.
    Use the gnomAD gene GraphQL endpoint's explicit ``variant.joint`` fields and
    leave any existing joint rows untouched.
    """

    if not _table_exists(conn, "main", "annotations_population"):
        return 0
    source_v4_rows = _count(
        conn,
        "annotations_population",
        "dataset IN (?, ?)",
        [GNOMAD_R4_EXOME_DATASET, GNOMAD_R4_GENOME_DATASET],
    )
    if not source_v4_rows:
        return 0

    before = _count(
        conn,
        "annotations_population",
        "dataset = ? AND pop = 'all'",
        [GNOMAD_R4_JOINT_DATASET],
    )
    variants = {
        f"{row['chromosome']}-{row['position']}-{row['ref']}-{row['alt']}": row["id"]
        for row in conn.execute(
            """
            SELECT id, chromosome, position, ref, alt
            FROM variants
            WHERE chromosome IS NOT NULL
              AND position IS NOT NULL
              AND ref IS NOT NULL
              AND alt IS NOT NULL
            """
        )
    }
    if not variants:
        return 0

    joint_rows = _fetch_gnomad_gene_joint_rows(gene)
    rows = []
    for variant_key, variant_id in variants.items():
        joint = joint_rows.get(variant_key)
        if not joint:
            continue
        an = joint.get("an")
        ac = joint.get("ac")
        if an is None or an <= 0 or ac is None:
            continue
        rows.append(
            (
                variant_id,
                GNOMAD_R4_JOINT_DATASET,
                "all",
                ac / an,
                ac,
                an,
                joint.get("n_homozygotes"),
                joint.get("filter_status"),
                GNOMAD_R4_JOINT_SOURCE,
            )
        )

    conn.executemany(
        """
        INSERT OR IGNORE INTO annotations_population
            (variant_id, dataset, pop, af, ac, an, n_homozygotes, filter_status, source)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        rows,
    )
    after = _count(
        conn,
        "annotations_population",
        "dataset = ? AND pop = 'all'",
        [GNOMAD_R4_JOINT_DATASET],
    )
    return after - before


def _fetch_gnomad_gene_joint_rows(gene: str) -> dict[str, dict[str, Any]]:
    """Fetch gnomAD r4/v4.1 gene variants and return explicit joint rows by key."""

    import requests

    query = """
    query GeneJoint($symbol: String!, $dataset: DatasetId!) {
      gene(gene_symbol: $symbol, reference_genome: GRCh38) {
        variants(dataset: $dataset) {
          variant_id
          exome { ac an }
          genome { ac an }
          joint {
            ac
            an
            homozygote_count
            filters
          }
        }
      }
    }
    """
    response = requests.post(
        GNOMAD_API,
        json={"query": query, "variables": {"symbol": gene.upper(), "dataset": "gnomad_r4"}},
        timeout=120,
    )
    response.raise_for_status()
    payload = response.json()
    if payload.get("errors"):
        messages = "; ".join(error.get("message", "?") for error in payload["errors"])
        raise PublishError(f"gnomAD API error while fetching {gene} joint rows: {messages}")

    gene_payload = (payload.get("data") or {}).get("gene")
    if not gene_payload:
        raise PublishError(f"gnomAD API returned no gene payload for {gene}")

    rows = {}
    mismatches = []
    for variant in gene_payload.get("variants") or []:
        variant_id = variant.get("variant_id")
        joint = variant.get("joint") or {}
        if not variant_id or joint.get("ac") is None or joint.get("an") is None:
            continue
        exome = variant.get("exome") or {}
        genome = variant.get("genome") or {}
        expected_ac = (exome.get("ac") or 0) + (genome.get("ac") or 0)
        if joint["ac"] != expected_ac:
            mismatches.append((variant_id, joint["ac"], expected_ac))
            continue
        filters = joint.get("filters") or []
        rows[variant_id] = {
            "ac": int(joint["ac"]),
            "an": int(joint["an"]),
            "n_homozygotes": joint.get("homozygote_count"),
            "filter_status": ",".join(filters) if filters else "PASS",
        }

    if mismatches:
        first = ", ".join(
            f"{vid}: joint={got} exome+genome={want}"
            for vid, got, want in mismatches[:5]
        )
        raise PublishError(
            f"gnomAD joint AC cross-check failed for {len(mismatches)} {gene} variants: {first}"
        )
    return rows


def _assert_slice_counts(
    conn: sqlite3.Connection,
    gene: str,
    tables: list[str],
) -> dict[str, int]:
    counts = {}
    for table in tables:
        if not _table_exists(conn, "src", table) or not _table_exists(conn, "main", table):
            continue
        where, params = _copy_where(table, gene)
        source_count = _count(conn, f"src.{_quote_ident(table)}", where, params)
        slice_count = _count(conn, _quote_ident(table), "1 = 1", [])
        if source_count != slice_count:
            raise AssertionError(
                f"Slice count mismatch for {table}: source={source_count}, slice={slice_count}"
            )
        counts[table] = slice_count
    return counts


def _count(conn: sqlite3.Connection, table_expr: str, where: str, params: list[Any]) -> int:
    row = conn.execute(f"SELECT COUNT(*) FROM {table_expr} WHERE {where}", params).fetchone()
    return int(row[0] or 0)


def _annotation_tables(conn: sqlite3.Connection) -> list[str]:
    rows = conn.execute(
        """
        SELECT name
        FROM main.sqlite_master
        WHERE type = 'table' AND name LIKE 'annotations\\_%' ESCAPE '\\'
        ORDER BY name
        """
    ).fetchall()
    return [row["name"] for row in rows]


def _table_exists(conn: sqlite3.Connection, schema: str, table: str) -> bool:
    row = conn.execute(
        f"SELECT 1 FROM {schema}.sqlite_master WHERE type = 'table' AND name = ?",
        [table],
    ).fetchone()
    return row is not None


def _table_columns(conn: sqlite3.Connection, schema: str, table: str) -> list[str]:
    return [row[1] for row in conn.execute(f"PRAGMA {schema}.table_info({_quote_ident(table)})")]


def _quote_ident(value: str) -> str:
    return '"' + value.replace('"', '""') + '"'


def _open_read_conn(db: Path | sqlite3.Connection | Any) -> tuple[sqlite3.Connection, bool]:
    if isinstance(db, sqlite3.Connection):
        db.row_factory = sqlite3.Row
        return db, False
    if hasattr(db, "conn"):
        db.conn.row_factory = sqlite3.Row
        return db.conn, False

    conn = sqlite3.connect(Path(db))
    conn.row_factory = sqlite3.Row
    return conn, True


def _gene_row_counts(conn: sqlite3.Connection, gene: str) -> dict[str, int]:
    tables = [
        table
        for table in [
            *DATA_TABLES,
            *_manifest_annotation_tables(conn),
        ]
        if _table_exists(conn, "main", table)
    ]
    return {table: _source_gene_count(conn, table, gene.upper()) for table in tables}


def _manifest_annotation_tables(conn: sqlite3.Connection) -> list[str]:
    rows = conn.execute(
        """
        SELECT name
        FROM sqlite_master
        WHERE type = 'table' AND name LIKE 'annotations\\_%' ESCAPE '\\'
        ORDER BY name
        """
    ).fetchall()
    return [row["name"] for row in rows]


def _source_gene_count(conn: sqlite3.Connection, table: str, gene: str) -> int:
    if table == "genes":
        return _count(conn, _quote_ident(table), "UPPER(symbol) = ?", [gene])
    if table == "gene_constraint":
        return _count(conn, _quote_ident(table), "UPPER(gene_symbol) = ?", [gene])
    if table == "transcripts":
        return _count(
            conn,
            _quote_ident(table),
            """
            UPPER(gene_symbol) = ?
            OR transcript_id IN (
                SELECT transcript_id
                FROM variant_consequences
                WHERE UPPER(gene_symbol) = ?
            )
            """,
            [gene, gene],
        )
    if table == "variant_consequences":
        return _count(conn, _quote_ident(table), "UPPER(gene_symbol) = ?", [gene])
    if table == "variants":
        return _count(
            conn,
            _quote_ident(table),
            """
            id IN (
                SELECT DISTINCT variant_id
                FROM variant_consequences
                WHERE UPPER(gene_symbol) = ?
            )
            """,
            [gene],
        )
    if table == "variant_aliases":
        return _count(conn, _quote_ident(table), _variant_id_where(), [gene])
    if table == "penetrance_estimates":
        return _count(conn, _quote_ident(table), _variant_id_where(), [gene])
    if table == "frameshift_nonsense_mappings":
        return _count(
            conn,
            _quote_ident(table),
            """
            UPPER(gene_symbol) = ?
            AND proxy_variant_id IN (
                SELECT DISTINCT variant_id
                FROM variant_consequences
                WHERE UPPER(gene_symbol) = ?
            )
            """,
            [gene, gene],
        )
    if table.startswith("annotations_"):
        return _count(conn, _quote_ident(table), _variant_id_where(), [gene])
    return 0


def _variant_id_where() -> str:
    return """
        variant_id IN (
            SELECT DISTINCT variant_id
            FROM variant_consequences
            WHERE UPPER(gene_symbol) = ?
        )
    """


def _read_schema_version(conn: sqlite3.Connection) -> int | str | None:
    user_version = conn.execute("PRAGMA user_version").fetchone()[0]
    if user_version:
        return int(user_version)

    for table, key_col, value_col in (
        ("metadata", "key", "value"),
        ("schema_metadata", "key", "value"),
    ):
        if not _table_exists(conn, "main", table):
            continue
        cols = set(_table_columns(conn, "main", table))
        if {key_col, value_col}.issubset(cols):
            row = conn.execute(
                f"""
                SELECT {_quote_ident(value_col)}
                FROM {_quote_ident(table)}
                WHERE {_quote_ident(key_col)} = 'schema_version'
                LIMIT 1
                """
            ).fetchone()
            if row and row[0] is not None:
                return row[0]

    if _table_exists(conn, "main", "schema_version"):
        cols = set(_table_columns(conn, "main", "schema_version"))
        for col in ("version", "schema_version"):
            if col in cols:
                row = conn.execute(
                    f"SELECT {_quote_ident(col)} FROM schema_version ORDER BY rowid DESC LIMIT 1"
                ).fetchone()
                if row and row[0] is not None:
                    return row[0]

    return None


def _source_versions(conn: sqlite3.Connection) -> dict[str, list[dict[str, Any]]]:
    sources: dict[str, list[dict[str, Any]]] = {}
    if _table_exists(conn, "main", "annotations_pathogenicity"):
        sources["pathogenicity"] = _grouped_rows(
            conn,
            """
            SELECT predictor, predictor_version AS version, source, COUNT(*) AS rows
            FROM annotations_pathogenicity
            GROUP BY predictor, predictor_version, source
            ORDER BY predictor, predictor_version, source
            """,
        )
    if _table_exists(conn, "main", "annotations_population"):
        sources["population"] = _grouped_rows(
            conn,
            """
            SELECT dataset, source, COUNT(*) AS rows
            FROM annotations_population
            GROUP BY dataset, source
            ORDER BY dataset, source
            """,
        )
    if _table_exists(conn, "main", "annotations_clinical"):
        sources["clinical"] = _grouped_rows(
            conn,
            """
            SELECT source, COUNT(*) AS rows
            FROM annotations_clinical
            GROUP BY source
            ORDER BY source
            """,
        )
    if _table_exists(conn, "main", "annotations_splice"):
        sources["splice"] = _grouped_rows(
            conn,
            """
            SELECT predictor, predictor_version AS version, source, COUNT(*) AS rows
            FROM annotations_splice
            GROUP BY predictor, predictor_version, source
            ORDER BY predictor, predictor_version, source
            """,
        )
    if _table_exists(conn, "main", "annotations_structure"):
        sources["structure"] = _grouped_rows(
            conn,
            """
            SELECT feature, feature_version AS version, source, COUNT(*) AS rows
            FROM annotations_structure
            GROUP BY feature, feature_version, source
            ORDER BY feature, feature_version, source
            """,
        )
    if _table_exists(conn, "main", "gene_constraint"):
        sources["gene_constraint"] = _grouped_rows(
            conn,
            """
            SELECT dataset, source, COUNT(*) AS rows
            FROM gene_constraint
            GROUP BY dataset, source
            ORDER BY dataset, source
            """,
        )
    return sources


def _grouped_rows(conn: sqlite3.Connection, sql: str) -> list[dict[str, Any]]:
    return [dict(row) for row in conn.execute(sql).fetchall()]


def _normalize_artifacts(
    slices: Mapping[str, Path] | list[PublishArtifact] | tuple[PublishArtifact, ...],
) -> list[PublishArtifact]:
    if isinstance(slices, Mapping):
        artifacts = []
        for key, value in slices.items():
            path = Path(value)
            if key in {"__full_db__", "full_db"}:
                artifacts.append(PublishArtifact(None, path, "full_db"))
            else:
                artifacts.append(PublishArtifact(str(key).upper(), path, "gene_slice"))
        return artifacts

    artifacts = []
    for item in slices:
        if isinstance(item, PublishArtifact):
            artifacts.append(item)
        elif isinstance(item, Mapping):
            gene = item.get("gene")
            kind = str(item.get("kind") or ("gene_slice" if gene else "full_db"))
            artifacts.append(PublishArtifact(str(gene).upper() if gene else None, Path(item["path"]), kind))
        else:
            gene, path = item
            artifacts.append(PublishArtifact(str(gene).upper(), Path(path), "gene_slice"))
    return artifacts


def _artifact_relative_path(artifact: PublishArtifact) -> str:
    if artifact.kind == "gene_slice":
        if not artifact.gene:
            raise PublishError("Gene slice artifact is missing a gene")
        return f"genes/{artifact.gene}.db"
    if artifact.kind == "full_db":
        return "variants.db"
    raise PublishError(f"Unknown artifact kind: {artifact.kind!r}")


def _planned_artifact(prefix: str, release: str, artifact: PublishArtifact) -> dict[str, Any]:
    relative = _artifact_relative_path(artifact)
    return {
        "gene": artifact.gene,
        "kind": artifact.kind,
        "local_path": str(artifact.path),
        "blob_path": f"{prefix}/{release}/{relative}",
        "sha256": sha256_file(artifact.path),
        "bytes": artifact.path.stat().st_size,
    }


def _latest_entries(
    prefix: str,
    release: str,
    built_at: str,
    manifest: dict[str, Any],
    artifacts: list[PublishArtifact],
) -> dict[str, dict[str, Any]]:
    schema_version = manifest.get("schema_version")
    entries = {}
    for artifact in artifacts:
        if artifact.kind != "gene_slice" or not artifact.gene:
            continue
        relative = _artifact_relative_path(artifact)
        entries[artifact.gene] = {
            "path": f"{prefix}/{release}/{relative}",
            "sha256": sha256_file(artifact.path),
            "built_at": built_at,
            "schema_version": schema_version,
        }
    return entries


def _latest_diff(
    before: dict[str, Any],
    after: dict[str, Any],
    changed: dict[str, Any],
) -> dict[str, dict[str, Any]]:
    return {
        gene: {
            "before": before.get(gene),
            "after": after.get(gene),
        }
        for gene in sorted(changed)
    }


def _manifest_bytes(manifest: dict[str, Any]) -> bytes:
    return (json.dumps(manifest, indent=2, sort_keys=True) + "\n").encode("utf-8")


def _clean_prefix(prefix: str) -> str:
    return prefix.strip("/")


def _load_azure() -> dict[str, Any]:
    global BlobServiceClient
    global ContentSettings
    global DefaultAzureCredential
    global MatchConditions
    global ResourceExistsError
    global ResourceModifiedError
    global ResourceNotFoundError

    if BlobServiceClient is None:
        try:
            from azure.core import MatchConditions as AzureMatchConditions
            from azure.core.exceptions import (
                ResourceExistsError as AzureResourceExistsError,
            )
            from azure.core.exceptions import (
                ResourceModifiedError as AzureResourceModifiedError,
            )
            from azure.core.exceptions import (
                ResourceNotFoundError as AzureResourceNotFoundError,
            )
            from azure.identity import DefaultAzureCredential as AzureDefaultAzureCredential
            from azure.storage.blob import BlobServiceClient as AzureBlobServiceClient
            from azure.storage.blob import ContentSettings as AzureContentSettings
        except ImportError as e:
            raise PublishError(
                "Azure publishing requires optional dependencies. "
                "Install with `pip install variantfeatures[azure]`."
            ) from e

        BlobServiceClient = AzureBlobServiceClient
        ContentSettings = AzureContentSettings
        DefaultAzureCredential = AzureDefaultAzureCredential
        MatchConditions = AzureMatchConditions
        ResourceExistsError = AzureResourceExistsError
        ResourceModifiedError = AzureResourceModifiedError
        ResourceNotFoundError = AzureResourceNotFoundError

    return {
        "BlobServiceClient": BlobServiceClient,
        "ContentSettings": ContentSettings,
        "DefaultAzureCredential": DefaultAzureCredential,
        "MatchConditions": MatchConditions,
        "ResourceExistsError": ResourceExistsError,
        "ResourceModifiedError": ResourceModifiedError,
        "ResourceNotFoundError": ResourceNotFoundError,
    }


def _blob_client(service: Any, container_client: Any, container: str, blob_path: str) -> Any:
    if hasattr(container_client, "get_blob_client"):
        return container_client.get_blob_client(blob_path)
    return service.get_blob_client(container=container, blob=blob_path)


def _blob_metadata_sha256(blob_client: Any, deps: dict[str, Any]) -> str | None:
    try:
        props = blob_client.get_blob_properties()
    except Exception as e:
        if _is_azure_exception(e, deps["ResourceNotFoundError"]):
            return None
        raise
    metadata = _prop(props, "metadata") or {}
    return metadata.get("sha256")


def _merge_latest_json(
    service: Any,
    container_client: Any,
    container: str,
    latest_blob_path: str,
    latest_entries: dict[str, dict[str, Any]],
    deps: dict[str, Any],
) -> dict[str, Any]:
    latest_client = _blob_client(service, container_client, container, latest_blob_path)
    for _attempt in range(3):
        before, etag = _read_latest(latest_client, deps)
        after = {**before, **latest_entries}
        latest_bytes = (json.dumps(after, indent=2, sort_keys=True) + "\n").encode("utf-8")
        kwargs = {
            "metadata": {
                "sha256": sha256_bytes(latest_bytes),
                "producer": "variantfeatures",
            },
            "content_settings": deps["ContentSettings"](content_type="application/json"),
        }
        if etag:
            kwargs.update(
                {
                    "overwrite": True,
                    "etag": etag,
                    "match_condition": deps["MatchConditions"].IfNotModified,
                }
            )
        else:
            kwargs.update({"overwrite": False})
        try:
            latest_client.upload_blob(latest_bytes, **kwargs)
        except Exception as e:
            if _is_azure_exception(e, deps["ResourceModifiedError"]) or _is_azure_exception(
                e, deps["ResourceExistsError"]
            ):
                continue
            raise
        return {
            "latest_before": before,
            "latest_after": after,
            "latest_diff": _latest_diff(before, after, latest_entries),
        }

    raise PublishError(
        "Could not update latest.json atomically after concurrent modifications; retry publish"
    )


def _read_latest(blob_client: Any, deps: dict[str, Any]) -> tuple[dict[str, Any], str | None]:
    try:
        props = blob_client.get_blob_properties()
        data = blob_client.download_blob().readall()
    except Exception as e:
        if _is_azure_exception(e, deps["ResourceNotFoundError"]):
            return {}, None
        raise
    if not data:
        return {}, _prop(props, "etag")
    return json.loads(data.decode("utf-8")), _prop(props, "etag")


def _prop(obj: Any, name: str) -> Any:
    if isinstance(obj, Mapping):
        return obj.get(name)
    return getattr(obj, name, None)


def _is_azure_exception(error: Exception, expected_type: Any) -> bool:
    if expected_type is not None and isinstance(error, expected_type):
        return True
    return error.__class__.__name__ == getattr(expected_type, "__name__", "")
