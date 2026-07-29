"""Command-line interface for VariantFeatures."""

from __future__ import annotations

import click
import csv
import json
import os
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path

import requests

from .database import DEFAULT_DB, SCHEMA_VERSION, SchemaVersionError, VariantDB
from .local_storage import (
    ALLOW_DB_CREATE_ENV,
    ALLOW_LOCAL_ENV,
    LocalStorageError,
    REPO_ROOT,
    external_link_target,
    external_path_state,
    external_volume_name,
)


class _StorageAwareGroup(click.Group):
    """Report storage and schema problems as a clean error, not a traceback.

    An unmounted volume, an absent database, or a database at the wrong schema
    version is an operator condition with a written remedy, so it should read
    like `Error: ...` and exit non-zero rather than dumping a stack through
    `VariantDB.__init__`.
    """

    def invoke(self, ctx):
        try:
            return super().invoke(ctx)
        except (LocalStorageError, SchemaVersionError) as e:
            raise click.ClickException(str(e)) from e


@click.group(cls=_StorageAwareGroup)
@click.version_option()
def main():
    """Aggregate predictive features for genetic variants."""
    pass


@main.command()
@click.option("--gene", "gene_opts", "-g", multiple=True, help="Gene symbol. Repeat or use --genes.")
@click.option("--genes", default=None, help="Comma-separated gene symbols (legacy-compatible)")
@click.option("--db", type=click.Path(), default=None, help="Database path")
@click.option(
    "--sources",
    default="core",
    show_default=True,
    help=(
        "Normalized sources or groups: core, all, pathogenicity, population, "
        "clinical, splice, expression, structure, alphamissense, revel, cadd, gnomad, ..."
    ),
)
@click.option(
    "--types",
    default="missense,nonsense,stop_lost,start_lost,synonymous",
    show_default=True,
    help="Variant types to enumerate before annotation.",
)
@click.option("--limit", type=int, default=None, help="Stop after N enumerated variants per transcript")
@click.option("--annotation-limit", type=int, default=None, help="Max jobs per runnable source")
@click.option("--no-run", is_flag=True, help="Only enumerate and queue jobs; do not run annotators")
@click.option("--strict", is_flag=True, help="Fail instead of skipping unavailable annotators")
@click.option(
    "--pext-bigwig-dir",
    type=click.Path(file_okay=False),
    default="data/pext/ucsc_hg38",
    show_default=True,
    help="Local gnomAD pext bigWig directory for the expression source.",
)
@click.option("--pext-dataset", default="ucsc_gnomad_pext_hg38", show_default=True)
@click.option(
    "--frameshift-max-steps",
    type=int,
    default=20,
    show_default=True,
    help="Amino-acid search radius for frameshift_proxy mappings.",
)
@click.option(
    "--isoforms",
    type=click.Choice(["canonical", "mane", "all"]),
    default="canonical",
    show_default=True,
    help="Transcripts to enumerate before annotation.",
)
@click.option(
    "--max-isoforms",
    type=int,
    default=None,
    help="Cap selected isoforms after MANE/canonical/APPRIS ranking.",
)
def build(
    gene_opts: tuple,
    genes: str,
    db: str,
    sources: str,
    types: str,
    limit: int,
    annotation_limit: int,
    no_run: bool,
    strict: bool,
    pext_bigwig_dir: str,
    pext_dataset: str,
    frameshift_max_steps: int,
    isoforms: str,
    max_isoforms: int,
):
    """Build/update the normalized database for one or more genes."""
    gene_list = _parse_gene_options(gene_opts, genes)
    if not gene_list:
        raise click.UsageError("Pass at least one --gene/-g or --genes value.")

    from .build import BuildError, build_gene, parse_sources
    from .transcripts import TranscriptError

    db_path = Path(db) if db else None
    selected_sources = parse_sources(sources)
    click.echo(f"Building normalized database for genes: {', '.join(gene_list)}")
    click.echo(f"Sources: {', '.join(sorted(selected_sources))}")
    click.echo(f"Isoforms: {isoforms}" + (f" (max {max_isoforms})" if max_isoforms else ""))
    if no_run:
        click.echo("Mode: enumerate + queue only")

    final_db_path = db_path
    for symbol in gene_list:
        click.echo(f"\n{'=' * 60}")
        click.echo(f"Processing {symbol}")
        click.echo(f"{'=' * 60}")
        try:
            result = build_gene(
                symbol,
                db_path=db_path,
                sources=sources,
                types=types,
                enumerate_limit=limit,
                run_annotations=not no_run,
                annotation_limit=annotation_limit,
                strict=strict,
                pext_bigwig_dir=pext_bigwig_dir,
                pext_dataset=pext_dataset,
                frameshift_max_steps=frameshift_max_steps,
                isoforms=isoforms,
                max_isoforms=max_isoforms,
            )
        except (BuildError, TranscriptError, ValueError) as e:
            click.echo(f"Error: {e}", err=True)
            sys.exit(1)

        final_db_path = result.db_path
        enum = result.enumeration
        if len(result.transcript_labels) == 1:
            click.echo(f"Transcript: {result.transcript_label}")
        else:
            click.echo(f"Transcripts: {len(result.transcript_labels)}")
            for label in result.transcript_labels:
                tx_summary = enum.get("by_transcript", {}).get(label, {})
                click.echo(
                    f"  {label:<24} variants={tx_summary.get('variants', 0)} "
                    f"consequences={tx_summary.get('consequences', 0)}"
                )
        click.echo(f"Variants inserted/updated: {enum['variants']}")
        click.echo(f"Consequences written:      {enum['consequences']}")
        click.echo(f"Annotation jobs queued:    {enum['jobs_queued']}")
        if enum["by_consequence"]:
            click.echo("By consequence:")
            for cons, n in sorted(enum["by_consequence"].items(), key=lambda x: -x[1]):
                click.echo(f"  {cons:<24} {n}")

        if result.sources:
            click.echo("Sources:")
            for run in result.sources:
                if run.status == "done":
                    summary = " ".join(f"{k}={v}" for k, v in run.summary.items())
                    click.echo(f"  {run.source:<16} done     {summary}")
                else:
                    click.echo(f"  {run.source:<16} skipped  {run.error}")

    click.echo(f"\nDatabase: {final_db_path or VariantDB().db_path}")


def _parse_gene_options(gene_opts: tuple, genes: str | None) -> list[str]:
    values: list[str] = []
    for item in gene_opts:
        values.extend(str(item).split(","))
    if genes:
        values.extend(genes.split(","))
    out = []
    seen = set()
    for value in values:
        symbol = value.strip().upper()
        if not symbol or symbol in seen:
            continue
        seen.add(symbol)
        out.append(symbol)
    return out


def _open_db(db: str | None, *, read_only: bool = False) -> VariantDB:
    db_path = Path(db) if db else None
    return VariantDB(db_path, initialize=not read_only, read_only=read_only)



@main.command()
@click.option("--gene", "-g", required=True, help="Gene symbol")
@click.option("--db", type=click.Path(), default=None, help="Database path")
@click.option("--format", "fmt", type=click.Choice(["csv", "json", "table"]), default="table")
def query(gene: str, db: str, fmt: str):
    """Query variants for a gene."""
    vdb = _open_db(db, read_only=True)

    from .normalized_export import build_wide_rows

    variants, fieldnames = build_wide_rows(vdb, gene.upper())

    if not variants:
        click.echo(f"No normalized variants found for {gene}")
        return

    if fmt == "json":
        click.echo(json.dumps(variants, indent=2, default=str))
    elif fmt == "csv":
        writer = csv.DictWriter(sys.stdout, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(variants)
    else:
        _print_normalized_query_table(gene, variants, fieldnames)



def _print_normalized_query_table(gene: str, variants: list[dict], fieldnames: list[str]) -> None:
    click.echo(f"Found {len(variants)} normalized variants for {gene}")
    feature_cols = [c for c in fieldnames if "." in c and not c.startswith("aliases.")]
    if feature_cols:
        click.echo(f"Feature columns available: {len(feature_cols)} (use --format csv/json for all)")
    click.echo()
    click.echo(f"{'VARIANT':<23} {'CONSEQUENCE':<22} {'ClinVar':<18} {'AM':>7} {'REVEL':>7} {'gnomAD AF':>12}")
    click.echo("-" * 95)

    for row in variants[:50]:
        variant = row.get("hgvs_p") or row.get("hgvs_c") or _format_variant_coord(row)
        consequence = row.get("consequence") or "N/A"
        clinvar = _first_feature_value(row, "clinical.", ".classification") or "N/A"
        am = _first_feature_value(row, "pathogenicity.alphamissense", ".score")
        revel = _first_feature_value(row, "pathogenicity.revel", ".score")
        gnomad = (
            row.get("population.gnomad_exomes_v4.all.af")
            if row.get("population.gnomad_exomes_v4.all.af") not in (None, "")
            else row.get("population.gnomad_genomes_v4.all.af")
        )
        if gnomad in (None, ""):
            gnomad = _first_feature_value(row, "population.", ".all.af")

        click.echo(
            f"{str(variant)[:22]:<23} "
            f"{str(consequence)[:21]:<22} "
            f"{str(clinvar)[:17]:<18} "
            f"{_format_score(am):>7} "
            f"{_format_score(revel):>7} "
            f"{_format_af(gnomad):>12}"
        )

    if len(variants) > 50:
        click.echo(f"\n... and {len(variants) - 50} more variants (use --format csv for full output)")


def _first_feature_value(row: dict, prefix: str, suffix: str):
    for key in sorted(row):
        value = row.get(key)
        if key.startswith(prefix) and key.endswith(suffix) and value not in (None, ""):
            return value
    return None


def _format_variant_coord(row: dict) -> str:
    return f"{row.get('chromosome')}:{row.get('position')}:{row.get('ref')}>{row.get('alt')}"


def _format_score(value) -> str:
    if value in (None, ""):
        return "N/A"
    try:
        return f"{float(value):.3f}"
    except (TypeError, ValueError):
        return str(value)[:7]


def _format_af(value) -> str:
    if value in (None, ""):
        return "N/A"
    try:
        return f"{float(value):.2e}"
    except (TypeError, ValueError):
        return str(value)[:12]


def _fmt_bytes(n: int) -> str:
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if n < 1024 or unit == "TiB":
            return f"{n:.0f} {unit}" if unit == "B" else f"{n:.1f} {unit}"
        n /= 1024
    return f"{n} B"


@main.command()
@click.option("--db", type=click.Path(), default=None, help="Database path to check")
def doctor(db: str):
    """Report whether storage and the database are ready to use.

    Read-only and never raises: the point is to say what state things are in,
    especially on a fresh clone where nothing is set up yet. Exits non-zero when
    something needed is missing, so it can gate a script.
    """
    problems: list[str] = []
    click.echo("VariantFeatures environment")
    click.echo("=" * 60)
    click.echo(f"Repo root:  {REPO_ROOT}")
    click.echo(f"Python:     {sys.version.split()[0]}")
    click.echo("")

    click.echo("Storage")
    for name in ("data", "annovar/humandb"):
        state = external_path_state(name)
        target = external_link_target(name)
        via = f" -> {target}" if target else ""
        if state == "linked":
            click.echo(f"  ok       {name}{via}")
        elif state == "dangling":
            volume = external_volume_name(name)
            detail = f"volume {volume!r} is not mounted" if volume else "target missing"
            click.echo(f"  MISSING  {name}{via}  ({detail})")
            problems.append(f"{name}: {detail}")
        else:
            click.echo(f"  absent   {name}  (not set up in this checkout)")
            # annovar/humandb is optional: without it only the ANNOVAR source is off.
            if name == "data":
                problems.append(f"{name}: not set up in this checkout")
    click.echo("")

    click.echo("Database")
    db_path = Path(db) if db else DEFAULT_DB
    click.echo(f"  path     {db_path}")
    if not db_path.is_file():
        click.echo("  MISSING  no database file at that path")
        problems.append("variants.db: not present")
    else:
        size = db_path.stat().st_size
        if size == 0:
            click.echo("  MISSING  file is empty (0 bytes)")
            problems.append("variants.db: empty file")
        else:
            click.echo(f"  size     {_fmt_bytes(size)}")
            try:
                vdb = VariantDB(db_path, initialize=False, read_only=True)
            except (LocalStorageError, SchemaVersionError) as e:
                click.echo(f"  UNUSABLE {e}".splitlines()[0])
                problems.append(f"variants.db: {str(e).splitlines()[0]}")
            else:
                try:
                    stamped = vdb.schema_version
                    click.echo(
                        f"  schema   version {stamped}"
                        if stamped
                        else f"  schema   unstamped (adopted as {SCHEMA_VERSION})"
                    )
                    n = vdb.conn.execute("SELECT COUNT(*) FROM variants").fetchone()[0]
                    click.echo(f"  variants {n:,}")
                    genes = vdb.conn.execute("SELECT COUNT(*) FROM genes").fetchone()[0]
                    click.echo(f"  genes    {genes:,}")
                finally:
                    vdb.conn.close()
    click.echo("")

    click.echo("Overrides in effect")
    any_override = False
    for env, effect in (
        (ALLOW_LOCAL_ENV, "plain local storage is permitted"),
        (ALLOW_DB_CREATE_ENV, "building a database from scratch is permitted"),
    ):
        raw = os.environ.get(env)
        if raw:
            click.echo(f"  {env}={raw}  ({effect})")
            any_override = True
    if not any_override:
        click.echo("  none")
    click.echo("")

    if not problems:
        click.echo("Ready.")
        return
    click.echo(f"Not ready — {len(problems)} thing(s) to fix:")
    for p in problems:
        click.echo(f"  - {p}")
    click.echo("")
    click.echo("See README.md > 'Local Data Storage' for setup, including the")
    click.echo("collaborator path when you have no external volume.")
    raise SystemExit(1)


@main.command()
@click.option("--db", type=click.Path(), default=None, help="Database path")
def stats(db: str):
    """Show database statistics."""
    vdb = _open_db(db, read_only=True)

    _stats_normalized(vdb)

    click.echo(f"\nDatabase: {vdb.db_path}")
    stamped = vdb.schema_version
    click.echo(
        f"Schema version: {stamped}"
        if stamped
        else f"Schema version: unstamped (adopted as {SCHEMA_VERSION})"
    )


def _stats_normalized(vdb: VariantDB) -> None:
    click.echo("VariantFeatures Normalized Database Statistics")
    click.echo("=" * 60)

    cur = vdb.conn.execute("""
        WITH gene_vars AS (
            SELECT UPPER(gene_symbol) AS gene, variant_id
            FROM variant_consequences
            WHERE source = 'enumerated'
              AND gene_symbol IS NOT NULL
              AND gene_symbol != ''
            GROUP BY UPPER(gene_symbol), variant_id
        )
        SELECT gene, COUNT(DISTINCT variant_id) AS total
        FROM gene_vars
        GROUP BY gene
        ORDER BY total DESC
    """)

    rows = {row["gene"]: {"gene": row["gene"], "total": row["total"]} for row in cur.fetchall()}
    click.echo("\nDenominator: variant_consequences.source='enumerated'")
    if not rows:
        click.echo("No normalized enumerated variants found.")
        return

    for field in (
        "pathogenicity", "population", "clinical", "splice",
        "structure", "expression", "conservation", "gene_constraint",
    ):
        for row in rows.values():
            row[field] = 0

    family_tables = {
        "pathogenicity": "annotations_pathogenicity",
        "population": "annotations_population",
        "clinical": "annotations_clinical",
        "splice": "annotations_splice",
        "structure": "annotations_structure",
        "expression": "annotations_expression",
        "conservation": "annotations_conservation",
    }
    for field, table in family_tables.items():
        for row in _normalized_family_counts(vdb, table):
            if row["gene"] in rows:
                rows[row["gene"]][field] = row["variants"]

    for row in _normalized_gene_constraint_counts(vdb):
        if row["gene"] in rows:
            rows[row["gene"]]["gene_constraint"] = row["variants"]

    click.echo(
        f"\n{'Gene':<10} {'Total':>8} {'Path':>8} {'Pop':>8} {'Clin':>8} "
        f"{'Splice':>8} {'Struct':>8} {'Expr':>8} {'Cons':>8} {'GeneC':>8}"
    )
    click.echo("-" * 92)
    sorted_rows = sorted(rows.values(), key=lambda r: (-r["total"], r["gene"]))
    for row in sorted_rows:
        click.echo(
            f"{row['gene']:<10} {row['total']:>8} {row['pathogenicity']:>8} "
            f"{row['population']:>8} {row['clinical']:>8} {row['splice']:>8} "
            f"{row['structure']:>8} {row['expression']:>8} "
            f"{row['conservation']:>8} {row['gene_constraint']:>8}"
        )


def _normalized_family_counts(vdb: VariantDB, table: str) -> list[dict]:
    sql = f"""
        WITH gene_vars AS (
            SELECT UPPER(gene_symbol) AS gene, variant_id
            FROM variant_consequences
            WHERE source = 'enumerated'
              AND gene_symbol IS NOT NULL
              AND gene_symbol != ''
            GROUP BY UPPER(gene_symbol), variant_id
        )
        SELECT gv.gene, COUNT(DISTINCT gv.variant_id) AS variants
        FROM gene_vars gv
        JOIN {table} a ON a.variant_id = gv.variant_id
        GROUP BY gv.gene
    """
    return [dict(row) for row in vdb.conn.execute(sql)]


def _normalized_gene_constraint_counts(vdb: VariantDB) -> list[dict]:
    sql = """
        WITH gene_vars AS (
            SELECT UPPER(gene_symbol) AS gene, variant_id
            FROM variant_consequences
            WHERE source = 'enumerated'
              AND gene_symbol IS NOT NULL
              AND gene_symbol != ''
            GROUP BY UPPER(gene_symbol), variant_id
        )
        SELECT gv.gene, COUNT(DISTINCT gv.variant_id) AS variants
        FROM gene_vars gv
        JOIN gene_constraint gc ON UPPER(gc.gene_symbol) = gv.gene
        GROUP BY gv.gene
    """
    return [dict(row) for row in vdb.conn.execute(sql)]



@main.command()
@click.option("--gene", "-g", required=True, help="Gene symbol")
@click.option("--db", type=click.Path(), default=None, help="Database path")
@click.option("--output", "-o", type=click.Path(), required=True, help="Output CSV file")
@click.option(
    "--layout",
    type=click.Choice(["wide", "transcript-wide", "long"]),
    default="wide",
    show_default=True,
)
@click.option("--groups", default="all", show_default=True,
              help="Feature groups: all,pathogenicity,population,clinical,splice,expression,structure,conservation,gene_constraint")
@click.option("--include-provenance", is_flag=True, help="Include source columns in wide export")
def export(gene: str, db: str, output: str, layout: str, groups: str, include_provenance: bool):
    """Export variants for downstream pipelines."""
    db_path = Path(db) if db else None
    vdb = VariantDB(db_path)

    from .normalized_export import ExportError, export_gene

    try:
        summary = export_gene(
            vdb,
            gene.upper(),
            output,
            layout=layout,
            groups=groups,
            include_provenance=include_provenance,
        )
    except ExportError as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)

    if layout in {"wide", "transcript-wide"}:
        if summary["variants"] == 0:
            click.echo(f"No normalized variants found for {gene}")
            return
        row_label = "variant/transcript rows" if layout == "transcript-wide" else "variants"
        click.echo(
            f"Exported {summary['variants']} normalized {row_label} "
            f"({summary['columns']} columns) to {summary['path']}"
        )
    else:
        if summary["rows"] == 0:
            click.echo(f"No normalized feature rows found for {gene}")
            return
        click.echo(
            f"Exported {summary['rows']} normalized feature rows "
            f"({summary['columns']} columns) to {summary['path']}"
        )


@main.command()
@click.option("--gene", "-g", "gene_opts", multiple=True, help="Gene symbol to publish. Repeatable.")
@click.option(
    "--all-built",
    is_flag=True,
    help="Publish every gene listed in the database's genes table.",
)
@click.option(
    "--db",
    type=click.Path(),
    default=str(DEFAULT_DB),
    show_default=True,
    help="Database path",
)
@click.option("--built-at", default=None, help="ISO-8601 UTC build timestamp")
@click.option("--prefix", default="variantfeatures", show_default=True, help="Azure blob prefix")
@click.option(
    "--include-full-db/--no-full-db",
    default=False,
    show_default=True,
    help="Include the complete SQLite database as an additional artifact.",
)
@click.option("--dry-run", is_flag=True, help="Write local artifacts and print the upload plan.")
def publish(
    gene_opts: tuple[str, ...],
    all_built: bool,
    db: str,
    built_at: str | None,
    prefix: str,
    include_full_db: bool,
    dry_run: bool,
):
    """Publish per-gene SQLite slices and provenance to Azure Blob Storage."""
    from .publish import (
        PublishError,
        build_manifest,
        export_gene_slice,
        format_built_at,
        get_git_sha,
        parse_built_at,
        release_name,
        upload,
    )

    db_path = Path(db)
    if all_built:
        genes = _built_gene_symbols(db_path)
    else:
        genes = _parse_gene_options(gene_opts, None)
        if not genes:
            raise click.UsageError("Pass at least one --gene/-g or use --all-built.")

    built_at_value = (
        format_built_at(parse_built_at(built_at))
        if built_at
        else format_built_at(datetime.now(timezone.utc))
    )
    git_sha = get_git_sha()
    release = release_name(built_at_value, git_sha)
    release_dir = Path("dist") / "publish" / release
    gene_dir = release_dir / "genes"
    gene_dir.mkdir(parents=True, exist_ok=True)

    slices: dict[str, Path] = {}
    for gene in genes:
        out_path = gene_dir / f"{gene}.db"
        summary = export_gene_slice(db_path, gene, out_path)
        slices[gene] = out_path
        click.echo(
            f"Wrote {gene} slice: {out_path} "
            f"({summary['row_counts'].get('variants', 0)} variants)"
        )

    artifacts: dict[str, Path] = dict(slices)
    if include_full_db:
        full_db_path = release_dir / "variants.db"
        shutil.copy2(db_path, full_db_path)
        artifacts["__full_db__"] = full_db_path

    manifest = build_manifest(db_path, genes, built_at_value, artifacts)
    manifest_path = release_dir / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    click.echo(f"Wrote manifest: {manifest_path}")

    account_url = os.environ.get("VF_BLOB_ACCOUNT_URL")
    container = os.environ.get("VF_BLOB_CONTAINER")
    effective_dry_run = dry_run or not account_url or not container
    if effective_dry_run and not dry_run:
        click.echo("Azure is not configured; running as a local dry run.")

    try:
        result = upload(
            account_url,
            container,
            prefix,
            built_at_value,
            git_sha,
            artifacts,
            manifest,
            dry_run=effective_dry_run,
        )
    except PublishError as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)

    action = "WOULD upload" if result["dry_run"] else "Uploaded"
    click.echo()
    for artifact in result["artifacts"]:
        status = artifact.get("status")
        suffix = f" [{status}]" if status else ""
        click.echo(
            f"{action}: {artifact['blob_path']} "
            f"sha256={artifact['sha256']} bytes={artifact['bytes']}{suffix}"
        )
    manifest_status = result["manifest"].get("status")
    manifest_suffix = f" [{manifest_status}]" if manifest_status else ""
    click.echo(
        f"{action}: {result['manifest']['blob_path']} "
        f"sha256={result['manifest']['sha256']} "
        f"bytes={result['manifest']['bytes']}{manifest_suffix}"
    )
    click.echo(f"latest.json: {result['latest_blob_path']}")
    click.echo("latest.json diff:")
    click.echo(json.dumps(result["latest_diff"], indent=2, sort_keys=True))
    click.echo(f"Local artifacts: {release_dir}")


def _built_gene_symbols(db_path: Path) -> list[str]:
    vdb = VariantDB(db_path, initialize=False, read_only=True)
    try:
        rows = vdb.conn.execute(
            "SELECT UPPER(symbol) AS symbol FROM genes ORDER BY UPPER(symbol)"
        ).fetchall()
    finally:
        vdb.close()
    genes = [row["symbol"] for row in rows if row["symbol"]]
    if not genes:
        raise click.UsageError("No genes found in the database's genes table.")
    return genes


@main.command()
@click.argument("identifier")
@click.option("--kind", default=None,
              type=click.Choice([
                  "ca_id", "rsid", "hgvs_g", "hgvs_c", "hgvs_p",
                  "clinvar_vcv", "clinvar_rcv", "clinvar_variation", "gnomad_id",
              ]),
              help="Force input kind (otherwise auto-detected)")
@click.option("--db", type=click.Path(), default=None, help="Database path")
@click.option("--no-persist", is_flag=True, help="Do not write the result to the database")
@click.option("--format", "fmt", type=click.Choice(["table", "json"]), default="table")
def identify(identifier: str, kind: str, db: str, no_persist: bool, fmt: str):
    """Resolve a variant identifier to canonical GRCh38 + all known aliases.

    Accepts: HGVS (g./c.), rsID, ClinGen CA-ID, ClinVar VCV/RCV/variation ID, gnomAD ID.
    """
    from .identity import resolve, persist, IdentityError

    try:
        rv = resolve(identifier, kind=kind)
    except IdentityError as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)
    except requests.RequestException as e:
        click.echo(f"Network error talking to ClinGen Allele Registry: {e}", err=True)
        sys.exit(2)

    if not no_persist:
        vdb = VariantDB(Path(db) if db else None)
        variant_id = persist(rv, vdb)
        if variant_id is None:
            click.echo("Warning: GRCh38 coordinates not available; aliases were not persisted.", err=True)

    if fmt == "json":
        click.echo(json.dumps(rv.to_dict(), indent=2, default=str))
        return

    click.echo(f"CA-ID:        {rv.ca_id or '(none)'}")
    if rv.chromosome and rv.position is not None:
        click.echo(f"GRCh38:       chr{rv.chromosome}:{rv.position} {rv.ref or '-'}>{rv.alt or '-'} ({rv.variant_type or '?'})")
    else:
        click.echo("GRCh38:       (not available)")
    click.echo(f"hgvs_g:       {rv.hgvs_g or '(none)'}")
    click.echo()
    click.echo(f"Aliases ({len(rv.aliases)}):")
    click.echo(f"  {'TYPE':<20} {'VALUE'}")
    click.echo("  " + "-" * 60)
    for a in rv.aliases:
        click.echo(f"  {a['alias_type']:<20} {a['alias_value']}")


@main.command(name="enumerate")
@click.option("--gene", "-g", required=True, help="Gene symbol (e.g. KCNH2)")
@click.option("--types", default="missense,nonsense",
              help="Comma-separated variant types: missense, nonsense, synonymous, stop_lost, start_lost")
@click.option("--db", type=click.Path(), default=None, help="Database path")
@click.option("--no-jobs", is_flag=True,
              help="Don't queue annotation jobs after inserting variants")
@click.option("--limit", type=int, default=None,
              help="Stop after N enumerated variants (for smoke testing)")
def enumerate_cmd(gene: str, types: str, db: str, no_jobs: bool, limit: int):
    """Pre-populate every possible missense/nonsense SNV for a gene's canonical CDS.

    Pulls the canonical (preferentially MANE Select) transcript from Ensembl REST,
    walks every codon, and writes one row per allele to `variants` + a per-transcript
    row to `variant_consequences`. Annotation jobs are queued for downstream sources
    (AlphaMissense, REVEL, CADD, ClinVar, gnomAD, ClinGen AR) unless --no-jobs.
    """
    from .transcripts import fetch_canonical_transcript, TranscriptError
    from .enumerate import populate_for_transcript

    try:
        click.echo(f"Fetching canonical transcript for {gene}...")
        transcript = fetch_canonical_transcript(gene)
    except TranscriptError as e:
        click.echo(f"Transcript lookup failed: {e}", err=True)
        sys.exit(1)
    except requests.RequestException as e:
        click.echo(f"Network error talking to Ensembl REST: {e}", err=True)
        sys.exit(2)

    label = transcript.refseq_match or transcript.transcript_id_versioned
    click.echo(
        f"  {gene} -> {label}  chr{transcript.chromosome}  strand={transcript.strand:+d}  "
        f"CDS={transcript.cds_length}bp  MANE_Select={transcript.is_mane_select}"
    )

    vdb = VariantDB(Path(db) if db else None)
    type_list = [t.strip() for t in types.split(",") if t.strip()]

    summary = populate_for_transcript(
        transcript, vdb, types=type_list, enqueue=not no_jobs, limit=limit
    )

    click.echo()
    click.echo(f"Variants inserted/updated: {summary['variants']}")
    click.echo(f"Consequences written:      {summary['consequences']}")
    if not no_jobs:
        click.echo(f"Annotation jobs queued:    {summary['jobs_queued']}")
    click.echo()
    click.echo("By consequence:")
    for cons, n in sorted(summary["by_consequence"].items(), key=lambda x: -x[1]):
        click.echo(f"  {cons:<24} {n}")
    click.echo()
    click.echo(f"Database: {vdb.db_path}")


@main.command(name="frameshift-map")
@click.option("--gene", "-g", required=True, help="Gene symbol")
@click.option(
    "--aa-position",
    "aa_positions",
    type=int,
    multiple=True,
    help="Observed frameshift amino-acid position. Repeatable. Omit to map all positions.",
)
@click.option(
    "--n-steps",
    type=int,
    default=20,
    show_default=True,
    help="Maximum amino-acid distance to search left/right for a stop-gained SNV proxy.",
)
@click.option("--db", type=click.Path(), default=None, help="Database path")
@click.option(
    "--queue-source",
    "queue_sources",
    multiple=True,
    help="Annotation source to queue for proxy stop-gained variants, e.g. cadd. Repeatable.",
)
@click.option(
    "--no-persist",
    is_flag=True,
    help="Compute mappings without writing proxy variants or mapping rows to SQLite.",
)
@click.option("--format", "fmt", type=click.Choice(["table", "csv", "json"]), default="table")
def frameshift_map(
    gene: str,
    aa_positions: tuple[int, ...],
    n_steps: int,
    db: str,
    queue_sources: tuple[str, ...],
    no_persist: bool,
    fmt: str,
):
    """Map frameshift positions to concrete stop-gained SNV proxy variants.

    The proxy variants are normal `variants` rows, so CADD/MyVariant/gnomAD and
    other predictor annotations stay queryable through the existing schema.
    """
    from . import frameshift
    from .transcripts import TranscriptError

    gene = gene.upper()
    try:
        transcript = frameshift.fetch_canonical_transcript(gene)
    except TranscriptError as e:
        click.echo(f"Transcript lookup failed: {e}", err=True)
        sys.exit(1)
    except requests.RequestException as e:
        click.echo(f"Network error talking to Ensembl REST: {e}", err=True)
        sys.exit(2)

    try:
        if no_persist:
            rows = _frameshift_preview_rows(transcript, aa_positions, n_steps)
        else:
            vdb = VariantDB(Path(db) if db else None)
            if aa_positions:
                frameshift.annotate_positions(
                    vdb,
                    transcript,
                    aa_positions,
                    max_steps=n_steps,
                    enqueue_sources=queue_sources,
                    gene_symbol=gene,
                )
                rows = []
                for aa_pos in aa_positions:
                    rows.extend(
                        vdb.get_frameshift_nonsense_mappings(
                            gene,
                            frameshift_aa_pos=aa_pos,
                            transcript_id=transcript.transcript_id_versioned,
                        )
                    )
            else:
                frameshift.annotate_gene(
                    vdb,
                    gene,
                    transcript=transcript,
                    max_steps=n_steps,
                    enqueue_sources=queue_sources,
                )
                rows = vdb.get_frameshift_nonsense_mappings(
                    gene,
                    transcript_id=transcript.transcript_id_versioned,
                )
            _attach_proxy_pathogenicity(vdb, rows)
    except ValueError as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)

    output_rows = [_frameshift_output_row(row, gene, transcript.transcript_id_versioned) for row in rows]
    _print_frameshift_rows(output_rows, fmt)


def _frameshift_preview_rows(transcript, aa_positions: tuple[int, ...], n_steps: int) -> list[dict]:
    from . import frameshift

    positions = aa_positions or tuple(range(1, frameshift.protein_codon_count(transcript) + 1))
    rows = []
    for aa_pos in positions:
        for proxy in frameshift.map_frameshift_position(transcript, aa_pos, max_steps=n_steps):
            rows.append(proxy.to_record())
    return rows


def _attach_proxy_pathogenicity(vdb: VariantDB, rows: list[dict]) -> None:
    if not rows:
        return
    ids = sorted({row.get("proxy_variant_id") for row in rows if row.get("proxy_variant_id")})
    if not ids:
        return
    placeholders = ", ".join("?" for _ in ids)
    cur = vdb.conn.execute(
        f"""
        SELECT variant_id, predictor, predictor_version, score, category
        FROM annotations_pathogenicity
        WHERE variant_id IN ({placeholders})
        ORDER BY predictor, predictor_version
        """,
        ids,
    )
    by_variant: dict[int, list[dict]] = {}
    for row in cur.fetchall():
        by_variant.setdefault(row["variant_id"], []).append(dict(row))

    for row in rows:
        annotations = by_variant.get(row.get("proxy_variant_id"), [])
        row["cadd_phred"] = _first_predictor_score(annotations, "cadd_phred")
        row["cadd_raw"] = _first_predictor_score(annotations, "cadd_raw")
        row["pathogenicity"] = "|".join(_format_pathogenicity_annotation(a) for a in annotations)


def _first_predictor_score(annotations: list[dict], predictor: str):
    for annotation in annotations:
        if annotation.get("predictor") == predictor:
            return annotation.get("score")
    return None


def _format_pathogenicity_annotation(annotation: dict) -> str:
    name = annotation.get("predictor") or ""
    version = annotation.get("predictor_version") or ""
    if version:
        name = f"{name}@{version}"
    value = annotation.get("score")
    if annotation.get("category"):
        return f"{name}={value}({annotation['category']})"
    return f"{name}={value}"


def _frameshift_output_row(row: dict, gene: str, transcript_id: str) -> dict:
    proxy_variant = ""
    if row.get("proxy_chromosome") and row.get("proxy_position"):
        proxy_variant = (
            f"chr{row['proxy_chromosome']}:{row['proxy_position']}:"
            f"{row.get('proxy_ref')}>{row.get('proxy_alt')}"
        )
    return {
        "gene": row.get("gene_symbol") or gene,
        "transcript_id": row.get("transcript_id") or transcript_id,
        "frameshift_aa_pos": row.get("frameshift_aa_pos"),
        "direction": row.get("direction"),
        "n_steps_wrt_frameshift": row.get("n_steps_wrt_frameshift"),
        "proxy_variant_id": row.get("proxy_variant_id"),
        "proxy_variant": proxy_variant,
        "proxy_hgvs_c": row.get("proxy_hgvs_c"),
        "proxy_hgvs_p": row.get("proxy_hgvs_p"),
        "proxy_aa_pos": row.get("proxy_aa_pos"),
        "proxy_aa_ref": row.get("proxy_aa_ref"),
        "proxy_codon_ref": row.get("proxy_codon_ref"),
        "proxy_codon_alt": row.get("proxy_codon_alt"),
        "cadd_phred": row.get("cadd_phred"),
        "cadd_raw": row.get("cadd_raw"),
        "pathogenicity": row.get("pathogenicity"),
    }


def _print_frameshift_rows(rows: list[dict], fmt: str) -> None:
    fieldnames = [
        "gene",
        "transcript_id",
        "frameshift_aa_pos",
        "direction",
        "n_steps_wrt_frameshift",
        "proxy_variant_id",
        "proxy_variant",
        "proxy_hgvs_c",
        "proxy_hgvs_p",
        "proxy_aa_pos",
        "proxy_aa_ref",
        "proxy_codon_ref",
        "proxy_codon_alt",
        "cadd_phred",
        "cadd_raw",
        "pathogenicity",
    ]
    if fmt == "json":
        click.echo(json.dumps(rows, indent=2, default=str))
        return
    if fmt == "csv":
        writer = csv.DictWriter(sys.stdout, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
        return
    if not rows:
        click.echo("No frameshift proxy mappings found.")
        return
    click.echo(
        f"{'FS AA':>6} {'DIR':<5} {'D':>3} {'PROXY':<24} "
        f"{'HGVS.p':<28} {'CADD':>7}"
    )
    click.echo("-" * 82)
    for row in rows[:80]:
        click.echo(
            f"{str(row.get('frameshift_aa_pos')):>6} "
            f"{str(row.get('direction'))[:5]:<5} "
            f"{str(row.get('n_steps_wrt_frameshift')):>3} "
            f"{str(row.get('proxy_variant'))[:23]:<24} "
            f"{str(row.get('proxy_hgvs_p') or '')[:27]:<28} "
            f"{_format_score(row.get('cadd_phred')):>7}"
        )
    if len(rows) > 80:
        click.echo(f"... and {len(rows) - 80} more (use --format csv/json for full output)")


@main.command(name="annotate-pending")
@click.option("--source", default=None,
              help="Restrict to one source (e.g. clingen_ar). Default: all registered sources.")
@click.option("--gene", "-g", default=None, help="Restrict to pending jobs for one enumerated gene")
@click.option("--db", type=click.Path(), default=None, help="Database path")
@click.option("--limit", type=int, default=None,
              help="Max jobs to process this run (default: drain everything pending)")
@click.option("--batch-size", type=int, default=100, help="Jobs claimed per round trip")
@click.option("--rate-limit", type=float, default=None,
              help="Seconds to sleep between handler calls. Overrides handler default.")
def annotate_pending(source: str, gene: str, db: str, limit: int, batch_size: int, rate_limit: float):
    """Process pending annotation jobs from the queue.

    Each job calls the handler registered for its source, marks the job done on
    success, or failed (with error) on any exception. Designed to be safe to run
    repeatedly (e.g. on a launchd/cron schedule) — already-done jobs are not retouched.
    """
    from .worker import run_pending
    from .handlers import HANDLERS

    if source and source not in HANDLERS:
        click.echo(f"Error: no handler registered for source {source!r}.", err=True)
        click.echo(f"Registered: {', '.join(sorted(HANDLERS))}", err=True)
        sys.exit(1)

    vdb = VariantDB(Path(db) if db else None)
    gene_filter = gene.upper() if gene else None

    # Show what we're about to do
    total_pending = _pending_job_count(vdb, source=source, gene_filter=gene_filter)
    if total_pending == 0:
        where = []
        if source:
            where.append(f"source {source!r}")
        if gene_filter:
            where.append(f"gene {gene_filter}")
        suffix = f" for {' and '.join(where)}" if where else ""
        click.echo(f"No pending jobs{suffix}.")
        return

    where = []
    if source:
        where.append(source)
    if gene_filter:
        where.append(gene_filter)
    click.echo(f"{total_pending} pending job(s)" + (f" for {', '.join(where)}" if where else "") + ".")
    if limit is not None:
        click.echo(f"Processing up to {limit} this run.")

    last = {"n": 0}

    def on_event(ev):
        last["n"] += 1
        if last["n"] % 25 == 0:
            click.echo(f"  processed {last['n']}...", err=True)

    summary = run_pending(
        vdb,
        source=source,
        gene_filter=gene_filter,
        batch_size=batch_size,
        max_jobs=limit,
        rate_limit_sec=rate_limit,
        progress_callback=on_event,
    )

    click.echo()
    click.echo(f"Claimed: {summary['claimed']}   done: {summary['done']}   failed: {summary['failed']}   skipped: {summary['skipped']}")
    if summary["by_source"]:
        for src, counts in sorted(summary["by_source"].items()):
            click.echo(f"  {src:<20} done={counts['done']:<6} failed={counts['failed']:<6} skipped={counts['skipped']}")


def _pending_job_count(vdb: VariantDB, *, source: str | None, gene_filter: str | None) -> int:
    params: list = []
    join = ""
    where = ["j.status = 'pending'"]
    if gene_filter:
        join = """
            JOIN variant_consequences c
              ON c.variant_id = j.variant_id
             AND c.source = 'enumerated'
             AND UPPER(c.gene_symbol) = UPPER(?)
        """
        params.append(gene_filter)
    if source:
        where.append("j.source = ?")
        params.append(source)
    row = vdb.conn.execute(
        f"""
        SELECT COUNT(DISTINCT j.id)
        FROM annotation_jobs j
        {join}
        WHERE {' AND '.join(where)}
        """,
        params,
    ).fetchone()
    return int(row[0] or 0)


@main.command(name="myvariant-batch-run")
@click.option("--db", type=click.Path(), default=None, help="Database path")
@click.option("--gene", "-g", default=None, help="Restrict to pending jobs for one enumerated gene")
@click.option("--limit", type=int, default=None, help="Max pending jobs to claim this run")
@click.option("--batch-size", type=int, default=500, show_default=True)
@click.option("--source-label", default=None,
              help="Stored annotation source label. Default: myvariant_batch_hg38_<today>")
def myvariant_batch_run(db: str, gene: str, limit: int, batch_size: int, source_label: str):
    """Drain pending MyVariant jobs using MyVariant.info's batch endpoint."""
    from .handlers import myvariant

    vdb = VariantDB(Path(db) if db else None)
    try:
        result = myvariant.run_batch(
            vdb,
            limit=limit,
            gene_filter=gene,
            batch_size=batch_size,
            source_label=source_label,
        )
    except myvariant.HandlerError as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)
    click.echo(
        f"Claimed: {result['claimed']}   done: {result['done']}   "
        f"failed: {result['failed']}   not found: {result['notfound']}   "
        f"batches: {result['batches']}"
    )


@main.command()
@click.option("--source", required=True, help="Source name to queue (must have a registered handler)")
@click.option("--gene", default=None, help="Filter to one gene symbol (otherwise all variants in DB)")
@click.option("--db", type=click.Path(), default=None, help="Database path")
def queue(source: str, gene: str, db: str):
    """Enqueue a pending annotation job per variant for one source.

    Useful when you've added a new annotator after enumerating variants. Idempotent —
    rows that already exist for (variant, source) are left alone.
    """
    from .handlers import HANDLERS, BATCH_HANDLERS

    if source not in HANDLERS and source not in BATCH_HANDLERS:
        click.echo(f"Warning: no handler registered for source {source!r}.", err=True)
        click.echo(f"Registered per-job handlers: {', '.join(sorted(HANDLERS))}", err=True)
        click.echo(f"Registered batch handlers:   {', '.join(sorted(BATCH_HANDLERS))}", err=True)
        click.echo("Queueing anyway (you may have a custom handler).", err=True)

    vdb = VariantDB(Path(db) if db else None)
    n = vdb.enqueue_source_for_all_variants(source, gene_filter=gene)
    where = f" for gene {gene}" if gene else ""
    click.echo(f"Queued {n} new pending job(s) for source {source!r}{where}.")


@main.command(name="alphamissense-run")
@click.option("--db", type=click.Path(), default=None, help="Database path")
@click.option("--file", "file_path", type=click.Path(), default=None,
              help="AlphaMissense TSV(.gz) path (env: ALPHAMISSENSE_FILE)")
@click.option("--gene", "-g", default=None, help="Restrict to pending jobs for one enumerated gene")
@click.option("--limit", type=int, default=None, help="Max pending jobs to claim this run")
def alphamissense_run(db: str, file_path: str, gene: str, limit: int):
    """Drain pending alphamissense jobs by single-pass over the AlphaMissense TSV."""
    from .handlers import alphamissense as am

    vdb = VariantDB(Path(db) if db else None)
    try:
        result = am.run_batch(vdb, file_path=file_path, limit=limit, gene_filter=gene)
    except am.HandlerError as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)
    click.echo(f"Claimed {result['claimed']}   matched: {result['matched']}   "
               f"failed: {result['failed']}   lines scanned: {result['lines_scanned']}")


@main.command(name="revel-run")
@click.option("--db", type=click.Path(), default=None, help="Database path")
@click.option("--file", "file_path", type=click.Path(), default=None,
              help="REVEL revel_with_transcript_ids file or zip path (env: REVEL_FILE)")
@click.option("--gene", "-g", default=None, help="Restrict to pending jobs for one enumerated gene")
@click.option("--limit", type=int, default=None, help="Max pending jobs to claim this run")
def revel_run(db: str, file_path: str, gene: str, limit: int):
    """Drain pending REVEL jobs by single-pass over the local REVEL file."""
    from .handlers import revel

    vdb = VariantDB(Path(db) if db else None)
    try:
        result = revel.run_batch(vdb, file_path=file_path, limit=limit, gene_filter=gene)
    except revel.HandlerError as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)
    click.echo(
        f"Claimed {result['claimed']}   matched: {result['matched']}   "
        f"failed: {result['failed']}   lines scanned: {result['lines_scanned']}"
    )


@main.command(name="annovar-run")
@click.option("--db", type=click.Path(), default=None, help="Database path")
@click.option("--gene", "-g", default=None, help="Restrict to pending jobs for one enumerated gene")
@click.option("--limit", type=int, default=None, help="Max pending jobs to claim this run")
@click.option("--build", default="hg38", show_default=True)
@click.option("--protocols", default=None,
              help="ANNOVAR -protocol value. Default uses refGeneWithVer,gnomad,clinvar,dbnsfp.")
@click.option("--operations", default=None,
              help="ANNOVAR -operation value matching --protocols.")
@click.option("--keep-outputs", is_flag=True, help="Don't delete the temp work directory.")
def annovar_run(db: str, gene: str, limit: int, build: str, protocols: str, operations: str, keep_outputs: bool):
    """Drain pending annovar jobs by invoking ANNOVAR's table_annovar.pl once.

    Requires ANNOVAR installed. Set ANNOVAR_HOME (and ANNOVAR_DB if humandb is elsewhere).
    """
    from .handlers import annovar

    vdb = VariantDB(Path(db) if db else None)
    kwargs = {"limit": limit, "gene_filter": gene, "build": build, "keep_outputs": keep_outputs}
    if protocols:
        kwargs["protocols"] = protocols
    if operations:
        kwargs["operations"] = operations
    try:
        result = annovar.run_batch(vdb, **kwargs)
    except annovar.HandlerError as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)
    click.echo(f"Claimed {result['claimed']}   matched: {result['matched']}   "
               f"missing: {result['missing']}   parsed lines: {result['lines_parsed']}")


@main.command(name="annovar-import")
@click.argument("multianno_path", type=click.Path(exists=True, dir_okay=False))
@click.option("--db", type=click.Path(), default=None, help="Database path")
def annovar_import(multianno_path: str, db: str):
    """Import an externally-produced *_multianno.txt without running ANNOVAR."""
    from .handlers import annovar

    vdb = VariantDB(Path(db) if db else None)
    counts = annovar.import_multianno(vdb, multianno_path)
    click.echo(f"Imported {counts['variants']} variants ({counts['annotated']} with annotation rows).")


@main.command(name="dbscsnv-import")
@click.argument("dbscsnv_path", type=click.Path(exists=True, dir_okay=False))
@click.option("--db", type=click.Path(), default=None, help="Database path")
@click.option("--gene", "-g", default=None, help="Restrict matching to one enumerated gene")
@click.option("--source-label", default="annovar_hg38_dbscsnv11", show_default=True)
def dbscsnv_import(dbscsnv_path: str, db: str, gene: str, source_label: str):
    """Import dbscSNV scores directly from an ANNOVAR humandb text file."""
    from .handlers import annovar

    vdb = VariantDB(Path(db) if db else None)
    counts = annovar.import_dbscsnv(
        vdb,
        dbscsnv_path,
        gene_filter=gene,
        source_label=source_label,
    )
    click.echo(
        f"Indexed {counts['indexed_variants']} variants   scanned rows: {counts['rows']}   "
        f"matched rows: {counts['matched_rows']}   splice annotations: {counts['splice']}"
    )


@main.command(name="clinvar-import")
@click.argument("variant_summary_path", required=False, type=click.Path(dir_okay=False))
@click.option("--db", type=click.Path(), default=None, help="Database path")
@click.option("--gene", "-g", required=True, help="Restrict matching to one enumerated gene")
@click.option("--assembly", default="GRCh38", show_default=True)
@click.option("--include-all-types", is_flag=True, help="Import non-SNV ClinVar rows too")
def clinvar_import(
    variant_summary_path: str,
    db: str,
    gene: str,
    assembly: str,
    include_all_types: bool,
):
    """Import ClinVar variant_summary rows into normalized clinical annotations."""
    from .handlers import clinvar

    vdb = VariantDB(Path(db) if db else None)
    try:
        counts = clinvar.import_gene(
            vdb,
            gene,
            data_file=variant_summary_path,
            assembly=assembly,
            include_all_types=include_all_types,
        )
    except clinvar.HandlerError as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)
    click.echo(
        f"Scanned ClinVar rows: {counts['rows']}   matched rows: {counts['matched_rows']}   "
        f"clinical annotations: {counts['annotations']}   aliases: {counts['aliases']}"
    )


@main.command(name="vep-run")
@click.option("--db", type=click.Path(), default=None, help="Database path")
@click.option("--gene", "-g", default=None, help="Restrict to pending jobs for one enumerated gene")
@click.option("--limit", type=int, default=None, help="Max pending jobs to claim this run")
@click.option("--plugin", "plugins", multiple=True,
              help="VEP --plugin argument; repeatable. Overrides VEP_PLUGINS env.")
@click.option("--keep-outputs", is_flag=True)
def vep_run(db: str, gene: str, limit: int, plugins: tuple, keep_outputs: bool):
    """Drain pending vep jobs by invoking VEP once on a generated VCF.

    Requires VEP installed (`brew install ensembl-vep`) and a cache (`~/.vep` by default).
    """
    from .handlers import vep

    vdb = VariantDB(Path(db) if db else None)
    try:
        result = vep.run_batch(
            vdb,
            limit=limit,
            gene_filter=gene,
            plugins=list(plugins) or None,
            keep_outputs=keep_outputs,
        )
    except vep.HandlerError as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)
    click.echo(f"Claimed {result['claimed']}   annotated: {result['annotated']}   "
               f"missing: {result['missing']}")


@main.command(name="alphafold-run")
@click.option("--db", type=click.Path(), default=None, help="Database path")
@click.option("--gene", "-g", default=None, help="Restrict to pending jobs for one enumerated gene")
@click.option("--limit", type=int, default=None, help="Max pending jobs to claim this run")
def alphafold_run(db: str, gene: str, limit: int):
    """Drain pending alphafold jobs via AlphaFold DB pLDDT JSON files."""
    from .handlers import alphafold

    vdb = VariantDB(Path(db) if db else None)
    result = alphafold.run_batch(vdb, limit=limit, gene_filter=gene)
    click.echo(
        f"Claimed {result['claimed']}   annotated: {result['annotated']}   "
        f"failed: {result['failed']}"
    )


@main.command(name="pext-import")
@click.argument("path", type=click.Path(exists=True, dir_okay=False))
@click.option("--gene", default=None, help="Filter to one gene symbol")
@click.option("--db", type=click.Path(), default=None, help="Database path")
@click.option("--dataset", default="gnomad_pext_v10", show_default=True)
def pext_import(path: str, gene: str, db: str, dataset: str):
    """Import local gnomAD pext CSV/TSV rows into annotations_expression."""
    from .handlers import pext

    vdb = VariantDB(Path(db) if db else None)
    try:
        result = pext.import_file(vdb, path, gene_filter=gene, dataset=dataset)
    except pext.HandlerError as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)
    click.echo(
        f"Rows: {result['rows']}   matched rows: {result['matched_rows']}   "
        f"annotations: {result['annotations']}"
    )


@main.command(name="pext-run")
@click.option("--file", "file_path", type=click.Path(exists=True, dir_okay=False), required=True,
              help="Local pext CSV/TSV(.gz) file")
@click.option("--gene", default=None, help="Filter to one gene symbol")
@click.option("--db", type=click.Path(), default=None, help="Database path")
@click.option("--limit", type=int, default=None, help="Max pending jobs to claim this run")
@click.option("--dataset", default="gnomad_pext_v10", show_default=True)
def pext_run(file_path: str, gene: str, db: str, limit: int, dataset: str):
    """Drain queued pext jobs from a local gnomAD pext export."""
    from .handlers import pext

    vdb = VariantDB(Path(db) if db else None)
    result = pext.run_batch(vdb, file_path=file_path, limit=limit, gene_filter=gene, dataset=dataset)
    click.echo(
        f"Claimed {result['claimed']}   annotated variants: {result['annotated']}   "
        f"failed: {result['failed']}   rows scanned: {result['rows']}"
    )


@main.command(name="pext-bigwig-run")
@click.option("--dir", "dir_path", type=click.Path(exists=True, file_okay=False),
              default="data/pext/ucsc_hg38", show_default=True,
              help="Directory containing one pext .bw file per tissue")
@click.option("--gene", default=None, help="Filter to one gene symbol")
@click.option("--db", type=click.Path(), default=None, help="Database path")
@click.option("--dataset", default="ucsc_gnomad_pext_hg38", show_default=True)
@click.option("--source-label", default="ucsc_hg38_gnomad_pext_bigwig", show_default=True)
@click.option("--bigwig-average-bin", default=None,
              help="Path to bigWigAverageOverBed (env: BIGWIG_AVERAGE_OVER_BED)")
def pext_bigwig_run(
    dir_path: str,
    gene: str,
    db: str,
    dataset: str,
    source_label: str,
    bigwig_average_bin: str,
):
    """Import UCSC/gnomAD pext bigWig tracks at variant positions."""
    from .handlers import pext

    vdb = VariantDB(Path(db) if db else None)
    try:
        result = pext.import_bigwig_dir(
            vdb,
            dir_path,
            gene_filter=gene,
            dataset=dataset,
            source_label=source_label,
            bigwig_average_bin=bigwig_average_bin,
        )
    except pext.HandlerError as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)
    click.echo(
        f"Variants: {result['variants']}   tissues: {result['tissues']}   "
        f"annotations: {result['annotations']}"
    )


@main.command(name="absplice-import")
@click.argument("path", type=click.Path(exists=True, dir_okay=False))
@click.option("--gene", default=None, help="Filter to one gene symbol")
@click.option("--db", type=click.Path(), default=None, help="Database path")
@click.option("--dataset", default="absplice", show_default=True)
@click.option("--source-label", default="absplice", show_default=True)
def absplice_import(path: str, gene: str, db: str, dataset: str, source_label: str):
    """Import AbSplice/AbExp tabular output by genomic coordinate."""
    from .handlers import absplice

    vdb = VariantDB(Path(db) if db else None)
    try:
        result = absplice.import_file(
            vdb,
            path,
            gene_filter=gene,
            dataset=dataset,
            source_label=source_label,
        )
    except absplice.HandlerError as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)
    click.echo(
        f"Rows: {result['rows']}   matched rows: {result['matched_rows']}   "
        f"splice annotations: {result['splice']}   expression annotations: {result['expression']}"
    )


@main.command(name="nmd-import")
@click.argument("path", type=click.Path(exists=True, dir_okay=False))
@click.option("--predictor", required=True, help="Predictor name, e.g. nmdep or nmdetective")
@click.option("--gene", default=None, help="Filter to one gene symbol")
@click.option("--score-column", default=None, help="Explicit score column")
@click.option("--category-column", default=None, help="Optional category/class column")
@click.option("--db", type=click.Path(), default=None, help="Database path")
def nmd_import(path: str, predictor: str, gene: str, score_column: str, category_column: str, db: str):
    """Import coordinate-keyed NMDEP/NMDetective-style NMD predictor output."""
    from .handlers import nmd_external

    vdb = VariantDB(Path(db) if db else None)
    try:
        result = nmd_external.import_file(
            vdb,
            path,
            predictor=predictor,
            gene_filter=gene,
            score_column=score_column,
            category_column=category_column,
        )
    except nmd_external.HandlerError as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)
    click.echo(
        f"Rows: {result['rows']}   matched rows: {result['matched_rows']}   "
        f"annotations: {result['annotations']}"
    )


@main.command(name="vep-plugin-config")
@click.option("--spliceai-snv", type=click.Path(), default=None, help="SpliceAI SNV VCF.gz")
@click.option("--spliceai-indel", type=click.Path(), default=None, help="SpliceAI indel VCF.gz")
@click.option("--dbscsnv", type=click.Path(), default=None, help="dbscSNV tabix-indexed TXT.gz")
@click.option("--maxentscan-dir", type=click.Path(), default=None, help="Unpacked MaxEntScan directory")
@click.option("--loftee-dir", type=click.Path(), default=None, help="LOFTEE plugin directory")
@click.option("--loftee-data-dir", type=click.Path(), default=None, help="LOFTEE data directory")
@click.option("--include-nmd", is_flag=True, help="Include VEP NMD plugin")
def vep_plugin_config(
    spliceai_snv: str,
    spliceai_indel: str,
    dbscsnv: str,
    maxentscan_dir: str,
    loftee_dir: str,
    loftee_data_dir: str,
    include_nmd: bool,
):
    """Print VEP_PLUGINS value for the deferred splice/NMD plugins."""
    plugins = []
    if spliceai_snv and spliceai_indel:
        plugins.append(f"SpliceAI,snv={spliceai_snv},indel={spliceai_indel},split_output=1")
    if dbscsnv:
        plugins.append(f"dbscSNV,{dbscsnv},GRCh38")
    if maxentscan_dir:
        plugins.append(f"MaxEntScan,{maxentscan_dir},SWA,NCSS")
    if loftee_dir:
        data_dir = Path(loftee_data_dir) if loftee_data_dir else Path("data") / "loftee_data"
        plugins.append(
            "LoF,"
            f"loftee_path:{loftee_dir},"
            f"human_ancestor_fa:{data_dir / 'human_ancestor.fa.gz'},"
            f"gerp_bigwig:{data_dir / 'gerp_conservation_scores.homo_sapiens.GRCh38.bw'},"
            f"conservation_file:{data_dir / 'loftee.sql'}"
        )
    if include_nmd:
        plugins.append("NMD")
    if not plugins:
        click.echo("No plugin paths provided.", err=True)
        return
    click.echo(";".join(plugins))


@main.command(name="nmd-rule")
@click.option("--gene", "-g", required=True, multiple=True,
              help="Gene symbol; pass multiple for multiple genes")
@click.option("--db", type=click.Path(), default=None, help="Database path")
def nmd_rule_cmd(gene: tuple, db: str):
    """Apply the last-exon + 50-nt-rule NMD trigger/escape classifier.

    Writes one annotations_pathogenicity row per stop_gained variant with
    predictor='nmd_rule', score 0.0 (escape) or 1.0 (triggers), and a category
    string for which rule fired.
    """
    from .handlers import nmd_rules

    vdb = VariantDB(Path(db) if db else None)
    for sym in gene:
        sym = sym.upper()
        try:
            result = nmd_rules.annotate_gene(vdb, sym)
        except nmd_rules.HandlerError as e:
            click.echo(f"{sym}: error - {e}", err=True)
            continue
        click.echo(f"{sym}: {result['considered']} stop_gained variants")
        click.echo(f"  triggers NMD: {result['triggers']}    escapes: {result['escapes']}")
        for cat, n in sorted(result["by_category"].items(), key=lambda x: -x[1]):
            click.echo(f"    {cat:<32} {n}")


@main.command(name="gene-constraint")
@click.option("--gene", "-g", required=True, multiple=True,
              help="Gene symbol; pass multiple times for multiple genes")
@click.option("--db", type=click.Path(), default=None, help="Database path")
def gene_constraint_cmd(gene: tuple, db: str):
    """Pull gnomAD pLI / LOEUF / mis_z / o-e LoF for one or more genes.

    Most predictive single feature for novel LoF variant interpretation.
    Stored in gene_constraint, joinable via variant_consequences.gene_symbol.
    """
    from .handlers import gnomad_constraint

    vdb = VariantDB(Path(db) if db else None)
    for sym in gene:
        sym = sym.upper()
        try:
            ok = gnomad_constraint.annotate_gene(vdb, sym)
            if ok:
                cur = vdb.conn.execute(
                    "SELECT pli, oe_lof_upper, mis_z, obs_lof, exp_lof FROM gene_constraint WHERE gene_symbol = ?",
                    [sym],
                )
                row = cur.fetchone()
                if row:
                    click.echo(
                        f"{sym}: pLI={row['pli']:.3f}  LOEUF={row['oe_lof_upper']:.3f}  "
                        f"mis_z={row['mis_z']:.2f}  obs_lof={row['obs_lof']}/{row['exp_lof']:.0f}"
                    )
            else:
                click.echo(f"{sym}: no gnomAD constraint data available")
        except gnomad_constraint.HandlerError as e:
            click.echo(f"{sym}: error - {e}", err=True)


@main.command()
@click.option("--db", type=click.Path(), default=None, help="Database path")
def jobs(db: str):
    """Show annotation job queue status by source."""
    vdb = VariantDB(Path(db) if db else None)
    rows = vdb.job_status_counts()
    if not rows:
        click.echo("(queue is empty)")
        return
    click.echo(f"{'STATUS':<12} {'SOURCE':<20} {'N':>10}")
    click.echo("-" * 44)
    for r in rows:
        click.echo(f"{r['status']:<12} {r['source']:<20} {r['n']:>10}")


@main.command(name="feature-schema")
@click.option("--format", "fmt", type=click.Choice(["table", "json"]), default="table", show_default=True)
def feature_schema_cmd(fmt: str):
    """Show where normalized feature families are stored and exported."""
    from .normalized_export import feature_schema, feature_schema_json

    if fmt == "json":
        click.echo(feature_schema_json())
        return

    rows = feature_schema()
    click.echo(f"{'GROUP':<18} {'TABLE':<28} EXPORT COLUMN PATTERN")
    click.echo("-" * 100)
    for group, info in rows.items():
        click.echo(f"{group:<18} {info['table']:<28} {info['wide_prefix']}")


@main.command(name="feature-coverage")
@click.option("--gene", "-g", required=True, help="Gene symbol")
@click.option("--db", type=click.Path(), default=None, help="Database path")
@click.option(
    "--kind",
    type=click.Choice([
        "all", "pathogenicity", "population", "clinical", "splice",
        "structure", "expression", "conservation",
    ]),
    default="all",
    show_default=True,
)
@click.option("--limit", type=int, default=40, show_default=True)
def feature_coverage(gene: str, db: str, kind: str, limit: int):
    """Show normalized annotation coverage for a gene."""
    vdb = _open_db(db, read_only=True)
    gene = gene.upper()
    total = vdb.conn.execute(
        """
        SELECT COUNT(DISTINCT variant_id)
        FROM variant_consequences
        WHERE UPPER(gene_symbol) = ? AND source = 'enumerated'
        """,
        [gene],
    ).fetchone()[0]
    if not total:
        click.echo(f"No normalized enumerated variants found for {gene}")
        return

    click.echo(f"{gene}: {total} normalized enumerated variant(s)")

    def show(title: str, sql: str, params: list):
        click.echo()
        click.echo(title)
        click.echo("-" * len(title))
        rows = vdb.conn.execute(sql, params).fetchall()
        if not rows:
            click.echo("(none)")
            return
        click.echo(f"{'feature':<48} {'variants':>9} {'rows':>9} {'pct':>7}")
        for r in rows[:limit]:
            pct = 100.0 * (r["variants"] or 0) / total
            click.echo(f"{r['feature']:<48} {r['variants']:>9} {r['rows']:>9} {pct:>6.1f}%")
        if len(rows) > limit:
            click.echo(f"... {len(rows) - limit} more")

    gene_cte = """
        WITH gene_vars AS (
            SELECT DISTINCT variant_id
            FROM variant_consequences
            WHERE UPPER(gene_symbol) = ? AND source = 'enumerated'
        )
    """
    sections = {
        "pathogenicity": (
            "Pathogenicity",
            gene_cte + """
            SELECT p.predictor AS feature, COUNT(DISTINCT p.variant_id) AS variants, COUNT(*) AS rows
            FROM annotations_pathogenicity p JOIN gene_vars g ON g.variant_id = p.variant_id
            GROUP BY p.predictor
            ORDER BY variants DESC, feature
            """,
        ),
        "population": (
            "Population",
            gene_cte + """
            SELECT p.dataset AS feature, COUNT(DISTINCT p.variant_id) AS variants, COUNT(*) AS rows
            FROM annotations_population p JOIN gene_vars g ON g.variant_id = p.variant_id
            GROUP BY p.dataset
            ORDER BY variants DESC, feature
            """,
        ),
        "clinical": (
            "Clinical",
            gene_cte + """
            SELECT c.source AS feature, COUNT(DISTINCT c.variant_id) AS variants, COUNT(*) AS rows
            FROM annotations_clinical c JOIN gene_vars g ON g.variant_id = c.variant_id
            GROUP BY c.source
            ORDER BY variants DESC, feature
            """,
        ),
        "splice": (
            "Splice",
            gene_cte + """
            SELECT s.predictor || ' [' || COALESCE(s.source, '') || ']' AS feature,
                   COUNT(DISTINCT s.variant_id) AS variants, COUNT(*) AS rows
            FROM annotations_splice s JOIN gene_vars g ON g.variant_id = s.variant_id
            GROUP BY s.predictor, s.source
            ORDER BY variants DESC, feature
            """,
        ),
        "structure": (
            "Structure",
            gene_cte + """
            SELECT s.feature AS feature, COUNT(DISTINCT s.variant_id) AS variants, COUNT(*) AS rows
            FROM annotations_structure s JOIN gene_vars g ON g.variant_id = s.variant_id
            GROUP BY s.feature
            ORDER BY variants DESC, feature
            """,
        ),
        "expression": (
            "Expression",
            gene_cte + """
            SELECT e.metric || ' [' || e.dataset || ']' AS feature,
                   COUNT(DISTINCT e.variant_id) AS variants, COUNT(*) AS rows
            FROM annotations_expression e JOIN gene_vars g ON g.variant_id = e.variant_id
            GROUP BY e.metric, e.dataset
            ORDER BY variants DESC, feature
            """,
        ),
        "conservation": (
            "Conservation",
            gene_cte + """
            SELECT c.metric AS feature, COUNT(DISTINCT c.variant_id) AS variants, COUNT(*) AS rows
            FROM annotations_conservation c JOIN gene_vars g ON g.variant_id = c.variant_id
            GROUP BY c.metric
            ORDER BY variants DESC, feature
            """,
        ),
    }

    for key, (title, sql) in sections.items():
        if kind in ("all", key):
            show(title, sql, [gene])


if __name__ == "__main__":
    main()
