"""Command-line interface for VariantFeatures."""

import click
import csv
import json
import sys
from pathlib import Path

import requests

from .database import VariantDB


@click.group()
@click.version_option()
def main():
    """Aggregate predictive features for genetic variants."""
    pass


@main.command()
@click.option("--genes", "-g", required=True, help="Comma-separated gene symbols")
@click.option("--db", type=click.Path(), default=None, help="Database path")
@click.option("--sources", default="all", help="Sources: alphamissense,gnomad,clinvar (comma-separated or 'all')")
@click.option("--skip-download", is_flag=True, help="Skip downloading data files")
def build(genes: str, db: str, sources: str, skip_download: bool):
    """Fetch features and build/update database."""
    gene_list = [g.strip().upper() for g in genes.split(",")]
    db_path = Path(db) if db else None
    
    # Parse sources
    if sources == "all":
        source_list = ["alphamissense", "gnomad", "clinvar"]
    else:
        source_list = [s.strip().lower() for s in sources.split(",")]
    
    click.echo(f"Building database for genes: {', '.join(gene_list)}")
    click.echo(f"Sources: {', '.join(source_list)}")
    
    vdb = VariantDB(db_path)
    
    for gene in gene_list:
        click.echo(f"\n{'='*60}")
        click.echo(f"Processing {gene}")
        click.echo(f"{'='*60}")
        
        # AlphaMissense
        if "alphamissense" in source_list:
            click.echo(f"\n[AlphaMissense] Fetching scores for {gene}...")
            try:
                from .fetchers.alphamissense import fetch_alphamissense
                count = 0
                for variant in fetch_alphamissense(gene):
                    vdb.upsert_missense(
                        gene=gene,
                        hgvs_p=variant['hgvs_p'],
                        alphamissense_score=variant['alphamissense_score'],
                        alphamissense_class=variant['alphamissense_class'],
                    )
                    count += 1
                click.echo(f"  Loaded {count} AlphaMissense scores")
            except FileNotFoundError as e:
                click.echo(f"  Skipped: AlphaMissense data not downloaded yet")
                click.echo(f"  Run: python -m variantfeatures.fetchers.alphamissense")
            except ValueError as e:
                click.echo(f"  Skipped: {e}")
            except Exception as e:
                click.echo(f"  Error: {e}")
        
        # gnomAD
        if "gnomad" in source_list:
            click.echo(f"\n[gnomAD] Fetching frequencies for {gene}...")
            try:
                from .fetchers.gnomad import fetch_gnomad
                count = 0
                for variant in fetch_gnomad(gene):
                    if variant.get('hgvs_p'):
                        vdb.upsert_missense(
                            gene=gene,
                            hgvs_p=variant['hgvs_p'],
                            hgvs_c=variant.get('hgvs_c'),
                            gnomad_af=variant.get('gnomad_af'),
                            gnomad_homozygotes=variant.get('gnomad_homozygotes'),
                        )
                        count += 1
                click.echo(f"  Loaded {count} gnomAD variants")
            except Exception as e:
                click.echo(f"  Error: {e}")
        
        # ClinVar - uses existing data from load_clinvar.py
        if "clinvar" in source_list:
            click.echo(f"\n[ClinVar] Note: ClinVar data loaded separately via scripts/load_clinvar.py")
    
    click.echo(f"\n{'='*60}")
    click.echo("Build complete!")
    click.echo(f"Database: {vdb.db_path}")


@main.command()
@click.option("--gene", "-g", required=True, help="Gene symbol")
@click.option("--db", type=click.Path(), default=None, help="Database path")
@click.option("--format", "fmt", type=click.Choice(["csv", "json", "table"]), default="table")
@click.option("--include-lof", is_flag=True, help="Include loss-of-function variants")
def query(gene: str, db: str, fmt: str, include_lof: bool):
    """Query variants for a gene."""
    db_path = Path(db) if db else None
    vdb = VariantDB(db_path)
    
    # Get missense variants (primary)
    variants = vdb.get_gene_missense(gene.upper())
    
    # Optionally include LOF
    if include_lof:
        lof_variants = vdb.get_gene_lof(gene.upper())
        # Add type marker
        for v in variants:
            v['variant_type'] = 'missense'
        for v in lof_variants:
            v['variant_type'] = 'lof'
        variants.extend(lof_variants)
    
    if not variants:
        click.echo(f"No variants found for {gene}")
        return
    
    if fmt == "json":
        # JSON output
        click.echo(json.dumps(variants, indent=2, default=str))
    
    elif fmt == "csv":
        # CSV output
        if variants:
            fieldnames = list(variants[0].keys())
            writer = csv.DictWriter(sys.stdout, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(variants)
    
    else:
        # Table output (default)
        click.echo(f"Found {len(variants)} variants for {gene}")
        click.echo()
        click.echo(f"{'HGVS':<20} {'ClinVar':<20} {'AM Score':<10} {'gnomAD AF':<12} {'REVEL':<8}")
        click.echo("-" * 80)
        
        for v in variants[:50]:  # Limit output
            hgvs = v.get('hgvs_p', 'N/A')[:19]
            clinvar = (v.get('clinvar_significance') or 'N/A')[:19]
            am = v.get('alphamissense_score')
            am_str = f"{am:.3f}" if am else 'N/A'
            gnomad = v.get('gnomad_af')
            gnomad_str = f"{gnomad:.2e}" if gnomad else 'N/A'
            revel = v.get('revel_score')
            revel_str = f"{revel:.3f}" if revel else 'N/A'
            
            click.echo(f"{hgvs:<20} {clinvar:<20} {am_str:<10} {gnomad_str:<12} {revel_str:<8}")
        
        if len(variants) > 50:
            click.echo(f"\n... and {len(variants) - 50} more variants (use --format csv for full output)")


@main.command()
@click.option("--db", type=click.Path(), default=None, help="Database path")
def stats(db: str):
    """Show database statistics."""
    db_path = Path(db) if db else None
    vdb = VariantDB(db_path)
    
    click.echo("VariantFeatures Database Statistics")
    click.echo("=" * 60)
    
    # Get counts per gene
    cur = vdb.conn.execute("""
        SELECT 
            gene,
            COUNT(*) as total,
            SUM(CASE WHEN clinvar_id IS NOT NULL THEN 1 ELSE 0 END) as clinvar,
            SUM(CASE WHEN alphamissense_score IS NOT NULL THEN 1 ELSE 0 END) as alphamissense,
            SUM(CASE WHEN gnomad_af IS NOT NULL THEN 1 ELSE 0 END) as gnomad,
            SUM(CASE WHEN revel_score IS NOT NULL THEN 1 ELSE 0 END) as revel
        FROM variants_missense 
        GROUP BY gene
        ORDER BY total DESC
    """)
    
    click.echo(f"\n{'Gene':<10} {'Total':<8} {'ClinVar':<10} {'AlphaMissense':<15} {'gnomAD':<10} {'REVEL':<10}")
    click.echo("-" * 70)
    
    for row in cur:
        click.echo(f"{row[0]:<10} {row[1]:<8} {row[2]:<10} {row[3]:<15} {row[4]:<10} {row[5]:<10}")
    
    # LOF stats
    lof_cur = vdb.conn.execute("SELECT gene, COUNT(*) FROM variants_lof GROUP BY gene ORDER BY gene")
    lof_data = list(lof_cur)
    
    if lof_data:
        click.echo(f"\nLOF Variants:")
        click.echo("-" * 30)
        for row in lof_data:
            click.echo(f"  {row[0]}: {row[1]}")
    
    click.echo(f"\nDatabase: {vdb.db_path}")


@main.command()
@click.option("--gene", "-g", required=True, help="Gene symbol")
@click.option("--db", type=click.Path(), default=None, help="Database path")
@click.option("--output", "-o", type=click.Path(), required=True, help="Output CSV file")
def export(gene: str, db: str, output: str):
    """Export variants for downstream pipelines."""
    db_path = Path(db) if db else None
    vdb = VariantDB(db_path)
    
    variants = vdb.get_gene_missense(gene.upper())
    
    if not variants:
        click.echo(f"No variants found for {gene}")
        return
    
    output_path = Path(output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_path, 'w', newline='') as f:
        fieldnames = list(variants[0].keys())
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(variants)
    
    click.echo(f"Exported {len(variants)} variants to {output_path}")


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


@main.command(name="annotate-pending")
@click.option("--source", default=None,
              help="Restrict to one source (e.g. clingen_ar). Default: all registered sources.")
@click.option("--db", type=click.Path(), default=None, help="Database path")
@click.option("--limit", type=int, default=None,
              help="Max jobs to process this run (default: drain everything pending)")
@click.option("--batch-size", type=int, default=100, help="Jobs claimed per round trip")
@click.option("--rate-limit", type=float, default=None,
              help="Seconds to sleep between handler calls. Overrides handler default.")
def annotate_pending(source: str, db: str, limit: int, batch_size: int, rate_limit: float):
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

    # Show what we're about to do
    pending_rows = [r for r in vdb.job_status_counts() if r["status"] == "pending"]
    if source:
        pending_rows = [r for r in pending_rows if r["source"] == source]
    total_pending = sum(r["n"] for r in pending_rows)
    if total_pending == 0:
        click.echo(f"No pending jobs{f' for source {source!r}' if source else ''}.")
        return

    click.echo(f"{total_pending} pending job(s)" + (f" for {source}" if source else "") + ".")
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
@click.option("--limit", type=int, default=None, help="Max pending jobs to claim this run")
def alphamissense_run(db: str, file_path: str, limit: int):
    """Drain pending alphamissense jobs by single-pass over the AlphaMissense TSV."""
    from .handlers import alphamissense as am

    vdb = VariantDB(Path(db) if db else None)
    try:
        result = am.run_batch(vdb, file_path=file_path, limit=limit)
    except am.HandlerError as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)
    click.echo(f"Claimed {result['claimed']}   matched: {result['matched']}   "
               f"failed: {result['failed']}   lines scanned: {result['lines_scanned']}")


@main.command(name="annovar-run")
@click.option("--db", type=click.Path(), default=None, help="Database path")
@click.option("--limit", type=int, default=None, help="Max pending jobs to claim this run")
@click.option("--build", default="hg38", show_default=True)
@click.option("--protocols", default=None,
              help="ANNOVAR -protocol value. Default uses refGene,gnomad,clinvar,dbnsfp.")
@click.option("--operations", default=None,
              help="ANNOVAR -operation value matching --protocols.")
@click.option("--keep-outputs", is_flag=True, help="Don't delete the temp work directory.")
def annovar_run(db: str, limit: int, build: str, protocols: str, operations: str, keep_outputs: bool):
    """Drain pending annovar jobs by invoking ANNOVAR's table_annovar.pl once.

    Requires ANNOVAR installed. Set ANNOVAR_HOME (and ANNOVAR_DB if humandb is elsewhere).
    """
    from .handlers import annovar

    vdb = VariantDB(Path(db) if db else None)
    kwargs = {"limit": limit, "build": build, "keep_outputs": keep_outputs}
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


@main.command(name="vep-run")
@click.option("--db", type=click.Path(), default=None, help="Database path")
@click.option("--limit", type=int, default=None, help="Max pending jobs to claim this run")
@click.option("--plugin", "plugins", multiple=True,
              help="VEP --plugin argument; repeatable. Overrides VEP_PLUGINS env.")
@click.option("--keep-outputs", is_flag=True)
def vep_run(db: str, limit: int, plugins: tuple, keep_outputs: bool):
    """Drain pending vep jobs by invoking VEP once on a generated VCF.

    Requires VEP installed (`brew install ensembl-vep`) and a cache (`~/.vep` by default).
    """
    from .handlers import vep

    vdb = VariantDB(Path(db) if db else None)
    try:
        result = vep.run_batch(vdb, limit=limit, plugins=list(plugins) or None, keep_outputs=keep_outputs)
    except vep.HandlerError as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)
    click.echo(f"Claimed {result['claimed']}   annotated: {result['annotated']}   "
               f"missing: {result['missing']}")


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


if __name__ == "__main__":
    main()
