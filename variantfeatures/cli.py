"""Command-line interface for VariantFeatures."""

import click
import csv
import json
import sys
from pathlib import Path

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


if __name__ == "__main__":
    main()
