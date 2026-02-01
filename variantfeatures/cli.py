"""Command-line interface for VariantFeatures."""

import click
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
@click.option("--sources", default="all", help="Sources: alphamissense,revel,cadd,clinvar,gnomad")
def build(genes: str, db: str, sources: str):
    """Fetch features and build/update database."""
    gene_list = [g.strip().upper() for g in genes.split(",")]
    db_path = Path(db) if db else None
    
    click.echo(f"Building database for genes: {', '.join(gene_list)}")
    
    # TODO: Implement fetchers
    click.echo("Fetchers not yet implemented. Structure ready.")


@main.command()
@click.option("--gene", "-g", required=True, help="Gene symbol")
@click.option("--db", type=click.Path(), default=None, help="Database path")
@click.option("--format", "fmt", type=click.Choice(["csv", "json", "table"]), default="table")
def query(gene: str, db: str, fmt: str):
    """Query variants for a gene."""
    db_path = Path(db) if db else None
    vdb = VariantDB(db_path)
    
    variants = vdb.get_gene_variants(gene.upper())
    
    if not variants:
        click.echo(f"No variants found for {gene}")
        return
    
    click.echo(f"Found {len(variants)} variants for {gene}")
    # TODO: Format output
    for v in variants[:10]:
        click.echo(f"  {v['hgvs_p']}: AM={v['alphamissense_score']}, REVEL={v['revel_score']}")


if __name__ == "__main__":
    main()
