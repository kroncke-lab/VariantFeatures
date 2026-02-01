"""Fetch gnomAD allele frequencies."""

from typing import Iterator
import requests

# gnomAD API: https://gnomad.broadinstitute.org/api
# GraphQL endpoint for queries

GNOMAD_API = "https://gnomad.broadinstitute.org/api"


def fetch_gnomad(gene: str) -> Iterator[dict]:
    """
    Fetch gnomAD frequencies for variants in a gene.
    
    Yields dicts with:
        - hgvs_p: e.g., "p.Arg528His"
        - gnomad_af: global allele frequency
        - gnomad_af_popmax: max population AF
        - gnomad_homozygotes: count of homozygotes
    """
    # gnomAD has a GraphQL API - example query:
    query = """
    query GeneVariants($geneSymbol: String!) {
      gene(gene_symbol: $geneSymbol, reference_genome: GRCh38) {
        variants(dataset: gnomad_r4) {
          variant_id
          hgvsp
          exome {
            af
            af_popmax
            homozygote_count
          }
          genome {
            af
            af_popmax
            homozygote_count
          }
        }
      }
    }
    """
    # TODO: Implement GraphQL call
    raise NotImplementedError("gnomAD fetcher not yet implemented")
