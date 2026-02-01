"""Fetch LOF-specific annotations (LOFTEE, NMD prediction)."""

from typing import Iterator

# LOFTEE is typically run as part of VEP
# For pre-computed annotations, gnomAD provides LOFTEE flags

GNOMAD_LOF_FIELDS = [
    "lof",           # HC, LC, or empty
    "lof_filter",    # LOFTEE filter reasons
    "lof_flags",     # Additional flags
    "lof_info",      # Free text info
]


def classify_lof_type(hgvs_c: str, hgvs_p: str, consequence: str) -> str:
    """
    Classify LOF variant type from HGVS or VEP consequence.
    
    Returns: 'nonsense', 'frameshift', 'splice_donor', 'splice_acceptor', or None
    """
    if hgvs_p and "Ter" in hgvs_p and "fs" not in hgvs_p:
        return "nonsense"
    if hgvs_p and "fs" in hgvs_p:
        return "frameshift"
    if consequence:
        if "splice_donor" in consequence:
            return "splice_donor"
        if "splice_acceptor" in consequence:
            return "splice_acceptor"
        if "stop_gained" in consequence:
            return "nonsense"
        if "frameshift" in consequence:
            return "frameshift"
    return None


def predict_nmd_escape(truncation_position: float, last_exon: bool) -> bool:
    """
    Predict whether a truncating variant escapes nonsense-mediated decay.
    
    NMD escape typically occurs when:
    - Variant is in the last exon
    - Variant is within 50-55bp of the last exon-exon junction
    - Truncation is near the end of the protein (>90% remaining)
    
    Args:
        truncation_position: Fraction of protein remaining (0-1)
        last_exon: Whether variant is in the last exon
        
    Returns:
        True if likely to escape NMD
    """
    if last_exon:
        return True
    if truncation_position > 0.90:
        return True
    return False


def calculate_truncation_position(variant_position: int, protein_length: int) -> float:
    """
    Calculate fraction of protein remaining after truncation.
    
    Args:
        variant_position: Amino acid position of truncation
        protein_length: Total protein length
        
    Returns:
        Fraction remaining (0-1)
    """
    if protein_length <= 0:
        return 0.0
    return variant_position / protein_length


def fetch_loftee_annotations(gene: str) -> Iterator[dict]:
    """
    Fetch LOFTEE annotations for LOF variants in a gene.
    
    These come from gnomAD's pre-computed VEP annotations.
    
    Yields dicts with:
        - hgvs_c: e.g., "c.1234C>T"
        - lof_type: nonsense, frameshift, splice_donor, splice_acceptor
        - loftee_confidence: HC or LC
        - loftee_flags: comma-separated flags
        - nmd_escape: boolean prediction
    """
    # TODO: Implement - parse gnomAD VCF/Hail table for LOFTEE fields
    raise NotImplementedError("LOFTEE fetcher not yet implemented")
