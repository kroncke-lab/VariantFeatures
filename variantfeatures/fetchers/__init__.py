"""Data fetchers for variant annotation sources."""

from .alphamissense import fetch_alphamissense
from .clinvar import fetch_clinvar
from .gnomad import fetch_gnomad
from .lof import fetch_loftee_annotations, classify_lof_type, predict_nmd_escape

__all__ = [
    "fetch_alphamissense", 
    "fetch_clinvar", 
    "fetch_gnomad",
    "fetch_loftee_annotations",
    "classify_lof_type",
    "predict_nmd_escape",
]
