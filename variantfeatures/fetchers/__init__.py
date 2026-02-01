"""Data fetchers for variant annotation sources."""

from .alphamissense import fetch_alphamissense
from .clinvar import fetch_clinvar
from .gnomad import fetch_gnomad

__all__ = ["fetch_alphamissense", "fetch_clinvar", "fetch_gnomad"]
