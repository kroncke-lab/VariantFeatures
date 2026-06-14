"""Compatibility import for the frameshift proxy mapper.

The implementation now lives in :mod:`variantfeatures.frameshift` so it can use
the repo's transcript, enumeration, SQLite, and annotation machinery.
"""

from variantfeatures.frameshift import frameshift_variant_processing

__all__ = ["frameshift_variant_processing"]
