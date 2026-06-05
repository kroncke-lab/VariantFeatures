from __future__ import annotations

from unittest.mock import patch

from variantfeatures.fetchers.gnomad import fetch_single_variant


def test_fetch_single_variant_delegates_to_normalized_handler():
    payload = {
        "variant_id": "7-140753336-A-T",
        "rsids": ["rs113488022"],
        "exome": {"ac": 2, "an": 1460618, "af": 1.369e-6, "homozygote_count": 0},
        "genome": None,
    }

    with patch("variantfeatures.fetchers.gnomad.gnomad_handler.fetch_gnomad", return_value=payload) as fetch:
        result = fetch_single_variant("chr7", 140753336, "A", "T")

    fetch.assert_called_once_with("chr7", 140753336, "A", "T", dataset="gnomad_r4")
    assert result == {
        "variant_id": "7-140753336-A-T",
        "rsids": ["rs113488022"],
        "hgvs_p": None,
        "hgvs_c": None,
        "gnomad_af": 1.369e-6,
        "gnomad_homozygotes": 0,
        "gnomad_an": 1460618,
    }
