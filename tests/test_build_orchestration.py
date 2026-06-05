from __future__ import annotations

from variantfeatures.build import parse_sources, sources_by_consequence


def test_parse_sources_expands_groups():
    assert {"myvariant", "gnomad", "alphamissense", "revel", "alphafold"} <= parse_sources("core")
    assert parse_sources("population,pext") == {"gnomad", "pext_bigwig"}


def test_sources_by_consequence_keeps_missense_only_predictors_on_missense():
    mapping = sources_by_consequence({"alphamissense", "revel", "gnomad", "cadd"})

    assert "alphamissense" in mapping["missense_variant"]
    assert "revel" in mapping["missense_variant"]
    assert "gnomad" in mapping["stop_gained"]
    assert "cadd" in mapping["synonymous_variant"]
    assert "alphamissense" not in mapping["stop_gained"]
    assert "revel" not in mapping["synonymous_variant"]
