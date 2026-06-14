from __future__ import annotations

from pathlib import Path

from variantfeatures.build import build_gene, parse_sources, sources_by_consequence
from variantfeatures.database import VariantDB
from variantfeatures.transcripts import CdsSegment, Transcript


def test_parse_sources_expands_groups():
    assert {
        "myvariant",
        "gnomad",
        "clinvar",
        "alphamissense",
        "revel",
        "alphafold",
        "frameshift_proxy",
    } <= parse_sources("core")
    assert parse_sources("clinvar") == {"clinvar", "myvariant"}
    assert parse_sources("population,pext") == {"gnomad", "pext_bigwig"}
    assert "frameshift_proxy" in parse_sources("lof")


def test_sources_by_consequence_keeps_missense_only_predictors_on_missense():
    mapping = sources_by_consequence({"alphamissense", "revel", "gnomad", "cadd"})

    assert "alphamissense" in mapping["missense_variant"]
    assert "revel" in mapping["missense_variant"]
    assert "gnomad" in mapping["stop_gained"]
    assert "cadd" in mapping["synonymous_variant"]
    assert "alphamissense" not in mapping["stop_gained"]
    assert "revel" not in mapping["synonymous_variant"]


def test_build_gene_enumerates_selected_isoforms(monkeypatch, tmp_path: Path):
    transcripts = [
        _test_transcript("ENST_CANON", genomic_start=100, is_canonical=True),
        _test_transcript("ENST_ALT", genomic_start=200),
    ]
    monkeypatch.setattr(
        "variantfeatures.build.fetch_gene_transcripts",
        lambda gene, include, max_transcripts: transcripts,
    )

    result = build_gene(
        "KCNH2",
        db_path=tmp_path / "test.db",
        sources="gnomad",
        types="missense",
        enumerate_limit=1,
        run_annotations=False,
        isoforms="all",
    )

    assert result.transcript_labels == ["ENST_CANON.1", "ENST_ALT.1"]
    assert set(result.enumeration["by_transcript"]) == {"ENST_CANON.1", "ENST_ALT.1"}
    db = VariantDB(tmp_path / "test.db")
    assert {row["transcript_id"] for row in db.get_transcripts("KCNH2")} == {
        "ENST_CANON.1",
        "ENST_ALT.1",
    }


def _test_transcript(transcript_id: str, *, genomic_start: int, is_canonical: bool = False) -> Transcript:
    cds = "ATGGAGTAA"
    return Transcript(
        transcript_id=transcript_id,
        transcript_version="1",
        gene_symbol="KCNH2",
        gene_ensembl="ENSG_TEST",
        chromosome="7",
        strand=1,
        cds_length=len(cds),
        cds_sequence=cds,
        cds_segments=[
            CdsSegment(
                genomic_start=genomic_start,
                genomic_end=genomic_start + len(cds) - 1,
                cds_start=1,
                cds_end=len(cds),
            )
        ],
        is_canonical=is_canonical,
    )
