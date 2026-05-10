"""Tests for the rule-based NMD trigger/escape classifier."""

from __future__ import annotations

from pathlib import Path

import pytest

from variantfeatures.database import VariantDB
from variantfeatures.handlers.nmd_rules import classify_ptc, annotate_gene
from variantfeatures.transcripts import CdsSegment, Transcript


def _transcript_with_segments(*segs: CdsSegment) -> Transcript:
    """Build a transcript with the given CDS segments and a sane cds_length."""
    cds_length = sum(s.genomic_end - s.genomic_start + 1 for s in segs)
    return Transcript(
        transcript_id="ENST_TEST",
        gene_symbol="TEST",
        chromosome="1",
        strand=1,
        cds_length=cds_length,
        cds_sequence="A" * cds_length,
        cds_segments=list(segs),
    )


def test_ptc_in_last_exon_escapes():
    # Two exons. Last exon is segments[-1] = cds positions 31-90.
    t = _transcript_with_segments(
        CdsSegment(genomic_start=100, genomic_end=129, cds_start=1, cds_end=30),
        CdsSegment(genomic_start=200, genomic_end=259, cds_start=31, cds_end=90),
    )
    # codon 12 -> cds_pos = 34, in last exon
    score, cat = classify_ptc(t, codon_pos=12)
    assert score == 0.0
    assert cat == "last_exon_escape"


def test_ptc_within_50nt_of_last_junction_escapes():
    # Last exon: 91-200; second-to-last: 31-90. 50nt rule applies to second-to-last.
    t = _transcript_with_segments(
        CdsSegment(genomic_start=100, genomic_end=129, cds_start=1, cds_end=30),
        CdsSegment(genomic_start=200, genomic_end=259, cds_start=31, cds_end=90),
        CdsSegment(genomic_start=300, genomic_end=409, cds_start=91, cds_end=200),
    )
    # cds_pos 70: in segments[-2] (31-90), distance to junction = 90 - 70 = 20 (<= 50) -> escape
    # codon 24 -> cds_pos = 70
    score, cat = classify_ptc(t, codon_pos=24)
    assert score == 0.0
    assert cat == "near_last_junction_escape"


def test_ptc_far_from_last_junction_triggers():
    t = _transcript_with_segments(
        CdsSegment(genomic_start=100, genomic_end=129, cds_start=1, cds_end=30),
        CdsSegment(genomic_start=200, genomic_end=259, cds_start=31, cds_end=90),
        CdsSegment(genomic_start=300, genomic_end=409, cds_start=91, cds_end=200),
    )
    # cds_pos 35: in segments[-2], distance to junction = 90 - 35 = 55 (> 50) -> triggers
    # codon 12 -> cds_pos = 34, codon 13 -> cds_pos = 37
    score, cat = classify_ptc(t, codon_pos=13)
    assert score == 1.0
    assert cat == "triggers_nmd"


def test_ptc_in_first_exon_triggers():
    t = _transcript_with_segments(
        CdsSegment(genomic_start=100, genomic_end=129, cds_start=1, cds_end=30),
        CdsSegment(genomic_start=200, genomic_end=259, cds_start=31, cds_end=90),
        CdsSegment(genomic_start=300, genomic_end=409, cds_start=91, cds_end=200),
    )
    # codon 3 -> cds_pos = 7, in segments[0] (first exon)
    score, cat = classify_ptc(t, codon_pos=3)
    assert score == 1.0
    assert cat == "triggers_nmd"


def test_single_exon_transcript_always_escapes():
    """Single-exon genes have no exon-exon junction so all PTCs escape NMD."""
    t = _transcript_with_segments(
        CdsSegment(genomic_start=100, genomic_end=399, cds_start=1, cds_end=300),
    )
    score, cat = classify_ptc(t, codon_pos=10)
    assert score == 0.0
    assert cat == "last_exon_escape"


def test_annotate_gene_writes_pathogenicity_rows(tmp_path):
    """End-to-end DB persistence."""
    db = VariantDB(tmp_path / "test.db")

    # Three-exon transcript so the first exon isn't trivially within 50nt of the last junction.
    # Segments: [1..150], [151..210], [211..300]
    last_exon_var  = db.upsert_variant(chromosome="1", position=999, ref="A", alt="T", variant_type="SNV")
    first_exon_var = db.upsert_variant(chromosome="1", position=110, ref="A", alt="T", variant_type="SNV")
    db.upsert_consequence(
        variant_id=last_exon_var, transcript_id="ENST_T", source="enumerated",
        gene_symbol="TEST", consequence="stop_gained", aa_pos=85,   # cds_pos 253, in last exon
    )
    db.upsert_consequence(
        variant_id=first_exon_var, transcript_id="ENST_T", source="enumerated",
        gene_symbol="TEST", consequence="stop_gained", aa_pos=10,   # cds_pos 28, in first exon, far from last junction
    )

    transcript = _transcript_with_segments(
        CdsSegment(genomic_start=100,  genomic_end=249,  cds_start=1,   cds_end=150),
        CdsSegment(genomic_start=400,  genomic_end=459,  cds_start=151, cds_end=210),
        CdsSegment(genomic_start=600,  genomic_end=689,  cds_start=211, cds_end=300),
    )
    transcript.gene_symbol = "TEST"

    result = annotate_gene(db, "TEST", transcript=transcript)
    assert result["considered"] == 2
    assert result["triggers"] == 1
    assert result["escapes"] == 1
    assert result["by_category"]["last_exon_escape"] == 1
    assert result["by_category"]["triggers_nmd"] == 1

    cur = db.conn.execute(
        "SELECT variant_id, score, category FROM annotations_pathogenicity "
        "WHERE predictor = 'nmd_rule' ORDER BY variant_id"
    )
    rows = [dict(r) for r in cur.fetchall()]
    assert len(rows) == 2
    assert {(r["score"], r["category"]) for r in rows} == {
        (0.0, "last_exon_escape"),
        (1.0, "triggers_nmd"),
    }
