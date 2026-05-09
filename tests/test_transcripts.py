"""Tests for transcripts.py — CDS-to-genomic mapping and Ensembl REST."""

from __future__ import annotations

import os

import pytest

from variantfeatures.transcripts import (
    CdsSegment,
    Transcript,
    TranscriptError,
    _build_cds_segments,
    _choose_canonical_transcript,
    fetch_canonical_transcript,
    reverse_complement,
)


# ---------------------------------------------------------------------------
# CDS segment construction
# ---------------------------------------------------------------------------

def test_build_cds_segments_plus_strand_single_exon():
    exons = [{"start": 100, "end": 200}]
    segs = _build_cds_segments(exons, cds_genomic_start=110, cds_genomic_end=190, strand=1)
    assert len(segs) == 1
    assert segs[0].genomic_start == 110
    assert segs[0].genomic_end == 190
    assert segs[0].cds_start == 1
    assert segs[0].cds_end == 81


def test_build_cds_segments_plus_strand_multi_exon():
    exons = [
        {"start": 100, "end": 110},   # exon 1
        {"start": 200, "end": 220},   # exon 2 (with intron)
        {"start": 300, "end": 305},   # exon 3
    ]
    # CDS spans from 105 to 303 -> last 6bp of exon1, all of exon2, first 4bp of exon3
    segs = _build_cds_segments(exons, cds_genomic_start=105, cds_genomic_end=303, strand=1)
    assert [(s.genomic_start, s.genomic_end, s.cds_start, s.cds_end) for s in segs] == [
        (105, 110, 1, 6),
        (200, 220, 7, 27),
        (300, 303, 28, 31),
    ]


def test_build_cds_segments_minus_strand():
    exons = [
        {"start": 100, "end": 110},
        {"start": 200, "end": 220},
        {"start": 300, "end": 305},
    ]
    # Same genomic CDS interval, but minus strand -> CDS walks high->low genomically.
    segs = _build_cds_segments(exons, cds_genomic_start=105, cds_genomic_end=303, strand=-1)
    assert [(s.genomic_start, s.genomic_end, s.cds_start, s.cds_end) for s in segs] == [
        (300, 303, 1, 4),
        (200, 220, 5, 25),
        (105, 110, 26, 31),
    ]


# ---------------------------------------------------------------------------
# Transcript.cds_to_genomic
# ---------------------------------------------------------------------------

def _make_transcript(strand: int, segments: list[CdsSegment], cds_seq: str = "") -> Transcript:
    cds_length = sum(s.genomic_end - s.genomic_start + 1 for s in segments)
    return Transcript(
        transcript_id="ENST_TEST",
        strand=strand,
        cds_length=cds_length,
        cds_sequence=cds_seq or ("A" * cds_length),
        cds_segments=segments,
    )


def test_cds_to_genomic_plus_strand_single_segment():
    t = _make_transcript(1, [CdsSegment(genomic_start=100, genomic_end=199, cds_start=1, cds_end=100)])
    assert t.cds_to_genomic(1) == 100
    assert t.cds_to_genomic(50) == 149
    assert t.cds_to_genomic(100) == 199


def test_cds_to_genomic_plus_strand_multi_segment_jumps_intron():
    t = _make_transcript(1, [
        CdsSegment(genomic_start=100, genomic_end=109, cds_start=1, cds_end=10),
        CdsSegment(genomic_start=200, genomic_end=219, cds_start=11, cds_end=30),
    ])
    assert t.cds_to_genomic(10) == 109
    assert t.cds_to_genomic(11) == 200  # crossed the intron
    assert t.cds_to_genomic(30) == 219


def test_cds_to_genomic_minus_strand_single_segment():
    # CDS pos 1 corresponds to the highest genomic coord of the segment.
    t = _make_transcript(-1, [CdsSegment(genomic_start=100, genomic_end=199, cds_start=1, cds_end=100)])
    assert t.cds_to_genomic(1) == 199
    assert t.cds_to_genomic(2) == 198
    assert t.cds_to_genomic(100) == 100


def test_cds_to_genomic_minus_strand_multi_segment():
    t = _make_transcript(-1, [
        CdsSegment(genomic_start=300, genomic_end=303, cds_start=1, cds_end=4),
        CdsSegment(genomic_start=200, genomic_end=210, cds_start=5, cds_end=15),
    ])
    assert t.cds_to_genomic(1) == 303
    assert t.cds_to_genomic(4) == 300
    assert t.cds_to_genomic(5) == 210  # jumped across intron going downward
    assert t.cds_to_genomic(15) == 200


def test_cds_to_genomic_out_of_range_raises():
    t = _make_transcript(1, [CdsSegment(genomic_start=100, genomic_end=109, cds_start=1, cds_end=10)])
    with pytest.raises(TranscriptError):
        t.cds_to_genomic(0)
    with pytest.raises(TranscriptError):
        t.cds_to_genomic(11)


# ---------------------------------------------------------------------------
# genomic_ref_for_cds_base
# ---------------------------------------------------------------------------

def test_genomic_ref_plus_strand_is_identity():
    t = _make_transcript(1, [CdsSegment(1, 1, 1, 1)])
    for b in "ACGT":
        assert t.genomic_ref_for_cds_base(b) == b


def test_genomic_ref_minus_strand_is_complement():
    t = _make_transcript(-1, [CdsSegment(1, 1, 1, 1)])
    assert t.genomic_ref_for_cds_base("A") == "T"
    assert t.genomic_ref_for_cds_base("C") == "G"
    assert t.genomic_ref_for_cds_base("G") == "C"
    assert t.genomic_ref_for_cds_base("T") == "A"


def test_reverse_complement():
    assert reverse_complement("ACGT") == "ACGT"
    assert reverse_complement("AAAA") == "TTTT"
    assert reverse_complement("ATG") == "CAT"


# ---------------------------------------------------------------------------
# _choose_canonical_transcript
# ---------------------------------------------------------------------------

def test_choose_canonical_prefers_mane_select():
    transcripts = [
        {"id": "T1", "is_canonical": True, "Translation": {"start": 1, "end": 100}},
        {"id": "T2", "MANE": [{"type": "MANE_Select"}], "Translation": {"start": 1, "end": 50}},
    ]
    assert _choose_canonical_transcript(transcripts)["id"] == "T2"


def test_choose_canonical_fallback_to_is_canonical():
    transcripts = [
        {"id": "T1", "is_canonical": True, "Translation": {"start": 1, "end": 100}},
        {"id": "T2", "Translation": {"start": 1, "end": 200}},
    ]
    assert _choose_canonical_transcript(transcripts)["id"] == "T1"


def test_choose_canonical_fallback_to_longest_coding():
    transcripts = [
        {"id": "T1", "biotype": "protein_coding", "Translation": {"start": 1, "end": 100}},
        {"id": "T2", "biotype": "protein_coding", "Translation": {"start": 1, "end": 200}},
        {"id": "T3", "biotype": "lincRNA"},
    ]
    assert _choose_canonical_transcript(transcripts)["id"] == "T2"


# ---------------------------------------------------------------------------
# Live integration (KCNH2)
# ---------------------------------------------------------------------------
@pytest.mark.skipif(
    os.environ.get("RUN_INTEGRATION") != "1",
    reason="Set RUN_INTEGRATION=1 to hit Ensembl REST.",
)
def test_fetch_canonical_transcript_kcnh2_live():
    t = fetch_canonical_transcript("KCNH2")
    assert t.gene_symbol == "KCNH2"
    assert t.chromosome == "7"
    assert t.strand == -1
    assert t.cds_length % 3 == 0
    # KCNH2 protein is 1159 aa -> CDS = 1159*3 + 3 (stop) = 3480 bp
    assert t.cds_length == 3480
    # CDS starts with ATG and ends with a stop
    assert t.cds_sequence.startswith("ATG")
    assert t.cds_sequence[-3:] in {"TAA", "TAG", "TGA"}
    # And the segments cover exactly cds_length
    assert sum(s.genomic_end - s.genomic_start + 1 for s in t.cds_segments) == t.cds_length
