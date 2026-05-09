"""Tests for saturation enumeration."""

from __future__ import annotations

from pathlib import Path

import pytest

from variantfeatures.database import VariantDB
from variantfeatures.enumerate import (
    CODON_TABLE,
    enumerate_snvs,
    normalize_types,
    populate_for_transcript,
)
from variantfeatures.transcripts import CdsSegment, Transcript


# ---------------------------------------------------------------------------
# normalize_types
# ---------------------------------------------------------------------------

def test_normalize_types_friendly_names():
    assert normalize_types(["missense", "nonsense"]) == {"missense_variant", "stop_gained"}


def test_normalize_types_passthrough_so_terms():
    assert normalize_types(["missense_variant"]) == {"missense_variant"}


def test_normalize_types_string_input():
    assert normalize_types("missense,nonsense") == {"missense_variant", "stop_gained"}


def test_normalize_types_unknown_raises():
    with pytest.raises(ValueError):
        normalize_types(["banana"])


# ---------------------------------------------------------------------------
# Codon table sanity
# ---------------------------------------------------------------------------

def test_codon_table_has_64_entries():
    assert len(CODON_TABLE) == 64


def test_codon_table_stops():
    assert CODON_TABLE["TAA"] == "*"
    assert CODON_TABLE["TAG"] == "*"
    assert CODON_TABLE["TGA"] == "*"


def test_codon_table_start():
    assert CODON_TABLE["ATG"] == "M"


# ---------------------------------------------------------------------------
# Enumeration on a hand-crafted minimal transcript
# ---------------------------------------------------------------------------

def _make_plus_strand_transcript(cds: str, genomic_start: int = 100) -> Transcript:
    """One-exon + strand transcript: CDS at genomic [start, start+len-1]."""
    n = len(cds)
    return Transcript(
        transcript_id="ENST_PLUS",
        transcript_version="1",
        gene_symbol="TESTG",
        chromosome="1",
        strand=1,
        cds_length=n,
        cds_sequence=cds,
        cds_segments=[CdsSegment(genomic_start=genomic_start, genomic_end=genomic_start + n - 1, cds_start=1, cds_end=n)],
        is_canonical=True,
        is_mane_select=True,
    )


def _make_minus_strand_transcript(cds: str, genomic_end: int = 200) -> Transcript:
    """One-exon - strand transcript: CDS at genomic [end-len+1, end]."""
    n = len(cds)
    return Transcript(
        transcript_id="ENST_MINUS",
        transcript_version="1",
        gene_symbol="TESTG",
        chromosome="1",
        strand=-1,
        cds_length=n,
        cds_sequence=cds,
        cds_segments=[CdsSegment(genomic_start=genomic_end - n + 1, genomic_end=genomic_end, cds_start=1, cds_end=n)],
        is_canonical=True,
        is_mane_select=True,
    )


def test_enumerate_three_codon_plus_strand_missense_and_nonsense():
    # ATG GAG TAA -> Met-Glu-Stop
    cds = "ATGGAGTAA"
    t = _make_plus_strand_transcript(cds, genomic_start=100)

    snvs = list(enumerate_snvs(t, types=("missense", "nonsense")))

    # Codon 1 (Met) at cds 1-3 (genomic 100-102) is the start codon.
    # ATG -> any other codon = start_lost (not missense), so excluded.
    # Codon 2 (Glu, GAG, cds 4-6, genomic 103-105) -> generates missense + 1 nonsense (E -> *).
    # Codon 3 (Stop, TAA, cds 7-9) -> only stop_lost variants, excluded.

    # GAG at codon 2:
    #   pos1 (G): G>A=AAG K (mis), G>C=CAG Q (mis), G>T=TAG * (nonsense)
    #   pos2 (A): A>C=GCG A (mis), A>G=GGG G (mis), A>T=GTG V (mis)
    #   pos3 (G): G>A=GAA E (synonymous, excluded), G>C=GAC D (mis), G>T=GAT D (mis)
    # Total kept: 7 missense + 1 nonsense = 8

    assert len(snvs) == 8

    by_consequence = {}
    for s in snvs:
        by_consequence.setdefault(s.consequence, []).append(s)
    assert len(by_consequence["missense_variant"]) == 7
    assert len(by_consequence["stop_gained"]) == 1

    nonsense = by_consequence["stop_gained"][0]
    assert nonsense.aa_pos == 2
    assert nonsense.aa_ref == "E"
    assert nonsense.aa_alt == "*"
    # GAG -> TAG is at codon position 1 (CDS pos 4), not pos 3.
    # plus strand, genomic_start=100 -> genomic = 100 + (4-1) = 103
    assert nonsense.cds_pos == 4
    assert nonsense.position == 103
    assert nonsense.ref == "G"
    assert nonsense.alt == "T"
    assert nonsense.codon_alt == "TAG"
    assert nonsense.hgvs_c == "ENST_PLUS.1:c.4G>T"
    assert nonsense.hgvs_p == "ENST_PLUS.1:p.Glu2Ter"


def test_enumerate_minus_strand_emits_complement_in_vcf():
    # CDS in transcript orientation, but genomic + strand has the complement.
    # ATG GAG TAA -> on - strand at genomic 200..192 (CDS 1..9 maps high->low)
    cds = "ATGGAGTAA"
    t = _make_minus_strand_transcript(cds, genomic_end=200)

    snvs = list(enumerate_snvs(t, types=("missense", "nonsense")))

    nonsense = next(s for s in snvs if s.consequence == "stop_gained")
    # GAG -> TAG is at codon position 1 (CDS pos 4).
    # minus strand, genomic_end=200 -> genomic = 200 - (4-1) = 197
    assert nonsense.cds_pos == 4
    assert nonsense.position == 197
    # CDS base G -> + strand ref C (complement); CDS alt T -> + strand alt A
    assert nonsense.ref == "C"
    assert nonsense.alt == "A"
    assert nonsense.hgvs_c == "ENST_MINUS.1:c.4G>T"
    assert nonsense.hgvs_p == "ENST_MINUS.1:p.Glu2Ter"


def test_enumerate_synonymous_excluded_by_default():
    cds = "ATGGAGTAA"
    t = _make_plus_strand_transcript(cds)
    snvs = list(enumerate_snvs(t, types=("missense", "nonsense")))
    assert all(s.consequence != "synonymous_variant" for s in snvs)


def test_enumerate_synonymous_included_when_requested():
    cds = "ATGGAGTAA"
    t = _make_plus_strand_transcript(cds)
    snvs = list(enumerate_snvs(t, types=("synonymous",)))
    # Only GAG p3 G>A (E -> E) is a true synonymous_variant.
    # TAA -> TGA / TAG are stop->stop, classified separately as stop_retained_variant.
    assert all(s.consequence == "synonymous_variant" for s in snvs)
    assert len(snvs) == 1
    assert snvs[0].cds_pos == 6
    assert snvs[0].codon_alt == "GAA"


def test_enumerate_start_lost_distinct_from_missense():
    cds = "ATGGAGTAA"
    t = _make_plus_strand_transcript(cds)
    snvs = list(enumerate_snvs(t, types=("start_lost",)))
    # ATG codon: 9 single-nt changes total, 6 yield non-Met -> 6 start_lost (3 are Met-self or syn)
    # Actually ATG -> ATA, ATC, ATT (Ile), ACG (Thr), AAG (Lys), AGG (Arg), CTG (Leu), GTG (Val), TTG (Leu)
    # That's 9 alts. None equal Met. So all 9 are start_lost.
    assert len(snvs) == 9
    assert all(s.consequence == "start_lost" for s in snvs)
    assert all(s.aa_pos == 1 for s in snvs)


def test_enumerate_skips_codon_not_in_genetic_code():
    """Pretend we have a CDS containing N — those codons are skipped, not crashed."""
    cds = "ATGNNN"
    t = _make_plus_strand_transcript(cds)
    snvs = list(enumerate_snvs(t, types=("missense", "nonsense")))
    # Codon 1 = ATG (start_lost only, excluded). Codon 2 = NNN -> aa_ref None, skipped entirely.
    assert snvs == []


# ---------------------------------------------------------------------------
# populate_for_transcript persistence + job queue
# ---------------------------------------------------------------------------

@pytest.fixture
def db(tmp_path: Path) -> VariantDB:
    return VariantDB(tmp_path / "test.db")


def test_populate_persists_variants_consequences_and_jobs(db):
    cds = "ATGGAGTAA"
    t = _make_plus_strand_transcript(cds)
    summary = populate_for_transcript(t, db, types=("missense", "nonsense"))

    assert summary["variants"] == 8
    assert summary["consequences"] == 8
    assert summary["jobs_queued"] > 0
    assert summary["by_consequence"]["missense_variant"] == 7
    assert summary["by_consequence"]["stop_gained"] == 1

    # variants table has 8 rows
    cur = db.conn.execute("SELECT COUNT(*) AS n FROM variants")
    assert cur.fetchone()["n"] == 8

    # consequences table has 8 rows tagged enumerated
    cur = db.conn.execute("SELECT COUNT(*) AS n FROM variant_consequences WHERE source = 'enumerated'")
    assert cur.fetchone()["n"] == 8

    # Aliases exist for each variant (hgvs_c at minimum, often hgvs_p too)
    cur = db.conn.execute("SELECT COUNT(*) AS n FROM variant_aliases WHERE alias_type = 'hgvs_c'")
    assert cur.fetchone()["n"] == 8

    # Annotation jobs were queued
    cur = db.conn.execute("SELECT COUNT(*) AS n FROM annotation_jobs WHERE status = 'pending'")
    assert cur.fetchone()["n"] == summary["jobs_queued"]


def test_populate_no_enqueue(db):
    cds = "ATGGAGTAA"
    t = _make_plus_strand_transcript(cds)
    summary = populate_for_transcript(t, db, types=("missense", "nonsense"), enqueue=False)
    assert summary["jobs_queued"] == 0
    cur = db.conn.execute("SELECT COUNT(*) AS n FROM annotation_jobs")
    assert cur.fetchone()["n"] == 0


def test_populate_idempotent(db):
    cds = "ATGGAGTAA"
    t = _make_plus_strand_transcript(cds)
    s1 = populate_for_transcript(t, db, types=("missense", "nonsense"))
    s2 = populate_for_transcript(t, db, types=("missense", "nonsense"))
    assert s1["variants"] == s2["variants"]
    cur = db.conn.execute("SELECT COUNT(*) AS n FROM variants")
    assert cur.fetchone()["n"] == 8  # not doubled


def test_populate_respects_limit(db):
    cds = "ATGGAGTAA"
    t = _make_plus_strand_transcript(cds)
    summary = populate_for_transcript(t, db, types=("missense", "nonsense"), limit=3)
    assert summary["variants"] == 3
