"""Transcript metadata + CDS-to-genomic mapping via Ensembl REST.

Used by the saturation enumerator to walk a CDS codon-by-codon and emit
canonical-strand (chrom, pos, ref, alt) tuples for every possible SNV.

Coordinates returned are GRCh38, 1-based inclusive (matching VCF / HGVS).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import requests


ENSEMBL_REST = "https://rest.ensembl.org"
DEFAULT_TIMEOUT = 30
JSON_HEADERS = {"Accept": "application/json"}

_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


class TranscriptError(Exception):
    """Raised when transcript metadata cannot be retrieved or interpreted."""


def reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


@dataclass
class CdsSegment:
    """One contiguous stretch of CDS that lives inside a single exon.

    `cds_start` and `cds_end` are 1-based positions in the transcript-oriented
    CDS sequence. `genomic_start` <= `genomic_end` always (genomic coords).
    """

    genomic_start: int
    genomic_end: int
    cds_start: int
    cds_end: int


@dataclass
class Transcript:
    transcript_id: str            # Ensembl, e.g. ENST00000262186
    transcript_version: Optional[str] = None
    refseq_match: Optional[str] = None  # MANE Select RefSeq, e.g. NM_000238.4
    protein_id: Optional[str] = None    # Ensembl ENSP, used for hgvs_p prefix
    gene_symbol: Optional[str] = None
    gene_ensembl: Optional[str] = None
    chromosome: Optional[str] = None
    strand: int = 1               # +1 or -1
    cds_length: int = 0
    cds_sequence: str = ""        # transcript orientation (5'->3')
    cds_segments: list[CdsSegment] = field(default_factory=list)
    is_canonical: bool = False
    is_mane_select: bool = False
    is_mane_plus_clinical: bool = False

    @property
    def transcript_id_versioned(self) -> str:
        if self.transcript_version:
            return f"{self.transcript_id}.{self.transcript_version}"
        return self.transcript_id

    def cds_to_genomic(self, cds_pos: int) -> int:
        """Map a 1-based CDS position to a 1-based genomic position.

        Raises TranscriptError if cds_pos is outside the CDS.
        """
        if cds_pos < 1 or cds_pos > self.cds_length:
            raise TranscriptError(f"CDS position {cds_pos} out of range [1, {self.cds_length}]")
        for seg in self.cds_segments:
            if seg.cds_start <= cds_pos <= seg.cds_end:
                offset = cds_pos - seg.cds_start
                if self.strand == 1:
                    return seg.genomic_start + offset
                return seg.genomic_end - offset
        raise TranscriptError(f"CDS position {cds_pos} not within any segment")

    def genomic_ref_for_cds_base(self, cds_base: str) -> str:
        """Translate a CDS base (transcript orientation) to the + strand reference base."""
        if self.strand == 1:
            return cds_base.upper()
        return cds_base.translate(_COMPLEMENT).upper()


# ---------------------------------------------------------------------------
# Ensembl REST fetch
# ---------------------------------------------------------------------------

def fetch_canonical_transcript(
    gene_symbol: str,
    *,
    species: str = "homo_sapiens",
    timeout: int = DEFAULT_TIMEOUT,
) -> Transcript:
    """Look up a gene's canonical transcript and build a complete Transcript object.

    Prefers the MANE Select transcript if Ensembl exposes one, falling back to
    `is_canonical`.
    """
    payload = _ensembl_get(
        f"{ENSEMBL_REST}/lookup/symbol/{species}/{gene_symbol}",
        params={"expand": 1, "mane": 1},
        timeout=timeout,
    )
    transcripts = payload.get("Transcript") or []
    if not transcripts:
        raise TranscriptError(f"No transcripts returned for gene {gene_symbol!r}")

    chosen = _choose_canonical_transcript(transcripts)
    return _build_transcript(chosen, gene_symbol=gene_symbol, gene_ensembl=payload.get("id"), timeout=timeout)


def _choose_canonical_transcript(transcripts: list[dict]) -> dict:
    # MANE Select wins if we have it
    for t in transcripts:
        mane = t.get("MANE")
        if isinstance(mane, list):
            for m in mane:
                if (m.get("type") or "").lower() == "mane_select":
                    return t
    # else: Ensembl-flagged canonical
    for t in transcripts:
        if t.get("is_canonical"):
            return t
    # else: longest CDS-bearing transcript
    coding = [t for t in transcripts if t.get("biotype") == "protein_coding" and t.get("Translation")]
    if coding:
        return max(coding, key=lambda t: (t.get("Translation", {}).get("end", 0) - t.get("Translation", {}).get("start", 0)))
    raise TranscriptError("No coding transcript found")


def _build_transcript(t: dict, *, gene_symbol: Optional[str], gene_ensembl: Optional[str], timeout: int) -> Transcript:
    translation = t.get("Translation")
    if not translation:
        raise TranscriptError(f"Transcript {t.get('id')} has no Translation section")
    cds_genomic_start = translation["start"]
    cds_genomic_end = translation["end"]
    if cds_genomic_start > cds_genomic_end:
        cds_genomic_start, cds_genomic_end = cds_genomic_end, cds_genomic_start

    exons = t.get("Exon") or []
    if not exons:
        raise TranscriptError(f"Transcript {t.get('id')} has no Exon section")

    strand = int(t.get("strand", 1))
    chromosome = t.get("seq_region_name") or t.get("chromosome") or ""

    cds_segments = _build_cds_segments(exons, cds_genomic_start, cds_genomic_end, strand)
    cds_length = sum(seg.genomic_end - seg.genomic_start + 1 for seg in cds_segments)

    transcript_id = t["id"]
    transcript_version = str(t.get("version")) if t.get("version") is not None else None

    is_mane_select = False
    is_mane_plus_clinical = False
    refseq_match = None
    for m in t.get("MANE") or []:
        if not isinstance(m, dict):
            continue
        kind = (m.get("type") or "").lower()
        if kind == "mane_select":
            is_mane_select = True
            refseq_match = m.get("refseq_match") or refseq_match
        elif kind == "mane_plus_clinical":
            is_mane_plus_clinical = True
            refseq_match = refseq_match or m.get("refseq_match")

    cds_sequence = _fetch_cds_sequence(transcript_id, timeout=timeout)
    if len(cds_sequence) != cds_length:
        # Coordinate math should match the sequence length; warn loudly via exception
        raise TranscriptError(
            f"CDS length mismatch for {transcript_id}: segments imply {cds_length} bp, "
            f"sequence is {len(cds_sequence)} bp"
        )

    protein_id = translation.get("id")

    return Transcript(
        transcript_id=transcript_id,
        transcript_version=transcript_version,
        refseq_match=refseq_match,
        protein_id=protein_id,
        gene_symbol=gene_symbol,
        gene_ensembl=gene_ensembl,
        chromosome=str(chromosome),
        strand=strand,
        cds_length=cds_length,
        cds_sequence=cds_sequence.upper(),
        cds_segments=cds_segments,
        is_canonical=bool(t.get("is_canonical")),
        is_mane_select=is_mane_select,
        is_mane_plus_clinical=is_mane_plus_clinical,
    )


def _build_cds_segments(
    exons: list[dict],
    cds_genomic_start: int,
    cds_genomic_end: int,
    strand: int,
) -> list[CdsSegment]:
    """Intersect each exon with the CDS interval, then assign CDS positions in transcript order."""
    intervals: list[tuple[int, int]] = []
    for ex in exons:
        ex_start = ex.get("start")
        ex_end = ex.get("end")
        if ex_start is None or ex_end is None:
            continue
        lo = max(ex_start, cds_genomic_start)
        hi = min(ex_end, cds_genomic_end)
        if lo <= hi:
            intervals.append((lo, hi))

    intervals.sort(key=lambda x: x[0])
    if strand == -1:
        intervals.reverse()

    segments: list[CdsSegment] = []
    cursor = 1
    for lo, hi in intervals:
        length = hi - lo + 1
        segments.append(CdsSegment(
            genomic_start=lo,
            genomic_end=hi,
            cds_start=cursor,
            cds_end=cursor + length - 1,
        ))
        cursor += length
    return segments


def _fetch_cds_sequence(transcript_id: str, *, timeout: int) -> str:
    payload = _ensembl_get(
        f"{ENSEMBL_REST}/sequence/id/{transcript_id}",
        params={"type": "cds"},
        timeout=timeout,
    )
    seq = payload.get("seq")
    if not seq:
        raise TranscriptError(f"No CDS sequence returned for transcript {transcript_id!r}")
    return seq


def _ensembl_get(url: str, *, params: Optional[dict] = None, timeout: int) -> dict:
    resp = requests.get(url, params=params or {}, headers=JSON_HEADERS, timeout=timeout)
    if resp.status_code == 429:
        raise TranscriptError("Ensembl REST rate-limited (HTTP 429); back off and retry")
    if resp.status_code == 404:
        raise TranscriptError(f"Ensembl REST returned 404 for {url}")
    resp.raise_for_status()
    data = resp.json()
    if isinstance(data, dict) and data.get("error"):
        raise TranscriptError(f"Ensembl REST error: {data['error']}")
    return data
