"""Rule-based NMD trigger/escape annotation.

Implements the canonical "last exon + 50-nt rule" used by aenmd and most
clinical NMD predictors:

    - Premature termination codons (PTCs) in the LAST CDS exon escape NMD.
    - PTCs within 50 nt upstream of the LAST exon-exon junction (i.e., the
      3' end of the second-to-last CDS exon) also escape.
    - All other PTCs trigger NMD.

The rule is implemented per (gene, transcript) by walking the transcript's
CDS segments. We re-fetch the transcript via Ensembl REST so this works
regardless of whether the saturation enumerator's metadata is still in memory;
the lookup is idempotent and cheap (one REST call per gene).

Output: writes one annotations_pathogenicity row per stop_gained variant with
predictor='nmd_rule', score=1.0 (triggers NMD) or 0.0 (escapes), and a category
string with the specific rule that fired.
"""

from __future__ import annotations

from typing import Iterable, Optional

from ..transcripts import Transcript, fetch_canonical_transcript


SOURCE = "nmd_rule"
PREDICTOR = "nmd_rule"
NMD_50NT_THRESHOLD = 50


class HandlerError(Exception):
    pass


def classify_ptc(transcript: Transcript, codon_pos: int) -> tuple[float, str]:
    """Return (score, category) for a PTC at amino-acid position `codon_pos`.

    score = 1.0  -> triggers NMD
    score = 0.0  -> escapes NMD

    Categories:
        triggers_nmd                  -- standard rule fires
        last_exon_escape              -- PTC in last CDS exon
        near_last_junction_escape     -- PTC within 50nt of last exon-exon junction
    """
    if not transcript.cds_segments:
        return 1.0, "triggers_nmd"

    ptc_cds_pos = (codon_pos - 1) * 3 + 1   # 1-based CDS position of the codon's 1st base
    segments = transcript.cds_segments      # in transcript 5'->3' order

    # Find which segment contains the PTC.
    ptc_segment_idx = None
    for i, seg in enumerate(segments):
        if seg.cds_start <= ptc_cds_pos <= seg.cds_end:
            ptc_segment_idx = i
            break

    last_idx = len(segments) - 1

    if ptc_segment_idx is None:
        # Off the CDS map; default to "triggers" (most conservative for pathogenic call).
        return 1.0, "triggers_nmd"

    # Last-exon rule: any PTC in the terminal CDS exon escapes NMD.
    if ptc_segment_idx == last_idx:
        return 0.0, "last_exon_escape"

    # 50-nt rule: PTC within 50 nt upstream of the last exon-exon junction
    # (i.e., near the 3' end of the segment immediately before the last) escapes.
    if ptc_segment_idx == last_idx - 1:
        seg = segments[ptc_segment_idx]
        distance_to_junction = seg.cds_end - ptc_cds_pos
        if distance_to_junction <= NMD_50NT_THRESHOLD:
            return 0.0, "near_last_junction_escape"

    return 1.0, "triggers_nmd"


def annotate_gene(db, gene_symbol: str, *, transcript: Optional[Transcript] = None) -> dict:
    """Apply the NMD rule to every stop_gained variant of `gene_symbol`.

    Returns: {"considered": int, "triggers": int, "escapes": int, "by_category": dict}
    """
    if transcript is None:
        try:
            transcript = fetch_canonical_transcript(gene_symbol)
        except Exception as e:
            raise HandlerError(f"Could not fetch transcript for {gene_symbol}: {e}")

    cur = db.conn.execute(
        """
        SELECT v.id AS variant_id, c.aa_pos AS aa_pos
        FROM variants v
        JOIN variant_consequences c ON c.variant_id = v.id
        WHERE c.source = 'enumerated' AND c.consequence = 'stop_gained' AND c.gene_symbol = ?
        ORDER BY c.aa_pos
        """,
        [gene_symbol.upper()],
    )
    by_category: dict[str, int] = {}
    triggers = 0
    escapes = 0
    n = 0
    for row in cur.fetchall():
        n += 1
        score, category = classify_ptc(transcript, int(row["aa_pos"]))
        by_category[category] = by_category.get(category, 0) + 1
        if score >= 1.0:
            triggers += 1
        else:
            escapes += 1
        db.upsert_pathogenicity(
            row["variant_id"],
            PREDICTOR,
            score=score,
            category=category,
            source=SOURCE,
        )
    return {
        "considered": n,
        "triggers": triggers,
        "escapes": escapes,
        "by_category": by_category,
    }
