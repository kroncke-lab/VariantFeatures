"""Map frameshift positions to finite stop-gained SNV proxy variants.

Frameshift variants can be represented at many inserted/deleted sequence
states, but most downstream predictors in this repo are keyed by concrete
SNVs. This module maps an amino-acid position affected by a frameshift to the
nearest codon whose reference sequence can become a stop codon by one SNV, then
persists the concrete stop-gained SNV(s) as proxy variants.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Optional

from .enumerate import CODON_TABLE, EnumeratedSNV, _format_hgvs_c, _format_hgvs_p
from .enumerate import populate_for_transcript
from .transcripts import Transcript, fetch_canonical_transcript


SOURCE = "frameshift_proxy"
MAPPING_SOURCE = "frameshift_nearest_stop_snv"
DEFAULT_MAX_STEPS = 20
STOP_CODONS = frozenset(codon for codon, aa in CODON_TABLE.items() if aa == "*")


@dataclass(frozen=True)
class FrameshiftNonsenseProxy:
    """One mapped frameshift position -> one concrete stop-gained SNV."""

    frameshift_aa_pos: int
    direction: str
    n_steps_wrt_frameshift: int
    proxy: EnumeratedSNV

    def to_record(self, *, variant_id: Optional[int] = None) -> dict:
        return {
            "frameshift_aa_pos": self.frameshift_aa_pos,
            "direction": self.direction,
            "n_steps_wrt_frameshift": self.n_steps_wrt_frameshift,
            "proxy_variant_id": variant_id,
            "proxy_chromosome": self.proxy.chromosome,
            "proxy_position": self.proxy.position,
            "proxy_ref": self.proxy.ref,
            "proxy_alt": self.proxy.alt,
            "proxy_hgvs_c": self.proxy.hgvs_c,
            "proxy_hgvs_p": self.proxy.hgvs_p,
            "proxy_aa_pos": self.proxy.aa_pos,
            "proxy_aa_ref": self.proxy.aa_ref,
            "proxy_codon_ref": self.proxy.codon_ref,
            "proxy_codon_alt": self.proxy.codon_alt,
        }


def protein_codon_count(transcript: Transcript) -> int:
    """Return the number of protein-coding codons, excluding terminal stop."""
    n_codons = len(transcript.cds_sequence) // 3
    if n_codons and CODON_TABLE.get(transcript.cds_sequence[(n_codons - 1) * 3 : n_codons * 3]) == "*":
        return n_codons - 1
    return n_codons


def codons_one_snv_from_stop() -> frozenset[str]:
    """Return non-stop codons that can become a stop codon via one SNV."""
    out = []
    for codon, aa in CODON_TABLE.items():
        if aa == "*":
            continue
        for offset, base in enumerate(codon):
            for alt_base in "ACGT":
                if alt_base == base:
                    continue
                alt_codon = codon[:offset] + alt_base + codon[offset + 1 :]
                if alt_codon in STOP_CODONS:
                    out.append(codon)
                    break
            else:
                continue
            break
    return frozenset(out)


CODONS_ONE_SNV_FROM_STOP = codons_one_snv_from_stop()


def stop_gained_snvs_for_aa_pos(transcript: Transcript, aa_position: int) -> tuple[EnumeratedSNV, ...]:
    """Enumerate concrete stop-gained SNVs available at one amino-acid position."""
    _validate_aa_position(transcript, aa_position)
    cds = transcript.cds_sequence
    codon_start = (aa_position - 1) * 3
    codon = cds[codon_start : codon_start + 3]
    aa_ref = CODON_TABLE.get(codon)
    if aa_ref in (None, "*"):
        return ()

    snvs: list[EnumeratedSNV] = []
    for offset, base in enumerate(codon):
        for alt_base in "ACGT":
            if alt_base == base:
                continue
            alt_codon = codon[:offset] + alt_base + codon[offset + 1 :]
            if CODON_TABLE.get(alt_codon) != "*":
                continue

            cds_pos = codon_start + offset + 1
            snvs.append(
                EnumeratedSNV(
                    consequence="stop_gained",
                    cds_pos=cds_pos,
                    aa_pos=aa_position,
                    aa_ref=aa_ref,
                    aa_alt="*",
                    codon_ref=codon,
                    codon_alt=alt_codon,
                    chromosome=transcript.chromosome or "",
                    position=transcript.cds_to_genomic(cds_pos),
                    ref=transcript.genomic_ref_for_cds_base(base),
                    alt=transcript.genomic_ref_for_cds_base(alt_base),
                    hgvs_c=_format_hgvs_c(transcript, cds_pos, base, alt_base),
                    hgvs_p=_format_hgvs_p(
                        transcript,
                        aa_position,
                        aa_ref,
                        "*",
                        "stop_gained",
                    ),
                )
            )
    return tuple(snvs)


def map_frameshift_position(
    transcript: Transcript,
    aa_position: int,
    *,
    max_steps: int = DEFAULT_MAX_STEPS,
) -> tuple[FrameshiftNonsenseProxy, ...]:
    """Map one frameshift amino-acid position to nearest stop-gained SNV proxies.

    Search walks symmetrically from the frameshift amino-acid position. It stops
    at the first distance with at least one stop-gained SNV proxy. Position zero
    distance is reported with direction ``self``; otherwise directions are
    ``left`` and ``right`` in protein-coordinate space.
    """
    if max_steps < 0:
        raise ValueError("max_steps must be >= 0")
    _validate_aa_position(transcript, aa_position)

    for step in range(max_steps + 1):
        candidates: list[FrameshiftNonsenseProxy] = []
        if step == 0:
            search_positions = (("self", aa_position),)
        else:
            search_positions = (
                ("left", aa_position - step),
                ("right", aa_position + step),
            )

        for direction, proxy_aa_pos in search_positions:
            if not _aa_position_in_range(transcript, proxy_aa_pos):
                continue
            for snv in stop_gained_snvs_for_aa_pos(transcript, proxy_aa_pos):
                candidates.append(
                    FrameshiftNonsenseProxy(
                        frameshift_aa_pos=aa_position,
                        direction=direction,
                        n_steps_wrt_frameshift=step,
                        proxy=snv,
                    )
                )
        if candidates:
            return tuple(candidates)

    return ()


def annotate_gene(
    db,
    gene_symbol: str,
    *,
    transcript: Optional[Transcript] = None,
    max_steps: int = DEFAULT_MAX_STEPS,
    enqueue_sources: Iterable[str] = (),
) -> dict:
    """Persist mappings for every protein position in a gene.

    Stop-gained proxy SNVs are first upserted into the normalized variant tables
    so queued or existing CADD/MyVariant/etc. annotations remain attached to the
    concrete proxy variant IDs.
    """
    if transcript is None:
        transcript = fetch_canonical_transcript(gene_symbol)
    gene = (transcript.gene_symbol or gene_symbol).upper()
    summary = _ensure_stop_gained_proxy_space(db, transcript)
    db.clear_frameshift_nonsense_mappings(
        gene,
        transcript.transcript_id_versioned,
        source=MAPPING_SOURCE,
    )
    return _persist_positions(
        db,
        transcript,
        range(1, protein_codon_count(transcript) + 1),
        max_steps=max_steps,
        gene_symbol=gene,
        proxy_summary=summary,
        enqueue_sources=enqueue_sources,
    )


def annotate_positions(
    db,
    transcript: Transcript,
    aa_positions: Iterable[int],
    *,
    max_steps: int = DEFAULT_MAX_STEPS,
    enqueue_sources: Iterable[str] = (),
    gene_symbol: Optional[str] = None,
) -> dict:
    """Persist mappings for selected frameshift amino-acid positions."""
    positions = list(dict.fromkeys(int(pos) for pos in aa_positions))
    if not positions:
        return {
            "positions": 0,
            "mapped_positions": 0,
            "unmapped_positions": 0,
            "mappings": 0,
            "nonsense_variants": 0,
            "jobs_queued": 0,
            "by_direction": {},
        }

    gene = (gene_symbol or transcript.gene_symbol or "").upper()
    if not gene:
        raise ValueError("gene_symbol is required when transcript.gene_symbol is missing")

    summary = _ensure_stop_gained_proxy_space(db, transcript)
    db.clear_frameshift_nonsense_mappings(
        gene,
        transcript.transcript_id_versioned,
        source=MAPPING_SOURCE,
        frameshift_aa_positions=positions,
    )
    return _persist_positions(
        db,
        transcript,
        positions,
        max_steps=max_steps,
        gene_symbol=gene,
        proxy_summary=summary,
        enqueue_sources=enqueue_sources,
    )


def frameshift_variant_processing(coding_sequence: str, aa_position: int, n_steps: int) -> dict:
    """Compatibility wrapper for the original notation_conversion helper.

    It returns the nearest left/right codon that can become a stop codon by one
    SNV using only a coding sequence. New code should prefer
    :func:`map_frameshift_position`, which returns concrete genomic SNV proxies.
    """
    try:
        if n_steps < 0:
            raise ValueError("n_steps must be >= 0")
        cds = coding_sequence.upper()
        n_codons = len(cds) // 3
        if aa_position < 1 or aa_position > n_codons:
            raise ValueError("aa_position out of range")

        for step in range(n_steps + 1):
            found: dict[str, dict] = {}
            if step == 0:
                positions = (("left", aa_position), ("right", aa_position))
            else:
                positions = (("left", aa_position - step), ("right", aa_position + step))

            for side, pos in positions:
                if pos < 1 or pos > n_codons:
                    continue
                codon = cds[(pos - 1) * 3 : pos * 3]
                if codon not in CODONS_ONE_SNV_FROM_STOP:
                    continue
                aa = CODON_TABLE.get(codon)
                found[side] = {
                    "ref_aa_position": pos,
                    "ref_aa": aa,
                    "n_steps_wrt_ref_aa": step,
                    "near_stop_codon": codon,
                    "lof_variant": f"{aa}{pos}X",
                }
            if found:
                return found
    except Exception:
        pass

    return {
        "left": _empty_legacy_result(aa_position),
        "right": _empty_legacy_result(aa_position),
    }


def _persist_positions(
    db,
    transcript: Transcript,
    aa_positions: Iterable[int],
    *,
    max_steps: int,
    gene_symbol: str,
    proxy_summary: dict,
    enqueue_sources: Iterable[str] = (),
) -> dict:
    mapped_positions: set[int] = set()
    by_direction: dict[str, int] = {}
    n_mappings = 0
    n_positions = 0
    n_jobs = 0
    sources = tuple(dict.fromkeys(s for s in enqueue_sources if s))

    for aa_pos in aa_positions:
        _validate_aa_position(transcript, aa_pos)
        n_positions += 1
        proxies = map_frameshift_position(transcript, aa_pos, max_steps=max_steps)
        if proxies:
            mapped_positions.add(aa_pos)
        for proxy in proxies:
            variant_id = _variant_id_for_snv(db, proxy.proxy)
            db.upsert_frameshift_nonsense_mapping(
                gene_symbol=gene_symbol,
                transcript_id=transcript.transcript_id_versioned,
                frameshift_aa_pos=aa_pos,
                direction=proxy.direction,
                n_steps_wrt_frameshift=proxy.n_steps_wrt_frameshift,
                max_search_steps=max_steps,
                proxy_variant_id=variant_id,
                proxy_aa_pos=proxy.proxy.aa_pos,
                proxy_aa_ref=proxy.proxy.aa_ref,
                proxy_codon_ref=proxy.proxy.codon_ref,
                proxy_codon_alt=proxy.proxy.codon_alt,
                proxy_hgvs_p=proxy.proxy.hgvs_p,
                source=MAPPING_SOURCE,
            )
            n_mappings += 1
            by_direction[proxy.direction] = by_direction.get(proxy.direction, 0) + 1
            for source in sources:
                db.enqueue_job(variant_id, source=source)
                n_jobs += 1

    return {
        "positions": n_positions,
        "mapped_positions": len(mapped_positions),
        "unmapped_positions": n_positions - len(mapped_positions),
        "mappings": n_mappings,
        "nonsense_variants": proxy_summary.get("variants", 0),
        "jobs_queued": proxy_summary.get("jobs_queued", 0) + n_jobs,
        "by_direction": by_direction,
    }


def _ensure_stop_gained_proxy_space(
    db,
    transcript: Transcript,
) -> dict:
    return populate_for_transcript(
        transcript,
        db,
        types=("nonsense",),
        sources_by_consequence={"stop_gained": ()},
        enqueue=False,
    )


def _variant_id_for_snv(db, snv: EnumeratedSNV) -> int:
    cur = db.conn.execute(
        """
        SELECT id
        FROM variants
        WHERE chromosome = ? AND position = ? AND ref = ? AND alt = ?
        """,
        [snv.chromosome, snv.position, snv.ref, snv.alt],
    )
    row = cur.fetchone()
    if row is not None:
        return int(row["id"])
    return db.upsert_variant(
        chromosome=snv.chromosome,
        position=snv.position,
        ref=snv.ref,
        alt=snv.alt,
        variant_type="SNV",
    )


def _validate_aa_position(transcript: Transcript, aa_position: int) -> None:
    if not _aa_position_in_range(transcript, aa_position):
        n_codons = protein_codon_count(transcript)
        raise ValueError(f"aa_position {aa_position} out of range [1, {n_codons}]")


def _aa_position_in_range(transcript: Transcript, aa_position: int) -> bool:
    return 1 <= aa_position <= protein_codon_count(transcript)


def _empty_legacy_result(aa_position: int) -> dict:
    return {
        "ref_aa_position": aa_position,
        "ref_aa": None,
        "n_steps_wrt_ref_aa": None,
        "near_stop_codon": None,
        "lof_variant": None,
    }
