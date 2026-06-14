"""Saturation enumeration of single-nucleotide variants in a CDS.

Walks every codon, generates all 3 possible substitutions at each position,
and classifies each as missense_variant / stop_gained / synonymous_variant /
stop_lost / start_lost using the standard genetic code.

For each kept variant we:
- compute the GRCh38 (chrom, pos, ref, alt) tuple, accounting for strand;
- compute hgvs_c (transcript orientation) and hgvs_p (3-letter, HGVS Ter);
- upsert the canonical row in `variants`;
- upsert a `variant_consequences` row tagged source='enumerated';
- enqueue annotation jobs for each requested external source.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Iterator, Optional

from .transcripts import Transcript


# Standard genetic code (DNA bases).
CODON_TABLE: dict[str, str] = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

AA_3LETTER: dict[str, str] = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys",
    "E": "Glu", "Q": "Gln", "G": "Gly", "H": "His", "I": "Ile",
    "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro",
    "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val",
    "*": "Ter",
}

# User-friendly type names map to SO terms.
TYPE_TO_CONSEQUENCES: dict[str, set[str]] = {
    "missense": {"missense_variant"},
    "nonsense": {"stop_gained"},
    "synonymous": {"synonymous_variant"},
    "stop_lost": {"stop_lost"},
    "start_lost": {"start_lost"},
}

# Default annotation sources to queue per consequence type.
# `myvariant` is a meta-source that brings back CADD + dbNSFP (~25 predictors)
# + gnomAD + ClinVar in a single REST call; it's queued for every variant.
# `clingen_ar` is queued for full alias resolution (every name a variant has).
# Per-source jobs (alphamissense, revel, ...) are also queued for missense so
# higher-quality dedicated fetchers can refine values later.
DEFAULT_SOURCES_BY_CONSEQUENCE: dict[str, tuple[str, ...]] = {
    "missense_variant": ("clingen_ar", "myvariant", "gnomad", "alphamissense", "revel", "alphafold", "annovar", "vep"),
    "stop_gained": ("clingen_ar", "myvariant", "gnomad", "alphafold", "annovar", "vep"),
    "stop_lost": ("clingen_ar", "myvariant", "gnomad", "alphafold", "annovar", "vep"),
    "start_lost": ("clingen_ar", "myvariant", "gnomad", "alphafold", "annovar", "vep"),
    "synonymous_variant": ("clingen_ar", "myvariant", "gnomad", "annovar", "vep"),
}


@dataclass
class EnumeratedSNV:
    """One single-nucleotide variant produced by saturation mutagenesis."""

    consequence: str
    cds_pos: int                        # 1-based CDS position
    aa_pos: int                         # 1-based amino acid (codon number)
    aa_ref: str                         # 1-letter (or '*' for stop)
    aa_alt: str
    codon_ref: str                      # 3-letter codon (transcript orientation)
    codon_alt: str
    chromosome: str
    position: int                       # 1-based GRCh38
    ref: str                            # + strand reference base
    alt: str
    hgvs_c: str
    hgvs_p: Optional[str]

    def to_dict(self) -> dict:
        return self.__dict__.copy()


# ---------------------------------------------------------------------------
# Pure enumeration
# ---------------------------------------------------------------------------

def normalize_types(types: Iterable[str]) -> set[str]:
    """Translate user-friendly type names ('missense', 'nonsense') into SO consequences."""
    if isinstance(types, str):
        types = [t.strip() for t in types.split(",") if t.strip()]
    out: set[str] = set()
    for t in types:
        key = t.strip().lower()
        if key in TYPE_TO_CONSEQUENCES:
            out |= TYPE_TO_CONSEQUENCES[key]
        elif key in {c for cs in TYPE_TO_CONSEQUENCES.values() for c in cs}:
            out.add(key)
        else:
            raise ValueError(f"Unknown variant type: {t!r}")
    return out


def enumerate_snvs(
    transcript: Transcript,
    types: Iterable[str] = ("missense", "nonsense"),
) -> Iterator[EnumeratedSNV]:
    """Yield every kept SNV produced by walking every codon * every alt base."""
    keep = normalize_types(types)
    cds = transcript.cds_sequence
    if not cds:
        return
    n_codons = len(cds) // 3

    for codon_idx in range(n_codons):
        codon_start_cds = codon_idx * 3 + 1     # 1-based CDS position of codon's 1st base
        codon = cds[codon_idx * 3 : codon_idx * 3 + 3]
        if len(codon) != 3:
            continue
        aa_ref = CODON_TABLE.get(codon)
        if aa_ref is None:
            continue
        for j in range(3):
            base = codon[j]
            for alt_base in "ACGT":
                if alt_base == base:
                    continue
                new_codon = codon[:j] + alt_base + codon[j + 1 :]
                aa_alt = CODON_TABLE.get(new_codon)
                if aa_alt is None:
                    continue
                consequence = _classify_codon_change(codon_idx, aa_ref, aa_alt)
                if consequence not in keep:
                    continue

                cds_pos = codon_start_cds + j
                genomic_pos = transcript.cds_to_genomic(cds_pos)
                vcf_ref = transcript.genomic_ref_for_cds_base(base)
                vcf_alt = transcript.genomic_ref_for_cds_base(alt_base)

                yield EnumeratedSNV(
                    consequence=consequence,
                    cds_pos=cds_pos,
                    aa_pos=codon_idx + 1,
                    aa_ref=aa_ref,
                    aa_alt=aa_alt,
                    codon_ref=codon,
                    codon_alt=new_codon,
                    chromosome=transcript.chromosome or "",
                    position=genomic_pos,
                    ref=vcf_ref,
                    alt=vcf_alt,
                    hgvs_c=_format_hgvs_c(transcript, cds_pos, base, alt_base),
                    hgvs_p=_format_hgvs_p(transcript, codon_idx + 1, aa_ref, aa_alt, consequence),
                )


def _classify_codon_change(codon_idx: int, aa_ref: str, aa_alt: str) -> str:
    if codon_idx == 0 and aa_ref == "M" and aa_alt != "M":
        return "start_lost"
    if aa_ref == "*" and aa_alt == "*":
        return "stop_retained_variant"
    if aa_ref == aa_alt:
        return "synonymous_variant"
    if aa_alt == "*":
        return "stop_gained"
    if aa_ref == "*":
        return "stop_lost"
    return "missense_variant"


def _format_hgvs_c(transcript: Transcript, cds_pos: int, ref: str, alt: str) -> str:
    return f"{transcript.transcript_id_versioned}:c.{cds_pos}{ref}>{alt}"


def _format_hgvs_p(transcript: Transcript, aa_pos: int, aa_ref: str, aa_alt: str, consequence: str) -> Optional[str]:
    prefix = transcript.protein_id or transcript.transcript_id_versioned
    aa_ref3 = AA_3LETTER.get(aa_ref, aa_ref)
    aa_alt3 = AA_3LETTER.get(aa_alt, aa_alt)
    if consequence == "synonymous_variant":
        return f"{prefix}:p.{aa_ref3}{aa_pos}="
    if consequence == "start_lost":
        return f"{prefix}:p.Met1?"
    if consequence == "stop_lost":
        return f"{prefix}:p.Ter{aa_pos}{aa_alt3}ext*?"
    return f"{prefix}:p.{aa_ref3}{aa_pos}{aa_alt3}"


# ---------------------------------------------------------------------------
# Persistence
# ---------------------------------------------------------------------------

def populate_for_transcript(
    transcript: Transcript,
    db,
    *,
    types: Iterable[str] = ("missense", "nonsense"),
    sources_by_consequence: Optional[dict[str, tuple[str, ...]]] = None,
    enqueue: bool = True,
    limit: Optional[int] = None,
) -> dict:
    """Enumerate, persist, and (optionally) queue annotation jobs.

    Returns a summary dict: {'variants': int, 'consequences': int, 'jobs_queued': int,
    'by_consequence': {consequence: count}}.
    """
    sources_map = sources_by_consequence or DEFAULT_SOURCES_BY_CONSEQUENCE

    counts_by_consequence: dict[str, int] = {}
    n_variants = 0
    n_consequences = 0
    n_jobs = 0

    for snv in enumerate_snvs(transcript, types=types):
        if limit is not None and n_variants >= limit:
            break
        counts_by_consequence[snv.consequence] = counts_by_consequence.get(snv.consequence, 0) + 1

        variant_id = db.upsert_variant(
            chromosome=snv.chromosome,
            position=snv.position,
            ref=snv.ref,
            alt=snv.alt,
            variant_type="SNV",
        )
        n_variants += 1

        # Add hgvs aliases for the canonical c. and p. notations
        aliases = [{"alias_type": "hgvs_c", "alias_value": snv.hgvs_c, "source": "enumerated"}]
        if snv.hgvs_p:
            aliases.append({"alias_type": "hgvs_p", "alias_value": snv.hgvs_p, "source": "enumerated"})
        db.add_aliases(variant_id, aliases)

        db.upsert_consequence(
            variant_id=variant_id,
            transcript_id=transcript.transcript_id_versioned,
            source="enumerated",
            gene_symbol=transcript.gene_symbol,
            gene_ensembl=transcript.gene_ensembl,
            consequence=snv.consequence,
            hgvs_c=snv.hgvs_c,
            hgvs_p=snv.hgvs_p,
            aa_pos=snv.aa_pos,
            aa_ref=snv.aa_ref,
            aa_alt=snv.aa_alt,
            codon_pos=snv.aa_pos,
            codon_ref=snv.codon_ref,
            codon_alt=snv.codon_alt,
            is_canonical=int(transcript.is_canonical),
            is_mane_select=int(transcript.is_mane_select),
            is_mane_plus_clinical=int(transcript.is_mane_plus_clinical),
        )
        n_consequences += 1

        if enqueue:
            for source in sources_map.get(snv.consequence, ()):  # default to none if unknown
                n_jobs += db.enqueue_job(variant_id, source=source)

    return {
        "variants": n_variants,
        "consequences": n_consequences,
        "jobs_queued": n_jobs,
        "by_consequence": counts_by_consequence,
    }
