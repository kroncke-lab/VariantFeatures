"""Variant identity resolution.

Given any common form of a variant identifier (HGVS, rsID, ClinVar ID,
ClinGen Allele Registry CA-ID, gnomAD ID, VCF tuple), resolve to a canonical
GRCh38 representation plus the full set of known aliases.

Backed primarily by the ClinGen Allele Registry REST API
(https://reg.clinicalgenome.org), which already aggregates aliases across
dbSNP, ClinVar, gnomAD, and MyVariant.info. Read access is unauthenticated.

Internally, GRCh38 is the canonical build; any GRCh37 input is resolved
through ClinGen AR which returns both builds.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Optional

import requests


CLINGEN_AR_BASE = "https://reg.clinicalgenome.org"
DEFAULT_TIMEOUT = 30
SOURCE_CLINGEN = "clingen_ar"

_RS_RE = re.compile(r"^rs\d+$", re.IGNORECASE)
_CA_RE = re.compile(r"^CA\d+$", re.IGNORECASE)
_VCV_RE = re.compile(r"^VCV0*(\d+)$", re.IGNORECASE)
_RCV_RE = re.compile(r"^RCV0*(\d+)$", re.IGNORECASE)
_NUMERIC_RE = re.compile(r"^\d+$")
_GNOMAD_RE = re.compile(r"^(?:chr)?[\dXYMT]+[-:]\d+[-:][ACGTN]+[-:][ACGTN]+$", re.IGNORECASE)
_HGVS_G_RE = re.compile(r":g\.")
_HGVS_C_RE = re.compile(r":c\.")
_HGVS_P_RE = re.compile(r":p\.")


class IdentityError(Exception):
    """Raised when a variant cannot be resolved."""


@dataclass
class ResolvedVariant:
    """Canonical variant identity + all known aliases.

    Coordinates are GRCh38, 1-based, left-aligned (VCF-normalized).
    """

    ca_id: Optional[str] = None
    chromosome: Optional[str] = None
    position: Optional[int] = None
    ref: Optional[str] = None
    alt: Optional[str] = None
    variant_type: Optional[str] = None
    hgvs_g: Optional[str] = None
    aliases: list[dict] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            "ca_id": self.ca_id,
            "chromosome": self.chromosome,
            "position": self.position,
            "ref": self.ref,
            "alt": self.alt,
            "variant_type": self.variant_type,
            "hgvs_g": self.hgvs_g,
            "aliases": list(self.aliases),
        }


def detect_kind(value: str) -> str:
    """Best-effort classification of an input identifier."""
    v = value.strip()
    if _CA_RE.match(v):
        return "ca_id"
    if _RS_RE.match(v):
        return "rsid"
    if _VCV_RE.match(v):
        return "clinvar_vcv"
    if _RCV_RE.match(v):
        return "clinvar_rcv"
    if _HGVS_G_RE.search(v):
        return "hgvs_g"
    if _HGVS_C_RE.search(v):
        return "hgvs_c"
    if _HGVS_P_RE.search(v):
        return "hgvs_p"
    if _GNOMAD_RE.match(v):
        return "gnomad_id"
    if _NUMERIC_RE.match(v):
        return "clinvar_variation"
    raise IdentityError(f"Cannot classify identifier: {value!r}")


def resolve(value: str, kind: Optional[str] = None, *, timeout: int = DEFAULT_TIMEOUT) -> ResolvedVariant:
    """Resolve any supported identifier to canonical form + aliases.

    p.HGVS is not supported directly by ClinGen AR; pass the corresponding
    c.HGVS or g.HGVS instead.
    """
    if kind is None:
        kind = detect_kind(value)
    payload = _fetch_clingen(value, kind, timeout=timeout)
    return _parse_clingen(payload)


def _fetch_clingen(value: str, kind: str, *, timeout: int) -> dict:
    """Route an identifier to the right ClinGen AR endpoint.

    The registry exposes two endpoints:
      - /allele              -- single-result lookups: by ?hgvs= or via path /allele/{CA_ID}
      - /alleles             -- index lookups (rsID, MyVariantInfo, gnomAD, ClinVar AlleleId);
                                returns a JSON list.
    """
    if kind == "ca_id":
        url = f"{CLINGEN_AR_BASE}/allele/{value.upper()}"
        params: dict = {}
    elif kind in {"hgvs_g", "hgvs_c"}:
        url = f"{CLINGEN_AR_BASE}/allele"
        params = {"hgvs": value}
    elif kind == "rsid":
        rs_num = value[2:] if value.lower().startswith("rs") else value
        url = f"{CLINGEN_AR_BASE}/alleles"
        params = {"dbSNP.rs": rs_num}
    elif kind == "gnomad_id":
        url = f"{CLINGEN_AR_BASE}/alleles"
        params = {"MyVariantInfo_hg38.id": value} if value.lower().startswith("chr") else {"gnomAD.id": value}
    elif kind in {"clinvar_vcv", "clinvar_rcv", "clinvar_variation"}:
        # ClinGen AR's /alleles endpoint indexes by ClinVar AlleleId, not variationId.
        # For VCV / RCV / variationId, the user must look up the AlleleId via ClinVar first.
        raise IdentityError(
            "ClinVar variation/RCV/VCV is not directly indexed by ClinGen Allele Registry. "
            "Pass the ClinVar Allele ID (numeric) via --kind=ca_id-style HGVS instead, "
            "or look it up via NCBI ClinVar's API and resolve the resulting HGVS."
        )
    elif kind == "hgvs_p":
        raise IdentityError(
            "Protein HGVS (p.) cannot be resolved directly by ClinGen Allele Registry. "
            "Pass the corresponding c.HGVS (e.g. NM_xxx:c.NNNT>C) instead."
        )
    else:
        raise IdentityError(f"Unsupported identifier kind: {kind!r}")

    resp = requests.get(url, params=params, timeout=timeout)
    if resp.status_code == 404:
        raise IdentityError(f"Variant not found in ClinGen Allele Registry: {value!r}")
    if resp.status_code == 400:
        # ClinGen AR returns 400 with a structured error body (e.g. reference mismatch)
        try:
            body = resp.json()
        except ValueError:
            body = {}
        msg = body.get("description") or body.get("message") or resp.text
        raise IdentityError(f"ClinGen AR rejected {value!r}: {msg}")
    resp.raise_for_status()
    data = resp.json()
    if isinstance(data, list):
        if not data:
            raise IdentityError(f"Empty result list from ClinGen AR for: {value!r}")
        data = data[0]
    if isinstance(data, dict) and data.get("errorType"):
        raise IdentityError(f"ClinGen AR error for {value!r}: {data.get('description', data)}")
    return data


def _parse_clingen(payload: dict) -> ResolvedVariant:
    rv = ResolvedVariant()

    # CA-ID is the last path segment of @id
    at_id = payload.get("@id", "")
    if at_id:
        ca = at_id.rstrip("/").rsplit("/", 1)[-1]
        if _CA_RE.match(ca):
            rv.ca_id = ca.upper()
            _add_alias(rv, "ca_id", rv.ca_id, SOURCE_CLINGEN)

    # Genomic alleles: prefer GRCh38; collect every g.HGVS as an alias
    grch38_entry = None
    for ga in payload.get("genomicAlleles", []) or []:
        for h in ga.get("hgvs", []) or []:
            _add_alias(rv, "hgvs_g", h, SOURCE_CLINGEN)
        if ga.get("referenceGenome") == "GRCh38" and grch38_entry is None:
            grch38_entry = ga

    if grch38_entry is not None:
        rv.chromosome = grch38_entry.get("chromosome")
        coords = grch38_entry.get("coordinates") or []
        if coords:
            c = coords[0]
            # ClinGen AR uses 0-based start, 1-based end (per HGVS doc).
            # For a SNV this is start=position-1, end=position. Convert to 1-based VCF position.
            start_0 = c.get("start")
            end_1 = c.get("end")
            ref = c.get("referenceAllele", "") or ""
            alt = c.get("allele", "") or ""
            if start_0 is not None and end_1 is not None:
                rv.position = start_0 + 1
                rv.ref = ref
                rv.alt = alt
                rv.variant_type = _classify_variant(ref, alt)
        # Pick the canonical NC_-style hgvs as hgvs_g if present
        for h in grch38_entry.get("hgvs", []) or []:
            if h.startswith("NC_"):
                rv.hgvs_g = h
                break
        if rv.hgvs_g is None and (grch38_entry.get("hgvs") or []):
            rv.hgvs_g = grch38_entry["hgvs"][0]

    # Transcript alleles: collect c.HGVS and the protein effect p.HGVS
    for ta in payload.get("transcriptAlleles", []) or []:
        for h in ta.get("hgvs", []) or []:
            _add_alias(rv, "hgvs_c", h, SOURCE_CLINGEN)
        protein = ta.get("proteinEffect") or {}
        p_hgvs = protein.get("hgvs")
        if p_hgvs:
            _add_alias(rv, "hgvs_p", p_hgvs, SOURCE_CLINGEN)

    # External records
    ext = payload.get("externalRecords") or {}

    for r in ext.get("dbSNP", []) or []:
        rs = r.get("rs")
        if rs:
            _add_alias(rv, "rsid", f"rs{rs}", SOURCE_CLINGEN)

    for r in ext.get("ClinVarVariations", []) or []:
        vid = r.get("variationId") or r.get("id")
        if vid:
            _add_alias(rv, "clinvar_vcv", f"VCV{int(vid):09d}", SOURCE_CLINGEN)
        for rcv in r.get("RCV", []) or []:
            _add_alias(rv, "clinvar_rcv", rcv, SOURCE_CLINGEN)

    for r in ext.get("ClinVarAlleles", []) or []:
        aid = r.get("alleleId") or r.get("id")
        if aid:
            _add_alias(rv, "clinvar_allele", str(aid), SOURCE_CLINGEN)

    for r in ext.get("gnomAD", []) or []:
        gid = r.get("id")
        if gid:
            _add_alias(rv, "gnomad_id", gid, SOURCE_CLINGEN)

    for key in ("MyVariantInfo_hg38", "MyVariantInfo_hg19"):
        for r in ext.get(key, []) or []:
            mid = r.get("id")
            if mid:
                _add_alias(rv, "myvariant_id", mid, SOURCE_CLINGEN)

    return rv


def _add_alias(rv: ResolvedVariant, alias_type: str, alias_value: str, source: str) -> None:
    if not alias_value:
        return
    entry = {"alias_type": alias_type, "alias_value": alias_value, "source": source}
    if entry not in rv.aliases:
        rv.aliases.append(entry)


def _classify_variant(ref: str, alt: str) -> str:
    if ref == "" and alt:
        return "INS"
    if alt == "" and ref:
        return "DEL"
    if len(ref) == 1 and len(alt) == 1:
        return "SNV"
    if len(ref) == len(alt):
        return "MNV"
    if len(ref) < len(alt) and alt.startswith(ref):
        return "INS"
    if len(ref) > len(alt) and ref.startswith(alt):
        return "DEL"
    return "DELINS"


def persist(rv: ResolvedVariant, db) -> Optional[int]:
    """Write a resolved variant + aliases to the database. Returns variant_id."""
    if rv.chromosome is None or rv.position is None or rv.ref is None or rv.alt is None:
        return None
    variant_id = db.upsert_variant(
        chromosome=rv.chromosome,
        position=rv.position,
        ref=rv.ref,
        alt=rv.alt,
        ca_id=rv.ca_id,
        variant_type=rv.variant_type,
        hgvs_g=rv.hgvs_g,
    )
    db.add_aliases(variant_id, rv.aliases)
    return variant_id
