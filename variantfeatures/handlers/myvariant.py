"""MyVariant.info annotation handler.

A single REST call to https://myvariant.info/v1/variant/<id> returns a deeply
aggregated JSON document covering CADD, dbNSFP (~25 in-silico predictors),
gnomAD exomes/genomes, ClinVar, dbSNP, ClinGen Allele Registry, and more.

This handler parses the response into the generic annotation tables. Each
sub-source is queried with a single round trip — much cheaper than calling
each underlying API individually.

Public docs: https://docs.myvariant.info/
"""

from __future__ import annotations

import json
from typing import Any, Iterable, Optional

import requests

from .clingen_ar import GRCH38_NC_ACCESSIONS, normalize_chromosome  # reuse the chr -> NC table


SOURCE = "myvariant"
MYVARIANT_BASE = "https://myvariant.info/v1/variant"
DEFAULT_RATE_LIMIT_SEC = 0.10
DEFAULT_TIMEOUT = 30


# Population codes we surface (gnomAD's typical subdivisions). Sex and sub-pop breakdowns are skipped here.
GNOMAD_POPS = ["all", "afr", "amr", "asj", "eas", "fin", "nfe", "sas", "oth", "popmax"]

# dbNSFP predictors with a uniform sub-shape: {score, rankscore, pred}
DBNSFP_FLAT_PREDICTORS: list[tuple[str, str]] = [
    # (json key in dbnsfp, our predictor name)
    ("alphamissense", "alphamissense"),
    ("revel", "revel"),
    ("metasvm", "metasvm"),
    ("metalr", "metalr"),
    ("m-cap", "mcap"),
    ("primateai", "primateai"),
    ("clinpred", "clinpred"),
    ("eve", "eve"),
    ("esm1b", "esm1b"),
    ("dann", "dann"),
    ("deogen2", "deogen2"),
    ("eigen", "eigen"),
    ("eigen-pc", "eigen_pc"),
    ("fathmm-mkl", "fathmm_mkl"),
    ("fathmm-xf", "fathmm_xf"),
    ("genocanyon", "genocanyon"),
    ("gmvp", "gmvp"),
    ("metarnn", "metarnn"),
    ("mutpred", "mutpred"),
    ("mutationassessor", "mutation_assessor"),
    ("mutationtaster", "mutation_taster"),
    ("provean", "provean"),
    ("sift", "sift"),
    ("vest4", "vest4"),
]


class HandlerError(Exception):
    pass


def _myvariant_id(chromosome: str, position: int, ref: str, alt: str) -> str:
    """Build a MyVariant.info HGVS-style ID. SNV form: chr<chr>:g.<pos><ref>><alt>."""
    chrom = normalize_chromosome(chromosome)
    if chrom not in GRCH38_NC_ACCESSIONS:
        raise HandlerError(f"Unrecognized chromosome {chromosome!r}")
    return f"chr{chrom}:g.{position}{ref}>{alt}"


def fetch_myvariant(myvariant_id: str, *, assembly: str = "hg38", timeout: int = DEFAULT_TIMEOUT) -> Optional[dict]:
    """GET MyVariant.info for one variant. Returns None if not found (404)."""
    url = f"{MYVARIANT_BASE}/{myvariant_id}"
    resp = requests.get(url, params={"assembly": assembly}, timeout=timeout)
    if resp.status_code == 404:
        return None
    resp.raise_for_status()
    data = resp.json()
    if isinstance(data, dict) and data.get("notfound"):
        return None
    return data


# ---------------------------------------------------------------------------
# Parsers
# ---------------------------------------------------------------------------

def _coerce_numeric(v) -> Optional[float]:
    """Force numeric-or-None. Lists collapse to max of their numeric entries."""
    if v is None:
        return None
    if isinstance(v, bool):
        # bools are ints in Python; explicitly drop them for score fields.
        return None
    if isinstance(v, (int, float)):
        return float(v)
    if isinstance(v, list):
        nums = [x for x in v if isinstance(x, (int, float)) and not isinstance(x, bool)]
        return float(max(nums)) if nums else None
    return None  # strings, dicts, etc. are not valid scores


def _coerce_category(v) -> Optional[str]:
    """Force string-or-None. Lists pick first; dicts get JSON-encoded; everything else stringified."""
    if v is None:
        return None
    if isinstance(v, list):
        v = v[0] if v else None
        if v is None:
            return None
    if isinstance(v, dict):
        return json.dumps(v, separators=(",", ":"), sort_keys=True)
    return str(v)


def _extract_score_triplet(node: Any) -> tuple[Optional[float], Optional[float], Optional[str]]:
    """For dbNSFP-style {score, rankscore, pred}, pull (score, rank_score, category)."""
    if not isinstance(node, dict):
        return None, None, None
    return (
        _coerce_numeric(node.get("score")),
        _coerce_numeric(node.get("rankscore")),
        _coerce_category(node.get("pred")),
    )


def parse_pathogenicity(payload: dict) -> Iterable[dict]:
    """Yield dicts of {predictor, score, rank_score, category} for every recognized predictor."""
    dbnsfp = payload.get("dbnsfp") or {}

    # CADD: prefer dbnsfp.cadd (it has rankscore alongside phred and raw)
    cadd = dbnsfp.get("cadd") or {}
    if cadd:
        phred = _coerce_numeric(cadd.get("phred"))
        if phred is not None:
            yield {"predictor": "cadd_phred", "score": phred, "rank_score": None, "category": None}
        raw = _coerce_numeric(cadd.get("raw_score"))
        if raw is not None:
            yield {
                "predictor": "cadd_raw",
                "score": raw,
                "rank_score": _coerce_numeric(cadd.get("raw_rankscore")),
                "category": None,
            }

    # Top-level CADD as fallback (sometimes present even when dbnsfp.cadd absent)
    top_cadd = payload.get("cadd") or {}
    if top_cadd and not cadd:
        phred = _coerce_numeric(top_cadd.get("phred"))
        if phred is not None:
            yield {"predictor": "cadd_phred", "score": phred, "rank_score": None, "category": None}
        raw = _coerce_numeric(top_cadd.get("raw_score"))
        if raw is not None:
            yield {"predictor": "cadd_raw", "score": raw, "rank_score": None, "category": None}

    # Flat dbNSFP predictors
    for key, name in DBNSFP_FLAT_PREDICTORS:
        sub = dbnsfp.get(key)
        score, rank_score, category = _extract_score_triplet(sub)
        if score is not None or rank_score is not None or category is not None:
            yield {"predictor": name, "score": score, "rank_score": rank_score, "category": category}

    # PolyPhen-2 has hdiv and hvar variants
    pp2 = dbnsfp.get("polyphen2") or {}
    for sub_name, our_name in (("hdiv", "polyphen2_hdiv"), ("hvar", "polyphen2_hvar")):
        s, r, c = _extract_score_triplet(pp2.get(sub_name))
        if s is not None or r is not None or c is not None:
            yield {"predictor": our_name, "score": s, "rank_score": r, "category": c}

    # BayesDel has add_af and no_af variants
    bd = dbnsfp.get("bayesdel") or {}
    for sub_name, our_name in (("add_af", "bayesdel_add_af"), ("no_af", "bayesdel_no_af")):
        s, r, c = _extract_score_triplet(bd.get(sub_name))
        if s is not None or r is not None or c is not None:
            yield {"predictor": our_name, "score": s, "rank_score": r, "category": c}


def parse_conservation(payload: dict) -> Iterable[dict]:
    """Yield dicts of {metric, score, rank_score} for cross-species conservation tracks."""
    dbnsfp = payload.get("dbnsfp") or {}

    phylop = (dbnsfp.get("phylop") or {})
    if isinstance(phylop, dict) and "phylop" in phylop:
        # Some MyVariant releases nest phylop one extra level
        phylop = phylop.get("phylop") or {}
    for sub_name, metric in (
        ("100way_vertebrate", "phylop100way_vertebrate"),
        ("470way_mammalian", "phylop470way_mammalian"),
        ("17way_primate", "phylop17way_primate"),
    ):
        node = phylop.get(sub_name) if isinstance(phylop, dict) else None
        if isinstance(node, dict):
            yield {"metric": metric, "score": _coerce_numeric(node.get("score")), "rank_score": _coerce_numeric(node.get("rankscore"))}

    phastcons = (dbnsfp.get("phastcons") or {})
    if isinstance(phastcons, dict) and "phastcons" in phastcons:
        phastcons = phastcons.get("phastcons") or {}
    for sub_name, metric in (
        ("100way_vertebrate", "phastcons100way_vertebrate"),
        ("470way_mammalian", "phastcons470way_mammalian"),
        ("17way_primate", "phastcons17way_primate"),
    ):
        node = phastcons.get(sub_name) if isinstance(phastcons, dict) else None
        if isinstance(node, dict):
            yield {"metric": metric, "score": _coerce_numeric(node.get("score")), "rank_score": _coerce_numeric(node.get("rankscore"))}

    gerp = dbnsfp.get("gerp++") or {}
    if "rs" in gerp or "nr" in gerp:
        yield {"metric": "gerp_pp_rs", "score": _coerce_numeric(gerp.get("rs")), "rank_score": _coerce_numeric(gerp.get("rs_rankscore"))}
        nr = _coerce_numeric(gerp.get("nr"))
        if nr is not None:
            yield {"metric": "gerp_pp_nr", "score": nr, "rank_score": None}

    siphy = dbnsfp.get("siphy_29way") or {}
    if "logodds" in siphy:
        yield {
            "metric": "siphy_29way_logodds",
            "score": _coerce_numeric(siphy.get("logodds")),
            "rank_score": _coerce_numeric(siphy.get("logodds_rankscore")),
        }


def parse_population(payload: dict) -> Iterable[dict]:
    """Yield dicts of {dataset, pop, af, ac, an, n_homozygotes, filter_status}."""
    for top_key, dataset in (("gnomad_exome", "gnomad_exome"), ("gnomad_genome", "gnomad_genome")):
        node = payload.get(top_key) or {}
        af_dict = node.get("af") or {}
        ac_dict = node.get("ac") or {}
        an_node = node.get("an")
        an_dict = an_node if isinstance(an_node, dict) else {}
        hom_dict = node.get("hom") or {}
        filter_status = node.get("filter")
        if isinstance(filter_status, list):
            filter_status = ",".join(str(f) for f in filter_status)

        for pop in GNOMAD_POPS:
            if pop == "all":
                af = af_dict.get("af")
                ac = ac_dict.get("ac")
                an = an_dict.get("an") if isinstance(an_node, dict) else an_node
                hom = hom_dict.get("hom")
            else:
                af = af_dict.get(f"af_{pop}")
                ac = ac_dict.get(f"ac_{pop}")
                an = an_dict.get(f"an_{pop}") if isinstance(an_node, dict) else None
                hom = hom_dict.get(f"hom_{pop}")
            # Skip empty rows (no data for this pop in this dataset)
            if af is None and ac is None and an is None and hom is None:
                continue
            yield {
                "dataset": dataset,
                "pop": pop,
                "af": af,
                "ac": ac,
                "an": an,
                "n_homozygotes": hom,
                "filter_status": filter_status,
            }


def parse_clinical(payload: dict) -> Iterable[dict]:
    """Yield ClinVar RCV-level dicts: {source, record_id, classification, review_status, last_evaluated, conditions}."""
    clinvar = payload.get("clinvar") or {}
    if not clinvar:
        return
    rcvs = clinvar.get("rcv") or []
    if isinstance(rcvs, dict):
        rcvs = [rcvs]
    if not rcvs and clinvar.get("clinical_significance"):
        # Fallback: at least record the variant-level summary if no RCV list
        yield {
            "source": "clinvar",
            "record_id": str(clinvar.get("variant_id") or clinvar.get("allele_id") or ""),
            "classification": clinvar.get("clinical_significance"),
            "review_status": clinvar.get("review_status"),
            "stars": clinvar.get("gold_stars"),
            "last_evaluated": clinvar.get("last_evaluated"),
            "conditions": None,
        }
        return
    for r in rcvs:
        if not isinstance(r, dict):
            continue
        cond = r.get("conditions")
        if isinstance(cond, dict):
            cond_str = cond.get("name")
        elif isinstance(cond, list):
            cond_str = ",".join((c.get("name") or "") for c in cond if isinstance(c, dict))
        else:
            cond_str = None
        yield {
            "source": "clinvar",
            "record_id": str(r.get("accession") or ""),
            "classification": _coerce_category(r.get("clinical_significance")),
            "review_status": _coerce_category(r.get("review_status")),
            "stars": r.get("gold_stars") if isinstance(r.get("gold_stars"), int) else None,
            "last_evaluated": _coerce_category(r.get("last_evaluated")),
            "conditions": _coerce_category(cond_str),
        }


def parse_aliases(payload: dict) -> Iterable[dict]:
    """Yield {alias_type, alias_value} pairs for the canonical alias table."""
    clingen = payload.get("clingen") or {}
    caid = clingen.get("caid")
    if caid:
        yield {"alias_type": "ca_id", "alias_value": caid}

    dbsnp = payload.get("dbsnp") or {}
    rsid = dbsnp.get("rsid")
    if rsid:
        yield {"alias_type": "rsid", "alias_value": str(rsid)}

    mv_id = payload.get("_id")
    if mv_id:
        yield {"alias_type": "myvariant_id", "alias_value": mv_id}

    clinvar = payload.get("clinvar") or {}
    vid = clinvar.get("variant_id")
    if vid:
        yield {"alias_type": "clinvar_vcv", "alias_value": f"VCV{int(vid):09d}"}
    aid = clinvar.get("allele_id")
    if aid:
        yield {"alias_type": "clinvar_allele", "alias_value": str(aid)}


# ---------------------------------------------------------------------------
# Persistence
# ---------------------------------------------------------------------------

def persist(db, variant_id: int, payload: dict) -> dict:
    """Write all parseable annotation rows from a MyVariant payload. Returns counts."""
    counts = {"pathogenicity": 0, "conservation": 0, "population": 0, "clinical": 0, "aliases": 0}

    for row in parse_pathogenicity(payload):
        db.upsert_pathogenicity(variant_id, source=SOURCE, **row)
        counts["pathogenicity"] += 1

    for row in parse_conservation(payload):
        db.upsert_conservation(variant_id, source=SOURCE, **row)
        counts["conservation"] += 1

    for row in parse_population(payload):
        db.upsert_population(variant_id, source=SOURCE, **row)
        counts["population"] += 1

    for row in parse_clinical(payload):
        db.upsert_clinical(variant_id, **row)
        counts["clinical"] += 1

    aliases = list(parse_aliases(payload))
    if aliases:
        counts["aliases"] = db.add_aliases(variant_id, aliases, source=SOURCE)

    # Update canonical CA-ID on the variants row if we learned one (and it's not already set)
    clingen = payload.get("clingen") or {}
    if clingen.get("caid"):
        cur = db.conn.execute("SELECT chromosome, position, ref, alt, ca_id FROM variants WHERE id = ?", [variant_id])
        row = cur.fetchone()
        if row and not row["ca_id"]:
            db.upsert_variant(
                chromosome=row["chromosome"], position=row["position"], ref=row["ref"], alt=row["alt"],
                ca_id=clingen["caid"],
            )

    return counts


def handle(db, variant_id: int, payload: Optional[str] = None, *, timeout: int = DEFAULT_TIMEOUT) -> None:
    """Worker entry point: fetch MyVariant.info for the variant and persist all annotations."""
    cur = db.conn.execute(
        "SELECT chromosome, position, ref, alt FROM variants WHERE id = ?",
        [variant_id],
    )
    row = cur.fetchone()
    if row is None:
        raise HandlerError(f"variant_id {variant_id} not found in variants table")
    if len(row["ref"]) != 1 or len(row["alt"]) != 1:
        raise HandlerError(
            f"myvariant handler currently only handles SNVs (got ref={row['ref']!r}, alt={row['alt']!r})"
        )

    mv_id = _myvariant_id(row["chromosome"], row["position"], row["ref"], row["alt"])
    data = fetch_myvariant(mv_id, assembly="hg38", timeout=timeout)
    if data is None:
        # Not found in MyVariant.info — record nothing but don't fail. The job is "done".
        return
    persist(db, variant_id, data)
