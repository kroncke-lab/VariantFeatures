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

from datetime import date
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
    ("mvp", "mvp"),
    ("mpc", "mpc"),
    ("sift4g", "sift4g"),
    ("lrt", "lrt"),
    ("list-s2", "list_s2"),
    ("bstatistic", "bstatistic"),
    ("dann", "dann"),
    ("deogen2", "deogen2"),
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


def fetch_myvariant_batch(
    myvariant_ids: list[str],
    *,
    assembly: str = "hg38",
    timeout: int = DEFAULT_TIMEOUT,
) -> list[dict]:
    """POST a batch of MyVariant IDs to the batch annotation endpoint."""
    if not myvariant_ids:
        return []
    resp = requests.post(
        MYVARIANT_BASE,
        data={"ids": ",".join(myvariant_ids), "assembly": assembly},
        timeout=timeout,
    )
    resp.raise_for_status()
    data = resp.json()
    if not isinstance(data, list):
        raise HandlerError("Unexpected MyVariant batch response shape")
    return data


def fetch_myvariant_batch_resilient(
    myvariant_ids: list[str],
    *,
    assembly: str = "hg38",
    timeout: int = DEFAULT_TIMEOUT,
) -> list[dict]:
    """Fetch a batch, recursively splitting chunks that trigger server errors."""
    try:
        return fetch_myvariant_batch(myvariant_ids, assembly=assembly, timeout=timeout)
    except requests.RequestException as e:
        if len(myvariant_ids) <= 1:
            return [{"query": myvariant_ids[0], "_batch_error": str(e)}]
        mid = len(myvariant_ids) // 2
        return (
            fetch_myvariant_batch_resilient(myvariant_ids[:mid], assembly=assembly, timeout=timeout)
            + fetch_myvariant_batch_resilient(myvariant_ids[mid:], assembly=assembly, timeout=timeout)
        )


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
        _coerce_numeric(node.get("score") if "score" in node else node.get("coding_score")),
        _coerce_numeric(
            node.get("rankscore")
            if "rankscore" in node
            else node.get("converted_rankscore")
            if "converted_rankscore" in node
            else node.get("coding_rankscore")
        ),
        _coerce_category(node.get("pred") if "pred" in node else node.get("coding_pred")),
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

    # VARITY has four related score families.
    varity = dbnsfp.get("varity") or {}
    for sub_name, our_name in (
        ("r", "varity_r"),
        ("er", "varity_er"),
        ("r_loo", "varity_r_loo"),
        ("er_loo", "varity_er_loo"),
    ):
        s, r, c = _extract_score_triplet(varity.get(sub_name))
        if s is not None or r is not None or c is not None:
            yield {"predictor": our_name, "score": s, "rank_score": r, "category": c}

    # Eigen/Eigen-PC and fitCons use predictor-specific field names.
    eigen = dbnsfp.get("eigen") or {}
    if isinstance(eigen, dict):
        s = _coerce_numeric(eigen.get("raw_coding"))
        if s is not None:
            yield {
                "predictor": "eigen_raw_coding",
                "score": s,
                "rank_score": _coerce_numeric(eigen.get("raw_coding_rankscore")),
                "category": None,
            }
    eigen_pc = dbnsfp.get("eigen-pc") or {}
    if isinstance(eigen_pc, dict):
        s = _coerce_numeric(eigen_pc.get("raw_coding"))
        if s is not None:
            yield {
                "predictor": "eigen_pc_raw_coding",
                "score": s,
                "rank_score": _coerce_numeric(eigen_pc.get("raw_coding_rankscore")),
                "category": None,
            }

    fitcons = dbnsfp.get("fitcons") or {}
    integrated = fitcons.get("integrated") if isinstance(fitcons, dict) else None
    if isinstance(integrated, dict):
        s, r, c = _extract_score_triplet(integrated)
        if s is not None or r is not None or c is not None:
            yield {"predictor": "integrated_fitcons", "score": s, "rank_score": r, "category": c}


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

def persist(db, variant_id: int, payload: dict, *, source_label: str = SOURCE) -> dict:
    """Write all parseable annotation rows from a MyVariant payload. Returns counts."""
    counts = {"pathogenicity": 0, "conservation": 0, "population": 0, "clinical": 0, "aliases": 0}

    for row in parse_pathogenicity(payload):
        db.upsert_pathogenicity(variant_id, source=source_label, **row)
        counts["pathogenicity"] += 1

    for row in parse_conservation(payload):
        db.upsert_conservation(variant_id, source=source_label, **row)
        counts["conservation"] += 1

    for row in parse_population(payload):
        db.upsert_population(variant_id, source=source_label, **row)
        counts["population"] += 1

    for row in parse_clinical(payload):
        db.upsert_clinical(variant_id, **row)
        counts["clinical"] += 1

    aliases = list(parse_aliases(payload))
    if aliases:
        counts["aliases"] = db.add_aliases(variant_id, aliases, source=source_label)

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


def run_batch(
    db,
    *,
    limit: Optional[int] = None,
    batch_size: int = 500,
    source_label: Optional[str] = None,
    timeout: int = DEFAULT_TIMEOUT,
) -> dict:
    """Drain pending MyVariant jobs through the POST batch endpoint."""
    jobs = db.claim_pending_jobs(source=SOURCE, limit=limit if limit is not None else 1_000_000)
    if not jobs:
        return {"claimed": 0, "done": 0, "failed": 0, "notfound": 0, "batches": 0}

    ids = [j["id"] for j in jobs]
    placeholders = ",".join("?" * len(ids))
    cur = db.conn.execute(
        f"""
        SELECT j.id AS job_id, j.variant_id AS variant_id,
               v.chromosome AS chromosome, v.position AS position, v.ref AS ref, v.alt AS alt
        FROM annotation_jobs j JOIN variants v ON v.id = j.variant_id
        WHERE j.id IN ({placeholders})
        """,
        ids,
    )
    rows = [dict(r) for r in cur.fetchall()]
    source_label = source_label or f"myvariant_batch_hg38_{date.today().isoformat()}"

    requests_to_make: list[dict] = []
    failed = 0
    for row in rows:
        if len(row["ref"]) != 1 or len(row["alt"]) != 1:
            db.mark_job_failed(row["job_id"], "myvariant batch handler currently only handles SNVs")
            failed += 1
            continue
        try:
            row["query"] = _myvariant_id(row["chromosome"], row["position"], row["ref"], row["alt"])
        except HandlerError as e:
            db.mark_job_failed(row["job_id"], str(e))
            failed += 1
            continue
        requests_to_make.append(row)

    done = 0
    notfound = 0
    batches = 0
    for i in range(0, len(requests_to_make), batch_size):
        chunk = requests_to_make[i : i + batch_size]
        payloads = fetch_myvariant_batch_resilient([r["query"] for r in chunk], timeout=timeout)
        by_query = {p.get("query") or p.get("_id"): p for p in payloads if isinstance(p, dict)}
        found_payloads = []
        done_job_ids = []
        for row in chunk:
            payload = by_query.get(row["query"])
            if payload and payload.get("_batch_error"):
                db.mark_job_failed(row["job_id"], payload["_batch_error"])
                failed += 1
                continue
            if not payload or payload.get("notfound"):
                done_job_ids.append(row["job_id"])
                done += 1
                notfound += 1
                continue
            found_payloads.append((row["variant_id"], payload))
            done_job_ids.append(row["job_id"])
            done += 1
        _persist_many(db, found_payloads, source_label=source_label)
        _mark_jobs_done_bulk(db, done_job_ids)
        batches += 1

    return {
        "claimed": len(jobs),
        "done": done,
        "failed": failed,
        "notfound": notfound,
        "batches": batches,
    }


def _persist_many(db, items: list[tuple[int, dict]], *, source_label: str) -> None:
    if not items:
        return

    pathogenicity_rows = []
    conservation_rows = []
    population_rows = []
    clinical_rows = []
    alias_rows = []
    ca_updates = []

    for variant_id, payload in items:
        for row in parse_pathogenicity(payload):
            pathogenicity_rows.append((
                variant_id,
                row["predictor"],
                row.get("predictor_version") or "",
                row.get("score"),
                row.get("rank_score"),
                row.get("category"),
                source_label,
            ))

        for row in parse_conservation(payload):
            conservation_rows.append((
                variant_id,
                row["metric"],
                row.get("score"),
                row.get("rank_score"),
                source_label,
            ))

        for row in parse_population(payload):
            population_rows.append((
                variant_id,
                row["dataset"],
                row.get("pop") or "all",
                row.get("af"),
                row.get("ac"),
                row.get("an"),
                row.get("n_homozygotes"),
                row.get("filter_status"),
                source_label,
            ))

        for row in parse_clinical(payload):
            clinical_rows.append((
                variant_id,
                row["source"],
                row.get("record_id") or "",
                row.get("classification"),
                row.get("review_status"),
                row.get("stars"),
                row.get("last_evaluated"),
                row.get("conditions"),
            ))

        for row in parse_aliases(payload):
            alias_rows.append((variant_id, row["alias_type"], row["alias_value"], source_label))

        caid = (payload.get("clingen") or {}).get("caid")
        if caid:
            ca_updates.append((caid, variant_id))

    if pathogenicity_rows:
        db.conn.executemany(
            """
            INSERT INTO annotations_pathogenicity
                (variant_id, predictor, predictor_version, score, rank_score, category, source)
            VALUES (?, ?, ?, ?, ?, ?, ?)
            ON CONFLICT(variant_id, predictor, predictor_version) DO UPDATE SET
                score = excluded.score,
                rank_score = excluded.rank_score,
                category = excluded.category,
                source = excluded.source,
                fetched_at = CURRENT_TIMESTAMP
            """,
            pathogenicity_rows,
        )

    if conservation_rows:
        db.conn.executemany(
            """
            INSERT INTO annotations_conservation
                (variant_id, metric, score, rank_score, source)
            VALUES (?, ?, ?, ?, ?)
            ON CONFLICT(variant_id, metric) DO UPDATE SET
                score = excluded.score,
                rank_score = excluded.rank_score,
                source = excluded.source,
                fetched_at = CURRENT_TIMESTAMP
            """,
            conservation_rows,
        )

    if population_rows:
        db.conn.executemany(
            """
            INSERT INTO annotations_population
                (variant_id, dataset, pop, af, ac, an, n_homozygotes, filter_status, source)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
            ON CONFLICT(variant_id, dataset, pop) DO UPDATE SET
                af = excluded.af,
                ac = excluded.ac,
                an = excluded.an,
                n_homozygotes = excluded.n_homozygotes,
                filter_status = excluded.filter_status,
                source = excluded.source,
                fetched_at = CURRENT_TIMESTAMP
            """,
            population_rows,
        )

    if clinical_rows:
        db.conn.executemany(
            """
            INSERT INTO annotations_clinical
                (variant_id, source, record_id, classification, review_status, stars, last_evaluated, conditions)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?)
            ON CONFLICT(variant_id, source, record_id) DO UPDATE SET
                classification = excluded.classification,
                review_status = excluded.review_status,
                stars = excluded.stars,
                last_evaluated = excluded.last_evaluated,
                conditions = excluded.conditions,
                fetched_at = CURRENT_TIMESTAMP
            """,
            clinical_rows,
        )

    if alias_rows:
        db.conn.executemany(
            """
            INSERT OR IGNORE INTO variant_aliases
                (variant_id, alias_type, alias_value, source)
            VALUES (?, ?, ?, ?)
            """,
            alias_rows,
        )

    for caid, variant_id in ca_updates:
        db.conn.execute(
            "UPDATE variants SET ca_id = COALESCE(ca_id, ?), updated_at = CURRENT_TIMESTAMP WHERE id = ?",
            [caid, variant_id],
        )

    db.conn.commit()


def _mark_jobs_done_bulk(db, job_ids: list[int]) -> None:
    if not job_ids:
        return
    placeholders = ",".join("?" * len(job_ids))
    db.conn.execute(
        f"""
        UPDATE annotation_jobs
        SET status = 'done', error = NULL, updated_at = CURRENT_TIMESTAMP
        WHERE id IN ({placeholders})
        """,
        job_ids,
    )
    db.conn.commit()


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
