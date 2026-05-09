"""Unit tests for the MyVariant.info handler."""

from __future__ import annotations

import os
from pathlib import Path
from unittest.mock import patch

import pytest

from variantfeatures.database import VariantDB
from variantfeatures.handlers import myvariant as mv


# Slimmed-down BRAF V600E payload that exercises every parser branch.
SAMPLE = {
    "_id": "chr7:g.140753336A>T",
    "chrom": "7",
    "clingen": {"caid": "CA123643"},
    "dbsnp": {"rsid": "rs113488022"},
    "clinvar": {
        "variant_id": 13961,
        "allele_id": 29000,
        "rcv": [
            {
                "accession": "RCV000014992",
                "clinical_significance": "Pathogenic",
                "review_status": "no assertion criteria provided",
                "last_evaluated": "2014-09-04",
                "conditions": {"name": "Carcinoma of colon"},
            },
            {
                "accession": "RCV000014993",
                "clinical_significance": "Pathogenic",
                "review_status": "no assertion criteria provided",
                "last_evaluated": "2014-09-04",
                "conditions": {"name": "Papillary thyroid carcinoma"},
            },
        ],
    },
    "cadd": {"phred": 29.8, "raw_score": 5.2},
    "dbnsfp": {
        "cadd": {"phred": 29.8, "raw_score": 5.298, "raw_rankscore": 0.89},
        "alphamissense": {
            "score": [0.985, 0.998, 0.993],  # multiple transcripts
            "rankscore": 0.972,
            "pred": ["P", "P", "P"],
        },
        "revel": {"score": 0.931, "rankscore": 0.984},
        "metasvm": {"score": -0.77, "rankscore": 0.57, "pred": "T"},
        "metalr": {"score": 0.64, "rankscore": 0.66, "pred": "D"},
        "m-cap": {"score": 0.10, "rankscore": 0.50, "pred": "D"},
        "primateai": {"score": 0.85, "rankscore": 0.92, "pred": "D"},
        "clinpred": {"score": 0.99, "rankscore": 0.95, "pred": "D"},
        "esm1b": {"score": -8.5, "rankscore": 0.95},
        "eve": {"score": 0.80, "rankscore": 0.88, "pred": "P"},
        "polyphen2": {
            "hdiv": {"score": 1.0, "rankscore": 0.95, "pred": "D"},
            "hvar": {"score": 1.0, "rankscore": 0.93, "pred": "D"},
        },
        "bayesdel": {
            "add_af": {"score": 0.40, "rankscore": 0.90, "pred": "D"},
            "no_af": {"score": 0.34, "rankscore": 0.90, "pred": "D"},
        },
        "phylop": {
            "100way_vertebrate": {"score": 9.236, "rankscore": 0.94},
            "470way_mammalian": {"score": 11.21, "rankscore": 0.90},
        },
        "phastcons": {
            "100way_vertebrate": {"score": 1.0, "rankscore": 0.71},
        },
        "gerp++": {"rs": 5.65, "nr": 5.65, "rs_rankscore": 0.87},
        "siphy_29way": {"logodds": 16.0, "logodds_rankscore": 0.85},
    },
    "gnomad_exome": {
        "af": {"af": 4e-06, "af_afr": 0.0, "af_sas": 3.27e-05, "af_nfe": 0.0},
        "ac": {"ac": 1, "ac_afr": 0, "ac_sas": 1, "ac_nfe": 0},
        "an": {"an": 251260, "an_afr": 16252, "an_sas": 30000, "an_nfe": 100000},
        "hom": {"hom": 0, "hom_afr": 0, "hom_sas": 0, "hom_nfe": 0},
    },
}


@pytest.fixture
def db(tmp_path: Path) -> VariantDB:
    return VariantDB(tmp_path / "test.db")


@pytest.fixture
def vid(db) -> int:
    return db.upsert_variant(chromosome="7", position=140753336, ref="A", alt="T", variant_type="SNV")


# ---------------------------------------------------------------------------
# myvariant_id
# ---------------------------------------------------------------------------

def test_myvariant_id_for_snv():
    assert mv._myvariant_id("7", 140753336, "A", "T") == "chr7:g.140753336A>T"


def test_myvariant_id_strips_chr_prefix():
    assert mv._myvariant_id("chr7", 1, "A", "G") == "chr7:g.1A>G"


def test_myvariant_id_unknown_chromosome():
    with pytest.raises(mv.HandlerError):
        mv._myvariant_id("99", 1, "A", "G")


# ---------------------------------------------------------------------------
# parse_pathogenicity
# ---------------------------------------------------------------------------

def test_parse_pathogenicity_includes_cadd_and_dbnsfp():
    rows = list(mv.parse_pathogenicity(SAMPLE))
    by_name = {r["predictor"]: r for r in rows}

    # CADD
    assert by_name["cadd_phred"]["score"] == 29.8
    assert by_name["cadd_raw"]["score"] == 5.298
    assert by_name["cadd_raw"]["rank_score"] == 0.89

    # AlphaMissense — score is a list, parser collapses to max
    assert by_name["alphamissense"]["score"] == 0.998
    assert by_name["alphamissense"]["rank_score"] == 0.972
    assert by_name["alphamissense"]["category"] == "P"

    # REVEL
    assert by_name["revel"]["score"] == 0.931

    # Nested predictors
    assert "polyphen2_hdiv" in by_name
    assert "polyphen2_hvar" in by_name
    assert by_name["polyphen2_hdiv"]["score"] == 1.0
    assert "bayesdel_add_af" in by_name
    assert by_name["bayesdel_add_af"]["score"] == 0.40

    # Sanity: at least 10 predictors covered
    assert len(rows) >= 10


def test_parse_pathogenicity_handles_missing_dbnsfp():
    rows = list(mv.parse_pathogenicity({"cadd": {"phred": 10}}))
    assert len(rows) == 1
    assert rows[0]["predictor"] == "cadd_phred"


def test_parse_pathogenicity_skips_predictor_with_only_nulls():
    rows = list(mv.parse_pathogenicity({"dbnsfp": {"revel": {}}}))
    assert rows == []


# ---------------------------------------------------------------------------
# parse_conservation
# ---------------------------------------------------------------------------

def test_parse_conservation_basic():
    rows = list(mv.parse_conservation(SAMPLE))
    by_metric = {r["metric"]: r for r in rows}
    assert "phylop100way_vertebrate" in by_metric
    assert by_metric["phylop100way_vertebrate"]["score"] == 9.236
    assert "phastcons100way_vertebrate" in by_metric
    assert "gerp_pp_rs" in by_metric
    assert by_metric["gerp_pp_rs"]["score"] == 5.65
    assert "siphy_29way_logodds" in by_metric


# ---------------------------------------------------------------------------
# parse_population
# ---------------------------------------------------------------------------

def test_parse_population_skips_empty_pops():
    rows = list(mv.parse_population(SAMPLE))
    by_pop = {(r["dataset"], r["pop"]): r for r in rows}

    assert ("gnomad_exome", "all") in by_pop
    assert by_pop[("gnomad_exome", "all")]["af"] == 4e-06
    assert by_pop[("gnomad_exome", "all")]["ac"] == 1
    assert by_pop[("gnomad_exome", "all")]["an"] == 251260
    assert by_pop[("gnomad_exome", "all")]["n_homozygotes"] == 0

    assert ("gnomad_exome", "sas") in by_pop
    assert by_pop[("gnomad_exome", "sas")]["af"] == 3.27e-05

    # popmax was not in fixture -> should be skipped
    assert ("gnomad_exome", "popmax") not in by_pop

    # gnomad_genome absent in fixture
    assert not any(r["dataset"] == "gnomad_genome" for r in rows)


# ---------------------------------------------------------------------------
# parse_clinical
# ---------------------------------------------------------------------------

def test_parse_clinical_emits_one_row_per_rcv():
    rows = list(mv.parse_clinical(SAMPLE))
    assert len(rows) == 2
    assert {r["record_id"] for r in rows} == {"RCV000014992", "RCV000014993"}
    for r in rows:
        assert r["source"] == "clinvar"
        assert r["classification"] == "Pathogenic"
        assert r["last_evaluated"] == "2014-09-04"


def test_parse_clinical_falls_back_to_summary_if_no_rcv():
    payload = {"clinvar": {
        "variant_id": 12345, "clinical_significance": "VUS",
        "review_status": "criteria provided", "gold_stars": 1,
    }}
    rows = list(mv.parse_clinical(payload))
    assert len(rows) == 1
    assert rows[0]["classification"] == "VUS"
    assert rows[0]["stars"] == 1


# ---------------------------------------------------------------------------
# parse_aliases
# ---------------------------------------------------------------------------

def test_parse_aliases_includes_clingen_dbsnp_clinvar():
    out = list(mv.parse_aliases(SAMPLE))
    pairs = {(a["alias_type"], a["alias_value"]) for a in out}
    assert ("ca_id", "CA123643") in pairs
    assert ("rsid", "rs113488022") in pairs
    assert ("clinvar_vcv", "VCV000013961") in pairs
    assert ("clinvar_allele", "29000") in pairs
    assert ("myvariant_id", "chr7:g.140753336A>T") in pairs


# ---------------------------------------------------------------------------
# persist + handle (DB-side)
# ---------------------------------------------------------------------------

def test_persist_writes_to_all_tables(db, vid):
    counts = mv.persist(db, vid, SAMPLE)

    assert counts["pathogenicity"] >= 10
    assert counts["conservation"] >= 3
    assert counts["population"] >= 4
    assert counts["clinical"] == 2
    assert counts["aliases"] >= 4

    # Spot check rows
    cur = db.conn.execute("SELECT predictor, score FROM annotations_pathogenicity WHERE variant_id = ? AND predictor = 'alphamissense'", [vid])
    assert dict(cur.fetchone())["score"] == 0.998

    cur = db.conn.execute("SELECT classification FROM annotations_clinical WHERE variant_id = ? AND record_id = 'RCV000014992'", [vid])
    assert dict(cur.fetchone())["classification"] == "Pathogenic"

    cur = db.conn.execute("SELECT af FROM annotations_population WHERE variant_id = ? AND dataset = 'gnomad_exome' AND pop = 'all'", [vid])
    assert dict(cur.fetchone())["af"] == 4e-06


def test_persist_updates_canonical_ca_id_when_missing(db, vid):
    mv.persist(db, vid, SAMPLE)
    cur = db.conn.execute("SELECT ca_id FROM variants WHERE id = ?", [vid])
    assert dict(cur.fetchone())["ca_id"] == "CA123643"


class _FakeResp:
    def __init__(self, status_code=200, payload=None):
        self.status_code = status_code
        self._payload = payload

    def json(self):
        return self._payload

    def raise_for_status(self):
        if self.status_code >= 400:
            from requests import HTTPError
            raise HTTPError(f"{self.status_code}")


def test_handle_404_is_not_an_error(db, vid):
    """If MyVariant.info doesn't have the variant, the job is treated as done with no rows written."""
    with patch("variantfeatures.handlers.myvariant.requests.get", return_value=_FakeResp(404)):
        mv.handle(db, vid)  # should not raise

    cur = db.conn.execute("SELECT COUNT(*) FROM annotations_pathogenicity WHERE variant_id = ?", [vid])
    assert cur.fetchone()[0] == 0


def test_handle_persists_when_response_present(db, vid):
    with patch("variantfeatures.handlers.myvariant.requests.get", return_value=_FakeResp(200, SAMPLE)):
        mv.handle(db, vid)
    cur = db.conn.execute("SELECT COUNT(*) FROM annotations_pathogenicity WHERE variant_id = ?", [vid])
    assert cur.fetchone()[0] >= 10


def test_handle_rejects_indels(db):
    vid = db.upsert_variant(chromosome="7", position=1, ref="A", alt="AT")
    with pytest.raises(mv.HandlerError):
        mv.handle(db, vid)


# ---------------------------------------------------------------------------
# Live integration (BRAF V600E)
# ---------------------------------------------------------------------------
@pytest.mark.skipif(
    os.environ.get("RUN_INTEGRATION") != "1",
    reason="Set RUN_INTEGRATION=1 to hit MyVariant.info live.",
)
def test_handle_braf_v600e_live(tmp_path):
    """Hit the real MyVariant.info; verify CADD/REVEL/AlphaMissense come back populated."""
    d = VariantDB(tmp_path / "live.db")
    vid = d.upsert_variant(chromosome="7", position=140753336, ref="A", alt="T", variant_type="SNV")
    mv.handle(d, vid, timeout=60)

    cur = d.conn.execute(
        "SELECT predictor, score FROM annotations_pathogenicity WHERE variant_id = ? AND predictor IN ('cadd_phred', 'revel', 'alphamissense')",
        [vid],
    )
    by = {r["predictor"]: r["score"] for r in cur.fetchall()}
    assert by.get("cadd_phred") is not None and by["cadd_phred"] > 20
    assert by.get("revel") is not None and by["revel"] > 0.5
    assert by.get("alphamissense") is not None and by["alphamissense"] > 0.5
