"""Unit tests for variantfeatures.identity (resolver + ClinGen AR parser).

The integration test against the live ClinGen Allele Registry is opt-in
(set RUN_INTEGRATION=1).
"""

from __future__ import annotations

import os
from unittest.mock import patch

import pytest

from variantfeatures import identity


# ---------------------------------------------------------------------------
# Hand-crafted ClinGen AR response fixture (KCNH2 c.1841T>C / p.M614T-ish)
# ---------------------------------------------------------------------------
SAMPLE_PAYLOAD = {
    "@id": "https://reg.genome.network/allele/CA12345678",
    "type": "nucleotide",
    "genomicAlleles": [
        {
            "referenceGenome": "GRCh37",
            "chromosome": "7",
            "hgvs": ["NC_000007.13:g.150951925A>G"],
            "coordinates": [
                {"start": 150951924, "end": 150951925, "referenceAllele": "A", "allele": "G"}
            ],
        },
        {
            "referenceGenome": "GRCh38",
            "chromosome": "7",
            "hgvs": ["NC_000007.14:g.151251182A>G"],
            "coordinates": [
                {"start": 151251181, "end": 151251182, "referenceAllele": "A", "allele": "G"}
            ],
        },
    ],
    "transcriptAlleles": [
        {
            "hgvs": ["NM_000238.4:c.1841T>C"],
            "geneSymbol": "KCNH2",
            "proteinEffect": {"hgvs": "NP_000229.1:p.Met614Thr"},
        }
    ],
    "externalRecords": {
        "dbSNP": [{"rs": 199472953}],
        "ClinVarVariations": [{"variationId": 67244, "RCV": ["RCV000058854"]}],
        "ClinVarAlleles": [{"alleleId": 82268}],
        "gnomAD": [{"id": "7-151251182-A-G"}],
        "MyVariantInfo_hg38": [{"id": "chr7:g.151251182A>G"}],
    },
}


# ---------------------------------------------------------------------------
# detect_kind
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    "value,expected",
    [
        ("CA12345678", "ca_id"),
        ("rs199472953", "rsid"),
        ("RS199472953", "rsid"),
        ("VCV000067244", "clinvar_vcv"),
        ("RCV000058854", "clinvar_rcv"),
        ("NM_000238.4:c.1841T>C", "hgvs_c"),
        ("NC_000007.14:g.151251182A>G", "hgvs_g"),
        ("NP_000229.1:p.Met614Thr", "hgvs_p"),
        ("7-151251182-A-G", "gnomad_id"),
        ("chr7:151251182:A:G", "gnomad_id"),
        ("67244", "clinvar_variation"),
    ],
)
def test_detect_kind(value, expected):
    assert identity.detect_kind(value) == expected


def test_detect_kind_unknown():
    with pytest.raises(identity.IdentityError):
        identity.detect_kind("not-a-real-identifier")


# ---------------------------------------------------------------------------
# _classify_variant
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    "ref,alt,expected",
    [
        ("A", "G", "SNV"),
        ("AC", "GT", "MNV"),
        ("A", "AT", "INS"),
        ("AT", "A", "DEL"),
        ("AT", "GC", "MNV"),
        ("ATG", "AC", "DELINS"),
    ],
)
def test_classify_variant(ref, alt, expected):
    assert identity._classify_variant(ref, alt) == expected


# ---------------------------------------------------------------------------
# _parse_clingen against fixture
# ---------------------------------------------------------------------------
def test_parse_clingen_canonical_fields():
    rv = identity._parse_clingen(SAMPLE_PAYLOAD)
    assert rv.ca_id == "CA12345678"
    assert rv.chromosome == "7"
    assert rv.position == 151251182  # 1-based, GRCh38
    assert rv.ref == "A"
    assert rv.alt == "G"
    assert rv.variant_type == "SNV"
    assert rv.hgvs_g == "NC_000007.14:g.151251182A>G"


def test_parse_clingen_alias_inventory():
    rv = identity._parse_clingen(SAMPLE_PAYLOAD)
    pairs = {(a["alias_type"], a["alias_value"]) for a in rv.aliases}

    # CA-ID and both build g.HGVS
    assert ("ca_id", "CA12345678") in pairs
    assert ("hgvs_g", "NC_000007.13:g.150951925A>G") in pairs
    assert ("hgvs_g", "NC_000007.14:g.151251182A>G") in pairs

    # Coding/protein
    assert ("hgvs_c", "NM_000238.4:c.1841T>C") in pairs
    assert ("hgvs_p", "NP_000229.1:p.Met614Thr") in pairs

    # External records
    assert ("rsid", "rs199472953") in pairs
    assert ("clinvar_vcv", "VCV000067244") in pairs
    assert ("clinvar_rcv", "RCV000058854") in pairs
    assert ("clinvar_allele", "82268") in pairs
    assert ("gnomad_id", "7-151251182-A-G") in pairs
    assert ("myvariant_id", "chr7:g.151251182A>G") in pairs


def test_parse_clingen_handles_missing_grch38():
    payload = {
        "@id": "https://reg.genome.network/allele/CA1",
        "genomicAlleles": [
            {
                "referenceGenome": "GRCh37",
                "chromosome": "1",
                "hgvs": ["NC_000001.10:g.100A>G"],
                "coordinates": [{"start": 99, "end": 100, "referenceAllele": "A", "allele": "G"}],
            }
        ],
    }
    rv = identity._parse_clingen(payload)
    assert rv.ca_id == "CA1"
    # No GRCh38 entry -> coordinates not set
    assert rv.chromosome is None
    assert rv.position is None
    # But the GRCh37 g.HGVS still ends up as an alias
    assert any(a["alias_type"] == "hgvs_g" for a in rv.aliases)


# ---------------------------------------------------------------------------
# resolve() with mocked HTTP layer
# ---------------------------------------------------------------------------
class _FakeResp:
    def __init__(self, status_code=200, payload=None):
        self.status_code = status_code
        self._payload = payload or {}

    def json(self):
        return self._payload

    def raise_for_status(self):
        if self.status_code >= 400:
            from requests import HTTPError

            raise HTTPError(f"{self.status_code}")


def test_resolve_routes_hgvs_c_to_hgvs_param():
    captured = {}

    def fake_get(url, params=None, timeout=None):
        captured["url"] = url
        captured["params"] = params
        return _FakeResp(200, SAMPLE_PAYLOAD)

    with patch("variantfeatures.identity.requests.get", side_effect=fake_get):
        rv = identity.resolve("NM_000238.4:c.1841T>C")

    assert captured["params"] == {"hgvs": "NM_000238.4:c.1841T>C"}
    assert rv.ca_id == "CA12345678"


def test_resolve_routes_rsid():
    captured = {}

    def fake_get(url, params=None, timeout=None):
        captured["url"] = url
        captured["params"] = params
        # /alleles (plural) returns a JSON list; route returns first
        return _FakeResp(200, [SAMPLE_PAYLOAD])

    with patch("variantfeatures.identity.requests.get", side_effect=fake_get):
        identity.resolve("rs199472953")

    assert captured["url"].endswith("/alleles")
    assert captured["params"] == {"dbSNP.rs": "199472953"}


def test_resolve_clinvar_variation_rejected_with_clear_message():
    """ClinGen AR doesn't index ClinVar variationId directly; we should fail loudly."""
    with pytest.raises(identity.IdentityError, match="ClinVar"):
        identity.resolve("VCV000067244")
    with pytest.raises(identity.IdentityError, match="ClinVar"):
        identity.resolve("67244")


def test_resolve_routes_ca_id_to_path():
    captured = {}

    def fake_get(url, params=None, timeout=None):
        captured["url"] = url
        captured["params"] = params
        return _FakeResp(200, SAMPLE_PAYLOAD)

    with patch("variantfeatures.identity.requests.get", side_effect=fake_get):
        identity.resolve("CA12345678")

    assert captured["url"].endswith("/allele/CA12345678")
    assert captured["params"] == {}


def test_resolve_protein_hgvs_rejected():
    with pytest.raises(identity.IdentityError):
        identity.resolve("NP_000229.1:p.Met614Thr")


def test_resolve_404():
    def fake_get(url, params=None, timeout=None):
        return _FakeResp(404)

    with patch("variantfeatures.identity.requests.get", side_effect=fake_get):
        with pytest.raises(identity.IdentityError):
            identity.resolve("NM_999999.1:c.1A>T")


# ---------------------------------------------------------------------------
# Live integration (opt-in)
# ---------------------------------------------------------------------------
@pytest.mark.skipif(
    os.environ.get("RUN_INTEGRATION") != "1",
    reason="Set RUN_INTEGRATION=1 to hit the live ClinGen Allele Registry.",
)
def test_resolve_braf_v600e_live():
    """BRAF V600E (rs113488022, CA123643) is a well-tested allele in the registry."""
    rv = identity.resolve("NC_000007.13:g.140453136A>T")
    assert rv.ca_id == "CA123643"
    assert rv.chromosome == "7"
    # GRCh38 BRAF V600E is at chr7:140753336
    assert rv.position == 140753336
    assert any(a["alias_type"] == "rsid" and a["alias_value"] == "rs113488022" for a in rv.aliases)
    # Alias inventory should at least include c. and p. forms
    types = {a["alias_type"] for a in rv.aliases}
    assert "hgvs_c" in types
    assert "hgvs_p" in types


@pytest.mark.skipif(
    os.environ.get("RUN_INTEGRATION") != "1",
    reason="Set RUN_INTEGRATION=1 to hit the live ClinGen Allele Registry.",
)
def test_resolve_braf_v600e_via_rsid_live():
    rv = identity.resolve("rs113488022")
    assert rv.ca_id == "CA123643"
    assert rv.chromosome == "7"
    assert rv.position == 140753336
