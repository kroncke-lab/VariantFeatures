from __future__ import annotations

import csv
from pathlib import Path

from click.testing import CliRunner

from variantfeatures.cli import main
from variantfeatures.database import VariantDB
from variantfeatures.normalized_export import build_wide_rows, iter_long_rows


def _seed(db: VariantDB) -> int:
    vid = db.upsert_variant(
        chromosome="7",
        position=150000001,
        ref="A",
        alt="G",
        ca_id="CA1",
        variant_type="SNV",
        hgvs_g="NC_000007.14:g.150000001A>G",
    )
    db.add_aliases(vid, [{"alias_type": "rsid", "alias_value": "rs123"}], source="test")
    db.upsert_consequence(
        vid,
        "ENST00000262186.10",
        "enumerated",
        gene_symbol="KCNH2",
        consequence="missense_variant",
        hgvs_c="ENST00000262186.10:c.1A>G",
        hgvs_p="ENSP00000262186:p.Lys1Arg",
        aa_pos=1,
        aa_ref="K",
        aa_alt="R",
        is_mane_select=1,
    )
    db.upsert_pathogenicity(vid, "revel", predictor_version="1.3", score=0.91, source="revel")
    db.upsert_population(vid, "gnomad_exomes_v4", "all", af=1e-5, ac=2, an=200000, source="gnomad")
    db.upsert_splice(vid, "spliceai", "donor_gain", score=0.22, distance=4, source="vep")
    db.upsert_expression(vid, "pext", dataset="ucsc_gnomad_pext_hg38", tissue="Heart_LeftVentricle", score=0.73)
    db.upsert_structure(vid, "alphafold_plddt", feature_version="AFDB_v6", residue_number=1, score=95.0)
    db.upsert_clinical(vid, "clinvar", "VCV000000001", classification="Pathogenic", stars=2)
    db.upsert_gene_constraint("KCNH2", dataset="gnomad_v4", pli=0.99, oe_lof_upper=0.2)
    return vid


def test_wide_export_namespaces_feature_families(tmp_path: Path):
    db = VariantDB(tmp_path / "test.db")
    _seed(db)

    rows, fieldnames = build_wide_rows(db, "KCNH2")

    assert len(rows) == 1
    row = rows[0]
    assert row["aliases.rsid"] == "rs123"
    assert "pathogenicity.revel.version_1_3.score" in fieldnames
    assert row["pathogenicity.revel.version_1_3.score"] == 0.91
    assert row["population.gnomad_exomes_v4.all.af"] == 1e-5
    assert row["splice.spliceai.donor_gain.score"] == 0.22
    assert row["expression.pext.ucsc_gnomad_pext_hg38.Heart_LeftVentricle.score"] == 0.73
    assert row["structure.alphafold_plddt.version_AFDB_v6.score"] == 95.0
    assert row["clinical.clinvar.classification"] == "Pathogenic"
    assert row["gene_constraint.gnomad_v4.pli"] == 0.99


def test_long_export_identifies_storage_group(tmp_path: Path):
    db = VariantDB(tmp_path / "test.db")
    _seed(db)

    rows = list(iter_long_rows(db, "KCNH2", groups={"population", "expression"}))

    groups = {(r["feature_group"], r["feature_name"], r["field"]) for r in rows}
    assert ("population", "af", "af") in groups
    assert ("expression", "pext", "score") in groups


def test_transcript_wide_export_keeps_isoform_rows_and_scopes_expression(tmp_path: Path):
    db = VariantDB(tmp_path / "test.db")
    vid = db.upsert_variant(chromosome="7", position=150000001, ref="A", alt="G")
    db.upsert_consequence(
        vid,
        "ENST_CANON.1",
        "enumerated",
        gene_symbol="KCNH2",
        consequence="missense_variant",
        hgvs_p="P_CANON:p.Gly1Arg",
        is_canonical=1,
    )
    db.upsert_consequence(
        vid,
        "ENST_ALT.1",
        "enumerated",
        gene_symbol="KCNH2",
        consequence="synonymous_variant",
        hgvs_p="P_ALT:p.Gly1=",
    )
    db.upsert_expression(
        vid,
        "pext",
        dataset="test_pext",
        tissue="heart",
        transcript_id="ENST_ALT.1",
        score=0.03,
    )
    db.upsert_expression(
        vid,
        "pext",
        dataset="test_pext",
        tissue="brain",
        score=0.50,
    )

    rows, fieldnames = build_wide_rows(db, "KCNH2", row_grain="transcript")

    assert len(rows) == 2
    assert "expression.pext.test_pext.heart.ENST_ALT_1.score" in fieldnames
    assert "expression.pext.test_pext.brain.score" in fieldnames
    by_transcript = {row["transcript_id"]: row for row in rows}
    assert by_transcript["ENST_ALT.1"]["expression.pext.test_pext.heart.ENST_ALT_1.score"] == 0.03
    assert "expression.pext.test_pext.heart.ENST_ALT_1.score" not in by_transcript["ENST_CANON.1"]
    assert by_transcript["ENST_CANON.1"]["expression.pext.test_pext.brain.score"] == 0.50


def test_cli_export_defaults_to_normalized_wide(tmp_path: Path):
    db_path = tmp_path / "test.db"
    db = VariantDB(db_path)
    _seed(db)
    output = tmp_path / "kcnh2.csv"

    result = CliRunner().invoke(
        main,
        ["export", "-g", "KCNH2", "--db", str(db_path), "--output", str(output)],
    )

    assert result.exit_code == 0, result.output
    with open(output, newline="") as f:
        rows = list(csv.DictReader(f))
    assert rows[0]["pathogenicity.revel.version_1_3.score"] == "0.91"
    assert rows[0]["population.gnomad_exomes_v4.all.af"] == "1e-05"
