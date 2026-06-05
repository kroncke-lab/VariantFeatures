from __future__ import annotations

import csv
from pathlib import Path

from click.testing import CliRunner

from variantfeatures.cli import main
from variantfeatures.database import VariantDB


def _seed_normalized(db: VariantDB) -> tuple[int, int]:
    enumerated = db.upsert_variant(chromosome="7", position=150000001, ref="A", alt="G")
    db.upsert_consequence(
        enumerated,
        "ENST00000262186.10",
        "enumerated",
        gene_symbol="KCNH2",
        gene_ensembl="ENSG00000055118",
        consequence="missense_variant",
        hgvs_c="ENST00000262186.10:c.1A>G",
        hgvs_p="ENSP00000262186:p.Lys1Arg",
        is_mane_select=1,
    )
    db.upsert_pathogenicity(enumerated, "revel", predictor_version="1.3", score=0.91)
    db.upsert_population(enumerated, "gnomad_exomes_v4", "all", af=1e-5)
    db.upsert_clinical(enumerated, "clinvar", "VCV1", classification="Pathogenic")

    artifact = db.upsert_variant(chromosome="7", position=150000002, ref="C", alt="T")
    db.upsert_consequence(
        artifact,
        "NM_ARTIFACT.1",
        "annovar",
        gene_symbol="KCNH2",
        consequence="missense_variant",
        hgvs_c="c.2C>T",
        hgvs_p="p.Pro1Leu",
    )
    db.upsert_population(artifact, "gnomad_exomes_v4", "all", af=2e-5)
    return enumerated, artifact


def test_query_defaults_to_normalized_csv(tmp_path: Path):
    db_path = tmp_path / "test.db"
    db = VariantDB(db_path)
    _seed_normalized(db)

    result = CliRunner().invoke(
        main,
        ["query", "-g", "KCNH2", "--db", str(db_path), "--format", "csv"],
    )

    assert result.exit_code == 0, result.output
    rows = list(csv.DictReader(result.output.splitlines()))
    assert rows
    assert "pathogenicity.revel.version_1_3.score" in rows[0]
    assert "population.gnomad_exomes_v4.all.af" in rows[0]


def test_stats_uses_enumerated_denominator(tmp_path: Path):
    db_path = tmp_path / "test.db"
    db = VariantDB(db_path)
    _seed_normalized(db)

    result = CliRunner().invoke(main, ["stats", "--db", str(db_path)])

    assert result.exit_code == 0, result.output
    assert "Normalized Database Statistics" in result.output
    assert "Denominator: variant_consequences.source='enumerated'" in result.output
    assert "KCNH2" in result.output
    # Total and population coverage are both 1; the annovar-only artifact is ignored.
    assert "KCNH2             1        1        1        1" in result.output


def test_feature_coverage_uses_enumerated_denominator(tmp_path: Path):
    db_path = tmp_path / "test.db"
    db = VariantDB(db_path)
    _seed_normalized(db)

    result = CliRunner().invoke(
        main,
        ["feature-coverage", "-g", "KCNH2", "--db", str(db_path)],
    )

    assert result.exit_code == 0, result.output
    assert "KCNH2: 1 normalized enumerated variant(s)" in result.output
    assert "gnomad_exomes_v4" in result.output
    assert "  100.0%" in result.output
