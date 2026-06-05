from __future__ import annotations

from pathlib import Path

import pytest

from variantfeatures.database import VariantDB
from variantfeatures.handlers import absplice, nmd_external, pext


@pytest.fixture
def db(tmp_path: Path) -> VariantDB:
    return VariantDB(tmp_path / "test.db")


@pytest.fixture
def variant(db) -> int:
    vid = db.upsert_variant(chromosome="7", position=140753336, ref="A", alt="T", variant_type="SNV")
    db.upsert_consequence(
        variant_id=vid,
        transcript_id="ENST1",
        source="enumerated",
        gene_symbol="BRAF",
        consequence="missense_variant",
        aa_pos=600,
    )
    return vid


def test_pext_import_long_form(db, variant, tmp_path):
    path = tmp_path / "pext.tsv"
    path.write_text("chrom\tpos\ttissue\tpext\n7\t140753336\theart\t0.76\n")

    result = pext.import_file(db, path, gene_filter="BRAF")
    assert result == {"rows": 1, "matched_rows": 1, "annotations": 1}

    cur = db.conn.execute(
        "SELECT metric, dataset, tissue, score FROM annotations_expression WHERE variant_id = ?",
        [variant],
    )
    assert dict(cur.fetchone()) == {
        "metric": "pext",
        "dataset": "gnomad_pext_v10",
        "tissue": "heart",
        "score": 0.76,
    }


def test_absplice_import(db, variant, tmp_path):
    path = tmp_path / "absplice.tsv"
    path.write_text(
        "variant\tgene_id\ttissue\tAbSplice_DNA\tAbExp\n"
        "chr7:140753336:A>T\tENSG00000157764\theart\t0.84\t-1.5\n"
    )

    result = absplice.import_file(db, path, gene_filter="BRAF")
    assert result == {"rows": 1, "matched_rows": 1, "splice": 1, "expression": 1}

    cur = db.conn.execute(
        "SELECT predictor, score_type, score FROM annotations_splice WHERE variant_id = ?",
        [variant],
    )
    assert dict(cur.fetchone()) == {"predictor": "absplice_dna", "score_type": "heart", "score": 0.84}

    cur = db.conn.execute(
        "SELECT metric, tissue, score FROM annotations_expression WHERE variant_id = ? AND metric = 'abexp'",
        [variant],
    )
    assert dict(cur.fetchone()) == {"metric": "abexp", "tissue": "heart", "score": -1.5}


def test_nmd_external_import(db, variant, tmp_path):
    path = tmp_path / "nmdep.csv"
    path.write_text("variant,nmdep_score,prediction\nchr7:140753336:A>T,0.42,escape\n")

    result = nmd_external.import_file(db, path, predictor="nmdep", gene_filter="BRAF")
    assert result == {"rows": 1, "matched_rows": 1, "annotations": 1}

    cur = db.conn.execute(
        "SELECT predictor, score, category, source FROM annotations_pathogenicity WHERE variant_id = ?",
        [variant],
    )
    assert dict(cur.fetchone()) == {
        "predictor": "nmdep",
        "score": 0.42,
        "category": "escape",
        "source": "nmd_external",
    }
