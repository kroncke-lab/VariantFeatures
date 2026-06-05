from __future__ import annotations

from pathlib import Path

from variantfeatures.database import VariantDB
from variantfeatures.handlers import cadd


def test_cadd_handler_persists_phred_and_raw(monkeypatch, tmp_path: Path):
    db = VariantDB(tmp_path / "test.db")
    vid = db.upsert_variant(chromosome="7", position=10, ref="A", alt="G", variant_type="SNV")

    def fake_fetch(chromosome, position, ref=None, alt=None, version=None):
        assert (chromosome, position, ref, alt) == ("7", 10, "A", "G")
        return {"cadd_phred": 23.4, "cadd_raw": 4.56}

    monkeypatch.setattr(cadd, "fetch_cadd_single", fake_fetch)

    cadd.handle(db, vid)

    rows = [
        dict(r)
        for r in db.conn.execute(
            "SELECT predictor, predictor_version, score, source FROM annotations_pathogenicity WHERE variant_id = ? ORDER BY predictor",
            [vid],
        )
    ]
    assert rows == [
        {
            "predictor": "cadd_phred",
            "predictor_version": "GRCh38-v1.6",
            "score": 23.4,
            "source": "cadd_api",
        },
        {
            "predictor": "cadd_raw",
            "predictor_version": "GRCh38-v1.6",
            "score": 4.56,
            "source": "cadd_api",
        },
    ]
