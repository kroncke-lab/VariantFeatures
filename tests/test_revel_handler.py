from __future__ import annotations

from pathlib import Path

from variantfeatures.database import VariantDB
from variantfeatures.handlers import revel


def test_revel_run_batch_matches_coordinate_and_transcript(tmp_path: Path):
    db = VariantDB(tmp_path / "test.db")
    vid = db.upsert_variant(chromosome="7", position=10, ref="A", alt="G", variant_type="SNV")
    db.upsert_consequence(
        vid,
        "ENST000001.2",
        "enumerated",
        gene_symbol="KCNH2",
        consequence="missense_variant",
    )
    db.enqueue_job(vid, "revel")

    path = tmp_path / "revel.csv"
    path.write_text(
        "chr,hg19_pos,grch38_pos,ref,alt,aaref,aaalt,REVEL,Ensembl_transcriptid\n"
        "7,10,10,A,G,K,R,0.84,ENST000001\n"
    )

    result = revel.run_batch(db, file_path=str(path))

    assert result["claimed"] == 1
    assert result["matched"] == 1
    row = db.conn.execute(
        "SELECT predictor, predictor_version, score, source FROM annotations_pathogenicity WHERE variant_id = ?",
        [vid],
    ).fetchone()
    assert dict(row) == {
        "predictor": "revel",
        "predictor_version": "1.3",
        "score": 0.84,
        "source": "revel_zenodo_7072866",
    }
