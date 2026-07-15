from __future__ import annotations

from pathlib import Path

import pytest

from variantfeatures.database import VariantDB
from variantfeatures.handlers import alphafold


@pytest.fixture
def db(tmp_path: Path) -> VariantDB:
    return VariantDB(tmp_path / "test.db")


def test_alphafold_run_batch_writes_plddt(db, monkeypatch):
    vid = db.upsert_variant(chromosome="7", position=1, ref="A", alt="G", variant_type="SNV")
    db.upsert_consequence(
        variant_id=vid,
        transcript_id="ENST1",
        source="enumerated",
        gene_symbol="KCNH2",
        consequence="stop_gained",
        aa_pos=3,
    )
    db.enqueue_job(vid, "alphafold")

    # Gene->UniProt is resolved dynamically (no hardcoded map); stub it offline.
    monkeypatch.setattr(
        alphafold,
        "resolve_uniprot_accession",
        lambda gene, timeout=30, extra=None: "Q12809",
    )
    monkeypatch.setattr(
        alphafold,
        "fetch_predictions",
        lambda accession, timeout=30: [{
            "sequenceStart": 1,
            "sequenceEnd": 5,
            "latestVersion": 6,
            "plddtDocUrl": "https://example.org/confidence.json",
        }],
    )
    monkeypatch.setattr(
        alphafold,
        "fetch_confidence",
        lambda url, timeout=30: {"residueNumber": [1, 2, 3, 4, 5], "confidenceScore": [55, 70, 91.2, 80, 60]},
    )

    result = alphafold.run_batch(db)
    assert result == {"claimed": 1, "annotated": 1, "failed": 0}

    cur = db.conn.execute(
        "SELECT feature, feature_version, protein_accession, residue_number, score, category FROM annotations_structure WHERE variant_id = ?",
        [vid],
    )
    row = dict(cur.fetchone())
    assert row == {
        "feature": "alphafold_plddt",
        "feature_version": "AFDB_v6",
        "protein_accession": "Q12809",
        "residue_number": 3,
        "score": 91.2,
        "category": "very_high_confidence",
    }

    cur = db.conn.execute("SELECT status FROM annotation_jobs WHERE variant_id = ? AND source = 'alphafold'", [vid])
    assert cur.fetchone()["status"] == "done"
