from __future__ import annotations

import gzip
from pathlib import Path

from variantfeatures.database import VariantDB
from variantfeatures.handlers import clinvar


def test_import_gene_matches_clinvar_summary_by_coordinate(tmp_path: Path):
    db = VariantDB(tmp_path / "test.db")
    variant_id = db.upsert_variant(
        chromosome="7",
        position=150977910,
        ref="G",
        alt="T",
        variant_type="SNV",
    )
    db.upsert_consequence(
        variant_id,
        "NM_000238.4",
        "enumerated",
        gene_symbol="KCNH2",
        consequence="missense_variant",
    )
    other_id = db.upsert_variant(chromosome="7", position=150977901, ref="T", alt="A")
    db.upsert_consequence(
        other_id,
        "NM_OTHER.1",
        "enumerated",
        gene_symbol="OTHER",
        consequence="missense_variant",
    )

    path = tmp_path / "variant_summary.txt.gz"
    _write_variant_summary(path)

    result = clinvar.import_gene(db, "KCNH2", data_file=path)

    assert result == {"rows": 1, "matched_rows": 1, "annotations": 1, "aliases": 2}
    row = db.conn.execute(
        """
        SELECT source, record_id, classification, review_status, stars, last_evaluated
        FROM annotations_clinical
        WHERE variant_id = ?
        """,
        [variant_id],
    ).fetchone()
    assert dict(row) == {
        "source": "clinvar",
        "record_id": "VCV000012345",
        "classification": "Likely pathogenic",
        "review_status": "criteria provided, single submitter",
        "stars": 1,
        "last_evaluated": "2024-01-02",
    }
    aliases = {(row["alias_type"], row["alias_value"]) for row in db.get_aliases(variant_id)}
    assert ("clinvar_vcv", "VCV000012345") in aliases
    assert ("clinvar_allele", "111") in aliases


def _write_variant_summary(path: Path) -> None:
    header = [f"col{i}" for i in range(34)]
    row = [""] * 34
    row[0] = "111"
    row[1] = "single nucleotide variant"
    row[2] = "NM_000238.4(KCNH2):c.4G>T (p.Pro2Thr)"
    row[4] = "KCNH2"
    row[6] = "Likely pathogenic"
    row[8] = "Jan 02, 2024"
    row[16] = "GRCh38"
    row[18] = "7"
    row[19] = "150977910"
    row[24] = "criteria provided, single submitter"
    row[30] = "12345"
    row[31] = "150977910"
    row[32] = "G"
    row[33] = "T"
    with gzip.open(path, "wt") as fh:
        fh.write("\t".join(header) + "\n")
        fh.write("\t".join(row) + "\n")
