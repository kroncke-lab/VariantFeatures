from __future__ import annotations

import csv
from pathlib import Path

from click.testing import CliRunner

from variantfeatures import frameshift
from variantfeatures.cli import main
from variantfeatures.database import VariantDB
from variantfeatures.transcripts import CdsSegment, Transcript


def _make_plus_strand_transcript(cds: str, genomic_start: int = 100) -> Transcript:
    n = len(cds)
    return Transcript(
        transcript_id="ENST_PLUS",
        transcript_version="1",
        gene_symbol="TESTG",
        chromosome="1",
        strand=1,
        cds_length=n,
        cds_sequence=cds,
        cds_segments=[
            CdsSegment(
                genomic_start=genomic_start,
                genomic_end=genomic_start + n - 1,
                cds_start=1,
                cds_end=n,
            )
        ],
        is_canonical=True,
        is_mane_select=True,
    )


def test_map_frameshift_position_self_to_stop_gained_snv():
    transcript = _make_plus_strand_transcript("ATGGAGTAA")

    [proxy] = frameshift.map_frameshift_position(transcript, 2, max_steps=0)

    assert proxy.direction == "self"
    assert proxy.n_steps_wrt_frameshift == 0
    assert proxy.proxy.aa_pos == 2
    assert proxy.proxy.aa_ref == "E"
    assert proxy.proxy.codon_ref == "GAG"
    assert proxy.proxy.codon_alt == "TAG"
    assert proxy.proxy.hgvs_c == "ENST_PLUS.1:c.4G>T"
    assert proxy.proxy.hgvs_p == "ENST_PLUS.1:p.Glu2Ter"


def test_map_frameshift_position_searches_left_and_right():
    transcript = _make_plus_strand_transcript("ATGGAGTAA")

    [proxy] = frameshift.map_frameshift_position(transcript, 1, max_steps=1)

    assert proxy.direction == "right"
    assert proxy.n_steps_wrt_frameshift == 1
    assert proxy.proxy.aa_pos == 2


def test_original_frameshift_variant_processing_import_shape():
    result = frameshift.frameshift_variant_processing("ATGGAGTAA", 1, 1)

    assert result["right"] == {
        "ref_aa_position": 2,
        "ref_aa": "E",
        "n_steps_wrt_ref_aa": 1,
        "near_stop_codon": "GAG",
        "lof_variant": "E2X",
    }


def test_annotate_positions_persists_proxy_mapping(tmp_path: Path):
    db = VariantDB(tmp_path / "test.db")
    transcript = _make_plus_strand_transcript("ATGGAGTAA")

    summary = frameshift.annotate_positions(
        db,
        transcript,
        [1],
        max_steps=1,
        gene_symbol="TESTG",
    )

    assert summary["positions"] == 1
    assert summary["mapped_positions"] == 1
    assert summary["mappings"] == 1

    [row] = db.get_frameshift_nonsense_mappings("TESTG", frameshift_aa_pos=1)
    assert row["direction"] == "right"
    assert row["proxy_variant_id"]
    assert row["proxy_hgvs_c"] == "ENST_PLUS.1:c.4G>T"
    assert row["proxy_hgvs_p"] == "ENST_PLUS.1:p.Glu2Ter"
    assert row["proxy_chromosome"] == "1"
    assert row["proxy_position"] == 103
    assert row["proxy_ref"] == "G"
    assert row["proxy_alt"] == "T"


def test_annotate_positions_queues_only_mapped_proxy_variants(tmp_path: Path):
    db = VariantDB(tmp_path / "test.db")
    transcript = _make_plus_strand_transcript("ATGGAGTAA")

    summary = frameshift.annotate_positions(
        db,
        transcript,
        [1],
        max_steps=1,
        enqueue_sources=("cadd",),
        gene_symbol="TESTG",
    )

    [mapping] = db.get_frameshift_nonsense_mappings("TESTG", frameshift_aa_pos=1)
    assert summary["jobs_queued"] == 1
    cur = db.conn.execute("SELECT variant_id, source FROM annotation_jobs")
    assert [dict(row) for row in cur.fetchall()] == [
        {"variant_id": mapping["proxy_variant_id"], "source": "cadd"}
    ]


def test_frameshift_map_cli_csv(monkeypatch, tmp_path: Path):
    transcript = _make_plus_strand_transcript("ATGGAGTAA")
    monkeypatch.setattr(frameshift, "fetch_canonical_transcript", lambda gene: transcript)

    result = CliRunner().invoke(
        main,
        [
            "frameshift-map",
            "-g",
            "TESTG",
            "--aa-position",
            "1",
            "--n-steps",
            "1",
            "--db",
            str(tmp_path / "test.db"),
            "--format",
            "csv",
        ],
    )

    assert result.exit_code == 0, result.output
    rows = list(csv.DictReader(result.output.splitlines()))
    assert len(rows) == 1
    assert rows[0]["frameshift_aa_pos"] == "1"
    assert rows[0]["direction"] == "right"
    assert rows[0]["proxy_hgvs_p"] == "ENST_PLUS.1:p.Glu2Ter"
