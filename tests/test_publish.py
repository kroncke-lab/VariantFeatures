from __future__ import annotations

import json
import sqlite3
from pathlib import Path
from types import SimpleNamespace

from click.testing import CliRunner

from variantfeatures import publish
from variantfeatures.cli import main
from variantfeatures.database import VariantDB
from variantfeatures.publish import build_manifest, export_gene_slice, upload


def _seed_publish_db(db: VariantDB) -> None:
    db.conn.execute(
        """
        INSERT INTO transcripts
            (transcript_id, gene_symbol, gene_ensembl, refseq_match, protein_id,
             biotype, cds_length, is_canonical, is_mane_select, source)
        VALUES
            ('ENST_KCNH2.1', 'KCNH2', 'ENSG_KCNH2', 'NM_KCNH2.1', 'P_KCNH2',
             'protein_coding', 1200, 1, 1, 'test'),
            ('ENST_BRCA1.1', 'BRCA1', 'ENSG_BRCA1', 'NM_BRCA1.1', 'P_BRCA1',
             'protein_coding', 1800, 1, 1, 'test')
        """
    )
    db.conn.commit()

    k1 = db.upsert_variant(
        chromosome="7",
        position=150000001,
        ref="A",
        alt="G",
        ca_id="CA_K1",
        variant_type="SNV",
        hgvs_g="NC_000007.14:g.150000001A>G",
    )
    db.add_aliases(k1, [{"alias_type": "rsid", "alias_value": "rsK1"}], source="test")
    db.upsert_consequence(
        k1,
        "ENST_KCNH2.1",
        "enumerated",
        gene_symbol="KCNH2",
        gene_ensembl="ENSG_KCNH2",
        consequence="missense_variant",
        hgvs_c="c.1A>G",
        hgvs_p="p.Lys1Arg",
        aa_pos=1,
        is_canonical=1,
        is_mane_select=1,
    )
    db.upsert_pathogenicity(k1, "revel", predictor_version="1.3", score=0.91, source="test")
    db.upsert_population(k1, "gnomad_exomes_v4", "all", af=1e-5, source="test")
    db.upsert_clinical(k1, "clinvar", "VCV1", classification="Pathogenic")
    db.upsert_splice(k1, "spliceai", "overall", score=0.1, source="test")
    db.upsert_expression(k1, "pext", dataset="test", tissue="heart", score=0.8, source="test")
    db.upsert_structure(k1, "alphafold_plddt", feature_version="AFDB_v6", score=90.0)
    db.upsert_conservation(k1, "phylop", score=2.2, source="test")
    db.store_penetrance(
        "missense",
        k1,
        "KCNH2",
        "p.Lys1Arg",
        0.4,
        0.35,
        0.2,
        0.6,
        model_version="test",
        n_carriers=5,
    )

    k2 = db.upsert_variant(chromosome="7", position=150000002, ref="C", alt="T")
    db.upsert_consequence(
        k2,
        "ENST_KCNH2.1",
        "annovar",
        gene_symbol="KCNH2",
        consequence="synonymous_variant",
        hgvs_c="c.2C>T",
        hgvs_p="p.Gly1=",
    )
    db.upsert_pathogenicity(k2, "cadd_phred", predictor_version="GRCh38-v1.6", score=12.3)

    b1 = db.upsert_variant(chromosome="17", position=43000001, ref="G", alt="A")
    db.upsert_consequence(
        b1,
        "ENST_BRCA1.1",
        "enumerated",
        gene_symbol="BRCA1",
        gene_ensembl="ENSG_BRCA1",
        consequence="missense_variant",
        hgvs_c="c.1G>A",
        hgvs_p="p.Val1Ile",
    )
    db.upsert_pathogenicity(b1, "revel", predictor_version="1.3", score=0.12, source="test")
    db.upsert_gene_constraint("KCNH2", dataset="gnomad_v4", pli=0.99, source="test")
    db.upsert_gene_constraint("BRCA1", dataset="gnomad_v4", pli=1.0, source="test")
    db.enqueue_job(k1, "cadd")


def _count(conn: sqlite3.Connection, table: str) -> int:
    return int(conn.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0])


def test_export_gene_slice_copies_only_one_gene(tmp_path: Path):
    db_path = tmp_path / "source.db"
    db = VariantDB(db_path)
    _seed_publish_db(db)
    db.close()

    out_path = tmp_path / "KCNH2.db"
    summary = export_gene_slice(db_path, "KCNH2", out_path)

    assert summary["row_counts"]["variants"] == 2
    assert summary["row_counts"]["variant_consequences"] == 2
    assert summary["row_counts"]["annotations_pathogenicity"] == 2

    conn = sqlite3.connect(out_path)
    conn.row_factory = sqlite3.Row
    try:
        assert _count(conn, "genes") == 1
        assert _count(conn, "gene_constraint") == 1
        assert _count(conn, "transcripts") == 1
        assert _count(conn, "variants") == 2
        assert _count(conn, "variant_consequences") == 2
        assert _count(conn, "variant_aliases") == 1
        assert _count(conn, "annotations_population") == 1
        assert _count(conn, "annotations_clinical") == 1
        assert _count(conn, "annotations_splice") == 1
        assert _count(conn, "annotations_expression") == 1
        assert _count(conn, "annotations_structure") == 1
        assert _count(conn, "annotations_conservation") == 1
        assert _count(conn, "penetrance_estimates") == 1
        assert _count(conn, "annotation_jobs") == 0
        genes = {
            row["gene_symbol"]
            for row in conn.execute("SELECT DISTINCT gene_symbol FROM variant_consequences")
        }
        assert genes == {"KCNH2"}
    finally:
        conn.close()


def test_build_manifest_includes_counts_and_artifact_hashes(monkeypatch, tmp_path: Path):
    db_path = tmp_path / "source.db"
    db = VariantDB(db_path)
    _seed_publish_db(db)
    db.close()
    slice_path = tmp_path / "KCNH2.db"
    export_gene_slice(db_path, "KCNH2", slice_path)
    monkeypatch.setattr(publish, "get_git_sha", lambda: "abcdef123456")

    manifest = build_manifest(
        db_path,
        ["KCNH2"],
        "2026-06-13T12:34:56Z",
        {"KCNH2": slice_path},
    )

    assert manifest["producer"] == "variantfeatures"
    assert manifest["git_sha"] == "abcdef123456"
    assert manifest["schema_version"] is None
    assert manifest["genes"]["KCNH2"]["row_counts"]["variants"] == 2
    assert manifest["source_versions"]["pathogenicity"]
    artifact = manifest["artifacts"]["genes/KCNH2.db"]
    assert artifact["sha256"] == publish.sha256_file(slice_path)
    assert artifact["bytes"] == slice_path.stat().st_size


class FakeResourceNotFoundError(Exception):
    pass


class FakeResourceExistsError(Exception):
    pass


class FakeResourceModifiedError(Exception):
    pass


class FakeContentSettings:
    def __init__(self, content_type: str):
        self.content_type = content_type


class FakeCredential:
    pass


class FakeBlobServiceClient:
    store: dict[str, dict] = {}
    calls: list[dict] = []
    instances: list["FakeBlobServiceClient"] = []

    def __init__(self, account_url: str, credential):
        self.account_url = account_url
        self.credential = credential
        self.instances.append(self)

    def get_container_client(self, container: str):
        return FakeContainerClient(container)


class FakeContainerClient:
    def __init__(self, container: str):
        self.container = container

    def get_blob_client(self, blob: str):
        return FakeBlobClient(blob)


class FakeDownloader:
    def __init__(self, data: bytes):
        self.data = data

    def readall(self) -> bytes:
        return self.data


class FakeBlobClient:
    def __init__(self, blob: str):
        self.blob = blob

    def get_blob_properties(self):
        if self.blob not in FakeBlobServiceClient.store:
            raise FakeResourceNotFoundError(self.blob)
        record = FakeBlobServiceClient.store[self.blob]
        return SimpleNamespace(metadata=record.get("metadata") or {}, etag=record["etag"])

    def download_blob(self):
        if self.blob not in FakeBlobServiceClient.store:
            raise FakeResourceNotFoundError(self.blob)
        return FakeDownloader(FakeBlobServiceClient.store[self.blob]["data"])

    def upload_blob(self, data, **kwargs):
        exists = self.blob in FakeBlobServiceClient.store
        if not kwargs.get("overwrite") and exists:
            raise FakeResourceExistsError(self.blob)
        etag = kwargs.get("etag")
        if etag and (not exists or FakeBlobServiceClient.store[self.blob]["etag"] != etag):
            raise FakeResourceModifiedError(self.blob)
        payload = data.read() if hasattr(data, "read") else data
        version = len(FakeBlobServiceClient.calls) + 1
        FakeBlobServiceClient.store[self.blob] = {
            "data": payload,
            "metadata": kwargs.get("metadata") or {},
            "etag": f"etag-{version}",
        }
        FakeBlobServiceClient.calls.append({"blob": self.blob, "kwargs": kwargs})


def _install_fake_azure(monkeypatch):
    FakeBlobServiceClient.store = {}
    FakeBlobServiceClient.calls = []
    FakeBlobServiceClient.instances = []
    monkeypatch.setattr(publish, "BlobServiceClient", FakeBlobServiceClient)
    monkeypatch.setattr(publish, "DefaultAzureCredential", FakeCredential)
    monkeypatch.setattr(publish, "ContentSettings", FakeContentSettings)
    monkeypatch.setattr(
        publish,
        "MatchConditions",
        SimpleNamespace(IfNotModified="IfNotModified"),
    )
    monkeypatch.setattr(publish, "ResourceNotFoundError", FakeResourceNotFoundError)
    monkeypatch.setattr(publish, "ResourceExistsError", FakeResourceExistsError)
    monkeypatch.setattr(publish, "ResourceModifiedError", FakeResourceModifiedError)


def test_upload_uses_versioned_layout_and_merges_latest(monkeypatch, tmp_path: Path):
    _install_fake_azure(monkeypatch)
    slice_path = tmp_path / "KCNH2.db"
    slice_path.write_bytes(b"sqlite slice")
    manifest = {
        "producer": "variantfeatures",
        "git_sha": "abcdef123456",
        "built_at": "2026-06-13T12:34:56Z",
        "schema_version": "test-schema",
        "genes": {"KCNH2": {"row_counts": {"variants": 1}}},
        "artifacts": {},
    }
    FakeBlobServiceClient.store["variantfeatures/latest.json"] = {
        "data": json.dumps({"KCNQ1": {"path": "old/path.db"}}).encode(),
        "metadata": {},
        "etag": "etag-existing",
    }

    result = upload(
        "https://acct.blob.core.windows.net",
        "variantfeatures",
        "variantfeatures",
        "2026-06-13T12:34:56Z",
        "abcdef123456",
        {"KCNH2": slice_path},
        manifest,
        dry_run=False,
    )

    expected = "variantfeatures/20260613-1234__abcdef1/genes/KCNH2.db"
    assert result["artifacts"][0]["blob_path"] == expected
    assert result["manifest"]["blob_path"] == (
        "variantfeatures/20260613-1234__abcdef1/manifest.json"
    )
    latest = json.loads(FakeBlobServiceClient.store["variantfeatures/latest.json"]["data"])
    assert latest["KCNQ1"]["path"] == "old/path.db"
    assert latest["KCNH2"] == {
        "path": expected,
        "sha256": publish.sha256_file(slice_path),
        "built_at": "2026-06-13T12:34:56Z",
        "schema_version": "test-schema",
    }
    latest_call = FakeBlobServiceClient.calls[-1]
    assert latest_call["blob"] == "variantfeatures/latest.json"
    assert latest_call["kwargs"]["etag"] == "etag-existing"
    assert latest_call["kwargs"]["match_condition"] == "IfNotModified"


def test_publish_cli_dry_run_writes_local_artifacts_without_azure(
    monkeypatch,
    tmp_path: Path,
):
    _install_fake_azure(monkeypatch)
    monkeypatch.setattr(publish, "get_git_sha", lambda: "abcdef123456")
    runner = CliRunner()

    with runner.isolated_filesystem(temp_dir=tmp_path):
        db_path = Path("source.db")
        db = VariantDB(db_path)
        _seed_publish_db(db)
        db.close()

        result = runner.invoke(
            main,
            [
                "publish",
                "--gene",
                "KCNH2",
                "--db",
                str(db_path),
                "--built-at",
                "2026-06-13T12:34:56Z",
                "--dry-run",
            ],
            env={"VF_BLOB_ACCOUNT_URL": "", "VF_BLOB_CONTAINER": ""},
        )

        assert result.exit_code == 0, result.output
        release_dir = Path("dist/publish/20260613-1234__abcdef1")
        assert (release_dir / "genes/KCNH2.db").exists()
        assert (release_dir / "manifest.json").exists()
        assert (
            "WOULD upload: variantfeatures/20260613-1234__abcdef1/genes/KCNH2.db"
            in result.output
        )
        assert "latest.json diff:" in result.output
        assert FakeBlobServiceClient.instances == []
