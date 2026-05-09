"""Ensembl VEP wrapper.

Batch-style handler: writes the pending VEP jobs to a VCF, calls `vep` with
`--cache` + plugins, and parses the JSON output. One VEP invocation handles
many variants.

VEP plugin output we recognize (when present in the JSON):
- `AlphaMissense` (am_pathogenicity, am_class)
- `REVEL` (REVEL)
- `CADD` (CADD_PHRED, CADD_RAW)
- `SpliceAI` (donor_gain/loss, acceptor_gain/loss + distances)
- `dbNSFP` (configurable subset; we look for common fields)
- `pLI`, `LoF`, `LoF_filter`, `LoF_flags` (LOFTEE / pLI plugins)

Configuration (env vars override CLI flags):
- `VEP_BIN`: path to the `vep` binary. Default: discovered on PATH.
- `VEP_CACHE_DIR`: cache directory. Default: `~/.vep`.
- `VEP_SPECIES`: default `homo_sapiens`.
- `VEP_ASSEMBLY`: default `GRCh38`.
- `VEP_PLUGINS`: comma-separated plugin invocations to pass through, e.g.
  `AlphaMissense,file=/path/AlphaMissense_hg38.tsv.gz;REVEL,file=...`. Each
  plugin string here is appended to a `--plugin` flag verbatim (semicolons
  separate plugin args within one plugin per VEP convention).

Install VEP via:
    brew install ensembl-vep
or
    cpanm --quiet --notest --installdeps Bio::EnsEMBL::VEP
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Any, Iterable, Optional

from .clingen_ar import normalize_chromosome


SOURCE = "vep"
DEFAULT_RATE_LIMIT_SEC = 0.0  # batch mode; per-job rate limit not relevant


class HandlerError(Exception):
    pass


# ---------------------------------------------------------------------------
# Config helpers
# ---------------------------------------------------------------------------

def vep_bin() -> Optional[str]:
    return os.environ.get("VEP_BIN") or shutil.which("vep")


def vep_cache_dir() -> Path:
    env = os.environ.get("VEP_CACHE_DIR")
    if env:
        return Path(env)
    return Path.home() / ".vep"


def is_installed() -> bool:
    return bool(vep_bin())


def configured_plugins() -> list[str]:
    raw = os.environ.get("VEP_PLUGINS", "")
    return [p for p in (s.strip() for s in raw.split(",")) if p]


# ---------------------------------------------------------------------------
# Batch run
# ---------------------------------------------------------------------------

def run_batch(
    db,
    *,
    limit: Optional[int] = None,
    species: Optional[str] = None,
    assembly: Optional[str] = None,
    plugins: Optional[list[str]] = None,
    extra_args: Optional[list[str]] = None,
    work_dir: Optional[str] = None,
    keep_outputs: bool = False,
) -> dict:
    """Drain pending vep jobs in one VEP invocation.

    Returns: {"claimed": int, "annotated": int, "missing": int}
    """
    binary = vep_bin()
    if not binary:
        raise HandlerError(
            "VEP binary not found. Install via Homebrew (`brew install ensembl-vep`) "
            "or set VEP_BIN to the absolute path of `vep`."
        )

    cache = vep_cache_dir()
    if not cache.exists():
        raise HandlerError(
            f"VEP cache directory {cache} does not exist. Set VEP_CACHE_DIR or run "
            f"`vep_install --AUTO cf --SPECIES homo_sapiens --ASSEMBLY GRCh38`."
        )

    species = species or os.environ.get("VEP_SPECIES", "homo_sapiens")
    assembly = assembly or os.environ.get("VEP_ASSEMBLY", "GRCh38")
    plugins = list(plugins or configured_plugins())
    extra_args = list(extra_args or [])

    jobs = db.claim_pending_jobs(source=SOURCE, limit=limit if limit is not None else 1_000_000)
    if not jobs:
        return {"claimed": 0, "annotated": 0, "missing": 0}

    ids = [j["id"] for j in jobs]
    placeholders = ",".join("?" * len(ids))
    cur = db.conn.execute(
        f"""
        SELECT j.id AS job_id, j.variant_id AS variant_id,
               v.chromosome AS chromosome, v.position AS position, v.ref AS ref, v.alt AS alt
        FROM annotation_jobs j JOIN variants v ON v.id = j.variant_id
        WHERE j.id IN ({placeholders})
        """,
        ids,
    )
    rows = [dict(r) for r in cur.fetchall()]

    work = Path(work_dir) if work_dir else Path(tempfile.mkdtemp(prefix="vep_run_"))
    vcf_path = work / "input.vcf"
    json_path = work / "output.json"

    _write_vcf(rows, vcf_path)
    _invoke_vep(binary, vcf_path, json_path, cache=cache, species=species, assembly=assembly,
                plugins=plugins, extra_args=extra_args)

    by_key: dict[tuple, list[dict]] = {}
    for r in rows:
        key = (normalize_chromosome(r["chromosome"]), int(r["position"]), r["ref"], r["alt"])
        by_key.setdefault(key, []).append(r)

    annotated = 0
    seen_keys: set[tuple] = set()
    with open(json_path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            record = json.loads(line)
            key = _key_from_vep(record)
            if key is None:
                continue
            jobs_here = by_key.get(key)
            if not jobs_here:
                continue
            for job in jobs_here:
                _persist_record(db, job["variant_id"], record)
                db.mark_job_done(job["id"])
                annotated += 1
            seen_keys.add(key)

    missing = 0
    for key, jobs_here in by_key.items():
        if key in seen_keys:
            continue
        for job in jobs_here:
            db.mark_job_failed(job["id"], "no VEP record returned for this variant")
            missing += 1

    if not keep_outputs:
        shutil.rmtree(work, ignore_errors=True)

    return {"claimed": len(rows), "annotated": annotated, "missing": missing}


def _write_vcf(rows: list[dict], path: Path) -> None:
    """Minimal VCF for VEP input. VEP only needs CHROM/POS/ID/REF/ALT."""
    with open(path, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for r in rows:
            chrom = normalize_chromosome(r["chromosome"])
            f.write(f"{chrom}\t{r['position']}\t.\t{r['ref']}\t{r['alt']}\t.\tPASS\t.\n")


def _invoke_vep(binary: str, vcf_path: Path, out_path: Path, *,
                cache: Path, species: str, assembly: str,
                plugins: list[str], extra_args: list[str]) -> None:
    cmd = [
        binary,
        "-i", str(vcf_path),
        "-o", str(out_path),
        "--cache",
        "--dir_cache", str(cache),
        "--species", species,
        "--assembly", assembly,
        "--json",
        "--no_stats",
        "--offline",
        "--force_overwrite",
        "--canonical",
        "--mane",
    ]
    for plugin in plugins:
        cmd += ["--plugin", plugin]
    cmd += extra_args

    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise HandlerError(
            f"vep exited with code {proc.returncode}.\nstderr:\n{proc.stderr.strip()}"
        )


def _key_from_vep(record: dict) -> Optional[tuple]:
    """Recover (chrom, pos, ref, alt) from a VEP JSON record."""
    if not record:
        return None
    chrom = record.get("seq_region_name") or record.get("chrom")
    pos = record.get("start") or record.get("pos")
    alleles = record.get("allele_string") or ""
    if not (chrom and pos and "/" in alleles):
        return None
    ref, alt = alleles.split("/", 1)
    return (normalize_chromosome(chrom), int(pos), ref, alt)


# ---------------------------------------------------------------------------
# Persistence: parse one VEP JSON record into our annotation tables
# ---------------------------------------------------------------------------

def _persist_record(db, variant_id: int, record: dict) -> None:
    transcript_csqs = record.get("transcript_consequences") or []
    if not isinstance(transcript_csqs, list):
        transcript_csqs = [transcript_csqs] if transcript_csqs else []

    # Pick the canonical / MANE Select consequence (preferred) or first transcript.
    primary = _pick_primary_csq(transcript_csqs)

    if primary:
        consequence_terms = primary.get("consequence_terms") or []
        consequence = consequence_terms[0] if consequence_terms else None
        db.upsert_consequence(
            variant_id=variant_id,
            transcript_id=str(primary.get("transcript_id") or ""),
            source=SOURCE,
            gene_symbol=primary.get("gene_symbol"),
            gene_ensembl=primary.get("gene_id"),
            consequence=consequence,
            hgvs_c=primary.get("hgvsc"),
            hgvs_p=primary.get("hgvsp"),
            aa_pos=primary.get("protein_start"),
            is_canonical=int(bool(primary.get("canonical"))),
            is_mane_select=int(bool(primary.get("mane_select"))),
            is_mane_plus_clinical=int(bool(primary.get("mane_plus_clinical"))),
        )

    # Pull plugin-emitted scores from the (per-transcript) consequence dicts. We collapse
    # across transcripts by taking the max numeric value (consistent with our earlier
    # convention for AlphaMissense list-shaped scores from MyVariant.info).
    for predictor, key in _VEP_PLUGIN_PATHOGENICITY:
        values = [c.get(key) for c in transcript_csqs if isinstance(c, dict)]
        score = _max_numeric(values)
        if score is not None:
            db.upsert_pathogenicity(variant_id, predictor, score=score, source=SOURCE)

    # AlphaMissense plugin emits am_pathogenicity (score) + am_class (category).
    am_scores = [c.get("am_pathogenicity") for c in transcript_csqs if isinstance(c, dict)]
    am_classes = [c.get("am_class") for c in transcript_csqs if isinstance(c, dict)]
    am_score = _max_numeric(am_scores)
    if am_score is not None:
        am_class = next((c for c in am_classes if isinstance(c, str)), None)
        category = {"likely_benign": "B", "ambiguous": "A", "likely_pathogenic": "P"}.get(am_class, am_class)
        db.upsert_pathogenicity(variant_id, "alphamissense", score=am_score, category=category, source=SOURCE)

    # SpliceAI plugin emits a delta-score block per transcript. Take the first
    # transcript that has any SpliceAI score and write all 4 directional scores
    # from that one consequence dict.
    for csq in transcript_csqs:
        if not isinstance(csq, dict):
            continue
        spliceai_keys = (
            ("acceptor_gain", "spliceai_pred_ds_ag"),
            ("acceptor_loss", "spliceai_pred_ds_al"),
            ("donor_gain", "spliceai_pred_ds_dg"),
            ("donor_loss", "spliceai_pred_ds_dl"),
        )
        if not any(csq.get(k) is not None for _, k in spliceai_keys):
            continue
        for score_type, key in spliceai_keys:
            score = csq.get(key)
            if score is None:
                continue
            try:
                score_val = float(score)
            except (TypeError, ValueError):
                continue
            distance = csq.get(key.replace("_ds_", "_dp_"))
            try:
                distance_val = int(distance) if distance is not None else None
            except (TypeError, ValueError):
                distance_val = None
            db.upsert_splice(
                variant_id, "spliceai", score_type=score_type,
                score=score_val, distance=distance_val, source=SOURCE,
            )
        break  # one transcript's worth of SpliceAI scores is enough per variant

    # Existing variation IDs (rsIDs etc.)
    for ev in record.get("colocated_variants") or []:
        if not isinstance(ev, dict):
            continue
        rsid = ev.get("id")
        if rsid and str(rsid).startswith("rs"):
            db.add_aliases(variant_id, [{"alias_type": "rsid", "alias_value": rsid, "source": SOURCE}])


_VEP_PLUGIN_PATHOGENICITY: list[tuple[str, str]] = [
    ("revel", "revel_score"),
    ("cadd_phred", "cadd_phred"),
    ("cadd_raw", "cadd_raw"),
    ("primateai", "primateai_score"),
    ("metasvm", "metasvm_score"),
    ("metalr", "metalr_score"),
    ("mcap", "mcap_score"),
    ("clinpred", "clinpred_score"),
    ("eve", "eve_score"),
    ("esm1b", "esm1b_score"),
    ("polyphen2_hvar", "polyphen2_hvar_score"),
    ("sift", "sift_score"),
    ("provean", "provean_score"),
    ("mutpred", "mutpred_score"),
    ("vest4", "vest4_score"),
    # LOFTEE / constraint
    ("loftee_lof", "lof"),  # 'HC' / 'LC' string -> won't pass numeric coercion, skipped
]


def _pick_primary_csq(csqs: Iterable[dict]) -> Optional[dict]:
    csqs = [c for c in csqs if isinstance(c, dict)]
    for c in csqs:
        if c.get("mane_select"):
            return c
    for c in csqs:
        if c.get("canonical"):
            return c
    return csqs[0] if csqs else None


def _max_numeric(values) -> Optional[float]:
    nums: list[float] = []
    for v in values or []:
        if isinstance(v, bool):
            continue
        if isinstance(v, (int, float)):
            nums.append(float(v))
        elif isinstance(v, str):
            try:
                nums.append(float(v))
            except ValueError:
                continue
    return max(nums) if nums else None
