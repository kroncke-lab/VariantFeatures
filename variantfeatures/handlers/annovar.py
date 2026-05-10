"""ANNOVAR wrapper + multianno.txt importer.

Two modes:
1. **Run ANNOVAR locally**: `run_batch(db, ...)` writes a temporary AVINPUT file
   for the pending ANNOVAR jobs, calls `table_annovar.pl` from the user's
   ANNOVAR install, and parses the resulting `*_multianno.txt`.
2. **Import a pre-existing multianno.txt** (collaborator-supplied):
   `import_multianno(db, path)` parses the file and matches rows back to the
   `variants` table by (chrom, pos, ref, alt). New variants are inserted; jobs
   are marked done if they exist.

Configure via:
- `ANNOVAR_HOME` env var: path to the ANNOVAR install (containing `table_annovar.pl`).
- `ANNOVAR_DB` env var: path to the humandb folder. Default: `$ANNOVAR_HOME/humandb`.
- `--protocols` / `--operations`: which ANNOVAR databases to call.

Default protocols cover the same broad covariates we get from MyVariant.info,
plus per-transcript consequence annotation (refGene/ensGene).

ANNOVAR site: https://annovar.openbioinformatics.org/
"""

from __future__ import annotations

import csv
import os
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Iterable, Optional

from .clingen_ar import normalize_chromosome


SOURCE = "annovar"

# Default protocols / operations. Conservative subset that's typically present in a
# fresh ANNOVAR humandb install. Users with bigger dbs can pass custom values.
DEFAULT_PROTOCOLS = "refGene,gnomad411_genome,clinvar_20240611,dbnsfp47a"
DEFAULT_OPERATIONS = "g,f,f,f"
DEFAULT_BUILD = "hg38"


class HandlerError(Exception):
    pass


# ---------------------------------------------------------------------------
# Path resolution
# ---------------------------------------------------------------------------

def annovar_home() -> Optional[Path]:
    """Return the ANNOVAR install directory if discoverable, else None."""
    env = os.environ.get("ANNOVAR_HOME")
    if env:
        return Path(env)
    # Try locating table_annovar.pl on PATH and walking back up.
    found = shutil.which("table_annovar.pl")
    if found:
        return Path(found).parent
    return None


def annovar_db_dir() -> Optional[Path]:
    env = os.environ.get("ANNOVAR_DB")
    if env:
        return Path(env)
    home = annovar_home()
    if home and (home / "humandb").exists():
        return home / "humandb"
    return None


def is_installed() -> bool:
    home = annovar_home()
    return bool(home and (home / "table_annovar.pl").exists())


# ---------------------------------------------------------------------------
# Run ANNOVAR
# ---------------------------------------------------------------------------

def run_batch(
    db,
    *,
    limit: Optional[int] = None,
    build: str = DEFAULT_BUILD,
    protocols: str = DEFAULT_PROTOCOLS,
    operations: str = DEFAULT_OPERATIONS,
    work_dir: Optional[str] = None,
    keep_outputs: bool = False,
) -> dict:
    """Run ANNOVAR over pending `annovar` jobs and load the results.

    Returns: {"claimed": int, "matched": int, "missing": int, "lines_parsed": int}
    """
    if not is_installed():
        raise HandlerError(
            "ANNOVAR not detected. Set ANNOVAR_HOME to the directory containing "
            "table_annovar.pl, or place table_annovar.pl on PATH."
        )

    db_dir = annovar_db_dir()
    if db_dir is None or not db_dir.exists():
        raise HandlerError(
            "ANNOVAR humandb directory not found. Set ANNOVAR_DB or ensure "
            "$ANNOVAR_HOME/humandb exists."
        )

    jobs = db.claim_pending_jobs(source=SOURCE, limit=limit if limit is not None else 1_000_000)
    if not jobs:
        return {"claimed": 0, "matched": 0, "missing": 0, "lines_parsed": 0}

    # Pull (chrom, pos, ref, alt) for each claimed job.
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

    work = Path(work_dir) if work_dir else Path(tempfile.mkdtemp(prefix="annovar_run_"))
    avinput = work / "input.avinput"
    out_prefix = work / "out"

    _write_avinput(rows, avinput)
    multianno = _invoke_table_annovar(
        avinput=avinput, out_prefix=out_prefix, build=build, db_dir=db_dir,
        protocols=protocols, operations=operations,
    )

    # Index pending jobs by (chrom, pos, ref, alt) for matching.
    by_key: dict[tuple, list[dict]] = {}
    for r in rows:
        key = (normalize_chromosome(r["chromosome"]), int(r["position"]), r["ref"], r["alt"])
        by_key.setdefault(key, []).append(r)

    matched = 0
    lines_parsed = 0
    seen_keys: set[tuple] = set()
    for record in parse_multianno(multianno):
        lines_parsed += 1
        key = (normalize_chromosome(record["chromosome"]), record["position"], record["ref"], record["alt"])
        jobs_here = by_key.get(key)
        if not jobs_here:
            continue
        for job in jobs_here:
            _persist_record(db, job["variant_id"], record)
            db.mark_job_done(job["job_id"])
            matched += 1
        seen_keys.add(key)

    # Mark anything that wasn't in the multianno output as failed.
    missing = 0
    for key, jobs_here in by_key.items():
        if key in seen_keys:
            continue
        for job in jobs_here:
            db.mark_job_failed(job["job_id"], "no ANNOVAR multianno row for this variant")
            missing += 1

    if not keep_outputs:
        try:
            shutil.rmtree(work)
        except OSError:
            pass

    return {"claimed": len(rows), "matched": matched, "missing": missing, "lines_parsed": lines_parsed}


def _write_avinput(rows: list[dict], path: Path) -> None:
    """ANNOVAR AVINPUT format: chrom start end ref alt (tab-delimited)."""
    with open(path, "w") as f:
        for r in rows:
            chrom = normalize_chromosome(r["chromosome"])
            pos = int(r["position"])
            ref = r["ref"]
            alt = r["alt"]
            # SNV: start == end. (Indel handling is more involved; defer.)
            f.write(f"{chrom}\t{pos}\t{pos}\t{ref}\t{alt}\n")


def _invoke_table_annovar(*, avinput: Path, out_prefix: Path, build: str,
                          db_dir: Path, protocols: str, operations: str) -> Path:
    home = annovar_home()
    cmd = [
        str(home / "table_annovar.pl"),
        str(avinput),
        str(db_dir),
        "-buildver", build,
        "-out", str(out_prefix),
        "-remove",
        "-protocol", protocols,
        "-operation", operations,
        "-nastring", ".",
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise HandlerError(
            f"table_annovar.pl exited with code {proc.returncode}.\nstderr:\n{proc.stderr.strip()}"
        )
    multianno = Path(f"{out_prefix}.{build}_multianno.txt")
    if not multianno.exists():
        raise HandlerError(f"Expected ANNOVAR output not produced: {multianno}")
    return multianno


# ---------------------------------------------------------------------------
# Import standalone multianno.txt
# ---------------------------------------------------------------------------

def import_multianno(db, path: str | Path, *, source_label: str = SOURCE) -> dict:
    """Load a pre-existing multianno.txt into the DB, creating canonical variants
    and writing annotations. Useful for collaborator-supplied ANNOVAR runs.

    Returns: {"variants": int, "annotated": int}
    """
    counts = {"variants": 0, "annotated": 0}
    for record in parse_multianno(path):
        variant_id = db.upsert_variant(
            chromosome=record["chromosome"], position=record["position"],
            ref=record["ref"], alt=record["alt"],
            variant_type=_classify(record["ref"], record["alt"]),
        )
        counts["variants"] += 1
        if _persist_record(db, variant_id, record, source_label=source_label):
            counts["annotated"] += 1
    return counts


def parse_multianno(path: str | Path) -> Iterable[dict]:
    """Yield record dicts from a *_multianno.txt file."""
    p = Path(path)
    with open(p, "r", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            yield _normalize_record(row)


def _normalize_record(row: dict) -> dict:
    """Pull the columns we know how to use into a canonical dict."""
    # ANNOVAR's gene-table columns are named after the protocol used —
    # `refGene`, `refGeneWithVer`, `ensGene`, `knownGene`, etc. Pick the first
    # match by scanning the row's keys.
    func = _first_present(row, "Func.refGene", "Func.refGeneWithVer", "Func.ensGene", "Func.knownGene")
    exonic_func = _first_present(row, "ExonicFunc.refGene", "ExonicFunc.refGeneWithVer", "ExonicFunc.ensGene", "ExonicFunc.knownGene")
    aa_change = _first_present(row, "AAChange.refGene", "AAChange.refGeneWithVer", "AAChange.ensGene", "AAChange.knownGene")
    gene = _first_present(row, "Gene.refGene", "Gene.refGeneWithVer", "Gene.ensGene", "Gene.knownGene")
    # gnomAD: detect which version's columns are present and pull per-pop AFs.
    # We support v4.1 (`gnomad41_genome_*`), v4.0 (`gnomad40_*`), v2.1 (`gnomad211_*`),
    # and the legacy ALL/AF column from older bundles.
    gnomad_pops = _gnomad_per_pop_from_annovar(row)
    return {
        "chromosome": row.get("Chr") or row.get("#Chr") or row.get("chrom"),
        "position": int(row.get("Start") or row.get("Pos") or 0),
        "ref": row.get("Ref"),
        "alt": row.get("Alt"),
        "func_refgene": func,
        "exonic_func_refgene": exonic_func,
        "aa_change_refgene": aa_change,
        "gene_refgene": gene,
        # Legacy single-AF accessors (kept for backwards-compat with the old _persist_record path)
        "gnomad_af": _f(row, "AF",
                        "gnomad41_genome_AF", "gnomad40_genome_AF",
                        "gnomad411_genome_AF", "gnomAD_genome_ALL", "gnomad_genome_AF"),
        "gnomad_af_popmax": _f(row, "AF_popmax",
                                "gnomad41_genome_AF_grpmax", "gnomad40_genome_AF_popmax",
                                "gnomad411_genome_AF_popmax"),
        "gnomad_pops": gnomad_pops,  # list of {dataset, pop, af} dicts
        # ClinVar
        "clinvar_clnsig": row.get("CLNSIG") or row.get("clinvar_20240611") or row.get("CLNSIG_clinvar"),
        "clinvar_clnrevstat": row.get("CLNREVSTAT") or row.get("clinvar_20240611_REVSTAT"),
        "clinvar_id": row.get("CLNALLELEID") or row.get("CLNVID") or row.get("clinvar_VARID"),
        # dbNSFP key fields (names vary by version: dbnsfp47a / dbnsfp42a)
        "revel": _f(row, "REVEL_score", "dbnsfp_REVEL_score"),
        "cadd_phred": _f(row, "CADD_phred", "dbnsfp_CADD_phred"),
        "alphamissense": _f(row, "AlphaMissense_score", "dbnsfp_AlphaMissense_score"),
        "metasvm": _f(row, "MetaSVM_score", "dbnsfp_MetaSVM_score"),
        "metalr": _f(row, "MetaLR_score", "dbnsfp_MetaLR_score"),
        "primateai": _f(row, "PrimateAI_score", "dbnsfp_PrimateAI_score"),
        "phylop100": _f(row, "phyloP100way_vertebrate", "phyloP100way_vertebrate_score"),
        "phastcons100": _f(row, "phastCons100way_vertebrate"),
        "gerp_rs": _f(row, "GERP++_RS", "GERP_pp_RS"),
        "siphy": _f(row, "SiPhy_29way_logOdds"),
        "_raw": row,
    }


# Population codes that gnomAD ships per-pop AFs for. `_grpmax`/`_popmax` is the
# popmax frequency; we surface it as pop="popmax".
_GNOMAD_POP_CODES = ("afr", "ami", "amr", "asj", "eas", "fin", "mid", "nfe", "sas", "remaining")


def _gnomad_per_pop_from_annovar(row: dict) -> list[dict]:
    """Detect ANNOVAR's gnomAD column scheme and pull a list of {dataset, pop, af} rows."""
    candidates = [
        ("gnomad41_genome", "gnomad41_genome"),
        ("gnomad40_genome", "gnomad40_genome"),
        ("gnomad41_exome",  "gnomad41_exome"),
        ("gnomad40_exome",  "gnomad40_exome"),
        ("gnomad211_genome", "gnomad211_genome"),
        ("gnomad211_exome",  "gnomad211_exome"),
    ]
    out: list[dict] = []
    for prefix, dataset in candidates:
        # Skip if the schema isn't present in this row.
        if f"{prefix}_AF" not in row:
            continue
        # Overall AF
        all_af = _f(row, f"{prefix}_AF")
        if all_af is not None:
            out.append({"dataset": dataset, "pop": "all", "af": all_af})
        # popmax / grpmax
        popmax = _f(row, f"{prefix}_AF_grpmax", f"{prefix}_AF_popmax")
        if popmax is not None:
            out.append({"dataset": dataset, "pop": "popmax", "af": popmax})
        # Per-population AFs
        for pop in _GNOMAD_POP_CODES:
            v = _f(row, f"{prefix}_AF_{pop}")
            if v is not None:
                out.append({"dataset": dataset, "pop": pop, "af": v})
    return out


def _f(row: dict, *keys: str) -> Optional[float]:
    """First non-null/non-NA value among the candidate keys, parsed as float."""
    for k in keys:
        v = row.get(k)
        if v is None or v in (".", "NA", "", "-"):
            continue
        try:
            return float(v)
        except (ValueError, TypeError):
            continue
    return None


def _first_present(row: dict, *keys: str) -> Optional[str]:
    """Return the first key's value that's present and non-NA."""
    for k in keys:
        v = row.get(k)
        if v not in (None, ".", "NA", "", "-"):
            return v
    return None


def _classify(ref: Optional[str], alt: Optional[str]) -> str:
    if not ref or not alt:
        return "UNKNOWN"
    if len(ref) == 1 and len(alt) == 1:
        return "SNV"
    if len(ref) == len(alt):
        return "MNV"
    if len(ref) < len(alt):
        return "INS"
    return "DEL"


# ---------------------------------------------------------------------------
# Persistence
# ---------------------------------------------------------------------------

def _persist_record(db, variant_id: int, record: dict, *, source_label: str = SOURCE) -> bool:
    """Write all annotation rows we can extract from one ANNOVAR multianno record. Returns True if any rows were written."""
    wrote_any = False

    # Pathogenicity / functional predictors
    for predictor, key in [
        ("revel", "revel"),
        ("cadd_phred", "cadd_phred"),
        ("alphamissense", "alphamissense"),
        ("metasvm", "metasvm"),
        ("metalr", "metalr"),
        ("primateai", "primateai"),
    ]:
        score = record.get(key)
        if score is not None:
            db.upsert_pathogenicity(variant_id, predictor, score=score, source=source_label)
            wrote_any = True

    # Conservation
    for metric, key in [
        ("phylop100way_vertebrate", "phylop100"),
        ("phastcons100way_vertebrate", "phastcons100"),
        ("gerp_pp_rs", "gerp_rs"),
        ("siphy_29way_logodds", "siphy"),
    ]:
        score = record.get(key)
        if score is not None:
            db.upsert_conservation(variant_id, metric, score=score, source=source_label)
            wrote_any = True

    # gnomAD population: full per-pop breakdown if the gnomad41/gnomad40 schema
    # is present; fall back to the single AF / popmax columns otherwise.
    pops = record.get("gnomad_pops") or []
    if pops:
        for r in pops:
            db.upsert_population(variant_id, r["dataset"], r["pop"], af=r["af"], source=source_label)
            wrote_any = True
    else:
        af = record.get("gnomad_af")
        if af is not None:
            db.upsert_population(variant_id, "gnomad_annovar", "all", af=af, source=source_label)
            wrote_any = True
        af_popmax = record.get("gnomad_af_popmax")
        if af_popmax is not None:
            db.upsert_population(variant_id, "gnomad_annovar", "popmax", af=af_popmax, source=source_label)
            wrote_any = True

    # ClinVar
    cln_sig = record.get("clinvar_clnsig")
    if cln_sig and cln_sig not in (".", "NA", ""):
        db.upsert_clinical(
            variant_id, "clinvar",
            record_id=str(record.get("clinvar_id") or ""),
            classification=cln_sig,
            review_status=record.get("clinvar_clnrevstat"),
        )
        wrote_any = True

    # Per-transcript consequence (one row, the refGene-canonical view).
    # Prefer the more specific ExonicFunc when Func is just "exonic".
    func = record.get("func_refgene")
    exonic = record.get("exonic_func_refgene")
    if (func == "exonic" or func is None) and exonic and exonic not in (".", "NA", ""):
        func = exonic
    aa = record.get("aa_change_refgene")
    gene = record.get("gene_refgene")
    if (func or aa) and aa not in (".", None):
        # AAChange.refGene format: GENE:NM_xxx:exonN:c.Nxxxx:p.Xxxx,GENE:NM_yyy:...
        first = (aa.split(",", 1)[0] if isinstance(aa, str) else "").split(":")
        hgvs_c = next((s for s in first if s.startswith("c.")), None)
        hgvs_p = next((s for s in first if s.startswith("p.")), None)
        transcript = next((s for s in first if s.startswith("NM_") or s.startswith("NR_") or s.startswith("ENST")), "annovar")
        db.upsert_consequence(
            variant_id=variant_id,
            transcript_id=transcript,
            source=source_label,
            gene_symbol=gene,
            consequence=_func_to_so(func),
            hgvs_c=hgvs_c,
            hgvs_p=hgvs_p,
        )
        wrote_any = True

    return wrote_any


def _func_to_so(func: Optional[str]) -> Optional[str]:
    """Map ANNOVAR Func.refGene / ExonicFunc.refGene labels to SO terms (best effort)."""
    if not func:
        return None
    mapping = {
        "exonic": "exonic",
        "intronic": "intron_variant",
        "splicing": "splice_region_variant",
        "UTR3": "3_prime_UTR_variant",
        "UTR5": "5_prime_UTR_variant",
        "intergenic": "intergenic_variant",
        "ncRNA_exonic": "non_coding_transcript_exon_variant",
        "ncRNA_intronic": "non_coding_transcript_intron_variant",
        "upstream": "upstream_gene_variant",
        "downstream": "downstream_gene_variant",
        # ExonicFunc values:
        "nonsynonymous_SNV": "missense_variant",
        "synonymous_SNV": "synonymous_variant",
        "stopgain": "stop_gained",
        "stoploss": "stop_lost",
        "frameshift_deletion": "frameshift_variant",
        "frameshift_insertion": "frameshift_variant",
        "nonframeshift_deletion": "inframe_deletion",
        "nonframeshift_insertion": "inframe_insertion",
        "startloss": "start_lost",
    }
    key = func.replace(" ", "_")
    return mapping.get(key, func)
