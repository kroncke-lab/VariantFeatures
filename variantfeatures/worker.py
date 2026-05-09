"""Generic annotation queue runner.

Drains pending rows from `annotation_jobs` and dispatches each to the handler
registered for the row's `source`. Handlers raise on failure; the worker marks
the job failed (with the error message) and continues.
"""

from __future__ import annotations

import time
import traceback
from typing import Callable, Optional

from . import identity
from .handlers import HANDLERS, BATCH_HANDLERS


class WorkerSummary(dict):
    """A plain dict with helpful keys: claimed, done, failed, skipped, by_source."""


def run_pending(
    db,
    source: Optional[str] = None,
    *,
    batch_size: int = 100,
    max_jobs: Optional[int] = None,
    rate_limit_sec: Optional[float] = None,
    progress_callback: Optional[Callable[[dict], None]] = None,
) -> WorkerSummary:
    """Drain pending jobs.

    `source` filters to one source (e.g. 'clingen_ar'). If None, runs every source
    that has a registered handler.

    `max_jobs` caps total processed jobs in this run (None = unbounded).
    `rate_limit_sec` overrides the handler's default delay between calls.
    `progress_callback(event)` is called after each job with a dict describing it.
    """
    summary = WorkerSummary(claimed=0, done=0, failed=0, skipped=0, by_source={})

    while True:
        if max_jobs is not None and summary["claimed"] >= max_jobs:
            break

        remaining = (max_jobs - summary["claimed"]) if max_jobs is not None else batch_size
        limit = min(batch_size, remaining)
        jobs = db.claim_pending_jobs(source=source, limit=limit)
        if not jobs:
            break

        for job in jobs:
            summary["claimed"] += 1
            src = job["source"]
            summary["by_source"].setdefault(src, {"done": 0, "failed": 0, "skipped": 0})

            handler = HANDLERS.get(src)
            if handler is None:
                summary["skipped"] += 1
                summary["by_source"][src]["skipped"] += 1
                if src in BATCH_HANDLERS:
                    # Re-set status back to pending: batch handlers process these via their own
                    # subcommand, not through this loop.
                    db.conn.execute(
                        "UPDATE annotation_jobs SET status = 'pending', attempts = attempts - 1, last_attempted_at = NULL, updated_at = CURRENT_TIMESTAMP WHERE id = ?",
                        [job["id"]],
                    )
                    db.conn.commit()
                else:
                    db.mark_job_failed(job["id"], f"no handler registered for source {src!r}")
                if progress_callback:
                    progress_callback({"event": "skipped", "job": job, "source": src})
                continue

            delay = rate_limit_sec if rate_limit_sec is not None else getattr(handler, "DEFAULT_RATE_LIMIT_SEC", 0.0)

            try:
                handler.handle(db, job["variant_id"], job.get("payload"))
                db.mark_job_done(job["id"])
                summary["done"] += 1
                summary["by_source"][src]["done"] += 1
                if progress_callback:
                    progress_callback({"event": "done", "job": job, "source": src})
            except identity.IdentityError as e:
                db.mark_job_failed(job["id"], f"identity error: {e}")
                summary["failed"] += 1
                summary["by_source"][src]["failed"] += 1
                if progress_callback:
                    progress_callback({"event": "failed", "job": job, "source": src, "error": str(e)})
            except Exception as e:
                # Catch-all so one bad row doesn't stop the loop.
                tb = traceback.format_exception_only(type(e), e)[-1].strip()
                db.mark_job_failed(job["id"], tb[:1000])
                summary["failed"] += 1
                summary["by_source"][src]["failed"] += 1
                if progress_callback:
                    progress_callback({"event": "failed", "job": job, "source": src, "error": tb})

            if delay > 0 and len(jobs) > 1:
                time.sleep(delay)

    return summary
