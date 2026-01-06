# back/ (local archive)

This directory is a **local-only archive** to keep the working tree clean.

We move old runtime artifacts here, e.g.:
- Slurm logs (`slurm_*.out`, `slurm_*.err`)
- One-off output logs from previous debugging runs

Notes:
- `back/` is ignored by git (`back/**`) so it will **not** be pushed.
- Keep current job logs in the main directory if a job is actively running.


