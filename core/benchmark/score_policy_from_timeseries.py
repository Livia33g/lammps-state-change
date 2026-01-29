#!/usr/bin/env python3
"""
Kaggle-style scoring from the timeseries produced by:
  sim_templates/state_change/analyze_trajectory_target_yield_and_work.py

This is intentionally "policy-agnostic": it does not re-run LAMMPS; it just scores
what happened (yield vs time) and the energy-kick proxy (ΔPE) and flip counts.

Why this exists
---------------
In the proposal framing ([file://Simons_Molecular_Computing.pdf](file://Simons_Molecular_Computing.pdf)),
protocol optimization is the key challenge: compare non-equilibrium protocols by
how efficiently they drive the system toward a desired readout.

We interpret flips as energy consumption and aim to maximize:
  - yield, quickly, with minimal energy/work proxy.

Current limitations
-------------------
The CSV contains ΔPotEng between thermo outputs, which includes background dynamics.
We therefore track both:
  - work_abs_total: sum |ΔPE| over all intervals
  - work_abs_flip_intervals: sum |ΔPE| restricted to intervals with >=1 flip

If/when fixes print instantaneous per-event dE, we should switch energy accounting to Σ dE.
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional


@dataclass(frozen=True)
class Row:
    timestep: int
    pe: float
    dpe: float
    yield_val: float
    n_flips_interval: int


def _to_float(x: str) -> float:
    x = x.strip()
    if x.lower() == "nan":
        return float("nan")
    return float(x)


def load_timeseries(path: Path) -> List[Row]:
    rows: List[Row] = []
    with path.open("r", newline="") as f:
        r = csv.DictReader(f)
        required = {"timestep", "pe", "dpe", "yield", "n_statechanges_interval"}
        missing = required - set(r.fieldnames or [])
        if missing:
            raise RuntimeError(f"{path} is missing columns: {sorted(missing)}")

        for rec in r:
            rows.append(
                Row(
                    timestep=int(float(rec["timestep"])),
                    pe=_to_float(rec["pe"]),
                    dpe=_to_float(rec["dpe"]),
                    yield_val=_to_float(rec["yield"]),
                    n_flips_interval=int(float(rec["n_statechanges_interval"])),
                )
            )
    rows.sort(key=lambda rr: rr.timestep)
    return rows


def first_timestep_reaching(rows: List[Row], threshold: float) -> Optional[int]:
    for rr in rows:
        if math.isnan(rr.yield_val):
            continue
        if rr.yield_val >= threshold:
            return rr.timestep
    return None


def compute_metrics(rows: List[Row], yield_threshold: float) -> Dict[str, float]:
    # filter to rows with defined yield
    valid = [rr for rr in rows if not math.isnan(rr.yield_val)]
    if not valid:
        raise RuntimeError("No finite yield values found in CSV (did you sample dump at thermo steps?)")

    final_y = valid[-1].yield_val
    t_reach = first_timestep_reaching(valid, yield_threshold)

    n_flips_total = sum(rr.n_flips_interval for rr in valid)
    work_abs_total = sum(abs(rr.dpe) for rr in valid[1:])  # skip t0
    work_abs_flip_intervals = sum(abs(rr.dpe) for rr in valid[1:] if rr.n_flips_interval > 0)

    # Guard against divide-by-zero; if final_y==0 then efficiency is undefined/bad.
    work_per_yield = work_abs_flip_intervals / final_y if final_y > 0 else float("inf")

    # A simple scalar score:
    #  - reward final yield (dominant)
    #  - penalize energy proxy and slowness to hit threshold
    #
    # score = final_y - alpha * work_abs_flip_intervals - beta * (t_reach / t_max)
    #
    # NOTE: This is a placeholder; for a real "competition" we'd likely use a Pareto frontier
    # (or report multiple metrics) instead of collapsing to one scalar.
    t_max = valid[-1].timestep if valid[-1].timestep > 0 else 1
    t_norm = (t_reach / t_max) if t_reach is not None else 1.0

    alpha = 0.25  # weight for energy proxy (tune per project)
    beta = 0.10   # weight for slowness (tune per project)

    score = float(final_y) - alpha * float(work_abs_flip_intervals) - beta * float(t_norm)

    out: Dict[str, float] = {
        "final_yield": float(final_y),
        "yield_threshold": float(yield_threshold),
        "t_reach_threshold": float(t_reach if t_reach is not None else float("nan")),
        "t_max": float(t_max),
        "n_flips_total": float(n_flips_total),
        "work_abs_total": float(work_abs_total),
        "work_abs_flip_intervals": float(work_abs_flip_intervals),
        "work_per_yield": float(work_per_yield),
        "score": float(score),
    }
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True, help="Timeseries CSV from analyze_trajectory_target_yield_and_work.py")
    ap.add_argument("--yield-threshold", type=float, default=0.6, help="Threshold for time-to-yield metric")
    ap.add_argument("--out", default=None, help="Optional output path for a one-row CSV (leaderboard format)")
    args = ap.parse_args()

    csv_path = Path(args.csv)
    rows = load_timeseries(csv_path)
    metrics = compute_metrics(rows, args.yield_threshold)

    # Print a readable summary
    print(f"csv: {csv_path}")
    for k in (
        "final_yield",
        "t_reach_threshold",
        "n_flips_total",
        "work_abs_flip_intervals",
        "work_per_yield",
        "score",
    ):
        print(f"{k}: {metrics[k]}")

    if args.out:
        out_path = Path(args.out)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with out_path.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(metrics.keys()))
            w.writeheader()
            w.writerow(metrics)
        print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()


