#!/usr/bin/env python3
"""
Run a single benchmark task described by a JSON file.

This script exists to make benchmarking many different statechange simulations easy:
  - you write a new task JSON (paths + target definition)
  - you run this script
  - it produces:
      1) a standardized timeseries CSV/PNG (yield + ΔPE + flips)
      2) a one-row leaderboard CSV (metrics + score)

It intentionally calls the existing analysis/scoring scripts as subprocesses so it works
in any activated environment (e.g. your jax env for plotting).
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional


def _req(d: Dict[str, Any], key: str) -> Any:
    if key not in d:
        raise RuntimeError(f"Task is missing required key: '{key}'")
    return d[key]


def _resolve_path(base: Path, p: Optional[str]) -> Optional[Path]:
    if p is None:
        return None
    pp = Path(p)
    return pp if pp.is_absolute() else (base / pp).resolve()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--task", required=True, help="Path to task JSON (see task_schema.md)")
    ap.add_argument(
        "--dry-run",
        action="store_true",
        help="Print commands but do not execute them",
    )
    args = ap.parse_args()

    task_path = Path(args.task).resolve()
    task_base = task_path.parent
    task = json.loads(task_path.read_text())

    name = str(_req(task, "name"))
    dump = _resolve_path(task_base, str(_req(task, "dump")))
    thermo = _resolve_path(task_base, str(_req(task, "thermo")))
    events = _resolve_path(task_base, task.get("events"))

    site_types = _req(task, "site_types")
    if not isinstance(site_types, list) or not site_types:
        raise RuntimeError("'site_types' must be a non-empty list of ints")
    site_types_str = [str(int(x)) for x in site_types]

    bond_cutoff = float(_req(task, "bond_cutoff"))
    yield_mode = str(task.get("yield_mode", "fraction_molecules"))
    label_mode = str(task.get("label_mode", "majority_site_type"))
    sample = str(task.get("sample", "thermo"))
    max_frames = task.get("max_frames", None)

    out_dir = _resolve_path(task_base, task.get("out_dir", "../analysis/"))
    if out_dir is None:
        raise RuntimeError("out_dir resolved to None")
    out_dir.mkdir(parents=True, exist_ok=True)

    # Where outputs go
    out_prefix = (out_dir / name).resolve()
    out_csv = out_prefix.with_suffix(".csv")
    out_png = out_prefix.with_suffix(".png")

    # Script locations (relative to this file)
    bench_dir = Path(__file__).resolve().parent
    state_change_dir = bench_dir.parent
    analyze_py = (state_change_dir / "analyze_trajectory_target_yield_and_work.py").resolve()
    score_py = (bench_dir / "score_policy_from_timeseries.py").resolve()

    # Build analysis command
    cmd: List[str] = [
        sys.executable,
        str(analyze_py),
        "--dump",
        str(dump),
        "--thermo",
        str(thermo),
        "--site-types",
        *site_types_str,
        "--bond-cutoff",
        str(bond_cutoff),
        "--label-mode",
        label_mode,
        "--sample",
        sample,
        "--out",
        str(out_prefix),
        "--yield-mode",
        yield_mode,
    ]
    if events is not None and events.exists():
        cmd += ["--events", str(events)]

    if yield_mode == "species_fraction":
        species_label = task.get("species_label", None)
        if species_label is None:
            raise RuntimeError("yield_mode=species_fraction requires 'species_label'")
        cmd += ["--species-label", str(int(species_label))]
    else:
        target_size = task.get("target_size", None)
        if target_size is None:
            raise RuntimeError("Cluster yield modes require 'target_size'")
        cmd += ["--target-size", str(int(target_size))]
        if task.get("target_composition") is not None:
            cmd += ["--target-composition", str(task["target_composition"])]

    if max_frames is not None:
        cmd += ["--max-frames", str(int(max_frames))]

    # Scoring command
    yield_threshold = float(task.get("yield_threshold", 0.6))
    leaderboard_row = (out_dir / f"{name}.leaderboard.csv").resolve()
    cmd_score = [
        sys.executable,
        str(score_py),
        "--csv",
        str(out_csv),
        "--yield-threshold",
        str(yield_threshold),
        "--out",
        str(leaderboard_row),
    ]

    print(f"=== Task: {name} ===")
    print(f"- dump:   {dump}")
    print(f"- thermo: {thermo}")
    if events is not None:
        print(f"- events: {events} ({'exists' if events.exists() else 'missing'})")
    print(f"- out_csv: {out_csv}")
    print(f"- out_png: {out_png}")
    print(f"- leaderboard: {leaderboard_row}")

    print("\n--- Running analysis ---")
    print(" ".join(cmd))
    if not args.dry_run:
        subprocess.run(cmd, check=True)

    print("\n--- Scoring ---")
    print(" ".join(cmd_score))
    if not args.dry_run:
        subprocess.run(cmd_score, check=True)

    print("\n✅ Done.")


if __name__ == "__main__":
    main()


