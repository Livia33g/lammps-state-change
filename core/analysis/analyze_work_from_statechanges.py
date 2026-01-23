#!/usr/bin/env python3
"""
Compute "work" injected by state-change events from LAMMPS logs.

Definition used:
- Each state-change event occurs at timestep t_event (parsed from stderr lines like:
    "STATECHANGE ...: step 3500 ...")
- Let PE(t) be the potential energy from the thermo table.
- Define event work as:
    W_event = PE(t_event) - PE(t_prev)
  where t_prev is the most recent thermo timestep strictly less than t_event.

Outputs:
- CSV with per-event work and cumulative work
- Optional PNG plot (if matplotlib is available)

This avoids reading huge trajectory dumps; it uses the thermo and statechange logs.
"""

from __future__ import annotations

import argparse
import bisect
import csv
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple


STATECHANGE_RE = re.compile(r"STATECHANGE\s+\S+:\s+step\s+(\d+)\s+(.+)$")
DENERGY_RE = re.compile(r"\bdE\s+([+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)\b")


def parse_statechange_steps(path: Path) -> List[int]:
    steps: List[int] = []
    for line in path.read_text(errors="ignore").splitlines():
        m = STATECHANGE_RE.search(line)
        if not m:
            continue
        steps.append(int(m.group(1)))
    steps.sort()
    return steps


def parse_statechange_step_dE(path: Path) -> List[Tuple[int, float]]:
    """
    Parse instantaneous energy change per flip event if the fix prints it, e.g.:
      STATECHANGE ...: step 3500 ... dE -0.0123
    Returns list of (timestep, dE) rows (one per printed flip line).
    """
    out: List[Tuple[int, float]] = []
    for line in path.read_text(errors="ignore").splitlines():
        m = STATECHANGE_RE.search(line)
        if not m:
            continue
        step = int(m.group(1))
        md = DENERGY_RE.search(line)
        if not md:
            continue
        out.append((step, float(md.group(1))))
    return out


@dataclass(frozen=True)
class ThermoSeries:
    steps: List[int]
    pe: List[float]


def _try_parse_thermo_table(lines: List[str]) -> Optional[ThermoSeries]:
    """
    Parse a LAMMPS thermo table of the form:
      Step Temp PotEng KinEng TotEng
      0  ...
      1000 ...
    Works for both old/new LAMMPS output as long as there is a header with Step and a PE column.
    """
    header_idx = None
    header_cols: List[str] = []
    for i, line in enumerate(lines):
        if "Step" in line and ("PotEng" in line or "pe" in line.lower()):
            # Use the first matching header
            header_idx = i
            header_cols = line.split()
            break
    if header_idx is None:
        return None

    # Find Step and PE column
    try:
        step_col = header_cols.index("Step")
    except ValueError:
        return None

    pe_col = None
    for cand in ("PotEng", "PE", "pe"):
        if cand in header_cols:
            pe_col = header_cols.index(cand)
            break
    if pe_col is None:
        # try case-insensitive match
        for j, c in enumerate(header_cols):
            if c.lower() in ("poteng", "pe"):
                pe_col = j
                break
    if pe_col is None:
        return None

    steps: List[int] = []
    pes: List[float] = []

    for line in lines[header_idx + 1 :]:
        cols = line.split()
        if len(cols) < max(step_col, pe_col) + 1:
            continue
        # Stop if we hit timing breakdown or other section
        if cols[0] in ("Loop", "Performance:", "MPI", "Section"):
            break
        try:
            s = int(float(cols[step_col]))
            pe = float(cols[pe_col])
        except ValueError:
            continue
        steps.append(s)
        pes.append(pe)

    if not steps:
        return None
    return ThermoSeries(steps=steps, pe=pes)


def parse_thermo_pe(path: Path) -> ThermoSeries:
    lines = path.read_text(errors="ignore").splitlines()
    series = _try_parse_thermo_table(lines)
    if series is None:
        raise RuntimeError(
            f"Could not find/parse thermo table with Step and PotEng/pe in {path}"
        )
    return series


def compute_event_work(series: ThermoSeries, event_steps: List[int]) -> List[Tuple[int, int, float, float, float]]:
    """Backward-compatible helper (kept), but prefer compute_step_work()."""
    step_list = series.steps
    pe_list = series.pe
    rows: List[Tuple[int, int, float, float, float]] = []

    step_to_pe: Dict[int, float] = {s: pe for s, pe in zip(step_list, pe_list)}
    for t in event_steps:
        if t not in step_to_pe:
            # event step missing in thermo; use nearest previous + nearest next if available
            idx = bisect.bisect_left(step_list, t)
            if idx == 0 or idx >= len(step_list):
                continue
            t_prev = step_list[idx - 1]
            t_next = step_list[idx]
            # Approximate pe_event as pe at next printed step (less ideal)
            pe_prev = step_to_pe[t_prev]
            pe_event = step_to_pe[t_next]
            work = pe_event - pe_prev
            rows.append((t, t_prev, pe_prev, pe_event, work))
            continue

        idx = bisect.bisect_left(step_list, t)
        if idx == 0:
            continue
        t_prev = step_list[idx - 1]
        pe_prev = step_to_pe[t_prev]
        pe_event = step_to_pe[t]
        work = pe_event - pe_prev
        rows.append((t, t_prev, pe_prev, pe_event, work))
    return rows


def compute_step_work(
    series: ThermoSeries, event_steps: List[int]
) -> List[Tuple[int, int, int, float, float, float]]:
    """
    Aggregate by timestep so multiple flips at the same timestep do NOT double-count energy.

    Returns list of rows:
      (t_event, n_events, t_prev, pe_prev, pe_event, work)
    """
    # count flips per step
    counts: Dict[int, int] = {}
    for t in event_steps:
        counts[t] = counts.get(t, 0) + 1
    unique_steps = sorted(counts.keys())

    step_list = series.steps
    step_to_pe: Dict[int, float] = {s: pe for s, pe in zip(series.steps, series.pe)}

    rows: List[Tuple[int, int, int, float, float, float]] = []
    for t in unique_steps:
        n_events = counts[t]

        idx = bisect.bisect_left(step_list, t)
        if idx == 0:
            continue

        t_prev = step_list[idx - 1]
        pe_prev = step_to_pe[t_prev]

        if t in step_to_pe:
            pe_event = step_to_pe[t]
        else:
            # If the flip step isn't printed in thermo, approximate with next printed step
            if idx >= len(step_list):
                continue
            pe_event = step_to_pe[step_list[idx]]

        work = pe_event - pe_prev
        rows.append((t, n_events, t_prev, pe_prev, pe_event, work))

    return rows


def maybe_plot(out_png: Path, rows: List[Tuple[int, int, float, float, float]]) -> None:
    try:
        import matplotlib.pyplot as plt  # type: ignore
    except Exception:
        return

    ts = [r[0] for r in rows]
    works = [r[4] for r in rows]
    cum = []
    c = 0.0
    for w in works:
        c += w
        cum.append(c)

    fig, ax = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
    ax[0].plot(ts, works, lw=1)
    ax[0].set_ylabel("Work per event (ΔPE)")
    ax[0].axhline(0.0, color="k", lw=0.5)
    ax[1].plot(ts, cum, lw=1)
    ax[1].set_ylabel("Cumulative work")
    ax[1].set_xlabel("Timestep")
    fig.tight_layout()
    fig.savefig(out_png, dpi=200)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--thermo", required=True, help="Path to LAMMPS thermo output (lammps_stdout.log or slurm .out)")
    ap.add_argument("--events", required=True, help="Path to statechange log (slurm .err with STATECHANGE lines)")
    ap.add_argument("--out", required=True, help="Output prefix (directory/filename prefix)")
    args = ap.parse_args()

    thermo_path = Path(args.thermo)
    events_path = Path(args.events)
    out_prefix = Path(args.out)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    # Prefer instantaneous dE if present in the events file; otherwise fall back to ΔPE from thermo.
    step_dEs = parse_statechange_step_dE(events_path)
    if step_dEs:
        # Aggregate dE by timestep
        counts: Dict[int, int] = {}
        sums: Dict[int, float] = {}
        for t, de in step_dEs:
            counts[t] = counts.get(t, 0) + 1
            sums[t] = sums.get(t, 0.0) + de
        steps = sorted(counts.keys())

        out_csv = out_prefix.with_suffix(".csv")
        with out_csv.open("w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["t_event", "n_flips", "work_dE_sum", "cum_work"])
            cum = 0.0
            for t in steps:
                cum += sums[t]
                w.writerow([t, counts[t], sums[t], cum])

        out_png = out_prefix.with_suffix(".png")
        try:
            import matplotlib.pyplot as plt  # type: ignore
        except Exception:
            print(f"Wrote {out_csv}")
            print("matplotlib not available; skipped plot (CSV still written)")
            return

        works = [sums[t] for t in steps]
        cumw = []
        c = 0.0
        for wv in works:
            c += wv
            cumw.append(c)
        fig, ax = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
        ax[0].plot(steps, works, lw=1)
        ax[0].set_ylabel("Instant work per timestep (Σ dE)")
        ax[0].axhline(0.0, color="k", lw=0.5)
        ax[1].plot(steps, cumw, lw=1)
        ax[1].set_ylabel("Cumulative work")
        ax[1].set_xlabel("Timestep")
        fig.tight_layout()
        fig.savefig(out_png, dpi=200)

        print(f"Wrote {out_csv}")
        print(f"Wrote {out_png}")
        return

    # Fallback: ΔPE from thermo output (coarser)
    series = parse_thermo_pe(thermo_path)
    event_steps = parse_statechange_steps(events_path)
    rows = compute_step_work(series, event_steps)

    out_csv = out_prefix.with_suffix(".csv")
    with out_csv.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["t_event", "n_flips", "t_prev", "pe_prev", "pe_event", "work_dpe", "cum_work"])
        cum = 0.0
        for (t_event, n_flips, t_prev, pe_prev, pe_event, work) in rows:
            cum += work
            w.writerow([t_event, n_flips, t_prev, pe_prev, pe_event, work, cum])

    out_png = out_prefix.with_suffix(".png")
    # Plot expects (t_event, ..., work) at index 5 with this row shape, so adapt:
    plot_rows = [(r[0], r[2], r[3], r[4], r[5]) for r in rows]
    maybe_plot(out_png, plot_rows)

    print(f"Wrote {out_csv}")
    if out_png.exists():
        print(f"Wrote {out_png}")
    else:
        print("matplotlib not available; skipped plot (CSV still written)")


if __name__ == "__main__":
    main()


