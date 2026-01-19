#!/usr/bin/env python3
"""
Calculate work from state changes using consecutive dump/thermo frames.

This script addresses thermal fluctuation noise by:
1. Only calculating work AT state change events (not continuously)
2. Using the smallest possible time window between measurements
3. Preferring instantaneous dE if printed by the fix
4. Otherwise using PE difference between closest consecutive frames

Scientific rationale:
- State changes are configuration changes (type flips)
- Work should reflect the potential energy change from the new interaction landscape
- Thermal KE fluctuations are maintained by thermostat (not relevant to work)
- We want ΔPE at the moment of state change, not integrated over many timesteps

Usage:
    python analyze_work_statechange_frames.py \\
        --thermo lammps_stdout.log \\
        --events slurm-12345.err \\
        --out work_analysis

Output:
    work_analysis.csv - Per-event work values
    work_analysis.png - Plots of work per event and cumulative work
"""

from __future__ import annotations

import argparse
import bisect
import csv
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# Regex patterns
STATECHANGE_RE = re.compile(r"STATECHANGE\s+\S+:\s+step\s+(\d+)\b")
DENERGY_RE = re.compile(r"\bdE\s+([+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)\b")


@dataclass(frozen=True)
class ThermoSeries:
    """Time series of potential energy from LAMMPS thermo output."""
    steps: List[int]
    pe: List[float]


def parse_statechange_events(path: Path) -> List[Tuple[int, Optional[float]]]:
    """
    Parse state change events from stderr/events file.

    Returns list of (timestep, dE) tuples.
    If dE is not printed in the line, returns (timestep, None).
    """
    events: List[Tuple[int, Optional[float]]] = []

    for line in path.read_text(errors="ignore").splitlines():
        m = STATECHANGE_RE.search(line)
        if not m:
            continue

        timestep = int(m.group(1))

        # Try to extract instantaneous dE if present
        dE_match = DENERGY_RE.search(line)
        dE = float(dE_match.group(1)) if dE_match else None

        events.append((timestep, dE))

    events.sort(key=lambda x: x[0])
    return events


def parse_thermo_pe(path: Path) -> ThermoSeries:
    """
    Parse LAMMPS thermo output to extract Step and PotEng columns.
    """
    lines = path.read_text(errors="ignore").splitlines()

    # Find thermo table header
    header_idx = None
    header_cols: List[str] = []

    for i, line in enumerate(lines):
        if "Step" in line and ("PotEng" in line or "pe" in line.lower()):
            header_idx = i
            header_cols = line.split()
            break

    if header_idx is None:
        raise RuntimeError(f"Could not find thermo table header with Step and PotEng in {path}")

    # Find column indices
    try:
        step_col = header_cols.index("Step")
    except ValueError:
        raise RuntimeError(f"No 'Step' column in header: {header_cols}")

    # Try multiple variations of PE column name
    pe_col = None
    for cand in ("PotEng", "PE", "pe"):
        if cand in header_cols:
            pe_col = header_cols.index(cand)
            break

    if pe_col is None:
        # Case-insensitive fallback
        for j, c in enumerate(header_cols):
            if c.lower() in ("poteng", "pe"):
                pe_col = j
                break

    if pe_col is None:
        raise RuntimeError(f"No PotEng/PE column in header: {header_cols}")

    # Parse data rows
    steps: List[int] = []
    pes: List[float] = []

    for line in lines[header_idx + 1:]:
        cols = line.split()

        # Stop at performance summary or other sections
        if not cols or cols[0] in ("Loop", "Performance:", "MPI", "Section"):
            break

        if len(cols) <= max(step_col, pe_col):
            continue

        try:
            step = int(float(cols[step_col]))
            pe = float(cols[pe_col])
            steps.append(step)
            pes.append(pe)
        except (ValueError, IndexError):
            continue

    if not steps:
        raise RuntimeError(f"Could not parse any thermo data from {path}")

    return ThermoSeries(steps=steps, pe=pes)


def compute_work_from_consecutive_frames(
    events: List[Tuple[int, Optional[float]]],
    thermo: ThermoSeries
) -> List[Dict[str, float]]:
    """
    Calculate work for each state change event.

    Strategy:
    1. If event has instantaneous dE, use that (most accurate)
    2. Otherwise, find PE at event timestep and previous thermo step
    3. Calculate work = PE(event) - PE(previous)

    Returns list of dictionaries with:
        - timestep: when state change occurred
        - pe_before: PE before state change
        - pe_after: PE at/after state change
        - work: PE(after) - PE(before)
        - source: 'instantaneous_dE' or 'frame_diff'
        - n_steps_between: timesteps between measurements (for quality assessment)
    """
    step_to_pe: Dict[int, float] = {s: pe for s, pe in zip(thermo.steps, thermo.pe)}
    results: List[Dict[str, float]] = []

    # Group events by timestep to avoid double-counting multiple flips at same time
    timestep_events: Dict[int, List[Optional[float]]] = {}
    for t, dE in events:
        if t not in timestep_events:
            timestep_events[t] = []
        timestep_events[t].append(dE)

    for event_step in sorted(timestep_events.keys()):
        dE_values = timestep_events[event_step]

        # If ALL events at this timestep have instantaneous dE, use sum of dE
        if all(dE is not None for dE in dE_values):
            total_dE = sum(dE for dE in dE_values if dE is not None)

            # Get PE at event step for reference (if available)
            if event_step in step_to_pe:
                pe_after = step_to_pe[event_step]
                pe_before = pe_after - total_dE
            else:
                # Approximate: find nearest PE
                idx = bisect.bisect_left(thermo.steps, event_step)
                if idx >= len(thermo.steps):
                    idx = len(thermo.steps) - 1
                pe_after = step_to_pe[thermo.steps[idx]]
                pe_before = pe_after - total_dE

            results.append({
                'timestep': event_step,
                'n_flips': len(dE_values),
                'pe_before': pe_before,
                'pe_after': pe_after,
                'work': total_dE,
                'source': 'instantaneous_dE',
                'n_steps_between': 0  # Instantaneous
            })
            continue

        # Otherwise, calculate from consecutive thermo frames
        idx = bisect.bisect_left(thermo.steps, event_step)

        if idx == 0:
            # Event before first thermo output - skip
            continue

        # Get previous thermo step
        step_before = thermo.steps[idx - 1]
        pe_before = step_to_pe[step_before]

        # Get PE at or after event
        if event_step in step_to_pe:
            step_after = event_step
            pe_after = step_to_pe[event_step]
        elif idx < len(thermo.steps):
            step_after = thermo.steps[idx]
            pe_after = step_to_pe[step_after]
        else:
            # Event after last thermo output - skip
            continue

        work = pe_after - pe_before
        n_steps = step_after - step_before

        results.append({
            'timestep': event_step,
            'n_flips': len(dE_values),
            'pe_before': pe_before,
            'pe_after': pe_after,
            'work': work,
            'source': 'frame_diff',
            'n_steps_between': n_steps
        })

    return results


def plot_work(results: List[Dict[str, float]], out_path: Path) -> None:
    """Create plots of work per event and cumulative work."""
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available; skipping plot")
        return

    if not results:
        print("No results to plot")
        return

    timesteps = [r['timestep'] for r in results]
    work_values = [r['work'] for r in results]

    # Calculate cumulative work
    cumulative = []
    total = 0.0
    for w in work_values:
        total += w
        cumulative.append(total)

    # Create figure with two subplots
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

    # Color by source
    colors = ['blue' if r['source'] == 'instantaneous_dE' else 'orange' for r in results]

    # Top: work per event
    ax1.scatter(timesteps, work_values, c=colors, alpha=0.6, s=20)
    ax1.axhline(0, color='black', linestyle='--', linewidth=0.5, alpha=0.5)
    ax1.set_ylabel('Work per state change event (ΔPE)', fontsize=11)
    ax1.set_title('Work from State Changes (avoiding thermal fluctuations)', fontsize=12)
    ax1.grid(True, alpha=0.3)

    # Add legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='blue', label='Instantaneous dE'),
        Patch(facecolor='orange', label='Frame difference')
    ]
    ax1.legend(handles=legend_elements, loc='best', fontsize=9)

    # Bottom: cumulative work
    ax2.plot(timesteps, cumulative, color='darkgreen', linewidth=1.5)
    ax2.set_xlabel('Timestep', fontsize=11)
    ax2.set_ylabel('Cumulative work', fontsize=11)
    ax2.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches='tight')
    print(f"Wrote plot to {out_path}")


def write_csv(results: List[Dict[str, float]], out_path: Path) -> None:
    """Write results to CSV file."""
    with out_path.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=[
            'timestep', 'n_flips', 'pe_before', 'pe_after', 'work',
            'source', 'n_steps_between', 'cumulative_work'
        ])
        writer.writeheader()

        cumulative = 0.0
        for r in results:
            cumulative += r['work']
            row = r.copy()
            row['cumulative_work'] = cumulative
            writer.writerow(row)

    print(f"Wrote {len(results)} events to {out_path}")


def print_summary(results: List[Dict[str, float]]) -> None:
    """Print summary statistics."""
    if not results:
        print("No results to summarize")
        return

    total_events = sum(r['n_flips'] for r in results)
    n_instant = sum(1 for r in results if r['source'] == 'instantaneous_dE')
    n_frame = len(results) - n_instant

    total_work = sum(r['work'] for r in results)
    mean_work = total_work / len(results) if results else 0.0

    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    print(f"Total state change events: {total_events}")
    print(f"Unique timesteps with events: {len(results)}")
    print(f"  - Using instantaneous dE: {n_instant}")
    print(f"  - Using frame difference: {n_frame}")
    print(f"\nTotal work: {total_work:.6f}")
    print(f"Mean work per event: {mean_work:.6f}")

    if n_frame > 0:
        frame_diff_steps = [r['n_steps_between'] for r in results if r['source'] == 'frame_diff']
        mean_steps = sum(frame_diff_steps) / len(frame_diff_steps)
        print(f"\nFor frame differences:")
        print(f"  Mean timesteps between measurements: {mean_steps:.1f}")
        print(f"  (smaller is better - less thermal fluctuation)")

    print("="*60 + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Calculate work from state changes using consecutive frames"
    )
    parser.add_argument(
        "--thermo",
        type=Path,
        required=True,
        help="LAMMPS thermo output file (log.lammps or stdout)"
    )
    parser.add_argument(
        "--events",
        type=Path,
        required=True,
        help="State change events file (stderr with STATECHANGE lines)"
    )
    parser.add_argument(
        "--out",
        type=Path,
        required=True,
        help="Output prefix (will create .csv and .png)"
    )

    args = parser.parse_args()

    # Create output directory if needed
    args.out.parent.mkdir(parents=True, exist_ok=True)

    # Parse inputs
    print(f"Parsing state change events from {args.events}...")
    events = parse_statechange_events(args.events)
    print(f"Found {len(events)} state change events")

    print(f"\nParsing thermo data from {args.thermo}...")
    thermo = parse_thermo_pe(args.thermo)
    print(f"Found {len(thermo.steps)} thermo outputs from step {thermo.steps[0]} to {thermo.steps[-1]}")

    # Compute work
    print("\nComputing work from consecutive frames...")
    results = compute_work_from_consecutive_frames(events, thermo)

    # Print summary
    print_summary(results)

    # Write outputs
    csv_path = args.out.with_suffix('.csv')
    write_csv(results, csv_path)

    png_path = args.out.with_suffix('.png')
    plot_work(results, png_path)


if __name__ == "__main__":
    main()
