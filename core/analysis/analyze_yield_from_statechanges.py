#!/usr/bin/env python3
"""
Plot "yield" (normalized concentrations) over time for A, C, E, F in the twins simulations.

Key idea:
- We do NOT parse the big trajectory dump.
- We parse initial molecule identities from the LAMMPS data file (by molecule id / atom types).
- We parse state-change events from SLURM stderr (STATECHANGE lines).
- Since the dynamics are irreversible (A->C and E->F), we can reconstruct counts over time.

Yield definition:
- For each timestep t, define counts A(t), C(t), E(t), F(t).
- Yield is normalized over ONLY these 4 species:
    yA = A / (A+C+E+F), etc
  so yA+yC+yE+yF = 1.0 always (up to float error).

Inputs:
  --data   path to data file (e.g. data.dimer_ksat_1core_twosideB_twins)
  --events path to slurm .err (or any file containing STATECHANGE lines)
  --out    output prefix (writes .csv and .png)
"""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple


# Matches both:
#  STATECHANGE ...: step 1234 molA 5 flipped 2->4 ...
#  STATECHANGE ...: step 1234 molE 7 flipped 8->10 ...
EVENT_RE = re.compile(
    r"STATECHANGE\s+\S+:\s+step\s+(?P<step>\d+)\s+(?P<tag>molA|molE)\s+(?P<mol>\d+)\s+flipped\s+(?P<from>\d+)->(?P<to>\d+)"
)


def parse_data_molecule_species(data_path: Path) -> Dict[int, str]:
    """
    Return mapping mol_id -> one of {"A","C","E","F","OTHER"} inferred from atom types present.

    We classify based on presence of patch atom types:
    - A: has any type 2
    - C: has any type 4
    - E: has any type 8
    - F: has any type 10
    If none of these patch types appear, classify as OTHER.
    """
    text = data_path.read_text(errors="ignore").splitlines()
    # Find "Atoms" section (we assume "Atoms # full")
    start = None
    for i, line in enumerate(text):
        if line.strip().startswith("Atoms"):
            start = i
            break
    if start is None:
        raise RuntimeError(f"Could not find Atoms section in {data_path}")

    # After "Atoms" header, there may be blank line(s); then atom lines until next section
    i = start + 1
    while i < len(text) and text[i].strip() == "":
        i += 1

    mol_patch_types: Dict[int, Set[int]] = defaultdict(set)
    while i < len(text):
        line = text[i].strip()
        i += 1
        if not line:
            continue
        # next section header begins with a word like Bonds/Angles/Velocities/etc
        if line[0].isalpha():
            break
        # atom line: id mol type q x y z
        parts = line.split()
        if len(parts) < 4:
            continue
        try:
            mol = int(parts[1])
            typ = int(parts[2])
        except ValueError:
            continue
        mol_patch_types[mol].add(typ)

    mol_species: Dict[int, str] = {}
    for mol, types in mol_patch_types.items():
        # Prefer exact mapping; if multiple markers appear (shouldn't), pick in this order:
        if 2 in types:
            mol_species[mol] = "A"
        elif 4 in types:
            mol_species[mol] = "C"
        elif 8 in types:
            mol_species[mol] = "E"
        elif 10 in types:
            mol_species[mol] = "F"
        else:
            mol_species[mol] = "OTHER"
    return mol_species


@dataclass
class Event:
    step: int
    mol: int
    tag: str  # molA or molE
    from_type: int
    to_type: int


def parse_events(events_path: Path) -> List[Event]:
    events: List[Event] = []
    for line in events_path.read_text(errors="ignore").splitlines():
        m = EVENT_RE.search(line)
        if not m:
            continue
        events.append(
            Event(
                step=int(m.group("step")),
                mol=int(m.group("mol")),
                tag=m.group("tag"),
                from_type=int(m.group("from")),
                to_type=int(m.group("to")),
            )
        )
    events.sort(key=lambda e: e.step)
    return events


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True, help="LAMMPS data file path")
    ap.add_argument("--events", required=True, help="Statechange log path (slurm .err)")
    ap.add_argument("--out", required=True, help="Output prefix (writes .csv and .png)")
    args = ap.parse_args()

    data_path = Path(args.data)
    events_path = Path(args.events)
    out_prefix = Path(args.out)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    mol_species = parse_data_molecule_species(data_path)

    # Initial counts
    counts: Dict[str, int] = {"A": 0, "C": 0, "E": 0, "F": 0}
    for sp in mol_species.values():
        if sp in counts:
            counts[sp] += 1

    total = sum(counts.values())
    if total == 0:
        raise RuntimeError("No A/C/E/F molecules found in data file; cannot compute yields.")

    # Track which molecules have already flipped (prevent double count)
    flipped_A: Set[int] = set()
    flipped_E: Set[int] = set()

    events = parse_events(events_path)

    # Aggregate events by step
    by_step: Dict[int, List[Event]] = defaultdict(list)
    for e in events:
        by_step[e.step].append(e)

    steps_sorted = sorted(by_step.keys())

    # Build time series: include t=0 row
    rows: List[Tuple[int, int, int, int, int]] = []
    rows.append((0, counts["A"], counts["C"], counts["E"], counts["F"]))

    for step in steps_sorted:
        for e in by_step[step]:
            if e.tag == "molA":
                if e.mol in flipped_A:
                    continue
                # Expect A->C
                flipped_A.add(e.mol)
                if counts["A"] > 0:
                    counts["A"] -= 1
                counts["C"] += 1
            elif e.tag == "molE":
                if e.mol in flipped_E:
                    continue
                flipped_E.add(e.mol)
                if counts["E"] > 0:
                    counts["E"] -= 1
                counts["F"] += 1

        rows.append((step, counts["A"], counts["C"], counts["E"], counts["F"]))

    out_csv = out_prefix.with_suffix(".csv")
    with out_csv.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["timestep", "A", "C", "E", "F", "yA", "yC", "yE", "yF"])
        for (t, a, c, e, ff) in rows:
            denom = a + c + e + ff
            if denom <= 0:
                yA = yC = yE = yF = 0.0
            else:
                yA = a / denom
                yC = c / denom
                yE = e / denom
                yF = ff / denom
            w.writerow([t, a, c, e, ff, yA, yC, yE, yF])

    # Plot (optional)
    out_png = out_prefix.with_suffix(".png")
    try:
        import matplotlib.pyplot as plt  # type: ignore
    except Exception:
        print(f"Wrote {out_csv}")
        print("matplotlib not available; skipped plot (CSV still written)")
        return

    ts = [r[0] for r in rows]
    ys = []
    for (t, a, c, e, ff) in rows:
        denom = a + c + e + ff
        if denom <= 0:
            ys.append((0.0, 0.0, 0.0, 0.0))
        else:
            ys.append((a / denom, c / denom, e / denom, ff / denom))

    yA = [u[0] for u in ys]
    yC = [u[1] for u in ys]
    yE = [u[2] for u in ys]
    yF = [u[3] for u in ys]

    fig, ax = plt.subplots(1, 1, figsize=(10, 4))
    ax.step(ts, yA, where="post", label="A")
    ax.step(ts, yC, where="post", label="C")
    ax.step(ts, yE, where="post", label="E")
    ax.step(ts, yF, where="post", label="F")
    ax.set_xlabel("Timestep")
    ax.set_ylabel("Yield (A+C+E+F normalized)")
    ax.set_ylim(-0.02, 1.02)
    ax.legend(ncol=4, fontsize=9)
    fig.tight_layout()
    fig.savefig(out_png, dpi=200)

    print(f"Wrote {out_csv}")
    print(f"Wrote {out_png}")


if __name__ == "__main__":
    main()


