#!/usr/bin/env python3
"""
Plot/compute yields (normalized concentrations) over time for monomers A, C, E, F.

We define yield as fraction among ONLY these four species:
  y_X(t) = N_X(t) / (N_A(t)+N_C(t)+N_E(t)+N_F(t))

This script avoids parsing huge dump trajectories by using:
1) Initial composition from a LAMMPS data file (counts by molecule ID and patch types)
2) State-change events from the SLURM .err (STATECHANGE lines)

Works for the "twins" setup where:
- A->C flips are logged as: "... flipped 2->4 ..."
- E->F flips are logged as: "... flipped 8->10 ..."

Outputs:
- CSV with counts and yields vs timestep
- PNG plot if matplotlib is available
"""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Tuple


STATECHANGE_RE = re.compile(r"STATECHANGE\s+\S+:\s+step\s+(\d+)\s+(.+)$")
FLIP_AC_RE = re.compile(r"\bflipped\s+2->4\b")
FLIP_EF_RE = re.compile(r"\bflipped\s+8->10\b")


@dataclass(frozen=True)
class InitialCounts:
    nA: int
    nC: int
    nE: int
    nF: int


def _iter_data_atoms_section(lines: List[str]) -> Iterable[Tuple[int, int, int]]:
    """
    Yield tuples: (mol_id, atom_type, atom_id) from the 'Atoms' section.
    Supports 'Atoms # full' format used by our generators.
    """
    # Find "Atoms" section
    start = None
    for i, line in enumerate(lines):
        if line.strip().startswith("Atoms"):
            start = i
            break
    if start is None:
        raise RuntimeError("Could not find 'Atoms' section in data file")

    # data files have a blank line after section header
    i = start + 1
    while i < len(lines) and lines[i].strip() == "":
        i += 1

    for line in lines[i:]:
        s = line.strip()
        if not s:
            continue
        # Stop when next section starts
        if s[0].isalpha():
            break
        # Expected: atom_id mol_id type q x y z
        cols = s.split()
        if len(cols) < 3:
            continue
        atom_id = int(cols[0])
        mol_id = int(cols[1])
        atom_type = int(cols[2])
        yield mol_id, atom_type, atom_id


def infer_initial_counts_from_data(data_path: Path) -> InitialCounts:
    """
    Infer monomer counts by molecule ID using patch types.
    Classification (by presence of patch types in the molecule):
      - A: has type 2 patches
      - C: has type 4 patches
      - E: has type 8 patches
      - F: has type 10 patches
    """
    lines = data_path.read_text(errors="ignore").splitlines()
    mol_patch_types: Dict[int, set] = defaultdict(set)

    for mol_id, atom_type, _atom_id in _iter_data_atoms_section(lines):
        # only track patch types of interest
        if atom_type in (2, 4, 8, 10):
            mol_patch_types[mol_id].add(atom_type)

    nA = nC = nE = nF = 0
    for mol_id, pts in mol_patch_types.items():
        # A/C/E/F should be mutually exclusive by construction
        if 2 in pts:
            nA += 1
        elif 4 in pts:
            nC += 1
        elif 8 in pts:
            nE += 1
        elif 10 in pts:
            nF += 1

    return InitialCounts(nA=nA, nC=nC, nE=nE, nF=nF)


def parse_flip_events(events_path: Path) -> Dict[int, Dict[str, int]]:
    """
    Returns dict: step -> { "AC": count, "EF": count }
    """
    per_step = defaultdict(lambda: {"AC": 0, "EF": 0})
    for line in events_path.read_text(errors="ignore").splitlines():
        m = STATECHANGE_RE.search(line)
        if not m:
            continue
        step = int(m.group(1))
        if FLIP_AC_RE.search(line):
            per_step[step]["AC"] += 1
        if FLIP_EF_RE.search(line):
            per_step[step]["EF"] += 1
    return dict(per_step)


def compute_yields_over_time(
    init: InitialCounts, flips_by_step: Dict[int, Dict[str, int]]
) -> List[Dict[str, float]]:
    """
    Produces rows with:
      step, N_A,N_C,N_E,N_F, y_A,y_C,y_E,y_F
    """
    nA, nC, nE, nF = init.nA, init.nC, init.nE, init.nF
    denom0 = nA + nC + nE + nF
    if denom0 <= 0:
        raise RuntimeError("Initial denominator N_A+N_C+N_E+N_F is 0; cannot define yields.")

    steps = sorted(flips_by_step.keys())

    rows: List[Dict[str, float]] = []
    # include step 0 snapshot
    rows.append(
        {
            "step": 0,
            "N_A": nA,
            "N_C": nC,
            "N_E": nE,
            "N_F": nF,
            "y_A": nA / denom0,
            "y_C": nC / denom0,
            "y_E": nE / denom0,
            "y_F": nF / denom0,
        }
    )

    for step in steps:
        dAC = flips_by_step[step]["AC"]
        dEF = flips_by_step[step]["EF"]

        nA = max(0, nA - dAC)
        nC = nC + dAC
        nE = max(0, nE - dEF)
        nF = nF + dEF

        denom = nA + nC + nE + nF
        if denom <= 0:
            denom = 1.0

        rows.append(
            {
                "step": step,
                "N_A": nA,
                "N_C": nC,
                "N_E": nE,
                "N_F": nF,
                "y_A": nA / denom,
                "y_C": nC / denom,
                "y_E": nE / denom,
                "y_F": nF / denom,
            }
        )

    return rows


def maybe_plot(out_png: Path, rows: List[Dict[str, float]]) -> None:
    try:
        import matplotlib.pyplot as plt  # type: ignore
    except Exception:
        return

    steps = [int(r["step"]) for r in rows]
    yA = [r["y_A"] for r in rows]
    yC = [r["y_C"] for r in rows]
    yE = [r["y_E"] for r in rows]
    yF = [r["y_F"] for r in rows]

    plt.figure(figsize=(10, 5))
    plt.plot(steps, yA, label="A")
    plt.plot(steps, yC, label="C")
    plt.plot(steps, yE, label="E")
    plt.plot(steps, yF, label="F")
    plt.xlabel("Timestep")
    plt.ylabel("Yield (normalized among A,C,E,F)")
    plt.ylim(-0.02, 1.02)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True, help="LAMMPS data file path (for initial counts)")
    ap.add_argument("--events", required=True, help="SLURM .err path containing STATECHANGE lines")
    ap.add_argument("--out", required=True, help="Output prefix (directory/filename prefix)")
    args = ap.parse_args()

    data_path = Path(args.data)
    events_path = Path(args.events)
    out_prefix = Path(args.out)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    init = infer_initial_counts_from_data(data_path)
    flips = parse_flip_events(events_path)
    rows = compute_yields_over_time(init, flips)

    out_csv = out_prefix.with_suffix(".csv")
    with out_csv.open("w", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=["step", "N_A", "N_C", "N_E", "N_F", "y_A", "y_C", "y_E", "y_F"],
        )
        w.writeheader()
        for r in rows:
            w.writerow(r)

    out_png = out_prefix.with_suffix(".png")
    maybe_plot(out_png, rows)

    print(f"Wrote {out_csv}")
    if out_png.exists():
        print(f"Wrote {out_png}")
    else:
        print("matplotlib not available; skipped plot (CSV still written)")


if __name__ == "__main__":
    main()


