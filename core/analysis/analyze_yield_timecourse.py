#!/usr/bin/env python3
"""
Compute and plot time-dependent "yields" (normalized concentrations) of monomer types
A, C, E, F from a LAMMPS dump (.lammpstrj).

Definition used:
- We classify each MOLECULE (by mol id) as A/C/E/F if it contains patch atoms of:
    A: patch type 2
    C: patch type 4
    E: patch type 8
    F: patch type 10
  (These correspond to the "twins" simulation conventions.)
- At each timestep, we count N_A, N_C, N_E, N_F as number of molecules of each type.
- Yield is normalized so that:
    y_A + y_C + y_E + y_F = 1
  i.e., divide by N_total = N_A + N_C + N_E + N_F.

Outputs:
- <out>.csv with counts and yields vs timestep
- <out>.png if matplotlib is available

This script streams the dump so it can handle large trajectories.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple


DEFAULT_MAP = {
    "A": 2,
    "C": 4,
    "E": 8,
    "F": 10,
}


@dataclass
class FrameYields:
    timestep: int
    n_A: int
    n_C: int
    n_E: int
    n_F: int

    @property
    def n_total(self) -> int:
        return self.n_A + self.n_C + self.n_E + self.n_F

    def yields(self) -> Tuple[float, float, float, float]:
        tot = self.n_total
        if tot <= 0:
            return (float("nan"),) * 4
        return (self.n_A / tot, self.n_C / tot, self.n_E / tot, self.n_F / tot)


def _parse_atoms_header(line: str) -> List[str]:
    # "ITEM: ATOMS id mol type x y z ..."
    parts = line.strip().split()
    if len(parts) < 3 or parts[0] != "ITEM:" or parts[1] != "ATOMS":
        raise ValueError(f"Unexpected ATOMS header: {line!r}")
    return parts[2:]


def iter_lammpstrj_frames(
    path: Path,
    stride: int = 1,
    max_frames: Optional[int] = None,
) -> Iterable[Tuple[int, List[str], List[List[str]]]]:
    """
    Yields tuples: (timestep, atom_columns, atom_rows)
    where atom_rows are tokenized strings (no numeric conversion here).
    """
    with path.open("r", errors="ignore") as f:
        frame_idx = 0
        yielded = 0
        while True:
            line = f.readline()
            if not line:
                return
            if not line.startswith("ITEM: TIMESTEP"):
                continue
            ts_line = f.readline()
            if not ts_line:
                return
            timestep = int(ts_line.strip())

            # ITEM: NUMBER OF ATOMS
            line = f.readline()
            if not line:
                return
            if not line.startswith("ITEM: NUMBER OF ATOMS"):
                raise ValueError(f"Expected NUMBER OF ATOMS, got: {line!r}")
            n_atoms = int(f.readline().strip())

            # ITEM: BOX BOUNDS
            line = f.readline()
            if not line:
                return
            if not line.startswith("ITEM: BOX BOUNDS"):
                raise ValueError(f"Expected BOX BOUNDS, got: {line!r}")
            # skip 3 box lines
            for _ in range(3):
                _ = f.readline()

            # ITEM: ATOMS ...
            line = f.readline()
            if not line:
                return
            cols = _parse_atoms_header(line)

            # decide whether to parse this frame
            take = (frame_idx % stride == 0)
            rows: List[List[str]] = []
            for _ in range(n_atoms):
                atom_line = f.readline()
                if not atom_line:
                    return
                if take:
                    rows.append(atom_line.split())

            if take:
                yield (timestep, cols, rows)
                yielded += 1
                if max_frames is not None and yielded >= max_frames:
                    return

            frame_idx += 1


def classify_molecules(
    cols: List[str],
    rows: List[List[str]],
    patch_types: Dict[str, int],
) -> FrameYields:
    """
    Determine molecule counts for A/C/E/F by scanning atom types per molecule.
    """
    try:
        mol_idx = cols.index("mol")
    except ValueError:
        mol_idx = cols.index("molecule") if "molecule" in cols else None  # type: ignore
    if mol_idx is None:
        raise ValueError("Could not find 'mol' column in ATOMS header")
    type_idx = cols.index("type")

    # Collect which patch types appear in each molecule.
    mol_seen: Dict[int, set] = {}
    target_types = set(patch_types.values())
    for r in rows:
        try:
            mol = int(float(r[mol_idx]))
            t = int(float(r[type_idx]))
        except (ValueError, IndexError):
            continue
        if mol <= 0:
            continue
        if t not in target_types:
            continue
        s = mol_seen.get(mol)
        if s is None:
            s = set()
            mol_seen[mol] = s
        s.add(t)

    # Classify each molecule. Prefer "A over C" and "E over F" if both appear (shouldn't happen
    # if flips update all patches, but dump timing could show mixed states).
    n_A = n_C = n_E = n_F = 0
    tA = patch_types["A"]
    tC = patch_types["C"]
    tE = patch_types["E"]
    tF = patch_types["F"]

    for types in mol_seen.values():
        if tA in types:
            n_A += 1
        elif tC in types:
            n_C += 1
        elif tE in types:
            n_E += 1
        elif tF in types:
            n_F += 1

    # timestep is assigned by caller
    return FrameYields(timestep=0, n_A=n_A, n_C=n_C, n_E=n_E, n_F=n_F)


def maybe_plot(out_png: Path, frames: List[FrameYields]) -> None:
    try:
        import matplotlib.pyplot as plt  # type: ignore
    except Exception:
        return

    ts = [fr.timestep for fr in frames]
    yA, yC, yE, yF = zip(*(fr.yields() for fr in frames))

    plt.figure(figsize=(10, 5))
    plt.plot(ts, yA, label="A")
    plt.plot(ts, yC, label="C")
    plt.plot(ts, yE, label="E")
    plt.plot(ts, yF, label="F")
    plt.xlabel("Timestep")
    plt.ylabel("Yield (normalized within A+C+E+F)")
    plt.ylim(0, 1)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump", required=True, help="Path to .lammpstrj dump file")
    ap.add_argument("--out", required=True, help="Output prefix (e.g. analysis/yield_twins)")
    ap.add_argument("--stride", type=int, default=1, help="Only analyze every Nth frame (default: 1)")
    ap.add_argument("--max_frames", type=int, default=None, help="Stop after this many analyzed frames")
    args = ap.parse_args()

    dump_path = Path(args.dump)
    out_prefix = Path(args.out)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    frames: List[FrameYields] = []
    for (timestep, cols, rows) in iter_lammpstrj_frames(
        dump_path, stride=args.stride, max_frames=args.max_frames
    ):
        fr = classify_molecules(cols, rows, DEFAULT_MAP)
        fr.timestep = timestep
        frames.append(fr)

    out_csv = out_prefix.with_suffix(".csv")
    with out_csv.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["timestep", "n_A", "n_C", "n_E", "n_F", "n_total", "y_A", "y_C", "y_E", "y_F"])
        for fr in frames:
            yA, yC, yE, yF = fr.yields()
            w.writerow([fr.timestep, fr.n_A, fr.n_C, fr.n_E, fr.n_F, fr.n_total, yA, yC, yE, yF])

    out_png = out_prefix.with_suffix(".png")
    maybe_plot(out_png, frames)

    print(f"Wrote {out_csv}")
    if out_png.exists():
        print(f"Wrote {out_png}")
    else:
        print("matplotlib not available; skipped plot (CSV still written)")


if __name__ == "__main__":
    main()


