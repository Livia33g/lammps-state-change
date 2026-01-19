#!/usr/bin/env python3
"""
Compute per-timestep "yields" (normalized concentrations) of A, C, E, F from a LAMMPS dump.

Definition:
- For each dump frame at timestep t, count molecules in each state:
    A: molecule has any atom of type 2  (A patches)
    C: molecule has any atom of type 4  (C patches)
    E: molecule has any atom of type 8  (E patches)
    F: molecule has any atom of type 10 (F patches)
  (These correspond to the twins simulation conventions.)
- Yield is normalized over ONLY these four:
    denom = N_A + N_C + N_E + N_F
    yA = N_A / denom, etc.
  So yA + yC + yE + yF == 1 (if denom > 0).

Outputs:
- CSV: timestep, N_A,N_C,N_E,N_F, yA,yC,yE,yF
- Optional PNG plot if matplotlib is available

Notes:
- This parser streams the dump file; it does not load all frames into memory.
- You can subsample frames with --stride_frames or stop early with --max_frames.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple


TYPE_A = 2
TYPE_C = 4
TYPE_E = 8
TYPE_F = 10


@dataclass
class FrameCounts:
    timestep: int
    nA: int
    nC: int
    nE: int
    nF: int

    @property
    def denom(self) -> int:
        return self.nA + self.nC + self.nE + self.nF

    def yields(self) -> Tuple[float, float, float, float]:
        d = self.denom
        if d <= 0:
            return (0.0, 0.0, 0.0, 0.0)
        return (self.nA / d, self.nC / d, self.nE / d, self.nF / d)


def _expect(line: str, prefix: str) -> None:
    if not line.startswith(prefix):
        raise RuntimeError(f"Expected line starting with {prefix!r}, got: {line[:80]!r}")


def iter_frames(
    dump_path: Path, stride_frames: int = 1, max_frames: Optional[int] = None
) -> Iterable[FrameCounts]:
    """
    Stream a LAMMPS custom dump in lammpstrj format.
    Assumes standard sections:
      ITEM: TIMESTEP
      <t>
      ITEM: NUMBER OF ATOMS
      <N>
      ITEM: BOX BOUNDS ...
      <3 lines>
      ITEM: ATOMS <cols...>
      <N atom lines>
    """
    stride_frames = max(1, stride_frames)
    n_yield_frames = 0
    frame_idx = 0

    with dump_path.open("r", errors="ignore") as f:
        while True:
            line = f.readline()
            if not line:
                break
            line = line.strip()
            if not line:
                continue

            _expect(line, "ITEM: TIMESTEP")
            t_line = f.readline()
            if not t_line:
                break
            timestep = int(t_line.strip())

            _expect(f.readline().strip(), "ITEM: NUMBER OF ATOMS")
            n_atoms = int(f.readline().strip())

            # box bounds header + 3 lines
            _expect(f.readline().strip(), "ITEM: BOX BOUNDS")
            f.readline()
            f.readline()
            f.readline()

            atoms_header = f.readline().strip()
            _expect(atoms_header, "ITEM: ATOMS")
            cols = atoms_header.split()[2:]

            # required columns
            try:
                mol_idx = cols.index("mol")
            except ValueError:
                mol_idx = cols.index("molecule") if "molecule" in cols else -1
            if mol_idx < 0:
                raise RuntimeError("Dump must include 'mol' column")
            type_idx = cols.index("type")

            # subsample frames
            frame_idx += 1
            keep = ((frame_idx - 1) % stride_frames == 0)

            # gather molecule -> flags (A/C/E/F)
            mol_flags: Dict[int, int] = {}
            # bitmask: 1=A, 2=C, 4=E, 8=F
            for _ in range(n_atoms):
                atom_line = f.readline()
                if not atom_line:
                    break
                if not keep:
                    continue
                parts = atom_line.split()
                if len(parts) < max(mol_idx, type_idx) + 1:
                    continue
                mol = int(parts[mol_idx])
                typ = int(parts[type_idx])
                if mol <= 0:
                    continue
                flag = mol_flags.get(mol, 0)
                if typ == TYPE_A:
                    flag |= 1
                elif typ == TYPE_C:
                    flag |= 2
                elif typ == TYPE_E:
                    flag |= 4
                elif typ == TYPE_F:
                    flag |= 8
                mol_flags[mol] = flag

            if not keep:
                continue

            # count unique molecules in each state
            nA = nC = nE = nF = 0
            for flag in mol_flags.values():
                # Prefer "post-flip" if both show up (shouldn't happen, but be safe)
                if flag & 2:
                    nC += 1
                elif flag & 1:
                    nA += 1

                if flag & 8:
                    nF += 1
                elif flag & 4:
                    nE += 1

            yield FrameCounts(timestep=timestep, nA=nA, nC=nC, nE=nE, nF=nF)
            n_yield_frames += 1
            if max_frames is not None and n_yield_frames >= max_frames:
                break


def maybe_plot_png(out_png: Path, frames: List[FrameCounts]) -> bool:
    try:
        import matplotlib.pyplot as plt  # type: ignore
    except Exception:
        return False

    ts = [fr.timestep for fr in frames]
    yA, yC, yE, yF = zip(*(fr.yields() for fr in frames)) if frames else ([], [], [], [])

    fig, ax = plt.subplots(1, 1, figsize=(10, 4))
    ax.plot(ts, yA, label="A", lw=1)
    ax.plot(ts, yC, label="C", lw=1)
    ax.plot(ts, yE, label="E", lw=1)
    ax.plot(ts, yF, label="F", lw=1)
    ax.set_xlabel("Timestep")
    ax.set_ylabel("Yield (normalized over A,C,E,F)")
    ax.set_ylim(-0.05, 1.05)
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_png, dpi=200)
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump", required=True, help="Path to .lammpstrj dump file")
    ap.add_argument("--out", required=True, help="Output prefix (e.g. analysis/yields_jobid)")
    ap.add_argument("--stride_frames", type=int, default=1, help="Keep every Nth frame (default: 1)")
    ap.add_argument("--max_frames", type=int, default=None, help="Stop after this many kept frames")
    args = ap.parse_args()

    dump_path = Path(args.dump)
    out_prefix = Path(args.out)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    frames = list(iter_frames(dump_path, stride_frames=args.stride_frames, max_frames=args.max_frames))

    out_csv = out_prefix.with_suffix(".csv")
    with out_csv.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["timestep", "N_A", "N_C", "N_E", "N_F", "yA", "yC", "yE", "yF"])
        for fr in frames:
            yA, yC, yE, yF = fr.yields()
            w.writerow([fr.timestep, fr.nA, fr.nC, fr.nE, fr.nF, yA, yC, yE, yF])

    out_png = out_prefix.with_suffix(".png")
    plotted = maybe_plot_png(out_png, frames)

    print(f"Wrote {out_csv}")
    if plotted:
        print(f"Wrote {out_png}")
    else:
        print("matplotlib not available; skipped plot (CSV still written)")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Compute and plot time-series yields (normalized concentrations) of monomers A, C, E, F
from a LAMMPS .lammpstrj dump.

Assumptions for the "twins" simulations:
- Molecules are identified by `mol` (molecule-ID) in the dump.
- Monomer identity is determined by patch atom types within a molecule:
    A: patch type 2
    C: patch type 4
    E: patch type 8
    F: patch type 10
  (B/D are ignored for this yield.)

Yield definition (per timestep):
  yA = nA / (nA + nC + nE + nF)
  yC = nC / (nA + nC + nE + nF)
  yE = nE / (nA + nC + nE + nF)
  yF = nF / (nA + nC + nE + nF)

Outputs:
- <out>.csv with per-timestep counts and yields
- <out>.png plot (if matplotlib is available)
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Tuple


PATCH_TO_LABEL = {
    2: "A",
    4: "C",
    8: "E",
    10: "F",
}


@dataclass
class FrameYields:
    step: int
    nA: int
    nC: int
    nE: int
    nF: int

    @property
    def total(self) -> int:
        return self.nA + self.nC + self.nE + self.nF

    def yields(self) -> Tuple[float, float, float, float]:
        tot = self.total
        if tot <= 0:
            return (0.0, 0.0, 0.0, 0.0)
        return (self.nA / tot, self.nC / tot, self.nE / tot, self.nF / tot)


def _readline(it: Iterator[str]) -> str:
    return next(it).rstrip("\n")


def iter_lammpstrj_frames(path: Path) -> Iterator[Tuple[int, List[Tuple[int, int]]]]:
    """
    Yield frames as (timestep, [(mol, type), ...]) pairs.
    Requires dump format containing at least: id mol type ...
    """
    with path.open("r", errors="ignore") as f:
        it = iter(f)
        while True:
            try:
                line = _readline(it)
            except StopIteration:
                return
            if not line:
                continue
            if not line.startswith("ITEM: TIMESTEP"):
                # Skip until next frame
                continue
            step = int(_readline(it).strip())

            # ITEM: NUMBER OF ATOMS
            _ = _readline(it)
            n_atoms = int(_readline(it).strip())

            # ITEM: BOX BOUNDS ...
            _ = _readline(it)
            _ = _readline(it)
            _ = _readline(it)
            _ = _readline(it)

            # ITEM: ATOMS ...
            atoms_header = _readline(it)
            if not atoms_header.startswith("ITEM: ATOMS"):
                raise RuntimeError(f"Unexpected ATOMS header: {atoms_header}")
            cols = atoms_header.split()[2:]
            try:
                mol_idx = cols.index("mol")
                type_idx = cols.index("type")
            except ValueError as e:
                raise RuntimeError(
                    f"Dump must contain 'mol' and 'type' columns. Found: {cols}"
                ) from e

            mol_type_pairs: List[Tuple[int, int]] = []
            mol_type_pairs.reserve if hasattr(mol_type_pairs, "reserve") else None  # no-op in py
            for _i in range(n_atoms):
                parts = _readline(it).split()
                mol = int(parts[mol_idx])
                typ = int(parts[type_idx])
                mol_type_pairs.append((mol, typ))

            yield step, mol_type_pairs


def compute_yields_timeseries(
    dump_path: Path,
    stride: int = 1,
    max_frames: Optional[int] = None,
) -> List[FrameYields]:
    out: List[FrameYields] = []
    frame_i = 0
    kept = 0
    for step, mol_type_pairs in iter_lammpstrj_frames(dump_path):
        if stride > 1 and (frame_i % stride != 0):
            frame_i += 1
            continue
        frame_i += 1

        # Determine label per molecule by seeing any patch type in PATCH_TO_LABEL
        mol_label: Dict[int, str] = {}
        for mol, typ in mol_type_pairs:
            lbl = PATCH_TO_LABEL.get(typ)
            if not lbl:
                continue
            # Each monomer should only ever match one of these labels; first match is enough.
            if mol not in mol_label:
                mol_label[mol] = lbl

        nA = sum(1 for v in mol_label.values() if v == "A")
        nC = sum(1 for v in mol_label.values() if v == "C")
        nE = sum(1 for v in mol_label.values() if v == "E")
        nF = sum(1 for v in mol_label.values() if v == "F")
        out.append(FrameYields(step=step, nA=nA, nC=nC, nE=nE, nF=nF))

        kept += 1
        if max_frames is not None and kept >= max_frames:
            break

    return out


def write_csv(rows: List[FrameYields], out_csv: Path) -> None:
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["timestep", "nA", "nC", "nE", "nF", "total", "yA", "yC", "yE", "yF"])
        for r in rows:
            yA, yC, yE, yF = r.yields()
            w.writerow([r.step, r.nA, r.nC, r.nE, r.nF, r.total, yA, yC, yE, yF])


def plot_png(rows: List[FrameYields], out_png: Path, title: str) -> bool:
    try:
        import matplotlib.pyplot as plt  # type: ignore
    except Exception:
        return False

    ts = [r.step for r in rows]
    yA = [r.yields()[0] for r in rows]
    yC = [r.yields()[1] for r in rows]
    yE = [r.yields()[2] for r in rows]
    yF = [r.yields()[3] for r in rows]

    plt.figure(figsize=(10, 5))
    plt.plot(ts, yA, label="A")
    plt.plot(ts, yC, label="C")
    plt.plot(ts, yE, label="E")
    plt.plot(ts, yF, label="F")
    plt.ylim(0.0, 1.0)
    plt.xlabel("Timestep")
    plt.ylabel("Yield (normalized over A,C,E,F)")
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_png, dpi=200)
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump", required=True, help="Path to .lammpstrj dump (must contain mol and type columns)")
    ap.add_argument("--out", required=True, help="Output prefix path (writes .csv and optionally .png)")
    ap.add_argument("--stride", type=int, default=1, help="Keep every Nth frame (default: 1 = keep all)")
    ap.add_argument("--max_frames", type=int, default=None, help="Optional cap on number of frames processed")
    ap.add_argument("--title", default="Yields vs time", help="Plot title")
    args = ap.parse_args()

    dump_path = Path(args.dump)
    out_prefix = Path(args.out)
    rows = compute_yields_timeseries(dump_path, stride=args.stride, max_frames=args.max_frames)

    out_csv = out_prefix.with_suffix(".csv")
    write_csv(rows, out_csv)
    print(f"Wrote {out_csv}")

    out_png = out_prefix.with_suffix(".png")
    if plot_png(rows, out_png, args.title):
        print(f"Wrote {out_png}")
    else:
        print("matplotlib not available; skipped plot (CSV still written)")


if __name__ == "__main__":
    main()


