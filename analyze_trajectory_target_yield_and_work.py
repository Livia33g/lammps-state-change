#!/usr/bin/env python3
"""
Trajectory-based monitoring of:
  1) "Yield" of a user-defined target structure (via cluster detection on molecules)
  2) "Work" as the energy kick proxy: ΔPotEng between consecutive thermo outputs

Motivation
----------
State-change simulations are out of equilibrium; each statechange can be viewed as an
energy injection/removal ("kick"). In practice, what you typically have on disk is:
  - A LAMMPS dump trajectory (positions + types)
  - A thermo stream in log/stdout (Step ... PotEng ...)
  - Optionally a statechange event stream in stderr (lines starting with "STATECHANGE ...")

This script joins those:
  - It computes yield from the DUMP at thermo timesteps (default), so you can plot yield vs time.
  - It computes work as ΔPotEng between consecutive thermo rows (your definition).
  - If an events file is provided, it counts statechange events in each thermo interval.

Target yield definition (cluster-based)
--------------------------------------
We build an undirected graph on molecules:
  - Each molecule is a node (molecule ids come from the dump's `mol` column).
  - An edge exists between two molecules if ANY pair of selected "site" atoms from the two
    molecules are within a distance cutoff (with PBC).
Then a "cluster" is a connected component in this molecule graph.

The "target structure" is defined by:
  - `--target-size K`: cluster size K (in molecules)
  - optional `--target-composition`: require cluster to have specific molecule labels
    (labels are derived from atom types in each molecule; see `--label-mode`)

Yield can be reported as:
  - fraction of molecules that are in target clusters (default)
  - number of target clusters

Species yield (common for state-change sims)
-------------------------------------------
If your "target" is simply "how many molecules are in state C now?", use:
  --yield-mode species_fraction --species-label <LABEL>

Where LABEL is the per-molecule label derived from site atom types (see --label-mode).
For example, in dimer_ksat_1core, A patches are type 2 and C patches are type 4, so:
  --yield-mode species_fraction --species-label 4

Example (dimer_ksat_1core):
--------------------------
python3 analyze_trajectory_target_yield_and_work.py \
  --dump dimer_ksat_1core_simulation_cpp/dump.dimer_ksat_1core.lammpstrj \
  --thermo dimer_ksat_1core_simulation_cpp/lammps_stdout.log \
  --events slurm_dimer_ksat_1core-15466594.err \
  --bond-cutoff 0.7 \
  --site-types 2 3 4 \
  --yield-mode species_fraction \
  --species-label 4 \
  --out analysis/dimer_ksat_1core_target_dimer
"""

from __future__ import annotations

import argparse
import bisect
import csv
import math
import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Set, Tuple


STATECHANGE_RE = re.compile(r"STATECHANGE\s+\S+:\s+step\s+(\d+)\b")


@dataclass(frozen=True)
class ThermoSeries:
    steps: List[int]
    pe: List[float]


def parse_statechange_steps(path: Path) -> List[int]:
    steps: List[int] = []
    for line in path.read_text(errors="ignore").splitlines():
        m = STATECHANGE_RE.search(line)
        if not m:
            continue
        steps.append(int(m.group(1)))
    steps.sort()
    return steps


def _try_parse_thermo_table(lines: List[str]) -> Optional[ThermoSeries]:
    """
    Parse a LAMMPS thermo table with a header containing Step and PotEng/pe.
    """
    header_idx = None
    header_cols: List[str] = []
    for i, line in enumerate(lines):
        if "Step" in line and ("PotEng" in line or "pe" in line.lower()):
            header_idx = i
            header_cols = line.split()
            break
    if header_idx is None:
        return None

    # Step column
    try:
        step_col = header_cols.index("Step")
    except ValueError:
        return None

    # Potential energy column (robust to case)
    pe_col = None
    for cand in ("PotEng", "PE", "pe"):
        if cand in header_cols:
            pe_col = header_cols.index(cand)
            break
    if pe_col is None:
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
        raise RuntimeError(f"Could not parse thermo table (Step + PotEng/pe) from: {path}")
    return series


@dataclass(frozen=True)
class DumpFrame:
    timestep: int
    box: Tuple[float, float, float]  # Lx, Ly, Lz
    # Atom records: (mol_id, type, x, y, z)
    atoms: List[Tuple[int, int, float, float, float]]


def iter_lammpstrj_frames(
    dump_path: Path, *, keep_types: Optional[Set[int]] = None
) -> Iterator[DumpFrame]:
    """
    Streaming reader for LAMMPS 'dump custom' style with sections:
      ITEM: TIMESTEP
      <t>
      ITEM: NUMBER OF ATOMS
      <n>
      ITEM: BOX BOUNDS ...
      <xlo xhi>
      <ylo yhi>
      <zlo zhi>
      ITEM: ATOMS id mol type x y z ...

    We only require that the ATOMS header includes: id mol type x y z (in that order).
    Extra columns are ignored.
    """
    with dump_path.open("r", errors="ignore") as f:
        while True:
            line = f.readline()
            if not line:
                return
            if not line.startswith("ITEM: TIMESTEP"):
                continue
            t_line = f.readline()
            if not t_line:
                return
            timestep = int(t_line.strip())

            # number of atoms
            if not f.readline().startswith("ITEM: NUMBER OF ATOMS"):
                raise RuntimeError("Unexpected dump format (missing NUMBER OF ATOMS)")
            n = int(f.readline().strip())

            # box bounds
            box_hdr = f.readline()
            if not box_hdr.startswith("ITEM: BOX BOUNDS"):
                raise RuntimeError("Unexpected dump format (missing BOX BOUNDS)")
            xlo, xhi = [float(x) for x in f.readline().split()[:2]]
            ylo, yhi = [float(x) for x in f.readline().split()[:2]]
            zlo, zhi = [float(x) for x in f.readline().split()[:2]]
            lx, ly, lz = (xhi - xlo, yhi - ylo, zhi - zlo)

            atoms_hdr = f.readline()
            if not atoms_hdr.startswith("ITEM: ATOMS"):
                raise RuntimeError("Unexpected dump format (missing ATOMS header)")
            cols = atoms_hdr.split()[2:]  # after "ITEM: ATOMS"
            try:
                id_idx = cols.index("id")
                mol_idx = cols.index("mol")
                typ_idx = cols.index("type")
                x_idx = cols.index("x")
                y_idx = cols.index("y")
                z_idx = cols.index("z")
            except ValueError as e:
                raise RuntimeError(f"Dump ATOMS header must include id mol type x y z. Got: {cols}") from e

            atoms: List[Tuple[int, int, float, float, float]] = []
            for _ in range(n):
                row = f.readline()
                if not row:
                    raise RuntimeError("Truncated dump file")
                parts = row.split()
                mol = int(parts[mol_idx])
                typ = int(parts[typ_idx])
                if keep_types is not None and typ not in keep_types:
                    continue
                x = float(parts[x_idx])
                y = float(parts[y_idx])
                z = float(parts[z_idx])
                atoms.append((mol, typ, x, y, z))

            yield DumpFrame(timestep=timestep, box=(lx, ly, lz), atoms=atoms)


def _min_image(d: float, box_l: float) -> float:
    if box_l <= 0.0:
        return d
    # round to nearest integer number of box lengths
    return d - box_l * round(d / box_l)


def build_molecule_labels(
    frame: DumpFrame, *, label_mode: str = "majority_site_type", site_types: Optional[Set[int]] = None
) -> Dict[int, int]:
    """
    Produce mol_id -> label (int) for this frame.

    label_mode:
      - "majority_site_type": most common atom type among site_types within the molecule
      - "min_site_type": minimum atom type among site_types within the molecule
      - "max_site_type": maximum atom type among site_types within the molecule

    If a molecule has no atoms matching site_types, it gets label 0.
    """
    if site_types is None:
        # if no site types given, label from all atom types
        site_types = set()

    mol_to_types: Dict[int, List[int]] = defaultdict(list)
    for mol, typ, *_ in frame.atoms:
        if site_types and typ not in site_types:
            continue
        mol_to_types[mol].append(typ)

    labels: Dict[int, int] = {}
    for mol, types in mol_to_types.items():
        if not types:
            labels[mol] = 0
            continue
        if label_mode == "majority_site_type":
            labels[mol] = Counter(types).most_common(1)[0][0]
        elif label_mode == "min_site_type":
            labels[mol] = min(types)
        elif label_mode == "max_site_type":
            labels[mol] = max(types)
        else:
            raise ValueError(f"Unknown label_mode: {label_mode}")

    # ensure all mols in frame are present
    for mol, *_ in frame.atoms:
        labels.setdefault(mol, 0)
    return labels


class UnionFind:
    def __init__(self, items: Iterable[int]):
        self.parent: Dict[int, int] = {x: x for x in items}
        self.rank: Dict[int, int] = {x: 0 for x in items}

    def find(self, x: int) -> int:
        p = self.parent.get(x, x)
        if p != x:
            self.parent[x] = self.find(p)
        return self.parent.get(x, x)

    def union(self, a: int, b: int) -> None:
        ra = self.find(a)
        rb = self.find(b)
        if ra == rb:
            return
        if self.rank[ra] < self.rank[rb]:
            ra, rb = rb, ra
        self.parent[rb] = ra
        if self.rank[ra] == self.rank[rb]:
            self.rank[ra] += 1


def _cell_index(x: float, box_l: float, cell: float) -> int:
    # map coordinate to [0, ncell) then to int cell index
    ncell = max(1, int(math.floor(box_l / cell)))
    # shift to [0, box_l)
    u = x % box_l if box_l > 0 else x
    return int(math.floor(u / box_l * ncell)) if box_l > 0 else int(math.floor(u / cell))


def molecule_clusters_from_sites(
    frame: DumpFrame, *, site_types: Set[int], cutoff: float
) -> Dict[int, Set[int]]:
    """
    Build molecule clusters using site-site distances with PBC.

    Returns dict root_mol -> set(mol_ids) for each connected component.
    """
    cutoff2 = cutoff * cutoff
    lx, ly, lz = frame.box

    # collect site atoms
    sites: List[Tuple[int, float, float, float]] = []
    mols: Set[int] = set()
    for mol, typ, x, y, z in frame.atoms:
        mols.add(mol)
        if typ in site_types:
            sites.append((mol, x, y, z))

    uf = UnionFind(mols)
    if not sites:
        return {m: {m} for m in mols}

    # cell list for speed
    cell = cutoff
    nx = max(1, int(math.floor(lx / cell))) if lx > 0 else 1
    ny = max(1, int(math.floor(ly / cell))) if ly > 0 else 1
    nz = max(1, int(math.floor(lz / cell))) if lz > 0 else 1

    grid: Dict[Tuple[int, int, int], List[int]] = defaultdict(list)
    xs = [s[1] for s in sites]
    ys = [s[2] for s in sites]
    zs = [s[3] for s in sites]

    for idx, (_, x, y, z) in enumerate(sites):
        ix = int(((x % lx) / lx) * nx) if lx > 0 else int(math.floor(x / cell))
        iy = int(((y % ly) / ly) * ny) if ly > 0 else int(math.floor(y / cell))
        iz = int(((z % lz) / lz) * nz) if lz > 0 else int(math.floor(z / cell))
        grid[(ix, iy, iz)].append(idx)

    def wrap(i: int, n: int) -> int:
        return i % n if n > 0 else i

    # neighbor offsets (-1,0,1)
    for i, (moli, xi, yi, zi) in enumerate(sites):
        ix = int(((xi % lx) / lx) * nx) if lx > 0 else int(math.floor(xi / cell))
        iy = int(((yi % ly) / ly) * ny) if ly > 0 else int(math.floor(yi / cell))
        iz = int(((zi % lz) / lz) * nz) if lz > 0 else int(math.floor(zi / cell))
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    key = (wrap(ix + dx, nx), wrap(iy + dy, ny), wrap(iz + dz, nz))
                    for j in grid.get(key, []):
                        if j <= i:
                            continue
                        molj, xj, yj, zj = sites[j]
                        if molj == moli:
                            continue
                        ddx = _min_image(xj - xi, lx)
                        ddy = _min_image(yj - yi, ly)
                        ddz = _min_image(zj - zi, lz)
                        rsq = ddx * ddx + ddy * ddy + ddz * ddz
                        if rsq <= cutoff2:
                            uf.union(moli, molj)

    comps: Dict[int, Set[int]] = defaultdict(set)
    for m in mols:
        comps[uf.find(m)].add(m)
    return dict(comps)


def parse_target_composition(spec: Optional[str]) -> Optional[Dict[int, int]]:
    """
    Parse like "2:1,3:1" meaning: label 2 count 1, label 3 count 1.
    """
    if not spec:
        return None
    out: Dict[int, int] = {}
    for part in spec.split(","):
        part = part.strip()
        if not part:
            continue
        if ":" not in part:
            raise ValueError(f"Bad --target-composition entry '{part}', expected 'label:count'")
        k, v = part.split(":", 1)
        out[int(k.strip())] = int(v.strip())
    return out or None


def cluster_matches_composition(cluster: Set[int], labels: Dict[int, int], target: Dict[int, int]) -> bool:
    c = Counter(labels.get(m, 0) for m in cluster)
    return all(c.get(lbl, 0) == n for lbl, n in target.items())


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump", required=True, help="LAMMPS dump (.lammpstrj) path")
    ap.add_argument("--thermo", required=True, help="LAMMPS thermo output (log.lammps or lammps_stdout.log)")
    ap.add_argument("--events", default=None, help="Optional file containing STATECHANGE lines (e.g. slurm .err)")
    ap.add_argument("--out", required=True, help="Output prefix (writes .csv and .png)")

    ap.add_argument("--site-types", nargs="+", type=int, required=True, help="Atom types to use as bonding sites")
    ap.add_argument("--bond-cutoff", type=float, required=True, help="Distance cutoff for site-site bond")

    ap.add_argument("--target-size", type=int, default=None, help="Target cluster size in molecules (cluster yield modes)")
    ap.add_argument(
        "--target-composition",
        default=None,
        help="Optional exact label composition for target clusters, e.g. '2:1,3:1'",
    )
    ap.add_argument(
        "--label-mode",
        choices=["majority_site_type", "min_site_type", "max_site_type"],
        default="majority_site_type",
        help="How to label each molecule (for composition constraints)",
    )
    ap.add_argument(
        "--yield-mode",
        choices=["fraction_molecules", "n_clusters", "species_fraction"],
        default="fraction_molecules",
        help="Yield type: cluster-based (fraction_molecules, n_clusters) or per-species (species_fraction)",
    )
    ap.add_argument(
        "--species-label",
        type=int,
        default=None,
        help="Required for --yield-mode species_fraction: molecule label to count (e.g. 4 for C)",
    )
    ap.add_argument(
        "--sample",
        choices=["thermo", "all"],
        default="thermo",
        help="Compute yield at thermo steps only (fast) or at every dump frame",
    )
    ap.add_argument(
        "--max-frames",
        type=int,
        default=None,
        help="Optional cap on number of dump frames processed (debugging/smoke tests)",
    )
    args = ap.parse_args()

    dump_path = Path(args.dump)
    thermo_path = Path(args.thermo)
    events_path = Path(args.events) if args.events else None
    out_prefix = Path(args.out)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    series = parse_thermo_pe(thermo_path)
    thermo_steps = series.steps
    pe_by_step: Dict[int, float] = {s: pe for s, pe in zip(series.steps, series.pe)}

    event_steps: List[int] = []
    if events_path is not None and events_path.exists():
        event_steps = parse_statechange_steps(events_path)

    # For each thermo row i>0, count number of events in (prev_step, step]
    events_in_interval: Dict[int, int] = {}
    if event_steps:
        idx = 0
        for i in range(1, len(thermo_steps)):
            prev_s = thermo_steps[i - 1]
            s = thermo_steps[i]
            count = 0
            while idx < len(event_steps) and event_steps[idx] <= s:
                if event_steps[idx] > prev_s:
                    count += 1
                idx += 1
            events_in_interval[s] = count

    # Decide which dump frames we care about
    sample_steps: Optional[Set[int]] = None
    if args.sample == "thermo":
        sample_steps = set(thermo_steps)

    keep_types = set(args.site_types)
    target_comp = parse_target_composition(args.target_composition)

    # We write rows keyed by timestep for later join with thermo
    yield_by_step: Dict[int, Tuple[float, int, int]] = {}  # step -> (yield, n_clusters, n_mols)

    n_processed = 0
    for frame in iter_lammpstrj_frames(dump_path, keep_types=keep_types):
        if sample_steps is not None and frame.timestep not in sample_steps:
            continue
        labels = build_molecule_labels(frame, label_mode=args.label_mode, site_types=keep_types)

        n_mols = len({mol for mol, *_ in frame.atoms})
        n_target_clusters = 0

        if args.yield_mode == "species_fraction":
            if args.species_label is None:
                raise RuntimeError("--yield-mode species_fraction requires --species-label")
            n_species = sum(1 for _, lbl in labels.items() if lbl == args.species_label)
            y = (n_species / n_mols) if n_mols > 0 else 0.0
        else:
            if args.target_size is None or args.target_size < 1:
                raise RuntimeError("Cluster yield modes require --target-size >= 1")
            comps = molecule_clusters_from_sites(frame, site_types=keep_types, cutoff=args.bond_cutoff)

            # Determine which clusters match target
            mols_in_target = 0
            for cluster in comps.values():
                if len(cluster) != args.target_size:
                    continue
                if target_comp is not None and not cluster_matches_composition(cluster, labels, target_comp):
                    continue
                n_target_clusters += 1
                mols_in_target += len(cluster)

            if args.yield_mode == "fraction_molecules":
                y = (mols_in_target / n_mols) if n_mols > 0 else 0.0
            else:
                y = float(n_target_clusters)

        yield_by_step[frame.timestep] = (y, n_target_clusters, n_mols)

        n_processed += 1
        if args.max_frames is not None and n_processed >= args.max_frames:
            break

    # Build output table at thermo steps (always), joining yield (if missing -> NaN)
    rows_out: List[Dict[str, object]] = []
    cum_work = 0.0
    for i, s in enumerate(thermo_steps):
        pe = pe_by_step[s]
        if i == 0:
            dpe = 0.0
        else:
            dpe = pe - pe_by_step[thermo_steps[i - 1]]
            cum_work += dpe

        y, ncl, nm = yield_by_step.get(s, (float("nan"), 0, 0))
        ev = events_in_interval.get(s, 0) if i > 0 else 0

        rows_out.append(
            {
                "timestep": s,
                "pe": pe,
                "dpe": dpe,
                "cum_work": cum_work,
                "yield": y,
                "n_target_clusters": ncl,
                "n_molecules": nm,
                "n_statechanges_interval": ev,
            }
        )

    out_csv = out_prefix.with_suffix(".csv")
    with out_csv.open("w", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=[
                "timestep",
                "pe",
                "dpe",
                "cum_work",
                "yield",
                "n_target_clusters",
                "n_molecules",
                "n_statechanges_interval",
            ],
        )
        w.writeheader()
        for r in rows_out:
            w.writerow(r)

    # Plot (optional)
    out_png = out_prefix.with_suffix(".png")
    try:
        import matplotlib.pyplot as plt  # type: ignore
    except Exception:
        print(f"Wrote {out_csv}")
        print("matplotlib not available; skipped plot (CSV still written)")
        return

    ts = [int(r["timestep"]) for r in rows_out]
    ys = [float(r["yield"]) for r in rows_out]
    dpes = [float(r["dpe"]) for r in rows_out]
    cumw = [float(r["cum_work"]) for r in rows_out]
    evs = [int(r["n_statechanges_interval"]) for r in rows_out]

    fig, ax = plt.subplots(3, 1, figsize=(11, 8), sharex=True)
    ax[0].plot(ts, ys, lw=1.5)
    ax[0].set_ylabel("Target yield")
    ax[0].set_title("Target yield and work (ΔPotEng) vs time")

    ax[1].plot(ts, dpes, lw=1.0)
    ax[1].axhline(0.0, color="k", lw=0.5)
    ax[1].set_ylabel("Work proxy: ΔPE\n(per thermo interval)")

    ax[2].plot(ts, cumw, lw=1.0)
    ax[2].set_ylabel("Cumulative ΔPE")
    ax[2].set_xlabel("Timestep")

    # If events exist, annotate with a light bar overlay on ax[1]
    if any(evs):
        ax1b = ax[1].twinx()
        ax1b.step(ts, evs, where="post", color="tab:orange", alpha=0.35)
        ax1b.set_ylabel("# statechanges\nin interval", color="tab:orange")
        ax1b.tick_params(axis="y", labelcolor="tab:orange")

    fig.tight_layout()
    fig.savefig(out_png, dpi=200)

    print(f"Wrote {out_csv}")
    print(f"Wrote {out_png}")


if __name__ == "__main__":
    main()


