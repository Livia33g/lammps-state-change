#!/usr/bin/env python3
import argparse, csv, math

BIND_TYPES = {"A_RIGHT": 3, "B_LEFT": 4, "B_RIGHT": 5, "C_LEFT": 6}

def min_image(dx, L):
    if L <= 0: return dx
    h = 0.5*L
    if dx >  h: dx -= L
    if dx < -h: dx += L
    return dx

def dist_pbc(a, b, box):
    (xlo,xhi),(ylo,yhi),(zlo,zhi) = box
    dx = min_image(a[0]-b[0], xhi-xlo)
    dy = min_image(a[1]-b[1], yhi-ylo)
    dz = min_image(a[2]-b[2], zhi-zlo)
    return math.sqrt(dx*dx + dy*dy + dz*dz)

def parse_dump(path):
    with open(path, "r") as f:
        while True:
            line = f.readline()
            if not line: break
            if not line.startswith("ITEM: TIMESTEP"): continue
            ts = int(f.readline().strip())

            # number of atoms
            assert f.readline().startswith("ITEM: NUMBER OF ATOMS")
            natoms = int(f.readline().strip())

            # box bounds (assumes periodic)
            assert f.readline().startswith("ITEM: BOX BOUNDS")
            xlo, xhi = map(float, f.readline().split()[:2])
            ylo, yhi = map(float, f.readline().split()[:2])
            zlo, zhi = map(float, f.readline().split()[:2])
            box = ((xlo,xhi),(ylo,yhi),(zlo,zhi))

            # atoms header (must include id mol type x y z)
            atoms_hdr = f.readline().strip()
            assert atoms_hdr.startswith("ITEM: ATOMS")
            cols = atoms_hdr.split()[2:]
            idx = {c:i for i,c in enumerate(cols)}
            required = ["id","mol","type","x","y","z"]
            if any(r not in idx for r in required):
                raise RuntimeError(f"Dump must include columns {required}; found {cols}")

            atoms = []
            for _ in range(natoms):
                p = f.readline().split()
                atoms.append({
                    "id":   int(p[idx["id"]]),
                    "mol":  int(p[idx["mol"]]),
                    "type": int(p[idx["type"]]),
                    "xyz":  (float(p[idx["x"]]), float(p[idx["y"]]), float(p[idx["z"]]))
                })
            yield ts, box, atoms

def count_frame(box, atoms, cutoff):
    # collect binding patch atoms
    A3 = [a for a in atoms if a["type"] == BIND_TYPES["A_RIGHT"]]
    B4 = [a for a in atoms if a["type"] == BIND_TYPES["B_LEFT" ]]
    B5 = [a for a in atoms if a["type"] == BIND_TYPES["B_RIGHT"]]
    C6 = [a for a in atoms if a["type"] == BIND_TYPES["C_LEFT" ]]

    # AB edges via 3–4 within cutoff; unique by (A_mol,B_mol)
    AB = set()
    for a in A3:
        for b in B4:
            if a["mol"] == b["mol"]: continue
            if dist_pbc(a["xyz"], b["xyz"], box) <= cutoff:
                AB.add((a["mol"], b["mol"]))

    # BC edges via 5–6 within cutoff; unique by (B_mol,C_mol)
    BC = set()
    for b in B5:
        for c in C6:
            if b["mol"] == c["mol"]: continue
            if dist_pbc(b["xyz"], c["xyz"], box) <= cutoff:
                BC.add((b["mol"], c["mol"]))

    # neighbor maps by B
    AB_by_B = {}
    for a,b in AB:
        AB_by_B.setdefault(b,set()).add(a)
    BC_by_B = {}
    for b,c in BC:
        BC_by_B.setdefault(b,set()).add(c)

    # reverse maps (for "pure" exclusivity)
    B_by_A = {}
    for a,b in AB:
        B_by_A.setdefault(a,set()).add(b)
    B_by_C = {}
    for b,c in BC:
        B_by_C.setdefault(c,set()).add(b)

    # all A–B–C triplets (cartesian product per shared B)
    triples_all = set()
    for b in set(AB_by_B) & set(BC_by_B):
        for a in AB_by_B[b]:
            for c in BC_by_B[b]:
                triples_all.add((a,b,c))

    # "clean": each B has exactly one A and one C neighbor
    triples_clean = set()
    for b in set(AB_by_B) & set(BC_by_B):
        if len(AB_by_B[b]) == 1 and len(BC_by_B[b]) == 1:
            a = next(iter(AB_by_B[b]))
            c = next(iter(BC_by_B[b]))
            triples_clean.add((a,b,c))

    # "pure": clean + A and C are not attached to any other B
    triples_pure = set()
    for (a,b,c) in triples_clean:
        if len(B_by_A.get(a,set())) == 1 and len(B_by_C.get(c,set())) == 1:
            triples_pure.add((a,b,c))

    return {
        "nAB_edges": len(AB),
        "nBC_edges": len(BC),
        "nABC_all": len(triples_all),
        "nABC_clean": len(triples_clean),
        "nABC_pure": len(triples_pure),
        "nApatch3": len(A3), "nBpatch4": len(B4),
        "nBpatch5": len(B5), "nCpatch6": len(C6)
    }

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dump")
    ap.add_argument("--cutoff", type=float, default=0.28)
    ap.add_argument("--out", default="abc_counts.csv")
    args = ap.parse_args()

    rows = []
    for ts, box, atoms in parse_dump(args.dump):
        row = count_frame(box, atoms, args.cutoff)
        row["timestep"] = ts
        rows.append(row)

    # sort by time and write CSV
    rows.sort(key=lambda r: r["timestep"])
    fields = ["timestep","nABC_pure","nABC_clean","nABC_all","nAB_edges","nBC_edges",
          "nApatch3","nBpatch4","nBpatch5","nCpatch6"]
    # fields = ["timestep","nAB_edges","nBC_edges","nABC_all","nABC_clean","nABC_pure",
    #           "nApatch3","nBpatch4","nBpatch5","nCpatch6"]
    with open(args.out, "w", newline="") as fp:
        w = csv.DictWriter(fp, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"[done] wrote {args.out} with {len(rows)} frames.")

if __name__ == "__main__":
    main()
