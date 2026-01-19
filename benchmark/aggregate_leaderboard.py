#!/usr/bin/env python3
"""
Aggregate many one-row leaderboard CSVs into a single table for ranking.

Expected input: files produced by run_task.py, named like: *.leaderboard.csv
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Dict, List


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True, help="Directory containing *.leaderboard.csv files")
    ap.add_argument("--out", required=True, help="Output CSV path")
    ap.add_argument("--sort-by", default="score", help="Column to sort by (default: score)")
    ap.add_argument("--desc", action="store_true", help="Sort descending (default: ascending)")
    args = ap.parse_args()

    d = Path(args.dir)
    files = sorted(d.glob("*.leaderboard.csv"))
    if not files:
        raise SystemExit(f"No *.leaderboard.csv files found in {d}")

    rows: List[Dict[str, str]] = []
    for f in files:
        with f.open("r", newline="") as fh:
            r = csv.DictReader(fh)
            for rec in r:
                rec = dict(rec)
                rec["source_file"] = str(f)
                rows.append(rec)
                break

    # determine fieldnames (union)
    fieldnames: List[str] = []
    seen = set()
    for rec in rows:
        for k in rec.keys():
            if k not in seen:
                seen.add(k)
                fieldnames.append(k)

    def to_float(x: str) -> float:
        try:
            return float(x)
        except Exception:
            return float("-inf")

    sort_key = args.sort_by
    rows.sort(key=lambda rr: to_float(rr.get(sort_key, "")), reverse=args.desc)

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        for rec in rows:
            w.writerow(rec)

    print(f"Wrote {out_path} ({len(rows)} rows)")


if __name__ == "__main__":
    main()


