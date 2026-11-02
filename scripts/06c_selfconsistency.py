#!/usr/bin/env python3
"""
Stage 06c - two-seed POSE agreement.  Hard gate.

Challenge 03 found two seeds agreeing on interface confidence to within 0.03
while placing the binder 37 A apart. Score agreement is not pose agreement, and
a design whose predicted pose is seed-dependent is not a prediction at all.

Rather than port comp-03's 42_selfconsistency.py (which maps poses through a
receptor-renumbering table built for a different target), this compares the two
seeds on the quantity we actually care about here: WHICH TARGET RESIDUES the
binder touches, in construct numbering.

For each design, both seeds give a contact set on the target. We report

  jaccard   |A n B| / |A u B|    overlap of the two contact sets
  centroid  distance between the binder centre of mass in the two poses,
            after superposing on the target (the target is template-pinned to
            the same structure in both runs, so its frame is already common)

A design passes when the two seeds agree on where it binds:
  jaccard >= --min-jaccard  AND  centroid shift <= --max-shift

  venv-phscan/bin/python -I scripts/06c_selfconsistency.py
"""
import argparse
import csv
import glob
import os
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser

ROOT = Path(__file__).resolve().parent.parent
OFFSET = 121
CONTACT = 5.0


def pose(path, binder="X", target="A"):
    m = next(PDBParser(QUIET=True).get_structure("x", path).get_models())
    try:
        B = [r for r in m[binder] if r.id[0] == " "]
        T = {r.id[1]: r for r in m[target] if r.id[0] == " "}
    except KeyError:
        return None, None
    bx = np.array([a.coord for r in B for a in r])
    contacts = set()
    for i, r in T.items():
        tx = np.array([a.coord for a in r])
        if np.linalg.norm(bx[:, None, :] - tx[None, :, :], axis=-1).min() < CONTACT:
            contacts.add(i + OFFSET)
    return contacts, bx.mean(axis=0)


def find(root, name):
    h = glob.glob(str(root / "**" / f"{name}_model_0.pdb"), recursive=True)
    return h[0] if h else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed1", type=Path, nargs="+",
                    default=[ROOT / "boltz/refold/pred", ROOT / "boltz/refold_r2/pred"])
    ap.add_argument("--seed2", type=Path, default=ROOT / "boltz/selfconsist/pred")
    ap.add_argument("--names", type=Path, default=ROOT / "gen/selfconsist_input.fasta")
    ap.add_argument("--min-jaccard", type=float, default=0.50)
    ap.add_argument("--max-shift", type=float, default=5.0)
    ap.add_argument("--out", type=Path, default=ROOT / "results/selfconsistency.csv")
    args = ap.parse_args()

    names = [l[1:].strip() for l in open(args.names) if l.startswith(">")]
    rows = []
    for n in names:
        p1 = next((find(d, n) for d in args.seed1 if find(d, n)), None)
        p2 = find(args.seed2, n)
        if not p1 or not p2:
            rows.append(dict(name=n, jaccard=None, shift=None, n1=None, n2=None,
                             passes=False, reason="missing a seed"))
            continue
        c1, m1 = pose(p1)
        c2, m2 = pose(p2)
        if c1 is None or c2 is None:
            rows.append(dict(name=n, jaccard=None, shift=None, n1=None, n2=None,
                             passes=False, reason="unparseable"))
            continue
        u = len(c1 | c2)
        j = (len(c1 & c2) / u) if u else 0.0
        d = float(np.linalg.norm(m1 - m2))
        ok = j >= args.min_jaccard and d <= args.max_shift
        rows.append(dict(name=n, jaccard=round(j, 3), shift=round(d, 2),
                         n1=len(c1), n2=len(c2), passes=ok,
                         reason="" if ok else
                         (f"jaccard {j:.2f}" if j < args.min_jaccard
                          else f"shift {d:.1f}A")))

    rows.sort(key=lambda r: (not r["passes"], -(r["jaccard"] or 0)))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    npass = sum(r["passes"] for r in rows)
    print(f"[06c] {len(rows)} designs, {npass} agree across seeds "
          f"(jaccard >= {args.min_jaccard}, shift <= {args.max_shift} A)")
    big = [r for r in rows if r["shift"] and r["shift"] > 10]
    print(f"[06c] designs whose binder moved >10 A between seeds: {len(big)}")
    for r in big[:6]:
        print(f"        {r['name'][:46]:<48} shift {r['shift']:.1f} A  "
              f"jaccard {r['jaccard']}")
    print(f"[out] {args.out}")


if __name__ == "__main__":
    main()
