#!/usr/bin/env python3
"""
Stage 36 - build a PER-POSITION histidine bias aimed at the conserved cations.

Why. Round 1 used ProteinMPNN's global His bias, which raises His everywhere on
the binder. Measured outcome across 153 designs with a usable switch:

    species_specific  142      (R208 / K371, strong: up to 2.21 kcal/mol = 42x)
    cross_species      11      (K374, weak: 0.29-1.15 kcal/mol = 2-7x)

and ten of those eleven are a SINGLE histidine. That is the signature of luck,
not design: only one conserved cation (K374) is available, a global bias
scatters His across the whole binder, and a His lands near K374 only
occasionally. Meanwhile R208 and K371 are more exposed and more abundant, so
the global bias preferentially builds the switch we cannot use cross-species.

Objectives 1 and 2 are the top two ranked criteria and only a conserved-cation
switch serves both. So bias His at exactly the binder positions that FACE a
conserved cation, and nowhere else.

Conserved, non-vetoed cations in the generation crop (construct numbering):
    K374  crop index 139   the anchor: conserved K->K376, no carboxylate <10 A
    K224  crop index  34   conserved K->K226, nearest carboxylate E272 @ 8.4 A
K344 is excluded: conserved, but salt-bridged to E214 at 2.8 A, which would
inverT the switch (stage 00). K231 is excluded: 8.7 A from the N251 sequon.

Numbering. RFdiffusion writes binder as chain A renumbered from 1 and the
cropped target as chain B starting immediately after, so chain B residue
= binder_length + crop_index. Rather than rely on that arithmetic, the target
residue is taken POSITIONALLY as the Nth residue of chain B, which holds
whatever the binder length is.

ProteinMPNN alphabet order is 'ACDEFGHIKLMNPQRSTVWYX'; His is index 6.

  venv-phscan/bin/python -I scripts/36_targeted_his_bias.py
"""
import argparse
import json
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser

ROOT = Path(__file__).resolve().parent.parent
ALPHABET = "ACDEFGHIKLMNPQRSTVWYX"
HIS = ALPHABET.index("H")

# (construct position, crop index) for conserved, non-vetoed cations
CONSERVED_CATIONS = [(374, 139), (224, 34)]
RADIUS = 10.0          # A, binder CB within this of the cation -> bias His
BIAS = 2.5             # logit bump; global pool used 1.5 everywhere


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--backbones", type=Path, default=ROOT / "gen/r2_backbones")
    ap.add_argument("--out", type=Path, default=ROOT / "gen/bias_targeted.jsonl")
    ap.add_argument("--radius", type=float, default=RADIUS)
    ap.add_argument("--bias", type=float, default=BIAS)
    ap.add_argument("--report", type=Path, default=ROOT / "gen/bias_targeted_report.tsv")
    args = ap.parse_args()

    p = PDBParser(QUIET=True)
    out, rows = {}, []
    for pdb in sorted(args.backbones.glob("*.pdb")):
        try:
            m = next(p.get_structure("x", pdb).get_models())
            A = [r for r in m["A"] if r.id[0] == " "]
            B = [r for r in m["B"] if r.id[0] == " "]
        except Exception as exc:
            print(f"  skip {pdb.name}: {exc}")
            continue
        if not A or len(B) < 140:
            print(f"  skip {pdb.name}: chains A={len(A)} B={len(B)}")
            continue

        # positional lookup of each conserved cation in the cropped target
        targets = []
        for construct, idx in CONSERVED_CATIONS:
            if idx - 1 < len(B):
                r = B[idx - 1]
                at = r["CB"] if "CB" in r else (r["CA"] if "CA" in r else None)
                if at is not None:
                    targets.append((construct, np.array(at.coord)))
        if not targets:
            continue

        bias = np.zeros([len(A), 21], dtype=float)
        hits = []
        for i, r in enumerate(A):
            at = r["CB"] if "CB" in r else (r["CA"] if "CA" in r else None)
            if at is None:
                continue
            c = np.array(at.coord)
            near = [(cn, float(np.linalg.norm(c - t))) for cn, t in targets]
            best = min(near, key=lambda x: x[1])
            if best[1] <= args.radius:
                bias[i, HIS] = args.bias
                hits.append((i + 1, best[0], round(best[1], 1)))

        out[pdb.stem] = {"A": bias.tolist(),
                         "B": np.zeros([len(B), 21]).tolist()}
        rows.append((pdb.stem, len(A), len(hits),
                     ";".join(f"{i}->{cn}@{d}" for i, cn, d in hits)))

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w") as f:
        f.write(json.dumps(out) + "\n")
    with open(args.report, "w") as f:
        f.write("backbone\tbinder_len\tn_biased\tpositions\n")
        for r in rows:
            f.write("\t".join(map(str, r)) + "\n")

    n = [r[2] for r in rows]
    print(f"[36] {len(rows)} backbones biased")
    print(f"[36] binder positions facing a conserved cation: "
          f"min {min(n)}  median {int(np.median(n))}  max {max(n)}  "
          f"total {sum(n)}")
    print(f"[36] backbones with ZERO facing positions: {sum(1 for x in n if x == 0)}")
    print(f"[out] {args.out}")
    print(f"[out] {args.report}")


if __name__ == "__main__":
    main()
