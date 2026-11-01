#!/usr/bin/env python3
"""
Stage 04e - epitope fidelity gate.  A HARD filter, not a ranking term.

Why this exists. Scoring against the full 639-residue monomer lets Boltz dock
the binder anywhere on the receptor, and it does. In the pilot the HIGHEST
scoring design (plain_d_1_s1, binder-target ipTM 0.775) hit zero of the five
targeted hotspots and instead parked on residues 317-320, 469-473, 636-648 and
680-689 -- i.e. on the N317 glycan sequon AND inside the transferrin site
(619-670). It would block receptor recycling, be glycan-occluded in the assay,
and forfeit the functional-epitope credit in objective 3. Its ipTM says none of
that.

So interface confidence is necessary and nowhere near sufficient. A design is
only a candidate if it binds WHERE WE AIMED IT.

Numbering: Boltz renumbers the target chain 1..639. Construct numbering is
resnum + 121. Everything reported here is CONSTRUCT numbering; the raw Boltz
index never leaves this script.

Gates (all must pass):
  - at least --min-hotspots of the targeted hotspots contacted
  - ZERO contacts inside the transferrin site 619-670  (objective-3 killer)
  - at most --max-glycan contacts in the glycan shadow (assay killer)
  - C-terminus not buried in the interface (the Twin-Strep capture tag
    attaches there and must not participate -- [OFFICIAL], and a silent
    disqualifier)

  venv-phscan/bin/python -I scripts/04e_epitope_qc.py --pred boltz/screen/pred
"""
import argparse
import csv
import glob
import json
import os
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser

ROOT = Path(__file__).resolve().parent.parent
OFFSET = 121                      # boltz target resi i  ==  construct i + OFFSET
TF_SITE = (619, 670)
CONTACT = 5.0                     # A, heavy-atom
CTERM_MIN = 8.0                   # A, C-term must sit this far from the target

HOTSPOT_SETS = {
    "shoulder": [209, 211, 348, 371, 374],
    "swap_locus": [211, 348, 371],
    "apex": [204, 208, 211, 212],
}


def batch_of(name: str) -> str:
    for key in ("shoulder", "apex", "swap"):
        if f"_{key}" in name or name.startswith(key):
            return key
    return "shoulder"


def hotspots_for(name: str):
    b = batch_of(name)
    if b == "apex":
        return HOTSPOT_SETS["apex"], "apex"
    if b == "swap":
        return HOTSPOT_SETS["swap_locus"], "swap_locus"
    return HOTSPOT_SETS["shoulder"], "shoulder"


def load_glycan_shadow():
    p = ROOT / "target/nondesignable.txt"
    if not p.exists():
        return set()
    return {int(t[1:]) for t in p.read_text().strip().split(",") if t.strip()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pred", type=Path, default=ROOT / "boltz/screen/pred")
    ap.add_argument("--out", type=Path, default=ROOT / "gen/epitope_qc.tsv")
    ap.add_argument("--min-hotspots", type=int, default=2)
    ap.add_argument("--max-glycan", type=int, default=2)
    ap.add_argument("--binder-chain", default="X")
    ap.add_argument("--target-chain", default="A")
    args = ap.parse_args()

    shadow = load_glycan_shadow()
    p = PDBParser(QUIET=True)
    files = sorted(glob.glob(str(args.pred / "**" / "*_model_0.pdb"),
                             recursive=True))
    print(f"[04e] {len(files)} predicted complexes; glycan shadow = "
          f"{len(shadow)} residues")

    rows = []
    for f in files:
        name = os.path.basename(f).replace("_model_0.pdb", "")
        try:
            m = next(p.get_structure("x", f).get_models())
            tgt = {r.id[1]: r for r in m[args.target_chain] if r.id[0] == " "}
            bres = [r for r in m[args.binder_chain] if r.id[0] == " "]
        except Exception as exc:
            print(f"  skip {name}: {exc}")
            continue
        if not bres or not tgt:
            continue
        bx = np.array([a.coord for r in bres for a in r])

        contacts = []
        for i, r in tgt.items():
            tx = np.array([a.coord for a in r])
            if np.linalg.norm(bx[:, None, :] - tx[None, :, :],
                              axis=-1).min() < CONTACT:
                contacts.append(i + OFFSET)
        contacts = sorted(contacts)

        hs, hs_name = hotspots_for(name)
        hit = sorted(set(contacts) & set(hs))
        tf = [c for c in contacts if TF_SITE[0] <= c <= TF_SITE[1]]
        gly = [c for c in contacts if c in shadow]

        # C-terminus exposure: the Twin-Strep tag attaches there
        cterm = bres[-1]
        cx = np.array([a.coord for a in cterm])
        allt = np.array([a.coord for r in tgt.values() for a in r])
        d_ct = float(np.linalg.norm(cx[:, None, :] - allt[None, :, :],
                                    axis=-1).min())

        reasons = []
        if len(hit) < args.min_hotspots:
            reasons.append(f"only {len(hit)}/{len(hs)} hotspots")
        if tf:
            reasons.append(f"transferrin site {tf[:4]}")
        if len(gly) > args.max_glycan:
            reasons.append(f"glycan shadow x{len(gly)}")
        if d_ct < CTERM_MIN:
            reasons.append(f"C-term {d_ct:.1f}A from target")

        rows.append(dict(
            name=name, hotspot_set=hs_name, n_contacts=len(contacts),
            n_hotspots=len(hit), hotspots=";".join(map(str, hit)),
            n_tf=len(tf), n_glycan=len(gly), cterm_dist=round(d_ct, 1),
            passes=(not reasons), reasons="; ".join(reasons),
            contacts=";".join(map(str, contacts))))

    rows.sort(key=lambda r: (-r["passes"], -r["n_hotspots"]))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    npass = sum(r["passes"] for r in rows)
    print(f"\n[04e] PASS {npass}/{len(rows)}")
    from collections import Counter
    why = Counter()
    for r in rows:
        if not r["passes"]:
            for part in r["reasons"].split("; "):
                why[part.split(" ")[0] if part.startswith("only")
                    else part.split(" ")[0]] += 1
    for k, v in why.most_common(8):
        print(f"       rejected, {k}: {v}")
    print(f"[out] {args.out}")


if __name__ == "__main__":
    main()
