#!/usr/bin/env python3
"""
Design self-consistency: does the predicted complex reproduce the designed one?

Boltz re-predicts the binder pose rather than templating it (only the receptor
is pinned). So a high interface confidence can mean either "this design binds as
intended" or "Boltz found some other pose it likes". The brief's own benchmark
warns that Boltz-2 confidence overestimates correctness for misplaced binders,
so confidence alone cannot distinguish those.

Self-consistency separates them. Superpose the predicted receptor onto the
design's receptor, then ask whether the binder ends up in the same place:

  binder_rmsd      Ca RMSD of the binder after receptor superposition. Low means
                   the prediction agrees with the design.
  contact_jaccard  overlap of the receptor residues contacted, designed vs
                   predicted. Robust to small pose shifts, and directly answers
                   "is it still on the intended epitope".

A design that scores well on Boltz but moves its binder tens of angstroms is not
a validated design; it is a different complex that happens to score well.

Usage:
  python 42_selfconsistency.py --boltz gen/boltz --designs gen/rfdiff \
      --scores gen/boltz/boltz_scores.csv --out gen/boltz/selfconsistency.csv
"""

import argparse
import csv
import glob
import os
import pickle
import sys

import numpy as np

CONTACT_A = 5.0


def chains(pdb):
    ch = {}
    for l in open(pdb):
        if l[:4] != "ATOM":
            continue
        ch.setdefault(l[21], {}).setdefault(int(l[22:26]), []).append(
            (l[12:16].strip(),
             np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])))
    return ch


def ca(resmap):
    out = {}
    for n, atoms in resmap.items():
        for an, xyz in atoms:
            if an == "CA":
                out[n] = xyz
                break
    return out


def kabsch(P, Q):
    """Rotation+translation mapping P onto Q."""
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, qc - R @ pc


def contacts(binder_atoms, target_resmap, cutoff=CONTACT_A):
    b = [np.array([x for _, x in a]) for a in binder_atoms]
    out = set()
    for tn, atoms in target_resmap.items():
        tv = np.array([x for _, x in atoms])
        if min(np.linalg.norm(x[:, None, :] - tv[None, :, :], axis=-1).min()
               for x in b) < cutoff:
            out.add(tn)
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--boltz", default="gen/boltz")
    p.add_argument("--designs", default="gen/rfdiff")
    p.add_argument("--scores", default="gen/boltz/boltz_scores.csv")
    p.add_argument("--out", default="gen/boltz/selfconsistency.csv")
    # PORTED from challenge 03. The target PDB was hardcoded there; it is an
    # argument here so the apical-crop and full-monomer scoring modes can both
    # be mapped back to author numbering.
    p.add_argument("--target-pdb", default="target/tfr1_A.pdb")
    args = p.parse_args()

    # boltz receptor is renumbered 1..N; map to original receptor numbering
    seqmap, i = {}, 0
    for l in open(args.target_pdb):
        if l[:4] == "ATOM" and l[12:16].strip() == "CA":
            i += 1
            seqmap[i] = int(l[22:26])

    rows = list(csv.DictReader(open(args.scores)))
    out = []
    for r in rows:
        name = r["design"]
        parent = name.rsplit("_s", 1)[0]
        hits = glob.glob(os.path.join(args.boltz, "**", f"{name}_model_0.pdb"),
                         recursive=True)
        dpdb = os.path.join(args.designs, parent + ".pdb")
        if not hits or not os.path.isfile(dpdb):
            continue

        pred = chains(hits[0])
        des = chains(dpdb)
        if "A" not in pred or "B" not in pred or "A" not in des or "B" not in des:
            continue

        # design target numbering -> original, via the .trb
        trb = dpdb[:-4] + ".trb"
        if not os.path.isfile(trb):
            continue
        d = pickle.load(open(trb, "rb"))
        tmap = {h[1]: rr[1] for rr, h in zip(d["receptor_con_ref_pdb_idx"],
                                             d["receptor_con_hal_pdb_idx"])}

        des_t_ca = {tmap[n]: c for n, c in ca(des["B"]).items() if n in tmap}
        pred_t_ca = {seqmap.get(n, n): c for n, c in ca(pred["B"]).items()}
        shared = sorted(set(des_t_ca) & set(pred_t_ca))
        if len(shared) < 20:
            continue

        P = np.array([pred_t_ca[n] for n in shared])
        Q = np.array([des_t_ca[n] for n in shared])
        R, t = kabsch(P, Q)

        des_b = ca(des["A"])
        pred_b = ca(pred["A"])
        common = sorted(set(des_b) & set(pred_b))
        if not common:
            continue
        pb = np.array([R @ pred_b[n] + t for n in common])
        db = np.array([des_b[n] for n in common])
        brmsd = float(np.sqrt(((pb - db) ** 2).sum(1).mean()))

        des_contacts = {tmap[n] for n in
                        contacts(list(des["A"].values()), des["B"]) if n in tmap}
        pred_contacts = {seqmap.get(n, n) for n in
                         contacts(list(pred["A"].values()), pred["B"])}
        inter = des_contacts & pred_contacts
        union = des_contacts | pred_contacts
        jac = len(inter) / len(union) if union else 0.0

        rr = dict(r)
        rr["binder_rmsd_vs_design"] = round(brmsd, 2)
        rr["contact_jaccard"] = round(jac, 3)
        rr["n_contacts_pred"] = len(pred_contacts)
        out.append(rr)

    if not out:
        sys.exit("no comparable pairs found")
    keys = list(out[0].keys())
    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
        w.writeheader()
        w.writerows(out)

    import statistics as st
    for c in "ABC":
        v = [o for o in out if o["design"].startswith("class" + c)]
        if not v:
            continue
        rm = [o["binder_rmsd_vs_design"] for o in v]
        jc = [o["contact_jaccard"] for o in v]
        print(f"  class {c}: n={len(v):>3}  binder RMSD median {st.median(rm):6.1f} A"
              f"   jaccard median {st.median(jc):.2f}"
              f"   self-consistent (<5 A) {sum(1 for x in rm if x < 5)}")
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
