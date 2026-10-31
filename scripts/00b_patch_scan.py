#!/usr/bin/env python3
"""
Stage 00b - systematic surface-patch scan over 11QC chain A.

Consumes target/conservation_table.tsv (written by 00_conservation.py) and
ranks every candidate binder footprint on the TfR1 monomer surface, so the
epitope choice is made from a table rather than from the literature.

For each surface residue, takes the 13 A shell around it as a candidate
footprint and reports: size, human/mouse conservation fraction, buried area,
conserved cations available as strategy-1 switch anchors, carboxylate count,
distance to the apical apex / transferrin site / nearest glycan sequon, and
how many of the nine VP1u apical-footprint residues it covers.

Also reports, for each candidate switch cation, the nearest target carboxylate
-- a target Asp/Glu within ~6 A inverts the switch (protonated His is then
STABILISED at low pH), which is a hard disqualifier.

Run:  venv-phscan/bin/python -I scripts/00b_patch_scan.py
"""
import os
import csv
import json

import numpy as np
from Bio.PDB import MMCIFParser
from Bio.PDB.Polypeptide import is_aa

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TGT = os.path.join(ROOT, "target")

SHELL = 13.0            # A, candidate binder footprint radius
MIN_PATCH = 10          # discard shells with fewer surface residues
MIN_AREA = 600          # A^2, brief section 4 step 6
TF_CLEAR = 12.0         # A, must be this far from the transferrin site
GLY_CLEAR = 8.0         # A, must be this far from a glycan sequon
ANION_VETO = 6.0        # A, target carboxylate this close inverts the switch

# the nine residues VP1u actually contacts on TfR1 in 11RN (computed in 00)
VP1U_FOOTPRINT = {202, 204, 208, 209, 210, 211, 212, 348, 371}
GLYCAN_SEQUONS = [251, 317, 727]
SWITCH_CANDIDATES = [344, 374, 231, 224, 193, 382, 371, 208]


def load():
    rows = {int(r["pos"]): r
            for r in csv.DictReader(open(os.path.join(TGT,
                                    "conservation_table.tsv")), delimiter="\t")}
    st = MMCIFParser(QUIET=True).get_structure(
        "x", os.path.join(TGT, "structures", "11QC.cif"))
    chain = next(st.get_models())["A"]
    res = {r.id[1]: r for r in chain
           if r.id[0] == " " and is_aa(r, standard=False)}
    xyz = {p: np.array([a.coord for a in r]) for p, r in res.items()}
    return rows, res, xyz


def tip(res, p):
    """Distal side-chain heavy atom - where the charge actually sits."""
    r = res[p]
    ca = r["CA"].coord
    sc = [a for a in r if a.get_id() not in ("N", "CA", "C", "O")]
    return (max(sc, key=lambda a: float(np.linalg.norm(a.coord - ca))).coord
            if sc else ca)


def main():
    rows, res, xyz = load()
    surf = sorted(p for p in rows
                  if rows[p]["surface"] == "True" and p in xyz)
    print(f"[scan] {len(surf)} surface residues on 11QC chain A")

    # pairwise min heavy-atom distance matrix
    n = len(surf)
    idx = {p: i for i, p in enumerate(surf)}
    D = np.zeros((n, n))
    for i, p in enumerate(surf):
        for j in range(i + 1, n):
            q = surf[j]
            d = float(np.linalg.norm(xyz[p][:, None, :] - xyz[q][None, :, :],
                                     axis=-1).min())
            D[i, j] = D[j, i] = d

    patches = []
    for p in surf:
        near = [surf[j] for j in range(n) if D[idx[p], j] <= SHELL]
        if len(near) < MIN_PATCH:
            continue
        ncon = sum(1 for q in near if rows[q]["conserved"] == "True")
        area = sum(float(rows[q]["sasa"]) for q in near)
        ccat = [q for q in near if rows[q]["charge"] == "cationic"
                and rows[q]["conserved"] == "True"]
        anion = [q for q in near if rows[q]["charge"] == "anionic"]
        dtf = min(float(rows[q]["d_tf"]) for q in near
                  if float(rows[q]["d_tf"]) >= 0)
        dgly = min(float(rows[q]["d_glycan"]) for q in near
                   if float(rows[q]["d_glycan"]) >= 0)
        dapex = min(float(rows[q]["d_apex"]) for q in near)
        patches.append(dict(
            centre=p, aa=rows[p]["aa"], domain=rows[p]["domain"],
            n=len(near), n_conserved=ncon, frac_conserved=round(ncon / len(near), 2),
            area=round(area), conserved_cations=ccat, n_cons_cation=len(ccat),
            n_anion=len(anion), d_tf=dtf, d_glycan=dgly, d_apex=dapex,
            n_vp1u=len(set(near) & VP1U_FOOTPRINT), members=sorted(near)))

    elig = [o for o in patches
            if o["d_tf"] > TF_CLEAR and o["d_glycan"] > GLY_CLEAR
            and o["n_cons_cation"] >= 1 and o["area"] >= MIN_AREA]
    print(f"[scan] {len(patches)} patches, {len(elig)} pass "
          f"Tf>{TF_CLEAR} / glycan>{GLY_CLEAR} / >=1 conserved cation / "
          f"area>={MIN_AREA}\n")

    def show(title, lst, k=12):
        print("=" * 104)
        print(title)
        print("=" * 104)
        print(f"{'ctr':>6} {'domain':>14} {'n':>3} {'fCons':>6} {'areaA2':>7} "
              f"{'nCat+':>6} {'nAnion':>7} {'dApex':>6} {'dTf':>6} {'dGly':>6} "
              f"{'nVP1u':>6}  conserved cations")
        for o in lst[:k]:
            print(f"{o['aa']}{o['centre']:<5} {o['domain']:>14} {o['n']:>3} "
                  f"{o['frac_conserved']:>6} {o['area']:>7} "
                  f"{o['n_cons_cation']:>6} {o['n_anion']:>7} {o['d_apex']:>6} "
                  f"{o['d_tf']:>6} {o['d_glycan']:>6} {o['n_vp1u']:>6}  "
                  f"{o['conserved_cations']}")
        print()

    by_cons = lambda o: (-o["frac_conserved"], -o["n_cons_cation"])
    show("A. Most conserved eligible patches, any location", sorted(elig, key=by_cons))
    show("B. Eligible patches overlapping the VP1u apical footprint "
         "(functional epitope)",
         sorted([o for o in elig if o["n_vp1u"] >= 1], key=by_cons))
    show("C. Eligible patches near the apical apex (<=15 A)",
         sorted([o for o in elig if o["d_apex"] <= 15], key=by_cons))

    # ---- switch-cation viability -------------------------------------------
    print("=" * 104)
    print("SWITCH-CATION VIABILITY - nearest target carboxylate to each "
          "candidate anchor")
    print(f"(a target Asp/Glu within {ANION_VETO} A INVERTS the switch: "
          "protonated His would be stabilised at pH 6)")
    print("=" * 104)
    anions = [p for p in rows
              if rows[p]["charge"] == "anionic" and p in xyz]
    verdicts = {}
    for c in SWITCH_CANDIDATES:
        if c not in xyz:
            continue
        near = sorted((round(float(np.linalg.norm(tip(res, c) - tip(res, a))), 1),
                       f"{rows[a]['aa']}{a}") for a in anions
                      if float(np.linalg.norm(tip(res, c) - tip(res, a))) <= 10)
        conserved = rows[c]["conserved"] == "True"
        vetoed = bool(near) and near[0][0] <= ANION_VETO
        ok = conserved and not vetoed
        verdicts[c] = dict(conserved=conserved, nearest_anion=near[0] if near else None,
                           vetoed=vetoed, usable=ok)
        tag = ("USABLE  " if ok else
               "veto:anion" if vetoed else "veto:not-conserved")
        print(f"  {rows[c]['aa']}{c:<5} -> mouse {rows[c]['mouse_aa']}"
              f"{rows[c]['mouse_pos']:<5} {'CONS' if conserved else 'diff'}  "
              f"{tag:18s} "
              + (", ".join(f"{a}@{d}A" for d, a in near[:3]) if near
                 else "no carboxylate within 10 A"))

    json.dump(dict(patches=patches, eligible=[o["centre"] for o in elig],
                   switch_cations=verdicts),
              open(os.path.join(TGT, "patch_scan.json"), "w"), indent=1)
    print(f"\n[out] target/patch_scan.json")


if __name__ == "__main__":
    main()
