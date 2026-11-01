#!/usr/bin/env python3
"""
Stage 05b - apply CONSERVATION POLICY to the raw pH-scan sites.

05_ph_scan.py is the physics engine and is deliberately target-agnostic: it
finds any interface position whose His substitution would be destabilised by a
nearby target cation. It knows nothing about which cations survive into mouse,
or which are chemically poisoned. Run bare against TfR1 it produced 497 sites,
and two thirds of them are unusable:

  target cation   construct   conserved?        sites   verdict
  A87ARG          R208        NO (-> N211)        236   species-specific
  A223LYS         K344        yes                  96   VETO: E214 at 2.8 A
  A250LYS         K371        NO (-> Q373)         78   species-specific
  A253LYS         K374        yes                   77   THE anchor

K344 is the trap stage 00 identified: it is conserved and exposed, but
salt-bridged to E214 at 2.8 A, so a protonated His beside it is STABILISED at
pH 6.0 and the design binds *tighter* when it should let go. Sites whose only
nearby cation is vetoed are therefore discarded, not merely downranked.

R208 and K371 are real switches in human and do nothing in mouse. They are not
discarded -- they are exactly what the class-C apex hedge is for -- but they
are labelled `species_specific` so objective 2 is never claimed for them.

Tiering, instead of challenge 02's single 2.8 kcal/mol pass/fail. Measured on
this target the maximum single-site ddG is 1.29 kcal/mol and no design has
three independent strong sites, so 2.8 (~100x) is unreachable here and a bare
pass/fail reports zero and hides the real landscape:

  strong    >= 2.0 kcal/mol   ~ 30x Kd shift
  moderate  >= 1.4            ~ 10x
  weak      >= 1.0            ~  5x

That is in the range of the Baker-lab strategy-1 results at the weaker end
(they report 122x on TNFR2, up to ~1000x on EphA2), and well short of the best
case. The methods report must say so rather than implying a 100x switch.

Additivity still only holds between electrostatically independent sites, so
the >8 A separation rule from stage 05 is re-applied here using real
coordinates from the predicted complex.

  venv-phscan/bin/python -I scripts/05b_switch_select.py
"""
import argparse
import csv
import glob
import itertools
import json
import os
from collections import defaultdict
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser

ROOT = Path(__file__).resolve().parent.parent
OFFSET = 121            # boltz target resnum + OFFSET = construct numbering
RT = 0.593              # kcal/mol at 298 K
MIN_SEP = 8.0           # A between His, for additivity to be defensible
MAX_SUBS = 3

# Policy, straight out of target/conservation_report.md.
VETOED = {344: "salt-bridged to E214 at 2.8 A - would invert the switch"}
CONSERVED_OK = {374, 224, 231, 193, 382, 161}   # conserved AND not vetoed
TIERS = [("strong", 2.0), ("moderate", 1.4), ("weak", 1.0)]


def fold_shift(ddg):
    return float(np.exp(ddg / RT))


def cation_positions(field):
    """'A253LYS@4.67;A87ARG@5.1' -> [(374,'LYS',4.67), ...] in construct numbering."""
    out = []
    for tok in (field or "").split(";"):
        tok = tok.strip()
        if not tok or "@" not in tok:
            continue
        label, dist = tok.split("@")
        digits = "".join(ch for ch in label if ch.isdigit())
        if not digits:
            continue
        out.append((int(digits) + OFFSET, label[-3:], float(dist)))
    return out


def classify(cats):
    """Return (usable, tier, detail) for a site given its nearby cations."""
    if not cats:
        return False, None, "no cation"
    live = [c for c in cats if c[0] not in VETOED]
    if not live:
        why = VETOED[cats[0][0]]
        return False, None, f"only vetoed cation {cats[0][0]} ({why})"
    conserved = [c for c in live if c[0] in CONSERVED_OK]
    if conserved:
        return True, "cross_species", \
            "conserved cation " + ",".join(str(c[0]) for c in conserved)
    return True, "species_specific", \
        "human-only cation " + ",".join(str(c[0]) for c in live)


def cb_coords(pdb_path, chain, positions):
    try:
        m = next(PDBParser(QUIET=True).get_structure("x", pdb_path).get_models())
        ch = m[chain]
    except Exception:
        return {}
    out = {}
    for r in ch:
        if r.id[0] != " " or r.id[1] not in positions:
            continue
        at = r["CB"] if "CB" in r else (r["CA"] if "CA" in r else None)
        if at is not None:
            out[r.id[1]] = at.coord
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scan", type=Path, default=ROOT / "results/ph_scan.tsv")
    ap.add_argument("--pred", type=Path, default=ROOT / "boltz/screen/pred")
    ap.add_argument("--binder-chain", default="X")
    ap.add_argument("--out", type=Path, default=ROOT / "results/switch_candidates.tsv")
    ap.add_argument("--max-subs", type=int, default=MAX_SUBS)
    args = ap.parse_args()

    rows = list(csv.DictReader(open(args.scan), delimiter="\t"))
    by_design = defaultdict(list)
    for r in rows:
        by_design[r["design"]].append(r)
    print(f"[05b] {len(rows)} raw sites across {len(by_design)} designs")

    # --- site-level policy ------------------------------------------------
    kept, dropped = defaultdict(list), defaultdict(int)
    for design, rs in by_design.items():
        for r in rs:
            cats = cation_positions(r.get("near_cationic", ""))
            ok, tier, detail = classify(cats)
            ddg = float(r["ddg_kcal"])
            if not ok or ddg <= 0:
                dropped[detail.split(" (")[0] if not ok else "ddg<=0"] += 1
                continue
            kept[design].append(dict(pos=int(r["position"]), ddg=ddg,
                                     wt=r["wt"], tier=tier, detail=detail,
                                     cats=cats))
    nk = sum(len(v) for v in kept.values())
    print(f"[05b] {nk} sites survive policy; dropped:")
    for k, v in sorted(dropped.items(), key=lambda x: -x[1])[:6]:
        print(f"        {k}: {v}")

    # --- combinations -----------------------------------------------------
    out = []
    for design, sites in kept.items():
        hits = glob.glob(str(args.pred / "**" / f"{design}.pdb"), recursive=True)
        coords = cb_coords(hits[0], args.binder_chain,
                           {s["pos"] for s in sites}) if hits else {}
        sites.sort(key=lambda s: -s["ddg"])
        sites = sites[:8]
        best = None
        for n in range(1, args.max_subs + 1):
            for combo in itertools.combinations(sites, n):
                if n > 1 and coords:
                    ok = True
                    for a, b in itertools.combinations(combo, 2):
                        ca, cb = coords.get(a["pos"]), coords.get(b["pos"])
                        if ca is None or cb is None:
                            continue
                        if float(np.linalg.norm(ca - cb)) < MIN_SEP:
                            ok = False
                            break
                    if not ok:
                        continue
                elif n > 1 and not coords:
                    continue        # cannot verify independence -> do not claim it
                tot = sum(c["ddg"] for c in combo)
                xs = "cross_species" if all(c["tier"] == "cross_species"
                                            for c in combo) else "species_specific"
                if best is None or tot > best["total_ddg"]:
                    best = dict(design=design, n_his=n, total_ddg=round(tot, 3),
                                fold_shift=round(fold_shift(tot), 1),
                                mutations=";".join(f"{c['wt']}{c['pos']}H"
                                                   for c in combo),
                                scope=xs,
                                cations=";".join(sorted({str(x[0]) for c in combo
                                                         for x in c["cats"]})))
        if best:
            tier = next((t for t, cut in TIERS if best["total_ddg"] >= cut), "sub")
            best["tier"] = tier
            out.append(best)

    out.sort(key=lambda r: -r["total_ddg"])
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(out)

    print(f"\n[05b] {len(out)} designs with at least one usable switch")
    from collections import Counter
    print("  by tier :", dict(Counter(r["tier"] for r in out)))
    print("  by scope:", dict(Counter(r["scope"] for r in out)))
    print(f"\n  {'design':<40}{'ddG':>6}{'fold':>7}{'n':>3}  {'scope':<16} mutations")
    for r in out[:20]:
        print(f"  {r['design'][:39]:<40}{r['total_ddg']:>6.2f}{r['fold_shift']:>7.0f}"
              f"{r['n_his']:>3}  {r['scope']:<16} {r['mutations']}")
    json.dump(out, open(args.out.with_suffix(".json"), "w"), indent=1)
    print(f"\n[out] {args.out}")


if __name__ == "__main__":
    main()
