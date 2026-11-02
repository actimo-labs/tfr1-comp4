#!/usr/bin/env python3
"""
Stage 80 - rank for the three objectives and build the submission package.

Ranking follows the challenge's own order, which is NOT a weighted sum:

  1. pH-selective binding   (human 7.4 vs human 6.0)
  2. mouse cross-reactivity (mouse 7.4 vs human 7.4), among designs whose
                             human KD is below 1 uM
  3. human affinity at 7.4, with a functional epitope counting in favour

Because objective 2 is explicitly conditioned on objective-3 performance, and
objective 1 is scored by a mechanism we can only estimate, the slate is built
as a PORTFOLIO rather than a single ranked cut - which is also what the
published guidance rewards, since in-silico scores have predicted real binding
at roughly chance in prior Adaptyv rounds.

Classes, in priority order of slot allocation:
  A  cross-species switch + survived refold   objectives 1+2+3
  B  strongest human switches                 objectives 1+3
  C  highest affinity, no deliberate switch   objective 3, and the internal
                                              control for the switch hypothesis
  D  apex-epitope affinity hedge              objective 3

Diversity caps are applied LAST and are allowed to displace higher-scoring
designs: no more than --max-per-backbone designs from one backbone, and a
sequence-identity ceiling. Twenty variations of one scaffold is one bet.

Writes results/submission.csv, results/final_ranked_full.csv,
results/rejected_candidates.csv, results/tfr1_metrics.csv and structures/.
"""
import argparse
import csv
import glob
import json
import os
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
N_SLOTS = 20


def load_scores(p):
    if not os.path.exists(p):
        return {}
    raw = json.load(open(p))
    rows = raw["scores"] if isinstance(raw, dict) else raw
    return {r["name"]: r for r in rows}


def ident(a, b):
    if not a or not b:
        return 0.0
    n = min(len(a), len(b))
    return sum(1 for i in range(n) if a[i] == b[i]) / max(len(a), len(b))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-per-backbone", type=int, default=2)
    ap.add_argument("--max-identity", type=float, default=0.70)
    ap.add_argument("--out", type=Path, default=ROOT / "results")
    ap.add_argument("--structures", type=Path, default=ROOT / "structures")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    args.structures.mkdir(parents=True, exist_ok=True)

    # ---- gather ---------------------------------------------------------
    refold = []
    for f, rnd in (("results/06_refold.tsv", "r1"),
                   ("results/06_refold_r2.tsv", "r2")):
        if os.path.exists(ROOT / f):
            for r in csv.DictReader(open(ROOT / f), delimiter="\t"):
                r["round"] = rnd
                refold.append(r)
    mouse = load_scores(ROOT / "boltz/mouse/scores.json")
    qc = {}
    for f in ("gen/epitope_qc.tsv", "gen/epitope_qc_r2.tsv"):
        if os.path.exists(ROOT / f):
            for r in csv.DictReader(open(ROOT / f), delimiter="\t"):
                qc[r["name"]] = r

    # binder sequences from every yaml we wrote
    seqs = {}
    for tag in ("refold", "refold_r2", "screen", "screen_r2", "mouse"):
        for y in glob.glob(str(ROOT / f"boltz/{tag}/yaml/*.yaml")):
            L = open(y).read().splitlines()
            for i, l in enumerate(L):
                if l.strip() == "id: X":
                    for n in L[i + 1:i + 3]:
                        if n.strip().startswith("sequence:"):
                            seqs[Path(y).stem] = n.split("sequence:", 1)[1].strip()
                    break

    cands = []
    for r in refold:
        if r["survives"] != "True":
            continue
        name = r["variant"]
        parent = r["parent"]
        seq = seqs.get(name)
        if not seq:
            continue
        m = mouse.get(name, {})
        hum = float(r["variant_iptm"])
        mou = m.get("binder_target_iptm")
        q = qc.get(parent, {})
        backbone = parent.split("_", 1)[1].rsplit("_s", 1)[0]
        cands.append(dict(
            name=name, parent=parent, backbone=backbone, round=r["round"],
            sequence=seq, length=len(seq),
            human_iptm=hum,
            mouse_iptm=round(mou, 4) if mou is not None else None,
            mouse_ratio=(round(mou / hum, 3) if mou and hum else None),
            ddg=float(r["total_ddg"] or 0),
            fold_shift=round(2.718281828 ** (float(r["total_ddg"] or 0) / 0.593), 1),
            scope=r["scope"], n_his=r["n_his"], mutations=r["mutations"],
            n_hotspots=int(q.get("n_hotspots", 0) or 0),
            hotspots=q.get("hotspots", ""),
            hotspot_set=q.get("hotspot_set", ""),
            cterm_dist=float(q.get("cterm_dist", 0) or 0),
            binder_ptm=float(r["binder_ptm"]), iplddt=float(r["iplddt"])))

    # ---- unmutated parents, as the no-switch control -------------------
    # §7 asks for designs with NO deliberate switch: pure objective-3
    # insurance AND the internal control for the switch hypothesis. Every
    # refold survivor carries His substitutions by construction, so without
    # this the slate has no control arm and the switch claim has nothing to be
    # compared against.
    parents_scored = {}
    for tag in ("screen", "screen_r2"):
        parents_scored.update(load_scores(ROOT / f"boltz/{tag}/scores.json"))
    chosen_parents = {c["parent"] for c in cands}
    for pname, ps in parents_scored.items():
        q = qc.get(pname)
        if not q or q.get("passes") != "True":
            continue
        seq = seqs.get(pname)
        if not seq or pname in chosen_parents:
            continue
        hum = ps.get("binder_target_iptm") or 0
        if hum < 0.45:
            continue
        cands.append(dict(
            name=pname, parent=pname,
            backbone=pname.split("_", 1)[1].rsplit("_s", 1)[0],
            round="r1" if "_r2" not in pname else "r2",
            sequence=seq, length=len(seq), human_iptm=hum,
            # controls were mouse-scored too, in the same boltz/mouse run
            mouse_iptm=(round(mouse[pname]["binder_target_iptm"], 4)
                        if pname in mouse else None),
            mouse_ratio=(round(mouse[pname]["binder_target_iptm"] / hum, 3)
                         if pname in mouse and hum else None),
            ddg=0.0, fold_shift=1.0,
            scope="none", n_his=0, mutations="",
            n_hotspots=int(q.get("n_hotspots", 0) or 0),
            hotspots=q.get("hotspots", ""), hotspot_set=q.get("hotspot_set", ""),
            cterm_dist=float(q.get("cterm_dist", 0) or 0),
            binder_ptm=float(ps.get("binder_ptm") or 0),
            iplddt=float(ps.get("complex_iplddt") or 0),
            is_control=True))
    # ---- self-consistency gate (HARD) ----------------------------------
    # Two seeds must agree on WHERE the binder sits, not merely on how
    # confident they are. Challenge 03 saw two seeds agree on interface
    # confidence to within 0.03 while placing the binder 37 A apart. Here
    # eight designs fail: three flip to an entirely different site (jaccard
    # 0.0, centroid shifts 41-76 A after superposing the targets), and one of
    # those had already been selected into class A on a mouse ratio of 1.14.
    sc_path = ROOT / "results/selfconsistency.csv"
    if sc_path.exists():
        sc = {r["name"]: r for r in csv.DictReader(open(sc_path))}
        before = len(cands)
        cands = [c for c in cands
                 if sc.get(c["name"], {}).get("passes", "True") == "True"]
        print(f"[80] self-consistency gate: {before} -> {len(cands)} "
              f"({before - len(cands)} dropped for seed-dependent poses)")
    else:
        print("[80] WARNING: no selfconsistency.csv - pose gate NOT applied")

    print(f"[80] {len(cands)} candidates "
          f"({sum(1 for c in cands if c.get('is_control'))} unmutated controls)")
    print(f"[80] mouse-scored: {sum(1 for c in cands if c['mouse_iptm'] is not None)}")

    # ---- classes --------------------------------------------------------
    # IMPORTANT, learned from the mouse screen: a conserved-cation switch does
    # NOT imply a cross-reactive binder. Those are independent properties. The
    # switch is one residue; the interface is 15-20, and most of them differ
    # between species. Measured: of the four surviving cross-species-switch
    # designs, mouse:human ipTM ratios are 1.34, 0.46, 0.20 and 0.19 - and the
    # best human binder of the whole campaign (0.755) sits at 0.19. Meanwhile
    # several SPECIES-SPECIFIC-switch designs cross-react well (1.34, 1.14,
    # 1.08).
    #
    # Objective 2 is scored by what the assay measures - mouse binding - not by
    # our mechanism label. So class A is defined by MEASURED mouse retention,
    # with the switch as a secondary ordering term, not the other way round.
    for c in cands:
        ratio = c["mouse_ratio"] or 0
        if ratio >= 0.75 and c["mouse_iptm"] and c["mouse_iptm"] >= 0.30:
            c["cls"] = "A"
        elif c["ddg"] >= 0.7:
            c["cls"] = "B"
        elif c.get("is_control"):
            c["cls"] = "D"
        else:
            c["cls"] = "C"

    # within-class ordering
    def key_A(c):   # objective 2: measured mouse retention first, then switch
        return (-(c["mouse_ratio"] or 0), -(c["ddg"]), -c["human_iptm"])

    def key_B(c):   # human switch: switch strength, then affinity
        return (-(c["ddg"]), -c["human_iptm"])

    def key_aff(c):  # affinity: ipTM, then epitope coverage
        return (-c["human_iptm"], -c["n_hotspots"])

    order = {"A": key_A, "B": key_B, "C": key_aff, "D": key_aff}
    quota = {"A": 8, "B": 5, "C": 4, "D": 3}

    chosen, used_bb, rejected = [], {}, []

    def try_add(c, why):
        if len(chosen) >= N_SLOTS:
            return False
        if used_bb.get(c["backbone"], 0) >= args.max_per_backbone:
            rejected.append((c, f"backbone cap ({c['backbone']})"))
            return False
        for s in chosen:
            if ident(c["sequence"], s["sequence"]) > args.max_identity:
                rejected.append((c, f">{args.max_identity:.0%} identity to {s['name']}"))
                return False
        c["selected_as"] = why
        chosen.append(c)
        used_bb[c["backbone"]] = used_bb.get(c["backbone"], 0) + 1
        return True

    # Pools, not partitions. A and B are defined by their mechanism, but C is
    # the AFFINITY slot and must draw from everything still unselected -
    # otherwise a design with both top affinity and a decent switch gets
    # labelled B, loses B's slots to stronger switches, and never reaches C.
    # That happened on the first pre-flight: the two best human binders in the
    # campaign (ipTM 0.772 and 0.755) were excluded from the slate entirely,
    # while objective 3 is a ranked criterion.
    for cls in ("A", "B", "C", "D"):
        if cls == "C":
            pool = sorted([c for c in cands if c not in chosen], key=key_aff)
        elif cls == "D":
            pool = sorted([c for c in cands
                           if c not in chosen and c.get("is_control")],
                          key=key_aff)
        else:
            pool = sorted([c for c in cands if c["cls"] == cls], key=order[cls])
        n = 0
        for c in pool:
            if n >= quota[cls]:
                break
            if try_add(c, cls):
                n += 1
        print(f"[80] class {cls}: {n}/{quota[cls]} filled from {len(pool)} candidates")

    # backfill any unfilled slots by affinity, respecting the same caps
    if len(chosen) < N_SLOTS:
        rest = sorted([c for c in cands if c not in chosen], key=key_aff)
        for c in rest:
            if len(chosen) >= N_SLOTS:
                break
            try_add(c, f"{c['cls']}-backfill")
    print(f"[80] slate: {len(chosen)}/{N_SLOTS}")

    # ---- outputs --------------------------------------------------------
    with open(args.out / "submission.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["name", "sequence", "molecule_class"])
        for c in chosen:
            w.writerow([c["name"], c["sequence"], "single_chain"])

    cols = ["name", "selected_as", "cls", "round", "length", "human_iptm",
            "mouse_iptm", "mouse_ratio", "ddg", "fold_shift", "scope", "n_his",
            "mutations", "n_hotspots", "hotspots", "hotspot_set", "cterm_dist",
            "binder_ptm", "iplddt", "backbone", "parent", "sequence"]
    with open(args.out / "tfr1_metrics.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(chosen)
    with open(args.out / "final_ranked_full.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for c in sorted(cands, key=key_aff):
            c.setdefault("selected_as", "")
            w.writerow(c)
    with open(args.out / "rejected_candidates.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["name", "human_iptm", "ddg", "scope", "reason"])
        for c, why in rejected:
            w.writerow([c["name"], c["human_iptm"], c["ddg"], c["scope"], why])

    # structures
    n = 0
    for c in chosen:
        hits = []
        for tag in ("refold", "refold_r2", "screen", "screen_r2"):
            hits += glob.glob(
                str(ROOT / f"boltz/{tag}/pred/**/{c['name']}_model_0.pdb"),
                recursive=True)
        if hits:
            shutil.copy(hits[0], args.structures / f"{c['name']}_complex.pdb")
            n += 1
    print(f"[80] structures copied: {n}/{len(chosen)}")

    from collections import Counter
    print(f"\n[80] slate composition: {dict(Counter(c['cls'] for c in chosen))}")
    print(f"[80] distinct backbones: {len(set(c['backbone'] for c in chosen))}")
    print(f"[80] length range: {min(c['length'] for c in chosen)}-"
          f"{max(c['length'] for c in chosen)}")
    print(f"\n  {'name':<42}{'cls':>4}{'hIPTM':>7}{'mIPTM':>7}{'ratio':>7}{'ddG':>6}{'fold':>6}")
    for c in chosen:
        print(f"  {c['name'][:41]:<42}{c['cls']:>4}{c['human_iptm']:>7.3f}"
              f"{(c['mouse_iptm'] or 0):>7.3f}{(c['mouse_ratio'] or 0):>7.2f}"
              f"{c['ddg']:>6.2f}{c['fold_shift']:>6.0f}")
    print(f"\n[out] {args.out}/submission.csv")


if __name__ == "__main__":
    main()
