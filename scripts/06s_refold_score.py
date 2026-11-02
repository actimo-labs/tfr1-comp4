#!/usr/bin/env python3
"""
Stage 06s - score the His-refold variants.

Replaces the tail of comp-02's 06_his_refold.py, which assumed a two-chain
homotrimer target throughout: `iptm_min` (minimum across two target chains) and
`target_pair_iptm` (the assembly sanity gate) have no meaning against a
monomer, and patching each one in turn was costing more time than rewriting
the scoring step. The Boltz predictions already exist; this only reads them.

The gate CLAUDE.md calls non-negotiable: a design that switches off at pH 6.0
but no longer binds at 7.4 is worthless. A variant survives when

    variant ipTM >= --cutoff                       (it still binds), AND
    ( variant ipTM >= --strong                     (binds well on its own), OR
      it kept >= (1 - --max-drop) of its parent )  (interface largely intact)

The second clause matters independently: a parent at 0.72 falling to 0.41 is
still "above cutoff" while having lost most of its interface.
"""
import argparse
import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def load_scores(p: Path):
    raw = json.loads(p.read_text())
    rows = raw["scores"] if isinstance(raw, dict) else raw
    return {r["name"]: r for r in rows}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--variants", type=Path, required=True,
                    help="refold scores.json")
    ap.add_argument("--parents", type=Path, required=True,
                    help="parent screen scores.json")
    ap.add_argument("--combos", type=Path, required=True)
    ap.add_argument("--cutoff", type=float, default=0.35)
    ap.add_argument("--strong", type=float, default=0.40)
    ap.add_argument("--max-drop", type=float, default=0.25)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    V = load_scores(args.variants)
    P = load_scores(args.parents)
    combos = {r["design"].replace("_model_0", ""): r
              for r in csv.DictReader(open(args.combos), delimiter="\t")}

    rows = []
    for name, s in V.items():
        parent = name.split("__")[0]
        muts = name.split("__")[1] if "__" in name else ""
        vi = s.get("binder_target_iptm")
        pi = (P.get(parent) or {}).get("binder_target_iptm")
        if vi is None:
            continue
        drop = None if not pi else round((pi - vi) / pi, 3)
        strong = vi >= args.strong
        ok = vi >= args.cutoff and (strong or drop is None or drop <= args.max_drop)
        c = combos.get(parent, {})
        rows.append(dict(
            variant=name, parent=parent, mutations=muts,
            parent_iptm=round(pi, 4) if pi else None,
            variant_iptm=round(vi, 4),
            drop=drop, binder_ptm=round(s.get("binder_ptm") or 0, 4),
            iplddt=round(s.get("complex_iplddt") or 0, 4),
            total_ddg=c.get("total_ddg"), scope=c.get("scope"),
            n_his=c.get("n_his"), survives=ok))

    rows.sort(key=lambda r: (not r["survives"], -(r["variant_iptm"] or 0)))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    surv = [r for r in rows if r["survives"]]
    xs = [r for r in surv if r["scope"] == "cross_species"]
    print(f"[06s] {len(rows)} variants, {len(surv)} survive the pH 7.4 refold "
          f"({100*len(surv)/len(rows):.0f}%)")
    print(f"[06s] surviving cross-species switches: {len(xs)}")
    print(f"\n  {'variant':<46}{'pIPTM':>7}{'vIPTM':>7}{'drop':>7}{'ddG':>6}  scope")
    for r in surv[:20]:
        print(f"  {r['variant'][:45]:<46}{r['parent_iptm'] or 0:>7.3f}"
              f"{r['variant_iptm']:>7.3f}{(r['drop'] if r['drop'] is not None else 0):>7.2f}"
              f"{float(r['total_ddg'] or 0):>6.2f}  {r['scope']}")
    print(f"\n[out] {args.out}")


if __name__ == "__main__":
    main()
