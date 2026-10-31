#!/usr/bin/env python3
"""
De novo compliance check, plus the intra-submission diversity rule.

Two separate requirements from the brief:

  [OFFICIAL] section 1  -- de novo designs only.
  section 5             -- do not emit His-Ala-Glu-Gly or any contiguous run
                           matching GLP-1, exendin-4, glucagon or GIP.
  section 4             -- no two submitted designs above 70% identity.

This matters most for class B. Those designs are built by scaffolding GLP-1's
message-module BACKBONE (a geometric constraint), with the sequence designed
from scratch by ProteinMPNN. That is defensible as de novo precisely because no
GLP-1 sequence survives -- so it has to actually be checked, not asserted.

Reports, per design:
  * longest common substring against each natural incretin peptide
  * whether that exceeds --max-run
  * max identity against every other design in the set

Exit code 1 if anything fails, so it can gate the submission build.
"""

import argparse
import csv
import sys
from itertools import combinations

# Natural incretin/glucagon-family peptides the design must not reproduce.
NATURAL = {
    "GLP-1(7-36)":  "HAEGTFTSDVSSYLEGQAAKEFIAWLVKGR",
    "GLP-1(1-37)":  "HDEFERHAEGTFTSDVSSYLEGQAAKEFIAWLVKGRG",
    "exendin-4":    "HGEGTFTSDLSKQMEEEAVRLFIEWLKNGGPSSGAPPPS",
    "glucagon":     "HSQGTFTSDYSKYLDSRRAQDFVQWLMNT",
    "GIP(1-42)":    "YAEGTFISDYSIAMDKIHQQDFVNWLLAQKGKKNDWKHNITQ",
    "oxyntomodulin":"HSQGTFTSDYSKYLDSRRAQDFVQWLMNTKRNRNNIA",
}


def longest_common_substring(a, b):
    """Return (length, substring). O(len(a)*len(b)) DP, fine at peptide scale."""
    best, best_end = 0, 0
    prev = [0] * (len(b) + 1)
    for i in range(1, len(a) + 1):
        cur = [0] * (len(b) + 1)
        for j in range(1, len(b) + 1):
            if a[i - 1] == b[j - 1]:
                cur[j] = prev[j - 1] + 1
                if cur[j] > best:
                    best, best_end = cur[j], i
        prev = cur
    return best, a[best_end - best:best_end]


def pct_identity(a, b, min_overlap_frac=0.8):
    """
    Ungapped identity, best offset, normalised by the SHORTER sequence.

    Normalising by the overlap instead (the obvious implementation) is wrong and
    silently catastrophic: a one-residue overhang that happens to match scores
    100%, so every pair of designs looks like a duplicate. Normalising by the
    shorter length means a partial overlap can never exceed the fraction of that
    sequence actually aligned. Offsets that align less than `min_overlap_frac`
    of the shorter sequence are skipped outright.
    """
    if not a or not b:
        return 0.0
    shorter = min(len(a), len(b))
    need = max(1, int(shorter * min_overlap_frac))
    best = 0
    for off in range(-(len(b) - 1), len(a)):
        m = n = 0
        for i in range(len(a)):
            j = i - off
            if 0 <= j < len(b):
                n += 1
                m += a[i] == b[j]
        if n >= need:
            best = max(best, m / shorter)
    return best * 100


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--manifest", required=True,
                   help="CSV with columns: design, sequence")
    p.add_argument("--max-run", type=int, default=5,
                   help="longest contiguous match to a natural peptide allowed")
    p.add_argument("--max-identity", type=float, default=70.0)
    p.add_argument("--out", default="")
    args = p.parse_args()

    rows = list(csv.DictReader(open(args.manifest)))
    if not rows:
        sys.exit(f"{args.manifest} is empty")
    designs = [(r["design"], r["sequence"].strip().upper()) for r in rows
               if r.get("sequence", "").strip()]
    print(f"checking {len(designs)} designs\n")

    out, n_fail = [], 0
    print(f"{'design':<22}{'worst natural match':<34}{'len':>4}  verdict")
    for name, seq in designs:
        worst_len, worst_sub, worst_src = 0, "", ""
        for src, nat in NATURAL.items():
            L, sub = longest_common_substring(seq, nat)
            if L > worst_len:
                worst_len, worst_sub, worst_src = L, sub, src
        bad = worst_len > args.max_run
        n_fail += bad
        print(f"{name:<22}{worst_src + ' ' + repr(worst_sub):<34}{worst_len:>4}  "
              f"{'FAIL' if bad else 'ok'}")
        out.append({"design": name, "sequence": seq,
                    "worst_natural_source": worst_src,
                    "worst_natural_match": worst_sub,
                    "worst_natural_len": worst_len,
                    "denovo_ok": not bad})

    # explicit check on the exact motif the brief calls out
    print()
    for name, seq in designs:
        if "HAEG" in seq:
            print(f"  !! {name} contains HAEG (His-Ala-Glu-Gly)")
            n_fail += 1

    # pairwise identity
    print(f"\npairwise identity (limit {args.max_identity:.0f}%):")
    worst = []
    for (n1, s1), (n2, s2) in combinations(designs, 2):
        idp = pct_identity(s1, s2)
        if idp > args.max_identity:
            worst.append((idp, n1, n2))
    if worst:
        for idp, n1, n2 in sorted(worst, reverse=True)[:20]:
            print(f"  FAIL {idp:5.1f}%  {n1}  vs  {n2}")
        n_fail += len(worst)
    else:
        mx = max((pct_identity(s1, s2) for (_, s1), (_, s2)
                  in combinations(designs, 2)), default=0.0)
        print(f"  all pairs under limit (max observed {mx:.1f}%)")

    if args.out:
        with open(args.out, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(out[0].keys()))
            w.writeheader()
            w.writerows(out)
        print(f"\nwrote {args.out}")

    print(f"\n{'FAILED' if n_fail else 'PASSED'}: {n_fail} violation(s)")
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
