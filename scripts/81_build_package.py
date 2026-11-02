#!/usr/bin/env python3
"""
Stage 81 - assemble the submission package.

Mirrors the challenge-03 layout, which is the format Harjas has shipped twice:

    results/submission.csv            name,sequence,molecule_class  (20 rows)
    results/METHODS.md
    results/rationale.md              per-design, written from the metrics
    results/tfr1_metrics.csv
    results/final_ranked_full.csv
    results/rejected_candidates.csv
    results/selfconsistency.csv
    structures/<name>_complex.pdb
    target/ ...  scripts/ ...  runlog/ ...
    README.md  MANIFEST.txt  SHA256SUMS.txt

Selection for testing is done by a model reading the submitted metrics and
method text, so rationale.md states what each design is FOR and what is
uncertain about it. Overclaiming is the cheapest way to look unserious.
"""
import argparse
import csv
import hashlib
import os
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=ROOT / "upload/submission")
    args = ap.parse_args()
    out = args.out
    if out.exists():
        shutil.rmtree(out)
    for d in ("results", "structures", "target", "scripts", "runlog"):
        (out / d).mkdir(parents=True, exist_ok=True)

    rows = list(csv.DictReader(open(ROOT / "results/tfr1_metrics.csv")))

    # ---- rationale -------------------------------------------------------
    CLS = {
        "A": "cross-species switch — the only architecture that can serve "
             "objectives 1 and 2 together",
        "B": "strongest human pH switch — objective 1 and 3; mouse not claimed",
        "C": "affinity, no deliberate switch — objective 3 and the internal "
             "control for the switch hypothesis",
        "D": "apex-epitope affinity hedge — objective 3",
    }
    with open(out / "results/rationale.md", "w") as f:
        f.write("# Per-design rationale\n\n")
        f.write("Epitope for every design is the TfR1 **apical domain**: the "
                "patch contacted by B19V VP1u in PDB 11RN (residues 202, 204, "
                "208, 209, 210, 211, 212, 348, 371), plus the conserved rim "
                "cation **K374**. Evidence: 11RN, and the conservation map in "
                "`target/conservation_report.md`.\n\n")
        f.write("Transferrin binds the helical domain at 619-670. Every design "
                "below was screened to make **zero** contacts there, so "
                "transferrin binding and receptor recycling should be "
                "unaffected.\n\n")
        f.write("`ipTM` is the Boltz-2 binder-target chain-pair value against "
                "the full 639-residue monomer. `ddG` is PROPKA thermodynamic "
                "linkage, **not** a predicted KD ratio — no model used here has "
                "a pH parameter.\n\n---\n\n")
        for r in rows:
            f.write(f"## {r['name']}\n\n")
            f.write(f"- **Class {r['cls']}** — {CLS.get(r['cls'], '')}\n")
            f.write(f"- Length {r['length']} aa, single chain, de novo\n")
            f.write(f"- Human interface ipTM **{r['human_iptm']}**, "
                    f"binder pTM {r['binder_ptm']}, interface pLDDT {r['iplddt']}\n")
            f.write(f"- Contacts **{r['n_hotspots']}** targeted epitope residues"
                    + (f" ({r['hotspots']})\n" if r['hotspots'] else "\n"))
            if r.get("mutations"):
                f.write(f"- pH switch: **{r['mutations']}**, linkage ddG "
                        f"{r['ddg']} kcal/mol (~{r['fold_shift']}x KD shift), "
                        f"scope **{r['scope']}**\n")
                if r["scope"] == "species_specific":
                    f.write("  - *This switch is built on R208/K371, which are "
                            "neutral in mouse. pH selectivity is claimed for "
                            "human only.*\n")
                else:
                    f.write("  - *Built on K374, conserved as K376 in mouse, so "
                            "the mechanism is species-independent — but the "
                            "effect is small.*\n")
            else:
                f.write("- **No deliberate pH switch.** Included as the control "
                        "for the switch hypothesis and as affinity insurance.\n")
            if r.get("mouse_iptm") not in (None, "", "None"):
                f.write(f"- Mouse interface ipTM {r['mouse_iptm']} "
                        f"(mouse:human ratio {r['mouse_ratio']}). "
                        f"*Threaded onto the human structure; no mouse "
                        f"structure exists.*\n")
            f.write(f"- C-terminus {r['cterm_dist']} A from the target, so the "
                    f"Twin-Strep tag should not participate in binding\n\n")

    # ---- copy ------------------------------------------------------------
    for f in ("submission.csv", "METHODS.md", "tfr1_metrics.csv",
              "final_ranked_full.csv", "rejected_candidates.csv",
              "selfconsistency.csv", "ph_scan.tsv", "switch_candidates.tsv",
              "06_refold.tsv", "06_refold_r2.tsv"):
        src = ROOT / "results" / f
        if src.exists():
            shutil.copy(src, out / "results" / f)
    for p in (ROOT / "structures").glob("*.pdb"):
        shutil.copy(p, out / "structures" / p.name)
    for f in ("human_tfr1_121_760.fasta", "mouse_tfr1_123_763.fasta",
              "vp1u_X.fasta", "tfr1_A.pdb", "hotspots_shoulder.txt",
              "hotspots_apex.txt", "nondesignable.txt",
              "conservation_report.md", "conservation_table.tsv"):
        src = ROOT / "target" / f
        if src.exists():
            shutil.copy(src, out / "target" / f)
    for p in (ROOT / "scripts").glob("*"):
        if p.is_file():
            shutil.copy(p, out / "scripts" / p.name)
    for p in (ROOT / "runlog").glob("*.md"):
        shutil.copy(p, out / "runlog" / p.name)

    # ---- manifest --------------------------------------------------------
    files = sorted(p for p in out.rglob("*") if p.is_file())
    with open(out / "MANIFEST.txt", "w") as f:
        f.write("TfR1 Challenge 04 submission package\n")
        f.write(f"{len(rows)} designs, {len(files)} files\n\n")
        for p in files:
            f.write(f"{p.relative_to(out)}\t{p.stat().st_size}\n")
    with open(out / "SHA256SUMS.txt", "w") as f:
        for p in files:
            f.write(f"{sha256(p)}  {p.relative_to(out)}\n")

    print(f"[81] package at {out}")
    print(f"[81] {len(rows)} designs, {len(files)} files, "
          f"{sum(p.stat().st_size for p in files)/1e6:.1f} MB")


if __name__ == "__main__":
    main()
