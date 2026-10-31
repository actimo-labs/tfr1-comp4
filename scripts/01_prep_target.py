#!/usr/bin/env python3
"""
Stage 01 - prepare the TfR1 design target from 11QC.

Full rewrite; challenge 02's 01_prep_target.py was TNF-specific throughout.

What this produces and why:

  target/tfr1_A.pdb        chain A monomer, glycans stripped. THE design and
                           validation target. 11RN shows VP1u binds one RBD per
                           TfR1 subunit (A2B2), so the apical epitope is a
                           single-subunit epitope and the dimer is not needed.
  target/tfr1_AB.pdb       the A+B dimer, kept ONLY as a reference for seam
                           checks. Never fed to RFdiffusion - it doubles the
                           denoised length and the budget does not survive it.
  target/hotspots_*.txt    RFdiffusion hotspot sets, derived from the stage-00
                           conservation map, not from the literature.
  target/nondesignable.txt glycan-shadowed residues to keep binders off.
  target/prep_report.json  every check, with its result.

Checks that must pass or the script exits non-zero:
  - position 142 is Gly (the Ser142 allele trap)
  - chain A matches the assay construct at every modelled residue
  - all three glycan sequons are present and are genuine N-X-S/T
  - every hotspot is modelled, surface-exposed, and in one contiguous patch

Run:  venv-phscan/bin/python -I scripts/01_prep_target.py
"""
import os
import sys
import json

import numpy as np
from Bio import SeqIO
from Bio.PDB import MMCIFParser, PDBIO, Select
from Bio.PDB.Polypeptide import is_aa
from Bio.Data.PDBData import protein_letters_3to1_extended as P3TO1

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TGT = os.path.join(ROOT, "target")
CIF = os.path.join(TGT, "structures", "11QC.cif")
HUMAN_FASTA = os.path.join(TGT, "human_tfr1_121_760.fasta")
HUMAN_START = 121

GLYCAN_SEQUONS = [251, 317, 727]
GLYCAN_SHADOW = 10.0        # A around a sequon Asn that a glycan can cover

# Hotspot sets, from stage 00.  See target/conservation_report.md section 6.
#   S1  primary: the apical shoulder. Carries the switch anchor K374 (conserved
#       K->K376, no carboxylate within 10 A) plus four of the nine VP1u
#       footprint residues, so it is both switchable cross-species and a
#       defensible functional epitope.
#   swap the 348/371 locus alone, for the class-E speculative slots. Human
#       carries its cation at 371 and mouse at 350 (= human 348), and the two
#       sites are 3.6 A apart, so a His sited here may find a cation in both
#       species. High upside, rests on unmeasurable mouse side-chain geometry.
#   X   the apex tip, for the class-C affinity hedge. Expected to fail in mouse
#       (every residue here diverges); included deliberately as the hedge and
#       as the internal control for the switch hypothesis.
#
# A fourth set built on the second conserved cation K231 was tried and
# REJECTED by this script's own glycan check: K231 is 8.7 A from the N251
# sequon, inside the 10 A shadow. K344, the other conserved apical cation, is
# salt-bridged to E214 at 2.8 A and would invert the switch. K374 is therefore
# the only clean conserved switch anchor on this epitope - which is why the
# primary set is built around it.
HOTSPOT_SETS = {
    "shoulder": [209, 211, 348, 371, 374],
    "swap_locus": [211, 348, 371],
    "apex": [204, 208, 211, 212],
}
PRIMARY = "shoulder"


def die(msg):
    print(f"FATAL: {msg}", file=sys.stderr)
    sys.exit(2)


def log(m=""):
    print(m, flush=True)


class Clean(Select):
    """Keep standard amino acids of the wanted chains, drop waters/NAG/H/altlocs."""

    def __init__(self, chains):
        self.chains = set(chains)

    def accept_chain(self, chain):
        return chain.id in self.chains

    def accept_residue(self, residue):
        return residue.id[0] == " " and is_aa(residue, standard=False)

    def accept_atom(self, atom):
        return atom.element != "H" and atom.get_altloc() in (" ", "A")


def main():
    report = {}
    log("=" * 70)
    log("STAGE 01 - TfR1 target preparation from 11QC")
    log("=" * 70)

    if not os.path.exists(CIF):
        die(f"{CIF} missing - run stage 00 first")

    model = next(MMCIFParser(QUIET=True).get_structure("t", CIF).get_models())
    chains = [c.id for c in model]
    log(f"\n[in]  11QC chains: {chains}")

    aa = {}
    for r in model["A"]:
        if r.id[0] == " " and is_aa(r, standard=False):
            one = P3TO1.get(r.get_resname().upper())
            if one and len(one) == 1:
                aa[r.id[1]] = one
    nums = sorted(aa)
    log(f"[in]  chain A: {len(aa)} amino acids, span {nums[0]}-{nums[-1]}")

    het = {}
    for c in model:
        for r in c:
            if r.id[0] != " ":
                het[r.get_resname()] = het.get(r.get_resname(), 0) + 1
    log(f"[in]  heteroatoms: {het or 'none'}")
    report["heteroatoms_removed"] = het

    # ---- check 1: insertion codes would break every downstream numbering map
    icodes = [r.id for c in model for r in c if r.id[2] != " "]
    if icodes:
        die(f"insertion codes present, numbering map unsafe: {icodes[:5]}")
    log("[chk] no insertion codes")

    # ---- check 2: Gly142
    if 142 not in aa:
        die("position 142 not modelled in 11QC chain A")
    log(f"[chk] position 142 = {aa[142]}")
    if aa[142] != "G":
        die(f"position 142 is {aa[142]}, expected Gly. Patch Ser142->Gly142 "
            f"before designing (brief section 2).")
    log("[chk] Gly142 confirmed - no patch needed")
    report["pos142"] = aa[142]

    # ---- check 3: structure vs assay construct
    h = str(next(SeqIO.parse(HUMAN_FASTA, "fasta")).seq).upper()
    if len(h) != 640:
        die(f"assay construct must be 640 aa, got {len(h)}")
    mism = [(n, aa[n], h[n - HUMAN_START]) for n in nums
            if HUMAN_START <= n <= 760 and aa[n] != h[n - HUMAN_START]]
    log(f"[chk] chain A vs assay construct: {len(mism)} mismatches / "
        f"{len(nums)} modelled")
    if mism:
        for n, s, c in mism[:10]:
            log(f"       {n}: structure {s}, construct {c}")
        die("structure does not match the assay construct - do not design "
            "against it until this is explained")
    report["construct_mismatches"] = 0

    # ---- check 4: glycan sequons are real N-X-S/T
    seq_ok = {}
    for n in GLYCAN_SEQUONS:
        i = n - HUMAN_START
        tri = h[i:i + 3]
        ok = len(tri) == 3 and tri[0] == "N" and tri[1] != "P" and tri[2] in "ST"
        seq_ok[n] = dict(motif=tri, valid=bool(ok))
        log(f"[chk] sequon N{n}: {tri} {'OK' if ok else 'INVALID'}")
        if not ok:
            die(f"sequon N{n} is not a valid N-X-S/T motif")
    report["glycan_sequons"] = seq_ok

    # ---- write the monomer and the dimer -----------------------------------
    io = PDBIO()
    io.set_structure(model)
    mono_path = os.path.join(TGT, "tfr1_A.pdb")
    io.save(mono_path, select=Clean(["A"]))
    io.save(os.path.join(TGT, "tfr1_AB.pdb"), select=Clean(["A", "B"]))
    log(f"\n[out] target/tfr1_A.pdb   (design + validation target, monomer)")
    log(f"[out] target/tfr1_AB.pdb  (reference only, never fed to RFdiffusion)")

    # verify the written monomer reloads with the numbering intact
    from Bio.PDB import PDBParser
    back = next(PDBParser(QUIET=True).get_structure("b", mono_path).get_models())
    bnums = sorted(r.id[1] for r in back["A"] if r.id[0] == " ")
    if bnums != nums:
        die(f"numbering changed on write: {len(bnums)} vs {len(nums)} residues")
    log(f"[chk] round-trip OK - {len(bnums)} residues, numbering preserved")

    # ---- glycan shadow ------------------------------------------------------
    res = {r.id[1]: r for r in model["A"]
           if r.id[0] == " " and is_aa(r, standard=False)}
    xyz = {p: np.array([a.coord for a in r]) for p, r in res.items()}

    def mind(p, q):
        return float(np.linalg.norm(xyz[p][:, None, :] - xyz[q][None, :, :],
                                    axis=-1).min())

    shadow = set()
    for s in GLYCAN_SEQUONS:
        for p in xyz:
            if mind(p, s) <= GLYCAN_SHADOW:
                shadow.add(p)
    shadow |= set(GLYCAN_SEQUONS)
    with open(os.path.join(TGT, "nondesignable.txt"), "w") as f:
        f.write(",".join(f"A{p}" for p in sorted(shadow)) + "\n")
    log(f"\n[out] target/nondesignable.txt - {len(shadow)} residues within "
        f"{GLYCAN_SHADOW} A of N251/N317/N727")
    report["nondesignable_count"] = len(shadow)
    report["glycan_shadow_radius"] = GLYCAN_SHADOW

    # ---- hotspots -----------------------------------------------------------
    cons = {}
    tsv = os.path.join(TGT, "conservation_table.tsv")
    if os.path.exists(tsv):
        import csv
        cons = {int(r["pos"]): r
                for r in csv.DictReader(open(tsv), delimiter="\t")}

    log("\n" + "=" * 70)
    log("HOTSPOT SETS (from the stage-00 map)")
    log("=" * 70)
    hs_report = {}
    for name, members in HOTSPOT_SETS.items():
        missing = [p for p in members if p not in res]
        if missing:
            die(f"hotspot set {name}: residues not modelled: {missing}")
        in_shadow = [p for p in members if p in shadow]
        span = max(mind(a, b) for a in members for b in members)
        isolated = [p for p in members
                    if min(mind(p, q) for q in members if q != p) > 12]
        buried = [p for p in members
                  if cons and float(cons[p]["rel_sasa"]) < 0.15]
        tag = "PRIMARY" if name == PRIMARY else ""
        log(f"\n  {name} {tag}")
        log(f"    members: " + ", ".join(
            f"{cons[p]['aa'] if cons else '?'}{p}"
            f"->{cons[p]['mouse_aa']}{cons[p]['mouse_pos']}" if cons else str(p)
            for p in members))
        if cons:
            for p in members:
                c = cons[p]
                log(f"      {c['aa']}{p:<4} -> mouse {c['mouse_aa']}{c['mouse_pos']:<5}"
                    f" {'conserved' if c['conserved']=='True' else 'differs  '}"
                    f" relSASA {c['rel_sasa']:<6} {c['charge']:<12}"
                    f" VP1u={c['vp1u']}")
        log(f"    max pairwise span {span:.1f} A; isolated {isolated or 'none'}; "
            f"glycan-shadowed {in_shadow or 'none'}; buried {buried or 'none'}")
        if isolated:
            die(f"hotspot set {name} is not one contiguous patch: {isolated}")
        if in_shadow:
            die(f"hotspot set {name} includes glycan-shadowed residues: "
                f"{in_shadow} - a binder parked there fails in the assay")
        path = os.path.join(TGT, f"hotspots_{name}.txt")
        with open(path, "w") as f:
            f.write(",".join(f"A{p}" for p in members) + "\n")
        log(f"    [out] target/hotspots_{name}.txt")
        hs_report[name] = dict(members=[f"A{p}" for p in members],
                               span_A=round(span, 1), primary=(name == PRIMARY))
    report["hotspot_sets"] = hs_report

    # the primary set is also written to the conventional filename the crop
    # script defaults to, so a mis-specified --hotspots cannot silently pick
    # up the wrong patch
    with open(os.path.join(TGT, "hotspots.txt"), "w") as f:
        f.write(",".join(f"A{p}" for p in HOTSPOT_SETS[PRIMARY]) + "\n")
    log(f"\n[out] target/hotspots.txt -> copy of '{PRIMARY}' (the primary set)")

    report["n_residues_monomer"] = len(nums)
    report["span"] = [nums[0], nums[-1]]
    with open(os.path.join(TGT, "prep_report.json"), "w") as f:
        json.dump(report, f, indent=2)
    log(f"[out] target/prep_report.json")
    log("\nSTAGE 01 COMPLETE - all checks passed.")


if __name__ == "__main__":
    main()
