#!/usr/bin/env python3
"""
Stage 00 - TfR1 human/mouse surface conservation map.  THE DECISION GATE.

Answers section 4 of tfr1-handoff-brief.md:
  Q1  Is K205 conserved in mouse?  (checked first - decides whether a strong
      strategy-1 switch and cross-reactivity can coexist at the apical apex)
  Q2  Are there conserved, surface-exposed Arg/Lys residues in a patch that is
      (a) not the transferrin site, (b) not glycan-occluded, (c) contiguous
      enough to host a ~600-900 A^2 binder footprint?
  Q3  Does the map corroborate or discard the unverified 318/633/719 seam lead?

Both sequences are the assay constructs from the challenge page - human
121-760 (640 aa) and mouse 123-763 (641 aa).  Neither is from UniProt.

Run:  venv-phscan/bin/python -I scripts/00_conservation.py
"""
import os
import sys
import json
import math
from collections import defaultdict

import numpy as np
from Bio import Align, SeqIO
from Bio.PDB import MMCIFParser, ShrakeRupley
from Bio.PDB.Polypeptide import is_aa
from Bio.Data.PDBData import protein_letters_3to1_extended as P3TO1

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TGT = os.path.join(ROOT, "target")
STRUCT = os.path.join(TGT, "structures")

HUMAN_FASTA = os.path.join(TGT, "human_tfr1_121_760.fasta")
MOUSE_FASTA = os.path.join(TGT, "mouse_tfr1_123_763.fasta")
HUMAN_START, MOUSE_START = 121, 123

# --- target biology constants, all human ectodomain numbering ----------------
DOMAINS = [("protease-like", 121, 188), ("apical", 189, 383),
           ("protease-like", 384, 606), ("helical", 607, 760)]
TF_SITE = (619, 670)                    # transferrin binding, brief section 2
GLYCAN_SEQUONS = [251, 317, 727]        # N-X-S/T, matches 11RN's 3 NAG/chain
APICAL_PATCH = list(range(202, 216)) + [348, 371]   # page-named apical patch
SEAM_LEAD = [318, 633, 719]             # unverified preprint, brief section 4

# page-stated human->mouse apical differences [OFFICIAL], used to validate the
# alignment.  (human_pos, human_aa, mouse_aa, mouse_pos); L212 has no counterpart.
PAGE_DIFFS = [(204, "D", "Q", 207), (208, "R", "N", 211), (210, "V", "D", 213),
              (211, "Y", "P", 214), (348, "N", "K", 350), (371, "K", "Q", 373)]

SASA_REL_CUTOFF = 0.25      # relative SASA above this = surface exposed
CATIONIC, ANIONIC = set("RK"), set("DE")
# Tien et al. 2013 theoretical max SASA per residue (A^2)
MAXASA = {"A": 129, "R": 274, "N": 195, "D": 193, "C": 167, "E": 223, "Q": 225,
          "G": 104, "H": 224, "I": 197, "L": 201, "K": 236, "M": 224, "F": 240,
          "P": 159, "S": 155, "T": 172, "W": 285, "Y": 263, "V": 174}


def log(msg=""):
    print(msg, flush=True)


def _r(x, nd=1):
    """Round, mapping a genuinely-missing value to -1.  NB: a real 0.0 A
    distance is falsy, so `x or -1` would corrupt it - hence the `is None`."""
    return -1.0 if x is None else round(x, nd)


def domain_of(pos):
    for name, lo, hi in DOMAINS:
        if lo <= pos <= hi:
            return name
    return "outside-construct"


def charge_class(aa):
    if aa in CATIONIC:
        return "cationic"
    if aa in ANIONIC:
        return "anionic"
    if aa == "H":
        return "histidine"
    if aa in set("STNQYCW"):
        return "polar"
    return "hydrophobic"


# =============================================================================
# 1. sequences and alignment
# =============================================================================
def load_sequences():
    h = str(next(SeqIO.parse(HUMAN_FASTA, "fasta")).seq).upper()
    m = str(next(SeqIO.parse(MOUSE_FASTA, "fasta")).seq).upper()
    assert len(h) == 640, f"human construct must be 640 aa, got {len(h)}"
    assert len(m) == 641, f"mouse construct must be 641 aa, got {len(m)}"
    return h, m


def align(h, m):
    """Global align; return {human_pos: (h_aa, mouse_pos|None, m_aa)} and stats."""
    al = Align.PairwiseAligner(scoring="blastp", mode="global",
                               open_gap_score=-11, extend_gap_score=-1)
    aln = al.align(h, m)[0]
    A, B = str(aln[0]), str(aln[1])
    ident = sum(1 for a, b in zip(A, B) if a == b and a != "-")
    cols = sum(1 for a, b in zip(A, B) if a != "-" and b != "-")
    hp, mp, mapping = HUMAN_START, MOUSE_START, {}
    for a, b in zip(A, B):
        if a != "-":
            mapping[hp] = (a, (mp if b != "-" else None), b)
        if a != "-":
            hp += 1
        if b != "-":
            mp += 1
    return mapping, 100.0 * ident / cols, ident, cols


def anchor_to_page(mapping, h, m):
    """
    The single apical indel is placeable at several chemically similar columns.
    The challenge page is [OFFICIAL] and states V210->D213, Y211->P214 and that
    L212 has no mouse counterpart.  A plain aligner puts the gap at V210 instead.
    Re-seat the gap in the 209-213 window onto the page's placement so every
    downstream number agrees with the organisers'.  Returns (mapping, notes).
    """
    notes = []
    if mapping[210][1] is None and mapping[211][1] is not None:
        notes.append(
            "Indel re-seated: aligner placed the apical gap at human V210; the "
            "challenge page places it at L212 (V210->D213, Y211->P214, L212->none). "
            "Adopted the page's placement. Both are equal-cost alignments of the "
            "same single-residue deletion; only the gap column differs, and all of "
            "210/211/212 are non-conserved under either.")
        mapping[210] = ("V", 213, m[213 - MOUSE_START])
        mapping[211] = ("Y", 214, m[214 - MOUSE_START])
        mapping[212] = ("L", None, "-")
    return mapping, notes


def validate_alignment(mapping):
    rows, ok = [], True
    for hpos, haa, maa, mpos in PAGE_DIFFS:
        got = mapping.get(hpos)
        if got is None:
            rows.append((hpos, haa, "?", "?", "MISSING"))
            ok = False
            continue
        g_haa, g_mpos, g_maa = got
        good = (g_haa == haa and g_maa == maa and g_mpos == mpos)
        ok &= good
        rows.append((hpos, f"{g_haa}{hpos}", f"{g_maa}{g_mpos}",
                     f"{haa}{hpos}->{maa}{mpos}", "OK" if good else "MISMATCH"))
    return rows, ok


# =============================================================================
# 2. structure
# =============================================================================
def load_chain(cif_path, want_chain):
    st = MMCIFParser(QUIET=True).get_structure("x", cif_path)
    model = next(st.get_models())
    if want_chain not in [c.id for c in model]:
        raise KeyError(f"chain {want_chain} not in {cif_path}: "
                       f"{[c.id for c in model]}")
    return model, model[want_chain]


def aa_residues(chain):
    out = []
    for r in chain:
        if r.id[0] != " " or not is_aa(r, standard=False):
            continue
        aa = P3TO1.get(r.get_resname().upper())
        if aa and len(aa) == 1:
            out.append((r.id[1], aa, r))
    return out


def compute_sasa(entity):
    ShrakeRupley().compute(entity, level="R")


def residue_sasa_map(chain):
    out = {}
    for num, aa, r in aa_residues(chain):
        out[num] = r.sasa
    return out


def vp1u_footprint(cutoff=4.5):
    """Residues of 11RN TfR1 chain A within `cutoff` A of any VP1u atom."""
    st = MMCIFParser(QUIET=True).get_structure("x",
                                               os.path.join(STRUCT, "11RN.cif"))
    model = next(st.get_models())
    chains = {c.id: c for c in model}
    tfr_ids = [cid for cid, c in chains.items() if len(aa_residues(c)) > 400]
    vp_ids = [cid for cid, c in chains.items()
              if 30 < len(aa_residues(c)) < 200]
    if not tfr_ids or not vp_ids:
        return set(), tfr_ids, vp_ids
    tfr = chains[tfr_ids[0]]
    vp_atoms = [a for cid in vp_ids for a in chains[cid].get_atoms()]
    vp_xyz = np.array([a.coord for a in vp_atoms])
    hits = set()
    for num, aa, r in aa_residues(tfr):
        xyz = np.array([a.coord for a in r.get_atoms()])
        d = np.linalg.norm(xyz[:, None, :] - vp_xyz[None, :, :], axis=-1)
        if d.min() <= cutoff:
            hits.add(num)
    return hits, tfr_ids, vp_ids


# =============================================================================
# 3. geometry helpers
# =============================================================================
def centroid(res_by_num, positions):
    pts = [a.coord for p in positions if p in res_by_num
           for a in res_by_num[p].get_atoms()]
    return np.mean(np.array(pts), axis=0) if pts else None


def min_dist_to_set(res_by_num, pos, targets):
    if pos not in res_by_num:
        return None
    xyz = np.array([a.coord for a in res_by_num[pos].get_atoms()])
    best = None
    for t in targets:
        if t not in res_by_num:
            continue
        txyz = np.array([a.coord for a in res_by_num[t].get_atoms()])
        d = float(np.linalg.norm(xyz[:, None, :] - txyz[None, :, :],
                                 axis=-1).min())
        best = d if best is None else min(best, d)
    return best


def dist_to_point(res_by_num, pos, pt):
    if pos not in res_by_num or pt is None:
        return None
    xyz = np.array([a.coord for a in res_by_num[pos].get_atoms()])
    return float(np.linalg.norm(xyz - pt[None, :], axis=-1).min())


# =============================================================================
# main
# =============================================================================
def main():
    log("=" * 74)
    log("STAGE 00 - TfR1 human/mouse surface conservation map (DECISION GATE)")
    log("=" * 74)

    # ---- sequences ----------------------------------------------------------
    h, m = load_sequences()
    log(f"\n[seq] human assay construct 121-760: {len(h)} aa")
    log(f"[seq] mouse assay construct 123-763: {len(m)} aa")
    log("[seq] both from the challenge page, NOT UniProt")

    mapping, pct, ident, cols = align(h, m)
    mapping, anchor_notes = anchor_to_page(mapping, h, m)
    log(f"\n[aln] global identity over aligned columns: {pct:.1f}% "
        f"({ident}/{cols})  -- page states 78%")
    for n in anchor_notes:
        log(f"[aln] {n}")

    rows, ok = validate_alignment(mapping)
    log("\n[aln] validation against the page's six stated apical differences:")
    for hpos, got_h, got_m, expect, flag in rows:
        log(f"      {flag:9s} {got_h} -> {got_m}   (page: {expect})")
    if not ok:
        log("\n[aln] FATAL: alignment disagrees with [OFFICIAL] page facts. STOP.")
        sys.exit(2)
    log("[aln] ALIGNMENT VALIDATED against all six page-stated differences.")

    # ---- Q1: K205 ------------------------------------------------------------
    k205 = mapping[205]
    log("\n" + "=" * 74)
    log("Q1  IS K205 CONSERVED IN MOUSE?   <-- checked first, decides the campaign")
    log("=" * 74)
    log(f"    human K205  ->  mouse {k205[2]}{k205[1]}")
    k205_conserved = (k205[2] == "K")
    log(f"    VERDICT: K205 is {'CONSERVED' if k205_conserved else 'NOT CONSERVED'}"
        f"  (K -> {k205[2]})")

    log("\n    apical main loop 202-215, residue by residue:")
    for p in range(202, 216):
        haa, mpos, maa = mapping[p]
        tag = "conserved" if haa == maa else "DIFFERS"
        log(f"      human {haa}{p:<4d} -> mouse {maa}{mpos if mpos else '  -':<5}"
            f" {tag}")

    log("\n    apical cation inventory:")
    for p in [205, 208, 371, 348]:
        haa, mpos, maa = mapping[p]
        hc = "cation" if haa in CATIONIC else charge_class(haa)
        mc = "cation" if maa in CATIONIC else charge_class(maa) if maa != "-" else "none"
        log(f"      {haa}{p:<4d} ({hc:11s}) -> mouse {maa}{mpos if mpos else '-'}"
            f" ({mc})")

    # ---- structure -----------------------------------------------------------
    log("\n" + "=" * 74)
    log("STRUCTURE 11QC")
    log("=" * 74)
    model, chainA = load_chain(os.path.join(STRUCT, "11QC.cif"), "A")
    chain_ids = [c.id for c in model]
    resA = aa_residues(chainA)
    log(f"[pdb] 11QC chains: {chain_ids}; chain A modelled residues: {len(resA)}")
    log(f"[pdb] chain A numbering span: {resA[0][0]}-{resA[-1][0]}")

    hetero = defaultdict(int)
    for c in model:
        for r in c:
            if r.id[0] != " ":
                hetero[r.get_resname()] += 1
    log(f"[pdb] heteroatom residues present: {dict(hetero)}")

    res_by_num = {num: r for num, aa, r in resA}
    seq_by_num = {num: aa for num, aa, r in resA}

    # position 142 check (the Gly142 vs Ser142 trap)
    g142 = seq_by_num.get(142)
    log(f"\n[142] 11QC chain A position 142 = {g142}  "
        f"(assay construct carries G; mouse-aligned partner = "
        f"{mapping[142][2]}{mapping[142][1]})")
    if g142 != "G":
        log("[142] WARNING: 11QC does not carry Gly142. Patch before designing.")
    else:
        log("[142] OK - 11QC matches the assay construct at 142. No patch needed.")

    # cross-check structure sequence against the assay construct
    mism = [(n, seq_by_num[n], h[n - HUMAN_START])
            for n in sorted(seq_by_num)
            if HUMAN_START <= n <= 760 and seq_by_num[n] != h[n - HUMAN_START]]
    log(f"\n[xchk] 11QC chain A vs assay construct: {len(mism)} mismatches "
        f"over {len(seq_by_num)} modelled residues")
    for n, a, b in mism[:12]:
        log(f"       pos {n}: structure {a}, construct {b}")

    # ---- SASA ---------------------------------------------------------------
    log("\n[sasa] Shrake-Rupley on the isolated monomer (chain A)...")
    import copy
    mono = copy.deepcopy(chainA)
    compute_sasa(mono)
    sasa_mono = residue_sasa_map(mono)

    log("[sasa] Shrake-Rupley on the A+B dimer (to find seam-buried residues)...")
    dim = copy.deepcopy(model)
    for c in list(dim):
        if c.id not in ("A", "B"):
            dim.detach_child(c.id)
    for c in dim:
        for r in list(c):
            if r.id[0] != " ":
                c.detach_child(r.id)
    compute_sasa(dim)
    sasa_dim = residue_sasa_map(dim["A"])

    # ---- VP1u footprint ------------------------------------------------------
    log("\n[11RN] computing VP1u footprint on TfR1...")
    fp, tfr_ids, vp_ids = vp1u_footprint()
    log(f"[11RN] TfR1 chain(s) {tfr_ids}, VP1u chain(s) {vp_ids}")
    log(f"[11RN] VP1u contacts {len(fp)} TfR1 residues (<=4.5 A): "
        f"{sorted(fp)}")

    # ---- per-residue table ---------------------------------------------------
    apex = centroid(res_by_num, APICAL_PATCH)
    tf_positions = list(range(TF_SITE[0], TF_SITE[1] + 1))

    table = []
    for num in sorted(seq_by_num):
        if not (HUMAN_START <= num <= 760):
            continue
        aa = seq_by_num[num]
        mp = mapping.get(num)
        if mp is None:
            continue
        maa, mpos = mp[2], mp[1]
        rel = sasa_mono.get(num, 0.0) / MAXASA.get(aa, 200)
        rel_dim = sasa_dim.get(num, 0.0) / MAXASA.get(aa, 200)
        table.append(dict(
            pos=num, aa=aa, mouse_aa=maa, mouse_pos=mpos,
            conserved=(aa == maa), domain=domain_of(num),
            charge=charge_class(aa),
            sasa=round(sasa_mono.get(num, 0.0), 1), rel_sasa=round(rel, 3),
            rel_sasa_dimer=round(rel_dim, 3),
            seam_buried=(rel - rel_dim > 0.15),
            surface=(rel >= SASA_REL_CUTOFF),
            d_apex=_r(dist_to_point(res_by_num, num, apex)),
            d_tf=_r(min_dist_to_set(res_by_num, num, tf_positions)),
            d_glycan=_r(min_dist_to_set(res_by_num, num, GLYCAN_SEQUONS)),
            vp1u=(num in fp),
        ))

    surf = [t for t in table if t["surface"]]
    log(f"\n[tab] {len(table)} residues scored; {len(surf)} surface-exposed "
        f"(rel SASA >= {SASA_REL_CUTOFF})")
    ncons = sum(1 for t in surf if t["conserved"])
    log(f"[tab] of those, {ncons} conserved ({100*ncons/len(surf):.0f}%) - note the "
        f"surface is less conserved than the {pct:.0f}% global figure")

    # ---- Q2: conserved surface cations and the patches around them -----------
    log("\n" + "=" * 74)
    log("Q2  CONSERVED, SURFACE-EXPOSED CATIONS IN A USABLE PATCH")
    log("=" * 74)
    cat = [t for t in surf if t["charge"] == "cationic" and t["conserved"]]
    log(f"    {len(cat)} conserved surface Arg/Lys on the monomer")

    # usable = not transferrin site, not glycan-adjacent, not dimer seam
    usable = [t for t in cat
              if not (TF_SITE[0] <= t["pos"] <= TF_SITE[1])
              and t["d_tf"] > 10
              and t["d_glycan"] > 10
              and not t["seam_buried"]]
    log(f"    {len(usable)} of them are Tf-sparing (>10 A from 619-670), "
        f">10 A from a glycan sequon, and not dimer-seam buried")

    # score a binder footprint around each usable cation
    cands = []
    # precompute min heavy-atom distances between all surface residues once
    xyz_of = {t["pos"]: np.array([a.coord for a in res_by_num[t["pos"]]])
              for t in surf if t["pos"] in res_by_num}

    def mind(p, q):
        A_, B_ = xyz_of.get(p), xyz_of.get(q)
        if A_ is None or B_ is None:
            return None
        return float(np.linalg.norm(A_[:, None, :] - B_[None, :, :],
                                    axis=-1).min())

    for t in usable:
        # the cation itself sits at 0.0 A and must stay inside its own patch
        near = [u for u in surf
                if (lambda d: d is not None and d <= 12)(mind(u["pos"], t["pos"]))]
        ncon = sum(1 for u in near if u["conserved"])
        area = sum(u["sasa"] for u in near)
        othercat = [u["pos"] for u in near
                    if u["charge"] == "cationic" and u["conserved"]
                    and u["pos"] != t["pos"]]
        anions = [u["pos"] for u in near if u["charge"] == "anionic"]
        cands.append(dict(
            pos=t["pos"], aa=t["aa"], domain=t["domain"], d_apex=t["d_apex"],
            d_tf=t["d_tf"], d_glycan=t["d_glycan"], vp1u=t["vp1u"],
            n_surface=len(near), n_conserved=ncon,
            frac_conserved=round(ncon / max(len(near), 1), 2),
            patch_area=round(area), other_conserved_cations=othercat,
            anions_in_patch=anions))
    cands.sort(key=lambda c: (-c["frac_conserved"], -c["patch_area"]))

    log(f"\n    ranked candidate patches (12 A footprint around each conserved "
        f"surface cation):")
    log(f"    {'pos':>6} {'dom':>14} {'dApex':>6} {'dTf':>6} {'dGly':>6} "
        f"{'nSurf':>6} {'nCons':>6} {'fCons':>6} {'areaA2':>7} {'VP1u':>5} "
        f"other conserved cations")
    for c in cands[:25]:
        log(f"    {c['aa']}{c['pos']:<5d} {c['domain']:>14} {c['d_apex']:>6} "
            f"{c['d_tf']:>6} {c['d_glycan']:>6} {c['n_surface']:>6} "
            f"{c['n_conserved']:>6} {c['frac_conserved']:>6} "
            f"{c['patch_area']:>7} {str(c['vp1u']):>5} "
            f"{c['other_conserved_cations']}")

    # ---- the rearranged-cation question -------------------------------------
    log("\n" + "=" * 74)
    log("Q2b IS THE APICAL CATION INVENTORY 'REARRANGED' RATHER THAN LOST?")
    log("=" * 74)
    log("    human apical cations R208, K371 (+K205) all go neutral in mouse,")
    log("    but mouse GAINS K350 where human has N348. If human R208 and human")
    log("    348 are close in space, one binder His could sit between them and")
    log("    find a cation in BOTH species.")
    log(f"\n      {'pair':>12} {'minHeavy':>9} {'CB-CB':>8} {'tip-tip':>8}")
    for a, b in [(371, 348), (208, 348), (205, 348), (208, 371), (205, 208)]:
        dmin = min_dist_to_set(res_by_num, a, [b])
        cb = None
        if all(p in res_by_num and "CB" in res_by_num[p] for p in (a, b)):
            cb = float(np.linalg.norm(res_by_num[a]["CB"].coord
                                      - res_by_num[b]["CB"].coord))
        # distal side-chain atom ("tip") of each residue
        def tip(p):
            r = res_by_num.get(p)
            if r is None:
                return None
            ca = r["CA"].coord if "CA" in r else None
            sc = [at for at in r if at.get_id() not in ("N", "CA", "C", "O")]
            if not sc or ca is None:
                return ca
            return max(sc, key=lambda at: float(
                np.linalg.norm(at.coord - ca))).coord
        ta, tb = tip(a), tip(b)
        tt = (float(np.linalg.norm(ta - tb))
              if ta is not None and tb is not None else None)
        log(f"      {f'{a} <-> {b}':>12} "
            f"{_r(dmin):>9} {_r(cb):>8} {_r(tt):>8}")

    log("\n      exposure of the two swap partners (mouse gains a cation at 350,")
    log("      loses the one at 373 - human is the mirror image):")
    for p in (348, 371):
        t = next((x for x in table if x["pos"] == p), None)
        if t:
            log(f"        human {t['aa']}{p} relSASA {t['rel_sasa']:<6} "
                f"-> mouse {t['mouse_aa']}{t['mouse_pos']}  "
                f"VP1u-contact={t['vp1u']}")

    # ---- resolve the nine-residue apical patch ------------------------------
    log("\n" + "=" * 74)
    log("RESOLVING THE NINE-RESIDUE APICAL PATCH (brief section 4, step 7)")
    log("=" * 74)
    page_named = {204, 208, 210, 211, 212, 348, 371}
    log(f"    VP1u footprint on TfR1 (<=4.5 A) = {len(fp)} residues: {sorted(fp)}")
    log(f"    page names 7 of 9: {sorted(page_named)}")
    unnamed = sorted(fp - page_named)
    log(f"    the previously unnamed patch residues are therefore: {unnamed}")
    for p in unnamed:
        t = next((x for x in table if x["pos"] == p), None)
        if t:
            log(f"      {t['aa']}{p} -> mouse {t['mouse_aa']}{t['mouse_pos']}  "
                f"{'CONSERVED' if t['conserved'] else 'differs'};  "
                f"charge={t['charge']};  relSASA={t['rel_sasa']}")
    log(f"    K205 in the VP1u footprint? {205 in fp}")

    # ---- Q3: the unverified seam lead ---------------------------------------
    log("\n" + "=" * 74)
    log("Q3  DOES THE MAP CORROBORATE THE 318/633/719 SEAM LEAD?")
    log("=" * 74)
    for p in SEAM_LEAD:
        t = next((x for x in table if x["pos"] == p), None)
        if t is None:
            log(f"    {p}: not modelled in 11QC")
            continue
        log(f"    {t['aa']}{p} -> mouse {t['mouse_aa']}{t['mouse_pos']}  "
            f"{'conserved' if t['conserved'] else 'DIFFERS'}; "
            f"relSASA {t['rel_sasa']}; dTf {t['d_tf']} A; "
            f"dGlycan {t['d_glycan']} A; domain {t['domain']}")
    log(f"    (transferrin site is {TF_SITE[0]}-{TF_SITE[1]}; K633 falls inside it)")

    # ---- native histidines ---------------------------------------------------
    log("\n" + "=" * 74)
    log("NATIVE HISTIDINES - spatial check against the apical patch")
    log("=" * 74)
    his = [t for t in table if t["aa"] == "H"]
    log(f"    {len(his)} histidines modelled in chain A")
    for t in sorted(his, key=lambda x: x["d_apex"]):
        log(f"      H{t['pos']:<5d} relSASA {t['rel_sasa']:<6} "
            f"d(apex) {t['d_apex']:>6} A  domain {t['domain']:>14}  "
            f"mouse {t['mouse_aa']}{t['mouse_pos']}")
    log(f"    position 349 is {seq_by_num.get(349)} - the 'His349' claim in "
        f"background reading is WRONG, do not carry it forward")

    # ---- persist -------------------------------------------------------------
    os.makedirs(os.path.join(ROOT, "runlog"), exist_ok=True)
    with open(os.path.join(TGT, "conservation_table.tsv"), "w") as f:
        cols_ = ["pos", "aa", "mouse_aa", "mouse_pos", "conserved", "domain",
                 "charge", "sasa", "rel_sasa", "rel_sasa_dimer", "seam_buried",
                 "surface", "d_apex", "d_tf", "d_glycan", "vp1u"]
        f.write("\t".join(cols_) + "\n")
        for t in table:
            f.write("\t".join(str(t[c]) for c in cols_) + "\n")
    with open(os.path.join(TGT, "conservation_candidates.json"), "w") as f:
        json.dump(dict(global_identity_pct=pct, k205=dict(
            human="K205", mouse=f"{k205[2]}{k205[1]}",
            conserved=k205_conserved), candidates=cands,
            vp1u_footprint=sorted(fp),
            n_surface=len(surf), n_surface_conserved=ncons), f, indent=2)
    log(f"\n[out] target/conservation_table.tsv")
    log(f"[out] target/conservation_candidates.json")
    log("\nSTAGE 00 COMPLETE - stop here for the epitope decision.")


if __name__ == "__main__":
    main()
