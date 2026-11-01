#!/usr/bin/env python
"""Stage 04 -- Boltz-2 validation of designed binders against TfR1.

Ported from challenge 02. Four substantive changes:

1. THE TARGET IS A MONOMER. 11RN is A2B2 with one VP1u RBD per TfR1 subunit,
   so the apical epitope is a single-subunit epitope. There is exactly one
   binder-target chain pair, which removes challenge 02's "minimum across
   pairs" rule and its target_pair_iptm assembly gate -- there is no
   target-internal pair to go wrong. The non-negotiable still holds in spirit:
   we read pair_chains_iptm for the binder-target pair, never global ipTM.

2. NO HARDCODED TARGET SEQUENCE. Challenge 02 pinned the target sequence as a
   module constant, which the brief's port table flags as the thing to swap.
   Instead both target sequences are DERIVED at run time -- human from the
   prepared structure itself, mouse by alignment onto it -- so the sequence fed
   to Boltz cannot drift from the structure used to template it.

3. --target-mode {full,apical}. `full` is the 639-residue monomer; `apical` is
   the contiguous 189-383 apical domain. The 18 A generation crop is NOT usable
   here: it has five discontinuous segments and feeding it as one chain would
   invent peptide bonds across the gaps.

4. --species {human,mouse}. Mouse is a SECOND-STAGE filter on human survivors,
   never a parallel arm (brief section 5). Mouse scores are against a sequence
   threaded onto the human structure and carry more uncertainty than human
   ones; the methods report must say so.

Composition pre-filter: designs that stage 08 would reject anyway (free Cys,
N-X-S/T sequon, homopolymer run >=5) are dropped BEFORE folding rather than
after. ProteinMPNN's poly-Ala/poly-Glu bias produces these readily and folding
them is wasted GPU.

  conda run -n boltz python scripts/04_boltz.py --fasta-dir pilot/mpnn/plain/seqs
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

BINDER_ID = "X"
TARGET_ID = "A"

APICAL = (189, 383)          # contiguous apical domain, all hotspots inside
HUMAN_START = 121            # assay construct numbering


# ---------------------------------------------------------------- sequences
def human_from_structure(pdb: Path) -> dict[int, str]:
    """Read the prepared monomer and return {author_resnum: one_letter}."""
    from Bio.PDB import PDBParser
    from Bio.PDB.Polypeptide import is_aa
    from Bio.Data.PDBData import protein_letters_3to1_extended as P3
    model = next(PDBParser(QUIET=True).get_structure("t", str(pdb)).get_models())
    out = {}
    for r in model[TARGET_ID]:
        if r.id[0] == " " and is_aa(r, standard=False):
            aa = P3.get(r.get_resname().upper())
            if aa and len(aa) == 1:
                out[r.id[1]] = aa
    return out


def mouse_by_alignment(human_seq_by_num: dict[int, str]) -> dict[int, str]:
    """Thread the mouse construct onto human numbering.

    Returns {human_resnum: mouse_one_letter}; positions with no mouse
    counterpart (the single apical indel) are omitted, which deletes that
    residue from the mouse target rather than inventing one.
    """
    from Bio import Align, SeqIO
    h = str(next(SeqIO.parse(ROOT / "target/human_tfr1_121_760.fasta",
                             "fasta")).seq).upper()
    m = str(next(SeqIO.parse(ROOT / "target/mouse_tfr1_123_763.fasta",
                             "fasta")).seq).upper()
    al = Align.PairwiseAligner(scoring="blastp", mode="global",
                               open_gap_score=-11, extend_gap_score=-1)
    A, B = (str(x) for x in al.align(h, m)[0])
    hp, out = HUMAN_START, {}
    for a, b in zip(A, B):
        if a != "-":
            if b != "-":
                out[hp] = b
            hp += 1
    return {n: out[n] for n in human_seq_by_num if n in out}


def target_sequence(mode: str, species: str) -> tuple[str, Path]:
    """Return (sequence, template_cif) for the requested target."""
    pdb = ROOT / ("target/tfr1_A_apical.pdb" if mode == "apical"
                  else "target/tfr1_A.pdb")
    cif = ROOT / ("target/tfr1_A_apical.cif" if mode == "apical"
                  else "target/tfr1_A.cif")
    by_num = human_from_structure(pdb)
    if species == "mouse":
        by_num = mouse_by_alignment(by_num)
    seq = "".join(by_num[n] for n in sorted(by_num))
    return seq, cif


# ---------------------------------------------------------------- filtering
SEQUON = re.compile(r"N[^P][ST]")
HOMOPOLYMER = re.compile(r"(.)\1{4,}")


def composition_reject(seq: str) -> str | None:
    """Reasons stage 08 would reject this design anyway. None = keep."""
    if "C" in seq:
        return "free Cys"
    if SEQUON.search(seq):
        return f"N-X-S/T sequon at {SEQUON.search(seq).start()}"
    m = HOMOPOLYMER.search(seq)
    if m:
        return f"homopolymer {m.group(1)}x{len(m.group(0))}"
    return None


def read_binders(fasta_dir: Path, per_backbone: int | None) -> dict[str, str]:
    """MPNN fastas: record 0 is the input backbone, records 1..n are designs."""
    out: dict[str, str] = {}
    for fa in sorted(fasta_dir.glob("*.fa")):
        recs, name, buf = [], None, []
        for line in fa.read_text().splitlines():
            if line.startswith(">"):
                if name is not None:
                    recs.append((name, "".join(buf)))
                name, buf = line[1:], []
            elif line.strip():
                buf.append(line.strip())
        if name is not None:
            recs.append((name, "".join(buf)))
        designs = recs[1:]
        if per_backbone:
            designs = designs[:per_backbone]
        for hdr, seq in designs:
            m = re.search(r"sample=(\d+)", hdr)
            idx = m.group(1) if m else str(len(out))
            out[f"{fa.stem}_s{idx}"] = seq.split("/")[0]
    return out


def write_yaml(path: Path, binder: str, target: str, template: Path | None):
    lines = ["version: 1", "sequences:"]
    for cid, seq in ((BINDER_ID, binder), (TARGET_ID, target)):
        lines += ["  - protein:", f"      id: {cid}",
                  f"      sequence: {seq}", "      msa: empty"]
    if template:
        lines += ["templates:", f"  - cif: {template}",
                  f"    chain_id: [{TARGET_ID}]",
                  f"    template_id: [{TARGET_ID}]"]
    path.write_text("\n".join(lines) + "\n")


def parse_confidence(pred_root: Path, name: str) -> dict | None:
    hits = sorted(pred_root.rglob(f"confidence_{name}_model_0.json"))
    if not hits:
        return None
    data = json.loads(hits[0].read_text())
    pair = data.get("pair_chains_iptm") or {}

    def pair_val(i, j):
        got = []
        for a, b in ((i, j), (j, i)):
            try:
                got.append(float(pair[str(a)][str(b)]))
            except (KeyError, TypeError, ValueError):
                continue
        return sum(got) / len(got) if got else None

    # binder is chain index 0, target is 1, in YAML order
    bt = pair_val(0, 1)
    if bt is None:
        return None
    return {
        "name": name,
        "binder_target_iptm": bt,          # THE score
        "binder_ptm": (data.get("chains_ptm") or {}).get("0"),
        "target_ptm": (data.get("chains_ptm") or {}).get("1"),
        "global_iptm": data.get("iptm"),   # recorded, never used to rank
        "complex_plddt": data.get("complex_plddt"),
        "complex_iplddt": data.get("complex_iplddt"),
        "ptm": data.get("ptm"),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta-dir", type=Path, nargs="*", default=None)
    ap.add_argument("--out", type=Path, default=ROOT / "boltz")
    ap.add_argument("--tag", default="run")
    ap.add_argument("--target-mode", choices=("full", "apical"), default="full")
    ap.add_argument("--species", choices=("human", "mouse"), default="human")
    ap.add_argument("--per-backbone", type=int, default=None)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--recycling-steps", type=int, default=3)
    ap.add_argument("--diffusion-samples", type=int, default=1)
    ap.add_argument("--sampling-steps", type=int, default=200)
    ap.add_argument("--extra-seqs", type=Path, default=None)
    ap.add_argument("--devices", type=int, default=1)
    ap.add_argument("--override", action="store_true")
    ap.add_argument("--no-template", action="store_true")
    ap.add_argument("--no-composition-filter", action="store_true")
    args = ap.parse_args()

    binders: dict[str, str] = {}
    for d in (args.fasta_dir or []):
        pool = d.parent.name if d.name == "seqs" else d.name
        for name, seq in read_binders(d, args.per_backbone).items():
            binders[f"{pool}_{name}"] = seq

    if args.extra_seqs:
        name, buf = None, []
        for line in args.extra_seqs.read_text().splitlines():
            if line.startswith(">"):
                if name:
                    binders[name] = "".join(buf)
                name, buf = re.sub(r"\W+", "_", line[1:].strip()), []
            elif line.strip():
                buf.append(line.strip())
        if name:
            binders[name] = "".join(buf)

    if not binders:
        print("no sequences to score", file=sys.stderr)
        return 1

    if not args.no_composition_filter:
        kept, rejected = {}, []
        for n, s in binders.items():
            why = composition_reject(s)
            (rejected.append((n, why)) if why else kept.update({n: s}))
        if rejected:
            print(f"composition pre-filter dropped {len(rejected)}/"
                  f"{len(binders)} before folding:")
            for n, why in rejected[:10]:
                print(f"  {n}: {why}")
            (args.out / args.tag).mkdir(parents=True, exist_ok=True)
            (args.out / args.tag / "composition_rejects.json").write_text(
                json.dumps([{"name": n, "reason": w} for n, w in rejected],
                           indent=2))
        binders = kept
    if args.limit:
        binders = dict(list(binders.items())[: args.limit])
    if not binders:
        print("every sequence failed the composition filter", file=sys.stderr)
        return 1

    target, cif = target_sequence(args.target_mode, args.species)
    tmpl = None if args.no_template else cif
    if tmpl and not tmpl.exists():
        print(f"template not found: {tmpl}", file=sys.stderr)
        return 1

    print(f"target : {args.species} / {args.target_mode} / {len(target)} aa")
    print(f"template: {tmpl or 'NONE'}")
    print(f"designs : {len(binders)}")

    work = args.out / args.tag
    yaml_dir = work / "yaml"
    yaml_dir.mkdir(parents=True, exist_ok=True)
    for name, seq in binders.items():
        write_yaml(yaml_dir / f"{name}.yaml", seq, target, tmpl)

    pred_dir = work / "pred"
    # Resolve the boltz binary next to the interpreter running us. Invoking a
    # bare "boltz" only works if the conda env happens to be ACTIVATED, which
    # is true in an interactive shell and false when the autopilot calls this
    # script by absolute interpreter path -- that cost a failed stage once.
    boltz_bin = Path(sys.executable).with_name("boltz")
    if not boltz_bin.exists():
        boltz_bin = Path("boltz")
    cmd = [str(boltz_bin), "predict", str(yaml_dir),
           "--out_dir", str(pred_dir), "--model", "boltz2",
           "--accelerator", "gpu", "--devices", str(args.devices),
           "--recycling_steps", str(args.recycling_steps),
           "--diffusion_samples", str(args.diffusion_samples),
           "--sampling_steps", str(args.sampling_steps),
           "--output_format", "pdb"]
    if args.override:
        cmd.append("--override")
    print("+", " ".join(cmd), flush=True)
    t0 = time.time()
    subprocess.run(cmd, check=True)
    dt = time.time() - t0

    rows = [r for r in (parse_confidence(pred_dir, n) for n in binders) if r]
    rows.sort(key=lambda r: r["binder_target_iptm"], reverse=True)
    meta = dict(species=args.species, target_mode=args.target_mode,
                target_len=len(target), n_designs=len(binders),
                n_scored=len(rows), seconds=round(dt, 1),
                seconds_per_design=round(dt / max(len(binders), 1), 2))
    (work / "scores.json").write_text(json.dumps(
        {"meta": meta, "scores": rows}, indent=2))

    print(f"\n{len(rows)}/{len(binders)} scored in {dt/60:.1f} min "
          f"({meta['seconds_per_design']} s/design)")
    for r in rows[:15]:
        print(f"  {r['name']:<34} binder-target ipTM={r['binder_target_iptm']:.3f} "
              f"binder_ptm={r['binder_ptm']:.3f} "
              f"iplddt={r['complex_iplddt']:.3f}")
    print(f"scores -> {work/'scores.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
