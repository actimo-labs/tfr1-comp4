#!/usr/bin/env python3
"""
Stage 35 - choose which designs are worth a Boltz call.

Why this stage exists. The pilot measured Boltz at 39.5 s/design on the full
639-residue monomer, so the screen has room for roughly 900 designs in the
time available. MPNN produced 4 sequences per backbone per pool across 540
backbones = 4320 candidates. Something has to choose, and the choice should
protect the two things the portfolio actually needs:

  1. BACKBONE DIVERSITY. §7 of the brief makes structural diversity a hard
     requirement -- challenge 02 shipped 6 clusters across 20 designs and
     called it its weakest axis. So the rule is best-N-per-backbone, never
     global top-N: a global cut would hand every slot to whichever handful of
     backbones MPNN happens to score well and collapse the slate onto one fold.

  2. EPITOPE BALANCE. Quotas are proportional to each batch's backbone count,
     so the apex and swap-locus arms cannot be squeezed out by the larger
     shoulder batches. The portfolio classes have to survive this stage.

MPNN `score` is used only to pick between sequences OF THE SAME BACKBONE, where
it is a reasonable discriminator. It is deliberately NOT used to rank across
backbones -- it is a sequence-likelihood measure, not a binding predictor, and
treating it as one is how you get twenty variations of a single bet.

Writes gen/selected.fasta and gen/selected.tsv.
"""
import argparse
import json
import re
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

SEQUON = re.compile(r"N[^P][ST]")
HOMOPOLYMER = re.compile(r"(.)\1{4,}")


def composition_reject(seq: str) -> str | None:
    if "C" in seq:
        return "free Cys"
    m = SEQUON.search(seq)
    if m:
        return f"N-X-S/T sequon@{m.start()}"
    m = HOMOPOLYMER.search(seq)
    if m:
        return f"homopolymer {m.group(1)}x{len(m.group(0))}"
    return None


def read_pool(seqs_dir: Path, pool: str):
    """Yield (backbone, batch, pool, sample, score, seq) for every design."""
    for fa in sorted(seqs_dir.glob("*.fa")):
        backbone = fa.stem                      # e.g. shoulder_m_d_17
        batch = backbone.rsplit("_d_", 1)[0]    # e.g. shoulder_m
        name, buf, recs = None, [], []
        for line in fa.read_text().splitlines():
            if line.startswith(">"):
                if name is not None:
                    recs.append((name, "".join(buf)))
                name, buf = line[1:], []
            elif line.strip():
                buf.append(line.strip())
        if name is not None:
            recs.append((name, "".join(buf)))
        for hdr, seq in recs[1:]:               # record 0 is the backbone
            s = re.search(r"score=([\d.]+)", hdr)
            n = re.search(r"sample=(\d+)", hdr)
            yield dict(backbone=backbone, batch=batch, pool=pool,
                       sample=(n.group(1) if n else "0"),
                       score=float(s.group(1)) if s else 99.0,
                       seq=seq.split("/")[0])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mpnn", type=Path, default=ROOT / "gen/mpnn")
    ap.add_argument("--pools", default="plain,hisb")
    ap.add_argument("--cap", type=int, default=900)
    ap.add_argument("--per-backbone", type=int, default=1,
                    help="best N sequences kept per backbone per pool")
    ap.add_argument("--out-fasta", type=Path, default=ROOT / "gen/selected.fasta")
    ap.add_argument("--out-tsv", type=Path, default=ROOT / "gen/selected.tsv")
    args = ap.parse_args()

    allrecs = []
    for pool in [p.strip() for p in args.pools.split(",") if p.strip()]:
        d = args.mpnn / pool / "seqs"
        if not d.exists():
            print(f"missing pool dir {d}")
            continue
        allrecs += list(read_pool(d, pool))
    print(f"[35] {len(allrecs)} raw MPNN designs")

    kept, rejected = [], []
    for r in allrecs:
        why = composition_reject(r["seq"])
        (rejected.append((r, why)) if why else kept.append(r))
    print(f"[35] composition filter dropped {len(rejected)}, kept {len(kept)}")
    rc = defaultdict(int)
    for _, why in rejected:
        rc[why.split("@")[0].split(" x")[0]] += 1
    for k, v in sorted(rc.items(), key=lambda x: -x[1])[:6]:
        print(f"       {k}: {v}")

    # best-N per (backbone, pool)
    by_bb = defaultdict(list)
    for r in kept:
        by_bb[(r["backbone"], r["pool"])].append(r)
    shortlist = []
    for k, rs in by_bb.items():
        rs.sort(key=lambda r: r["score"])
        shortlist += rs[: args.per_backbone]
    print(f"[35] best-{args.per_backbone}-per-backbone-per-pool "
          f"-> {len(shortlist)} designs")

    # proportional per-batch quotas so no arm is squeezed out
    per_batch = defaultdict(list)
    for r in shortlist:
        per_batch[r["batch"]].append(r)
    total = len(shortlist)
    selected = []
    if total <= args.cap:
        selected = shortlist
        print(f"[35] under cap ({args.cap}) - keeping all")
    else:
        for batch, rs in sorted(per_batch.items()):
            quota = max(1, round(args.cap * len(rs) / total))
            rs.sort(key=lambda r: r["score"])
            selected += rs[:quota]
        selected = selected[: args.cap]
        print(f"[35] capped to {len(selected)} with proportional batch quotas")

    counts = defaultdict(int)
    for r in selected:
        counts[(r["batch"], r["pool"])] += 1
    print("\n[35] selected by batch / pool:")
    for (b, p), n in sorted(counts.items()):
        print(f"       {b:<12} {p:<6} {n}")

    args.out_fasta.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out_fasta, "w") as f:
        for r in selected:
            f.write(f">{r['pool']}_{r['backbone']}_s{r['sample']}\n{r['seq']}\n")
    with open(args.out_tsv, "w") as f:
        f.write("name\tbackbone\tbatch\tpool\tsample\tmpnn_score\tlength\n")
        for r in selected:
            f.write(f"{r['pool']}_{r['backbone']}_s{r['sample']}\t{r['backbone']}\t"
                    f"{r['batch']}\t{r['pool']}\t{r['sample']}\t{r['score']:.4f}\t"
                    f"{len(r['seq'])}\n")
    json.dump([{"name": f"{r['pool']}_{r['backbone']}_s{r['sample']}",
                "reason": w} for r, w in rejected],
              open(args.out_tsv.with_name("composition_rejects.json"), "w"),
              indent=1)
    print(f"\n[out] {args.out_fasta}  ({len(selected)} designs)")
    print(f"[out] {args.out_tsv}")


if __name__ == "__main__":
    main()
