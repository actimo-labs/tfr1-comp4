# TfR1 Challenge 04 — Anthropic × Adaptyv

De novo single-chain binders to **human transferrin receptor 1** (TfR1), designed for:

1. **pH-selective binding** — bind at pH 7.4, release at pH 6.0
2. **Mouse cross-reactivity** — bind mouse TfR1 as well as human
3. **Affinity** — lowest KD against human TfR1 at pH 7.4

Ranked in that order. 20 designs submitted, 61–118 aa, all `single_chain`.

**Status: complete.** Submission package at `upload/submission/`.

---

## The submission

| | |
|---|---|
| `upload/submission/results/submission.csv` | the 20 designs: `name, sequence, molecule_class` |
| `upload/submission/results/METHODS.md` | full methods, filters with provenance, deviations |
| `upload/submission/results/methods_summary.md` | ~640-word summary |
| `upload/submission/results/rationale.md` | per-design rationale and caveats |
| `upload/submission/structures/` | 20 predicted complexes |
| `upload/submission/SHA256SUMS.txt` | verify with `sha256sum -c` |

72 files, 10.2 MB. **Never pushed to remote, never submitted** — both require explicit approval.

## Pipeline

```
RFdiffusion → ProteinMPNN → Boltz-2        (structure + affinity)
PROPKA thermodynamic linkage               (pH selectivity — separate, on CPU)
```

```
710 backbones
  → ~6,300 MPNN sequences (3 pools)
  → 1,372 Boltz-scored against the FULL 639-residue monomer
  →   254 pass the epitope-fidelity gate        (19%)
  →   221 carry a usable pH switch
  →    57 survive the pH 7.4 His refold
  →    49 agree across two random seeds
  →    20 submitted
```

Stage scripts are in `scripts/`, numbered in execution order (`00_` → `81_`).

## The epitope decision

Chosen by measurement, not from the literature — see `target/conservation_report.md`.

The pH switch needs a target Arg/Lys beside an interface histidine: protonation
at pH 6.0 creates charge repulsion and the binder releases. The problem is that
the apical cations **R208 and K371 are both neutral in mouse**, so a switch
built on them cannot be cross-reactive. K205, the candidate the brief hoped for,
is a **serine** in mouse and is not in the epitope at all.

**K374 is the one usable conserved anchor** (K→K376, exposed, no carboxylate
within 10 Å). Two apparent alternatives were vetoed by the map:

- **K344** — conserved and exposed, but salt-bridged to **E214 at 2.8 Å**. A
  protonated His there is *stabilised* at pH 6.0, so the design would bind
  **tighter** when it should release. A conservation-only screen picks this.
- **K231** — conserved, but 8.7 Å from the N251 glycan sequon.

Final epitope: the apical shoulder around **S370/K374**, overlapping 5 of the 9
VP1u contacts (functional epitope, with 11RN as a real positive control) and
33 Å clear of the transferrin site, so receptor recycling is untouched.

## Why there are four hard gates beyond confidence

**Interface confidence does not tell you where a binder binds.** Scoring against
a 639-residue monomer lets Boltz dock the binder anywhere, and it does — only
19% land on target. The **highest-scoring design of all 900** (ipTM 0.878)
contacts the transferrin site and a glycan sequon with **zero** targeted
hotspots. Ranking on confidence alone would have put it at the top of the slate.

| gate | what it caught |
|---|---|
| `04e` epitope fidelity | off-target docking, transferrin-site and glycan binders, buried C-termini |
| `05b` conservation policy | two thirds of raw switch sites sat on cations already ruled out, 96 on the K344 trap |
| `06` His refold | 164 of 221 variants no longer bind at pH 7.4 after His substitution |
| `06c` two-seed pose agreement | 8 designs whose pose is seed-dependent — one was already in the slate |

## Results

| class | n | role |
|---|---|---|
| **A** | 6 | mouse cross-reactive — mouse:human ipTM ratios 1.34, 1.33, 1.07, 0.90, 0.86, 0.75 |
| **B** | 7 | human pH switches — up to 1.65 kcal/mol (**16×** KD shift) |
| **C** | 2 | affinity |
| **D** | 5 | unmutated controls — the internal comparison for the switch hypothesis, to ipTM **0.813** |

**20 distinct backbones**, worst pairwise identity 29% against a 70% ceiling.
Challenge 02 shipped 6 structural clusters across 20 designs and called that its
weakest axis.

## Honest summary

**Objective 3 (affinity) is the strongest result** — 61 on-epitope designs at
ipTM ≥0.5, best 0.813, at a functional epitope.

**Objective 1 (pH selectivity) is modest** — best 16×, and that switch is
human-only. The ~100× bar was unreachable here: maximum single-site ddG measured
is 1.29 kcal/mol and no design has three independent strong sites.

**Objective 2 (cross-reactivity) is better than the mechanism work suggested,
for an instructive reason.** A conserved-cation switch does **not** imply a
cross-reactive binder — the switch is one residue, the interface is 15–20, and
most of those differ between species. Cross-reactive *binding* is achievable
(8 designs at ratio ≥0.75); cross-reactive *switching* is what is hard, because
K374 sits at median **11.6 Å** from the nearest binder atom while Y211/L209 sit
at 5.1 Å. RFdiffusion centres binders on the protruding apical tip regardless of
where it is pointed — round 2 removed 209/211 from the hotspot set entirely and
169 of 170 designs still gripped them.

**No design is claimed to be both strongly pH-selective and cross-reactive.**
Class A cross-reacts but switches weakly; class B switches well but is human-only.

**Nothing here is an experimental measurement.** No structure-prediction model
used has a pH parameter — pH selectivity is a mechanism hypothesis from
thermodynamic linkage, not a predicted KD ratio. No mouse structure exists; all
mouse numbers come from threading the mouse sequence onto the human structure.

## Layout

```
scripts/     pipeline, stages 00 → 81
target/      assay-construct sequences (human, mouse, VP1u control),
             prepared target, crops, hotspot sets, glycan mask,
             conservation_report.md  ← the epitope decision
results/     submission.csv, METHODS.md, rationale, metrics, rejected log
structures/  predicted complexes for the slate
runlog/      environment audit, schedule (incl. a CORRECTION section), open items
gen/         backbones, MPNN output, selection and QC tables
boltz/       predictions and scores (gitignored — large, regenerable)
upload/      staged packages: batch1, batch2 (interim), submission (final)
```

Reproduce the epitope decision in a few minutes of CPU:

```bash
venv-phscan/bin/python -I scripts/00_conservation.py
venv-phscan/bin/python -I scripts/00b_patch_scan.py
```

## Environment

One A100-80GB, 45.5 h wall-clock. RFdiffusion via Docker (Hydra overrides only,
container-absolute paths under `/work`). ProteinMPNN and Boltz-2 in conda envs
`mpnn` / `boltz`; PROPKA and Biopython in `venv-phscan`. See
`runlog/00-environment-audit.md` for the verified layout.
