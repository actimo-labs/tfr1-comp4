# TfR1 Challenge 04 — conditional, cross-species binder

Read `tfr1-handoff-brief.md` in full before acting. It is the complete context:
challenge rules, the epitope decision gate, filter stack, portfolio and working
agreement.

Deadline: **Mon 2 Nov 22:59 Melbourne (AEDT)** = Sun 1 Nov 23:59 AoE. Hard.

## What this is

Anthropic × Adaptyv Challenge 04. Design 20 de novo single-chain binders
(10–250 aa) that:

1. bind human TfR1 at pH 7.4 and show **no detectable binding at pH 6.0**
2. bind **mouse** TfR1 at pH 7.4 as well as human
3. bind human TfR1 as tightly as possible

Ranked in that order. Unlike challenge 02, cross-reactivity is **required**,
not dropped.

## Environment

- One A100-80GB, budget-constrained. Do not launch speculative GPU jobs.
- Verify the layout before assuming it. First action is the audit in §0 of the
  brief — `ls`, `nvidia-smi`, `git rev-parse`, record what is actually on disk.
- Reuse the challenge-02 code (`tnf-upload/scripts/`). It is a source of
  working *code*, not of *strategy* — the epitope logic inverts here.
- RFdiffusion runs via the Docker wrapper. Hydra overrides only (the image has
  an ENTRYPOINT), container-absolute paths under `/work`, include
  `hydra.run.dir=/work/.hydra`.

## Non-negotiables

- **The epitope is not chosen yet.** Stage 00 computes the human/mouse surface
  conservation map and the decision follows from it. Do not hardcode apical
  apex hotspots from the literature and start generating.
- **`target/human_tfr1_121_760.fasta` is the target sequence.** It is from the
  challenge page, 640 aa, span 121–760, position 142 = **Gly**. Verified. Do not
  substitute UniProt. The mouse sequence must be pulled from the page's Mouse
  tab in a browser — it is JS-gated and did not survive a plain fetch.
- **Verify `git config user.name`/`user.email` before the first commit.** They
  must be Harjas's, with an email verified on his GitHub. Report a mismatch;
  do not set them yourself.
- **pH selectivity cannot be scored by a folding model.** AF2, Boltz and Chai
  ignore protonation. It is scored by thermodynamic linkage on CPU in
  `05_ph_scan.py`. Do not substitute a confidence score for it.
- **A switch histidine needs a target Arg/Lys, never a target Asp/Glu.** A
  protonated His beside a carboxylate strengthens binding at low pH — backwards.
- **The apical cations R208 and K371 are not conserved in mouse** (→ N211,
  Q373). An interface-repulsion switch built on them cannot be cross-reactive.
  See §3 of the brief — this is the central conflict of this challenge.
- **Crop before generating.** The TfR1 ectodomain is 639 aa per chain. Uncropped
  generation does not fit the deadline by roughly 3×. See §5.
- **Design against a monomer, not the dimer**, unless the conservation map sends
  the epitope to the dimer seam. 11RN is A2B2 with one VP1u per subunit — the
  apical epitope is a single-subunit epitope.
- **Every His variant must be re-folded** to confirm the pH 7.4 interface
  survives. A design that switches but no longer binds at 7.4 is worthless.
- **Mouse scoring is a second-stage filter on human survivors**, not a parallel
  arm. Scoring both species on everything doubles the Boltz bill for nothing.
- **Do not introduce a new generative model.** RFdiffusion → ProteinMPNN →
  Boltz-2 is measured and debugged. See §6.

## Working style

- Commit after every meaningful unit. Commits are pre-authorised and local.
- A failed commit is logged and worked around, never a reason to stall.
- **Never push to remote. Never submit to the competition site.**
- Never modify or delete challenge-02 or challenge-03 artifacts.
- Long GPU jobs run in tmux. Never in a bare foreground shell.
- Check in before anything costing more than ~30 GPU-minutes.
- Say when a number comes from a run versus an estimate.
- Flag it when the brief stops matching reality rather than improvising around
  it quietly.