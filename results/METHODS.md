# TfR1 Challenge 04 — Methods

De novo single-chain binders to human transferrin receptor 1, designed for
pH-selective binding (7.4 on, 6.0 off), mouse cross-reactivity, and affinity,
in that ranked order.

**Nothing here is an experimental measurement.** Every number is a
computational estimate from the pipeline below. No design in this submission is
a validated binder until Adaptyv's sensorgrams say so. §9 lists what was not
computed.

---

## 1. Stack and versions

| Tool | Version / commit | Role |
|---|---|---|
| RFdiffusion | `86507b6`, Docker `rosettacommons/rfdiffusion:latest` (`33d504109509`) | backbone generation |
| ProteinMPNN | `8907e66`, soluble weights `v_48_020` | sequence design |
| Boltz-2 | `boltz2_conf.ckpt` | complex prediction and scoring |
| PROPKA | 3.x (`propka` via venv) | pKa shifts for the pH linkage |
| Biopython | 1.88 | structure handling, SASA, alignment |

Hardware: one A100-80GB, shared. Campaign wall-clock 45.5 h.

## 2. Target preparation

Target is **PDB 11QC chain A**, the human TfR1 ectodomain monomer.

- **Monomer, not dimer.** 11RN is A2B2 with one B19V VP1u RBD bound per TfR1
  subunit, so the apical epitope is a single-subunit epitope. Designing and
  scoring against the dimer would have roughly doubled the denoised length for
  no epitope that isn't already present.
- **Sequence verified against the assay construct**, not UniProt. The
  challenge page's 640-residue human construct (121–760) matches 11QC chain A
  at **all 639 modelled residues, zero mismatches**.
- **Gly142 confirmed present.** The Ser142 allele carried by 1CX8/3KAS/3S9L/6WRV
  would have been a silent corruption; no patch was needed.
- **Glycans stripped** (6 NAG) and the three sequons **N251, N317, N727** plus
  everything within 10 Å (73 residues) masked as non-designable. A binder
  parked on a glycan fails in the assay for a reason no score reveals.

The **mouse** construct (641 aa, 123–763) was recovered from the challenge
page's raw HTML — the Mouse tab is JS-gated in rendered text but both sequences
are server-rendered in the payload. **UniProt Q62351 was not used.** The
sequence was cross-validated against all six [OFFICIAL] page-stated apical
differences, which reproduce exactly (D204→Q207, R208→N211, V210→D213,
Y211→P214, N348→K350, K371→Q373), at 78.1% global identity against the page's
stated 78%.

## 3. Epitope selection — measured, not assumed

The epitope was chosen by a conservation map (`00_conservation.py`,
`00b_patch_scan.py`), not from the literature. Key results:

**K205 is not conserved** (K→S208) **and is not in the epitope.** The VP1u
footprint on 11RN is exactly nine residues — 202, 204, 208, 209, 210, 211, 212,
348, 371 — and the two the page does not name are **I202 and L209**, both
conserved and both hydrophobic. There is **no conserved cation in the epitope
core**.

This matters because the pH switch mechanism needs a target cation: a
protonated interface His is repelled by Arg/Lys at pH 6.0. The apical cations
R208 and K371 are both neutral in mouse, so a switch built on them cannot be
cross-reactive.

**K374 is the one usable conserved anchor.** Conserved (K→K376), exposed
(rel. SASA 0.282), **no carboxylate within 10 Å**, 6.7–7.9 Å from R208/K371/Y211,
40.6 Å from the transferrin site, 15.2 Å from the nearest glycan sequon.

Two cations were explicitly **vetoed**:
- **K344** — conserved and exposed, but salt-bridged to **E214 at 2.8 Å**. A
  protonated His beside it is *stabilised* at pH 6.0, so the design would bind
  **tighter** when it should release. A conservation-only screen would have
  selected it.
- **K231** — conserved, but 8.7 Å from the N251 sequon, inside the glycan shadow.

The unverified literature lead at **H318 / K633 / R719** was **discarded**:
H318 sits 1.3 Å from the N317 sequon, K633 is inside the transferrin site
(619–670) and not conserved (K→R636), R719 abuts it.

**Chosen epitope:** the apical shoulder centred on S370/K374, overlapping 5 of
the 9 VP1u contacts (functional-epitope credit, and 11RN provides a real
positive control) and 33 Å clear of the transferrin site, so receptor recycling
is untouched.

## 4. Generation

| round | hotspots | backbones | lengths |
|---|---|---|---|
| 1 | A209, A211, A348, A371, A374 (shoulder) | 340 | 55–120 |
| 1 | A204, A208, A211, A212 (apex) | 135 | 70–120 |
| 1 | A211, A348, A371 (swap locus) | 65 | 70–95 |
| 2 | **A348, A371, A374** (K374-centric) | 170 | 70–120 |

Generation used an **18 Å crop** (145 target residues). At 15 Å the shell kept
only 109 residues — the apical domain protrudes, so there is less surrounding
protein than TNF had — which risked binders clashing with cropped-away atoms.
The crop is for generation only; **all scoring used the uncropped 639-residue
monomer**.

Sequences: ProteinMPNN soluble weights, `--omit_AAs CX`, three pools — plain
(control), global His bias (1.5), and a **targeted per-position His bias** (2.5)
applied only at binder positions within 10 Å of a conserved cation.

## 5. Filters, with provenance

| stage | filter | threshold | provenance |
|---|---|---|---|
| 35 | composition | no free Cys, no N-X-S/T, no homopolymer ≥5 | §7 of the brief |
| 04 | Boltz-2 binder–target pair ipTM | ranked, not cut | monomer target → one chain pair |
| **04e** | **epitope fidelity** | ≥2 hotspots, **0** transferrin contacts, ≤2 glycan-shadow contacts, C-term ≥8 Å | **this campaign** |
| 05 | PROPKA linkage ddG | sites with ddG > 0.3, cation ≤6 Å, no anion ≤6 Å, His–His >8 Å | comp-02 |
| **05b** | **conservation policy on switch sites** | vetoed-cation sites discarded; scope labelled | **this campaign** |
| 06 | His refold | variant ipTM ≥0.35 and (≥0.40 or ≤25% drop) | comp-02 |
| 06b | mouse screen | second-stage, on human survivors only | §5 of the brief |
| 06c | two-seed pose agreement | Jaccard ≥0.50, centroid shift ≤5 Å **after superposing the targets** | comp-03 |

### 5a. Why epitope fidelity is a hard gate

Scoring against a 639-residue monomer lets Boltz dock the binder anywhere, and
it does. **172 of 900 round-1 designs (19%) land on the targeted epitope.** The
distribution is bimodal — 374 of 494 sampled hit *zero* targeted hotspots —
so designs either grip the patch or are somewhere else entirely.

Confidence does not catch this. The **highest-scoring design of all 900**
(ipTM 0.878) fails the gate, and a pilot design scoring **0.775** bound the
transferrin site and a glycan sequon with 0/5 hotspots. Ranking on ipTM alone
would have put an off-epitope binder at the top of the slate.

ipTM and fidelity do correlate (10% pass below 0.3; 68% pass at 0.5–0.6) but 44
designs scoring ≥0.5 still dock off-target.

### 5b. Why the pH scan needed a conservation policy

`05_ph_scan.py` is target-agnostic: it finds any interface position whose His
substitution is destabilised by a nearby cation. Run bare it produced 497 sites,
of which **two thirds were unusable**: 236 on R208 and 78 on K371 (neither
conserved), and **96 on K344** (the E214 salt-bridge trap). Only 77 sat on K374.

`05b` discards vetoed-cation sites and labels the rest `cross_species` or
`species_specific`, so objective 2 is never claimed for a human-only switch.

### 5c. Scoring target — full monomer, not a crop

Timed on the pilot: apical crop (195 aa) 27.3 s/design, full monomer (639 aa)
39.5 s/design. Only 1.4× slower, and the two **disagree**. The full monomer was
chosen because it *exposes* off-target docking that a crop hides — a crop
cannot show a binder reaching the transferrin site because the transferrin site
is not in it.

## 6. Results

| | |
|---|---|
| backbones generated | 710 (540 round 1 + 170 round 2) |
| MPNN sequences | ~6,300 across three pools |
| Boltz-scored | 1,372 (900 round 1 + 472 round 2) |
| pass epitope gate | 254 (172 + 82) |
| designs with a usable switch | 221 |
| **survive the pH 7.4 His refold** | **57** |
| of those, cross-species switch | 4 |
| **agree across two seeds (pose gate)** | **49 / 57** |
| mouse-scored | 62 (57 survivors + 5 controls) |
| **submitted** | **20** |

Final slate: 6 cross-reactive (class A), 7 human pH switches (B), 2 affinity
(C), 5 unmutated controls (D). **20 distinct backbones**, 61-118 aa.

## 7. Honest statement on each objective

**Objective 1 — pH selectivity.** Achievable but **modest**. Best total linkage
ddG is **2.21 kcal/mol ≈ 42× KD shift**, and that design's switch is built on
R208/K371 and so is human-only. The 2.8 kcal/mol (~100×) target from challenge
02 is **unreachable on this target**: maximum single-site ddG measured is 1.29
kcal/mol and no design has three independent strong sites at >8 Å separation.
Results are reported in tiers (strong ≥2.0, moderate ≥1.4, weak ≥1.0) with the
fold-shift stated, not as pass/fail. This sits at the weak end of the published
strategy-1 range (122× on TNFR2, up to ~1000× on EphA2).

**Objective 2 — mouse cross-reactivity.** Better than the mechanism work
suggested, for a reason worth stating plainly.

Measured across the 57 refold survivors scored against the mouse target:
median mouse:human ipTM ratio **0.33**, but **8 designs at ratio ≥0.75** and
**5 at ≥0.90**, the best at **1.34** (`hisb_shoulder_m_d_67_s4__E72H`, human
0.451 / mouse 0.606 — it binds mouse *better* than human). Six survive the pose
gate and are submitted as class A.

**A conserved-cation switch does NOT imply a cross-reactive binder.** These are
independent properties and we initially conflated them. Of the four surviving
cross-species-switch designs the mouse:human ratios are 1.34, 0.46, 0.20 and
0.19 — and the best human binder of the campaign (ipTM 0.755) sits at 0.19.
Meanwhile several designs whose switch is *human-only* cross-react well (1.34,
1.33, 1.07). In hindsight this is obvious: the switch is **one residue**, the
interface is 15–20, and most of those differ between species. Objective 2 is
therefore ranked on **measured mouse retention**, not on our mechanism label.

So cross-reactive *binding* is achievable here; what is hard is cross-reactive
*switching*. On that, the cause was measured rather than guessed: across the
123 on-epitope backbones the closest binder atom sits at median **5.1 Å from
Y211/L209 but 11.6 Å from K374**, the only clean conserved cation. RFdiffusion
centres binders on the protruding apical tip regardless of where it is pointed
— round 2 removed 209/211 from the hotspot set entirely and 169 of 170 designs
still gripped them at 5.1 Å. A histidine cannot bridge 11.6 Å.

The targeted His bias did work at its specific job — cross-species switches by
pool were plain 0, global-His 3, **targeted 5** — but K374-centric designs bind
measurably worse (top-five ipTM 0.148–0.270 vs 0.813 for the best round-1
design) and survive the refold at 10% vs 33%. **The conserved anchor sits where
binding is hard.** That trade-off is the central finding of this campaign.

**No design in this slate is claimed to be both strongly pH-selective and
cross-reactive.** The class-A designs cross-react but switch weakly (0.16–1.39
kcal/mol); the class-B designs switch up to 16× but are human-only.

**Objective 3 — human affinity.** **The strongest result.** 61 on-epitope
designs at ipTM ≥0.5, best 0.813, at a functional epitope overlapping 5 of the 9
VP1u contacts and 33 Å clear of the transferrin site, so transferrin binding and
receptor recycling should be unaffected.

## 8. Deviations and corrections

- **Scoring target** changed from crop to full monomer after the pilot timing
  (§5c). The original rationale recorded for this was **wrong** and is corrected
  in `runlog/00-schedule.md`: the crop scored the off-target design low because
  it lacks the site that design binds, not because it mis-ranks good designs.
- **swap_locus arm abandoned** after round 1 (12/115 = 10% epitope fidelity);
  three tightly clustered hotspots give a binder too little to grip.
- **No new generative model** was introduced. The secondary-structure census
  showed round-1 backbones are **93.3% helix / 1.6% strand**, i.e. one fold
  class. Fold-conditioned generation was available in-engine
  (`scaffoldguided` + `Complex_Fold_base_ckpt.pt`) but was **not run** — time
  went to the K374 problem instead. This is a real diversity limitation of the
  slate and is stated rather than hidden.
- **The pose gate was initially broken and is worth recording.** Its first
  version compared binder centroids between the two seeds without superposing
  the targets, and reported 0/57 agreeing with 55-87 A shifts for designs whose
  contact sets were *identical* (Jaccard 1.0). Boltz writes each prediction in
  its own global frame; template pinning constrains conformation, not
  placement. After superposing on target CA atoms: **49/57 agree**, 3 flip to a
  genuinely different site. One of those three had already been selected into
  class A on a mouse ratio of 1.14 and was removed.
- **Four conditions / three conditions.** The challenge page's ranking intro
  says "four conditions" then names three. Every ranking criterion is scoped to
  a named condition, so no stated criterion consumes a fourth. Proceeded on
  three. Unresolved at submission time.

## 9. Not computed

- **No experimental validation of anything.**
- **No pH was modelled by any structure predictor.** Boltz-2 has no pH
  parameter. pH selectivity is thermodynamic linkage from PROPKA pKa shifts —
  a mechanism hypothesis, **not a predicted KD ratio**.
- **No mouse structure exists.** Every mouse number comes from threading the
  mouse sequence onto the *human* structure. Whether mouse K350 projects into
  the interface where human has the buried N348 is unmodelled and unknowable
  from this data.
- Additivity across His sites assumes electrostatic independence; pairs closer
  than 8 Å were rejected rather than scored, because a coupled pKa calculation
  (constant-pH MD or Poisson–Boltzmann) did not fit the time budget.
- No expression, aggregation, immunogenicity or stability prediction.
- Boltz ran in single-sequence mode (no MSA); binders are de novo and the
  target MSA is identical across designs.

## 10. De novo compliance

No design contains a contiguous run matching VP1u (11RN), Machupo GP1 (3KAS),
the 6WRV minibinders, transferrin, or the published brain-shuttle modules.
All sequences originate from RFdiffusion backbones with ProteinMPNN sequences;
no natural scaffold was grafted.

`molecule_class` is `single_chain` throughout. A de novo nanobody or scFv would
have required RFantibody, which was out of scope (§6 of the brief) and would
not have been de novo if built on a real VHH framework.
