# Stage 00 — TfR1 human/mouse surface conservation map

**Decision gate. Generated Sat 31 Oct 2026, 14:30–14:40 UTC (Sun 1 Nov 01:30 AEDT).**

Inputs: `target/human_tfr1_121_760.fasta` (640 aa) and
`target/mouse_tfr1_123_763.fasta` (641 aa) — **both the assay constructs from
the challenge page**, plus PDB **11QC** (human TfR1 ectodomain, 2.40 Å) and
**11RN** (B19V VP1u RBD + TfR1).

Reproduce with:

```bash
venv-phscan/bin/python -I scripts/00_conservation.py   # → conservation_table.tsv
venv-phscan/bin/python -I scripts/00b_patch_scan.py    # → patch_scan.json
```

---

## Headline

| Question | Answer |
|---|---|
| **Is K205 conserved in mouse?** | **No. K205 → Ser208.** And K205 is **not in the epitope at all.** |
| Does the §3 conflict dissolve? | **The switch half does. The affinity half does not.** |
| Is there a conserved cation usable for a strategy-1 switch? | **Yes — K374 (→K376), on the rim of the apical epitope, with no carboxylate within 10 Å.** |
| Recommended epitope | **Apical shoulder, centred on the S370/K374 locus**, overlapping 5 of the 9 VP1u footprint residues |
| Which §4 decision rule fires | **Rule 1, in modified form** — not via K205, via K374 |
| The 318/633/719 seam lead | **Discarded.** See Q3. |

The brief's best case was "K205 is conserved, strategy 1 works cross-species at
the recommended epitope, §3 dissolves." That specific case is **dead**. A
different and almost-as-good one is **alive**: the conserved cation is K374
rather than K205, it sits on the rim rather than the core of the apical
epitope, and it is cleaner than K205 would have been because nothing anionic is
near it.

---

## 1. Alignment, and why it can be trusted

Global pairwise alignment (BLOSUM62, gap −11/−1) of the two **page** constructs:

- **78.1 % identity** over 639 aligned columns (499/639). The page states 78 %. ✔
- Human 121–760 ↔ mouse 123–763, one single-residue deletion in the apical loop.

**The mouse sequence was recovered from the challenge page, not UniProt.** The
Mouse tab is JS-gated in rendered text, but the raw HTML payload carries both
sequences server-rendered; the 641-residue mouse string was extracted directly
and matches the stated 123–763 span exactly.

It is then **cross-validated against the six [OFFICIAL] human→mouse
differences the page states** — all six reproduce, at the exact mouse
positions:

| Page states | Computed | |
|---|---|---|
| D204 → Q207 | D204 → Q207 | ✔ |
| R208 → N211 | R208 → N211 | ✔ |
| V210 → D213 | V210 → D213 | ✔ |
| Y211 → P214 | Y211 → P214 | ✔ |
| N348 → K350 | N348 → K350 | ✔ |
| K371 → Q373 | K371 → Q373 | ✔ |

Six independent page-stated facts reproducing from an independently-obtained
sequence is strong evidence the sequence *and* the numbering map are right.

*One honest caveat.* A plain aligner puts the apical gap at **V210**; the page
puts it at **L212**. Both are equal-cost alignments of the same single deletion
— only the gap column differs. The page is [OFFICIAL], so its placement was
adopted (`anchor_to_page()` in `00_conservation.py`, which logs the override).
Under either placement all of 210/211/212 are non-conserved, so nothing
downstream turns on it.

### 11QC is clean

- 639 residues modelled in chain A, span 122–760; chains A + B; **6 NAG** (3 per chain) ✔
- **Position 142 = Gly.** The Ser142 trap does **not** apply. No patch needed.
- **Zero mismatches** between 11QC chain A and the assay construct across all
  639 modelled residues. The structure is the construct.

---

## 2. Q1 — K205

```
human K205  →  mouse S208        NOT CONSERVED
```

Worse for the hypothesis: **K205 is not in the epitope.** Computing the VP1u
footprint on 11RN (TfR1 chain A heavy atoms within 4.5 Å of any VP1u atom)
returns **exactly nine residues** — the number the page implies:

```
202, 204, 208, 209, 210, 211, 212, 348, 371
```

The page names seven (204, 208, 210, 211, 212, 348, 371). **The two unnamed
residues are I202 and L209** — both conserved, both hydrophobic, and I202 is
largely buried (rel. SASA 0.025) while L209 is well exposed (0.579).

So the nine-residue patch resolves as: **six divergent, three conserved, and
the three conserved ones are I202, L209 and nothing else.** There is **no
conserved cation inside the epitope core.** K205 is not a member of it.

### The apical cation inventory

| Human | rel.SASA | Mouse | Charge change |
|---|---|---|---|
| K205 | 0.563 | S208 | cation → polar |
| R208 | 0.553 | N211 | cation → polar |
| K371 | 0.268 | Q373 | cation → polar |
| N348 | 0.122 | **K350** | polar → **cation** |

Every human apical cation is neutralised in mouse, and mouse gains one where
human has none. §3 of the brief is confirmed exactly as written — *for the
epitope core*.

---

## 3. The cation is rearranged, not deleted — and how close the near-miss was

The brief hypothesised that because mouse gains K350 where human has N348, a
binder His sited between human R208 and position 348 might find a cation in
both species. **Measured on 11QC:**

| Pair | min heavy-atom | Cβ–Cβ | tip–tip |
|---|---|---|---|
| **371 ↔ 348** | **3.6 Å** | 5.2 Å | 4.6 Å |
| 208 ↔ 371 | 8.4 Å | 11.1 Å | 14.1 Å |
| 208 ↔ 348 | 10.4 Å | 14.5 Å | 17.2 Å |
| 205 ↔ 348 | 11.0 Å | 12.7 Å | 13.3 Å |

**Positions 348 and 371 are 3.6 Å apart and both are VP1u contacts.** Human
carries its cation at 371; mouse carries one at 350 (= human 348). The cation
does not leave the epitope — it **moves about 5 Å within it**.

This is tantalising and I want to be careful not to oversell it. A His placed
at that locus would be ~4–6 Å from K371 in human and ~4–6 Å from K350 in mouse,
which is the geometry strategy 1 wants, in both species, from one design.

**But it rests on an inference I cannot measure:** human N348 has rel. SASA
0.122 — it is substantially buried. Whether a *lysine* substituted at 348 in
mouse projects its charged tip out into the interface, or tucks along the
surface, is a property of a structure that does not exist. Lys is longer than
Asn, so it plausibly reaches further, but that is a model-level guess and the
mouse structure will be a homology model, not a measurement.

**So: a real opportunity, worth a minority of design slots, and not the
foundation of the campaign.** The foundation is the next section, which does
not require guessing about a side chain nobody has seen.

---

## 4. Q2 — conserved surface cations, and the switch anchors

Shrake–Rupley SASA on the isolated monomer, surface = rel. SASA ≥ 0.25:

- **244 surface residues**, of which **163 conserved (67 %)**. Note the surface
  is markedly *less* conserved than the 78 % global figure — the divergence is
  concentrated where binders bind.
- **32 conserved surface Arg/Lys**; **23** survive the Tf-sparing (>10 Å from
  619–670), glycan-clear (>10 Å from N251/N317/N727) and
  not-dimer-seam-buried filters.

### Switch-anchor viability

A switch cation is only usable if it is **conserved** *and* has **no target
carboxylate within ~6 Å** — a nearby Asp/Glu inverts the mechanism, stabilising
the protonated His at pH 6.0 and strengthening binding exactly when we want it
to fail.

| Cation | Mouse | Conserved | Nearest carboxylate | Verdict |
|---|---|---|---|---|
| **K374** | K376 | ✔ | **none within 10 Å** | **USABLE — best available** |
| K231 | K233 | ✔ | E369 @ 6.7 Å | usable, marginal |
| K224 | K226 | ✔ | E272 @ 8.4 Å | usable |
| K382 | K384 | ✔ | E156 @ 8.2 Å | usable, peripheral |
| K193 | K195 | ✔ | E156 @ 9.5 Å | usable, peripheral |
| K344 | K346 | ✔ | **E214 @ 2.8 Å** | **VETO — salt-bridged** |
| K371 | Q373 | ✘ | E350 @ 8.1 Å | VETO — not conserved |
| R208 | N211 | ✘ | D204 @ 5.2 Å | VETO — both counts |

Two results here matter beyond bookkeeping:

**K344 is a trap.** It is a conserved, surface-exposed, apical lysine sitting
3.2 Å from Y211 and 3.8 Å from N348 — on paper the perfect answer, and it is
what a conservation-only analysis would have selected. It is salt-bridged to
**E214 at 2.8 Å**. Building the switch there would have produced designs that
bind *better* at pH 6.0. This is precisely the failure mode the brief's
non-negotiable warns about, and it is invisible unless carboxylate proximity is
checked explicitly.

**R208 is already compensated even in human** — D204 at 5.2 Å. The apex switch
the brief earmarked for class C is weaker than §3 assumes, independently of the
mouse problem.

### K374 — the recommendation's load-bearing residue

```
K374 → mouse K376        conserved
rel. SASA 0.282          solvent-exposed
nearest carboxylate      >10 Å — nothing
distance to R208          6.7 Å
distance to K371          7.9 Å
distance to Y211          7.9 Å
distance to Tf site      40.6 Å
distance to N251 sequon  15.2 Å
```

K374 is **within 8 Å of three of the nine epitope residues**. A binder with a
typical 600–900 Å² footprint that engages the apical apex will cover it
whether or not it is nominated. It is conserved, exposed, carboxylate-free and
far from both transferrin and every glycan.

**This is the residue §3 said did not exist.** §3's claim — "the only cations
in the apical patch are R208 and K371" — is true of the strict nine-residue
VP1u footprint and false of the binder-sized neighbourhood around it.

---

## 5. Patch scan — where the epitope should go

Every surface residue was taken as the centre of a 13 Å candidate footprint and
scored (`00b_patch_scan.py`). 244 patches; 55 pass Tf/glycan/cation/area
filters. The three that matter:

| Centre | Domain | n | %cons | Area Å² | Conserved cations | Anions | d(apex) | d(Tf) | d(glycan) | VP1u overlap |
|---|---|---|---|---|---|---|---|---|---|---|
| **S370** | apical | 17 | **53 %** | 1493 | **K231, K344, K374** | 2 | 3.1 Å | 33.4 Å | 8.7 Å | **5 of 9** |
| A225 | apical | 15 | **60 %** | 1272 | K193, K224, K374, K382 | 1 | 6.4 Å | 29.3 Å | 9.5 Å | **0 of 9** |
| Q520 | protease-like | 16 | 62 % | 1362 | K495, K508, K514 | 1 | 39.8 Å | 12.5 Å | 32.4 Å | 0 of 9 |

The trade is explicit and three-way:

- **Q520 and friends** are the most conserved surface on the protein (60–62 %)
  and are useless for objective 3's functional-epitope credit — 40 Å from the
  apex, on the protease-like domain, nowhere near anything the receptor does.
  They also sit 12 Å from the transferrin site.
- **A225** is the best-conserved patch *near* the apex (60 %, four conserved
  cations) but touches **none** of the nine epitope residues. It is the
  cross-reactivity-maximising choice that forfeits the functional-epitope
  claim and the 11RN positive control.
- **S370** covers **five of the nine** epitope residues, is 3.1 Å from the apex
  centroid, carries **K374** plus K231, is the furthest of the three from
  transferrin (33 Å), and costs 7 points of conservation against A225.

---

## 6. Recommended epitope

> **The apical shoulder, centred on the S370/K374 locus** — the face spanning
> roughly **209–215, 343–351, 362–374**, deliberately overlapping the apical
> apex at **L209 / Y211 / N348 / K371** while extending onto the conserved
> shoulder at **K231 / K374**.
>
> **Switch anchor: K374** (conserved K→K376, no carboxylate within 10 Å).
> **Secondary: K231.** **Explicitly excluded from switch geometry: R208, K371
> (not conserved) and K344 (E214 at 2.8 Å).**

Reasoning, in the order the objectives are ranked:

1. **pH selectivity (obj. 1).** K374 gives a genuine strategy-1 interface
   switch — a real target cation, correctly charged, with no carboxylate to
   invert it. This is the strong mechanism, not the weak intrinsic one.
2. **Cross-reactivity (obj. 2).** K374 is conserved, so *the switch itself is
   species-independent* — the thing that was supposed to be impossible here.
   The surrounding patch is 53 % conserved, which is honest-to-middling: better
   than the 50 % apex core, worse than the 60–62 % available off-epitope.
3. **Affinity + functional epitope (obj. 3).** Five of nine VP1u contacts are
   in footprint, so the "name the epitope and link your evidence" claim is
   defensible and **11RN's VP1u RBD becomes a directly relevant positive
   control** for calibration (brief §8, 04c). 33 Å from the transferrin site —
   receptor recycling untouched.

### Why not simply take A225 and maximise conservation

Because objective 2 only ranks designs **whose human KD is below 1 µM**. A
binder on the protease-like domain or the off-epitope shoulder must first bind
tightly, and we have no positive control there, no natural binder to calibrate
against, and no evidence the surface is ligandable at all. The apical face is
the one surface on this protein with multiple independent natural and designed
binders — VP1u, Machupo GP1, the 6WRV minibinders. **S370 keeps the ligandable,
controllable surface and pays 7 points of conservation for it.** A225 is the
right hedge, not the right primary.

### Honest statement of what this does and does not solve

**Solved:** the switch no longer depends on a species-specific residue. §3's
central conflict — "the residues that power objective 1 are the residues that
break objective 2" — is **false as stated**, because K374 powers objective 1
and is conserved.

**Not solved:** the epitope core is still ~50 % conserved. The residues that
determine *affinity* — the shape and chemistry a binder reads at the apex —
still differ at six of nine positions. Cross-reactivity remains the genuine
long shot the brief calls it in §12. Stage 00 has decoupled the pH mechanism
from the species problem; it has not made the species problem go away.

---

## 7. Q3 — the 318/633/719 seam lead: discard

The brief flagged an unverified preprint claiming a conserved,
transferrin-sparing VHH epitope at H318 / K633 / R719, and asked the map to
corroborate or discard it. **The map discards it**, on all three residues:

| Residue | Mouse | Conserved | rel.SASA | d(Tf) | d(glycan) | Problem |
|---|---|---|---|---|---|---|
| H318 | H320 | ✔ | 0.798 | 32.7 Å | **1.3 Å** | **sits on the N317 sequon — glycan-occluded in the assay.** 49 Å from the apex. |
| K633 | **R636** | ✘ | 0.671 | **0.0 Å** | 14.0 Å | **inside the transferrin site.** Would block recycling and forfeit the functional-epitope credit. |
| R719 | R722 | ✔ | 0.221 | **2.8 Å** | 6.8 Å | immediately adjacent to the transferrin site; marginally exposed. |

The three residues do not form a usable patch: one is under a glycan, one is in
the transferrin footprint, one abuts it. Note K633→R636 is a *conservative*
substitution that keeps the positive charge — so a naive conservation screen
scoring charge rather than identity would have passed it. **No design should be
built on this seam. O5 in §11 is closed; the preprint needs no further
verification because nothing now depends on it.**

---

## 8. Native histidines

Twelve histidines in chain A: **186, 234, 300, 302, 318, 401, 475, 515, 603,
684, 699, 707.** (The construct also contains these twelve — consistent.)

**The favourable prior holds, and now spatially, not just in sequence.** The
nearest His to the apical apex centroid is **H234 at 16.9 Å — and it is
completely buried (rel. SASA 0.000).** The next are H300 (26.6 Å) and H302
(28.4 Å). Nothing titratable sits in or near the recommended epitope.

- **H318 is the only strongly exposed one** (rel. SASA 0.798) and it is 49.2 Å
  from the apex and 1.3 Å from the N317 glycan sequon. Irrelevant to this
  epitope and glycan-shadowed anyway.
- **Position 349 is methionine.** The "His318 and His349" claim from background
  reading is confirmed wrong — carry only H318 forward.
- **H401 → mouse R403** is the only non-conserved histidine. Buried
  (rel. SASA 0.048), 60 Å from the apex, irrelevant here, but worth noting
  because it is a His→Arg swap and would be a species-divergent pH sensor for
  anyone targeting that face.

`native_his_ddg` in stage 05 must still be computed, not skipped — the prior is
favourable, which is a reason to expect a small number, not to omit it.

---

## 9. Portfolio implication

§7's pre-committed shape survives with class A and B both sited on the apical
shoulder rather than split between a "conserved patch" and the apex:

| Class | Proposed n | Site | Switch |
|---|---|---|---|
| A | 7 | apical shoulder (S370/K374) | **strategy 1 on K374** |
| B | 4 | apical shoulder (S370/K374) | **strategy 2**, buried binder-core His network — species-independent insurance |
| C | 4 | apical apex proper (R208/K371) | strategy 1 — affinity hedge, mouse expected to fail |
| D | 3 | apical shoulder / apex | **no deliberate switch** — the control for the switch hypothesis |
| E | 2 | the **348/371 swap locus** (§3 above) | strategy 1 sited to reach K371 in human and K350 in mouse |

Class E is new and is the speculative slot — two designs on the rearranged-cation
idea, which is high-upside and rests on an unmeasurable side-chain inference.

**This is a proposal for the gate, not a decision.** No hotspots have been
fixed, no backbones generated, no GPU job started.

---

## 10. Not computed — read before relying on any of this

- **No pH was modelled anywhere in this stage.** Everything above is sequence
  and geometry. The switch claims are *mechanism hypotheses* to be scored by
  thermodynamic linkage in stage 05, never by a folding model.
- **The mouse structure does not exist.** Every mouse statement is from the
  sequence alignment onto the *human* structure. Mouse side-chain geometry —
  including, critically, whether K350 projects into the interface — is
  unmodelled. Stage 06b's mouse target will be a homology model and must be
  labelled as one.
- SASA is from the **isolated monomer**. Dimer-seam burial was computed and
  used as a filter but the recommended epitope is well away from the seam.
- **Glycans were stripped** (6 NAG removed); N251/N317/N727 and neighbours are
  masked as non-designable. The recommended patch's closest approach to a
  sequon is **8.7 Å (to N251, via K231)** — clear of the Asn itself, but a real
  N-glycan extends well beyond its attachment point, so **the K231 end of the
  patch carries genuine glycan-shadow risk that no number here captures.**
  K374, at 15.2 Å from N251, does not.
- The 13 Å footprint radius and the 0.25 rel.SASA surface cutoff are
  conventions chosen to approximate a 600–900 Å² binder footprint, not
  measurements.
- Conservation here is **identity**. K633→R636 shows why charge-class
  conservation and identity conservation can disagree; both are in the table,
  neither is the whole story.
