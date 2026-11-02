# Methodology summary — TfR1 Challenge 04

**Pipeline:** RFdiffusion → ProteinMPNN → Boltz-2, with pH selectivity scored
separately by thermodynamic linkage (PROPKA). 20 de novo single-chain binders,
61–118 aa.

**Epitope chosen by measurement, not literature.** We built a human/mouse
surface conservation map on PDB 11QC and computed the B19V VP1u footprint from
11RN (exactly nine residues: 202, 204, 208, 209, 210, 211, 212, 348, 371). The
pH switch mechanism needs a target Arg/Lys next to an interface histidine —
protonation at pH 6.0 creates charge repulsion and the binder releases. The
apical cations R208 and K371 are both neutral in mouse, so a switch there
cannot be cross-reactive. K205, the obvious candidate, is a serine in mouse and
is not in the epitope at all. **K374 is the one usable conserved anchor**
(K→K376, exposed, no carboxylate within 10 Å). Two apparent alternatives were
vetoed: K344 is conserved but salt-bridged to E214 at 2.8 Å, which would make
binding *stronger* at pH 6.0; K231 sits inside the N251 glycan shadow. Target
epitope is the apical shoulder around S370/K374 — overlapping 5 of the 9 VP1u
contacts, 33 Å clear of the transferrin site so receptor recycling is
unaffected.

**Target prep.** Monomer, not dimer (11RN shows independent A2B2 binding).
Verified zero sequence mismatches against the assay construct across all 639
modelled residues, confirmed Gly142, stripped glycans and masked the three
sequons plus neighbours. The mouse construct was taken from the challenge page,
not UniProt, and cross-validated against all six page-stated apical differences.

**Generation and screening.** 710 backbones across three epitope variants and
three length ranges; ~6,300 ProteinMPNN sequences in three pools (plain,
global His-bias, and a per-position His bias targeted at residues facing the
conserved cation); 1,372 scored with Boltz-2 against the **full 639-residue
monomer**, not a crop.

**Filters.** Four hard gates beyond confidence:

1. **Epitope fidelity.** Scoring against a full monomer lets the binder dock
   anywhere, and it does — only 19% land on target. Confidence does not catch
   this: the highest-scoring design of all 900 contacts the transferrin site
   and a glycan sequon with zero targeted hotspots. We require ≥2 hotspot
   contacts, **zero** transferrin-site contacts, minimal glycan-shadow contact,
   and the C-terminus ≥8 Å from the target so the Twin-Strep tag cannot
   participate in binding.
2. **pH linkage with a conservation policy.** PROPKA pKa shifts give the
   ddG between pH 7.4 and 6.0. Run bare it found 497 sites, two thirds of them
   on cations we had already ruled out — so sites on vetoed cations are
   discarded and the rest are labelled human-only or cross-species.
3. **His refold.** Every histidine variant is re-folded to confirm the pH 7.4
   interface survives. 57 of 221 did.
4. **Two-seed pose agreement.** Scores agreeing is not poses agreeing. 49 of 57
   placed the binder in the same site under a second seed; the rest were cut.

**Portfolio.** Ranked in the challenge's own objective order, as four classes
rather than one list: 6 cross-reactive, 7 human pH switches, 2 affinity, and
**5 unmutated controls** that carry no deliberate switch and serve as the
internal comparison for the switch hypothesis. Diversity caps were applied last
and allowed to displace higher-scoring designs: **20 distinct backbones**,
worst pairwise identity 29%.

**Honest limits.** No structure-prediction model used here has a pH parameter —
pH selectivity is a mechanism hypothesis from thermodynamic linkage, not a
predicted KD ratio. No mouse structure exists; all mouse numbers come from
threading the mouse sequence onto the human structure. Best switch is ~16×
(human-only); the 100× target was unreachable on this target. A conserved-cation
switch turned out **not** to imply a cross-reactive binder — the switch is one
residue, the interface is 15–20 — so cross-reactivity is ranked on measured
mouse binding. **No design is claimed to be both strongly pH-selective and
cross-reactive:** class A cross-reacts but switches weakly, class B switches
well but is human-only.
