# Challenge 04 — TfR1 conditional cross-species binder: Claude Code handoff brief

Written 2026-10-31 ~20:00 Melbourne. Source of truth:
https://proteinbase.com/competitions/anthropic-adaptyv-2026/challenges/tfr1

**Deadline: Sun 1 Nov 23:59 AoE = Mon 2 Nov 22:59 Melbourne (AEDT). 51 hours
from the time this was written.**

Labels used throughout: **[OFFICIAL]** = competition rule, taken from the
challenge page; **[EVIDENCE]** = published experimental finding;
**[HYPOTHESIS]** = structural proxy with no validation against the actual
experimental readout; **[UNVERIFIED]** = claim I could not check and which must
be checked before it is relied on.

---

## 0. Working agreement for the agent

**First action, before anything else.** Run the audit and write it into the run
log before touching any pipeline stage:

```
ls ~/                       # find the challenge-02/03 trees
nvidia-smi                  # GPU state, driver, free VRAM
git -C <repo> rev-parse --show-toplevel
conda env list
ls ~/tnf-track3/RFdiffusion/models/ ~/tnf-track3/ProteinMPNN/
docker images | grep -i rfdiff
```

Record: actual directory layout, which tool checkouts exist and at what commit,
which weights are already on disk, which conda/venv environments work, GPU
state. **Do not assume any path or checkpoint from this brief — verify.** If a
model or weight set is missing, say so immediately rather than pulling a large
checkpoint unprompted.

**Git identity.** Commits must be authored by Harjas, not by the agent. The
identity is already configured in the repo — **verify it before the first
commit and stop if it is wrong**:

```bash
git -C <repo> config user.name    # expect: harjas2001
git -C <repo> config user.email   # must be an email VERIFIED on his GitHub
```

GitHub attributes a commit to a profile by matching the **author email** to a
verified address on the account. If the email does not match, the commits land
as an unattributed author and no amount of fixing later is cheap. Do not change
these values yourself — report a mismatch and let Harjas set it.

**Git.** New tree for challenge 04. Work on `main` there, or on branch
`tfr1-challenge04` if it is a subdirectory of an existing repo. Never commit
onto the challenge-02 or challenge-03 branches; never modify or delete their
artifacts.

**Long jobs run in their own detached tmux session, not in the foreground.**
The agent's own session is already inside tmux; a foreground GPU job blocks the
agent for hours and dies with the pane. Pattern:

```bash
tmux new-session -d -s rfdiff 'bash -lc "cd ~/<repo> && source env.sh && \
    ./rfdiffusion.sh ... 2>&1 | tee logs/02_rfdiffusion.log"'
tmux capture-pane -pt rfdiff | tail -20     # poll
tmux list-sessions                           # still running?
```

Poll the log rather than blocking. Every long batch checkpoints so it resumes
rather than restarts — with 51 hours there is no budget for a lost run.

- **Commit constantly and atomically.** Pre-authorised. Message format
  `tfr1: <what changed>` with counts where relevant
  (`tfr1: RFdiffusion batch 1, 600 backbones, seed 11`).
- Commits are local and need no credentials. If one fails, fix it once if the
  fix is obvious, otherwise log it and carry on. Never stall the science on git.
- Commit code, configs, manifests, metrics tables, run logs. Do **not** commit
  weights or raw structure dumps over ~50 MB — gitignore them, record paths in
  the manifest.
- **No pushes, no force-push, no history rewriting, no submission to the
  competition site without explicit approval from Harjas.**
- Checkpoint every long batch so it resumes rather than restarts.

**Reuse, don't rewrite.** Challenge 02 (`tnf-upload/scripts/`) is the base —
it is the same problem shape. Port status:

| Script | Action |
|---|---|
| `02_crop_target.py` | reuse as-is, already generic (`--pdb`, `--hotspots`, `--radius`) |
| `05_ph_scan.py` | reuse as-is — linkage physics is target-agnostic. The key asset. |
| `06_his_refold.py` | reuse as-is |
| `03_mpnn.py` | chain-ID constants only |
| `04_boltz.py` | swap the hardcoded target sequence (line ~35); add the mouse arm |
| `07_submit.py` | re-rank for three objectives instead of one |
| `01_prep_target.py` | **full rewrite** — TNF-specific throughout |

From challenge 03 (`glp1r-comp3/scripts/`), port three modules challenge 02
lacked:

- `af2ig/` — PyRosetta-free AF2 initial-guess. Restores the Bennett 2023
  `pae_interaction` filter, the best-validated binder filter available, which
  challenge 02 had to substitute with a calibrated Boltz cutoff.
- `42_selfconsistency.py` — two-seed pose agreement. In challenge 03 this
  removed 22 of 46 candidates and caught a case where two seeds agreed on
  interface confidence to within 0.03 while placing the binder 37 Å apart.
  **Single best filter in either repo. Do not skip it.**
- `20_denovo_check.py` — novelty screen, with a TfR1-specific reference set
  (see §8).

Leave the rest of challenge 03 alone; it is GPCR-specific.

---

## 1. Verified challenge spec [OFFICIAL]

| Item | Value |
|---|---|
| Target | Human TfR1 ectodomain, UniProt **P02786-1, residues 121–760** (640 aa) |
| Assay construct | Ectodomain 121–760, **tag-free** (Gly-Pro scar before the first residue, not modelled) |
| Mouse construct | Ectodomain **123–763**, tag-free |
| Oligomeric state | **Homodimer in solution** — two copies of the epitope per target |
| Structure reference | PDB **11QC** chains A–B |
| Apical complex | PDB **11RN** chain A with VP1u |
| Recommended epitope | **Apical domain** (a hint, not a rule) |
| Identity | Human/mouse ectodomains **78% identical**; apical patch differs at **6 of 9** positions |
| Length / format | 10–250 aa; single-chain protein, nanobody, scFv or Fab; **de novo only** |
| Cap | **20** designs (Tracks 2 and 3), 40 (Track 1) — see open item O3 |
| Capture | Design immobilised via **C-terminal Twin-Strep tag**, which **must not be part of the binding** |
| Binding gate | Fittable sensorgram giving a KD, **or** association shift **>300%** over negative control. Otherwise "no detectable binding" |
| Conditions | Human pH 7.4, human pH 6.0, mouse pH 7.4 — design immobilised, target in solution (see open item O1) |

**Ranking [OFFICIAL], in this order:**

1. **pH-selective binding** — measurable KD against human TfR1 at pH 7.4 and no
   detectable binding at pH 6.0. Designs binding at both rank after these, by
   the ratio KD(6.0)/KD(7.4). *"Designs that rely on a conformational change or
   a designed pH mechanism are not penalised when we select designs for
   testing; the measurement decides."*
2. **Mouse cross-reactivity** — ranked by how close the mouse:human KD ratio is
   to 1, **among designs whose human KD is below 1 µM**.
3. **Affinity** — lowest KD against human TfR1 at pH 7.4. Binding at a
   functional epitope counts in its favour: *"name the epitope you target and
   link your evidence, or we read it from a predicted complex."*

The winner must clear all three. Each step also awards an honourable mention to
the best design that did not win.

**Selection for testing is separate from winning.** Track 1 teams get their
ranked designs tested; for Tracks 2 and 3, Claude selects from all submissions
using a published prompt. **The submitted metrics, structures and method
descriptions are part of what it reads.** Write the per-design dossiers
accordingly — named epitope, linked evidence, declared pH mechanism, honest
not-computed list. Overclaiming is the cheapest way to look unserious.

### Structures verified against RCSB, 2026-10-31

Both accessions in the challenge page are real and were confirmed directly —
they are **not** OCR errors for 1IQC/1IRN. Both were released 2026-05-27, which
is past the model training cutoff, so they will not be in prior knowledge.

| PDB | Content | Method | Resolution | Chains |
|---|---|---|---|---|
| **11QC** | Human TfR1 ectodomain | cryo-EM | 2.40 Å | A, B — 639 res each, *H. sapiens* |
| **11RN** | B19V VP1u RBD + human TfR1 ectodomain | cryo-EM | 2.40 Å | A, B = TfR1 (637 res); C, D (author X, Y) = VP1u (62 res) |

**11RN is A2B2 — one 62-residue VP1u RBD per TfR1 subunit, bound
independently.** This is load-bearing: unlike challenge 02, where the receptor
groove genuinely spanned two protomers, **the apical epitope here is a
single-subunit epitope.** Design against a monomer. That is the difference
between a feasible and an infeasible compute budget (§5).

11RN also carries **6 NAG instances, 3 per TfR1 chain** — the N251/N317/N727
glycosylation sites. Handle per §4.

### The assay sequence is authoritative — use it, not UniProt

`target/human_tfr1_121_760.fasta` ships with this brief. It is the sequence
**from the challenge page**, not from UniProt, and it is what the assay measures.
Verified 2026-10-31: **640 residues, spanning 121–760 exactly.**

Checks already run against it, so you do not need to repeat them:

| Check | Result |
|---|---|
| Length / span | 640 aa, 121–760 ✓ |
| **Position 142** | **G142 — the construct carries Gly**, as the page states |
| Apical patch | D204, R208, V210, Y211, L212, N348, K371 all confirmed at the stated positions ✓ |
| Apical main loop 202–215 | `IVDKNGRLVYLVEN` |
| Glycosylation sequons (N-X-S/T) | exactly three: **N251, N317, N727** — independently matching 11RN's 3 NAG per chain ✓ |
| Transferrin site 619–670 | `LSFVRDLNQYRADIKEMGLSLQWLYSARGDFFRATSRLTTDFGNAEKTDRFV` |

**G142 is confirmed in the construct.** 11QC should match, but verify it on the
downloaded coordinates. If you model from 1CX8, 3KAS, 3S9L or 6WRV, patch
Ser142 → Gly142.

**The mouse sequence is not in this brief.** The page's Mouse TfR1 tab is
JavaScript-gated and did not render to a plain fetch. **Pull it from the page in
a browser** — do not substitute UniProt Q62351 without checking it against the
page first, since the construct is what gets assayed. Mouse construct is stated
as ectodomain **123–763**, tag-free.

---

## 2. Target biology, the parts that matter

**Architecture [EVIDENCE].** Each TfR1 monomer has three domains: a
protease-like domain, an **apical domain** (the tip, where viruses and
brain-shuttle antibodies bind), and a **helical domain** which mediates the
dimer contacts and carries the transferrin site. Disulfide-linked homodimer.

**The apical epitope [OFFICIAL, from the challenge page].** Residues
**202–215, 348, 371**, the patch around **Tyr211**. Parvovirus B19 (11RN),
Machupo virus (3KAS) and the designed binders of 6WRV all land there.

**Transferrin binds the helical domain, residues 619–670**, and only grazes the
apical domain of the partner subunit. Targeting the apical domain therefore
leaves transferrin docking and receptor recycling intact — which is the stated
intent of the challenge and part of what "functional epitope" means in
objective 3.

**The human→mouse apical differences [OFFICIAL]:**

| Human | Mouse |
|---|---|
| D204 | Q207 |
| R208 | N211 |
| V210 | D213 |
| Y211 | P214 |
| N348 | K350 |
| K371 | Q373 |
| L212 | *no counterpart* |

Six of nine positions, plus an indel. This is why TfR1 binders are
"notoriously species-specific" and why brain-shuttle antibodies are normally
made twice.

**Gly142 vs Ser142 [OFFICIAL].** UniProt and the assay construct carry
**Gly142**; most older crystal structures (1CX8, 3KAS, 3S9L, 6WRV) carry
Ser142, the other common allele. The page states 11QC and 11RN match UniProt at
every modelled residue. **Verify this on the downloaded coordinates rather than
trusting it** — this is the direct analogue of challenge 02's Leu143→Asp patch,
which was real and would have silently corrupted the target. If you model from
an older entry, correct position 142.

---

## 3. The central conflict — read this before choosing hotspots

**Objectives 1 and 2 fight each other at the apical apex, mechanically.**

Challenge 02's switch worked by placing interface histidines next to a target
Arg/Lys: protonation at pH 6.0 creates charge–charge repulsion and the complex
falls apart. That epitope was chosen *because* B31/B32 were an exposed arginine
pair. The counter-charge was the whole mechanism.

Apply the same logic here. The only cations in the apical patch are **R208** and
**K371**. In mouse they become **N211** and **Q373** — both neutral.

So an interface-repulsion switch built at the apical apex:

- works in human,
- has **no counter-charge at all in mouse**, and
- sits on the single most divergent surface on the receptor.

The residues that power objective 1 are the residues that break objective 2.
Running challenge 02's strategy unchanged optimises one objective by sacrificing
the other.

### The way out: move the switch into the binder

Baker-lab work on designed pH-sensitive binders [EVIDENCE, preprint] describes
two mechanisms, not one:

- **Strategy 1 — interface destabilisation.** Interface His as H-bond acceptor
  beside a target cation. Largest effects reported (up to ~1000× on EphA2,
  122× on TNFR2), requires ≥2 protonated-His–cation contacts. **Species-
  dependent: needs the target to supply the cation.**
- **Strategy 2 — monomer destabilisation.** A buried **His–His / His–Arg /
  His–Lys network in the binder's own core**. Protonation destabilises the
  binder's fold and speeds dissociation. Smaller effects reported (best ~79× on
  TNFα, 6× IL-6, 3.5× PCSK9). **Species-independent — the target supplies
  nothing, so the switch behaves identically against human and mouse.**

Strategy 2 decouples objective 1 from objective 2 entirely. It is weaker, so the
portfolio hedges across both (§7).

**Also relevant:** random histidine incorporation does **not** produce pH
sensitivity in that work — simply raising His count gave little effect. Donor/
acceptor status and electrostatic context are what matter. Challenge 02's
`05_ph_scan.py` already encodes this correctly (cation within 6 Å good, anion
within 6 Å disqualifying, >8 Å His–His separation enforced). Keep those rules.

**Native target histidines — corrected against the actual sequence.**
Challenge 02 learned the hard way that pre-existing His must be counted, with
the correct unbound reference, and that the contribution can be negative.

Background reading described TfR1's native pH sensors as "His318 and His349".
**Checked against the assay construct: H318 is real, but position 349 is
methionine, not histidine.** Do not carry the His349 claim forward.

The construct contains **12 histidines total**: 186, 234, 300, 302, **318**,
401, 475, 515, 603, 684, 699, 707. Four fall in the apical domain (234, 300,
302, 318). **None is within 15 sequence positions of the apical patch**, which
is a favourable prior — native-His interference at the apical epitope looks
low. Confirm it spatially on 11QC, since sequence distance is not structural
distance. The `native_his_ddg` term in stage 05 must still be computed, not
skipped.

---

## 4. Stage 00 — the conservation map (decision gate, do this first)

**The epitope is not chosen in this brief. It is measured.** This is a few
minutes of CPU and it converts the biggest decision in the campaign from a
literature argument into a number.

Build `scripts/00_conservation.py`:

1. Fetch **11QC** (and **11RN** for the VP1u footprint). Fetch from
   `https://files.rcsb.org/download/<ID>.cif` directly — Biopython `PDBList`
   uses a retired endpoint and 403s on everything.
2. Pull human **P02786** and mouse **Q62351** sequences; align the ectodomains
   (human 121–760, mouse 123–763). Report global identity and confirm it lands
   near the stated 78% — if it does not, the alignment or the isoform is wrong.
3. Map identity onto 11QC chain A per residue.
4. Compute per-residue SASA (Shrake-Rupley, as in challenge 02's
   `01_prep_target.py`).
5. Emit, for every **surface-exposed** residue (SASA above threshold), a table:
   residue, domain, conserved yes/no, charge class, distance to the apical apex
   centroid, distance to the transferrin site (619–670), distance to the nearest
   glycan sequon (N251/N317/N727).
6. **Then answer the one question that decides the campaign:** are there
   **conserved, surface-exposed Arg/Lys residues** in a patch that is (a) not
   the transferrin site, (b) not glycan-occluded, and (c) contiguous enough to
   host a binder footprint of ~600–900 Å²?

7. **Resolve the nine-residue patch.** The page names only seven of the nine
   apical patch residues — the six that differ plus L212. **Two or three are
   conserved and unnamed, and their identity may decide the whole campaign.**
   **K205 is the specific one to check.** It sits in the apical main loop
   (`IVDKNGRLVYLVEN`, 202–215), it is a cation, and the page does not list it
   as differing. **If K205 is conserved in mouse, strategy 1 becomes viable
   cross-species at the organiser-recommended epitope** — a switch built on
   K205 rather than R208/K371 would repel in both species, and the conflict in
   §3 dissolves. That is the single best outcome available and it costs one
   alignment to find out. Check it first.
   Note also that mouse *gains* a lysine where human has N348 (→ K350), so the
   cation inventory is not simply depleted in mouse — it is rearranged.

Output `target/conservation_report.md` with the candidate patches ranked, and
**stop there for a decision from Harjas.**

Decision rules, pre-committed so the gate is fast:

- **K205 (or another apical-patch cation) is conserved in mouse** → stay at the
  apical apex, build strategy-1 switches on the *conserved* cation only, and
  explicitly exclude R208 and K371 from the switch geometry. Best case: the
  recommended epitope, a strong switch, and cross-reactivity all at once.
- **A conserved cationic patch exists elsewhere, Tf-sparing** → that patch is the
  primary epitope. It supports strategy 1 *and* cross-reactivity.
- **Conserved patches exist but with no cations** → design there, use strategy 2
  for the switch, and keep a minority of apical-apex slots running strategy 1 as
  the high-affinity hedge.
- **No usable conserved patch** → apical apex primary, strategy 1, and accept
  that objective 2 is likely forfeit. Say so explicitly in the methods report
  rather than implying cross-reactivity was attempted.

[UNVERIFIED] My background research surfaced a 2025 preprint reporting a
camelid VHH family binding a conserved, transferrin-sparing epitope at the
apical–helical seam around **His318 / Lys633 / Arg719**, with genuine
human+mouse reactivity. **I did not verify this preprint, and Lys633 falls
inside the transferrin site (619–670)**, so taking it at face value risks
blocking transferrin and losing the functional-epitope credit in objective 3.
Treat it as a hypothesis the conservation map should test, not as a result.
If the map independently puts a conserved cationic patch near 318/719, that is
corroboration; if it does not, discard the lead.

**Glycans.** Strip NAG from the coordinates, but mask N251, N317 and N727 and
their immediate neighbours as non-designable. In silico they look like ordinary
surface; in the assay they are covered. A binder that parks on a glycan site
will fail for a reason no score will show.

---

## 5. Compute budget — this is what kills the run

Measured rates from challenge 02 on the same A100-80GB:

| Operation | Measured |
|---|---|
| RFdiffusion, 382 total residues (binder + 2-protomer target) | 1.84 min/backbone |
| RFdiffusion, 263 total residues (15 Å crop) | 1.09 min/backbone |
| Boltz-2, templated complex, batched | 20.7 s/design |
| PROPKA linkage scan, 10 processes / 12 cores | ~6 s/design |

TfR1 is much larger. Monomer 639 aa + an 80 aa binder = ~719 residues.
Extrapolating the two measured RFdiffusion points (which scale somewhere
between ^1.4 and ^2):

| Configuration | Total residues | Est. min/backbone | 1000 backbones |
|---|---|---|---|
| Dimer, uncropped | ~1358 | ~20–45 | infeasible |
| Monomer, uncropped | ~719 | ~4.4–8.1 | **73–135 GPU-h — infeasible** |
| Monomer, ~180-residue crop | ~260 | ~1.1 | **~18 GPU-h — fits** |

These are **estimates, not measurements**. The pilot in §9 replaces them with
real numbers before anything scales.

**Consequences, all mandatory:**

1. **Crop.** Reuse `02_crop_target.py` with a 15 Å shell, targeting ~180
   residues kept. The crop is for generation only; the uncropped monomer stays
   the validation target, exactly as in challenge 02.
2. **Monomer, not dimer** — justified by 11RN's independent A2B2 binding.
   Revisit only if the conservation map sends the epitope to the dimer seam, in
   which case the budget must be recut before generating.
3. **Scoring also needs a crop decision.** Challenge 02 could validate against
   the full 314-residue target at 20 s/design. The full TfR1 dimer plus binder
   is ~1358 residues — likely to be an order of magnitude slower or to OOM.
   **Time one Boltz call on the full monomer and one on the crop during the
   pilot, before committing.**
4. **Mouse is a second-stage filter on human survivors.** Roughly 80% of the
   mouse scoring bill disappears and nothing is lost — a design that fails
   against human is not going to be submitted on mouse strength.

---

## 6. On the newer models — do not swap the engine

Harjas raised RFD4-Proteina, Mol-1, OpenDDE, chai-lab, RFantibody, OpenFold3,
Genie3 and FreeBindCraft. With 51 hours the verdict is: **keep
RFdiffusion → ProteinMPNN → Boltz-2.**

- It is measured, debugged, and has four documented failure modes already fixed
  (template pinning, the PROPKA imidazole bug, native-His accounting, per-pair
  scoring).
- `chai-lab` is Chai-1 — prediction only. Chai-2's design capability, the one
  with the reported hit rates, is not in the open repo.
- RFD4-Proteina, BoltzGen and Mol-1 are days-to-weeks old. Setup risk with no
  recovery time.
- AlphaProteo has no public weights.
- The real wins this round are the conservation map, the mouse arm and
  challenge 03's AF2-IG — none of which need a new model.

**One exception.** If the pilot shows the generation stage coming in
comfortably under budget, FreeBindCraft is worth adding as a *parallel backbone
source* for fold diversity — it is PyRosetta-free, which is why it is the only
one of the eight that fits the existing constraints. Secondary, never primary,
and only with time already banked.

---

## 7. Portfolio (20 designs)

Final allocation is set after the stage-00 gate. Pre-committed shape:

| Class | n | Switch mechanism | Rationale |
|---|---|---|---|
| **A. Conserved-patch, intrinsic switch** | 8 | Strategy 2 (buried His network in binder core) | The only architecture that can satisfy objectives 1 and 2 simultaneously |
| **B. Conserved-patch, interface switch** | 4 | Strategy 1, only if the map finds a conserved cation | Stronger switch where the target supplies the counter-charge in both species |
| **C. Apical apex, interface switch** | 5 | Strategy 1 on R208/K371 | Maximises objectives 1 and 3; expect mouse to fail. The affinity hedge. |
| **D. Apical apex, no deliberate switch** | 3 | none | Pure objective-3 insurance and the internal control for the switch hypothesis |

If stage 00 finds no usable conserved patch, collapse A and B into C and D and
state the forfeit plainly.

Constraints on the slate:

- No two designs above **70% sequence identity**.
- Cluster at **TM-score 0.6**, cap per cluster. Challenge 02 shipped 6 clusters
  across 20 and called that its weakest axis — beat it.
- Prior Adaptyv rounds show in-silico scores predict real binding at roughly
  chance (ROC-AUC ≈ 0.5–0.64). **Twenty variations on one scaffold is one bet,
  not twenty.** The diversity cap costs predicted quality on purpose.
- **C-terminus must point away from the interface** on every design — the
  Twin-Strep capture tag attaches there and must not participate in binding.
  Check this on the predicted complex; it is an [OFFICIAL] requirement and a
  silent disqualifier.
- Composition: no free Cys, no N-X-S/T sequons, no homopolymer runs ≥5, watch
  the ProteinMPNN poly-Ala/poly-Glu bias that produced two awkward designs in
  challenge 02.

---

## 8. Filter stack

| Stage | Script | Notes |
|---|---|---|
| 00 | `00_conservation.py` | **new** — the decision gate, §4 |
| 01 | `01_prep_target.py` | **rewrite** — 11QC, Gly142 check, glycan masking, numbering map, hotspots from the map not the literature |
| 02 | `02_crop_target.py` | reuse — 15 Å shell, ~180 residues |
| 02b | `rfdiffusion.sh` | reuse wrapper; Hydra overrides only, `hydra.run.dir=/work/.hydra` |
| 03 | `03_mpnn.py` | reuse — soluble weights, two pools (plain + His-biased) |
| 04 | `04_boltz.py` | swap target; template-pin the target (challenge 02 correction #1); score **binder↔target pairs only** |
| 04b/c | `04b_make_controls.py`, `04c_calibrate.py` | reuse — shuffled-sequence negatives. **Add positives this time** via `--positives`: VP1u RBD from 11RN is an ideal apical positive control and anchors the top end, which challenge 02 never did. |
| 04d | `af2ig/predict_pdb.py` | **port from challenge 03** — `pae_interaction` triage |
| 05 | `05_ph_scan.py` | reuse — add strategy-2 scoring for buried binder-core networks; `native_his_ddg` must include His318/His349 |
| 06 | `06_his_refold.py` | reuse — pH 7.4 interface must survive |
| 06b | `06b_mouse.py` | **new** — score survivors against the mouse target |
| 06c | `42_selfconsistency.py` | **port from challenge 03** — two-seed pose agreement. Hard gate. |
| 07 | `20_denovo_check.py` | **port** — TfR1 reference set below |
| 08 | `07_submit.py` | re-rank for three objectives, dossiers, rejected-candidate log |

**Mouse target construction.** Build a mouse TfR1 model (Q62351, 123–763) by
homology on 11QC, or predict it. Record which, and treat it as a model, not a
structure — mouse scores carry more uncertainty than human ones and the methods
report should say so.

**De novo reference set for the novelty check:** VP1u (11RN), Machupo GP1
(3KAS), the 6WRV designed minibinders, transferrin itself, and the published
brain-shuttle antibodies (8D3, OX26, the Roche and Denali TfR1 modules). No
contiguous run matching any of them.

---

## 9. 51-hour schedule

| Window | Work |
|---|---|
| 0–2 h | Audit, repo setup, fetch 11QC/11RN, **stage 00 conservation map** |
| 2–3 h | **Epitope decision gate with Harjas.** Do not generate before this closes. |
| 3–6 h | Target prep, crop, **pilot: ~10 designs end to end through every stage.** Record wall-clock per design per tool, including a Boltz timing on full monomer vs crop. **Do not scale before this passes.** |
| 6–22 h | Main RFdiffusion generation + ProteinMPNN (both pools) |
| 22–32 h | Boltz human screen, calibration with positives, AF2-IG triage |
| 32–40 h | pH scan (CPU, runs parallel with GPU work), His refold |
| 40–45 h | Mouse screen on survivors, self-consistency gate |
| 45–49 h | Ranking, portfolio selection, submission package, methods report, pre-submission audit |
| 49–51 h | Buffer |

If generation runs slower than the pilot predicts, **cut sample counts, not
validation stages.** Challenge 02's funnel guidance holds: under ~150 survivors
into the His scan is too thin, because the scan keeps only 10–20%. If the funnel
is thin, raise `num_designs` — do not loosen the confidence cutoff, since a
design failing confidence will not survive His substitution either.

---

## 10. Traps already hit in challenges 02 and 03 — do not rediscover

- **Biopython `PDBList` uses a retired RCSB endpoint and 403s.** Fetch from
  `https://files.rcsb.org/download/<ID>.cif` directly.
- A real PDB file has hundreds of header lines before the first ATOM — do not
  validate a download by scanning its first few KB.
- **The RFdiffusion image has an ENTRYPOINT.** Hydra overrides only — no
  `python`, no script path. Container-absolute paths under `/work`. Include
  `hydra.run.dir=/work/.hydra` or it fails on a permission error writing to
  root-owned `/app`.
- **Boltz does not assemble a multi-chain target on its own.** In challenge 02
  the target-internal pair scored 0.291 until both chains were template-pinned,
  at which point it rose to 0.708 and binder scores roughly doubled. If you use
  the dimer at any point, template-pin it and gate on a `target_pair_iptm`
  sanity check.
- **A histidine without its imidazole is invisible to PROPKA.** Challenge 02's
  scan silently returned zero candidates for every design because the mutation
  truncated the side chain to CB. The fixed version grafts a real imidazole on
  the N/CA/CB frame with χ1 chosen by clash. Use the fixed version.
- **Score binder↔target pairs only**, and take the minimum across pairs where
  the target has multiple chains.
- Check the `boltz` env's torch build against the driver — challenge 02 lost
  time to `torch 2.14.1+cu130` on a CUDA 12.4 driver making the GPU unusable.
- **Score agreement is not pose agreement.** Challenge 03 found two seeds
  agreeing to within 0.03 on interface confidence with the binder placed 37 Å
  apart. This is why the self-consistency gate exists.

---

## 11. Open items to close before submitting

- **O1 — the fourth condition. CONFIRMED as a real discrepancy, 2026-10-31.**
  The live page says designs are measured "in the same four conditions" and
  then lists three: human TfR1 at pH 7.4, human at pH 6.0, mouse at pH 7.4.
  Either a fourth exists (mouse pH 6.0 is the obvious candidate) or the page
  miscounts. **Ask in Slack.** If mouse pH 6.0 is measured, the
  species-independent strategy-2 switch becomes more valuable still, since it
  would switch in both species.
- **O2 — FAQ.** Pull the full FAQ with the page's *Copy Markdown* button; a
  plain fetch truncates it. Confirm submission field names and `molecule_class`
  values before building the CSV.
- **O3 — Track and cap.** Challenge 02 was entered as "Track 3". On the
  challenge-04 page all three objectives apply to every design and the ranking
  is a single ordered list, while the cap is stated as 40 for Track 1 and 20 for
  Tracks 2 and 3 — which reads as a participant category, not an objective
  split. **Confirm which track applies and whether the cap is 20 or 40.** Plan
  for 20; if it is 40, the portfolio in §7 scales by doubling each class rather
  than adding scaffolds to one.
- **O4 — selection prompt.** Read the published selection prompt before writing
  the dossiers. It is tab-gated and needs a browser.
- **O5 — the preprint in §4.** Verify or discard before any design is built on
  the 318/633/719 seam.

---

## 12. Honest expectation

Binding to TfR1 at the apical domain is tractable — it is a well-characterised,
accessible epitope with multiple natural and designed binders, and there is a
clean positive control in 11RN.

The other two objectives are harder, and in opposite directions. **Cross-
reactivity is the hard one**: the literature position is that TfR1 binders are
species-specific almost by default, and the apical patch differs at 6 of 9
positions, which is why brain-shuttle antibodies are normally made twice. **pH
selectivity is achievable but unmeasurable in advance** — no structure
prediction model has a pH parameter, so everything stage 05 produces is a
mechanism hypothesis scored by thermodynamic linkage, not a predicted KD ratio.

A realistic outcome for a 51-hour, single-GPU, 20-design campaign is **several
confirmed human binders, one or two with a measurable pH shift, and
cross-reactivity as the genuine long shot.** Design and document for that. No
design in this submission is a pH-selective cross-species binder until Adaptyv's
sensorgrams say so, and the methods report should state exactly that, with the
not-computed list in full.