# Open items — status

## O1 — the "four conditions" discrepancy. RESOLVED AS DOCUMENTATION, 31 Oct.

The ranking intro claims "four conditions" and then names three:
human TfR1 pH 7.4, human pH 6.0, mouse pH 7.4.

Every other mention of mouse on the page is pH 7.4 only:
- Objective 2: "we will test binding against mouse TfR1 at pH 7.4"
- Objective 3: "we will measure binding against human TfR1 in both in vitro
  conditions (pH 7.4 and 6.0)"

No sentence anywhere states mouse is tested at pH 6.0. No FAQ answer touches
assay conditions.

**Why it does not affect design.** Every ranking criterion is explicitly scoped
to a named condition:

| # | Criterion | Conditions consumed |
|---|---|---|
| 1 | pH-selective | human 7.4 vs human 6.0 |
| 2 | cross-reactivity | mouse 7.4 vs human 7.4 |
| 3 | affinity | human 7.4 |

Even if mouse pH 6.0 is measured, **no stated criterion consumes it** — it
would be unranked data. The design objective is identical under both readings.

**Decision: proceed on three ranked conditions.** Carry the discrepancy into
METHODS.md as an unresolved ambiguity, with the quotes above. Harjas is asking
in Slack in parallel; an answer changes documentation, not design.

Note this *does* still favour the strategy-2 (binder-core) switch slightly, for
the reason the brief gives: a species-independent switch would also switch in
mouse if mouse pH 6.0 ever were measured. That is already hedged in the
portfolio and needs no change.

## O3 — track and cap. RESOLVED, 31 Oct.

**Cap is 20 designs.** Confirmed by Harjas. The §7 portfolio stands as a
20-slot allocation; no doubling.

## O2 — full FAQ via *Copy Markdown*. OPEN, non-blocking.
Needed before the submission CSV is built, to confirm field names and
`molecule_class` values. Does not block generation or scoring.

## O4 — the published selection prompt. OPEN, non-blocking.
Needed before writing the per-design dossiers. Tab-gated, needs a browser.

## O5 — the 318/633/719 preprint. CLOSED by stage 00.
The conservation map discarded the seam lead on all three residues (H318 is
1.3 Å from the N317 sequon; see target/conservation_report.md §7). Nothing in the
campaign depends on it, so the preprint needs no verification.
