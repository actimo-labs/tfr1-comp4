# Rebased schedule — 45.5 h, not 51 h

Written at **Sat 31 Oct 14:38 UTC = Sun 1 Nov 01:38 AEDT**.
Deadline **Mon 2 Nov 22:59 AEDT = Mon 2 Nov 11:59 UTC**. Measured remaining at
time of writing: **45.35 h**.

The machine clock is **UTC**; Melbourne is **UTC+11 (AEDT)**. Every time below
is given in both. All `date` output in logs will be UTC — do not read it as
local.

§9 of the brief was written against 51 h from Sat 31 Oct ~20:00 AEDT. It has
lost **5.5 h**. Per the brief's own rule, the cut comes out of **generation
volume**, not validation stages.

## Verdict up front: the plan still fits, but not for the reason §9 assumes

The 5.5 h cut is survivable, and in fact cheaper than it looks, because
**§9 mis-sites the bottleneck**. §9 gives generation 16 h and the Boltz screen
10 h. On challenge-02's measured rates that is backwards for this target:

- RFdiffusion on a ~180-residue crop + 80 aa binder ≈ 260 total residues,
  ~1.1 min/backbone (estimate extrapolated from challenge 02, not yet measured
  here) → 16 h would buy ~870 backbones.
- Boltz-2 at challenge 02's measured 20.7 s/design would score only ~1700
  designs in 10 h — and that was against a **314-residue** target. TfR1's
  cropped target is comparable, but the **uncropped monomer is 639 aa**, so the
  validation-stage Boltz call is likely several times slower. At 60 s/design,
  10 h buys ~600 designs.

So Boltz, not RFdiffusion, sets the funnel width. Cutting backbones from ~870
to ~570 costs nothing we were going to use. **The 5.5 h comes off generation
and the funnel target stays where the brief wants it: ≥150 survivors into the
pH scan.**

The pilot in block C must replace both estimates with real numbers. If the
cropped-target Boltz call comes in under ~30 s/design, widen generation again
with the banked time.

## Blocks

`t` = hours from 14:38 UTC Sat.

| # | Block | Δh | UTC | AEDT |
|---|---|---|---|---|
| A | Audit, schedule, **stage 00 conservation map** | 0.5 | Sat 14:08–14:38 | Sun 01:08–01:38 |
| B | **EPITOPE DECISION GATE — Harjas** | 1.0 | Sat 14:38–15:38 | Sun 01:38–02:38 |
| C | Target prep (01), crop (02), **pilot 10 designs end-to-end** incl. Boltz timing crop-vs-monomer | 3.0 | Sat 15:38–18:38 | Sun 02:38–05:38 |
| D | **Generation** — RFdiffusion + ProteinMPNN, both pools | **10.5** | Sat 18:38 – Sun 05:08 | Sun 05:38–16:08 |
| E | Boltz human screen, calibration with VP1u positives, AF2-IG triage | 10.0 | Sun 05:08–15:08 | Sun 16:08 – Mon 02:08 |
| F | pH scan (CPU, overlaps E/G), His refold | 8.0 | Sun 15:08–23:08 | Mon 02:08–10:08 |
| G | Mouse screen on survivors, self-consistency gate | 5.0 | Sun 23:08 – Mon 04:08 | Mon 10:08–15:08 |
| H | Ranking, portfolio, submission package, methods report, pre-submission audit | 4.0 | Mon 04:08–08:08 | Mon 15:08–19:08 |
| I | **Buffer** | 3.5 | Mon 08:08–11:38 | Mon 19:08–22:38 |
| | **Total** | **45.5** | | deadline Mon 11:59 UTC / 22:59 AEDT |

Change against §9: **generation 16 h → 10.5 h**. Blocks E, F, G, H are
unchanged. Buffer grows 2 h → 3.5 h, because block A came in at 0.5 h against
the 2 h §9 allowed.

## What is already done

Block A is complete and came in **1.5 h under budget**. Stage 00 ran, the
mouse sequence was recovered from the challenge page (not substituted from
UniProt), and the epitope question is answered. That surplus is already banked
into the buffer above.

## Standing risks on this schedule

1. **Block B is the only blocking human dependency.** It is budgeted at 1 h.
   Every hour the gate stays open comes straight out of block I. If the gate is
   still open at Sun 04:38 AEDT (t+3), the fallback is the pre-committed rule
   in §4 of the brief — which stage 00 has already resolved to a single
   recommendation — and generation starts against it.
2. **Disk, not GPU, is the hard resource.** 23 GB free on a filesystem shared
   with another user who holds 1.6 TB. See the audit. Challenge 02's Boltz
   output alone was 4.0 GB. Block E must write only scoring summaries plus the
   top-N structures, and prune as it goes, or it will fill the disk around
   Sun 20:00 AEDT with no warning.
3. **Boltz throughput on a 639-residue target is unmeasured.** If the pilot
   shows >90 s/design on the uncropped monomer, validation moves to the crop
   and the uncropped monomer is used only on the final slate. That decision
   belongs in block C, not block E.
4. **Mouse scoring is a second-stage filter**, not a parallel arm. If it is run
   on everything, block G overruns by roughly 4×.
5. **No FreeBindCraft.** §6's exception requires time already banked. 5.5 h
   short is not banked time.

## Open items that need a browser or Slack, not compute

These do not block generation but must close before block H:

- **O1** — the "four conditions" discrepancy (mouse pH 6.0?). Ask in Slack now
  so the answer arrives before ranking.
- **O2** — full FAQ via *Copy Markdown*; confirm submission field names.
- **O3** — track and cap (20 vs 40 designs). Plan for 20.
- **O4** — the published selection prompt, before writing dossiers.
- **O5** — resolved by stage 00. The 318/633/719 seam lead is **discarded**;
  see the conservation report. No design will be built on it, so the preprint
  no longer needs verifying.
