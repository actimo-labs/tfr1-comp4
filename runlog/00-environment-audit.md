# Stage 00 environment audit

Run **Sat 31 Oct 14:23–14:38 UTC** (Sun 1 Nov 01:23–01:38 AEDT). Every path
below was verified on disk, not taken from the brief. Nothing was downloaded
to fix a gap; gaps are flagged.

## Verdict

**The environment is ready. Nothing is missing.** Two things in the brief are
wrong and one resource is genuinely scarce:

- The challenge-02 scripts are at `tnf-a-comp2/tnf-upload/scripts/`, not
  `tnf-upload/scripts/`.
- `dl_binder_design` **is** installed (challenge-02's CLAUDE.md says it is not —
  that is stale; challenge 03 installed it on 25 Oct). The AF2 weights the
  AF2-IG port needs are inside it.
- **Disk: 23 GB free.** This is the binding constraint on the campaign, not GPU.

## Host

| Item | Value |
|---|---|
| Kernel | Linux 5.10.0-46-cloud-amd64 |
| Clock | **UTC** (Melbourne = UTC+11 AEDT) |
| CPU | 12 cores |
| RAM | 167 GB total, 164 GB available |
| tmux | 3.1c, session `tfr1` active |
| Docker | working, non-root OK |

## GPU

```
NVIDIA A100-SXM4-80GB   driver 550.90.07   CUDA 12.4
1 MiB / 81920 MiB used,  0% util,  31 C,  no running processes
```

Idle and fully available. **Driver CUDA is 12.4** — see the torch note below,
this was a documented challenge-02 failure.

## Disk — the real constraint

```
/dev/sda1   2.0T   1.9T used   23G avail   99%
```

The filesystem is **shared with other users**:

| Path | Size |
|---|---|
| `/home/jaberjaber` | **1.6 T** |
| `/home/harjas50` | 116 G |
| `/home/HarjasG` | 55 G |
| `/var/lib/docker` | 12 G |

**1.6 TB belongs to another user. It is not ours and was not touched.** We
cannot free it, and if that user grows, we run out mid-campaign with no notice.

What we hold, for reference if Harjas decides to prune (nothing was deleted —
challenge-02/03 artifacts are read-only per the working agreement):

| Path | Size | Note |
|---|---|---|
| `~/miniforge3` | 23 G | conda envs, all four in use |
| `~/tnf-track3/rfdiff` | 4.3 G | challenge-02 backbone dump |
| `~/tnf-track3/boltz` | 4.0 G | challenge-02 Boltz output |
| `~/tnf-track3/glp1r-comp3` | 3.0 G | challenge-03 tree |
| `~/.boltz` | 7.6 G | Boltz-2 checkpoints + `mols.tar` (1.9 G, extractable, the extracted `mols/` is already present) |

`~/.boltz/mols.tar` is 1.9 GB and its contents are already extracted into
`~/.boltz/mols/`. That is the one obviously redundant file in our own space.
**Not deleted — flagging it for a decision.**

## Repository

| Item | Value |
|---|---|
| Toplevel | `/home/HarjasG/tnf-track3/tfr1-comp4` |
| Branch | `main` |
| Remote | `origin https://github.com/actimo-labs/tfr1-comp4.git` — **never pushed to** |
| HEAD | `3d2a1c8 spec brief` |
| Status at audit | clean |

### Git identity — MISMATCH, first commit withheld

```
git config user.name     harjas2001                 OK
git config user.email    harjas.g2001@gmail.com     NOT VERIFIABLE — see below
```

The repo-local email is an **outlier against every other piece of evidence on
this machine**:

| Source | Email |
|---|---|
| `~/.gitconfig` (global) | `116927987+harjas2001@users.noreply.github.com` |
| both existing commits in this repo | `116927987+harjas2001@users.noreply.github.com` |
| `glp1r-comp3` local config + all commits | `116927987+harjas2001@users.noreply.github.com` |
| `tnf-a-comp2` local config + all commits | `116927987+harjas2001@users.noreply.github.com` |
| **`tfr1-comp4` local config** | **`harjas.g2001@gmail.com`** |

The `116927987+harjas2001@users.noreply.github.com` form is GitHub's own
per-account noreply address — it is attributable to the account by
construction. The gmail address may well be verified on the account, but **I
cannot check another user's verified addresses**, and the account is evidently
configured for private-email mode, which also means GitHub will *reject a push*
of commits authored from a real address if push-privacy protection is on.

Per the brief: reported, **not changed by me**, and **no commit made**. See the
handover note at the end of this file.

## Tool checkouts

| Tool | Path | Commit | State |
|---|---|---|---|
| RFdiffusion | `~/tnf-track3/RFdiffusion` | `86507b6` (15 Jul 2026) | present |
| ProteinMPNN | `~/tnf-track3/ProteinMPNN` | `8907e66` (27 Jun 2023) | present |
| dl_binder_design | `~/tnf-track3/dl_binder_design` | — | present, 382 M |

### Weights on disk

| Weight set | Path | Size |
|---|---|---|
| RFdiffusion `Base_ckpt.pt` | `RFdiffusion/models/` | 484 M |
| RFdiffusion `Complex_base_ckpt.pt` | `RFdiffusion/models/` | 484 M |
| RFdiffusion `Complex_beta_ckpt.pt` | `RFdiffusion/models/` | 483 M |
| ProteinMPNN soluble | `ProteinMPNN/soluble_model_weights/` | `v_48_002/010/020/030.pt` |
| ProteinMPNN vanilla | `ProteinMPNN/vanilla_model_weights/` | `v_48_002/010/020/030.pt` |
| ProteinMPNN Cα | `ProteinMPNN/ca_model_weights/` | `v_48_002/010/020.pt` |
| Boltz-2 confidence | `~/.boltz/boltz2_conf.ckpt` | 2.3 G |
| Boltz-2 affinity | `~/.boltz/boltz2_aff.ckpt` | 2.1 G |
| **AF2 `model_1_ptm`** | `dl_binder_design/af2_initial_guess/model_weights/params/params_model_1_ptm.npz` | 373 M |

The AF2-IG port hardcodes
`DBD = ~/tnf-track3/dl_binder_design/af2_initial_guess` and reads
`model_weights/params` under it — **that path resolves correctly**. No AF2
download needed. (`~/tnf-track3/af2/` is an empty directory and is a red
herring.)

RFdiffusion Docker image: `rosettacommons/rfdiffusion:latest`,
`33d504109509`, 12.2 GB, 16 months old.

## Environments

| Env | Interpreter | Key versions | GPU |
|---|---|---|---|
| `venv-phscan` | `~/tnf-track3/venv-phscan/bin/python` (3.10) | biopython **1.88**, propka 3, numpy 2.2.6 | CPU only — correct |
| `boltz` | `~/miniforge3/envs/boltz/bin/python` | **torch 2.14.1+cu126** | `cuda.is_available() == True` |
| `mpnn` | `~/miniforge3/envs/mpnn/bin/python` | torch 2.5.1+cu121 | True |
| `af2ig` | `~/miniforge3/envs/af2ig/bin/python` | **jax 0.4.20** | `[cuda(id=0)]` |
| `viz` | `~/miniforge3/envs/viz` | — | not exercised |

**On the challenge-02 torch trap:** that failure was `torch 2.14.1+cu130` on a
CUDA 12.4 driver. The env now carries **`+cu126`**, and
`torch.cuda.is_available()` returns `True` against the 12.4 driver (minor-version
compatibility holds). **Already fixed — do not re-fix it.**

**biopython is 1.88**, where `Bio.pairwise2` is gone. Stage 00 uses
`Bio.Align.PairwiseAligner`. Any ported script still importing `pairwise2` will
fail — check before running, not after.

## Network

| Endpoint | Result |
|---|---|
| `files.rcsb.org/download/11QC.cif` | HTTP 200, 1 044 962 B |
| `files.rcsb.org/download/11RN.cif` | HTTP 200, 1 142 073 B |
| `rest.uniprot.org/.../P02786.fasta` | HTTP 200 |
| `rest.uniprot.org/.../Q62351.fasta` | HTTP 200 |
| challenge page (plain `curl`, browser UA) | HTTP 200, 260 834 B — **carries both sequences, see below** |

Biopython `PDBList` was not used (retired endpoint, 403s — brief §10).

## Port sources — verified present

Challenge 02, at `~/tnf-track3/tnf-a-comp2/tnf-upload/scripts/`
(**not** `tnf-upload/scripts/` as the brief states):

`01_prep_target.py` `02_crop_target.py` `03_mpnn.py` `04_boltz.py`
`04b_make_controls.py` `04c_calibrate.py` `05_ph_scan.py` `06_his_refold.py`
`07_submit.py` `make_template_cif.py` `render_complex.py`

Challenge 03, at `~/tnf-track3/glp1r-comp3/scripts/` — all three named modules
present: `af2ig/` (`predict_pdb.py`, `af2_util_nopyro.py`),
`42_selfconsistency.py`, `20_denovo_check.py`.

`rfdiffusion.sh`: challenge 02's copy is byte-identical to the one at
`~/tnf-track3/rfdiffusion.sh` (`4ac104e8…`); challenge 03's differs
(`cb54e543…`) and is the newer one. Diff them before adopting either.

**Both sibling trees were read only. Nothing in them was modified, deleted or
committed to.**

## Target files

| File | State |
|---|---|
| `target/human_tfr1_121_760.fasta` | shipped with brief; **640 aa verified**, spans 121–760 |
| `target/mouse_tfr1_123_763.fasta` | **newly created this session — 641 aa, recovered from the challenge page, NOT UniProt** |
| `target/structures/11QC.cif` | downloaded, 10 198 atom records, chains A+B, 639 res/chain, 6 NAG |
| `target/structures/11RN.cif` | downloaded, 11 200 atom records, TfR1 A+B, VP1u **X+Y** (not C/D as the brief says) |

### The mouse sequence was retrieved from the page, not substituted

The brief says the Mouse TfR1 tab is JS-gated and did not survive a plain
fetch. That is true of the *rendered* text — a markdown-converting fetch
returns the human tab only. But the raw HTML payload (260 KB) carries **both**
sequences; the mouse one is server-rendered into the document and only hidden
by the tab UI. Extracted with a plain regex for runs of ≥200 residue letters:

- 640 aa — human, byte-identical to the shipped human FASTA
- **641 aa — mouse**, `RLYWADLKTLLSEKLNSIEFADTIKQLSQNTY…TIQGVANALSGDIWNIDNEF`

641 aa is exactly the stated 123–763 construct span. It is saved to
`target/mouse_tfr1_123_763.fasta` and is **cross-validated against six
[OFFICIAL] page-stated human→mouse differences** — all six reproduce exactly
(see the conservation report). **UniProt Q62351 was never used.**

## Stage 00 artifacts produced

| File | Content |
|---|---|
| `scripts/00_conservation.py` | the conservation map |
| `scripts/00b_patch_scan.py` | systematic surface-patch ranking |
| `target/conservation_table.tsv` | 639 residues × 16 columns |
| `target/conservation_candidates.json` | ranked conserved-cation patches |
| `target/patch_scan.json` | all 244 patches + switch-cation verdicts |
| `target/conservation_report.md` | **the decision document** |
| `logs/00_conservation.log`, `logs/00b_patch_scan.log` | full run output |

## Two bugs found and fixed in my own stage-00 code

Recording these because both are the silent kind:

1. `round(x or -1, 1)` — a genuine **0.0 Å distance is falsy**, so every
   residue *inside* the transferrin site was reported as `d_tf = -1`
   (i.e. "missing") instead of 0. Replaced with an explicit `is None` test.
2. The same idiom in the patch-footprint filter (`or 99`) **excluded each
   cation from its own patch**, undercounting every footprint by one residue
   and skewing the conservation fractions. Replaced with a precomputed
   distance matrix.

Both were fixed before the numbers in the conservation report were generated.

## Handover — what Harjas needs to do

The audit and schedule are **written but uncommitted**, because the identity
check failed. One of:

```bash
# if the gmail IS verified on the GitHub account, keep it and commit:
cd ~/tnf-track3/tfr1-comp4 && git add -A && git commit

# if it is not — match the other two challenge repos:
git -C ~/tnf-track3/tfr1-comp4 config user.email \
    116927987+harjas2001@users.noreply.github.com
```

Then the stage-00 commit can go in. I have not run either.
