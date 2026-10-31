#!/usr/bin/env bash
# RFdiffusion via the maintained Docker image -- GLP1R (challenge 03) wrapper.
#
# Derived from ~/tnf-track3/rfdiffusion.sh (challenge 02), which is left
# untouched per the brief. ONE deliberate difference, see WEIGHTS below.
#
# Usage: pass ONLY Hydra overrides. The image's ENTRYPOINT already invokes
# run_inference.py, so a leading `python` or script path makes Hydra try to
# parse it as a config key. Do not override the container workdir either; the
# entrypoint resolves its script path relative to /app/RFdiffusion.
#
# $PWD is mounted at /work -- pass container-absolute paths (/work/target/...)
# and include hydra.run.dir=/work/.hydra.
#
# Mounts that must be writable because we run as the host UID while /app is
# root-owned in the image:
#   /app/RFdiffusion/schedules - IGSO3 diffusion schedule cache, mkdir'd on
#                                first run (model_runners.py:129). Without it
#                                the run dies with PermissionError before any
#                                sampling. Persisted across runs; saves several
#                                minutes of "Calculating IGSO3" each time.
#
# WEIGHTS -- why this differs from the challenge-02 wrapper
# --------------------------------------------------------
# The comp-02 wrapper bind-mounts the host dir
#   $HOME/tnf-track3/RFdiffusion/models  ->  /app/RFdiffusion/models:ro
# That host dir holds only 3 checkpoints (Base, Complex_base, Complex_beta),
# whereas the image ships 7 -- including ActiveSite_ckpt.pt and
# Complex_Fold_base_ckpt.pt. The mount therefore MASKS the image's own weights.
#
# GLP1R needs the masked ones:
#   ActiveSite_ckpt.pt        - small-motif scaffolding, for the class-B
#                               N-terminal message module
#   Complex_Fold_base_ckpt.pt - fold-conditioned binder design, for the
#                               "classes A and B span >=3 backbone topologies"
#                               requirement in brief section 4
# So we drop that mount and use the image's complete set. Costs no disk and no
# download. Verified present: see runlog/00-environment-audit.md finding F1.
set -euo pipefail

mkdir -p "$PWD/schedules"

exec docker run --rm --gpus all \
  --user "$(id -u):$(id -g)" \
  -v "$PWD:/work" \
  -v "$PWD/schedules:/app/RFdiffusion/schedules" \
  rosettacommons/rfdiffusion:latest "$@"
