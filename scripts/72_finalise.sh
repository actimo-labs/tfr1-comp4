#!/usr/bin/env bash
# Final chain: seed2 -> mouse-score the controls -> pose agreement -> rank ->
# package. Every stage here has been pre-flighted on real data.
set -uo pipefail
cd "$(dirname "$0")/.."
B=/home/HarjasG/miniforge3/envs/boltz/bin/python
P=/home/HarjasG/tnf-track3/venv-phscan/bin/python
note(){ echo "[$(date -u '+%F %T')] $*" | tee -a logs/72_finalise.log; }

note "waiting for seed 2"
while ! [ -f boltz/selfconsist/scores.json ]; do
  if ! tmux has-session -t sc 2>/dev/null && ! [ -f boltz/selfconsist/scores.json ]; then
    note "seed2 session gone before completion - restarting"
    tmux new-session -d -s sc "bash -lc 'cd ~/tnf-track3/tfr1-comp4 && $B scripts/04_boltz.py --extra-seqs gen/selfconsist_input.fasta --tag selfconsist --target-mode full --species human --seed 2024 2>&1 | tee -a logs/06c_seed2.log'"
  fi
  sleep 30
done
note "seed 2 complete"

note "mouse-scoring the 5 unmutated controls"
$B scripts/04_boltz.py --extra-seqs gen/mouse_controls.fasta --tag mouse \
    --target-mode full --species mouse >>logs/72_mouse_ctrl.log 2>&1 \
    && note "control mouse scoring ok" || note "control mouse scoring FAILED"

note "two-seed pose agreement"
$P -I scripts/06c_selfconsistency.py >>logs/72_selfconsist.log 2>&1 \
    && note "pose agreement ok" || note "pose agreement FAILED"

note "ranking"
$P -I scripts/80_rank_and_package.py >>logs/72_rank.log 2>&1 \
    && note "ranking ok" || note "ranking FAILED"

note "packaging"
$P -I scripts/81_build_package.py >>logs/72_package.log 2>&1 \
    && note "packaging ok" || note "packaging FAILED"

note "=== FINALISE COMPLETE ==="
echo FINALISE_DONE
