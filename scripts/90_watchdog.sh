#!/usr/bin/env bash
# Watchdog for the final stages. Twice this run a stage finished and the next
# one failed to start, leaving the A100 idle for hours (3h37m after block D,
# ~9h after the refold parse crashed). Both were silent. This detects the
# pattern - work outstanding, no tmux session, GPU idle - and restarts it.
cd "$(dirname "$0")/.."
B=/home/HarjasG/miniforge3/envs/boltz/bin/python
want=$(grep -c '^>' gen/selfconsist_input.fasta)
while true; do
  have=$(find boltz/selfconsist/pred -name 'confidence_*' 2>/dev/null | wc -l)
  if [ "$have" -ge "$want" ] || [ -f boltz/selfconsist/scores.json ]; then
    echo "WATCHDOG: seed2 complete ($have/$want)"; exit 0
  fi
  if ! tmux has-session -t sc 2>/dev/null; then
    echo "WATCHDOG: seed2 session gone at $have/$want - RESTARTING"
    tmux new-session -d -s sc "bash -lc 'cd ~/tnf-track3/tfr1-comp4 && $B scripts/04_boltz.py --extra-seqs gen/selfconsist_input.fasta --tag selfconsist --target-mode full --species human --seed 2024 2>&1 | tee -a logs/06c_seed2.log; echo SEED2_DONE_\$?'"
  fi
  sleep 60
done
