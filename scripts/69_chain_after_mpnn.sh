#!/usr/bin/env bash
# Wait for the MPNN tmux session to exit, then run the autopilot immediately.
# This is the structural fix for the 3h37m idle gap after block D: the handoff
# lives in the shell, not in anyone remembering to start the next stage.
set -uo pipefail
cd "$(dirname "$0")/.."
echo "[$(date -u '+%F %T')] waiting for mpnn session to exit..."
while tmux has-session -t mpnn 2>/dev/null; do sleep 30; done
echo "[$(date -u '+%F %T')] mpnn session gone; plain=$(find gen/mpnn/plain/seqs -name '*.fa' 2>/dev/null | wc -l) hisb=$(find gen/mpnn/hisb/seqs -name '*.fa' 2>/dev/null | wc -l)"
exec ./scripts/70_autopilot.sh
