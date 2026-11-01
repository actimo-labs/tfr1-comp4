#!/usr/bin/env bash
# Autopilot - chains the GPU stages so the card never sits idle between them.
#
# Written after block D ended at 22:32 UTC and the next stage did not start
# until 02:10: 3h37m of idle A100 because generation ran in a bare tmux session
# with nothing watching for its exit. This script removes the gap structurally
# rather than relying on someone noticing.
#
# Properties that matter:
#   RESUMABLE   every stage skips work already on disk, so re-running after a
#               crash costs only the unfinished part.
#   GUARDED     a failing stage is logged and marked, and the chain continues
#               to stages that do not depend on it rather than dying silently.
#   OBSERVABLE  one line per stage transition into logs/70_autopilot.log, and
#               a machine-readable gen/autopilot_state.json.
#   DISK-SAFE   halts cleanly above the free-space floor instead of dying
#               mid-write.
#
# Usage:  ./scripts/70_autopilot.sh [first_stage]
set -uo pipefail          # NB: deliberately no -e; stage guards handle failure
cd "$(dirname "$0")/.."

STATE="gen/autopilot_state.json"
MIN_FREE_GB="${MIN_FREE_GB:-5}"
BOLTZ_PY="/home/HarjasG/miniforge3/envs/boltz/bin/python"
MPNN_PY="/home/HarjasG/miniforge3/envs/mpnn/bin/python"
PH_PY="/home/HarjasG/tnf-track3/venv-phscan/bin/python"
mkdir -p gen logs

note() { echo "[$(date -u '+%F %T')] $*" | tee -a logs/70_autopilot.log; }

check_disk() {
  local g; g=$(df -BG --output=avail / | tail -1 | tr -dc '0-9')
  if [ "${g}" -lt "${MIN_FREE_GB}" ]; then
    note "HALT: ${g} GB free, below ${MIN_FREE_GB} GB floor. Nothing corrupted; free space and re-run."
    exit 3
  fi
}

declare -A STATUS
mark() { STATUS[$1]=$2; python3 - "$STATE" "$1" "$2" <<'PY' 2>/dev/null || true
import json,sys,os
p,k,v=sys.argv[1],sys.argv[2],sys.argv[3]
d=json.load(open(p)) if os.path.exists(p) else {}
d[k]=v; json.dump(d,open(p,"w"),indent=1)
PY
}

stage() {                      # stage <name> <command...>
  local name="$1"; shift
  check_disk
  note "START ${name}"
  local t0; t0=$(date +%s)
  if "$@" >>"logs/${name}.log" 2>&1; then
    note "OK    ${name}  ($(( $(date +%s) - t0 ))s)"
    mark "$name" ok
    return 0
  else
    note "FAIL  ${name}  ($(( $(date +%s) - t0 ))s) - see logs/${name}.log"
    mark "$name" failed
    return 1
  fi
}

note "=== autopilot start ==="

# ---- 35: choose which designs get a Boltz call --------------------------
if [ ! -s gen/selected.fasta ]; then
  stage 35_select "$PH_PY" -I scripts/35_select_for_boltz.py --cap 900
else
  note "SKIP  35_select (gen/selected.fasta exists)"
fi

if [ ! -s gen/selected.fasta ]; then
  note "ABORT: no gen/selected.fasta - nothing downstream can run"
  exit 1
fi
note "selected $(grep -c '^>' gen/selected.fasta) designs for the human screen"

# ---- 04: Boltz human screen, FULL 639-residue monomer -------------------
# full monomer not the apical crop. The crop is 1.4x faster but hides
# off-target docking: the pilot's top design by ipTM (0.775) binds the
# transferrin site and a glycan sequon and hits 0/5 hotspots, and the crop
# scored it low only because the crop lacks the site it wanted. 04e is the
# real filter. See runlog/00-schedule.md and its CORRECTION section.
stage 04_boltz_human "$BOLTZ_PY" scripts/04_boltz.py \
    --extra-seqs gen/selected.fasta \
    --tag screen --target-mode full --species human

# ---- 04e: epitope fidelity gate (HARD) ---------------------------------
# Must run before anything downstream spends time on a design. The pilot's
# top-scoring design by ipTM binds the transferrin site and a glycan sequon
# and hits zero targeted hotspots; confidence alone cannot catch that.
stage 04e_epitope_qc "$PH_PY" -I scripts/04e_epitope_qc.py \
    --pred boltz/screen/pred --out gen/epitope_qc.tsv

# ---- 04b/04c: controls and calibration ----------------------------------
# Positives matter this round. 11RN's VP1u RBD is a real apical binder and
# anchors the TOP of the scale, which challenge 02 never had - it calibrated
# against shuffled negatives only and had no idea what "good" looked like.
stage 04b_controls "$PH_PY" -I scripts/04b_make_controls.py || true
stage 04c_calibrate "$PH_PY" -I scripts/04c_calibrate.py || true

note "=== autopilot: GPU screen chain complete ==="
note "state: $(cat $STATE 2>/dev/null)"
echo AUTOPILOT_A_DONE
