#!/usr/bin/env bash
# Block D - main RFdiffusion generation across the portfolio.
#
# Budget: 10.5 h of the 45.5 h (see runlog/00-schedule.md). The brief's §9 gave
# generation 16 h, but Boltz is the real throttle, not RFdiffusion, so the 5.5 h
# schedule cut lands here and costs nothing we would have used.
#
# RESUME, not restart. Each batch counts the backbones already on disk and asks
# RFdiffusion only for the shortfall. Re-running this script after a crash, an
# OOM or a killed tmux picks up where it stopped. There is no budget for a lost
# run.
#
# TRAJECTORIES ARE OFF. inference.write_trajectory=False. With them on,
# RFdiffusion writes ~9.1 MB of trajectory per design against ~0.13 MB of
# actual output -- challenge 02's 4.3 GB rfdiff/ directory is almost entirely
# trajectory files nobody read. At 540 designs that is 4.9 GB saved on a
# filesystem with 23 GB free.
#
# Batches vary BOTH epitope and binder length. Length variation is the cheapest
# source of fold diversity, and §7 of the brief makes diversity a hard
# requirement: challenge 02 shipped 6 structural clusters across 20 designs and
# called it its weakest axis. Twenty variations on one scaffold is one bet.
set -euo pipefail
cd "$(dirname "$0")/.."

OUTROOT="gen/rfdiff"
mkdir -p "$OUTROOT" logs

# batch name | hotspot set | binder length | target count
#
# shoulder   = the primary epitope: conserved switch anchor K374 plus four of
#              the nine VP1u footprint residues. Classes A, B and D.
# apex       = R208/K371 tip. Class C affinity hedge; mouse expected to fail.
# swap_locus = the 348/371 cation-swap site. Class E, speculative.
BATCHES=(
  "shoulder_s|shoulder|55-75|110"
  "shoulder_m|shoulder|70-95|140"
  "shoulder_l|shoulder|90-120|90"
  "apex_m|apex|70-95|90"
  "apex_l|apex|90-120|45"
  "swap_m|swap_locus|70-95|65"
)

# Free-space guard. The filesystem is shared and another user holds 1.6 TB of
# the 2.0 TB; our 23 GB of headroom is not reliably ours. Halt cleanly between
# batches rather than dying mid-write with a half-written PDB. Everything here
# resumes, so stopping is cheap and corruption is not.
MIN_FREE_GB="${MIN_FREE_GB:-5}"
check_disk() {
  local free_gb
  free_gb=$(df -BG --output=avail / | tail -1 | tr -dc '0-9')
  if [ "${free_gb}" -lt "${MIN_FREE_GB}" ]; then
    echo "HALT: only ${free_gb} GB free, below the ${MIN_FREE_GB} GB floor."
    echo "Nothing is corrupted. Free space and re-run this script; it resumes."
    exit 3
  fi
}
check_disk

CONTIG=$(cat target/crop_contigs.txt)
# strip the trailing binder-length token; each batch supplies its own
CONTIG_BASE=$(echo "$CONTIG" | sed -E 's/ [0-9]+-[0-9]+\]$//' | sed 's/^\[//')

TOTAL_START=$(date +%s)
echo "=== block D generation started $(date -u '+%F %T UTC') ==="
echo "contig base: ${CONTIG_BASE}"
echo

for spec in "${BATCHES[@]}"; do
  IFS='|' read -r NAME SET LEN WANT <<< "$spec"
  OUT="${OUTROOT}/${NAME}"
  mkdir -p "$OUT"

  HAVE=$(ls -1 "${OUT}"/d_*.pdb 2>/dev/null | wc -l)
  NEED=$(( WANT - HAVE ))
  if [ "$NEED" -le 0 ]; then
    echo "[${NAME}] complete: ${HAVE}/${WANT} - skipping"
    continue
  fi

  check_disk
  HOT=$(tr -d '\n' < "target/hotspots_${SET}.txt")
  echo "[${NAME}] have ${HAVE}/${WANT}, generating ${NEED}"
  echo "         set=${SET} len=${LEN} hotspots=${HOT}"

  B_START=$(date +%s)
  # RFdiffusion numbers outputs from 0 and would overwrite on resume, so offset
  # the prefix by what already exists.
  ./rfdiffusion.sh \
    inference.output_prefix="/work/${OUT}/d" \
    inference.input_pdb=/work/target/tfr1_A_crop.pdb \
    "contigmap.contigs=[${CONTIG_BASE} ${LEN}]" \
    "ppi.hotspot_res=[${HOT}]" \
    inference.num_designs="${WANT}" \
    inference.design_startnum="${HAVE}" \
    denoiser.noise_scale_ca=0 \
    denoiser.noise_scale_frame=0 \
    inference.write_trajectory=False \
    hydra.run.dir=/work/.hydra 2>&1 | grep -vE "Timestep|FALLBACK|UserWarning|^ *sh = |codegen|PYTORCH" || true
  B_END=$(date +%s)

  NOW=$(ls -1 "${OUT}"/d_*.pdb 2>/dev/null | wc -l)
  echo "[${NAME}] done: ${NOW}/${WANT} in $(( B_END - B_START ))s"
  echo "         disk: $(du -sh "${OUT}" | cut -f1)   free: $(df -h / | awk 'NR==2{print $4}')"
  echo
done

echo "=== block D generation finished $(date -u '+%F %T UTC') ==="
echo "total elapsed: $(( $(date +%s) - TOTAL_START ))s"
echo "backbones: $(ls -1 ${OUTROOT}/*/d_*.pdb 2>/dev/null | wc -l)"
echo "disk used: $(du -sh ${OUTROOT} | cut -f1)"
echo "free space: $(df -h / | awk 'NR==2{print $4}')"
