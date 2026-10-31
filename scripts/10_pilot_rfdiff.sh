#!/usr/bin/env bash
# Block C pilot - RFdiffusion timing on the 18 A crop.
#
# Purpose: replace the brief's ESTIMATED 1.1 min/backbone with a measurement
# before block D commits 10.5 h of GPU. Section 5 of the brief is explicit that
# the numbers there are extrapolations, not measurements.
#
# 10 designs. At the estimated rate this is ~11 GPU-minutes, under the 30-min
# check-in threshold in CLAUDE.md.
#
# Hydra overrides ONLY - the image has an ENTRYPOINT. Container-absolute paths
# under /work. hydra.run.dir must be set or the run dies writing to root-owned
# /app.
set -euo pipefail
cd "$(dirname "$0")/.."

N="${1:-10}"
SET="${2:-shoulder}"
OUT="pilot/rfdiff_${SET}"
mkdir -p "$OUT" logs

CONTIG=$(cat target/crop_contigs.txt)
HOT=$(tr -d '\n' < "target/hotspots_${SET}.txt")

echo "=== RFdiffusion pilot ==="
echo "set      : ${SET}"
echo "hotspots : ${HOT}"
echo "contig   : ${CONTIG}"
echo "designs  : ${N}"
echo "start    : $(date -u '+%F %T UTC')"
echo

START=$(date +%s)

./rfdiffusion.sh \
  inference.output_prefix="/work/${OUT}/d" \
  inference.input_pdb=/work/target/tfr1_A_crop.pdb \
  "contigmap.contigs=${CONTIG}" \
  "ppi.hotspot_res=[${HOT}]" \
  inference.num_designs="${N}" \
  denoiser.noise_scale_ca=0 \
  denoiser.noise_scale_frame=0 \
  inference.write_trajectory=False \
  hydra.run.dir=/work/.hydra

END=$(date +%s)
ELAPSED=$((END-START))
MADE=$(ls -1 "${OUT}"/d_*.pdb 2>/dev/null | wc -l)

echo
echo "=== pilot result ==="
echo "elapsed      : ${ELAPSED}s"
echo "backbones    : ${MADE}"
if [ "${MADE}" -gt 0 ]; then
  echo "per backbone : $(awk -v e="${ELAPSED}" -v m="${MADE}" 'BEGIN{printf "%.2f min", e/m/60}')"
  echo "=> 10.5 h of GPU would buy $(awk -v e="${ELAPSED}" -v m="${MADE}" \
       'BEGIN{printf "%d", 10.5*3600/(e/m)}') backbones"
fi
echo "disk used by this run: $(du -sh "${OUT}" | cut -f1)"
echo "finish       : $(date -u '+%F %T UTC')"
