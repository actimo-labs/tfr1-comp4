#!/usr/bin/env bash
# Round 2 generation - K374-CENTRIC.
#
# Round 1 measured the problem precisely. Across the 123 backbones that docked
# on-epitope, the closest binder CB sits at median 5.1 A from Y211 and L209 but
# median 11.6 A from K374:
#
#     target   median closest binder CB    within 10 A
#     209           5.1 A                   123/123
#     211           5.1 A                   123/123
#     374          11.6 A                    22/123
#     371          13.1 A                    11/123
#
# RFdiffusion centred every binder on the two most exposed hotspots and left the
# conserved cation at the rim. That is why only 11 of 153 switches are
# cross-species and ten of those are a single lucky histidine: a His cannot
# reach a cation 11.6 A away.
#
# Objectives 1 and 2 are the top two ranked criteria and ONLY a conserved-cation
# switch serves both, so round 2 drops 209/211 from the hotspot set entirely and
# aims at the K374/K371/N348 cluster. Same crop (it already contains all three);
# only ppi.hotspot_res changes.
set -uo pipefail
cd "$(dirname "$0")/.."
OUTROOT="gen/rfdiff_r2"; mkdir -p "$OUTROOT" logs

MIN_FREE_GB="${MIN_FREE_GB:-5}"
check_disk(){ local g; g=$(df -BG --output=avail / | tail -1 | tr -dc '0-9')
  [ "$g" -lt "$MIN_FREE_GB" ] && { echo "HALT: ${g} GB free"; exit 3; }; return 0; }

CONTIG=$(cat target/crop_contigs.txt)
BASE=$(echo "$CONTIG" | sed -E 's/ [0-9]+-[0-9]+\]$//' | sed 's/^\[//')
HOT=$(tr -d '\n' < target/hotspots_k374.txt)

BATCHES=("k374_m|70-95|100" "k374_l|90-120|70")

echo "=== round 2 (K374-centric) started $(date -u '+%F %T UTC') ==="
echo "hotspots: ${HOT}"
for spec in "${BATCHES[@]}"; do
  IFS='|' read -r NAME LEN WANT <<< "$spec"
  OUT="${OUTROOT}/${NAME}"; mkdir -p "$OUT"; check_disk
  HAVE=$(find "$OUT" -maxdepth 1 -name 'd_*.pdb' 2>/dev/null | wc -l)
  [ "$HAVE" -ge "$WANT" ] && { echo "[${NAME}] complete ${HAVE}/${WANT}"; continue; }
  echo "[${NAME}] have ${HAVE}/${WANT} len=${LEN}"
  T0=$(date +%s)
  ./rfdiffusion.sh \
    inference.output_prefix="/work/${OUT}/d" \
    inference.input_pdb=/work/target/tfr1_A_crop.pdb \
    "contigmap.contigs=[${BASE} ${LEN}]" \
    "ppi.hotspot_res=[${HOT}]" \
    inference.num_designs="${WANT}" \
    inference.design_startnum="${HAVE}" \
    denoiser.noise_scale_ca=0 denoiser.noise_scale_frame=0 \
    inference.write_trajectory=False \
    hydra.run.dir=/work/.hydra 2>&1 | grep -vE "Timestep|FALLBACK|UserWarning|^ *sh = |codegen|PYTORCH" || true
  NOW=$(find "$OUT" -maxdepth 1 -name 'd_*.pdb' 2>/dev/null | wc -l)
  echo "[${NAME}] done ${NOW}/${WANT} in $(( $(date +%s) - T0 ))s  free $(df -h / | awk 'NR==2{print $4}')"
done
echo "=== round 2 generation finished $(date -u '+%F %T UTC') ==="
echo "backbones: $(find ${OUTROOT} -name 'd_*.pdb' | wc -l)"
echo R2GEN_DONE
