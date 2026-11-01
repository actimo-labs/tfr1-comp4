#!/usr/bin/env bash
# Overnight chain, round 2. ONLY stages already proven on real data in this
# campaign: MPNN -> select -> Boltz -> epitope QC -> pH scan -> switch policy.
#
# His refold, the mouse screen and the self-consistency gate are deliberately
# NOT in here. They have not yet been exercised on this target's data, and an
# untested stage failing silently at 03:00 costs more than running it attended.
# Waiters are armed so they get driven as soon as this chain lands.
set -uo pipefail
cd "$(dirname "$0")/.."
BOLTZ_PY=/home/HarjasG/miniforge3/envs/boltz/bin/python
MPNN_PY=/home/HarjasG/miniforge3/envs/mpnn/bin/python
PH_PY=/home/HarjasG/tnf-track3/venv-phscan/bin/python
note(){ echo "[$(date -u '+%F %T')] $*" | tee -a logs/71_overnight.log; }
guard(){ local g; g=$(df -BG --output=avail / | tail -1 | tr -dc '0-9')
  if [ "$g" -lt 5 ]; then note "HALT: ${g} GB free"; exit 3; fi; }

note "=== overnight round-2 chain start ==="

# 1. MPNN, three pools. 'targeted' biases His only at binder positions facing a
#    conserved cation; plain/hisb are the controls that let us measure whether
#    targeting beat scattering.
guard
if [ ! -d gen/r2b_mpnn/plain/seqs ]; then
  note "MPNN plain+hisb on 170 round-2 backbones"
  $MPNN_PY scripts/03_mpnn.py --rfdiff gen/r2b_backbones --out gen/r2b_mpnn \
      --num-seq 4 --seed 77 >>logs/71_mpnn.log 2>&1 && note "MPNN ok" || note "MPNN FAILED"
fi
if [ ! -d gen/r2b_mpnn_t/targeted/seqs ]; then
  note "MPNN targeted-bias pool"
  mkdir -p gen/r2b_mpnn_t
  $MPNN_PY /home/HarjasG/tnf-track3/ProteinMPNN/protein_mpnn_run.py \
      --jsonl_path gen/r2b_mpnn/parsed.jsonl \
      --chain_id_jsonl gen/r2b_mpnn/assigned.jsonl \
      --bias_by_res_jsonl gen/bias_r2.jsonl \
      --out_folder gen/r2b_mpnn_t/targeted \
      --num_seq_per_target 4 --sampling_temp 0.1 --seed 77 --batch_size 1 \
      --omit_AAs CX --use_soluble_model >>logs/71_mpnn_t.log 2>&1 \
      && note "MPNN targeted ok ($(find gen/r2b_mpnn_t/targeted/seqs -name '*.fa' 2>/dev/null | wc -l) fa)" \
      || note "MPNN targeted FAILED"
fi

# 2. select
guard
note "select"
$PH_PY -I scripts/35_select_for_boltz.py --mpnn gen/r2b_mpnn --cap 420 \
    --out-fasta gen/r2_selected.fasta --out-tsv gen/r2_selected.tsv \
    >>logs/71_select.log 2>&1 && note "select ok" || note "select FAILED"
if [ -d gen/r2b_mpnn_t/targeted/seqs ]; then
  $PH_PY -I scripts/35_select_for_boltz.py --mpnn gen/r2b_mpnn_t --pools targeted \
      --cap 180 --out-fasta gen/r2_selected_t.fasta --out-tsv gen/r2_selected_t.tsv \
      >>logs/71_select.log 2>&1 && cat gen/r2_selected_t.fasta >> gen/r2_selected.fasta
fi
note "round-2 designs selected: $(grep -c '^>' gen/r2_selected.fasta)"

# 3. Boltz, full monomer
guard
note "Boltz round-2 human screen"
$BOLTZ_PY scripts/04_boltz.py --extra-seqs gen/r2_selected.fasta \
    --tag screen_r2 --target-mode full --species human \
    >>logs/71_boltz_r2.log 2>&1 && note "Boltz r2 ok" || note "Boltz r2 FAILED"

# 4. epitope gate
guard
note "epitope QC round 2"
$PH_PY -I scripts/04e_epitope_qc.py --pred boltz/screen_r2/pred \
    --out gen/epitope_qc_r2.tsv >>logs/71_qc_r2.log 2>&1 \
    && note "QC ok: $(grep -c 'True' gen/epitope_qc_r2.tsv 2>/dev/null) pass" || note "QC FAILED"

# 5. pH scan on round-2 passers
guard
note "pH scan round 2"
$PH_PY -I - <<'PY' >>logs/71_phlist.log 2>&1
import csv,glob,os
rows=[r for r in csv.DictReader(open("gen/epitope_qc_r2.tsv"),delimiter="\t") if r["passes"]=="True"]
out=[]
for r in rows:
    h=glob.glob(f"boltz/screen_r2/pred/**/{r['name']}_model_0.pdb",recursive=True)
    if h: out.append(h[0])
open("gen/r2_passing_complexes.txt","w").write("\n".join(out)+"\n")
print(len(out))
PY
if [ -s gen/r2_passing_complexes.txt ]; then
  $PH_PY -I scripts/05_ph_scan.py --complex $(cat gen/r2_passing_complexes.txt | tr '\n' ' ') \
      --binder-chain X --out results/ph_scan_r2.tsv --combos-out results/ph_combos_r2.tsv \
      --jobs 10 >>logs/71_phscan_r2.log 2>&1 && note "pH scan r2 ok" || note "pH scan r2 FAILED"
  $PH_PY -I scripts/05b_switch_select.py --scan results/ph_scan_r2.tsv \
      --pred boltz/screen_r2/pred --out results/switch_candidates_r2.tsv \
      >>logs/71_switch_r2.log 2>&1 && note "switch select r2 ok" || note "switch r2 FAILED"
fi

note "=== overnight round-2 chain complete ==="
note "free: $(df -h / | awk 'NR==2{print $4}')"
echo OVERNIGHT_DONE
