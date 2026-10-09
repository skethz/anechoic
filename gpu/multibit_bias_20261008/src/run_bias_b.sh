#!/usr/bin/env bash
# Driver for PROTOCOL.md Amendment 1 (binary sca_gpu_bias_tg) and Amendment 2 (graph partitioning / TSP from the ReAIM study's
# hw_export). Usage (on gpu-host, from the project root): bash src/run_bias_b.sh PHASE   (BV | BT | BX | BV2 | B | BW)
set -euo pipefail
R=/scratch/USER/anechoic_gpu_multibit_bias_20261008
source $R/environment.sh
BIN=build/sca_gpu_bias_tg; CHK=build/check_ref_bias; SEED=20261004
MB=/scratch/USER/anechoic_gpu_multibit_20261007
UUID=GPU-xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx
K2=data/K2000.bin; D=data/inst; X=data/hw_export
gpulog() { echo "$(date -u +%FT%TZ) $1 $(nvidia-smi --id=$UUID --query-gpu=utilization.gpu,memory.used,power.draw,module.power.draw.average,clocks.sm --format=csv,noheader) | load $(cut -d' ' -f1-3 /proc/loadavg)" >> logs/gpu_status.log; }
run() { local tag=$1; shift; gpulog "start $tag"; local rc=0; $BIN "$@" --seed $SEED --output $tag > $tag.out 2>&1 < /dev/null || rc=$?; gpulog "end $tag"; tail -1 $tag.out; return $rc; }
declare -A VAR=(
  [r2g1]="--kbits 2 --r 2 --preg 2 --cs 8 --groups 1" [r2g2]="--kbits 2 --r 2 --preg 2 --cs 8 --groups 2"
  [r2g3]="--kbits 2 --r 2 --preg 2 --cs 8 --groups 3" [r2g4]="--kbits 2 --r 2 --preg 2 --cs 8 --groups 4"
  [r1g1]="--kbits 2 --r 1 --preg 1 --cs 8 --groups 1" [r1g2]="--kbits 2 --r 1 --preg 1 --cs 8 --groups 2"
  [k4]="--kbits 4 --r 2 --preg 2 --cs 8 --groups 1"
)
VK2=(r2g1 r2g2 r2g3 r2g4 r1g1 r1g2 k4)
declare -A SELV=([1_const]=r2g1 [1_ons]=r1g1 [60_const]=r2g1 [60_ons]=r2g1 [120_const]=r1g1 [120_ons]=r1g1 [240_const]=r1g2 [240_ons]=r2g4)
PL="--q 4 --lambda 0 --t0 30 --steps 40"; ON="--q 8 --lambda 1.05 --ramp --t0 12 --steps 960"
TT="--q 8 --tecT 1.75 --ramp --t0 15 --steps 280"; TE="--q 8 --tec-jv -4 --t0 30 --steps 100"
phase=${1:?phase}
mkdir -p logs
case $phase in
BV) mkdir -p results/BV
    vchk() { local tag=$1; shift; $CHK "$@" --seed $SEED --prefix $tag --threads 64 | tee -a results/BV/check.jsonl; }
    # (a, b) K2000, b = 0: one multibit phase-V1 run per mode x 7 variants, tables in shared memory and in global memory
    for spec in "plain_c8g1 4 0 PL" "ons_c8g2 8 900100 ON" "tect_c8g1 4 2000000 TT" "tec_c8g2 8 999990 TE"; do read v2tag ch off sv <<< "$spec"; sched=${!sv}
      for vn in "${VK2[@]}"; do for tb in smem global; do
        tag=results/BV/k2000_${v2tag}_${vn}_$tb
        run $tag --graph $K2 --tables $tb ${VAR[$vn]} --chains $ch --batches 1 --trace 1 --check-fields 1 --offset $off $sched
        vchk $tag --graph $K2 $sched
        python3 src/compare_runs.py $tag $MB/results/V/k2000_${v2tag}_$vn | tee -a results/BV/mb_identity.jsonl
      done; done; done
    # (c) biased instances with global tables (schedules of phase V1)
    off=8000000
    for spec in "rb2_n1000|--q 4 --lambda 0.5 --ramp --t0 8 --t1 2.5 --steps 300|${VK2[*]}" \
                "rb4_n800|--q 21 --lambda 0 --t0 41 --t1 13 --steps 300|k4" "rb4_n800|--q 21 --lambda 18.67 --ramp --t0 41 --t1 13 --steps 300|k4" \
                "rb4_n800|--q 21 --tecT 32.7 --ramp --t0 41 --t1 13 --steps 300|k4" "rb4_n800|--q 21 --tec-jv -10.5 --t0 41 --t1 13 --steps 150|k4" \
                "tsp40|--q 40 --lambda 0 --t0 300 --t1 3 --steps 300|k4" "tsp40|--q 40 --lambda 1.9 --ramp --t0 300 --t1 3 --steps 300|k4" \
                "tsp40|--q 40 --tecT 0.1 --ramp --t0 300 --t1 3 --steps 300|k4" "tsp40|--q 40 --tec-jv -4 --t0 300 --t1 3 --steps 150|k4" \
                "gp1000|--q 900 --lambda 9 --ramp --t0 400 --t1 5 --steps 300|k4" "wide4_n1200|--q 26 --lambda 18.6 --ramp --t0 50 --t1 16 --steps 300|k4"; do
      IFS='|' read inst sched vl <<< "$spec"; src="--jint8 $D/$inst.jint8 --bias $D/$inst.bias"
      for vn in $vl; do tag=results/BV/${inst}_tg_$(echo "$sched" | grep -oE 'lambda 0 |lambda [0-9.]+ --ramp|tecT|tec-jv' | head -1 | tr -cd 'a-zA-Z0-9')_$vn
        run $tag $src --target-energy 0 --tables global ${VAR[$vn]} --chains 8 --batches 1 --trace 1 --check-fields 1 --offset $off $sched
        vchk $tag $src $sched; off=$((off + 1000)); done; done
    # (d) long schedules beyond the shared-memory limit (tables in global memory automatically)
    for spec in "tsp40|--q 40 --lambda 1.9 --ramp --t0 300 --t1 3 --steps 8192|k4|8" "tsp40|--q 40 --tecT 0.1 --ramp --t0 300 --t1 3 --steps 8192|k4|60" \
                "rb2_n1000|--q 4 --lambda 0.5 --ramp --t0 8 --t1 2.5 --steps 6000|r2g1|8" "rb2_n1000|--q 4 --lambda 0 --t0 8 --t1 2.5 --steps 5000|r1g2|240"; do
      IFS='|' read inst sched vn ch <<< "$spec"; src="--jint8 $D/$inst.jint8 --bias $D/$inst.bias"
      tag=results/BV/${inst}_long_${vn}_$ch
      run $tag $src --target-energy 0 ${VAR[$vn]} --chains $ch --batches 1 --trace 1 --check-fields 1 --offset $off $sched
      vchk $tag $src $sched; off=$((off + 1000)); done
    tag=results/BV/k2000_long_r1g1; S8="--q 8 --tecT 1.75 --ramp --t0 15 --steps 8192"
    run $tag --graph $K2 ${VAR[r1g1]} --chains 8 --batches 1 --trace 1 --check-fields 1 --offset $off $S8
    vchk $tag --graph $K2 $S8
    python3 - <<'PY'
import json, glob
rows=[json.loads(l) for l in open('results/BV/check.jsonl')]
ids=[json.loads(l) for l in open('results/BV/mb_identity.jsonl')]
dm=0; fm=0; tg=0
for f in glob.glob('results/BV/*.summary.json'):
    s=json.load(open(f)); dm+=s['mismatches']; tg+=s['tables_global']
for f in glob.glob('results/BV/*.trials.jsonl'):
    for l in open(f): fm+=json.loads(l).get('field_mismatches',0)
tot=sum(r['trials'] for r in rows); ex=sum(r['exact'] for r in rows); idt=sum(r['trials'] for r in ids); idok=sum(r['trials'] for r in ids if r['identical'])
ok = ex==tot and idok==idt and dm==0 and fm==0
json.dump(dict(runs=len(rows),trials=tot,exact=ex,mb_identity_trials=idt,mb_identical_trials=idok,device_host_mismatches=dm,field_mismatches=fm,tables_global_runs=tg),open('results/BV/summary.json','w'))
print('BINARY VERIFICATION', 'PASS' if ok else 'FAIL', ex, '/', tot, 'exact;', idok, '/', idt, 'identical to multibit;', 'global-table runs', tg)
assert ok
PY
    ;;
BT) mkdir -p results/BT
    declare -A TC=([X5]="$TT" [O1]="$ON")
    for rep in 1 2 3; do
      for c in X5 O1; do cls=$([ $c = O1 ] && echo ons || echo const)
        for B in 1 60 120 240; do vn=${SELV[${B}_$cls]}
          for tb in smem global; do run results/BT/k2000_${c}_b${B}_${tb}_r$rep --graph $K2 --tables $tb ${VAR[$vn]} --chains $B --batches 16 --offset 5300000 ${TC[$c]}; done
        done; done
      for sched in "--q 40 --lambda 1.9 --ramp --t0 300 --t1 3 --steps 300" "--q 40 --lambda 0 --t0 300 --t1 3 --steps 300"; do md=$(echo "$sched" | grep -q ramp && echo on || echo pl)
        for B in 1 60; do for tb in smem global; do
          run results/BT/tsp40_${md}_b${B}_${tb}_r$rep --jint8 $D/tsp40.jint8 --bias $D/tsp40.bias --target-energy 0 --tables $tb ${VAR[k4]} --chains $B --batches 16 --offset 5400000 $sched
        done; done; done
    done
    ;;
*) echo "unknown phase $phase (BX, BV2, B, BW are defined by Amendment 2)"; exit 1;;
esac
