#!/usr/bin/env bash
# Driver for PROTOCOL.md Amendment 5 (GPU-selected G-set schedules). Usage (on gpu-host, from the project root):
#   bash src/run_gsel.sh STEP   (GSP plan | GSV verification | GSR counted runs | GSW power)
# All processes are pinned to 16 cores (CORES, default 128-143), the CPU limit for this phase.
set -euo pipefail
R=/scratch/USER/anechoic_gpu_multibit_bias_20261008
source $R/environment.sh
CORES=${CORES:-128-143}
TS="taskset -c $CORES"
BIN=build/sca_gpu_bias; CHK=build/check_ref_bias; SEED=20261004
UUID=GPU-xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx
gpulog() { echo "$(date -u +%FT%TZ) $1 $(nvidia-smi --id=$UUID --query-gpu=utilization.gpu,memory.used,power.draw,module.power.draw.average,clocks.sm --format=csv,noheader) | load $(cut -d' ' -f1-3 /proc/loadavg)" >> logs/gpu_status.log; }
run() { local tag=$1; shift; gpulog "start $tag"; local rc=0; $TS $BIN "$@" --seed $SEED --output $tag > $tag.out 2>&1 < /dev/null || rc=$?; gpulog "end $tag"; tail -1 $tag.out; return $rc; }
declare -A VAR=(
  [r2g1]="--kbits 2 --r 2 --preg 2 --cs 8 --groups 1" [r2g2]="--kbits 2 --r 2 --preg 2 --cs 8 --groups 2"
  [r2g4]="--kbits 2 --r 2 --preg 2 --cs 8 --groups 4" [r1g1]="--kbits 2 --r 1 --preg 1 --cs 8 --groups 1"
  [r1g2]="--kbits 2 --r 1 --preg 1 --cs 8 --groups 2"
)
# the variant rule of phase G (results/K selection of the multibit study, register fallback when the selection cannot hold B chains
# in one wave at the schedule length), so that every cohort runs the variant its cost-model cell was measured with
declare -A SELV=([60_const]=r2g1 [60_ons]=r2g1 [120_const]=r1g1 [120_ons]=r1g1 [240_const]=r1g2 [240_ons]=r2g4)
declare -A FALLBACK=([60]=r2g1 [120]=r2g2 [240]=r2g4)
pick() {
  local B=$1 cls=$2; shift 2; local vn=${SELV[${B}_${cls}]}
  local cpw=$($TS $BIN --info --gset data/gset/G1 --target 0 ${VAR[$vn]} "$@" < /dev/null | python3 -c "import json,sys;print(json.load(sys.stdin)['chains_per_wave'])")
  if [ "$cpw" -lt "$B" ]; then echo ${FALLBACK[$B]}; else echo $vn; fi; }
step=${1:?step}
mkdir -p logs results/GS
case $step in
GSP) $TS python3 src/gsel.py data/gset_study_results ;;
GSV) mkdir -p results/GSV
    k=0
    while IFS=$'\t' read -r key inst N target rule cls S B flags base labels <&3; do
      vn=$(pick $B $cls $flags); t=results/GSV/$key
      run $t --gset data/gset/$inst --target $target ${VAR[$vn]} --chains 4 --batches 1 --trace 1 --check-fields 1 --offset $((160000000 + 1000 * k)) $flags
      $TS $CHK --gset data/gset/$inst --prefix $t $flags --seed $SEED --threads 16 >> results/GSV/check.jsonl
      k=$((k + 1))
    done 3< results/GS/plan_cohorts.tsv
    python3 - <<'PY'
import json, glob
r=[json.loads(l) for l in open('results/GSV/check.jsonl')]
dm=sum(json.load(open(f))['mismatches'] for f in glob.glob('results/GSV/*.summary.json'))
fm=sum(json.loads(l).get('field_mismatches',0) for f in glob.glob('results/GSV/*.trials.jsonl') for l in open(f))
ok = len(r)==sum(1 for _ in open('results/GS/plan_cohorts.tsv')) and all(x['mismatch']==0 and x['exact']==x['trials'] for x in r) and dm==0 and fm==0
json.dump(dict(runs=len(r),trials=sum(x['trials'] for x in r),exact=sum(x['exact'] for x in r),device_host_mismatches=dm,field_mismatches=fm),open('results/GSV/summary.json','w'))
print('GS VERIFICATION', 'PASS' if ok else 'FAIL', sum(x['exact'] for x in r), '/', sum(x['trials'] for x in r), 'in', len(r), 'runs; device/host', dm, 'fields', fm)
assert ok
PY
    ;;
GSR) mkdir -p results/GSR
    while IFS=$'\t' read -r key inst N target rule cls S B flags base labels <&3; do
      vn=$(pick $B $cls $flags)
      run results/GSR/$key --gset data/gset/$inst --target $target ${VAR[$vn]} --chains $B --batches 128 --offset $base $flags
    done 3< results/GS/plan_cohorts.tsv
    ;;
GSW) mkdir -p results/GSW
    # one 40 s burn per (B, variant, class) combination used by the counted cohorts: the first such cohort in plan order
    python3 - > results/GSW/burns.tsv <<'PY'
import json
seen = {}
for line in open('results/GS/plan_cohorts.tsv'):
    f = line.rstrip('\n').split('\t'); key = f[0]
    s = json.load(open(f'results/GSR/{key}.summary.json'))
    combo = f"B{s['chains']}_R{s['r']}G{s['groups']}_{f[5]}"
    if combo not in seen: seen[combo] = f
for i, (combo, f) in enumerate(seen.items()):
    print('\t'.join([combo, str(i)] + f[:10]))
PY
    nvidia-smi --id=$UUID --query-gpu=timestamp,utilization.gpu,memory.used,power.draw,power.draw.instant,module.power.draw.average,module.power.draw.instant,clocks.sm --format=csv,noheader -lms 200 > results/GSW/power_samples.csv &
    SMI=$!
    trap "kill $SMI 2>/dev/null" EXIT
    echo "$(date -u +%FT%T.%3NZ) idle_start" >> results/GSW/marks.txt; sleep 30
    while IFS=$'\t' read -r combo i key inst N target rule cls S B flags base <&3; do
      vn=$(pick $B $cls $flags)
      echo "$(date -u +%FT%T.%3NZ) burn_${combo}_start" >> results/GSW/marks.txt
      $TS $BIN --gset data/gset/$inst --target $target ${VAR[$vn]} --chains $B --burn 40 --offset $((170000000 + 1000000 * i)) --seed $SEED $flags --output results/GSW/burn_$combo > results/GSW/burn_$combo.out 2>&1 < /dev/null
      echo "$(date -u +%FT%T.%3NZ) burn_${combo}_end" >> results/GSW/marks.txt; sleep 20
    done 3< results/GSW/burns.tsv
    echo "$(date -u +%FT%T.%3NZ) idle_end" >> results/GSW/marks.txt
    ;;
*) echo "unknown step $step"; exit 1;;
esac
