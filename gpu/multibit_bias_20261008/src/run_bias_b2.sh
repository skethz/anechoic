#!/usr/bin/env bash
# Driver for PROTOCOL.md Amendment 3: graph partitioning (bisection) and TSP from the ReAIM study's hw_export, K = 4.
# Usage (on gpu-host, from the project root): bash src/run_bias_b2.sh STEP   (BX | BV2 | BR | B | BW)
#   BX  import: hashes, bias conversion, objective spec + mapping check, table identity, range/residency checks, plan
#   BV2 validation traces against run_trial_bias (must pass before BR / B)
#   BR  replication of the export's precision-study trials (seed of the entry, trial ids 0..1023)
#   B   counted runs: B = 60 (128 batches) and B = 1 (256 launches), fresh trial ids
#   BW  power
set -euo pipefail
R=/scratch/USER/anechoic_gpu_multibit_bias_20261008
source $R/environment.sh
BIN=build/sca_gpu_bias_tg2; CHK=build/check_ref_bias
UUID=GPU-xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx
X=data/hw_export; K4="--kbits 4 --r 2 --preg 2 --cs 8 --groups 1"
gpulog() { echo "$(date -u +%FT%TZ) $1 $(nvidia-smi --id=$UUID --query-gpu=utilization.gpu,memory.used,power.draw,module.power.draw.average,clocks.sm --format=csv,noheader) | load $(cut -d' ' -f1-3 /proc/loadavg)" >> logs/gpu_status.log; }
run() { local tag=$1 seed=$2; shift 2; gpulog "start $tag"; local rc=0; $BIN "$@" --seed $seed --output $tag > $tag.out 2>&1 < /dev/null || rc=$?; gpulog "end $tag"; tail -1 $tag.out; return $rc; }
step=${1:?step}
mkdir -p logs results/B
case $step in
BX) python3 src/hw_import.py   # writes data/hw_export/plan_b.tsv, data/hw_export/import_check.json; fails on any check
    ;;
BV2) mkdir -p results/BV2
    k=0
    while IFS=$'\t' read -r tag inst kind N rule S seed jf bf spec flags <&3; do
      t=results/BV2/$tag
      run $t $seed --jint8 $jf --bias $bf --target-energy 0 $K4 --chains 4 --batches 1 --trace 1 --check-fields 1 --offset $((30000000 + 1000 * k)) $flags
      $CHK --jint8 $jf --bias $bf --prefix $t $flags --seed $seed --threads 64 >> results/BV2/check.jsonl
      k=$((k + 1))
    done 3< $X/plan_b.tsv
    python3 - <<'PY'
import json, glob
r=[json.loads(l) for l in open('results/BV2/check.jsonl')]
dm=sum(json.load(open(f))['mismatches'] for f in glob.glob('results/BV2/*.summary.json'))
ok = r and all(x['mismatch']==0 and x['exact']==x['trials'] for x in r) and dm==0
json.dump(dict(runs=len(r),trials=sum(x['trials'] for x in r),exact=sum(x['exact'] for x in r),device_host_mismatches=dm),open('results/BV2/summary.json','w'))
print('GP/TSP VALIDATION', 'PASS' if ok else 'FAIL', sum(x['exact'] for x in r), '/', sum(x['trials'] for x in r), 'in', len(r), 'runs')
assert ok
PY
    ;;
BR) mkdir -p results/BR
    while IFS=$'\t' read -r tag inst kind N rule S seed jf bf spec flags <&3; do
      t=results/BR/$tag
      run $t $seed --jint8 $jf --bias $bf --target-energy 0 $K4 --chains 32 --batches 32 --offset 0 $flags
      python3 src/eval_problem.py $spec $t > /dev/null
    done 3< $X/plan_b.tsv
    ;;
B)  while IFS=$'\t' read -r tag inst kind N rule S seed jf bf spec flags <&3; do
      run results/B/${tag}_b60 $seed --jint8 $jf --bias $bf --target-energy 0 $K4 --chains 60 --batches 128 --offset 10000000 $flags
      python3 src/eval_problem.py $spec results/B/${tag}_b60 > /dev/null
      run results/B/${tag}_b1 $seed --jint8 $jf --bias $bf --target-energy 0 $K4 --chains 1 --batches 256 --offset 20000000 $flags
      python3 src/eval_problem.py $spec results/B/${tag}_b1 > /dev/null
    done 3< $X/plan_b.tsv
    ;;
BW) mkdir -p results/BW
    # the TSP and the GPP cell with the largest S x N (first in plan order on ties), B = 60, 60 s each
    read -r ttag tseed tjf tbf tflags < <(python3 src/hw_import.py --power tsp)
    read -r gtag gseed gjf gbf gflags < <(python3 src/hw_import.py --power gpp)
    nvidia-smi --id=$UUID --query-gpu=timestamp,utilization.gpu,memory.used,power.draw,power.draw.instant,module.power.draw.average,module.power.draw.instant,clocks.sm --format=csv,noheader -lms 200 > results/BW/power_samples.csv &
    SMI=$!
    trap "kill $SMI 2>/dev/null" EXIT
    echo "$(date -u +%FT%T.%3NZ) idle_start" >> results/BW/marks.txt; sleep 30
    echo "$(date -u +%FT%T.%3NZ) burn_g_start" >> results/BW/marks.txt
    $BIN --jint8 $tjf --bias $tbf --target-energy 0 $K4 --chains 60 --burn 60 --offset 40000000 --seed $tseed $tflags --output results/BW/burn_$ttag > results/BW/burn_tsp.out 2>&1
    echo "$(date -u +%FT%T.%3NZ) burn_g_end" >> results/BW/marks.txt; sleep 30
    echo "$(date -u +%FT%T.%3NZ) burn_b_start" >> results/BW/marks.txt
    $BIN --jint8 $gjf --bias $gbf --target-energy 0 $K4 --chains 60 --burn 60 --offset 40000000 --seed $gseed $gflags --output results/BW/burn_$gtag > results/BW/burn_gpp.out 2>&1
    echo "$(date -u +%FT%T.%3NZ) burn_b_end" >> results/BW/marks.txt; sleep 15
    echo "$(date -u +%FT%T.%3NZ) idle_end" >> results/BW/marks.txt
    ;;
*) echo "unknown step $step"; exit 1;;
esac
