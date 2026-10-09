#!/usr/bin/env bash
# Driver for PROTOCOL.md. Usage (on gpu-host, from the project root): bash src/run_protocol.sh PHASE
#   V validation | K kernel-variant selection | M V80 configurations | P pilot | H held-out | W power
# Every run refuses to overwrite its outputs. Phases must be run in the order V K M P H W; V must pass before anything else.
set -euo pipefail
R=/scratch/USER/snowball_gpu_v2_20261007
source $R/environment.sh
BIN=build/sca_gpu; CHK=build/check_ref; SEED=20261004
UUID=GPU-xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx
gpulog() { echo "$(date -u +%FT%TZ) $1 $(nvidia-smi --id=$UUID --query-gpu=utilization.gpu,memory.used,power.draw,module.power.draw.average,clocks.sm --format=csv,noheader)" >> logs/gpu_status.log; }
declare -A CFG=(
  [X5]="--q 8 --tecT 1.75 --ramp --t0 15 --steps 280"
  [X1]="--q 8 --tecT 1.75 --ramp --t0 15 --steps 320"
  [X4]="--q 8 --tecT 2.0 --ramp --t0 15 --steps 760"
  [X2]="--q 6 --lambda 0.9 --ramp --t0 15 --steps 320"
  [O1]="--q 8 --lambda 1.05 --ramp --t0 12 --steps 960"
  [O4]="--q 6 --lambda 0.9 --t0 12 --steps 360"
  [P4]="--q 8 --lambda 0 --t0 30 --steps 960"
)
ORDER=(X5 X1 X4 X2 O1 O4 P4)
variant() {  # kernel variant chosen in phase K for a batch size: prints "--cs C --groups G"
  python3 -c "import json;d=json.load(open('results/K/selected_variants.json'));print(d['$1'])"
}
run() {  # run TAG ARGS...  (one GPU process; status logged before/after)
  local tag=$1; shift
  gpulog "start $tag"; $BIN "$@" --seed $SEED --output $tag > $tag.out 2>&1; local rc=$?; gpulog "end $tag"
  tail -1 $tag.out; return $rc
}
phase=${1:?phase}
case $phase in
V)  mkdir -p results/V
    v() { local tag=results/V/$1 cs=$2 g=$3 ch=$4 off=$5; shift 5
          run $tag --cs $cs --groups $g --chains $ch --batches 1 --trace 1 --offset $off "$@"
          $CHK "$@" --seed $SEED --prefix $tag --threads 64 | tee -a results/V/check.jsonl; }
    PL="--q 4 --lambda 0 --t0 30 --steps 40"; ON="--q 8 --lambda 1.05 --ramp --t0 12 --steps 960"
    TT="--q 8 --tecT 1.75 --ramp --t0 15 --steps 280"; TE="--q 8 --tec-jv -4 --t0 30 --steps 100"
    v plain_c8g1 8 1 4 0 $PL;          v plain_c4g2 4 2 8 900000 $PL;     v plain_l2 1 1 2 77 $PL
    v ons_c8g2 8 2 8 900100 $ON;       v ons_c4g1 4 1 4 123456 $ON;      v ons_l2 1 1 2 901234 $ON
    v ons_c8g1_wave 8 1 120 940000 $ON
    v tect_c8g1 8 1 4 2000000 $TT;     v tect_c4g2 4 2 8 950000 $TT;     v tect_c8g2 8 2 8 31337 $TT; v tect_l2 1 1 2 999 $TT
    v tect_c8g2_wave 8 2 240 960000 $TT
    v tec_c8g2 8 2 8 999990 $TE;       v tec_c4g1 4 1 4 42 $TE;          v tec_l2 1 1 2 905000 $TE
    python3 - <<'PY'
import json
rows=[json.loads(l) for l in open('results/V/check.jsonl')]
tot=sum(r['trials'] for r in rows); ex=sum(r['exact'] for r in rows); bad=sum(r['mismatch'] for r in rows)
json.dump(dict(runs=len(rows),trials=tot,exact=ex,mismatch=bad,rows=rows),open('results/V/validation_summary.json','w'),indent=1)
print('VALIDATION', 'PASS' if bad==0 and ex==tot else 'FAIL', ex, '/', tot)
assert bad==0 and ex==tot
PY
    ;;
K)  mkdir -p results/K; X5=${CFG[X5]}
    run results/K/b1_c8g1 --cs 8 --groups 1 --chains 1 --batches 16 --offset 5000000 $X5
    run results/K/b1_c4g1 --cs 4 --groups 1 --chains 1 --batches 16 --offset 5000100 $X5
    run results/K/b1_l2 --cs 1 --groups 1 --chains 1 --batches 4 --offset 5000200 $X5
    run results/K/b120_c8g1 --cs 8 --groups 1 --chains 120 --batches 16 --offset 5010000 $X5
    run results/K/b120_c4g1 --cs 4 --groups 1 --chains 120 --batches 16 --offset 5020000 $X5
    run results/K/b240_c8g2 --cs 8 --groups 2 --chains 240 --batches 16 --offset 5030000 $X5
    run results/K/b240_c4g2 --cs 4 --groups 2 --chains 240 --batches 16 --offset 5040000 $X5
    run results/K/b264_l2 --cs 1 --groups 1 --chains 264 --batches 4 --offset 5050000 $X5
    python3 - <<'PY'
import json
def t(tag): return json.load(open(f'results/K/{tag}.summary.json'))['device_ms_mean']
sel={}
for B,c in [('1',['b1_c8g1','b1_c4g1','b1_l2']),('120',['b120_c8g1','b120_c4g1']),('240',['b240_c8g2','b240_c4g2'])]:
    best=min(c,key=t); cs={'c8':8,'c4':4,'l2':1}[best.split('_')[1][:2]]; g=2 if best.endswith('g2') else 1
    sel[B]=f'--cs {cs} --groups {g}'; print(B,{x:t(x) for x in c},'->',best)
sel['480']=sel['240']
json.dump(sel,open('results/K/selected_variants.json','w'),indent=1)
PY
    ;;
M)  mkdir -p results/M
    for i in "${!ORDER[@]}"; do c=${ORDER[$i]}; base=$((10000000 + i * 1000000))
      run results/M/${c}_b1 $(variant 1) --chains 1 --batches 256 --offset $base ${CFG[$c]}
      run results/M/${c}_b120 $(variant 120) --chains 120 --batches 128 --offset $((base + 100000)) ${CFG[$c]}
      run results/M/${c}_b240 $(variant 240) --chains 240 --batches 128 --offset $((base + 200000)) ${CFG[$c]}
      run results/M/${c}_b480 $(variant 480) --chains 480 --batches 64 --offset $((base + 300000)) ${CFG[$c]}
    done
    ;;
P)  mkdir -p results/P; idx=0
    for mode in "--tecT 1.5" "--tecT 1.75" "--tecT 2.0" "--lambda 0.9" "--lambda 1.05"; do for q in 6 8; do for t0 in 12 15; do
      for S in 80 120 160 200 240 280; do
        base=$((3000000 + idx * 8000)); tag=$(printf "p%03d" $idx)
        args="$mode --ramp --q $q --t0 $t0 --steps $S"
        run results/P/${tag}_b240 $(variant 240) --chains 240 --batches 32 --offset $base $args
        run results/P/${tag}_b120 $(variant 120) --chains 120 --batches 32 --offset $base $args
        echo "$tag $args" >> results/P/index.txt; idx=$((idx + 1))
      done; done; done; done
    python3 - <<'PY'
import json,glob,sys,math
sys.path.insert(0,'src'); import analyze as A
cells=[A.estimate(f[:-13]) for f in sorted(glob.glob('results/P/p*_b*.summary.json'))]
for w in ['primary','secondary']:
    for c in cells: c['score_'+w]=A.selection_score(c,w)
    best=min(cells,key=lambda c:(c['score_'+w],c['tts_'+w+'_ms']))
    json.dump(best,open(f'results/P/selected_{w}.json','w'),indent=1); print(w,best['prefix'],best['score_'+w],best['tts_'+w+'_ms'])
json.dump(cells,open('results/P/pilot_cells.json','w'),indent=0)
PY
    ;;
H)  mkdir -p results/H
    args_of() { python3 -c "
import json; c=json.load(open('results/P/selected_$1.json'))
m=f\"--tecT {c['kappa']:g}\" if c['kappa']>0 else f\"--lambda {c['lam']:g}\"
print(f\"{m} --ramp --q {c['q']:g} --t0 {c['t0']:g} --steps {c['steps']}\")"; }
    B_of() { python3 -c "import json;print(json.load(open('results/P/selected_$1.json'))['B'])"; }
    a1=$(args_of primary); b1=$(B_of primary); a2=$(args_of secondary); b2=$(B_of secondary)
    echo "primary: $a1 B=$b1 ; secondary: $a2 B=$b2" | tee results/H/selected.txt
    run results/H/primary_b${b1} $(variant $b1) --chains $b1 --batches 128 --offset 4000000 $a1
    run results/H/primary_b1 $(variant 1) --chains 1 --batches 256 --offset 4100000 $a1
    ob=$([ "$b1" = 120 ] && echo 240 || echo 120)
    run results/H/primary_b${ob} $(variant $ob) --chains $ob --batches 128 --offset 4200000 $a1
    if [ "$a2 $b2" != "$a1 $b1" ]; then
      run results/H/secondary_b${b2} $(variant $b2) --chains $b2 --batches 128 --offset 4500000 $a2
    fi
    ;;
W)  mkdir -p results/W
    a1=$(python3 -c "
import json; c=json.load(open('results/P/selected_primary.json'))
m=f\"--tecT {c['kappa']:g}\" if c['kappa']>0 else f\"--lambda {c['lam']:g}\"
print(f\"{m} --ramp --q {c['q']:g} --t0 {c['t0']:g} --steps {c['steps']}\")")
    b1=$(python3 -c "import json;print(json.load(open('results/P/selected_primary.json'))['B'])")
    nvidia-smi --id=$UUID --query-gpu=timestamp,utilization.gpu,memory.used,power.draw,power.draw.instant,module.power.draw.average,module.power.draw.instant,clocks.sm --format=csv,noheader -lms 200 > results/W/power_samples.csv &
    SMI=$!
    echo "$(date -u +%FT%T.%3NZ) idle_start" >> results/W/marks.txt; sleep 30
    echo "$(date -u +%FT%T.%3NZ) burn_start" >> results/W/marks.txt
    $BIN $(variant $b1) --chains $b1 --burn 75 --offset 6000000 --seed $SEED --output results/W/burn $a1 > results/W/burn.out 2>&1
    echo "$(date -u +%FT%T.%3NZ) burn_end" >> results/W/marks.txt; sleep 15
    echo "$(date -u +%FT%T.%3NZ) idle_end" >> results/W/marks.txt
    kill $SMI
    ;;
*) echo "unknown phase $phase"; exit 1;;
esac
