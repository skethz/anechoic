#!/usr/bin/env bash
# Driver for PROTOCOL.md (bias + runtime-n study). Usage (on gpu-host, from the project root): bash src/run_bias.sh PHASE
#   V verification | T timing | GV G-set validation | G G-set runs | W power
# Every run refuses to overwrite its outputs. V must pass before anything else; GV before G.
set -euo pipefail
R=/scratch/USER/anechoic_gpu_multibit_bias_20261008
source $R/environment.sh
BIN=build/sca_gpu_bias; CHK=build/check_ref_bias; SEED=20261004
MB=/scratch/USER/anechoic_gpu_multibit_20261007; MBBIN=$MB/build/sca_gpu_mb
UUID=GPU-xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx
K2=data/K2000.bin; D=data/inst
gpulog() { echo "$(date -u +%FT%TZ) $1 $(nvidia-smi --id=$UUID --query-gpu=utilization.gpu,memory.used,power.draw,module.power.draw.average,clocks.sm --format=csv,noheader) | load $(cut -d' ' -f1-3 /proc/loadavg)" >> logs/gpu_status.log; }
run() {  # run TAG ARGS...  (one GPU process; status logged before/after)
  local tag=$1; shift
  gpulog "start $tag"; local rc=0; $BIN "$@" --seed $SEED --output $tag > $tag.out 2>&1 || rc=$?; gpulog "end $tag"
  tail -1 $tag.out; return $rc
}
runbin() {  # runbin BINARY TAG ARGS...  (timing runs with either binary)
  local bin=$1 tag=$2; shift 2
  gpulog "start $tag"; local rc=0; $bin "$@" --seed $SEED --output $tag > $tag.out 2>&1 || rc=$?; gpulog "end $tag"
  tail -1 $tag.out; return $rc
}
declare -A VAR=(
  [r2g1]="--kbits 2 --r 2 --preg 2 --cs 8 --groups 1" [r2g2]="--kbits 2 --r 2 --preg 2 --cs 8 --groups 2"
  [r2g3]="--kbits 2 --r 2 --preg 2 --cs 8 --groups 3" [r2g4]="--kbits 2 --r 2 --preg 2 --cs 8 --groups 4"
  [r1g1]="--kbits 2 --r 1 --preg 1 --cs 8 --groups 1" [r1g2]="--kbits 2 --r 1 --preg 1 --cs 8 --groups 2"
  [k4]="--kbits 4 --r 2 --preg 2 --cs 8 --groups 1"
)
VK2=(r2g1 r2g2 r2g3 r2g4 r1g1 r1g2 k4)
declare -A WAVE=([r2g1]=60 [r2g2]=120 [r2g3]=180 [r2g4]=240 [r1g1]=120 [r1g2]=240 [k4]=60)
# phase-K selection of the multibit study (results/K/selected_variants.json there): variant per (B, class)
declare -A SELV=([1_const]=r2g1 [1_ons]=r1g1 [60_const]=r2g1 [60_ons]=r2g1 [120_const]=r1g1 [120_ons]=r1g1 [240_const]=r1g2 [240_ons]=r2g4)
declare -A FALLBACK=([1]=r2g1 [60]=r2g1 [120]=r2g2 [240]=r2g4)
pick() {  # pick B CLASS FLAGS... -> variant (phase-K selection; register variant if the selection cannot hold B chains in one wave at
          # this schedule length -- the step tables live in shared memory)
  local B=$1 cls=$2; shift 2; local vn=${SELV[${B}_${cls}]}
  local cpw=$($BIN --info --gset data/gset/G1 --target 0 ${VAR[$vn]} "$@" < /dev/null | python3 -c "import json,sys;print(json.load(sys.stdin)['chains_per_wave'])")
  if [ "$cpw" -lt "$B" ]; then echo ${FALLBACK[$B]}; else echo $vn; fi; }
phase=${1:?phase}
mkdir -p logs
case $phase in
V)  mkdir -p results/V
    vchk() { local tag=$1; shift; $CHK "$@" --seed $SEED --prefix $tag --threads 64 | tee -a results/V/check.jsonl; }
    # ---- V0: K2000, b = 0: every multibit phase-V1 K2000 run (same schedule, chains, first trial id) with every variant ----
    PL="--q 4 --lambda 0 --t0 30 --steps 40"; ON="--q 8 --lambda 1.05 --ramp --t0 12 --steps 960"
    TT="--q 8 --tecT 1.75 --ramp --t0 15 --steps 280"; TE="--q 8 --tec-jv -4 --t0 30 --steps 100"
    V2RUNS=("plain_c8g1 4 0 PL" "plain_c4g2 8 900000 PL" "plain_l2 2 77 PL" "ons_c8g2 8 900100 ON" "ons_c4g1 4 123456 ON"
            "ons_l2 2 901234 ON" "ons_c8g1_wave 120 940000 ON" "tect_c8g1 4 2000000 TT" "tect_c4g2 8 950000 TT"
            "tect_c8g2 8 31337 TT" "tect_l2 2 999 TT" "tect_c8g2_wave 240 960000 TT" "tec_c8g2 8 999990 TE"
            "tec_c4g1 4 42 TE" "tec_l2 2 905000 TE")
    [ -f $D/zero2000.bias ] || build/gen_ising zero --n 2000 --out $D/zero2000.bias > /dev/null
    for spec in "${V2RUNS[@]}"; do read v2tag ch off sv <<< "$spec"; sched=${!sv}
      for vn in "${VK2[@]}"; do
        tag=results/V/k2000_${v2tag}_$vn
        run $tag --graph $K2 ${VAR[$vn]} --chains $ch --batches 1 --trace 1 --check-fields 1 --offset $off $sched
        vchk $tag --graph $K2 $sched
        python3 src/compare_runs.py $tag $MB/results/V/k2000_${v2tag}_$vn | tee -a results/V/mb_identity.jsonl
        case $v2tag in plain_c8g1|ons_c8g2|tect_c8g1|tec_c8g2)
          tag=results/V/k2000z_${v2tag}_$vn   # explicit all-zero bias file
          run $tag --graph $K2 --bias $D/zero2000.bias --target 33000 ${VAR[$vn]} --chains $ch --batches 1 --trace 1 --check-fields 1 --offset $off $sched
          vchk $tag --graph $K2 --bias $D/zero2000.bias $sched
          python3 src/compare_runs.py $tag $MB/results/V/k2000_${v2tag}_$vn | tee -a results/V/mb_identity.jsonl
          tag=results/V/k2000w_${v2tag}_$vn   # clamped (WIDE) kernel forced
          run $tag --graph $K2 --wide 1 ${VAR[$vn]} --chains $ch --batches 1 --trace 1 --check-fields 1 --offset $off $sched
          vchk $tag --graph $K2 $sched
          python3 src/compare_runs.py $tag $MB/results/V/k2000_${v2tag}_$vn | tee -a results/V/mb_identity.jsonl ;;
        esac
      done; done
    # ---- V1: biased instances (random K = 2 / K = 4 with n < 2000 and n = 1, 2048; WIDE instances; GP; TSP) ----
    declare -A BASE=(  # q t0 t1 | lambda kappa jv
      [rb2_n1000]="4 8 2.5|0.5 0.875 -2"      [rb2_n257]="2.5 5 1.5|0.69 1.21 -1.25"   [rb2_n17]="0.7 3 0.5|0.84 1.47 -0.35"
      [wide2_n1500]="5 10 3|0.5 0.875 -2.5"   [rb4_n800]="21 41 13|18.67 32.7 -10.5"  [rb4_n1999]="34 65 21|18.65 32.6 -17"
      [rb4_n2048]="34 66 21|18.65 32.6 -17"   [rb4_n1]="0.5 2 0.5|1 1 -0.25"           [wide4_n1200]="26 50 16|18.6 32.6 -13"
      [gp22]="40 60 2|1 0.2 -4"               [gp1000]="900 400 5|9 0.5 -90"           [tsp6]="10 60 2|12 0.2 -1"
      [tsp40]="40 300 3|1.9 0.1 -4" )
    sched_of() { local inst=$1 md=$2; local b=${BASE[$inst]}; local qtt=${b%%|*} mp=${b#*|}; read q t0 t1 <<< "$qtt"; read lam kap jv <<< "$mp"
      case $md in pl) echo "--q $q --lambda 0 --t0 $t0 --t1 $t1 --steps 300";; on) echo "--q $q --lambda $lam --ramp --t0 $t0 --t1 $t1 --steps 300";;
                  tt) echo "--q $q --tecT $kap --ramp --t0 $t0 --t1 $t1 --steps 300";; te) echo "--q $q --tec-jv $jv --t0 $t0 --t1 $t1 --steps 150";; esac; }
    off=7000000
    for inst in rb2_n1000 rb2_n257 rb2_n17 wide2_n1500 rb4_n800 rb4_n1999 rb4_n2048 rb4_n1 wide4_n1200 gp22 gp1000 tsp6 tsp40; do
      case $inst in rb2_*|wide2_*) vl=("${VK2[@]}");; *) vl=(k4);; esac
      src="--jint8 $D/$inst.jint8 --bias $D/$inst.bias"
      for md in pl on tt te; do sched=$(sched_of $inst $md)
        for vn in "${vl[@]}"; do tag=results/V/${inst}_${md}_$vn
          run $tag $src --target-energy 0 ${VAR[$vn]} --chains 8 --batches 1 --trace 1 --check-fields 1 --offset $off $sched
          vchk $tag $src $sched; off=$((off + 1000)); done; done
      # full waves
      case $inst in rb2_n1000) wl="r1g2:on r2g1:pl";; rb2_n257) wl="r2g3:pl r2g2:tt";; rb2_n17) wl="r1g1:te";; wide2_n1500) wl="r2g4:tt r1g2:on";;
                    *) wl="k4:on";; esac
      for w in $wl; do vn=${w%%:*}; md=${w#*:}; sched=$(sched_of $inst $md); tag=results/V/${inst}_wave_${md}_$vn
        run $tag $src --target-energy 0 ${VAR[$vn]} --chains ${WAVE[$vn]} --batches 1 --trace 1 --check-fields 1 --offset $off $sched
        vchk $tag $src $sched; off=$((off + 1000)); done
    done
    # clamped (WIDE) kernel forced on biased instances that do not need it
    for spec in "rb2_n1000 on" "rb4_n800 on" "tsp40 on" "gp1000 tt"; do read inst md <<< "$spec"; sched=$(sched_of $inst $md)
      case $inst in rb2_*) vl=("${VK2[@]}");; *) vl=(k4);; esac
      src="--jint8 $D/$inst.jint8 --bias $D/$inst.bias"
      for vn in "${vl[@]}"; do tag=results/V/${inst}_forcedwide_${md}_$vn
        run $tag $src --target-energy 0 --wide 1 ${VAR[$vn]} --chains 8 --batches 1 --trace 1 --check-fields 1 --offset $off $sched
        vchk $tag $src $sched; off=$((off + 1000)); done; done
    python3 - <<'PY'
import json, glob
rows=[json.loads(l) for l in open('results/V/check.jsonl')]
ids=[json.loads(l) for l in open('results/V/mb_identity.jsonl')]
tot=sum(r['trials'] for r in rows); ex=sum(r['exact'] for r in rows); bad=sum(r['mismatch'] for r in rows)
idt=sum(r['trials'] for r in ids); idok=sum(r['trials'] for r in ids if r['identical'])
fm=0; dm=0; wide_runs=0
for f in glob.glob('results/V/*.summary.json'):
    s=json.load(open(f)); dm+=s['mismatches']; wide_runs+=s['wide']
for f in glob.glob('results/V/*.trials.jsonl'):
    for l in open(f): fm+=json.loads(l).get('field_mismatches',0)
out=dict(runs=len(rows),trials=tot,exact=ex,mismatch=bad,mb_identity_runs=len(ids),mb_identity_trials=idt,mb_identical_trials=idok,
         device_host_mismatches=dm,field_mismatches=fm,wide_runs=wide_runs,rows=rows,mb_rows=ids)
json.dump(out,open('results/V/validation_summary.json','w'),indent=1)
ok = bad==0 and ex==tot and idok==idt and dm==0 and fm==0
print('VALIDATION', 'PASS' if ok else 'FAIL', ex, '/', tot, 'reference-exact;', idok, '/', idt, 'identical to the multibit kernel;',
      'device/host mismatches', dm, '; field mismatches', fm, '; WIDE runs', wide_runs)
assert ok
PY
    ;;
T)  mkdir -p results/T
    # T1: bias binary vs multibit binary (K2000, b = 0), interleaved A B A B A B, same trial ids
    declare -A TC=([X5]="--q 8 --tecT 1.75 --ramp --t0 15 --steps 280" [X4]="--q 8 --tecT 2.0 --ramp --t0 15 --steps 760"
                   [O4]="--q 6 --lambda 0.9 --t0 12 --steps 360" [O1]="--q 8 --lambda 1.05 --ramp --t0 12 --steps 960")
    declare -A TCL=([X5]=const [X4]=const [O4]=ons [O1]=ons)
    for rep in 1 2 3; do for c in X5 X4 O4 O1; do for B in 1 60 120 240; do vn=${SELV[${B}_${TCL[$c]}]}
      runbin $BIN results/T/t1_bias_${c}_b${B}_r$rep --graph $K2 ${VAR[$vn]} --chains $B --batches 16 --offset 5000000 ${TC[$c]}
      runbin $MBBIN results/T/t1_mb_${c}_b${B}_r$rep --graph $K2 ${VAR[$vn]} --chains $B --batches 16 --offset 5000000 ${TC[$c]}
    done; done; done
    # T2: with bias vs all-zero bias file (bias binary, same J, same trial ids)
    for rep in 1 2 3; do
      for spec in "rb2_n1000 120 r1g1 --q 4 --lambda 0.5 --ramp --t0 8 --t1 2.5 --steps 300" "rb2_n1000 120 r1g1 --q 4 --lambda 0 --t0 8 --t1 2.5 --steps 300" \
                  "tsp40 60 k4 --q 40 --lambda 1.9 --ramp --t0 300 --t1 3 --steps 300" "tsp40 60 k4 --q 40 --lambda 0 --t0 300 --t1 3 --steps 300"; do
        read inst B vn sched <<< "$spec"; md=$(echo "$sched" | grep -q ramp && echo on || echo pl)
        n=$(python3 -c "import json;print(json.load(open('$D/$inst.json'))['N'])")
        [ -f $D/zero_$inst.bias ] || build/gen_ising zero --n $n --out $D/zero_$inst.bias > /dev/null
        runbin $BIN results/T/t2_bias_${inst}_${md}_r$rep --jint8 $D/$inst.jint8 --bias $D/$inst.bias --target-energy 0 ${VAR[$vn]} --chains $B --batches 16 --offset 5100000 $sched
        runbin $BIN results/T/t2_zero_${inst}_${md}_r$rep --jint8 $D/$inst.jint8 --bias $D/zero_$inst.bias --target-energy 0 ${VAR[$vn]} --chains $B --batches 16 --offset 5100000 $sched
      done; done
    # T3: clamped (WIDE) kernel vs plain kernel (bias binary, K2000)
    for rep in 1 2 3; do for c in X5 O1; do for B in 1 120 240; do vn=${SELV[${B}_${TCL[$c]}]}
      runbin $BIN results/T/t3_plain_${c}_b${B}_r$rep --graph $K2 --wide 0 ${VAR[$vn]} --chains $B --batches 16 --offset 5200000 ${TC[$c]}
      runbin $BIN results/T/t3_wide_${c}_b${B}_r$rep --graph $K2 --wide 1 ${VAR[$vn]} --chains $B --batches 16 --offset 5200000 ${TC[$c]}
    done; done; done
    ;;
GV|G) mkdir -p results/G results/GV
    [ -f results/G/plan.tsv ] || python3 src/gset_plan.py
    if [ $phase = GV ]; then
      k=0
      while IFS=$'\t' read -r tag inst N target kind rule cls S flags base <&3; do
        if [ $kind = E12 ]; then bl="60 120 240"; else bl="1"; fi
        for vn in $(for B in $bl; do pick $B $cls $flags; done | sort -u); do
          t=results/GV/${tag}_$vn
          run $t --gset data/gset/$inst --target $target ${VAR[$vn]} --chains 4 --batches 1 --trace 1 --check-fields 1 --offset $((39000000 + 1000 * k)) $flags
          $CHK --gset data/gset/$inst --prefix $t $flags --seed $SEED --threads 64 >> results/GV/check.jsonl
          k=$((k + 1))
        done
      done 3< results/G/plan.tsv
      python3 - <<'PY'
import json, glob
r=[json.loads(l) for l in open('results/GV/check.jsonl')]
dm=sum(json.load(open(f))['mismatches'] for f in glob.glob('results/GV/*.summary.json'))
ok = r and all(x['mismatch']==0 and x['exact']==x['trials'] for x in r) and dm==0
json.dump(dict(runs=len(r),trials=sum(x['trials'] for x in r),exact=sum(x['exact'] for x in r),device_host_mismatches=dm),open('results/GV/summary.json','w'))
print('G-SET VALIDATION', 'PASS' if ok else 'FAIL', sum(x['exact'] for x in r), '/', sum(x['trials'] for x in r), 'in', len(r), 'runs; device/host mismatches', dm)
assert ok
PY
    else
      while IFS=$'\t' read -r tag inst N target kind rule cls S flags base <&3; do
        if [ $kind = E12 ]; then
          for spec in "60 0" "120 20000" "240 50000"; do read B d <<< "$spec"; vn=$(pick $B $cls $flags)
            run results/G/${tag}_b$B --gset data/gset/$inst --target $target ${VAR[$vn]} --chains $B --batches 128 --offset $((base + d)) $flags; done
        else
          vn=$(pick 1 $cls $flags)
          run results/G/${tag}_b1 --gset data/gset/$inst --target $target ${VAR[$vn]} --chains 1 --batches 256 --offset $((base + 90000)) $flags
        fi
      done 3< results/G/plan.tsv
    fi
    ;;
W)  mkdir -p results/W
    # G22, Onsager-online E12 configuration at B = 120; then tsp40 (biased, K = 4) at B = 60 with its V1 Onsager schedule
    IFS=$'\t' read -r tag inst N target kind rule cls S flags base < <(awk -F'\t' '$1=="G22_Onsageronline_E12"' results/G/plan.tsv)
    vn=$(pick 120 $cls $flags)
    nvidia-smi --id=$UUID --query-gpu=timestamp,utilization.gpu,memory.used,power.draw,power.draw.instant,module.power.draw.average,module.power.draw.instant,clocks.sm --format=csv,noheader -lms 200 > results/W/power_samples.csv &
    SMI=$!
    echo "$(date -u +%FT%T.%3NZ) idle_start" >> results/W/marks.txt; sleep 30
    echo "$(date -u +%FT%T.%3NZ) burn_g_start" >> results/W/marks.txt
    $BIN --gset data/gset/$inst --target $target ${VAR[$vn]} --chains 120 --burn 60 --offset 46000000 --seed $SEED $flags --output results/W/burn_g22 > results/W/burn_g22.out 2>&1
    echo "$(date -u +%FT%T.%3NZ) burn_g_end" >> results/W/marks.txt; sleep 30
    echo "$(date -u +%FT%T.%3NZ) burn_b_start" >> results/W/marks.txt
    $BIN --jint8 $D/tsp40.jint8 --bias $D/tsp40.bias --target-energy 0 ${VAR[k4]} --chains 60 --burn 60 --offset 46500000 --seed $SEED --q 40 --lambda 1.9 --ramp --t0 300 --t1 3 --steps 300 --output results/W/burn_tsp40 > results/W/burn_tsp40.out 2>&1
    echo "$(date -u +%FT%T.%3NZ) burn_b_end" >> results/W/marks.txt; sleep 15
    echo "$(date -u +%FT%T.%3NZ) idle_end" >> results/W/marks.txt
    kill $SMI
    ;;
*) echo "unknown phase $phase"; exit 1;;
esac
