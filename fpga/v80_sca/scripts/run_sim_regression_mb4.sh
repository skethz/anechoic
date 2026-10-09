#!/usr/bin/env bash
# r4 RTL regression (sca_core_mb_r4.v: revision 3 with 8192-entry schedule tables): a revision-3 subset (identity) plus the
# exported problems of research/reaim_benchmarks_20261008/hw_export at K = 4 (TSP with bias at S = 8192, GPP at S <= 2048).
# Each run detached in its own fresh directory under build/sim_mb4_<tag>/. Usage: bash scripts/run_sim_regression_mb4.sh <tag>
set -uo pipefail
source /scratch/USER/sca_v80_20261003/scripts/environment.sh > /dev/null
tag=$1; base=$TASK_ROOT/build/sim_mb4_$tag; C=$TASK_ROOT/data/r4_cases; D=$TASK_ROOT/data; A=$D/a4_problems
mkdir -p "$base" "$TASK_ROOT/logs/sim_mb4_$tag"
L() { local envs=$1 name=$2; shift 2
  env $envs setsid -f nohup bash "$TASK_ROOT/scripts/sim_mb4.sh" "$base/$name" "$@" > "$TASK_ROOT/logs/sim_mb4_$tag/$name.log" 2>&1 < /dev/null; }
K2000=k2000:$D/K2000.bin
L "BIAS=0 HBX=0" k2_k2000_short_n2000 2 $K2000 $C/cases_k2000_short.txt 0 2000
L "BIAS=1 HBX=17" k4b_k2000_short_bias0 4 $K2000 $C/cases_k2000_short.txt 1 2000
L "BIAS=1 HBX=17" k4b_rand4_n777 4 rand:20261008:0.5:7:777:300 $C/cases_rand4_bias.txt 1
L "BIAS=1 HBX=17" k4b_tsp25 4 ising:$D/ising_20261008/tsp5.ising $C/cases_bias_small.txt 1
L "BIAS=1 HBX=17" k4b_gr17_S8192 4 jint8:$A/gr17_K4.jint8:$A/gr17_K4.bias.bin $C/cases_gr17_K4.txt 1
L "BIAS=1 HBX=17" k4b_bays29_S8192 4 jint8:$A/bays29_K4.jint8:$A/bays29_K4.bias.bin $C/cases_bays29_K4.txt 1
L "BIAS=1 HBX=17" k4b_G1gpp_bm1 4 jint8:$A/G1_K4.jint8:$A/G1_K4.bias.bin $C/cases_G1gpp_K4.txt 1
L "BIAS=1 HBX=17" k4b_G1gpp_bm0 4 jint8:$A/G1_K4.jint8:$A/G1_K4.bias.bin $C/cases_G1gpp_K4.txt 0
echo "launched into $base"
