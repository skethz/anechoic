#!/usr/bin/env bash
# r3 RTL regression (sca_core_mb_r3.v): K = 2 without bias (runtime n; the board image for Max-Cut) and K = 4 with bias
# (HBX = 17; graph partitioning, TSP, random biases). Each run detached in its own fresh directory under build/sim_mb3_<tag>/.
# Usage: bash scripts/run_sim_regression_mb3.sh <tag> [k2|k4|all]
set -uo pipefail
source /scratch/USER/sca_v80_20261003/scripts/environment.sh > /dev/null
tag=$1; set=${2:-all}; base=$TASK_ROOT/build/sim_mb3_$tag; C=$TASK_ROOT/data/r3_cases; D=$TASK_ROOT/data
mkdir -p "$base" "$TASK_ROOT/logs/sim_mb3_$tag"
L() {  # env, name, sim_mb3.sh arguments after the run dir
  local envs=$1 name=$2; shift 2
  env $envs setsid -f nohup bash "$TASK_ROOT/scripts/sim_mb3.sh" "$base/$name" "$@" > "$TASK_ROOT/logs/sim_mb3_$tag/$name.log" 2>&1 < /dev/null
}
K2000=k2000:$D/K2000.bin
if [ "$set" = all ] || [ "$set" = k2 ]; then
  E="BIAS=0 HBX=0"
  L "$E" k2_k2000_short_n2000 2 $K2000 $C/cases_k2000_short.txt 0 2000
  L "$E" k2_k2000_short_n0 2 $K2000 $C/cases_k2000_short.txt 0 0
  L "$E" k2_k2000_proto 2 $K2000 $C/cases_k2000_proto.txt 0 2000
  L "$E" k2_g22_neg 2 gset:$D/gset/G22 $C/cases_neg.txt 0 2000
  for g in G1 G11 G14; do L "$E" k2_${g}_short 2 gset:$D/gset/$g $C/cases_ternary_short.txt 0; done
  for g in G1 G11 G14 G43 G51; do L "$E" k2_${g}_study 2 gset:$D/gset/$g $C/cases_${g}_study.txt 0; done
fi
if [ "$set" = all ] || [ "$set" = k4 ]; then
  E="BIAS=1 HBX=17"
  L "$E" k4b_k2000_short_nobias 4 $K2000 $C/cases_k2000_short.txt 0 2000
  L "$E" k4b_k2000_short_bias0 4 $K2000 $C/cases_k2000_short.txt 1 2000
  L "$E" k4b_rand4_n777 4 rand:20261008:0.5:7:777:300 $C/cases_rand4_bias.txt 1
  L "$E" k4b_rand4_n2048 4 rand:20261009:0.5:7:2048:300 $C/cases_rand4_bias.txt 1
  L "$E" k4b_gp20 4 ising:$D/ising_20261008/gp20_A2.ising $C/cases_bias_small.txt 1
  L "$E" k4b_tsp25 4 ising:$D/ising_20261008/tsp5.ising $C/cases_bias_small.txt 1
  L "$E" k4b_tsp16 4 ising:$D/ising_20261008/tsp4.ising $C/cases_bias_small.txt 1
  L "$E" k4b_g1_short 4 gset:$D/gset/G1 $C/cases_ternary_short.txt 0
fi
echo "launched $set into $base"
