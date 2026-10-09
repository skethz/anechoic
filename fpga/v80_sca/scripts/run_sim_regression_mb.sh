#!/usr/bin/env bash
# Launch the multi-bit RTL regression on fpga-host (each run detached, its own fresh directory under build/sim_mb_<tag>/).
# Usage: bash scripts/run_sim_regression_mb.sh <tag> [set]   (set: all (default) | k2000 | gset | rand | v64)
set -uo pipefail
source /scratch/USER/sca_v80_20261003/scripts/environment.sh > /dev/null
tag=$1; set=${2:-all}; base=$TASK_ROOT/build/sim_mb_$tag
mkdir -p "$base" "$TASK_ROOT/logs/sim_mb_$tag"
K2=k2000:$TASK_ROOT/data/K2000.bin; G22=gset:$TASK_ROOT/data/gset/G22; G32=gset:$TASK_ROOT/data/gset/G32; R4=rand:20261007:0.5:7
V64="SIM_CORE=$TASK_ROOT/build/mbcheck_20261007/sca_core_v64b_xsimidx.v SIM_TB=$TASK_ROOT/build/mbcheck_20261007/tb_core_v64x.v"
L() {  # env assignments (space-separated, may be empty), run name, then sim_mb.sh arguments after the run dir
  local envs=$1 name=$2; shift 2
  env $envs setsid -f nohup bash "$TASK_ROOT/scripts/sim_mb.sh" "$base/$name" "$@" > "$TASK_ROOT/logs/sim_mb_$tag/$name.log" 2>&1 < /dev/null
}
if [ "$set" = all ] || [ "$set" = k2000 ]; then
  L "" k2000_short_k2 2 $K2 short
  for c in 0 1 2 3 4 5 6; do L "" k2000_proto_c${c}_k2 2 $K2 proto $c; done
  L "DRAIN=1" k2000_short_k2_drain1 2 $K2 short
fi
if [ "$set" = all ] || [ "$set" = gset ]; then
  L "" g22_short_k2 2 $G22 short
  L "" g22_long_k2 2 $G22 long
  L "" g32_short_k2 2 $G32 short
  L "" g32_long_k2 2 $G32 long
fi
if [ "$set" = all ] || [ "$set" = rand ]; then
  L "" rand4_short_k4 4 $R4 short
  L "" rand4_long_k4 4 $R4 long
fi
if [ "$set" = all ] || [ "$set" = v64 ]; then
  # v6.4b (sim-only copy with the xsim index rewrite) on the same K2000 cases with K = 1 vectors: cycle reference
  L "$V64" v64x_k2000_short_k1 1 $K2 short
  for c in 0 3 6; do L "$V64" v64x_k2000_proto_c${c}_k1 1 $K2 proto $c; done
fi
echo "launched $set into $base"
