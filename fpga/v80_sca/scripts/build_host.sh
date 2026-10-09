#!/usr/bin/env bash
# Build the physical runner against the generated register map and timer offsets (after prepare_ip.py).
set -euo pipefail
source /scratch/USER/sca_v80_20261003/scripts/environment.sh
sfx=${SCA_VARIANT:+_$SCA_VARIANT}
DRV=$(dirname "$(find "$TASK_ROOT/build/hls_sca${sfx}/solution/impl/misc/drivers" -name xsca_v80_hw.h | head -1)")
test -f "$TASK_ROOT/build/sca_timer_offsets${sfx}.h"
"$SNOWBALL_ROOT/build/scratch_exec" /scratch/USER g++ -O3 -std=c++17 -Wall -Wextra \
  ${SCA_CFLAGS:-} -include "$TASK_ROOT/build/sca_timer_offsets${sfx}.h" -I"$DRV" -I"$TASK_ROOT/src" -I"$AVED_ROOT/sw/AMI/api/include" \
  src/${SCA_HOST_SRC:-host_sca.cpp} /usr/local/lib/libami.a -pthread -o build/host_sca${sfx}
# Read-only UUID probe; recorded but not fatal (the board may hold a different image during unattended builds).
"$SNOWBALL_ROOT/build/scratch_exec" /scratch/USER --v80-ami build/host_sca${sfx} --probe > results/board_probe${sfx}.json \
  || echo '{"probe":"failed"}' > results/board_probe${sfx}.json
cat results/board_probe${sfx}.json
