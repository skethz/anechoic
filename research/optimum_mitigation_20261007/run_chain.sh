#!/bin/bash
# Runs the pre-registered stages in order after stage a1, then the exploratory additions.
cd "$(cd "$(dirname "$0")/." && pwd)"
export NUMBA_NUM_THREADS=6 OMP_NUM_THREADS=6 VECLIB_MAXIMUM_THREADS=6 OPENBLAS_NUM_THREADS=6
while pgrep -f "study.py a1" >/dev/null; do sleep 5; done
for st in a2 c hold b hold_b b48; do
  python3 study.py $st > logs/$st.log 2>&1 || { echo "stage $st failed" >> logs/chain.log; exit 1; }
  echo "$(date '+%H:%M:%S') stage $st done" >> logs/chain.log
done
for e in e1 e2; do
  python3 explore.py $e > logs/$e.log 2>&1 || { echo "explore $e failed" >> logs/chain.log; exit 1; }
  echo "$(date '+%H:%M:%S') explore $e done" >> logs/chain.log
done
echo "$(date '+%H:%M:%S') chain complete" >> logs/chain.log
