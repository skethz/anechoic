#!/usr/bin/env bash
# Start ABL2 only after the ABL held-out stage has finished (avoids CPU contention).
cd "$(cd "$(dirname "$0")/." && pwd)"
while pgrep -f "abl.py all" >/dev/null; do sleep 20; done
python3 abl2.py all > abl2.log 2>&1
