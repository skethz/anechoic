"""Operational helper for Addendum 5: runs run_a5's frozen job list in reverse order (cheapest first) beside the forward
runner, to use idle cores. Each job's seed and output path depend only on the job, and finished jobs are skipped, so
results are identical whichever runner computes a job. Usage: python3 run_a5_rev.py [W]"""
import sys
from pathlib import Path

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_a2 as RA  # noqa: E402
import run_a5 as R5  # noqa: E402

if __name__ == '__main__':
    W = int(sys.argv[1]) if len(sys.argv) > 1 else 32
    RA.run_pool(list(reversed(R5.jobs())), R5.job, W, 'A5 reverse-order helper')
