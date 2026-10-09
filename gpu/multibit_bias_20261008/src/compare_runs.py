#!/usr/bin/env python3
"""Identity check of a run against an earlier kernel's raw files (same schedule and trial ids).
Usage: compare_runs.py NEW_PREFIX OLD_PREFIX. Compares per-trial (trial_id, cut, flips), the final-spin file byte for byte and,
if both runs traced, the n_lin trace file. Prints one JSON line; exits 3 if anything differs."""
import json, sys, os, hashlib

def trials(prefix):
    return [(d['trial_id'], d['cut'], d['flips']) for d in map(json.loads, open(prefix + '.trials.jsonl'))]

def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()

new, ref = sys.argv[1], sys.argv[2]
a, b = trials(new), trials(ref)
same_trials = a == b
sa, sb = digest(new + '.spins.bin'), digest(ref + '.spins.bin')
tr = None
if os.path.exists(new + '.trace.bin') and os.path.exists(ref + '.trace.bin'):
    tr = digest(new + '.trace.bin') == digest(ref + '.trace.bin')
ok = same_trials and sa == sb and tr is not False
print(json.dumps(dict(prefix=os.path.basename(new), v2_prefix=ref, trials=len(a), trials_identical=same_trials,
                      spins_sha256=sa, spins_identical=sa == sb, trace_identical=tr, identical=ok)))
sys.exit(0 if ok else 3)
