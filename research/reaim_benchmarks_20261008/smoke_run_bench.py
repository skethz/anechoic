"""Smoke test of run_bench.job_budget on SYNTHETIC instances (temp output dir, reduced runs); no benchmark instance."""
import os, sys, tempfile, time
from pathlib import Path
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import run_bench as RB
import calib as C
RB.OUT = Path(tempfile.mkdtemp(dir=os.environ.get('TMPDIR')))
RB.B_PILOT, RB.B_FINAL = 2, 4
RB._PROBS[('gpp', 1)] = C.syn_gpp('rand', 7); RB._PROBS[('gpp', 1)].meta.update(ref=7600, target=7676)
RB._PROBS[('tsp', 'gr17')] = C.syn_tsp(9, 8)
t0 = time.time()
for prob, inst, S in (('gpp', 1, 256), ('tsp', 'gr17', 512)):
    for m in RB.METHODS:
        print(RB.job_budget((prob, inst, m, S)), flush=True)
print('smoke ok', time.time() - t0)
