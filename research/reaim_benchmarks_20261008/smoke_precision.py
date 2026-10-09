"""Smoke test of precision.py job() and export_hw paths on fabricated schedules (temp dir, 4 trials); no study output."""
import os, sys, tempfile, json
from pathlib import Path
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import precision as PR
PR.OUT = Path(tempfile.mkdtemp(dir=os.environ.get('TMPDIR'))); PR.TRIALS = 4
P = PR.problem('tsp', 'gr17'); a = P.sigma_T; st = P.stats()
cfgs = [dict(family='tec', T0=16 * a, tfin=0.6 * a, q=16 * a, jv=-8 * a, S=512),
        dict(family='tecT', T0=16 * a, tfin=0.6 * a, q=16 * a, kappa=2 ** -8 * st['kappa_factor'], ramp=True, S=512),
        dict(family='onsager', T0=16 * a, tfin=0.6 * a, q=16 * a, lam=2 ** -6 * st['lam_factor'], ramp=True, S=512)]
for rule, c in zip(('TEC', 'Onsager-kT', 'Onsager-online'), cfgs):
    s = dict(problem='tsp', instance='gr17', rule=rule, S=512, cfg=c)
    for K in (None, 2, 8):
        print(PR.job((0, s, K)))
G = PR.problem('gpp', 14); a = G.sigma_T
s = dict(problem='gpp', instance=14, rule='SCA', S=256, cfg=dict(family='plain', T0=4 * a, tfin=0.15 * a, q=16 * a, S=256))
for K in (None, 2, 3):
    print(PR.job((0, s, K)))
M = PR.problem('mcp', 1); a = M.sigma_T
s = dict(problem='mcp', instance=1, rule='Onsager-online', S=250, cfg=dict(family='onsager', T0=1.35 * a, tfin=0.672 * a, q=2.16 * a, lam=0.5 * M.stats()['lam_factor'], ramp=True, S=250))
for K in (None, 2):
    print(PR.job((0, s, K)))
f = sorted((PR.OUT / 'runs').glob('*.json'))[0]; d = json.loads(f.read_text()); print(f.name, d['schedule'], d['checks'])
print('smoke precision ok')
