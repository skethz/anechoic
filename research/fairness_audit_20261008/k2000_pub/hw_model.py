"""Table 2 support (PROTOCOL_PUB.md, part H): model predictions of the 12-engine V80 TTS99 for schedules at the
published settings of STATICA and TEC, which the board never measured, and for the measured schedules (validation).
Model: float engine model abl.run (bit-faithful arithmetic of the engine), 2,052 trials per schedule (171 rounds of 12),
fresh seeds SeedSequence([20261008, 401, index]); cycles per trial = 0.1424 flips + 18.09 S + 886 at 250 MHz
(paper Eq. cycles, v6.4); t_round = E[max of 12 trial times] from the empirical per-trial distribution;
primary = t_round if P_round = 1-(1-p)^12 >= 0.995 else t_round ln0.01/ln(1-P_round); secondary = t_round ln0.01/(12 ln(1-p)).
No board action. Usage: python3 hw_model.py run [workers] | summarize"""
import json
import math
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import pub_methods as P  # noqa: E402
from common import M, ROOT  # noqa: E402

HERE = Path(__file__).resolve().parent
OUT = HERE / 'results_hw_model'
E_ENG, MHZ, NTRIAL = 12, 250.0, 2052
SEED = 20261008
CONFIGS = [   # key, label, cfg (abl.run families); measured keys refer to results/hwv6_board_v64_e12_250mhz
    ('SP_long', 'STATICA published long point (q 4, T 40->5, S 1560)', dict(family='plain', q=4.0, T0=40.0, S=1560)),
    ('SP_short', 'STATICA published short point (q 4, T 30->5, S 560)', dict(family='plain', q=4.0, T0=30.0, S=560)),
    ('TP_long', 'TEC J_v +30 on STATICA long', dict(family='tec', jv=30.0, q=4.0, T0=40.0, S=1560)),
    ('TP_short', 'TEC J_v +30 on STATICA short', dict(family='tec', jv=30.0, q=4.0, T0=30.0, S=560)),
    ('P1', 'measured P1', dict(family='plain', q=4.0, T0=30.0, S=1560)),
    ('P2', 'measured P2', dict(family='plain', q=6.0, T0=30.0, S=1560)),
    ('P3', 'measured P3 (paper: original plain SCA)', dict(family='plain', q=8.0, T0=30.0, S=1560)),
    ('P4', 'measured P4', dict(family='plain', q=8.0, T0=30.0, S=960)),
    ('T1', 'measured T1 = B4', dict(family='tec', jv=-4.0, q=8.0, T0=30.0, S=1560)),
    ('T2', 'measured T2 = B3', dict(family='tec', jv=-4.0, q=8.0, T0=30.0, S=960)),
    ('T3', 'measured T3', dict(family='tec', jv=4.0, q=8.0, T0=30.0, S=1560)),
    ('B1', 'measured B1', dict(family='plain', q=8.0, T0=40.0, S=960)),
    ('B2', 'measured B2', dict(family='plain', q=6.0, T0=40.0, S=1560)),
    ('X5', 'measured X5 = paper O1', dict(family='tecT', kappa=1.75, q=8.0, T0=15.0, ramp=True, S=280)),
    ('X4', 'measured X4 = paper O5', dict(family='tecT', kappa=2.0, q=8.0, T0=15.0, ramp=True, S=760)),
]


def cycles(flips, S):
    return 0.1424 * flips + 18.09 * S + 886


def e_max(x, k):
    x = np.sort(np.asarray(x, np.float64)); n = len(x); i = np.arange(1, n + 1)
    return float(np.sum(x * ((i / n) ** k - ((i - 1) / n) ** k)))


def tts(t_round, p):
    P_r = 1 - (1 - p) ** E_ENG
    prim = t_round if P_r >= 0.995 else (math.inf if p <= 0 else t_round * math.log(0.01) / math.log1p(-P_r))
    sec = math.inf if p <= 0 else (t_round if p >= 1 else t_round * math.log(0.01) / (E_ENG * math.log1p(-p)))
    return prim, sec, P_r


_J = None


def _init():
    global _J
    _J = M.m.load()


def _run(arg):
    idx, key, label, cfg = arg
    J, sumw = _J; t0 = time.time()
    rng = np.random.default_rng(np.random.SeedSequence([SEED, 401, idx]))
    s0 = rng.choice(np.array([-1.0, 1.0], dtype=np.float32), size=(NTRIAL, len(J)))
    s, flips = P.abl.run(J, s0, dict(cfg), rng)
    cuts = P.m.cut(J, sumw, s)
    return dict(key=key, label=label, cfg=cfg, cuts=cuts.tolist(), flips=flips.tolist(), sec=time.time() - t0)


def run(workers):
    OUT.mkdir(exist_ok=True)
    todo = [(i, k, l, c) for i, (k, l, c) in enumerate(CONFIGS) if not (OUT / f'{k}.json').exists()]
    with ProcessPoolExecutor(max_workers=workers, initializer=_init) as ex:
        for r in ex.map(_run, todo):
            (OUT / f"{r['key']}.json").write_text(json.dumps(r) + '\n'); print(r['key'], f"{r['sec']:.0f}s", flush=True)
    summarize()


def summarize():
    meas = json.loads((ROOT / 'fpga/v80_sca/results/hwv6_board_v64_e12_250mhz/hwm_summary.json').read_text())['rounds'] \
        if (ROOT / 'fpga/v80_sca/results/hwv6_board_v64_e12_250mhz/hwm_summary.json').exists() else {}
    rows = []
    for key, label, cfg in CONFIGS:
        f = OUT / f'{key}.json'
        if not f.exists():
            continue
        r = json.loads(f.read_text())
        cuts = np.array(r['cuts']); flips = np.array(r['flips']); k = int((cuts >= 33000).sum()); p = k / NTRIAL
        c = cycles(flips, cfg['S'])
        t_round = e_max(c, E_ENG) / (MHZ * 1e3)
        prim, sec, P_r = tts(t_round, p)
        row = dict(key=key, label=label, cfg=cfg, p=p, successes=k, mean_flips=float(flips.mean()),
                   t_trial_ms=float(c.mean() / (MHZ * 1e3)), t_round_ms=t_round, P_round=P_r, primary_ms=prim, secondary_ms=sec)
        if key in meas:
            mm = meas[key]
            row['measured'] = dict(p=mm['p'], P_round=mm['P_round'], t_round_ms=mm['t_round_ms'], primary_ms=mm['tts99_ms'],
                                   secondary_ms=mm['tts99_ind_ms'])
        rows.append(row)
    (OUT / 'summary.json').write_text(json.dumps(rows, indent=1) + '\n')
    for r in rows:
        mm = r.get('measured')
        print(f"{r['key']:9s} p={r['p']:.3f} t_round={r['t_round_ms']:.4f} prim={r['primary_ms']:.4f} sec={r['secondary_ms']:.4f}"
              + (f" | board p={mm['p']:.3f} t_round={mm['t_round_ms']:.4f} prim={mm['primary_ms']:.4f} sec={mm['secondary_ms']:.4f}" if mm else '')
              + f"  {r['label']}")


if __name__ == '__main__':
    if sys.argv[1] == 'run':
        run(int(sys.argv[2]) if len(sys.argv) > 2 else 16)
    else:
        summarize()
