"""Ablation sweep (PROTOCOL_ABL.md): Onsager SCA vs TEC-T, TEC-sched, APC-SCA (and re-costed plain / TEC-const)."""
import itertools
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parent / '../theory_ideas_20261003'))
import sweeps_abc as m  # noqa: E402

ROOT = Path(__file__).resolve().parent
J_RESULTS = Path(str(__import__('pathlib').Path(__file__).resolve().parent / '../theory_ideas_20261003/J_results.json'))
E_LIST = (1, 4, 16)
MHZ = 250.0


def t_ms(flips, S):
    return (0.2754 * flips + 31.06 * S + 753) / (MHZ * 1e3)


def r99(p):
    if p >= 1: return 1.0
    if p <= 0: return math.inf
    return math.log(0.01) / math.log1p(-p)


def dev_tts(t, p, E):
    R = r99(p)
    return math.inf if not math.isfinite(R) else t * math.ceil(math.ceil(R) / E)


def wilson_lower(k, n, z=1.96):
    p = k / n; d = 1 + z * z / n
    return max(0.0, (p + z * z / (2 * n) - z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))) / d)


def ramp_factor(t, S, ramp):
    return 1.0 if (not ramp or t < 0.7 * S) else (S - 1 - t) / max(1.0, S - 1 - 0.7 * S)


def run(J, s0, cfg, rng, ctable=None):
    """Clipped parallel SCA with the family's temporal term. Returns final spins and flips per run."""
    fam = cfg['family']; S = cfg['S']; T = m.sched(cfg['T0'], S)
    s = s0.astype(np.float32).copy(); h = s @ J
    B, N = s.shape
    flips = np.zeros(B); s_prev = None; c_prev = None; T_prev = None
    q = np.full((B, N), float(cfg.get('q', 0.0)), np.float32)
    if fam == 'apc':
        q[:] = cfg['q_reset']
    for t, Tt in enumerate(T):
        field = h
        if s_prev is not None:
            if fam == 'onsager':
                field = h - cfg['lam'] * ramp_factor(t, S, cfg['ramp']) * c_prev * s_prev
            elif fam == 'onsager_table':
                field = h - cfg['lam'] * ramp_factor(t, S, cfg['ramp']) * ctable[t] * s_prev
            elif fam == 'tecT':
                field = h - cfg['kappa'] * T_prev * ramp_factor(t, S, cfg['ramp']) * s_prev
            elif fam == 'tecS':
                jv = cfg['j_lo'] + (cfg['j_hi'] - cfg['j_lo']) * (Tt - m.T_FIN) / max(1e-9, cfg['T0'] - m.T_FIN)
                field = h + jv * s_prev
            elif fam == 'tec':
                field = h + cfg['jv'] * s_prev
        z = s * field + q
        flip = np.clip(z / (4 * Tt) + 0.5, 0, 1) < rng.random(s.shape, dtype=np.float32)
        if fam in ('onsager',):
            c_prev = ((z > -2 * Tt) & (z < 2 * Tt)).sum(1, keepdims=True) / (2 * Tt)
        if fam == 'apc':
            q = np.where(flip, np.float32(cfg['q_reset']), np.maximum(q * np.float32(cfg['r_q']), np.float32(cfg['q_lim'])))
        d = np.where(flip, -2 * s, 0).astype(np.float32)
        s_prev = s.copy(); T_prev = Tt
        s += d; h += d @ J; flips += flip.sum(1)
    return s, flips


def evaluate(J, sumw, cfgs, B, seeds, ctables=None, tag=''):
    out = []
    for i, cfg in enumerate(cfgs):
        rng = np.random.default_rng(seeds[i]); t1 = time.perf_counter()
        s0 = rng.choice(np.array([-1.0, 1.0], dtype=np.float32), size=(B, len(J)))
        s, fl = run(J, s0, cfg, rng, None if ctables is None else ctables.get(i))
        c = m.cut(J, sumw, s); k = int((c >= m.TARGET).sum())
        tr = t_ms(fl.mean(), cfg['S'])
        r = dict(cfg, trials=B, successes=k, p=k / B, mean_cut=float(c.mean()), mean_flips=float(fl.mean()), t_ms=tr,
                 sec=time.perf_counter() - t1)
        out.append(r)
        print(tag, json.dumps({kk: (round(v, 4) if isinstance(v, float) else v) for kk, v in r.items() if kk not in ('sec',)}), flush=True)
    return out


def grid_new():
    g = []
    for q, kap, ramp, T0, S in itertools.product((6.0, 8.0), (0.3, 0.45, 0.6, 0.75), (False, True), (12.0, 15.0, 20.0), (560, 960)):
        g.append(dict(family='tecT', q=q, kappa=kap, ramp=ramp, T0=T0, S=S))
    for jl, jh, q, T0, S in itertools.product((-8.0, -4.0, -2.0), (0.0, 4.0, 8.0, 16.0), (4.0, 8.0), (20.0, 30.0), (960, 1560)):
        g.append(dict(family='tecS', j_lo=jl, j_hi=jh, q=q, T0=T0, S=S))
    for qr, rq, ql, T0, S in itertools.product((8.0, 16.0, 32.0), (0.9, 0.97), (2.0, 4.0), (15.0, 20.0, 30.0), (560, 960)):
        g.append(dict(family='apc', q_reset=qr, r_q=rq, q_lim=ql, T0=T0, S=S))
    return g


def j_pilot_as_cfgs():
    """Reuse the J pilot (256 runs each), re-costed with the v5 model."""
    d = json.loads(J_RESULTS.read_text())
    out = []
    for r in d['pilot']:
        fam = {'plain': 'plain', 'onsager': 'onsager', 'tec': 'tec'}[r['mode']]
        cfg = dict(family=fam, q=r['q'], T0=r['T0'], S=r['S'])
        if fam == 'onsager': cfg.update(lam=r['par'], ramp=r['ramp'])
        if fam == 'tec': cfg.update(jv=r['par'])
        out.append(dict(cfg, trials=r['trials'], successes=r['successes'], p=r['p'], mean_flips=r['mean_flips'],
                        t_ms=t_ms(r['mean_flips'], r['S'])))
    return out


def select(rows):
    sel = {}
    for fam in sorted({r['family'] for r in rows}):
        for E in E_LIST:
            best = min((r for r in rows if r['family'] == fam),
                       key=lambda r: dev_tts(r['t_ms'], wilson_lower(r['successes'], r['trials']), E))
            sel[f'{fam}_E{E}'] = {k: v for k, v in best.items() if k in ('family', 'q', 'T0', 'S', 'lam', 'ramp', 'jv', 'kappa',
                                                                          'j_lo', 'j_hi', 'q_reset', 'r_q', 'q_lim')}
    return sel


def main():
    J, sumw = m.load()
    stage = sys.argv[1] if len(sys.argv) > 1 else 'all'
    pilot_path = ROOT / 'pilot_new.json'
    if stage in ('pilot', 'all'):
        cfgs = grid_new()
        seeds = np.random.SeedSequence(50001).spawn(len(cfgs))
        pilot = evaluate(J, sumw, cfgs, 256, seeds, tag='PILOT')
        pilot_path.write_text(json.dumps(pilot, indent=1) + '\n')
    if stage in ('holdout', 'all'):
        pilot = json.loads(pilot_path.read_text()) + j_pilot_as_cfgs()
        sel = select(pilot)
        distinct = []
        for v in sel.values():
            if v not in distinct: distinct.append(v)
        # Onsager-table arm: mean online c(t) from a separate 256-run pilot (seed 50003) for each Onsager-selected config
        ons = [v for v in distinct if v['family'] == 'onsager']
        tabs = []
        for v in ons:
            rng = np.random.default_rng(50003)
            s0 = rng.choice(np.array([-1.0, 1.0], dtype=np.float32), size=(256, len(J)))
            ct = np.zeros(v['S']); S = v['S']; T = m.sched(v['T0'], S)
            # record mean c(t-1) used at step t by replaying the online run
            s = s0.copy(); h = s @ J; sp = None; cp = None
            for t, Tt in enumerate(T):
                field = h if sp is None else h - v['lam'] * ramp_factor(t, S, v['ramp']) * cp * sp
                if cp is not None: ct[t] = float(cp.mean())
                z = s * field + v['q']
                flip = np.clip(z / (4 * Tt) + 0.5, 0, 1) < rng.random(s.shape, dtype=np.float32)
                cp = ((z > -2 * Tt) & (z < 2 * Tt)).sum(1, keepdims=True) / (2 * Tt)
                d = np.where(flip, -2 * s, 0).astype(np.float32); sp = s.copy(); s = s + d; h = h + d @ J
            tabs.append(ct)
        table_cfgs = [dict(v, family='onsager_table') for v in ons]
        hcfgs = distinct + table_cfgs
        ctables = {len(distinct) + i: tabs[i] for i in range(len(tabs))}
        seeds = [np.random.SeedSequence(50002)] * len(hcfgs)          # paired: identical starts and streams
        hold = evaluate(J, sumw, hcfgs, 1024, seeds, ctables=ctables, tag='HOLD')
        (ROOT / 'holdout.json').write_text(json.dumps(dict(selection=sel, holdout=hold), indent=1) + '\n')


if __name__ == '__main__':
    main()
