"""EXPLORATORY, pre-protocol, SYNTHETIC graphs only (no G-set instance): (a) wider plain-SCA scan of q, T0, T_fin in
units of sigma (calib.py found the optimum at its grid edges); (b) placement of each rule's own parameter (TEC J_v,
Onsager-kT kappa, Onsager-online lambda * d/N, APC) at plain's best shared parameters; (c) the online coefficient
lambda_eff * n_lin / (2T) expressed as kappa_equiv = coefficient / T. Output: calib2.json, calib2.log."""
import itertools
import json
import multiprocessing as mp

import numpy as np

import calib
import engine as E

PLUS = ('R800+', 'P800+', 'R2000+', 'P2000+', 'R1000+')
PM = ('R800+-', 'T800+-', 'P800+-', 'T2000+-')


def task(args):
    key, cfg_rel, B, seed = args
    g = calib.graph(key); a = g.sigma
    cfg = dict(cfg_rel)
    tfin = cfg.pop('tfin') * a
    for k in ('q', 'T0', 'jv', 'q_reset', 'q_lim'):
        if k in cfg:
            cfg[k] = cfg[k] * a
    if 'lam' in cfg:
        cfg['lam'] = cfg['lam'] * g.lam_factor
    rng = np.random.default_rng(seed)
    s0 = rng.choice(np.array([-1.0, 1.0], dtype=np.float32), size=(B, g.N))
    if cfg['family'] == 'onsager':
        s, fl, coef, band, T = E.run_traced(g, s0, cfg, rng, tfin)
        kq = [float(coef[t] / T[t - 1]) for t in (int(0.1 * len(T)), int(0.3 * len(T)), int(0.5 * len(T)), int(0.69 * len(T)))]
    else:
        s, fl = E.run(g, s0, cfg, rng, tfin); kq = None
    c = g.cut(s)
    return dict(graph=key, cfg=cfg_rel, mean_cut=float(c.mean()), sd=float(c.std()), max_cut=int(c.max()),
                mean_flips=float(fl.mean()), kappa_equiv_at_10_30_50_69pct=kq)


def main():
    S, B = 1000, 32
    jobs = []
    for k in PLUS:
        for i, (q, t0, tf) in enumerate(itertools.product((0.72, 1.08, 1.44, 2.16, 2.88, 4.32), (0.45, 0.9, 1.8), (0.112, 0.224, 0.448))):
            jobs.append((k, dict(family='plain', q=q, T0=t0, tfin=tf, S=S), B, 2000 + i))
    for k in PM:
        for i, (q, t0, tf) in enumerate(itertools.product((0.09, 0.18, 0.36, 0.54), (0.45, 0.9, 1.35, 1.8), (0.056, 0.112, 0.224))):
            jobs.append((k, dict(family='plain', q=q, T0=t0, tfin=tf, S=S), B, 2000 + i))
    with mp.Pool(6) as pool:
        rows = list(pool.imap_unordered(task, jobs, chunksize=2))
    best = {}
    for k in PLUS + PM:
        rr = sorted((r for r in rows if r['graph'] == k), key=lambda r: -r['mean_cut'])
        best[k] = rr[0]['cfg']
        mx = max(r['max_cut'] for r in rr)
        print(f"{k:8s} best-found {mx}", flush=True)
        for r in rr[:5]:
            c = r['cfg']
            print(f"   q/s {c['q']:.2f} T0/s {c['T0']:.2f} Tf/s {c['tfin']:.3f}: mean {r['mean_cut']:.1f} sd {r['sd']:.1f} "
                  f"({100 * r['mean_cut'] / mx:.2f}%) flips {r['mean_flips']:.0f}", flush=True)
    # (b) rule parameters at plain's best shared parameters
    jobs = []
    for k in PLUS + PM:
        b = best[k]; base = dict(q=b['q'], T0=b['T0'], tfin=b['tfin'], S=S)
        for x in (0.125, 0.25, 0.5, 1.0, 2.0, 4.0):
            jobs.append((k, dict(base, family='tecT', kappa=x, ramp=True), B, 3000))
            jobs.append((k, dict(base, family='onsager', lam=x, ramp=True), B, 3000))
        for x in (-0.36, -0.18, -0.09, -0.045, 0.045, 0.09, 0.18):
            jobs.append((k, dict(base, family='tec', jv=x), B, 3000))
        jobs.append((k, dict(base, family='plain'), B, 3000))
    with mp.Pool(6) as pool:
        rows2 = list(pool.imap_unordered(task, jobs, chunksize=1))
    for k in PLUS + PM:
        rr = [r for r in rows2 if r['graph'] == k]
        p0 = [r for r in rr if r['cfg']['family'] == 'plain'][0]['mean_cut']
        line = f"{k:8s} plain {p0:.1f} |"
        for fam, par in (('tecT', 'kappa'), ('onsager', 'lam'), ('tec', 'jv')):
            xs = sorted((r for r in rr if r['cfg']['family'] == fam), key=lambda r: r['cfg'][par])
            line += f" {fam}: " + " ".join(f"{r['cfg'][par]:g}:{r['mean_cut'] - p0:+.1f}" for r in xs) + " |"
        print(line, flush=True)
        ons = sorted((r for r in rr if r['cfg']['family'] == 'onsager'), key=lambda r: r['cfg']['lam'])
        print(f"          onsager kappa_equiv (coef/T at 10/30/50/69% of run): " +
              "; ".join(f"lam {r['cfg']['lam']:g}: " + ",".join(f"{v:.2f}" for v in r['kappa_equiv_at_10_30_50_69pct']) for r in ons),
              flush=True)
    (E.ROOT / 'calib2.json').write_text(json.dumps(dict(plain_scan=rows, best_plain=best, rule_params=rows2), indent=1) + '\n')


if __name__ == '__main__':
    main()
