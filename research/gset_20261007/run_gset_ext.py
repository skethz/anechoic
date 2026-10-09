"""PROTOCOL_AMENDMENT1.md runner (grid extension; run_gset.py is unchanged and imported).
Usage: python3 run_gset_ext.py A2 | C2
  A2: pilot the 16 extension points per (instance, engine rule, S) (64 runs each, seed [SEED, 4, g, method, S index, 16+i]),
      select over the union (original 16 + extension 16) by pilot mean cut (ties -> lower union index, original first);
      final: reuse the original 256-run final if the union selection is an original point, else 256 new runs
      (seed [SEED, 5, g, method, S index]).                       -> results/ext_budget/G<g>/<rule>_S<S>.json
  C2: TTS confirmations from the union finals: S*_1, S*_12 as in PROTOCOL.md, 1024 runs, seed [SEED, 6, g, objective],
      shared by the four rules (paired).                            -> results/ext_conf/G<g>/<rule>_<E1|E12>.json
"""
import os

for _v in ('NUMBA_NUM_THREADS', 'OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_v] = '1'

import itertools  # noqa: E402
import json  # noqa: E402
import math  # noqa: E402
import multiprocessing as mp  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402

import numpy as np  # noqa: E402

import run_gset as R  # noqa: E402

OUT = R.OUT
EXT = {   # extension blocks, values in units of sigma (kappa dimensionless, lambda in K2000 units)
    'P+': dict(plain_q=(1.8, 2.52, 3.6, 4.32), plain_T0=(1.125, 1.575, 2.25, 3.6), q=(2.16, 2.88), T0=(0.9, 1.125),
               jv=(-0.72, -0.36, -0.18, -0.09)),
    'Q+': dict(plain_q=(2.88, 3.6, 4.32, 5.76), plain_T0=(1.35, 1.8, 2.7, 3.6), q=(2.16, 2.88), T0=(0.9, 1.8),
               jv=(-0.36, -0.18, -0.09, 0.09)),
    'M+-': dict(plain_q=(0.135, 0.18, 0.27, 0.36), plain_T0=(1.125, 1.575, 1.8, 2.7), q=(0.18, 0.36), T0=(1.35, 1.8),
                jv=(-0.18, -0.09, -0.045, 0.045)),
}
KAPPA_EXT = (0.25, 0.5, 1.0, 2.0)
LAMBDA_EXT = (0.25, 0.5, 1.0, 2.0)


def grid_ext(method, G, g, S):
    c = EXT[R.wclass(g)]; tfr = R.SHARED[R.wclass(g)]['tfin']; a = G.sigma; tf = tfr * a; P = itertools.product
    out = []
    if method == 'SCA':
        for q, t in P(c['plain_q'], c['plain_T0']):
            out.append(dict(family='plain', q=q * a, T0=t * a, tfin=tf, S=S, rel=dict(q=q, T0=t, tfin=tfr)))
    elif method == 'TEC':
        for j, q, t in P(c['jv'], c['q'], c['T0']):
            out.append(dict(family='tec', jv=j * a, q=q * a, T0=t * a, tfin=tf, S=S, rel=dict(jv=j, q=q, T0=t, tfin=tfr)))
    elif method == 'Onsager-kT':
        for k, q, t in P(KAPPA_EXT, c['q'], c['T0']):
            out.append(dict(family='tecT', kappa=k, q=q * a, T0=t * a, ramp=True, tfin=tf, S=S, rel=dict(kappa=k, q=q, T0=t, tfin=tfr)))
    elif method == 'Onsager-online':
        for lam, q, t in P(LAMBDA_EXT, c['q'], c['T0']):
            out.append(dict(family='onsager', lam=lam * G.lam_factor, q=q * a, T0=t * a, ramp=True, tfin=tf, S=S,
                            rel=dict(lam_k2000=lam, lam_factor=G.lam_factor, q=q, T0=t, tfin=tfr)))
    assert len(out) == 16
    old = R.grid(method, G, g, S)
    assert not any(o['rel'] == n['rel'] for o in old for n in out), 'extension duplicates an original point'
    return out


def ext_path(g, method, S):
    return OUT / 'ext_budget' / f'G{g}' / f"{method.replace(' ', '_')}_S{S}.json"


def job_ext(args):
    g, method, S = args
    path = ext_path(g, method, S)
    if path.exists():
        return None
    t0 = time.time()
    G = R.graph(g); mi = R.METHODS.index(method); si = R.S_LIST.index(S); tgt = R.BKV[f'G{g}']['target']
    orig = json.loads(R.budget_path(g, method, S).read_text())
    gr = grid_ext(method, G, g, S)
    pilot = []
    for i, cfg in enumerate(gr):
        cuts, flips, finite = R.run_method(G, cfg, R.B_PILOT, np.random.SeedSequence([R.SEED, 4, g, mi, si, 16 + i]))
        pilot.append(dict(cfg=cfg, valid=True, mean_cut=float(cuts.mean()), p_target=float((cuts >= tgt).mean()),
                          mean_flips=float(np.mean(flips))))
    union = orig['pilot'] + pilot
    sel = max(range(len(union)), key=lambda i: (union[i]['mean_cut'], -i))
    rec = dict(instance=f'G{g}', method=method, S=S, wclass=R.wclass(g), ext_pilot=pilot, union_selected=sel,
               union_selected_cfg=union[sel]['cfg'], from_extension=bool(sel >= 16))
    if sel < 16:
        assert sel == orig['selected']
        rec.update(final=orig['final'], cuts=orig['cuts'], flips=orig['flips'], final_source='original final (seed [SEED,2,...])')
    else:
        cuts, flips, finite = R.run_method(G, union[sel]['cfg'], R.B_FINAL, np.random.SeedSequence([R.SEED, 5, g, mi, si]))
        rec.update(final=R.stats(cuts, flips, finite, g, S), cuts=cuts.tolist(), flips=np.asarray(flips).tolist(),
                   final_source='extension final (seed [SEED,5,...])')
    rec['sec'] = time.time() - t0
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix('.tmp'); tmp.write_text(json.dumps(rec) + '\n'); tmp.rename(path)
    f = rec['final']
    return (f"EXT G{g:<3d} {method:15s} S={S:<5d} union sel={sel}{'*' if sel >= 16 else ' '} pilot={union[sel]['mean_cut']:.1f} | "
            f"final mean={f['mean_cut']:.1f} p={f['p_target']:.3f} mcs99={f['mcs99']:.0f} tts1={f['tts1_ms']:.4g}ms "
            f"tts12={f['tts12_ms']:.4g}ms ({rec['sec']:.0f}s)")


def ext_conf_path(g, method, obj):
    return OUT / 'ext_conf' / f'G{g}' / f"{method.replace(' ', '_')}_{obj}.json"


def job_ext_conf(args):
    g, method, obj = args
    path = ext_conf_path(g, method, obj)
    if path.exists():
        return None
    t0 = time.time()
    best = None
    for S in R.S_LIST:
        d = json.loads(ext_path(g, method, S).read_text())
        v = d['final']['tts1_ms' if obj == 'E1' else 'tts12_ms']
        if math.isfinite(v) and (best is None or v < best[0]):
            best = (v, S, d['union_selected_cfg'])
    rec = dict(instance=f'G{g}', method=method, objective=obj)
    if best is None:
        rec.update(S=None, cfg=None, final=None, note='TTS infinite at every budget in the union finals: not confirmed')
    else:
        _, S, cfg = best
        cuts, flips, finite = R.run_method(R.graph(g), cfg, R.B_CONF, np.random.SeedSequence([R.SEED, 6, g, R.OBJECTIVES.index(obj)]))
        rec.update(S=S, cfg=cfg, tts_selection_256=best[0], final=R.stats(cuts, flips, finite, g, S), cuts=cuts.tolist(),
                   flips=np.asarray(flips).tolist())
    rec['sec'] = time.time() - t0
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix('.tmp'); tmp.write_text(json.dumps(rec) + '\n'); tmp.rename(path)
    f = rec['final']
    return (f"EXTCONF G{g:<3d} {method:15s} {obj:3s} S={rec['S']} p={f['p_target']:.3f} tts1={f['tts1_ms']:.4g}ms "
            f"tts12={f['tts12_ms']:.4g}ms ({rec['sec']:.0f}s)" if f else f"EXTCONF G{g} {method} {obj}: not confirmed")


def main():
    part = sys.argv[1]
    order_g = R.SUBSET + tuple(x for x in R.ALL if x not in R.SUBSET)
    if part == 'A2':
        jobs = [(g, m, S) for g in order_g for m in R.ENGINE_RULES for S in R.S_LIST]
        sub = [j for j in jobs if j[0] in R.SUBSET]; rest = [j for j in jobs if j[0] not in R.SUBSET]
        jobs = sorted(sub, key=lambda j: -R.graph_n(j[0]) * j[2]) + sorted(rest, key=lambda j: -R.graph_n(j[0]) * j[2])
        fn = job_ext
    elif part == 'C2':
        jobs = [(g, m, o) for g in order_g for m in R.ENGINE_RULES for o in R.OBJECTIVES]; fn = job_ext_conf
    else:
        raise SystemExit(__doc__)
    print(f'part {part}: {len(jobs)} jobs, {R.WORKERS} workers, start {time.strftime("%Y-%m-%d %H:%M:%S")}', flush=True)
    t0 = time.time(); done = 0
    with mp.Pool(R.WORKERS) as pool:
        for msg in pool.imap_unordered(fn, jobs, chunksize=1):
            done += 1
            if msg:
                print(f'[{done}/{len(jobs)} {time.time() - t0:.0f}s] {msg}', flush=True)
    print(f'part {part} complete {time.strftime("%Y-%m-%d %H:%M:%S")} ({time.time() - t0:.0f}s)', flush=True)


if __name__ == '__main__':
    main()
