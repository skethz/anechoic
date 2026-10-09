"""Addendum 4 validation runs: every baseline at its paper's settings against a result its paper reports
(PROTOCOL_ADDENDUM4.md lists the targets and criteria). Usage: python3 validate_a4.py run [W]; python3 validate_a4.py summary
Results: validation_a4/<job>.json. Seeds: SeedSequence([20261008, 60, job code, chunk]); real Neal: integer seeds
crc32(job) mod 2^31 (dwave-samplers rejects seeds >= 2^31 despite its message)."""
import os

for _v in ('NUMBA_NUM_THREADS', 'OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_v] = '1'

import json  # noqa: E402
import math  # noqa: E402
import subprocess  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
import zlib  # noqa: E402
from pathlib import Path  # noqa: E402

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
RES = HERE.parent
os.environ.setdefault('NUMBA_CACHE_DIR', str(HERE / '.numba_cache'))
sys.path.insert(0, str(HERE))
import numpy as np  # noqa: E402

OUT = HERE / 'validation_a4'
NEAL_PY = Path('/scratch/USER/anechoic_cpu_20261008/venv_neal/bin/python3')
GSET = RES / 'gset_20261007' / 'data'
K2000_RUD = RES.parent / 'fpga' / 'v80_snowball' / 'data' / 'WK2000_1.rud'
TABLE4 = ('G1', 'G2', 'G6', 'G7', 'G10', 'G11', 'G12', 'G13', 'G14', 'G19', 'G20')
APC_T3 = {'G22': (31, 1, -6545.5), 'G30': (2, 2, -6638.9), 'G32': (26, 3, -2726.3), 'G35': (27, 1, -3430.1)}


def seed(job, chunk=0):
    return np.random.SeedSequence([20261008, 60, zlib.crc32(job.encode()), chunk])


def mcp_problem(name):
    import k2000
    import run_a2 as RA
    if name == 'K2000':
        return k2000.problem()
    g = int(name[1:])
    if g <= 20:
        return RA.problem('mcp', g)
    import problems as PB
    lines = (GSET / name).read_text().split('\n')
    n, m = (int(x) for x in lines[0].split())
    a = np.array([ln.split() for ln in lines[1:] if ln.strip()], dtype=np.int64)
    i, j, w = a[:, 0] - 1, a[:, 1] - 1, a[:, 2]
    rows = np.concatenate([i, j]); cols = np.concatenate([j, i]); vals = -np.concatenate([w, w]).astype(np.float64)
    return PB.Problem(name, 'mcp', n, rows, cols, vals, 0.0, np.zeros(n), math.sqrt(2.0 * len(i) / n),
                      dict(ei=i, ej=j, w=w, W=int(w.sum()), ref=1, target=1, m=len(i)))


def cut_of(P, s):
    s = np.asarray(s)
    return ((s[:, P.meta['ei']] != s[:, P.meta['ej']]) * P.meta['w'][None, :]).sum(1)


def ising_energy(P, s):
    """APC convention H = -1/2 sum J s s - sum h s (Max-Cut: H = sum_{edges} w s_i s_j)."""
    s = np.asarray(s, np.float64)
    return (s[:, P.meta['ei']] * s[:, P.meta['ej']] * P.meta['w'][None, :]).sum(1)


def neal_defaults(P):
    import settings_a3 as ST
    C = ST.instance_constants(P)
    return 1.0 / C['neal']['hot_beta'], 1.0 / C['neal']['cold_beta'], C


# ------------------------------------------------------------------------------------------------ jobs
def job(args):
    name, kind, params = args
    path = OUT / f'{name}.json'
    if path.exists():
        return None
    t0 = time.time()
    import solvers as SV
    import solvers_a3 as SV3
    import solvers_a3b as SV3B
    import solvers_a4 as SV4
    rec = dict(job=name, kind=kind, params=params)
    if kind == 'sa_vs_neal':
        P = mcp_problem(params['inst']); S = params['S']; B = params['B']
        T0, T1, _ = neal_defaults(P)
        s = SV.sa(P, B, dict(family='SA', T0=T0, T1=T1, S=S), np.random.default_rng(seed(name)))
        mine = cut_of(P, s)
        src = str(K2000_RUD if params['inst'] == 'K2000' else GSET / params['inst'])
        out = subprocess.run([str(NEAL_PY), str(HERE / 'neal_ref.py'), src, str(S), str(B), str(zlib.crc32(name.encode()) % 2 ** 31)],
                             capture_output=True, text=True, check=True).stdout
        neal = json.loads(out)
        rec.update(T0=T0, T1=T1, ours=mine.tolist(), neal=neal['cuts'], neal_beta_range=neal['beta_range'],
                   ours_beta_range=[1 / T0, 1 / T1])
    elif kind == 'statica':
        P = mcp_problem('K2000'); B = params['B']
        cfg = dict(family='plain', q=params['q'], T0=params['T0'], tfin=params['tfin'], S=params['S'])
        s, flips, _ = SV.run_method(P, cfg, B, seed(name))
        rec.update(cuts=cut_of(P, s).tolist(), mean_flips=float(np.mean(flips)))
    elif kind == 'tec':
        P = mcp_problem('K2000'); B = params['B']
        cfg = dict(family='tec_seq', jv=params['jv'], T0=100.0, tfin=0.1, S=params['S'])
        s, flips, Etr = SV4.tec_seq(P, B, cfg, np.random.default_rng(seed(name, params['chunk'])), trace=True)
        W = int(P.meta['w'].sum())
        cut_tr = (W - Etr) / 2.0                        # spatial energy E = sum_edges w s s = W - 2 cut
        rec.update(final_cuts=cut_of(P, s).tolist(), cut_trace_mean=cut_tr.mean(0).tolist(),
                   first_pass=[int(np.argmax(c >= 31670)) + 1 if (c >= 31670).any() else None for c in cut_tr],
                   check_trace_final=bool(np.allclose(cut_tr[:, -1], cut_of(P, s))), mean_flips=float(np.mean(flips)))
    elif kind == 'tec_sync':
        P = mcp_problem('K2000'); B = params['B']
        cfg = dict(family='tec_sig', jv=params['jv'], q=0.0, T0=100.0, tfin=0.1, cb=2.0, S=params['S'])
        s, flips, _ = SV3.run_method(P, cfg, B, seed(name))
        rec.update(final_cuts=cut_of(P, s).tolist(), mean_flips=float(np.mean(flips)))
    elif kind == 'apc_t3':
        P = mcp_problem(params['inst']); B = params['B']
        lam = float(np.linalg.eigvalsh(-P.dense())[-1])
        cfg = dict(q_init=params['qi'] / 64 * lam, q_final=params['qf'] / 64 * lam, T0=10.0, tfin=0.1, S=1000)
        best, final, flips = SV4.sca_sig_q(P, B, cfg, np.random.default_rng(seed(name)))
        rec.update(lam=lam, energy_best=ising_energy(P, best).tolist(), energy_final=ising_energy(P, final).tolist())
    elif kind == 'apc_fig3':
        P = mcp_problem(params['inst']); B = params['B']
        lam = float(np.linalg.eigvalsh(-P.dense())[-1])
        cfg = dict(family='apc_sig', q_reset=lam / 2, r_q=params['r_q'], q_lim=0.0, T0=10.0, tfin=0.1, cb=1.0, S=params['S'])
        best, flips, _ = SV3B.run_method(P, cfg, B, seed(name))
        rec.update(lam=lam, energy_best=ising_energy(P, best).tolist())
    elif kind == 'reaim':
        import run_a2 as RA
        import settings_a3 as ST
        prob, inst = params['prob'], params['inst']
        P = RA.problem(prob, inst)
        cfg = ST.baseline_cfgs(prob, P, 'ReAIM ASA', params['S'], ST.instance_constants(P))[0]
        s, _, _ = SV3B.run_method(P, cfg, params['B'], seed(name))
        ev = P.evaluate(s)
        rec.update(values=ev['value'].tolist(), feasible=ev['feasible'].tolist(), cfg=cfg)
    elif kind == 'asb':
        P = mcp_problem('K2000'); B = params['B']
        import settings_a3 as ST
        C = ST.instance_constants(P)
        cfg = dict(family='aSB', dt=params['dt'], M=params['M'], xi=C['rms_J'] / C['sd_J'], S=params['S'])
        s, finite = SV.asb(P, B, cfg, np.random.default_rng(seed(name)))
        rec.update(cuts=cut_of(P, s).tolist(), finite=np.asarray(finite).tolist(), xi0=0.7 / (C['sd_J'] * math.sqrt(P.N)))
    elif kind == 'sb21':
        P = mcp_problem('K2000'); B = params['B']
        cfg = dict(family=params['family'], dt=params['dt'], xi=1.0, S=params['S'])
        s = SV.sb(P, B, cfg, np.random.default_rng(seed(name, params['chunk'])))
        rec.update(cuts=cut_of(P, s).tolist())
    rec['sec'] = time.time() - t0
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix('.tmp'); tmp.write_text(json.dumps(rec) + '\n'); tmp.rename(path)
    return f'{name} ({rec["sec"]:.0f}s)'


def jobs():
    J = []
    for inst in ('G1', 'G14', 'G22', 'K2000'):
        J.append((f'sa_vs_neal_{inst}', 'sa_vs_neal', dict(inst=inst, S=1000, B=256)))
    for inst in TABLE4:
        J.append((f'neal_t4_{inst}', 'sa_vs_neal', dict(inst=inst, S=1000, B=240)))
    J.append(('statica_long', 'statica', dict(S=1560, q=4.0, T0=40.0, tfin=5.0, B=1024)))
    J.append(('statica_short', 'statica', dict(S=560, q=4.0, T0=30.0, tfin=5.0, B=1024)))
    for jv in (0.0, 30.0, 6.0, -6.0):
        for c in range(8):
            J.append((f'tec_jv{jv:g}_c{c}', 'tec', dict(jv=jv, S=3000, B=8, chunk=c)))
    for jv in (0.0, 30.0):
        J.append((f'tec_sync_jv{jv:g}', 'tec_sync', dict(jv=jv, S=3000, B=64)))
    for inst, (qi, qf, _) in APC_T3.items():
        J.append((f'apc_t3_{inst}', 'apc_t3', dict(inst=inst, qi=qi, qf=qf, B=128)))
        J.append((f'apc_fig3_{inst}', 'apc_fig3', dict(inst=inst, r_q=0.45, S=1000, B=128)))
    for inst in TABLE4:
        J.append((f'reaim_t4_{inst}', 'reaim', dict(prob='mcp', inst=int(inst[1:]), S=4096, B=240)))
    for g in (1, 2, 3, 4, 5, 14, 15, 16, 17):
        J.append((f'reaim_t5_G{g}', 'reaim', dict(prob='gpp', inst=g, S=4096, B=240)))
    for t in ('gr17', 'gr21', 'gr24', 'fri26', 'bayg29', 'bays29'):
        J.append((f'reaim_t6_{t}', 'reaim', dict(prob='tsp', inst=t, S=8192, B=240)))
    J.append(('asb_k2000_fig2b', 'asb', dict(dt=0.9, M=2, S=186, B=256)))
    for fam in ('bSB', 'dSB'):
        for S, nch, B in ((100, 1, 256), (1000, 4, 64), (10000, 16, 16)):
            for c in range(nch):
                J.append((f'sb21_{fam}_S{S}_c{c}', 'sb21', dict(family=fam, dt=1.0, S=S, B=B, chunk=c)))
    return J


def run(W):
    import multiprocessing as mp
    J = jobs()
    order = {'sb21': 0, 'tec': 1, 'reaim': 2, 'sa_vs_neal': 3}
    J.sort(key=lambda j: (order.get(j[1], 9), -j[2].get('S', 0)))
    print(f'validation: {len(J)} jobs, {W} workers', flush=True)
    t0 = time.time()
    with mp.Pool(W) as pool:
        for k, msg in enumerate(pool.imap_unordered(job, J, chunksize=1)):
            if msg:
                print(f'[{k + 1}/{len(J)} {time.time() - t0:.0f}s] {msg}', flush=True)
    print(f'validation complete ({time.time() - t0:.0f}s)', flush=True)


if __name__ == '__main__':
    if sys.argv[1] == 'run':
        run(int(sys.argv[2]) if len(sys.argv) > 2 else 64)
