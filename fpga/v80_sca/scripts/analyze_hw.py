#!/usr/bin/env python3
"""Pre-registered analysis of PROTOCOL_HW.md. Usage: analyze_hw.py <hw_result_dir> <J_results.json>"""
import json
import math
import sys
from pathlib import Path

TARGET = 33000
KEYS = {  # protocol key -> J selection tuple (mode, q, par, ramp, T0, S)
    'P1': ('plain', 4.0, 0.0, False, 30.0, 1560), 'P2': ('plain', 6.0, 0.0, False, 30.0, 1560),
    'P3': ('plain', 8.0, 0.0, False, 30.0, 1560), 'P4': ('plain', 8.0, 0.0, False, 30.0, 960),
    'T1': ('tec', 8.0, -4.0, False, 30.0, 1560), 'T2': ('tec', 8.0, -4.0, False, 30.0, 960),
    'T3': ('tec', 8.0, 4.0, False, 30.0, 1560),
    'O1': ('onsager', 8.0, 1.05, True, 12.0, 960), 'O2': ('onsager', 6.0, 0.9, True, 15.0, 560),
    'O3': ('onsager', 6.0, 0.9, True, 12.0, 560), 'O4': ('onsager', 6.0, 0.9, False, 12.0, 360),
}


def binom_cdf(k, n, p):
    if p <= 0: return 1.0
    if p >= 1: return 1.0 if k >= n else 0.0
    lp, lq = math.log(p), math.log1p(-p)
    return min(1.0, sum(math.exp(math.lgamma(n + 1) - math.lgamma(i + 1) - math.lgamma(n - i + 1) + i * lp + (n - i) * lq)
                        for i in range(0, k + 1)))


def clopper_pearson(k, n, a=0.05):
    def bisect(f):
        lo, hi = 0.0, 1.0
        for _ in range(60):
            mid = (lo + hi) / 2
            if f(mid): hi = mid
            else: lo = mid
        return (lo + hi) / 2
    lower = 0.0 if k == 0 else bisect(lambda p: 1 - binom_cdf(k - 1, n, p) >= a / 2)
    upper = 1.0 if k == n else bisect(lambda p: binom_cdf(k, n, p) <= a / 2)
    return lower, upper


def r99(p):
    if p >= 1: return 1.0
    if p <= 0: return math.inf
    return math.log(0.01) / math.log1p(-p)


def model_cycles(flips, S, corr):
    return 0.434 * flips - 0.98 * S + 1989 + (4 * S if corr else 0)


def load(path):
    return [json.loads(line) for line in Path(path).read_text().splitlines() if line.strip()]


def launches(rows):
    """Group trials by launch (consecutive rows sharing launch_cycles and launch_trials)."""
    out, i = [], 0
    while i < len(rows):
        n = rows[i]['launch_trials']
        out.append(rows[i:i + n]); i += n
    return out


def tts(p, t_run):
    return t_run * r99(p), t_run * math.ceil(r99(p))


def main():
    d, jpath = Path(sys.argv[1]), Path(sys.argv[2])
    manifest = json.loads((d / 'image_manifest.json').read_text())
    mhz = manifest['clock']['actual_mhz']
    jhold = {(r['mode'], r['q'], r['par'], r['ramp'], r['T0'], r['S']): r for r in json.loads(jpath.read_text())['holdout']}
    summary = {'clock_mhz': mhz, 'logic_uuid': manifest['logic_uuid'], 'verify': {}, 'cohort': {}, 'single': {}}

    for f in sorted((d / 'verify').glob('*.jsonl')):
        rows = load(f)
        summary['verify'][f.stem] = {'trials': len(rows), 'exact_pass': sum(r['exact_reference'] == 'PASS' for r in rows),
                                     'scores_agree': sum(r['scores_agree'] for r in rows)}

    per = {}
    for key, cfg in KEYS.items():
        f = d / 'cohort' / f'{key}.jsonl'
        if not f.exists(): continue
        rows = load(f)
        ok = [r['success'] and r['scores_agree'] for r in rows]
        k, n = sum(ok), len(rows)
        L = launches(rows)
        t_launch = [g[0]['launch_cycles'] / len(g) / (mhz * 1e3) for g in L]  # ms per trial, amortized per launch
        t_run = sum(t_launch) / len(t_launch)
        corr = cfg[0] == 'onsager'
        ratio = sum(g[0]['launch_cycles'] for g in L) / sum(model_cycles(r['flips'], r['steps'], corr) for r in rows)
        p = k / n; lo, hi = clopper_pearson(k, n)
        jr = jhold[cfg]; pj = jr['p']
        se = math.sqrt(p * (1 - p) / n + pj * (1 - pj) / jr['trials'])
        diff = (p - pj - 1.96 * se, p - pj + 1.96 * se)
        t99, t99c = tts(p, t_run)
        per[key] = {'ok': ok, 't_launch': t_launch}
        summary['cohort'][key] = {
            'config': cfg, 'trials': n, 'successes': k, 'p': p, 'p_ci95': [lo, hi],
            'integrity_failures': sum(not r['scores_agree'] for r in rows),
            'mean_cut': sum(r['cut'] for r in rows) / n, 'mean_flips': sum(r['flips'] for r in rows) / n,
            't_run_ms': t_run, 'tts99_ms': t99, 'tts99_whole_run_ms': t99c,
            'projected_device_tts_ms': {E: t_run * math.ceil(math.ceil(r99(p)) / E) if p > 0 else math.inf for E in (1, 4, 8, 16, 32)},
            'j_model': {'p': pj, 't_ms_model': jr['t_ms'], 'mean_flips': jr['mean_flips']},
            'p_diff_ci95': list(diff), 'reproduced': abs(p - pj) <= 0.05 or diff[0] <= 0 <= diff[1],
            'measured_over_model_cycles': ratio}

    # Single-run latency and a descriptive per-trial cycle fit (cycles = a*flips + b*S + c).
    X, y = [], []
    for f in sorted((d / 'single').glob('*.jsonl')):
        rows = load(f)
        lat = sorted(r['launch_device_ms'] for r in rows)
        summary['single'][f.stem] = {'trials': len(rows), 'mean_ms': sum(lat) / len(lat), 'median_ms': lat[len(lat) // 2],
                                     'max_ms': lat[-1], 'host_mean_ms': sum(r['launch_host_ms'] for r in rows) / len(rows)}
        for r in rows:
            if not r['load_J']: X.append((r['flips'], r['steps'], 1.0)); y.append(r['launch_cycles'])
    if len(X) >= 3:
        # Normal equations for 3 parameters.
        A = [[sum(a[i] * a[j] for a in X) for j in range(3)] for i in range(3)]
        b = [sum(a[i] * v for a, v in zip(X, y)) for i in range(3)]
        for i in range(3):  # Gaussian elimination
            piv = A[i][i]
            for j in range(i + 1, 3):
                m = A[j][i] / piv
                A[j] = [aj - m * ai for aj, ai in zip(A[j], A[i])]; b[j] -= m * b[i]
        c = [0.0] * 3
        for i in (2, 1, 0):
            c[i] = (b[i] - sum(A[i][j] * c[j] for j in range(i + 1, 3))) / A[i][i]
        summary['cycle_fit'] = {'cycles_per_flip': c[0], 'cycles_per_step': c[1], 'per_run_constant': c[2], 'points': len(X)}

    # Primary comparison: O1 vs post-hoc best plain and best TEC (E=1, Eq. 9), bootstrap over trials and launches.
    def best(prefix, src):
        vals = {k: v for k, v in src.items() if k.startswith(prefix)}
        return min(vals.items(), key=lambda kv: kv[1]) if vals else (None, math.inf)
    point = {k: v['tts99_ms'] for k, v in summary['cohort'].items()}
    if 'O1' in point:
        import numpy as np
        rng = np.random.default_rng(20261004)
        B = 10000
        bt = {}
        for k, v in per.items():  # resampling 0/1 trials with replacement == binomial draw at the observed p
            n = len(v['ok']); s = rng.binomial(n, sum(v['ok']) / n, size=B)
            tl = np.array(v['t_launch']); t = tl[rng.integers(0, len(tl), size=(B, len(tl)))].mean(1)
            pb = s / n
            with np.errstate(divide='ignore'):
                r = np.where(pb >= 1, 1.0, np.where(pb <= 0, np.inf, np.log(0.01) / np.log1p(-np.minimum(pb, 1 - 1e-12))))
            bt[k] = t * r
        boots = {}
        for fam, pre in (('plain', 'P'), ('tec', 'T')):
            fam_min = np.min(np.stack([bt[k] for k in bt if k.startswith(pre)]), axis=0)
            boots[fam] = list(fam_min / bt['O1'])
        comp = {}
        for fam, pre in (('plain', 'P'), ('tec', 'T')):
            bk, bv = best(pre, point)
            xs = sorted(boots[fam])
            comp[fam] = {'best_key': bk, 'best_tts99_ms': bv, 'speedup_O1': bv / point['O1'],
                         'speedup_ci95': [xs[249], xs[9749]]}
        summary['primary'] = comp

    (d / 'hw_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(f"clock {mhz:.3f} MHz  uuid {summary['logic_uuid']}")
    for k, v in summary['verify'].items():
        print(f"verify {k}: exact {v['exact_pass']}/{v['trials']}")
    print('key  p_hw (95% CI)          p_J    repro  t_run ms  TTS99 ms  meas/model cycles')
    for k, v in summary['cohort'].items():
        print(f"{k:3s}  {v['p']:.3f} ({v['p_ci95'][0]:.3f}-{v['p_ci95'][1]:.3f})  {v['j_model']['p']:.3f}  "
              f"{'yes' if v['reproduced'] else 'NO ':3s}    {v['t_run_ms']:.4f}   {v['tts99_ms']:.4f}   {v['measured_over_model_cycles']:.3f}")
    for k, v in summary['single'].items():
        print(f"single {k}: mean {v['mean_ms']:.4f} ms, median {v['median_ms']:.4f}, max {v['max_ms']:.4f}")
    if 'cycle_fit' in summary: print('cycle fit', summary['cycle_fit'])
    if 'primary' in summary: print('primary', json.dumps(summary['primary']))


if __name__ == '__main__':
    main()
