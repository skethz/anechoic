"""Paper figure: DMFT solver (gpu_v2_M1000000_it32_S560_lam0_S560_lam0.7.json) against direct simulation on K2000 for the
same schedules (T 30 -> 5 geometric, S = 560, q = 4, lambda = 0 and 0.7 constant). Simulation: 256 runs each, float model.
(a) E[s u](t) with u = h/sqrt(N);  (b) echo coefficient c(t) = n_lin/(2T) (DMFT: N * P_lin / (2T)).
DMFT curves are drawn to step 500: closer to freezing the sampled solver (M = 1e6) shows artifacts.
Usage: python3 plot_theory_vs_sim.py <out.pdf>"""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42   # TrueType, not Type 3 (ACM PDF requirements)
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],  # Arial/Helvetica
                            'mathtext.fontset': 'custom', 'mathtext.rm': 'Arial', 'mathtext.it': 'Arial:italic',
                            'mathtext.bf': 'Arial:bold', 'mathtext.sf': 'Arial', 'mathtext.fallback': 'stixsans'})
import matplotlib.pyplot as plt
import numpy as np

# Locate the research folders relative to this file (the project folder was renamed; hashed modules keep their old
# absolute paths, so put the real folders first on sys.path and point abl.J_RESULTS at the real file).
_RESEARCH = Path(__file__).resolve().parent.parent
for _d in ('theory_ideas_20261003', 'reaim_reproduction_20261003', 'ablation_20261005', 'algorithm_compare_20261007'):
    sys.path.insert(0, str(_RESEARCH / _d))
import methods as M  # noqa: E402
import abl as _abl  # noqa: E402
_abl.J_RESULTS = _RESEARCH / 'theory_ideas_20261003' / 'J_results.json'

ROOT = Path(__file__).resolve().parent


def simulate(J, lam, B=256, seed=96000):
    S, T0, q = 560, 30.0, 4.0
    T = M.m.sched(T0, S); N = len(J)
    rng = np.random.default_rng(seed)
    s = rng.choice(np.array([-1.0, 1.0], dtype=np.float32), size=(B, N)); h = s @ J
    su = np.empty(S + 1); c = np.empty(S); s_prev = None; c_prev = None
    su[0] = np.mean(np.einsum('bi,bi->b', s, h)) / N / np.sqrt(N)
    for t, Tt in enumerate(T):
        field = h if (lam == 0 or s_prev is None) else h - lam * c_prev * s_prev
        z = s * field + q
        flip = np.clip(z / (4 * Tt) + 0.5, 0, 1) < rng.random(s.shape, dtype=np.float32)
        nlin = ((z > -2 * Tt) & (z < 2 * Tt)).sum(1, keepdims=True)
        c_prev = nlin / (2 * Tt); c[t] = float(c_prev.mean())
        d = np.where(flip, -2 * s, 0).astype(np.float32)
        s_prev = s.copy(); s += d; h += d @ J
        su[t + 1] = np.mean(np.einsum('bi,bi->b', s, h)) / N / np.sqrt(N)
    return T, su, c


def main():
    out = sys.argv[1]
    J, sumw = M.m.load(); N = len(J)
    D = json.loads((ROOT / 'gpu_v2_M1000000_it32_S560_lam0_S560_lam0.7.json').read_text())
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(3.45, 1.85))
    summary = {}
    for lam, key, col in ((0.0, 'S560_lam0', '#1f77b4'), (0.7, 'S560_lam0.7', '#d62728')):
        T, su, c = simulate(J, lam)
        d = D[key]; su_th = np.array(d['su']); c_th = N * np.array(d['plin']) / (2 * T)
        lab = 'Plain' if lam == 0 else '$\\lambda=0.7$'
        ax.plot(su, color=col, lw=1.2, label=f'{lab}, sim.'); ax.plot(np.arange(501), su_th[:501], color='k', lw=0.6, ls='--')
        bx.plot(c, color=col, lw=1.2); bx.plot(np.arange(500), c_th[:500], color='k', lw=0.6, ls='--')
        summary[key] = dict(max_abs_su_diff=float(np.max(np.abs(su - su_th))), final_su_sim=float(su[-1]), final_su_dmft=float(su_th[-1]),
                            max_rel_c_diff=float(np.max(np.abs(c - c_th) / np.maximum(c_th, 1e-9))))
    ax.plot([], [], color='k', lw=0.6, ls='--', label='Theory')
    ax.set_xlabel('Step', fontsize=7.5); ax.set_ylabel('$\\mathbb{E}[s\\,u]$', fontsize=7.5); ax.tick_params(labelsize=6.5)
    ax.legend(fontsize=6.2, frameon=False, loc='lower right', handlelength=1.6); ax.set_title('(a) Energy per Spin', fontsize=7.5, loc='left')
    bx.set_yscale('log'); bx.set_xlabel('Step', fontsize=7.5); bx.set_ylabel('$c=n_{\\mathrm{lin}}/2T$', fontsize=7.5)
    bx.tick_params(labelsize=6.5); bx.set_title('(b) Echo Coefficient', fontsize=7.5, loc='left')
    fig.tight_layout(pad=0.25); fig.savefig(out); fig.savefig(str(Path(out).with_suffix('.png')), dpi=220)
    (ROOT / 'theory_vs_sim.json').write_text(json.dumps(summary, indent=1) + '\n'); print(json.dumps(summary, indent=1))


if __name__ == '__main__':
    main()
