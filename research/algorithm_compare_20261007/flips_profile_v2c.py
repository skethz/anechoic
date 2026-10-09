"""v2c copy of flips_profile_v2.py (taller figure only; same seeds and data). Where the engine's cycles go: mean flips per step over the anneal for the hardware configurations, on a time axis from the
measured v6.4 cycle model (cycles/step = 18.09 + 0.1424*flips, +886 per trial, 250 MHz); the first version used the v6.2 model (16.91, 0.2714, 1148, 275 MHz). Float model (abl.py); 128 runs each.
Usage: python3 flips_profile_v2.py <out.pdf>"""
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
import matplotlib.ticker
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
# v2 (9 Oct 2026 KST): adds STATICA's two published K2000 schedules (B5, B7; q 4, geometric T to 5), so that the figure
# explains the headline comparison of Section 6.2. Seeds are tied to each configuration (SEED), so the five original
# curves are bit-identical to v1. Board p and round times for B5/B7: Table 2 (A5 run, hwv6_board_v64_e12_250mhz_a5).
CFGS = [  # key, label, cfg (identical to scripts/run_hw_v6.sh), on-board p, on-board mean round time [ms]
    ('S7', 'STATICA (B7)', dict(family='plain', q=4.0, T0=40.0, S=1560), 0.824, 0.302),
    ('S5', 'STATICA (B5)', dict(family='plain', q=4.0, T0=30.0, S=560), 0.134, 0.094),  # 0.094473 ms (a5_summary.json); was 0.095
    ('P3', 'Plain SCA (B3)', dict(family='plain', q=8.0, T0=30.0, S=1560), 0.439, 0.1993),
    ('B1', 'Plain SCA (B2)', dict(family='plain', q=8.0, T0=40.0, S=960), 0.461, 0.1568),
    ('T2', 'TEC-style (B1)', dict(family='tec', q=8.0, jv=-4.0, T0=30.0, S=960), 0.392, 0.1366),
    ('O4', 'Onsager-online (O2)', dict(family='onsager', q=6.0, lam=0.9, ramp=False, T0=12.0, S=360), 0.326, 0.0590),
    ('X5', 'Onsager-$\\kappa T$ (O1)', dict(family='tecT', q=8.0, kappa=1.75, ramp=True, T0=15.0, S=280), 0.338, 0.0503),
]
SEED = {'T2': 95000, 'B1': 95001, 'P3': 95002, 'X5': 95003, 'O4': 95004, 'S5': 95005, 'S7': 95006}  # v1 used 95000 + index
COL = {'S7': '#08306b', 'S5': '#4292c6', 'P3': '#1f77b4', 'B1': '#9ecae1', 'T2': '#9467bd', 'O4': '#d62728', 'X5': '#ff7f0e'}


def run(J, cfg, B, seed):
    """sca_traj with flips per step recorded (same arithmetic)."""
    rng = np.random.default_rng(seed)
    s0 = rng.choice(np.array([-1.0, 1.0], dtype=np.float32), size=(B, len(J)))
    fam = cfg['family']; S = cfg['S']; T = M.m.sched(cfg['T0'], S)
    s = s0.copy(); h = s @ J; s_prev = None; c_prev = None; T_prev = None
    q = np.full(s.shape, float(cfg.get('q', 0.0)), np.float32)
    flips = np.zeros((S, B))
    for t, Tt in enumerate(T):
        field = h
        if s_prev is not None:
            if fam == 'onsager':
                field = h - cfg['lam'] * M.abl.ramp_factor(t, S, cfg['ramp']) * c_prev * s_prev
            elif fam == 'tecT':
                field = h - cfg['kappa'] * T_prev * M.abl.ramp_factor(t, S, cfg['ramp']) * s_prev
            elif fam == 'tec':
                field = h + cfg['jv'] * s_prev
        z = s * field + q
        flip = np.clip(z / (4 * Tt) + 0.5, 0, 1) < rng.random(s.shape, dtype=np.float32)
        if fam == 'onsager':
            c_prev = ((z > -2 * Tt) & (z < 2 * Tt)).sum(1, keepdims=True) / (2 * Tt)
        d = np.where(flip, -2 * s, 0).astype(np.float32)
        s_prev = s.copy(); T_prev = Tt
        s += d; h += d @ J
        flips[t] = flip.sum(1)
    return flips, s, h


def main():
    out = sys.argv[1]
    J, sumw = M.m.load()
    fig, ax = plt.subplots(figsize=(3.4, 1.86))   # v2c (9 Oct, 19:30 KST): 0.18 in taller to fill a 1.4-line column gap on page 8
    summary = {}
    for i, (key, lab, cfg, p_hw, t_hw) in enumerate(CFGS):
        flips, s, h = run(J, cfg, 128, SEED[key])
        cuts = (sumw + 0.5 * np.einsum('bi,bi->b', s, h)) / 2
        mf = flips.mean(1); cyc = 18.09 + 0.1424 * mf
        t_us = np.cumsum(cyc) / 250.0
        tot_ms = (cyc.sum() + 886) / 250e3
        summary[key] = dict(mean_flips=float(flips.sum(0).mean()), model_trial_ms=tot_ms, p_model=float((cuts >= 33000).mean()),
                            p_board=p_hw, t_round_board_ms=t_hw)
        ax.plot(t_us, mf, color=COL[key], lw=1.4 if key in ('O4', 'X5') else 1.0,
                label=f"{lab}: {p_hw:.2f}, {t_hw * 1e3:.0f} µs")
        ax.plot([t_us[-1]], [mf[-1]], marker='o', ms=3, color=COL[key])
    ax.set_yscale('log'); ax.set_ylim(0.5, 2.5e5); ax.set_yticks([1, 10, 100, 1000]); ax.yaxis.set_minor_locator(matplotlib.ticker.FixedLocator([m * 10.0 ** e for e in range(-1, 4) for m in range(2, 10)]))   # v2: headroom for a two-column legend; ticks only where the data are
    ax.set_xlabel('Time on One Engine at 250 MHz [µs]', fontsize=8); ax.set_ylabel('Mean Flips per Step', fontsize=8)
    ax.tick_params(labelsize=7); ax.legend(fontsize=6.0, frameon=False, loc='upper right', handlelength=1.4, labelspacing=0.22, ncol=2, columnspacing=0.8)
    fig.tight_layout(pad=0.3); fig.savefig(out); fig.savefig(str(Path(out).with_suffix('.png')), dpi=200)
    (ROOT / 'flips_profile_v2c.json').write_text(json.dumps(summary, indent=1) + '\n')
    print(json.dumps(summary, indent=1))


if __name__ == '__main__':
    main()
