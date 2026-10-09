"""[Copy with TEC removed from ORDER: TEC as published could not be reproduced; see ../results_validation and HANDOVER, 8 Oct 2026]
Full-width (figure*) figure for the paper from the PROTOCOL.md sweep (results/S*.json): (a) mean H(t) at S = 1000, (b) P(cut >= 33,000) vs S.
Usage: python3 plot_alg.py <out.pdf>"""
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
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

ROOT = Path(__file__).resolve().parent
ORDER = ['SA', 'SCA', 'APC-SCA', 'ReAIM ASA', 'aSB', 'bSB', 'dSB', 'TEC-T (ours)', 'Onsager SCA (ours)']  # TEC excluded: not reproduced (fairness audit, 8 Oct 2026)
STY = {'SA': ('#7f7f7f', '-', 'o'), 'SCA': ('#1f77b4', '--', 's'), 'TEC': ('#9467bd', '--', 'v'),
       'APC-SCA': ('#8c564b', '--', '^'), 'ReAIM ASA': ('#e377c2', ':', 'D'), 'aSB': ('#bcbd22', ':', 'P'),
       'bSB': ('#17becf', ':', 'X'), 'dSB': ('#2ca02c', ':', '*'), 'TEC-T (ours)': ('#ff7f0e', '-', 'h'),
       'Onsager SCA (ours)': ('#d62728', '-', 'o')}
H_TARGET, H_BEST = -67040, -67714


def main():
    out = sys.argv[1]
    res = {int(p.stem[1:]): json.loads(p.read_text()) for p in sorted((ROOT / 'results').glob('S*.json'))}
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(7.0, 2.55), gridspec_kw=dict(width_ratios=[1.15, 1]))
    r = res.get(1000, {})
    DISPLAY = {'TEC-T (ours)': 'Onsager-$\\kappa T$ (this work)', 'Onsager SCA (ours)': 'Onsager-online (this work)'}
    hmap = {}
    for n in ORDER:
        c, ls, mk = STY[n]; lw = 1.1 if 'ours' not in n else 1.7
        hmap[n] = Line2D([], [], color=c, ls=ls, marker=mk, ms=3.5, lw=lw, label=DISPLAY.get(n, n))
        if n in r and 'final' in r[n]:
            H = np.array(r[n]['final']['H_mean'])
            ax.plot(np.arange(len(H)), H, color=c, ls=ls, lw=lw)
    handles = [hmap[n] for n in ORDER[:-2] + ['Onsager SCA (ours)', 'TEC-T (ours)']]
    handles = handles[:-2] + [Line2D([], [], color='none', label=' ')] + handles[-2:]  # blank slot keeps both 'this work' entries in the last column
    ax.axhline(H_TARGET, color='k', lw=0.7, ls='-.')
    ax.set_ylim(-68000, -36000); ax.set_xlim(0, 1000)   # headroom above the curves for the inset (no overlap)
    ax.set_xlabel('Monte Carlo Step', fontsize=8); ax.set_ylabel('Mean Ising Energy $H$', fontsize=8)
    ax.tick_params(labelsize=7); ax.set_yticks(range(-65000, -35000, 5000)); ax.set_title('(a) Mean Energy, $S=1000$', fontsize=8, loc='left')
    ins = inset_axes(ax, width='50%', height='44%', loc='upper right', borderpad=0.9)
    for n in ORDER:
        if n not in r or 'final' not in r[n]:
            continue
        c, ls, _ = STY[n]; H = np.array(r[n]['final']['H_mean'])
        ins.plot(np.arange(700, len(H)), H[700:], color=c, ls=ls, lw=0.9 if 'ours' not in n else 1.5)
    ins.axhline(H_TARGET, color='k', lw=0.6, ls='-.'); ins.set_ylim(-67550, -66300); ins.tick_params(labelsize=6.5); ins.set_yticks([-67500, -67000, -66500])
    for n in ORDER:
        xs, ys, lo, hi = [], [], [], []
        for S in sorted(res):
            if n in res[S] and 'final' in res[S][n]:
                f = res[S][n]['final']; xs.append(S); ys.append(f['p33000'])
                lo.append(max(0.0, f['p33000'] - f['p33000_wilson95'][0])); hi.append(max(0.0, f['p33000_wilson95'][1] - f['p33000']))
        if xs:
            c, ls, mk = STY[n]
            bx.errorbar(xs, ys, yerr=[lo, hi], color=c, ls=ls, marker=mk, ms=3.5, lw=1.1 if 'ours' not in n else 1.7, capsize=1.8)
    bx.set_xscale('log'); bx.set_xticks(sorted(res)); bx.set_xticklabels([str(s) for s in sorted(res)])
    bx.set_ylim(-0.02, 1.02); bx.set_xlabel('Step Budget $S$', fontsize=8); bx.set_ylabel('P(Cut $\\geq$ 33,000)', fontsize=8)
    bx.tick_params(labelsize=7); bx.set_title('(b) Success Probability vs. Step Budget', fontsize=8, loc='left')
    bx.minorticks_off()
    fig.legend(handles=handles, loc='upper center', ncol=5, fontsize=7, frameon=False, bbox_to_anchor=(0.5, 1.0),
               columnspacing=1.6, handlelength=2.8)
    fig.tight_layout(pad=0.4, rect=(0, 0, 1, 0.85))
    fig.savefig(out)
    fig.savefig(str(Path(out).with_suffix('.png')), dpi=200)
    print('saved', out, 'budgets', sorted(res))


if __name__ == '__main__':
    main()
