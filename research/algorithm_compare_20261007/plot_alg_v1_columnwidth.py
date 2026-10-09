"""Figure for the paper from the PROTOCOL.md sweep (results/S*.json): (a) mean H(t) at S = 1000, (b) P(cut >= 33,000) vs S.
Usage: python3 plot_alg.py <out.pdf>"""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42   # TrueType, not Type 3 (ACM PDF requirements)
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

ROOT = Path(__file__).resolve().parent
ORDER = ['SA', 'SCA', 'TEC', 'APC-SCA', 'ReAIM ASA', 'aSB', 'bSB', 'dSB', 'TEC-T (ours)', 'Onsager SCA (ours)']
STY = {'SA': ('#7f7f7f', '-', 'o'), 'SCA': ('#1f77b4', '--', 's'), 'TEC': ('#9467bd', '--', 'v'),
       'APC-SCA': ('#8c564b', '--', '^'), 'ReAIM ASA': ('#e377c2', ':', 'D'), 'aSB': ('#bcbd22', ':', 'P'),
       'bSB': ('#17becf', ':', 'X'), 'dSB': ('#2ca02c', ':', '*'), 'TEC-T (ours)': ('#ff7f0e', '-', 'h'),
       'Onsager SCA (ours)': ('#d62728', '-', 'o')}
H_TARGET, H_BEST = -67040, -67714


def main():
    out = sys.argv[1]
    res = {int(p.stem[1:]): json.loads(p.read_text()) for p in sorted((ROOT / 'results').glob('S*.json'))}
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(7.0, 2.6), gridspec_kw=dict(width_ratios=[1.25, 1]))
    r = res.get(1000, {})
    for n in ORDER:
        if n not in r or 'final' not in r[n]:
            continue
        c, ls, _ = STY[n]; H = np.array(r[n]['final']['H_mean'])
        ax.plot(np.arange(len(H)), H, color=c, ls=ls, lw=1.1 if 'ours' not in n else 1.6, label=n.replace(' (ours)', '*'))
    ax.axhline(H_TARGET, color='k', lw=0.6, ls='-.')
    ax.set_ylim(-68000, -50000); ax.set_xlim(0, 1000)
    ax.set_xlabel('Monte Carlo step', fontsize=7); ax.set_ylabel('mean Ising energy $H$', fontsize=7)
    ax.tick_params(labelsize=6); ax.set_title('(a) $S=1000$', fontsize=7, loc='left')
    ax.legend(fontsize=5.2, ncol=2, loc='upper right', frameon=False)
    ins = inset_axes(ax, width='42%', height='38%', loc='center right', borderpad=0.6)
    for n in ORDER:
        if n not in r or 'final' not in r[n]:
            continue
        c, ls, _ = STY[n]; H = np.array(r[n]['final']['H_mean'])
        ins.plot(np.arange(700, len(H)), H[700:], color=c, ls=ls, lw=0.9 if 'ours' not in n else 1.4)
    ins.axhline(H_TARGET, color='k', lw=0.5, ls='-.'); ins.set_ylim(-67550, -66300); ins.tick_params(labelsize=5)
    for n in ORDER:
        xs, ys, lo, hi = [], [], [], []
        for S in sorted(res):
            if n in res[S] and 'final' in res[S][n]:
                f = res[S][n]['final']; xs.append(S); ys.append(f['p33000'])
                lo.append(max(0.0, f['p33000'] - f['p33000_wilson95'][0])); hi.append(max(0.0, f['p33000_wilson95'][1] - f['p33000']))
        if xs:
            c, ls, mk = STY[n]
            bx.errorbar(xs, ys, yerr=[lo, hi], color=c, ls=ls, marker=mk, ms=3, lw=1.0 if 'ours' not in n else 1.5, capsize=1.5)
    bx.set_xscale('log'); bx.set_xticks(sorted(res)); bx.set_xticklabels([str(s) for s in sorted(res)])
    bx.set_ylim(-0.02, 1.02); bx.set_xlabel('step budget $S$', fontsize=7); bx.set_ylabel('P(cut $\\geq$ 33,000)', fontsize=7)
    bx.tick_params(labelsize=6); bx.set_title('(b) success vs budget (* = ours)', fontsize=7, loc='left')
    bx.minorticks_off()
    fig.tight_layout(pad=0.4)
    fig.savefig(out)
    fig.savefig(str(Path(out).with_suffix('.png')), dpi=200)
    print('saved', out, 'budgets', sorted(res))


if __name__ == '__main__':
    main()
