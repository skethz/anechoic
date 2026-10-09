"""Figure for the G-set speedups measured on the V80 (replaces the old Table 8a; 9 Oct 2026 KST).

Data: the board's Amendment 3 summary (a3_summary.json), extended grids ('extended_amendment1'), read through
board_variants.load(); no new data. Per instance: speedup = TTS99(reference rule) / TTS99(corrected rule), primary
12-engine estimator. Before plotting, the script checks that the geometric means and win counts equal the published
Table 8a values (board_variants.json), so the figure and the text cannot disagree.
Usage: python3 plot_gset_speedup.py <out.pdf>"""
import json
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42   # TrueType, not Type 3 (ACM PDF requirements)
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
                            'mathtext.fontset': 'custom', 'mathtext.rm': 'Arial', 'mathtext.it': 'Arial:italic',
                            'mathtext.bf': 'Arial:bold', 'mathtext.sf': 'Arial', 'mathtext.fallback': 'stixsans'})
import matplotlib.pyplot as plt
import matplotlib.ticker

sys.path.insert(0, str(Path(__file__).resolve().parent))
import board_variants as BV  # noqa: E402

VARIANT = 'extended_amendment1'
EXPECT = {  # Table 8a as published (board_variants.json, extended grids): geomean, wins of 51
    ('Onsager-kT', 'SCA'): (2.17, 49), ('Onsager-kT', 'TEC'): (1.68, 42),
    ('Onsager-online', 'SCA'): (1.26, 35), ('Onsager-online', 'TEC'): (0.98, 28)}
LABEL = {'Random, +1': 'Random\n$+1$', 'Random, +-1': 'Random\n$\\pm1$', 'Toroidal, +-1': 'Toroidal\n$\\pm1$',
         'Planar-like, +1': 'Planar-like\n$+1$', 'Planar-like, +-1': 'Planar-like\n$\\pm1$'}
REF = {'SCA': ('plain SCA', '#1f77b4', 'o'), 'TEC': ('TEC', '#9467bd', '^')}
GAP = 1.6   # horizontal gap between classes, in instance widths


def gmean(xs):
    return math.exp(sum(map(math.log, xs)) / len(xs))


def main():
    out = sys.argv[1]
    tv = BV.load()[VARIANT]
    assert len(tv) == 51
    # x positions: classes in table order, instances by G number within each class
    xs, x, spans = {}, 0.0, []
    for cls in BV.CLASSES:
        inst = sorted((g for g in tv if BV.wclass(int(g[1:])) == cls), key=lambda g: int(g[1:]))
        x0 = x
        for g in inst:
            xs[g] = x; x += 1.0
        spans.append((cls, inst, x0, x - 1.0)); x += GAP
    summary = {}
    for (a, b), (gm_pub, wins_pub) in EXPECT.items():
        r = [tv[g][b][0] / tv[g][a][0] for g in tv]
        gm, wins = gmean(r), sum(v > 1 for v in r)
        assert round(gm, 2) == gm_pub and wins == wins_pub, (a, b, gm, wins)
        summary[f'{a} vs {b}'] = dict(geomean=gm, wins=wins, min=min(r), max=max(r))
    fig, axes = plt.subplots(2, 1, figsize=(3.4, 2.2), sharex=True)
    for ax, rule, title, loc in zip(axes, ('Onsager-kT', 'Onsager-online'), ('(a) Onsager-$\\kappa T$', '(b) Onsager-online'),
                                    ('lower right', 'upper right')):   # legends where each panel has no markers
        ax.axhline(1.0, color='0.35', lw=0.6, ls=(0, (3, 2)), zorder=1)
        for ref in ('SCA', 'TEC'):
            name, col, mk = REF[ref]
            gm, wins = summary[f'{rule} vs {ref}']['geomean'], summary[f'{rule} vs {ref}']['wins']
            for cls, inst, x0, x1 in spans:
                ax.plot([x0 - 0.4, x1 + 0.4], [gmean([tv[g][ref][0] / tv[g][rule][0] for g in inst])] * 2,
                        color=col, lw=1.1, alpha=0.85, zorder=2, solid_capstyle='butt')
            ax.scatter([xs[g] for g in tv], [tv[g][ref][0] / tv[g][rule][0] for g in tv], s=7, marker=mk,
                       facecolors=col if ref == 'SCA' else 'white', edgecolors=col, linewidths=0.6, zorder=3,
                       label=f'vs {name}: {gm:.2f}$\\times$, faster on {wins}/51')
        ax.set_yscale('log', base=2); ax.set_ylim(0.14, 20)
        ax.yaxis.set_major_locator(matplotlib.ticker.FixedLocator([0.25, 0.5, 1, 2, 4, 8, 16]))
        ax.yaxis.set_major_formatter(matplotlib.ticker.FixedFormatter(['0.25', '0.5', '1', '2', '4', '8', '16']))
        ax.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
        ax.tick_params(labelsize=6.5, length=2); ax.set_ylabel('Speedup', fontsize=7.5, labelpad=1)
        ax.text(0.01, 0.97, title, transform=ax.transAxes, ha='left', va='top', fontsize=7, fontweight='bold')
        ax.legend(fontsize=6.0, frameon=False, loc=loc, handletextpad=0.2, borderaxespad=0.2,
                  labelspacing=0.15, markerscale=1.3)
        for cls, inst, x0, x1 in spans[1:]:
            ax.axvline(x0 - GAP / 2 - 0.5, color='0.8', lw=0.5, zorder=0)
    axes[1].set_xticks([(x0 + x1) / 2 for _, _, x0, x1 in spans])
    axes[1].set_xticklabels([f'{LABEL[c]} ({len(i)})' for c, i, _, _ in spans], fontsize=6.2, linespacing=0.95)
    axes[1].tick_params(axis='x', length=0, pad=2)
    axes[0].set_xlim(-1.0, x - GAP + 0.0)
    fig.tight_layout(pad=0.25, h_pad=0.3)
    fig.savefig(out); fig.savefig(str(Path(out).with_suffix('.png')), dpi=220)
    Path(__file__).with_name('plot_gset_speedup.json').write_text(json.dumps(summary, indent=1) + '\n')
    print(json.dumps(summary, indent=1))


if __name__ == '__main__':
    main()
