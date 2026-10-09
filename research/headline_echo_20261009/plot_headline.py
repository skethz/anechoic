"""Figure 1 of the paper (v3, 9 Oct 2026): (a) the echo of synchronous SCA, (b) Anechoic subtracts it, (c) its effect.
v3: the echo path starts at the pre-flip state s_i(t-1) and returns in h_i(t) as c(t-1) s_i(t-1) (Eq. 5);
echo shares exclude the first step, where no flip can be undone (same definition as Section 4.2: B5 67.1%, O1 37.3%).
(c) uses the engine's float model on K2000 (1,024 runs per schedule, echo_stats.json); V80 values are from Table 2.
Usage: python3 plot_headline.py <out.pdf>"""
import json, sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'pdf.fonttype': 42, 'ps.fonttype': 42, 'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'], 'mathtext.fontset': 'custom', 'mathtext.rm': 'Arial',
    'mathtext.it': 'Arial:italic', 'mathtext.bf': 'Arial:bold', 'mathtext.fallback': 'stixsans'})
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch
H = Path(__file__).resolve().parent
D = json.loads((H / 'echo_stats.json').read_text())
UP, DN, ECHO, OK, GREY, TXT, SUB = '#2b6cb0', '#dd6b20', '#c53030', '#2f855a', '#a0aec0', '#1a202c', '#4a5568'
R = 0.23
XS, Y0, NB = [0.42, 1.55, 2.68], 0.0, (1.55, 0.86)

def echo_share(k):
    f, b = D[k]['flips_step'], D[k]['back_step']
    return sum(b[1:]) / sum(f[1:])

def spin(ax, x, y, s):
    ax.add_patch(Circle((x, y), R, fc=UP if s > 0 else DN, ec='none', zorder=3))
    ax.text(x, y - 0.012, '+' if s > 0 else '−', ha='center', va='center', color='white', fontsize=11.5,
            fontweight='bold', zorder=4)

def arrow(ax, p, q, color, rad=0.0, lw=1.2, ls='-', alpha=1.0):
    ax.add_patch(FancyArrowPatch(p, q, connectionstyle=f'arc3,rad={rad}', arrowstyle='-|>', mutation_scale=9,
                                 lw=lw, color=color, ls=ls, alpha=alpha, zorder=2, shrinkA=2, shrinkB=2))

def schematic(ax, cancelled):
    ax.set_xlim(0.0, 3.62); ax.set_ylim(-0.50, 1.20); ax.set_aspect('equal'); ax.axis('off'); ax.set_anchor('N')
    for x, s in zip(XS, [+1, -1, -1 if cancelled else +1]): spin(ax, x, Y0, s)
    for x, lab in zip(XS, ['$s_i(t{-}1)$', '$s_i(t)$', '$s_i(t{+}1)$']):
        ax.text(x, Y0 - R - 0.07, lab, ha='center', va='top', fontsize=7.2, color=TXT)
    arrow(ax, (XS[0] + R, Y0), (XS[1] - R, Y0), GREY); arrow(ax, (XS[1] + R, Y0), (XS[2] - R, Y0), GREY)
    ax.text((XS[0] + XS[1]) / 2, Y0 + 0.07, 'flip', ha='center', va='bottom', fontsize=6.8, color=SUB)
    for dx in (-0.15, 0.0, 0.15): ax.add_patch(Circle((NB[0] + dx, NB[1]), 0.065, fc='#cbd5e0', ec='none', zorder=3))
    ax.text(NB[0], NB[1] + 0.11, 'neighbors', ha='center', va='bottom', fontsize=6.6, color=SUB)
    # the pre-flip state s_i(t-1) enters the neighbors' fields; their response returns into the decision at step t
    p0, p1 = (XS[0] + 0.06, Y0 + R), (NB[0] - 0.24, NB[1] - 0.02)
    arrow(ax, p0, p1, GREY, rad=-0.28, lw=1.0)
    q0, q1, rad = (NB[0] + 0.24, NB[1] - 0.02), (XS[2] - 0.06, Y0 + R), -0.28
    if not cancelled:
        arrow(ax, q0, q1, ECHO, rad=rad, lw=1.6)
        ax.text(2.47, 0.93, 'echo\n$c(t{-}1)\\,s_i(t{-}1)$', ha='left', va='center', fontsize=6.8, color=ECHO,
                linespacing=1.15)
        ax.text(XS[2] + R + 0.07, Y0, 'flip\nundone', ha='left', va='center', fontsize=7, color=ECHO,
                fontweight='bold', linespacing=1.05)
    else:
        arrow(ax, q0, q1, ECHO, rad=rad, lw=1.3, ls=(0, (2.5, 1.5)), alpha=0.45)
        dx, dy = q1[0] - q0[0], q1[1] - q0[1]
        mx, my = (q0[0] + q1[0]) / 2 + rad * dy / 2, (q0[1] + q1[1]) / 2 - rad * dx / 2
        ax.text(mx, my, '×', ha='center', va='center', fontsize=13, color=OK, fontweight='bold', zorder=5)
        ax.text(2.47, 0.93, 'echo\nsubtracted', ha='left', va='center', fontsize=6.8, color=OK, linespacing=1.15)
        ax.text(XS[2] + R + 0.07, Y0, 'flip\nkept', ha='left', va='center', fontsize=7, color=OK,
                fontweight='bold', linespacing=1.05)

def frac(k, w=0.04):
    f = np.array(D[k]['flips_step']); bk = np.array(D[k]['back_step']); S = len(f); n = max(3, int(w * S))
    ker = np.ones(n); num = np.convolve(bk, ker, 'same'); den = np.convolve(f, ker, 'same')
    x = np.arange(S) + 1; m = (den > 0.5) & (x > n // 2) & (x <= S - n // 2)
    return x[m], num[m] / den[m]

def main():
    out = Path(sys.argv[1])
    e_b5, e_o1 = round(100 * echo_share('B5')), round(100 * echo_share('O1'))
    fig = plt.figure(figsize=(7.0, 2.0))
    outer = fig.add_gridspec(1, 2, width_ratios=[2.0, 1.02], wspace=0.17, left=0.0, right=0.995, top=0.885, bottom=0.165)
    inner = outer[0].subgridspec(1, 2, wspace=0.0)
    a, b, c = fig.add_subplot(inner[0]), fig.add_subplot(inner[1]), fig.add_subplot(outer[1])
    schematic(a, cancelled=False); schematic(b, cancelled=True)
    A_TXT = ('1. All spins update at once, so the neighbors\n'
             "    respond to spin $i$'s old state $s_i(t{-}1)$ as it flips.\n"
             "2. Their response returns in the spin's local field\n"
             '    $h_i(t)$ as the echo $c(t{-}1)\\,s_i(t{-}1)$ of strength $c$\n'
             f'    and pushes the spin back: {e_b5}% of all flips undo\n'
             '    a flip one step earlier (c).')
    B_TXT = ('1. Each decision uses the corrected field\n'
             '    $h_i(t)-\\lambda\\,c(t{-}1)\\,s_i(t{-}1)$, with a gain $\\lambda$.\n'
             '2. $c$ needs one population count per step, or none\n'
             '    if set from the temperature (Onsager-$\\kappa T$).\n'
             f'3. The flip is kept: only {e_o1}% of flips are undone (c).')
    for ax, txt in ((a, A_TXT), (b, B_TXT)):
        ax.text(0.10, -0.56, txt, ha='left', va='top', fontsize=6.6, color=TXT, linespacing=1.3, clip_on=False)
    for k, lab, col, xy, ha in (('B5', "Plain SCA, STATICA's schedule\nTTS on the V80: 0.231 ms", ECHO, (598, 63.5), 'right'),
                                ('O1', 'Anechoic, Onsager-$\\kappa T$\nTTS on the V80: 0.050 ms', OK, (22, 43.5), 'left')):
        x, y = frac(k); c.plot(x, 100 * y, color=col, lw=1.7)
        c.plot([x[-1]], [100 * y[-1]], marker='o', ms=3.5, color=col)
        c.text(*xy, lab, color=col, fontsize=6.6, ha=ha, va='bottom', multialignment=ha, linespacing=1.15)
    c.set_xlim(0, 600); c.set_ylim(0, 80); c.set_xticks([0, 140, 280, 420, 560]); c.set_yticks([0, 20, 40, 60, 80])
    c.set_xlabel('Annealing Step', fontsize=7.5, labelpad=1)
    c.set_ylabel('Flips Undoing the\nPrevious Step [%]', fontsize=7.2, labelpad=2)
    c.tick_params(labelsize=6.8, length=2); c.spines[['top', 'right']].set_visible(False)
    for ax in (a, b): ax.apply_aspect()
    for ax, ttl, dx in ((a, '(a) Synchronous SCA: the echo', 0.0), (b, '(b) Anechoic: subtract the echo', 0.0),
                        (c, '(c) Effect on K2000', -0.062)):
        fig.text(ax.get_position().x0 + dx + 0.005, 0.985, ttl, fontsize=8.6, fontweight='bold', va='top', color=TXT)
    fig.savefig(out); fig.savefig(out.with_suffix('.png'), dpi=250)
    print('saved', out, 'echo', e_b5, e_o1)

if __name__ == '__main__':
    main()
