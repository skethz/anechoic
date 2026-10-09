"""Compare reproduced SCA statistics with STATICA's published values; fit the DDSS timing constant."""
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from statica_repro import GRID_S, GRID_TINIT, PAPER, clopper_pearson, tts

ROOT = Path(__file__).resolve().parent
CLOCK_PER_MS = 300e3  # 300 MHz
N = 2000

# Approximate digitization of Fig. 22 (mean energy, annealing time) and Fig. 24 (TTS), 400-dpi render.
# Labelled values (-66540, -67186, 0.13, 0.48) are exact; the rest are read from pixel positions (+-~10 energy, +-~0.01 ms).
FIG22_H = {50.0: [-66510, -66824, -67010, -67048, -67097, -67182, -67208],
           40.0: [-66494, -66769, -66907, -67060, -67113, -67155, -67186],
           30.0: [-66178, -66540, -66664, -66716, -66836, -66928, -66954]}
FIG22_T = {50.0: [0.163, 0.243, 0.320, 0.396, 0.473, 0.553, 0.629],
           40.0: [0.132, 0.189, 0.248, 0.306, 0.364, 0.420, 0.480],
           30.0: [0.094, 0.130, 0.166, 0.201, 0.237, 0.271, 0.309]}
FIG24_TTS = {50.0: [None, 3.53, 1.99, 2.48, 2.07, 1.89, 1.74],
             40.0: [None, 5.00, 2.85, 2.08, 1.65, 1.39, 1.50],
             30.0: [None, 8.23, 7.24, 5.33, 3.44, 2.52, 2.78]}


def implied_p(t, tts_ms):
    return 1 - math.exp(t * math.log(0.01) / tts_ms)


def main():
    raw = json.loads((ROOT / "raw_results.json").read_text())
    R = raw["results"]
    out = {"self_test": raw["self_test"], "instance": raw["instance"], "target_H": raw["target_H"]}

    # Timing constant c0: cycles = flips-only cycles + S * c0, least squares on the two Table II times.
    pts = [(R[k]["S"], R[k]["mean_cycles_flips_only"], PAPER[k]["t_ms"] * CLOCK_PER_MS) for k in ("short", "long")]
    c0 = sum(S * (c - f) for S, f, c in pts) / sum(S * S for S, _, _ in pts)
    out["timing"] = {"c0_cycles_per_step": c0,
                     "c0_from_each_point": {k: (PAPER[k]["t_ms"] * CLOCK_PER_MS - R[k]["mean_cycles_flips_only"]) / R[k]["S"]
                                            for k in ("short", "long")}}

    def t_model(r):
        return (r["mean_cycles_flips_only"] + r["S"] * c0) / CLOCK_PER_MS

    table2 = {}
    for k in ("short", "long"):
        r, p = R[k], PAPER[k]
        lo, hi = r["p_ci95"]
        paper_ci = clopper_pearson(round(p["p"] * 100), 100)
        table2[k] = dict(S=r["S"], T_init=r["T_init"], trials=r["trials"], successes=r["successes"], p=r["p"], p_ci95=r["p_ci95"],
                         paper_p=p["p"], paper_p_ci95_100runs=paper_ci, p_inside_paper_ci=paper_ci[0] <= r["p"] <= paper_ci[1],
                         mean_cut=r["mean_cut"], paper_mean_cut=p["mean_cut"], sd_cut=r["sd_H"] / 2,
                         z_mean_cut_vs_paper=(r["mean_cut"] - p["mean_cut"]) / (r["sd_H"] / 2 / math.sqrt(100)),
                         t_model_ms=t_model(r), paper_t_ms=p["t_ms"],
                         tts_paper_t_our_p=tts(p["t_ms"], r["p"]), tts_paper_t_our_p_range=(tts(p["t_ms"], hi), tts(p["t_ms"], lo)),
                         tts_model_t_our_p=tts(t_model(r), r["p"]), paper_tts=p["tts_ms"],
                         edge_rescore=r["edge_rescore"])
    out["table2"] = table2

    grid = []
    for T in GRID_TINIT:
        for i, S in enumerate(GRID_S):
            r = R[f"grid_S{S}_T{T:g}"]
            tm = t_model(r)
            ptts = FIG24_TTS[T][i]
            grid.append(dict(S=S, T_init=T, successes=r["successes"], trials=r["trials"], p=r["p"], p_ci95=r["p_ci95"],
                             mean_H=r["mean_H"], fig22_H=FIG22_H[T][i], dH=r["mean_H"] - FIG22_H[T][i],
                             se_H_100runs=r["sd_H"] / math.sqrt(100),
                             t_model_ms=tm, fig22_t_ms=FIG22_T[T][i],
                             tts_model_t=tts(tm, r["p"]), tts_fig22_t=tts(FIG22_T[T][i], r["p"]), fig24_tts=ptts,
                             fig24_implied_p=None if ptts is None else implied_p(FIG22_T[T][i], ptts)))
    out["grid"] = grid
    (ROOT / "results.json").write_text(json.dumps(out, indent=1) + "\n")

    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    colors = {50.0: "tab:blue", 40.0: "tab:green", 30.0: "goldenrod"}
    for T in GRID_TINIT:
        g = [x for x in grid if x["T_init"] == T]
        S = [x["S"] for x in g]
        c = colors[T]
        ax[0].plot(S, [x["mean_H"] for x in g], "o-", color=c, label=f"ours, T_init={T:g}")
        ax[0].plot(S, [x["fig22_H"] for x in g], "x--", color=c, alpha=.6, label=f"paper Fig. 22, T_init={T:g}")
        ax[1].plot(S, [x["t_model_ms"] for x in g], "o-", color=c)
        ax[1].plot(S, [x["fig22_t_ms"] for x in g], "x--", color=c, alpha=.6)
        ax[2].plot(S[1:], [x["tts_model_t"] for x in g][1:], "o-", color=c)
        ax[2].plot(S[1:], [x["fig24_tts"] for x in g][1:], "x--", color=c, alpha=.6)
    ax[0].axhline(raw["target_H"], color="k", lw=.8, ls=":")
    ax[0].set(xlabel="MC steps S", ylabel="mean final Ising energy", title="Energy (circles: ours, x: paper)")
    ax[1].set(xlabel="MC steps S", ylabel="annealing time [ms]", title=f"Time (ours: DDSS model, c0={c0:.1f} cycles/step)")
    ax[2].set(xlabel="MC steps S", ylabel="tts(0.99) [ms]", title="TTS, Eq. (9)", ylim=(0, 10))
    ax[0].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(ROOT / "reproduction.png", dpi=130)
    print(json.dumps({"timing": out["timing"], "table2": table2}, indent=1))
    for x in grid:
        print(f"S={x['S']:5d} T0={x['T_init']:4.0f}  p={x['p']:.3f} [{x['p_ci95'][0]:.3f},{x['p_ci95'][1]:.3f}]"
              f"  H={x['mean_H']:9.1f} vs {x['fig22_H']:7d} (d={x['dH']:+6.1f}, se100={x['se_H_100runs']:4.1f})"
              f"  t={x['t_model_ms']:.3f} vs {x['fig22_t_ms']:.3f}  TTS={x['tts_model_t']:.2f} vs {x['fig24_tts']}"
              f"  implied paper p={x['fig24_implied_p'] if x['fig24_implied_p'] is None else round(x['fig24_implied_p'], 3)}")


if __name__ == "__main__":
    main()
