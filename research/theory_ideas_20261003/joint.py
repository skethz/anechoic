"""Frozen joint sweep + TEC comparison (PROTOCOL_J.md)."""
import itertools
import json
import time
from pathlib import Path

import numpy as np

import sweeps_abc as m

ROOT = Path(__file__).resolve().parent
TS, SS = (12.0, 15.0, 20.0, 30.0), (360, 560, 960, 1560)


def run_sca(J, s0, S, t0, q, mode, par, ramp, rng):
    """mode: 'plain' | 'onsager' (par=lambda) | 'tec' (par=J_v)."""
    s = s0.astype(np.float32).copy(); h = s @ J
    flips = np.zeros(len(s)); s_prev = None; c_prev = None
    for t, T in enumerate(m.sched(t0, S)):
        if mode == "onsager" and s_prev is not None:
            lam = par if (not ramp or t < 0.7 * S) else par * (S - 1 - t) / max(1.0, S - 1 - 0.7 * S)
            field = h - lam * c_prev[:, None] * s_prev
        elif mode == "tec" and s_prev is not None:
            field = h + par * s_prev
        else:
            field = h
        z = s * field + q
        flip = np.clip(z / (4 * T) + 0.5, 0, 1) < rng.random(s.shape, dtype=np.float32)
        if mode == "onsager":
            c_prev = ((z > -2 * T) & (z < 2 * T)).sum(1) / (2 * T)
        d = np.where(flip, -2 * s, 0).astype(np.float32)
        s_prev = s.copy(); s += d; h += d @ J; flips += flip.sum(1)
    return s, flips


def grid():
    g = [("plain", q, 0.0, False, t0, S) for q, t0, S in itertools.product((4.0, 6.0, 8.0), TS, SS)]
    g += [("onsager", q, l, r, t0, S) for q, l, r, t0, S in itertools.product((4.0, 6.0, 8.0), (0.5, 0.7, 0.9, 1.05), (False, True), TS, SS)]
    g += [("tec", q, jv, False, t0, S) for q, jv, t0, S in itertools.product((4.0, 8.0), (-16.0, -8.0, -4.0, 4.0, 8.0, 16.0, 30.0), (12.0, 20.0, 30.0), (360, 960, 1560))]
    return g


def evaluate(J, sumw, configs, B, root):
    seeds = np.random.SeedSequence(root).spawn(len(configs)); out = []
    for cfg, ss in zip(configs, seeds):
        mode, q, par, ramp, t0, S = cfg
        rng = np.random.default_rng(ss); t = time.perf_counter()
        s0 = rng.choice(np.array([-1.0, 1.0], dtype=np.float32), size=(B, len(J)))
        s, fl = run_sca(J, s0, S, t0, q, mode, par, ramp, rng)
        c = m.cut(J, sumw, s); k = int((c >= m.TARGET).sum())
        tm = m.cycles(fl.mean(), S, 1.0 if mode == "onsager" else 0.0) / 300e3
        out.append(dict(mode=mode, q=q, par=par, ramp=ramp, T0=t0, S=S, trials=B, successes=k, p=k / B,
                        p_ci95=m.sr.clopper_pearson(k, B), mean_cut=float(c.mean()), mean_flips=float(fl.mean()), t_ms=tm,
                        score={E: m.tts(tm, m.wilson_lower(k, B), E) for E in m.E_LIST},
                        tts={E: m.tts(tm, k / B, E) for E in m.E_LIST}, sec=time.perf_counter() - t))
        print(json.dumps({k2: out[-1][k2] for k2 in ("mode", "q", "par", "ramp", "T0", "S", "p", "t_ms")}), flush=True)
    return out


def main():
    J, sumw = m.load()
    pilot = evaluate(J, sumw, grid(), 256, 40001)
    sel = {}
    for E in m.E_LIST:
        for fam in ("plain", "onsager", "tec"):
            best = min((r for r in pilot if r["mode"] == fam), key=lambda r: r["score"][E])
            sel[f"{fam}_E{E}"] = (best["mode"], best["q"], best["par"], best["ramp"], best["T0"], best["S"])
    uniq = sorted(set(sel.values()))
    hold = evaluate(J, sumw, uniq, 1024, 40002)
    h = {(r["mode"], r["q"], r["par"], r["ramp"], r["T0"], r["S"]): r for r in hold}
    summary = {}
    for E in m.E_LIST:
        row = {fam: dict(cfg=sel[f"{fam}_E{E}"], p=h[sel[f"{fam}_E{E}"]]["p"], p_ci95=h[sel[f"{fam}_E{E}"]]["p_ci95"],
                         t_ms=h[sel[f"{fam}_E{E}"]]["t_ms"], tts_ms=h[sel[f"{fam}_E{E}"]]["tts"][E]) for fam in ("plain", "onsager", "tec")}
        row["onsager_vs_plain"] = row["plain"]["tts_ms"] / row["onsager"]["tts_ms"]
        row["onsager_vs_tec"] = row["tec"]["tts_ms"] / row["onsager"]["tts_ms"]
        summary[E] = row
        print("SUMMARY", E, json.dumps(row), flush=True)
    (ROOT / "J_results.json").write_text(json.dumps(dict(pilot=pilot, selection={k: list(v) for k, v in sel.items()},
                                                         holdout=hold, summary=summary), indent=1) + "\n")


if __name__ == "__main__":
    main()
