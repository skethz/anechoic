"""A2: confirmatory test of theory-D predictions at N=2000 (PROTOCOL_A2.md)."""
import json
import time
from pathlib import Path

import numpy as np

import sweeps_abc as m

ROOT = Path(__file__).resolve().parent
SCHED = {"short": (560, 30.0), "long": (1560, 40.0)}


def sca_logged(J, s0, S, t0, q, lam_fn, rng):
    s = s0.astype(np.float32).copy(); h = s @ J
    flips = np.zeros(len(s)); back = 0; total = 0; prev_flip = None
    s_prev = None; c_prev = None; log = []
    for t, T in enumerate(m.sched(t0, S)):
        lam = lam_fn(t, S)
        field = h if (lam == 0 or s_prev is None) else h - lam * c_prev[:, None] * s_prev
        z = s * field + q
        n_lin = ((z > -2 * T) & (z < 2 * T)).sum(1)
        flip = np.clip(z / (4 * T) + 0.5, 0, 1) < rng.random(s.shape, dtype=np.float32)
        if prev_flip is not None:
            back += int((flip & prev_flip).sum()); total += int(flip.sum())
        d = np.where(flip, -2 * s, 0).astype(np.float32)
        s_prev, c_prev, prev_flip = s.copy(), n_lin / (2 * T), flip
        s += d; h += d @ J; flips += flip.sum(1)
        if t % 20 == 0 or t == S - 1:
            x = float(np.mean(n_lin / (4 * T * T)))
            log.append(dict(t=t, T=float(T), x=x, b=lam * x, flips=float(flip.sum(1).mean())))
    return s, flips, back / max(1, total), log


def main():
    J, sumw = m.load()
    out = {}
    for sname, (S, t0) in SCHED.items():
        init = np.random.default_rng(np.random.SeedSequence([30003, S])).choice(np.array([-1.0, 1.0], dtype=np.float32), size=(1024, len(J)))
        const = lambda lam: (lambda t, S_: lam)
        ramp = lambda lam: (lambda t, S_: lam if t < 0.7 * S_ else lam * (S_ - 1 - t) / max(1, (S_ - 1 - 0.7 * S_)))
        configs = [(f"q4_lam{l}", 4.0, const(l), l) for l in (0.0, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)]
        configs += [(f"q4_lam{l}_ramp", 4.0, ramp(l), l) for l in (0.5, 0.7)]
        configs += [("q8_lam0.0", 8.0, const(0.0), 0.0), ("q8_lam0.9", 8.0, const(0.9), 0.9)]
        res = {}
        for i, (name, q, fn, lam_nom) in enumerate(configs + [("rule", 4.0, None, None)]):
            if name == "rule":  # P6: lambda = 0.75/(4 x0), x0 = median x over T in [8,25] from the lambda=0 run
                lg = res["q4_lam0.0"]["log"]
                x0 = float(np.median([e["x"] for e in lg if 8 <= e["T"] <= 25]))
                lam_nom = 0.75 / (4 * x0); fn = ramp(lam_nom)
            rng = np.random.default_rng(np.random.SeedSequence([30003, S, i])); t = time.perf_counter()
            s, fl, back, log = sca_logged(J, init, S, t0, q, fn, rng)
            c = m.cut(J, sumw, s); k = int((c >= m.TARGET).sum())
            tm = m.cycles(fl.mean(), S, lam_nom) / 300e3
            res[name] = dict(q=q, lam=lam_nom, successes=k, p=k / 1024, p_ci95=m.sr.clopper_pearson(k, 1024),
                             mean_cut=float(c.mean()), mean_flips=float(fl.mean()), flip_back=back, t_ms=tm,
                             tts_ms={E: m.tts(tm, k / 1024, E) for E in m.E_LIST}, log=log, sec=time.perf_counter() - t)
            if name == "rule":
                res[name]["x0"] = x0
            print(sname, name, json.dumps({k2: res[name][k2] for k2 in ("lam", "p", "mean_cut", "mean_flips", "flip_back", "t_ms")}), flush=True)
        out[sname] = res
    (ROOT / "A2_results.json").write_text(json.dumps(out, indent=1) + "\n")


if __name__ == "__main__":
    main()
