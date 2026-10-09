"""Frozen A-C sweeps (PROTOCOL_ABC.md). Usage: python3 sweeps_abc.py {selftest|A|C|B}"""
import itertools
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
PROJ = ROOT.parents[1]
sys.path.insert(0, str(PROJ / "research/statica_reproduction_20261003"))
sys.path.insert(0, str(PROJ / "research/iamp_evaluation_20261003"))
import statica_repro as sr

Q, T_FIN, TARGET, N_DEFAULT = 4.0, 5.0, 33000, 2000
PILOT_SEED, HOLD_SEED = 30001, 30002
E_LIST = (1, 4, 8, 16, 32)


def load():
    Jf, edges = sr.load_graph()
    return Jf.astype(np.float32), float(edges[2].sum())


# ---------- cost model and statistics ----------
def cycles(flips, S, lam=0.0):
    return 0.434 * flips - 0.98 * S + 1989 + (4 * S if lam else 0)


def r99(p):
    return math.inf if p <= 0 else (1.0 if p >= 1 else math.log(0.01) / math.log(1 - p))


def tts(t_ms, p, E=1):
    r = r99(p)
    if not math.isfinite(r):
        return math.inf
    return t_ms * r if E == 1 else t_ms * math.ceil(math.ceil(r) / E)


def wilson_lower(k, n, z=1.959963984540054):
    if n == 0:
        return 0.0
    p = k / n
    return max(0.0, (p + z * z / (2 * n) - z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))) / (1 + z * z / n))


def sched(t0, S):
    return t0 * (T_FIN / t0) ** (np.arange(S) / max(1, S - 1))


def cut(J, sumw, s):
    return (sumw + 0.5 * np.einsum("bi,bi->b", s, s @ J)) / 2


# ---------- A: clipped SCA with Onsager echo correction ----------
def sca(J, s0, S, t0, lam, rng):
    s = s0.astype(np.float32).copy(); h = s @ J
    flips = np.zeros(len(s)); s_prev = None; c_prev = None
    for T in sched(t0, S):
        field = h if (lam == 0 or s_prev is None) else h - lam * c_prev[:, None] * s_prev
        z = s * field + Q
        n_lin = ((z > -2 * T) & (z < 2 * T)).sum(1)
        flip = np.clip(z / (4 * T) + 0.5, 0, 1) < rng.random(s.shape, dtype=np.float32)
        d = np.where(flip, -2 * s, 0).astype(np.float32)
        s_prev, c_prev = s.copy(), n_lin / (2 * T)
        s += d; h += d @ J; flips += flip.sum(1)
    return s, flips


def run_A(J, sumw, configs, B, root):
    seeds = np.random.SeedSequence(root).spawn(len(configs))
    res = []
    for (lam, t0, S), ss in zip(configs, seeds):
        rng = np.random.default_rng(ss); t = time.perf_counter()
        s0 = rng.choice(np.array([-1.0, 1.0], dtype=np.float32), size=(B, len(J)))
        s, fl = sca(J, s0, S, t0, lam, rng)
        c = cut(J, sumw, s); k = int((c >= TARGET).sum())
        tm = cycles(fl.mean(), S, lam) / 300e3
        res.append(dict(lam=lam, T0=t0, S=S, trials=B, successes=k, p=k / B, p_ci95=sr.clopper_pearson(k, B),
                        mean_cut=float(c.mean()), mean_flips=float(fl.mean()), t_ms=tm,
                        score={E: tts(tm, wilson_lower(k, B), E) for E in E_LIST},
                        tts={E: tts(tm, k / B, E) for E in E_LIST}, sec=time.perf_counter() - t))
        print("A", json.dumps({k2: res[-1][k2] for k2 in ("lam", "T0", "S", "p", "mean_flips", "t_ms")}), flush=True)
    return res


def phase_A():
    J, sumw = load()
    grid = list(itertools.product((0.0, 0.25, 0.4, 0.5, 0.6, 0.75), (15.0, 20.0, 30.0, 40.0), (360, 560, 960, 1560)))
    pilot = run_A(J, sumw, grid, 256, PILOT_SEED)
    sel = {}
    for E in E_LIST:
        for fam, ok in (("baseline", lambda r: r["lam"] == 0), ("onsager", lambda r: r["lam"] > 0)):
            best = min((r for r in pilot if ok(r)), key=lambda r: r["score"][E])
            sel[f"{fam}_E{E}"] = (best["lam"], best["T0"], best["S"])
    uniq = sorted(set(sel.values()))
    hold = run_A(J, sumw, uniq, 1024, HOLD_SEED)
    hmap = {(r["lam"], r["T0"], r["S"]): r for r in hold}
    summary = {}
    for E in E_LIST:
        b, o = hmap[sel[f"baseline_E{E}"]], hmap[sel[f"onsager_E{E}"]]
        summary[E] = dict(baseline=dict(cfg=sel[f"baseline_E{E}"], p=b["p"], p_ci95=b["p_ci95"], t_ms=b["t_ms"], tts_ms=b["tts"][E]),
                          onsager=dict(cfg=sel[f"onsager_E{E}"], p=o["p"], p_ci95=o["p_ci95"], t_ms=o["t_ms"], tts_ms=o["tts"][E]),
                          speedup=b["tts"][E] / o["tts"][E] if o["tts"][E] > 0 else None)
        print("A-summary E", E, json.dumps(summary[E]), flush=True)
    (ROOT / "A_results.json").write_text(json.dumps(dict(pilot=pilot, selection={k: list(v) for k, v in sel.items()},
                                                         holdout=hold, summary=summary), indent=1) + "\n")


# ---------- C: population annealing with exact PCA weights ----------
def logf(s, h, T):
    return np.logaddexp(0.0, -(Q + s * h) / T).sum(1)


def energy(s, h):
    return -0.5 * (s * h).sum(1)


def pa(J, sumw, npop, E, S, t0, variant, rng):
    R = npop * E; Nn = len(J)
    s = rng.choice(np.array([-1.0, 1.0], dtype=np.float32), size=(R, Nn)); h = s @ J
    logw = np.zeros(R); flips = np.zeros(R); n_resample = 0; T_prev = None
    for T in sched(t0, S):
        if variant != "IND" and T_prev is not None:
            inc = -energy(s, h) * (1 / T - 1 / T_prev)
            if variant == "PA-exact":
                inc = inc + logf(s, h, T) - logf(s, h, T_prev)
            logw += inc
            lw = logw.reshape(npop, E); w = np.exp(lw - lw.max(1, keepdims=True)); w /= w.sum(1, keepdims=True)
            ess = 1.0 / (w * w).sum(1)
            for p_ in np.nonzero(ess < E / 2)[0]:
                u = (rng.random() + np.arange(E)) / E
                idx = np.minimum(np.searchsorted(np.cumsum(w[p_]), u), E - 1) + p_ * E
                sl = slice(p_ * E, (p_ + 1) * E)
                s[sl] = s[idx]; h[sl] = h[idx]; flips[sl] = flips[idx]; logw[sl] = 0.0; n_resample += 1
        z = s * h + Q
        flip = rng.random(s.shape) < 1.0 / (1.0 + np.exp(np.clip(z / T, -60, 60)))
        d = np.where(flip, -2 * s, 0).astype(np.float32)
        s += d; h += d @ J; flips += flip.sum(1); T_prev = T
    c = cut(J, sumw, s).reshape(npop, E)
    success = (c >= TARGET).any(1)
    wave_cycles = np.array([cycles(f, S) for f in flips.reshape(npop, E).max(1)])
    return success, c, wave_cycles, n_resample


def run_C(J, sumw, configs, npop, root):
    seeds = np.random.SeedSequence(root).spawn(len(configs))
    res = []
    for (E, S, t0, var), ss in zip(configs, seeds):
        rng = np.random.default_rng(ss); t = time.perf_counter()
        succ, c, wc, nres = pa(J, sumw, npop, E, S, t0, var, rng)
        k = int(succ.sum()); tm = float(wc.mean()) / 300e3
        r = r99(k / npop)
        res.append(dict(E=E, S=S, T0=t0, variant=var, populations=npop, successes=k, P_pop=k / npop,
                        P_pop_ci95=sr.clopper_pearson(k, npop), single_replica_p=float((c >= TARGET).mean()),
                        wave_t_ms=tm, waves=math.ceil(r) if math.isfinite(r) else None,
                        tts_ms=tm * math.ceil(r) if math.isfinite(r) else math.inf,
                        score=tm * math.ceil(r99(wilson_lower(k, npop))) if wilson_lower(k, npop) > 0 else math.inf,
                        resamplings_per_pop=nres / npop, sec=time.perf_counter() - t))
        print("C", json.dumps({k2: res[-1][k2] for k2 in ("E", "S", "T0", "variant", "P_pop", "single_replica_p", "wave_t_ms", "resamplings_per_pop")}), flush=True)
    return res


def phase_C():
    J, sumw = load()
    grid = list(itertools.product((16, 32), (200, 360, 560), (20.0, 30.0), ("IND", "PA-exact", "PA-Gibbs")))
    pilot = run_C(J, sumw, grid, 64, PILOT_SEED + 100)
    sel = {}
    for E in (16, 32):
        sel[f"IND_E{E}"] = min((r for r in pilot if r["E"] == E and r["variant"] == "IND"), key=lambda r: r["score"])
        sel[f"PA_E{E}"] = min((r for r in pilot if r["E"] == E and r["variant"] != "IND"), key=lambda r: r["score"])
    uniq = sorted({(r["E"], r["S"], r["T0"], r["variant"]) for r in sel.values()})
    # Also hold out all three variants at each selected PA schedule (paired comparison).
    for E in (16, 32):
        r = sel[f"PA_E{E}"]
        uniq += [(E, r["S"], r["T0"], v) for v in ("IND", "PA-exact", "PA-Gibbs")]
    uniq = sorted(set(uniq))
    hold = run_C(J, sumw, uniq, 256, HOLD_SEED + 100)
    (ROOT / "C_results.json").write_text(json.dumps(dict(pilot=pilot, selection={k: [v["E"], v["S"], v["T0"], v["variant"]] for k, v in sel.items()},
                                                         holdout=hold), indent=1) + "\n")


# ---------- B: IAMP(t*) + SCA finish ----------
def phase_B():
    import quick_tests as qt  # IAMP recurrence (identical to iamp_eval) with early stop
    import iamp_eval as ie
    J, sumw = load(); Nn = len(J); table = ie.load_table()

    def run(configs, B, root):
        out = []
        seeds = np.random.SeedSequence(root).spawn(64)
        cache = {}
        for i, (K, frac, Ta, Sf, lam) in enumerate(configs):
            key = (K, frac)
            if key not in cache:
                g = np.random.default_rng(seeds[len(cache)]).normal(size=(B, Nn))
                kstop = int(round(frac * K))
                cache[key] = (np.where(qt.iamp_stop(g, K, kstop, table) >= 0, 1.0, -1.0).astype(np.float32), kstop)
            s0, kstop = cache[key]
            rng = np.random.default_rng(np.random.SeedSequence([root, i]))
            s, fl = sca(J, s0, Sf, Ta, lam, rng)
            c = cut(J, sumw, s); k = int((c >= TARGET).sum())
            cyc = (kstop - 1) * (0.434 * Nn + 100) + cycles(fl.mean(), Sf, lam)
            tm = cyc / 300e3
            out.append(dict(K=K, t_star=frac, Ta=Ta, Sf=Sf, lam=lam, trials=B, successes=k, p=k / B, p_ci95=sr.clopper_pearson(k, B),
                            mean_cut=float(c.mean()), iamp_steps=kstop, finish_flips=float(fl.mean()), t_ms=tm,
                            score={E: tts(tm, wilson_lower(k, B), E) for E in (1, 16)}, tts={E: tts(tm, k / B, E) for E in (1, 16)}))
            print("B", json.dumps({k2: out[-1][k2] for k2 in ("K", "t_star", "Ta", "Sf", "lam", "p", "t_ms")}), flush=True)
        return out

    grid = list(itertools.product((64, 128), (0.8, 0.9, 1.0), (6.0, 8.0, 12.0), (100, 200, 400), (0.0, 0.5)))
    pilot = run(grid, 256, PILOT_SEED + 200)
    sel = {E: min(pilot, key=lambda r: r["score"][E]) for E in (1, 16)}
    uniq = sorted({(r["K"], r["t_star"], r["Ta"], r["Sf"], r["lam"]) for r in sel.values()})
    hold = run(uniq, 1024, HOLD_SEED + 200)
    (ROOT / "B_results.json").write_text(json.dumps(dict(pilot=pilot, selection={E: [v["K"], v["t_star"], v["Ta"], v["Sf"], v["lam"]] for E, v in sel.items()},
                                                         holdout=hold), indent=1) + "\n")


def selftest():
    """Exact enumeration on a 10-spin graph: PA-exact weights equal log pi_PCA ratios; logistic kernel leaves pi_PCA invariant."""
    rng = np.random.default_rng(5); n = 10
    W = np.triu(rng.choice([-1.0, 1.0], size=(n, n)), 1); Js = -(W + W.T)
    S = np.array(list(itertools.product([-1.0, 1.0], repeat=n)))
    H = S @ Js; E_ = -0.5 * (S * H).sum(1)
    for Ta, Tb in ((7.0, 5.5), (3.0, 2.2)):
        def logpi(T):
            lp = -E_ / T + logf(S, H, T)
            return lp - np.logaddexp.reduce(lp)
        inc = -E_ * (1 / Tb - 1 / Ta) + logf(S, H, Tb) - logf(S, H, Ta)
        diff = (logpi(Tb) - logpi(Ta)) - inc
        assert np.ptp(diff) < 1e-9, "PA-exact weight mismatch"
        pf = 1 / (1 + np.exp((S * H + Q) / Tb))  # logistic kernel at Tb
        P = np.ones((len(S), len(S)))
        for a in range(len(S)):
            fl = S[a] != S
            P[a] = np.prod(np.where(fl, pf[a], 1 - pf[a]), 1)
        pi = np.exp(logpi(Tb))
        assert np.max(np.abs(pi @ P - pi)) < 1e-12, "logistic PCA not stationary for pi_PCA"
    return "PASS: PA-exact incremental weights equal exact log pi_PCA ratios; logistic kernel leaves pi_PCA invariant (10 spins, 1024 states)"


if __name__ == "__main__":
    {"selftest": lambda: print(selftest()), "A": phase_A, "B": phase_B, "C": phase_C}[sys.argv[1]]()
