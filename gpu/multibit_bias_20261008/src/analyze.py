#!/usr/bin/env python3
"""Estimators for the GPU v2 study (PROTOCOL.md section 6), computed from the raw per-run files.

primary   : P_batch = fraction of batches with >= 1 success; TTS99 = t_batch * ln(0.01) / ln(1 - P_batch) (P_batch = 1 -> t_batch)
secondary : p = per-trial success fraction; P = 1 - (1 - p)^B; same formula
t_batch   : mean CUDA-event device time of one kernel launch (one batch of B chains)
Intervals : Wilson 95% on P_batch and p, mapped through the (monotone) TTS formula.
"""
import json, math, sys, os, glob

Z = 1.959963984540054
LN01 = math.log(0.01)


def wilson(k, n):
    if n == 0:
        return (0.0, 1.0)
    ph = k / n
    den = 1 + Z * Z / n
    cen = (ph + Z * Z / (2 * n)) / den
    half = Z * math.sqrt(ph * (1 - ph) / n + Z * Z / (4 * n * n)) / den
    return (max(0.0, cen - half), min(1.0, cen + half))


def tts(P, t):
    if P >= 1.0:
        return t
    if P <= 0.0:
        return math.inf
    return t * LN01 / math.log(1.0 - P)


def tts_sec(p, B, t):
    # Amendment 1: secondary TTS99 = t * ln(0.01) / (B * ln(1 - p)) evaluated directly (no P = 1 floor, no 1-(1-p)^B rounding).
    if p <= 0.0:
        return math.inf
    if p >= 1.0:
        return 0.0
    return t * LN01 / (B * math.log(1.0 - p))


def tts_whole(P, t):
    if P >= 1.0:
        return t
    if P <= 0.0:
        return math.inf
    return t * max(1, math.ceil(LN01 / math.log(1.0 - P)))


def load_run(prefix):
    s = json.load(open(prefix + ".summary.json"))
    bat = [json.loads(l) for l in open(prefix + ".batches.jsonl")]
    return s, bat


def estimate(prefix, split=None):
    """split: evaluate P_batch on sub-batches of `split` chains (contiguous trial ids) -- used only for documentation."""
    s, bat = load_run(prefix)
    B, nb = s["chains"], s["batches"]
    t = sum(b["device_ms"] for b in bat) / nb
    th = sum(b["host_ms"] for b in bat) / nb
    k_b = sum(1 for b in bat if b["successes"] > 0)
    k_t = s["successes"]
    n_t = B * nb
    Pb = k_b / nb
    p = k_t / n_t
    Pb_lo, Pb_hi = wilson(k_b, nb)
    p_lo, p_hi = wilson(k_t, n_t)
    sec = lambda pp: tts_sec(pp, B, t)
    out = dict(prefix=os.path.basename(prefix), cs=s["cs"], groups=s["groups"], B=B, batches=nb, steps=s["steps"], t0=s["t0"], q=s["q"],
               lam=s["lambda"], ramp=s["ramp"], tec_jv=s["tec_jv"], kappa=s["tecT_kappa"], trial_offset=s["trial_offset"],
               successes=k_t, trials=n_t, p=p, p_wilson95=[p_lo, p_hi], batches_with_success=k_b, P_batch=Pb, P_batch_wilson95=[Pb_lo, Pb_hi],
               t_batch_ms=t, t_batch_min_ms=s["device_ms_min"], t_batch_max_ms=s["device_ms_max"], host_ms=th,
               us_per_step=1000.0 * t / s["steps"],
               tts_primary_ms=tts(Pb, t), tts_primary_ci_ms=[tts(Pb_hi, t), tts(Pb_lo, t)], tts_primary_whole_ms=tts_whole(Pb, t),
               tts_secondary_ms=sec(p), tts_secondary_ci_ms=[sec(p_hi), sec(p_lo)],
               cut_mismatches=s["cut_mismatches"], gpu_before=s["gpu_before"], gpu_after=s["gpu_after"])
    return out


def selection_score(e, which):
    """Pre-registered pilot selection score: TTS at the Wilson 95% lower bound of P_batch (primary) or p (secondary)."""
    if which == "primary":
        return tts(e["P_batch_wilson95"][0], e["t_batch_ms"])
    return tts_sec(e["p_wilson95"][0], e["B"], e["t_batch_ms"])


if __name__ == "__main__":
    rows = [estimate(p[: -len(".summary.json")]) for p in sorted(sys.argv[1:]) if p.endswith(".summary.json")]
    for r in rows:
        print(json.dumps(r))
