"""POST-HOC diagnostic, added after the first reproduction result exceeded the published short-point quality.

It tests whether STATICA's documented random-number circuit (JSSC 2021, Sec. IV-C2, Fig. 11) changes success.
The circuit works as follows:
- Take a 16-bit random value and AND it with T_MASK, which is all ones from the top set bit of 2T downward.
- If the result exceeds 2T, output its XOR with T_MASK (folding). Then attach a random sign, giving a threshold in [-2T, 2T].
- A spin flips if (h*sigma + q) < r, which is Eq. (7) with the threshold drawn from [-2T, 2T].

Assumptions not stated in the paper, so they are stated here:
- T and the threshold are in Q8.8 fixed point (resolution 1/256).
- The sign bit is independent of the 16 masked bits.
- Ideal uniform 16-bit inputs are used instead of XorShift32, so this isolates folding and quantization only.
- T is rounded to Q8.8 after every multiplication by r_T.
"""
import json
import sys
from pathlib import Path

import numpy as np

import statica_repro as sr

ROOT = Path(__file__).resolve().parent
FRAC = 8  # Q8.8


def folded_threshold(rng, shape, T):
    two_t = int(round(2 * T * (1 << FRAC)))
    mask = (1 << two_t.bit_length()) - 1
    v = rng.integers(0, 1 << 16, size=shape, dtype=np.int64) & mask
    v = np.where(v > two_t, v ^ mask, v)
    sign = np.where(rng.integers(0, 2, size=shape, dtype=np.int8) == 1, 1, -1)
    return (sign * v).astype(np.float64) / (1 << FRAC)


def sca_hw(J, S, t_init, t_fin, q, sigma0, rng, folded=True):
    sigma = sigma0.astype(np.float32).copy()
    h = sigma @ J
    r_T = (t_fin / t_init) ** (1.0 / (S - 1))
    T = round(t_init * (1 << FRAC)) / (1 << FRAC)
    for s in range(S):
        z = h * sigma + q
        if folded:
            flip = z < folded_threshold(rng, sigma.shape, T)
        else:  # quantized uniform threshold only (no folding), to separate the two effects
            r = np.floor(rng.uniform(-2 * T, 2 * T, size=sigma.shape) * (1 << FRAC)) / (1 << FRAC)
            flip = z < r
        d = np.where(flip, -2.0 * sigma, 0.0).astype(np.float32)
        sigma += d
        h += d @ J
        T = round(T * r_T * (1 << FRAC)) / (1 << FRAC)
    return sigma


def main():
    J, edges = sr.load_graph()
    sumw = float(edges[2].sum())
    seeds = np.random.SeedSequence(20261004).spawn(4)
    out, c = {}, 0
    for name, (S, t_init) in sr.TABLE2.items():
        for folded in (True, False):
            rng = np.random.default_rng(seeds[c]); c += 1
            sig = np.concatenate([sca_hw(J, S, t_init, sr.T_FIN, sr.Q,
                                         rng.choice(np.array([-1.0, 1.0], dtype=np.float32), size=(1024, 2000)),
                                         rng, folded) for _ in range(2)])
            H = sr.energies(J, sig)
            cut = (sumw - H) / 2
            k = int((cut >= sr.TARGET_CUT).sum())
            key = f"{name}_{'folded_Q8.8' if folded else 'uniform_Q8.8'}"
            out[key] = dict(S=S, T_init=t_init, trials=len(sig), successes=k, p=k / len(sig),
                            p_ci95=sr.clopper_pearson(k, len(sig)), mean_cut=float(cut.mean()), sd_cut=float(cut.std(ddof=1)))
            print(key, json.dumps(out[key]), flush=True)
    (ROOT / "hw_rng_sensitivity.json").write_text(json.dumps(out, indent=1) + "\n")


def grid():
    """Folded-threshold model on the full experiment-3 grid (1,024 trials each), for comparison with Fig. 22/24."""
    J, edges = sr.load_graph()
    sumw = float(edges[2].sum())
    seeds = np.random.SeedSequence(20261006).spawn(len(sr.GRID_S) * len(sr.GRID_TINIT))
    out, c = {}, 0
    for t_init in sr.GRID_TINIT:
        for S in sr.GRID_S:
            rng = np.random.default_rng(seeds[c]); c += 1
            sig = sca_hw(J, S, t_init, sr.T_FIN, sr.Q,
                         rng.choice(np.array([-1.0, 1.0], dtype=np.float32), size=(1024, 2000)), rng, True)
            H = sr.energies(J, sig)
            cut = (sumw - H) / 2
            k = int((cut >= sr.TARGET_CUT).sum())
            key = f"grid_S{S}_T{t_init:g}"
            out[key] = dict(S=S, T_init=t_init, trials=len(sig), successes=k, p=k / len(sig),
                            p_ci95=sr.clopper_pearson(k, len(sig)), mean_H=float(H.mean()), sd_H=float(H.std(ddof=1)))
            print(key, json.dumps(out[key]), flush=True)
    (ROOT / "hw_rng_grid.json").write_text(json.dumps(out, indent=1) + "\n")


if __name__ == "__main__":
    grid() if "--grid" in sys.argv else main()
