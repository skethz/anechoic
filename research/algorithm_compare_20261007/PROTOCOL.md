# Algorithm-level comparison on K2000: energy vs Monte Carlo step (frozen 7 October 2026, before any pilot of this protocol)

## Purpose

Ising-machine papers compare algorithms by the Ising energy H(s) = −½ sᵀJs reached after a given number of Monte Carlo steps (MCS). This study makes that comparison for our two methods and the published methods in the paper's comparison set, on the dense K2000 instance (WK2000_1, N = 2000) used for every TTS result.

- H = Σw − 2·cut, so cut ≥ 33,000 ⇔ H ≤ −67,040, and the best-known cut 33,337 ⇔ H = −67,714.
- The runs before this protocol (`tune.py`, `tune_sa_ext.json`, `final.py` at S = 960) were exploratory. They are kept, but nothing below is selected from them.

## Step definition (one step = each spin gets one update opportunity; all O(N²) work)

| Method | Source | One step |
|---|---|---|
| SA | Kirkpatrick et al. 1983 | one sequential sweep: N single-spin Metropolis updates in index order, geometric T0 → T1 |
| SCA | Yamamoto et al., JSSC 2021 (STATICA) | one synchronous update of all spins: clipped flip probability, constant self-coupling q, geometric T0 → 5 |
| TEC | Du et al., arXiv:2608.21753 (2026) | SCA step with field h + J_v s(t−1) |
| APC-SCA | Okonogi et al., IEICE E106-D (2023) | SCA step with per-spin q_i: reset to q_reset after a flip, otherwise q_i ← max(r_q·q_i, q_lim) |
| ReAIM ASA | Chiang et al., ISCA 2024 (Algs. 2–3) | one iteration: thresholded selection of at most k flips. Trial phases run one replica per candidate k (counted once, as the paper counts iterations). Implementation: `research/reaim_reproduction_20261003/pilot.py`, noise-free, ITER_trial 32, ITER_run 96 |
| aSB | Goto et al., Sci. Adv. 2019 | one integration step: M = 5 Kerr sub-steps, then one J·x product. Pump p = t/S, x(0) = 0, y(0) = 0.1·u·s₀. Same integrator as `gpu/prior_compare_20260929/sb_gpu.py` |
| bSB, dSB | Goto et al., Sci. Adv. 2021 | one symplectic-Euler step with inelastic walls; dSB uses sign(x) in the coupling term; a(t) linear 0 → 1; c0 = ξ·0.5/(σ√N) |
| TEC-T (ours) | this work | SCA step with field h − κ·T(t−1)·ramp(t)·s(t−1) |
| Onsager SCA (ours) | this work | SCA step with field h − λ·ramp(t)·c(t−1)·s(t−1), with c = n_lin/(2T) |

- The SCA family uses the float model `research/ablation_20261005/abl.py`. That is the same arithmetic used for every schedule run on the V80, where the engine is bit-exact to the fixed-point version.
- `ramp` means the coefficient falls linearly to 0 over the last 30% of steps, as on the hardware.

## Hyperparameter grids (each method gets its own grid; ours are not larger than the largest baseline grid)

| Method | Grid | Size |
|---|---|---|
| SA | T0 ∈ {20, 30, 50, 80} × T1 ∈ {0.5, 1, 2, 3} | 16 |
| SCA | q ∈ {2, 4, 6, 8} × T0 ∈ {15, 20, 30, 40} | 16 |
| TEC | J_v ∈ {−4, −2, +2, +4} × q ∈ {6, 8} × T0 ∈ {20, 30} | 16 |
| APC-SCA | q_reset ∈ {8, 16, 32} × r_q ∈ {0.9, 0.97} × q_lim ∈ {2, 4} × T0 ∈ {15, 20, 30} | 36 |
| ReAIM ASA | k set ∈ {{64,128,256,512}, {128,256,512,1024}, {256,512,1024,2000}} × T 1 → {0.1, 0.05} | 6 |
| aSB | dt ∈ {0.5, 0.9, 1.25} × ξ = {0.5, 1, 1.5} × 0.7/√N | 9 |
| bSB, dSB | dt ∈ {0.5, 0.75, 1.0, 1.25} × ξ ∈ {0.5, 1, 2} | 12 each |
| TEC-T (ours) | κ ∈ {1.25, 1.75, 2.0} × q ∈ {6, 8} × T0 ∈ {12, 15, 20}, ramp | 18 |
| Onsager (ours) | λ ∈ {0.7, 0.9, 1.05} × q ∈ {6, 8} × T0 ∈ {12, 15, 20}, ramp | 18 |

## Budgets, selection and final runs

- **Budgets.** S ∈ {250, 500, 1000, 2000, 4000} steps.
- **Pilot.** For each method, S and grid point: 64 runs, seed 80000 + 1000·(S index) + grid index (method-specific offset 100000·method index).
- **Selection rule (fixed).** Highest pilot mean final cut, with ties going to the lower grid index. This is one rule for every method; P(cut ≥ 33,000) is reported but not used to select.
- **Final (held out).** The selected configuration per (method, S), 256 runs, seed 90000 + 1000·(S index) + 100000·(method index), disjoint from all pilot seeds.
- **Reported, final runs only.**
  - Mean energy trajectory H(t) at every step, with 10th and 90th percentiles; the figure uses S = 1000.
  - Final mean cut ± sd.
  - P(cut ≥ 33,000) with a Wilson 95% interval.
  - P(cut ≥ 33,200).
  - Steps-to-solution MCS99 = S·ln(0.01)/ln(1 − p), or S if p = 1; the minimum over S is reported per method.
- **Final state.** Success uses the final state, as on the hardware. The best-visited state is not used.

## What this does not show

- **Steps are not time.** An SA sweep is N sequential updates, while an SCA step is one parallel update. bSB, dSB, aSB and ReAIM steps each need a dense real-valued (or analog) J·x product. On the V80, our SCA-family engine takes about 18 cycles per step.
- **Hardware time is reported separately.** The paper reports measured wall-clock TTS separately, and this study makes no time claim.
