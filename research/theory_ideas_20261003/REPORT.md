# Theory-motivated improvements to SCA: first exploratory tests

3 October 2026. Laptop CPU, WK2000_1, target cut ≥ 33,000, uniform random thresholds, q = 4. Each setting uses 1,024 runs (idea A) or 512 runs (idea B), from a single seed set. **Exploratory: no frozen protocol and no held-out cohort.** Script: [quick_tests.py](quick_tests.py); results: [quick_tests.json](quick_tests.json).

## A. Onsager-corrected SCA (dynamical TAP for synchronous updates)

- **Mechanism.** In parallel updates, spin i's value two steps ago echoes back into its own field through every coupled spin j. The linear-response size of that echo is `c_i = Σ_j J_ij² χ_j`, where `χ_j = 1/(2T)` inside Eq. 7's random band and 0 outside. For ±1 couplings, `c = n_lin/(2T)`.
- **Test.** Make decisions with `h_i − λ·c·s_i(t−1)`. Hardware cost: one stored spin vector, one popcount and one add per spin.

| Setting | λ | Success | Flips per run | Run time (cycle model) | TTS, 1 engine | Runs for 99% |
|---|---:|---:|---:|---:|---:|---:|
| Short (S=560, T 30→5) | 0 (SCA) | 0.133 | 86,693 | 0.130 ms | 4.21 ms | 33 |
|  | **0.5** | **0.498** | 129,304 | 0.192 ms | **1.28 ms** | **7** |
|  | 1.0 | 0 (diverges) | 459,354 | — | — | — |
|  | −0.5 | 0.017 | 72,388 | 0.110 ms | 30.1 ms | 276 |
| Long (S=1560, T 40→5) | 0 (SCA) | 0.830 | 329,526 | 0.478 ms | 1.24 ms | 3 |
|  | **0.5** | **0.913** | 472,598 | 0.685 ms | 1.29 ms | 2 |

- **Success.** Partial correction (λ = 0.5) raises success sharply at the same number of steps: 3.7× at the short setting. It is not only extra noise. At matched flip cost, short λ = 0.5 (129k flips, p = 0.50) beats plain SCA with T_init = 40 (127k flips, p = 0.27).
- **Sign.** Adding the echo instead (λ < 0) hurts. Full subtraction (λ = 1) is unstable.
- **TTS so far.** With the flip-based cycle model, the extra flips roughly cancel the gain, because the schedules were not retuned for λ.
- **Open questions:**
  - Why λ* ≈ 1/2? A discrete-spin (non-linear-response) derivation should predict it.
  - Can a cooler or shorter schedule keep the success while cutting flips?

## B. IAMP stopped at t* < 1, then a short SCA finish

IAMP's late collapse is removed. At K = 128, success goes from 0.137 with sign rounding to 0.58–0.68 with a 200–400-step SCA finish at T 8–12 → 5. At K = 256 with t* = 0.9 it reaches 0.67–0.76.

The cost is 115–230 dense multi-bit matrix-vector products. On the flip-sparse STATICA-class cycle model that is worse than SCA. On GPU it is estimated at roughly 1.5× better than SCA's 61 ms single-chain TTS (not measured).

## Further directions not yet tested

- **C. Exact sequential-Monte-Carlo weights for parallel engines.** Use the identity `π_PCA ∝ π_T / P_hold`: the correct resampling weight between temperatures is `exp(−ΔβH) · P_hold,T(s) / P_hold,T'(s)`. This enables population annealing across E hardware engines with weights valid for parallel kernels, instead of independent restarts.
- **D. DMFT-optimal schedules.** For dense couplings, the N→∞ dynamics of synchronous SCA is described by a single-site effective process (discrete-time generating functional). Optimizing T(t), q(t) and λ(t) offline costs nothing on hardware, and could also explain λ*. A 2026 preprint already does this for analog solvers, so the SCA-specific version needs a precise novelty statement.

## Next steps

1. Run a frozen λ × (T_init, S) sweep on pilot seeds, with held-out confirmation. Report TTS under the cycle model and on GPU.
2. Derive the Onsager coefficient for binary stochastic parallel dynamics (direction D) and test whether it predicts λ* and the schedule.
