# Theory-based improvements to SCA: results for A–D

3 October 2026. Laptop CPU, WK2000_1, target cut ≥ 33,000, final-state success.

- Protocols were frozen before the runs: [PROTOCOL_ABC.md](PROTOCOL_ABC.md), and [PROTOCOL_A2.md](PROTOCOL_A2.md), whose predictions were registered after the theory but before the N=2000 test.
- TTS uses the STATICA-class cycle model, `cycles = 0.434·flips − 0.98·S + 1989` at 300 MHz, plus 4 cycles per step when the correction is on. Device TTS with E engines is `t_run·⌈⌈R99⌉/E⌉`.
- These are model-based hardware times with ideal uniform random numbers, not measured chips.

## Summary

| Idea | Theory basis | Held-out result | Status |
|---|---|---|---|
| **A. Onsager-corrected SCA** | Echo in synchronous dynamics; the correction's strength b = λx is the hard-tanh AT/AMP parameter | **1.34–1.85× lower TTS than plain SCA at every engine count.** Best single engine: 0.69 ms vs 1.27 ms. | Strong; novelty must be distinguished from TEC (same functional term) |
| **C. Population annealing with exact PCA weights** | The identity `π_PCA ∝ π_T / P_hold` | With the same logistic kernel: P(population succeeds) 0.60 → 0.90 (16 engines) and 0.90 → 0.98 (32 engines); **1.7–3.3× lower TTS** than independent engines. Standard Gibbs weights perform *worse* than independent engines. | Strong mechanism. The logistic kernel is weaker than clipped SCA, so absolute TTS does not yet beat A. |
| **B. IAMP(t*) + SCA finish** | Structure of the El Alaoui–Montanari–Sellke theorem | 1 engine: 1.56 ms. 16 engines: 0.175 ms. | Not competitive on flip-sparse hardware |
| **D. Theory of A** | Linear stability, runaway through n_lin | Predictions P1–P5 confirmed at N=2000; the P6 rule ties the best constant λ but does not beat it | Explains λ_c and the q–λ coupling; does not yet give the optimal schedule |

## A. Onsager-corrected SCA, frozen sweep

96 pilot configurations (256 runs each). Selection used the Wilson lower bound on p. 1,024 held-out runs for each selected configuration.

| Engines | Plain SCA: config (λ, T_init, S), p, TTS | Onsager SCA: config, p, TTS | Speed-up |
|---:|---|---|---:|
| 1 | (0, 30, 1560), 0.671, 1.274 ms | (0.75, 20, 960), 0.866, **0.688 ms** | 1.85× |
| 4 | (0, 40, 1560), 0.824, 0.478 ms | (0.4, 30, 960), 0.736, 0.289 ms | 1.65× |
| 8 | (0, 40, 960), 0.615, 0.305 ms | (0.75, 20, 560), 0.620, 0.191 ms | 1.60× |
| 16 | (0, 30, 960), 0.396, 0.202 ms | (0.75, 15, 560), 0.397, 0.133 ms | 1.51× |
| 32 | (0, 30, 560), 0.143, 0.130 ms | (0.75, 15, 360), 0.240, 0.097 ms | 1.34× |

**Mechanism.** The correction allows a much cooler start (T 15–20 instead of 30–40), which removes the expensive flip-heavy early phase. It also reduces immediate flip-backs: 67% of flips at λ=0, 54% at λ=0.7.

## A2. Pre-registered test of theory-D predictions (1,024 runs per configuration, paired starts)

| ID | Prediction | Short (S=560, T 30→5) | Long (S=1560, T 40→5) | Verdict |
|---|---|---|---|---|
| P1 | x ≈ 0.3, flat over T ∈ [8, 25] at λ=0 | 0.329 | 0.309 | **Pass** |
| P2 | λ=1 never orders; b crosses 1 at T ≈ 16 | T=15.8; cut 9,721 | T=15.7; cut 9,971 | **Pass** |
| P3 | λ_c ≈ 0.8–0.9; λ=0.7 ≥ 0.5 | p: 0.5→0.517, 0.7→0.606, 0.8→0.612, **0.9→0.235**, 1→0 | 0.922, 0.940, 0.928, **0.757**, 0 | **Pass** (threshold in (0.8, 0.9]) |
| P4 | q8/λ0.9 > q4/λ0.5; q8/λ0 < q4/λ0 | 0.635 > 0.517; 0.057 < 0.118 | 0.970 > 0.922; 0.790 vs 0.824 (overlap) | **Pass** (one comparison unresolved) |
| P5 | Ramping λ to 0 over the last 30% helps | λ0.7: 0.699 vs 0.606 | λ0.7: 0.978 vs 0.940 | **Pass** at 0.7; unresolved at 0.5 |
| P6 | Rule λ = 0.75/(4x₀) + ramp beats the best constant λ | 0.619 vs 0.612 (overlap) | 0.957 vs 0.940 (overlap) | **Not confirmed** (ties) |

**Best single-engine TTS at these two fixed schedules**
- Short: λ0.7 with ramp, 1.000 ms (plain 4.76 ms).
- Long: q=8 with λ=0, 0.968 ms (plain 1.27 ms).

Both higher q and the correction are levers, and the frozen A grid fixed q=4. **A joint (q, λ, ramp, T_init, S) optimization is the obvious next sweep.**

## C. Population annealing with exact PCA weights (logistic kernel, 64 pilot / 256 held-out populations)

| Engines | Schedule | Independent: P_pop, TTS | Gibbs-weight PA | **Exact-weight PA** |
|---:|---|---|---|---|
| 16 | S=560, T 30→5 | 0.602, 1.60 ms | 0.523, 1.60 ms | **0.902, 0.478 ms** |
| 32 | S=560, T 30→5 | 0.895, 0.80 ms | 0.742, 0.90 ms | **0.977, 0.475 ms** |

- The correct stationary law of the *parallel* kernel is essential: naive Gibbs weights hurt.
- Resampling happens about 26–30 times per run. The cost of copying states between engines is not charged, which is optimistic.
- **Limitation:** exact weights require the logistic PCA kernel, which is weaker per run than clipped SCA at these schedules. Clipped SCA and Onsager SCA have no closed-form stationary law. Extending C needs either a better-tuned logistic schedule or approximate weights for corrected kernels.

## B. IAMP(t*) + SCA finish (frozen, 108 pilot configurations, 1,024 held-out runs)

| Target | Selected configuration | p | t | TTS |
|---|---|---|---|---|
| 1 engine | K=128, t*=0.9, finish T 12→5 over 400 steps, λ=0.5 | 0.696 | 0.402 ms | 1.56 ms |
| 16 engines | K=64, t*=0.8, finish T 12→5 over 200 steps | 0.259 | 0.175 ms | 0.175 ms |

These times assume optimistic hypothetical multi-bit multiply-accumulate hardware. Dense IAMP steps dominate the cost.

## D. Theory

See [REPORT_D.md](REPORT_D.md).

- The echo coefficient is `c = n_lin/(2T)`.
- The correction's strength is b = λx, with x the hard-tanh AT/AMP parameter, so λ=1 is a stochastic AMP iteration.
- The linear-stability window `x ∈ (x_ord(λ), 1/λ)` closes at λ=1.
- Runaway through n_lin gives λ_c ≈ 1/(4x₀).
- λ ≈ ½ is not principled; the invariant is b.

## Novelty and next steps

**Prior art.** TEC (arXiv:2608.21753) uses the same σ(t)σ(t−1) term, and APC-SCA uses functionally similar pinning. A defensible contribution is:
1. The derivation of the coefficient and the AMP/AT identity.
2. The predicted and confirmed threshold λ_c and q–λ coupling.
3. Exact-weight population annealing for parallel kernels.
4. Hardware TTS gains under a calibrated cycle model.

**Next steps**
- Run a joint (q, λ, ramp, T_init, S) frozen sweep.
- Compare against TEC's published rule on the same instance.
- Extend C to corrected and clipped kernels.
- Test fresh graphs and other N.
- Implement Onsager SCA on the V80 engine; the hardware cost is one stored spin vector, one popcount and one add per spin.

## J. Joint (q, λ, ramp, T_init, S) sweep and TEC comparison ([PROTOCOL_J.md](PROTOCOL_J.md))

558 pilot configurations (256 runs each); 1,024 held-out runs for each selected configuration. TEC (arXiv:2608.21753) is embedded in the same parallel SCA with a constant temporal coupling J_v, in both signs.

| Engines | Plain SCA (best q, T_init, S) | TEC (constant J_v) | **Onsager SCA** | vs plain | vs TEC |
|---:|---|---|---|---:|---:|
| 1 | q4 T30 S1560, p=0.675, 1.261 ms | q8 J_v=−4 T30 S1560, p=0.675, 1.010 ms | **q8 λ1.05+ramp T12 S960, p=0.937, 0.348 ms** | 3.62× | 2.90× |
| 4 | q6 T30 S1560, 0.478 ms | 0.480 ms | **q6 λ0.9+ramp T15 S560, 0.139 ms** | 3.44× | 3.46× |
| 8 | q6 T30 S1560, 0.239 ms | 0.246 ms | **q6 λ0.9+ramp T12 S560, 0.097 ms** | 2.48× | 2.55× |
| 16 | q8 T30 S1560, 0.190 ms | 0.160 ms | **q6 λ0.9 T12 S360, 0.076 ms** | 2.50× | 2.10× |
| 32 | q8 T30 S960, 0.126 ms | 0.148 ms | **q6 λ0.9 T12 S360, 0.076 ms** | 1.65× | 1.95× |

- **TEC.** Its best is a small *antiferromagnetic* constant coupling (J_v = −4). The ferromagnetic J_v = +4 to +30 recommended by the TEC paper is worse than plain SCA for the 33,000-cut target.
- **Onsager vs TEC.** The theory-scaled, time-varying correction with an end ramp beats every constant coupling by 2–3.5×.
- **Q–λ coupling.** The q–λ trade-off predicted by D (λ_c grows with q) is visible: q = 8 tolerates λ = 1.05.
- **Against STATICA.** The best single engine (0.348 ms) is about 4.3× better than STATICA's published 1.50 ms. Both are modelled times, and ours uses ideal uniform thresholds.
- **Caveat.** Several winners sit on grid edges (λ = 1.05, T_init = 12, S = 360). J2 extends the grid ([PROTOCOL_J2.md](PROTOCOL_J2.md)).

**J2 extension** ([PROTOCOL_J2.md](PROTOCOL_J2.md), 120 pilot configurations, held-out confirmation). Extending beyond J's edges (λ up to 1.4, T_init down to 8, q up to 12) did not improve TTS. Held-out Onsager TTS was 0.358 / 0.133 / 0.109 / 0.085 / 0.070 ms at E = 1 / 4 / 8 / 16 / 32, all within ±15% of J's values (unresolved by the pre-declared 1.2× rule), and λ = 1.05 was again selected. Raising q to 10–12 did not help plain SCA either. J's optimum is therefore essentially interior: **q ≈ 8, λ ≈ 1.05 with end ramp, T_init ≈ 10–12**.
