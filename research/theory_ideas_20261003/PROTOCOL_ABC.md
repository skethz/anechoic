# Frozen protocol for ideas A–C

3 October 2026. Written before any A–C sweep result. Earlier quick tests ([REPORT.md](REPORT.md)) were exploratory, used other seeds, and are not reused.

## Common settings

- Instance: WK2000_1 (raw-edge SHA256 verified). J = −w, q = 4, T_fin = 5, geometric T schedule.
- Success: **final** state cut ≥ 33,000.
- Seeds: pilot root 30001, held-out root 30002. Held-out runs never inform selection.
- Kernels:
  - Clipped Eq. 7 with uniform random thresholds, for A and B.
  - Logistic PCA rule `p_flip = 1/(1 + exp(z/T))`, for C. It has the same slope at z = 0, and its stationary law is known.
- Cost model (STATICA-class engine, 300 MHz): `cycles = 0.434 × flips − 0.98 × S + 1989`. For λ ≠ 0 add 4 cycles per step for the n_lin popcount.
- Device TTS with E engines: `t_run × ceil(ceil(R99)/E)`, where R99 = ln 0.01 / ln(1−p). For E = 1, use Eq. 9 (unrounded).
- Selection score: the same TTS evaluated at the pilot's **Wilson 95% lower bound** on p (conservative).
- Report held-out p with Clopper–Pearson 95% intervals.

## A. Onsager-corrected SCA

- Decision field: `h_i − λ·c·s_i(t−1)`, with c = n_lin(t−1)/(2T_{t−1}).
- Pilot grid: λ ∈ {0, 0.25, 0.4, 0.5, 0.6, 0.75} × T_init ∈ {15, 20, 30, 40} × S ∈ {360, 560, 960, 1560}. That is 96 configurations, 256 runs each.
- Selection: for each E ∈ {1, 4, 8, 16, 32}, pick the lowest-score configuration among λ = 0 (baseline) and among λ > 0.
- Held-out: 1,024 runs for each distinct selected configuration.
- Report held-out TTS for both families and their ratio.

## B. IAMP(t*) + SCA finish

- IAMP exactly as in `research/iamp_evaluation_20261003`: K ∈ {64, 128}, t* ∈ {0.8, 0.9, 1.0}, start from sign(m).
- Finish grid: SCA with T_a ∈ {6, 8, 12}, S_f ∈ {100, 200, 400}, λ_f ∈ {0, 0.5}. That is 108 configurations, 256 runs each.
- Hypothetical cost model: each IAMP step costs `0.434·N + 100` cycles, and IAMP's matrix-vector products are multi-bit multiply-accumulates. This favours IAMP.
- Selection: for E = 1 and E = 16. Held-out: 1,024 runs per selected configuration.

## C. Population annealing with exact PCA weights

- Kernel: logistic PCA, λ = 0, a population of E replicas with a shared schedule.
- Incremental log weight per step, using `π_PCA,T ∝ exp(−H/T) · f_T`, `log f_T = Σ_i softplus(−(q + s_i h_i)/T)`:
  `Δlog w = −H·(1/T_{t+1} − 1/T_t) + log f_{T_{t+1}} − log f_{T_t}`
- Systematic resampling whenever ESS < E/2.
- Variants:
  - **IND:** independent replicas, no resampling.
  - **PA-exact:** the weights above.
  - **PA-Gibbs:** drops the log f terms, to test whether the identity matters.
- Population success: any replica's final cut ≥ 33,000. Estimated from 64 populations per configuration in the pilot, and 256 held-out populations for the selected configurations.
- Pilot grid: E ∈ {16, 32} × S ∈ {200, 360, 560} × T_init ∈ {20, 30} × {IND, PA-exact, PA-Gibbs}.
- Device TTS: one population is one wave of E engines, so TTS = `t_run × ceil(ln 0.01 / ln(1 − P_pop))`. The resampling copy cost is reported separately and not charged (an optimistic assumption).
- Selection: for each E, pick the lowest-score configuration for IND and for the best PA variant.

All runs are on the laptop CPU. GPU timing is a separate later step.
