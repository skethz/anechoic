# Dynamical mean-field theory of synchronous stochastic annealers (SCA family)

Status, 4 October 2026: the derivation and a validated numerical solver exist. The optimal-control step (deriving schedules) has not started. Code: [dmft.py](dmft.py).

## 1. Model

**Couplings.** Symmetric, with zero diagonal. K2000 has J_ij = ±1, N = 2000. Rescaled units:

A = J/√N,  T̃ = T/√N,  q̃ = q/√N,  β̃ = β/√N,  u_i(t) = Σ_j A_ij s_j(t) = h_i(t)/√N.

**Synchronous update.** All spins update at once:

g_i(t) = u_i(t) + q̃_t s_i(t) − β̃_t s_i(t−1),  P(s_i(t+1) = +1) = (1 + f_t(g_i(t)))/2,

- Clipped rule (STATICA Eq. 7, rewritten): f_t(g) = clip(g/(2T̃_t), −1, 1).
- Logistic rule: f_t(g) = tanh(g/(2T̃_t)).

**Special cases** of the lag-1 term β:

| Case | β |
|---|---|
| Plain SCA | 0 |
| TEC | −J_v (constant) |
| Onsager SCA | λ_t · c(t−1), with c = n_lin/(2T) |

The stay-probability form STATICA uses, clamp((s·h + q)/(4T) + ½), is the same rule written for s(t+1) = s(t).

## 2. Effective single-spin process (N → ∞)

**Generating functional.**
- Write the path probability with u_i(t) enforced by auxiliary fields û_i(t).
- Average over symmetric couplings with variance 1/N. For ±1 entries the higher cumulants vanish as N → ∞ (universality).
- The quadratic term is N Σ_{t,t′}[Q̂(t,t′)C(t,t′) + K(t,t′)K(t′,t)]/2, with order parameters C = N⁻¹Σ s s, Q̂ = N⁻¹Σ û û and K = N⁻¹Σ û s.
- The saddle point gives Q̂ = 0 and iK(t,t′) = −G(t′,t), which leaves one effective spin:

  **u(t) = η(t) + Σ_{t′<t} G(t,t′) s(t′),  η Gaussian with E[η(t)η(t′)] = C(t,t′),**

  with self-consistency C(t,t′) = E[s(t)s(t′)] and G(t,t′) = ∂E[s(t)]/∂θ(t′). Here θ(t′) is a field perturbation added to g(t′).

This is the discrete-time, parallel-update counterpart of the Sompolinsky–Zippelius and Eissfeller–Opper equations for SK dynamics; Coolen 2001 treats parallel dynamics. Our numerical validation (section 4) is an independent check.

**The retarded self-interaction Σ G(t,t′)s(t′) is the "echo".** It is the field each spin receives from its own past through the rest of the network. It is absent for asymmetric couplings (Mézard & Sakellariou 2011).

## 3. Exact consequences

1. **The hardware's Onsager coefficient is the lag-1 response.**
   - Under the clipped rule, G(t,t−1) = E[f′(g(t−1))] = P_lin(t−1)/(2T̃_{t−1}). In raw units this is exactly c(t−1) = n_lin/(2T).
   - The online n_lin popcount is therefore a finite-N estimate of a self-averaging deterministic quantity. Onsager SCA with λ = 1 cancels the lag-1 echo term, leaving (1 − λ)G(t,t−1)s(t−1) + Σ_{k≥2} G(t,t−k)s(t−k).
2. **Flip-backs come from the lag-1 echo.**
   - The term G(t,t−1)s(t−1) enters the decision that produces s(t+1) and biases it toward s(t−1), which is exactly a flip-back.
   - This explains the measured drop in immediate flip-backs from 67% to 54% at λ = 0.7.
3. **A constant coupling cannot track the echo.**
   - G(t,t−1) ∝ P_lin/T̃ changes over the anneal, so the best constant TEC coupling must depend on temperature. TEC reports this empirically: "the optimal Jv shifts monotonically ... with increasing temperature".
   - The theory predicts the sign at low T (J_v = −β < 0) and gives the coefficient.
4. **Energy and cost are observables of the effective process.**
   - E[s u](t) = 2 Σ_{t′<t} G(t,t′)C(t,t′) (Gaussian integration by parts).
   - cut = (Σ_{i<j} w_ij + ½ N^{3/2} E[s u])/2.
   - Flips per step = N·P(s(t+1) ≠ s(t)), and n_lin = N·P_lin.
   - The cut ≥ 33,000 target corresponds to E[s u] ≈ 1.499. The Parisi ground-state value of E[s u] is 2 × 0.7633 = 1.527.
5. **The coefficients can be precomputed.**
   - Because c(t) is deterministic as N → ∞, β_t can come from a table, which would remove the hardware popcount.
   - This is testable: fixed table vs online n_lin.

## 4. Numerical solution and validation

**Method.** Sample M paths of the effective process with C and G held fixed. Re-estimate C empirically, and estimate G with the exact score identity, Rao-Blackwellized:
- Lag 1 uses the exact E[f′].
- Lags ≥ 2 use E[(f(g(t−1)) − mean)·score(t′)].

**Solvers tried, in order.** All are in [dmft.py](dmft.py); logs are in this folder.

| Solver | Result |
|---|---|
| Single forward pass (`solve`) | Unstable: final E[s u] spread 1.32–1.54 across seeds at M = 160k. Noise in many G entries adds up coherently for persistent spins and acts like a random error in q. |
| Gaussian integration by parts for G | Fails at low T: C is ill-conditioned. |
| Global damped iteration (`solve_iter`) | Good at S = 200 with 120 iterations. At S = 560 it converges very slowly, and intermediate iterates exceed the ground-state bound. |
| **Causal block iteration (`solve_causal`)** | Uses the fact that rows up to t don't depend on later times: windows of 40 steps, 1/n running averages, and a PSD projection of C. **Accuracy is set by passes per window:** with 64 passes, S = 200 matches within 0.005. |

**Validation:** plain SCA, S = 200, T 30→5, q = 4, metric E[s u](t).

| t | 25 | 50 | 100 | 150 | 200 |
|---|---|---|---|---|---|
| Simulation, random instances, N = 1000 | 0.793 | 1.159 | 1.418 | 1.468 | 1.476 |
| Simulation, N = 2000 | 0.797 | 1.168 | 1.422 | 1.470 | 1.478 |
| Simulation, N = 4000 (sd across instances ≤ 0.0002 late) | 0.793 | 1.168 | 1.424 | 1.470 | 1.478 |
| Simulation, K2000 | 0.782 | 1.159 | 1.418 | 1.467 | 1.475 |
| DMFT, causal, 64 passes per window | 0.802 | 1.175 | 1.421 | 1.475 | 1.473 |

- **Finite-N effects are negligible** on these time scales, so every discrepancy is solver error.
- **Flip cost** (flips per run) is predicted within 1–5% at all λ, even by the cruder solvers.
- **The cut target is tight for the solver.** cut ≥ 33,000 corresponds to E[s u] = 1.499, and 0.005 in E[s u] is about 110 cuts. DMFT can rank schedules and locate thresholds, but predicting success probabilities would also need the finite-N spread of the final cut (about 70–190 cuts across runs).

## 5. Onsager threshold test

Short schedule: S = 560, T 30→5, q = 4. Files: `lambda_scan.log`, `conv_causal64.log`, `finite_n_lam.log`.

- **What the theory gets right.**
  - **The failure at λ = 1:** E[s u] = 0.40 from DMFT vs 0.43 simulated; the system never orders.
  - **Flip cost** at every λ (global iteration): within 1–5%, from 86k flips per run at λ = 0 to 464k at λ = 1.
- **Trajectory accuracy**, causal solver with 64 passes per window:

| E[s u] at t | 140 | 280 | 420 | 560 | Flips per run |
|---|---|---|---|---|---|
| λ = 0, DMFT (2 seeds) | 1.297 / 1.297 | 1.460 / 1.460 | 1.497 / 1.488 | 1.514 / 1.526 | 87.5k / 88.0k |
| λ = 0, simulation (K2000) | 1.288 | 1.453 | 1.486 | 1.492 | 86.8k |
| λ = 0.7, DMFT (2 seeds) | 0.869 / 0.879 | 1.345 / 1.352 | 1.490 / 1.443 | 1.449 / 1.487 | 186.5k / 185.3k |
| λ = 0.7, simulation (K2000) | 0.819 | 1.317 | 1.474 | 1.500 | 173.6k |
| λ = 0.7, simulation, random instances N = 1k / 2k / 4k / 8k | 0.836 / 0.836 / 0.849 / 0.844 | 1.322 / 1.327 / 1.330 / 1.332 | 1.471 / 1.479 / 1.481 / 1.484 | 1.501 / 1.504 / 1.506 / 1.508 | 173k / 171k / 170k / 169k |

- **Plain SCA is accurate within about 0.01 until the final cold stretch,** which is biased upward by 0.02–0.03.
- **The Onsager mode shows about 0.03–0.05 error and 7% excess flips.** Finite-N effects are larger than for λ = 0 but explain only part of it. Most of the gap is solver error.
- **This is not yet accurate enough** to derive schedules near the cut ≥ 33,000 target, where 0.005 is about 110 cuts.
- **Likely sources of the solver error:**
  - Score-estimator noise at low T, where f′ = 1/(2T̃) is large and the clipped score is heavy-tailed.
  - The PSD projection of a nearly singular C.
  - Slow convergence of the fixed point in the cold, aging regime.


## 6. Next steps

1. **Optimal control on the DMFT.** Choose T_t, q_t and β_t (and possibly lag-2 terms) to maximize the final E[s u] at a fixed flip and step budget, using the measured hardware cost model. Then pre-register the predicted schedule and its predicted ranking.
2. **Precomputed vs online coefficient** (a hardware simplification).
3. **Solver precision:** larger M, variance reduction, and checking the late-time regime at S = 960–1,560.

## 7. Track B: multi-lag echo correction (exploratory, negative)

Files: `multilag.py`, `multilag_measure.log`, `multilag_test.py`, `multilag_test.log`, `multilag_test.json`.

**Measurement.** Echo coefficients ĉ_k(t) = N·∂E[s(t)]/∂h(t−k) were measured exactly in the N = 2,000 system by noise injection (ε = 1). An exact control variate makes this feasible: s(t−k) is independent of the injected r(t−k).

| | Lag 1 | Lags 3 and 5 | Lags 2 and 4 |
|---|---|---|---|
| Plain SCA (P1) | 16 → 1.8 over the anneal | Notable: 4.5 and 2.2 early, roughly 0.3× lag 1 (period-2 structure) | Small |

- **Onsager SCA (O1):** lag 1 falls from 26 to 0.4. Lags 2–6 persist (0.3–0.9 each late), and late in the run their sum exceeds lag 1 (2.6 vs 1.45 at t = 720).

**Test.** O1 plus `−λ_x·ramp(t)·ĉ_k(t)·s(t−k)` for k = 2..K, with the measured tables smoothed over 21 steps. Pilot: 512 paired runs per configuration.
- With λ₁ = 0.9: the extra terms restore success from 0.13 to 0.90 (K = 6, λ_x = 0.5).
- With λ₁ = 1.05: the extra terms over-correct (flips 340k–470k) and success drops to 0.03–0.68.
- **Confirmation (2,048 fresh paired runs):** O1 TTS99 0.59 ms (p = 0.942) vs best multi-lag 0.77 ms (p = 0.877). **No gain.**

**Explanation.**
- The residual long-lag echo acts mostly on persistent spins, where s(t−k) ≈ s(t).
- So Σ_k ĉ_k s(t−k) ≈ (Σ_k ĉ_k)·s(t), which is a negative inertia term duplicating the existing q.
- The multi-lag terms therefore trade off against λ₁ and q rather than adding a new degree of freedom. This is consistent with the q–λ coupling found in theory D.

## 8. Precomputed coefficient table vs online popcount (prediction 5: not confirmed)

File: `table_vs_online.log`. O1 configuration, 2,048 paired runs per variant; the table is the mean online c(t) from a separate 256-run pilot.

| Variant | Success | Mean cut | Flips per run |
|---|---|---|---|
| Online popcount | 0.927 | 33,132 | 133,261 |
| Precomputed table | 0.901 | 33,120 | 129,290 |

- **The table is slightly worse.**
- **The reason is collective, run-level fluctuation of n_lin at N = 2,000 in Onsager mode.** Across single runs, the coefficient's spread (sd/mean) is 7% at t = 100, 39% at t = 300, 22% at t = 600 and 66% at t = 900.
- **The online coefficient feeds each run's own fluctuation back into its dynamics, and that helps.** Keep the hardware popcount.
- **This also suggests why Onsager-mode DMFT (N → ∞) differs more from N = 2,000 than plain mode does:** near the correction's instability, finite-N fluctuations are large.

## 9. GPU solver (Track A): accuracy reached

`dmft_torch.py` runs on gpu-host GPU1 (UUID-pinned, everything under `/scratch/USER/snowball_dmft_20261004`), with M = 1–2 million paths per pass.

**Key fix (`solve_causal2`).** After a short warm-up, each window replaces the whole leading block of C and G with the common-weight mean of that window's estimates.
- C then stays an exact average of Gram matrices, so it is consistent and positive semidefinite.
- Per-row averaging had broken this. It distorted the tiny low-temperature innovations of η and biased the cold end upward (S = 200 final: +0.026 → +0.006).

| Case (E[s u] at the listed t) | DMFT, GPU v2 | Simulation |
|---|---|---|
| Plain, S = 200, t = 25/50/100/150/200 | 0.799 / 1.175 / 1.429 / 1.472 / 1.484 | 0.793 / 1.168 / 1.424 / 1.470 / 1.478 (N = 4,000 random instances) |
| Plain, S = 560, t = 140/280/420/560 | 1.303 / 1.463 / 1.486 / 1.495; flips 85.8k | 1.288 / 1.453 / 1.486 / 1.492; flips 86.8k (K2000) |
| λ = 0.7, S = 560, mean of 4 seeds | 0.866 / 1.340 / 1.487 / 1.495 (final spread ±0.012); flips 170–173k | 0.844 / 1.332 / 1.484 / 1.508 (N = 8,000); flips about 171k |

## 10. Assessment of the theory program so far

- **The DMFT works as an explanatory theory.**
  - It reproduces the dynamics within about 0.01 in E[s u] and flip cost within a few percent.
  - It identifies the Onsager coefficient exactly as the lag-1 response.
  - It explains flip-backs, the λ = 1 failure, and TEC's temperature-dependent optimum.
  - It explains why multi-lag correction duplicates q (section 7).
- **It is not a better schedule optimizer for N = 2,000 hardware.**
  - About 0.01 in E[s u] is roughly 200 cuts, coarser than the differences between good schedules (50–150 cuts).
  - In Onsager mode, the N = 2,000 system itself deviates from N → ∞ by 0.01–0.02 (sections 5 and 8). The online coefficient's finite-N fluctuations even help (section 8).
  - Direct finite-N simulation (seconds per schedule) remains the better tool for numerical tuning.
- **Both algorithmic consequences tested so far are negative:**
  - Multi-lag correction (section 7).
  - A precomputed coefficient table (section 8).
- **The surviving algorithm is one-step, online Onsager SCA.** The theory now explains why it is the right minimal correction.
