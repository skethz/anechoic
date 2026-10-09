# Narrowing the gap at K2000's best-known cut (33,337): pre-registered protocol

Written 7 October 2026, before any pilot or held-out run of this study. CPU only (Apple M4 Max, 6 threads).
Instance WK2000_1 (N = 2000, ±1 weights), loaded by `abl.m.load()` (SHA-256-checked .rud file).

## 0. Model, code and what was run before this protocol

- **Model.** `sca_fast.py`: the clipped synchronous SCA of `research/ablation_20261005/abl.py::run` (families `plain`,
  `onsager` = Onsager-online with coefficient λ·n_lin/2T, `tecT` = Onsager-κT), with integer state, a sparse field
  update and one splitmix64 stream per trial. Optional final temperature T1 (abl.py fixes T1 = 5) and an optional
  finishing phase appended to the schedule (§3).
- **Development runs (disclosed, not used for selection or estimates).** `test_sca_fast.py`:
  (i) exact replay against `abl.run` with injected uniforms: identical final spins and flips for plain, Onsager (ramp and
  no ramp) and Onsager-κT, 6 trials each (seed [777, i]); (ii) a speed check (S = 8000, seeds 1–2; 24 + 96 trials,
  max cut 33,336); (iii) a statistical check at S = 2000 against `abl.run` (seeds 4242/4243, 256 trials each: mean cut
  33,196.1 vs 33,195.6, SE 5.5). Logs in `logs/test_*.log`. No other runs of this study precede this protocol. The prior
  exploratory probe (`research/optimum_target_20261006`, 1/256 at 33,337 for Onsager q8 λ1.05 T0=12 S8000) is known.
- **Seeds.** `np.random.SeedSequence([20261007, stage, index])`: stage 1 = (a) screen, 2 = (a) confirm-pilot,
  3 = (c) pilot, 4 = (b) pilot, 9 = held-out (a)+(c), 10 = held-out (b) R = 12, 11 = held-out (b) R = 48.
  Held-out runs of different configurations share initial spins and RNG seeds (paired starts).

## 1. Success, metrics, statistics

- **Success.** FINAL-state cut ≥ 33,337 (the engine returns the final state). Best-visited cut (max over all states
  of the trajectory) is reported separately and is never used for TTS. Also reported: final cut ≥ 33,300 / 33,320,
  mean final cut, and the cut after a host-side steepest single-flip descent from the final state (reference only).
- **Intervals.** Wilson 95% intervals on every probability; TTS intervals map the p interval through the TTS formula.
- **Cost model (v6.4 engine, measured).** cycles/trial = 0.1424·flips + 18.09·S_total + 886, S_total = anneal +
  finishing steps. V80: 250 MHz. ASIC: the same cycles at 1.0 GHz and at 1.25 GHz. Twelve engines in all cases.
- **Round time.** t_R = E[max of 12 trial cycles] / f, estimated by bootstrap (4,000 draws of 12) from the per-trial
  cycle distribution of the same runs.
- **TTS (Eq. 9).** TTS99 = t_R · ln(0.01) / ln(1 − P_R), P_R = 1 − (1 − p)^12; P_R ≥ 0.99 gives one round.
  Single-engine TTS99 = t_trial · ln(0.01)/ln(1 − p) is also reported. Ratio to GbSB = TTS99 / 9.61 ms (Goto, Hidaka,
  Tatsumura, Phys. Rev. Applied 25, 044011, 2026; Agilex 7 FPGA).
- **Engine limit.** The v6 engine's schedule table has 4,096 rows (`tab_mem[0:4095]`, host rejects S > 4096).
  Configurations with S_total ≤ 4096 run on the existing engine with new table contents only; S_total > 4096 needs a
  deeper table or a per-row repeat count (new hardware). Results are reported in both classes.

## 2. Family (a): target-specific schedules

**Screen grid (stage 1, 256 runs per configuration).** T(t) = T0·(T1/T0)^(t/(S−1)); end ramp of the correction over the
last 30% for both corrected families. T1 ∈ {5, (q+1)/2, q/2}: at T1 ≤ (q+1)/2 a spin with s·h = +1 can no longer flip
(z = q+1 ≥ 2T), so these end points quench; at T1 = q/2 a spin with s·h = −1 still flips with probability 1/(2q).
- Onsager-online: S ∈ {2000, 3600, 8000, 16000} × T0 ∈ {10, 12, 15} × q ∈ {6, 8} × λ ∈ {0.9, 1.05} × T1 (3) = 144.
- Onsager-κT: S (same 4) × T0 ∈ {12, 15, 20} × q ∈ {6, 8} × κ ∈ {1.5, 2.0} × T1 (3) = 144.
- Plain SCA (reference): S ∈ {3600, 8000, 16000} × T0 ∈ {30, 40} × q ∈ {4, 6, 8} × T1 (3) = 54.

**Stage-1 proxy.** Within each (family, S): rank by 12-engine V80 TTS99 at the 33,300 threshold with p(final ≥ 33,300)
at its Wilson lower bound; ties and all-zero cases by mean final cut. The top 2 per (family, S) go to stage 2.

**Stage 2 (2,048 fresh pilot runs per candidate).** Rank by 12-engine V80 TTS99 at 33,337 with p at its Wilson lower
bound. Per family select (i) the best overall and (ii) the best with S_total ≤ 4096 (existing engine). If a family has
no candidate with k ≥ 1 at 33,337, rank that family at 33,320, then 33,300.

**Held-out (8,192 runs each, stage 9 seeds, paired).** The selected configurations (up to 2 per corrected family, 1 for
plain) plus the fixed reference = the prior probe configuration (Onsager q8 λ1.05 T0=12 S8000 T1=5, ramp).

## 3. Family (c): cold finishing phase and host greedy descent

**Finishing phase.** L extra steps at constant (q_f, T_f), correction off (λ = κ = 0), appended to the anneal: only
extra rows of the existing per-step tables (fourT, q, kcorr = 0, kconst = 0).
Grid: L ∈ {50, 100, 200, 400} × (q_f, T_f) ∈ {(8, 4.5), (6, 3.5), (4, 2.5), (2, 1.5), (1, 1.0), (12, 6.5),
(q, T1) of the anneal} = 28 options; for existing-engine configurations only options with S + L ≤ 4096.
In (q_f, T_f) = (8, 4.5), (6, 3.5), (4, 2.5), (2, 1.5), (12, 6.5) a spin with s·h ≥ +1 never flips and a spin with
s·h = −1 flips with probability 1/9, 1/7, 1/5, 1/3, 1/13 (a stochastic parallel descent); (1, 1.0) flips such spins
with probability 1/2 and those with s·h ≤ −3 always.
**Pilot (stage 3, 4,096 fresh anneals per selected (a) configuration).** Each option is applied to the same final
anneal states. Select per configuration the option with the lowest 12-engine V80 TTS99 at 33,337 (Wilson lower bound),
charging the finishing steps and flips. **Held-out:** the selected option is applied to the same 8,192 held-out (a)
anneals (paired with and without finishing); all other options are also run there and labelled exploratory.
**Host greedy descent.** Steepest single-flip descent to a 1-flip local optimum, from the final state (with and without
finishing); host cost reported separately, not in engine cycles.

## 4. Family (b): population annealing (PA) across R replicas

- **Kernel and schedule.** Onsager-online with (T0, q, λ, T1, ramp) of the (a)-selected overall Onsager-online
  configuration, S ∈ {2000, 4000, 8000}.
- **Resampling.** K stages; at t_k = round(k·S/K), k = 1..K−1, every population computes the replicas' energies and
  resamples systematically with weights w ∝ exp(Δlog w).
  - `gibbs-α`: Δlog w = −α·(β(t_{k+1}) − β(t_k))·E, β = 1/T (Boltzmann reweighting; approximate, because the clipped,
    corrected parallel kernel has no closed-form stationary law; α is a selection-strength factor).
  - `pca-1`: Δlog w = −(β(t_{k+1}) − β(t_k))·E + Σ_i [ℓ(T(t_{k+1})) − ℓ(T(t_k))], ℓ(T) = log(1 + exp(−(q + s_i h_i)/T)):
    the exact PA weight of the LOGISTIC parallel kernel (π_PCA ∝ π_T / P_hold, verified in
    `theory_ideas_20261003/sweeps_abc.py::selftest`), used here only as an approximation for the clipped kernel.
  - A copied replica takes s, s_prev, h and n_lin of its parent and keeps its own RNG stream.
- **Grid.** S (3) × K ∈ {4, 8, 16} × α ∈ {0.5, 1, 2} (gibbs) + `pca-1` at K = 8 for each S = 30 PA configurations,
  plus independent runs (IND, no resampling) at each S.
- **Cost.** Replicas run on 12 engines in lockstep: t_pop = 886 + Σ_k max_r(0.1424·flips_{r,k} + 18.09·L_k) +
  (K−1)·C_rs. The current engine has no state-copy feature; assumed new hardware: C_rs = 886 cycles per resampling
  event (copy of s and s_prev between engines and recomputation of h by the existing initial-field path, all copies in
  parallel), with an on-chip resampler. Sensitivity: C_rs = 0, and 4·886 cycles + 5 µs of host round-trip per event.
  IND: t_R = 886 + max_r(0.1424·flips_r + 18.09·S), same lockstep convention.
- **Metric.** P_pop = P(at least one replica's FINAL cut ≥ 33,337); TTS99 = t_pop·ln(0.01)/ln(1 − P_pop).
- **Pilot (stage 4).** 256 populations of R = 12 per configuration; select the PA configuration with the lowest V80
  TTS99 at 33,337 at the Wilson lower bound of P_pop (fallback 33,320 if no PA configuration has a success).
- **Held-out (stage 10).** 1,024 populations (12,288 runs) of the selected PA configuration and of IND at the same S on
  the same initial replicas and RNG seeds (paired).
- **R = 48 (stage 11, if time allows).** The R = 12 selection without re-tuning, 256 populations, with paired IND-48.
  Costed as 48 replicas time-multiplexed on 12 engines (4× the lockstep stage time, switching cost ignored,
  i.e. optimistic) and, alternatively, on 48 engines.

## 5. Pre-registered expectations (honest guesses, not commitments)

- **H-a1.** The selected Onsager-online schedule has held-out p(33,337) at least 1.5× the probe reference per unit
  cost, i.e. 12-engine V80 TTS99 at most ~60 ms (the reference is about 90 ms under the v6 model if p ≈ 0.004).
- **H-a2.** A quenching end point (T1 < 5) is selected for at least one corrected family.
- **H-a3.** Onsager-online or Onsager-κT beats plain SCA at 33,337 by at least 1.5× in TTS99.
- **H-c1.** The selected finishing phase raises held-out p(33,337) by at least 1.2× on the same anneals at ≤ 5% added
  cycles; host greedy descent gives a gain of similar size. If (a) already selects T1 < 5, the gain is smaller.
- **H-b1.** PA at R = 12 changes TTS99 relative to paired IND by a factor between 0.7 and 1.3; a ≥ 1.5× improvement
  would be a surprise. R = 48 helps more than R = 12.
- **H-overall.** At 250 MHz with 12 engines the best combination remains at least 2× slower than GbSB (9.61 ms); the
  ASIC (same cycles at 1.0–1.25 GHz) lands within a factor of 0.5–2 of GbSB.

## 6. Reporting

Pre-registered: the selections and held-out estimates above, the paired with/without-finishing comparison of the
selected option, and PA vs paired IND. Everything else (other finishing options on held-out runs, S > 16000, piecewise
constant tables, R = 48 costing variants) is labelled exploratory. All numbers are model-based (CPU simulation of the
engine algorithm with ideal uniform random numbers and the measured cycle model), not board measurements.

## Amendment 1 (7 October 2026, 16:24 CEST; before any stage-b run)

**Why.** In a development smoke test of the PA code (S = 300, 8 populations of 12, dev seeds [999, 3–4], not study data),
gibbs-α = 1 at K = 4 overwrote on average 10.2 of 12 replicas per resampling event, i.e. near-greedy selection. The
log-weight spread is roughly α·Δβ·sd(E) with Δβ per stage ≈ (1/T1 − 1/T0)/K, so α ∈ {0.5, 1, 2} covers only strong
selection. **Change.** Stage b's gibbs grid becomes α ∈ {0.125, 0.25, 0.5, 1, 2} × K ∈ {4, 8, 16} (45 gibbs + 3 pca + 3
IND = 51 configurations, 256 populations each); everything else unchanged. Written while stage a1 was running (rows
0–39 of a1 seen; they carry no information on PA).
