# Narrowing the gap at K2000's best-known cut (33,337): schedules, cold finish, population annealing

7 October 2026. CPU only (M4 Max, 6 threads). Instance WK2000_1. Success means the final-state cut is ≥ 33,337. All times are model-based: a CPU simulation of the engine algorithm with ideal uniform random numbers, costed with the measured v6.4 cycle model. Nothing here was run on the board.

## Bottom line

- **Best held-out result on the V80 model:** TTS99 = **53.8 ms [44.4, 65.3]** with 12 engines, **5.6× GbSB's 9.61 ms**.
  - Schedule: Onsager-online, S = 16,000, T0 = 12, T1 = 4, q = 8, λ = 1.05.
  - p = 103/8192 = 0.0126 [0.0104, 0.0152].
  - It needs more than the engine's 4,096 schedule-table rows, so it is new (small) hardware.
- **Same cycles on the ASIC:** 13.5 ms [11.1, 16.3] at 1.0 GHz (1.4× GbSB) and 10.8 ms [8.9, 13.1] at 1.25 GHz (1.12×; the interval includes parity).
- **Existing engine (≤ 4,096 rows, table contents only):**
  - Pre-registered best: Onsager-κT S3600 plus a 50-step finish, **104.6 ms [67.0, 163.4]**, 10.9× GbSB.
  - Exploratory: Onsager-online at S = 4096, **81.5 ms [52.8, 125.9]**, 8.5× GbSB.
- **Against the prior probe configuration** (S8000, T1 = 5): held-out p is 22/8192 = 0.0027, i.e. 143.8 ms. The best schedule is therefore 2.7× better.
- **What worked: ending colder.** The paper's end temperature T1 = 5 throws away most ground-state visits. In the screen, only 14–25% of runs that ever visit 33,337 still hold it at the end; with a quenched end point it is 53–75%.
  - Two equivalent ways to fix it: end the geometric schedule at T1 ≤ (q+1)/2, or append a 50–100-step zero-temperature finishing phase (table rows only).
  - On the T1 = 5 reference, the finish raises p 2.45× (22 → 54 of 8,192; paired 32 gained, 0 lost) for +0.8% cycles. This is the same as host greedy descent (22 → 53).
  - The two fixes are not additive: on already-quenched schedules the finish does nothing.
- **What did not work: population annealing.**
  - R = 12: 1.13× in TTS, McNemar p = 0.12.
  - R = 48: 1.02×, p = 0.85.
  - PA copies the replica that found the ground state; it does not find it more often.
- **Longer anneals help only slowly** (same Onsager-online parameters):

| S | 4096 | 8000 | 16,000 | 32,000 |
|---|---:|---:|---:|---:|
| p | 0.0024 | 0.0057 | 0.0126 | 0.0264 |
| TTS99, V80, 12 engines | 81.5 ms | 63.6 ms | 53.8 ms | 49.2 ms |

## 1. Method

- **Pre-registration.** PROTOCOL.md was hashed at 14:19:37Z, before any pilot.
  - Amendment 1 (14:24:59Z, before any PA run) only added gibbs α ∈ {0.125, 0.25}.
  - Its stated rationale was "10.2 of 12 replicas overwritten per event". That count includes slots that systematic resampling merely shifts, so it overstates lineage loss. The amendment only added configurations.
- **Model (`sca_fast.py`).** A numba re-implementation of `abl.run`: integer state, sparse field update, one splitmix64 stream per trial; about 500× faster than the BLAS model.
  - With identical injected uniforms it reproduces `abl.run`'s final spins and flip counts **exactly** for plain, Onsager (with and without ramp) and Onsager-κT.
  - Statistical check at S = 2000: mean cut 33,196.1 vs 33,195.6, SE 5.5.
- **Cost.** cycles/trial = 0.1424·flips + 18.09·S_total + 886. V80 at 250 MHz; ASIC = the same cycles at 1.0 and 1.25 GHz; 12 engines throughout.
  - Round time t_R = bootstrap mean of the maximum of 12 trial cycle counts.
  - Caveat: the model was fitted on S ≤ 1560 and is extrapolated here.
- **TTS99 (Eq. 9).** TTS99 = t_R·ln(0.01)/ln(1−P_R), with P_R = 1−(1−p)^12. Wilson 95% intervals on p, mapped through to TTS.
- **Engine limit found in the RTL.** `tab_mem[0:4095]`, and the host rejects S > 4096.
- **Stages.**
  - (a) 342-configuration screen × 256 runs, ranked within each (family, S) on a proxy: TTS at 33,300 using the Wilson lower bound. The top 2 per group go to a 2,048-run confirm-pilot, which selects by TTS at 33,337 (Wilson lower bound). Selected configurations plus the prior probe configuration then get 8,192 held-out runs with paired starts.
  - (c) 28 finishing options, selected on 4,096 pilot anneals, then applied to the same held-out anneals (paired).
  - (b) 51 PA configurations × 256 populations of 12; the selection then runs on 1,024 held-out populations against paired independent runs (IND); R = 48 on 256 populations.

## 2. (a) Target-specific schedules

**Screen findings**
- For Onsager-online, q = 8 with λ = 1.05 is again best; q = 6 with λ = 1.05 runs away (about 1–4 M flips per run).
- Lowering T1 inside the geometric schedule lowers the mean final cut slightly (by 0–17 cut units) but keeps ground-state visits:

| Family | Final / best-visited at T1 = 5 | At T1 = (q+1)/2 | At T1 = q/2 |
|---|---:|---:|---:|
| Onsager-online | 10/42 | 24/39 | 33/44 |
| Onsager-κT | 11/78 | 50/82 | 47/77 |
| Plain | 4/16 | 8/15 | 12/17 |

**Selections** (all have T1 < 5)
- Onsager-online, overall: S16000, T0 = 12, T1 = 4, q8, λ1.05.
- Onsager-online, existing engine: S2000, T0 = 15, T1 = 3, q6, λ0.9. This pick rested on 3/2048 and failed on held-out.
- Onsager-κT, overall: S8000, T0 = 15, T1 = 3, q6, κ2.
- Onsager-κT, existing engine: S3600, T0 = 15, T1 = 4.5, q8, κ2.
- Plain: S16000, T0 = 30, T1 = 2.5, q4.

**Held-out results** (8,192 runs each, final cut ≥ 33,337, TTS99 in ms with 12 engines)

| Configuration | Rows needed | p [95% CI] | Best-visited k | t_R V80 | V80 | ASIC 1.0 GHz | ASIC 1.25 GHz | V80 / GbSB |
|---|---:|---|---:|---:|---|---|---|---:|
| Onsager-online S16000 T1=4 | 16,000 | 0.0126 [0.0104, 0.0152] | 131 | 1.775 | **53.8 [44.4, 65.3]** | 13.5 [11.1, 16.3] | 10.8 [8.9, 13.1] | 5.6× |
| Onsager-κT S8000 T1=3 | 8,000 | 0.0066 [0.0051, 0.0086] | 75 | 1.344 | 78.0 [59.8, 101.7] | 19.5 [14.9, 25.4] | 15.6 [12.0, 20.4] | 8.1× |
| Onsager-κT S3600 T1=4.5 | 3,600 | 0.0022 [0.0014, 0.0035] | 26 | 0.629 | 109.8 [69.4, 173.5] | 27.4 [17.4, 43.4] | 22.0 [13.9, 34.7] | 11.4× |
| … + finish L50 (q12, T6.5) | 3,650 | 0.0023 [0.0015, 0.0036] | — | 0.633 | 104.6 [67.0, 163.4] | 26.1 [16.7, 40.8] | 20.9 [13.4, 32.7] | 10.9× |
| Onsager-online S2000 T1=3 | 2,000 | 0.0002 (2/8192) | 2 | 0.260 | 409 [112, 1492] | 102 | 82 | 42.6× |
| Plain S16000 T1=2.5 | 16,000 | 0.0042 [0.0030, 0.0058] | 44 | 1.993 | 183.9 [131.6, 257.0] | 46.0 [32.9, 64.3] | 36.8 [26.3, 51.4] | 19.1× |
| Reference: Onsager-online S8000 T1=5 | 8,000 | 0.0027 [0.0018, 0.0041] | 83 | 1.008 | 143.8 [95.0, 217.8] | 36.0 [23.8, 54.5] | 28.8 [19.0, 43.6] | 15.0× |
| Reference + finish L100 (q12, T6.5) | 8,100 | 0.0066 [0.0051, 0.0086] | — | 1.015 | 58.9 [45.2, 76.9] | 14.7 [11.3, 19.2] | 11.8 [9.0, 15.4] | 6.1× |

Single-engine TTS99 on the V80 for the best schedule: 611 ms.

## 3. (c) Cold finishing phase and host greedy descent

- **Mechanism.** With z = q + s·h and 2T_f ≤ q_f + 1, a spin with s·h ≥ +1 never flips. A spin with s·h = −1 flips with probability 1/9, 1/7, 1/5 or 1/3 for (q_f, T_f) = (8, 4.5), (6, 3.5), (4, 2.5), (2, 1.5). So the phase is a stochastic parallel descent.
- **On quenched anneals the finish is a no-op.** The final state is already a 1-flip local optimum: zero extra flips, identical success count, +0.2–1.4% cycles.
- **On the T1 = 5 reference** (held-out, paired):

| | Final-state hits / 8,192 | TTS99 V80 |
|---|---:|---:|
| No finish | 22 | 143.8 ms |
| Selected finish L100 (q12, T6.5), +1,809 cycles (+0.8%) | 54 (32 gained, 0 lost) | 58.9 ms |
| Any descent option, L = 50 | 52–54 | 58.7–61.0 ms |
| Continuing at the anneal's own end, (8, 5.0) | 19–25 | 130–167 ms |
| Host greedy descent | 53 | — |

- **Host greedy cost.** 0.14 ms per final state on one core; recomputing the fields dominates, the descent itself takes about 3 µs. It can overlap with the next round given about 2 host cores, and would be nearly free if the engine exported h.

## 4. (b) Population annealing (R replicas, approximate Boltzmann weights)

**Weights**
- `gibbs-α`: Δlog w = −α·Δβ·E, with systematic resampling at K − 1 fixed stage boundaries.
- Exact population-annealing weights exist only for the logistic kernel. The clipped, Onsager-corrected engine kernel has no closed-form stationary law, so the logistic formula applied here (`pca-1`) is only an approximation. It was worse in the pilot: 0, 3 and 6 of 256 at S = 2000, 4000 and 8000, against IND 1, 10 and 12.

**Copy-cost assumption (new hardware)**
- 886 cycles per resampling event: copy s and s_prev between engines, rebuild h through the existing initial-field path, with an on-chip resampler.
- Sensitivity: 0, and 4 × 886 cycles + 5 µs host round trip, gives 58.7–62.3 ms for the selected PA configuration.

**Pilot (populations reaching 33,337, of 256)**

| S | IND | Gibbs cells: median (range) |
|---:|---:|---|
| 2000 | 1 | 1 (0–8) |
| 4000 | 10 | 7 (2–14) |
| 8000 | 12 | 14 (5–26) |

Selected: S8000, K = 4, α = 0.25 (26/256 in the pilot).

**Held-out**

| | Populations | P_pop [95% CI] | Populations that visited 33,337 | TTS99 V80 | ASIC 1.0 GHz | ASIC 1.25 GHz |
|---|---:|---|---:|---|---|---|
| PA, R = 12 | 1,024 | 0.071 [0.057, 0.089] | 73 | 59.3 [47.2, 74.7] | 14.8 | 11.9 |
| IND, R = 12 | 1,024 | 0.063 [0.049, 0.079] | 78 | 67.1 [52.6, 85.8] | 16.8 | 13.4 |
| PA, R = 48 | 256 | 0.285 [0.233, 0.343] | 73 | 54.0* | 13.5 | 10.8 |
| IND, R = 48 | 256 | 0.277 [0.226, 0.335] | 83 | 55.0* | 13.7 | 11.0 |

\*12 engines time-multiplexed, switching cost ignored. On 48 engines: 13.5 vs 13.7 ms.

- **Paired comparison.** R = 12: 55 both / 18 PA only / 9 IND only, McNemar p = 0.12. R = 48: 59 / 14 / 12, p = 0.85.
- **Mechanism.** In each successful PA population, on average 7.4 replicas hold the ground state; about 7.8 of 12 parents survive each resampling. PA keeps a find rather than making more finds.

## 5. Exploratory additions (not pre-registered)

- **e1, S = 32,000:** 108/4096 = 0.0264 [0.0219, 0.0317]. TTS99 49.2 ms [40.8, 59.4] on the V80; 12.3 ms at 1.0 GHz; 9.8 ms [8.2, 11.9] at 1.25 GHz (1.02× GbSB).
- **e2, 4,096-row tables with a per-row repeat count:**
  - Onsager-online S16000: 85/8192 with repeats vs 99/8192 smooth (65.4 vs 56.1 ms).
  - Onsager-κT S8000: 54 vs 58.
  - Mean cuts equal within 0.5; intervals overlap. A possible loss of about 10–15% is not resolved.
- **e3, existing engine at its maximum length:**
  - Onsager-online S4096 T1 = 4: 0.0024 [0.0016, 0.0038], 81.5 ms (20.4 ms at 1.0 GHz, 16.3 ms at 1.25 GHz).
  - T1 = 5 at S3996 plus a 100-step finish: 5 → 14 of 8,192, 124.3 ms.
- Also exploratory: the descriptive screen aggregates, all non-selected finishing options on held-out runs, the McNemar tests and the distinct-parents diagnostic.

## 6. Pre-registered expectations: outcome

- **H-a1 — pass.** 53.8 ms, 2.7× better than the reference. The protocol had guessed the reference at ~90 ms; it is 143.8 ms because its p is 0.0027, not 0.0039.
- **H-a2 — pass.** Every selected configuration has T1 < 5.
- **H-a3 — pass.** Plain / Onsager-online = 3.4×.
- **H-c1 — partial.** On the T1 = 5 reference: 2.45× for +0.8% cycles, and greedy matches it. On the quenched selections: no gain, as the protocol's conditional clause anticipated.
- **H-b1 — first part passes, second fails.** The R = 12 ratio of 1.13 lies within [0.7, 1.3]. "R = 48 helps more" fails (1.02×).
- **H-overall — pass.** The V80 stays at least 2× slower (5.6×). The ASIC lands at 1.12–1.4× GbSB, inside the expected 0.5–2 range.

## 7. Disclosures

- study.py was edited twice after the protocol hash, without changing any random stream or result:
  1. the α grid line from Amendment 1, made while a1 ran (the a1 code path is unchanged);
  2. a PA diagnostic added before the held-out PA run.
- Every output records the file hashes at run time.
- Stage c also included the reference configuration; this is an addition to the protocol.
- The existing-engine Onsager-online selection was selection noise at small counts (3/2048 in the confirm-pilot, 2/8192 held-out).
- The comparison is lopsided: our numbers are modelled (ideal random numbers, cycle model extrapolated to long schedules) and come from one instance, while GbSB's are measured.

## 8. Recommended wording for the Limitations section

"At K2000's best-known cut (33,337), Anechoic is not competitive with the state of the art: in our CPU model of the engine (8,192 held-out runs), the best schedule we found (Onsager-online, 16,000 steps, ending below the paper's end temperature) reaches a TTS99 of about 54 ms (95% CI 44–65 ms) on twelve V80 engines, 5.6× slower than GbSB's measured 9.6 ms, and it needs a schedule table deeper than the engine's current 4,096 rows (about 80–105 ms, 8–11×, within it). Longer anneals (49 ms at 32,000 steps), a short zero-temperature finishing phase and population annealing over 12–48 replicas did not close the gap. The same cycle counts at 1.0–1.25 GHz in an ASIC would give roughly 11–14 ms, at best parity with GbSB."

## 9. Files

All in `research/optimum_mitigation_20261007/`:
- `PROTOCOL.md`, `PROTOCOL.sha256` (original hash 5502300a…, after Amendment 1 47405a36…)
- `sca_fast.py`, `test_sca_fast.py`, `study.py`, `explore.py`, `summarize.py`, `run_chain.sh`
- `results/`:
  - `a1_screen.json`, `a2_confirm.json`, `c_pilot.json`
  - `holdout_ac.json`, `hold_0.npz` … `hold_5.npz`
  - `b_pilot.json`, `holdout_b.json`, `holdout_b48.json` and their `*_gibbs.npz` / `*_IND.npz`
  - `explore_e1.json`, `explore_e2.json`, `explore_e3.json`
- `logs/`: one log per stage, `test_*.log`, `host_greedy_cost.log`, `chain.log`, and `summarize.log` (all tables, regenerable with `python3 summarize.py`)
