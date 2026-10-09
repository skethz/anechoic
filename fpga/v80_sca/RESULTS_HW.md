# Physical V80 results: Onsager-SCA engine

4 October 2026, fpga-host. The board run followed [PROTOCOL_HW.md](PROTOCOL_HW.md), which was frozen before any board data; its SHA256 is in `results/protocol_hw_frozen.sha256` on fpga-host.
- Raw per-trial records and final states: [results/hw_board_300mhz](results/hw_board_300mhz).
- Analysis: [scripts/analyze_hw.py](scripts/analyze_hw.py) → `hw_summary.json`.

## Bottom line

- **The physical engine runs the reference algorithm exactly.** All 24 bring-up trials matched the fixed-point reference bit for bit: every spin, the cut, the flip count, and the n_lin value at every step. All 11,456 cohort and latency trials passed host rescoring.
- **Measured single-engine TTS99 of Onsager SCA (O1) is 0.634 ms** (p = 0.930; 0.366 ms per run at 300 MHz).
  - That is **4.19× lower than the best plain SCA configuration** (bootstrap 95% CI 3.74–4.68×).
  - It is **3.72× lower than the best TEC configuration** (3.33–4.20×).
  - Both baselines are post-hoc minima over the measured configurations, which favours them.
- **The correction costs no measurable time.** The fitted per-run cycle term for the correction is −533 cycles, zero within noise. Its benefit comes entirely from needing fewer flips and steps for the same success.
- **This engine is about 1.7–2.6× slower per run than the STATICA-class cycle model used in J.** The cause is per-step overhead (about 40 cycles per step), not the flip path (0.52 cycles per flip). The Onsager configurations are the least affected because they use fewer steps, so J's ranking holds on hardware.

## Image

| Item | Value |
|---|---|
| Kernel clock | 300 MHz; routed WNS 0.000 ns with 0 failing endpoints out of 567,389; hold +0.010 ns |
| Logic UUID | `e818b91b2fb4d907f2b0d59397626a3c` (`build/board_300mhz`, PDI SHA256 `7dc5ae36…`) |
| Whole image, including shell | 200k LUT (7.8%), 227k FF (4.4%), 140.5 BRAM tiles (3.8%), 4 URAM, 283 DSP (2.6%) |
| Verification chain | Native testbench, HLS C simulation and RTL co-simulation: bit-exact. Board: 24/24 exact (phase 1). |

A 250 MHz image was also built and is timing-clean (WNS +0.030 ns). Per the protocol, it was not run.

## Phase 2 cohorts (1,024 paired trials each, 64 per launch)

| Key | Configuration | p_hw [95% CI] | p_J (float model) | Reproduced | t_run (ms) | TTS99 (ms) | Measured/model cycles |
|---|---|---|---|---|---:|---:|---:|
| P1 | Plain q4, T 30→5, S 1560 | 0.613 [0.583, 0.643] | 0.675 | **No** | 0.584 | 2.83 | 1.90 |
| P2 | Plain q6, T30, S1560 | 0.580 [0.549, 0.611] | 0.567 | Yes | 0.501 | 2.66 | 2.10 |
| P3 | Plain q8, T30, S1560 | 0.438 [0.408, 0.469] | 0.443 | Yes | 0.443 | 3.53 | 2.33 |
| P4 | Plain q8, T30, S960 | 0.180 [0.157, 0.205] | 0.205 | Yes | 0.282 | 6.55 | 2.24 |
| T1 | TEC J_v=−4, q8, T30, S1560 | 0.631 [0.600, 0.660] | 0.675 | Yes | 0.511 | 2.36 | 2.07 |
| T2 | TEC J_v=−4, q8, T30, S960 | 0.392 [0.362, 0.422] | 0.404 | Yes | 0.323 | 2.99 | 2.02 |
| T3 | TEC J_v=+4, q8, T30, S1560 | 0.200 [0.176, 0.226] | 0.198 | Yes | 0.392 | 8.07 | 2.64 |
| **O1** | **Onsager λ1.05 + ramp, q8, T 12→5, S960** | **0.930 [0.912, 0.945]** | 0.937 | Yes | **0.366** | **0.634** | 1.74 |
| O2 | Onsager λ0.9 + ramp, q6, T15, S560 | 0.744 [0.716, 0.771] | 0.737 | Yes | 0.232 | 0.784 | 1.67 |
| O3 | Onsager λ0.9 + ramp, q6, T12, S560 | 0.506 [0.475, 0.537] | 0.524 | Yes | 0.182 | 1.19 | 1.88 |
| O4 | Onsager λ0.9, q6, T12, S360 | 0.335 [0.306, 0.365] | 0.315 | Yes | 0.132 | 1.49 | 1.74 |

- **t_run** is device time per run, with table loading amortized over 64 trials.
- **TTS99** uses Eq. 9.
- **Mean flips per run** match the float model within 1% for every configuration (for example, O1: 134,111 vs 132,887).

### Replicated-engine projections

These use E identical engines at the measured per-run time. They are projections, not measurements; replication is not built.

| Engines | Best plain | Best TEC | Best Onsager | Onsager vs plain | Onsager vs TEC |
|---:|---|---|---|---:|---:|
| 1 (whole runs) | P1, 2.92 ms | T1, 2.56 ms | O1, 0.731 ms | 3.99× | 3.50× |
| 4 | P3, 0.885 ms | T2, 0.969 ms | O2, 0.232 ms | 3.82× | 4.18× |
| 8 | P3, 0.443 ms | T1, 0.511 ms | O3, 0.182 ms | 2.43× | 2.80× |
| 16 | P3, 0.443 ms | T2, 0.323 ms | O4, 0.132 ms | 3.34× | 2.44× |
| 32 | P4, 0.282 ms | T2, 0.323 ms | O4, 0.132 ms | 2.13× | 2.44× |

## Measured timing model

| Quantity | Measured | Source |
|---|---|---|
| Cycles per run, couplings resident | **0.522·flips + 40.2·S + 1,884** | 165 launches, max error 0.44% |
| Correction term | −533 cycles per run (none) | Same fit |
| Table loading, per launch | 35.3 cycles per step | Single vs batched launches |
| Coupling load (first launch only) | 61,000–68,000 cycles | Single-launch cohorts |
| Single-run device latency, including table loading | O1 0.478 ms; P1 0.771 ms; T1 0.699 ms (means of 64) | Phase 2 single-trial launches |

- The flip path achieves its design rate: two flip extractors give about 2 flips per cycle, minus imbalance.
- Per-step overhead is the cost: the 16-cycle lane pass, pipeline fill and drain, and the extractor start.
  - It is 35% of O1's run time and 36% of P1's.
- Table loading appears to use single-beat HBM reads, since about 35 cycles per step is consistent with bounded outstanding reads. It affects only the first run of a launch.

## Reproduction of the float model

10 of 11 configurations meet the pre-registered rule. **P1 does not** (difference −0.062, 95% CI [−0.103, −0.020]).

**Exploratory follow-up** (not part of the protocol; the files are `results/explore_*.txt`):
- The board is bit-exact with the fixed-point reference. P1, P2 and T1 were therefore rerun in that reference on 2,048 fresh seeds, and in the float model (the J sweep's own kernel) on 4,096 fresh runs each.

| Config | Fixed point: board / fresh seeds / pooled | Float: J held-out / fresh / pooled | Pooled gap |
|---|---|---|---|
| P1 | 0.613 / 0.636 / **0.628** (n=3,072) | 0.675 (and 0.671 in A) / 0.647 / **0.655** (n=6,144) | −0.027 (about 2.6σ) |
| T1 | 0.631 / 0.642 / **0.638** (n=3,072) | 0.675 / 0.667 / **0.668** (n=5,120) | −0.030 (about 2.8σ) |
| P2 | 0.580 / 0.548 / **0.559** (n=3,072) | 0.567 / 0.561 / **0.562** (n=5,120) | −0.003 |

- **The P1 non-reproduction is mostly chance in both cohorts.** The board cohort was low and the original float cohort was high.
- **A small fixed-point deficit of 0.02–0.03 for the T_init=30, S=1560 baselines (P1, T1) remains possible but is not established.** It is about 2.7σ, with three configurations examined after the fact. No mechanism has been identified: the arithmetic is exact for integer fields, and thresholds are quantized at 2⁻¹⁶.
- **Sensitivity of the primary result.** Suppose the baselines are credited with their pooled float p (P1 0.655, T1 0.668) at their measured per-run times.
  - P1's TTS99 becomes 2.53 ms and T1's becomes 2.13 ms.
  - O1's advantage shrinks from 4.19× to **3.99×** against plain, and from 3.72× to **3.37×** against TEC. The conclusion is unchanged.
- O1 itself reproduces the float model (0.930 vs 0.937).

## Context, with caveats

- **STATICA** (JSSC 2021) published TTS of 1.50 ms. That figure is from a simulated chip at 300 MHz, with its own RNG and on the same instance and target. This engine's measured O1 (0.634 ms) is 2.4× lower, but the technology, measurement method and algorithm all differ.
- **GPU** (GH200, [research/tts_platform_20261003](../../research/tts_platform_20261003/REPORT.md)): the best measured batch TTS was 8.69 ms, for plain SCA with 128 runs per launch. Onsager SCA has not been run on the GPU, so the 13.7× ratio mixes algorithm and platform effects.
- **Earlier Snowball V80 engines** (random scan, roulette): 1.81 s and 2.51 s. That is a different algorithm class; the gap is not a like-for-like hardware comparison.

## Engineering headroom (estimates, not measured)

| Change | Effect on O1 per-run cycles (109.5k measured) |
|---|---|
| Reduce per-step overhead from 40 to about 12 cycles (256 lanes, overlapping drain/fill) | About 83k → 0.28 ms per run |
| Plus 4 flips per cycle (four extractors, duplicated coupling ports) | About 48k → 0.16 ms per run |
| Burst table loading | Removes about 34k cycles from single-run latency only |
| Replicate engines (image uses 7.8% LUT / 3.8% BRAM including shell) | Makes the E>1 projections real |

## Board state

- The SCA image (UUID `e818b91b…`) is loaded and idle.
- Recovery: re-stage the corrected Snowball image `/scratch/USER/snowball_v80_20260928/build/board_300mhz_counter_fixed` (UUID `27fadf27…`) with its own `stage_program_image.py` and run the same scoped programming command.
- The programming log is in `results/programming_board_300mhz_*.log` on fpga-host. BAR2 was allocated without needing the rescan helper.

No external review. Single instance (WK2000_1), single target (cut ≥ 33,000).

# v5 engine (5 October 2026), pre-registered in [PROTOCOL_HW_v5.md](PROTOCOL_HW_v5.md)

**Image.**
- `build/board_v5_250mhz`, UUID `ce847edb4aacea1ad5fa78862e3b8960`, 250 MHz, routed setup slack +0.001 ns.
- The 300 MHz v5 build failed timing; the HLS estimate was 206 MHz. v5.1 (HLS estimate 297 MHz) is building at 300 and 275 MHz.

**Exactness and integrity.** Phase 1: 24/24 bit-exact against the 256-lane reference. All 11,456 trials passed host rescoring. Raw data: [results/hw_board_v5_250mhz](results/hw_board_v5_250mhz).

| Prediction | Measured | Verdict |
|---|---|---|
| H1: cycles per flip ≤ 0.33 | **0.275** ± 0.001 (v4: 0.520) | Pass |
| H2: cycles per step ≤ 34 | **31.1** ± 0.1 (v4: 40.8) | Pass |
| H3: table loading ≤ 6 cycles per step per launch | 25.8 (v4: 35.3) | **Fail.** HLS inferred a burst only for the 64-bit kcorr table; the three 32-bit tables missed it. Fix: pack the tables. This affects single-run latency only. |
| H4: O1 t_run ≤ 0.27 ms × 300/250 = 0.324 ms | **0.266 ms** (66.6k cycles vs 109.6k on v4; 1.65× fewer) | Pass |
| H5: quality unchanged vs v4 | All 11 configurations pass (largest gap: P1 +0.049, against v4's low draw) | Pass |

**Measured, one engine at 250 MHz.**
- **O1 TTS99: 0.464 ms** (p = 0.929; v4 at 300 MHz: 0.634 ms).
- **Against the best baselines:**
  - **3.94×** lower than the best plain configuration (P1, 1.83 ms; CI 3.52–4.44).
  - **3.58×** lower than the best TEC configuration (T1, 1.66 ms; CI 3.19–4.03).
- **Cycles per trial:** 1.53–1.65× fewer than v4 across all configurations.

The cycle counts are clock-independent, so at 300 MHz (v5.1, same cycles) O1 would take about 0.222 ms per run, giving a TTS99 of about 0.39 ms.

# Popcount-free rule on the board (TEC-T), pre-registered in [PROTOCOL_HW_TECT.md](PROTOCOL_HW_TECT.md)

**Setup.** v5 image at 250 MHz. Host `host_sca_v5t` adds `--tecT κ` (kconst[t] = κ·T_{t−1}·ramp(t), kcorr = 0); there is no hardware change. Raw data: [results/hw_tect_board_v5_250mhz](results/hw_tect_board_v5_250mhz).

**Exact gate.** 4/4 TEC-T trials (S = 960) are bit-exact against `sca_ref.hpp` with the TEC-T tables.

| Paired cohort, 1,024 trials (same starts as v5 O1) | TEC-T κ1.25 + ramp, q8, T20, S960 | Onsager O1 (λ1.05 + ramp, q8, T12, S960) |
|---|---|---|
| Success p | 0.929 | 0.929 |
| Discordant pairs | 67 TEC-T-only successes | 67 O1-only successes (McNemar p = 1.0) |
| Flips per run | 155.3k | 133.1k |
| Cycles per run | 73.3k | 66.4k |
| TTS99 (Eq. 9) | 0.512 ms | **0.464 ms** |

| Prediction | Result | Verdict |
|---|---|---|
| R1: p in 0.89–0.93 | 0.929 | Pass |
| R2: cycle ratio 1.0–1.2 | 1.103 | Pass |
| R3: TTS ratio 0.9–1.3 | 1.102 | Pass |

**Reading.** The theory's anti-echo rule works on hardware in both forms.
- **Success is identical.** The popcount-free fixed schedule matches the online correction.
- **The online coefficient wins on cost.** It needs about 14% fewer flips per run, giving about 10% lower TTS. This matches the CPU ablation (ABL2/ABL3).


# Three engines: measured multi-engine TTS (6 October 2026), pre-registered in [PROTOCOL_HW_MULTI.md](PROTOCOL_HW_MULTI.md), Amendment 3

## Image

- `build/board_v53_e3_225mhz`, UUID `730d57ff5762ef0970d64b4b60fab8f8`, 225 MHz.
- **Engines.** 3 × v5.3, each confined by a pblock to its own SLR (`scripts/pblocks_e3.xdc`).
- **Timing.** Routing missed by 0.014 ns. A post-route phys_opt (AggressiveExplore) closed it to setup/hold slack 0.000 ns with 0 failing endpoints.
- **Resources.** 1.05M LUT (40.6%), 1.11M FF (21.6%), 890 BRAM tiles (23.8%).
- **Timers.** The timer-structure check found 3 independent start-gated 64-bit counters.
- **Clock record.** The fractional MMCM (94 + 32/64) was verified against the clock-wizard multiplier and the routed period.
- **Earlier attempts** are in DESIGN.md: 4 engines unconstrained; 2 engines unconstrained; 3 floorplanned engines at 250 MHz (−0.137 ns).

## Phases

**Phase 1.** 33/33 trials bit-exact on all three engines (W1–W4).

**Phase 2a (concurrency).** Per-engine cycles per trial are 0.984× v5's measured fit (range 0.984–0.985). **M1 passes**: there is no slowdown from running three engines at once. The ratio is below 1 because v5.3's burst table loading more than offsets its extra pipeline drain.

**Phase 2b.** 1,026 trials per configuration, as 342 rounds × 3 engines, 1 trial per engine per round. Round time is the maximum of the three engines' cycle counters.

| Key | p (trial) | P_round, measured (predicted 1−(1−p)³) | Round time (ms) | TTS99, Eq. 9 on rounds (ms) | Whole-round TTS (ms) |
|---|---:|---:|---:|---:|---:|
| P1 | 0.662 | 0.956 (0.961) | 0.5153 | 0.759 | 1.031 |
| P2 | 0.552 | 0.912 (0.910) | 0.4567 | 0.864 | 0.913 |
| P3 | 0.434 | 0.833 (0.818) | 0.4148 | 1.066 | 1.245 |
| P4 | 0.209 | 0.506 (0.504) | 0.2613 | 1.707 | 1.829 |
| T1 | 0.656 | 0.956 (0.959) | 0.4639 | 0.683 | 0.928 |
| T2 | 0.376 | 0.766 (0.757) | 0.2911 | 0.923 | 1.164 |
| T3 | 0.211 | 0.523 (0.508) | 0.3778 | 2.348 | 2.645 |
| O1 | 0.929 | 1.000 (1.000) | 0.3374 | 0.337 | 0.337 |
| O2 | 0.757 | 0.988 (0.986) | 0.2042 | 0.211 | 0.408 |
| O3 | 0.519 | 0.904 (0.889) | 0.1686 | 0.332 | 0.337 |
| O4 | 0.327 | 0.699 (0.695) | 0.1198 | 0.460 | 0.479 |

**M2 (independence) passes** for all 11 configurations.

## Primary result

- **Best measured 3-engine Onsager TTS99 = 0.211 ms** (O2). **M3″ passes** (≤ 0.278 ms).
- O1 succeeded in every one of the 342 rounds, so its TTS equals its round time, 0.337 ms.
- **M4″ passes.** Onsager vs best plain (P1): **3.59×** [95% CI 2.78, 4.81]. Onsager vs best TEC (T1): **3.23×** [2.48, 4.34]. Both are post-hoc minima over the measured configurations, which favours the baselines.
- **Against one engine:** the single-engine v5 at 250 MHz measured 0.464 ms (O1). Three engines at 225 MHz give **0.211 ms**, a measured 2.2× reduction despite the lower clock.


# v6 RTL engine × 9 at 300 MHz (6 October 2026), pre-registered in [PROTOCOL_HW_V6.md](PROTOCOL_HW_V6.md) with Amendments 1 and 2

## Image

- `board_v6_e9_300mhz`, UUID `33da866c81e6bf357f128e58ec48dd9f`.
- **Engines.** 9 × (RTL core v6.1 + HLS front end), 3 per SLR.
- **Timing.** Closed at 300 MHz with setup +0.030 ns and hold +0.010 ns after post-route phys_opt.
- **Resources.** 1.27M LUT (49%), 1.29M FF (25%), 2,488.5 BRAM tiles (66.5%).

## Frozen phases

**Phase 1.** 45/45 trials bit-exact on all 9 engines (W1–W4).

**M1 (cycle model).** Concurrency ratio 1.025 against the RTL-simulation cycle model. **Pass.**

**Rounds (primary TTS = Eq. 9 on measured rounds; secondary = per-trial p with tested independence).**

| Key | p | P_round (pred) | t_round ms | TTS99 primary | TTS99 secondary |
|---|---:|---:|---:|---:|---:|
| P1 | 0.652 | 1.000 (1.000) | 0.4216 | 0.4216 | 0.2046 |
| P2 | 0.556 | 1.000 (0.999) | 0.3817 | 0.3817 | 0.2405 |
| P3 | 0.439 | 1.000 (0.994) | 0.3542 | 0.3542 | 0.3140 |
| P4 | 0.207 | 0.877 (0.876) | 0.2221 | 0.4878 | 0.4897 |
| T1 | 0.667 | 1.000 (1.000) | 0.3838 | 0.3838 | 0.1785 |
| T2 | 0.392 | 0.974 (0.989) | 0.2401 | 0.3039 | 0.2466 |
| T3 | 0.202 | 0.855 (0.869) | 0.3319 | 0.7908 | 0.7517 |
| O1 | 0.933 | 1.000 (1.000) | 0.2883 | 0.2883 | 0.0545 |
| O2 | 0.764 | 1.000 (1.000) | 0.1667 | 0.1667 | 0.0591 |
| O3 | 0.515 | 1.000 (0.999) | 0.1409 | 0.1409 | 0.0996 |
| O4 | 0.326 | 0.965 (0.971) | 0.0981 | 0.1349 | 0.1275 |
| X1 | 0.428 | 0.996 (0.993) | 0.0958 | 0.0813 | 0.0877 |
| X2 | 0.455 | 0.996 (0.996) | 0.1036 | 0.0879 | 0.0873 |
| X3 | 0.970 | 1.000 (1.000) | 0.3575 | 0.3575 | 0.0520 |
| X4 | 0.950 | 1.000 (1.000) | 0.2457 | 0.2457 | 0.0419 |
| B1 | 0.461 | 0.991 (0.996) | 0.2724 | 0.2649 | 0.2258 |
| B2 | 0.819 | 1.000 (1.000) | 0.4727 | 0.4727 | 0.1416 |
| B3 | 0.392 | 0.974 (0.989) | 0.2400 | 0.3039 | 0.2466 |
| B4 | 0.667 | 1.000 (1.000) | 0.3839 | 0.3839 | 0.1786 |
| X5 | 0.338 | 0.982 (0.976) | 0.0846 | 0.0963 | 0.1048 |

**M2 (independence).** Passes for all 19 configurations.

## Outcome against the frozen predictions

| ID | Result | Value | Verdict |
|---|---|---|---|
| M3 | Best Onsager (primary) | O4 0.135 ms against threshold 0.105 ms | **FAIL** |
| M4 | Onsager vs plain / vs TEC (primary) | 2.63× / 2.25× | Pass |
| M5 | ≤ 0.13 ms | 0.135 ms | **FAIL**, by 0.005 ms |
| M6 | Secondary ≤ 0.05 ms, ratios ≥ 2.5 | O1 0.0545 ms; 3.75× / 3.28× | **FAIL** (absolute only) |
| A1 | min(X1, X2) primary ≤ 0.08 ms | 0.0813 ms | **FAIL**, by 0.0013 ms |
| A2 | min(X3, X4) secondary ≤ 0.04 ms | 0.0419 ms | **FAIL**, by 0.0019 ms |
| A3 | Ratios against re-selected baselines | 3.26× / 3.38× | Pass |

## Why the absolute predictions failed

- **The cause.** The predictions modelled per-launch table streaming as S + 300 cycles. On the board, each one-trial launch costs about **20.8 extra cycles per step**: excess over the core model fits about 20.8·S across all 19 configurations. In the 64-trials-per-launch concurrency phase, that excess falls to 1.4 cycles per step.
- **Explanation.** All 9 engines re-read the same tables every round through two NoC ports configured for 2 GB/s.
- **Not the core.** The core's own cycle model holds (M1).
- **Fix: v6.2.** v6.2 keeps the tables resident after round 0. It is verified in simulation and its build is the next item.

## Headline

- **Measured best TTS99:**
  - primary: **0.081 ms** (X1, TEC-T S320);
  - secondary: **0.042 ms** (X4, TEC-T S760).
- **Against our 3-engine v5.3 result (0.211 ms):** 2.6× lower.
- **Against Toshiba bSBM's 0.26 ms TTT99** (≈33,004, measured FPGA): 3.2× lower on the primary estimator. The targets differ slightly.
- **The frozen O1–O4 schedules alone:** 0.135 ms primary.


# v6.2 RTL engine × 12 at 275 MHz with resident tables (6 October 2026), pre-registered as PROTOCOL_HW_V6 Amendment 3

## Image

- `board_v62_e12_275mhz`, UUID `c52770e9474879f38c65c7cd709d498c`.
- **Engines.** 12 × (RTL core v6.2 + HLS front end), 4 per SLR.
- **Timing.** Setup +0.045 ns and hold +0.010 ns after post-route phys_opt.
- **Resources.** 1.68M LUT (65%), 1.71M FF (33%), 3,318 BRAM tiles (88.7%).
- **Tables.** Resident after each invocation's first round (`--keep-tables`).

## Results

**Exactness.** Phase 1: 60/60 bit-exact on all 12 engines, with resident tables.

**R1 (cross-image reproducibility).** 41,040 trials compared against the 9-engine v6.1 image, trial id by trial id. **Zero** mismatches in cut or flips.

**Cycle model (M1).** Ratio 1.020. **Independence (M2).** Passes for all configurations.

| Key | p | P_round | t_round ms | TTS99 primary | TTS99 secondary |
|---|---:|---:|---:|---:|---:|
| P1 | 0.652 | 1.000 | 0.3146 | 0.3146 | 0.1145 |
| P2 | 0.556 | 1.000 | 0.2700 | 0.2700 | 0.1276 |
| P3 | 0.439 | 1.000 | 0.2396 | 0.2396 | 0.1593 |
| P4 | 0.207 | 0.947 | 0.1533 | 0.2397 | 0.2534 |
| T1 | 0.667 | 1.000 | 0.2722 | 0.2722 | 0.0950 |
| T2 | 0.392 | 1.000 | 0.1730 | 0.1730 | 0.1333 |
| T3 | 0.202 | 0.906 | 0.2151 | 0.4181 | 0.3653 |
| O1 | 0.933 | 1.000 | 0.2352 | 0.2352 | 0.0333 |
| O2 | 0.764 | 1.000 | 0.1341 | 0.1341 | 0.0356 |
| O3 | 0.515 | 1.000 | 0.1055 | 0.1055 | 0.0559 |
| O4 | 0.326 | 0.994 | 0.0775 | 0.0695 | 0.0756 |
| X1 | 0.428 | 1.000 | 0.0771 | 0.0771 | 0.0529 |
| X2 | 0.455 | 1.000 | 0.0874 | 0.0874 | 0.0553 |
| X3 | 0.970 | 1.000 | 0.2851 | 0.2851 | 0.0311 |
| X4 | 0.950 | 1.000 | 0.1984 | 0.1984 | 0.0254 |
| B1 | 0.461 | 0.994 | 0.2082 | 0.1865 | 0.1295 |
| B2 | 0.819 | 1.000 | 0.3697 | 0.3697 | 0.0831 |
| B3 | 0.392 | 1.000 | 0.1730 | 0.1730 | 0.1333 |
| B4 | 0.667 | 1.000 | 0.2722 | 0.2722 | 0.0950 |
| X5 | 0.338 | 1.000 | 0.0686 | 0.0686 | 0.0638 |

## Outcome against the frozen predictions (Amendment 3): all pass

| ID | Value |
|---|---|
| M3, M5 | Best Onsager, O4 = 0.0695 ms (predicted 0.071; thresholds 0.08 and 0.13) |
| M4 | 3.45× vs plain (P3), 2.49× vs TEC (T2) |
| M6 | O1 secondary 0.0333 ms |
| A1 | min(X1, X2, X5) primary = 0.0686 ms (X5; predicted 0.070) |
| A2 | min(X3, X4) secondary = 0.0254 ms (X4; predicted 0.027) |
| A3 | 2.52× / 3.27× against re-selected baselines |

## Headline (measured on the V80; target cut ≥ 33,000)

- **Primary TTS99 (Eq. 9 on measured 12-engine rounds): 0.069 ms** (X5 0.0686; frozen O4 0.0695).
- **Secondary TTS99 (per-trial p, independence tested): 0.025 ms** (X4).
- **Against Toshiba bSBM** (TTT99 0.26 ms, 99% of 33,337, measured FPGA): 3.8× faster on the primary estimator and 10× on the secondary.
- **Against our earlier results:** 3.1× faster than the 3-engine v5.3 (0.211 ms) and 6.7× faster than the single-engine v5 (0.464 ms).

# Board power and energy, 12 × v6.2 at 275 MHz (7 October 2026)

Raw data:
- `results/power_board_v62_e12_275mhz_20261006T235156Z/` (`hwmon.txt`, `load_*.json`, `power_summary.json`)
- Script: `scripts/measure_power_v6.sh`, analysed by `scripts/analyze_power_v6.py`
- Host: `src/host_sca_multi_v6_power.cpp` rev 2. This is a separate binary that adds `--power-run`. The frozen `host_sca_multi_v6.cpp` and `build/host_sca_v62` are unchanged.

**Method.**
- The board ran each configuration back to back on all 12 engines for 150 s: 512 trials per engine per launch, tables resident, no readback.
- `Total_Power` (board input, the sum of the 12 V PEX, AUX1 and AUX2 rails and the 3.3 V rail) and the VCCINT current were sampled every 0.25 s from the AMI driver's hwmon node.
- Before and after every load there was 60 s of idle with the image loaded. The first 10 s of each phase were dropped.

**Measurement pitfall found and fixed.**
- The first run sampled with `ami_tool sensors`, which stalls the host's BAR accesses.
- That left the engines idle between launches: the device was busy only 71–73% of the time, so that run understates load power. It is kept as a record (`..._20261006T233706Z`).
- A 60 s test (`throttle_test_v62_*`) showed host submit-to-done time equal to cycles/275 MHz (ratio 1.0008) once sampling stopped. **The kernel clock is not throttled under sustained load.**
- In the hwmon run the devices were busy 99.5–99.8% of the time (host/device 1.0003–1.0008).

| Phase | Board power (W) | VCCINT (W) |
|---|---:|---:|
| Idle, before / between / after | 62.4, 63.9, 64.4, 64.4 (rises as the board warms) | 36.0–37.5 |
| X5 load (S = 280, TEC-T) | 85.15 ± 0.83 | 54.6 |
| X4 load (S = 760, TEC-T) | 85.84 ± 0.65 | 55.5 |
| O1 load (S = 960, Onsager) | 85.45 ± 0.82 | 55.5 |

**Dynamic power** (load minus the mean of the two neighbouring idle phases) is 21.1–22.0 W for 12 engines, about 1.8 W per engine. VCCINT accounts for 18.0 W of it.

Vivado's vectorless routed estimate (`route_report_power_0.rpt`, low confidence, Tj = 100 °C) is 101.9 W on chip, 3.65–3.80 W per engine. It overestimates the measured engine power by about 2×, as expected for default toggle rates.

| Key | Primary TTS99 | Energy to solution, board (85 W × TTS) | Dynamic part | Board energy per trial | Sustained trials/s |
|---|---:|---:|---:|---:|---:|
| **X5** | **0.0686 ms** | **5.85 mJ** | 1.51 mJ | 477 µJ | 178,700 |
| X4 | 0.198 ms | 17.0 mJ | 4.30 mJ | 1,399 µJ | 61,400 |
| O1 | 0.235 ms | 20.1 mJ | 4.95 mJ | 1,428 µJ | 59,900 |

Board energy uses whole-board input power, including the 62–64 W the shell, HBM, DDR, PS and static power draw at idle.

# v6.4 RTL engine × 12 at 250 MHz (7 October 2026), PROTOCOL_HW_V6 Amendment 4: all predictions pass

**Image.**
- `board_v64_e12_250mhz`, UUID `8669ad38a456bcdf17ad8485fa0cec7a`, PDI SHA256 `67e591f3…`.
- Timing: routed setup +0.022 ns, hold +0.010 ns, after post-route phys_opt.
- Programmed 02:25 UTC by the unattended chain `scripts/bringup_v64.sh` (log `logs/bringup_v64.log` on fpga-host): UUID match and BAR2 OK.
- Run with `KEEP_TABLES=1`, same trial ids as the 12 × v6.2 run.
- Raw data: `results/hwv6_board_v64_e12_250mhz/` (`hwm_summary.json`, `amendment4.json`).

**Engine change.** 8 flip extractors and 4 coupling copies (2 BRAM, 2 URAM via XPM). Cycle model 0.1424·flips + 18.09·S + 886.

**Exactness.**
- Phase 1: W1–W4 exact on all 12 engines.
- R1: all 41,040 common trials have identical cut and flips to the 12 × v6.2 image. Same algorithm, bit for bit.

| ID | Result | Verdict |
|---|---|---|
| M1 | Concurrency ratio 1.024 against the v6.4 model (within 3%) | Pass |
| M2 | Independence holds for every configuration | Pass |
| R1 | 41,040 trials identical | Pass |
| M3 | Best Onsager primary: O4 0.0528 ms (≤ 0.06; predicted 0.053) | Pass |
| M4 | 3.74× vs plain (P4), 2.58× vs TEC (T2) | Pass |
| M5 | ≤ 0.13 ms; 4.9× under SOTA 0.26 ms | Pass |
| M6 | O1 secondary 0.0246 ms (≤ 0.03); 3.55× / 3.07× | Pass |
| A1 | min(X1, X2, X5) primary = **0.0503 ms** (X5; predicted 0.051) | Pass |
| A2 | min(X3, X4) secondary = **0.0181 ms** (X4; predicted 0.018) | Pass |
| A3 | 2.71× / 3.36× against re-selected baselines | Pass |

| Key | p | P_round | t_round ms | TTS99 primary | TTS99 secondary |
|---|---:|---:|---:|---:|---:|
| X5 | 0.338 | 1.000 | 0.0503 | **0.0503** | 0.0468 |
| O4 | 0.326 | 0.994 | 0.0590 | 0.0528 | 0.0575 |
| X1 | 0.428 | 1.000 | 0.0567 | 0.0567 | 0.0389 |
| X2 | 0.455 | 1.000 | 0.0628 | 0.0628 | 0.0397 |
| X4 | 0.950 | 1.000 | 0.1419 | 0.1419 | **0.0181** |
| O1 | 0.933 | 1.000 | 0.1733 | 0.1733 | 0.0246 |
| B1 (plain) | 0.461 | 0.994 | 0.1568 | 0.1405 | 0.0975 |
| B2 (plain) | 0.819 | 1.000 | 0.2713 | 0.2713 | 0.0610 |
| T2 (TEC) | 0.392 | 1.000 | 0.1366 | 0.1366 | 0.1052 |
| T1 (TEC) | 0.667 | 1.000 | 0.2164 | 0.2164 | 0.0755 |

**Cycles per step on the board.** X5 13,047 cycles per trial (46.6 per step, 186 ns); O4 14,501 (40.3 per step, 161 ns). This is 0.68–0.72× of v6.2 with identical flips.

**Board power** (`results/power_board_v64_e12_250mhz_20261007T031120Z/`; hwmon, 512 trials per launch, 150 s loads; devices 99.4–99.8% busy, host/device ≤ 1.001):

| Phase | Board W | VCCINT W |
|---|---:|---:|
| Idle (before, then between and after loads) | 66.2, 71.3, 73.3, 73.0 | 39.2–44.5 |
| X5 load | 108.2 ± 2.3 | 73.8 |
| X4 load | 113.1 ± 1.8 | 77.3 |
| O1 load | 109.0 ± 1.9 | 73.8 |

- Dynamic power (load minus bracketing idle): 35.8–40.8 W.
- X5 energy to solution: **5.45 mJ** board (108.18 W × 0.0503 ms), 1.98 mJ dynamic.
- X5 board energy per trial: 446 µJ. Sustained: 242,700 trials/s.

**Resources** (hierarchical, routed; `results/util_20261007/`):

| | v6.2 per engine | v6.4 per engine | v6.4 total |
|---|---:|---:|---:|
| LUT | 130,473 | 149,512 | 1,906,208 (74.0%) |
| FF | 126,456 | 131,555 | 1,767,564 (34.3%) |
| RAMB36 | 258 | 258 | 3,312 (+12 RAMB18) |
| URAM | 2 | 66 | 792 (41.1%) |
| DSP | 282 | 282 | 3,384 (31.2%) |

**Headline (12 × v6.4, measured).** Primary TTS99 **0.050 ms**, 5.2× faster than bSBM's 0.26 ms. Secondary 0.018 ms.
