# v6 hardware measurement protocol: RTL engine × 9 (frozen 6 October 2026, before any v6 board data)

## Design under test

- **Engine v6.1.** Hand-written RTL core (`src/v6/sca_core.v`) plus a small HLS front end (`src/v6/sca_mover.cpp`). The front end keeps the v5.2/5.3 AXI-Lite register map, timers (offsets 0x18/0x24/0x30) and HBM layout.
  - **Why.** v5.3 was 342k LUT and 356k FF per engine and was limited by HLS loop-outlining (duplicated state, high-fanout clock enables).
  - **v6.1 size.** 128k LUT, 126k FF, 258 BRAM and 282 DSP per engine.
  - **Standalone timing.** The engine alone, constrained to 3×3 clock regions of SLR1, meets 300 MHz with +0.323 ns slack and 0 of 272k endpoints failing.
- **Replication.** 9 engines, 3 per SLR (`scripts/pblocks_v6_e9.xdc`).
  - Fully registered AXI register slices on every control and memory port.
  - Two memory masters per engine, on two SmartConnects into NoC S04/S05.
  - Host: `src/host_sca_multi_v6.cpp`, which is `host_sca_multi.cpp` with up to 12 engines and a 1 MiB trace stride.
- **Image selection (fixed now).** The first timing-clean image in this order:
  1. 9 engines at 300 MHz (build `board_v6_e9_300mhz`, started 13:24 UTC);
  2. otherwise 9 engines at 275 MHz;
  3. otherwise 9 engines at 250 MHz.

  Any other engine count needs an amendment written before its data.
- **Verification before hardware.**
  - **RTL simulation (iverilog), bit-exact against `sca_ref.hpp`.** Spins, cut, flips and n_lin trace all match:
    - 8 short trials (plain, Onsager+ramp, TEC; random stream stalls; couplings resident across launches);
    - one full trial each of O1, O2, O3, O4, P1 and T1.
  - **Drain margin.** DRAIN=2 is exact and DRAIN=1 fails, confirming the field-update pipeline analysis.
  - **HLS front end.** C simulation and C/RTL co-simulation pass.
  - **Block design.** The 9-engine design validates.

## Phases (`scripts/run_hw_v6.sh`; trial counts derived from E = 9; seed 20261004)

**Phase 1: exact gate on all engines.** Every trial must match `sca_ref.hpp` exactly (spins, cut, flips, n_lin trace). Any mismatch stops the protocol.

| ID | Configuration | Trials (trials per launch) |
|---|---|---|
| W1 | Plain q4, T30, S40 | 9 (1) |
| W2 | Onsager λ0.7 + ramp, q4, T20, S560 | 18 (2) |
| W3 | TEC −8, S40 | 9 (1) |
| W4 | O1, S960 | 9 (1) |

**Phase 2a: concurrency.** O1 with 9 engines × 128 trials, 64 trials per launch per engine.

**Phase 2b: measured 9-engine rounds.**
- **Configurations.** All 11 frozen configurations (P1–P4, T1–T3, O1–O4, as in `run_hw_multi.sh`).
- **Trials.** 2,052 per configuration: 228 rounds × 9 engines, 1 trial per engine per round, trial ids 0–2051 (round r, engine e → id 9r + e).
- **Round time.** The maximum of the 9 engines' start-to-done cycles. It includes the per-launch header and table streaming. Couplings load in round 0, which is excluded from timing.

## Analysis (fixed in advance; `scripts/analyze_hw_v6.py`, cycle model `results/v6_cycle_model_sim.json`)

- **Primary.** Measured 9-engine TTS99 = t_R · ln(0.01)/ln(1 − P_round). P_round is measured over the 228 rounds; P_round = 1 gives one round time.
  - Comparison: best of O1–O4 against best of P1–P4 and best of T1–T3 (post-hoc minima, which favours the baselines).
  - Uncertainty: bootstrap 95% intervals over rounds (10,000 resamples).
- **Secondary (pre-registered).** TTS99 = t_R · ln(0.01)/(9 · ln(1 − p)), with p the per-trial success over all 2,052 trials.
  - This is valid only where M2 holds.
  - It removes the P_round = 1 floor of the primary form. With 9 engines the primary form floors O1, O2 and most baselines at one round, which compresses the differences between rules.
- Whole-round TTS (t_R · ⌈R99⌉) is also reported.

## Predictions

**Basis.** The algorithm is bit-identical to v5.3, so trial ids 0–1025 have known outcomes and flips from the 3-engine run (`results/hwm_board_v53_e3_225mhz`). This means the success side of these predictions is not an independent test. What the run tests is:
- exactness on hardware;
- cycles and clock;
- concurrency;
- the 1,026 new trials (ids 1026–2051).

**Inputs used.**
- Cycle model from RTL simulation: cycles/trial = 0.2714·flips + 16.91·S + 1148.
- Per-launch streaming: S + 300 cycles.
- Unseen trials resampled.
- Output: `results/v6_predictions_e9.json`.
- At 300 MHz the predicted primary best is O4 at 0.088 ms (90% range 0.071–0.104 ms).

| ID | Prediction | Pass condition |
|---|---|---|
| M1 | Cycle model holds on hardware with 9 concurrent engines | Mean of (launch cycles / Σ model cycles) in 2a within 3% of 1 |
| M2 | Engines are independent | For every configuration, P_round inside the 95% interval implied by p |
| M3 | 9-engine gain | Best measured Onsager TTS99 (primary) ≤ 0.105 ms × (300/f) |
| M4 | Onsager advantage persists (primary) | Best Onsager / best plain and best Onsager / best TEC both ≥ 2× (point estimates; predicted 2.5× and 2.2×) |
| M5 | Faster than the measured state of the art at the 33,000 level | Best measured Onsager TTS99 (primary) ≤ 0.13 ms, i.e. at least 2× below Toshiba bSBM's 0.26 ms TTT99 (Goto et al., Sci. Adv. 2021; target 99% of 33,337, i.e. about 33,004) |
| M6 | Secondary estimator | Best Onsager ≤ 0.05 ms × (300/f), and Onsager / plain and Onsager / TEC both ≥ 2.5× (predicted 0.042 ms, 3.4× and 2.9× at 300 MHz) |

## Limits stated in advance

- Single instance (WK2000_1); target cut ≥ 33,000.
- At the best-known cut, 33,337, this engine is far from the state of the art. Toshiba's GbSBM reports TTS99 9.6 ms; SCA-type schedules here reach 33,337 in under 1% of runs even at S = 8000 (CPU model).
- The Toshiba comparison mixes targets (≥ 33,000 here, about 33,004 there) and metrics (TTS99 here, TTT99 there).
- The project's earlier preprint figures (0.085 and 0.128 ms) are not used. The forensic review in HANDOVER.md found them unsupported.
- No external review.

## Amendment 1 (6 October 2026, before any v6 board data)

**What it adds.** Eight configurations, run in phase 2b after O4 (`run_hw_v6.sh`), each with 2,052 trials. They come from a separate pre-specified re-selection for the v6 cost model with 9 engines (`research/v6_schedule_20261006/PROTOCOL.md`: pilot of 256 runs per configuration, Wilson-bound selection, held-out of 2,048 fresh runs).

| Key | Schedule | Selected for | Held-out p | Held-out TTS (300 MHz) |
|---|---|---|---|---|
| X1 | TEC-T κ1.75 + ramp, q8, T15, S320 | primary | 0.441 | 0.066 ms |
| X2 | Onsager λ0.9 + ramp, q6, T15, S320 | primary | 0.445 | 0.069 ms |
| X3 | Onsager λ1.05 + ramp, q8, T12, S1260 | secondary | 0.975 | 0.032 ms |
| X4 | TEC-T κ2.0 + ramp, q8, T15, S760 | secondary | 0.950 | 0.033 ms |
| B1 | Plain q8, T40, S960 | primary (baseline) | 0.466 | 0.198 ms |
| B2 | Plain q6, T40, S1560 | secondary (baseline) | 0.814 | 0.108 ms |
| B3 | TEC J_v −4, q8, T30, S960 | primary (baseline) | 0.410 | 0.161 ms |
| B4 | TEC J_v −4, q8, T30, S1560 | secondary (baseline) | 0.653 | 0.126 ms |

**Caveats on the selection.** X1–X4 were chosen for this cost model and engine count, and the baselines B1–B4 were re-selected the same way. Several selections sit on a grid edge: TEC-T κ, q8, S1260, and plain T40/S1560.

**Host change.** `host_sca_multi_v6.cpp` gains `--tecT κ`, ported unchanged from `host_sca.cpp` and measured bit-exact on v5 in PROTOCOL_HW_TECT. The frozen phases do not use it. New hashes are in `PROTOCOL_HW_V6_A1.sha256`.

**Analysis.** The frozen analysis above is unchanged. The extension is reported separately and never pooled with M1–M6.

| ID | Prediction | Pass condition |
|---|---|---|
| A1 | Fastest re-selected schedule, primary | min(X1, X2) TTS99 (primary) ≤ 0.08 ms × (300/f) |
| A2 | Fastest re-selected schedule, secondary | min(X3, X4) TTS99 (secondary) ≤ 0.04 ms × (300/f) |
| A3 | Advantage over equally re-selected baselines | min(B1, B3) / min(X1, X2) ≥ 2 (primary) and min(B2, B4) / min(X3, X4) ≥ 2.5 (secondary) |

## Amendment 2 (6 October 2026, before any v6 board data): 12 engines at 275 MHz

**Second image.** `board_v6_e12_275mhz`: 12 engines (4 per SLR, `pblocks_v6_e12.xdc`), same RTL, started 16:18 UTC.
- It is measured with the same phases (E = 12: W1/W2/W3/W4 = 12/24/12/12, NC 1,536, NR 2,052 = 171 rounds) into its own result directory.
- If it does not close timing, there is no 12-engine result. The 9-engine image remains the primary measurement.

**Added configuration.** X5 = TEC-T κ1.75 + ramp, q8, T15, S280. It was selected for 12 engines (`research/v6_schedule_20261006` addendum 2): held-out p 0.343, primary TTS 0.071 ms at 275 MHz.
- On the 9-engine image X5 runs as exploratory, with no prediction.

**Correction to Amendment 1's sweep (before data).**
- The 9-engine round-maximum factor is about 1.08 (from known trial flips), not 1.035. Amendment 1's held-out TTS values are about 4% optimistic.
- Its thresholds (0.08 and 0.04 ms) are unchanged.

**12-engine predictions.**
- Known trials ids 0–1019 and cycle model, at 275 MHz (`results/v6_predictions_e12_275.json`).
- Primary (O1–O4):
  - best O4, 0.080 ms (90% range 0.072–0.091);
  - best plain P3, 0.239 ms;
  - best TEC T2, 0.178 ms.
- Secondary: O1 0.035 ms; plain P1 0.114 ms; TEC T1 0.100 ms.

| ID (12 engines) | Pass condition |
|---|---|
| M1, M2 | As above |
| M3 | Best measured Onsager TTS99 (primary, O1–O4) ≤ 0.095 ms |
| M4 | Ratios against plain and TEC ≥ 2× (primary) |
| M5 | ≤ 0.13 ms |
| M6 | Best Onsager (secondary) ≤ 0.042 ms, with ratios ≥ 2.5× |
| A1 | min(X1, X2, X5) primary ≤ 0.085 ms |
| A2 | min(X3, X4) secondary ≤ 0.035 ms |
| A3 | min(B1, B3) / min(X1, X2, X5) ≥ 2 and min(B2, B4) / min(X3, X4) ≥ 2.5 |

## Amendment 3 (6 October 2026, after the 9-engine v6.1 run and before any v6.2 data): v6.2 × 12 engines at 275 MHz

**Why.** The 9-engine run showed that per-launch table streaming costs about 21.9 cycles per step on the board (NoC-bound). That cost, not the core, made M3, M5, M6, A1 and A2 fail.

**Change (v6.2).** Tables stay resident after each invocation's first round:
- `load_J` bit 1 tells the front end to skip the table stream, and header bit 34 tells the core to skip LOADT.
- The host gains `--keep-tables`, and the runner passes it with KEEP_TABLES=1.
- Verification:
  - RTL simulation is bit-exact, including a launch that reuses resident tables (660 words, 0 errors);
  - the front-end C sim and co-sim pass with `load_J` = 0, 1, 2 and 3.
- The v6.1 sources are kept as `src/v6/sca_core_v61.v` and `src/v6/sca_mover_v61.cpp`.

**Image.** `board_v62_e12_275mhz` (12 engines, 4 per SLR), started 17:10 UTC. It replaces the stopped v6.1 12-engine build.

**Run.** `KEEP_TABLES=1 SCA_VARIANT=v62 NUM_ENGINES=12`, with the same phases, configurations and trial ids. The exact gate therefore also covers resident tables on hardware.

**Predictions.** These come from the 9-engine run's per-trial outcomes and launch cycles (same ids, bit-exact), minus the fitted table stream 21.87·S − 1293 (residual sd 870 cycles), regrouped into 12-engine rounds at 275 MHz (`results/v62_predictions_e12_275.json`):

| Estimator | Ours | Plain | TEC |
|---|---|---|---|
| Primary (O1–O4) | O4 0.071 ms | P4 0.260 ms | T2 0.186 ms |
| Primary (X configurations) | X5 0.070 ms | | |
| Secondary | X4 0.027 ms; O1 0.034 ms | P1 0.123 ms | T1 0.103 ms |

| ID | Pass condition |
|---|---|
| Phase 1 | All exact (with resident tables) |
| M1, M2 | As before |
| R1 (new) | Every trial id of every configuration has the same cut and flips as in the 9-engine v6.1 run (deterministic, bit-exact algorithm) |
| M3 | Best Onsager (primary, O1–O4) ≤ 0.08 ms |
| M4 | Ratios ≥ 2× |
| M5 | ≤ 0.13 ms |
| M6 | Secondary ≤ 0.04 ms with ratios ≥ 2.5× |
| A1 | min(X1, X2, X5) primary ≤ 0.08 ms |
| A2 | min(X3, X4) secondary ≤ 0.031 ms |
| A3 | min(B1, B3) / min(X1, X2, X5) ≥ 2, and min(B2, B4) / min(X3, X4) ≥ 2.5 |

**Commands.**
- `analyze_hw_v6.py <dir> results/v6_cycle_model_sim.json 0.08 0.26 0.04`
- `analyze_v62_amend3.py <dir> results/hwv6_board_v6_e9_300mhz`

## Amendment 4 (6 October 2026, after the 12-engine v6.2 run and before any v6.4 data): v6.4 × 12 engines at 250 MHz

**Why.** On the 12-engine v6.2 image all Amendment 3 predictions passed. The cycle breakdown (RTL simulation, O4) showed that 4-way flip extraction is 64% of step cycles.

**Change (v6.4).** 8 extractors and 4 coupling copies: 2 BRAM, plus 2 URAM through XPM `xpm_memory_tdpram` under `SYNTHESIS` (URAM inference was rejected for its write mode).
- **Verification.** RTL simulation is bit-exact on the 660-word suite and on full O1–O4, P1 and T1 trials.
- **Cycle model.** 0.1424·flips + 18.09·S + 886 (`results/v64_cycle_model_sim.json`): 0.68–0.73× of v6.2 per trial.
- **Timing fixes (v6.4b).** An explicit two-level tree for the accumulate-mode bit, and one more stage on the seed-shift strobe. Standalone, v6.4a met 275 MHz with +0.119 ns.
- **Previous versions** are kept as `src/v6/sca_core_v64a.v` and `src/v6/sca_core_v62.v`.

**Image.** `board_v64_e12_250mhz` (12 engines, 4 per SLR, resident tables), started 22:18 UTC. 250 MHz was chosen because v6.4a's standalone Fmax was about 284 MHz and board-level results run about 15% lower.

**Run.** `KEEP_TABLES=1 SCA_VARIANT=v64 NUM_ENGINES=12`, same phases and trial ids. The trial ids and rounds equal the 12-engine v6.2 run, so success patterns are known exactly.

**Predictions.** The v6.2 12-engine launch cycles per trial, scaled by the v6.4/v6.2 model ratio, at 250 MHz (`results/v64_predictions_e12_250.json`):

| Estimator | Ours | Plain | TEC |
|---|---|---|---|
| Primary (O1–O4) | O4 0.053 ms | P4 0.197 ms | T2 0.137 ms |
| Primary (X configurations) | X5 0.051 ms | | |
| Secondary | X4 0.018 ms; O1 0.025 ms | | |

| ID | Pass condition |
|---|---|
| Phase 1 | All exact |
| M1 | Against the v6.4 cycle model, within 3% |
| M2 | Independence holds |
| R1 | Identical cut and flips per trial id against the 12-engine v6.2 run |
| M3 | ≤ 0.06 ms |
| M4 | Ratios ≥ 2× |
| M5 | ≤ 0.13 ms |
| M6 | Secondary ≤ 0.03 ms, with ratios ≥ 2.5× |
| A1 | min(X1, X2, X5) primary ≤ 0.06 ms |
| A2 | min(X3, X4) secondary ≤ 0.022 ms |
| A3 | Ratios ≥ 2 and ≥ 2.5 |

**Commands.**
- `analyze_hw_v6.py <dir> results/v64_cycle_model_sim.json 0.06 0.26 0.03`
- `analyze_v64_amend4.py <dir> results/hwv6_board_v62_e12_275mhz`

**If the image does not close timing at 250 MHz,** there is no v6.4 result. The 12-engine v6.2 image remains the measured best.
