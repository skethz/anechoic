# Multi-bit engine hardware protocol (PROTOCOL_HW_MB): K = 2 couplings, 12 engines

Frozen on 7 October 2026, before any multi-bit board data. Hashes are in `PROTOCOL_HW_MB.sha256`, and the bring-up script (`scripts/bringup_mb2.sh`) refuses to program the board if any frozen file has changed. It extends PROTOCOL_HW_V6 (Amendments 1–4), which is not modified, to a new engine.

## Design under test

- **Engine.** `src/v6/sca_core_mb.v` with K = 2: a sign plane and one magnitude plane, so J ∈ {−1, 0, +1}. It follows `MULTIBIT_SPEC.md` and is derived from v6.4b (`src/v6/sca_core.v`, SHA-256 32e1c2f8…, which is unchanged).
  - **Single-copy residue banks.** Extractor e only issues spins with x mod 8 = e: o_x[2:0] = E in the RTL, and this is asserted in every simulation. Each port therefore has its own bank holding rows x ≡ e (mod 8) once.
    - Residues 0–3 are block RAM (one 256 × 64 bank per plane). Residues 4–7 are UltraRAM (one URAM per plane pair, via XPM).
    - Per engine: 258 RAMB36 (as in v6.4b) and 130 URAM (v6.4b: 66), from out-of-context synthesis.
  - **Field update.** Two 8-input population counts per magnitude plane.
  - **Diagonal and fields.** J_xx = 0, so the stored field is the true field (v6.4b stored h − s).
  - **Unchanged.** Step pipeline, latencies, DRAIN = 2 and the stream protocol, except that the coupling load is K × 32768 beats.
- **Front end.** `src/v6/sca_mover_mb.cpp` (MB_K = 2). Same register map and timers as v6.4. HLS C simulation and C/RTL co-simulation pass.
- **Host.** `src/host_sca_multi_mb.cpp` (MB_K = 2), with matrices from `src/mb_common.hpp`.
  - The golden model is `src/sca_ref.hpp` `run_trial(..., dense)`, unchanged.
  - Every returned state is rescored independently: the K2000 packed graph, or the G-set raw edge list.
  - `--corr-scale auto` folds MULTIBIT_SPEC's γ factor, mean_i Σ_j J_ij² / N, into λ (Onsager) or κ (TEC-T).
  - The default scale of 1 reproduces v6.4's tables exactly.
- **Image.** `board_mb2_e12_250mhz`: 12 engines, 4 per SLR (`scripts/pblocks_v6_e12.xdc`), REGSLICE=1, POST_ROUTE_OPT=1, integration `scripts/integrate_kernel_mb.tcl`. The build started 15:11 UTC.
  - The logic UUID and PDI hash are recorded in the image manifest at build time.
  - The bring-up requires the timing sign-off to be PASS.
  - **Fallback, fixed now.** If 250 MHz does not close, the same design is rebuilt at 225 MHz and measured under this protocol: absolute-time thresholds scale by 250/f, and R1/C1 are unaffected. If 225 MHz also fails, there is no board result for this protocol.
- **Verification before hardware** (xsim, Vivado 2025.1; vectors from `src/v6/gen_vectors_mb.cpp`; `scripts/run_sim_regression_mb.sh`). Results are in `results/multibit_20261007/sim/`.
  - **K2000 (±1).** The v6.4 short suite (10 trials, including resident tables) and full O1, O2, O3, O4, P1, T1 and X5 trials.
    - Spins, cut, flips and n_lin trace are bit-exact against the reference.
    - Expected outputs are byte-identical to v6.4's generator.
    - Per-trial cycles are identical to v6.4b on the same vectors (O1 41,300; O4 14,206; X5 12,409; all 10 short trials).
  - **G22 and G32 (ternary, N = 2000).** Short suite (10 trials) and full-length S = 600 Onsager and TEC-T trials: bit-exact.
  - **Random K = 4 matrix** (|J| ≤ 7, density 0.5; sim-only K = 4 build): short suite and S = 300/400 trials, bit-exact.
  - **Drain margin.** DRAIN = 1 fails (652/660 words), DRAIN = 2 passes, as for v6.4b.
  - **xsim workaround.** The simulator evaluates v6.4b's index `{g_k, m[2:0]}` (a concatenation containing a loop variable) as `m`. sca_core_mb.v writes `8 * g_k + m` instead; this is the same in Verilog and in synthesis.
    - A sim-only copy of v6.4b with only that rewrite passes the v6.4 vectors, which validates the flow.
    - v6.4b itself was verified with iverilog and on the board.

## Phases (`scripts/run_hw_mb.sh`; E = 12; seed 20261004; tables resident; one trial per launch in cohorts)

**Phase 1: exact gate on all engines.** Every trial must be bit-exact with the reference (spins, cut, flips, n_lin trace). Any failure stops the protocol.

| ID | Matrix | Configuration | Trials (per launch) | Trial ids |
|---|---|---|---|---|
| W1–W4 | K2000 | As PROTOCOL_HW_V6 (plain S40; Onsager λ0.7 + ramp S560; TEC −8 S40; O1 S960) | 12, 24 (2), 12, 12 | 900000+ (the v6.4 ids) |
| G22a | G22 | Onsager λ0.9 + ramp, q7, T 4 → 0.5, S600, scale auto | 24 (2) | 910000+ |
| G22b | G22 | TEC-T κ1.5 + ramp, q7, T 4 → 0.5, S600, scale auto | 12 | 910100+ |
| G32a | G32 | Onsager λ0.9 + ramp, q2, T 4 → 0.5, S600, scale auto | 24 (2) | 920000+ |
| G32b | G32 | TEC J_v −1, q2, T 4 → 0.5, S40 | 12 | 920100+ |

**Phase 2: K2000 regression.** X5 (TEC-T κ1.75 + ramp, q8, T15, S280) and O4 (Onsager λ0.9, q6, T12, S360).
- Each has 2,052 trials in 171 rounds of 12 engines, trial ids 0–2051 (round r, engine e → id 12r + e).
- Configuration strings, seed and launch shape are exactly those of the 12 × v6.4 run (`results/hwv6_board_v64_e12_250mhz`).

**Phase 3: G-set cohorts.** The configurations are in the frozen list `PROTOCOL_HW_MB_gset.txt`: 2,052 trials each, trial ids 0–2051, same seed and launch shape. Their selection is described below.

**Phase 4: board power** (`scripts/measure_power_mb.sh`; frozen list `PROTOCOL_HW_MB_power.txt`).
- Sustained 150 s loads on all 12 engines: 512 trials per launch, tables resident, no readback.
- 60 s idle phases bracket each load.
- Power is sampled from `/sys/class/hwmon/hwmon4` every 0.25 s. `ami_tool sensors` is never used: it stalls BAR access.

## Analysis (fixed now; `scripts/analyze_hw_mb.py <dir> results/hwv6_board_v64_e12_250mhz --predictions PROTOCOL_HW_MB_predictions.json`)

- **TTS99 primary.** t_R · ln(0.01)/ln(1 − P_round), with P_round measured over rounds and round 0 (coupling load) excluded from t_R.
- **TTS99 secondary.** t_R · ln(0.01)/(12 · ln(1 − p)).
- **Uncertainty.** Bootstrap 95% intervals over rounds.
- **Success.** The independently rescored cut is at least the configuration's target.
  - K2000: 33,000.
  - G-set: ⌈0.99 × best-known⌉, with best-known G22 = 13,359 and G32 = 1,410 (Benlic and Hao, *Breakout local search for the max-cut problem*, Eng. Appl. Artif. Intell. 26 (2013)). Targets: G22 ≥ 13,226; G32 ≥ 1,396.

## Predictions (pass conditions)

| ID | Prediction | Pass condition |
|---|---|---|
| E1 | Exact on hardware | Every phase-1 trial exact, on all 12 engines |
| R1 | Same algorithm as v6.4, bit for bit | Every X5 and O4 trial id has the same cut and flips as the 12 × v6.4 run (2 × 2,052 trials); W1–W4 likewise against v6.4's verify sets |
| C1 | Same cycles as v6.4 (identical in RTL simulation) | Mean per-launch cycle ratio (MB / v6.4) over common launches of X5 and O4 within 1 ± 0.01 |
| T1 | Same K2000 TTS | X5 primary TTS99 within ±2% of v6.4's 0.0503 ms × (250/f); O4 within ±5% of 0.0528 ms × (250/f) |
| M2 | Engines independent | For every cohort, P_round inside the interval implied by the per-trial p (as `analyze_hw_v6.py`) |
| G1 | Held-out success equals the reference model's | For each G-set configuration, \|p_board − p_validation\| ≤ 1.96 · sqrt(p̄(1 − p̄)(1/n_board + 1/n_val)) + 1/n_board |
| G2 | Cycle model | For each G-set configuration, measured mean round time within 1 ± 0.10 of the selection model t_R = (1.03 (max(22.6 S, 0.1424 F + 18.09 S) + 886) + 600) / f, evaluated at the measured mean flips F |
| G3 | G-set TTS | Measured primary TTS99 inside the band listed for its configuration below (validation p interval × round time ± 10%) |

Board power is reported without a pass condition; it is a first measurement for this image. The verdicts are computed mechanically by `analyze_hw_mb.py` from `PROTOCOL_HW_MB_predictions.json` (validation counts, predicted round times, G3 bands, tolerances).

## G-set configuration selection

**Why these configurations.** The software G-set study's selection (`research/gset_20261007/selected_configs.json`) did not exist when this protocol was frozen (checked 15:51 UTC). So the configurations are **pre-declared** here by a contingency procedure, `scripts/gset_select_mb.py`.
- The procedure uses only the reference model (`build/explore_mb`, which calls `sca_ref.hpp` `run_trial(..., dense)`).
- It ran before any multi-bit board data.
- Records are in `results/multibit_20261007/gset_select/`.
- If the study's selection arrives later, it can only be added by a further amendment, hashed before its board data.

**Procedure.**
1. **Pilot** (seed 0x5EED): 1,413 configurations, 48 trials each for G22 and 32 for G32.
   - Rules: plain SCA; Onsager online (λ, ramp, correction scale auto); TEC-T (κ, ramp, scale auto).
   - G22 grid: T0 ∈ {5, 6, 8, 10}, T1 ∈ {0.3, 0.5, 0.7}, q ∈ {5, 6, 7}, S ∈ {512, 1024, 2048, 4096}, λ ∈ {0.5, 0.9, 1.3}, κ ∈ {1, 1.5, 2}.
   - G32 grid: T0 ∈ {1.5, 2, 3}, T1 ∈ {0.3, 0.5, 0.7}, q ∈ {0.5, 1, 1.5}, S ∈ {1024, 2048, 4096}, λ ∈ {0.5, 0.9}, κ ∈ {1, 2}.
2. **Selection**, per instance and rule family: minimum predicted 12-engine TTS99 at the Wilson 95% lower bound of the pilot p, with at least one round, t_R · max(1, r99(P_round)); ties go to the shorter S.
   - Cycle model: c = max(22.6 S, 0.1424 F + 18.09 S) + 886, fitted on the MB RTL simulations. Round time t_R = (1.03 c + 600) / 250 MHz.
   - **Rule change, recorded.** The first run (rule v1) used the continuous form without the one-round floor. As P_round → 1 that form falls below one round, so its Wilson-bound estimate came out faster than its point estimate. The floor was added before validation, and the v1 output is kept as `selected_v1.json`.
3. **Validation**: the selected configurations re-run with a fresh seed, 0xBA11D, 512 trials each. This gives the unbiased p used below.
4. **Board**: seed 20261004, trial ids 0–2051, which are disjoint from the pilot and the validation.

**Frozen configurations** (`PROTOCOL_HW_MB_gset.txt`) **and predictions.** The validation p is given with its Wilson 95% interval. Predicted TTS is the primary form, t_R · r99(1 − (1 − p)^12), at 250 MHz. The G3 band runs from the upper p bound with 0.9 t_R to the lower p bound with 1.1 t_R (at least one round).

| Key | Schedule (S, T0 → T1, q, rule) | Validation p | t_R (ms) | Predicted TTS99 (ms) | G3 band (ms) |
|---|---|---|---|---|---|
| G22_ons | 1024, 5 → 0.5, q5, Onsager λ0.5 + ramp, scale auto | 0.430 [0.388, 0.473] (220/512) | 0.1014 | 0.069 | 0.055–0.112 |
| G22_plain | 1024, 5 → 0.7, q5, plain | 0.412 [0.370, 0.455] (211/512) | 0.1014 | 0.073 | 0.058–0.112 |
| G22_tecT | 1024, 5 → 0.7, q5, TEC-T κ1.5 + ramp, scale auto | 0.408 [0.367, 0.451] (209/512) | 0.1014 | 0.074 | 0.058–0.112 |
| G32_ons | 4096, 1.5 → 0.5, q1.5, Onsager λ0.9 + ramp, scale auto | 0.002 [0.000, 0.011] (1/512) | 0.387 | 76 | 12–474 |
| G32_plain | 4096, 2 → 0.5, q1, plain | 0.010 [0.004, 0.023] (5/512) | 0.552 | 21.6 | 8.3–55.6 |
| G32_tecT | 4096, 2 → 0.5, q1, TEC-T κ2 + ramp, scale auto | 0.012 [0.005, 0.025] (6/512) | 0.553 | 18.0 | 7.5–43.3 |

- **What the hardware result tests.** The board is bit-exact with the model (E1), so the G-set cohorts test the protocol (held-out p, G1), the cycle model (G2) and measured end-to-end time (G3). They do not test the arithmetic again.
- **G22 is cheap at this target** with 12 engines: one round of S = 1024.
- **G32 is hard** within the 4096-step limit: about 4–24 expected hits in 2,052 trials, so its TTS has wide intervals.
- **Power configurations** (`PROTOCOL_HW_MB_power.txt`): X5 (K2000), G22_ons and G32_tecT.

## Limits stated in advance

- **N is fixed at 2000.** The engine's N is a build parameter and the golden model's `sca::N` is 2000, so only N = 2000 G-set instances (G22–G42) run unpadded. G22 (random, +1 weights) and G32 (toroidal, ±1) are the spec's random and toroidal instances.
- **S ≤ 4096.** The step-table memory holds 4,096 entries. In the reference model, the coarse pilots needed hotter starts to reach the 99% target within that limit.
- **K = 4 and K = 8** are not built. K = 4 is simulation-verified; for K = 4 and K = 8, resource estimates come from out-of-context synthesis.
- **Single-copy banks halve the stored bits** (8 Mbit per engine at K = 2, against v6.4b's 16 Mbit). On the FPGA, however, the coupling store is bound by read ports, not capacity: K × 8 plane reads of 64 bits per group per cycle. The engine therefore keeps v6.4b's block-RAM count and doubles its UltraRAM.
- **No external review.**
