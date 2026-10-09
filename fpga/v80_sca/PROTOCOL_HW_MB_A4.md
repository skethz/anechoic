# PROTOCOL_HW_MB, Amendment 4: board validation of the bias path (K = 4 + bias image, revision 4)

Frozen on 8 October 2026, before any data from this image. At freezing time the image is not built yet: it is queued after the K = 2 image of Amendment 3. Hashes are in `PROTOCOL_HW_MB_A4.sha256`. Earlier protocols are unchanged.

## Purpose and scope

- **Purpose.** Functional board validation of the per-spin bias path and of K = 4 couplings: bit-exact against the shared golden model `src/sca_ref_bias.hpp` `run_trial_bias`.
- **Not in scope.** No solution-quality or TTS claims are made for graph partitioning or TSP. The software ReAIM study found that the engine rules do not reach the targets on these problems. The exported instances serve only as test inputs.
- **Own check of the export.** The export's TSP results reproduce exactly under the golden model. For all six TSP instances with the SCA schedule, the number of valid tours among trials 0–1023 equals the export's p_feasible × 1024: 993, 595, 974, 844, 158 and 678.
  - The export's GPP `precision_study` values are wrong. They hold the Max-Cut values of the same-named instances; see `research/reaim_benchmarks_20261008/ERRATUM_hw_export_precision_study.md`, with corrected values in `precision/precision_summary.json`.
  - Consistent with the corrected values (failure; for example G1 K = 4 SCA, p_feasible 0.077), the golden model gave no feasible GPP state in any checked trial (4 trials each of five schedules).
  - This amendment makes no comparison with GPP precision values. The GPP sets are dense-K = 4 exactness inputs only, with the update path at its maximum flip rate.

## Image and software

- **Core.** `src/v6/sca_core_mb_r4.v`: revision 3 with the schedule-table and trace memories parameterised (TLOG = 13, so S ≤ 8192; revision 3 allows 4096). The exported TSP schedules need S = 8192.
  - Parameters: K = 4, BIAS = 1, HBX = 17. HB = 17 covers 2047·7 plus any 16-bit bias; the lane high part is 19 bits.
  - Simulation sign-off before the build is in `build/sim_mb4_*`.
- **Front end.** `sca_mover_mb3.cpp` (variant mb4b; co-simulation PASS).
- **Host.** `src/host_sca_multi_mb4.cpp` with `src/mb_common4.hpp`, which loads SCAJINT8 matrices and int32 biases. The export's files are copied with SHA-256 to `data/a4_problems`; all 30 K = 4 files match the export's own hashes.
- **Build.** 6 engines, 2 per SLR. This is the most the URAM read ports allow: 258 URAM per engine, and the SLR pblocks hold 650–710 each. The build is queued at 250 MHz, with a 225 MHz rebuild if that fails timing (`scripts/chain_builds_mb3.sh`).
  - The bring-up uses whichever image passes timing, preferring 250 MHz. If neither passes, the board is not touched.

## Procedure (`scripts/bringup_mb4b.sh`, unattended)

1. **Preconditions.** Wait for the image. Wait until no Amendment-3 chain is running; every such chain ends with v6.4 restored. Check this amendment's hashes. No programming after 02:00 UTC on 9 October (the board stop is 04:00 UTC).
2. **Program.** Back up the staged PDI and program the image with the scoped command.
3. **Run the amendment** (`scripts/run_hw_mb4b.sh`), 6 engines, tables resident:
   - **Phase 1, exact gate.** Every trial is verified bit for bit against `run_trial_bias`: spins, score word (Σ s·h in bias mode, the cut otherwise), flips and n_lin trace. The sets are:
     - K2000 W1–W4 with bias mode off, with the trial ids of PROTOCOL_HW_V6;
     - K2000 W1 and W4 with bias mode on and b = 0;
     - random K = 4 matrices (|J| ≤ 7, density 0.5) with random biases in [−300, 300], at n = 777 and n = 2048, two schedules each;
     - `PROTOCOL_HW_MB_A4_verify.txt`, 6 trials each (one per engine): all 6 exported TSP instances (bias, n = 289–841, S = 8192), and GPP G1, G14 and G17 (dense K = 4, b = 0), each with its exported SCA schedule and seed.
   - **Phase 2, K2000 regression with bias mode off.** X5 and O4, 2052 trials each, trial ids 0–2051, seed 20261004.
   - **Phase 3, TSP cohorts.** `PROTOCOL_HW_MB_A4_cohorts.txt`: every exported TSP instance with its SCA schedule and the export's seed, 1026 trials (trial ids 0–1025).
4. **Power.** `scripts/measure_power_mb3.sh` with `PROTOCOL_HW_MB_A4_power.txt`: X5; bays29 (bias, n = 841, S = 8192); and G1 GPP, the maximum flip rate. 150 s per load.
5. **Restore and analyse.** Restore the v6.4 headline image regardless of the outcome. Then run `scripts/analyze_hw_mb_a4.py <dir> results/hwv6_board_v64_e12_250mhz PROTOCOL_HW_MB_A4_predictions.json`, followed by the power analysis.

## Predictions (pass conditions)

| ID | Prediction | Pass condition |
|---|---|---|
| A4-E | Bit-exact on every set | Every phase-1 trial is exact, padding bits are 0, and every set covers all 6 engines |
| A4-R | K2000 with bias mode off is identical to v6.4 | X5 and O4: the same cut and flips per trial id as the v6.4 12-engine run (2052 each); W1–W4 identical to v6.4's sets |
| A4-C | Same cycle behaviour as v6.4 | Launch-cycle ratio against v6.4 within 1% (X5, O4; launches after the coupling load) |
| A4-F | The engine reproduces the export trial for trial | Per TSP cohort: the number of valid tours (the c² one-hot spins form a permutation matrix) among trial ids 0–1023 equals the export's count exactly: gr17 993, gr21 595, gr24 974, fri26 844, bayg29 158, bays29 678 |

## Reporting

- Verdicts and phase-1 counts.
- K2000 X5/O4 TTS99 with 6 engines, against v6.4's 12.
- Valid-tour fractions and best energies, as functional by-products only.
- Board power per load.
- Resources and timing of the K = 4 + bias image against the K = 2 image.

## Limits

- **Quantized TSP.** Tour quality at K = 4 is limited by the quantization; see the ReAIM study's summary. Nothing here measures solution quality.
- **Collapsing GPP sets.** The GPP schedules collapse: all spins flip every step. They test the update path at its maximum flip rate, not the problem.
