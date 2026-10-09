> **Note (9 October 2026).** Some values here come from revision 1 of the K = 2 image (for example 124 W and 88.6% LUT); the paper uses revision 3, which has no bias. The G-set counts in the Amendment 3 section (48 and 41 of 51) use the pre-registered grids. The paper reports these and, in Table 8, the extended grids of Amendment 1 (49 and 42).

# Multi-bit couplings on the V80 engine (K = 2): design, verification and board results

**Status (8 October 2026, 19:00 UTC).**
- **Revision-1 board results.** The 12 × K = 2 image is measured: main protocol, Amendment 1 and power, all below.
- **Revision 3 (runtime n, optional per-spin bias).** Simulation-verified: 20 of 20 xsim runs and both front-end co-simulations pass.
- **K = 2 image.** The 12-engine revision-3 K = 2 image (`board_mb2n_e12_250mhz`) is measured: Amendment 3, all verdicts PASS, below. v6.4 was restored at 09:08 UTC.
- **K = 4 + bias image.** Six engines failed routing at 250 MHz and passed timing at 225 MHz. Amendment 4 measured them, all verdicts PASS, and v6.4 was restored at 18:49 UTC.

**Specification.** `MULTIBIT_SPEC.md`. **Pre-registered protocol.** `PROTOCOL_HW_MB.md`, frozen 15:55 UTC, before any board data; hashes in `PROTOCOL_HW_MB.sha256`.

**Records.** Local: `results/multibit_20261007/`. FPGA host: `/scratch/USER/sca_v80_20261003`.

**Board state.** At 22:27 UTC on 7 October the v6.4 headline image was restored by the unattended chain: UUID 8669ad38… matched, BAR2 OK, staged PDI SHA-256 67e591f3…. Every later chain restores v6.4 when it finishes.

## What was built

| File | Role |
|---|---|
| `src/v6/sca_core_mb.v` | Multi-bit core, parameter K (built at K = 2), derived from v6.4b (`sca_core.v`, SHA-256 32e1c2f8…, unchanged) |
| `src/v6/sca_mover_mb.cpp`, `tb_mover_mb.cpp` | HLS front end with MB_K planes in the coupling stream (same register map and timers) |
| `src/mb_common.hpp` | Shared helpers: integer matrices (K2000 / G-set edge list / random), bit planes, correction scaling, range checks, independent scorer |
| `src/host_sca_multi_mb.cpp` | Board host: `--matrix`, `--target`, `--corr-scale`, new HBM layout, `--verify` / `--cohort` / `--power-run` |
| `src/v6/gen_vectors_mb.cpp`, `tb_core_mb.v` | Stream vectors from `sca_ref.hpp` `run_trial(..., dense)` (unchanged golden model) and the xsim testbench |
| `src/v6/explore_mb.cpp`, `scripts/gset_select_mb.py` | Reference-model schedule exploration and the pre-declared G-set selection |
| `scripts/sim_mb.sh`, `run_sim_regression_mb.sh`, `ooc_mb.tcl`, `integrate_kernel_mb.tcl`, `run_hw_mb.sh`, `bringup_mb2.sh`, `measure_power_mb.sh`, `analyze_hw_mb.py`, `analyze_power_mb.py` | Flow |

### Arithmetic (as `MULTIBIT_SPEC.md`)

- **Planes.** Plane 0 is the sign (1 iff J < 0). Planes 1..K−1 hold bits 0..K−2 of |J|. The layout is the existing row/bank layout.
- **Field update.** For the up to 8 flips (x_e, σ_e) of a cycle, with τ_e = σ_e (−1)^sg:
  - per magnitude plane b, a_b = #{e: m_b ∧ τ_e = +1} and c_b = #{e: m_b};
  - the field changes by Σ_b 2^b · mode · (2 a_b − c_b), with mode 2 in a decision step and 1 at initialisation.
  - This takes two 8-input population counts per plane and no multipliers. K = 1 reduces to v6.4b's 2·pc − nvalid.
- **Widths.** Fields are 12 / 15 / 19 bits at K = 2 / 4 / 8, and the field increment is K + 4 bits.
  - The lane z keeps v6.4b's 16-bit high part for K ≤ 4. This is exact while |q| + |corr| < 2^29, which the host checks. At K = 8 the high part is 21 bits.
- **Diagonal.** J_xx = 0 as in the spec, so the stored field is the true field. v6.4b stored h − s via a −1 diagonal and added the 1 back in the lane; the decision values are the same.
- **Correction scaling.** Done on the host only: `--corr-scale auto` multiplies λ (Onsager) or κ (TEC-T) by mean_i Σ_j J_ij² / N, the γ factor of the spec. Scale 1, the default, gives v6.4's tables bit for bit.
  - K2000 itself would get 1999/2000 under the spec formula, so the K2000 regression uses scale 1.

### Residue banks

- **The residue property, checked.** In v6.4b, extractor e emits `o_x <= {curk, enc32(lb), E[2:0]}`, so x mod 8 = e always, and coupling port e is driven only by extractor e. Every simulation also asserts this at run time: 0 violations.
- **Layout.** Each port has one single-port bank holding rows x ≡ e (mod 8) at address x >> 3. The store is K · 2048² bits per engine: 8 Mbit at K = 2, against v6.4b's four 1-bit copies at 16 Mbit.
- **FPGA finding.** On the V80 the coupling store is bound by read ports, not by bits.
  - Each group needs 8 · K plane reads of 64 bits per cycle. A BRAM36 gives one such read, and a URAM gives two.
  - At K = 2 a group therefore needs 16 read ports: 8 RAMB36 (residues 0–3, one 256 × 64 bank per plane) plus 4 URAM (residues 4–7, one URAM per plane pair).
  - Per engine that is 258 RAMB36 (as v6.4b) and 130 URAM (v6.4b: 66), confirmed by out-of-context synthesis.
  - So the single-copy banks keep the BRAM count flat while the read width doubles. They do not halve BRAM at K = 2.
  - With 12 engines, URAM use is about 81% of the device (650–710 per SLR pblock, 520 used).

## Simulation verification (xsim, Vivado 2025.1, on fpga-host)

All runs are bit-exact against `sca_ref.hpp` `run_trial(..., dense)`: final spins, cut, flips and the per-step n_lin trace. Random input gaps and output back-pressure were applied. Records are in `results/multibit_20261007/sim/simpack_20261007a/` (first core) and `simpack_20261007b/` (r2 regression and negative-score suite).

| Run | Matrix | Cases | Result |
|---|---|---|---|
| k2000_short_k2 | K2000 ±1 | v6.4 short suite: plain, Onsager + ramp, TEC, TEC with resident tables (10 trials, 660 words) | PASS |
| k2000_proto_c0..c6_k2 | K2000 ±1 | Full O1, O2, O3, O4, P1, T1, X5 | PASS (7/7) |
| g22_short_k2, g22_long_k2 | G22 (random, +1) | Short suite; S = 600 Onsager and TEC-T (correction scale auto) | PASS |
| g32_short_k2, g32_long_k2 | G32 (toroidal, ±1) | Same | PASS |
| rand4_short_k4, rand4_long_k4 | Random K = 4 (\|J\| ≤ 7, density 0.5, seed 20261007) | Short suite; S = 400 Onsager and S = 300 TEC-T | PASS |
| k2000_short_k2_drain1 | K2000 | DRAIN = 1 | FAIL, as expected (DRAIN = 2 is the minimum, as for v6.4b) |
| rand8_short_k8 (extra) | Random K = 8 (\|J\| ≤ 127, density 0.01, seed 20261008) | Short suite (21-bit lane high part) | PASS |

- **K2000 regression.** The expected outputs from the multi-bit generator are byte-identical to v6.4's generator (short suite and O4).
  - v6.4b, simulated on the same cases (K = 1 vectors), has identical per-trial cycles: all 10 short-suite trials, O1 41,300, O4 14,206 and X5 12,409.
  - The only cycle difference is the one-off coupling load: K × 32,768 beats.
- **Cycle model.** Since the cycles equal v6.4b's, v6.4's linear model (0.1424·flips + 18.09·S + 886) holds for K2000-like flip rates.
  - At the low flip rates of sparse G-set runs, about 2 flips per step, a per-step floor of about 22.6 cycles dominates. The linear model then underestimates by about 21% (G22, S = 600).
  - The selection model c = max(22.6 S, 0.1424 F + 18.09 S) + 886 fits the simulated full-length G-set trials within 1% (G22 1.000; G32 0.993–0.998).
  - It overestimates the dense K = 4 trials by about 5%.
- **xsim pitfall.** Vivado 2025.1's simulator evaluates an array index built by concatenation with a loop variable, such as v6.4b's `h[{g_k, m[2:0]}]`, as `h[m]`, at every optimisation level.
  - Unmodified v6.4b therefore fails under xsim, although it is bit-exact with iverilog and on the board.
  - sca_core_mb.v uses `8 * g_k + m`, which is identical in Verilog semantics and in synthesis.
  - A sim-only copy of v6.4b with just that rewrite passes all 660 v6.4 words, which validates the xsim flow.
  - The testbench also avoids printing a string-valued conditional, which xsim prints as a number.

## Resources (out-of-context synthesis, one engine, xcv80)

| | v6.4b (synth / routed) | MB K = 2 (synth) | MB K = 4 (synth) | MB K = 8 (synth) |
|---|---:|---:|---:|---:|
| LUT | 171,376 / 149,113 | 203,287 | 369,105 | 785,967 |
| FF | | 138,633 | 167,254 | 220,809 |
| RAMB36 | 258 | 258 | 514 | 1,026 |
| URAM | 66 | 130 | 258 | 514 |
| DSP | 282 | 282 | 282 | 282 |
| Engines that fit the V80 | 12 (built) | 12 (built) | 6 (BRAM- and LUT-bound) | 3 |

The +32k LUT per engine (+8 per field) comes from the second population count per field. Projected at 12 engines: about 87% LUT (v6.4: 74%), so timing closure at 250 MHz is less certain than for v6.4.

## Board image (12 × K = 2 at 250 MHz)

- **Build.** `board_mb2_e12_250mhz`: 12 engines, 4 per SLR, the same pblocks, register slices and post-route phys_opt as v6.4. Started 15:11 UTC.
- **Post-placement** (all 12 engines):

| Resource | Used | Device | Per SLR | v6.4 |
|---|---:|---:|---:|---:|
| LUT | 2,279,831 | 88.6% | 86.3 / 90.5 / 89.1% | 73.9% |
| FF | 1,857,475 | 36.1% | | |
| RAMB36 | 3,318 | 88.7% | | same as v6.4 |
| URAM | 1,560 | 81.0% | 76.8 / 83.3 / 83.3% | 41.1% |
| DSP | 3,384 | 31.2% | | |

- **Timing.**
  - Post-placement estimate: WNS +0.331 ns (v6.4: +0.311 ns).
  - Routing took 77 minutes at congestion level 5; peak memory was 60 GB. After the first iteration the setup WNS was −0.221 ns (TNS −0.959 ns).
  - Routed: setup **+0.007 ns**, hold +0.010 ns, 0 failing endpoints.
  - **Final sign-off** (after post-route phys_opt, which changed nothing): PASS, setup +0.007 ns, hold +0.010 ns, pulse width 0.000 ns, 0 failing endpoints.
- **Routed utilisation per engine** (`results/util_mb_board_mb2_e12_250mhz/mb2_hier.rpt`, post-route phys_opt):

| Per engine | v6.4 (routed) | MB K = 2 (routed) |
|---|---:|---:|
| LUT | 149,512 | 180,502–181,853 (+21%) |
| FF | 131,555 | 138,818–139,890 |
| RAMB36 | 258 | 258 |
| URAM | 66 | 130 |
| DSP | 282 | 282 |

  The whole design uses 2,289,586 LUT (88.9%), 1,857,475 FF, 3,312 RAMB36 (+12 RAMB18), 1,560 URAM and 3,384 DSP.
- **Image.** Logic UUID `257f202a432b4f4914f0c933fe037434`; PDI SHA-256 `331d3349a82dfaf7e16aaca598a62c1edd53e083f9d67f07c09835e460301640` (111,366,480 bytes).
  - Kernel clock 249.99999975 MHz (nominal, from the MMCM dividers).
  - Pipeline finished 20:37 UTC: 5 h 26 min from start.

## Board run, main protocol (`PROTOCOL_HW_MB.md`, frozen 15:55 UTC)

**Image.** Programmed 20:37–21:28 UTC by `scripts/bringup_mb2_a1.sh`, which wraps the frozen `bringup_mb2.sh`.
- The staged v6.4 headline image was backed up first (SHA-256 67e591f3…, `results/staged_image_backup_20261007T203744Z`).
- UUID matched and BAR2 was OK; no rescan was needed.
- All phases ran in 13 s.
- Records: `results/multibit_20261007/board/board_mb2_e12_250mhz/` (`hwmb_board_mb2_e12_250mhz/analysis.txt`, `mb_summary.json`).

| ID | Result | Verdict |
|---|---|---|
| E1 | Exact gate 132/132 bit-exact on all 12 engines: K2000 W1–W4 (60 trials, the v6.4 ids), G22a/b (36), G32a/b (36) | **PASS** |
| R1 | X5 and O4 each 2,052/2,052 trial ids with identical cut and flips to the 12 × v6.4 run; W1–W4 identical too | **PASS** |
| C1 | Per-launch cycles over common launches, MB / v6.4: X5 1.00024, O4 1.00022; identical in 1,020 of 2,040 launches (the rest is NoC jitter) | **PASS** |
| T1 | X5 primary TTS99 0.0503 ms (v6.4: 0.0503); O4 0.0529 ms (v6.4: 0.0528) | **PASS** |
| M2 | Engine independence holds for every cohort | **PASS** |
| G1 | All 6 G-set cohorts' p agree with the reference-model validation | **PASS** |
| G2 | Measured round time against the selection cycle model: 5 of 6 within ±10%; G22_ons 1.20 | **FAIL** |
| G3 | Measured TTS inside its band: 5 of 6; G22_ons 0.1218 ms against a band of 0.055–0.112 ms | **FAIL** |

**Why G2 and G3 fail for G22_ons.** In 3 of its 2,052 trials the run fell into the all-equal period-2 state: about 2.0M flips, each such trial costing 266k cycles against a median round of 26.2k. Those three rounds raise the mean round time from about 0.105 to 0.1218 ms. The cycle model, evaluated at the mean flips, cannot represent that tail.

**Contingency G-set results** (pre-declared configurations; target ⌈0.99 BKV⌉; 12 engines at 250 MHz):

| Key | S | p (board) | p (validation) | t_round (ms) | TTS99 primary (ms) [95%] | TTS99 secondary (ms) |
|---|---:|---:|---:|---:|---:|---:|
| G22_ons | 1024 | 0.451 | 0.430 | 0.1218 | 0.1218 [0.105, 0.144] | 0.078 |
| G22_plain | 1024 | 0.375 | 0.412 | 0.1033 | 0.1033 | 0.084 |
| G22_tecT | 1024 | 0.424 | 0.408 | 0.1036 | 0.1036 | 0.072 |
| G32_ons | 4096 | 0.006 | 0.002 | 0.403 | 23.5 [14.9, 51.9] | 24.3 |
| G32_plain | 4096 | 0.015 | 0.010 | 0.530 | 13.6 [9.7, 20.7] | 13.8 |
| G32_tecT | 4096 | 0.014 | 0.012 | 0.532 | 13.7 [10.0, 20.8] | 14.3 |

### Finding: the device cut word is wrong when Σ s·h < 0, a latent bug inherited from v6.4b

The 3 oscillating G22_ons trials were flagged as integrity failures: the device cut word was 2^62, while the host's independent rescoring gave 0.
- **Spins and flips are exact.** They match `run_trial` on the CPU for all three trials.
- **Cause.** v6.4b's halving `(x + {63'd0, x[63]}) >>> 1` mixes signed and unsigned operands, so `>>>` shifts logically. The result is wrong exactly when Σ s·h < 0 (cut < sumw/2).
  - With Σ s·h = −39,980 this gives cut_half = 0x7FFF…B1EA, cut_sum = 2^63 and cut = 2^62, as observed.
- **Impact.**
  - K2000 (sumw = −1040): impossible, so v6.4b and every K2000 result are unaffected.
  - +1-weight G-set graphs: only failed low-cut states are affected.
  - Success always uses the host rescoring, so p and TTS are unaffected.
  - The integrity flags were 3 trials here, and 4, 2 and 177 in Amendment 1's G35/G37/G38 TEC cohorts. All 177 G38 flags have cut ≤ 5,616 < sumw/2 = 5,889.5.
- **Simulation.** It was missed because every simulated state had Σ s·h > 0. A new negative-score suite (G22 with q = 1; 8 trials ending at cut 0) reproduces the board exactly on `sca_core_mb.v`: 8 errors, all cut words.
- **Fix.** `src/v6/sca_core_mb_r2.v` makes the sign operand `$signed`. It passes that suite and the regression with identical cycles.
- **Re-validation on the board.** A revision-2 image is pre-registered as Amendment 2; its result is below.

## Amendment 1: the G-set study's 12-engine selections (`PROTOCOL_HW_MB_A1.md`, frozen 18:16 UTC)

The study's `selected_configs.json` (SHA-256 291401004efe27e4…) arrived after the main protocol was frozen. It was added as a hashed amendment before any board data.
- **Scope.** All 21 N = 2000 instances (G22–G42); the four engine-implemented rules; the `tts_E12` selection in both grid variants.
- **Size.** 100 cohorts × 2,052 trials, plus a 12-trial exact gate per instance.
- **Run time.** It ran in 82 s.

| ID | Result | Verdict |
|---|---|---|
| A1-E | 21 × 12/12 exact-gate trials bit-exact | **PASS** |
| A1-P | 100/100 cohorts' board p agree with the study's float-model p (within 1.96σ + 0.02) | **PASS** |
| A1-T | Median board TTS99 / study `tts_E12` = 1.000 (range 0.27–1.87) | **PASS** |

**Measured 12-engine TTS99 (ms, primary)**, original-grid selections. In brackets: board p at the target ⌈0.99 BKV⌉. "Best" is the post-hoc minimum over all 4 rules and both grid variants.

| Instance | Target | SCA | TEC | Onsager-kT | Onsager-online | Best (rule, grid) |
|---|---:|---:|---:|---:|---:|---:|
| G22 | 13226 | 0.080 (0.66) | 0.038 (0.28) | 0.037 (0.33) | 0.064 (0.85) | 0.037 (tecT_o) |
| G23 | 13211 | 0.069 (0.21) | 0.038 (0.71) | 0.036 (0.63) | 0.035 (0.39) | 0.035 (ons_o) |
| G24 | 13204 | 0.057 (0.87) | 0.034 (0.58) | 0.044 (0.56) | 0.031 (0.35) | 0.031 (ons_o) |
| G25 | 13207 | 0.063 (0.16) | 0.038 (0.67) | 0.036 (0.62) | 0.037 (0.26) | 0.036 (tecT_o) |
| G26 | 13195 | 0.057 (0.91) | 0.034 (0.64) | 0.036 (0.62) | 0.035 (0.39) | 0.034 (tec_o) |
| G27 | 3308 | 0.237 (0.55) | 0.141 (0.23) | 0.099 (0.45) | 0.137 (0.51) | 0.099 (tecT_o) |
| G28 | 3266 | 0.240 (0.70) | 0.115 (0.33) | 0.105 (0.20) | 0.139 (0.54) | 0.105 (tecT_o) |
| G29 | 3371 | 0.238 (0.38) | 0.219 (0.44) | 0.088 (0.33) | 0.207 (0.24) | 0.088 (tecT_o) |
| G30 | 3379 | 0.201 (0.20) | 0.115 (0.35) | 0.101 (0.18) | 0.136 (0.46) | 0.101 (tecT_o) |
| G31 | 3277 | 0.260 (0.17) | 0.103 (0.31) | 0.099 (0.57) | 0.140 (0.48) | 0.099 (tecT_o) |
| G32 | 1396 | 20.706 (0.01) | 34.877 (0.00) | 2.491 (0.07) | 17.859 (0.01) | 2.491 (tecT_o) |
| G33 | 1369 | 22.820 (0.01) | 39.354 (0.01) | 3.825 (0.05) | 14.302 (0.02) | 3.825 (tecT_o) |
| G34 | 1371 | 19.952 (0.01) | 13.824 (0.02) | 2.570 (0.07) | 7.449 (0.04) | 2.570 (tecT_o) |
| G35 | 7611 | 0.454 (0.75) | 0.321 (0.67) | 0.199 (0.28) | 0.587 (0.63) | 0.092 (tecT_x) |
| G36 | 7604 | 0.396 (0.18) | 0.250 (0.69) | 0.189 (0.25) | 0.469 (0.59) | 0.092 (tecT_x) |
| G37 | 7615 | 0.454 (0.71) | 0.314 (0.65) | 0.333 (0.64) | 0.993 (0.10) | 0.134 (tecT_x) |
| G38 | 7612 | 0.453 (0.76) | 0.249 (0.24) | 0.176 (0.28) | 0.549 (0.59) | 0.092 (tecT_x) |
| G39 | 2384 | 2.349 (0.14) | 1.192 (0.19) | 1.778 (0.22) | 1.543 (0.06) | 1.192 (tec_o) |
| G40 | 2376 | 2.042 (0.12) | 2.149 (0.13) | 2.451 (0.07) | 1.836 (0.11) | 1.836 (ons_o) |
| G41 | 2381 | 3.262 (0.08) | 3.307 (0.05) | 1.496 (0.09) | 1.981 (0.09) | 1.496 (tecT_o) |
| G42 | 2457 | 2.303 (0.06) | 2.296 (0.13) | 2.634 (0.17) | 1.476 (0.12) | 1.476 (ons_o) |

Geometric mean over the 21 instances (original grid): SCA 0.578 ms, TEC 0.407 ms, **Onsager-kT 0.263 ms**, Onsager-online 0.425 ms.

## Board power (12 × K = 2 at 250 MHz; hwmon4 every 0.25 s; 150 s loads, 512 trials per launch, devices 99.4–99.8% busy, host/device time 1.000–1.001)

| Phase | Board W | VCCINT W | Dynamic W (vs bracketing idle) | Energy to solution, board |
|---|---:|---:|---:|---:|
| Idle (before; between loads, warming) | 73.1; 81.8, 77.7, 77.2 | 45.1–51.4 | | |
| X5 (K2000) | 124.3 ± 3.5 | 86.4 | 46.9 | 6.26 mJ (0.0503 ms); 512 µJ per trial |
| G22_ons | 95.3 ± 1.1 | 62.6 | 15.6 | 11.6 mJ (0.1218 ms) |
| G32_tecT | 92.6 ± 0.4 | 61.0 | 15.2 | 1.27 J (13.69 ms) |

For comparison, v6.4 under X5 drew 108.2 W (35.8–40.8 W dynamic), giving 5.45 mJ to solution and 446 µJ per trial. The multi-bit engine draws about 16 W more at K2000 flip rates: twice the coupling read width and URAM count, plus the second population count. At sparse G-set flip rates the dynamic power is about a third of that.

## Revision-2 image (Amendment 2)

- **What happened.** The revision-2 image (`sca_core_mb_r2.v`, cut-word fix only) started building at 21:43 UTC on 7 October. Its unattended board validation (Amendment 2) was cancelled at 01:06 UTC on 8 October, before any board action.
- **Why.** Revision 3, below, supersedes revision 2: it contains the same fix and adds runtime n.
- **The build.** It was stopped at 01:17 UTC, during routing. The intermediate WNS was −0.399 ns; revision 1 passed through −0.221 ns at the same stage before closing at +0.007 ns. The memory was needed for the revision-3 simulations and build.
- **What this leaves.** Revision 2 has no board data. Its fix is board-validated through revision 3 (Amendment 3: the negative-score sets in the exact gate).

## Revision 3: runtime n and per-spin bias (`src/v6/sca_core_mb_r3.v`)

Golden model: the shared `src/sca_ref_bias.hpp` `run_trial_bias(n, dense, bias, tables, seed, trial, trace)`, unchanged. Its own test (`test_ref_bias`) passes 16 of 16 checks on fpga-host.

- **Runtime n.** Header H0[47:36] = n, 1..2048; 0 means the build-time N = 2000.
  - Each group compares its lane's spin index with a registered copy of n. Spins i ≥ n never flip and do not count in n_lin or in flips.
  - Output bits of spins i ≥ n are masked to 0. The random-number streams advance as before.
  - The cycles do not depend on n, because the engine still sweeps 2048 slots.
- **Bias** (parameter `BIAS = 1`). Each group keeps its 64 biases as 16-bit registers, 1,024 FF per group.
  - **Loading.** In bias mode (H0[35]), a 256-beat bias segment follows the coupling planes in the load stream. Beat i holds the 8 biases of group i >> 3, fields 8(i & 7) + q.
  - **Trial start.** The trial-start clear loads b into the field, so the INIT batch builds h(0) = b + J s(0). The per-flip update is unchanged.
  - **Score word.** Output word 32 becomes Σ s_i h_i, and the host computes H = −(Σ s h + Σ b s)/2.
  - **Width.** `HBX` overrides the field width. HB = 17 covers 2047·7 + 32767 at K = 4 with any 16-bit bias. The lane high part is then HB + 2 bits.
  - **Load time.** The one-off load grows by 256 beats: 0.8% at K = 2 and 0.4% at K = 4.
- **Front end.** `src/v6/sca_mover_mb3.cpp` passes n and the bias flag through the existing `load_J` argument: bit 2 is the bias flag, bits 27..16 are n. The register map is unchanged.
  - HLS C simulation and C/RTL co-simulation PASS for both variants, mb2n (K = 2) and mb4b (K = 4). The 8 cases include n = 800, 1000, 17 and 2048, with and without bias.
- **Host.** `src/host_sca_multi_mb3.cpp` with `src/mb_common3.hpp`, a copy of the frozen `mb_common.hpp` with runtime n, a text Ising format, the bias segment and an independent Σ s h / Σ b s scorer. The host:
  - writes n on every launch;
  - checks the field bound max Σ|J| + max|b| < 2^(HB−1), and that padding bits are 0;
  - rescores every result independently;
  - verifies traced trials against `run_trial_bias`.

**xsim regression** (`build/sim_mb3_20261008b`; records in `results/multibit_20261007/sim/simpack_20261008_r3/`). 20 of 20 runs pass bit-exact (spins, score word, flips, n_lin trace), with 0 residue violations. The first launch, `sim_mb3_20261008a`, failed to compile because a wire was used before its declaration; that was fixed and the suite rerun.

| Run | Image parameters | Matrix (n) | Result |
|---|---|---|---|
| k2_k2000_short_n2000, _n0 | K = 2, BIAS = 0 | K2000 (2000; header n = 2000, or 0 = default) | PASS. Expected words byte-identical to revision 2 / v6.4; with n = 0 the stimulus is byte-identical to revision 2's; identical per-trial cycles (suite end at cycle 105,956, as revision 2) |
| k2_k2000_proto | K = 2 | K2000 O4 and X5 | PASS (O4 trial 14,206 cycles, as v6.4b) |
| k2_g22_neg | K = 2 | G22, q = 1 (negative Σ s h) | PASS (cut-word fix) |
| k2_G1/G11/G14_short | K = 2 | G1 random, G11 toroidal, G14 planar (n = 800) | PASS |
| k2_G1/G11/G14_study, k2_G43/G51_study | K = 2 | The study's E12 configuration of every rule, one traced trial each (n = 800, 1000) | PASS |
| k4b_k2000_short_nobias, _bias0 | K = 4, BIAS = 1, HB = 17 | K2000 with bias mode off, and on with b = 0 | PASS (the same expected words as v6.4) |
| k4b_rand4_n777, _n2048 | K = 4 + bias | Random \|J\| ≤ 7, density 0.5, biases uniform in [−300, 300] (n = 777, 2048) | PASS |
| k4b_gp20 | K = 4 + bias | Graph partition, n = 20 (A = 2, B = 2; \|J\| ≤ 4, b = 16) | PASS |
| k4b_tsp16, k4b_tsp25 | K = 4 + bias | TSP 4 and 5 cities (n = 16, 25; \|J\| ≤ 6, \|b\| ≤ 54) | PASS |
| k4b_g1_short | K = 4 + bias | G1 Max-Cut (n = 800), bias mode off | PASS |

- **The GP and TSP instances** were generated by `src/v6/gen_gp_tsp.cpp` and checked by brute force over all 2^n states. In each, the Ising ground state is a feasible partition or tour with the optimal cut or length:
  - graph partition: optimum cut 21;
  - TSP: tour lengths 7 and 6.
- **Ground state reached in simulation.** The engine reached the TSP-16 ground state (H = −252) in several simulated trials.
- **Instance files.** In `simpack_20261008_r3/*.ising`.

## Amendment 3: revision-3 K = 2 image, runtime n, all 51 G-set instances (`PROTOCOL_HW_MB_A3.md`, frozen 01:35 UTC on 8 October)

**Image `board_mb2n_e12_250mhz`.** 12 engines of `sca_core_mb_r3.v` at K = 2 without bias, 250 MHz. This is the pre-registered clock; no fallback was needed.
- **Routed timing:** setup slack +0.024 ns, hold +0.010 ns, 0 failing endpoints; signed off as PASS.
- **Utilisation:** 2,302,342 LUT (89.4%), 1,858,655 FF (36.1%), 3,312 RAMB36, 1,560 URAM (81.0%), 3,384 DSP. Revision 1 used 2,289,586 LUT and 1,857,475 FF.
- **Image identity:** logic UUID 14f910d903d81b622198e745a92a49af; PDI SHA-256 a34ed2ec….
- **Timeline.**
  - The pipeline started 01:25 UTC and finished 07:07 UTC.
  - Programming ran 07:07–07:59. The staged v6.4 PDI was backed up first: `staged_image_backup_20261008T070730Z`, SHA-256 67e591f3….
  - The amendment ran 07:59–08:05 and power 08:05–08:23.
  - v6.4 was restored 08:23–09:08: UUID matched, BAR2 OK.
- **Records:** `results/multibit_20261007/board/board_mb2n_e12_250mhz/`. They contain `hwmb_board_mb2n_e12_250mhz_a3/analysis.txt`, `a3_summary.json`, `a3_report.txt` and `a3_tables.md`. The cohort records are in `gset_cohorts_jsonl.tgz`; the spin states stay on fpga-host, with hashes in `gset_spins_on_tempo.sha256`.

| ID | Result | Verdict |
|---|---|---|
| A3-E | 61 sets, 768/768 trials bit-exact against `run_trial_bias`, all 12 engines, padding bits 0. The sets: K2000 W1–W4, G22/G32 sets, the two negative-score sets, and one set per instance for all 51 instances (n = 800, 1000, 2000) | **PASS** |
| A3-R | X5 and O4 each 2,052/2,052 trial ids with the same cut and flips as the 12 × v6.4 run; W1–W4 identical | **PASS** |
| A3-C | Launch cycles against v6.4: X5 1.00036, O4 1.00032. 1,019 and 1,020 of 2,040 launches are cycle-identical; the rest is NoC jitter | **PASS** |
| A3-I | All 100 N = 2000 cohorts (21 instances) identical to revision 1 (Amendment 1): 205,200 trials compared, 0 mismatches in cut or flips | **PASS** |
| A3-P | 264/264 cohorts' board p agree with the study's p | **PASS** |
| A3-T | Median board TTS99 / study `tts_E12` = 0.986 (range 0.26–1.48) | **PASS** |

- **K2000 TTS99:** X5 0.0504 ms (v6.4: 0.0503), O4 0.0529 ms (v6.4: 0.0528).
- **Integrity:** 0 device/host score disagreements in 1,625,184 cohort trials. The cut-word fix holds.
- **M2** (round independence, reported only): holds for 262 of 264 cohorts.

**Per-class results** (12 engines, primary TTS99, target ⌈0.99 BKV⌉, 6,156 trials per cohort; classes from the G-set weight signs):

**Geometric-mean primary TTS99 (ms), pre-registered (original) grid, 12 engines.** In brackets: the number of instances.

| Rule | random +1 | random ±1 | toroidal ±1 | planar +1 | planar ±1 | all |
|---|---:|---:|---:|---:|---:|---:|
| SCA | 0.055 (15) | 0.121 (10) | 6.070 (6) | 0.249 (12) | 0.767 (8) | 0.241 (51) |
| TEC | 0.036 (15) | 0.096 (10) | 6.119 (6) | 0.203 (12) | 0.654 (8) | 0.188 (51) |
| Onsager-kT | 0.029 (15) | 0.054 (10) | 1.409 (6) | 0.162 (12) | 0.489 (8) | 0.121 (51) |
| Onsager-online | 0.040 (15) | 0.072 (10) | 4.311 (6) | 0.364 (12) | 0.440 (8) | 0.191 (51) |

**Geometric-mean primary TTS99 (ms), extended grid (extended_amendment1), 12 engines.** In brackets: the number of instances.

| Rule | random +1 | random ±1 | toroidal ±1 | planar +1 | planar ±1 | all |
|---|---:|---:|---:|---:|---:|---:|
| SCA | 0.053 (15) | 0.119 (10) | 7.409 (6) | 0.249 (12) | 0.707 (8) | 0.240 (51) |
| TEC | 0.030 (15) | 0.102 (10) | 6.216 (6) | 0.197 (12) | 0.774 (8) | 0.186 (51) |
| Onsager-kT | 0.031 (15) | 0.054 (10) | 1.624 (6) | 0.092 (12) | 0.508 (8) | 0.111 (51) |
| Onsager-online | 0.040 (15) | 0.072 (10) | 5.114 (6) | 0.306 (12) | 0.488 (8) | 0.190 (51) |

Onsager-kT TTS99 below plain SCA on 48 of 51 instances, below TEC on 41, below both on 41 (original grid).
  not below SCA: ['G15', 'G40', 'G47']
  not below TEC: ['G14', 'G15', 'G22', 'G24', 'G26', 'G37', 'G39', 'G40', 'G47', 'G53']

**Per instance, pre-registered grid** (TTS99 in ms; board p in brackets). "Best" is the post-hoc minimum over the 4 rules and both grids; "Study best E12" is the study's own minimum `tts_E12`.

| Instance | N | Class | Target | SCA | TEC | Onsager-kT | Onsager-online | Best (key) | Study best E12 |
|---|---:|---|---:|---:|---:|---:|---:|---:|---:|
| G1 | 800 | random | 11508 | 0.053 (0.67) | 0.063 (0.77) | 0.019 (0.46) | 0.068 (0.69) | 0.019 (G1_tecT_o) | 0.023 |
| G2 | 800 | random | 11504 | 0.065 (0.15) | 0.037 (0.24) | 0.026 (0.62) | 0.064 (0.15) | 0.026 (G2_tec_x) | 0.023 |
| G3 | 800 | random | 11506 | 0.058 (0.75) | 0.041 (0.29) | 0.026 (0.53) | 0.082 (0.14) | 0.026 (G3_tec_x) | 0.023 |
| G4 | 800 | random | 11530 | 0.055 (0.72) | 0.045 (0.21) | 0.019 (0.50) | 0.057 (0.73) | 0.019 (G4_tecT_o) | 0.023 |
| G5 | 800 | random | 11515 | 0.081 (0.23) | 0.030 (0.36) | 0.026 (0.59) | 0.054 (0.17) | 0.026 (G5_tec_x) | 0.023 |
| G6 | 800 | random | 2157 | 0.049 (0.39) | 0.065 (0.40) | 0.025 (0.35) | 0.038 (0.67) | 0.025 (G6_tecT_o) | 0.032 |
| G7 | 800 | random | 1986 | 0.061 (0.30) | 0.081 (0.26) | 0.031 (0.30) | 0.041 (0.32) | 0.031 (G7_tecT_o) | 0.036 |
| G8 | 800 | random | 1985 | 0.070 (0.36) | 0.066 (0.17) | 0.025 (0.43) | 0.032 (0.41) | 0.025 (G8_tecT_o) | 0.035 |
| G9 | 800 | random | 2034 | 0.084 (0.24) | 0.073 (0.29) | 0.034 (0.31) | 0.038 (0.38) | 0.034 (G9_tecT_o) | 0.032 |
| G10 | 800 | random | 1980 | 0.081 (0.24) | 0.062 (0.31) | 0.035 (0.31) | 0.030 (0.41) | 0.030 (G10_ons_o) | 0.034 |
| G11 | 800 | toroidal | 559 | 1.347 (0.11) | 1.190 (0.14) | 0.446 (0.15) | 1.462 (0.04) | 0.446 (G11_tecT_o) | 0.471 |
| G12 | 800 | toroidal | 551 | 1.219 (0.07) | 1.085 (0.14) | 0.502 (0.27) | 0.867 (0.20) | 0.502 (G12_tecT_o) | 0.585 |
| G13 | 800 | toroidal | 577 | 2.961 (0.05) | 2.821 (0.06) | 1.345 (0.06) | 2.222 (0.08) | 1.345 (G13_tecT_o) | 1.364 |
| G14 | 800 | planar | 3034 | 0.162 (0.36) | 0.110 (0.45) | 0.120 (0.52) | 0.153 (0.34) | 0.095 (G14_tecT_x) | 0.083 |
| G15 | 800 | planar | 3020 | 0.174 (0.31) | 0.172 (0.22) | 0.177 (0.55) | 0.375 (0.16) | 0.079 (G15_tecT_x) | 0.072 |
| G16 | 800 | planar | 3022 | 0.209 (0.29) | 0.124 (0.28) | 0.115 (0.35) | 0.260 (0.21) | 0.070 (G16_tecT_x) | 0.083 |
| G17 | 800 | planar | 3017 | 0.174 (0.32) | 0.110 (0.34) | 0.094 (0.41) | 0.188 (0.27) | 0.094 (G17_tecT_o) | 0.083 |
| G18 | 800 | planar | 983 | 0.201 (0.23) | 0.196 (0.23) | 0.118 (0.22) | 0.110 (0.22) | 0.110 (G18_ons_o) | 0.129 |
| G19 | 800 | planar | 897 | 0.228 (0.22) | 0.145 (0.27) | 0.115 (0.35) | 0.101 (0.23) | 0.101 (G19_ons_o) | 0.114 |
| G20 | 800 | planar | 932 | 0.128 (0.37) | 0.134 (0.30) | 0.083 (0.18) | 0.067 (0.34) | 0.067 (G20_ons_o) | 0.072 |
| G21 | 800 | planar | 922 | 0.481 (0.21) | 0.367 (0.24) | 0.255 (0.21) | 0.218 (0.21) | 0.218 (G21_ons_o) | 0.245 |
| G22 | 2000 | random | 13226 | 0.080 (0.67) | 0.037 (0.30) | 0.040 (0.33) | 0.064 (0.85) | 0.037 (G22_tec_o) | 0.033 |
| G23 | 2000 | random | 13211 | 0.070 (0.21) | 0.038 (0.71) | 0.036 (0.63) | 0.029 (0.39) | 0.029 (G23_ons_o) | 0.034 |
| G24 | 2000 | random | 13204 | 0.057 (0.88) | 0.034 (0.57) | 0.044 (0.55) | 0.031 (0.35) | 0.031 (G24_ons_o) | 0.034 |
| G25 | 2000 | random | 13207 | 0.068 (0.16) | 0.038 (0.67) | 0.036 (0.62) | 0.036 (0.27) | 0.036 (G25_ons_o) | 0.036 |
| G26 | 2000 | random | 13195 | 0.057 (0.92) | 0.034 (0.65) | 0.036 (0.64) | 0.035 (0.40) | 0.034 (G26_tec_o) | 0.034 |
| G27 | 2000 | random | 3308 | 0.236 (0.56) | 0.155 (0.23) | 0.099 (0.45) | 0.137 (0.51) | 0.099 (G27_tecT_o) | 0.102 |
| G28 | 2000 | random | 3266 | 0.239 (0.70) | 0.103 (0.35) | 0.090 (0.20) | 0.139 (0.55) | 0.090 (G28_tecT_o) | 0.099 |
| G29 | 2000 | random | 3371 | 0.176 (0.37) | 0.219 (0.45) | 0.102 (0.32) | 0.183 (0.24) | 0.102 (G29_tecT_o) | 0.099 |
| G30 | 2000 | random | 3379 | 0.195 (0.19) | 0.114 (0.34) | 0.096 (0.19) | 0.136 (0.46) | 0.096 (G30_tecT_o) | 0.099 |
| G31 | 2000 | random | 3277 | 0.237 (0.17) | 0.110 (0.31) | 0.099 (0.58) | 0.140 (0.49) | 0.099 (G31_tecT_o) | 0.102 |
| G32 | 2000 | toroidal | 1396 | 20.08 (0.01) | 33.22 (0.00) | 2.586 (0.07) | 17.87 (0.01) | 2.586 (G32_tecT_o) | 2.274 |
| G33 | 2000 | toroidal | 1369 | 31.12 (0.01) | 24.77 (0.01) | 3.633 (0.05) | 15.14 (0.02) | 3.633 (G33_tecT_o) | 4.706 |
| G34 | 2000 | toroidal | 1371 | 16.45 (0.02) | 17.50 (0.01) | 2.767 (0.07) | 8.419 (0.04) | 2.767 (G34_tecT_o) | 3.194 |
| G35 | 2000 | planar | 7611 | 0.453 (0.75) | 0.317 (0.68) | 0.235 (0.27) | 0.561 (0.61) | 0.106 (G35_tecT_x) | 0.093 |
| G36 | 2000 | planar | 7604 | 0.392 (0.18) | 0.245 (0.68) | 0.179 (0.26) | 0.456 (0.60) | 0.092 (G36_tecT_x) | 0.081 |
| G37 | 2000 | planar | 7615 | 0.454 (0.70) | 0.303 (0.63) | 0.334 (0.64) | 1.024 (0.10) | 0.123 (G37_tecT_x) | 0.099 |
| G38 | 2000 | planar | 7612 | 0.453 (0.76) | 0.234 (0.24) | 0.210 (0.27) | 0.540 (0.60) | 0.097 (G38_tecT_x) | 0.088 |
| G39 | 2000 | planar | 2384 | 2.350 (0.14) | 1.389 (0.18) | 1.448 (0.23) | 1.477 (0.06) | 1.389 (G39_tec_o) | 1.374 |
| G40 | 2000 | planar | 2376 | 2.181 (0.11) | 2.299 (0.12) | 2.373 (0.07) | 1.673 (0.11) | 1.673 (G40_ons_o) | 1.965 |
| G41 | 2000 | planar | 2381 | 3.275 (0.08) | 3.346 (0.05) | 1.528 (0.08) | 2.169 (0.09) | 1.528 (G41_tecT_o) | 1.489 |
| G42 | 2000 | planar | 2457 | 2.532 (0.05) | 2.241 (0.13) | 2.175 (0.18) | 1.636 (0.11) | 1.636 (G42_ons_o) | 1.716 |
| G43 | 1000 | random | 6594 | 0.053 (0.21) | 0.028 (0.52) | 0.024 (0.50) | 0.027 (0.36) | 0.024 (G43_tecT_o) | 0.025 |
| G44 | 1000 | random | 6584 | 0.043 (0.26) | 0.030 (0.65) | 0.029 (0.61) | 0.029 (0.45) | 0.029 (G44_ons_o) | 0.028 |
| G45 | 1000 | random | 6588 | 0.046 (0.26) | 0.030 (0.63) | 0.029 (0.57) | 0.021 (0.41) | 0.021 (G45_ons_o) | 0.028 |
| G46 | 1000 | random | 6583 | 0.037 (0.28) | 0.030 (0.67) | 0.029 (0.65) | 0.029 (0.49) | 0.029 (G46_ons_o) | 0.028 |
| G47 | 1000 | random | 6591 | 0.032 (0.27) | 0.030 (0.64) | 0.033 (0.58) | 0.029 (0.45) | 0.029 (G47_ons_o) | 0.028 |
| G51 | 1000 | planar | 3810 | 0.166 (0.35) | 0.296 (0.29) | 0.126 (0.41) | 0.249 (0.26) | 0.096 (G51_tecT_x) | 0.085 |
| G52 | 1000 | planar | 3813 | 0.148 (0.36) | 0.306 (0.09) | 0.147 (0.48) | 0.729 (0.89) | 0.096 (G52_tecT_x) | 0.085 |
| G53 | 1000 | planar | 3812 | 0.190 (0.39) | 0.144 (0.42) | 0.150 (0.49) | 0.183 (0.32) | 0.096 (G53_tecT_x) | 0.085 |
| G54 | 1000 | planar | 3814 | 0.329 (0.22) | 0.257 (0.22) | 0.174 (0.28) | 0.406 (0.16) | 0.071 (G54_tecT_x) | 0.078 |

**Board power** (hwmon4; 150 s loads, 512 trials per launch, devices 98.9–99.8% busy):

| Load | Board W | Dynamic W | Energy to solution |
|---|---:|---:|---:|
| Idle | 74.4 (first); 77.1–82.2 between loads | | |
| X5 (K2000) | 126.3 | 48.0 | 6.36 mJ (0.0504 ms); 521 µJ per trial |
| G22_ons | 95.8 | 15.6 | |
| G32_tecT | 93.0 | 15.2 | |
| G1_ons_o (N = 800) | 93.2 | 15.9 | 6.29 mJ (0.0675 ms) |
| G43_ons_o (N = 1000) | 93.8 | 16.6 | 2.54 mJ (0.0271 ms) |

## PROTOCOL_HW_V6 Amendment 5: STATICA's published K2000 points on the 12 × v6.4 image (`PROTOCOL_HW_V6_A5.md`, frozen 13:10 UTC on 8 October)

**Run.** The programmed v6.4 headline image (UUID 8669ad38…, no reprogramming) and the frozen v6.4 host `build/host_sca_v64` (SHA-256 3935b5a1…). Seed 20261004, trial ids 0–2051, one trial per launch, tables resident.
- **Duration.** 13:10:55–13:11:03 UTC. The Amendment-4 chain was paused with SIGSTOP for those 8 s and resumed afterwards (`a5_a4_pause.log`).
- **Records.** `results/multibit_20261007/board/hwv6_board_v64_e12_250mhz_a5/` (`analysis.txt`, `a5_summary.json`, `rounds/`).
- **Verdict.** A5-C passes: the same-session X5 control has the same cut and flips as the original X5 run for 2052 of 2052 trial ids.

| Cohort | Schedule | p | P_round | t_round (ms) | TTS99 primary (ms) | TTS99 secondary (ms) | Model (k2000_pub): p, primary, secondary |
|---|---|---:|---:|---:|---:|---:|---|
| X5 = O1 (control) | q 8, TEC-T κ 1.75 + ramp, T 15→5, S 280 | 0.338 | 1.000 | 0.0503 | 0.0503 | 0.0468 | |
| SP_long | plain SCA, q 4, T 40→5, S 1560 | 0.824 | 1.000 | 0.3015 | 0.3015 | 0.0666 | 0.829, 0.308, 0.067 (p within band) |
| SP_short | plain SCA, q 4, T 30→5, S 560 | 0.134 | 0.848 | 0.0945 | 0.2310 | 0.2520 | 0.138, 0.249, 0.249 (p within band) |

**Ratios.** O5 = X4, from the original records: secondary 0.0181 ms.
- SP_long / O1, primary: **5.99** (O1 / SP_long = 0.167).
- SP_long / O5, secondary: **3.67** (O5 / SP_long = 0.273).
- SP_short / O1, primary: 4.59. SP_short / O5, secondary: 13.9.
- For reference, the previous "original plain SCA" baseline P3 (q 8, T0 30, S 1560) measured 0.1993 ms primary and 0.1325 ms secondary.

## Revision 4: 8192-entry schedule tables for the K = 4 + bias image (`src/v6/sca_core_mb_r4.v`)

- **Why.** The ReAIM export (`research/reaim_benchmarks_20261008/hw_export`) uses S = 8192 for every TSP schedule. Revision 3's table and trace memories hold 4096 steps.
- **What changed.** Revision 4 parameterises both memories (TLOG = 13), and nothing else. Revision 3 stays as built for the K = 2 image, which is frozen in Amendment 3.
  - The K = 4 + bias build uses revision 4 (`scripts/integrate_kernel_mb4.tcl`, host `src/host_sca_multi_mb4.cpp`).
  - `src/mb_common4.hpp` adds the SCAJINT8 / int32-bias loader.
- **xsim** (`build/sim_mb4_20261008a/b`, bit-exact against `run_trial_bias`):
  - K2000 with bias mode on and b = 0 (same suite end cycle as revision 3, 176,090);
  - random K = 4 + bias at n = 777 (same end cycle, 185,770);
  - gr17 TSP (n = 289, bias, S = 8192): traced Onsager-online plus TEC;
  - bays29 TSP (n = 841, S = 8192);
  - G1 GPP (n = 800, dense K = 4) with bias mode on.
- **Export check with the golden model.** For all six TSP instances, the SCA schedule's valid-tour count over trials 0–1023 matches the export's precision study exactly: 993, 595, 974, 844, 158 and 678.
  - The export's GPP precision values were wrong: they are the Max-Cut values of the same-named instances. See `ERRATUM_hw_export_precision_study.md` in that folder.
  - The golden model gives failure on the GPP schedules (the all-equal oscillation, or frozen states), consistent with the corrected values.
  - No GPP or TSP quality claims are made. On the board these problems serve only as bias-path and K = 4 exactness checks (Amendment 4).
- **Amendment 4** (`PROTOCOL_HW_MB_A4.md`, frozen 03:22 UTC on 8 October; hash file SHA-256 7eb65515…, 52 files including the 30 copied K = 4 export files) pre-registers the board validation of the bias path on the 6-engine K = 4 + bias image:
  - exact gate: K2000 with bias mode off and on; random K = 4 + bias at n = 777 and 2048; all 6 TSP and 3 GPP exports;
  - K2000 X5/O4 identical to v6.4;
  - for each TSP cohort, valid-tour counts over trials 0–1023 equal to the export's;
  - power.
  - Its unattended chain (`scripts/bringup_mb4b.sh`) waits for the build and for the Amendment-3 chain, and does not program after 02:00 UTC on 9 October.
- **Build at 250 MHz failed in routing** (`board_mb4b_e6_250mhz`, pipeline exit 1 at 12:06 UTC).
  - Post-placement physopt reached WNS +0.001 ns.
  - Then `route_design` failed: 4,679 signals could not be routed because of congestion (4,054 node overlaps).
  - Six K = 4 engines put about 740k LUT into each SLR pblock.
- **225 MHz rebuild.** The pre-declared rebuild started at 12:06 UTC (`board_mb4b_e6_225mhz`). The congestion is structural, so it may recur at the lower clock.

## Amendment 4: K = 4 + bias image, 6 engines, bias-path validation (`PROTOCOL_HW_MB_A4.md`, frozen 03:22 UTC on 8 October)

**Builds.** Revision 4 (`sca_core_mb_r4.v`), K = 4, BIAS = 1, HB = 17, 6 engines (2 per SLR).
- **250 MHz** (`board_mb4b_e6_250mhz`): failed in routing. Post-placement physopt WNS was +0.001 ns. Estimated congestion was level 5, and `route_design` failed with 4,679 unroutable signals and 4,054 node overlaps.
- **225 MHz** (`board_mb4b_e6_225mhz`, the pre-declared fallback): routed and passed.
  - Timing: setup +0.002 ns, hold +0.010 ns, 0 failing endpoints, signed off as PASS. The router passed through −0.049 ns before converging.
  - Utilisation: 2,213,218 LUT (86.0%), 1,370,759 FF (26.6%), 3,180 RAMB36 (+6 RAMB18; 85.1% of Block RAM tiles), 1,566 URAM (81.4%), 1,692 DSP58.
  - Image: UUID 015c977e…; PDI SHA-256 2dd9b530….
  - Report file: `results/multibit_20261007/A4_K4_BIAS_BUILD_REPORT_20261008.md`.

**Board run** (unattended chain).
- Programming 17:02–17:52 UTC. The staged v6.4 PDI was backed up first: `staged_image_backup_20261008T170203Z`, SHA-256 67e591f3….
- Amendment 17:52, power 17:52–18:04, v6.4 restored 18:49 UTC (UUID matched, BAR2 OK).
- Records: `results/multibit_20261007/board/board_mb4b_e6_225mhz/`.

| ID | Result | Verdict |
|---|---|---|
| A4-E | 156/156 trials bit-exact on all 6 engines. Bias mode off: 42/42 (K2000 W1–W4). Bias mode on: 114/114, made up of K2000 with b = 0 24/24, random K = 4 + bias at n = 777 and 2048 36/36, GPP exports 18/18, and TSP exports with bias at S = 8192 36/36 | **PASS** |
| A4-R | X5 and O4 (bias mode off) each 2,052/2,052 identical to v6.4; W1–W4 identical | **PASS** |
| A4-C | Launch cycles against v6.4: X5 0.99947, O4 0.99953 | **PASS** |
| A4-F | TSP valid-tour counts over trial ids 0–1023 equal the export exactly: 993 / 595 / 974 / 844 / 158 / 678 | **PASS** |

- **K2000 on this image** (6 engines at 225 MHz): X5 TTS99 0.0981 ms and O4 0.1171 ms. The 12 × v6.4 image at 250 MHz gives 0.0503 and 0.0528 ms.
- **Power:**

| Load | Board W | Dynamic W | Notes |
|---|---:|---:|---|
| Idle | 67.3–70.5 | | |
| X5 | 93.5 | 24.7 | 9.17 mJ to solution |
| bays29 TSP | 78.2 | 8.4 | |
| G1 GPP | 77.0 | 8.1 | |

- **Engines per device:** K = 2 runs 12 engines at 250 MHz; K = 4 with bias runs 6 engines at 225 MHz, limited by URAM read ports and routing congestion; K = 8 would fit 3 (out-of-context synthesis estimate).

## Deviations from the plan, and limits

1. **BRAM is not halved on the FPGA.** Single-copy residue banks halve the stored bits (8 Mbit per engine at K = 2, against 16). But the V80 coupling store is bound by read ports: 16 plane reads per group per cycle at K = 2. So the engine keeps v6.4's 258 RAMB36 and doubles its URAM (130). It also needs +21% LUT. The 12 engines closed at 250 MHz with only +0.007 ns of slack (88.6% LUT, 81% URAM).
2. **Simulator.** The heavy work ran on fpga-host with xsim, as instructed; v6.4b had been verified with iverilog on the laptop. xsim mis-evaluates concatenated loop-variable indices, so the multi-bit RTL uses the equivalent arithmetic form (identical in Verilog and in synthesis). A sim-only copy of v6.4b with the same rewrite passes the v6.4 vectors.
3. **N was fixed at 2000 in revisions 1 and 2**, so only G22–G42 ran there. Revision 3 adds runtime n with the shared `run_trial_bias` semantics, and Amendment 3 measured all 51 instances (N = 800, 1000, 2000).
4. **S ≤ 4096** (table memory) in revisions 1–3. The study's selections stay within S ≤ 4000. Revision 4 (the K = 4 + bias image) allows S ≤ 8192, because the exported TSP schedules use S = 8192.
5. **G-set configurations.**
   - The study's `selected_configs.json` was not available when the main protocol was frozen. The main protocol therefore pre-declares contingency configurations: a pilot, a Wilson-bound selection, and a fresh-seed validation, all with the reference model.
   - The selection rule was changed once before validation (a one-round floor), and this is recorded.
   - When the study's file arrived, it was added as Amendment 1 (hashed before data). The waiting bring-up was stopped before any board action and relaunched through a wrapper that runs the frozen bring-up unchanged.
6. **Device cut word.** The cut word was wrong for Σ s·h < 0, a bug inherited from v6.4b and latent there. It was found on the board, reproduced and fixed in simulation (`sca_core_mb_r2.v`), and re-validated as Amendment 2 on a rebuilt image (see that section). On the first image, success counts and TTS were unaffected, because success always uses the host's independent rescoring.
7. **Pre-registered failures.** G2 and G3 fail for G22_ons alone, because three oscillating trials create long-tail rounds. The other five G-set configurations pass G2 and G3. Board power had no pass condition; it is higher than v6.4's at K2000 flip rates (124 W against 108 W for X5).
8. **Not implemented.** A loader for raw dense-int8 matrix files. Edge lists (G-set), the K2000 packed file and seeded random matrices are supported in `mb_common.hpp`.
9. **K = 8 was not built.** K = 4 + bias (revision 4, 6 engines) is built for Amendment 4. The out-of-context estimates for K = 4 without bias, and for K = 8, follow.
   - K = 4: simulation-verified (|J| ≤ 7); synthesis 369k LUT, 514 RAMB36, 258 URAM per engine, so about 6 engines fit.
   - K = 8: simulation-verified on a sparse |J| ≤ 127 matrix; synthesis 786k LUT, 1,026 RAMB36, 514 URAM per engine, so about 3 engines fit.
   - Also, at K = 8 the lane temperature format (4T < 2048 in a 27-bit Q16.16 constant) limits T to below 512, so dense |J| ≤ 127 problems would need a rescaled temperature format.
10. **Engine count.** The board image uses N = 2000 and 12 engines. No engine-count or clock fallback was needed.

## Files

| What | Where |
|---|---|
| RTL | `src/v6/sca_core_mb.v` (first image), `src/v6/sca_core_mb_r2.v` (cut fix), `src/v6/sca_core_mb_r3.v` (runtime n, bias; parameters K, BIAS, HBX), `src/v6/tb_core_mb.v`, `src/v6/tb_core_mb3.v` |
| Front end | `src/v6/sca_mover_mb.cpp`, `src/v6/tb_mover_mb.cpp`; revision 3: `src/v6/sca_mover_mb3.cpp`, `src/v6/tb_mover_mb3.cpp` |
| Host and helpers | `src/host_sca_multi_mb.cpp`, `src/mb_common.hpp`; revision 3: `src/host_sca_multi_mb3.cpp`, `src/mb_common3.hpp`; references `src/sca_ref.hpp`, `src/sca_ref_bias.hpp` (shared, unchanged) |
| Vectors and exploration | `src/v6/gen_vectors_mb.cpp` (suites short/proto/long/neg, K = 1 v6.4 format), `src/v6/explore_mb.cpp`; revision 3: `src/v6/gen_vectors_mb3.cpp` (cases files, `run_trial_bias`), `src/v6/gen_gp_tsp.cpp` (brute-force-checked GP/TSP Ising instances) |
| Protocols | `PROTOCOL_HW_MB.md` (+ `_gset.txt`, `_power.txt`, `_predictions.json`, `.sha256`); `PROTOCOL_HW_MB_A1.md` (+ lists, predictions, `.sha256`); `PROTOCOL_HW_MB_A2.md` (+ `.sha256`); `PROTOCOL_HW_MB_A3.md` (+ `_gset.txt`, `_verify.txt`, `_power.txt`, `_predictions.json`, `.sha256`, hash-file SHA-256 9dec6c76…); `PROTOCOL_HW_MB_A4.md` (+ `_verify.txt`, `_cohorts.txt`, `_power.txt`, `_predictions.json`, `.sha256`, hash-file SHA-256 7eb65515…) |
| Scripts (revision 3) | `scripts/{sim_mb3.sh, run_sim_regression_mb3.sh, build_mb3.sh, chain_builds_mb3.sh, integrate_kernel_mb3.tcl, run_hw_mb3.sh, bringup_mb3.sh, measure_power_mb3.sh, gset_amend_mb3.py, analyze_hw_mb_a3.py, a3_tables.py, run_hw_mb4b.sh, bringup_mb4b.sh, analyze_hw_mb_a4.py, ooc_group_mb3.tcl}` |
| Scripts | `scripts/{sim_mb.sh, run_sim_regression_mb.sh, ooc_mb.tcl, integrate_kernel_mb{,_r2}.tcl, run_hw_mb{,_a1,_a2neg}.sh, bringup_mb2{,_a1,r2}.sh, build_mb2r2.sh, measure_power_mb.sh, analyze_hw_mb{,_a1}.py, analyze_power_mb.py, gset_select_mb.py, gset_amend_mb.py, compare_r1_r2.py, collect_mb_results.sh, post_build_mb.sh, rebuild_mb2_fallback.sh}` |
| Results (local) | `results/multibit_20261007/`: `sim/` (simulation records), `ooc/` (K = 2/4/8 synthesis reports), `data/` (21 G-set files + hashes), `gset_select/` (pilot, selection, validation), `gset_study/` (copy of the study's selection + hash), `board/board_mb2_e12_250mhz/` (main protocol, Amendment 1 (cohort files in `gset_cohorts.tgz`), power, bring-up and programming logs, build sign-off and utilisation), `SHA256SUMS_sources.txt` |
| FPGA host | `/scratch/USER/sca_v80_20261003`: `build/board_mb2_e12_250mhz` (image), `results/hwmb_board_mb2_e12_250mhz{,_a1}`, `results/power_board_mb2_e12_250mhz_20261007T212851Z`, `results/staged_image_backup_20261007T203744Z` (v6.4 PDI copy, SHA-256 67e591f3…) |
