# Addendum 5: Table 8b under the authors' decisions after the fairness audit (8 October 2026, ~15:20 CEST, before any A5 run)

## Source

- **What these decisions are.** The authors' note of 8 October (after the validation in `research/fairness_audit_20261008/k2000_pub`, PROTOCOL_PUB_V1/V2/A3) gave decisions "so both studies agree". This addendum applies them as given.
- **What had already been seen when it was written:**
  - all Addendum 4 validation results (`validation_a4/`, `validation_a4_summary.json`);
  - the A4 and A4.1 re-runs;
  - the exploratory `diag_tec.py`.
- **What is not seen yet:** the exploratory `diag_statica.py` is still running.
- **Who made the choices.** Nothing below is chosen from our results. The k sets, the Δt rule, the aSB setting and the TEC/STATICA treatment are the authors' decisions.

## A5.1 ReAIM ASA: ReAIM Table I's per-problem k

- **Table I, "Best flips k per iteration":** MCP ≥ 16, GPP 6, TSP 1.
  - **MCP:** k ∈ {128, 256, 512, 1024}. These are all ≥ 16, and with this set the audit reproduces ReAIM's K2000 Table VII: P_a 0.54/0.74 against 0.47/0.80. Each k is capped at N, so G-set (N = 800) uses {128, 256, 512, 800}.
  - **GPP:** k ∈ {6}.
  - **TSP:** k ∈ {1}.
  - **Replaces:** A3's union {1, 2, 6, 16}, which the audit found gives P_a = 0.00 on K2000.
- **Unchanged:**
  - Table II T_init/T_final;
  - F per Sec. IV-C (max for MCP and GPP, min for TSP);
  - ITER_trial/ITER_run 32/96 (not given by the paper; the audit uses the same);
  - FIFO initialised with N;
  - noise-free;
  - **x_best output:** Algorithm 3, the lowest-energy run-phase end state (`solvers_a3b.reaim_best`).
- **Code:** `run_a5.reaim_cfg`.
- **Finals:** every instance and S of the grid, 256 runs, `SeedSequence([20261008, 80, p, i, m, s])`, written to `results_a5/final/reaim_t1/`.

## A5.2 bSB, dSB: one Δt per instance by time-to-target (Goto 2021)

- **What Goto 2021 says:**
  - "the time step Δt for bSBM and dSBM is set to the best value for each problem among five values (0.25, 0.5, 0.75, 1, and 1.25)";
  - "The number of time steps, N_step, is also optimized for TTT or TTS separately";
  - TTT = T_com·log(1 − 0.99)/log(1 − P_S), with target = "99% of the optimal or best known value".

  This is our MCS99 with our targets (MCP ⌈0.99·BKV⌉, GPP ⌊1.01·R⌋, TSP ⌊1.01·L*⌋), with T_com ∝ S.
- **Pilots:** for each instance, method, Δt ∈ {0.25, 0.5, 0.75, 1, 1.25} and S of the grid: 256 runs, `SeedSequence([20261008, 82, p, i, m, s, k])`.
  - TTT_pilot(Δt) = min over S of MCS99.
- **Selection:**
  - the Δt with the smallest TTT_pilot;
  - ties (including all infinite) go to the higher pilot mean quality at the S attaining the minimum (the largest S when all are infinite), then to the smaller Δt.

  This is the audit's rule ("one Δt ... chosen by the minimum over S of the steps to solution on pilot runs"), applied per instance as Goto 2021 does per problem.
- **Finals:** all five Δt at every S, 256 runs, `SeedSequence([20261008, 83, p, i, m, s, k])`. Table 8b uses the finals at the selected Δt; all five are kept in the supplement.
- **TSP fields:** the ancillary spin of Addendum 4.1.
- **Unchanged:** a0 = 1, c0 = 0.5/(⟨J⟩√N), and initial conditions as in A3.
- **Replaces:** A3's per-budget choice by pilot mean quality.

## A5.3 aSB: the 2019 setting only

- **Setting:** Δt = 0.9, M = 2, ξ0 = 0.7/(SD(J)√N), p linear from 0 to 1. This reproduces the paper's own number at N_step = 186:
  - **ours (V6):** 32,761.8 over the 255/256 finite runs;
  - **the audit:** 32,756;
  - **the paper:** 32,768.
- **Sources for Table 8b:** A3 finals for MCP and GPP; the A4.1 `aSB_dt09` finals (ancillary spin) for TSP.
- **Divergence.** Runs with non-finite dynamics count as failures (quality 0). The fraction per problem is reported as a property of this setting in floating point, not as a failure of the method: the 2019 FPGA used saturating fixed point.
- **Supplement only:** Δt 0.5 and the stability-margin rule (A4 (b) and (c)).

## A5.4 TEC: not reproduced

- **Pre-registered V3 failed:**

  | | Ours | Paper (Fig. 3b) |
  |---|---|---|
  | J_v = 0 | 739 cycles | 1,392 |
  | J_v = 30J | 577 cycles | 740 |
  | Ratio | 1.28 | 1.88 |

  The Fig. 3(a) ordering is reproduced: J_v = 6 is faster and J_v = −6 slower than J_v = 0.
- **Exploration also failed.** The exploratory `diag_tec.py` (index-order, random-site and random-permutation updates × geometric or linear schedule) brackets the paper's J_v = 0 value but reproduces neither value nor the ratio (geometric: 745–752 and 575–606; linear: 2,447–2,460 and 2,134–2,155).
- **The audit agrees.** Its independent kernel gives 736/577.
- **Consequence.** TEC appears in Table 8b as "not reproduced", with no numbers. The A3 synchronous reading and the A4 sequential σ-transfer files are kept on disk but not reported.

## A5.5 Settings not specified by the paper

- **In Table 8b:**
  - STATICA on G-set Max-Cut, GPP and TSP: "settings not specified by the paper" (its paper gives K2000 settings only).
  - APC-SCA on GPP: also "settings not specified" (its paper covers Max-Cut and TSP).
- **Supplementary record, labelled:**
  - STATICA with the σ-transfer (A3, OUR RULE);
  - APC-SCA on GPP with its TSP-column settings (A3.1, OUR CHOICE).

## A5.6 Final Table 8b composition

| Method | Max-Cut G1–G20 | GPP (9) | TSP (6) |
|---|---|---|---|
| SA (Neal defaults) | A3 | A3 | A3 |
| STATICA | not specified | not specified | not specified |
| TEC | not reproduced | not reproduced | not reproduced |
| APC-SCA | A3.1 (paper's Max-Cut settings) | not specified | A3.1 (paper's TSP settings; its TSP results cannot be reproduced: unpublished Ising model) |
| ReAIM ASA | A5.1 | A5.1 | A5.1 |
| aSB | A3 (Δt 0.9, M 2) | A3 | A4.1 aSB_dt09 |
| bSB, dSB | A5.2 | A5.2 | A5.2 |
| Onsager-κT, Onsager-online (ours) | A3 | A3 | A3 |

- **Statistics:** as A3. Mean quality at the ReAIM-matched S (MCP 4000, GPP 4096, TSP 8192), feasibility, instances solved (finite MCS99 at some S), and the geometric-mean MCS99 over the Max-Cut instances that every reported row solves.

## A5.7 Validation V8 (`validate_a5.py`; `SeedSequence([20261008, 81, crc32(job), chunk])`)

- **V8a: K2000, k {128, 256, 512, 1024}, 4,096 and 6,400 iterations, 512 runs each.**
  - Measure: P_a = P(cut ≥ 33,000) of x_best, against ReAIM Table VII's 0.47 and 0.80.
  - Pass: each paper value is inside our 99% Wilson interval or within 0.05 (the audit's criterion).
  - If V8a fails, ReAIM is listed as "not validated" and its numbers are not used.
- **V8b: Tables IV–VI with the per-problem k sets,** 240 runs (12 blocks of 20). Deviations are reported, with no pass/fail.

## A5.8 Notes on the A4 validation, recorded before the A5 runs

- **V6 (aSB): two frozen criteria disagree.**
  - The frozen protocol text says "within 3 standard errors of a 100-run mean". The finite-run mean passes: z = −0.46.
  - The frozen summary code also required every run to be finite, and 1/256 diverged. That condition is not in the protocol text.
  - Both are reported. The verdict follows the protocol text, with the divergence stated.
- **V2 (STATICA long point): fails its pre-declared P_a criterion.**
  - The mean is within the criterion: 33,095.2, z = +2.42 against a 100-run SE.
  - P_a = 0.867 is above the paper's 95% interval [0.675, 0.848].
  - This also differs from `research/statica_reproduction_20261003` (0.827, n = 4,096) and the audit (0.823), which run the same algorithm on the same graph. `diag_statica.py` (exploratory, running) tests whether the gap is statistical or systematic.
  - STATICA has no Table 8b numbers either way (A5.5).

## Files

- `run_a5.py`: finals and pilots; `select` writes `results_a5/sb_ttt_selection.json`.
- `validate_a5.py`: V8.
- `analyze_a5.py`: Table 8b. Written after the runs; its logic is fixed by A5.6.
- **Directories:** `results_a5/`, `validation_a5/`.
- **Seed families:** 80–83 (new).
