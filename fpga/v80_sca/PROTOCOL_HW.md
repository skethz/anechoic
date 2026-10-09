# Physical V80 measurement protocol: Onsager-SCA engine (frozen 4 October 2026, before any board data)

## Purpose

1. Show that the physical engine runs the reference algorithm exactly.
2. Measure success probability and device time on the board for the configurations selected in the joint sweep J ([REPORT_ABCD.md](../../research/theory_ideas_20261003/REPORT_ABCD.md)).

J ranked configurations with the STATICA-class cycle model. This protocol re-ranks them with **measured** cycles of this engine.

## Image selection

- Two images are built from the same co-simulated RTL: 250 MHz and 300 MHz kernel clocks.
- Use the highest-clock image that passes routed timing (setup, hold and pulse-width slack ≥ 0) and image-manifest checks.
- Results at the other clock are not collected. Only the timing scale would differ, because the arithmetic is clock-independent.

## Programming

- Programming uses only the scoped commands from `fpga/v80_snowball/scripts/admin/README.md`: stage the image to the one allowed path, then `ami_tool cfgmem_program` to partition 1. If BAR2 is missing after re-enumeration, run the argument-free rescan helper.
- No reboot, no general PCI removal.
- Recovery image: the corrected Snowball image `build/board_300mhz_counter_fixed` (UUID `27fadf27c7874b94b07c60df1337bc7e`), whose manifest is unchanged.

## Phase 1: exact bring-up (gate)

`host_sca --verify` with per-step n_lin trace. Every trial must match `sca_ref.hpp` exactly: all 2,000 spins, cut, flip count, and every n_lin value.

| ID | Mode | q | λ / J_v | Ramp | T_init→T_fin | S | Trials | Trials per launch |
|---|---|---|---|---|---|---|---|---|
| V1 | Plain | 4 | 0 | No | 30→5 | 40 | 8 | 4 |
| V2 | Onsager | 4 | 0.7 | Yes | 20→5 | 560 | 8 | 8 |
| V3 | TEC | 8 | J_v=−8 | No | 20→5 | 40 | 4 | 1 |
| V4 | Onsager (E1 selection) | 8 | 1.05 | Yes | 12→5 | 960 | 4 | 2 |

Seed 20261004, trial offsets starting at 900000. Each host invocation loads the couplings on its first launch only. V2 runs one launch of 8 trials; V1 and V4 run two launches, the second with `load_J=0`, so reuse of the BRAM-resident couplings is checked; V3 runs four single-trial launches. **Any mismatch stops the protocol.** Phase 2 does not run.

## Phase 2: physical cohorts

- **Configurations.** The 11 distinct configurations selected in J (fixed q, T_fin=5, geometric T, ramp = linear λ to 0 over the last 30%).
- **Trials.** 1,024 per configuration, seed 20261004, trial offsets 0–1023 for every configuration. The same seeds give the same initial spins across configurations, so comparisons are paired.
- **Launches.** 64 trials per launch (16 launches). The couplings are loaded on the first launch only.
- **Single-run latency.** An extra 64 trials (offsets 100000–100063), one trial per launch, for the J E=1 selections (plain, TEC, Onsager).

| Key | Mode | q | λ or J_v | Ramp | T_init | S |
|---|---|---|---|---|---|---|
| P1 | Plain | 4 | 0 | No | 30 | 1560 |
| P2 | Plain | 6 | 0 | No | 30 | 1560 |
| P3 | Plain | 8 | 0 | No | 30 | 1560 |
| P4 | Plain | 8 | 0 | No | 30 | 960 |
| T1 | TEC | 8 | −4 | No | 30 | 1560 |
| T2 | TEC | 8 | −4 | No | 30 | 960 |
| T3 | TEC | 8 | +4 | No | 30 | 1560 |
| O1 | Onsager | 8 | 1.05 | Yes | 12 | 960 |
| O2 | Onsager | 6 | 0.9 | Yes | 15 | 560 |
| O3 | Onsager | 6 | 0.9 | Yes | 12 | 560 |
| O4 | Onsager | 6 | 0.9 | No | 12 | 360 |

**Integrity checks (every trial).**
- The host rescores every returned state with an independent scalar edge loop and requires `device_cut == host_cut`.
- Padding bits must be zero.
- Any disagreement is reported. Such a trial counts as a failure, not a success.

## Analysis (fixed in advance)

- **Success.** A trial succeeds if its final-state cut is ≥ 33,000 (host score). p is reported with a Clopper–Pearson 95% interval.
- **Device time per run.**
  - `t_run = launch_cycles / trials_in_launch / f_clk`, averaged over the 16 launches. This includes table loading once per launch, amortized over 64 trials.
  - The single-run cohort reports the full per-launch device latency, including table loading, as mean, median and maximum.
- **TTS.**
  - `TTS99 = t_run · ln(0.01)/ln(1−p)` (Eq. 9, unrounded).
  - Whole-run form: `t_run · ⌈R99⌉`.
- **E>1 projection.** Replicated-engine values are labelled projections, not measurements: `t_run · ⌈⌈R99⌉/E⌉`, assuming E identical engines at the measured per-run time.
- **Consistency with the J model run.** For each configuration, the 95% interval of p_hw − p_J is reported, where p_J is the Python held-out p (1,024 runs, different RNG). The run counts as "reproduced" if |p_hw − p_J| ≤ 0.05 or the interval contains 0.
- **Cycle-model check.** Measured cycles per trial are compared with the J model evaluated at the trial's own flip count, `0.434·flips − 0.98·S + 1989 + 4S·[correction on]`. The ratio is reported. A per-step fit `cycles ≈ a·flips + b·S + c` is reported as descriptive only.
- **Primary comparison.** Measured E=1 TTS99 of the J E=1 Onsager selection (O1) against the best measured plain configuration (P1–P4) and the best TEC configuration (T1–T3).
  - This is a post-hoc minimum over the measured configurations, so the comparison favours the baselines.
  - The speed-up is reported with a bootstrap 95% interval (10,000 resamples of trials within each configuration).

## Not claimed

- No external review.
- No claims for other graphs or N.
- Engine replication (E>1) is not built in this protocol.
- Results use the 16-bit xoshiro128** thresholds of this engine, not STATICA's RNG.
