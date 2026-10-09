# GH200 K-bit SCA engine with bias (external field) and runtime n: verification, step time, G-set Max-Cut

Frozen 8 October 2026 (UTC), before any verification, timing or measurement run of the binaries below. Changes after freezing go into dated amendments at the end, each hashed before its data. Development runs under `tmp/` are not results.

Builds on `gpu/multibit_20261007` (protocol `e8acc882…` and amendments). Estimators, timing boundary and intervals are as there (`src/analyze.py`, unchanged).

## 1. Platform and rules

- Host `gpu-host`, GPU `GPU-xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx` (GH200) only, selected by UUID; the binary refuses any other device.
- Project directory `/scratch/USER/anechoic_gpu_multibit_bias_20261008`; caches and temporary files stay there. The multibit study's directory and binary are only read (identity checks, timing reference).
- GPU1 is shared (another process holds about 37 GB at 0 % utilisation at freeze time). GPU state and host load are logged before and after every run (`logs/gpu_status.log`).
- Seed 20261004 throughout.

## 2. Implementation under test

**Spec.** `fpga/v80_sca/MULTIBIT_SPEC.md`, including "Bias (external field)".

**Golden reference.** `sca::run_trial_bias(n, dense, bias, tables, seed, trial, trace)` in `src/sca_ref_bias.hpp`, an unchanged copy of `fpga/v80_sca/src/sca_ref_bias.hpp`. Its own test passes on gpu-host: `test_ref_bias`, 16 of 16 checks.

**Kernel changes against `sca_gpu_mb.cu`.**
1. The fields start at the bias, h(0) = b + J s(0). The per-step update is unchanged.
2. A WIDE instantiation clamps s·h to ±HC inside the decision, with 65536·HC > xm_max and xm_max = max_t(|q| + corr_max + 4T).
   - For |s·h| ≥ HC, both the flip test and the n_lin band test are decided by the sign of s·h, so clamping is exact.
   - It is used only when (max_y Σ_x |J_xy| + max|b|)·65536 + xm_max ≥ 2^31, i.e. when int32 decisions could overflow (`--wide auto`). `--wide 1` forces it for tests.
3. Host side:
   - bias files (SCABIAS1, |b| ≤ 32767, zero padding);
   - the range plan;
   - energy E = −(Σ s h + Σ b s)/2 from the device fields, checked against an independent host energy −Σ_{i<j} J s s − Σ b s;
   - `--check-fields 1` compares every final field with b_i + Σ_j J_ij s_j recomputed from the final spins.

**Runtime n.** The kernel already took n at runtime (padding slots are never decided, every lane draws 8 numbers per step); it is verified here for n ≠ 2000.

| File | sha256 |
|---|---|
| `src/sca_gpu_bias.cu` | `0f906503c496e4cb056421ec25b7b54972d5df0646c71940c1d4ad48a105dd52` |
| `src/jmat.hpp` | `e8682a889fc9184cf31f6bdbd2514878d7350489336e2778d63e6d3f0e31832a` |
| `src/check_ref_bias.cpp` | `b73b9af8f7739622e0e8d2f0beb7fa011a4991ab21a3d5c8253c057b2b079d2d` |
| `src/gen_ising.cpp` | `e3841a4f1c54a896a0017e2fc6b85c500ab86440eeca8692eb0584a28c2a1f32` |
| `src/gset_plan.py` | `dbcc24885c7ee9939ce98e4d5062483cdb77e439f74a70673e1f9f3f622d9486` |
| `src/run_bias.sh` | `d5191312640ca5bba6f7e4f2d84e5d50eda80273740a587c7d7d74ff8b7e5146` |
| `src/compare_runs.py` | `e509c88952de0c99f4f8a3365e7d95c134b15b01f7f5cfc7e465998ce17586eb` |
| `src/sca_ref.hpp`, `src/sca_ref_bias.hpp`, `src/test_ref_bias.cpp` (unchanged copies) | `24d04a41…`, `d063b77c…`, `0f77457e…` |
| `src/tables.hpp`, `src/analyze.py` (unchanged) | `f1805772…`, `d6f2ba71…` |
| `src/build.sh`, `src/environment.sh` | `90341f2d…`, `c8718884…` |

**Build.** `bash src/build.sh` (nvcc 13.3 `-O3 -arch=sm_90 -DSCA_LANES=256`, g++ 13.3). Binary hashes:
- `build/sca_gpu_bias`: `b6cd781e426f08456dddb47ff574e9aac518ba331f71873a1e5799a01c7a5d31`;
- `build/check_ref_bias`: `746f0a53…`;
- `build/gen_ising`: `51fe3024…`;
- `build/test_ref_bias`: `f905e089…`.

Spills: 34 of 42 instantiations are spill-free. The exceptions are r1g2 (all modes) and r2g4 (n_lin modes), as in the multibit study.

**Instances** (`data/inst`, `gen_ising`; hashes in `data/inst/SHA256SUMS`, sha256 `376c8ad1…`):

| Instance | Construction | n | Bias |
|---|---|---:|---|
| rb2_n1000 | random ternary | 1,000 | ±40 |
| rb2_n257 | random ternary | 257 | ±20 |
| rb2_n17 | random ternary | 17 | ±5 |
| rb4_n800 | random K = 4 | 800 | ±300 |
| rb4_n1999 | random K = 4 | 1,999 | ±300 |
| rb4_n2048 | random K = 4 | 2,048 | ±300 |
| rb4_n1 | single spin | 1 | ±3 |
| wide2_n1500 | ternary | 1,500 | half of the spins up to ±32,767, half ±30 |
| wide4_n1200 | K = 4 | 1,200 | half up to ±32,767, half ±200 |
| gp22 | graph partitioning: 3-regular graph, partition sizes k = 8 and 14, A = 3, B = 1; brute-force checked: the Ising ground states are exactly the 19 optimal partitions | 22 | yes |
| gp1000 | graph partitioning: 3-regular graph, partition sizes k = 400 and 600 | 1,000 | yes |
| tsp6 | TSP, 6 cities, distances in {1, 2}, A = 3, B = 1; brute-force checked: the ground states are exactly the 8 optimal tours | 25 spins | yes |
| tsp40 | TSP, 40 cities | 1,521 spins | yes |

The two WIDE instances need the clamped decision. The GP and TSP instances need K = 4.

G-set files (51 instances, G1–G47 and G51–G54) were downloaded from Stanford on gpu-host. Their hashes are in `data/gset/SHA256SUMS` (sha256 `aa40ffe2…`); G22 and G32 equal the earlier downloads.

## 3. Phase V: verification (must pass before any other phase)

`bash src/run_bias.sh V`. Every trial is re-run with `run_trial_bias(trace = true)` (`check_ref_bias`). The following must be equal: final spins, Σ s h, Σ b s, energy (host and device), total flips, and every per-step n_lin. Every run also uses `--check-fields 1`.

- **V0, b = 0 regression on K2000.** Each of the 15 multibit phase-V1 K2000 runs (plain, Onsager, TEC-T, TEC; 2 to 240 chains) is repeated with each of the 7 variants, using the same schedule and trial ids. Each run must also be identical to the multibit kernel's raw files (`compare_runs.py`: per-trial cut and flips, final-spin file, n_lin trace). Four of the runs, one per mode, are repeated with an explicit all-zero bias file and with `--wide 1`. In total: 161 runs.
- **V1, biased instances.**
  - Every instance in all four modes, 8 chains each, with 300 steps (TEC: 150); schedules are in `run_bias.sh`.
  - K = 2 instances run with all 7 variants; K = 4 instances with k4.
  - Full waves, so that every variant runs a full wave on some instance.
  - `--wide 1` forced on rb2_n1000 (all 7 variants), rb4_n800, tsp40 and gp1000.
- **Pass criterion.** Any mismatch, any identity difference, or any device/host or field mismatch stops the protocol.

## 4. Phase T: step time (timing only; outcomes are identical by construction)

`bash src/run_bias.sh T`. 16 batches per run, three interleaved repetitions.

- **T1.** The bias binary (no bias) against the multibit binary on K2000, with the same trial ids. Configurations:
  - X5 (TEC-T, S = 280) and X4 (TEC-T, S = 760);
  - O4 (Onsager, S = 360) and O1 (Onsager, S = 960);
  - at B = 1, 60, 120 and 240 with the multibit phase-K variants.

  Reported: mean t_batch per binary and their ratio; the per-step slope per binary from the two S values per class, and the slope ratio. The expectation is a ratio of 1.00 within run-to-run spread.
- **T2.** With bias against an all-zero bias file (same binary, J and trial ids): rb2_n1000 at B = 120 (r1g1, Onsager and plain) and tsp40 at B = 60 (k4, Onsager and plain).
- **T3.** WIDE against the plain kernel on K2000: X5 and O1 at B = 1, 120 and 240.
- **Static check (not a timing).** Instruction counts of the step loop for the bias and the multibit kernels, from `cuobjdump -sass`.

## 5. Phases GV and G: G-set Max-Cut, K = 2, b = 0, all 51 instances (N = 800, 1000, 2000)

**Configurations.**
- Source: `research/gset_20261007/selected_configs.json` (software study), copied unchanged to `data/gset_selected_configs.json` (sha256 `291401004efe27e45aed38c7f402b9c85d24f1df769e6cceadfb45abf6f8174d`).
- Version `extended_amendment1` (the study's later, wider-grid selection).
- Rules: SCA (plain), TEC, Onsager-kT and Onsager-online. APC-SCA is not an engine mode and is omitted.
- Per instance and rule: the study's `tts_E1` (single engine) and `tts_E12` (12 engines) configurations. This gives 408 cells.
- `src/gset_plan.py` maps the study's rule definitions one-to-one onto the V80 table flags: absolute q, T0 and tfin; κ and λ with the end ramp; jv for TEC.

**Target.** The study's per-instance target ⌈0.99 × BKV⌉. Success means a final-state cut ≥ target, rescored on the host from the edge list.

**Variants.** Variants follow the multibit phase-K selection per (B, class). If the selected variant cannot hold B chains in one wave at the cell's S (shared-memory variants at S = 4,000), the register variant takes its place: r2g1, r2g2 or r2g4 for B = 60, 120 or 240.

**GV, validation (must pass before G).**
- Each E12 configuration is traced with 4 chains on every variant it uses at B = 60, 120 and 240.
- Each E1 configuration is traced on its B = 1 variant.
- All traces are checked against `run_trial_bias(n = N, b = 0)`, with trial ids from 39,000,000. Any mismatch stops phase G.

**G, counted runs** (cell index c = 4 × instance order + rule index; id base 40,000,000 + 100,000·c):

| Configuration | B | Runs | First trial id (offset from the base) |
|---|---:|---|---:|
| E12 | 60 | 128 batches | +0 |
| E12 | 120 | 128 batches | +20,000 |
| E12 | 240 | 128 batches | +50,000 |
| E1 | 1 | 256 single-chain launches | +90,000 |

**Reported.**
- Per cell: p, P_batch, t_batch, µs per step, primary and secondary TTS99 with Wilson 95 % intervals, mean cut and best cut.
- Beside them, the study's modelled V80 times (`tts_ms`, `p` from the same file), labelled as the study's numbers.
- Aggregates by graph class (random, toroidal, planar; by N, and by weights +1 or ±1).
- The best GPU cell per instance, over rules and B, is reported as such; it is a selection among 12 held-out cells, not a pre-registered single estimate.

## 6. Phase W: power

nvidia-smi on the UUID every 200 ms. Sequence:
1. 30 s idle.
2. 60 s of back-to-back batches of the G22 Onsager-online E12 configuration at B = 120.
3. 30 s idle.
4. 60 s of tsp40 (K = 4, biased) at B = 60 with its V1 Onsager schedule.
5. 15 s idle.

## 7. Graph partitioning and TSP from `research/reaim_benchmarks_20261008/hw_export`

The export did not exist at freeze time. When it exists, an amendment hashed before its data will declare the instances, schedules, targets, runs, validation and power, using K = 4.

K = 8 is not built:
- Decisions would be exact through the clamp (the WIDE kernel accepts any |h| < 2^31).
- But eight 2048 × 2048 bit planes (4 MB) do not fit the registers and shared memory of an 8-SM cluster.

## 8. Known limits stated in advance

- One shared GPU; clocks are not locked.
- G-set schedules are the software study's float-model selections. They are used unchanged (a transfer test), not re-tuned for the GPU's batch sizes.
- The validation instances test arithmetic, not solution quality.

## Amendment 1 (8 October 2026, during phase G, before any data of the phases below)

- **Reason.** The ReAIM study's planned TSP budgets reach S = 8192 (its `PROTOCOL.md` and `run_bench.py`; read only). Its export rule takes the largest S when no S reaches the target. With K = 4, two coupling planes (128 KB) plus 20 bytes of step tables per step exceed the 227 KB of shared memory beyond S ≈ 4,600, and the frozen binary accepts at most 4,096 steps.
- **New binary** `build/sca_gpu_bias_tg`, from the new file `src/sca_gpu_bias_tg.cu`. `src/sca_gpu_bias.cu` and `build/sca_gpu_bias` are unchanged and remain the binaries of phases V, T, GV, G and W.
  - The source adds a template flag TG: the per-step tables are read from global memory through the read-only path, and the decision's table values are loaded before the message wait.
  - TG = false instantiations are the unchanged source; their ptxas spill pattern is identical to the frozen build.
  - `--tables auto` selects TG iff S > 4096 or the tables do not fit beside the planes. `--tables smem|global` force either.
  - The step limit is 65,536.
  - Sources: `src/sca_gpu_bias_tg.cu` (`00e2f1e853ac6820e21778c24f51afc5d885427dbb5623f20eccba690f94397a`), `src/run_bias_b.sh` (`1f646b85dbb6a6f87d8b93f222b4dc61f82e30f3bef36769763429e7b3586180`).
  - Build: nvcc as `build.sh`; `build/sca_gpu_bias_tg` sha256 `816bc560289fd197ca0d1d1d872707a182edeba193d3010a3995d605a4a9eeae`.
- **Helpers for phase B.**
  - `src/raw2bias.py` (`134e79da…`): raw int32 bias to SCABIAS1, values unchanged.
  - `src/tables_equal.cpp` (`42e160fa…`) with an unchanged copy of `fpga/v80_sca/src/mb_common.hpp` (`2934a3a3…`). It checks, per schedule, that the engine's table code equals `mb::tables` (the study's `hw_check` path).
- **Phase BV, verification of the new binary** (`bash src/run_bias_b.sh BV`; must pass before BT and B). Every run uses trace and `--check-fields 1` and is checked against `run_trial_bias`. Any mismatch stops.
  - (a)/(b): one multibit K2000 run per mode × 7 variants, each with `--tables smem` and with `--tables global`, plus identity with the multibit raw files (56 runs).
  - (c): biased instances with global tables: rb2_n1000 (7 variants), rb4_n800 and tsp40 (all four modes), gp1000, and wide4_n1200 (WIDE together with TG).
  - (d): schedules beyond the shared-memory limit: tsp40 at S = 8192 (8 and 60 chains), rb2_n1000 at S = 6000 and 5000 (including 240 chains on r1g2), and K2000 at S = 8192.
- **Phase BT, timing of global against shared-memory tables.** Three interleaved repetitions of 16 batches:
  - K2000 X5 and O1 at B = 1, 60, 120 and 240 with the phase-K variants;
  - tsp40 (k4, Onsager and plain, S = 300) at B = 1 and 60.

  Reported as ratios.
- **Phase B, graph partitioning and TSP.** A second amendment will fix the instance files, their hashes, the schedules, the targets, the evaluation and the runs, once `research/reaim_benchmarks_20261008/hw_export/` exists. It will be hashed before any phase-B data. K = 4, with the k4 variant; B = 60 (one full wave) and B = 1.

## Amendment 2 (8 October 2026, 02:1x UTC, after a failed start of phase W, before any power data)

- **What failed.** Phase W read the plan line of G22_Onsageronline_E12 without the tab separator (the GV and G loops use it), so the schedule flags were cut after `--lambda`.
  - The burn binary stopped at once ("FAIL: stod") after the 30 s idle window, so no load was applied.
  - Under `set -e` the script stopped before stopping its nvidia-smi sampler; the sampler was killed by hand.
  - The attempt (idle samples only) is kept unchanged as `results/W_attempt1` and `logs/phase_W_attempt1.log`, and is not used.
  - The chain therefore stopped before BV and BT, as designed.
- **Fix.** `src/run_bias.sh` reads that line with `IFS=$'\t'`; nothing else changes. New sha256 `c840d63a31c9e7e87f4971c3d801e0db510a575ee0d2961cdc0590a89be62c2d`.
- **Then.** Phase W is rerun as specified in section 6, followed by BV and BT of Amendment 1.
- **Phase B.** The amendment for graph partitioning and TSP (Amendment 1's "second amendment") becomes Amendment 3.

## Amendment 3 (8 October 2026, ~02:45 UTC): graph partitioning and TSP (phase B), before any phase-B data

- **Source.** `research/reaim_benchmarks_20261008/hw_export/`, read only, copied unchanged to `data/hw_export/` on gpu-host. 91 files; `selected_configs.json` sha256 `776d05cbb8399b26f90cf6fcb2e97a4135b5e1a6ec1b12ecf553476bbe791eca`; the list of file hashes is `data/hw_export_SHA256SUMS` (sha256 `7c5c8064…`). The study's own `export_check.json` did not exist at copy time. This study's checks (below) do not depend on it.
- **Instances and precision.** All 15 exported instances: GPP bisection on G1–G5 and G14–G17 (N = 800), and TSP gr17, gr21, gr24, fri26, bayg29 and bays29 (N = n² = 289–841).
  - **K = 4 files only** (`<inst>_K4.jint8`, `<inst>_K4.bias.bin`).
  - K = 2 is not part of the request.
  - K = 8 is not run: the engine is not built for eight planes. Its int32 decisions with clamping would be exact, but 4 MB of planes do not fit an 8-SM cluster.
- **Cells.** The four rules exported per instance (SCA, TEC, Onsager-kT, Onsager-online) give 60 cells.
- **Schedules.** Each rule's entry (t0, t1, S, q, lam, jv, kappa, ramp; seed) maps to `--t0 --t1 --steps --q --lambda --tec-jv --tecT [--ramp]`, with tables identical to `mb_common.hpp::tables` (scale 1), checked per cell (`tables_equal`).
- **Engine.** `build/sca_gpu_bias_tg` (Amendment 1), variant k4 (60 chains per wave), `--tables auto`: global tables for S > 4096 (all TSP cells, S = 8192), shared memory otherwise. The range plan is recorded per cell, and the import step requires 60 chains per wave.
- **Objective and targets.** As in the export's objective mapping, evaluated on the true, unquantized problem from the final spins by `src/eval_problem.py` (independent code).
  - GPP: G-set graph (`data/gset`); feasible iff Σ s = 0; cut ≤ target; quality R/cut.
  - TSP: own TSPLIB parser, which reproduces the published optimal tour lengths of gr24, fri26, bayg29 and bays29; feasible iff permutation matrix; L ≤ target; quality L*/L.
  - Targets and references are the export's ("feasible and cut/L ≤ T"; reference_R, optimum_Lstar).
- **Step BX, import.** Hashes; raw bias to SCABIAS1; a full mapping check per instance (the exported matrix and bias must equal round(α·J) and round(α·b) rebuilt from the G-set edges or TSPLIB distances with P = the export's full-precision max|J|); table identity; range plan; residency. Any failure stops.
- **Step BV2, validation.** Each cell is traced with 4 chains (trial ids from 30,000,000) and checked against `run_trial_bias` (spins, Σ s h, Σ b s, energy, flips, n_lin). Any mismatch stops.
- **Step BR, replication.** Each cell is run on the entry's seed with trial ids 0–1023 (32 batches of 32 chains), exactly the trials of the export's `precision_study` (the study's bit-exact model). The GPU's feasibility, P(target) and mean quality are compared with the exported values; they are expected to be identical. This is a check, not a TTS measurement.
- **Step B, counted runs** (entry's seed):
  - B = 60, 128 batches, trial ids from 10,000,000;
  - B = 1, 256 single-chain launches, trial ids from 20,000,000.

  Reported per cell: feasibility, P(target) and mean quality (all runs and feasible runs), best value against the reference, optimum hits, TTS99 primary (P_batch at B = 60) and secondary (p) with Wilson intervals, and single-chain TTS99 (B = 1). With no target hit, TTS is ∞ and not reported as a number.
- **Step BW, power.** The TSP cell and the GPP cell with the largest S·N, each 60 s back-to-back at B = 60, with idle windows.
- **New files.** `src/eval_problem.py` (`4b8b20f7…`), `src/hw_import.py` (`31372e39…`), `src/run_bias_b2.sh` (`f3b19756…`); collector `src/make_tables_bias.py` (`847e2f8f…`); TSPLIB copies `data/tsplib` (hash list sha256 `e200c395…`).

## Amendment 4 (8 October 2026, ~02:50 UTC, after a failed step BX, before any phase-B data)

- **What failed.** Step BX stopped at G17 Onsager-online (λ = 97.2, t1 = 0.90). The host range plan of `sca_gpu_bias_tg` refused the schedule because the clamped (WIDE) fallback would not be exact (65536·HC + xm_max ≥ 2^31). But the plain int32 kernel is exact for it: H·65536 + xm_max = 1.98·10^9 < 2^31. The refusal was a host-side logic error: it should refuse only if neither kernel is exact.
  - A replica of the plan over all 60 cells finds this the only affected cell; no cell needs WIDE.
  - No GPU run had been made. The failed log is kept as `logs/phase_BX_attempt1.log`. Derived files (converted biases, specs) are regenerated.
- **Fix.** `src/sca_gpu_bias_tg2.cu` (`a4a78887…`) is `sca_gpu_bias_tg.cu` with only the host plan changed:
  - refuse iff neither the plain nor the clamped kernel is exact;
  - `--wide 1` refuses if clamping is not exact.

  `build/sca_gpu_bias_tg2` sha256 `bab96a1005de8a18e6ebc7abcaf5c8dec5920a72f10326797857e0b70eeab5f6`.
  - Its device code is identical to `build/sca_gpu_bias_tg`: `cuobjdump -sass`, 136,608 instruction lines, byte-identical. So the BV verification and BT timing of Amendment 1 apply unchanged.
  - `src/run_bias_b2.sh` (`c9e3689d…`) and `src/hw_import.py` (`a466b63a…`) now use `build/sca_gpu_bias_tg2`; nothing else changes.
- **Then.** Phase B runs as in Amendment 3.

## Amendment 5 (8 October 2026, ~11:35 UTC): GPU-selected G-set schedules, before any data of this amendment

- **Why.** Phase G ran the study's E12 configurations. The study selected them for the V80: S*₁₂ = argmin over S of its V80 12-engine round estimator (cycle model 0.1424·flips + 18.09·S + 886 at 250 MHz, E[max of 12 trials]), evaluated on the 256-run finals of the per-budget mean-cut selections. This was checked for all 408 original and extended E12 choices, and likewise for E1. The same S is not optimal for the GPU's batch costs, so the earlier G-set comparison (identical schedules) is unfair to the GPU. Here each platform gets its own selection by the same procedure. The V80 side is its measured A3 board results (`data/v80_a3/a3_summary.json`, an unchanged copy of `fpga/v80_sca/results/multibit_20261007/board/board_mb2n_e12_250mhz/hwmb_board_mb2n_e12_250mhz_a3/a3_summary.json`, sha256 `8edaa7e6…`).
- **GPU cost model.** t(S, B, class) is the median measured device time per batch of the phase-G G-set runs, for B ∈ {60, 120, 240}, class ∈ {const (SCA, TEC, Onsager-kT), ons (Onsager-online)} and S ∈ {250, 500, 1000, 2000, 4000}.
  - All 30 cells are measured, with 7–47 runs each and a spread of about ±1%. The table is `results/GS/cost_model.json` (sha256 `0a9b3f5f…`).
  - Each value is the measured batch-time floor of the variant the driver uses for that cell (phase G's rule, including the register fallback at S = 4000). Examples (ms): B = 60 const 0.182, 0.346, 0.672, 1.326, 2.642; B = 240 const 0.503 … 9.377.
- **Selection** (`src/gsel.py`, sha256 `f1bb4a181b36328394341620f7631ae809deb8f00d836ee42a1e6c90fbe4d96b`). Per instance, rule (SCA, TEC, Onsager-kT, Onsager-online) and grid (original; extended = Amendment 1 union):
  - **Candidates.** The 15 pairs (S, B). The configuration at S is the study's per-budget selection, with its 256-run final (k of 256). These are exactly the data and candidates the V80's selection used.
  - **Primary.** P_B = 1 − (1 − p)^B; TTS = t if P_B ≥ 0.995, else t·ln 0.01/ln(1 − P_B). This is the study's round convention with B in place of 12.
  - **Secondary.** TTS = t·ln 0.01/(B·ln(1 − p)).
  - **Estimate of p.** p = k/256, with k = 256 counted as 255.5.
  - **Choice.** Argmin of each prediction separately; ties go to smaller S, then smaller B.
  - **Sufficiency.** All 408 (instance, rule, grid) have a candidate with k > 0. No GPU pilot is run: a pilot over all 32 pilot configurations per budget would give the GPU a wider selection than the V80 had.
- **Plan** (generated before any run; identical on the laptop and gpu-host). `results/GS/plan_labels.json` (sha256 `1f72b35194fccaf1aaf49fdef7b42e8493aac251a98346fc24df6b10b56f4408`) and `results/GS/plan_cohorts.tsv` (`d1e6d8c404857c9c58198b2f48cf1755168c6c5c11591c98c1b7e4dccd02d96f`).
  - 816 labels (51 instances × 4 rules × 2 grids × 2 estimators) in 499 unique cohorts; identical (configuration, S, B) cohorts are run once and shared.
  - Original grid, primary: the GPU's S equals the V80's E12 S in 102 of 204 cells. Mostly S = 250 at B = 60 (95), then B = 120 at larger S.
  - Original grid, secondary: mostly B = 240 at S = 1000–2000, and B = 120 at S = 4000.
- **Inputs.** The study's results, copied unchanged to `data/gset_study_results` (budget, ext_budget, conf, ext_conf; hash list sha256 `3b7c4a90…`). The V80 A3 power summary is copied unchanged (`data/v80_a3/power_summary.json`, `97900616…`).
- **Binary.** `build/sca_gpu_bias` (phases V–W, unchanged) with phase G's variant rule. All S ≤ 4000, so the step tables stay in shared memory.
- **CPU.** Every process is pinned to 16 cores (`taskset -c 128-143`); the checker uses 16 threads.
- **Runs** (`src/run_gsel.sh`, sha256 `5cd5f44e11af2a980b68754ef4a8846c75106e2ca3572132f6dead4c593638ca`), in order; seed 20261004.
  - **GSV, verification.** Every unique cohort is traced with 4 chains and `--check-fields 1`, then checked against `run_trial_bias` (spins, Σ s h, energy, flips, n_lin). Trial ids are 160,000,000 + 1000·k. Any mismatch stops the amendment.
  - **GSR, held-out.** Each unique cohort runs 128 batches of B chains (as phase G's E12 cells). Fresh trial ids: id_base = 110,000,000 + 200,000·c + 50,000·j, where c = 4·instance order + rule index and j = 2·grid + estimator, for the first label that uses the cohort.
  - **GSW, power.** One 40 s burn per (B, variant, class) combination present in GSR, using the first such cohort in plan order. There is 30 s of idle first and 20 s idle between burns. Trial ids are 170,000,000 + 1,000,000·i; outcomes are not read.
- **Analysis** (`src/make_gsel_tables.py`, sha256 `d8beb0d0def7333e0fde9b12eea6a2c69cdfad5d04852f99a70f85663cf9d835`; estimators `src/analyze.py`, unchanged). Original grid is compared with the V80's original-grid cohorts (the paper's Table 8a); extended with extended.
  - **GPU values.** GPU primary = primary TTS99 of the primary-selected cohort; GPU secondary = secondary TTS99 of the secondary-selected cohort.
  - **V80 values.** V80 primary = A3 `tts99_ms`; V80 secondary = A3 `tts99_ind_ms`.
  - **Aggregates.** Geometric means over the 51 instances per rule, V80/GPU ratios, and per-instance win counts (point estimates). Also the same per class (random G1–10, G22–31, G43–47; toroidal G11–13, G32–34; planar G14–21, G35–42, G51–54), and the best of the four rules per instance on each platform, labelled post hoc.
- **Stated in advance (fairness of the secondary).** A3 measured the V80 only at its E12 (primary-selected) schedules. The V80's own secondary-optimal schedule is the argmin over S of t_round12·ln 0.01/(12·ln(1 − p)) on the same finals. It is reported as a labelled **modelled** column: the study's float-model p (1,024-run confirmation where that S is an E1 or E12 confirmation, else the 256-run final) with the study's t_round12, not board-measured. The earlier identical-schedule numbers are reported again for reference.
- **Energy.** GPU: board power of the cohort's (B, variant, class) burn × TTS. V80: the mean A3 board power of its G-set loads × TTS.
- **Correction before data.** The first freeze of this amendment (sha256 `37521d6b…`, 11:32:30 UTC) lost the A3 file path from the "Why" item, because of a shell substitution while writing. Only that path is restored here. No run of this amendment existed at either freeze.

## Amendment 6 (8 October 2026, ~11:45 UTC, during GSV, before any analysis output)

- **What.** Amendment 5's reference table of the earlier identical-schedule comparison paired the phase-G cells with the V80's original-grid A3 cohorts. Phase G ran the study's extended-grid (extended_amendment1) E12 configurations, so the identical schedules are the V80's extended-grid cohorts. A pre-analysis check found the mismatch: only the extended pairing reproduces the paper's figures (4.0 and 1.3).
- **Change.** `src/make_gsel_tables.py` (new sha256 `316ba4cd10d9162e1944934b1e8f656d57b381678be82f3eb2867379d49f521c`; the previous version is kept as `src/make_gsel_tables_v1.py`, sha256 `d8beb0d0…`) pairs that table with the extended-grid cohorts. It shows two B conventions: the per-instance best of B, and the fixed B = 60 (primary) / B = 240 (secondary) that reproduces the paper's figures.
- **Unchanged.** The main comparison (original against original, extended against extended), the selection, the plan and all runs.
