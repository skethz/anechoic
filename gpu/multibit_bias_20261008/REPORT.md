> **Note (9 October 2026).** The final GPU G-set tables are `results/GS_tables.md`. They supersede the numbers that this report calls final.

## GPU multi-bit kernel with bias and runtime n (GH200): final report

Bias (external field) and runtime-n support are done on the K-bit GH200 kernel and match `run_trial_bias` bit for bit. The G-set Max-Cut measurements (K = 2, all 51 instances, the software study's schedules) are done with TTS99, power and energy. Everything finished by about 05:05 CEST.

Folders:
- Local: `gpu/multibit_bias_20261008/`
- Remote: `/scratch/USER/anechoic_gpu_multibit_bias_20261008`

The full numbers are in `results/tables.md` and `results/tables.json`.

## What changed in the kernel

The new kernel is `src/sca_gpu_bias.cu`; the multibit kernel is untouched.
- **Bias:** fields start at the bias, h(0) = b + J s(0). The per-step update is unchanged.
- **Runtime n:** the kernel already took n at runtime. It is now verified for n = 1, 17, 257, 800, 1000, 1999 and 2048.
- **Overflow guard:** it now includes max|b|. Bias files are limited to |b| ≤ 32767 (the spec's 16-bit bias).
- **WIDE variant:** used only when int32 decisions could overflow. It clamps s·h inside the decision, which is exact because large |s·h| already fixes both the flip test and the n_lin band test. Decisions are therefore exact for any |h| < 2³¹; no instance needed 64-bit arithmetic.
- **Energy:** the host reports Σ s·h, E = −(Σ s·h + Σ b·s)/2, and checks it against an independent host energy. An option compares every final field with a recomputation from the final spins.
- **Long schedules:** a second binary reads the step tables from global memory, needed for the TSP budgets of S = 8192 at K = 4. Its last build (tg2) fixed a host-only range-check error; its device code is byte-identical to the first build.
- **Tools:** an instance generator (random K-bit with bias, graph partitioning, TSP, brute force), the bias checker, a G-set plan builder from the study's file, and an independent GP/TSP evaluator with its own TSPLIB parser. The parser reproduces the published optimal tour lengths of gr24, fri26, bayg29 and bays29.
- **K = 8:** not built. The clamped decision would be exact, but 4 MB of coupling planes don't fit an 8-SM cluster.

## Bit-exactness against `run_trial_bias`

Checked: final spins, Σ s·h, Σ b·s, energy, flips and every n_lin value. Device/host and final-field mismatches are 0 in every run.

| Check | Trials exact |
|---|---|
| Bias study verification (335 runs) | 6,308 / 6,308 |
| — of which K2000 with b = 0, identical to the multibit kernel's raw files (no bias, all-zero bias file, WIDE forced; all 7 variants) | 3,304 / 3,304 |
| G-set validation before counting, n = 800, 1000 and 2000 (816 runs) | 3,264 / 3,264 |
| Global-table binary, including S = 5000–8192 | 796 / 796 |
| — of which K2000 identical to the multibit kernel | 336 / 336 |
| GP/TSP export schedules | 240 / 240 |

The biased instances cover random K = 2 and K = 4 matrices with n < 2000, two WIDE instances, and generated graph-partitioning and TSP cases. The small GP and TSP cases were brute-force checked: their ground states are exactly the 19 optimal partitions and the 8 optimal tours.

## Step time

- **The bias itself costs nothing per step.** Same binary, same J and trial ids, with bias against an all-zero bias file: ratio 0.996–1.001.
- **Rebuilt kernel against the multibit build:** not exactly zero. Per-step slopes differ by −3.0% to +3.6% depending on the variant:
  - +2.5% to +3.6% for the register variant at B = 1 and 60;
  - −3.0% for the shared-memory variant at B = 120;
  - 1.000 at B = 240.

  The instruction mix is the same (identical tensor-core counts; the 8–48 extra instructions are start-up bias loads), so this is the compiler rescheduling code, not per-step work.
- **WIDE variant:** +1.4% to +6.1%.
- **Global step tables:** −3% to +5% for K = 4; −15% to +37% for K = 2.

## G-set Max-Cut, K = 2

Setup:
- 51 instances; four rules (SCA, TEC, Onsager-kT, Onsager-online); APC-SCA omitted because the engine has no such mode.
- The study's `extended_amendment1` single-engine (E1) and 12-engine (E12) configurations.
- E12 run at B = 60, 120 and 240 (128 batches each); E1 at B = 1 (256 launches).
- About 11.0 million trials; device and host cuts agree in every run.

Results:
- **Coverage:** every one of the 204 E12 cells reaches the target at B = 120.
- **Agreement with the study's float model:** median |Δp| is 0.009 (E12) and 0.011 (E1), and the mean shift is about zero. 26 of 408 cells (6.4%) fall outside a two-proportion test at 95%, against 5% expected. This assumes the study used 1,024 trials per cell (every p × 1024 is an integer).
- **Quality:** the best cut equals the best-known value on 16 of 51 instances.

Geometric means (ms) of the per-instance best over rules and B. This best is a selection among 12 held-out cells. The V80 columns are the software study's modelled times, not measurements.

| Class | GPU primary | GPU secondary | V80 E12 model | GPU single chain | V80 E1 model |
|---|---:|---:|---:|---:|---:|
| Random (N = 800/1000/2000) | 0.18–0.27 | 0.005–0.034 | 0.023–0.10 | 0.32–1.7 | 0.04–0.49 |
| Planar | 0.34–1.52 | 0.057–0.84 | 0.080–1.78 | 2.2–66 | 0.34–19.8 |
| Toroidal | 1.58 / 3.69 | 0.68 / 3.23 | 0.69 / 4.42 | 56 / 208 | 7.4 / 58.8 |
| **All 51** | **0.42** | **0.052** | **0.108** | **2.9** | **0.52** |

The GPU's secondary TTS beats the V80 12-engine model on 48 of 51 instances; its primary TTS beats it on only 6. The primary TTS is set by the GPU's batch-time floor: about 0.18 ms for 60 chains at S = 250.

**Power and energy** (GPU1, 200 ms sampling):
- G22 at B = 120 draws 368 W on the board; idle is 103 W.
- Energy per solution:

  | Basis | Energy |
  |---|---:|
  | Primary TTS99 | 155 mJ (112 mJ above idle) |
  | Secondary TTS99 | 19 mJ |
  | Single chain | 1.06 J |

- Power was measured on one configuration; other full waves drew 316–368 W.
- Module power (531 W) was inflated by another study's CPU load on the shared Grace CPU, so I used board power.

## GP/TSP (functional checks only, as you asked)

- **Import:** 90 files hash-checked. For all 15 instances, my rebuild of the K = 4 matrix and bias from the original G-set graph or TSPLIB distances equals the export exactly. All 60 schedules give tables identical to `mb_common.hpp`'s.
- **Counted runs:** no cell reaches its target. GPP partitions are balanced in 3–9% of runs; TSP tours are valid in 14–98% of runs, with quality 0.06–0.58. This matches your software result.
- **Export inconsistency.** Running the GPU on the export's own seeds and trials 0–1023 reproduces its `precision_study` numbers exactly in 44 of 60 cells. The other 16 are GPP cells where the export claims 100% feasibility (P(target) up to 0.996). On trials 0–31 of each, the GPU is bit-exact with `run_trial_bias` (512/512), and the reference itself gives 3–10% feasibility. Those 16 exported entries don't match the golden reference; the study's own `export_check.json` didn't exist when I copied the export.

## Protocol and deviations

- **Protocol:** frozen before any data (`a47bf701…`, 01:21:51 UTC). Amendments were each hashed before the data they govern:

  | Amendment | Hash | Time (UTC) | Content |
  |---|---|---|---|
  | 1 | `5cde8e9b…` | 01:53:29 | global-table binary |
  | 2 | `b9dd85d3…` | 02:11:25 | power-phase fix |
  | 3 | `2bd86f5a…` | 02:40:10 | GP/TSP |
  | 4 | `2cfd3e59…` | 02:43:15 | host range-check fix |

  The approximate times written inside amendments 3 and 4 are slightly off; the hash files are authoritative.
- **Power phase:** its first start failed on my bug (the schedule line was read without the tab separator, so no load ran). The attempt is kept as `results/W_attempt1`, and the phase was rerun.
- **GP/TSP import:** its first start failed on a host range-check error for one cell. No GPU run had been made; the fixed binary has byte-identical device code.
- **Timing:** the rebuilt kernel's per-step time is within −3% to +4% of the multibit build, not exactly zero (see Step time).
- **GP/TSP runs:** they ran completely, as pre-registered, before your message arrived; they are reported only as functional checks.
- **Cache writes in the home directory:** `/home/USER/.ccache` was touched at 03:18–03:30 CEST, but not by this study. ccache isn't on its PATH, its compilers resolve directly to gcc-13, and its build times don't coincide.
- **Rules compliance:**
  - all remote work stayed in scratch, with GPU1 selected by UUID;
  - the laptop ran at most one process at a time;
  - the manuscript was not edited;
  - from `research/gset_20261007` I read only `selected_configs.json`;
  - `research/reaim_benchmarks_20261008` was only read.

## Files

All in `gpu/multibit_bias_20261008/`:
- `PROTOCOL.md`, plus `PROTOCOL.sha256` and `PROTOCOL_A1`–`A4.sha256`;
- `src/`: the kernels (`sca_gpu_bias.cu`, plus the global-table versions tg and tg2), the checker, the generator, the G-set and GP/TSP drivers, the evaluator and the analysis script;
- `results/`: `tables.md` and `tables.json`, per-phase summaries, per-trial CSVs, and `raw_manifest.sha256`. Raw spins and traces (4.9 GB) stay on gpu-host;
- `remote_meta/`: logs, build hashes, instance and export hashes.
