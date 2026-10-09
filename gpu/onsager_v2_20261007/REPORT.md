# GPU v2 (GH200): bit-exact V80 SCA engine and measured TTS99

7 October 2026. gpu-host, GPU `GPU-xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx` (GH200 144G HBM3e, 132 SMs, SM clock 1980 MHz during runs). Remote directory: `/scratch/USER/snowball_gpu_v2_20261007`.

The protocol was frozen before any counted run: [PROTOCOL.md](PROTOCOL.md), sha256 `feb0785b…` (`PROTOCOL.sha256`), plus Amendment 1 (`PROTOCOL_A1.sha256`, an estimator rounding fix made before any analysis). All numbers come from `src/analyze.py` / `src/make_tables.py` applied to the raw files (`results/tables.md`).

This replaces nothing. The earlier, weaker GPU kernel ([../onsager_20261005/REPORT.md](../onsager_20261005/REPORT.md): 10–12 µs per step, O1 single-chain TTS99 14.9 ms) is the previous version.

## Summary

- **Bit-exact.** 424/424 GPU trials equal `sca::run_trial(trace=true)` (`sca_ref.hpp`, SCA_LANES = 256): final spins, cut, total flips and every per-step n_lin. This covers plain, Onsager, TEC-T and TEC; every kernel variant used for timing; trial ids up to 2,000,000+; and two full waves. In all 299 counted runs (2,013,448 trials), the device-side cut equals the independent host rescoring.
- **Per step.** About 0.63 µs (constant correction) or 0.72 µs (Onsager) for a single chain. 0.84 / 0.95 µs per step for a 120-chain wave, and 1.53 / 1.78 µs for a 240-chain wave. Batch time does not depend on the number of flips. This is about 15× faster per step than the previous GPU kernel.
- **Chains per wave.** 30 co-resident clusters of 8 CTAs: 120 chains with one 4-chain group per cluster, and **240 chains (the maximum co-resident)** with two groups.
- **Best held-out GPU TTS99, primary** (GPU-selected TEC-T κ2.0 + ramp, q8, T15, S120, B = 120): **0.114 ms**, 95% interval [0.114, 0.148]. The V80 (12 engines) measured 0.0686 ms, so **the V80 is 1.66× lower on the primary estimator**.
- **Best GPU TTS99, secondary.**
  - V80 configuration X4 at B = 240: **0.0078 ms** [0.0077, 0.0079], against the V80's best secondary of 0.0254 ms. **The GPU is 3.3× lower.**
  - The held-out GPU-selected secondary schedule gives 0.0112 ms.
- **Single chain (B = 1).** X5 TTS99 **2.36 ms** [1.89, 2.96]; O1 **1.31 ms** [1.12, 1.55]. The previous GPU kernel gave O1 14.9 ms.
- **Power.** 382 W board (power.draw) and 457 W GH200 module under sustained load, against 96 W / 172 W idle. Energy to solution at the held-out primary is about 43 mJ (board) or 52 mJ (module).

## 1. Design (`src/sca_gpu.cu`)

**Choice of field-update engine.** The suggested design (J columns in shared memory, popcount over flip words) was costed first and rejected.
- On sm_90, `POPC` issues at 16/clk/SM. With 100–900 flips per step, nearly all 64 flip words are non-zero for most of the schedule, so it would cost about 2,000 × 64 POPCs per chain-step.
- Instead the field update uses the **binary tensor core**. `mma.sync.m16n8k256.b1.and.popc` (SASS `BMMA.168256.AND.POPC`) is native on Hopper. A microbenchmark measured 0.54 instructions/clk/SM, i.e. 17.7k AND-popc bit-pairs per clock per SM. The fragment layout was verified against a CPU product.

**Data layout.**
- One cluster of CS = 8 CTAs × 4 warps (CS = 4 × 8 warps also built) holds the whole J (2048 × 2048 bits) in **registers**: 128 registers per thread as BMMA A-fragments.
- Each warp owns 8 RNG lanes × 8 rounds, the V80 slots k·256 + l.
- Lane L of a warp owns chain L mod 4 and RNG lane 8·warp + L/4. The 8 BMMA columns are 4 chains × (F∧σ′, F), where F is the flip mask and σ′ the new spins. So each lane receives popc(J_y ∧ F∧σ′) and popc(J_y ∧ F) for exactly the slots it decides.

**Exact field update.**
- The kernel keeps the shifted field ĥ = h − s. Then ĥ += 2·(4·cA − 2·cB + |F| − 2·|F∧σ′|), with no diagonal term.
- The initial field uses the same path with F := valid slots (factor 1).
- The reference's z = s·h·65536 + q ∓ corr becomes s·ĥ·65536 + (q + 65536) ∓ corr.
- A host guard proves for every step that the int32 evaluation cannot overflow (`check_ranges`).

**Decisions and random numbers.**
- Thresholds r = ((u16 − 32768)·4T) >> 16 are computed exactly, in 64-bit, from the per-lane xoshiro128** stream (splitmix64(seed, trial, lane) seeding).
- Draws for padding slots are consumed. The next step's 8 draws per thread are generated while messages are in flight.
- For TEC-T, TEC and plain (kcorr = 0), corr = kconst[t] is folded into the thresholds. n_lin is then only counted when tracing, since it does not affect the trajectory.

**Exchange.**
- Every warp packs its per-chain masks with one shuffle level. It sends 64 B of masks plus a 4-byte count word per chain to every CTA by `st.async … mbarrier::complete_tx`.
- Each CTA waits on a double-buffered mbarrier (2560 B per message). There is no CTA or cluster barrier inside the step loop. A buffer can only be overwritten after every warp of the cluster has sent the next message, so all reads of it have finished.

**Single-CTA comparison kernel (`--cs 1`).** One CTA per chain, J^T bit-packed in global memory (L2), with the brief's popcount update dh_y = Σ_w [F_w ≠ 0](2·popc(F_w ∧ m) − 4·popc((J^T_y,w ⊕ S′_w) ∧ F_w ∧ m)).

**Tables.** Per-step tables are built by `src/tables.hpp`, which is the V80 host code (`host_sca_multi_v6.cpp`: geometric T, λ ramp, `make_tables`, TEC-T override) copied verbatim.

**Development record.** `src/sca_gpu_v0_dev.cu`, `_v1_dev`, `_v2_dev` are preserved.
- v0 (CTA barrier and ballot packing) ran X5 in about 1.2 µs per step for a single cluster.
- The final kernel removes the barrier, packs with one shuffle level, uses 8 independent BMMA accumulator chains and overlaps the RNG with the message flight. No kernel spills registers (192–246 registers per thread).

## 2. Validation (phase V, before any timing)

| Configuration | Kernel variants (CS/G: chains) | Trials | Exact |
|---|---|---:|---:|
| plain q4 λ0 T30 S40 | 8/1: 4, 4/2: 8, L2: 2 | 14 | 14 |
| Onsager q8 λ1.05 ramp T12 S960 | 8/2: 8, 4/1: 4, L2: 2, **8/1 full wave 120** | 134 | 134 |
| TEC-T q8 κ1.75 ramp T15 S280 | 8/1: 4, 4/2: 8, 8/2: 8, L2: 2, **8/2 full wave 240** | 262 | 262 |
| TEC q8 J_v −4 T30 S100 | 8/2: 8, 4/1: 4, L2: 2 | 14 | 14 |
| **Total** | 15 runs, seed 20261004, trial ids 0 … 2,000,003 | **424** | **424** |

The success probabilities agree with the V80 for every configuration. For example:

| Config | p on GPU (B = 240) | p on V80 |
|---|---:|---:|
| X5 | 0.341 | 0.338 |
| X4 | 0.945 | 0.950 |
| O1 | 0.934 | 0.933 |
| O4 | 0.315 | 0.326 |

## 3. Kernel variants and chains per wave (phase K, X5)

| B | Variant | t_batch (ms) | µs/step |
|---:|---|---:|---:|
| 1 | **CS 8, G 1** | **0.186** | 0.66 |
| 1 | CS 4, G 1 | 0.235 | 0.84 |
| 1 | single-CTA, J in L2 | 13.06 | 46.6 |
| 120 | **CS 8, G 1** | **0.249** | 0.89 |
| 120 | CS 4, G 1 | 0.279 | 1.00 |
| 240 | **CS 8, G 2** | **0.444** | 1.58 |
| 240 | CS 4, G 2 | 0.466 | 1.66 |
| 264 | single-CTA, J in L2 (its own wave) | 20.18 | 72.1 |

- **Co-residency.** `cudaOccupancyMaxActiveClusters` gives **30 clusters** for CS = 4 and CS = 8, i.e. 120 SMs.
- **Waves.** One wave is 120 chains (G = 1) or **240 chains** (G = 2, the maximum co-resident). G = 3 does not fit in registers.
- **Choice.** CS = 8 was selected for every batch size. The single-CTA L2 kernel is about 70× slower at B = 1.

## 4. Per-step time and dependence on flips

**Fits over the 240 pilot cells** (S = 80–280; t_batch = a + b·S, and with mean flips per trial as an extra regressor):

| Cells | a (ms) | b (µs/step) | Flip coefficient | max \|residual\| |
|---|---:|---:|---:|---:|
| B = 120, constant correction (TEC-T) | 0.0141 | 0.839 | −0.025 ns per flip | 2.8 µs |
| B = 120, Onsager | 0.0149 | 0.952 | +0.001 ns | 3.1 µs |
| B = 240, constant correction | 0.0146 | 1.532 | +0.020 ns | 2.2 µs |
| B = 240, Onsager | 0.0152 | 1.779 | +0.001 ns | 5.2 µs |

- **Fixed cost.** The intercept, about 14–15 µs, is in-kernel setup: loading the 512 KB of J fragments, seeding and the initial field.
- **Single chain** (phase M): 0.633 µs/step with constant correction (X5, X1, X4) and 0.718 µs/step for Onsager (X2, O1).
- **No flip dependence.** Time does not depend on flips (less than 1% of batch time), because the tensor-core product is dense. This is unlike the V80 (about 0.27 cycles per flip) and the previous GPU kernel.

**Where a 120-chain step goes.** Diagnostic clock64 build, thread 0, `results/diag`, not a protocol result: 1,645 cycles (0.83 µs).

| Part | Cycles |
|---|---:|
| Receive, 32 BMMAs and field update | 619 |
| Waiting for messages | 319 |
| Next-step RNG (overlapping the flight) | 342 |
| Pack and st.async | 172 |
| Decisions | 104 |
| Loop | 90 |

- The cross-SM message round trip alone is about 250 ns (microbenchmarked; the cluster barrier is 175–330 ns).
- Nsight Compute shows issue slots 45% busy and about 2 active warps per scheduler. Dependency (wait) and scoreboard stalls dominate.

## 5. TTS99

Estimators as on the FPGA:
- **Primary:** t̄_batch · ln 0.01 / ln(1 − P_batch), where P_batch is the fraction of batches with at least one success, and t̄_batch if P_batch = 1.
- **Secondary:** t̄_batch · ln 0.01 / (B · ln(1 − p)).
- **Intervals:** Wilson 95%.
- **Time:** t_batch is the CUDA-event time of one kernel launch of all B chains.

### 5.1 V80 configurations on the GPU (phase M; trial ids from 10,000,000, distinct per configuration and batch size)

**Single chain, B = 1** (256 runs each):

| Config | p | t_run (ms) | TTS99 (ms) [95%] |
|---|---:|---:|---|
| X5 | 0.301 | 0.183 | **2.359** [1.894, 2.963] |
| X1 | 0.477 | 0.209 | 1.486 [1.247, 1.787] |
| X4 | 0.918 | 0.487 | 0.897 [0.770, 1.067] |
| X2 | 0.496 | 0.237 | 1.591 [1.340, 1.908] |
| O1 | 0.914 | 0.697 | **1.307** [1.123, 1.553] |
| O4 | 0.309 | 0.265 | 3.306 [2.662, 4.141] |
| P4 | 0.211 | 0.613 | 11.917 [9.170, 15.613] |

**Full waves.** 128 batches each. Every batch at B ≥ 120 had at least one success, so the primary TTS99 equals one batch time. Primary 95% upper bounds come from P_batch ≥ 0.971 (about +30%).

| Config | B = 120: t_batch = primary (ms) | B = 120: secondary (ms) | B = 240: t_batch = primary (ms) | B = 240: secondary (ms) | V80 primary / secondary (ms) |
|---|---:|---:|---:|---:|---:|
| X5 | 0.248 | 0.0231 | 0.443 | 0.0204 | 0.0686 / 0.0638 |
| X1 | 0.282 | 0.0195 | 0.504 | 0.0172 | 0.0771 / 0.0529 |
| X4 | 0.651 | 0.0085 | 1.180 | **0.0078** | 0.1984 / **0.0254** |
| X2 | 0.318 | 0.0193 | 0.584 | 0.0178 | — |
| O1 | 0.927 | 0.0130 | 1.723 | 0.0121 | 0.2352 / 0.0333 |
| O4 | 0.357 | 0.0354 | 0.655 | 0.0332 | 0.0695 / 0.0756 |
| P4 | 0.817 | 0.1332 | 1.487 | 0.1271 | — |

**Two-wave check, B = 480** (64 batches): batch times are 2.03–2.06× those at B = 240 (for example X5 0.901 ms, O1 3.519 ms). Secondary TTS is unchanged within its interval (X4 0.0081 ms).

### 5.2 GPU-selected schedule (pilot, then held-out)

- **Pilot (exploratory, trial ids from 3,000,000).** 120 schedules: TEC-T κ ∈ {1.5, 1.75, 2.0} and Onsager λ ∈ {0.9, 1.05}, all with end ramp; q ∈ {6, 8}; T_init ∈ {12, 15}; S ∈ {80, 120, 160, 200, 240, 280}. Each schedule ran at B = 240 and B = 120, 32 batches each.
- **Selection rule (pre-registered).** Minimum TTS99 at the lower Wilson bound of P_batch (primary) or of p (secondary).
- **Selected.**
  - Primary: **TEC-T κ2.0 + ramp, q8, T15, S120 at B = 120** (pilot point 0.114 ms).
  - Secondary: **TEC-T κ2.0 + ramp, q6, T15, S280 at B = 240** (pilot 0.0114 ms).
  - The next primary cells were TEC-T κ2 q6 S120 (0.115 ms) and Onsager λ0.9 q6 T15 S120 (0.128 ms).

**Held-out runs** (trial ids from 4,000,000, 128 batches):

| Run | B | p [95%] | P_batch [95%] | t_batch (ms) | TTS99 primary (ms) [95%] | TTS99 secondary (ms) [95%] |
|---|---:|---|---|---:|---|---|
| **Primary-selected** κ2 q8 T15 S120 | **120** | 0.041 [0.038, 0.044] | 1.000 [0.971, 1] | 0.1136 | **0.1136** [0.1136, 0.1479] | 0.1049 [0.0970, 0.1135] |
| same schedule | 240 | 0.044 [0.042, 0.047] | 1.000 [0.971, 1] | 0.1976 | 0.1976 [0.1976, 0.2573] | 0.0837 |
| same schedule | 1 (256 runs) | 0.035 [0.019, 0.065] | — | 0.0825 | 10.62 [5.61, 20.23] | — |
| **Secondary-selected** κ2 q6 T15 S280 | **240** | 0.531 [0.525, 0.536] | 1.000 [0.971, 1] | 0.4422 | 0.4422 | **0.0112** [0.0110, 0.0114] |

Both selections sit on the grid edge (κ = 2.0, and S = 280 for the secondary). The pre-specified V80 configuration X4 (κ2.0, S760) has a lower secondary (0.0078 ms) than the selected schedule, so the secondary optimum lies outside the pilot grid.

### 5.3 Comparison

| | V80 (12 × v6.2, 275 MHz, measured) | GPU v2 (this study) | Ratio |
|---|---|---|---|
| Primary TTS99, best | 0.0686 ms (X5) | 0.1136 ms (held-out, κ2 S120, B = 120) | **V80 1.66× lower** |
| Primary, same configuration X5 | 0.0686 ms | 0.248 ms (B = 120) | V80 3.6× lower |
| Secondary TTS99, best | 0.0254 ms (X4) | 0.0078 ms (X4, B = 240) | **GPU 3.3× lower** |
| Secondary, GPU-selected held-out | — | 0.0112 ms | GPU 2.3× lower than the V80's best |
| Per step, one chain | about 0.245 µs (X5 12-engine round / S) | 0.63–0.72 µs | V80 2.6–2.9× lower latency |
| Throughput | 12 chains per round, about 20 ns per chain-step (X5) | 240 chains, 6.4 ns per chain-step | GPU about 3.2× |

**Against the previous GPU kernel.** Different protocols and cohorts, so these are indicative only.
- Per step: 10–12 µs before, 0.63–0.72 µs now (about 15×).
- O1 single-chain TTS99: 14.9 ms before, 1.31 ms now.
- Best batch TTS: 4.8 ms frozen and 1.68 ms exploratory at B = 128 before, against 0.114 ms held-out now.

**Interpretation.** The GPU wins whenever the metric rewards aggregate throughput (secondary estimator, many chains per wave). The V80 wins on the primary (one-batch reliability) estimator. There the GPU's per-step latency floor, about 0.6–0.9 µs including a roughly 250 ns cross-SM message round trip per step, cannot be offset by its 10–20× larger wave.

## 6. Power and energy (phase W)

**Sampling.** nvidia-smi on the UUID every 200 ms. 30 s idle, then 75 s of back-to-back batches of the held-out primary configuration (710,208 batches, 0.1055 ms each), then 15 s idle.

| Window | Board power.draw (W) | Module power (W) | Utilisation | SM clock |
|---|---:|---:|---:|---:|
| Idle | 96.0 | 171.9 | 0% | 345 MHz |
| Load (primary configuration) | **381.6** | **457.4** | 100% | 1980 MHz |
| Idle after | 99.7 | 177.0 | 0% | 347 MHz |

The GH200 module figure includes the Grace CPU and memory. Energy to solution at the held-out primary TTS99 (0.1136 ms):

| Basis | Energy |
|---|---:|
| Board | 43.3 mJ |
| Board, idle-subtracted | 32.4 mJ |
| Module | 52.0 mJ |

Power was not measured separately for X4. At the same load level its secondary TTS would correspond to roughly 3 mJ (board), but that is an estimate, not a measurement.

## 7. Caveats and deviations from the brief

1. **Shared GPU.** About 21.7 GB on GPU1 was held by another process throughout, at 0% utilisation. Every run logged GPU state before and after (`results/remote_meta/logs/gpu_status.log`, `summary.json`). The clocks were not locked; the SM clock was 1980 MHz under load.
2. **Design.** The register-resident tensor-core design replaces the suggested smem-column popcount design (reasons in §1).
   - A cluster always computes 4 chains per group. **B = 1 is one 8-CTA cluster launch that reports one chain.** Its time equals that of 4 chains.
   - "Full wave" is 240 chains (30 clusters × 2 groups × 4). The 120-chain G = 1 wave is also reported because it is faster per batch.
3. **Mode-specific kernels.** For non-Onsager configurations the timed kernel does not count n_lin, which has no effect on the trajectory. The validation runs counted it, and it matched.
4. **Protocol details.**
   - The two-wave check uses 64 batches, not 128.
   - The pilot's B = 120 cells replay the first half of each B = 240 cell's trial ids (identical outcomes; a timing replay).
   - Amendment 1 corrected a double-precision rounding error in the frozen secondary estimator code. The formula itself was already frozen in §6. The error was fixed before any analysis.
5. **Selection on the grid edge.** Both selections sit at κ = 2.0, and the secondary at S = 280 as well (§5.2). The primary at P_batch = 1 is limited by statistics: with 128/128 successes, the 95% interval for P_batch is [0.971, 1].
6. **Comparison scope.** V80 numbers are copied from `fpga/v80_sca/RESULTS_HW.md`, not re-measured. Both systems use the same algorithm, instance, target, estimators and seed, but different trial ids and platforms. GPU timing includes in-kernel setup of about 14 µs per launch. It excludes host-to-device upload of the tables and graph (once per process).
7. **No energy comparison.** No V80 power figure was used, so no energy comparison is made.

## 8. Files

**Source** (`src/`):
- `sca_gpu.cu` (kernels and runner), `check_ref.cpp` (CPU checker);
- `tables.hpp` (V80 host table code, verbatim), `sca_ref.hpp` (unchanged);
- `analyze.py` (Amendment 1), `make_tables.py`, `run_protocol.sh`, `build.sh`;
- development versions `sca_gpu_v{0,1,2}_dev.cu`.

**Protocol.** `PROTOCOL.md`, `PROTOCOL.sha256`, `PROTOCOL_A1.sha256`.

**Results** (`results/`):
- `tables.md` and `tables.json` (all estimates);
- per-run `*.summary.json` and `*.batches.jsonl` in `V/`, `K/`, `M/`, `P/` and `H/`;
- `P/pilot_cells.json`, `P/selected_*.json`;
- `W/` (power samples, marks, burn record);
- `diag/` (Nsight metrics);
- `remote_meta/` (remote logs, build hashes, ptxas log, environment).

**Per-trial data.**
- `V/K/M/H_trials.csv.gz` hold the per-trial cut, device cut, flips and success.
- `raw_manifest.sha256` lists the sha256 of every raw file on gpu-host. That includes the final spin states (`*.spins.bin`, 492 MB) and the n_lin traces of the validation runs, which stay at `/scratch/USER/snowball_gpu_v2_20261007/results`.
