# GPU v2 (GH200) bit-exact Onsager/TEC-T SCA: measurement protocol

Frozen on 7 October 2026 (UTC), before any timed or counted run of the final binary. Changes after freezing go into dated amendments at the end. Earlier development runs (`tmp/dev*`, `tmp/prof*`) are not results and are not used.

## 1. Platform and rules

- Host `gpu-host`, GPU `GPU-xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx` (GH200 144G HBM3e, sm_90, 132 SMs) only.
  - `CUDA_VISIBLE_DEVICES` is set to the UUID (`environment.sh`).
  - The binary refuses to run unless exactly one device is visible and its `cudaDeviceProp.uuid` equals this UUID.
- Project directory: `/scratch/USER/snowball_gpu_v2_20261007`. Caches and temporary files are kept there; nothing is written under `/home/USER`.
- GPU1 may be shared. Utilisation, memory, board and module power, and SM clock are logged before and after every run (`logs/gpu_status.log` and each run's `summary.json`).
- Graph: `data/K2000.bin` (SNOWGPU1), SHA256 `a3058db6cffdbb632bc46223958cc7f8dce17da22eb2e4bbed6f60f811e6c0e9`. The binary checks `sumw == -1040`.
- Seed: 20261004. Success means a returned final-state cut of at least 33,000, rescored on the host from the final spins (`scalar_cut`, as in the V80 host). The device-side cut must also equal the host rescoring.

## 2. Implementation under test

Source hashes (sha256):

| File | sha256 |
|---|---|
| `src/sca_ref.hpp` (unchanged copy of `fpga/v80_sca/src/sca_ref.hpp`) | `24d04a41097ab21a4ef20b86217f48a56794cc470c2d7ab963566192ea582cab` |
| `src/tables.hpp` (V80 host table code, verbatim) | `f180577219494dcf26a7606e83a2aad281f51e57f5d8dca1bd95076cd4e3370a` |
| `src/sca_gpu.cu` | `9c79619c28ce8500f38b24a5960a343941bf333b540e0b39686fcfb320a467d9` |
| `src/check_ref.cpp` | `b41db0e3a6315917505f7bc36b53776da377832e0914fea8b1fba9220f8c9f39` |
| `src/analyze.py` | `32214f1e416702c779c38e0228bf2e1b976080dae2b82d7d24826e7e9ce9ffa0` |
| `src/run_protocol.sh` | `4476b5faae6c4ad05ea50d826bda96ebfa46c24290751d0e0a44ffc6f2d57b16` |
| `src/build.sh` | `d61b547a0dd84ec46d97372110da45b4822425845c5fcdb821b49f952d3677d0` |

The build is `bash src/build.sh` (nvcc 13.3, `-O3 -arch=sm_90 -DSCA_LANES=256`; g++ 13.3 for the checker). It produced `build/sca_gpu` with sha256 `5c9935afa123a478028731ea325ee7f31533548fc1c2beb291f0974a5c9a6862`. No kernel spills registers.

**Kernels.**

- **Cluster kernel**, template parameters CS ∈ {4, 8} (CTAs per cluster) and G ∈ {1, 2} (groups of 4 chains per cluster).
  - The cluster holds J as register-resident binary tensor-core fragments (`mma.m16n8k256.b1.and.popc`). Each warp owns 8 RNG lanes × 8 rounds: the V80 slots k·256 + l.
  - Per step, each thread makes its 8 decisions with the reference's Q16.16 integer arithmetic. Flip and flip∧new-spin masks plus counts go to every CTA by `st.async` with mbarrier `complete_tx`.
  - Every warp then computes popc(J_y ∧ F∧σ′) and popc(J_y ∧ F) for its 64 field slots and 4 chains, and updates the shifted field h − s exactly.
  - The random stream is generated ahead, while messages are in flight.
  - Mode variants: Onsager (correction from the previous n_lin), constant correction (TEC-T, TEC, plain: corr = kconst[t] folded into the thresholds), and n_lin counting (on in Onsager mode or when tracing).
- **Single-CTA comparison kernel** (`--cs 1`): one CTA of 256 threads per chain, J^T bit-packed in global memory (L2). It uses the brief's popcount update over non-zero flip words.
- **Chains per launch.** A cluster always computes 4·G chains. A launch reports `--chains B` of them (chain index < B). For B = 1, one cluster computes 4 chains (G = 1) and reports chain 0.

## 3. Phase V: validation (must pass before any other phase)

Run `bash src/run_protocol.sh V`. Every trial is re-run on the CPU with `sca::run_trial(trace=true)` from `sca_ref.hpp` (`check_ref`). Final spins, cut (device and host), total flips and every per-step n_lin must be equal. The runs:

| Config | Variant (CS/G, chains, first trial id) |
|---|---|
| plain q4 λ0 T30 S40 | 8/1, 4, 0; 4/2, 8, 900000; L2, 2, 77 |
| Onsager q8 λ1.05 ramp T12 S960 | 8/2, 8, 900100; 4/1, 4, 123456; L2, 2, 901234; 8/1 full wave 120, 940000 |
| TEC-T q8 κ1.75 ramp T15 S280 | 8/1, 4, 2000000; 4/2, 8, 950000; 8/2, 8, 31337; L2, 2, 999; 8/2 full wave 240, 960000 |
| TEC q8 J_v −4 T30 S100 | 8/2, 8, 999990; 4/1, 4, 42; L2, 2, 905000 |

This is 424 trials. **Any mismatch stops the protocol.**

## 4. Phase K: kernel variant per batch size (timing only; the variants are bit-identical)

- **Runs.** X5, trial ids from 5,000,000 (distinct ranges), 16 batches each (4 for the slow single-CTA L2 kernel):
  - B = 1: 8/1, 4/1 and L2;
  - B = 120: 8/1 and 4/1;
  - B = 240: 8/2 and 4/2;
  - L2 at B = 264, its own wave.
- **Selection.** The variant with the lowest mean device time per batch is used for that batch size in every later phase. B = 480 uses the B = 240 choice (two waves).

## 5. Timing

- **t_batch.** CUDA events recorded on the launch stream immediately before and after one kernel launch of all B chains. All chains start together.
  - This includes in-kernel setup: loading J fragments from global memory, seeding, the initial field and all S steps.
  - It excludes host table and graph upload, which happen once per process.
- **Host wall time.** Launch to event synchronisation, recorded per batch (`host_ms`).
- **Warm-up.** One warm-up launch per process, with trial ids ≥ 0xF0000000, is discarded.
- **Batch sizes.** Co-residency comes from `cudaOccupancyMaxActiveClusters`, measured during development: 30 clusters for CS = 4 and CS = 8, i.e. 120 chains at G = 1 and 240 chains at G = 2.
  - B = 1: single chain.
  - B = 120: a full wave of G = 1 clusters.
  - B = 240: a full wave at the maximum co-resident chain count, G = 2.
  - B = 480: two waves, as a check.

## 6. Estimators (`src/analyze.py`, identical to the FPGA ones)

- **Primary.** P_batch is the fraction of batches with at least one success. TTS99 = t̄_batch · ln(0.01) / ln(1 − P_batch); if P_batch = 1, TTS99 = t̄_batch.
- **Secondary.** p is the per-trial success fraction. P = 1 − (1 − p)^B, and the same formula applies, i.e. t̄_batch · ln(0.01) / (B · ln(1 − p)).
- **Intervals.** Wilson 95% intervals on P_batch and p, mapped through the (monotone) formula. Whole-batch TTS (t̄_batch · ⌈ln 0.01 / ln(1 − P_batch)⌉) is also reported.
- **Per-step time.** Fitted as t_batch = a + b·S across the pilot cells for each kernel variant. Dependence on flips is tested by adding mean flips per trial as a regressor.

## 7. Phase M: the V80 configurations (same flags as the V80 host; t1 = 5)

X5 `--q 8 --tecT 1.75 --ramp --t0 15 --steps 280`; X1 `--q 8 --tecT 1.75 --ramp --t0 15 --steps 320`; X4 `--q 8 --tecT 2.0 --ramp --t0 15 --steps 760`; X2 `--q 6 --lambda 0.9 --ramp --t0 15 --steps 320`; O1 `--q 8 --lambda 1.05 --ramp --t0 12 --steps 960`; O4 `--q 6 --lambda 0.9 --t0 12 --steps 360`; P4 `--q 8 --lambda 0 --t0 30 --steps 960`.

For configuration index i (X5 = 0, X1, X4, X2, O1, O4, P4 = 6) the trial ids are:

| Batch size | Runs | First trial id |
|---|---|---|
| B = 1 | 256 single-chain launches | 10,000,000 + 1,000,000·i |
| B = 120 | 128 batches | 10,000,000 + 1,000,000·i + 100,000 |
| B = 240 | 128 batches | 10,000,000 + 1,000,000·i + 200,000 |
| B = 480 | 64 batches (two-wave check) | 10,000,000 + 1,000,000·i + 300,000 |

## 8. Phase P: GPU schedule search (pilot, exploratory)

- **Grid** (120 schedules, index order as in `run_protocol.sh`):
  - mode ∈ {TEC-T κ 1.5, 1.75, 2.0; Onsager λ 0.9, 1.05}, all with the end ramp;
  - q ∈ {6, 8};
  - T_init ∈ {12, 15};
  - S ∈ {80, 120, 160, 200, 240, 280}.
- **Cells.** Each schedule is a cell at B = 240 (32 batches, trial ids 3,000,000 + 8,000·index onward) and a cell at B = 120 (32 batches).
  - The B = 120 cell uses the first 3,840 trial ids of the B = 240 range. It is therefore a timing replay with identical per-trial outcomes. Its P_batch is evaluated on its own 32 batches of 120.
- **Selection rule (pre-specified).**
  - Primary: the cell with the smallest primary TTS99 evaluated at the **lower** Wilson 95% bound of P_batch. This is an upper confidence bound that penalises small-sample P_batch = 1.
  - Secondary: the smallest secondary TTS99 at the lower Wilson bound of p.
  - Ties go to the smaller point estimate.
- Pilot numbers are exploratory and are not reported as results.

## 9. Phase H: held-out confirmation

- **Selected primary cell** (schedule, B):
  - 128 batches at B, trial ids from 4,000,000;
  - 256 single-chain runs, ids from 4,100,000;
  - 128 batches at the other wave size, ids from 4,200,000.
- **Selected secondary cell**, if different: 128 batches, ids from 4,500,000.
- Only held-out numbers are reported for the GPU-selected schedules.

## 10. Phase W: power

- **Sampling.** `nvidia-smi --id=<UUID> -lms 200` logs power.draw, power.draw.instant, module.power.draw.average, module.power.draw.instant, utilisation and SM clock.
- **Sequence.**
  1. 30 s idle;
  2. 75 s of back-to-back batches of the held-out primary configuration (`--burn 75`, trial ids 6,000,000+, outcomes not read);
  3. 15 s idle.
- **Energy to solution.** Mean power during the burn × held-out primary TTS99. Both whole-module (GH200 module, including the Grace CPU) and GPU board power are reported, together with idle-subtracted values.

## 11. Known limits stated in advance

- **Single instance and target.** Final-state success only. A shared GPU (another user holds memory on GPU1, about 21 GB at freeze time) may add timing noise.
- **Not a matched-hardware comparison.** The V80 comparison numbers come from `fpga/v80_sca/RESULTS_HW.md` (12 × v6.2 engines at 275 MHz). Platforms, clocks and power envelopes differ.
- **B = 1 is a 4-chain cluster launch reporting one chain.** Its device time equals that of 4 chains.

## Amendment 1 (7 October 2026, after phases V, K and M had run, before any analysis or pilot run)

- **Error.** The frozen `analyze.py` evaluated the secondary estimator as `tts(1 - (1-p)**B, t)`. In double precision 1 − (1 − p)^B rounds to exactly 1 for large B (for example p ≈ 0.2, B = 240), and `tts` then applies the P = 1 floor and returns t̄_batch.
- **Why this is an error.** It contradicts the formula frozen in §6, t̄_batch · ln(0.01) / (B · ln(1 − p)), and the FPGA's secondary estimator, which has no floor (for example the V80's X4 secondary is 0.0254 ms against a 0.198 ms round).
- **Fix.** `analyze.py` now evaluates t̄ · ln(0.01) / (B · ln(1 − p)) directly (`tts_sec`), in the estimates and in the secondary pilot selection score.
- **Not changed.** The primary estimator and all run definitions.
- **Data.** No analysis output existed when this was fixed. Phases V, K and M produced raw files only.
- **Summary fields.** The per-run `summary.json` fields `tts99_primary_ms` and `tts99_secondary_ms` written by the binary are informational, and the secondary field has the same rounding issue. Reported numbers come from `analyze.py` only.
