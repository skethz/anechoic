# ReAIM K2000 TTS: partial, exploratory reproduction

3 October 2026. Source: Chiang et al., ISCA 2024, DOI 10.1109/ISCA59077.2024.00015 (Algorithms 2–3, Tables II, III, VII, Sec. VI-D). Laptop CPU only. **These are exploratory pilot runs, not a frozen-protocol reproduction.** The k sets were extended after seeing early results.

## Bottom line

- **The published success probabilities are reproducible by the noise-free algorithm**, under the paper's own Max-Cut settings (T 1 → 0.1, F = max, Ising form), but only if the flip caps are large (k ≥ 64–128). At ~4,100 iterations noise-free ASA gives p ≈ 0.47; at ~6,300–6,500 iterations it gives p ≈ 0.80.
- **The published pair is internally consistent.** The iteration ratio between the two points (about 1.55–1.6) matches the published time ratio 0.23 / 0.15 = 1.53.
- **The published times cannot be checked.** They imply about **37 ns per iteration**, about 44 cycles at 1.2 GHz, each including a full 2,000 × 2,000 crossbar matrix-vector product. That is about 6× faster than the paper's own nurse-scheduling rate of 225 ns per iteration (2,048 iterations in 0.46 ms at WL = 32). It would need about 100 or more word lines active at once, which raises ReRAM noise. The paper does not give its K2000 MNSIM configuration, so this time is an unverified hardware projection.
- **ReRAM noise was not modelled.** Its parameters are not reported.

## What the paper specifies and what it does not

| Specified | Not reported for K2000 |
|---|---|
| Algorithm 2: thresholds, \|N\| vs q, random selection of ≤ k eligible spins, T ← αT; q = max of a 20-entry FIFO of \|N\| | Iterations per run (0.15 / 0.23 ms) |
| Algorithm 3: parallel k trial replicas, then a run phase with the best k; H_best updated after run phases | Candidate k set, ITER_trial, ITER_run, stopping rule |
| Max-Cut settings (Table II): T_init / T_final = 1 / 0.1; Ising 1-bit form; F = max for MCP (Sec. IV-B) | ReRAM noise parameters, WL_on |
| Hardware (Table III): 1.2 GHz, 256 × 256 1-bit crossbars, 16 × 8-bit 1.28 GS/s ADCs per IMA; K2000 scaled to 2,048 spins | Per-iteration latency from modified MNSIM; replicas per success event |

Our assumptions: ITER_trial = 32, ITER_run = 96, the FIFO starts at N (q starts large), and T follows the continuing trajectory.

## Results (noise-free, WK2000_1, best cut ≥ 33,000, 128–256 runs per point)

| Setting | 512 | 1,024 | 2,048 | 4,096 | 6,144 | 8,192 iterations |
|---|---:|---:|---:|---:|---:|---:|
| k ∈ {1, 2, 4, 8} (256 runs) | 0 | 0 | 0 | 0 | | |
| k ∈ {16, 32, 64, 128} | | 0 | | 0.188 | | |
| k ∈ {64, 128, 256, 512} | 0 | 0.055 | 0.188 | **0.469** | | |
| k ∈ {128, 256, 512, 1024} | 0.008 | 0.070 | 0.195 | **0.492** | **0.750** [0.67, 0.82] | **0.820** [0.74, 0.88] |
| k ∈ {256, …, 2000} | 0.016 | 0.047 | 0.148 | 0.406 | | |

The sweep also tried T 0.1 → 0.01 and T 0.3 → 0.003 with F = max or min. Every one of those settings had zero successes at 1,024 and 4,096 iterations, and F = min with T 1 → 0.1 also had zero.

**TTS check.** With the paper's times and our probabilities: 0.15 ms × r99(0.469) = **1.09 ms** (published 1.11), and 0.23 ms × r99(0.80) = **0.66 ms** (published 0.68).

## Interpretation

- **The ReAIM result is algorithmically plausible.** With large k, ASA behaves like a many-spin parallel update, and roughly 4–6 k iterations achieve STATICA-like success.
- **Its speed claim is a claim about analog matrix-vector throughput.** It needs about 4,100 full dense products per search. STATICA, by contrast, does about 166 product-equivalents of flip-sparse ±1 updates per search.
- Whether a ReRAM crossbar does a full 2,000 × 2,000 product with ADC readout in about 37 ns, at noise low enough to keep p, is the unverified part.
- Treat ReAIM's Table VII as an unverified hardware projection whose success rates are consistent with its algorithm. Confirming it needs the authors' K2000 settings, noise model and MNSIM configuration.

## Files

- [pilot.py](pilot.py): vectorized ASA. Fields were checked exact on a 40-spin graph; flip caps verified.
- [sensitivity.py](sensitivity.py), [sensitivity_bigk.py](sensitivity_bigk.py), and the `*_results.json` and `*.log` files.
- [ratio_check.json](ratio_check.json).
