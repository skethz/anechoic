# STATICA K2000 TTS: reproduced

3 October 2026. Source: Yamamoto et al., IEEE JSSC 56(1), 2021, DOI 10.1109/JSSC.2020.3027702 (Algorithm 1, Eq. 7–9, Table I–II, Figs. 13, 19, 22, 24, 25). The protocol in [PROTOCOL.md](PROTOCOL.md) was frozen before any run. Laptop CPU only; no remote server was used.

## Bottom line

- **The published long-run point is reproduced** by the plain algorithm: S = 1,560, T_init = 40, 0.48 ms, P_a = 0.77, TTS = 1.50 ms.
- **The published short-run point is reproduced only once STATICA's documented random-number circuit is modelled.** This was a post-hoc diagnostic. The point is S = 560, T_init = 30, 0.13 ms, P_a = 0.07, TTS = 8.23 ms.
- With uniform random thresholds, the same SCA rule performs better than the chip simulator at low T_init.
- **The published annealing times are linear in our simulated flip counts.** All 21 digitized Fig. 22 times are fitted to within 1.7 µs.

**Reference for later comparisons:** the published long point, `t = 0.48 ms, P_a = 0.77, tts(0.99) = 1.50 ms`. It is backed by this reproduction. Its own 100-run sampling interval is [1.17, 1.97] ms.

## Table II points

Each point uses 4,096 independent trials. Every final state was rescored from the raw edge list and all agreed.

| Point | Source | P_a [95% CI] | Mean cut | tts(0.99), paper's t |
|---|---|---:|---:|---:|
| Long, S=1560, T_init=40 | **Paper** (100 runs) | **0.77** [0.675, 0.848] | **33,073** | **1.50 ms** |
|  | SCA, uniform thresholds | 0.827 [0.815, 0.839] | 33,089.9 (+1.8σ) | 1.26 [1.21, 1.31] ms |
|  | SCA + STATICA folded threshold (post hoc, 2,048 runs) | 0.806 [0.788, 0.823] | 33,084.4 (+1.2σ) | 1.35 ms |
| Short, S=560, T_init=30 | **Paper** (100 runs) | **0.07** [0.029, 0.139] | **32,750** | **8.23 ms** |
|  | SCA, uniform thresholds | 0.138 [0.128, 0.149] | 32,831.6 (+5.3σ) | 4.02 [3.70, 4.37] ms |
|  | SCA + STATICA folded threshold (post hoc, 2,048 runs) | 0.081 [0.070, 0.094] | 32,766.8 (+1.0σ) | 7.08 ms |
|  | SCA, quantized uniform thresholds, no folding (post hoc) | 0.120 [0.106, 0.135] | 32,827.9 | 4.68 ms |

σ is the standard error of a 100-run mean, estimated from our trial spread.

**Why the short point needs the circuit model.** Fig. 11 / Sec. IV-C2 describe how each spin's random threshold is drawn:
- A 16-bit random value is ANDed with `T_MASK`.
- If it exceeds 2T, it is XOR-folded back into range, and then a random sign is attached.

Folding doubles the density of small thresholds, which suppresses small uphill moves. Quantization alone, with no folding, leaves quality at the plain-algorithm level. Folding brings both Table II points into agreement with the paper.

The diagnostic's unstated assumptions are:
- Q8.8 fixed point for T and the threshold.
- A sign bit independent of the masked bits.
- Ideal 16-bit inputs instead of XorShift32.

## Experiment-3 grid (Fig. 22 and 24)

There are 21 configurations with 1,024 trials each. Uniform thresholds were used, as in the frozen protocol.

| T_init | Mean energy vs Fig. 22 | Success vs paper's implied P_a |
|---|---|---|
| 50 | Agrees: within about ±60, roughly 2σ, with mixed signs | Close at S ≥ 960; ours lower at 560 and 760 |
| 40 | Ours lower by about 30–85 (1–3σ) | Ours slightly higher |
| 30 | **Ours lower by 167–244 (5–8σ)** | **Ours much higher** |

The paper's implied P_a is back-computed from the digitized Fig. 24 TTS and Fig. 22 times. The low-T_init gap is the same effect as at the short point. See [reproduction.png](reproduction.png) and [results.json](results.json).

**Folded-threshold model on the full grid** (post hoc; 1,024 trials per point; [hw_rng_grid.json](hw_rng_grid.json)). Mean energy minus Fig. 22:

| T_init | Uniform thresholds | Folded threshold | Standard error of a 100-run mean |
|---|---|---|---|
| 50 | −33 to +60 | −55 to +31 (all within 2σ) | 19–32 |
| 40 | −86 to +3 | −61 to −6 (within 2.3σ) | 21–33 |
| 30 | −244 to −167 | **−121 to −30** | 25–37 |

The circuit model brings T_init = 50 and 40 into agreement and halves the T_init = 30 gap. A 2–4σ residual remains at T_init = 30 for every S except 560. Other unmodelled chip details are plausible causes: XorShift32 output shared between adjacent spins, the exact fixed-point formats, and where the sign bit comes from. Not established.

## Annealing-time model

The chip's time comes from the authors' cycle-accurate 2K-spin simulator, which is unavailable. The DDSS architecture makes time depend on the number of flipped spins. Fitting the digitized Fig. 22 times against our mean flips per run gives

`cycles ≈ 0.434 × total flips − 0.98 × S + 1989`  (300 MHz)

with a maximum residual of 1.7 µs over all 21 points, which is within digitization error. That's about 2.3 flipped spins per cycle and essentially no per-step overhead. A simpler pre-registered model, `N/2 + Σ ceil(F/2) + c0·S`, fitted to the two Table II times, predicts the other 19 within 7%. Its fitted c0 is negative, so the true per-flip cost is below half a cycle.

For the long point the model gives **0.478 ms**, against the published 0.48 ms. Its work is 331,504 coupling-row operations per run: N for the initial field plus 329,504 flipped-spin row updates, all ±1 additions.

## Checks and files

- The vectorized SCA matches a literal transcription of Algorithm 1 exactly: 6 random 12-spin graphs, 40 steps, shared uniforms.
- The graph SHA256 is verified. Σw = −1,040, so the target H = −67,040, matching the paper.
- Final states of both Table II points are saved bit-packed in `final_states_*.npy`.
- Files: [statica_repro.py](statica_repro.py), [analyze.py](analyze.py), [hw_rng_sensitivity.py](hw_rng_sensitivity.py), [raw_results.json](raw_results.json), [hw_rng_sensitivity.json](hw_rng_sensitivity.json), plus logs.
- The anaconda SciPy install in this environment is broken, so exact binomial intervals are computed by bisection in `statica_repro.py`.
