# STATICA K2000 TTS reproduction: frozen protocol

3 October 2026. Written before any result was produced. The goal is to reproduce Yamamoto et al., IEEE JSSC 56(1), 2021, DOI 10.1109/JSSC.2020.3027702: Table II and Figs. 22, 24 and 25. It then becomes the reference for evaluating IAMP. No parameter is tuned.

## Algorithm, from the paper

- Algorithm 1 with the clipped-linear stay probability of Eq. (7): `P = clamp((h·σ + q)/(4T) + 1/2, 0, 1)`, where `h = Jσ`. A spin flips if `P < u` with `u ~ U[0,1)`. All N decisions use the old configuration, and all are committed together.
- `q` is fixed at 4.0 (experiment 3, Table I). `T_1 = T_init`, `T_{s+1} = r_T T_s`, `r_T = (T_fin/T_init)^(1/(S−1))`, so step S uses `T_fin`.
- Random initial spins and no external field. The output is the **final** configuration after S steps.

## Instance and target

Canonical WK2000_1, `fpga/v80_snowball/data/WK2000_1.rud`, raw-file SHA256 `9ed615e5…252ba7`. `J = −w` with zero diagonal. `H = −σᵀJσ/2`, and cut `C = (Σw − H)/2`. Success means final cut ≥ 33,000, equivalently `H ≤ −67,040`.

## Runs

- Experiment-3 grid: `S ∈ {360, 560, 760, 960, 1160, 1360, 1560}` × `T_init ∈ {50, 40, 30}`, `T_fin = 5`, with 1,024 independent trials each.
- The two Table II points get 4,096 independent trials each: short is `S=560, T_init=30` and long is `S=1560, T_init=40`.
- Each configuration has its own seed derived from a fixed root seed, and each trial has an independent initial state and random stream.
- Arithmetic is double-precision uniforms in float32 matrix arithmetic, which is exact for these integer fields. STATICA's 16-bit fixed-point T/q and its XorShift32 + T_MASK random-number circuit are not modelled. This is a recorded deviation, so the reproduction is of the SCA algorithm, not of the bit-level chip.

## Reference values, read from the paper

| Point | Mean H | Mean cut | P_a | t | tts(0.99) |
|---|---:|---:|---:|---:|---:|
| Short (Fig. 22/25, Table II) | −66,540 | 32,750 | 0.07 | 0.13 ms | 8.23 ms |
| Long (Fig. 22/25, Table II) | −67,186 | 33,073 | 0.77 | 0.48 ms | 1.50 ms |

Each paper value is the mean over 100 annealing runs. The exact 95% binomial intervals of the published success probabilities are [0.029, 0.139] for 7/100 and [0.675, 0.848] for 77/100. A point counts as **reproduced** if our success-probability estimate lies inside the paper's interval and our mean cut is within the sampling uncertainty of a 100-run mean. Values read off Figs. 22 and 24 are approximate digitizations and are used only for qualitative comparison.

## Timing model, for t

The paper's t comes from a cycle-accurate 2K-spin simulator at 300 MHz. That simulator is unavailable, so t is modelled from the published architecture:
- Step 1 computes all fields at two spins per cycle: N/2 cycles.
- Each later step accumulates the previous step's flipped spins at two per cycle (DDSS with two SRAMs): `ceil(F_{s−1}/2)` cycles.
- Each step also adds a constant pipeline overhead `c0`.

`c0` is fitted to the two Table II times and then used unchanged to predict the rest of the grid. The flips-only version (`c0 = 0`) is also reported. This is a model, not the authors' simulator.

## TTS

Eq. (9): `tts(0.99) = t · ln(0.01)/ln(1 − P_a)`, unrounded, as in the paper. Report it three ways: with the paper's t and our P_a, with modelled t and our P_a, and with exact binomial intervals for P_a.

## Checks

- A literal scalar implementation of Algorithm 1 must match the vectorized one exactly on small random graphs with shared uniforms.
- Every final state of the two Table II points must be rescored from the raw edge list, independently of the matrix.
