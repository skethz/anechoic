# G-set transfer test: Onsager-corrected SCA vs plain SCA, TEC and the K2000 comparison set (frozen 7 October 2026, before any G-set run)

## Question

All hardware results are on one dense instance (K2000, WK2000_1). Does the Onsager correction's advantage hold on the standard sparse G-set Max-Cut benchmark?

- **Model.** These are results of the bit-faithful software model of the engine, not hardware results.
- **Why not on hardware.** The V80 engine stores ±1 couplings only and cannot represent the zeros of a sparse graph.

## Instances and data

- **Source.** Y. Ye's G-set page, https://web.stanford.edu/~yyye/yyye/Gset/, downloaded 7 Oct 2026.
- **Instances.** All 51 instances with N ≤ 2048 (the engine limit): G1–G47 and G51–G54. G48–G50 (N = 3000) are excluded.
- **Verification.** See [manifest.json](manifest.json) and [verify_data.py](verify_data.py).
  - For every file: the edge count equals the header, weights are ±1, indices are in range, and there are no self-loops or duplicate edges.
  - SHA-256 is recorded, and two independent downloads are byte-identical.
  - The ±1 classes have exactly the edge sets of their +1 sibling classes, including G23, whose server date is 2023.
- **Mapping.** J_ij = −w_ij, H(s) = −Σ_{i<j} J_ij s_i s_j, cut(s) = (W − H(s))/2 with W = Σ w_ij. The cut is computed exactly in integers.
- **Classes** (12):

| Class | Instances | N | Weights |
|---|---|---|---|
| Random 6% | G1–G5 | 800 | +1 |
| Random 6% | G6–G10 | 800 | ±1 |
| Toroidal | G11–G13 | 800 | ±1 |
| Planar-like | G14–G17 | 800 | +1 |
| Planar-like | G18–G21 | 800 | ±1 |
| Random 1% | G22–G26 | 2000 | +1 |
| Random 1% | G27–G31 | 2000 | ±1 |
| Toroidal | G32–G34 | 2000 | ±1 |
| Planar-like | G35–G38 | 2000 | +1 |
| Planar-like | G39–G42 | 2000 | ±1 |
| Random 2% | G43–G47 | 1000 | +1 |
| Planar-like | G51–G54 | 1000 | +1 |

- **Representative subset** (pre-declared, 12 instances, the first instance of every class): G1, G6, G11, G14, G18, G22, G27, G32, G35, G39, G43, G51.

## Best-known cuts and target

Values are parsed programmatically from two published tables ([build_bkv.py](build_bkv.py), [bkv.json](bkv.json)); none are typed by hand. Local copies are in `sources/`.

- **A.** U. Benlic and J.-K. Hao, "Breakout local search for the Max-Cut problem", *Engineering Applications of Artificial Intelligence* 26(3):1162–1173, 2013, doi:10.1016/j.engappai.2012.09.001. Table 2, read from the authors' accepted manuscript.
- **B.** F. Ma and J.-K. Hao, "A multiple search operator heuristic for the max-k-cut problem", *Annals of Operations Research* 248:365–403, 2017, doi:10.1007/s10479-016-2234-0. Table 5 (best known in the literature, including parallel GES, and MOH), read from arXiv:1510.09156v1.

**Rule.**
- BKV = max(A.f_best, B.f_pre, B.MOH).
- **Target = ⌈0.99·BKV⌉**, the same convention as bSB's TTT99 and close to K2000's 33,000 = 98.99% of 33,337.

**Flags.**
- **G23.** A's "previous best" column lists 13,354. No algorithm in either table reports it, and B's best-known compilation gives 13,344. BKV is taken as 13,344 (target 13,211); success at 13,221 (99% of 13,354) is also reported.
- **Values improved after BLS (from B):** G30 3413, G31 3310, G35 7687, G36 7680, G37 7691, G38 7688.
- **Preprints.** Both tables were read from preprints; the published versions were not compared.

## Methods

These are the ten methods of the K2000 algorithm comparison (`research/algorithm_compare_20261007/PROTOCOL.md`), with the same step definitions.

| Index | Method | Implementation (this folder; the original modules are not edited) |
|---|---|---|
| 0 | SA | `others.sa`: `methods._sa_kernel` over CSR rows. Bit-identical decisions (verify_others.json). |
| 1 | SCA (plain, STATICA) | `engine.run`, family plain |
| 2 | TEC | `engine.run`, family tec: field h + J_v s(t−1) |
| 3 | APC-SCA | `engine.run`, family apc |
| 4 | ReAIM ASA (noise-free) | `others.reaim`: `pilot.Batch` with an exact sparse field update. Bit-identical. |
| 5 | aSB | `others.asb`: the Goto 2019 integrator of `methods.asb_traj` (M = 5) |
| 6 | bSB | `others.sb`: the integrator of `methods.sb_traj` |
| 7 | dSB | `others.sb`, family dSB. Bit-identical to `methods.sb_traj`. |
| 8 | Onsager-kT (ours) | `engine.run`, family tecT: field h − κ·T(t−1)·ramp(t)·s(t−1) |
| 9 | Onsager-online (ours) | `engine.run`, family onsager: field h − λ_eff·ramp(t)·c(t−1)·s(t−1), c = n_lin/(2T) counted online over all N spins |

**The engine model.**
- `engine.run` is `research/ablation_20261005/abl.py::run` (unedited) with two changes:
  - the final temperature T_fin is a parameter (abl.run hard-wires the K2000 value 5);
  - the field update is an exact sparse scatter.
- On 45 cases (5 sparse synthetic graphs × 9 configurations covering every family) the final spins and flip counts are **bit-identical to abl.run itself** ([verify_engine.json](verify_engine.json)).
- Schedule: T(t) = T0·(T_fin/T0)^(t/(S−1)).
- `ramp`: the coefficient falls linearly to 0 over the last 30% of steps, as on the hardware. It is on for both Onsager forms.

**Other methods.**
- **Normalisation.** bSB, dSB and aSB use Goto's normalisation with σ_J = sqrt(Σ J²/(N(N−1))): c0 = ξ·0.5/(σ_J√N) and ξ0 = ξ·0.7/(σ_J√N).
  - For K2000, σ_J = 1, which is what methods.py hard-wires.
  - Hard-wiring σ_J = 1 on sparse graphs would weaken the coupling by about sqrt(N/d).
- **Numerics.** bSB and aSB sum J·x in a different floating-point order than dense BLAS. They are statistically equivalent to methods.py (verify_others.json), not bit-identical.

## Scaling to sparse instances (stated before any G-set run)

Per instance:
- **Field scale.** σ = sqrt(d̄), where d̄ = (1/N)Σ_ij J_ij² = 2m/N is the mean degree. σ is the RMS local field under random spins. K2000 has σ_K = sqrt(1999) = 44.7.
- **Temperatures and field-valued parameters** are given in units of σ. This covers T0, T_fin, q, J_v, q_reset, q_lim, and SA's T0 and T1.
- **Onsager-online coefficient.** The lag-one echo coefficient of spin i is Σ_j J_ij² f′(g_j). For near-uniform degree it is about (d̄/N)·n_lin/(2T). The engine's global popcount n_lin/(2T) therefore over-states it by N/d̄ on sparse graphs.
  - Hence **λ_eff = λ·(d̄/N)/(1999/2000)**, with λ given in K2000 units; the grid is below.
- **Onsager-kT coefficient.** κ is dimensionless.
  - With T scaled by σ, the term κ·T·s already scales with the field.
  - The echo coefficient also scales as σ when T does: (d/N)·N·x/(2T) with T ∝ sqrt(d) gives ∝ sqrt(d).
  - So κ carries no extra d/N factor. Against the risk of over-correction, the κ grid extends 16× below the K2000 optimum (κ ≈ 2).

**Exploratory calibration.** This was done before the protocol, on **synthetic graphs of each class only, never a G-set instance** ([calib.py](calib.py), [calib2.py](calib2.py), [calib3.py](calib3.py), with logs and JSON).
- It located the shared SCA parameters for plain SCA, the baseline.
- **Key finding.** On all-positive-weight (+1) graphs, synchronous SCA with the K2000-scaled self-coupling (q ≤ 0.36σ) oscillates between the all-up and all-down states and ends at cut 0. Plain SCA needs q ≈ 1.1–2.2σ and T_fin ≈ 0.67σ there.
  - Any anti-echo term also damps this oscillation. So that the corrected rules cannot "win" only by curing an under-sized q, the shared q range is set per weight class so that plain SCA works.
- The ±1 classes behave like K2000 scaled: q ≈ 0.09–0.36σ, T0 ≈ 0.9σ, T_fin ≈ 0.11σ.
- Rule parameters were only placed: the λ and κ optima on synthetic graphs were 0.125–1, and TEC's preferred J_v sign depends on the weight class. Nothing below is selected from these runs.

## Grids (equal size: 16 points for every method)

Shared-parameter classes:
- **P+**: +1 random (G1–5, G22–26, G43–47).
- **Q+**: +1 planar-like (G14–17, G35–38, G51–54).
- **M±**: every ±1 instance (G6–13, G18–21, G27–34, G39–42).

T_fin/σ is fixed per class, as T_fin = 5 was fixed on K2000: **0.672** for P+ and Q+, and **0.112** (the K2000 value 5/44.7) for M±.

| Method | Grid (values in units of σ unless stated) |
|---|---|
| SCA (plain) | q × T0, 4 × 4. P+: q {1.08, 1.44, 2.16, 2.88}, T0 {0.9, 1.35, 1.8, 2.7}. Q+: q {0.72, 1.08, 1.44, 2.16}, T0 {0.45, 0.9, 1.35, 1.8}. M±: q {0.045, 0.09, 0.18, 0.36}, T0 {0.45, 0.67, 0.9, 1.35} |
| TEC | J_v × q × T0, 4 × 2 × 2. J_v: P+ and Q+ {−0.36, −0.18, −0.09, +0.09}; M± {−0.09, −0.045, +0.045, +0.09} (K2000's {−4, −2, 2, 4}/σ_K) |
| Onsager-kT | κ {0.125, 0.25, 0.5, 1.0} (dimensionless) × q × T0, 4 × 2 × 2 |
| Onsager-online | λ {0.125, 0.25, 0.5, 1.0} (K2000 units; λ_eff = λ·(d̄/N)/(1999/2000)) × q × T0, 4 × 2 × 2 |
| APC-SCA | q_reset × r_q {0.9, 0.97} × q_lim × T0, 2 × 2 × 2 × 2. P+: q_reset {2.88, 4.32}, q_lim {1.44, 2.16}. Q+: q_reset {2.16, 2.88}, q_lim {1.08, 1.44}. M±: q_reset {0.18, 0.36}, q_lim {0.045, 0.09} |
| SA | T0 {0.447, 0.671, 1.118, 1.789} × T1 {0.0224, 0.0447, 0.0671, 0.1118}, i.e. K2000's {20, 30, 50, 80} × {1, 2, 3, 5}/σ_K. T1 is shifted one step up from K2000's {0.5, 1, 2, 3} because calib3 found the upper edge preferred |
| ReAIM ASA | k set × T1 {0.2, 0.1, 0.05, 0.025}, 4 × 4. k sets {32..256}, {64..512}, {128..1024}, {256..2000} (doubling), scaled by N/2000 and rounded (min 1, max N); T from 1 |
| aSB | dt {0.25, 0.5, 0.9, 1.25} × ξ {0.5, 0.75, 1.0, 1.5} |
| bSB, dSB | dt {0.5, 0.75, 1.0, 1.25} × ξ {0.5, 0.75, 1.0, 2.0} |

The shared (q, T0) values for TEC, the two Onsager forms and APC-SCA:
- **P+**: q {1.44, 2.16}, T0 {1.35, 1.8}.
- **Q+**: q {1.08, 1.44}, T0 {0.9, 1.35}.
- **M±**: q {0.09, 0.18}, T0 {0.45, 0.9}.

These are subsets of plain SCA's 4 × 4 grid, so plain SCA has the densest search over the shared parameters, as on K2000.

## Procedure (identical for all methods)

- **Budgets.** S ∈ {250, 500, 1000, 2000, 4000} steps. The step definitions are those of the K2000 comparison: an SA sweep, an SCA parallel update, an SB integration step, a ReAIM iteration.
- **Pilot.** For each instance, method, S and grid point: 64 runs, seed SeedSequence([20261007, 1, g, method index, S index, grid index]), where g is the G-set number.
- **Selection (fixed).** Highest pilot mean final cut among valid configurations, ties to the lower grid index. This is the K2000 comparison's rule. A configuration is invalid if any run has non-finite dynamics (aSB).
- **Final, held-out.** The selected configuration, 256 runs, seed SeedSequence([20261007, 2, g, method index, S index]).
- **Success.** Final-state cut ≥ target, as on the hardware. The best-visited state is not used.
  - In a final or confirmation run, a run with non-finite dynamics counts as a failure and is excluded from mean-cut statistics; the count is reported.
  - If no pilot configuration is valid, that (method, S) has no final and p = 0.
- **TTS confirmation** (the four engine rules: plain, TEC, Onsager-kT, Onsager-online).
  1. For each rule, choose S*₁ = argmin over S of single-engine TTS99 on its 256-run finals. Choose S*₁₂ the same way on the 12-engine round TTS.
  2. Rerun each distinct (selected configuration, S*) with **1,024 fresh runs**, seed SeedSequence([20261007, 3, g, objective]). This seed is identical for all four rules, so runs are paired.
  3. Reported TTS values come from these confirmation runs only.

## Measures

- **Statistics.** p = successes/runs with a Wilson 95% interval. P(cut ≥ BKV) is also reported.
- **Steps to solution** (hardware-independent, all methods): MCS99(S) = S·ln(0.01)/ln(1 − p), or S if p = 1 and ∞ if p = 0, on the 256-run finals. The minimum over S is reported, as on K2000.
- **Engine time** (the four engine rules only; v6.4 cycle model).
  - Cycles per trial = 0.1424·flips + 18.09·S + 886, at 250 MHz. Flips are the mean measured flips of the runs.
  - **Assumption.** The engine stores dense coupling rows and processes 2,048 spin slots per step whatever the instance, so the model is unchanged for sparse graphs and for N < 2048 (padding slots are inactive and excluded from n_lin).
  - A ternary-coupling engine would need two population counts per field update but has the same cycle structure.
- **Single-engine TTS99.** t_trial·ln(0.01)/ln(1 − p) (Eq. 9).
- **12-engine round estimator** (primary-style, model-based; the convention of `research/v6_schedule_20261006/sweep_e12.py`):
  - t_round = E[maximum of 12 trial times], from the empirical per-run cycle distribution;
  - P_round = 1 − (1 − p)^12;
  - TTS = t_round if P_round ≥ 0.995, else t_round·ln(0.01)/ln(1 − P_round).
- **Ratios.**
  - Per instance: TTS(plain)/TTS(rule) and TTS(TEC)/TTS(rule), with a paired bootstrap 95% interval over the 1,024 confirmation runs.
  - Summary: geometric mean over instances where both rules are finite, with min, median and max; win counts per class; and counts of instances where a rule never reaches the target.

## Run order and scope

1. **Part A.** The four engine rules on all 51 instances; the 12-instance subset runs first.
2. **Part C.** TTS confirmations.
3. **Part B.** The six other methods (SA, APC-SCA, ReAIM, aSB, bSB, dSB) on the 12-instance subset.

Further notes:
- If time runs short, Part B is limited to the subset, as stated here. Any extension of Part B is declared by an amendment before it runs.
- Threads: 6 worker processes, each single-threaded. CPU only.

## Hypotheses (Parts A and C; engine TTS from the confirmation runs)

| ID | Hypothesis | Pass condition (all 51 instances) |
|---|---|---|
| G1 | Onsager-online beats plain SCA | Lower single-engine TTS99 on more than half of the instances. An instance where both are infinite counts as not lower |
| G2 | Onsager-kT beats plain SCA | Same |
| G3 | Onsager-online beats TEC | Same |
| G4 | Onsager-kT beats TEC | Same |
| G5 | The K2000-size gain transfers | Geometric-mean TTS ratio plain/min(online, kT) ≥ 2 on single-engine TTS and on the 12-engine estimator (K2000: 3.7× on hardware, 2.6–5.0× on fresh dense instances on CPU) |
| G1′–G4′ | The same as G1–G4 for the 12-engine estimator | Same |

Parts B and C are descriptive: steps-to-solution and final-cut statistics of all ten methods on the subset, and the frequency of reaching the BKV.

**Expectation, stated honestly.** The synthetic calibration showed only small mean-cut gains for the corrections on sparse graphs (at most about 0.1% of the cut), so G5 may well fail. Results are reported as measured.

## What this does not show

- No hardware runs. The V80 engine cannot store zero couplings.
- Steps are not time across methods. An SA sweep is N sequential updates; an SCA step is one parallel update.
- **Weight classes.** The tuning ranges differ by weight class: per-class shared q and T_fin, chosen on synthetic graphs. The +1 classes are a regime that the K2000 engine never faced.
