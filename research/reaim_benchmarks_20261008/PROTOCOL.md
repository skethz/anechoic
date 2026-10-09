# ReAIM benchmark suite at the software level: ten methods on Max-Cut, graph partitioning and TSP (frozen 8 October 2026, before any held-out run)

## Question

The paper's algorithm comparison (ten methods) is on Max-Cut. ReAIM (Chiang et al., ISCA 2024, IEEE 10609617) evaluates three problem classes: Max-Cut (MCP), balanced graph partitioning (GPP) and TSP. How do the ten methods compare on the same suite, under one identical tuning procedure?

- **Model.** Software models only. The SCA family runs the float engine model (`research/gset_20261007/engine.py`, bit-identical to `abl.run`); the coupling-precision study runs a bit-exact port of the hardware golden reference `fpga/v80_sca/src/sca_ref_bias.hpp`. No hardware run.
- **Not on today's hardware.** GPP and TSP need weighted (multi-bit) couplings and, for TSP, a bias (external field); the V80 engine stores ±1 couplings only.

## What ReAIM evaluated (verified in `research/benchmark_conventions_20260929/sources/reaim.txt`)

| Item | MCP | GPP | TSP |
|---|---|---|---|
| Benchmark (Table II) | G-set | G-set | TSPLIB [Reinelt 1991] |
| Instances (Table II) | G1–G20 (Table IV shows G1, G2, G6, G7, G10–G14, G19, G20) | G1–G5, G14–G17 | gr17, gr21, gr24, fri26, bayg29, bays29 |
| Size | 800 nodes | 800 nodes | 17/21/24/26/29/29 cities |
| Iterations (Table II) | 4096 | 4096 | 8192 |
| T_init/T_final (Table II) | 1/0.1 | 1/0.01 | 0.5/0.1 |
| Bit width QUBO/Ising (Table II) | 4/1 | 13/4 | 10/11 |
| Constraints m (Table I) | 0 | 1 | 2 |
| Metric | cut value, higher better; best of 20 runs (Sec. VI-A, following [4]) | cut value, lower better; best of 20 runs | ARPD (%), lower better: mean over 20 runs of 100·(L − L*)/L* (Sec. VI-A, following [4]) |

- Fig. 3(a) plots a "Normalized Cut Value (Higher is better)" for MCP; the normalization is not stated. Tables IV–VI report raw values.
- ReAIM does **not** state: the GPP or TSP Ising/QUBO formulas, the penalty weights, or how infeasible solutions (unbalanced partitions, invalid tours) enter the GPP cut or the ARPD.
- ReAIM's [4] is M. Ayodele, "Penalty weights in QUBO formulations: permutation problems", EvoCOP 2022 (arXiv:2206.11040, `sources/ayodele2022*.{pdf,txt}`). It uses the same six TSPLIB instances (Table 1 there), 20 runs and ARPD, and finds that MQC (penalty = maximum distance) gives the best ARPD for TSP with all runs feasible (its Tables 3–5).
- ReAIM's best GPP values (Table V) are kept for comparison (`gpp_reference.py`).

## Data

- **G-set:** the files and SHA-256 of `research/gset_20261007/data` (manifest there), unchanged. Best-known Max-Cut values and targets: `research/gset_20261007/bkv.json`.
- **TSPLIB:** `https://comopt.ifi.uni-heidelberg.de/software/TSPLIB95/tsp/` (`<name>.tsp.gz`, and `.opt.tour.gz` for gr24, fri26, bayg29, bays29), downloaded 8 Oct 2026, two byte-identical downloads. SHA-256, formats and checks: `data/tsplib_manifest.json` (`verify_tsplib.py`). All six are `EXPLICIT` (LOWER_DIAG_ROW: gr17, gr21, gr24, fri26; UPPER_ROW: bayg29; FULL_MATRIX: bays29), parsed per the TSPLIB 95 documentation; display coordinates are not used. GEO is implemented per the TSPLIB FAQ but not needed.
- **Optima verified three ways:** TSPLIB's optimum list (`sources/TSP-BEST.html`): gr17 2085, gr21 2707, gr24 1272, fri26 937, bayg29 1610, bays29 2020; the shipped optimal tours (4 of 6) evaluate to these lengths; and an independent exact integer program (HiGHS, subtour elimination) returns the same six values.

## Formulations (Lucas 2014, Front. Phys. 2:5, Secs. 2.2 and 7.1; ReAIM states none, so this is the pre-declared standard choice)

Engine convention: E(s) = −Σ_{i<j} J_ij s_i s_j − Σ_i b_i s_i, local field h_i = Σ_j J_ij s_j + b_i, flip cost 2 s_i h_i. All J, b are integers (`problems.py`).

- **MCP:** J = −w, b = 0 (exactly the G-set study).
- **GPP (min-cut balanced bisection):** H = A(Σ s_i)² + B Σ_{(ij)∈E} (1 − s_i s_j)/2, with **A = B = 1**. Integer form (×2): J_ij = [ij ∈ E] − P with **P = 4A = 4**, dense on all pairs; b = 0 (N even). E = 2·cut + (P/2)·M² + const, M = Σ s.
  - **Feasible** iff M = 0 (an exact 400/400 bisection).
  - **Penalty choice.** A = B = 1 is the only setting consistent with ReAIM's reported bit widths for GPP (Table II: 13-bit QUBO, 4-bit Ising): the QUBO diagonal is 4A(1 − N) + deg_i ≈ −3,150, which needs 13 bits only for A ≈ 1, and the Ising couplings {−4, −3}/2 then need 4 bits. This is an inference, not a stated value. Lucas's sufficient bound (A/B ≥ min(2Δ, N)/8, P ≥ max degree = 67–153) is not used: it would make every balanced state a single-flip local minimum and the couplings 8 bits wide. On 14-node cases the brute-force check finds balanced ground states at P = 4.
- **TSP:** x_{v,j} = (1 + s_{v,j})/2 (city v at position j, positions cyclic, spin index v·n + j, n² spins = 289–841). H = A Σ_v (1 − Σ_j x_vj)² + A Σ_j (1 − Σ_v x_vj)² + B Σ_j Σ_{u≠v} W_uv x_{u,j} x_{v,j+1}, **B = 1, A = max_uv W_uv** (Ayodele's MQC, the best TSP rule of ReAIM's reference [4]; Lucas's condition is A > B·max W, so MQC is its boundary).
  - Integer form (×4): J = −2A on same-city and same-position pairs; J = −W_uv on (u, j)–(v, j ± 1), u ≠ v; b_{v,j} = −4A(n − 2) − 2 Σ_u W_uv. E = 4H + const.
  - **Feasible** iff x is a permutation matrix. The tour length is computed from the TSPLIB matrix.
- **Verification** (`verify.py`, `verify_hope.json`, all passed on gpu-host):
  - by brute force over all 2^N states (GPP N ≤ 14; TSP n = 3, 4, i.e. 9 and 16 spins), the energy equals the Lucas objective plus the derived constant; the fields equal J s + b;
  - GPP ground states are balanced with the minimum bisection; for TSP the feasible minimum energy is 4·L* + const, and no infeasible state is lower (MQC may tie).
- **Coupling statistics** (`instance_stats.py`, `data/instance_stats.json`):
  - MCP: one |J| level (ternary).
  - GPP: dense, |J| ∈ {3, 4}, 4 signed bits.
  - TSP: 119–257 distinct |J| levels, max |J| = 2A = 560–1,730 (11–12 signed bits), max |b| = 36,546–88,060 (17–18 signed bits).

## Field (bias) handling: identical for all ten methods

Every method uses h_i = Σ_j J_ij s_j + b_i. This is the authors' option (b), which equals option (a): a frozen auxiliary spin s_0 = +1 with J_i0 = b_i that never updates and is counted nowhere (not in n_lin, not in any statistic).

- **Incremental methods** (SCA family, SA, ReAIM): b enters the fields at initialization and persists through the per-flip updates. This is the hardware scheme of `MULTIBIT_SPEC.md` (bias section).
- **SB methods:** b_i is added to every J·x product (the frozen x_0 = 1).
- **Dense balance penalty:** carried exactly as a uniform coupling γ on all pairs: field += γ(M − s_i). Integer arithmetic makes this identical to a dense matrix.
- **Verified:**
  - with b = 0 and γ = 0, every kernel is bit-identical to `engine.run` / `others.*` (G1, G14, synthetic graphs, and `run_gset.run_method` on G2);
  - with b ≠ 0 or γ ≠ 0, the SCA family is bit-identical to a dense copy of `abl.run` with h(0) = sJ + b; SA, ReAIM and dSB are bit-identical to dense references.

## Methods (the ten of the paper's comparison; `solvers.py`)

| Index | Method | Implementation |
|---|---|---|
| 0 | SA | sequential Metropolis sweeps, geometric T0 → T1 |
| 1 | SCA (plain) | engine family plain |
| 2 | TEC | field h + J_v s(t−1) |
| 3 | APC-SCA | per-spin q (reset/decay) |
| 4 | ReAIM ASA, noise-free | Algorithms 2–3, ITER_trial 32, ITER_run 96, T 1 → T1 |
| 5 | aSB | Goto 2019 |
| 6 | bSB | Goto 2021 |
| 7 | dSB | Goto 2021 |
| 8 | Onsager-kT | h − κ_eff·T(t−1)·ramp(t)·s(t−1) |
| 9 | Onsager-online | h − λ_eff·ramp(t)·c(t−1)·s(t−1), c = n_lin/(2T) |

- `ramp` is the hardware ramp (the coefficient goes to 0 over the last 30%).
- Step definitions are those of the K2000 comparison: an SA sweep, an SCA parallel update, an SB step, a ReAIM iteration.

## Scaling (pre-declared)

- **Temperature unit σ_T.** All temperatures, q, J_v, q_reset and q_lim are in units of σ_T:
  - MCP: sqrt(mean degree), as in the G-set study;
  - GPP: sqrt(mean graph degree), the RMS cut-part field; the penalty part vanishes at balance;
  - TSP: A = max W.
- **Onsager-online:** λ_eff = λ·lam_factor, with lam_factor = (mean_i Σ_j J_ij²/N)/(1999/2000), λ in K2000 units, and J the full matrix including penalties (the G-set rule and `MULTIBIT_SPEC.md`'s γ factor).
- **Onsager-kT:** κ_eff = κ·kappa_factor, with kappa_factor = mean_i Σ_j J_ij²/σ_T². For MCP this is exactly 1 (the G-set rule); with absolute K2000 temperatures it is `MULTIBIT_SPEC.md`'s κ scaling.
- **Limits, stated honestly.** The lag-one echo coefficient of spin i is Σ_j J_ij² f′(g_j); both forms use one global, uniform approximation of it.
  - With dense uniform penalties (GPP) or row/column penalties (TSP), Σ J² is dominated by the penalty couplings (benchmark instances: kappa_factor 260 (G1–G5), ≈ 1,085 (G14–G17), 134–236 (TSP); lam_factor 15.6–15.9 (GPP), 2.4·10⁴–2.9·10⁵ (TSP)).
  - The approximation is derived for zero-mean, disordered couplings. The collective (uniform-mode) feedback of a dense antiferromagnetic penalty is not an echo, and no global coefficient can represent it. The grids therefore span 2^−10 to 2^1 in these units.
- **SB:** Goto's normalisation with σ_J = sqrt(Σ J²/(N(N−1))) of the full J. b is not in σ_J.

## Exploratory calibration (pre-protocol; synthetic instances only, never a benchmark instance)

`calib.py`, `calib_analyze.py`, `calib/`, `logs/calib_*.log`.

- **Instances:** GPP on random G(800, 19176) and planar-like 800-node graphs (G-set study's `synth.py`, +1 weights, P = 4); TSP on random Euclidean instances with n = 17 and 29 (exact optima by integer programming).
- **Runs:** 16 runs per point at S = 1024 (GPP) and 2048 (TSP), with broad log-spaced sweeps (stage 1, 1b, 2b, 2c).
- **Findings** (best mean normalized quality, infeasible = 0):

| Method | GPP | TSP |
|---|---|---|
| SA | 0.98–1.00 | 0.51–0.60 |
| bSB | 0.99–1.00 (ξ ≈ 16, dt ≈ 0.06) | 0.96–1.00 |
| ReAIM | 0.86–0.94 (smallest k sets) | 0.52–0.65 |
| APC-SCA | 0.47–0.81 (needs q_reset ≥ 128σ) | 0.56–0.57 |
| TEC | ≤ 0.25 | 0.44–0.49 |
| Onsager-kT | ≤ 0.25 | 0.43–0.47 |
| Onsager-online | ≤ 0.25 | 0.41–0.46 |
| plain SCA | ≤ 0.15 (feasibility ≤ 31%) | 0.38–0.46 |
| dSB | ≤ 0.11 | 0.50–0.56 |
| aSB | ≤ 0.25 | 0 |

- **Mechanism (GPP).** The dense balance penalty makes the uniform mode of synchronous updates unstable (gain ≈ P·n_lin/(2T) ≫ 1): from a random start every majority spin flips at once.
- **TSP.** The SCA family needs q ≈ 16σ_T and T0 ≈ 16–32σ_T.
- **Grids.** They are placed around these regions, and every method has the same number of points.

## Grids (`grids.py`; 32 points for every method and problem; values in units of σ_T)

| Method | GPP | TSP |
|---|---|---|
| SA | T0 {0.25, 0.5, 1, 2} × T1 {0.01, 0.02, 0.04, 0.08, 0.16, 0.32, 0.48, 0.64} | T0 {0.5, 1, 2, 4} × T1 {0.04, 0.08, 0.16, 0.24, 0.32, 0.48, 0.64, 0.96} |
| SCA | q {4, 8, 16, 32} × T0 {0.5, 1, 2, 4} × T_fin {0.15, 0.6} | q {8, 12, 16, 24} × T0 {8, 16, 32, 64} × T_fin {0.15, 0.6} |
| TEC | J_v {−16, −8, −4, −1} × shared | J_v {−16, −8, −4, −2} × shared |
| Onsager-kT | κ {2^−10, 2^−8, 2^−6, 2^−4} × shared | same κ × shared |
| Onsager-online | λ {2^−8, 2^−5, 2^−2, 2} × shared | λ {2^−8, 2^−6, 2^−4, 2^−2} × shared |
| (shared, ⊂ SCA grid) | q {8, 32} × T0 {1, 4} × T_fin {0.15, 0.6} | q {12, 16} × T0 {16, 32} × T_fin {0.15, 0.6} |
| APC-SCA | q_reset {64, 128, 256, 512} × r_q {0.97, 0.99} × q_lim {2, 8} × T0 {1, 4}, T_fin 0.15 | q_reset {32, 64, 256, 512} × r_q {0.97, 0.99} × q_lim {0.125, 0.5} × T0 {1, 4}, T_fin 0.15 |
| ReAIM ASA | k set {(1,1,2,2), (1,2,3,4), (1,2,4,8), (2,4,8,16)} × T1 {0.4, 0.2, 0.1, 0.05, 0.025, 0.0125, 0.00625, 0.003125} | same |
| bSB | dt {1/32, 1/16, 1/8, 1/4} × ξ {6, 8, 11, 16, 22, 32, 45, 64} | dt {0.125, 0.25, 0.5, 0.75} × ξ {0.5, 0.71, 1, 1.41, 2, 2.83, 4, 8} |
| dSB | dt {1/32, 1/16, 1/8, 1/4} × ξ {0.5, 1, 2, 4, 8, 16, 32, 64} | dt {1/16, 1/8, 1/4, 1/2} × ξ {1, 1.41, 2, 2.83, 4, 5.66, 8, 16} |
| aSB | dt {1/32, 1/16, 1/8, 1/4} × ξ {1, 2, 4, 8, 12, 16, 24, 32} | dt {1/32, 1/16, 1/8, 1/4} × ξ {0.5, 1, 2, 4, 8, 16, 32, 64} |

## Procedure (identical for all methods; `run_bench.py`)

- **Budgets.** GPP S ∈ {256, 512, 1024, 2048, 4096}; TSP S ∈ {512, 1024, 2048, 4096, 8192}. The maxima are ReAIM's iteration counts.
- **Pilot.** For each instance, method, S and grid point: 64 runs, seed SeedSequence([20261008, 1, p, i, method index, S index, grid index]).
  - p = 1 (GPP) or 2 (TSP);
  - i = the G-set number (GPP) or 100 + the TSPLIB index (gr17 = 100 … bays29 = 105).
- **Selection (fixed).** Highest pilot mean normalized quality (infeasible runs count 0), ties to the lower grid index. A configuration is invalid if any run is non-finite (aSB).
- **Held-out final.** The selected configuration, 256 runs, seed SeedSequence([20261008, 2, p, i, method index, S index]).
- **Final state** only, as on hardware. No repair or post-processing of any run.
- **Starts.** All methods start from random spins (SB from small random amplitudes), as in the G-set study.
- **Max-Cut G1–G20.**
  - The four engine rules are taken unchanged from `research/gset_20261007` (pre-registered grids, 256-run finals at S ∈ {250, …, 4000}).
  - The six other methods (SA, APC-SCA, ReAIM, aSB, bSB, dSB) already exist there for G1, G6, G11, G14 and G18 (Part B).
  - **Extension, amendment-style.** The same six methods on the other 15 instances (G2–G5, G7–G10, G12, G13, G15–G17, G19, G20), using `run_gset.py`'s own `grid()`, `run_method()` and `stats()`, its 16-point grids, its budgets and its seed formulas (seed 20261007, pilot [.., 1, g, mi, si, gi], final [.., 2, g, mi, si]). Results are written here (`results/mcp_ext/`).

## Measures

- **Normalized quality per run** (ReAIM-style, higher is better; 0 if infeasible):
  - MCP: cut/BKV;
  - GPP: R/cut for a balanced partition;
  - TSP: L*/L for a valid tour.
- **GPP references R** (`gpp_reference.py`, `data/gpp_reference.json`, computed before any held-out run; exploratory).
  - Method: the best balanced cut from 256 SA runs of 65,536 sweeps each, then exact-balance repair and a best-improvement balanced swap descent.
  - Values: G1 7590, G2 7581, G3 7579, G4 7588, G5 7583, G14 1089, G15 1091, G16 1071, G17 1058.
  - Each is at or below ReAIM's Table V minimum (7597, 7611, 7588, 7620, 7612, 1115, 1121, 1085, 1075). It is an upper bound on the true minimum bisection.
- **Targets** (all within 1% of the optimum or best-known):
  - MCP: cut ≥ ⌈0.99·BKV⌉;
  - GPP: balanced and cut ≤ ⌊1.01·R⌋;
  - TSP: valid and L ≤ ⌊1.01·L*⌋.
  - Also reported: the optimum or best-known itself (MCP cut ≥ BKV, GPP cut ≤ R, TSP L = L*).
- **Per (instance, method, S), from the 256-run finals:**
  - feasibility rate (Wilson 95%);
  - mean normalized quality, overall and over feasible runs;
  - P(target) (Wilson 95%) and P(optimum);
  - best value;
  - TSP ARPD over feasible runs;
  - ReAIM's best-of-20: the median over the 12 disjoint blocks of 20 runs (runs 0–239), with the number of blocks that have no feasible run;
  - steps to target: MCS99(S) = S·ln 0.01/ln(1 − p) (S if p = 1, ∞ if p = 0), minimum over S.
- **Paper summary.** ReAIM-matched budgets are MCP S = 4000 (the closest G-set budget to 4096), GPP 4096 and TSP 8192. Steps-to-target uses the minimum over S.

## Hypotheses (stated before the held-out runs; calibration-based, honest)

| ID | Expectation | Pass condition |
|---|---|---|
| R1 | On GPP the SCA engine rules (plain, TEC, kT, online) are far from SA and bSB | Each engine rule's mean quality at S = 4096 is below SA's on at least 8 of 9 instances |
| R2 | The Onsager forms do not rescue SCA on dense-penalty GPP | Neither form's feasibility exceeds 50% at S = 4096 on more than 4 of 9 instances |
| R3 | On TSP bSB is the best of the ten | Highest mean quality at S = 8192 on at least 4 of 6 instances |
| R4 | On TSP the Onsager forms are on par with plain SCA and TEC | A difference in mean quality at S = 8192 of at most 0.05 on at least 4 of 6 instances |
| R5 | MCP G1–G20 ordering as on the 12-instance subset | SA and dSB need fewer steps (MCS99) than every engine rule on at least 15 of 20 instances |

These are descriptive expectations, and results are reported as measured whatever the outcome. The key comparison for the paper is ReAIM-style quality and feasibility per method and problem.

## Coupling-precision study (four engine rules; specified now, run after the held-out runs)

- **Instances:** MCP G1–G20 (the reference case), GPP (9) and TSP (6).
- **Rules:** plain SCA, TEC, Onsager-kT and Onsager-online.
- **Schedule per (instance, rule).**
  - The held-out-selected configuration at S_exp = argmin over S of MCS99, ties to the smaller S.
  - If MCS99 is ∞ at every S: the S with the highest final mean quality, ties to the larger S.
  - MCP schedules come from the G-set study's pre-registered finals.
- **Precision K ∈ {2, 3, 4, 6, 8} and full** (the unquantized integer J).
  - α_K = (2^(K−1) − 1)/max|J|.
  - J^K = round(α_K·J) entrywise, rounding half away from zero, on the full matrix (sparse part plus the uniform penalty).
  - b^K = round(α_K·b), stored as int32. The bias is not limited to K bits; its bit width is reported (the hardware's B_BITS = 16).
  - The schedule maps as follows: T0, T_fin, q and J_v are multiplied by α_K; λ_eff and κ_eff are recomputed from J^K (lam_factor and kappa_factor of J^K, with σ_T^K = α_K·σ_T); the same S and ramp.
  - The objective is always evaluated on the true, unquantized problem from the final spins.
- **Arithmetic:** `hwmodel.py`, a port of `sca_ref_bias.hpp::run_trial_bias` with `sca_ref.hpp`'s Rng and `mb_common.hpp`'s tables, LANES = 256.
  - For K ≤ 8 it is checked bit for bit against the C++ reference (tables, spins, flips, n_lin trace, Σ s·h, Σ b·s).
  - Status: `verify_hw.py`, 60 of 60 cases identical on gpu-host, built with c++ −O2 −std=c++17 −DSCA_LANES=256.
  - Full precision uses the same integer arithmetic with 64-bit tables.
  - Trials 0–1023, seed 20261008 + 1000·p + 10·(instance code mod 100) + rule index (p = 0 MCP, 1 GPP, 2 TSP; rule index 0 plain, 1 TEC, 2 kT, 3 online), identical across K, so that runs are paired.
- **Measures per K:**
  - mean normalized quality (infeasible = 0), feasibility rate and P(target);
  - the paired difference from full precision, with a 95% interval.
- **Tolerance (pre-declared).** K is within tolerance if the mean quality difference (K − full) is ≥ −0.005.
  - K_min = the smallest tested K such that it and every larger tested K are within tolerance, per instance and rule.
  - Summary: the largest K_min over instances, and the fraction of instances within tolerance at each K.
- **Hardware range checks** are reported for each exported schedule (`hwmodel.table_checks`): int32 tables; the FPGA engine's 27-bit 4T lane constant (T < 512) and the |q| + |corr| < 2^29 bound (K ≤ 4).

## Hardware export (`hw_export/`)

- **Files** for each GPP and TSP instance at K ∈ {2, 4, 8}:
  - `<inst>_K<K>.jint8`: SCAJINT8, an 8-byte magic, uint32 N, uint32 0, then N×N int8 row-major (`gpu/multibit_20261007/src/jmat.hpp`);
  - `<inst>_K<K>.bias.bin`: int32 little-endian.
- **`selected_configs.json`**, per rule:
  - t0, t1, S, q, J_v, κ_eff, λ_eff and ramp in that K's units;
  - the range checks;
  - the measured precision-study quality at that K;
  - the objective mapping (spin layout, decoding, feasibility, true objective, E = −(Σ s·h + Σ b·s)/2);
  - SHA-256 of every file.

## Platform

- **Remote run:** gpu-host (gpu-host, aarch64, 144 Grace cores), Python 3.12.3, numpy 1.26.4, numba 0.60.0, scipy 1.13.1, in a venv under `/scratch/USER/anechoic_cpu_20261008`.
- **Laptop:** arm64, Python 3.12.2, numpy 1.26.4, numba 0.60.0.
- **Reproducibility:** laptop and gpu-host are bit-identical on a 33-case fingerprint of all ten methods and the hardware model (`repro.py`, `repro_laptop.json`, `repro_hope.json`).
- **Workers:** up to 64 single-threaded processes on gpu-host.

## What this does not show

- No hardware runs, and no timing. Steps are not time across methods (an SA sweep is N sequential updates; an SCA step is one parallel update).
- The penalty weights are fixed by the rules above, not tuned per method. Other penalties change all results.
- The GPP references are our own heuristic values, not proven optima.
- The calibration used two synthetic instances per problem. Grid edges will be reported.
