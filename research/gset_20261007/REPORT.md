## G-set transfer test (software model of the engine): final report

`selected_configs.json` (for the FPGA and GPU ternary engines) is written and holds both the frozen-grid and the extended selections.

### Bottom line

The answer splits by form of the correction.

- **Onsager-kT (κ·T coefficient, the form of the hardware's fastest schedule, X5) transfers to G-set.**
  - It reaches the target faster than plain SCA on **48 of 51** instances, geometric mean **1.91×** (12-engine estimator, pre-registered grids). With the pre-registered grid extension: 50 of 51, **2.13×**.
  - It beats TEC on 41 of 51 instances (1.51×); extended, 40 of 51 (1.62×).
  - It needs about 2× fewer steps than plain SCA, and never more steps on any instance.
- **Onsager-online (n_lin/2T popcount) helps much less.**
  - About 1.2× over plain SCA, and roughly on par with TEC (0.93–1.07×).
  - It loses to TEC on 11–12 of the 12 planar-like +1 instances. Those graphs have broad degree distributions (degree CV ≈ 1.0–1.2), where one global coefficient cannot match the per-spin echo Σ_j J_ij² f′.
- **Against the other methods (12-instance subset),** sequential SA and dSB need far fewer steps than any SCA rule on these sparse graphs. Any engine advantage there would come from per-step cost, not from fewer steps.

All of these are software-model results; the ±1-only hardware cannot run G-set.

### What was run

- **Data.**
  - G1–G47 and G51–G54 (51 instances, N ≤ 2048) from https://web.stanford.edu/~yyye/yyye/Gset/.
  - All checks pass: edge count equals the header, weights ±1, no loops or duplicates.
  - SHA-256 is in `manifest.json`. Two independent downloads are byte-identical; the second was made after the 16:12 VPN drop, during which G25 had first timed out and been re-fetched.
  - Every ±1 class has exactly the edge set of its +1 sibling, including G23, whose server file is dated 2023.
- **Best-known cuts.** Parsed programmatically from two tables (`bkv.json`, `build_bkv.py`):
  - Benlic & Hao, EAAI 26(3):1162–1173, 2013, doi:10.1016/j.engappai.2012.09.001, Table 2.
  - Ma & Hao, Ann. Oper. Res. 248:365–403, 2017, doi:10.1007/s10479-016-2234-0, Table 5.
  - BKV = max(BLS f_best, Ma–Hao f_pre, Ma–Hao MOH). Target = ⌈0.99·BKV⌉, computed exactly.
  - Both tables were read from preprints (the authors' accepted manuscript; arXiv:1510.09156v1). The published versions were not compared.
  - **Flagged:** G23. BLS's "previous best" of 13,354 appears nowhere else, so BKV = 13,344 (target 13,211). At 13,221 conclusions are unchanged; online confirmation success goes from 976 to 906 of 1,024.
  - **Values improved after BLS:** G30, G31 and G35–G38.
  - G22 is 13,359 (target 13,226) and G32 is 1,410 (target 1,396), matching the GPU study's values.
- **Model.** `engine.py` is bit-identical to `abl.run` (45 sparse cases: identical spins and flips). The only changes are a T_fin parameter and an exact sparse field update. `others.py` holds SA, ReAIM and dSB (bit-identical to `methods.py`) and bSB and aSB (statistically equivalent).
- **Cost model.** v6.4, 0.1424·flips + 18.09·S + 886 cycles at 250 MHz.
  - **Stated assumption:** dense rows and 2,048 slots per step mean the same model applies to sparse graphs; a ternary engine needs two popcounts per update but has the same cycle structure.
  - Single-engine TTS99 is from Eq. 9.
  - The 12-engine round estimator (model-based) uses t_round = E[max of 12 trial times], P_round = 1 − (1 − p)^12, and counts one round if P_round ≥ 0.995.
- **Scaling.**
  - T and q are in units of σ = sqrt(mean degree), the RMS local field.
  - λ_eff = λ·(d̄/N)/(1999/2000).
  - κ carries no d/N factor: with T ∝ σ, κ·T already scales with the field.
  - Tuning was identical for all methods: S ∈ {250, 500, 1000, 2000, 4000}, 16-point grids, 64-run pilots, selection by pilot mean cut, 256-run held-out finals. For the engine rules, S* is chosen per estimator and confirmed with 1,024 fresh paired runs.
- **Exploratory calibration (pre-protocol, synthetic graphs only).** On +1-weight graphs, synchronous SCA with K2000-scaled q oscillates between the all-up and all-down states and ends at cut 0. The shared q range was therefore set per weight class (q ≈ 1–3σ on +1 classes).
- **Protocol history.**
  - PROTOCOL.md was frozen at 14:36:48Z, before any G-set run.
  - The subset pilots showed the three-parameter rules' 2×2 shared (q, T0) subset was edge-limited in every class. Amendment 1 (a grid extension: 16 more points for every rule, union selection, fresh seeds) was frozen at 14:59:01Z, before its own runs but after the subset's original finals had been logged. Both analyses are reported.
- **Compute.** 6 single-threaded workers. Part A 83 min, C 23 min, A2 45 min, C2 20 min, B 33 min.

### Results: four engine rules on all 51 instances

TTS ratio = TTS(reference)/TTS(rule), so > 1 means the rule is faster. Geometric mean over 51 instances, with the number of instances where the rule is faster. Every rule reached the target on every instance.

| Ratio | Pre-registered, 1 engine | Pre-registered, 12 engines | Extended, 1 engine | Extended, 12 engines |
|---|---|---|---|---|
| plain / Onsager-kT | 1.89× (47) | **1.91× (48)** | 2.21× (49) | **2.13× (50)** |
| plain / Onsager-online | 1.37× (44) | 1.20× (37) | 1.35× (39) | 1.22× (33) |
| TEC / Onsager-kT | 1.48× (41) | 1.51× (41) | 1.62× (45) | 1.62× (40) |
| TEC / Onsager-online | 1.07× (29) | 0.95× (25) | 0.99× (22) | 0.93× (24) |
| plain / TEC | 1.28× (38) | 1.26× (40) | 1.36× (40) | 1.31× (37) |
| plain / best Onsager form | 2.02× | 1.99× | 2.33× | 2.23× |

- **Spread (12 engines).** plain/kT ranges 0.85–14.2× (pre-registered) and 0.87–18.4× (extended).
- **Paired-bootstrap 95% intervals.**
  - plain/kT is entirely above 1 on 46/51 (pre-registered) and 45–48/51 (extended), and below 1 on at most one instance.
  - TEC/online is split, about 16–20 instances each way.

**Hypotheses.**

| Grids | Single engine | 12 engines |
|---|---|---|
| Pre-registered | G1–G4 pass; G5 passes at 2.02 | G1′, G2′, G4′ pass; **G3′ fails** (online beats TEC on 25/51, 26 needed); **G5′ fails** narrowly (1.99 < 2) |
| Extended | G1, G2, G4, G5 (2.33) pass; **G3 fails** (22/51) | G1′, G2′, G4′, G5′ (2.23) pass; **G3′ fails** (24/51) |

**Steps to solution** (MCS99 = S·ln 0.01/ln(1 − p), minimum over S; hardware-independent):

| Ratio | Pre-registered | Extended |
|---|---|---|
| plain / Onsager-kT | 2.06× (kT fewer steps on 50, more on 0) | 2.21× (50/0) |
| plain / Onsager-online | 1.61× (45/2) | 1.51× (41/6) |
| TEC / Onsager-kT | 1.60× (43/4) | 1.67× (44/3) |
| TEC / Onsager-online | 1.25× (25/25) | 1.15× (24/27) |

**Best-known cut reached** (5×256 finals plus 2×1024 confirmations, pre-registered grids): plain 19/51 instances, TEC 21, Onsager-kT 21, Onsager-online 19. No rule ever reached it on G14–G18, G22–G26, G30–G42 or G51–G54.

### Ten methods on the 12-instance subset

Subset: G1, G6, G11, G14, G18, G22, G27, G32, G35, G39, G43, G51. Pre-registered grids, all methods at 16 points. MCS99 is the geometric mean relative to Onsager-online over solved instances (< 1 means fewer steps).

| Method | Solved | MCS99 vs online | BKV hits (runs, instances) |
|---|---|---|---|
| SA (sweeps) | 12/12 | 0.16× | 840, 6 |
| dSB | 12/12 | 0.34× | 295, 6 |
| bSB | 8/12 | 0.53× | 0 |
| APC-SCA | 12/12 | 0.76× | 146, 5 |
| Onsager-kT | 12/12 | 0.83× | 158, 5 |
| Onsager-online | 12/12 | 1.00× | 108, 4 |
| aSB | 6/12 | 1.41× | 0 |
| TEC | 12/12 | 1.48× | 101, 5 |
| plain SCA | 12/12 | 1.67× | 106, 5 |
| ReAIM ASA | 9/12 | 9.2× | 2, 2 |

- **SA and dSB** are strongest on G-set, as in the literature. On G32, SA needs 6,719 sweeps against 195,675 steps for Onsager-kT.
- **APC-SCA** edges out both Onsager forms on steps. Its per-spin pinning directly damps the +1 oscillation, but it needs per-spin q registers our engine lacks, so it has no engine TTS.
- Steps are not time: an SA sweep is N sequential updates, an SCA step one parallel update.

### Per-instance table (12-engine model TTS99 in µs; S* and p from the 1,024-run confirmations; extended grids, with pre-registered ratios at right)

| Inst | Class | BKV | Target | plain | TEC | kT | online | plain/kT | plain/online | TEC/kT | TEC/online | pre-reg plain/kT | plain/online | TEC/kT | TEC/online |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| G1 | R800+ | 11624 | 11508 | 51.7 (S500, p0.68) | 23.0 (S250, p0.47) | 26.2 (S250, p0.81) | 64.5 (S500, p0.68) | 1.97 | 0.80 | 0.88 | 0.36 | 1.86 | 0.69 | 2.60 | 0.97 |
| G2 | R800+ | 11620 | 11504 | 45.8 (S500, p0.50) | 23.0 (S250, p0.61) | 25.0 (S250, p0.91) | 60.4 (S250, p0.14) | 1.83 | 0.76 | 0.92 | 0.38 | 2.97 | 1.22 | 1.76 | 0.72 |
| G3 | R800+ | 11622 | 11506 | 46.7 (S500, p0.75) | 23.0 (S250, p0.51) | 28.5 (S250, p0.79) | 60.8 (S250, p0.16) | 1.64 | 0.77 | 0.81 | 0.38 | 2.06 | 0.65 | 1.63 | 0.51 |
| G4 | R800+ | 11646 | 11530 | 49.0 (S500, p0.71) | 23.0 (S250, p0.43) | 27.5 (S250, p0.85) | 54.1 (S500, p0.70) | 1.78 | 0.91 | 0.84 | 0.43 | 2.06 | 1.04 | 1.96 | 0.99 |
| G5 | R800+ | 11631 | 11515 | 45.8 (S500, p0.54) | 23.0 (S250, p0.53) | 28.7 (S250, p0.90) | 53.2 (S250, p0.19) | 1.60 | 0.86 | 0.80 | 0.43 | 2.92 | 1.22 | 1.45 | 0.61 |
| G6 | R800± | 2178 | 2157 | 67.6 (S500, p0.39) | 69.5 (S500, p0.43) | 35.2 (S250, p0.37) | 39.8 (S250, p0.66) | 1.92 | 1.70 | 1.97 | 1.75 | 2.09 | 1.70 | 2.15 | 1.75 |
| G7 | R800± | 2006 | 1986 | 72.2 (S500, p0.31) | 91.5 (S500, p0.27) | 38.0 (S250, p0.30) | 41.1 (S250, p0.33) | 1.90 | 1.76 | 2.41 | 2.22 | 2.11 | 1.86 | 2.33 | 2.05 |
| G8 | R800± | 2005 | 1985 | 70.0 (S500, p0.36) | 74.5 (S250, p0.16) | 35.3 (S250, p0.44) | 40.7 (S250, p0.41) | 1.98 | 1.72 | 2.11 | 1.83 | 2.15 | 1.86 | 2.02 | 1.74 |
| G9 | R800± | 2054 | 2034 | 83.3 (S500, p0.28) | 75.9 (S500, p0.31) | 34.4 (S250, p0.33) | 42.0 (S250, p0.40) | 2.42 | 1.98 | 2.21 | 1.81 | 2.63 | 2.00 | 2.73 | 2.07 |
| G10 | R800± | 2000 | 1980 | 92.1 (S500, p0.22) | 62.9 (S500, p0.32) | 38.3 (S250, p0.30) | 43.2 (S250, p0.41) | 2.41 | 2.14 | 1.64 | 1.46 | 2.59 | 2.06 | 1.89 | 1.51 |
| G11 | T800± | 564 | 559 | 1720 (S1000, p0.03) | 1274 (S4000, p0.14) | 520 (S1000, p0.08) | 1631 (S1000, p0.03) | 3.31 | 1.05 | 2.45 | 0.78 | 2.64 | 0.85 | 2.56 | 0.83 |
| G12 | T800± | 556 | 551 | 1326 (S4000, p0.12) | 1116 (S4000, p0.14) | 480 (S4000, p0.28) | 1323 (S4000, p0.17) | 2.76 | 1.00 | 2.33 | 0.84 | 2.47 | 1.50 | 2.02 | 1.22 |
| G13 | T800± | 582 | 577 | 3891 (S4000, p0.05) | 3156 (S4000, p0.05) | 1342 (S2000, p0.06) | 2482 (S4000, p0.08) | 2.90 | 1.57 | 2.35 | 1.27 | 2.00 | 1.06 | 1.59 | 0.84 |
| G14 | P800+ | 3064 | 3034 | 172 (S2000, p0.37) | 111 (S1000, p0.44) | 82.8 (S1000, p0.62) | 153 (S1000, p0.32) | 2.08 | 1.12 | 1.34 | 0.72 | 1.44 | 1.13 | 0.90 | 0.71 |
| G15 | P800+ | 3050 | 3020 | 168 (S2000, p0.32) | 188 (S1000, p0.21) | 72.4 (S1000, p0.36) | 202 (S2000, p0.29) | 2.33 | 0.84 | 2.60 | 0.93 | 0.96 | 0.46 | 0.98 | 0.46 |
| G16 | P800+ | 3052 | 3022 | 225 (S2000, p0.25) | 138 (S1000, p0.27) | 82.8 (S1000, p0.42) | 268 (S1000, p0.21) | 2.71 | 0.84 | 1.67 | 0.52 | 1.73 | 0.66 | 1.12 | 0.43 |
| G17 | P800+ | 3047 | 3017 | 165 (S2000, p0.33) | 106 (S1000, p0.35) | 82.9 (S1000, p0.51) | 209 (S1000, p0.26) | 1.99 | 0.79 | 1.28 | 0.51 | 1.30 | 0.78 | 0.97 | 0.59 |
| G18 | P800± | 992 | 983 | 226 (S1000, p0.22) | 245 (S500, p0.12) | 136 (S500, p0.21) | 134 (S500, p0.20) | 1.66 | 1.69 | 1.80 | 1.83 | 1.35 | 1.49 | 1.50 | 1.66 |
| G19 | P800± | 906 | 897 | 190 (S1000, p0.24) | 203 (S1000, p0.25) | 150 (S500, p0.19) | 150 (S1000, p0.43) | 1.27 | 1.27 | 1.35 | 1.35 | 1.58 | 1.96 | 1.07 | 1.33 |
| G20 | P800± | 941 | 932 | 135 (S1000, p0.38) | 157 (S250, p0.10) | 88.7 (S250, p0.17) | 149 (S250, p0.13) | 1.52 | 0.91 | 1.77 | 1.06 | 1.59 | 2.02 | 1.54 | 1.96 |
| G21 | P800± | 931 | 922 | 424 (S2000, p0.21) | 473 (S2000, p0.24) | 323 (S1000, p0.18) | 208 (S1000, p0.24) | 1.31 | 2.04 | 1.46 | 2.27 | 1.62 | 1.85 | 1.31 | 1.50 |
| G22 | R2000+ | 13359 | 13226 | 79.8 (S500, p0.67) | 36.1 (S250, p0.30) | 34.0 (S250, p0.33) | 62.4 (S500, p0.86) | 2.35 | 1.28 | 1.06 | 0.58 | 2.43 | 1.28 | 1.13 | 0.60 |
| G23 | R2000+ | 13344 | 13211 | 78.5 (S250, p0.19) | 39.0 (S250, p0.67) | 36.1 (S250, p0.63) | 33.9 (S250, p0.40) | 2.17 | 2.31 | 1.08 | 1.15 | 2.17 | 2.31 | 1.08 | 1.15 |
| G24 | R2000+ | 13337 | 13204 | 54.4 (S500, p0.88) | 33.6 (S250, p0.58) | 45.3 (S250, p0.56) | 30.0 (S250, p0.35) | 1.20 | 1.81 | 0.74 | 1.12 | 1.20 | 1.60 | 0.74 | 0.99 |
| G25 | R2000+ | 13340 | 13207 | 72.2 (S250, p0.15) | 39.1 (S250, p0.66) | 36.1 (S250, p0.59) | 34.4 (S250, p0.30) | 2.00 | 2.10 | 1.08 | 1.14 | 1.72 | 1.60 | 1.08 | 1.01 |
| G26 | R2000+ | 13328 | 13195 | 54.4 (S500, p0.92) | 33.7 (S250, p0.62) | 36.2 (S250, p0.65) | 34.0 (S250, p0.42) | 1.50 | 1.60 | 0.93 | 0.99 | 1.50 | 1.60 | 0.93 | 0.99 |
| G27 | R2000± | 3341 | 3308 | 247 (S1000, p0.60) | 229 (S1000, p0.65) | 102 (S500, p0.48) | 142 (S500, p0.50) | 2.42 | 1.74 | 2.24 | 1.62 | 2.42 | 1.74 | 1.52 | 1.09 |
| G28 | R2000± | 3298 | 3266 | 249 (S1000, p0.71) | 106 (S500, p0.35) | 108 (S250, p0.18) | 144 (S500, p0.54) | 2.32 | 1.73 | 0.99 | 0.74 | 2.53 | 1.73 | 1.19 | 0.81 |
| G29 | R2000± | 3405 | 3371 | 248 (S1000, p0.37) | 229 (S1000, p0.44) | 96.6 (S500, p0.33) | 213 (S500, p0.23) | 2.57 | 1.16 | 2.37 | 1.07 | 2.51 | 1.19 | 2.31 | 1.10 |
| G30 | R2000± | 3413 | 3379 | 205 (S500, p0.19) | 124 (S500, p0.31) | 111 (S250, p0.18) | 140 (S500, p0.50) | 1.85 | 1.46 | 1.12 | 0.88 | 2.12 | 1.50 | 1.23 | 0.87 |
| G31 | R2000± | 3310 | 3277 | 213 (S500, p0.18) | 122 (S500, p0.32) | 102 (S500, p0.60) | 145 (S500, p0.50) | 2.08 | 1.47 | 1.19 | 0.84 | 1.84 | 1.30 | 1.12 | 0.79 |
| G32 | T2000± | 1410 | 1396 | 22328 (S4000, p0.01) | 26965 (S4000, p0.01) | 5473 (S2000, p0.02) | 29889 (S2000, p0.01) | 4.08 | 0.75 | 4.93 | 0.90 | 10.71 | 1.09 | 57.12 | 5.80 |
| G33 | T2000± | 1382 | 1369 | 90855 (S4000, p0.00) | 22991 (S4000, p0.01) | 4931 (S4000, p0.04) | 33813 (S4000, p0.01) | 18.43 | 2.69 | 4.66 | 0.68 | 14.23 | 3.25 | 4.48 | 1.02 |
| G34 | T2000± | 1384 | 1371 | 10490 (S4000, p0.02) | 23673 (S4000, p0.01) | 3190 (S4000, p0.06) | 11427 (S4000, p0.03) | 3.29 | 0.92 | 7.42 | 2.07 | 4.86 | 1.77 | 3.69 | 1.35 |
| G35 | P2000+ | 7687 | 7611 | 411 (S4000, p0.75) | 297 (S2000, p0.71) | 93.2 (S1000, p0.36) | 436 (S4000, p0.68) | 4.41 | 0.94 | 3.19 | 0.68 | 2.17 | 0.72 | 1.38 | 0.46 |
| G36 | P2000+ | 7680 | 7604 | 410 (S2000, p0.18) | 243 (S2000, p0.69) | 81.4 (S1000, p0.36) | 473 (S2000, p0.60) | 5.03 | 0.87 | 2.98 | 0.51 | 2.26 | 0.86 | 1.30 | 0.49 |
| G37 | P2000+ | 7691 | 7615 | 411 (S4000, p0.70) | 273 (S2000, p0.61) | 98.9 (S1000, p0.30) | 1189 (S1000, p0.08) | 4.16 | 0.35 | 2.76 | 0.23 | 1.16 | 0.42 | 0.81 | 0.29 |
| G38 | P2000+ | 7688 | 7612 | 411 (S4000, p0.75) | 236 (S500, p0.22) | 87.6 (S1000, p0.33) | 423 (S4000, p0.76) | 4.69 | 0.97 | 2.70 | 0.56 | 1.96 | 0.71 | 1.01 | 0.37 |
| G39 | P2000± | 2408 | 2384 | 2086 (S4000, p0.15) | 2246 (S4000, p0.16) | 1551 (S4000, p0.21) | 1415 (S1000, p0.07) | 1.35 | 1.47 | 1.45 | 1.59 | 1.44 | 1.46 | 1.26 | 1.29 |
| G40 | P2000± | 2400 | 2376 | 2418 (S4000, p0.11) | 2333 (S2000, p0.07) | 2789 (S2000, p0.06) | 2087 (S2000, p0.09) | 0.87 | 1.16 | 0.84 | 1.12 | 1.13 | 1.35 | 0.79 | 0.94 |
| G41 | P2000± | 2405 | 2381 | 2642 (S4000, p0.10) | 2280 (S2000, p0.07) | 1925 (S2000, p0.07) | 2526 (S2000, p0.08) | 1.37 | 1.05 | 1.18 | 0.90 | 2.19 | 1.22 | 1.90 | 1.06 |
| G42 | P2000± | 2481 | 2457 | 2331 (S4000, p0.11) | 3562 (S4000, p0.11) | 2307 (S4000, p0.17) | 1784 (S2000, p0.10) | 1.01 | 1.31 | 1.54 | 2.00 | 1.06 | 1.32 | 1.13 | 1.40 |
| G43 | R1000+ | 6660 | 6594 | 62.8 (S250, p0.18) | 27.7 (S250, p0.54) | 33.7 (S250, p0.46) | 28.0 (S250, p0.38) | 1.87 | 2.25 | 0.82 | 0.99 | 1.60 | 2.17 | 0.82 | 1.12 |
| G44 | R1000+ | 6650 | 6584 | 41.0 (S250, p0.26) | 30.6 (S250, p0.66) | 29.1 (S250, p0.58) | 28.0 (S250, p0.45) | 1.41 | 1.46 | 1.05 | 1.09 | 1.42 | 1.48 | 1.05 | 1.09 |
| G45 | R1000+ | 6654 | 6588 | 42.1 (S250, p0.25) | 30.6 (S250, p0.61) | 29.1 (S250, p0.59) | 28.0 (S250, p0.40) | 1.45 | 1.50 | 1.05 | 1.09 | 1.52 | 1.58 | 1.05 | 1.09 |
| G46 | R1000+ | 6649 | 6583 | 39.2 (S250, p0.27) | 30.5 (S250, p0.70) | 29.0 (S250, p0.68) | 28.0 (S250, p0.47) | 1.35 | 1.40 | 1.05 | 1.09 | 1.20 | 1.25 | 1.05 | 1.09 |
| G47 | R1000+ | 6657 | 6591 | 34.1 (S250, p0.25) | 30.5 (S250, p0.63) | 33.6 (S250, p0.63) | 28.0 (S250, p0.41) | 1.01 | 1.22 | 0.91 | 1.09 | 0.85 | 1.03 | 0.91 | 1.09 |
| G51 | P1000+ | 3848 | 3810 | 177 (S2000, p0.32) | 204 (S2000, p0.53) | 84.7 (S1000, p0.48) | 216 (S2000, p0.28) | 2.09 | 0.82 | 2.41 | 0.95 | 0.92 | 0.64 | 1.41 | 0.99 |
| G52 | P1000+ | 3851 | 3813 | 178 (S2000, p0.40) | 305 (S500, p0.08) | 84.7 (S1000, p0.57) | 163 (S2000, p0.35) | 2.10 | 1.09 | 3.60 | 1.87 | 1.25 | 0.22 | 2.22 | 0.39 |
| G53 | P1000+ | 3850 | 3812 | 178 (S2000, p0.42) | 141 (S1000, p0.42) | 84.6 (S1000, p0.60) | 196 (S1000, p0.32) | 2.11 | 0.91 | 1.66 | 0.72 | 1.15 | 0.82 | 0.95 | 0.68 |
| G54 | P1000+ | 3852 | 3814 | 273 (S2000, p0.22) | 215 (S1000, p0.24) | 78.2 (S1000, p0.34) | 365 (S4000, p0.58) | 3.49 | 0.75 | 2.75 | 0.59 | 1.49 | 0.60 | 1.36 | 0.54 |

More detail:
- Single-engine tables, per-instance bootstrap intervals and per-class win counts: `analysis.txt` (pre-registered), `analysis_ext.txt` (extended), and `results_summary*.json`.
- Full selected configurations, for both estimators and both variants: `selected_configs.json`. It gives q, T0, T_fin, S, ramp, κ, λ_eff (with the d/N factor already applied), J_v and APC parameters.
- Markdown versions: `python3 report_tables.py results_summary[_ext].json`.

### Recommended paper paragraph

> **Beyond K2000 (software model).** All hardware measurements use one dense instance. To test whether the correction's advantage depends on it, we ran the engine's bit-exact software model on all 51 G-set Max-Cut instances that fit the engine (N ≤ 2048: random, toroidal and planar-like graphs with +1 or ±1 weights). The target was 99% of the best-known cut [Benlic & Hao 2013; Ma & Hao 2017]. Every rule received the same pre-registered per-budget tuning: 16-point grids, 64-run pilots, 256-run held-out finals and 1,024-run paired confirmations. Temperatures and self-coupling were scaled to each instance's RMS local field, and the online coefficient by d/N. Under the v6.4 cycle model with 12 engines, the temperature-proportional correction (κ·T, the form of our fastest hardware schedule) reached the target faster than plain SCA on 48 of 51 instances (geometric mean 1.9×; 2.1× with a pre-registered grid extension). It was faster than TEC on 41 of 51 (1.5×), with the largest gains on toroidal graphs (about 4–5×), and needed about 2× fewer steps than plain SCA. The online-popcount form gained less (1.2× over plain SCA) and was on par with TEC overall. It was slower than TEC on planar-like +1 graphs, whose broad degree distributions a single global coefficient cannot represent. On these sparse graphs, sequential SA and discrete simulated bifurcation needed fewer steps than any SCA rule, so the engine's advantage there would come from its per-step cost, not from fewer steps. The current hardware stores only ±1 couplings and cannot yet run these instances; all G-set numbers are model results.

**Compact table.** Model TTS99 ratio with 12-engine rounds, rule faster on (k/n) instances; pre-registered grids unless stated.

| G-set group (instances) | Onsager-kT vs plain SCA | Onsager-kT vs TEC | Onsager-online vs plain SCA | Onsager-online vs TEC |
|---|---|---|---|---|
| Random, +1 (15) | 1.74× (14/15) | 1.21× (11/15) | 1.31× (13/15) | 0.91× (7/15) |
| Random, ±1 (10) | 2.28× (10/10) | 1.77× (10/10) | 1.67× (10/10) | 1.29× (7/10) |
| Toroidal, ±1 (6) | 4.61× (6/6) | 4.45× (6/6) | 1.43× (5/6) | 1.38× (4/6) |
| Planar-like, +1 (12) | 1.42× (10/12) | 1.16× (7/12) | 0.62× (1/12) | 0.51× (0/12) |
| Planar-like, ±1 (8) | 1.46× (8/8) | 1.27× (7/8) | 1.56× (8/8) | 1.36× (7/8) |
| All 51 | **1.91× (48/51)** | 1.51× (41/51) | 1.20× (37/51) | 0.95× (25/51) |
| All 51, extended grids (Amendment 1) | **2.13× (50/51)** | 1.62× (40/51) | 1.22× (33/51) | 0.93× (24/51) |

Footnote for comparison: K2000 on the V80 (measured, 12× v6.4, primary estimator; `fpga/v80_sca/RESULTS_HW.md`, Amendment 4) gave 3.74× vs plain SCA and 2.58× vs TEC for the best frozen online-form schedule, and 2.71× vs re-selected baselines for the kT-form X5.

### Caveats

1. **Model, not hardware.** The cycle model is extrapolated to sparse graphs under the stated assumption, and the 12-engine estimator is model-based.
2. **Shared parameters are per weight class.** They were placed on synthetic graphs (not G-set). The +1 classes need large q (synchronous oscillation), a regime K2000 never posed.
3. **Grid edges remain after Amendment 1, mostly for the corrected rules** (`edge_report_ext.txt`):
   - Onsager-kT on planar +1 sits at κ = 2, q = 2.88σ, T0 = 0.9σ in 56 of 60 selections.
   - All three-parameter rules on random +1 sit at the lowest q (1.44σ) in about 50 of 75.
   - TEC on ±1 sits at its most negative J_v in 47 of 120.
   - Plain SCA is interior. The bias is therefore against the corrected rules, which makes the kT gains conservative.
4. **Amendment 1 is post hoc.** It was frozen before its own runs, but after the subset's original finals had been logged. Both analyses are reported.
5. **Hard instances.** Where p ≈ 0.01–0.1 (the 2000-node tori, planar ±1 2000), choosing S* on the 256-run finals is noisy. On G32, kT's extended confirmation at S* = 2000 was worse than the pre-registered one at 4000.
6. **Notes on the GPU study's heads-up** (received after the held-out runs had started; post hoc):
   - The two-cycle on all-J = −1 instances matches my calibration; my +1 grids use q ≥ 1.08σ for every rule.
   - Its q = 4 (0.89σ) with T 4 → 0.5 is a different working region from mine (q ≈ 1.4σ, T_fin ≈ 0.67σ); both work on G22.
   - Its κ = 1.75·d/N is about 50–500× weaker than my κ, which explains why it saw near-identical behaviour. My data show strong κ effects in the range 0.125–2.
   - λ_eff uses the same d/N scaling in both studies; my grids reach λ = 1 (pre-registered) and 2 (extension).
7. **Folders not edited by me.** The manuscript, other research folders and hashed files were not edited. All four imported modules match their recorded hashes, and the old Snowball path was never created. `fpga27/` changed timestamps during the session, but not through any of my commands.

### Files (all in research/gset_20261007/)

- **Protocol:** `PROTOCOL.md`, `PROTOCOL_AMENDMENT1.md`, `PROTOCOL.sha256` (every entry verifies).
- **Data and best-known values:** `manifest.json`, `data/G*`, `sources/` (BLS and Ma–Hao PDFs and text, the G-set index page), `bkv.json`, `build_bkv.py`, `verify_data.py`.
- **Code and verification:** `engine.py`, `others.py`, `verify_engine.json`, `verify_others.json`, `run_gset.py`, `run_gset_ext.py`.
- **Logs:** `logs/run_{A,C,A2,C2,B}.log`.
- **Raw per-run results:** `results/{budget,conf,ext_budget,ext_conf}/`.
- **Analysis:** `analyze_gset.py`, `results_summary.json`, `results_summary_ext.json`, `analysis.txt`, `analysis_ext.txt`, `selected_configs.json`, `degree_stats.json`, `edge_report_{subset,rest,ext}.txt`.
- **Exploratory calibration (synthetic graphs only):** `calib{,2,3}.{py,json,log}`.
