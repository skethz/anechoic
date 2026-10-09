# PROTOCOL_PUB Amendment A6: K2000 bSB and dSB with one Δt chosen by time-to-target (8 October 2026, ~15:50 CEST, before any A6 run)

## Why

- **The authors' request:** re-run K2000 bSB and dSB with the Δt rule used for the COP study's Table 8b (`research/reaim_benchmarks_20261008/PROTOCOL_ADDENDUM5.md`, A5.2), so that Figure 2 and Table 8b agree.
- **The two rules:**
  - **A2** (current Figure 2 data) picks Δt per budget by the highest pilot mean cut.
  - **Table 8b** picks one Δt per problem by time-to-target, as Goto et al. 2021 prescribe: Δt "set to the best value for each problem among five values (0.25, 0.5, 0.75, 1, and 1.25)"; "N_step is also optimized for TTT or TTS separately".

## Rule (identical to PROTOCOL_PUB's frozen `run_pub.sb_select`, applied to fresh pilots)

For each of bSB and dSB:
- **Pilots:** 256-run pilots at every S ∈ {250, 500, 1000, 2000, 4000} and every Δt ∈ {0.25, 0.5, 0.75, 1, 1.25}.
- **Selection:** TTT(Δt) = min over S of S·ln 0.01/ln(1 − p), with p = P(cut ≥ 33,000), Figure 2's success target for every method. The selected Δt minimises TTT; ties go to the higher pilot mean cut at the S attaining the minimum, then to the smaller Δt.
- **Relation to Table 8b:** this is Table 8b's A5.2 rule with Figure 2's target in place of the G-set target ⌈0.99·BKV⌉. The MCS99 convention is that of both studies' code: S if p = 1, uncapped for 0.99 < p < 1.
- **Sensitivity** (reported, never used for selection):
  - target 33,004 = ⌈0.99 × 33,337⌉ (Goto's "99% of the best known");
  - Goto's TTT convention, which caps the step count at T_com = S when P_S > 0.99.

## Runs (`sb_ttt_A6.py`; gpu-host, ≤ 72 processes, scratch only)

- **Kernel and format:** PROTOCOL_PUB's kernel and format, unchanged: `pub_methods.run_vectorized` → frozen `methods.sb_traj` (a0 = 1, a(t) = t/S, c0 = 0.5/(⟨J⟩√N), x, y uniform in (−0.1, 0.1), inelastic walls, ξ = 1), and `pub_methods.summarize_final`.
- **Pilots:** 256 runs per (method, Δt, S), `SeedSequence([20261008, 701, code, S index, k])`.
- **Finals:** 256 runs per (method, Δt, S), all five Δt, `SeedSequence([20261008, 702, code, S index, k])`. Figure 2 uses the finals at the selected Δt; the rest are kept.
- **Codes:** bSB 8, dSB 9 (`pub_methods.CODE`).
- **Seeds:** the 7xx families are unused anywhere in `research/`.

## Output

- `results_pub_A6/S<S>.json`: `results_pub_A2/S<S>.json` with the bSB and dSB entries replaced by the A6 finals at the selected Δt (same keys: `final` with H_mean/H_p10/H_p90, cuts, p33000, mean_cut, mcs99, best_visited, job, sec; and `setting`). Every other entry is copied unchanged, and this is asserted.
- `selection_A6.json`, `claims_A6.txt`, `compare_A6_vs_A2.json`.
- Raw runs: `results_pub_A6/raw/`.
- Nothing else in this folder is modified.

## Checks (pre-declared)

- **Kernel equivalence with the COP study's validated SB** (V7, `research/reaim_benchmarks_20261008/validation_a4/sb21_*`: `solvers.sb`, K2000, Δt = 1, S = 1000, 256 runs):

  | | V7 mean cut | SE |
  |---|---|---|
  | bSB | 33,217.2 | 0.36 |
  | dSB | 33,116.3 | 5.83 |

  The A6 finals at Δt = 1, S = 1000 must agree within 3 combined SE. Against Goto 2021 Fig. 2A (N_step = 1000, read ≈33,215 for bSB and ≈33,110 for dSB), they must fall within 40 + 3 SE.

## Disclosure

- **Already known when this was written:** the PROTOCOL_PUB selection from its seed-301 pilots, Δt = 1.25 for both bSB and dSB (`results_pub/summary.json`).
- **Why fresh runs anyway:** the authors asked for fresh seeds and a fresh run. The rule is fixed above and nothing is chosen by us.
