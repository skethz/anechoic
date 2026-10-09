# Addendum 2: fair re-tuning by automatic grid-edge extension (frozen 8 October 2026, before any A2 run)

## Why

**The finding.** The authors required that every comparison be fair, and that any comparison found unfair be done again. Their review of `edge_report.txt` found that many selected configurations of the first run (A1) sit at a grid end. Examples:

- **Engine rules:** T0 at its minimum and T_fin at its maximum for the SCA-family rules on TSP. On GPP: T0, q and J_v ends.
- **APC-SCA:** q_reset, r_q, T0 and q_lim at their ends.
- **ReAIM:** the k set at its minimum on TSP.
- **SB:** aSB's dt and ξ (all-infeasible ties) and dSB's dt.

**ReAIM's own options.** ReAIM's F ∈ {min, max} (its Step 4) was fixed to max, and its published TSP temperature schedule was not used.

**What this addendum does.** It re-tunes all ten methods on all three problems with one procedure. It uses fresh seeds and a fresh results directory, and leaves every A1 result untouched.

## Unchanged from PROTOCOL.md

- **Problems and formulations:** problems, instances, formulations, penalties (GPP P = 4; TSP A = max distance), the bias handling, and the methods' implementations. ReAIM is the exception (options below).
- **Budgets:**
  - Max-Cut S ∈ {250, 500, 1000, 2000, 4000};
  - GPP S ∈ {256, …, 4096};
  - TSP S ∈ {512, …, 8192}.
- **Runs:**
  - pilot: 64 runs per grid point and budget;
  - selection: highest pilot mean normalized quality (infeasible = 0; for Max-Cut, cut/BKV, which ranks exactly like the G-set study's mean cut), ties to the lower point index;
  - held-out final: 256 runs;
  - measures and targets: as before.
- **Max-Cut G1–G20:** all ten methods are now re-tuned here with the same procedure.
  - Base: the G-set study's pre-registered 16-point grids per weight class: P+ = G1–G5, Q+ = G14–G17, M± = G6–G13 and G18–G20.
  - These are the grids Table 8b used.
  - My solvers run the G-set study's code path bit for bit (verify.py part B).

## The procedure (identical for every method; `grids_a2.py`, `run_a2.py`)

### Grid families and base grid

1. **Grid families:**
   - `gpp`, `tsp`: relative units of `grids.py`;
   - `mcp:P+`, `mcp:Q+`, `mcp:M+-`: relative units of `research/gset_20261007/run_gset.py`.
2. **A1 grid:** the frozen grid of the first run. All A1 points are kept, in their A1 order (verified: `verify_a2.py`).
3. **A2 base additions.**
   - **ReAIM ASA:**
     - F ∈ {max, min}: the |N|-FIFO reduction of Algorithm 2, line 22 (ReAIM Step 4). F = max is the A1 value.
     - Its published temperature schedule (Table II), added for every k set and F: GPP T 1 → 0.01; TSP T 0.5 → 0.1. The Max-Cut schedule 1 → 0.1 is already on its grid.
     - `solvers_a2.reaim_a2` adds F and T0. With F = max and T0 = 1 it is bit-identical to the frozen ReAIM.
   - **APC-SCA on GPP and TSP:** a T_fin axis {0.15, 0.6}. Every other SCA-family rule there already had these values; APC's T_fin was fixed at 0.15.

### Axes and edge counts

4. **Axis:** a numeric parameter with at least two values. Categorical F and fixed parameters are not axes and are not extended.
   - Fixed by design and stated as such:
     - the penalty weights;
     - the SCA-family T_fin on Max-Cut (fixed per class for all five SCA-family rules in the G-set study);
     - the SB internal constants;
     - ReAIM's ITER_trial and ITER_run (32/96; not published);
     - the hardware ramp.
5. **Edge count:** for each family, method and axis, the fraction of the family's selections (instances × 5 budgets) at the axis minimum and at its maximum.
   - A selection whose pilot means are identical at every grid point (an all-tied grid, for example all runs infeasible) carries no direction. It counts at both ends of every axis.

### Extension rule

6. **Trigger:** an end with fraction > 0.20.
7. **Extension:** one new value beyond that end with the existing spacing.
   - **Geometric**, v_new = v_end²/v_neighbour, when the end value and its neighbour have the same sign. This covers positive axes and the all-negative J_v axes.
   - **Linear**, v_new = 2·v_end − v_neighbour, when they differ in sign.
   - **APC decay factor r_q:** geometric in (1 − r_q).
   - **ReAIM k set:** halve every flip cap (low end, floor 1; it stops at (1,1,1,1)), or double it (high end, capped at N).
8. **Constraints on new points only:**
   - T1 ≤ T0 (SA, ReAIM);
   - T_fin ≤ T0 (SCA family);
   - q_lim < q_reset (APC).
   - New points that break a constraint are not added. A1 points are kept even where the A1 grid had a rising schedule.
9. **Cap:** a (method, family) grid holds at most A1 size + **224** points, the same cap for every method and family.
   - Every A2 addition counts towards it, including ReAIM's F = min copies and published-schedule points and APC's T_fin values.
   - Within a round, triggers are applied in decreasing order of their fraction (ties: axis order, low end first).
   - An extension that would exceed the cap is skipped and logged.

### Rounds

10. **Rounds:** up to three.
    - Round 1 is triggered by the A1 selections.
    - Rounds 2 and 3 are triggered by the A2 selections after the previous round (the union of all A2 pilots so far).
    - **Round 1 pilots every point of its grid with fresh seeds**, including the A1 points. Rounds 2 and 3 pilot only their new points.
    - Residual edges after round 3 are reported.
11. **Held-out finals:** after round 3, 256 runs per (instance, method, budget) of the final union selection.

## Seeds and files (fresh)

- **Pilots:** SeedSequence([20261008, 20 + R, p, i, m, s, point index]), round R ∈ {1, 2, 3}.
- **Finals:** SeedSequence([20261008, 31, p, i, m, s]).
- **Indices:**
  - p = 0 (Max-Cut), 1 (GPP), 2 (TSP);
  - i = G-set number, or 100 + TSPLIB index;
  - m = method index;
  - s = budget index.
- **Disjointness:** earlier runs used, as second entry, 1, 2 and 7–11 (this study); 77 and 99 (checks and references); and seed 20261007 in the G-set study.
- **Files:**
  - `results_a2/grids_r<R>.json`: grids, extension logs and edge counts;
  - `results_a2/pilot/r<R>/...`;
  - `results_a2/final/...`.
- A1 results (`results/`, `research/gset_20261007/results/`) are not modified.

## Measures (Table 8b, recomputed; old and new reported side by side)

- **Quality and feasibility:** mean normalized quality and feasibility at the ReAIM-matched budget (Max-Cut S = 4000, GPP 4096, TSP 8192).
- **Steps to the 1% target:** MCS99, the minimum over S per instance.
  - Max-Cut: the geometric mean over the instances that **every** method solves (finite MCS99), plus solved counts out of 20.
  - GPP and TSP: the same; if the common set is empty, per method over its own solved instances, with counts.
- **Edge counts** after re-tuning.
- **Conclusions re-checked:**
  - **C1:** the engine rules fail on GPP (3–6% balanced at S = 4096).
  - **C2:** engine-rule TSP quality ≈ 0.48, and bSB is the only strong TSP method.
  - **C3:** SA and dSB need 2.6× and 1.5× fewer steps than Onsager-online on the common Max-Cut set (A1: 2.57× and 1.53× over the 11 common instances).
  - **C4:** TSP needs K = 6.
- **Precision study.** If the A2 selections change an engine rule's precision-study schedule (the S_exp rule of `precision.py`), the pre-registered precision study is re-run with the A2 schedules. The same code, tolerance and K_min apply. The hardware seed is 30261008 + 1000·p + 10·(instance code mod 100) + rule index, and the output goes to the fresh `precision_a2/`. It is run for all three problems.

## Platform

- **Host:** gpu-host, Python 3.12.3, numpy 1.26.4, numba 0.60.0 (venv under `/scratch/USER/anechoic_cpu_20261008`), at most 72 worker processes because the machine is shared.
- **Verification:** `verify_a2.py` (on gpu-host before the first A2 run) checks:
  - ReAIM option bit-identity;
  - A1-point equality and order against `grids.py` and `run_gset.py`;
  - a deterministic plan.
