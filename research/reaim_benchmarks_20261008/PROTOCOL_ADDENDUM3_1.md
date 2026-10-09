# Addendum 3.1: the papers' own output rules for APC-SCA and ReAIM (8 October 2026, 14:05 CEST)

## What changes

Addendum 3 kept PROTOCOL.md's rule "final state only" for every method. Two baseline papers define a different output, and Addendum 3 asks for strict adherence to the papers. **Their stated outputs are therefore adopted:**

| Method | Output rule | Source |
|---|---|---|
| APC-SCA | "Output: optimized spin states: argmin_{1≤s≤S+1} {H(σ(s))}": the lowest-energy state visited, including the initial state | Okonogi et al. 2023, Algorithms 1 and 2 |
| ReAIM ASA | "Output: best spin state x_best, best Hamiltonian H_best", updated after every run phase (lines 6–8) | Chiang et al. 2024, Algorithm 3 |

- All other methods return their final state, as their papers describe:
  - Neal: its samples are final states.
  - STATICA (Algorithm 1): "spin configuration".
  - SB: the sign of the final positions.
  - TEC does not state an output; the final state is used.
  - Our Onsager forms: final state, as on the hardware.
- Settings, budgets and everything else are unchanged from Addendum 3.

## Disclosure

- **Timing.** This addendum was written while the Addendum 3 finals were running.
- **What I had seen.** I had seen the Addendum 3 log lines of ReAIM ASA on TSP: valid tours in only 1–7% of final states.
- **What prompted it.** Re-reading the two papers' Output lines then showed that their algorithms return a best state, not the final state. ReAIM's x_best and APC's argmin are part of the algorithms as published.
- **What is reported.** The final-state results of Addendum 3 (seed entry 41) are kept and reported as well.

## Runs

- **Runner and seeds:** `run_a3b.py`, fresh seeds SeedSequence([20261008, 42, p, i, m, s]), output in `results_a3/final_bv/`.
- **Kernels:** `solvers_a3b.py` performs the identical dynamics (identical final states and flips; `verify_a3b.py`) and returns the paper's output state.
- **Table 8b** uses these outputs for APC-SCA and ReAIM ASA.
