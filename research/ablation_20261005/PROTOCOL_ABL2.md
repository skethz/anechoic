# Ablation extension ABL2: TEC-T grid edges (frozen 5 October 2026, after the ABL pilot and before any TEC-T held-out result)

## Why

In the ABL pilot, TEC-T (−κ·T·s(t−1)) was best at the upper edge of its grid: κ = 0.75 and T_init = 20. Best p rose monotonically with κ:

| κ | 0.3 | 0.45 | 0.6 | 0.75 |
|---|---:|---:|---:|---:|
| Best p | 0.145 | 0.242 | 0.461 | 0.672 |

A fair comparison needs the interior of TEC-T's grid, just as J2 extended Onsager's.

## Grid (TEC-T2)

κ ∈ {0.75, 1.0, 1.25, 1.5, 2.0} × ramp ∈ {off, on} × q ∈ {4, 6, 8} × T_init ∈ {20, 25, 30} × S ∈ {560, 960, 1560}: 270 configurations.

- Everything else is as in PROTOCOL_ABL.md: pilot of 256 runs per configuration (seed root 50011, spawned per configuration) and the same cost model.
- **Selection.** For E ∈ {1, 4, 16}, minimise TTS at the Wilson lower bound over the union of the TEC-T and TEC-T2 pilots.
- **Held-out.** 1,024 runs with the same paired seed as ABL (root 50002), so the results are directly paired with the ABL held-out runs.

## Decision

A1 is re-evaluated against the best TEC-T held-out result after the extension. Both the original A1 verdict and the extended one are reported.
