# A2: confirmatory test of theory-D predictions at N=2000 — pre-registered

3 October 2026. Written after [REPORT_D.md](REPORT_D.md) and **before** running any of the configurations below at N=2000. The frozen A sweep ([PROTOCOL_ABC.md](PROTOCOL_ABC.md)) was already running when this was written; A2 does not use its results for selection.

## Predictions under test (from REPORT_D §4)

| ID | Prediction |
|---|---|
| P1 | At λ=0, x(t) = n_lin/(4T²) ≈ 0.3, roughly flat over T ∈ [8, 25] (both STATICA schedules). |
| P2 | At λ=1, the run never orders: b = λx crosses 1 near T ≈ 16. |
| P3 | Threshold λ_c ≈ 0.8–0.9 at q=4. λ=0.7 succeeds at least as often as λ=0.5. λ=0.9 collapses. |
| P4 | q=8 with λ=0.9 beats q=4 with λ=0.5, and q=8 with λ=0 is worse than q=4 with λ=0. |
| P5 | Ramping λ to 0 over the final 30% of steps improves success over a constant λ. |
| P6 | **The rule** λ = κ/(4x₀) with κ=0.75 and the 30% ramp beats the best constant-λ setting at the same schedule. x₀ is measured from the λ=0 pilot runs of this protocol. |

## Design

- **Schedules:** short (S=560, T 30→5) and long (S=1560, T 40→5); clipped Eq. 7; uniform random thresholds.
- **Seeds:** root 30003. **Runs:** 1,024 per configuration; initial states paired across configurations within a schedule.
- **Configurations per schedule:**
  - q=4 with λ ∈ {0, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0}
  - q=4 with λ ∈ {0.5, 0.7} plus the ramp
  - q=8 with λ ∈ {0, 0.9}
  - the rule (P6)
- **Logged:** x(t) and b(t) means every 20 steps; flip-back fraction; final cut; flips.
- **Pass/fail:**
  - P1 holds if the median x over T ∈ [8, 25] is in [0.2, 0.4].
  - P2 holds if b crosses 1 at a temperature in [12, 20].
  - P3–P6 are judged on held-out success rates. A difference counts if the 95% intervals of the two success rates do not overlap; otherwise it is reported as "not resolved".
- **TTS:** cycle model as in PROTOCOL_ABC, plus 4 cycles per step when λ ≠ 0.
