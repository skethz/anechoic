# V80 check of the popcount-free rule (TEC-T) vs online Onsager (frozen 5 October 2026, before any run)

## Question

The CPU ablation (ABL2/ABL3) found that a fixed temperature-proportional anti-coupling with end ramp (TEC-T) comes within 4–9% of the online Onsager correction. Does this hold on the board?

## Image and host

- **Image:** v5 at 250 MHz (UUID `ce847edb…`), already programmed. It is the same image as the v5 cohorts in RESULTS_HW.md.
- **Host:** `host_sca_v5t`. It is the v5 host plus `--tecT κ`, which sets kconst[t] = κ·T_{t−1}·ramp(t) in Q16.16 with kcorr = 0. No hardware change is involved.

## Runs

| Phase | Configuration | Trials | Notes |
|---|---|---|---|
| Exact gate | TT: TEC-T κ1.25 + ramp, q8, T 20→5, S960 | 4, two per launch, offsets 900400+ | Every trial must match `sca_ref.hpp` exactly with the TEC-T tables (spins, cut, flips, n_lin trace) |
| Cohort | TT | 1,024, 64 per launch, seed 20261004, offsets 0–1023 | Paired with the existing v5 O1 cohort (same seed and offsets, so same initial spins) |

## Analysis and predictions (from CPU ABL2 held-out)

The same analysis as PROTOCOL_HW is used: t_run is cycles per trial amortised over the launch, and device TTS = t_run·⌈R99⌉.

| ID | Prediction | Pass condition |
|---|---|---|
| R1 | TT success matches the CPU estimate | p_TT within the 95% interval of CPU p = 0.914 (1,024 runs), i.e. 0.89–0.93 |
| R2 | TT per-run cycles are slightly above O1 | Cycle ratio TT/O1 within 1.0–1.2 (CPU flips ratio 1.16) |
| R3 | The two rules are equivalent at E = 1 | Measured TTS ratio TT/O1 within 0.9–1.3 (CPU 1.09) |
