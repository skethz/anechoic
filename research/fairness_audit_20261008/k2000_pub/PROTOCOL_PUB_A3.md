# PROTOCOL_PUB Amendment A3: aSB at Goto et al. 2021's K2000 aSB time step (8 October 2026; before this run)

## Why

At the 2019 K2000 setting (Δt = 0.9, M = 2), our aSB reproduces the paper's own number: N_step = 186 gives a mean cut of 32,756 against 32,768 and P_a = 0.020 against 0.04 (validation V1, `asb_186`). In floating point, however, it diverges in 1%, 15%, 63%, 82% and 100% of runs at S = 250–4000 (Amendment A1). The paper's FPGA used saturating 16-bit fixed point (its supplement), and our float model has no saturation. The COP study diagnoses this as an explicit-kick instability and adds Δt = 0.5, M = 2 as its option (b) (`research/reaim_benchmarks_20261008/PROTOCOL_ADDENDUM4.md`).

Δt = 0.5 is also the authors' own aSB time step on K2000 in Goto et al. 2021, Fig. 2A ("Δt is set to 0.5 (aSB)"). That figure supplies a reported result to validate against.

## Runs (`asb_dt05.py`; gpu-host; fresh seeds)

- **Finals.** K2000, Δt = 0.5, M = 2, ξ0 = 0.7/√N, p linear from 0 to 1, x(0) = 0, y(0) uniform in (−0.1, 0.1). S ∈ {250, 500, 1000, 2000, 4000}, 256 runs each, seed `SeedSequence([20261008, 304, 7, S index])`. A run with non-finite dynamics counts as a failure (A1 rule), and the number of such runs is reported.
- **Validation.** N_step ∈ {100, 1000, 10000}, 256 runs, seed `[20261008, 612, index]`. The targets are Goto 2021 Fig. 2A's aSB line, read from the figure: 32,030, 32,950 and 32,975 (±20). The pass condition (the COP study's rule) is |ours − read value| ≤ 40 + 3 standard errors.

## Use

The A3 aSB replaces the provisional aSB in the Figure 2 data only if it passes. The 2019-setting A1 numbers are reported alongside it.
