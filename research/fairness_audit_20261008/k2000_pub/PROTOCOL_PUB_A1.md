# PROTOCOL_PUB Amendment A1: per-run divergence accounting for aSB (8 October 2026; written after the PROTOCOL_PUB finals, before this computation)

## Observation that prompted it

Under PROTOCOL_PUB, aSB at its published K2000 setting (Goto et al. 2019: Δt = 0.9, M = 2, ξ0 = 0.7/√N) was flagged invalid at every budget. The frozen `methods.run_method` rule is applied per batch: if any of the 256 runs produces a non-finite x or y at any step, every cut of the batch is set to −∞.

The paper itself reports the modified symplectic Euler method as stable "when Δt ≤ 1 (for large M)". With M = 2 the float32 integration overflows in some runs.

## Change (analysis only; the PROTOCOL_PUB results are kept as pre-registered)

- **Re-simulation.** The aSB finals are re-simulated with exactly the same seeds (`SeedSequence([20261008, 302, 7, S index, 0])`, 256 runs) and the same arithmetic (`methods.asb_traj`). The only addition is a per-run flag that records whether that run's x or y ever became non-finite.
- **Per-run rule.** This is the rule of `research/gset_20261007/PROTOCOL.md`: a run with non-finite dynamics counts as a failure and is excluded from the mean-cut statistics. The number of such runs is reported.
- **Identity check.** The non-diverged runs must give exactly the same final cuts and energy trajectories as the PROTOCOL_PUB finals. Runs are independent rows of the batch.
- **Reporting.** Both the pre-registered (batch-invalid) and the per-run results are reported. The Figure 2 data use the per-run result, marked as Amendment A1.
