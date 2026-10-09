# PROTOCOL_PUB validation amendment V2: our TEC kernel against TEC Fig. 3(b) (8 October 2026, before this run)

The authors asked that TEC be taken from the COP study's validated implementation. That implementation (`research/reaim_benchmarks_20261008/solvers_a4.tec_seq`) has the same semantics as ours (`pub_methods._tec_glauber_kernel`):
- sequential Glauber updates in index order, one sweep per cycle;
- a current spatial field;
- the temporal field from the configuration at the end of the previous cycle;
- a geometric T.

To confirm this quantitatively, our kernel is run on the COP study's pre-declared TEC check (its PROTOCOL_ADDENDUM4.md V3), with the same targets and criteria.

## Run

- K2000.
- J_v ∈ {0, 30, 6, −6}.
- k_BT from 100 to 0.1, geometric over 3,000 cycles.
- 64 runs each.
- Seeds `SeedSequence([20261008, 611, J_v index])`.

## Measure and pass condition

- **Measure:** the first cycle at which the 64-run mean cut reaches 31,670.
- **Pass:**
  - J_v = 0 within ±15% of 1,392;
  - J_v = 30 within ±15% of 740;
  - the ratio within ±20% of 1.88.
- **Also reported:** the Fig. 3(a) ordering, and P(cut ≥ 33,000) at the end.
