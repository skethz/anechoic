# Formal verification (Lean 4 + Mathlib, and Rocq/Coq 9.3), 7 October 2026

The statements are in [SPEC.md](SPEC.md). Each one is proved **independently** in both systems:

| Statement | Lean (`SnowballFormal/SnowballFormal/Basic.lean`) | Coq (`coq/Snowball.v`) |
|---|---|---|
| F1 rule equivalence | `F1_rule_equiv_plus/minus` | `F1_rule_equiv_plus/minus` |
| F2 linear-band response (derivative 1/2T inside, 0 outside; summed response n_lin/2T; E[s′] = f) | `F2a_deriv_inside`, `F2b_deriv_outside`, `F2_dresp_is_derivative`, `F2c_summed_response`, `F2_expected_spin` (Mathlib `HasDerivAt`) | same names (Stdlib `derivable_pt_lim`, ε–δ) |
| F3 units √N·G = n_lin/2T | `F3_units` | `F3_units` |
| F4 q_eff identity and cases | `F4_qeff`, `F4_qeff_cases` | `F4_qeff`, `F4_qeff_cases` |
| F5 popcount field update (RTL P/Q selection) | `F5_popcount_decision/init` | `F5_popcount_decision/init` |
| F6 residue banking (v6.2 mod 4, v6.4 mod 8, unique representation) | `F6_*` | `F6_*` |
| F7 cut–energy identity and per-edge contribution | `F7_cut_energy`, `F7_edge_contribution` | `F7_cut_energy`, `F7_edge_contribution` |

**Results.**
- **No unproved steps.** Neither development contains `sorry`, `admit`, `Admitted` or a user axiom.
- **Lean:** `lake build` succeeds. `#print axioms` (`SnowballFormal/Axioms.lean`) lists only Lean's standard axioms `propext`, `Quot.sound` and `Classical.choice`. F5 needs no `Classical.choice`.
- **Coq:** `coqc` and the independent kernel checker `coqchk` pass. `Print Assumptions` lists only the Stdlib axioms of the classical real numbers (`ClassicalDedekindReals.sig_forall_dec`/`sig_not_dec`) and functional extensionality. F5 and F6 are closed: they depend on no axiom.

**Not formalised.** The large-N dynamical mean-field derivation (generating functional and saddle point) is not formalised. That derivation identifies the echo coefficient with the lag-one response G(t, t−1). What is formalised:
- the response of the clamped rule is exactly `1/(2T)` in the band and `0` outside;
- the summed response is `n_lin/(2T)`;
- the conversion to raw units;
- the hardware identities.

**How to re-check.**
- Coq: `cd coq && coqc Snowball.v && coqchk -silent -o -R . "" Snowball`.
- Lean: `cd SnowballFormal && lake exe cache get && lake build && lake env lean Axioms.lean`.
- Toolchains: Rocq 9.3.0 (Homebrew), Lean `v4.35.0-rc4` with Mathlib (from `lean-toolchain`/`lake-manifest.json`).
