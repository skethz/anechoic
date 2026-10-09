# Addendum 5.1: Table 8b includes TEC and STATICA, and every baseline is reported with its validation status (8 October 2026, ~15:30 CEST)

## Source

- **The change.** The authors made a change that overrides A5.4 and A5.5:
  - **TEC:** use the sequential p-bit implementation (`solvers_a4.tec_seq`: sequential single-spin Glauber updates within a time step, temporal field J_v σ_i(t−1)) at the paper's K2000 settings, J_v = 30J and T from 100J to 0.1J.
  - **TEC validation note:** "this reading reaches the paper's 95% target faster than reported (736/577 cycles against 1,392/740), with a smaller J_v gain" (the audit's numbers; ours, V3, are 739/577).
  - **Synchronous reading:** do not use it anywhere.
  - **Beyond K2000 (G-set Max-Cut, GP, TSP):** give TEC and STATICA numbers with their K2000 settings scaled to each instance by our σ rule. Label the rule as ours, and report no "settings not specified" cells.
  - **Unchanged:** ReAIM Table I k, Goto's Δt rule for bSB/dSB, aSB at its 2019 setting with its divergence reported.
- **What had been seen when this was written:**
  - all A4 validations, plus V8a: ReAIM K2000 P_a 0.523 at 4,096 iterations against 0.47 (pass), and 0.748 at 6,400 against 0.80, 0.002 outside the tolerance (fail);
  - the A3 STATICA and A4 TEC σ-transfer finals;
  - the exploratory `diag_statica.py` result: P_a 0.833 (n = 4,096) and 0.832 for the reference kernel (n = 2,048).
- **Still running:** the A5 ReAIM and SB finals.

## Changes to A5.6 (Table 8b composition)

| Method | Max-Cut G1–G20 | GPP (9) | TSP (6) |
|---|---|---|---|
| SA (Neal defaults) | A3 | A3 | A3 |
| STATICA, K2000 setting (q 4, T 40 → 5) × σ_inst/σ_K2000 **[our transfer rule]** | A3 `SCA` | A3 | A3 |
| TEC, sequential p-bit, K2000 setting (J_v 30J, T 100J → 0.1J, geometric) × σ_inst/σ_K2000 **[our transfer rule]** | A4 `tec_seq` | A4 | A4 |
| APC-SCA | A3.1 (paper's Max-Cut settings) | A3.1, TSP-column settings **[our choice]** | A3.1 (paper's TSP settings) |
| ReAIM ASA (Table I k) | A5.1 | A5.1 | A5.1 |
| aSB (Δt 0.9, M 2; non-finite fraction reported) | A3 | A3 | A4.1 `aSB_dt09` |
| bSB, dSB (Δt by TTT) | A5.2 | A5.2 | A5.2 |
| Onsager-κT, Onsager-online (ours) | A3 | A3 | A3 |

- **σ rule:** σ_inst = sqrt(mean_i Σ_j J_ij²) of the instance's full Ising matrix (penalty couplings included), σ_K2000 = √1999.
- **Seeds and freshness.** The TEC and STATICA numbers come from finals that were already run with fresh seeds under their own pre-registration: A4 seeds [.., 51, ..] and A3 seeds [.., 41, ..]. No setting is selected from results, so they are used as they are.

## Change to the usage rule (A4 §4 and A5.7)

- **Old rule:** "a baseline whose validation fails is shown as 'not validated', and its numbers are not used".
- **New rule:** the authors' decision to include TEC, whose Fig. 3b check failed, sets the policy that every baseline is reported with its validation status. It is applied to all baselines alike:
  - **TEC:** V3 failed (faster than reported, smaller J_v gain).
  - **ReAIM:** V8a reproduced the 0.15 ms point; the 0.23 ms point mapped to 6,400 iterations is 0.002 outside the tolerance.
  - **SA:** V1a matches the real library; V1b is within 0.5% on 9/11 instances, and the real library shows the same two deviations.
  - **STATICA:** V2's pre-registered n = 1,024 sample failed on P_a (0.867); the exploratory n = 4,096 sample is inside the paper's interval (0.833).
- **Disclosure.** This rule change was made after the V8a result was seen.
- **Strict version also reported.** Under the old rule, ReAIM, SA, STATICA and TEC would have no numbers.

## Unchanged

A5.1–A5.3 and A5.7; `run_a5.py`, `validate_a5.py`.

## Files

- **`analyze_a5.py`:** written after the runs; this addendum fixes its composition. It no longer has a "not specified" or "not reproduced" row, and the synchronous TEC reading appears nowhere.
- **`VALIDATED_BASELINES.md`:** per method, the settings, source, reproduction (ours against the paper's) and verdict. The σ transfer is labelled as our rule.
