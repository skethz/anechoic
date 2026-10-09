# Amendment 1: grid extension for the four engine rules (frozen 7 October 2026, before any extension run)

## Timing

This amendment was written about 16:55–17:00 CEST (14:55–15:00 UTC) on 7 October 2026. At that point:

- **Part A was complete for the 12 subset instances:** pilots and the 256-run finals, all of which had been printed to `logs/run_A.log`.
- **Part A was still running** for the other 39 instances.
- **No extension run had started.**

The design below uses **pilot data only**: [edge_report_subset.txt](edge_report_subset.txt), produced by [edge_report.py](edge_report.py). The frozen [PROTOCOL.md](PROTOCOL.md) and its results are unchanged and stay the pre-registered primary analysis. This follows the ABL2 precedent in `research/ablation_20261005`.

## Why

The frozen grids gave TEC, Onsager-kT and Onsager-online only a 2 × 2 subset of plain SCA's (q, T0) values. On G-set that subset is edge-limited in every weight class (subset pilots, all five budgets):

| Class | Plain SCA (4 × 4 shared grid) | TEC, Onsager-kT, Onsager-online (2 × 2 shared subset) |
|---|---|---|
| Q+ (planar +1) | At the (q_max, T0_max) corner in 15/15 selections | q_max and T0_min in 15/15 for each rule. q = 1.08σ makes the synchronous dynamics oscillate: best pilot cut −2,100 to −2,800 |
| P+ (random +1) | Interior | q = 1.44σ oscillates (−3,800), so only q = 2.16σ is usable. kT is at T0_min in 15/15; TEC prefers J_v_min |
| M± (±1) | Mostly interior (T0_max in 10/30) | TEC: T0_max 30/30, q_max 24/30, J_v_min 15/30. kT: κ_max 27/30, q_max 20/30. Online: q_max 28/30, T0_max 22/30 |

**Why the calibration missed this.** The synthetic planar proxy has no high-degree hubs. Real G-set planar graphs have degrees up to 326 and need a larger self-coupling q to stop the oscillation.

**Consequence.** The frozen comparison is biased against the three-parameter rules relative to plain SCA. It is unbiased between TEC and the Onsager forms, since they share identical subsets.

## Extension (16 new points per rule, instance and budget; union = 32 points for every rule)

Values are in units of σ; κ is dimensionless and λ is in K2000 units (λ_eff = λ·(d̄/N)/(1999/2000)). T_fin per class is unchanged.

| Rule | P+ | Q+ | M± |
|---|---|---|---|
| Plain SCA: q × T0 (4 × 4) | q {1.8, 2.52, 3.6, 4.32} × T0 {1.125, 1.575, 2.25, 3.6} (interleaved and outward) | q {2.88, 3.6, 4.32, 5.76} × T0 {1.35, 1.8, 2.7, 3.6} (beyond the corner) | q {0.135, 0.18, 0.27, 0.36} × T0 {1.125, 1.575, 1.8, 2.7} |
| Shared (q, T0) block for TEC, kT and online (2 × 2) | q {2.16, 2.88} × T0 {0.9, 1.125} | q {2.16, 2.88} × T0 {0.9, 1.8} | q {0.18, 0.36} × T0 {1.35, 1.8} |
| TEC J_v (× block) | {−0.72, −0.36, −0.18, −0.09} | {−0.36, −0.18, −0.09, +0.09} (interior, unchanged) | {−0.18, −0.09, −0.045, +0.045} |
| Onsager-kT κ (× block) | {0.25, 0.5, 1.0, 2.0} | same | same |
| Onsager-online λ (× block) | {0.25, 0.5, 1.0, 2.0} | same | same |

Design notes:
- Each rule-parameter range moves one doubling in that rule's own edge direction: TEC towards more negative J_v where J_v_min was hit, κ up (κ_max hit on M±), and λ up (λ_max hit on P+ and Q+ in 8/15 and 9/15).
- Plain SCA gets 16 new points as well, so that every rule has 32. On Q+ they extend beyond its corner; on P+ and M± they interleave and extend its range.
- No new point duplicates an original one. `run_gset_ext.py` asserts this.

## Procedure

- **Extension pilot.** For each instance (all 51), engine rule and S: 64 runs per new point, seed SeedSequence([20261007, 4, g, method index, S index, 16 + i]).
- **Selection.** Highest pilot mean cut over the union of the original and extension pilots. Ties go to the lower union index, so original points come first.
- **Final.**
  - If the union selection is an original point, its original 256-run final is the held-out estimate. It is the same configuration and seed as in PROTOCOL.md.
  - Otherwise 256 new runs, seed SeedSequence([20261007, 5, g, method index, S index]).
- **TTS confirmation.**
  - S*₁ and S*₁₂ are chosen as in PROTOCOL.md, on the union finals.
  - 1,024 runs, seed SeedSequence([20261007, 6, g, objective]), shared by the four rules so that runs are paired.
- **Order.** Part A, then Part C (the original confirmations), then A2 (extension), then C2 (extension confirmations), then Part B.

## Reporting

- **Both analyses are reported.**
  - **Original:** the frozen PROTOCOL.md grids and confirmations, i.e. the pre-registered verdicts on G1–G5 and G1′–G5′.
  - **Extended (Amendment 1):** the same measures and hypotheses on the union grids.
- **Labelling.** The extended analysis is labelled as a pre-registered post-hoc extension.
- **Disagreement.** Wherever the two disagree, both are stated.
