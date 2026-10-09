# PROTOCOL_HW_MB, Amendment 2: revision-2 image (cut-word fix)

Frozen on 7 October 2026, before any board data from the revision-2 image (its build had not started). Hashes are in `PROTOCOL_HW_MB_A2.sha256`. The main protocol and Amendment 1 are unchanged, and their results on the first image (`board_mb2_e12_250mhz`) stand as recorded.

## Finding on the first image

In the main protocol's G22_ons cohort, 3 of the 2,052 trials ended in the all-equal period-2 state: cut 0, and Σ_i s_i h_i = −2·sumw < 0.
- The device's cut word was 2^62, while the host's independent rescoring gave 0. The protocol counts such trials as integrity failures, and they are not successes in any case.
- Spins and flips matched the reference model exactly; this was rechecked on the CPU for all three trials.

**Cause.** `sca_core_mb.v` inherits the cut halving from v6.4b, `(x + {63'd0, x[63]}) >>> 1`. The concatenation is unsigned, which makes the whole expression unsigned, so `>>>` shifts logically.
- The result is wrong exactly when Σ s·h < 0, which means cut < sumw/2.
- On K2000, sumw = −1040 and cut ≥ 0, so this never happens: v6.4b and every K2000 result are unaffected.
- On +1-weight G-set graphs it can only happen in failed, low-cut states.

**Simulation.** The negative-score suite (`gen_vectors_mb neg`: G22 with q = 1; all 8 trials end at cut 0) has exactly 8 errors on `sca_core_mb.v`, all of them the cut words (0x4000000000000000 against 0).

**Fix.** `src/v6/sca_core_mb_r2.v` makes the sign operand `$signed`. Nothing else changes: same module, interface, pipeline and cycles.
- It passes the negative-score suite.
- It also passes the regression: K2000 short suite, O4, X5, G22 short/long, G32 short and random K = 4, with cycle counts identical to the first image's simulation.

## Revision-2 image and run

**Build.** `board_mb2r2_e12_250mhz` (`scripts/build_mb2r2.sh`), with the same front end (variant mb2r2, rebuilt from the unchanged source), integration (`scripts/integrate_kernel_mb_r2.tcl`) and settings as the first image.
- If it does not close timing at 250 MHz, there is no revision-2 board result, and the first image's results stand with the caveat above.

**Run.** `scripts/bringup_mb2r2.sh`, unattended.
1. It waits for the build and for the first image's chain to end.
2. It requires timing PASS and all three hash files, backs up the staged image, and programs the revision-2 image.
3. It then runs, unchanged:
   - the main protocol's phases 1–3 (`run_hw_mb.sh`: exact gate W1–W4 and G22a/b, G32a/b; K2000 X5 and O4; the six contingency G-set cohorts);
   - Amendment 1 (`run_hw_mb_a1.sh`: 21 exact-gate sets and 100 cohorts);
   - two negative-score exact-gate sets (`run_hw_mb_a2neg.sh`: G22, q = 1, 12 trials each, ids 940000+).
4. Finally it restores the v6.4 headline image. Board power is not re-measured, because only the cut logic changed.

## Predictions (pass conditions)

| ID | Prediction | Pass condition |
|---|---|---|
| A2-E | Exact on hardware, including negative scores | Every exact-gate trial of the three runs is bit-exact, including the 24 negative-score trials (cut 0) that the first image's cut word fails in simulation |
| A2-R | Same results as the first image | `compare_r1_r2.py`: every trial id of X5, O4, the six main G-set cohorts and the 100 Amendment 1 cohorts has identical spins, flips and independently rescored cut on both images |
| A2-C | Device cut word now correct | No integrity failure: device cut = rescored cut for every revision-2 trial |
| R1 | Same algorithm as v6.4 | As in the main protocol: X5 and O4 identical to the 12 × v6.4 run |
