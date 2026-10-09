# v5 engine measurement: addendum to PROTOCOL_HW.md (frozen 5 October 2026, before any v5 board data)

v5 changes the engine only:
- 256 decision lanes (8 rounds per step).
- 4 interleaved flip extractors with a duplicated coupling memory.
- Burst-loaded parameter tables.

The algorithm is unchanged except for the random-number mapping. Lane l now serves spins l, l+256, and so on. In the fixed-point reference, success on O1 is 0.932 at 128 lanes vs 0.933 at 256 lanes (2,048 trials each, `results/quality_lanes.txt`).

## Procedure

Identical to PROTOCOL_HW.md, run with `SCA_VARIANT=v5`, the host built with `-DSCA_LANES=256`, and a fresh result directory `results/hw_board_v5_<f>mhz`:
- **Image selection:** the highest-clock image that passes timing.
- **Phase 1 exact gate:** V1–V4, compared against the 256-lane reference.
- **Phase 2:** the same 11 configurations × 1,024 trials, plus single-run sets for P1, T1 and O1.
- **Analysis:** the same script.

## Pre-registered predictions

These come from the v4 measurements and the v5 design. v4 measured `cycles = 0.522·flips + 40.2·S + 1884` per run, and table loading at 35.3 cycles per step per launch.

| ID | Prediction | Pass condition |
|---|---|---|
| H1 | The flip coefficient drops (4 extractors) | Fitted cycles per flip ≤ 0.33 |
| H2 | The per-step coefficient drops (8 rounds instead of 16) | Fitted cycles per step ≤ 34 |
| H3 | Table loading becomes cheap | Per-launch overhead ≤ 6 cycles per step (single vs batched) |
| H4 | O1 runs faster | Measured O1 t_run ≤ 0.27 ms at 300 MHz (v4: 0.366 ms), i.e. at least 1.35× faster |
| H5 | Quality is unchanged | For every configuration, |p_v5 − p_v4| ≤ 0.05 or the 95% CI of the difference contains 0 |

If the image closes timing only below 300 MHz, H4 is evaluated at the achieved clock: t_run ≤ 0.27 ms × 300/f.

## Addendum: v5.1 (5 October 2026, before any v5 or v5.1 board data)

- **What changes.** v5.1 changes only the extractor's refill implementation: it uses a precomputed queue to close timing (HLS estimate 297 MHz vs 206 MHz for v5). Extraction order and per-cycle behaviour are identical, so cycle counts are identical, and bit-exactness to the 256-lane reference is unchanged.
- **Measurement order.**
  - v5 is measured at its best timing-clean clock (250 MHz; its 300 MHz build failed timing). It tests H1–H3 and H5.
  - v5.1 is then measured with the same procedure at its best timing-clean clock, in `results/hw_board_v51_<f>mhz`. It tests H4 at that clock and H5 again.
- **Additional prediction H6.** v5.1's fitted cycle coefficients equal v5's within the fit uncertainty, since the cycles are the same.
