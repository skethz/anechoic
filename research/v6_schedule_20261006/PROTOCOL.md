# v6 schedule re-selection for 9 engines (frozen 6 October 2026, before running)

**Question.** Do the v6 cost model and 9-engine rounds favour schedules other than the frozen O1–O4?

## Models

- **Cost.**
  - Cycles per trial: 0.2714·flips + 16.91·S + 1148.
  - Plus S + 300 cycles per launch.
  - Rounds: × 1.035 for the maximum over 9 engines, at 300 MHz.
- **Primary objective.** Eq. 9 on rounds, with P_round = 1 − (1 − p)^9. If P_round ≥ 0.995, TTS is one round (the measured P_round would be 1).
- **Secondary objective.** t_R · ln(0.01) / (9 · ln(1 − p)).
- **Simulator.** The ablation CPU model (`abl.run`).

## Procedure (`sweep.py`)

1. **Pilot.** 423 configurations, 256 runs each, seed 61001.
   - Onsager and TEC-T, ramped.
   - Short schedules: S 200–360.
   - Long schedules: S 560–1260.
2. **Selection.** For each family (Onsager, TEC-T) and each objective, take the top 3 by TTS at the Wilson 95% lower bound of p.
3. **Held-out.** The selected configurations, 2,048 fresh runs each, seed 61002 (shared starts, i.e. paired).

## Decision rule

A held-out configuration is proposed for the hardware run (as an amendment to PROTOCOL_HW_V6 written before any v6 board data) only if its held-out TTS beats the frozen-configuration prediction for the same objective by at least 10%:
- primary: 0.088 ms (O4);
- secondary: 0.042 ms (O1).

Otherwise the hardware run uses O1–O4 only.

## Addendum (written after the Onsager/TEC-T held-out results, before any baseline re-selection)

The plain and constant-TEC baselines get the same re-selection as Onsager and TEC-T: same cost model, objectives, Wilson-bound top 3 and 2,048-run held-out (`baselines.py`).

- **Plain grid:** q ∈ {2, 4, 6, 8}, T0 ∈ {15, 20, 30, 40}, S ∈ {200, 280, 360, 560, 960, 1560}.
- **TEC grid:** J_v ∈ {−8, −4, −2, 2, 4}, q ∈ {4, 8}, T0 ∈ {20, 30}, S ∈ {280, 360, 560, 960, 1560}.
- **Seeds:** pilot 61003, held-out 61004.

Any re-selected Onsager/TEC-T configuration added to the hardware run is compared only against these re-selected baselines. The frozen O1–O4 vs P1–P4 / T1–T3 analysis of PROTOCOL_HW_V6 is unchanged.

## Addendum 2 (before any v6 board data): selection for 12 engines

**Why.** The 12-engine image (275 MHz) changes the primary objective. With 12 engines, configurations with p ≳ 0.44 succeed in essentially every round, so the primary estimator floors at one round.

**Procedure.** The same pilot (`pilot.json`) is re-ranked with E = 12 at 275 MHz.
- **Round-maximum factor.** 1.112, computed from the known O1–O4 trial flips; it replaces the guessed 1.035.
- **Selection.** Top 3 per family (Onsager, TEC-T) by primary TTS at the Wilson lower bound.
- **Held-out.** 2,048 fresh runs each (seed 61005), in `holdout_e12.json`.

**Decision rule.** The best held-out configuration per family is added to the 12-engine hardware run if its held-out primary TTS is at least 10% below the best frozen O1–O4 12-engine prediction (O4: 0.080 ms).
