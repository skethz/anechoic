# PROTOCOL_HW_V6, Amendment 5: STATICA's published K2000 schedules on the 12 × v6.4 image

Frozen on 8 October 2026, before any data from these schedules. At freezing time the board holds the v6.4 headline image, idle; it was restored at 09:08 UTC, UUID 8669ad38a456bcdf17ad8485fa0cec7a. Hashes are in `PROTOCOL_HW_V6_A5.sha256`. `PROTOCOL_HW_V6.md` and Amendments 1–4 are unchanged.

## Why

- The paper's "original plain SCA" baseline is the v6 protocol's P3: q 8, T0 30, S 1560. That is not STATICA's published setting.
- STATICA (JSSC 2021, Fig. 22, Eq. 8, geometric r_T) fixes q = 4.0 for K2000. The fairness audit (`research/fairness_audit_20261008/k2000_pub/hw_model.py`) found that the published point had never been measured on the board.

## What

- **Image and host.** The programmed `board_v64_e12_250mhz` image (12 engines, 250 MHz) and the frozen v6.4 host `build/host_sca_v64` (SHA-256 3935b5a1…, built 02:24 UTC on 7 Oct, before the original v6.4 runs).
  - No reprogramming. The run stops if the board does not report the v6.4 UUID.
  - Tables come from the host's own generator (`sca_ref.hpp` `make_tables`, geometric T from T0 to 5).
- **Launch shape.** As `scripts/run_hw_v6.sh` phase 2b: seed 20261004, `--t1 5`, tables resident, one trial per launch, 2052 trials (171 rounds × 12 engines), trial ids 0–2051.
- **Cohorts, in this order** (`scripts/run_hw_v6_a5.sh`):
  - **X5** (paper O1), the same-session control: `--q 8 --tecT 1.75 --ramp --t0 15 --steps 280`.
  - **SP_long** (STATICA published long point, plain SCA): `--q 4 --lambda 0 --t0 40 --steps 1560`.
  - **SP_short** (STATICA published short point, plain SCA): `--q 4 --lambda 0 --t0 30 --steps 560`.
- **Estimators** (`scripts/analyze_hw_v6_a5.py`, identical to `analyze_hw_v6.py`). Success means cut ≥ 33000 and the device and host scores agree.
  - Primary TTS99 = t_R · ln(0.01)/ln(1 − P_round), where t_R is the mean round time after round 0.
  - Secondary = t_R · ln(0.01)/(12 · ln(1 − p)).
- **Isolation.** No overlap with Amendment 4: its chain (`bringup_mb4b.sh`) is paused with SIGSTOP for the duration and resumed afterwards. No `ami_tool sensors` calls.

## Predictions

| ID | Prediction | Pass condition |
|---|---|---|
| A5-C | The control reproduces the original X5 run | X5: the same cut and flips per trial id as `results/hwv6_board_v64_e12_250mhz/rounds/X5.jsonl` (2052 of 2052) |
| A5-M (reported) | Board p matches the audit's reference model (k2000_pub) | \|p_board − p_model\| ≤ 1.96·sqrt(p̄(1 − p̄)(2/2052)) + 0.02, against SP_long 0.829 and SP_short 0.138. The model predicts SP_long at 0.309 ms primary and 0.067 ms secondary, and SP_short at 0.249 ms for both |

## Reporting

- For SP_long, SP_short and the X5 control: primary and secondary TTS99, p, P_round and round time.
- Primary ratios against O1 (X5): this session's run and the original run.
- Secondary ratios against O5 (X4, from the original records: secondary 0.0181 ms).
