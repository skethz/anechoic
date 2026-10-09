# G-set Table 8a baselines at their papers' settings: engine-model study and board pre-registration (8 October 2026, before any run)

Companion to `../k2000_pub/PROTOCOL_PUB.md` (same rule: baselines are not tuned by us). It covers Table 8a (plain SCA and TEC against the two Onsager forms, measured on the V80 in PROTOCOL_HW_MB Amendment 3). `research/gset_20261007` and the board results are not changed.

## What the baseline papers prescribe for G-set

- **STATICA** (Yamamoto et al., JSSC 2021). Its K2000 study is as read in `research/statica_reproduction_20261003`; the full JSSC text was not re-read today (paywall). It gives no G-set setting and no rule for scaling T or q to other weight scales or degrees.
  - Its K2000 recipe: clipped Eq. (7), q = 4, T_fin = 5, geometric T, T_init ∈ {30, 40, 50}, S ∈ {360, 560, …, 1560}, choice by TTS.
  - G-set and K2000 both have |J| = 1, so this recipe is applied unchanged. It is the "general recipe or default" case of the rule.
- **TEC** (Du et al., arXiv:2608.21753) evaluates K2000 and graphs of 4–10 nodes only. It gives no G-set setting.
  - Its K2000 values: J_v = +30 J, k_BT from 100 J to 0.1 J, in its own sequential Glauber dynamics. That dynamics cannot run on the engine.
  - On the engine, TEC therefore takes STATICA's recipe plus the published J_v = +30 J (field h + J_v s(t−1)).
- **What was run in Table 8a** differs from both. The G-set study (`research/gset_20261007`) tuned every rule per instance:
  - temperatures and q are scaled to σ = √(mean degree);
  - 16-point grids (32 with Amendment 1) per budget S ∈ {250, …, 4000}.
- **Context only.** The APC-SCA paper (Okonogi et al., IEICE 2023; by the STATICA group) runs SCA on G22/G30/G32/G35 with these settings:
  - q = λ/2 (λ = the largest eigenvalue of −J);
  - exponential T from 10 to 0.1;
  - S = 1,000 (its Fig. 1, Table 2).

  It is modelled as a third, context row ("SCA-apc"), with the engine's clipped rule in place of the paper's logistic one.

## Model study (no board action)

- **Engine model.** `research/gset_20261007/engine.py` (unedited), which is bit-identical to `abl.run`. Instances, targets (⌈0.99·BKV⌉) and the 12-engine primary estimator are those of the G-set study:
  - cycles = 0.1424·flips + 18.09·S + 886 at 250 MHz;
  - t_round = E[max of 12 trial times];
  - P_round = 1 − (1 − p)^12;
  - one round if P_round ≥ 0.995.
- **SCA-pub and TEC-pub.** For each of the 51 instances and both rules:
  - **Pilots:** the 21 configurations of STATICA's grid (S × T_init), 256 runs each, seed `SeedSequence([20261008, 501, g, rule, k])`.
  - **Selection:** the minimum 12-engine primary TTS99. If that is infinite everywhere, the highest pilot mean cut. Ties go to the lower grid index.
  - **Confirmation:** the selected configuration, 1,024 fresh runs, seed `[20261008, 502, g, rule]`.
- **SCA-apc.** One configuration per instance, 1,024 runs, seed `[20261008, 503, g]`.
- **Output.**
  - Per-instance predictions.
  - The model's Table 8a with these baselines: Onsager-κT and Onsager-online from the G-set study's 1,024-run confirmations, both grid variants (`selected_configs.json`).
  - Geometric means are taken over instances where both TTS values are finite. Counts of instances where a baseline never reaches the target are reported separately.
- **Execution.** gpu-host, at most 48 single-threaded processes, scratch only.

## Board pre-registration (after the model runs; then STOP)

- The selected SCA-pub and TEC-pub configurations are written as host-argument lines in the format of `fpga/v80_sca/PROTOCOL_HW_MB_A3_gset.txt` (6,156 trials, seed 20261004).
- They are written together with the model predictions to `PROTOCOL_GSET_PUB_board.txt` and `PROTOCOL_GSET_PUB_predictions.json`, and hashed before any board data.
- No board run is made here. A board run would be coordinated separately.
