# Multi-engine measurement protocol: v5.2 × 4 engines (frozen 5 October 2026, before any multi-engine data)

## Design under test

- **Engine.** v5.2 = v5.1 plus one packed, burst-read table array. It has the same cycles per step and per flip as v5/v5.1, is bit-exact against the 256-lane reference, and uses host `host_sca_multi.cpp`.
- **Replication.** `NUM_ENGINES=4` copies of the same verified IP in one image.
  - Each copy has its own control window at BAR0 0x01100000 + e·0x1000, its own passive cycle counter, and its own output and trace regions.
  - All copies share the coupling and table inputs in HBM.
- **Image selection.** The highest timing-clean clock among the builds attempted (300 and 275 MHz).

## Phase 1: exact gate on all engines

`--verify --engines 4`. Every trial on every engine must match `sca_ref.hpp` exactly (spins, cut, flips, n_lin trace). Seed 20261004.

| ID | Configuration | Trials | Trials per launch |
|---|---|---|---|
| W1 | Plain q4, T 30→5, S40 | 8 | 1 |
| W2 | Onsager λ0.7 + ramp, q4, T20, S560 | 8 | 2 |
| W3 | TEC J_v=−8, q8, T20, S40 | 4 | 1 |
| W4 | O1: Onsager λ1.05 + ramp, q8, T12, S960 | 8 | 1 |

Each host invocation loads the couplings on its first round only.

**Any mismatch stops the protocol.**

## Phase 2

**2a. Concurrency.** O1 with 4 engines × 256 trials, 64 trials per launch per engine. Per-engine cycles per trial are compared with v5's measured fit, 0.2754·flips + 31.06·S + 753.

**2b. Measured 4-engine rounds.**
- **Configurations.** All 11 PROTOCOL_HW configurations (P1–P4, T1–T3, O1–O4).
- **Trials.** 1,024 per configuration (256 rounds × 4 engines), 1 trial per engine per round, trial offsets 0–1023.
- **Round time.** The maximum of the 4 engines' start-to-done cycles. This includes per-launch table loading. Couplings are loaded in round 0 only, and round 0 is excluded from timing.

## Analysis (fixed in advance)

- **Trial level.** p and Clopper–Pearson 95% interval, as before.
- **Round level.**
  - P_round = fraction of rounds with at least one success, from 256 rounds.
  - Mean round time t_R over rounds 1–255.
  - **Measured 4-engine TTS99 = t_R · ln(0.01)/ln(1 − P_round)** (Eq. 9 applied to rounds).
  - Whole-round form: t_R · ⌈R99⌉.
- **Independence check.** P_round compared with 1 − (1 − p)^4.
- **Primary.** The best measured 4-engine TTS among O1–O4, compared with the best among P1–P4 and among T1–T3. These are post-hoc minima, which favours the baselines. Bootstrap 95% intervals use 10,000 resamples of rounds.

## Predictions

| ID | Prediction | Pass condition |
|---|---|---|
| M1 | No slowdown from concurrency | Per-engine cycles per trial within 2% of the v5 fit (2a) |
| M2 | Engines are independent | For every configuration, P_round is inside the 95% interval implied by p |
| M3 | Large multi-engine gain | Best measured 4-engine Onsager TTS99 ≤ 0.20 ms × (300/f) |
| M4 | Onsager advantage persists at 4 engines | Best Onsager / best plain and best Onsager / best TEC both ≥ 2× (point estimate) |

**Basis for M3.** O2 is predicted at about 44k cycles per round (0.148 ms at 300 MHz) with P_round ≈ 0.997, giving TTS99 ≈ 0.12–0.15 ms. This falls back to M4-only if timing forces a lower clock.

No external review. Single instance (WK2000_1), single target (cut ≥ 33,000).

## Amendment (5 October 2026, before any multi-engine data)

The single-engine v5.1 builds failed routed timing:

| Clock | WNS |
|---|---|
| 300 MHz | −0.269 ns |
| 275 MHz | −0.060 ns |

The paths are route-dominated: the per-flip ±d broadcast to the 2,048 field adders. v5 closed at 250 MHz (+0.001 ns).

The 4-engine v5.2 image is therefore built at **250 MHz**, with 225 MHz as fallback if 250 fails. Everything else is unchanged, and M3 is evaluated with its clock scaling (≤ 0.20 ms × 300/f).

## Amendment 2 (5 October 2026, before any multi-engine data): 2 engines of v5.3

- **Why.** The 4-engine v5.2 build at 250 MHz could not close: post-placement WNS was about −1.8 ns, on engine-internal LANE and ACCUM paths, so the run was stopped. The image under test becomes **v5.3 × 2 engines at 250 MHz**.
  - v5.3 = v5.2 plus a registered per-field partial sum and a spin-word queue. It has the same cycles and is bit-exact; its HLS estimate is 331 MHz and its RTL co-simulation passed.
- **Procedure.** Phases 1, 2a and 2b are unchanged, with `--engines 2`. All trial counts are divisible by 2 × trials-per-launch.
- **Predictions for 2 engines.** These come from v5's measured cycles plus burst table loading.
  - **M1** and **M2**: unchanged.
  - **M3′**: the best measured 2-engine Onsager TTS99 is ≤ 0.30 ms × (250/f). Predicted: O1 has P_round ≈ 0.995 and a round of about 0.28 ms, giving TTS99 ≈ 0.24 ms.
  - **M4′**: best Onsager vs best plain and vs best TEC are each ≥ 2× (point estimate). Predicted about 3.5–4×.

## Amendment 3 (5 October 2026, before any multi-engine data): 3 floorplanned engines

**Why.** The unconstrained 2-engine v5.3 build also failed at placement (WNS −1.68 ns) and was stopped.
- The v5 routed checkpoint shows one engine filling SLR0, rows 1–3, columns X1–X8 (about 21 clock regions).
- The xcv80 has three SLRs: rows 0–4, 5–7 and 8–11. Unconstrained replicated engines straddle die boundaries.

**Image under test.** v5.3 × **3 engines at 250 MHz**, each confined by a pblock to its own SLR (`scripts/pblocks_e3.xdc`). If it fails, 225 MHz is the fallback. Phases 1, 2a and 2b are unchanged, with `--engines 3`.

**Trial counts.** 1,024 is not divisible by 3, so the counts are adjusted, all divisible by 3 × trials-per-launch:
- 2b: 1,026 trials per configuration (342 rounds), offsets 0–1025.
- 2a: 3 × 384 = 1,152 trials, 64 per launch.
- Phase 1: W1 = 9 trials, W2 = 12 (per-launch 2), W3 = 3, W4 = 9.

**Predictions** (v5 measured cycles plus burst table loading):
- **M1** and **M2**: unchanged.
- **M3″**: best measured 3-engine Onsager TTS99 (Eq. 9 on rounds) ≤ 0.25 ms × (250/f). Predicted: O1 has P_round ≈ 0.9997 and t_R ≈ 0.28 ms, giving about 0.16 ms.
- **M4″**: best Onsager vs best plain and vs best TEC are each ≥ 2× (point estimate).
