# Ablation: is the theory-derived Onsager coefficient better than cheaper temporal-coupling rules? (frozen 5 October 2026, before any run)

## Question

Reviewers will argue that our method is "TEC with a schedule" or "APC-SCA" (see `research/literature_20261004/PRIOR_ART.md`). This study tests whether the derived, online coefficient c(t) = n_lin/(2T) beats the obvious cheaper alternatives. Each alternative gets its own tuning grid.

## Common settings

As in the J sweep:
- WK2000_1, clipped Eq. 7 with uniform thresholds, T_fin = 5, geometric schedule, final-state cut ≥ 33,000.

**Cost model.** The measured v5 V80 engine, `cycles = 0.2754·flips + 31.06·S + 753` at 250 MHz (RESULTS_HW.md, v5). The same model applies to every family; this favours the baselines, since APC needs per-spin q registers. Device TTS with E engines is `t_run·⌈⌈R99⌉/E⌉`. This is a projection for E > 1.

## Families

The decision field is h + (temporal term); the stay test is clamp((s·field + q_i)/(4T) + ½).

| Family | Temporal term or pinning | Pilot grid |
|---|---|---|
| Plain | none | From J: q ∈ {4,6,8} × T_init ∈ {12,15,20,30} × S ∈ {360,560,960,1560} (J pilot data reused, re-costed) |
| TEC-const | +J_v·s(t−1), constant | From J: q ∈ {4,8} × J_v ∈ {±4,±8,±16,30} × T_init ∈ {12,20,30} × S ∈ {360,960,1560} (reused) |
| **Onsager (ours)** | −λ_t·c(t−1)·s(t−1), c = n_lin/(2T) online | From J: q ∈ {4,6,8} × λ ∈ {0.5,0.7,0.9,1.05} × ramp × T_init × S (reused) |
| **TEC-T** (new) | −κ·T_{t−1}·s(t−1): a temperature-proportional anti-coupling, the theory's coefficient without the popcount (c ≈ 2x·T when x is flat); optional ramp to 0 over the last 30% | q ∈ {6,8} × κ ∈ {0.3,0.45,0.6,0.75} × ramp ∈ {off,on} × T_init ∈ {12,15,20} × S ∈ {560,960} (96) |
| **TEC-sched** (new) | +J_v(T)·s(t−1), J_v linear in T from J_lo at T_fin to J_hi at T_init (TEC's empirical trend: negative when cold, positive when hot) | J_lo ∈ {−8,−4,−2} × J_hi ∈ {0,4,8,16} × q ∈ {4,8} × T_init ∈ {20,30} × S ∈ {960,1560} (96) |
| **APC-SCA** (new) | Per-spin pinning q_i: after a flip, q_i ← q_reset; otherwise q_i ← max(q_i·r_q, q_lim); start q_i = q_reset | q_reset ∈ {8,16,32} × r_q ∈ {0.9,0.97} × q_lim ∈ {2,4} × T_init ∈ {15,20,30} × S ∈ {560,960} (72) |

New pilots use 256 runs per configuration, with seed root 50001 spawned per configuration.

## Selection and held-out confirmation

1. **Selection.** For each family and each E ∈ {1, 4, 16}, select the configuration minimising device TTS at the pilot's Wilson 95% lower bound on p.
2. **Held-out runs.** Every distinct selected configuration gets 1,024 held-out runs. All of them use the same seed (root 50002), so starts and threshold streams are paired across configurations.
3. **Extra held-out arm, Onsager-table.** The Onsager-selected configurations are rerun with c(t) replaced by a fixed table: the mean online c(t) from a separate 256-run pilot, seed 50003.

## Hypotheses (held-out, device TTS)

| ID | Hypothesis | Pass condition |
|---|---|---|
| A1 | Onsager beats TEC-T | Onsager TTS ≤ TEC-T TTS / 1.2 at E = 1, 4 and 16 |
| A2 | Onsager beats TEC-sched | Same rule |
| A3 | Onsager beats APC-SCA | Same rule |
| A4 | The online coefficient is at least as good as a fixed table | Online p ≥ table p for the selected configurations (descriptive; prior O1 result 0.927 vs 0.901) |

Outcomes are reported as observed. If A1 fails, the online popcount is not needed: the derived coefficient reduces to "anti-coupling ∝ T", and the paper must say so.
