# Ablation: the derived Onsager coefficient vs cheaper temporal-coupling rules

5 October 2026. Protocols were frozen before the runs: [PROTOCOL_ABL.md](PROTOCOL_ABL.md), and [PROTOCOL_ABL2.md](PROTOCOL_ABL2.md), an extension of TEC-T's grid frozen before its held-out results (hashes in `protocol_frozen.sha256`).
- **Setup:** laptop CPU, WK2000_1, final-state cut ≥ 33,000.
- **Cost model:** the measured v5 V80 engine, 0.2754·flips + 31.06·S + 753 cycles at 250 MHz. It is applied to every family, which favours APC-SCA (per-spin registers) and the popcount-free rules.
- **Held-out runs:** 1,024 per selected configuration, **paired** (identical starts and threshold streams).
- **Files:** [abl.py](abl.py), [abl2.py](abl2.py), [analyze_abl.py](analyze_abl.py), `pilot_*.json`, `holdout*.json`, `abl*.log`.

## Families

- **Onsager (ours):** −λ_t·c(t−1)·s(t−1), with c = n_lin/(2T) computed online.
- **TEC, constant:** +J_v·s(t−1) (arXiv:2608.21753).
- **TEC-sched:** J_v linear in T, negative when cold and positive when hot (TEC's empirical trend).
- **TEC-T:** −κ·T·s(t−1). This is the theory's coefficient approximated without the popcount; with the end ramp it is exactly the hardware-free version of our rule.
- **APC-SCA:** per-spin adaptive pinning (IEICE 2023).
- **Plain SCA.**

## Held-out device TTS (ABL)

| Family | E = 1 | E = 4 | E = 16 |
|---|---|---|---|
| Plain SCA | p 0.665, 2.149 ms | p 0.561, 0.755 ms | p 0.381, 0.273 ms |
| TEC, constant | 0.659, 1.917 | 0.395, 0.724 | 0.395, 0.241 |
| TEC-sched | 0.625, 2.218 | 0.468, 0.718 | 0.363, 0.282 |
| APC-SCA | 0.312, 2.697 | 0.312, 0.830 | 0.312, 0.207 |
| TEC-T (grid edges; see ABL2) | 0.610, 1.244 | 0.610, 0.498 | 0.372, 0.150 |
| **Onsager** | **0.940, 0.540** | **0.760, 0.169** | **0.337, 0.098** |

| Hypothesis | E = 1 | E = 4 | E = 16 | Verdict |
|---|---|---|---|---|
| A1: Onsager vs TEC-T ≥ 1.2× | 2.31× | 2.95× | 1.54× | Pass (re-evaluated after ABL2) |
| A2: Onsager vs TEC-sched ≥ 1.2× | 4.11× | 4.25× | 2.88× | **Pass** |
| A3: Onsager vs APC-SCA ≥ 1.2× | 5.00× | 4.92× | 2.12× | **Pass** |

**A4: online coefficient vs fixed table.** Same configurations; the table is the mean online c(t) from a separate pilot. The online coefficient is better in all three:

| Configuration | Online p | Table p |
|---|---:|---:|
| T12, S960, λ1.05 | 0.940 | 0.923 |
| T15, S560, λ0.9 | 0.760 | 0.697 |
| T12, S360, λ0.9 | 0.337 | 0.281 |

## TEC-T grid extension (ABL2)

| E | Selected TEC-T (ABL ∪ ABL2) | p | t_run (ms) | TEC-T TTS | Onsager TTS | Onsager advantage |
|---|---|---|---|---|---|---|
| 1 | κ1.25 +ramp, q8, T20, S960 | 0.914 | 0.2936 | 0.587 ms | 0.540 ms | 1.09× |
| 4 | κ1.25 +ramp, q8, T20, S560 | 0.729 | 0.1753 | 0.175 ms | 0.169 ms | 1.04× |
| 16 | κ0.75 +ramp, q6, T20, S560 | 0.372 | 0.1505 | 0.150 ms | 0.098 ms | 1.54× |

- **A1 after the extension: FAIL** at E = 1 and E = 4: 1.09× and 1.04×, below the pre-registered 1.2×. The original ABL verdict (pass) is kept on record, as registered.
- **E = 16 is not equal-opportunity.** No TEC-T grid offered S < 560, while Onsager's E = 16 selection uses S = 360. ABL3 (frozen) gives both families the same short schedules.
- **What this means.**
  - The derived coefficient λ·n_lin/(2T) is about 2λx·T when the band fraction x is roughly constant.
  - A *fixed* temperature-proportional anti-coupling with the end ramp (κ = 1.25 ≈ 2λx with x ≈ 0.6) therefore captures most of the gain. The online popcount adds a few percent (consistent with A4).
  - The theory supplies the sign and the T-scaling. TEC's empirical trend has the opposite direction: its coupling becomes more positive when hot.
  - The popcount-free rule needs no hardware beyond a per-step table, which the V80 engine already has (`kconst`).

## ABL3: equal short schedules (final A1)

Frozen [PROTOCOL_ABL3.md](PROTOCOL_ABL3.md). Both families get S ∈ {200, 280, 360}. Selection is over all pilots; held-out runs are paired.

| E | Onsager (online) | TEC-T (fixed −κT + ramp) | Onsager advantage |
|---|---|---|---|
| 1 | λ1.05 + ramp, q8, T12, S960: p 0.940, **0.540 ms** | κ1.25 + ramp, q8, T20, S960: p 0.914, **0.587 ms** | 1.09× |
| 4 | λ0.9 + ramp, q6, T15, S560: p 0.760, **0.169 ms** | κ1.25 + ramp, q8, T20, S560: p 0.729, **0.175 ms** | 1.04× |
| 16 | λ0.9 + ramp, q6, T15, S280: p 0.393, **0.097 ms** | κ1.5 + ramp, q8, T15, S360: p 0.316, **0.100 ms** | 1.04× |

- **Final A1 (equal opportunity): FAIL at every E.** The online coefficient is better by only 4–9%, below the pre-registered 1.2×.
- TEC-T's E = 16 selection sits at ABL3's κ = 1.5 edge, so its margin could shrink slightly further.
- **Bottom line.** The 2–5× gains over plain SCA, TEC's recipes and APC-SCA come from a **temperature-scaled anti-echo term with an end ramp**. That is the theory's predicted form. Whether it is implemented online (popcount) or as a fixed κT schedule matters only at the 5–10% level.

## Interpretation

- **TEC's own recipe adds nothing at this target.** Its best schedules keep J_v ≤ 0, and the ferromagnetic high-T coupling it recommends for a 95%-of-optimum target does not help at cut ≥ 33,000.
- **Per-spin adaptive pinning is the weakest rule.** APC-SCA's mechanism matches ours in direction (q_eff = q − λc·s(t)s(t−1)), but its heuristic memory is far less effective than the one-step derived coefficient.
- **The online coefficient beats every fixed schedule tested, but only narrowly against the best one.**
  - A strong temperature-proportional anti-coupling with the end ramp (TEC-T, κ ≈ 1.25) comes within 4–9% at E = 1 and E = 4 (ABL2).
  - The exact mean c(t) used as a fixed table is 0.017–0.063 lower in p than online (A4). Each run's own n_lin carries information beyond the mean, because the coefficient fluctuates from run to run (DMFT report, section 8).
  - The measured echo coefficient is not proportional to T, so a fixed κT only approximates it: about 2.4·T early in an Onsager run and about 0.07·T late (DMFT report, section 7). The end ramp absorbs most of the late mismatch.
- **What the theory contributes is the form:** an anti-coupling whose strength scales with n_lin/(2T) ≈ 2x·T. Both the online and the fixed-κT implementations follow from it. TEC's constant or hot-positive couplings do not.
