# PROTOCOL_GSET_PUB Amendment A1: the COP study's sigma transfer of the published K2000 settings (8 October 2026; written after the unscaled runs, before this run)

## Why

Neither STATICA nor TEC gives a G-set setting. Under PROTOCOL_GSET_PUB, the unscaled K2000 values (q 4, T_fin 5, J_v +30) reach the 99% target on none of the 51 instances in the engine model. The cause is the field scale: G-set has σ = 2.0–6.9 against σ = 44.7 for K2000.

The COP study (`research/reaim_benchmarks_20261008/PROTOCOL_ADDENDUM3.md`) transfers K2000 settings to other instances by multiplying them by α = σ/σ_K2000, a rule it labels as its own. For consistency with that study, the same transfer is modelled here as a second reading.

## Configurations

- **SCA-pub-σ.** STATICA's K2000 long point transferred: q = 4α, T from 40α to 5α (clipped Eq. 7, geometric).
- **TEC-pub-σ.** The same plus TEC's J_v = 30α, in the engine's TEC adaptation: field h + J_v s(t−1) inside synchronous SCA with pinning q. TEC's own dynamics (Glauber, no q) cannot run on the engine.
- **S.** Chosen like S*₁₂ for every Table 8a rule: from {250, 500, 1000, 2000, 4000}, by the 12-engine primary TTS99 on 256-run pilots. Ties, or infinite TTS everywhere, go to the highest mean cut, then the smaller S.
- **Confirmation.** 1,024 runs.
- **Seeds.** Pilots `[20261008, 511, g, rule, S index]`; confirmations `[20261008, 512, g, rule]`.

## Reporting

- The model Table 8a against these baselines, both grid variants of the Onsager forms.
- The board pre-registration (`gset_pub_prereg.py`) includes these cohorts; it is marked as the COP-consistent reading.
- No board action.
