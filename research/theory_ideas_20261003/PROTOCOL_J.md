# J: joint (q, λ, ramp, T_init, S) sweep and TEC comparison: frozen protocol

3 October 2026. Written before any run of this sweep. The common settings follow [PROTOCOL_ABC.md](PROTOCOL_ABC.md):
- WK2000_1, final-state cut ≥ 33,000, clipped Eq. 7 with uniform thresholds, T_fin = 5, geometric schedule.
- Cycle model `0.434·flips − 0.98·S + 1989` at 300 MHz, plus 4 cycles per step for the Onsager popcount.
- Seeds: pilot root 40001 and held-out root 40002. Pilot runs are 256 per configuration; held-out runs are 1,024.

## Families (all inside parallel SCA)

| Family | Decision field | Grid |
|---|---|---|
| **Plain** | `h` | q ∈ {4, 6, 8} × T_init ∈ {12, 15, 20, 30} × S ∈ {360, 560, 960, 1560} |
| **Onsager** (ours) | `h − λ(t)·c(t)·s(t−1)`, c = n_lin/(2T); λ(t) constant or ramped linearly to 0 over the last 30% | q ∈ {4, 6, 8} × λ ∈ {0.5, 0.7, 0.9, 1.05} × ramp ∈ {off, on} × T_init × S as above |
| **TEC** (arXiv:2608.21753, constant temporal coupling) | `h + J_v·s(t−1)`, constant over the anneal | q ∈ {4, 8} × J_v ∈ {−16, −8, −4, 4, 8, 16, 30} × T_init ∈ {12, 20, 30} × S ∈ {360, 960, 1560} |

TEC's sign convention: J_v > 0 is ferromagnetic (keep the previous state). Its paper's best on K2000 was J_v = +30 J, for a 95%-of-optimum target. Our correction corresponds to J_v = −λc < 0, time-varying. TEC was originally specified for p-bit Glauber dynamics; embedding it in parallel SCA is our adaptation, so that only the temporal term differs between families.

## Selection and held-out confirmation

- For each E ∈ {1, 4, 8, 16, 32} and each family, select the configuration that minimises device TTS evaluated at the pilot's Wilson 95% lower bound on p.
- Run 1,024 held-out runs per distinct selected configuration.
- Report held-out p (Clopper–Pearson), t_run, device TTS, and the ratios Onsager/Plain and Onsager/TEC.
