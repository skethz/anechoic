# J2: extension beyond the J grid edges — frozen protocol

3–4 October 2026, written before any J2 run. J's held-out winners sat on grid edges (λ = 1.05, T_init = 12, S = 360). J2 extends the grid in those directions, following theory D (λ_c increases with q). Everything else is unchanged from [PROTOCOL_J.md](PROTOCOL_J.md): same instance, target, kernel, cycle model, and selection by device TTS at the pilot's Wilson lower bound.

## Grid

| Family | Values |
|---|---|
| Onsager (ramp on) | q ∈ {8, 10, 12} × λ ∈ {1.05, 1.2, 1.4} × T_init ∈ {8, 10, 12} × S ∈ {240, 360, 560, 960} |
| Plain | q ∈ {10, 12} × T_init ∈ {12, 20, 30} × S ∈ {560, 960, 1560} |

- Seeds: pilot root 50001, held-out root 50002. Pilot runs are 256 per configuration; held-out runs are 1,024.
- **Held-out:** for each E ∈ {1, 4, 8, 16, 32}, select the best J2 Onsager and best J2 plain configuration.
- **Report:** whether J2 improves on J's held-out TTS for each family. The comparison is between separate held-out cohorts, so differences within about 1.2× are treated as unresolved.
