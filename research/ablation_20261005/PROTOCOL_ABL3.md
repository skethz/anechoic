# ABL3: equal short schedules for TEC-T and Onsager (frozen 5 October 2026, after ABL2 and before any ABL3 run)

## Why

- After the ABL2 extension, TEC-T comes within 1.09× (E = 1) and 1.04× (E = 4) of Onsager.
- At E = 16, Onsager's selection uses S = 360, but no TEC-T grid ever offered S < 560. That comparison is therefore not equal-opportunity.

## Grids

The same settings as PROTOCOL_ABL apply. Pilot of 256 runs per configuration, seed root 50021.

| Family | Grid | Count |
|---|---|---|
| TEC-T | S ∈ {200, 280, 360} × κ ∈ {0.75, 1.0, 1.25, 1.5} × ramp ∈ {off, on} × q ∈ {6, 8} × T_init ∈ {12, 15, 20} | 144 |
| Onsager | S ∈ {200, 280} × λ ∈ {0.9, 1.05} × ramp ∈ {off, on} × q ∈ {6, 8} × T_init ∈ {12, 15} | 32 (S = 360 is already covered by J) |

## Selection and held-out runs

- **Selection.** For E ∈ {1, 4, 16}, use the Wilson-lower-bound rule over all pilots of each family:
  - TEC-T: ABL + ABL2 + ABL3.
  - Onsager: J + ABL3.
- **Held-out runs.** 1,024 per newly selected configuration, paired with seed root 50002.

## Final verdict (both reported)

A1 at E ∈ {1, 4, 16}: Onsager TTS ≤ TEC-T TTS / 1.2.
