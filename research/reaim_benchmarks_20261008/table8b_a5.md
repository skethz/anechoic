## Table 8b (final, Addendum 5)

| Method | Max-Cut quality (S 4000) | Max-Cut feasible | Max-Cut solved/20 | GPP quality (S 4096) | GPP feasible | GPP solved/9 | TSP quality (S 8192) | TSP valid | TSP solved/6 |
|---|---|---|---|---|---|---|---|---|---|
| SA (Neal defaults) | 0.995 | 1.00 | 20/20 | 0.979 | 1.00 | 5/9 | 0.528 | 1.00 | 0/6 |
| STATICA, K2000 setting x σ-transfer [our rule] | 0.542 | 1.00 | 11/20 | 0.000 | 0.00 | 0/9 | 0.000 | 0.00 | 0/6 |
| TEC sequential p-bit, K2000 setting x σ-transfer [our rule] | 0.986 | 1.00 | 19/20 | 0.433 | 0.62 | 0/9 | 0.052 | 0.11 | 0/6 |
| APC-SCA (GPP: TSP-column settings [our choice]) | 0.960 | 1.00 | 15/20 | 0.166 | 0.26 | 0/9 | 0.542 | 0.98 | 0/6 |
| ReAIM ASA (Table I k) | 0.962 | 1.00 | 11/20 | 0.567 | 0.58 | 9/9 | 0.416 | 0.52 | 1/6 |
| aSB (Δt 0.9, M 2) | 0.148 (non-finite 85%) | 0.15 | 7/20 | 0.000 (non-finite 100%) | 0.00 | 0/9 | 0.000 (non-finite 100%) | 0.00 | 0/6 |
| bSB (Δt by TTT) | 0.988 | 1.00 | 17/20 | 0.640 | 0.66 | 5/9 | 0.317 | 0.38 | 1/6 |
| dSB (Δt by TTT) | 0.995 | 1.00 | 20/20 | 0.000 | 0.00 | 0/9 | 0.288 | 0.58 | 0/6 |
| Onsager-κT (ours) | 0.994 | 1.00 | 20/20 | 0.026 | 0.05 | 0/9 | 0.476 | 1.00 | 0/6 |
| Onsager-online (ours) | 0.990 | 1.00 | 20/20 | 0.026 | 0.04 | 0/9 | 0.474 | 1.00 | 0/6 |

## Max-Cut steps to the 1% target (MCS99, min over S): geometric means over common solved sets

- all reported rows (5 instances): SA (Neal defaults) 4,726; STATICA, K2000 setting x σ-transfer [our rule] 4,594; TEC sequential p-bit, K2000 setting x σ-transfer [our rule] 92,387; APC-SCA (GPP: TSP-column settings [our choice]) 5,875; ReAIM ASA (Table I k) 43,187; aSB (Δt 0.9, M 2) 23,386; bSB (Δt by TTT) 879; dSB (Δt by TTT) 2,861; Onsager-κT (ours) 2,407; Onsager-online (ours) 1,729
- all reported rows except aSB (6 instances): SA (Neal defaults) 5,181; STATICA, K2000 setting x σ-transfer [our rule] 6,147; TEC sequential p-bit, K2000 setting x σ-transfer [our rule] 86,651; APC-SCA (GPP: TSP-column settings [our choice]) 6,210; ReAIM ASA (Table I k) 74,915; bSB (Δt by TTT) 1,925; dSB (Δt by TTT) 3,170; Onsager-κT (ours) 2,728; Onsager-online (ours) 2,130
- SA, dSB, Onsager-online (C3) (20 instances): SA (Neal defaults) 2,261; dSB (Δt by TTT) 2,481; Onsager-online (ours) 5,150

## Max-Cut paired step ratios, method / Onsager form, over instances both solve (n; Onsager fewer steps on k)

| Method | / Onsager-κT | / Onsager-online |
|---|---|---|
| SA (Neal defaults) | 0.56× (20, 8) | 0.44× (20, 8) |
| STATICA, K2000 setting x σ-transfer [our rule] | 3.60× (11, 11) | 3.75× (11, 11) |
| TEC sequential p-bit, K2000 setting x σ-transfer [our rule] | 10.83× (19, 19) | 8.65× (19, 19) |
| APC-SCA (GPP: TSP-column settings [our choice]) | 1.51× (15, 12) | 1.32× (15, 9) |
| ReAIM ASA (Table I k) | 38.64× (11, 11) | 34.05× (11, 11) |
| aSB (Δt 0.9, M 2) | 5.32× (7, 7) | 6.26× (7, 5) |
| bSB (Δt by TTT) | 0.59× (17, 3) | 0.46× (17, 2) |
| dSB (Δt by TTT) | 0.62× (20, 8) | 0.48× (20, 7) |
| Onsager-κT (ours) | 1.00× (20, 0) | 0.78× (20, 6) |
| Onsager-online (ours) | 1.28× (20, 14) | 1.00× (20, 0) |

C3 (instances solved by SA, dSB and Onsager-online: 20): Onsager-online/SA = 2.278, Onsager-online/dSB = 2.076 (ratio of geometric-mean steps; > 1 means more steps for ours)

## bSB/dSB: Δt selected by pilot TTT (count of instances)

- gpp/bSB: Δt 0.5: 9
- gpp/dSB: Δt 0.25: 9
- mcp/bSB: Δt 0.25: 1, Δt 0.75: 9, Δt 1.0: 3, Δt 1.25: 7
- mcp/dSB: Δt 0.25: 4, Δt 0.5: 1, Δt 0.75: 3, Δt 1.0: 9, Δt 1.25: 3
- tsp/bSB: Δt 1.0: 2, Δt 1.25: 4
- tsp/dSB: Δt 1.25: 6

## Supplementary record (labelled; not Table 8b)

| Method | Max-Cut quality (S 4000) | Max-Cut feasible | Max-Cut solved/20 | GPP quality (S 4096) | GPP feasible | GPP solved/9 | TSP quality (S 8192) | TSP valid | TSP solved/6 |
|---|---|---|---|---|---|---|---|---|---|
| aSB Δt 0.5, M 2 [our choice] | 0.935 (non-finite 0%) | 1.00 | 11/20 | 0.000 (non-finite 100%) | 0.00 | 0/9 | 0.000 (non-finite 100%) | 0.00 | 0/6 |
| aSB stability-margin Δt [our rule] (TSP: frozen ancilla) | 0.586 (non-finite 40%) | 0.60 | 12/20 | 0.091 | 0.09 | 5/9 | 0.000 (non-finite 100%) | 0.00 | 0/6 |
| ReAIM k {1,2,6,16} (A3.1) [superseded] | 0.979 | 1.00 | 16/20 | 0.697 | 0.73 | 5/9 | 0.140 | 0.27 | 0/6 |
| bSB per-budget Δt by pilot mean quality (A3) [superseded] | 0.989 | 1.00 | 17/20 | 0.635 | 0.65 | 5/9 | 0.942 | 1.00 | 3/6 |
| dSB per-budget Δt by pilot mean quality (A3) [superseded] | 0.995 | 1.00 | 20/20 | 0.000 | 0.00 | 0/9 | 0.540 | 1.00 | 0/6 |
| bSB fixed Δt 0.25 (A5 finals) | 0.987 | 1.00 | 17/20 | 0.295 | 0.31 | 0/9 | 0.000 | 0.00 | 0/6 |
| bSB fixed Δt 0.5 (A5 finals) | 0.988 | 1.00 | 16/20 | 0.640 | 0.66 | 5/9 | 0.070 | 0.08 | 0/6 |
| bSB fixed Δt 0.75 (A5 finals) | 0.988 | 1.00 | 17/20 | 0.000 | 0.00 | 0/9 | 0.205 | 0.24 | 1/6 |
| bSB fixed Δt 1 (A5 finals) | 0.738 | 1.00 | 13/20 | 0.000 | 0.00 | 0/9 | 0.278 | 0.33 | 1/6 |
| bSB fixed Δt 1.25 (A5 finals) | 0.737 | 1.00 | 8/20 | 0.000 | 0.00 | 0/9 | 0.228 | 0.27 | 0/6 |
| dSB fixed Δt 0.25 (A5 finals) | 0.992 | 1.00 | 20/20 | 0.000 | 0.00 | 0/9 | 0.000 | 0.00 | 0/6 |
| dSB fixed Δt 0.5 (A5 finals) | 0.994 | 1.00 | 20/20 | 0.000 | 0.00 | 0/9 | 0.002 | 0.00 | 0/6 |
| dSB fixed Δt 0.75 (A5 finals) | 0.994 | 1.00 | 20/20 | 0.000 | 0.00 | 0/9 | 0.035 | 0.07 | 0/6 |
| dSB fixed Δt 1 (A5 finals) | 0.746 | 1.00 | 15/20 | 0.000 | 0.00 | 0/9 | 0.160 | 0.31 | 0/6 |
| dSB fixed Δt 1.25 (A5 finals) | 0.743 | 1.00 | 13/20 | 0.000 | 0.00 | 0/9 | 0.288 | 0.58 | 0/6 |
