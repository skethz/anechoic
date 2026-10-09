## Max-Cut steps to the 1% target (MCS99, min over S), geometric means over common solved sets

| Set (n) | SA (Neal) | SCA (STATICA) | TEC | APC-SCA | ReAIM ASA | aSB | bSB | dSB | Onsager-κT | Onsager-online |
|---|---|---|---|---|---|---|---|---|---|---|
| all ten methods (0) | empty | empty | empty | empty | empty | empty | empty | empty | empty | empty |
| all methods except TEC (0/20) and aSB (6) | 5,181 | 6,147 | – | 6,210 | 198,522 | – | 2,957 | 3,356 | 2,728 | 2,130 |
| SA, dSB, Onsager-online (C3) (20) | 2,261 | – | – | – | – | – | – | 2,432 | 4,031 | 5,150 |
| A1 common set (11) (11) | 1,096 | – | – | – | 18,152 | – | 594 | 1,136 | 1,397 | 1,681 |
| solved / 20 | 20 | 11 | 0 | 15 | 16 | 7 | 17 | 20 | 20 | 20 |

## Max-Cut paired step ratios, method / Onsager form, over instances both solve (n; Onsager fewer steps on k)

| Method | vs Onsager-κT: geo ratio (n, k) | vs Onsager-online: geo ratio (n, k) |
|---|---|---|
| SA (Neal) | 0.56× (20, 8) | 0.44× (20, 8) |
| SCA (STATICA) | 3.60× (11, 11) | 3.75× (11, 11) |
| TEC | – (0 solved by TEC) | – (0 solved by TEC) |
| APC-SCA | 1.51× (15, 12) | 1.32× (15, 9) |
| ReAIM ASA | 25.38× (16, 16) | 19.91× (16, 16) |
| aSB | 5.32× (7, 7) | 6.26× (7, 5) |
| bSB | 0.68× (17, 2) | 0.53× (17, 1) |
| dSB | 0.60× (20, 7) | 0.47× (20, 7) |
| Onsager-κT | 1.00× (20, 0) | 0.78× (20, 6) |
| Onsager-online | 1.28× (20, 14) | 1.00× (20, 0) |

C3: {"online_over_SA_all20": 2.278, "online_over_dSB_all20": 2.118, "online_over_SA_A1set": 1.534, "online_over_dSB_A1set": 1.48, "kT_over_SA_all20": 1.783, "kT_over_dSB_all20": 1.657}

## Compact Table 8b (A3.1)

| | SA (Neal) | SCA (STATICA) | TEC | APC-SCA | ReAIM ASA | aSB | bSB | dSB | Onsager-κT | Onsager-online |
|---|---|---|---|---|---|---|---|---|---|---|
| Max-Cut G1–G20: quality (S=4000) | 0.995 | 0.542 | -0.004 | 0.960 | 0.979 | 0.148 | 0.989 | 0.995 | 0.994 | 0.990 |
| Max-Cut: steps, geo-mean over own solved (solved/20) | 2,261 (20/20) | 29,960 (11/20) | – (0/20) | 11,323 (15/20) | 56,200 (16/20) | 35,719 (7/20) | 2,154 (17/20) | 2,432 (20/20) | 4,031 (20/20) | 5,150 (20/20) |
| GPP (9): quality (S=4096) | 0.979 | 0.000 | 0.000 | 0.166 | 0.697 | 0.000 | 0.635 | 0.000 | 0.026 | 0.026 |
| GPP: feasible | 1.00 | 0.00 | 0.00 | 0.26 | 0.73 | 0.00 | 0.65 | 0.00 | 0.05 | 0.04 |
| GPP: steps, own solved (solved/9) | 3,851 (5/9) | – (0/9) | – (0/9) | – (0/9) | 178,476 (5/9) | – (0/9) | 60,250 (5/9) | – (0/9) | – (0/9) | – (0/9) |
| TSP (6): quality (S=8192) | 0.528 | 0.000 | 0.000 | 0.542 | 0.140 | 0.000 | 0.942 | 0.540 | 0.476 | 0.474 |
| TSP: feasible | 1.00 | 0.00 | 0.00 | 0.98 | 0.27 | 0.00 | 1.00 | 1.00 | 1.00 | 1.00 |
| TSP: steps, own solved (solved/6) | – (0/6) | – (0/6) | – (0/6) | – (0/6) | – (0/6) | – (0/6) | 133,134 (3/6) | – (0/6) | – (0/6) | – (0/6) |

## A1 (tuned baselines) vs A3.1 (papers' settings): quality (feasibility)

| Method | MCP A1 | MCP A3.1 | GPP A1 | GPP A3.1 | TSP A1 | TSP A3.1 |
|---|---|---|---|---|---|---|
| SA (Neal) | 0.997 (1.00) | 0.995 (1.00) | 0.985 (1.00) | 0.979 (1.00) | 0.567 (1.00) | 0.528 (1.00) |
| SCA (STATICA) | 0.991 (1.00) | 0.542 (1.00) | 0.039 (0.06) | 0.000 (0.00) | 0.476 (1.00) | 0.000 (0.00) |
| TEC | 0.993 (1.00) | -0.004 (1.00) | 0.034 (0.06) | 0.000 (0.00) | 0.482 (1.00) | 0.000 (0.00) |
| APC-SCA | 0.994 (1.00) | 0.960 (1.00) | 0.616 (0.91) | 0.166 (0.26) | 0.560 (1.00) | 0.542 (0.98) |
| ReAIM ASA | 0.976 (1.00) | 0.979 (1.00) | 0.856 (0.88) | 0.697 (0.73) | 0.586 (0.92) | 0.140 (0.27) |
| aSB | 0.984 (1.00) | 0.148 (0.15) | 0.108 (0.11) | 0.000 (0.00) | 0.000 (0.00) | 0.000 (0.00) |
| bSB | 0.989 (1.00) | 0.989 (1.00) | 0.980 (0.99) | 0.635 (0.65) | 0.976 (1.00) | 0.942 (1.00) |
| dSB | 0.995 (1.00) | 0.995 (1.00) | 0.023 (0.03) | 0.000 (0.00) | 0.563 (1.00) | 0.540 (1.00) |
| Onsager-κT | 0.994 (1.00) | 0.994 (1.00) | 0.031 (0.06) | 0.026 (0.05) | 0.477 (1.00) | 0.476 (1.00) |
| Onsager-online | 0.991 (1.00) | 0.990 (1.00) | 0.017 (0.03) | 0.026 (0.04) | 0.475 (1.00) | 0.474 (1.00) |

## APC-SCA and ReAIM: final state (A3, seed entry 41) vs the papers' output rule (A3.1, seed entry 42)

| Method | Problem | final-state quality (feasible) | paper-output quality (feasible) |
|---|---|---|---|
| APC-SCA | MCP | 0.744 (1.00) | 0.960 (1.00) |
| APC-SCA | GPP | 0.000 (0.00) | 0.166 (0.26) |
| APC-SCA | TSP | 0.522 (0.99) | 0.542 (0.98) |
| ReAIM ASA | MCP | 0.978 (1.00) | 0.979 (1.00) |
| ReAIM ASA | GPP | 0.697 (0.73) | 0.697 (0.73) |
| ReAIM ASA | TSP | 0.022 (0.04) | 0.140 (0.27) |
