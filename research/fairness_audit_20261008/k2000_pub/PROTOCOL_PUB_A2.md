# PROTOCOL_PUB Amendment A2: alignment with the COP study's readings (8 October 2026; written after the PROTOCOL_PUB finals, before this computation)

## Why

The authors asked that Figure 2 use the same settings and readings as the COP study (`research/reaim_benchmarks_20261008/PROTOCOL_ADDENDUM3.md` and `PROTOCOL_ADDENDUM3_1.md`) unless an error is found. Four readings differ from PROTOCOL_PUB:

| Item | PROTOCOL_PUB | COP study | Handling |
|---|---|---|---|
| STATICA T_init | Per budget, from STATICA's {30, 40, 50} by STATICA's TTS | 40 (the Table II long point) at every budget | COP reading adopted: the T_init = 40 finals already exist (k = 1) |
| bSB/dSB Δt | One Δt for K2000, chosen from Goto's five values by steps to solution (TTT-style, N_step optimised) | Per budget, highest pilot mean quality, ties to the smaller Δt | Both reported, from the existing pilots and finals of all five Δt. Goto 2021 states "the best value for each problem" with "N_step … optimized for TTT or TTS", which supports PROTOCOL_PUB; this is reported to the authors as a disagreement |
| Output rule, APC-SCA and ReAIM | Final state (primary); best visited as a secondary measure | The papers' own outputs: APC argmin_s H(σ(s)) including s = 1; ReAIM x_best updated after each run phase (Alg. 3) | COP reading adopted. APC: from the stored best-visited statistics (same definition). ReAIM: x_best is extracted by re-simulating the PROTOCOL_PUB ReAIM finals with identical seeds (`reaim_xbest.py`; the final cuts are asserted identical) |
| ReAIM k set | {128, 256, 512, 1024}: reproduces ReAIM's K2000 P_a (3 Oct reproduction; re-checked in V1) | {1, 2, 6, 16}: the k values reported in ReAIM's Table I/Fig. 3 | Decided by validation V1 (`reaim_cop_6400`): a k set that cannot reproduce ReAIM's own K2000 numbers is not used for K2000. The outcome goes to the authors |

TEC and aSB are taken from the COP study's VALIDATED_BASELINES.md once it exists. Until then the PROTOCOL_PUB TEC (sequential Glauber) and the A1 aSB (per-run divergence) are shown, labelled provisional.

## Output

`results_pub_A2/S<S>.json` in the format of `plot_alg.py`, and a claims table. Nothing is re-selected after seeing results beyond the fixed rules above.
