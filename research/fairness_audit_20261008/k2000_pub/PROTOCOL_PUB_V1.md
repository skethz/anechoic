# PROTOCOL_PUB validation amendment V1: each baseline must reproduce a number its own paper reports (8 October 2026, before these runs)

## Why

On 8 October the authors added this rule: no baseline number may be used unless the implementation, at its paper's settings, reproduces a result that its own paper reports. The PROTOCOL_PUB finals had already run when the rule arrived. They are kept, and each method's use is now conditional on the checks below.

## Checks (`validate.py`; fresh seeds `SeedSequence([20261008, 601, case])`; gpu-host, at most 11 processes)

| Method | Our run | Paper's number (where) | Pass condition |
|---|---|---|---|
| SCA (STATICA) | K2000, q 4, T 40→5, S 1560, 1,024 runs | Table II, Fig. 22/25: P_a = 0.77 (100 runs), mean cut 33,073 | P_a inside our 99% Wilson interval, and mean cut within 2.6·sd/√100 of 33,073 |
| SCA (STATICA) | same, T 30→5, S 560 | P_a = 0.07, mean cut 32,750 | same rule (known risk: STATICA's chip RNG, see `research/statica_reproduction_20261003`) |
| SA (Neal) | done: `results_neal_check` (real dwave-samplers vs our kernel, S = 250 and 1000) and the PROTOCOL_PUB finals | ReAIM Table VII, Neal on K2000: P_a = 0.38 (4,610 ms) and 0.77 (5,646 ms); sweeps not stated | kernel = real Neal (z-tests already passed); P_a values compared with our 1,000- and 2,000-sweep finals (the sweep counts are our inference from the time difference) |
| ReAIM ASA | K2000, k {128, 256, 512, 1024}, F = max, T 1→0.1, ITER 32/96, x_best output, 4,096 and 6,400 iterations, 512 runs each | Table VII: P_a = 0.47 (0.15 ms) and 0.8 (0.23 ms); time ratio 1.53 vs iteration ratio 1.56 | each published P_a inside our 99% Wilson interval or within 0.05 |
| ReAIM ASA, COP study's k set | k {1, 2, 6, 16}, 6,400 iterations, 256 runs | same 0.8 | same; a failure means this k set cannot reproduce ReAIM on K2000 |
| aSB | K2000, Δt 0.9, M 2, N_step 186, 1,024 runs, per-run divergence counted as failure | Goto et al. 2019 FPGA (0.5 ms = 186 steps), quoted in STATICA Table II and ReAIM Table VII: P_a = 0.04 (100 runs) | our p inside the exact 95% interval of 4/100, [0.011, 0.099] |
| APC-SCA, shared machinery | G22, G30, G32, G35: the paper's "fine-tuned SCA" (Algorithm 1, logistic, q(s) exponential with its Table 3 pair, T 10→0.1, S 1000, argmin output), 128 rounds | APC Table 3: −6,545.5, −6,638.9, −2,726.3, −3,430.1 | within 2.6 two-sample standard errors (our sd, 128 rounds each side) |
| APC-SCA | G22, Algorithm 2, r_q 0.45, q_limit 0, S 1000, argmin output, 128 rounds | Sect. 4.2 and Figs. 3–4: lower average energy than fine-tuned SCA | our mean below −6,545.5 |
| bSB, dSB | no new run | Goto et al. 2021 gives K2000 TTT/TTS only with its Table S1 and Sect. S8 (supplement), which we could not obtain | reported as not validated quantitatively; qualitative checks from the PROTOCOL_PUB finals (dSB reaches 33,337 and bSB does not; bSB reaches 99% of best fastest) |
| TEC | none here | to be taken from the COP study's VALIDATED_BASELINES.md | — |

## Use

- A baseline that fails is flagged in every table and figure, and the cause is investigated. The investigation is either an implementation error, which is fixed and re-run under a new amendment, or a documented property of the paper's own system.
- Nothing is re-selected after seeing these results.
