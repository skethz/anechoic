# K2000 algorithm comparison with every baseline at its published settings (pre-registered 8 October 2026, before any run)

Amendment to `research/algorithm_compare_20261007/PROTOCOL.md` (Figure 2 of the paper). That protocol and its results are unchanged.

## Why

The frozen comparison tuned every baseline on a grid chosen by us (6–36 points per budget). Most baseline selections sat on a grid edge (plain SCA at T0 = 40 in 5/5 budgets, TEC at T0 = 30 in 5/5, APC-SCA at q_reset = 8 and T0 = 30 in 5/5, aSB at Δt = 1.25 in 5/5). On 8 October 2026 the authors ruled: *"Do not tune — you should strictly follow what they describe. Tuning means we create our own methods."* Baselines therefore take the settings their papers describe: parameter values, schedules, scaling recipes and options. Where a paper prescribes its own search procedure, that procedure is followed. Where a paper gives no setting for K2000, its general recipe or default is used and recorded below. The two forms of the correction (ours) keep the procedure the paper describes for them (the frozen per-budget tuning).

## Settings (no value below is chosen by us from data)

Sources are stored in `../sources/` with SHA-256 in `PROTOCOL_PUB.sha256`.

| Method | Setting used | Source and quoted basis |
|---|---|---|
| **SA (Neal)** | dwave-samplers 1.5.0 `SimulatedAnnealingSampler` defaults: β range = `_default_ising_beta_range(h, J)` = [ln 2/(2·1999), ln(2000/0.01)/2] = [1.7337e-4, 6.1030] on K2000 (T from 5767.9 to 0.16385); geometric β schedule with one β per sweep (`num_sweeps_per_beta` = 1); Metropolis; sequential variable order (`randomize_order` = False); random initial state; `num_sweeps` = S | Package source `dwave/samplers/sa/sampler.py` (hot β "so that all spins are able to flip with probability at least 50%"; cold β so that "the probability to excite any spin must be bounded by 0.01 on the final sweep"); the paper's Table 3 quotes Neal. Implemented with the frozen kernel `methods._sa_kernel` (sequential Metropolis sweeps in index order). Validated against the real package (below). |
| **SCA (STATICA)** | Clipped Eq. (7), q = 4, T_fin = 5, geometric T, random start, final state; T_init ∈ {30, 40, 50}, selected per budget by STATICA's figure of merit, TTS = t·ln 0.01/ln(1 − p), with t from STATICA's time model (cycles = 0.434·flips − 0.98·S + 1989 at 300 MHz, fitted to the published Fig. 22 times) | Yamamoto et al., JSSC 56(1) 2021, as read in `research/statica_reproduction_20261003/PROTOCOL.md`: "q is fixed at 4.0 (experiment 3, Table I)", T_fin = 5, the Experiment-3 grid T_init ∈ {50, 40, 30}; Table II points (S = 560, T_init = 30) and (S = 1560, T_init = 40). Our budgets are not STATICA's, so STATICA's own grid and figure of merit choose T_init on pilot runs. Finals are run for all three T_init values and reported. |
| **TEC (published)** | Simulated annealing of a p-bit network with Glauber updates and the temporal field: flip probability 1/(1 + exp{2β σ_i(t)[h_i + J_v σ_i(t−1)]}); J_v = +30 J; k_BT from 100 J to 0.1 J | Du et al., arXiv:2608.21753, Eq. (4) and Fig. 3: "benchmark results for a 2000-node instance, obtained using simulated annealing with TEC … over a temperature range from k_BT = 100 J down to 0.1 J"; "The optimal value, J_v = 30 J". **Declared (not stated in the paper):** one annealing cycle = one sequential sweep in index order; σ(t−1) = the configuration at the end of the previous cycle (the random start for the first cycle); geometric temperature schedule; one uniform per update attempt. One step = one cycle (N sequential updates), as for SA. |
| **TEC on STATICA** (supplementary) | STATICA's SCA as above plus TEC's published J_v = +30 in the field (h + J_v s(t−1)); T_init chosen as for STATICA | The paper presents TEC as a correction of synchronous SCA; this row combines the two published settings. Reported separately; not a Figure 2 curve. |
| **APC-SCA (published)** | Algorithm 2: flip probability sigmoid(−(h_i σ_i + q_i)/T); q_i(1) = λ/2 with λ = largest eigenvalue of −J (λ = 88.9058 on K2000); q_i ← λ/2 after a flip, else q_i ← max(q_i·r_q, q_limit); r_q = 0.45, q_limit = 0; exponential T from T_init = 10 to T_final = 0.1 | Okonogi et al., IEICE Trans. E106-D(12):1969, 2023: Algorithm 2; Table 2 ("Common parameter settings", max-cut: T_init 10, T_final 0.1); Sect. 4.2 ("We assumed the lower limit of pinning parameter q_limit = 0"; "we solved G22 by APC-SCA with r_q = 0.45 and q_limit = 0"; "the range from 0.4 to 0.9 looks suitable"). The paper's max-cut instances have |J| = 1, as K2000; no K2000-specific setting exists. The paper's output is the best visited state; the primary estimator here is the final state for every method (see Estimators). |
| **ReAIM ASA** | Algorithms 2–3, noise-free; Max-Cut T 1 → 0.1, F = max. **Declared (not reported for K2000):** k set {128, 256, 512, 1024}, ITER_trial 32, ITER_run 96, FIFO initialised with N | Chiang et al., ISCA 2024, Table II and Sect. IV-B as read in `research/reaim_reproduction_20261003/REPORT.md`. The undeclared knobs take the values with which that reproduction recovers ReAIM's published K2000 success probabilities (0.49 vs 0.47 at about 4,100 iterations; 0.75–0.82 vs 0.80 at 6,100–8,200). Implementation: frozen `methods.reaim_traj`. |
| **aSB** | K = Δ = 1; ξ0 = 0.7Δ/(σ√N), σ = 1 (SD of J on K2000); p(t) linear from 0 to 1; modified explicit symplectic Euler with Δt = 0.9 and M = 2; x(0) = 0, y(0) uniform in (−0.1, 0.1); spins = signs of the final x | Goto, Tatsumura and Dixon, Sci. Adv. 5, eaav2372 (2019), Methods and Fig. 2 caption (Europe PMC full text PMC6474767): "Constants in SB are set as K = Δ = 1 and ξ0 = 0.7Δ/(σ√N)"; "p(t) is linearly increased from 0 to 1"; "we set Δt = 0.9 and M = 2 in Fig. 2" (the K2000 FPGA experiment). Implementation: frozen `methods.asb_traj` with M = 2. |
| **bSB, dSB** | a0 = 1; a(t) linear from 0 to a0; c0 = 0.5/(⟨J⟩√N), ⟨J⟩ = (ΣJ²/(N(N−1)))^½ = 1; symplectic Euler with inelastic walls at ±1; x, y uniform in (−0.1, 0.1); Δt = the best of {0.25, 0.5, 0.75, 1, 1.25} for the problem, chosen together with N_step for the target | Goto et al., Sci. Adv. 7, eabe7953 (2021), Methods "Parameter setting" (Europe PMC PMC11323291): "we set the constants as a0 = 1 and c0 = 0.5/(⟨J⟩√N)"; "the setting of Δt is more sensitive to performance. We therefore selected the best value among five values: 0.25, 0.5, 0.75, 1, and 1.25"; "The number of time steps, N_step, is also optimized for TTT or TTS separately". Here one Δt is chosen for K2000 by the minimum over S of the steps to solution on pilot runs, and used at every S. Finals are run for all five Δt values and reported. Implementation: frozen `methods.sb_traj` with ξ = 1. |
| **Onsager-κT, Onsager-online (ours)** | Unchanged: the frozen per-budget selections and 256-run finals of `research/algorithm_compare_20261007/results/S*.json` | Paper Sect. 4.4 procedure. Additionally re-run at the same selections on fresh seeds (256 runs) only to obtain the best-visited secondary estimator and as a reproducibility check. |

## Procedure

- **Instance and target.** WK2000_1 (SHA-256 of the .rud file checked by the loader); success = cut ≥ 33,000 (H ≤ −67,040).
- **Budgets.** S ∈ {250, 500, 1000, 2000, 4000}, the frozen budgets. One step is as in the frozen protocol: an SA or TEC sweep (N sequential updates), an SCA/APC parallel update, an SB integration step, a ReAIM iteration.
- **Pilots (only where the paper prescribes a choice among its own values).** STATICA and TEC on STATICA: the three T_init values at every S; bSB and dSB: the five Δt values at every S. 256 runs each, seed `SeedSequence([20261008, 301, method code, S index, value index])`.
- **Selection rules (pilot data only, fixed now).**
  - STATICA (and TEC on STATICA), per budget: minimum TTS as defined above; if p = 0 for all three, the highest mean cut; remaining ties to the lower T_init.
  - bSB and dSB, per method: the Δt minimising min over S of S·ln 0.01/ln(1 − p); ties to the higher pilot mean cut at that S, then to the smaller Δt.
- **Finals (held out).** Every configuration of every baseline at every budget: 256 runs, seed `SeedSequence([20261008, 302, method code, S index, value index])`. The pilot-selected configuration is the reported one; the others are reported as sensitivity.
- **Ours.** Primary values: the frozen finals. Re-run for the secondary estimator: 256 runs, seed `SeedSequence([20261008, 303, method code, S index])`.
- **Method codes.** SA (Neal) 1, SCA (STATICA) 2, TEC (published) 3, TEC on STATICA 4, APC-SCA (published) 5, ReAIM ASA 6, aSB 7, bSB 8, dSB 9, Onsager-κT 10, Onsager-online 11. No seed of the frozen protocol (integer seeds 70001–994000) is reused.
- **Estimators.**
  - **Primary (unchanged from the frozen protocol, all methods):** the final state. Mean energy trajectory with 10th/90th percentiles, final mean cut ± sd, P(cut ≥ 33,000) with Wilson 95% interval, P(cut ≥ 33,200), steps to solution MCS99 = S·ln 0.01/ln(1 − p), minimised over S.
  - **Secondary (all methods):** the best state visited along the recorded trajectory (APC-SCA and ReAIM define their output as the best visited state).
- **Validation of the Neal implementation.** The real `dwave.samplers.SimulatedAnnealingSampler` (defaults, `num_sweeps` = S) on WK2000_1 at S = 250 and S = 1000, 256 reads each, seeds 20261008 + read index. Pass if both p and the mean cut agree with the kernel's finals within 2.6 standard errors (two-sample, 99%). A failure is reported; the kernel results are then flagged.
- **Execution.** gpu-host, at most 48 single-threaded worker processes, scratch only (`/scratch/USER/anechoic_fairness_20261008`). Sequential methods run in chunks of 32 runs whose inputs are drawn once per 256-run job, so chunking does not change results (checked in `pub_methods.selftest`).

## Part H: Table 2 baselines at published settings (model only, no board action)

The V80 measured the "original" plain-SCA and TEC schedules P1–P4 and T1–T3. These come from the 3 October joint sweep J (our tuning; `fpga/v80_sca/PROTOCOL_HW.md`), not from the STATICA or TEC papers. STATICA's published K2000 points (q = 4, T_init = 40, S = 1,560 and q = 4, T_init = 30, S = 560) and TEC's published J_v = +30 were never measured on the board.

`hw_model.py` predicts their 12-engine TTS99 with the paper's cycle model (cycles = 0.1424·flips + 18.09·S + 886 at 250 MHz; t_round = E[max of 12 trial times] from the empirical per-trial distribution; primary = t_round if P_round ≥ 0.995, else t_round·ln 0.01/ln(1 − P_round); secondary = t_round·ln 0.01/(12·ln(1 − p))). It uses 2,052 model trials per schedule, seed `SeedSequence([20261008, 401, index])`.

The same model is run on the measured schedules P1–P4, T1–T3, B1, B2, X5 (paper O1) and X4 (paper O5) to validate it against the board. Predictions for unmeasured schedules are reported with the validation error. No board run is made or requested by this protocol.

## Reporting

- The Figure 2 data as `results_pub/S<S>.json` in the format read by `plot_alg.py`. The baselines are at their published settings and the two forms of ours are the frozen finals. A regenerated figure goes to a new file name; the paper's `fpga27/images/alg_compare.pdf` is not overwritten.
- The claims checked:
  - "the two forms of the correction reach the lowest mean energy at every budget" among discrete methods;
  - the S = 500 success rates;
  - the steps to solution (Onsager-κT 830, Onsager-online 1,358, SA 1,661, SCA 4,095, TEC/APC-SCA 5,778–5,919, ASA 26,278, bSB 307);
  - "most step-efficient of the evaluated discrete annealing rules".
- Every number is reported as measured, including failures of our methods to lead.

## What this does not show

- Settings published for other targets or instances are used where the paper gives no K2000 setting; such a baseline may be far from its best possible performance on K2000, which is the authors' chosen trade-off (no tuning by us).
- Steps are not time.
