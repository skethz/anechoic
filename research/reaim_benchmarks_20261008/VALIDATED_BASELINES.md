# Validated baselines: COP study, 8 October 2026

This file records the settings for re-running the K2000 Figure 2 comparison, and it backs Table 8b.

- **The authors' rule:** no baseline number may be used unless its implementation, at its paper's settings, reproduces a result that the paper itself reports.
- **The authors' change (Addendum 5.1):** every baseline is reported with its validation status, TEC included.
- **Sources:**
  - protocols: `PROTOCOL_ADDENDUM4.md` (V1–V7), `PROTOCOL_ADDENDUM5.md` (V8) and `PROTOCOL_ADDENDUM5_1.md`;
  - results: `validation_a4/`, `validation_a5/`, `validation_a4_summary.json`, `validation_a5_summary.json`, `validation_v8b.json`;
  - exploratory: `diag_a4_asb.json`, `diag_tec.json`, `diag_statica.json`;
  - fairness audit: `research/fairness_audit_20261008/k2000_pub/results_validation/`, cited as "audit".
- **Conventions:**
  - E(s) = −Σ_{i<j} J_ij s_i s_j − Σ_i b_i s_i, with local field h = Js + b.
  - Max-Cut: J = −w, and cut = (W − Σ_edges w s_i s_j)/2.
  - K2000 = `fpga/v80_snowball/data/WK2000_1.rud` (SHA-256 9ed615e5…): N = 2000, ±1 weights, SD(J) = 0.99975, ⟨J⟩ = 1, λ_max(−J) = 88.906, best cut 33,337.
  - P_a = P(cut ≥ 33,000).

## One row per method

| Method | Settings (paper's; K2000 values) | Source | Reproduction: ours against the paper's | Verdict |
|---|---|---|---|---|
| **SA (D-Wave Neal)** | dwave-samplers 1.2.0 `SimulatedAnnealingSampler` defaults: β range from `_default_ising_beta_range` (K2000: 1.734e-4 → 6.103), geometric in β, sequential Metropolis sweeps in index order, random start, final state. Sweeps = budget (default 1,000). | Neal software; ReAIM ISCA 2024 [15] Tables IV and VII (its CPU baseline) | **V1a** (our numba SA against the real library, 1,000 sweeps, 256 reads; mean cut): G1 11,603.0 / 11,602.7; G14 3,045.0 / 3,045.2; G22 13,322.9 / 13,323.9; K2000 32,960.0 / 32,958.0 (SE 11.4). β ranges are identical. **V1b** (ReAIM Table IV Neal column, best of 20, median of 12 blocks; ours / paper): 9 of 11 within 0.5%. G13 579 / 576 (+0.5%) and G19 902.5 / 893 (+1.1%); the real library gives 580 and 904, the same deviations. | **Validated**: identical to the real library. The strict V1b criterion (all 11 within 0.5%) fails on G13 and G19, where the library itself beats ReAIM's reported Neal (sweeps not stated by ReAIM). |
| **STATICA (SCA)** | Synchronous SCA. Stay probability clip((σ_i h_i + q)/(4T) + ½, 0, 1) (Eq. 7); q = 4.0 constant; geometric T from T_init to T_fin = 5 over S steps (step S at T_fin); random start; final state. Long point S = 1,560, T_init = 40. Short point S = 560, T_init = 30. | Yamamoto et al., JSSC 56(1) 2021, Table I (q = 4, experiment 3), Table II, Figs. 22 and 25 | **Long point.** Pre-registered V2 (n = 1,024): mean 33,095.2 vs 33,073 (+2.4 SE₁₀₀), P_a 0.867 vs 0.77 (paper's 100-run 95% interval [0.675, 0.848]). Exploratory, fresh seeds: n = 4,096 gives 33,091.6 and P_a 0.833; a transcription of the earlier reproduction's kernel (n = 2,048) gives 33,093.0 and 0.832. Earlier reproduction (`research/statica_reproduction_20261003`): 0.827. Audit: 0.823. **Short point:** 32,830.9 vs 32,750, P_a 0.136 vs 0.07. Known: this point is reproduced only with STATICA's RNG circuit modelled. | **Long point reproduced.** The pre-registered n = 1,024 sample missed the P_a criterion; the 4,096-run sample and two independent kernels agree at about 0.83. Short point not reproduced without the chip RNG. |
| **TEC (temporal coupling)** | Sequential p-bit Glauber updates, one sweep in index order per time step. P_flip,i = 1/(1 + exp{2βσ_i(t)[h_i(t) + J_v σ_i(t−1)]}): current spatial field, temporal field from the previous step's configuration. J_v = 30J (Fig. 3b optimum); k_BT from 100J to 0.1J, geometric (shape not stated); 3,000 cycles; random start; final state. | Du et al. 2026, arXiv:2608.21753, Eq. (4), Fig. 3 | **V3** (K2000, 64 runs, first cycle at which the mean cut reaches 31,670 = 95% of 33,337): J_v = 0: **739** (median first passage 708) vs **1,392**; J_v = 30J: **577** (551) vs **740**; ratio 1.28 vs 1.88. Fig. 3(a) ordering reproduced: J_v = 6 gives 680 and J_v = −6 gives 838. Final mean cut at 3,000 cycles: 33,107.5 (J_v = 0), 32,642.0 (J_v = 30J). Exploratory, other update orders or schedule shapes: geometric 745–752 / 575–606; linear 2,447–2,460 / 2,134–2,155. Audit: 736/577. | **Not reproduced numerically.** This reading reaches the paper's 95% target faster than reported (739/577 against 1,392/740), with a smaller J_v gain. Used in Table 8b with this note (user's decision). The synchronous reading is withdrawn and used nowhere. |
| **APC-SCA** | Algorithm 2: SCA with logistic flip probability 1/(1 + exp((σ_i h_i + q_i)/T)); q_i(1) = λ/2 with λ = λ_max(−J); q_i reset to λ/2 when spin i flips, else q_i ← max(r_q·q_i, q_limit); exponential T; output = lowest-energy visited state. **Max-Cut:** T 10 → 0.1, r_q = 0.45, q_limit = 0. **TSP:** T 0.2·N·min\|J\| → 0.1·max\|J\|, r_q = 0.8, q_limit = 0. | Okonogi et al. 2023: Algorithm 2, Table 2, Sect. 4.2 and 4.4, Fig. 3, Table 3 | **V4a**, the paper's Algorithm 1 "fine-tuned SCA" (same machinery, Table 3's exponential q pairs, S = 1,000, T 10 → 0.1, 128 rounds, mean energy): G22 −6,540.2 vs −6,545.5 (z +1.3); G30 −6,634.6 vs −6,638.9 (+1.6); G32 −2,724.9 vs −2,726.3 (+1.2); G35 −3,428.9 vs −3,430.1 (+0.6). **V4b**, APC itself (Fig. 3, r_q = 0.45, S = 1,000): G22 98.33% of optimum vs 98.3%; G30 97.50% vs 97.5%; G32 98.07% vs 98.05%; G35 83.1% (paper: below the 95% axis). Audit: Table 3 within 1.5σ; G22 98.34%. | **Validated on Max-Cut.** TSP uses the paper's TSP settings, but its TSP results cannot be reproduced (they use the authors' unpublished Ising model). GPP: no paper setting; the TSP column is used (**our choice**). |
| **ReAIM ASA** | Algorithms 2–3, noise-free. Table II: Max-Cut T 1 → 0.1; GPP 1 → 0.01; TSP 0.5 → 0.1. F = max (Max-Cut, GPP), min (TSP) (Sect. IV-C). Candidate k per Table I: Max-Cut ≥ 16 → {128, 256, 512, 1024} (capped at N); GPP {6}; TSP {1}. Not given by the paper: ITER_trial 32, ITER_run 96 (same as audit); FIFO of 20 initialised with N. Output x_best (Algorithm 3: the best run-phase end state). | Chiang et al., ISCA 2024: Algorithms 2–3, Tables I, II, IV–VII | **V8a** (K2000 Table VII, 512 runs, x_best): 4,096 iterations: P_a **0.523** (99% CI 0.467–0.580) vs **0.47**, pass; 6,400 iterations: **0.748** (0.696–0.794) vs **0.80**, 0.002 outside the pre-declared tolerance. The iteration counts are the audit's mapping of 0.15 and 0.23 ms; the 3 Oct reproduction found 0.75–0.82 over 6,100–8,200 iterations. Final-state output gives 0.381 and 0.525, so x_best matters. Audit: 0.543/0.738; k {1, 2, 6, 16}: 0.00. **V8b** (240 runs, best-of-20 medians): Table VI TSP ARPD with k = 1: 15.7, 29.7, 20.7, 36.5, 30.6, 24.9 vs 16.4, 32.1, 21.8, 44.0, 31.0, 28.5. Table V GPP with k = 6: −1.9% to +0.2%. Table IV G-set Max-Cut with k {128..800}: −3.7% to +0.2% (G19 −3.7%, G20 −3.4%, G11–G14 −1.6% to −1.8%). | **Validated on K2000 at 0.15 ms and on TSP (Table VI)**; the 0.23 ms point is a near miss (0.002). G-set Max-Cut sits up to 3.7% below Table IV with this k set, which the paper does not specify. |
| **aSB (adiabatic SB)** | K = Δ = 1; ξ0 = 0.7Δ/(SD(J)·√N) (K2000: 0.01566); p(t) linear from 0 to 1; modified explicit symplectic Euler: M = 2 sub-steps of Δt/M for the x–y oscillator, then one kick ξ0·Δt·Jx (Eqs. 14–17); **Δt = 0.9, M = 2** (Fig. 2); x(0) = 0, y(0) uniform in (−0.1, 0.1); output sign(x). | Goto, Tatsumura and Dixon, Sci. Adv. 5, eaav2372 (2019): Methods, Fig. 2 | **V6** (Fig. 2B, K2000, N_step = 186, 256 runs): mean cut **32,761.8** over 255 finite runs vs **32,768** (z −0.46); 1 run of 256 non-finite. Audit: 32,756.4, 6/1,024 non-finite, P_a 0.020 vs 0.04. | **Reproduced at the paper's N_step.** In floating point the explicit kick at Δt = 0.9 diverges in a growing fraction of runs at longer budgets and on matrices whose most negative eigenvalue is large relative to SD(J)·√N: G-set +1-weight and planar graphs, penalty matrices (`diag_a4_asb.json`). The 2019 FPGA used saturating fixed point. This is reported as a property of the setting. |
| **bSB** | a0 = 1; a(t) linear from 0 to a0; c0 = 0.5/(⟨J⟩√N), ⟨J⟩ = (ΣJ²/(N(N−1)))^½ (K2000: 0.01118); symplectic Euler with perfectly inelastic walls at \|x\| = 1; x, y uniform in (−0.1, 0.1); output sign(x). **Δt chosen per problem from {0.25, 0.5, 0.75, 1, 1.25}, with N_step, by time-to-target** (target 99% of the best known). Local fields via an ancillary spin. | Goto et al., Sci. Adv. 7, eabe7953 (2021): Methods, Figs. 2A and 3, Sect. S1 | **V7** (Fig. 2A, K2000, Δt = 1, 256 runs, mean cut vs values read from the figure): N_step 100: 32,853.0 vs ≈32,860; 1,000: 33,217.2 vs ≈33,215; 10⁴: 33,209.0 vs ≈33,210. | **Validated.** |
| **dSB** | As bSB, with the discrete (sign) coupling (Eq. 19). | same | **V7:** N_step 100: 32,366.2 vs ≈32,375; 1,000: 33,116.3 vs ≈33,110; 10⁴: 33,244.8 vs ≈33,245, best 33,337 (the paper also reaches 33,337). | **Validated.** |

## Rules that are ours (label them as such)

- **σ transfer of K2000 settings (STATICA, TEC).**
  - Each temperature- or field-scale constant is multiplied by r = σ_inst/σ_K2000, with σ = sqrt(mean_i Σ_j J_ij²) of the full Ising matrix (penalties included) and σ_K2000 = √1999.
  - STATICA: q = 4r, T 40r → 5r. TEC: J_v = 30r, T 100r → 0.1r. The schedule length is the budget S.
  - **Neither paper gives settings for G-set, GPP or TSP.**
- **APC-SCA on GPP:** the paper's TSP-column settings.
- **ReAIM:**
  - k capped at N (G-set N = 800: {128, 256, 512, 800});
  - ITER_trial/ITER_run 32/96;
  - Max-Cut k set {128, 256, 512, 1024}: the authors' reading of Table I's "≥ 16", validated on K2000.
- **Reading choices** (the paper is silent; not tuning):
  - TEC: geometric T schedule.
  - aSB: the non-finite-run rule (a non-finite run counts as a failure).

**Not ours:** local fields via an ancillary spin for aSB, bSB and dSB on TSP. This is Goto 2021, Sect. S1: J'_{i,N+1} = b_i, decode s_i = s'_i·s'_{N+1}; exactly checked in `run_a4b.py check`.

## For the K2000 Figure 2 re-run

- **Settings:** the K2000 values above apply directly (r = 1). The audit's K2000 choices agree with these:
  - STATICA's T_init from {30, 40, 50} by its figure of merit;
  - one bSB/dSB Δt for K2000 by minimum TTT over S;
  - ReAIM k {128, 256, 512, 1024}, x_best;
  - aSB Δt 0.9, M 2.
- **Our checks the audit could not run:**
  - **V7 validates bSB/dSB on K2000 against Fig. 2A**, using the main-text figure rather than Table S1.
  - **V1a shows the SA kernel equals the real Neal library on K2000.**
- **TEC:** if used, carry the validation note: faster than reported (739/577 against 1,392/740 cycles), smaller J_v gain.
- **Code:**
  - `solvers_a4.tec_seq` (TEC);
  - `solvers.engine_run` family `plain` (STATICA);
  - `solvers_a3b.engine_sig_best` (APC);
  - `solvers_a3b.reaim_best` (ReAIM);
  - `solvers.asb`, `solvers.sb` (SB);
  - `solvers.sa` with `settings_a3.instance_constants(P)['neal']` (SA).
