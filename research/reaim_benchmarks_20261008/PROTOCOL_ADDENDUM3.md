# Addendum 3: baselines at their papers' own settings, no tuning (frozen 8 October 2026, before any A3 run)

## Why

**The authors' decision** at about 13:38 CEST: "Do not tune — you should strictly follow what they describe. Tuning means we create our own methods."

- Addendum 2 (grid-extension re-tuning) was stopped and abandoned; see `PROTOCOL_ADDENDUM2_ABANDONED.md` and `results_a2/ABANDONED.md`.
- **Baselines:** the eight baselines now run at the settings their own papers describe. Where a paper prescribes a selection procedure, that procedure is followed exactly.
- **Our two forms:** Onsager-κT and Onsager-online keep their unchanged A1 selection procedure (PROTOCOL.md) and A1 selections.
- **Finals:** all ten methods get fresh held-out finals in a fresh directory, `results_a3/`.

## Unchanged

- Problems, instances, formulations, penalties and the bias handling (PROTOCOL.md).
- Budgets:
  - Max-Cut S ∈ {250, 500, 1000, 2000, 4000};
  - GPP S ∈ {256, …, 4096};
  - TSP S ∈ {512, …, 8192}.
- Finals: 256 runs.
- Measures and targets as in PROTOCOL.md.

## Settings (resolved per instance by `settings_a3.py`; written to `settings_a3.json`)

**Sources:** local copies of the papers are in `sources/papers/`, with the PDF SHA-256 values in `sources/papers/pdf_sha256.txt`. "Not covered" means the paper gives no setting for that problem. The rule used is stated in each case.

| Method | Setting used | Source in the paper | Problems not covered by the paper, and what is used |
|---|---|---|---|
| SA (Neal; ReAIM's CPU baseline [15]) | Sequential sweeps in variable order, Metropolis, random initial states. Geometric β schedule over S sweeps between Neal's default β_hot = ln2/(2·max_i Σ|bias_i|) and β_cold = ln(n_min/0.01)/(2·min_i min nonzero |bias_i|) | dwave-samplers 1.2.0 `SimulatedAnnealingSampler` (`_default_ising_beta_range`, `cpu_sa.cpp`), the implementation behind `neal` 0.6.0. Neal's default num_sweeps (1000) is replaced by the budget S | None: the defaults are problem-independent formulas |
| SCA (STATICA, yamamoto2021statica) | Synchronous update with the piecewise-linear transition probability (Eq. 7). Constant q = 4.0, exponential T from 40 to 5 | Sec. V-C, Table I, Fig. 24, Fig. 25 (K2000: S = 1560, q = 4.0, T_init = 40, T_fin = 5) | Every instance here (STATICA reports K2000 only). The K2000 values are multiplied by σ/σ_K2000, σ = sqrt(mean_i Σ_j J_ij²), σ_K2000 = sqrt(1999). **This transfer rule is ours**: it is the pre-registered G-set study's rule for K2000 settings |
| TEC (du2026tec) | Synchronous p-bit update with the Glauber probability of Eq. (4): P = 1/(1 + exp{2σ_i(t)[h_i(t) + J_v σ_i(t−1)]/T}). No SCA pinning. J_v = 30J, T from 100J to 0.1J | Eqs. (2)–(4); Fig. 3 text (K2000: "temperature range from k_BT = 100 J down to 0.1 J"; "the optimal value, J_v = 30 J"). The schedule shape is not stated; geometric is used | Every instance here (K2000 only). Same σ transfer as STATICA (ours) |
| APC-SCA (okonogi2023apc) | Algorithm 2 with the exact sigmoid. q_i(1) = λ/2, reset to λ/2 after a flip, else q_i ← max(q_i·r_q, q_limit); λ = largest eigenvalue of −J. Exponential T (Eq. 3) | Max-Cut: T 10 → 0.1 (Table 2); r_q = 0.45, q_limit = 0 (Sec. 4.2). TSP: T 0.2·N·min_{J≠0}\|J\| → 0.1·max\|J\| (Table 2); r_q = 0.8 (Sec. 4.4, Table 5); q_limit = 0 (Sec. 4.2's value; not stated for TSP) | GPP: not covered. The TSP column (the paper's setting for a weighted, penalty-constrained Ising model) is used |
| ReAIM ASA (10609617) | Algorithms 2–3, noise-free. T_init/T_final: Max-Cut 1/0.1, GPP 1/0.01, TSP 0.5/0.1. F = max (Max-Cut, GPP), F = min (TSP). q from a 20-entry FIFO of \|N\| (initialised to N, "q large initially"). Candidate k set {1, 2, 6, 16} | Table II; Sec. IV-C ("MCP and GPP favor the maximum value, whereas TSP leans towards the minimum. These are tested and selected in Step 4"); Sec. V-A (20-entry FIFO); Table I and Fig. 3 (k = ≥16, 6, 1, and 16, 6, 2) | All three problems are covered. **Not given by the paper:** the candidate k set (we use every k value the paper reports), ITER_trial/ITER_run (32/96, as in the earlier ReAIM reproduction), and the noise model (WL_on, Steps 2–3: not applicable noise-free) |
| aSB (goto2019combinatorial) | K = Δ = 1; ξ0 = 0.7Δ/(σ_SD√N), σ_SD = standard deviation of the elements of J. p(t) linear 0 → 1. Modified symplectic Euler with Δt = 0.9 and M = 2. x(0) = 0, y(0) small random | Fig. 2 caption; Methods ("we set Δt = 0.9 and M = 2 in Fig. 2", the K2000 setting); Eqs. (14)–(17) | None: the normalisation is a formula. M = 2 is the paper's K2000 value (M = 5 is its 100,000-spin value) |
| bSB, dSB (goto2021high) | a0 = 1; c0 = 0.5/(⟨J⟩√N), ⟨J⟩ = sqrt(Σ J²/(N(N−1))). a(t) linear 0 → a0. Inelastic walls. x, y ~ U(−0.1, 0.1). **Δt is selected as the best of {0.25, 0.5, 0.75, 1, 1.25}, the paper's own procedure**: 64-run pilot per Δt, highest pilot mean quality, ties to the smaller Δt, per (instance, budget) | Methods, "Parameter setting" | None (the procedure is problem-independent) |
| Onsager-κT, Onsager-online (ours) | Unchanged A1 selections: the G-set study's pre-registered selections for Max-Cut; PROTOCOL.md's 32-point grid selections for GPP and TSP | This study's PROTOCOL.md and research/gset_20261007/PROTOCOL.md | — |

**Field (bias) and dense penalty:** handled identically for every method, as in PROTOCOL.md.

## Implementation

- **`solvers_a3.py`:** the exact-logistic synchronous kernel for APC-SCA and TEC. It is verified bit for bit against a numpy reference with a dense J, a bias and the uniform GPP coupling.
- **`solvers_a2.py`:** ReAIM with F and T0. With F = max and T0 = 1 it is bit-identical to the frozen ReAIM; with F = min and T0 = 0.5 it is identical to an independent transcription.
- **Frozen `solvers.py`:** SA, plain SCA (the STATICA piecewise-linear engine), aSB, bSB and dSB.
- **Neal's default range:** checked against a dictionary transcription of dwave-samplers 1.2.0.
- **Checks:** `verify_a3.py`, run on gpu-host before any A3 run.

## Seeds and files (fresh)

- **Finals:** SeedSequence([20261008, 41, p, i, m, s]).
- **bSB/dSB Δt pilots:** SeedSequence([20261008, 40, p, i, m, s, k]).
- **Indices:** p = 0 (Max-Cut), 1 (GPP), 2 (TSP); i = G-set number or 100 + TSPLIB index; m = method index; s = budget index.
- **Disjointness:** these seeds are disjoint from every earlier seed of this study and of the G-set study.
- **Output:** `results_a3/final/<problem>/<instance>/<method>_S<S>.json`.
- **Runner:** `run_a3.py`, at most 70 worker processes on gpu-host.

## Report (Table 8b, recomputed)

- **Quality and feasibility:** mean normalized quality and feasibility at Max-Cut S = 4000, GPP S = 4096 and TSP S = 8192.
- **Steps to the 1% target:** MCS99, the minimum over S.
  - Max-Cut: the geometric mean over the instances that every method solves, with solved counts out of 20.
  - GPP and TSP: the same; if the common set is empty, per method with counts.
- **Conclusions:** which earlier conclusions change.
- **Comparison with ReAIM's own tables:**
  - ReAIM's best-of-20 cut values (its Tables IV and V);
  - its ARPD over 20 runs (Table VI), for ReAIM ASA and Neal: the median over the 12 disjoint blocks of 20 runs of our 256-run finals, at S = 4000 (Max-Cut, ReAIM 4096), 4096 (GPP) and 8192 (TSP).
- **Labelling:** everything is software-model output; the transfer rule for STATICA and TEC is ours and is labelled as such.
