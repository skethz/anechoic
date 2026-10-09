# Addendum 4: validate every baseline against its own paper; fix TEC and diagnose aSB (frozen 8 October 2026, before any A4 validation or re-run)

## Why

- **The authors' review:** "TEC failing everywhere ... and aSB blowing up ... look like re-implementation errors, not the methods' behaviour."
- **New rule:** **no baseline number may be used unless its implementation, at its paper's settings, reproduces a result that the paper itself reports.**

## 1. TEC: re-implemented as the paper's p-bit network

- **The paper.** Du et al. (arXiv:2608.21753) describe a single p-bit network with Glauber dynamics. Eq. (4) is P_flip,i(t) = 1/(1 + exp{2β σ_i(t)[h_i^s(t) + J_v σ_i(t−1)]}), with the current spatial field h_i^s = Σ_j J_ij σ_j(t). The steps are long compared with spin relaxation.
- **Implementation:** `solvers_a4.tec_seq`.
  - Sequential single-spin Glauber updates in index order, one sweep per cycle.
  - The spatial field is always current.
  - The temporal field J_v σ_i(t−1) comes from the configuration at the end of the previous cycle.
  - It is verified against a pure-Python transcription (`verify_a4.py`, 12/12 identical).
- **K2000 setting:** J_v = 30J (the optimum of Fig. 3b), k_BT from 100J to 0.1J, 3,000 cycles (the x-axis of Figs. 3a and 3c).
  - **Schedule shape:** not stated; **geometric** is pre-declared. A linear schedule from 100 would keep T above 45 (the K2000 field scale) until about cycle 1,650. That contradicts the steep early rise of the cut in Fig. 3(a) and the end of flips near cycle 2,800 in Fig. 3(c).
- **Labelled as our misreadings:**
  - A1's `tec` family (SCA engine with the TEC term);
  - A3's synchronous Glauber `tec_sig`.

  Both are reported separately under these labels. A3's TEC numbers are withdrawn from Table 8b.

## 2. aSB: diagnosis (exploratory, before this freeze: `diag_a4.py`, `diag_a4_asb.json`, 16 runs per case)

- **K2000 reproduction.** At the 2019 paper's K2000 setting (Δt = 0.9, M = 2, ξ0 = 0.7/(SD(J)·√N), p linear 0 → 1, Nstep = 186), our implementation gives a mean cut of 32,777.5 against the paper's 32,768 (Fig. 2B). The integrator matches the paper's Eqs. (14)–(17).
- **The divergence is not a coding error. It is a linear instability of the explicit kick at Δt = 0.9 on these matrices.**
  - For a restoring mode with eigenvalue λ < 0 of J, the update is stable only when roughly Δt²·(1 − p + ξ0|λ|) < 4.
  - The paper's ξ0 normalisation assumes zero-mean random couplings, whose spectral edge is ≈ 2σ√N, giving ξ0|λ| ≈ 1.4 (K2000: λ_min = −88.9).
  - Graphs with +1 weights have a Perron mode, λ_min ≈ −(mean degree): G1 −48.8, so ξ0|λ_min|Δt² = 4.1, and it diverges by step 12.
  - Planar graphs with hubs: G14 −22.4, giving 3.7, and it diverges.
  - Dense GPP penalties: λ_min ≈ −P·N.
  - Random ±1 graphs (G6) are near the edge: they are stable at S = 1000 in 16 runs but diverge in 16–100% of runs at S ≥ 1000 in A3.
  - At Δt = 0.5, G1, G14 and G6 do not diverge.
- **Table 8b reports aSB three ways:**
  - **(a)** the paper's K2000 values as-is (Δt = 0.9, M = 2; the A3 numbers);
  - **(b)** Δt = 0.5, M = 2. This is the paper's stated stability bound for the explicit method and Goto 2021's aSB time step. **OUR CHOICE.**
  - **(c)** Δt = min(0.9, 0.9·sqrt((1 + ξ0_K|λ_min,K|)/(1 + ξ0|λ_min|))), M = 2, which keeps K2000's stability margin. **OUR RULE.**

## 3. Validation targets (the papers' own reported results) and pre-declared criteria (`validate_a4.py`)

### V1. SA (Neal)

- **V1a. Implementation:** our numba SA against the real Neal, dwave-samplers 1.2.0 `SimulatedAnnealingSampler` (installed in a separate scratch venv), at Neal's defaults with 1,000 sweeps, on G1, G14, G22 and K2000, 256 reads each.
  - **Pass:** the mean cuts differ by at most 3 combined standard errors.
- **V1b. A reported result:** ReAIM's Table IV Neal column (Max-Cut, best of 20; Neal's sweep count unstated, so the default of 1,000 is used).
  - Real Neal and our SA each run 240 reads (12 blocks of 20) on the 11 instances; the median best-of-20 is compared.
  - **Pass:** within 0.5% of the paper's value on every instance.

### V2. STATICA

- **Target:** K2000, S = 1,560, q = 4.0, T 40 → 5 (Table II, Fig. 25): mean cut 33,073 and P(cut ≥ 33,000) = 0.77 (100 runs).
- **Runs:** 1,024.
- **Pass:** the mean is within 3 standard errors of a 100-run mean, **and** P lies inside the paper's 100-run 95% interval [0.675, 0.848].
- **Short point** (S = 560, T 30 → 5, paper 32,750): reported. The known +5σ gap from STATICA's RNG circuit is documented in `research/statica_reproduction_20261003`.

### V3. TEC

- **Target:** Fig. 3(b), cycles to reach 95% of the optimal cut (31,670) on K2000, read from the figure: **J_v = 0: 1,392**; **J_v = 30J: 740** (±10); ratio 1.88.
- **Runs:** 64 per J_v, 3,000 cycles.
- **Measure:** the first cycle at which the 64-run mean cut reaches 31,670 (the median first-passage cycle is also reported).
- **Pass:** each value within ±15%, **and** the ratio within ±20%.
- **Also reported:** the ordering of Fig. 3(a) (J_v = 6 faster and −6 slower than 0); our synchronous misreading at the same settings.

### V4. APC-SCA

- **V4a. Table 3, the paper's Algorithm 1** (the same machinery as APC, except the q control): G22, G30, G32, G35.
  - Settings: S = 1,000, T 10 → 0.1, 128 rounds, argmin output, q from q_init to q_final as tabulated (×λ/64).
  - Targets (average Ising energy): −6,545.5, −6,638.9, −2,726.3, −3,430.1.
  - **Pass:** within 3 standard errors of a 128-run mean.
- **V4b. Fig. 3, APC itself:** r_q = 0.45, q_limit = 0, S = 1,000, q_i reset to λ/2.
  - Targets (fraction of optimal energy, read from the figure): G22 98.3%, G30 97.5%, G32 98.05% (±0.3 points); G35 off-scale (worse than 95%).

### V5. ReAIM ASA

- **Targets:** its Tables IV (ASA best of 20, Max-Cut, 4,096 iterations), V (GPP, 4,096) and VI (TSP ARPD, 8,192).
- **Runs:** our A3.1 settings, 240 runs (12 blocks of 20).
- **No pass/fail:** the paper's candidate k set, ITER_trial/ITER_run and noise model are not given, and its numbers include ReRAM noise. Deviations are reported.

### V6. aSB

- **Target:** Fig. 2B of the 2019 paper: K2000, Δt = 0.9, M = 2, Nstep = 186, mean cut 32,768 (100 trials).
- **Runs:** 256.
- **Pass:** within 3 standard errors of a 100-run mean.

### V7. bSB and dSB

- **Target:** Fig. 2A of the 2021 paper: K2000, Δt = 1, a0 = 1, c0 = 0.5/(⟨J⟩√N); average cut over 1,000 trials, read from the figure (±20):

  | | Nstep = 100 | 1,000 | 10,000 |
  |---|---|---|---|
  | bSB | 32,860 | 33,215 | 33,210 |
  | dSB | 32,375 | 33,110 | 33,245 |

- **Runs:** 256 (Nstep = 100 and 1,000) and 256 (Nstep = 10,000; 16 chunks of 16).
- **Pass:** |ours − read value| ≤ 40 + 3 standard errors.

**Seeds:** validation SeedSequence([20261008, 60, crc32(job), chunk]); real-Neal integer seeds crc32(job). All are fresh.

## 4. Table 8b re-runs (`run_a4.py`; fresh seeds; `results_a4/`)

- **TEC (sequential):**
  - the K2000 settings transferred with the σ rule (J_v = 30r, T 100r → 0.1r, r = σ/σ_K2000; **OUR RULE**);
  - seeds SeedSequence([20261008, 51, p, i, m, s]).
- **aSB (b):** seeds [.., 52, ..].
- **aSB (c):** seeds [.., 53, ..].
- **Settings not specified by a paper:** Table 8b shows STATICA and TEC (and APC-SCA on GPP) both as "settings not specified by the paper" (no number) and as the transfer, labelled with its rule.
- **Usage rule:** a baseline whose validation fails is shown as "not validated", and its numbers are not used.

## 5. `VALIDATED_BASELINES.md`

The validated per-method settings, implementation notes and each method's reproduction result (our number against the paper's, with source), for re-running the K2000 Figure 2 comparison.
