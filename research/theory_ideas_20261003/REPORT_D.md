# D. Onsager-corrected SCA: echo, stability, and why λ≈½ works

3 October 2026. This report was produced by a theory sub-study and saved here verbatim in substance. Tags: **[D]** derived (linear response or exact algebra); **[C]** conjecture or heuristic; **[E]** numerics at N=400 with T and q scaled by √(N/2000). Scripts and saved results are in [D_theory/](D_theory/). The experiments are `exp1`–`exp9` and the simulator is `sca_onsager.py`.

## 1. Literature (verified by the sub-study)

**Mean-field theory of kinetic Ising dynamics**
- **Roudi & Hertz**, arXiv:1103.1044 (J. Stat. Mech. 2011). Synchronous TAP: `m_i(t+1) = tanh[h_i + Σ_j J_ij m_j(t) − m_i(t+1) Σ_j J_ij²(1−m_j(t)²)]`. It has no lag-2 term, because it conditions on the actual m(t).
- Inference-oriented work, none of which proposes a lag-2 correction for sampling or annealing:
  - Kappen & Spanjers, PRE 61, 5658 (2000).
  - Mézard & Sakellariou, arXiv:1103.3433.
  - Aurell & Mahmoudi, arXiv:1109.3399.
  - Aguilera et al., arXiv:2002.04309.

**Parallel dynamics and AMP**
- Coolen, cond-mat/0006011, eq. (98): the single-site process for parallel dynamics has a retarded self-interaction whenever the couplings are symmetric.
- Eissfeller & Opper, PRL 68, 2094 (1992).
- Gardner, Derrida & Mottishaw, J. Phys. 48, 741 (1987).
- Bolthausen, arXiv:1201.2891: AMP/TAP with the lag-1 Onsager term converges up to the AT line.

**Ising machines with related terms**
- **APC-SCA** (Okonogi et al., IPDPSW 2022; IEICE 2023): per-spin pinning that resets after a flip. This is the closest functional analogue.
- **TEC** (arXiv:2608.21753): adds `−J_v Σ σ_i(t)σ_i(t−1)`. **This is the same functional form as our correction.** Our distinguishing content is the coefficient `c = n_lin/(2T)` and the stability theory.
- PIMI (arXiv:2604.17109).
- Chen & Aihara, chao-dyn/9701021.
- Momentum annealing (Okuyama et al., PRE 2019).

**Not found** (not proof of novelty): an "Onsager/TAP-corrected SCA"; population annealing with PCA kernels.

**Population annealing (PA)**
- Theory: arXiv:1508.05647.
- PA weights need the kernel's stationary law.
  - Heat-bath PCA has one in closed form: `π ∝ Π_i cosh((h_i + q σ_i)/2T)`. This is the `π_T / P_hold` identity in this project.
  - Clipped-linear SCA and the corrected chain do not.
- DMFT schedule design for the coherent Ising machine: Zhou, Wong et al., arXiv:2603.13778. Not applied to SCA.

## 2. Derivation

1. **SCA as kinetic Ising [D].** `E[s_i(t+1) | s(t)] = φ((h_i + q s_i)/2T)`, where φ is the hard tanh. SCA is synchronous Glauber dynamics with couplings `W = J + qI` and gain `χ = 1/(2T)` inside the band. The heat-bath version is block Gibbs sampling of `exp(σᵀWτ/2T)`.
2. **Echo [D].** `E h_i(t) = h_i^cav(t) + Σ_{k≥1} K_i(k) s_i(t−k)`, with `K_i(k) = Σ_j J_ij J_ji R_j(t, t−k)`.
   - Lag 2: `K_i(1) = Σ_j J_ij² χ_j = n_lin/(2T) = c`.
   - Longer lags enter through inertia and network echoes.
   - For asymmetric J the echo vanishes.
3. **Key identity [D].** The correction's dimensionless strength is `b ≡ λ c χ = λ n_lin/(4T²) = λ x`, where `x = Σ_j J_ij² χ_j²` is the hard-tanh analogue of the Plefka/AT (AMP-convergence) parameter. **At λ=1 the update is a stochastic AMP iteration.**
4. **Linear stability [D].** Take a homogeneous in-band block with a Wigner spectrum. Then `y(t+1) = χ(J+q) y(t) − b y(t−1)`, and each eigenmode obeys `r² − a r + b = 0` with `a ∈ χq ± 2√x`.
   - **b > 1 → oscillatory growth**, `|r| = √b`, giving a period-≈4 "flip, stay" lock.
   - **Ordering** requires `2√x + χq > 1 + λx`.
     - At λ=0 this is `x ≳ ¼`.
     - At λ=1 it is `(1−√x)² < χq`. The correction cancels the ordering drive until x≈1, exactly where all bulk modes reach |r|=1, so there is **no usable window**.
     - For λ<1 the window `x ∈ (x_ord(λ), 1/λ)` is finite.
   - Equivalently, the corrected SCA is a PCA with history-dependent inertia:
     - q − 2Tb for spins that did not just flip.
     - q + 2Tb for spins that just flipped.
5. **Runaway via n_lin [C+E].** At λ=0, `n_lin ≈ 4 x₀ T²` with x₀ roughly constant (pseudogap scaling; field histograms collapse when scaled by 1/T²).
   - Empirically `x ≈ x₀(1 + λx)²`.
   - Solutions exist only for **λ ≤ λ_c = 1/(4x₀)**, and at λ_c, b = 1.
   - Flips scale as (1+b)².
6. **Is ½ principled? Mostly no.**
   - A time-symmetric split gives ½ heuristically, but the causal derivation gives 1.
   - The invariant quantity is `b = λx`. At q=4, x₀ ≈ ⅓, so λ_c ≈ 0.77 and ½ ≈ 0.65 λ_c, inside a flat optimum.
   - The benefit is dynamical, not stationary.
   - At λ=0, **61% of flips are undone at the next step**; at λ=0.5, 55%.

## 3. Numerical checks [E] (N=400, scaled schedule, 560 steps, q=4)

| λ | x_mid | max b | flips ÷ flips at λ=0 | flip-back fraction | mean e | p (≤0.4% gap) |
|---|---|---|---|---|---|---|
| 0 | 0.32 | 0 | 1.00 | 0.61 | −0.7368 | 0.15 |
| 0.5 | 0.50 (pred 0.51) | 0.28 | 1.52 | 0.55 | −0.7388 | 0.24 |
| 0.7 | 0.67 (pred 0.76) | 0.53 | 2.06 | 0.52 | −0.7388 | 0.28 |
| 0.8 | 0.85 | 0.81 | 2.6 | 0.50 | −0.7353 | 0.25 |
| 0.9 | 1.2 | 1.8 | 3.8 | 0.48 | −0.66 | ≈0.05 |
| 1.0 | 1.5 | 3.4 | 5.1 | 0.43 | −0.25 | 0 |

**Reproducing the N=2000 failure.** λ=1 reproduces the observed collapse: it never orders, and e ≈ −0.25 corresponds to a cut of about 9,700. b crosses 1 at T₂₀₀₀ ≈ 16.

**Forcing b directly.** Fixing b above 1 produces period-4 locking, with ⟨s(t)s(t−2)⟩ = −0.96 at b=2. Without the n_lin feedback there is no runaway.

**λ_c against q (predicted vs. observed)**

| Setting | Predicted λ_c = 1/(4x₀) | Observed failure between |
|---|---|---|
| q=4 | 0.77 | 0.8 and 0.9 |
| q=8 | 1.09 | 1.0 and 1.2 |
| q=12 | 1.65 | 1.2 and 1.4 (overestimate) |
| Gaussian SK, q=4 | 0.76 | 0.8 and 0.9 |

**Schedule variants** (p at 0.2% / 0.4% gap)

| Variant | p |
|---|---|
| q4, λ=0 | 0.10 / 0.19 |
| q4, λ=0.6 | 0.11 / 0.26 |
| **q4, λ=0.7, ramped to 0 over the last 30%** | **0.21 / 0.39** |
| q8, λ=0 | 0.07 / 0.13 |
| **q8, λ=0.9** | **0.20 / 0.35** |

**Negative results**
- A multi-lag memory kernel is worse, and diverges at λ=0.7.
- A fixed-b controller is worse.
- Re-tuning T or q at λ=0 does not recover the gain.

## 4. Predictions for WK2000 (q=4, unscaled T)

1. At λ=0, `n_lin/(4T²) ≈ 0.3`, roughly flat for T ∈ [8, 25]. At λ=0.5, b ≈ 0.2–0.25 and the flip ratio is ≈ (1+b)².
2. The λ=1 failure is non-ordering from the start. b crosses 1 at T ≈ 16 (about step 195 on the short schedule, 680 on the long one).
3. **λ_c ≈ 0.8–0.9.** λ=0.7 matches or beats λ=0.5. λ=0.9 fails, with onset near T ≈ 14.
4. q=8 with λ ≈ 0.9 beats q=4 with λ=0.5. q=8 with λ=0 is worse than q=4 with λ=0.
5. Ramping λ → 0 over the final 30% helps. Multi-lag corrections hurt.

**Proposed rule**
- From a λ=0 pilot, take x₀ = median `n_lin/(4T²)`.
- Set `λ = κ/(4x₀)` with κ ≈ 0.7–0.8.
- Ramp λ to 0 over the last 30% of steps.
- Either keep q=4, or raise q to about 8 with λ ≈ 0.9.

## 5. Novelty and uncertainty

**Known**
- The lag-2 retarded self-interaction (1987–2001).
- The lag-1 Onsager subtraction in AMP.
- History-dependent pinning (APC-SCA).
- The σ(t)σ(t−1) coupling term (TEC 2026).
- Decaying negative self-feedback (Chen–Aihara).

**Apparently new**
- `c = n_lin/(2T)` as the SCA echo.
- The identity `b = λx`.
- The result "λ=1 = stochastic AMP, so no ordering window".
- `λ_c = 1/(4x₀(q))` and the joint q–λ rule.

**Uncertain**
- The pseudogap self-consistency is phenomenological and fails at large q.
- All numerics are at N=400.
