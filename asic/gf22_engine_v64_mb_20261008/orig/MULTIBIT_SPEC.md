# Multi-bit couplings: shared specification (FPGA, GPU, ASIC)

7 October 2026. One arithmetic for all three platforms, so their results stay bit-identical to one reference and to each other.

## Couplings
- J_xy is a signed integer in {-M, ..., +M} with M = 2^(K-1) - 1. K is a build parameter; build and measure K = 2 first (ternary {-1, 0, +1}, which covers every G-set instance), and report K = 4 and K = 8 at least as synthesis or resource estimates.
- J is symmetric, J_xx = 0, and padding rows and columns (N..2047) are 0. N <= 2048.
- The existing ±1 mode (no zeros) is the special case in which every coupling has magnitude 1. A K = 2 engine given a ±1 matrix must reproduce the v6.4 results bit for bit; use this as the first regression on K2000.

## Storage: sign-magnitude bit planes
- One sign plane sg (1 iff J_xy < 0) and K-1 magnitude planes m_b (bit b of |J_xy|), each a 2048 x 2048 bit matrix in the existing row/bank layout (row x = coupling row J_x,·; bank b holds columns b + 32j).
- Store the couplings ONCE, in residue banks: extractor e (0..7) produces only flipped spins with x mod 8 = e, so each extractor port gets its own single-port bank holding rows x ≡ e (mod 8) at address x >> 3 (verify the residue property in the v6.4 RTL first). Memory per engine is K · 2048² bits: 8 Mbit at K = 2, half of today's four 1-bit copies (16 Mbit); 16 Mbit at K = 4.

## Field update
- For the up to 8 flips (x_e, σ_e) applied in a cycle, each field y needs Δh_y = Σ_e 2·σ_e·J_{x_e,y}.
- With τ_e = σ_e·(-1)^sg_{x_e,y} in {±1}: Δh_y = 2·Σ_b 2^b·(2·popcount_e(m_b ∧ [τ_e = +1]) − popcount_e(m_b)), summed over the K−1 magnitude planes.
- For K = 1, this reduces to today's identity Σ_e 2σ_eJ_e = 4·#{e: J_e = σ_e} − 2·#E (statement F5).
- No multipliers: per plane, two 8-input population counts and shifts.

## Widths
- Fields h_y: H_BITS = 1 + ceil(log2(2047·M + 1)), that is 12 bits at K = 2, 15 bits at K = 4 and 19 bits at K = 8. The decision z = s·h·2^16 + q_t ∓ corr (Q16.16) widens accordingly.
- Tables (4T_t, q_t, γ_t, constant term) keep their formats. For weighted problems, the schedules must be re-tuned, since the field scale grows with M.

## Correction coefficient (theory note)
- With general couplings, the lag-one echo coefficient for spin i is Σ_j J_ij²·f'(g_j). For ±1 couplings this is n_lin/2T; for near-uniform degree d it is about (mean J² · d / N) · n_lin / 2T.
- The scale factor is folded into the per-step table constant γ_t = λ_t·(mean_i Σ_j J_ij²)/(N·2T_{t-1}), computed on the host from J. No hardware change is needed. Onsager-κT scales κ in the same way.

## Reference model
- The golden model is `fpga/v80_sca/src/sca_ref.hpp` `run_trial(..., dense)`: the dense NP×NP int8 J path already carries integer couplings through exactly the same arithmetic (fields int32, z int64).
- Extend the host and test tooling to load integer J (dense int8 binary, or an edge list converted to dense) and to emit the plane format for each platform; leave `run_trial`'s arithmetic unchanged. Keep the cut formula cut = (sumw + Σ s·h / 2) / 2, with sumw = Σ_{i<j} w_ij and w = −J.
- Every platform must match the reference bit for bit (final spins, cut, flips, n_lin trace) on: K2000 (±1), at least two G-set instances (ternary, one random and one toroidal), and one random K = 4 matrix.

## Bias (external field), added 8 October 2026
Needed for constrained problems such as graph partitioning and TSP. Their Ising forms have linear terms after x → s = 2x − 1.
- **Model.** H(s) = −Σ_{i<j} J_ij s_i s_j − Σ_i b_i s_i, with integer biases b_i, |b_i| ≤ B_MAX (signed B_BITS-bit; build B_BITS = 16) and b_i = 0 for padding spins. The local field is h_i = Σ_{j≠i} J_ij s_j + b_i, and the flip energy keeps the form ΔE_i = 2 s_i h_i. With b = 0 everything reduces to the coupling-only engine.
- **Initialization only.** Fields start at h_i(0) = b_i + Σ_j J_ij s_j(0). The per-flip update Δh_y = Σ_e 2 σ_e J_{x_e,y} is unchanged, so the bias persists exactly through the incremental updates. No per-step datapath change and no extra cycles per step; only the load grows by the bias stream.
- **Storage.** One bias memory of 2048 × B_BITS per engine (32 kbit at 16 bits), loaded with the couplings as its own stream segment and read only while the initial fields are built.
- **Widths.** H_BITS = 1 + ceil(log2(max_y Σ_x |J_xy| + max_i |b_i| + 1)). The decision z widens accordingly. The host must check this bound before every launch and refuse instances that would overflow.
- **Score.** The device reports Σ_i s_i h_i with the biased fields, as now. The host computes the problem objective from the returned spins (energy E = −(Σ s h + Σ b s)/2) and checks the device sum. For Max-Cut (b = 0), the cut formula is unchanged.
- **Correction.** Unchanged: n_lin counts linear decisions of the biased z. The coefficient scaling of the "Correction coefficient" section still applies to J.
- **Reference.** Add a golden function with a bias argument in a new header, for example `src/sca_ref_bias.hpp` with `run_trial_bias(..., dense, bias)`. Leave `run_trial` unchanged. With b = 0 it must equal `run_trial` bit for bit.
- **Verification**, bit-exact (final spins, Σ s h, flips, n_lin trace):
  - b = 0 on K2000 (identical to the coupling-only engine and to v6.4);
  - random biases on random K-bit matrices;
  - at least one graph-partitioning and one TSP Ising instance (from `research/reaim_benchmarks_20261008` when available, otherwise small generated cases checked by brute force).
