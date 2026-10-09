# Machine-checked statements behind the Snowball paper (Lean 4/Mathlib and Coq)

7 October 2026. Every statement below is proved twice, independently, in `lean/` and in `coq/`.

**Scope.** These are the exact, finite statements that the paper's algorithm section and engine rely on. The large-N dynamical mean-field derivation itself is **not** formalised; the paper cites Sompolinsky–Zippelius, Eissfeller–Opper and Coolen for it. That derivation identifies the echo coefficient with the lag-one response of the expected spin.

What is formalised:
- that the response is exactly `1/(2T)` inside the linear band and `0` outside;
- that the summed response is exactly `n_lin/(2T)`;
- the unit conversion between the theory's rescaled field and the hardware's raw field.

| ID | Statement | Used in the paper |
|---|---|---|
| F1 | **Rule equivalence.** For T > 0 and s ∈ {−1, +1}, STATICA's stay probability `clamp((s·h + q)/(4T) + 1/2, 0, 1)` equals the expected-spin form `P(s' = +1) = (1 + clamp((h + q·s)/(2T), −1, 1))/2` when s = +1, and `1 − P(s' = +1)` when s = −1. Lemmas: `clamp(y/2 + 1/2, 0, 1) = (1 + clamp(y, −1, 1))/2` and `clamp(−y, −1, 1) = −clamp(y, −1, 1)`. | Eq. (sca); theory's f(g) = clip(g/2T) |
| F2 | **Linear-band response.** `f(g) = clamp(g/(2T), −1, 1)`, T > 0. (a) If \|g\| < 2T, f has derivative `1/(2T)` at g. (b) If \|g\| > 2T, f has derivative `0` at g. (c) For finitely many fields g_i with \|g_i\| ≠ 2T, the sum of derivatives equals `n_lin/(2T)`, with n_lin = #{i : \|g_i\| < 2T}. (d) `E[s'] = 2·P(s' = +1) − 1 = f(h + q·s)`. | Eq. (c): c = n_lin/2T is the exact summed response |
| F3 | **Units.** With u = h/√N, T̃ = T/√N, P_lin = n_lin/N and G = P_lin/(2T̃), the echo in raw field units is `√N·G = n_lin/(2T)`. | §4.2, raw vs rescaled units |
| F4 | **q_eff identity.** `s·(h − λc·s′) + q = s·h + (q − λc·s·s′)`. If s, s′ ∈ {±1}, then q_eff = q − λc when s′ = s (no flip) and q + λc when s′ = −s (just flipped). | §4.3, relation to TEC and APC-SCA |
| F5 | **Popcount field update.** For a finite set of flips with coupling bits b_e (J_e = +1 iff b_e) and new-sign bits a_e (σ_e = +1 iff a_e), `Σ 2·σ_e·J_e = 4·#{e : b_e = a_e} − 2·#E`, and in initialisation mode `Σ σ_e·J_e = 2·#{e : b_e = a_e} − #E`. These are exactly the RTL's `{pc, 2'b00} − (nvalid << 1)` and `{pc, 1'b0} − nvalid`. | §5.3, the multiplier-free field update |
| F6 | **Residue banking.** v6.2: x = 256k + 32m + 4q + e with m, q < 8 and e < 4 gives x mod 4 = e. v6.4: x = 256k + 8p + e with p < 32 and e < 8 gives x mod 8 = e. Every x < 2048 has exactly one such representation. Extractor e therefore only reads coupling rows x ≡ e, and single-copy residue banks partition the rows. | §5, §6.5 (ASIC), future FPGA revision |
| F7 | **Cut–energy.** For any finite edge list with weights w_e, J_e = −w_e and s ∈ {±1}^N: `cut(s) = (W − H(s))/2`, where cut = Σ w_e(1 − s_i s_j)/2, H = −Σ J_e s_i s_j and W = Σ w_e. Each edge contributes w_e iff s_i ≠ s_j. | §4.1, target cut ↔ energy |
