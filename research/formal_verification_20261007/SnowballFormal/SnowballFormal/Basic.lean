import Mathlib

/-!
# Machine-checked statements behind the Snowball paper

See `../../SPEC.md`. F1 rule equivalence, F2 linear-band response, F3 units,
F4 q_eff identity, F5 popcount field update, F6 residue banking, F7 cut-energy identity.
-/

open Filter Topology

namespace Snowball

/-- clamp `x` to `[a, b]` -/
noncomputable def clamp (x a b : ℝ) : ℝ := max a (min x b)

lemma clamp_low {y : ℝ} (h : y ≤ -1) : clamp y (-1) 1 = -1 := by
  unfold clamp; rw [min_eq_left (by linarith), max_eq_left h]

lemma clamp_mid {y : ℝ} (h1 : -1 ≤ y) (h2 : y ≤ 1) : clamp y (-1) 1 = y := by
  unfold clamp; rw [min_eq_left h2, max_eq_right h1]

lemma clamp_high {y : ℝ} (h : 1 ≤ y) : clamp y (-1) 1 = 1 := by
  unfold clamp; rw [min_eq_right h, max_eq_right (by norm_num)]

/-! ## F1: rule equivalence -/

lemma clamp_affine (y : ℝ) : clamp (y / 2 + 1 / 2) 0 1 = (1 + clamp y (-1) 1) / 2 := by
  rcases le_total y (-1) with h | h
  · rw [clamp_low h]; unfold clamp
    rw [min_eq_left (by linarith), max_eq_left (by linarith)]; ring
  · rcases le_total y 1 with h' | h'
    · rw [clamp_mid h h']; unfold clamp
      rw [min_eq_left (by linarith), max_eq_right (by linarith)]; ring
    · rw [clamp_high h']; unfold clamp
      rw [min_eq_right (by linarith), max_eq_right (by norm_num)]; ring

lemma clamp_neg (y : ℝ) : clamp (-y) (-1) 1 = -clamp y (-1) 1 := by
  rcases le_total y (-1) with h | h
  · rw [clamp_low h, clamp_high (by linarith)]; ring
  · rcases le_total y 1 with h' | h'
    · rw [clamp_mid h h', clamp_mid (by linarith) (by linarith)]
    · rw [clamp_high h', clamp_low (by linarith)]

/-- STATICA: probability that spin `s` keeps its value -/
noncomputable def pStay (T q h s : ℝ) : ℝ := clamp ((s * h + q) / (4 * T) + 1 / 2) 0 1
/-- expected-spin form of the clamped rule -/
noncomputable def f (T g : ℝ) : ℝ := clamp (g / (2 * T)) (-1) 1
/-- probability that the next spin is `+1` -/
noncomputable def pUp (T q h s : ℝ) : ℝ := (1 + f T (h + q * s)) / 2

theorem F1_rule_equiv_plus {T : ℝ} (hT : 0 < T) (q h : ℝ) : pUp T q h 1 = pStay T q h 1 := by
  unfold pUp pStay f
  have : (1 * h + q) / (4 * T) + 1 / 2 = ((h + q * 1) / (2 * T)) / 2 + 1 / 2 := by
    field_simp; ring
  rw [this, clamp_affine]

theorem F1_rule_equiv_minus {T : ℝ} (hT : 0 < T) (q h : ℝ) :
    pUp T q h (-1) = 1 - pStay T q h (-1) := by
  unfold pUp pStay f
  have : (-1 * h + q) / (4 * T) + 1 / 2 = (-((h + q * (-1)) / (2 * T))) / 2 + 1 / 2 := by
    field_simp; ring
  rw [this, clamp_affine, clamp_neg]; ring

/-- F2(d): the expected next spin is `f` of the decision field `h + q s` -/
theorem F2_expected_spin (T q h s : ℝ) : 2 * pUp T q h s - 1 = f T (h + q * s) := by
  unfold pUp; ring

/-! ## F2: linear-band response -/

lemma f_lin {T g : ℝ} (hT : 0 < T) (hg : |g| ≤ 2 * T) : f T g = g / (2 * T) := by
  have h2 : 0 < 2 * T := by linarith
  obtain ⟨hl, hu⟩ := abs_le.mp hg
  unfold f
  apply clamp_mid
  · rw [le_div_iff₀ h2]; linarith
  · rw [div_le_iff₀ h2]; linarith

lemma f_high {T g : ℝ} (hT : 0 < T) (hg : 2 * T ≤ g) : f T g = 1 := by
  have h2 : 0 < 2 * T := by linarith
  unfold f; apply clamp_high; rw [le_div_iff₀ h2]; linarith

lemma f_low {T g : ℝ} (hT : 0 < T) (hg : g ≤ -(2 * T)) : f T g = -1 := by
  have h2 : 0 < 2 * T := by linarith
  unfold f; apply clamp_low; rw [div_le_iff₀ h2]; linarith

/-- (a) inside the linear band the response is `1/(2T)` -/
theorem F2a_deriv_inside {T g : ℝ} (hT : 0 < T) (hg : |g| < 2 * T) :
    HasDerivAt (f T) (1 / (2 * T)) g := by
  obtain ⟨hl, hu⟩ := abs_lt.mp hg
  have hlin : HasDerivAt (fun x : ℝ => x / (2 * T)) (1 / (2 * T)) g := by
    simpa using (hasDerivAt_id g).div_const (2 * T)
  apply hlin.congr_of_eventuallyEq
  filter_upwards [Ioo_mem_nhds hl hu] with x hx
  exact f_lin hT (abs_le.mpr ⟨hx.1.le, hx.2.le⟩)

/-- (b) outside the linear band the response is `0` -/
theorem F2b_deriv_outside {T g : ℝ} (hT : 0 < T) (hg : 2 * T < |g|) : HasDerivAt (f T) 0 g := by
  rcases lt_abs.mp hg with h | h
  · apply (hasDerivAt_const g (1 : ℝ)).congr_of_eventuallyEq
    filter_upwards [Ioi_mem_nhds h] with x hx
    exact f_high hT (le_of_lt hx)
  · apply (hasDerivAt_const g (-1 : ℝ)).congr_of_eventuallyEq
    filter_upwards [Iio_mem_nhds (show g < -(2 * T) by linarith)] with x hx
    exact f_low hT (le_of_lt hx)

/-- the response of one spin: `1/(2T)` if its decision is linear, `0` otherwise -/
noncomputable def dresp (T g : ℝ) : ℝ := if |g| < 2 * T then 1 / (2 * T) else 0

theorem F2_dresp_is_derivative {T g : ℝ} (hT : 0 < T) (hne : |g| ≠ 2 * T) :
    HasDerivAt (f T) (dresp T g) g := by
  unfold dresp
  split_ifs with h
  · exact F2a_deriv_inside hT h
  · exact F2b_deriv_outside hT (lt_of_le_of_ne (not_lt.mp h) (Ne.symm hne))

/-- (c) the summed response over all spins is `n_lin / (2T)` -/
theorem F2c_summed_response {ι : Type*} (s : Finset ι) (T : ℝ) (g : ι → ℝ) :
    ∑ i ∈ s, dresp T (g i) = ((s.filter (fun i => |g i| < 2 * T)).card : ℝ) / (2 * T) := by
  unfold dresp
  rw [Finset.sum_ite, Finset.sum_const_zero, add_zero, Finset.sum_const, nsmul_eq_mul]
  ring

/-! ## F3: units (DMFT rescaled field u = h/√N, T̃ = T/√N, P_lin = n_lin/N) -/

theorem F3_units {N T : ℝ} (hN : 0 < N) (hT : 0 < T) (n : ℝ) :
    Real.sqrt N * ((n / N) / (2 * (T / Real.sqrt N))) = n / (2 * T) := by
  have hs : 0 < Real.sqrt N := Real.sqrt_pos.mpr hN
  have hss : Real.sqrt N * Real.sqrt N = N := Real.mul_self_sqrt hN.le
  rw [show n / N = n / (Real.sqrt N * Real.sqrt N) by rw [hss]]
  field_simp

/-! ## F4: q_eff identity -/

theorem F4_qeff (s s' h q lc : ℝ) : s * (h - lc * s') + q = s * h + (q - lc * s * s') := by ring

theorem F4_qeff_cases {s : ℝ} (hs : s = 1 ∨ s = -1) (h q lc : ℝ) :
    s * (h - lc * s) + q = s * h + (q - lc) ∧ s * (h - lc * (-s)) + q = s * h + (q + lc) := by
  rcases hs with rfl | rfl <;> constructor <;> ring

/-! ## F5: popcount field update (one entry per broadcast flip: valid, coupling bit b, new-sign bit a) -/

def Jv (b : Bool) : ℤ := if b then 1 else -1

/-- RTL: `P = v & a`, `Q = v & ~a`; the coupling bit selects `P` (b = 1) or `Q` (b = 0) -/
def sel (e : Bool × Bool × Bool) : Bool := if e.2.1 then (e.1 && e.2.2) else (e.1 && !e.2.2)

def pc (l : List (Bool × Bool × Bool)) : ℤ := ((l.filter sel).length : ℤ)
def nvalid (l : List (Bool × Bool × Bool)) : ℤ := ((l.filter (fun e => e.1)).length : ℤ)
def sumDec (l : List (Bool × Bool × Bool)) : ℤ :=
  (l.map (fun e => if e.1 then 2 * Jv e.2.2 * Jv e.2.1 else 0)).sum
def sumInit (l : List (Bool × Bool × Bool)) : ℤ :=
  (l.map (fun e => if e.1 then Jv e.2.2 * Jv e.2.1 else 0)).sum

theorem F5_popcount_decision (l : List (Bool × Bool × Bool)) :
    sumDec l = 4 * pc l - 2 * nvalid l := by
  induction l with
  | nil => simp [sumDec, pc, nvalid]
  | cons e r ih =>
    obtain ⟨v, b, a⟩ := e
    simp only [sumDec, pc, nvalid, List.map_cons, List.sum_cons, List.filter_cons] at ih ⊢
    cases v <;> cases b <;> cases a <;> simp [sel, Jv] at ih ⊢ <;> omega

theorem F5_popcount_init (l : List (Bool × Bool × Bool)) : sumInit l = 2 * pc l - nvalid l := by
  induction l with
  | nil => simp [sumInit, pc, nvalid]
  | cons e r ih =>
    obtain ⟨v, b, a⟩ := e
    simp only [sumInit, pc, nvalid, List.map_cons, List.sum_cons, List.filter_cons] at ih ⊢
    cases v <;> cases b <;> cases a <;> simp [sel, Jv] at ih ⊢ <;> omega

/-! ## F6: residue banking -/

theorem F6_residue_v62 (k m q e : ℕ) (_hm : m < 8) (_hq : q < 8) (he : e < 4) :
    (256 * k + 32 * m + 4 * q + e) % 4 = e := by omega

theorem F6_residue_v64 (k p e : ℕ) (_hp : p < 32) (he : e < 8) :
    (256 * k + 8 * p + e) % 8 = e := by omega

theorem F6_repr_exists (x : ℕ) (hx : x < 2048) :
    x = 256 * (x / 256) + 32 * ((x / 32) % 8) + 4 * ((x / 4) % 8) + x % 4 := by omega

theorem F6_repr_unique (k m q e k' m' q' e' : ℕ) (hm : m < 8) (hm' : m' < 8) (hq : q < 8) (hq' : q' < 8)
    (he : e < 4) (he' : e' < 4) (h : 256 * k + 32 * m + 4 * q + e = 256 * k' + 32 * m' + 4 * q' + e') :
    k = k' ∧ m = m' ∧ q = q' ∧ e = e' := by omega

/-! ## F7: cut–energy identity (edge list `(i, j, w)`, `J_e = -w_e`) -/

noncomputable def cutv (s : ℕ → ℝ) (es : List (ℕ × ℕ × ℝ)) : ℝ :=
  (es.map (fun e => e.2.2 * (1 - s e.1 * s e.2.1) / 2)).sum
noncomputable def energy (s : ℕ → ℝ) (es : List (ℕ × ℕ × ℝ)) : ℝ :=
  (es.map (fun e => -((-e.2.2) * s e.1 * s e.2.1))).sum
noncomputable def wsum (es : List (ℕ × ℕ × ℝ)) : ℝ := (es.map (fun e => e.2.2)).sum

theorem F7_cut_energy (s : ℕ → ℝ) (es : List (ℕ × ℕ × ℝ)) :
    cutv s es = (wsum es - energy s es) / 2 := by
  induction es with
  | nil => simp [cutv, energy, wsum]
  | cons e r ih =>
    simp only [cutv, energy, wsum, List.map_cons, List.sum_cons] at ih ⊢
    rw [ih]; ring

theorem F7_edge_contribution {si sj : ℝ} (hi : si = 1 ∨ si = -1) (hj : sj = 1 ∨ sj = -1) :
    (1 - si * sj) / 2 = if si = sj then 0 else 1 := by
  rcases hi with rfl | rfl <;> rcases hj with rfl | rfl <;> norm_num

end Snowball
