(* Machine-checked statements behind the Snowball paper (see ../SPEC.md). Rocq 9.3 + Stdlib.
   F1 rule equivalence, F2 linear-band response, F3 units, F4 q_eff identity,
   F5 popcount field update, F6 residue banking, F7 cut-energy identity. *)
From Stdlib Require Import Reals Lra Psatz ZArith Lia List.
Import ListNotations.
Open Scope R_scope.

(* ------------------------------------------------------------------ helpers *)
Definition clamp (x a b : R) : R := Rmax a (Rmin x b).

Lemma Rabs_le_bounds x a : Rabs x <= a -> -a <= x <= a.
Proof. unfold Rabs; destruct (Rcase_abs x); intros; lra. Qed.

Lemma Rabs_lt_bounds x a : Rabs x < a -> -a < x < a.
Proof. unfold Rabs; destruct (Rcase_abs x); intros; lra. Qed.

Lemma Rabs_gt_cases x a : a < Rabs x -> x < -a \/ a < x.
Proof. unfold Rabs; destruct (Rcase_abs x); intros; lra. Qed.

(* clamp to [-1,1] in closed form *)
Lemma clamp_low y : y <= -1 -> clamp y (-1) 1 = -1.
Proof. intro H; unfold clamp; rewrite (Rmin_left y 1) by lra; rewrite (Rmax_left (-1) y) by lra; reflexivity. Qed.
Lemma clamp_mid y : -1 <= y <= 1 -> clamp y (-1) 1 = y.
Proof. intro H; unfold clamp; rewrite (Rmin_left y 1) by lra; rewrite (Rmax_right (-1) y) by lra; reflexivity. Qed.
Lemma clamp_high y : 1 <= y -> clamp y (-1) 1 = 1.
Proof. intro H; unfold clamp; rewrite (Rmin_right y 1) by lra; rewrite (Rmax_right (-1) 1) by lra; reflexivity. Qed.

(* ------------------------------------------------------------------ F1: rule equivalence *)
Lemma clamp_affine y : clamp (y/2 + 1/2) 0 1 = (1 + clamp y (-1) 1) / 2.
Proof.
  destruct (Rle_or_lt y (-1)) as [H|H].
  - rewrite clamp_low by lra. unfold clamp.
    rewrite (Rmin_left (y/2+1/2) 1) by lra. rewrite (Rmax_left 0 (y/2+1/2)) by lra. lra.
  - destruct (Rle_or_lt y 1) as [H'|H'].
    + rewrite clamp_mid by lra. unfold clamp.
      rewrite (Rmin_left (y/2+1/2) 1) by lra. rewrite (Rmax_right 0 (y/2+1/2)) by lra. lra.
    + rewrite clamp_high by lra. unfold clamp.
      rewrite (Rmin_right (y/2+1/2) 1) by lra. rewrite (Rmax_right 0 1) by lra. lra.
Qed.

Lemma clamp_neg y : clamp (- y) (-1) 1 = - clamp y (-1) 1.
Proof.
  destruct (Rle_or_lt y (-1)) as [H|H].
  - rewrite (clamp_low y) by lra. rewrite (clamp_high (-y)) by lra. lra.
  - destruct (Rle_or_lt y 1) as [H'|H'].
    + rewrite (clamp_mid y) by lra. rewrite (clamp_mid (-y)) by lra. reflexivity.
    + rewrite (clamp_high y) by lra. rewrite (clamp_low (-y)) by lra. reflexivity.
Qed.

(* STATICA: probability that spin s keeps its value; expected-spin form: probability of +1 *)
Definition p_stay (T q h s : R) : R := clamp ((s*h + q)/(4*T) + 1/2) 0 1.
Definition f (T g : R) : R := clamp (g/(2*T)) (-1) 1.
Definition p_up (T q h s : R) : R := (1 + f T (h + q*s)) / 2.

Theorem F1_rule_equiv_plus T q h : 0 < T -> p_up T q h 1 = p_stay T q h 1.
Proof.
  intro HT. unfold p_up, p_stay, f.
  replace ((1*h + q)/(4*T) + 1/2) with (((h + q*1)/(2*T))/2 + 1/2) by (field; lra).
  rewrite clamp_affine. reflexivity.
Qed.

Theorem F1_rule_equiv_minus T q h : 0 < T -> p_up T q h (-1) = 1 - p_stay T q h (-1).
Proof.
  intro HT. unfold p_up, p_stay, f.
  replace ((-1*h + q)/(4*T) + 1/2) with ((- ((h + q*(-1))/(2*T)))/2 + 1/2) by (field; lra).
  rewrite clamp_affine, clamp_neg. lra.
Qed.

(* F2(d): expected next spin is f of the decision field g = h + q s *)
Theorem F2_expected_spin T q h s : 2 * p_up T q h s - 1 = f T (h + q*s).
Proof. unfold p_up. lra. Qed.

(* ------------------------------------------------------------------ F2: linear-band response *)
Lemma f_lin T g : 0 < T -> Rabs g <= 2*T -> f T g = g/(2*T).
Proof.
  intros HT Hg. apply Rabs_le_bounds in Hg. unfold f. apply clamp_mid.
  set (k := / (2*T)). assert (Hk : 0 < k) by (unfold k; apply Rinv_0_lt_compat; lra).
  assert (Hk2 : 2*T*k = 1) by (unfold k; field; lra).
  unfold Rdiv; fold k. split; nra.
Qed.

Lemma f_high T g : 0 < T -> 2*T <= g -> f T g = 1.
Proof.
  intros HT Hg. unfold f. apply clamp_high.
  set (k := / (2*T)). assert (Hk : 0 < k) by (unfold k; apply Rinv_0_lt_compat; lra).
  assert (Hk2 : 2*T*k = 1) by (unfold k; field; lra).
  unfold Rdiv; fold k. nra.
Qed.

Lemma f_low T g : 0 < T -> g <= -(2*T) -> f T g = -1.
Proof.
  intros HT Hg. unfold f. apply clamp_low.
  set (k := / (2*T)). assert (Hk : 0 < k) by (unfold k; apply Rinv_0_lt_compat; lra).
  assert (Hk2 : 2*T*k = 1) by (unfold k; field; lra).
  unfold Rdiv; fold k. nra.
Qed.

(* (a) inside the band the derivative is 1/(2T) *)
Theorem F2a_deriv_inside T g : 0 < T -> Rabs g < 2*T -> derivable_pt_lim (f T) g (/ (2*T)).
Proof.
  intros HT Hg eps Heps.
  assert (Hd : 0 < 2*T - Rabs g) by lra.
  exists (mkposreal _ Hd). intros dh Hdh Hsmall. simpl in Hsmall.
  assert (Hin : Rabs (g + dh) <= 2*T).
  { pose proof (Rabs_triang g dh). lra. }
  rewrite (f_lin T (g + dh)) by (assumption || lra).
  rewrite (f_lin T g) by (assumption || lra).
  replace (((g + dh)/(2*T) - g/(2*T))/dh - / (2*T)) with 0 by (field; split; lra).
  rewrite Rabs_R0. exact Heps.
Qed.

(* (b) outside the band the derivative is 0 *)
Theorem F2b_deriv_outside T g : 0 < T -> 2*T < Rabs g -> derivable_pt_lim (f T) g 0.
Proof.
  intros HT Hg eps Heps.
  assert (Hd : 0 < Rabs g - 2*T) by lra.
  exists (mkposreal _ Hd). intros dh Hdh Hsmall. simpl in Hsmall.
  apply Rabs_lt_bounds in Hsmall.
  destruct (Rabs_gt_cases g (2*T) Hg) as [Hneg|Hpos].
  - assert (Hag : Rabs g = - g) by (apply Rabs_left; lra).
    rewrite (f_low T (g + dh)) by lra. rewrite (f_low T g) by lra.
    replace ((-1 - -1)/dh - 0) with 0 by (field; exact Hdh). rewrite Rabs_R0; exact Heps.
  - assert (Hag : Rabs g = g) by (apply Rabs_right; lra).
    rewrite (f_high T (g + dh)) by lra. rewrite (f_high T g) by lra.
    replace ((1 - 1)/dh - 0) with 0 by (field; exact Hdh). rewrite Rabs_R0; exact Heps.
Qed.

(* (c) the summed response over spins equals n_lin / (2T) *)
Definition dresp (T g : R) : R := if Rlt_dec (Rabs g) (2*T) then / (2*T) else 0.
Fixpoint n_lin (T : R) (gs : list R) : nat :=
  match gs with
  | [] => 0%nat
  | g :: r => ((if Rlt_dec (Rabs g) (2*T) then 1 else 0) + n_lin T r)%nat
  end.
Fixpoint sumR (xs : list R) : R := match xs with [] => 0 | x :: r => x + sumR r end.

Theorem F2c_summed_response T gs : sumR (map (dresp T) gs) = INR (n_lin T gs) * / (2*T).
Proof.
  induction gs as [|g r IH]; simpl.
  - lra.
  - rewrite IH. unfold dresp. destruct (Rlt_dec (Rabs g) (2*T)); rewrite plus_INR; simpl; lra.
Qed.

(* dresp is the actual derivative wherever the field is not exactly on the band edge *)
Theorem F2_dresp_is_derivative T g : 0 < T -> Rabs g <> 2*T -> derivable_pt_lim (f T) g (dresp T g).
Proof.
  intros HT Hne. unfold dresp. destruct (Rlt_dec (Rabs g) (2*T)) as [Hlt|Hnlt].
  - apply F2a_deriv_inside; assumption.
  - apply F2b_deriv_outside; [assumption|]. lra.
Qed.

(* ------------------------------------------------------------------ F3: units *)
Theorem F3_units N T n : 0 < N -> 0 < T -> sqrt N * ((n / N) / (2 * (T / sqrt N))) = n / (2*T).
Proof.
  intros HN HT. assert (Hs : 0 < sqrt N) by (apply sqrt_lt_R0; exact HN).
  assert (Hss : sqrt N * sqrt N = N) by (apply sqrt_sqrt; lra).
  replace (n / N) with (n / (sqrt N * sqrt N)) by (rewrite Hss; reflexivity).
  field. split; lra.
Qed.

(* ------------------------------------------------------------------ F4: q_eff identity *)
Theorem F4_qeff s s' h q lc : s*(h - lc*s') + q = s*h + (q - lc*s*s').
Proof. ring. Qed.

Lemma spin_sq s : s = 1 \/ s = -1 -> s*s = 1.
Proof. intros [H|H]; subst; ring. Qed.

Theorem F4_qeff_cases s h q lc : s = 1 \/ s = -1 ->
  (s*(h - lc*s) + q = s*h + (q - lc)) /\ (s*(h - lc*(-s)) + q = s*h + (q + lc)).
Proof. intro Hs. pose proof (spin_sq s Hs). split; nra. Qed.

(* ------------------------------------------------------------------ F5: popcount field update (integers) *)
Open Scope Z_scope.
Definition Jv (b : bool) : Z := if b then 1 else -1.

(* entries: (valid, coupling bit b, new-sign bit a) of the flips broadcast in one cycle *)
Fixpoint sum_dec (l : list (bool * bool * bool)) : Z :=
  match l with [] => 0 | (v, b, a) :: r => (if v then 2 * Jv a * Jv b else 0) + sum_dec r end.
Fixpoint sum_init (l : list (bool * bool * bool)) : Z :=
  match l with [] => 0 | (v, b, a) :: r => (if v then Jv a * Jv b else 0) + sum_init r end.
(* RTL: P = v & a, Q = v & ~a; pc counts entries whose coupling bit selects P (b = 1) or Q (b = 0) *)
Fixpoint pc (l : list (bool * bool * bool)) : Z :=
  match l with [] => 0 | (v, b, a) :: r => (if (if b then andb v a else andb v (negb a)) then 1 else 0) + pc r end.
Fixpoint nvalid (l : list (bool * bool * bool)) : Z :=
  match l with [] => 0 | (v, _, _) :: r => (if v then 1 else 0) + nvalid r end.

Theorem F5_popcount_decision l : sum_dec l = 4 * pc l - 2 * nvalid l.
Proof.
  induction l as [|[[v b] a] r IH]; [reflexivity|].
  cbn [sum_dec pc nvalid]. destruct v, b, a; unfold Jv; cbn [andb negb]; lia.
Qed.

Theorem F5_popcount_init l : sum_init l = 2 * pc l - nvalid l.
Proof.
  induction l as [|[[v b] a] r IH]; [reflexivity|].
  cbn [sum_init pc nvalid]. destruct v, b, a; unfold Jv; cbn [andb negb]; lia.
Qed.

(* ------------------------------------------------------------------ F6: residue banking *)
Theorem F6_residue_v62 k m q e : 0 <= k -> 0 <= m < 8 -> 0 <= q < 8 -> 0 <= e < 4 ->
  (256*k + 32*m + 4*q + e) mod 4 = e.
Proof. intros. Z.div_mod_to_equations. lia. Qed.

Theorem F6_residue_v64 k p e : 0 <= k -> 0 <= p < 32 -> 0 <= e < 8 -> (256*k + 8*p + e) mod 8 = e.
Proof. intros. Z.div_mod_to_equations. lia. Qed.

(* every row x < 2048 has exactly one (k, m, q, e) representation, so the residue banks partition the rows *)
Theorem F6_repr_exists x : 0 <= x < 2048 ->
  x = 256*(x/256) + 32*((x/32) mod 8) + 4*((x/4) mod 8) + x mod 4.
Proof. intros. Z.div_mod_to_equations. lia. Qed.

Theorem F6_repr_unique k m q e k' m' q' e' :
  0 <= k -> 0 <= k' -> 0 <= m < 8 -> 0 <= m' < 8 -> 0 <= q < 8 -> 0 <= q' < 8 -> 0 <= e < 4 -> 0 <= e' < 4 ->
  256*k + 32*m + 4*q + e = 256*k' + 32*m' + 4*q' + e' -> k = k' /\ m = m' /\ q = q' /\ e = e'.
Proof. intros. lia. Qed.
Close Scope Z_scope.

(* ------------------------------------------------------------------ F7: cut-energy identity *)
(* edge list (i, j, w); J_e = -w_e; H = - sum J_e s_i s_j; cut = sum w_e (1 - s_i s_j)/2; W = sum w_e *)
Fixpoint cutv (s : nat -> R) (es : list (nat * nat * R)) : R :=
  match es with [] => 0 | (i, j, w) :: r => w * (1 - s i * s j) / 2 + cutv s r end.
Fixpoint energy (s : nat -> R) (es : list (nat * nat * R)) : R :=
  match es with [] => 0 | (i, j, w) :: r => - ((- w) * s i * s j) + energy s r end.
Fixpoint wsum (es : list (nat * nat * R)) : R :=
  match es with [] => 0 | (_, _, w) :: r => w + wsum r end.

Theorem F7_cut_energy s es : cutv s es = (wsum es - energy s es) / 2.
Proof. induction es as [|[[i j] w] r IH]; simpl; [field|]. rewrite IH. field. Qed.

Theorem F7_edge_contribution si sj : (si = 1 \/ si = -1) -> (sj = 1 \/ sj = -1) ->
  (1 - si * sj) / 2 = if Req_EM_T si sj then 0 else 1.
Proof.
  intros [Hi|Hi] [Hj|Hj]; subst; destruct (Req_EM_T _ _) as [E|E]; try lra.
Qed.
