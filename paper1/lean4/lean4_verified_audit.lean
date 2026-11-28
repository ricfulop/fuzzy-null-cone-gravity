/-
  Mechanically Verified Lean 4.36.0 Formal Audit of arXiv:2611.34924v1
  ("An algebraic quantization of gravity from the null cone", V. Chandrasekaran)

  Zero external dependencies, zero `sorry`, verified directly by the Lean 4 kernel.
  Run:  ~/.elan/bin/lean lean4_verified_audit.lean

  Formalizes and proves 6 core theorems covering Tests 1, 2, and 5:
    1. `uniform_lower_bound_finite_trace_implies_dim_bound`
    2. `uniform_lower_bound_infinite_trace_unbounded`
    3. `ads_544_unconditional_bounded_iff_weights_bounded`
    4. `cycle_partition_left_inv`, `cycle_partition_right_inv`,
       `section52_counterexample_product_diverges_pos`,
       `section52_counterexample_cycle_diverges_neg`
    5. `clock_continuum_incompatibility_triad`
    6. `sec52_does_not_entail_finite_dim`, `circular_loop_admits_all_false_model`,
       `appendix_E2_valid_conditionally`
-/

namespace Audit2611_34924

/-! ## Helper lemmas: prefix sums, list sums, and `Nodup` bounds -/

def prefixSum (f : Nat → Int) : Nat → Int
  | 0     => 0
  | n + 1 => prefixSum f n + f n

theorem prefixSum_nonneg (f : Nat → Int) (hf : ∀ i, 0 ≤ f i) : ∀ n : Nat, 0 ≤ prefixSum f n
  | 0     => Int.le_refl 0
  | n + 1 => by
    unfold prefixSum
    have ih := prefixSum_nonneg f hf n
    have hn := hf n
    omega

theorem natCast_succ_mul (n : Nat) (c : Int) :
    ((n + 1 : Nat) : Int) * c = (n : Int) * c + c := by
  have h1 : ((n + 1 : Nat) : Int) = (n : Int) + 1 := by omega
  rw [h1, Int.add_mul, Int.one_mul]

theorem prefixSum_ge_mul (f : Nat → Int) (c : Int) (hf : ∀ i, c ≤ f i) :
    ∀ n : Nat, (n : Int) * c ≤ prefixSum f n
  | 0     => by simp [prefixSum]
  | n + 1 => by
    unfold prefixSum
    have ih := prefixSum_ge_mul f c hf n
    have hn := hf n
    rw [natCast_succ_mul]
    omega

def listSum : List Int → Int
  | []      => 0
  | x :: xs => x + listSum xs

theorem listSum_append (l₁ l₂ : List Int) :
    listSum (l₁ ++ l₂) = listSum l₁ + listSum l₂ := by
  induction l₁ with
  | nil => simp [listSum]
  | cons x xs ih =>
    simp [listSum, ih]
    omega

theorem listSum_map_const_mul (xs : List Nat) (f : Nat → Int) (k : Int) :
    listSum (xs.map (fun i => f i * k)) = listSum (xs.map f) * k := by
  induction xs with
  | nil => simp [listSum]
  | cons x xs ih =>
    simp [listSum, ih, Int.add_mul]

theorem listSum_range_eq_prefixSum (f : Nat → Int) :
    ∀ n : Nat, listSum ((List.range n).map f) = prefixSum f n
  | 0     => by simp [listSum, prefixSum]
  | n + 1 => by
    rw [List.range_succ, List.map_append, listSum_append, listSum_range_eq_prefixSum f n]
    simp [listSum, prefixSum]

theorem exists_nat_gt_int (B : Int) : ∃ N : Nat, 0 < N ∧ B < (N : Int) := by
  let Bpos : Int := max B 0
  let N : Nat := Bpos.toNat + 1
  refine ⟨N, by omega, ?_⟩
  have h1 : B ≤ Bpos := Int.le_max_left B 0
  have h2 : 0 ≤ Bpos := Int.le_max_right B 0
  have h3 : (Bpos.toNat : Int) = Bpos := Int.toNat_of_nonneg h2
  have h4 : (N : Int) = Bpos + 1 := by
    dsimp [N]
    omega
  omega

/-! ## Theorems 1 & 2: Uniform lower bound vs. trace-class finiteness (Section 5.2, eq. 5.30c) -/

/-- Theorem 1: If `c ≤ a i` on a `d`-dimensional basis and `∑_{i < d} a i ≤ T`, then `d * c ≤ T`
    directly—without any graviton two-point functions, boosts, or cycle permutations. -/
theorem uniform_lower_bound_finite_trace_implies_dim_bound
    (a : Nat → Int) (c T : Int) (d : Nat)
    (hbound : ∀ i, c ≤ a i) (htrace : prefixSum a d ≤ T) :
    (d : Int) * c ≤ T :=
  Int.le_trans (prefixSum_ge_mul a c hbound d) htrace

/-- Theorem 2: On an infinite basis `ℕ`, a uniform lower bound `0 < c ≤ a i` forces the partial
    trace `∑_{i < N} a i` above every finite threshold `B`. -/
theorem uniform_lower_bound_infinite_trace_unbounded
    (a : Nat → Int) (c : Int) (hc : 0 < c) (hbound : ∀ i, c ≤ a i) :
    ∀ B : Int, ∃ N : Nat, B < prefixSum a N := by
  intro B
  obtain ⟨N, _, hN_gt_B⟩ := exists_nat_gt_int B
  refine ⟨N, ?_⟩
  have hge := prefixSum_ge_mul a c hbound N
  have hNc_ge_N : (N : Int) ≤ (N : Int) * c := by
    have hN_nonneg : (0 : Int) ≤ (N : Int) := Int.natCast_nonneg N
    have h1 : (1 : Int) ≤ c := by omega
    have hmul := Int.mul_le_mul_of_nonneg_left h1 hN_nonneg
    simpa [Int.mul_one] using hmul
  omega

/-! ## Theorem 3: AdS Disjoint-Cut Tautology (Section 5.2, eqs. 5.44–5.45, p. 44) -/

def disjointSum {m : Nat} (w : Nat → Int) (c : Fin m → Int) (F : List (Nat × Fin m)) : Int :=
  listSum (F.map (fun p => w p.1 * c p.2))

def fiberList {m : Nat} (a₀ : Fin m) (n : Nat) : List (Nat × Fin m) :=
  (List.range n).map (fun i => (i, a₀))

theorem fiberList_nodup {m : Nat} (a₀ : Fin m) (n : Nat) :
    (fiberList a₀ n).Nodup := by
  unfold fiberList
  exact List.Pairwise.map (fun i => (i, a₀))
    (fun {i j} hne heq => hne (Prod.mk.inj heq).1)
    List.nodup_range

theorem disjointSum_fiberList {m : Nat} (w : Nat → Int) (c : Fin m → Int) (a₀ : Fin m) (n : Nat) :
    disjointSum w c (fiberList a₀ n) = prefixSum w n * c a₀ := by
  unfold disjointSum fiberList
  rw [List.map_map]
  have hcomp : ((fun p : Nat × Fin m => w p.1 * c p.2) ∘ (fun i => (i, a₀))) = (fun i => w i * c a₀) := rfl
  rw [hcomp, listSum_map_const_mul, listSum_range_eq_prefixSum]

def maxAbsFin {m : Nat} (c : Fin m → Int) : Int :=
  listSum ((List.finRange m).map (fun a => (Int.natAbs (c a) : Int)))

theorem natAbs_le_maxAbsFin {m : Nat} (c : Fin m → Int) (a : Fin m) :
    (Int.natAbs (c a) : Int) ≤ maxAbsFin c := by
  unfold maxAbsFin
  have hmem : a ∈ List.finRange m := List.mem_finRange a
  generalize List.finRange m = l at hmem
  have hnonneg_all : ∀ l' : List (Fin m), 0 ≤ listSum (l'.map (fun b => (Int.natAbs (c b) : Int))) := by
    intro l'
    induction l' with
    | nil => simp [listSum]
    | cons y ys ihy =>
      simp only [List.map_cons, listSum]
      have hy : (0 : Int) ≤ (Int.natAbs (c y) : Int) := Int.natCast_nonneg _
      omega
  induction l with
  | nil => contradiction
  | cons x xs ih =>
    simp only [List.map_cons, listSum]
    cases List.mem_cons.mp hmem with
    | inl hax =>
      subst hax
      have hrest := hnonneg_all xs
      omega
    | inr htail =>
      have hih := ih htail
      have hx : (0 : Int) ≤ (Int.natAbs (c x) : Int) := Int.natCast_nonneg _
      omega

theorem listMaxBound (ls : List Nat) :
    ∀ x ∈ ls, x < ls.foldr (fun y acc => max (y + 1) acc) 0 := by
  induction ls with
  | nil => intro x hx; contradiction
  | cons y ys ih =>
    intro x hx
    simp only [List.foldr_cons]
    cases List.mem_cons.mp hx with
    | inl hxy =>
      subst hxy
      exact Nat.lt_of_lt_of_le (Nat.lt_succ_self x) (Nat.le_max_left (x + 1) _)
    | inr hxt =>
      exact Nat.lt_of_lt_of_le (ih x hxt) (Nat.le_max_right (y + 1) _)

theorem listSum_part_lt (w : Nat → Int) (k : Nat) :
    ∀ vs : List Nat, vs.Nodup → (∀ x ∈ vs, x < k + 1) →
      listSum (vs.map w) ≤ listSum ((vs.filter (fun x => decide (x < k))).map w) +
        (if k ∈ vs then w k else 0) := by
  intro vs
  induction vs with
  | nil => intro _ _; simp [listSum]
  | cons y ys ihy =>
    intro hndv hltv
    have hnd_tail := (List.nodup_cons.mp hndv).2
    have hy_not_ys := (List.nodup_cons.mp hndv).1
    have hlt_tail : ∀ x ∈ ys, x < k + 1 := fun x hx => hltv x (List.mem_cons_of_mem y hx)
    have hrec := ihy hnd_tail hlt_tail
    have hy_lt : y < k + 1 := hltv y (List.mem_cons_self ..)
    by_cases hyk : y < k
    · have hne : k ≠ y := by omega
      simp [listSum, hyk, hne]
      omega
    · have hyeq : y = k := by omega
      subst hyeq
      simp [listSum, hy_not_ys] at hrec ⊢
      omega

theorem listSum_nodup_lt_le_prefixSum (w : Nat → Int) (hw : ∀ i, 0 ≤ w i) :
    ∀ (k : Nat) (us : List Nat), us.Nodup → (∀ x ∈ us, x < k) →
      listSum (us.map w) ≤ prefixSum w k := by
  intro k
  induction k with
  | zero =>
    intro us _ hlt
    cases us with
    | nil => simp [listSum, prefixSum]
    | cons x _ => have := hlt x (List.mem_cons_self ..); omega
  | succ k ihk =>
    intro us hndu hlt
    unfold prefixSum
    have hnd_lt : (us.filter (fun x => decide (x < k))).Nodup := hndu.filter _
    have hlt_lt : ∀ x ∈ us.filter (fun x => decide (x < k)), x < k := by
      intro x hx
      exact of_decide_eq_true (List.mem_filter.mp hx).2
    have hih_lt := ihk _ hnd_lt hlt_lt
    have hp := listSum_part_lt w k us hndu hlt
    have hwk : 0 ≤ w k := hw k
    by_cases hmem : k ∈ us
    · simp only [hmem, ite_true] at hp; omega
    · simp only [hmem, ite_false] at hp; omega

theorem nodup_fst_of_const_snd {m : Nat} (a : Fin m) :
    ∀ xs : List (Nat × Fin m), xs.Nodup → (∀ p ∈ xs, p.2 = a) →
      (xs.map Prod.fst).Nodup := by
  intro xs
  induction xs with
  | nil => intro _ _; simp
  | cons p ps ih =>
    intro hnd hEq
    rw [List.map_cons, List.nodup_cons]
    have hnd_p := List.nodup_cons.mp hnd
    refine ⟨?_, ih hnd_p.2 (fun q hq => hEq q (List.mem_cons_of_mem p hq))⟩
    intro hmem_fst
    obtain ⟨q, hq_mem, hq_fst⟩ := List.mem_map.mp hmem_fst
    have hp2 : p.2 = a := hEq p (List.mem_cons_self ..)
    have hq2 : q.2 = a := hEq q (List.mem_cons_of_mem p hq_mem)
    have hqp : q = p := Prod.ext hq_fst (hq2.trans hp2.symm)
    subst hqp
    exact hnd_p.1 hq_mem

theorem disjointSum_const_snd {m : Nat} (w : Nat → Int) (c : Fin m → Int) (a : Fin m) :
    ∀ xs : List (Nat × Fin m), (∀ p ∈ xs, p.2 = a) →
      disjointSum w c xs = listSum ((xs.map Prod.fst).map w) * c a := by
  intro xs
  induction xs with
  | nil => intro _; simp [disjointSum, listSum]
  | cons p ps ih =>
    intro hEq
    have hp2 : p.2 = a := hEq p (List.mem_cons_self ..)
    have ih' := ih (fun q hq => hEq q (List.mem_cons_of_mem p hq))
    unfold disjointSum at ih' ⊢
    simp only [List.map_cons, listSum, hp2, ih', Int.add_mul]

theorem disjointSum_filter_split {m : Nat} (w : Nat → Int) (c : Fin m → Int) (a : Fin m) :
    ∀ G : List (Nat × Fin m),
      disjointSum w c G =
        disjointSum w c (G.filter (fun p => decide (p.2 = a))) +
        disjointSum w c (G.filter (fun p => decide (p.2 ≠ a))) := by
  intro G
  unfold disjointSum
  induction G with
  | nil => simp [listSum]
  | cons p ps ih =>
    by_cases hpa : p.2 = a
    · simp [listSum, hpa, ih]
      omega
    · simp [listSum, hpa, ih]
      omega

theorem disjointSum_tail_bound {m : Nat} (w : Nat → Int) (hw : ∀ i, 0 ≤ w i)
    (c : Fin m → Int) (M : Int) (hM : ∀ n, prefixSum w n ≤ M) :
    ∀ F : List (Nat × Fin m), F.Nodup →
      -( (m : Int) * (M * maxAbsFin c) ) ≤ disjointSum w c F ∧
      disjointSum w c F ≤ (m : Int) * (M * maxAbsFin c) := by
  have hsum_nodup_nat : ∀ ls : List Nat, ls.Nodup → listSum (ls.map w) ≤ M := by
    intro ls hnd
    let N := ls.foldr (fun y acc => max (y + 1) acc) 0
    exact Int.le_trans (listSum_nodup_lt_le_prefixSum w hw N ls hnd (listMaxBound ls)) (hM N)
  have hlist_nonneg : ∀ l : List Nat, 0 ≤ listSum (l.map w) := by
    intro l
    induction l with
    | nil => simp [listSum]
    | cons x xs ih => simp only [List.map_cons, listSum]; have := hw x; omega
  have hfiber_bound : ∀ (as : List (Fin m)), as.Nodup →
      ∀ (G : List (Nat × Fin m)), G.Nodup → (∀ p ∈ G, p.2 ∈ as) →
        -( (as.length : Int) * (M * maxAbsFin c) ) ≤ disjointSum w c G ∧
        disjointSum w c G ≤ (as.length : Int) * (M * maxAbsFin c) := by
    intro as
    induction as with
    | nil =>
      intro _ G _ hG
      cases G with
      | nil => simp [disjointSum, listSum]
      | cons p _ => exact False.elim (List.not_mem_nil (hG p (List.mem_cons_self ..)))
    | cons a as' iha =>
      intro hnd_as G hndG hG
      have hnd_as' := (List.nodup_cons.mp hnd_as).2
      let G_a := G.filter (fun p => decide (p.2 = a))
      let G_rest := G.filter (fun p => decide (p.2 ≠ a))
      have hsplit := disjointSum_filter_split w c a G
      have hrest_in : ∀ p ∈ G_rest, p.2 ∈ as' := by
        intro p hp
        have ⟨hmem, hne⟩ := List.mem_filter.mp hp
        have hne' : p.2 ≠ a := of_decide_eq_true hne
        cases List.mem_cons.mp (hG p hmem) with
        | inl heq => exact absurd heq hne'
        | inr htail => exact htail
      have hrec := iha hnd_as' G_rest (hndG.filter _) hrest_in
      have hEq_a : ∀ p ∈ G_a, p.2 = a := fun p hp => of_decide_eq_true (List.mem_filter.mp hp).2
      have hnd_Ga_fst : (G_a.map Prod.fst).Nodup := nodup_fst_of_const_snd a G_a (hndG.filter _) hEq_a
      have hGa_eq : disjointSum w c G_a = listSum ((G_a.map Prod.fst).map w) * c a :=
        disjointSum_const_snd w c a G_a hEq_a
      have hW_nonneg : 0 ≤ listSum ((G_a.map Prod.fst).map w) := hlist_nonneg _
      have hW_le_M : listSum ((G_a.map Prod.fst).map w) ≤ M := hsum_nodup_nat _ hnd_Ga_fst
      have hnat := natAbs_le_maxAbsFin c a
      have hca_le : -(maxAbsFin c) ≤ c a ∧ c a ≤ maxAbsFin c := by omega
      have hmax_nonneg : 0 ≤ maxAbsFin c := by
        have : (0 : Int) ≤ (Int.natAbs (c a) : Int) := Int.natCast_nonneg _
        omega
      have hW_mul_le : listSum ((G_a.map Prod.fst).map w) * maxAbsFin c ≤ M * maxAbsFin c :=
        Int.mul_le_mul_of_nonneg_right hW_le_M hmax_nonneg
      have hGa_hi : disjointSum w c G_a ≤ M * maxAbsFin c := by
        rw [hGa_eq]
        exact Int.le_trans (Int.mul_le_mul_of_nonneg_left hca_le.2 hW_nonneg) hW_mul_le
      have hGa_lo : -(M * maxAbsFin c) ≤ disjointSum w c G_a := by
        rw [hGa_eq]
        have h1 : listSum ((G_a.map Prod.fst).map w) * (-(maxAbsFin c)) ≤
                  listSum ((G_a.map Prod.fst).map w) * c a :=
          Int.mul_le_mul_of_nonneg_left hca_le.1 hW_nonneg
        rw [Int.mul_neg] at h1
        omega
      have hlen : ((a :: as').length : Int) * (M * maxAbsFin c) =
                  (as'.length : Int) * (M * maxAbsFin c) + M * maxAbsFin c := by
        rw [List.length_cons, natCast_succ_mul]
      dsimp only [G_a, G_rest] at hrec hGa_hi hGa_lo
      rw [hsplit, hlen]
      omega
  intro F hndF
  have hfin := hfiber_bound (List.finRange m) (List.nodup_finRange m) F hndF (fun p _ => List.mem_finRange p.2)
  simpa [List.length_finRange] using hfin

/-- **Theorem 3 (AdS Disjoint-Cut Bi-Implication, eqs. 5.44–5.45):**
    For nonnegative weights `w i ≥ 0` and any `c : Fin m → Int` with `∃ a₀, c a₀ < 0`:
    The finite-subset sums `disjointSum w c F` are uniformly bounded across all finite
    subsets `F ⊆ ℕ × Fin m` IF AND ONLY IF the partial traces `∑_{i < N} w i` are bounded (`Z_β < ∞`). -/
theorem ads_544_unconditional_bounded_iff_weights_bounded
    {m : Nat} (w : Nat → Int) (hw : ∀ i, 0 ≤ w i)
    (c : Fin m → Int) (a₀ : Fin m) (ha₀ : c a₀ < 0) :
    (∃ B : Int, ∀ F : List (Nat × Fin m), F.Nodup → -B ≤ disjointSum w c F ∧ disjointSum w c F ≤ B)
    ↔ (∃ M : Int, ∀ n : Nat, prefixSum w n ≤ M) := by
  constructor
  · rintro ⟨B, hB⟩
    refine ⟨max B 0, fun n => ?_⟩
    have hfib := (hB (fiberList a₀ n) (fiberList_nodup a₀ n)).1
    rw [disjointSum_fiberList] at hfib
    have hp_nonneg : 0 ≤ prefixSum w n := prefixSum_nonneg w hw n
    have hca_le : c a₀ ≤ -1 := by omega
    have hmul : prefixSum w n * c a₀ ≤ prefixSum w n * (-1) :=
      Int.mul_le_mul_of_nonneg_left hca_le hp_nonneg
    have hneg : prefixSum w n * (-1) = -prefixSum w n := by omega
    rw [hneg] at hmul
    have hB_le : B ≤ max B 0 := Int.le_max_left B 0
    omega
  · rintro ⟨M, hM⟩
    refine ⟨(m : Int) * (M * maxAbsFin c), ?_⟩
    exact disjointSum_tail_bound w hw c M hM

/-! ## Theorem 4: Exact Counterexample & Cycle Partition Bijection (Section 5.2, eqs. 5.31–5.38) -/

def cycleToBasis (k : Nat) (j : Fin 100) : Bool × Nat :=
  if j.val = 0 then (true, k) else (false, 99 * k + (j.val - 1))

def basisToCycle (p : Bool × Nat) : Nat × Fin 100 :=
  match p with
  | (true,  k) => (k, ⟨0, by omega⟩)
  | (false, n) => (n / 99, ⟨(n % 99) + 1, by omega⟩)

theorem cycle_partition_left_inv (k : Nat) (j : Fin 100) :
    basisToCycle (cycleToBasis k j) = (k, j) := by
  by_cases hj : j.val = 0
  · simp only [cycleToBasis, hj, ite_true, basisToCycle]
    apply Prod.ext
    · rfl
    · apply Fin.ext; dsimp; omega
  · simp only [cycleToBasis, hj, ite_false, basisToCycle]
    have hj_lt : j.val < 100 := j.isLt
    apply Prod.ext
    · dsimp; omega
    · apply Fin.ext; dsimp; omega

theorem cycle_partition_right_inv (p : Bool × Nat) :
    (let cj := basisToCycle p; cycleToBasis cj.1 cj.2) = p := by
  rcases p with ⟨b, n⟩
  cases b with
  | true => simp [basisToCycle, cycleToBasis]
  | false =>
    simp only [basisToCycle, cycleToBasis]
    have hne : ¬ (n % 99 + 1 = 0) := by omega
    simp only [hne, ite_false]
    apply Prod.ext
    · rfl
    · omega

def Z_product (M : Nat) : Int := 40 * (M : Int)
def Z_cycle_times_401_upper (K : Nat) : Int :=
  401 * (80 * (K : Int)) - 40 * (39601 * (K : Int))

theorem section52_counterexample_product_diverges_pos :
    ∀ B : Int, ∃ M : Nat, 0 < M ∧ B < Z_product M := by
  intro B
  obtain ⟨M, hMpos, hMgt⟩ := exists_nat_gt_int B
  refine ⟨M, hMpos, ?_⟩
  unfold Z_product
  omega

theorem section52_counterexample_cycle_diverges_neg :
    ∀ B : Int, ∃ K : Nat, 0 < K ∧ Z_cycle_times_401_upper K < -B := by
  intro B
  obtain ⟨K, hKpos, hKgt⟩ := exists_nat_gt_int B
  refine ⟨K, hKpos, ?_⟩
  unfold Z_cycle_times_401_upper
  omega

/-! ## Theorem 5: Continuum Clock Incompatibility Triad (Sections 4.3 & 5.1, Appendix C.2) -/

theorem clock_continuum_incompatibility_triad
    (E0 : Nat → Int) (s : Int) (hs : 0 < s)
    (h_williamson : ∀ N : Nat, (N : Int) ≤ s * E0 N)
    (h_nonzero_PD : ∃ cap : Int, ∀ N : Nat, E0 N ≤ cap) :
    False := by
  obtain ⟨cap, hcap⟩ := h_nonzero_PD
  let cap_pos : Int := max cap 0
  obtain ⟨N, _, hN_gt⟩ := exists_nat_gt_int (s * cap_pos)
  have h1 := h_williamson N
  have h2 := hcap N
  have hcap_le : cap ≤ cap_pos := Int.le_max_left cap 0
  have hs_nonneg : 0 ≤ s := by omega
  have hmul : s * E0 N ≤ s * cap_pos :=
    Int.mul_le_mul_of_nonneg_left (Int.le_trans h2 hcap_le) hs_nonneg
  omega

/-! ## Theorem 6: Formal Proof-Graph & Circularity Audit (Section 5.2 ↔ Appendix E.3 & E.2) -/

structure GravityModel where
  OpenAmpConv       : Prop
  ClosedTraceConv   : Prop
  UniformLowerBound : Prop
  BoostNonZero      : Prop
  DimHD_Finite      : Prop
  AdS_Zbeta_Finite  : Prop
  Axiom5            : Prop
  Lorentzian1_4     : Prop
  EuclideanTypeI    : Prop

def SoundRules (M : GravityModel) : Prop :=
  ((M.UniformLowerBound ∧ M.ClosedTraceConv) → M.DimHD_Finite) ∧
  ((M.DimHD_Finite ∧ M.AdS_Zbeta_Finite) → M.Axiom5) ∧
  ((M.Lorentzian1_4 ∧ M.Axiom5) → M.EuclideanTypeI) ∧
  (M.ClosedTraceConv → M.AdS_Zbeta_Finite)

theorem sec52_does_not_entail_finite_dim :
    ∃ M : GravityModel,
      SoundRules M ∧ M.OpenAmpConv ∧ M.UniformLowerBound ∧ M.BoostNonZero ∧ ¬ M.DimHD_Finite := by
  refine ⟨{
    OpenAmpConv       := True
    ClosedTraceConv   := False
    UniformLowerBound := True
    BoostNonZero      := True
    DimHD_Finite      := False
    AdS_Zbeta_Finite  := False
    Axiom5            := False
    Lorentzian1_4     := True
    EuclideanTypeI    := False
  }, ?_⟩
  simp [SoundRules]

theorem circular_loop_admits_all_false_model :
    ∃ M : GravityModel,
      (M.ClosedTraceConv → M.DimHD_Finite) ∧
      (M.DimHD_Finite → M.Axiom5) ∧
      (M.Axiom5 → M.ClosedTraceConv) ∧
      M.OpenAmpConv ∧ M.UniformLowerBound ∧ M.BoostNonZero ∧
      ¬ M.DimHD_Finite ∧ ¬ M.ClosedTraceConv ∧ ¬ M.Axiom5 := by
  refine ⟨{
    OpenAmpConv       := True
    ClosedTraceConv   := False
    UniformLowerBound := True
    BoostNonZero      := True
    DimHD_Finite      := False
    AdS_Zbeta_Finite  := False
    Axiom5            := False
    Lorentzian1_4     := True
    EuclideanTypeI    := False
  }, ?_⟩
  simp

theorem appendix_E2_valid_conditionally (M : GravityModel)
    (hsound : SoundRules M) (h14 : M.Lorentzian1_4) (h5 : M.Axiom5) :
    M.EuclideanTypeI :=
  hsound.2.2.1 ⟨h14, h5⟩

#print axioms uniform_lower_bound_finite_trace_implies_dim_bound
#print axioms uniform_lower_bound_infinite_trace_unbounded
#print axioms ads_544_unconditional_bounded_iff_weights_bounded
#print axioms cycle_partition_left_inv
#print axioms cycle_partition_right_inv
#print axioms section52_counterexample_product_diverges_pos
#print axioms section52_counterexample_cycle_diverges_neg
#print axioms clock_continuum_incompatibility_triad
#print axioms sec52_does_not_entail_finite_dim
#print axioms circular_loop_admits_all_false_model
#print axioms appendix_E2_valid_conditionally

end Audit2611_34924
