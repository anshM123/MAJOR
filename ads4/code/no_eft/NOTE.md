# The cubic bound without a bulk effective theory

Note, 27 Nov 2026. Unrefereed, not committed. Nothing outside `ads4/code/no_eft/` was edited.
Paper: `prl/letter.tex` (assumption (ii), lines 86–90), `prl/supplement.tex` (Table `tab:inputs`, rows A2, H_K, L1–L7, LS1–LS6;
"Order of limits", lines 677–686; the c₁ = 0 chain, lines 698–806 and Table `tab:c1`; the angular route, lines 2073–2168),
`prl/sec_fullpos.tex` (Theorem FP, Eqs. (fp-target), (fp-H)).
Code here: `nl1_light_action.py`, `nl1b_sanity.py`, `nl1c_light_st.py` (outputs in `cache/`).

Grades as in the Supplement: [P] proved, [D] derived, [D-sketch], [DV] derived and checked numerically, [NE] numerical evidence,
[Q] quoted, [A] assumption, [I] open.

---

## 0. Summary

**Reading confirmed, with three corrections.** At k = 2, 3 the clause "light data = tree-level EFT + o(1/c_T) uniformly in spin" in
A2 is used for one thing only. The light side of the exact finite-c_T sum rule must be
c_T⁻¹L₁ + o(c_T⁻¹), where c_T⁻¹L₁ is the action of the T block with the exact ⟨TTT⟩. Your (a) and (b) are what that needs.
No other link of the cubic chain needs a bulk effective theory. The three corrections:

1. **L4 is not an assumption about the CFT at k = 2, 3.** It is a mathematical identity, used on one explicit function (the
   graviton-exchange Witten diagram) to evaluate a CFT-kinematic quantity: the action of the single T block. This quantity can also
   be computed directly from the light-ray factorisation, as in I3 and Steps 1–3 of the c₁ chain. Either way no bulk EFT enters.
   (The `used` column of row L4 says "Q", yet L3/c₁ = 0, which is item D, is graded "given … L4". That is a bookkeeping
   inconsistency.)
2. **The angular route τ_\* also uses the uniformity in spin.** Its step 𝓛_g → 0 ("[D] from planar factorisation") is a sum
   over all light states of spin ≥ 6. It needs the same uniform statement as (b), at relative order o(1) instead of o(1/c_T).
   The same replacement assumption covers it (§3.4).
3. **The remainder's size is not "a power Δ_g^p/c_T²".** On light double-twist states the stripped action of Ω^EE reaches
   |ĥ| ≈ 10⁶⁹¹ at Δ_g = 1000 (10²¹⁴ at Δ_g = 300). It scales like T_Δ·e^{πΛ/2} [NE, §4]. So the term-wise control of the remainder
   needs c_T ≫ e^{πΛ/2} ≈ 10⁶⁸¹. This holds with or without the bulk-EFT clause. The asymptotic statements are unaffected.

**What is derived.**

* **Reduction (Lemma 2, [P]).** The finite-c_T sum rule plus Theorem FP give, at every c_T,
  τ(c_T)² ≤ τ₀² + (57600/KΛ⁸)·|c_T R(c_T)|/n₀. Here R is the light remainder (light states other than 1, T and light single
  traces). So "leading order in 1/c_T" means exactly c_T R → 0, and a rate for R gives a finite-c_T bound.
* **Fixed spin (Lemma 3, [D given (F_J)]).** For each fixed J the double zeros (L1) make the spin-J part of R O(c_T⁻²) under
  fixed-spin factorisation. Only finitely many states are involved.
* **Spin summation (Theorem A, [D-sketch]; [P] for the step given its input).** (b) follows from one functional-independent
  spectral statement, (ii′_E) "detector decoupling". It requires only a bound on the stripped light action of Ω^EE (Lemma 5,
  [NE] for J ≤ 1000, [D-sketch] beyond):
  * the light states other than T carry a vanishing share of the paper's exact angular endpoint sum rule;
  * the share is weighted by (1+J)^{5/2}log, a mild weight in spin.
* **The same statement gives τ_\*.** (ii′_E) also gives 𝓛_g → 0 for τ_\* (Lemma 6, [D]).
* **Physicist form (ii′) ⇒ (ii′_E) (Corollary A′).** This uses the paper's own reflection-positive two-scale crossing bound, made
  uniform in c_T by Euclidean factorisation (Lemma 4, [D]; extension to twist below 10 [D-sketch]). (ii′) states: below the gap,
  the T̂ × T̂ OPE (T̂ = T/√c_T) is generalised-free up to corrections that vanish uniformly in spin:
  * every double twist [TT]_{n,J} has c_T γ_{n,J}² → 0 uniformly in (n, J);
  * every other light state and wrong-parity structure has squared OPE coefficient o(1/c_T), uniformly in spin shells.

  Factorisation gives γ = O(1/c_T) at each fixed spin, so (ii′) asks only for o(c_T^{-1/2}), uniformly.
* **Replacement assumption.** (ii′_E), or the stronger but more transparent (ii′), [A]. Both are pure CFT statements.

**What is not derived.** The uniform-in-spin part of (ii′_E)/(ii′) at spins above any fixed J₁.
* **Failing step.** The s-channel twist projection (z → 0 in the Lorentzian inversion integral) does not commute with the
  crossed-channel spin sums. This is the non-rigorous step of large-spin perturbation theory.
* **What the rigorous results give.** Pal–Qiao–Rychkov, van Rees, and the project's `LIGHTCONE_THEOREM.md` give leading-order,
  spin-averaged, rate-free statements at a fixed CFT, with constants that depend on the correlator. They give nothing uniform in
  c_T, and nothing at o(1/c_T) precision.
* **Status.** (ii′_E)/(ii′) stays [A]. It is much weaker than the bulk-EFT clause: no bulk locality, no derivative expansion, no structure
  of the O(1/c_T) light data, nothing about the heavy sector. The eikonal regime is unreachable below the gap (§5.3).

**Quartic bound.** It genuinely needs L7(b). The coefficient the sum rule measures, g̃₄, is a CFT quantity. But even its existence
as a limit needs first-order (single-zero) uniform control of the light data, which is stronger than (ii′). Its identification with
the bulk R⁴ coefficient g₄ needs the light double-trace data to be those of a tree-level bulk EFT with no conspiracy between the
heavy sum and the contact tower (CHMRS's eq. (5.3)). g₄ is a bulk-Lagrangian parameter, so no CFT statement can replace this.
Confirmed, not attempted.

---

## 1. Setting and the precise meaning of "leading order"

### 1.1 Light states and the remainder

* **Light states.** States of twist < Δ_g in T × T.
* **The remainder set.** 𝓛′ is the light set minus {1, T} and minus the light single traces of spin 0, 2 (rows LS).
* **The sum rule.** At every finite c_T, H_K [D] and Rindler absolute convergence (R3) give, for the cut-off pieces of Ω^EE:

  `Σ_{heavy} p_O ω[O] + p_T ω[T] + R(c_T) + Σ_{light ST} p_O ω[O] = 0,   R(c_T) := Σ_{O∈𝓛′} p_O ω[O]`   (1.1)

  The identity drops because L[T]|Ω⟩ = 0. The sum is absolutely convergent over all blocks (R3).

### 1.2 Exactness of the T term (Lemma 1)

**Lemma 1 [D, given L1's kinematics].** p_T ω[T] is an exact quadratic form in (t₄, α_T) at every c_T:

`p_T ω[T] = n₀ G [ K τ²Λ⁸/57600 − 16 L_A ]`,   with 8πG = 256/c_T.

Here n₀ > 0 is a fixed normalisation and K is the |ĝ₃|² coefficient (L2, I3). There is no linear term (c₁ = 0, every ν).

*Proof.*
* The action of a single conformal block is fixed by conformal kinematics and the block's OPE data.
* For T those data are (c_T, t₄, α_T) exactly, because T is unique (A1′).
* The c₁ = 0 chain, Steps 1–3, is a statement about the T-block term ⟨Ψ₄|A₃P_TA₁|Ψ₂⟩. The one-point matrices are affine in t₄
  (and α_T) at fixed c_T.
* K (K₀, δK₀, H17) is computed from the same light-ray factorisation (I3).

No 1/c_T expansion enters: G is the exact 1/c_T. ∎

### 1.3 The reduction (Lemma 2)

**Lemma 2 [P].** Assume (1.1), (Pos) on heavy states, and no light single trace couples to TT. Then at every c_T

`τ(c_T)² ≤ τ₀² + (57600 / (K Λ⁸)) · |c_T R(c_T)| / (32 n₀/π)`,   τ₀ = (240/Λ⁴)(16L_A/K)^{1/2}.   (1.2)

Hence c_T R → 0 implies limsup_{c_T→∞} τ(c_T) ≤ τ₀. A bound |R| ≤ C_R/c_T² gives the explicit finite-c_T statement
τ² ≤ τ₀² + O(C_R/(c_TΛ⁸)).

*Proof.*
* (Pos) holds on every heavy state at every c_T, by Theorem FP at Δ_g = 1000: it is a statement about (Δ, J) and parity only.
* So Σ_heavy ≥ 0 in (1.1), and p_T ω[T] + R ≤ 0.
* Insert Lemma 1 (n₀G = 32 n₀/(πc_T)) and divide. ∎

**What "at leading order in 1/c_T" should mean.** The Letter's phrase should be read as limsup_{c_T→∞} τ(c_T) ≤ τ₀, the limit taken
on the exact finite-c_T inequality. That is what the Supplement's "c_T → ∞ limit of the finite-c_T inequality" says. The bound
applies to the exact t₄, α_T of each member of the family only up to the remainder term of (1.2).

### 1.4 What (b) is

(b) is **exactly** the statement c_T R(c_T) → 0. Everything else in the order-of-limits paragraph of the Supplement is [P]/[D]
without bulk input.

---

## 2. Audit (Task 1): every use of A2

"EFT?" asks whether the bulk-EFT clause is needed, or only large-c_T factorisation (a)/(b) or less.

### 2.1 Cubic results

| # | place | what is used | EFT? | replacement |
|---|---|---|---|---|
| 1 | H_K order of limits (Supp. 677–686; row H_K "order of limits [D given A2]") | c_T → ∞ at fixed Δ_g; light side = c_T⁻¹L₁ + o(c_T⁻¹) | **no**: only (b) | Lemma 2 + Theorem A |
| 2 | Row L1: "[TT] at double zeros, contacts drop" | double zeros: CFT kinematics of light-ray amplitudes [D/P]. "Contacts drop" is block language for "O(1/c_T) double-trace data are killed by double zeros" | **no**. "Contacts" need not exist as objects: at k = 2, 3 the O(1/c_T) contact and exchange double-trace data enter only through γ², δP·γ. Their summed effect is part of R | (b) |
| 3 | Row L1: −Ĉ′_{EE,ν}\|_E = 16πG/(ν²+9/4); "H7's value is the Regge arc at n = k = 2" | value of the T-block action | **no** | Lemma 1. Regge-arc route = L4 on the exchange diagram (item 5) |
| 4 | Row L2 / I3: K = K₀(1 + δ₃), δK₀ [certified], H17 closed form | T-block light-ray factorisation | **no** | none needed. The H17 asymmetry band is [M], unrelated to EFT |
| 5 | Row L3 / Table c1 link (i): "Lemma-R value of a tree term = term-wise s-channel value (L4)" | identity on explicit Witten diagrams | **no physical EFT** | (i) is needed only if the light side is *defined* as the Regge-arc value of a tree correlator. With the light side defined block-wise (as in H_K and Sec. `sec:arc`), it drops out at k = 2, 3. See item 5a |
| 5a | the same, if kept | that the Einstein exchange diagram is (T block) + (its [TT] tower), and that the tower is killed by the double zeros | no | a true identity for that explicit function: its double-trace data are explicit and dispersive, with finitely many Mellin residues. It is a computing device and says nothing about the CFT |
| 6 | Table c1 link (ii): "at O(1/c_T) only the T block survives" | (a) + (b) | **no**: only (b) | Theorem A under (ii′) |
| 7 | L1: "parity-violating admixtures enter at O(c_T⁻²) [D]" | at fixed J, factorisation. Summed over J, uniformity | **no** | included in (ii′), second clause |
| 8 | Remark of Sec. `sec:allgap`.A: "states of twist < Δ_g are among those of twist < Δ_g′, whose data are tree level" | monotonicity of the assumption in the gap | **no** | (ii′) at gap Δ_g′ implies (ii′) at every Δ_g ≤ Δ_g′ [P]. It is a sup over a subset; states of twist in [Δ_g, Δ_g′) become heavy and are covered by (Pos) at Δ_g |
| 9 | Theorem FP; T0; Theorems W, L′ ("Assume A1′, A2") | per-state statements about (Δ, J, parity) | **no**. A2 is idle there | drop "A2" from W, L′ or leave it (harmless) |
| 10 | Energy-flux corollary, \|t₄\| ≤ τ₀, α_T bound, \|n_s − n_f\| ≤ … | inherit item (1) | no | as item (1) |
| 11 | Einstein-rigidity theorem | items (1) and (3) | no | (ii′) for both (§3.4) |
| 12 | Angular route τ_\* (Supp. 106, 2112–2114): "𝓛_g → 0 [D] from planar factorisation"; input (3) "fixed-window planar regularity" | 𝓛_g = Σ_{light, J ≥ 6} \|a_{O,4,2}\|² → 0: a sum over infinitely many spins | **no**, but the [D] is really [D given uniformity in spin] | §3.4: (ii′) ⇒ 𝓛_g → 0 [D given two-scale crossing] |
| 13 | All-gap rigidity (Theorems RP, MP; sec_fullpos 427–441) | same light-side logic, other weights | no | (ii′) plus the functional-side bound of Lemma 5 for those weights. **Not checked**: their weights have partial-fraction coefficients ~10⁵⁷⁰⁰, so the constant of Lemma 5 is expected to be larger still |
| 14 | Conditional remark Δ_g ≥ 300 (Sec. `sec:cor`) | as item 1 | no | (ii′) |
| 15 | LS rows (only if light single traces are present). LS4 "(E′): CHMRS's light-side formula carried to spinning k = 2, 3" | treats a light single trace plus its double-trace dressing as a tree exchange Witten diagram | **bulk-flavoured, removable at k = 2, 3** | the dressing sits at double zeros and goes into R. The single-trace block's action can be evaluated directly. Side check: §6 |
| 16 | ABJM: "A1′ and A2 hold at leading order (by the spin count alone)" | planar IIA | unchanged status [D-sketch] | (ii′) holds in planar ABJM-IIA if its light double-twist anomalous dimensions are O(Δ_g^p/c_T) uniformly. That is the supergravity statement; same grade |

### 2.2 Quartic results (brief)

* **k = 4, 5.** The first double-trace families have only **single** effective zeros. They enter at O(c_T⁻¹) through P⁽⁰⁾γ⁽¹⁾ω′ and
  δP⁽¹⁾ω (Supp. 680).
* **What a CFT-only light side would need.** Their Σ_J must converge to a c_T-independent limit, i.e. c_T γ_{n,J} → γ⁽¹⁾_{n,J}
  uniformly with summability against single-zero weights. That is a *first-order* uniform statement, strictly stronger than (ii′).
* **What the sum rule measures.** The quantity g̃₄ = −½𝒜₄[light] is CFT-defined, given that statement.
* **From g̃₄ to g₄.** The step needs:
  * L4–L6: the arc evaluation on Witten diagrams; Lemma A; the calibration by a heavy scalar;
  * L7(a): the term-wise value of a heavy exchange's double traces [D];
  * L7(b): the light double-trace data equal those of a tree-level bulk EFT, with no conspiracy between the heavy sum and the contact
    tower [A]. This is CHMRS's own technical assumption, their eq. (5.3), as cited in the Supplement.
* **Why L7(b) cannot be removed.** g₄ is defined as a coefficient of a bulk Lagrangian. A pure CFT statement can at most define g̃₄.
  So L7(b) is genuinely needed for "g̃₄ = g₄", and for nothing else.
* **Confirmed, not attempted.**

---

## 3. The spin summation

### 3.1 Notation

* **Detector amplitudes.** a_{O,b,k} are the normalised rest-frame detector amplitudes of the angular route (Sec. `sec:angular`).
  Their normalisation is the one in which the endpoint rule reads

  `Σ_{O≠T, b} w_{O,b} |a_{O,b,0}|² = 3 + 𝒟,  w_{O,b} = [T_Δ − b(b−1)]/3 + 3/2 ≥ 3/2,  T_Δ = (Δ−1)(Δ−2)`.   (3.1)

  𝒟 = 0 [P given the angular inputs].
* **Double zeros.** |a_{O,b,0}|² carries the zero factor of the light-ray amplitude:
  * parity-even structures: |K|² ∝ sin²(πτ/2);
  * parity-odd structures: ∝ cos²(πτ/2)Γ(τ/2 − ½)², vanishing at τ = 3, 5, … (L1; τ ≥ 7 for J = 0, 2).

  Write S_b(τ_O) for the factor and |a_{O,b,0}|² = S_b(τ_O)·α_{O,b}, so that α_{O,b} is the "stripped" amplitude.
* **Per-state action (Eq. (fp-target), with B(1−r)³ = (2πGΛ⁴/9)T_Δ taken out).**

  `ω[O] = (2πGΛ⁴/9) Σ_b |a_{O,b,0}|² ĥ_b(Δ, J)`,

  `ĥ_b = T_Δ Σ_{k=0}^{J−b} D_k² U_{b,k} [ Z_k(1/(k+3) − r/(k+1)) + (2Y_k + V_{b,k}Z_k + σ_bΛ²Q_k)/T_Δ ]`.   (3.2)

  ĥ is finite at r = 1, where B and (1−r)³ are separately singular (the n = 0 double twist at γ = 0).
* **Dyadic shells.** 𝒮_S = {O ∈ 𝓛′ : J_O ∈ [S, 2S)}, with light shell mass λ_S(c_T) := Σ_{O∈𝒮_S} Σ_b |a_{O,b,0}|².

### 3.2 The fixed-spin step (Lemma 3)

**(F_J) Fixed-spin factorisation [A; it is what large-N factorisation of ⟨TTTT⟩ on Euclidean compacta gives, spin by spin, for a
discrete spectrum with finitely many light single traces].** For each fixed J:
* the spin-J operators of 𝓛′ are double twists [TT]_{n,J}, finitely many n since τ < Δ_g, with γ_{n,J} = O(1/c_T) and
  P_{n,J} = P⁽⁰⁾_{n,J} + O(1/c_T);
* plus operators and wrong-parity structures with total squared coefficient O(c_T⁻²).

**Lemma 3 [D given (F_J)].** For each fixed J₁, c_T Σ_{O∈𝓛′, J_O<J₁} |p_O ω[O]| = O(1/c_T).

*Proof.*
* Finitely many states are involved.
* On double twists, ω = P·S(τ)·W with S = O(γ²) = O(c_T⁻²) and W smooth.
* On the rest, p = O(c_T⁻²). ∎

(Low spins, J ≤ 2 included, need only this.)

### 3.3 Uniform envelopes (Lemmas 4, 5) and Theorem A

**Lemma 4 (shell envelope) [D given the two-scale crossing bound of Sec. `sec:angcross` and Euclidean factorisation; twist < 10:
D-sketch].** There is a c_T-independent 𝒦̂ such that, for all shells S ≥ 4 and all large c_T,

`λ_S(c_T) ≤ 𝒦̂ Δ_g² S⁻⁵ · c_T · sup_{O∈𝒮_S} S_b(τ_O) + 𝒦̂ Δ_g² S⁻⁵ · c_T · μ^B_S(c_T)`.   (3.3)

Here μ^B_S is the relative weight (fraction of the GFF shell weight) of the states in 𝒮_S that are not double twists near a zero.

*Proof.*
* **The bound at a single c_T.** The two-scale crossing bound gives
  Σ_{shell} 4^Δ2^{−J}(|λ_{++}|² + |λ_{+−}|²) ≲ 𝒦 S^{3/2}U^{3/2}, with 𝒦 depending on the CFT only through T-component norms and two
  Euclidean four-point values. It then gives Σ_shell Σ_b |a_{O,b,0}|² ≲ 𝒦 S⁻⁵U², obtained by bounding |K|² with sin² ≤ 1 (Wendel).
  Keeping the zero factor instead of bounding it by 1 gives the sup.
* **Uniformity in c_T.** In the normalisation of (3.1), the T norms and the Euclidean values are O(c_T²) at leading order, the
  generalised-free part. So 𝒦 = c_T·𝒦̂(1 + o(1)) with 𝒦̂ fixed by GFF data [D given factorisation at the two Euclidean points;
  this normalisation bookkeeping was done here once and not independently re-derived].
* **Light twist.** U ≤ Δ_g/2 there. The proof in the Supplement assumes h = τ/2 ≥ 5. For 1 ≤ h < 5 the same inequalities hold
  with explicit Gamma-function constants [D-sketch]. ∎

**Lemma 5 (functional side) [NE for J ≤ 1000 at Δ_g = 1000 and J ≤ 300 at Δ_g = 300; D-sketch beyond; the continuation of
(fp-target) to light states is D-sketch].** On light states (Δ − J < Δ_g),

`|ĥ_b(Δ, J)| ≤ T_Δ [ C_h(Δ_g) + C′ (1+J)^{5/2} log(2+J) ]`,   C_h(1000) ≈ 10^{685.6},  C_h(300) ≈ 10^{209},   (3.4)

and C_h ≈ e^{πΛ/2}·10^{4.5±0.5}. Evidence and the asymptotic argument: §4.

**The minimal assumption: detector decoupling.**

**(ii′_E)** `𝓔^{(q)}_light(c_T) := Σ_{O∈𝓛′} (1+J_O)^q Σ_b w_{O,b} |a_{O,b,0}|² → 0` as c_T → ∞, for q = 5/2 with a log, i.e.
weight (1+J)^{5/2}log(2+J).

Without the weight factor, 𝓔⁽⁰⁾_light is the light states' share of the exact endpoint sum rule (3.1), and it is ≤ 3 at every c_T.
(ii′_E) says that, as c_T → ∞, the whole of the "3" is carried by states above the gap, with a little room in spin. Physically:
**energy detectors acting on the graviton state T|Ω⟩ leak into light states only by an amount that vanishes as c_T → ∞, uniformly
in spin.** This is a statement about the light-ray algebra of the CFT and nothing else. It does not mention the functional, and at
each fixed spin it follows from (F_J).

**Theorem A [D-sketch through Lemma 5; given (3.4) the step is P].** (ii′_E) implies c_T R(c_T) → 0. Hence (Lemma 2)
limsup τ ≤ τ₀ at Δ_g = 1000. The same holds for every theory with gap Δ_g′ ≥ 1000, by row 8 of §2.1.

*Proof.* From (3.2) and (3.4), on light states, using w_{O,b} ≥ T_Δ/3 − 4 and T_Δ ≥ J(J+1):

|ω[O]| ≤ (2πGΛ⁴/9)·3(C_h + C″(1+J)^{5/2}log(2+J))·Σ_b w_{O,b}|a_{O,b,0}|²   (up to O(1) constants).

So c_T|R| ≤ (2π(c_TG)Λ⁴/3)·max(C_h, C″)·𝓔^{(5/2,log)}_light(c_T) → 0, since c_TG = 32/π is fixed. ∎

**A sufficient condition in terms of OPE data (physicist form).**

**(ii′)** (F_J) holds, and for some fixed S₁
`lim_{c_T→∞} sup_{S ≥ S₁} [ c_T sup_{O∈𝒮_S, DT} S_b(τ_O) + c_T μ^B_S(c_T) ] = 0`.

In words: on light double twists c_Tγ² → 0 uniformly in spin. All other light states and wrong-parity structures carry o(1/c_T) of
the GFF shell weight, uniformly in spin.

**Corollary A′ [D given Lemma 4].** (ii′) ⇒ (ii′_E).

*Proof.*
* Shells below S₁: (F_J) and Lemma 3.
* Shells S ≥ S₁: by (3.3), the shell contributes at most (2S+Δ_g)^{9/2}log S · 𝒦̂Δ_g²S⁻⁵·ε(c_T) to 𝓔, with ε(c_T) → 0 by (ii′).
* The dyadic sum converges (margin S^{−1/2}), so the total → 0. ∎

The converse fails: (ii′_E) does not force a uniform sup over individual states. So (ii′_E) is the weaker, and minimal, assumption.
(ii′) is the form a physicist would check.

### 3.4 The angular route τ_\*

**Lemma 6 [D given (ii′_E); no functional-side input].** 𝓛_g := Σ_{O light, J≥6} |a_{O,4,2}|² → 0 as c_T → ∞.

*Proof.*
* By the ladder ratio of the angular route, |a_2|²/|a_0|² = 100(C−20)(C−30)/((T−20)(T−30)). This is ≤ 100(1 + O(1/J)) on light
  states, since T_Δ ≥ C there.
* Also w_{O,4} ≥ 3/2. So 𝓛_g ≤ 70·𝓔⁽⁰⁾_light ≤ 70·𝓔^{(q)}_light → 0. ∎

The light single traces of spin 0, 2 do not enter 𝓛_g (J ≥ 6). T is excluded by construction. So τ ≤ τ_\*(Δ_g) holds in the limit
under (ii′_E) and the other inputs of item (3). Even the unweighted q = 0 form of (ii′_E) suffices here. This replaces "𝓛_g → 0 [D] from planar factorisation" (Supp. 106, 2112). The
planar-factorisation argument there is correct state by state, but it does not control the spin sum.

---

## 4. Numerics: the stripped action of Ω^EE on light states (Lemma 5)

`nl1_light_action.py` evaluates (3.2) at light (Δ, J). It uses the certified moment balls of `fullpos/cache/moments_<Δg>_K<K>.json`
(read-only), exact Q_k, and ball arithmetic. The relative ball radius is ≤ 10⁻¹²⁴ in every entry quoted below.

**Check (i) [DV].** At heavy states (`nl1b_sanity.py`), ĥ(1−r)³/T equals the FP evaluator `fp_poly.Action.H_ang` to all printed
digits (four states).

**The table [NE]** gives sign and log₁₀|ĥ_b| at Δ_g = 1000 (`cache/nl1_1000_1000.txt`).

| b | twist | J = 4 | J = 8 | J = 64 | J = 256 | J = 512 | J = 1000 |
|---|---|---|---|---|---|---|---|
| 0 | 2 | +26.0 | +42.3 | +202.2 | +515.3 | +688.1 | −691.4 |
| 4 | 2 | −7.2 | −26.3 | +195.2 | +510.7 | +688.3 | +691.6 |
| 0 | 4 | +24.4 | +40.1 | +198.2 | +510.1 | +681.0 | −685.0 |
| 0 | 12 | +21.4 | +35.6 | +187.1 | +494.4 | −662.4 | +665.1 |
| 0 | 42 | +17.4 | +28.7 | −163.4 | +452.8 | −604.9 | +607.0 |
| 0 | 202 | +12.1 | +18.4 | −108.6 | +324.1 | −402.4 | +404.3 |
| 0 | 998 | +5.7 | +5.2 | +2.8 | +1.7 | +1.8 | +2.2 |

**Observations.**

* **Growth.** For J ≲ Λ, |ĥ| grows like (Λ²/m²)^J-type factors (local log–log slopes 50–680). It then plateaus at
  |ĥ|/T_Δ ≈ 10^{685.4} (Δ_g = 1000) and ≈ 10^{208.7} (Δ_g = 300, `cache/nl1_300_300.txt`). Compare e^{πΛ/2} = 10^{681.2} and
  10^{203.6}: the plateau is e^{πΛ/2} times 10^{4–5}.
* **Mechanism.** This is the one the Supplement describes for single light blocks ("grows like e^{πν/2} above its mass"). The
  weight extends to ν ≈ Λ + 100, and a light state's per-ν image grows from ν ≈ m up to there. The partial ladder sums inherit the
  10⁶⁸⁵ cancellation scale of Theorem FP (its "fact (a)") without the cancellation, because for light states U_{b,k} ≈ 1 up to the
  ladder end instead of ≈ r^k.
* **No sign.** The sign changes with J and twist, so light double twists cannot be moved to the heavy side.
* **Continuity.** In twist (`nl1b_sanity.py` (iii), J = 64): log₁₀|ĥ| falls smoothly 202 → 136 → 62 → 10 → 2.8 from twist 2 to 900
  to 998, and joins the heavy values (2.83 at twist 1000).
* **Beyond J = 1000 (D-sketch).**
  * For n = 0, ĥ/T is a partial sum Σ_{k≤J}P_k of the single-sequence coefficients.
  * Their generating function is F(z) = (2/(Λ(1−z)³))∫E_θ P_{−1/2+iν}((1+z)/(1−z)) dν (sec_fullpos, [D]). Near z = 1 it behaves like
    (1−z)^{−5/2}(log), so P_k = O(k^{3/2} log k) for k ≫ Λ.
  * This gives the C′J^{5/2} term in (3.4). It is irrelevant numerically until J ~ 10²⁷⁴.
  * For n ≥ 1, U_{b,k} → r^k smoothly. The giant piece then decays like n·Λ²/J³ relative to C_h.
* **Consequence for "how large c_T must be".**
  * With (3.4) and γ ~ n³/(c_T J), the term-wise bound on R is of order C_h·c_T⁻²·(polynomial in Λ).
  * The cubic inequality (1.2) is therefore numerically meaningful only for c_T ≳ 10⁶⁹⁰ at Δ_g = 1000. Term-wise, without
    cancellations between spins, this holds **with or without** the bulk-EFT clause: the EFT would fix γ but not ĥ.
  * The Supplement's parenthesis "(a power Δ_g^p/c_T²)" should be replaced; see PROPOSED_EDITS.md, E6.
  * Cancellations across J in the true one-loop light value cannot be excluded, but nothing in hand shows them.

---

## 5. Can (ii′) be derived? (Task 2)

### 5.1 What follows from the listed axioms

* **Fixed spin.** It follows from Euclidean factorisation, with a finite light single-trace spectrum (needed in any case: a spectrum
  accumulating below the gap is already [I] in the Supplement) [D, (F_J)].
* **J-summability at each fixed c_T.** It follows from Rindler positivity (R3) and Regge boundedness of the exact correlator. The
  sum converges, but with no uniformity [D].
* **Uniform OPE envelopes in spin shells.** They follow from the reflection-positive two-scale crossing bound (Lemma 4) [D]. This is
  what makes the uniform-in-spin question one about anomalous dimensions only, not about OPE coefficients.
* **The tree-level shape of the large-spin data.** Lorentzian inversion holds for J > J₀ (Caron-Huot; spinning form:
  Kravchuk–Simmons-Duffin). Since J₀ ≤ 1 at finite c_T (R1; KKSDZ Sec. 4.2.4), it holds for every J ≥ 2.
  * Inverting the single t-channel T block gives an explicit γ^{(T)}_{n,J} = O(p_T/J^{τ_T}) = O(1/(c_T J)) for bounded n,
    uniformly in J [D for a single block].
  * Light single traces give O(1/(c_T J^{τ_L})).
  * Light t-channel double twists enter dDisc_t through sin²(πγ/2). Their contribution is quadratic in the unknown, so a bootstrap
    estimate is self-consistent.
  * Heavy t-channel operators contribute O(J^{−τ_h}), τ_h ≥ Δ_g, each with OPE weight O(1/c_T).
  * The paper's own positivity could bound their Ω-weighted total (Theorem FP gives Σ_heavy p N_s ≤ ε⁻¹|light side| = O(1/c_T)).
    That would give a non-circular handle on the heavy part. It is not the route rejected in Supp. 693, which used the tree-level
    heavy Regge tail.

### 5.2 The step that fails

To extract γ_{n,J} for a twist family n ≥ 0 from the inversion integral one needs the small-z expansion of C^t(z, β) up to
z^{Δ_g/2}. Two things go wrong:

* **Term-wise expansion diverges.** The t-channel expansion of dDisc_t G is at the edge of its convergence region there (1−z → 1).
  Its spin sums at fixed t-channel twist cannot be expanded term-wise. For GFF-like OPE coefficients the coefficient of
  z^m log z behaves like Σ_{J′} P_{J′}4^{h̄′}h̄′^{2m}, which diverges. The resummed spin sums are what produce the s-channel images.
* **The rigorous replacements are not uniform.** They are Hardy–Littlewood/Tauberian (Pal–Qiao–Rychkov, arXiv:2212.04893; van Rees,
  arXiv:2412.06907, both read in full in `ads4/DLC_ESTIMATES.md`; the averaged first-moment theorem of `ads4/LIGHTCONE_THEOREM.md`).
  * They give twist accumulation and a spin-averaged leading coefficient, with **no rate**.
  * They hold at a **fixed CFT**, with constants depending on the correlator.
  * They give a second moment of (τ − 2Δ) "of smaller order than the first", not O(c_T⁻²).
  * `DLC_ESTIMATES.md` found no rigorous quantitative large-spin correction in d ≥ 3 through Nov 2026.

So even a version of `LIGHTCONE_THEOREM` with constants uniform in c_T would give c_T·(second moment) = O(1), not o(1). That bounds R
by O(1/c_T), which only shifts the bound by an uncontrolled O(1) amount. (ii′) needs o(1).

**Minimal replacement.** Detector decoupling (ii′_E): 𝓔^{(5/2,log)}_light → 0. Its sufficient OPE-data form is (ii′) of §3.3:
c_Tγ² → 0 uniformly on light double twists, and o(1/c_T) weight elsewhere, in spin shells.

**Possible partial upgrade (not done).** A Nachtmann-type monotonicity theorem for the leading-twist trajectory of T × T would give
the n = 0 part. If |γ_{0,J}| is non-increasing for J ≥ J_\*, uniformity follows from (F_J) at J_\*. Komargodski–Zhiboedov discuss
convexity of the leading twist for scalar external operators. We know of no such theorem for spinning external T, nor any
monotonicity for subleading families. So this stays [A].

### 5.3 The eikonal question

* **Bounded twist means no eikonalisation.** Light double twists have h = τ/2 < Δ_g/2, bounded independently of c_T.
* **Bulk heuristic.** In the AdS₄ eikonal (Cornalba–Costa–Penedones–Schiappa), graviton exchange gives
  |γ_{n,J}| ≈ c·G h³/h̄ for h̄ ≫ h, and at most c·G Δ_g^p for J ≲ n. Hence
  sup_{light} c_T γ² ≲ c²(c_TG)² Δ_g^{2p}/c_T → 0.
* **When the eikonal phase is O(1).** That needs h ≳ (c_T J)^{1/3}, i.e. twist ≳ c_T^{1/3} ≫ Δ_g. Those states are heavy and
  covered by (Pos).
* **Stringy exchanges** add O(1/c_T)·(h/Δ_g)^{#}J^{−Δ_g} to the light γ's.
* **CFT version.** It is the single-block inversion statement of §5.1, which is uniform in J [D]. The gap is the sum over the other
  t-channel operators (§5.2).
* **Order of limits.** With the limits as used (finite c_T first, c_T → ∞ last), the eikonal regime lies entirely on the heavy side.
* **Grade.** The uniformity (ii′) is physically expected with a large margin, c_Tγ² ~ Δ_g⁶/c_T. But it is [A].

### 5.4 Comparison with the literature

* **CHMRS.** They assume the light data are those of a tree-level AdS EFT with a derivative expansion (their Sec. 1), plus the
  no-conspiracy assumption (5.3). Their order-of-limits logic (Sec. 1.1.3) is the one used here. (ii′) replaces their tree-level EFT
  assumption **at k = 2, 3 only**, where double zeros make the O(1/c_T) light double-trace data invisible. It cannot replace it at
  k = 4, 5 (§2.2).
* **AHKT** (Einstein gravity 3-point functions from CFT). As we understand it, they derive the parametric vanishing of t₄-type
  structures from CFT causality at large c_T with a sparse spectrum. Stress-tensor exchange dominates the relevant (lightcone/Regge)
  limit, which plays the role of (ii′). Their statement is parametric in the gap; ours carries explicit constants and needs
  uniformity only at relative order o(c_T^{-1/2}) in γ.
* **FKPS / KZ.** Large-spin asymptotics at a fixed theory. They give the shape of γ_{n,J} but no uniformity in c_T.
* **Caron-Huot.** Lorentzian inversion, valid for J > J₀. It is used only for the single-block statement.
* **KKSDZ.** J₀ ≤ 1 (R1) and the light-ray algebra (detector amplitudes, the endpoint sum rule's framework).
* **Mukhametzhanov–Zhiboedov.** The rigorous **lightcone** Tauberian results are Pal–Qiao–Rychkov and van Rees.
  Mukhametzhanov–Zhiboedov's Tauberian work (arXiv:1808.03212) concerns the Euclidean OPE limit (per `ads4/DLC_ESTIMATES.md`).
  Citing them for the lightcone statement would be inaccurate.

---

## 6. Side check: light single traces of spin 0, 2 without (E′)

Not needed for the headline. `nl1c_light_st.py`.

**Block-wise route.** At k = 2, 3 the double-trace dressing of a light single trace sits at the double zeros. Its contribution goes
into R. So the light single trace's total contribution at O(1/c_T) is its block's action. LS4's Witten-diagram language is not
needed. The block action follows from (3.2) continued to (Δ, J ≤ 2) [D-sketch, the same continuation as Lemma 5]. The ladder has at
most three terms.

**Results.**
* **J = 0.** Exactly ĥ/Λ² = 20.771893·T_Δ/Λ² + 11.317149 (Δ_g = 1000; 11.318397 at 300). This is the flat polynomial
  F_D⁽⁰⁾ = 11.317027 + 20.771893 m² of the Supplement, with m² = T_Δ/Λ² = (Δ−1)(Δ−2)/Λ². In this representation the AdS mass shift
  δ(Δ_L) of LS4 is zero.
* **Positivity at J = 0.** It holds for every unitary Δ ≥ ½ (T_Δ ≥ −¼) [P given the continuation].
* **J = 2.** min ĥ/Λ² = 1.94 on a Δ grid (3.05 … Δ_g + 2), attained at Δ ≈ 0.71Λ. Positive [NE].

**Caveat.** The Supplement's LS paragraph says light states "cannot be treated as heavy-type images". That warning concerns
leading-Regge images, not the exact ladder formula. But the continuation of (fp-target) is itself unproved for light states, so this
is only suggestive. It would replace (E′) [D-sketch] by "the ladder formula at light J ≤ 2" [D-sketch], with no change of grade. It
deserves a separate check against a direct evaluation of the block action.

---

## 7. What remains assumed for the cubic results

1. A1′ (higher-spin gap; unique T) [A]: unchanged.
2. c_T → ∞ at fixed Δ_g, along a family with a discrete spectrum and finitely many light single traces [A].
3. **(ii′_E)** detector decoupling (§3.3) [A]. This is the only replacement for the bulk-EFT clause. At fixed spin it follows from
   (F_J), standard large-N factorisation on Euclidean compacta. Its uniform-in-spin part is the genuine assumption. It is implied by
   the OPE-data form (ii′) through Lemma 4 (h < 5 extension D-sketch).
4. The inputs of Table `tab:inputs` other than A2: R1, R1′, R2, R3, DS, LRC, H_K, L1–L3 (L4 now optional at k = 2, 3), and for
   τ_\* the angular inputs. All are unchanged.
5. The functional-side step of Theorem A: Lemma 5 (continuation of (fp-target) to light states; [NE] to J = 1000; asymptotics
   D-sketch). τ_\* does not need it (Lemma 6).

   **Net grade of item (1), stated honestly:** [D given A1′, (ii′_E), inputs] for the logic, and [D-sketch] through Lemma 5.
   **Net grade of item (3) (τ_\*):** [D given (ii′_E) and the angular inputs]. This no longer needs the "[D] from planar
   factorisation" of the spin sum.

   **Alternative.** One may state the assumption directly in functional-weighted form: (ii′_Ω) c_T Σ_{O∈𝓛′}|p_Oω[O]| → 0. Then
   item (1) is [P/D given (ii′_Ω)]. (ii′_Ω) is still a pure CFT statement, since Ω^EE is a fixed explicit functional. Theorem A
   shows that (ii′_E), and hence (ii′), suffices for it.

   **Note on the old grade.** The old grade had the same D-sketch hidden "inside A2" (Supp. 685).

For the quartic results, nothing changes: L7(b) [A] is needed for g̃₄ = g₄, and the k = 4, 5 light side needs first-order uniform
control, which is stronger than (ii′).
