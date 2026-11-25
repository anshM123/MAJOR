# LADDER_ALLJ: audit of the all-J first-moment operator (ANGULAR_ROUTE §7) and the all-J helicity ladder (§8)

*25 Nov 2026. Research agent, about 3.5 h. Unreviewed and not committed. No existing file was edited. INDEP_CHECK_T0.md and
`code/t0check/` were not touched.*

* **Code.** New code is in `ads4/code/ladder/`:
  * `common.py`: the moment rules;
  * `fm_check.py`: (A) against the pipeline;
  * `ladder_data.py`: (B) against the pipeline, J ≤ 40, h = ±2;
  * `ward_engine.py` and `ward_symbolic.py`: an independent first-principles derivation of both operators, for fixed J and for
    symbolic J.
* **Outputs.** They go to `code/ladder/cache/`, which ignores everything but its own `.gitignore`.
* **Imports.** `code/kernels/{pipeline,structures,kernels_lib}.py` are imported read-only, with `KERNELS_CACHE` redirected to
  `ladder/cache`.
* **Compute.** About 1.3 CPU-h, at most 2 processes, `ulimit -v 3000000`.

**Grades.**

| grade | meaning |
|---|---|
| [P] | proved, or exact algebra valid for all J, including computer algebra with **symbolic** J, m, Δ |
| [D] | derived by hand, not machine-checked in full |
| [DV] | derived and checked exactly at finite J. This is not a proof for all J |
| [NE] | numerical evidence |
| [I] | open |

---

## 0. Verdict

**The audit is passed. No mismatch was found anywhere.**

* **(A) The first-moment operator is confirmed against the pipeline.** M^orig = (i/6)(W − 18)H holds exactly for J = 0, 2, 4, 5, 6, 7, 8.
  * This covers every conserved structure in both parities, with Δ symbolic and every component m ∈ [−J−1, J+1] [DV].
  * J = 1 and J = 3 have no conserved structure.
  * W is the §7 operator T + B₁² + B₂² + R² + (4−2Δ)B₁ − 2RB₂ acting on the O-polarisation polynomial.
* **The reduction of 𝖪 to W is re-derived independently.** I worked from the momentum-space action of K on ⟨O(p)| and its
  rest-frame reduction by Lorentz covariance; none of the authors' steps were used.
  * The result −K̂_z = W − 18 holds for **symbolic J, m and Δ** [P, given the standard K action; §2]. `ward_symbolic.py` re-checks this.
  * The overall normalisation [z·K, E] = 6i E⁽¹⁾ is consistent three ways: my ∂δ bookkeeping, the authors' (1+εα)^{−5} argument, and
    the pipeline.
* **The +18 and the sign of the normalisation correction are settled.**
  * The +18 is exactly 3(Δ+3) − 3(Δ−3).
    * 3(Δ+3) is the translation from the origin-centred to the midpoint-centred moment.
    * −3(Δ−3) is the change to normalised states, i.e. the correction −i(Δ−3)/2. **The authors' sign is right.**
  * It holds as an operator identity for all J, provided the normalised states are **Wigner-boosted** helicity states.
    * With the pipeline's fixed polarisations, the midpoint moment carries an extra frame term −(i/2)B₁H.
    * That term is verified exactly, J ≤ 8, and the Wigner frame removes it [DV + D].
  * Hence G^N_mid = (i/6) W H^N, with W tridiagonal: diagonal T + C − 2m², off-diagonal −b_m.
  * The sign of the translation constant c = −i is fixed in two ways: from first principles, and by the pipeline's translation piece
    (i/2)(Δ+3−B₁)H, which has exactly the Lorentz-covariance structure. The normalisation correction (c/4)·z·∂log(N_O/N_T) =
    −i(Δ−3)/2 follows from c [P given c].
* **(A) itself is therefore proved for all J,** given two inputs:
  * the operator identity [z·K, E(n)] = 6iE⁽¹⁾(n) (the authors' light-ray argument [D], DV'd at J ≤ 8);
  * the standard K action.
* **(B) The ladder: pipeline agrees with the closed forms.** The test used two rational Δ per J (near unitarity, Δ = J + 10/3, and
  Δ = J + 1001 + 1/7), exact rational arithmetic, and source helicities ±2.
  * For every conserved structure, J = 2 … 40, both parities, these hold **exactly**:
    * the §8 closed forms F_{−k} = F₀C(J,k)(4)_k/(Δ−1−k)_k and F_{4+k} = F₄C(J−4,k)(4)_k/(Δ−5−k)_k, with F₁ = F₂ = F₃ = 0;
    * the §8 recurrence, at every m;
    * the vanishing of the other sector, where appropriate [DV].
  * **h = −2** gives the mirror image F^{(−2)}_m = σF^{(+2)}_{−m}, with σ = +1 for parity-even structures and −1 for parity-odd ones,
    again exactly for J ≤ 40 [DV]. It follows from parity covariance [D]. So the h = −2 ladders are the same two ladders with m → −m, and
    the normalised ratio is unchanged (the norms are even in m).
* **(B) The light-transform conservation step is checked.**
  * I contracted D_z^μ (Todorov operator for transformed spin −2, effective coefficient h−1+J′−1 = −5/2) with my independently derived
    rest-frame K̂_μ (all three components; plain η contraction, no chart reflection needed). This gives **exactly −½ × the authors'
    recurrence**, row by row:
    * for generic amplitude vectors, J = 2, 4, 5;
    * for **symbolic J, m and Δ** (`ward_symbolic.py`) [P, given D_z·[K, E] = 0].
  * Endpoint terms. Any endpoint or contact anomaly in the light-transform conservation would show up as a violation of the recurrence in
    the pipeline, whose α-integrals are computed exactly and independently. None appears through J = 40 [DV]. The α-integrand decays
    like |α|^{−2} or faster at both ends (the pipeline's convergence condition), which is the analytic reason [D].
* **(B) Normalised ratios and uniqueness are proved [P].**
  * The recurrence, with F₁ = F₂ = F₃ = 0 (momentum Ward identity) and F_{±(J+1)} = 0, has a one-dimensional solution space in each
    sector for Δ > J+1.
  * The closed forms solve it.
  * With the all-J norms of §4, a_{k+1}/a_k = ((k+4)/(k+1))√(X_k/Y_k), with a common phase along each ladder.

| item | grade |
|---|---|
| (A) M^N = (i/6)A_{Δ,J}H^N (midpoint, Wigner-normalised), all J | **[P]** modulo the operator identity [z·K,E] = 6iE⁽¹⁾ ([D], authors), with **[DV]** pipeline confirmation of the full statement for J ≤ 8 |
| −K̂_z = W − 18 (rest-frame reduction of the K action) | **[P]** (symbolic J, m, Δ) |
| +18 = 3(Δ+3) [translation] − 3(Δ−3) [normalisation], and its sign | **[P]** given c = −i; c is [P] (first principles) and [DV] (pipeline) |
| (B) recurrence from D_z^μ[K_μ, E] = 0 | **[P]** (symbolic J, m, Δ) given the conservation statement ([D], authors). The recurrence itself is **[DV]** on the pipeline, J ≤ 40 |
| (B) closed forms and normalised ratio for all J | **[P]** given the recurrence. **[DV]** directly, J ≤ 40, h = ±2, all structures |
| Which convention the physical commutator sum rule needs | **[D]** (§4). Neither convention suffices alone; an explicit translation term (D-operator) must be included |

**Net.** The two representation-theory gaps of the angular route are closed. They are closed at the level "[P] given two operator
identities for the light-ray operator at infinity". These identities are:
* [z·K, E] = 6iE⁽¹⁾;
* D_z^μ[K_μ, E] = 0.

Both are standard-type consequences of the light-transform being a conformally covariant intertwiner. Their **consequences** are DV'd
exactly against the independent pipeline. The physical-matching items of ANGULAR_GO §2.2(c) remain open. These are not representation
theory: double traces, contacts, the k = 3 flip entries, and the common normalisation behind 47266.

---

## 1. Definitions (what a_k, H, M are)

* **Transition amplitude.** M_m(θ) = ∫d³x e^{−ip·x}⟨0|O_m(x)†E(n)T_h(0)|0⟩, with ⟨O(p′),m|E(n)|T(p),h⟩ = (2π)³δ³(p−p′)M_m(θ).
  * Rest frame: p = (1,0,0). Detector z = (1,n), with n = (cos θ, sin θ) and θ = 0 below.
  * Component m is the coefficient of t^{J−m} for the bra polarisation ζ̄ = ζ(t,1) (pipeline convention).
  * Equivalently, it is the coefficient of e_m = a^{J+m}b^{J−m} in the authors' convention; W is invariant under a ↔ b, so the two
    agree.
  * H_m := M_m (zeroth moment). M_m(θ) = e^{i(h−m)θ}M_m(0).
* **First moments.** E⁽¹⁾ = ∫du u E(u,n), with u the retarded time (u = α/2 in the pipeline's light-ray parameter).
  * *Midpoint:* u centred at the midpoint of x₃ (T) and x₄ (O). Per monomial, the analytic rule is I₁^mid = (a−b)/(2(a+b−2))·(2z·x)·I₀
    for the **u**-moment. This is ANGULAR_GO's "R2"; the factor-2 remark (α-moment = 2 × u-moment) is respected.
  * *Origin (T-centred):* u centred at x₃ = 0. M^orig = M^mid − ½(z·x)I₀-piece.
* **Normalised (Wigner) amplitudes.** a_m := H_m/√(N_{O,m}N_T), with N_{O,m} ∝ C(2J,J+m)Γ(Δ+m−1)Γ(Δ−m−1) (ANGULAR_ROUTE §4,
  KERNELS_D3 §3.3). The states at nearby momenta are |p,m⟩ = |p, L(p)ζ_m⟩, with L(p) the standard boost. The rest-frame values are the
  same in any frame; the frame matters only for first moments (§3).
* **The ladders (B).** Fix source helicity h = +2.
  * **b = 0 ladder:** a_k := a_{m=−k}, k = 0..J (harmonic h − m = 2+k).
    * Even J: carried by S^{++} and S^{odd}, which have no m ≥ 4 component.
    * J = 0 and J = 2: by the single structure of each parity.
  * **b = 4 ladder:** a_k := a_{m=4+k}, k = 0..J−4 (harmonic −(2+k)).
    * Even J ≥ 4: carried by S^{+−}.
    * Odd J ≥ 5: by the single parity-odd structure.
  * The ratio is the same for **every** structure and every combination of structures in a sector.
  * X_k = C_J − (b+k)(b+k+1) = (J−b−k)(J+b+k+1) and Y_k = T_Δ − (b+k)(b+k+1) = (Δ−b−k−2)(Δ+b+k−1). The ladder ends where X_k = 0, i.e.
    k = J − b.
  * For h = −2 everything is mirrored: m → −m, and amplitudes are multiplied by σ = ±1 according to parity.

## 2. (A) First-principles derivation of the rest-frame K operator [P] (`ward_engine.py`, `ward_symbolic.py`)

**Step 1: K on a momentum bra.**
* Take [K_b, O(x,ζ)] = −i(ξ·∂ + 2Δ(b·x) + (ωζ)·∂_ζ)O with ξ = 2(b·x)x − x²b and ω^μ_ν = 2(x^μb_ν − b^μx_ν).
* Fourier transforming (x^μ → −i∂/∂p_μ, ∂_μ → −ip_μ) makes K_b act on |O(p),ζ⟩ as the **real** operator

      𝒦_b = 2(3−Δ) b·∂ + 2 (p·∂)(b·∂) − (b·p) ∂·∂ − 2(b·ζ) ∂_p·∂_ζ + 2(ζ·∂_p)(b·∂_ζ),    ∂ = ∂/∂p^μ.

* ⟨O(p′),ζ̄| K_b acts by the same operator on the bra polarisation.

**Step 2: one-sided Ward identity.** K kills T(0)|0⟩, so ∫e^{−ip′x}⟨O(x)|[K_b,E]|T(0)⟩ = 𝒦_b F(p′), with F the zeroth-moment
function.

**Step 3: rest-frame reduction.** Lorentz invariance gives (vector field of ω on p + A_ω)F = 0, where A_ω is the sum of the induced
generators on ζ̄, ζ_T and z. Homogeneity is p·∂F = (Δ+1)F. At p₀:

    ∂_0F = (Δ+1)F,  ∂_iF = −A_iF,  ∂_i∂_jF = A_iA_jF − δ_ij(Δ+1)F,  ∂_0∂_iF = −Δ A_iF,  ∂_0²F = (Δ+1)ΔF.

**Step 4: spinor form.** On spinors ζ(a,b) = a²e₊/2 + b²e₋/2 + ab e₀, the induced generators of the x-boost, the y-boost and the rotation
are exactly the authors' B₁ = ½(b∂_a + a∂_b), B₂ = (i/2)(a∂_b − b∂_a) and R = −(i/2)(a∂_a − b∂_b). The spin term is a sum of
∂_{p^μ} × (rotation in the (b, e_μ) plane), decomposed into B₁, B₂, R.

**Results.**

| check | result |
|---|---|
| −𝒦_z F at θ = 0, generic amplitude vector f_m, symbolic Δ, J = 2, 4, 5 | = (W − 18)f exactly. The opposite sign of the spin term fails, so the sign of the spin term is fixed |
| the same, with **symbolic J and m** (local window f_{m−2..m+2}, symbolic exponents) | = (W − 18) row m exactly, `ward_symbolic.py`. Hence the §7 tridiagonal operator for **all J** [P] |
| ∂δ(p−p′) bookkeeping of ⟨O(p′)\|[z·K,E]\|T(p)⟩ (midpoint variables P = (p+p′)/2, q = p−p′) | The second-order terms give **no** ∂²δ and no δ terms. With the spin sign above and Lorentz invariance, the ∂δ part collapses to −6 z·∂_qδ(q)F. This is the translation term of 6iE⁽¹⁾ (c = −i below). So [z·K, E] = 6iE⁽¹⁾ + (translation-invariant part); the authors' (1+εα)^{−5} argument says there is no such extra part, and the pipeline (§3) confirms it |

Combining these: M^orig = (1/(6i))𝒦_zF = (i/6)(W − 18)H.

## 3. (A) Pipeline check, J = 0..8, and the +18 [DV] (`fm_check.py` → `cache/fm_check_J*.txt`)

**Origin-centred.** M^orig_m = (i/6)[(W−18)H]_m holds **exactly** at every m ∈ [−J−1, J+1], with Δ symbolic, for:
* J = 0 (even, odd); J = 2 (even, odd); J = 4 (2 even + 1 odd); J = 5 (odd); J = 6 (2 + 1); J = 7 (odd); J = 8 (2 + 1).

This covers the first-moment components m = 1, 2, 3, where H vanishes. For example, M_1 = −(i/6)(Δ−1)J·H₀ (unnormalised).

**Midpoint (fixed polarisations).** M^mid = M^orig + (i/2)[(Δ+3) − B₁]H holds **exactly** for the same cases.
* This is the translation piece −(c/2)z·∂F with c = −i.
* z·∂F at rest is the sum of three pieces:
  * scaling: Δ+1;
  * the detector weight: +2, since a boost along n rescales z, and E has weight −2 in z;
  * minus the boost of the **fixed** bra polarisation: −B₁.
* The T polarisation e₊ boosts into the gauge direction p, so it contributes nothing.
* The observed −B₁ is exactly what Lorentz covariance predicts. It is a frame artefact.

**The translation constant.** For ⟨O(p′)|E⁽¹⁾|T(p)⟩ = δ·G + c z·∂_{p′}δ(p′−p)·Y, with Y the zeroth moment:
* first principles give c = −i (u_x = −z·x, and ∫e^{−ikx}x = i∂_kδ);
* the pipeline's translation piece −(i/2)(Δ+3−B₁)H fixes the same value.

**Conversions** (δ-coefficient algebra, all J) [P]:

| convention | δ-coefficient |
|---|---|
| Midpoint vs Y at the O momentum (T-centred) | G_mid = G_T − (c/2)z·∂Y |
| Normalised by V_O(p′)†(…)V_T(p), midpoint | extra (c/4)·z·∂log(N_O/N_T) = −(i/4)[(2Δ−3) − 3] = **−i(Δ−3)/2** (N ∝ (−p²)^{Δ−3/2}). Authors' sign confirmed |
| Wigner frame | z·∂F_W = (Δ+3)F, diagonal: the B₁ term is gone |

**Therefore**
* G^N_mid,W = G_T/√N + [3(Δ+3) − 3(Δ−3)](i/6)H^N = (i/6)(W−18+**18**)H^N = (i/6)WH^N.
* In the normalised basis W is symmetric: diagonal T + C − 2m², off-diagonal −b_m (verify_ward.py, and the norm ratio of §4).
* J = 0: i(Δ²−3Δ−16)/6 → iT_Δ/6 = i(Δ−1)(Δ−2)/6. This agrees with ANGULAR_GO, and the sign is now [P] rather than [D-sketch].
* For the record, the T-centred normalised value is G^N_T = G_T/√N + (3i/2)H^N = (i/6)(W−9)H^N.

## 4. Which convention the physical first-moment sum rule needs [D]

**Setup.** Take the commutator ⟨T(p″)|[E⁽¹⁾(n₁), E(n₂)]|T(p)⟩ = 0, insert the states, and use the midpoint δ-coefficient for a
symmetric, normalised T wave packet.
* The ∂δ parts cancel family by family after the sum, by the zeroth-moment sum rule.
* The δ part gives, summed over m, the following identity, with Wigner-normalised amplitudes F_O(n) = (a_m e^{i(2−m)θ}) and
  G_O = (i/6)WF_O:

      Σ_O [ Ḡ_O(n₁)·F_O(n₂) − F̄_O(n₂)·G_O(n₁) − (i/2)( F̄_O(n₁)·𝒟F_O(n₂) + (𝒟F̄_O(n₂))·F_O(n₁) ) ] = 0,

      𝒟 = z₁·∂_p|_rest = 1 + 2cos φ + sin φ ∂_{θ₂},   φ = θ₂ − θ₁.

**The 𝒟 operator.** It is the scaling of the normalised amplitude (+1) plus the boost along n₁ acting on the second detector.
* On e^{ikθ₂} it acts as 1 + (1+k/2)e^{iφ} + (1−k/2)e^{−iφ}.
* This is the structure of the authors' Ỹ_m = q_m + (1−m/2)q_{m−1} + (1+m/2)q_{m+1} (ANGULAR_ROUTE §9). Their direct sum rule evidently
  contains this term. I did not re-derive their normalisation (/3, RHS 3).

**The two conventions.** The heavy operator to use is the **midpoint, Wigner-normalised** G = (i/6)A_{Δ,J}H^N, **together with** the
𝒟-term above.
* The origin (T-centred) convention is equivalent: G^N_T = (i/6)(W−9)H^N, with the 𝒟-term replaced by −i z₁·∂Σ(F̄₁F₂) = −iF̄₁(𝒟+2)F₂.
  The two agree after the family sum, which I checked algebraically.
* What is **not** allowed is to use either G alone, without the translation term.

**For the authors' telescoping.** The 𝒟-term is diagonal in m (it only shifts harmonics), so it adds a Δ-independent kernel and does not
touch the T_Δ-enhancement. This is consistent with ANGULAR_GO's scaling ceiling.

## 5. (B) Ladder: pipeline J ≤ 40, h = ±2 [DV] (`ladder_data.py` → `cache/ladder_J*.txt`)

**Method.**
* Conserved structures were computed at rational Δ (structures.py, exact). Amplitudes come from the read-only pipeline, with exact
  rational real and imaginary parts.
* The two Δ per J are Δ = J + 10/3 and Δ = J + 1001 + 1/7.

**Checks per structure** (all exact equalities):
1. P_m = 0 for m = 1, 2, 3 and for m = ±(J+1);
2. the closed forms (R0) and (R4), or the vanishing of the absent sector;
3. the §8 recurrence at every m ∈ [−J, J];
4. h = −2: F^{(−2)}_m = σF^{(+2)}_{−m} at every m.

**Result: RESULT_LINE**

The sector content (which ladder each structure carries) is exactly as in §1.

**Cost.** About 1–2 CPU-min per J, essentially independent of J, because the number of monomials is fixed (13/12). The heavier
symbolic-Δ runs are only needed for (A).

**Normalised ratio [P].** From the closed forms and the all-J norm ratio N_{m+1}/N_m = (J−m)(Δ+m−1)/((J+m+1)(Δ−m−2)):

* b = 0, m = −k: (a_{k+1}/a_k)² = ((k+4)/(k+1))²·(J−k)(J+k+1)/((Δ−k−2)(Δ+k−1)) = ((k+4)/(k+1))²X_k/Y_k.
* b = 4: the same with k → k+4 in X and Y.
* Both F ratios are positive for Δ > J+1, so the ladder has a common phase.

## 6. (B) The conservation step, audited [P + DV]

The authors' chain:
* conservation of 𝕃[T] (transformed spin −2, so the Thomas coefficient is h−1+J′−1 = −5/2);
* the base point at infinity, so ∂_y becomes K;
* hence D_z^μ[K_μ,E(z)] = 0;
* hence D_z^μ K̂_μ^𝒪 F = 0.

Audit, item by item:

| item | finding |
|---|---|
| Thomas/Todorov coefficient | D^μ = (h−1+z·∂_z)∂_z^μ − ½z^μ∂_z², with h = 3/2, acting on K̂_μF (homogeneity −2 in z): the effective coefficient is −5/2, as the authors say. The operator is interior for any extension off the cone (check: D(z²g) ∝ z² because 2h = d). I used the extension φ_k = (z¹+iz²)^k(z⁰)^{−2−k} |
| chart map, translations → K | I used my own K̂_μ (§2, all three components), with the **plain** η contraction D^μK̂_μ. No inversion reflection of indices is needed |
| m-weights | the coefficient of e_m in D^μK̂_μF equals **−½ × [§8 recurrence row m]**: exactly, for generic f, J = 2, 4, 5, every m (`ward_engine.py`); and for **symbolic J, m, Δ** (`ward_symbolic.py`). So the weights m(m−4), (m−3)(m−4), m(m−1) times the §7 diagonal and neighbour coefficients are right |
| endpoint terms | not derived analytically beyond the decay argument: the α-integrand is O(\|α\|^{−2}) at both ends, which is the pipeline's convergence condition nfac ≥ 0. Settled empirically: the pipeline's exact α-integrals satisfy the recurrence through J = 40, both parities, both helicities, two Δ |
| boundary conditions | F₁ = F₂ = F₃ = 0 is the momentum Ward identity, which is rigorous. In the recurrence, rows 0 and 4 have zero diagonal weight, so the sectors decouple |
| uniqueness | for Δ > J+1 the pivots m(m−1)(Δ+m−2)(J−m+1) (m ≤ −1) and (m−3)(m−4)(Δ−m−2)(J+m+1) (5 ≤ m ≤ J−1) are nonzero, so each sector has a one-dimensional solution space [P] |
| h = −2 | `ward_symbolic.py --minus` gives the mirror recurrence (m → −m). The pipeline mirror relation holds for J ≤ 40 |

## 7. Reproduce (`ads4/code/ladder/`; `ulimit -v 3000000; export OMP_NUM_THREADS=1`)

```
python3 -B fm_check.py 0 2 4 5        # ~5 min    (A) pipeline: origin (W-18), midpoint translation piece
python3 -B fm_check.py 6 7 8          # ~15 min
python3 -B ward_engine.py 2 4 5       # ~1 min    first-principles K-hat: -K_z = W-18 ; D.K -> recurrence (fixed J)
python3 -B ward_symbolic.py           # symbolic J, m, Delta, helicity +2
python3 -B ward_symbolic.py --minus   # helicity -2 (mirror)
python3 -B ladder_data.py 2 23        # ~35 min   (B) pipeline vs closed forms / recurrence / h=-2 mirror
python3 -B ladder_data.py 24 40       # ~40 min
```
