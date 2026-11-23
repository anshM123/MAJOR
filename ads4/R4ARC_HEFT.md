# L6′ and H_EFT(a): the collinear integrability of the R⁴ light function, and term-wise EFT matching for heavy exchanges on spinning externals

*23 Nov 2026, R4-arc/H_EFT agent. Unreviewed, not committed. No existing file was edited.*
* *New code:* `ads4/code/r4heft/` (`rh_common.py`, `rh_explore.py`, `rh_degree.py`, `rh_lemma_count.py`, `rh_t1_collinear.py`,
  `rh_t1_frozen.py`, `rh_t1_laurent.py`, `rh_t1_tail.py`, `rh_t2_heft.py`, `rh_t2_report.py`, `run_all.sh`). Outputs in `ads4/code/r4heft/cache/` (git-ignored except its `.gitignore`).
* *Read-only imports:* `resultq/lsd_core.py`, `lsd_mellin.py`, `lsd_pos.py`, `og_free.py` (LIGHTSIDE_DIST's Mellin machinery on the
  product family); `exactskel/exactmom.py` (exact d = 3 Mack polynomials and weights).
* *Literature read from source:* CLS-D 2311.04271 (§2.1 "dispersive functional", footnote on subtractions; §4.2 on subtracted
  spinning sum rules), CHMRS 2106.10274 (§4 "Light contributions to holographic sum rules", §5 "technical assumption", §6.4
  "Comments on truncation errors", App. "Computing ω|_light"); Paulos 1107.1504 (§2 and §7: spin-1/spin-2 bulk-to-boundary
  propagators as differential operators on scalar ones; abstract and source checked); Costa–Gonçalves–Penedones 1404.5625 (abstract
  checked). Penedones 1011.1485 (abstract checked). The last three are **not** in `prl/refs_verified.bib`.
* *Compute:* about 2 CPU-h, never more than 2 processes, one thread each, < 1 GB each.

Grades: **[P]** proved; **[D]** derived; **[D-sketch]**; **[DV]** derived and checked numerically; **[Q]** quoted; **[NE]**; **[M]**;
**[A]** assumption; **[I]** open. A chain carries the grade of its weakest link.

---

## 0. Summary

| item | result | grade |
|---|---|---|
| **L6′** (a = b integrability of the borderline R⁴ light function against `C_ν`) | **Closed.** In the og/CLS-D gauge the R⁴ light function obeys `|L_{R⁴}(a,b)| ≤ C Δφ^{−2}(1 + |log Δφ|)` uniformly in the mean angle; the criterion needs `p < 7/2`, so `p = 2` holds with **margin 3/2**. Free-standing: it uses neither H_EFT nor LRC\*. The same bound holds for every contact term (any number of derivatives, parity even or odd) and for the double-trace part of every heavy exchange. | **[I] → [D]** given two standard representation facts [Q]; counting lemma [P] (exhaustive check); mechanism [DV] on the product family |
| Lemma Y (counting lemma behind L6′) | every term `monomial × [one ε] × D_{Δ′}` of a spin-2 contact diagram is `O(y^0)` (up to logs) in the light-ray configuration, `y = 2 sin²(Δφ/2)` | **[P]**; exhaustive enumeration of 997 836 structures, worst margin `0` |
| Twist rule (dispersive terms) | for `f_k ∝ y^{−1}` (k = 2, 4) the light function of a **dispersive** tree term X is `y^{τ_min/2 − 2}(A + B log y)`, `τ_min` = lowest twist in X's (13) channel: `p = 4 − τ_min` (consistent with the free-field `Δφ^{−3}` law, τ = 1, of RESULT_Q_OG_CHECK §3.6) | [D]; [DV] at frozen angles: `M(8y) = A + B log 8y` to `3·10⁻¹³` over `y = 10⁻³…10⁻⁹`, remainder `∝ y^{≥0.4}` (k = 2 and 4) |
| Per-term dispersiveness in the og gauge | for any term whose s-channel Mellin poles are at integers, only finitely many residues `δ₁₂ = −n` are non-zero (the residue integrand is analytic inside the contour disc) | [D]; [DV] (residues `≤ 10⁻¹⁶` beyond `n = 4/5` in classes s, t; `n = 0` in 1324; all zero in u) |
| **Non-dispersive borderline terms (why dispersiveness is the crux)** | in the toy's 1243 class (half-integer free-field (12)-links) the residue tail `r_n ∝ n^{−p(δ₁₃)}` has `p` falling through 1 near `Re δ₁₃ ≈ 0.4`, which creates a δ₁₃-singularity in the strip and **extra collinear growth**: local slope 4.7 → 3.5 at `Δφ = 1.2 → 0.04`, above the twist-rule value 3, i.e. at the edge of the `p < 7/2` criterion. LEMMAS's concern was therefore real for non-holographic terms; it is the dispersiveness of holographic terms that removes it | [DV]/[NE] (exponent not converged) |
| **L7(a) = H_EFT(a), spinning** | **Holds pointwise in celestial angle, hence ν-paired, for every dispersive term** (all tree terms of a holographic ⟨TTTT⟩). The term-wise EFT expansion of a heavy exchange's double-trace part converges to its actual (OPE-side) value as an asymptotic series in `1/m²`. | **[D-sketch] → [D]** given the same [Q] inputs plus CLS-D's dispersiveness [Q]; **[DV]** on the dispersive spinning sector of the product family (heavy exchanges `J_h = 0, 2, 4`, `τ_h = 41, 61, 81`, k = 2 and 4): `|TP_R − DT|/|DT|` falls geometrically to `10⁻¹⁸–10⁻³⁷` |
| F_NORM's "H_EFT(a) is false pointwise in angle [NE]" | **Re-read:** the failure is a property of the non-dispersive classes of the toy (twist-1 free-field content in the (12) channel, not double traces of T′). In the dispersive classes of the same family the pointwise statement holds to `10⁻¹⁸–10⁻³⁷`; in the non-dispersive pair the term-wise series diverges (`|TP_R − DT|/|DT|` up to `10²³`), reproducing B3 | [DV] |
| Leading-order calibration (route II step 5) in the toy, scalar heavy exchange | `DT/EFT_lead − 1 = +0.123, +0.026, +0.0086` (k = 2) and `+0.31, +0.075, +0.031` (k = 4) at `τ_h = 41, 61, 81`: → 0, locally like `τ_h^{−3…−3.6}` | [DV] (dispersive pair) |
| same, spinning heavy exchange | `J_h = 4`: `−0.24, −0.16, −0.10` (k = 2), `−0.26, −0.20, −0.14` (k = 4), falling; `J_h = 2`: `≈ −0.03 ± 0.01` at both k, not yet falling at `τ_h ≤ 81` | [NE]; not used (route II needs only a scalar χ) |
| H_EFT(b) | unchanged **[A]**; §4 gives the literature's statement and how our use compares | [A] |
| **Result Q, upper end** | the one [I] link (L6′) is gone. The implication "(Pos) ⇒ (2)" is now **[D-sketch]** (weakest: H_K through the EE step, and `g₅`-(P)), plus **[A]** for H_EFT(b) in (2′) | see §5 |

**Residual (stated exactly).** Nothing in this note evaluates the actual AdS₄ R⁴ (or χ-exchange) ⟨TTTT⟩. The derivations use two
quoted representation facts: (R1) every tree-level spin-2 contact diagram is a finite sum `Σ (monomial in Z_i·Z_j, Z_i·P_j, P_i·P_j)
[× one ε] × D_{Δ′}(P)` with `Δ′_i = 3 + deg_{P_i}` (Paulos 1107.1504 §2, §7: spinning bulk-to-boundary propagators as polynomial
differential operators on scalar ones; CGP 1404.5625), and (R2) the Symanzik/Mellin representation of `D_{Δ′}` with poles only from
`Γ(δ_ij)` [Q]. For H_EFT(a) also (R3) the spinning exchange diagram has a Mellin representation with polynomial residues and integer-
spaced double-trace poles [Q, same sources], and (R4) CLS-D's dispersiveness of `Ψ_k` (double zeros on all but finitely many double
twists) [Q]. The numerical checks exercise every step of the mechanism on spin-2 externals (`T′ = σ ⊗ T`), but not the graviton's own
tensor structures.

---

## 1. Setting

* **Gauge** (og_free / EIN / LEMMAS §1.1): `x₀ = 0`, `x₂ = −e`, `x₄ = e`, `x₃ = t₃z(b)`, `x₁⁺ = −t₁z(a)`, light-ray polarisations
  `z(a)`, `z(b)` at 1, 3; measure `t₁^w t₃^w` (`w = 3` for T, 5 for `T′ = σ ⊗ T`); subtraction `f_k(u′, v′)`.
  `X₁₂ = 2t₁ − 1`, `X₃₄ = 2t₃ − 1`, `X₁₃ = 2y t₁t₃`, `X₁₄ = 2t₁ + 1`, `X₂₃ = 2t₃ + 1`, `X₂₄ = 4`, `y = 1 − cos Δφ = 2 sin²(Δφ/2)`.
* **Light function** of a tree term X: `L_X(a,b) := Ω_k[X](a,b) = E⁺(a,b) + (−1)^k E⁺(b,a)`, the term-wise s-channel swapped
  functional at fixed angles (LIGHTSIDE_DIST §1.1). For a borderline X it equals the Regge arc `𝒜_k[X](a,b)` (LSD §1.2).
* **Criterion** (LEMMAS §4.1, [D]): `|C_ν| ≤ A₁₃^{5/4} K̄(Δφ)`, `K̄ = 2 log(1/Δφ) + 8.69 + o(1)` uniformly in ν, so `L` pairs absolutely
  with `C_ν` (uniformly in ν, then against any `w_θ ∈ L¹`) iff `|L| ≲ Δφ^{−p}` with `p < 7/2`.
* **Where y enters [P].** In this gauge `y` appears only through `X₁₃ = 2yt₁t₃` and through contractions linking the 1- and 3-legs:
  `Z₁·Z₃ = z(a)·z(b) = −y`, `Z₁·P₃ = −z(a)·x₁₃ = −t₃ y`, `Z₃·P₁ = t₁ y`, `P₁·P₃ = −½X₁₃ = −yt₁t₃` (embedding conventions
  `P·Q = −½x_{PQ}²`, `Z_i·P_j = −z_i·x_{ij}`; `z_i·x_i = 0` because `x_i ∥ z_i`). All other invariants are `y`-independent (smooth in
  the mean angle). The subtraction factors: `f₂ = v′ − u′ = y^{−1}F₂(t)`, `f₄ = y^{−1}F₄(t) + yF₄′(t)`, `f₃ = −1`, `f₅ = O(y⁰)`
  (`u′, v′ ∝ 1/s₁₃ = −1/(2yt₁t₃)`). **k = 2 and 4 carry an explicit `y^{−1}`; k = 3, 5 none.**

---

## 2. Target 1: L6′

### 2.1 Structure of the light function in Mellin space [D]
Write a tree term as a finite sum of terms `S_α(Z,P) × M_α`, where `S_α` is a monomial in the invariants (times at most one ε, §2.3)
and `M_α` a scalar Mellin integral `∫[dδ] ∏Γ(δ_ij) X_ij^{−δ_ij} 𝓜_α(δ)` (for a contact: `𝓜_α = const`, i.e. a D-function `D_{Δ′}`).
At fixed angles, the term-wise s-channel functional closes the `δ₁₂` contour to the left (LSD §1.2). Because `y` enters only as in §1,

    L_X(y) = Σ_α y^{m_α} ∫_{Re δ₁₃ = c_α} dδ₁₃/(2πi) (8y)^{−δ₁₃} Φ_α(δ₁₃),     Φ_α(δ₁₃) = Γ(δ₁₃)Γ(δ₂₄) Σ_n Res_{δ₁₂ = pole_n}[…]         (2.1)

with `Φ_α` independent of y, `m_α` = (y-order of `S_α` and of `f_k`), and the mean angle entering `Φ_α` smoothly. Its small-y behaviour
is read off by moving the `δ₁₃` contour to the left: each pole `δ₁₃ = δ*` of `Φ_α` gives `y^{m_α − δ*}(log y)^j`.

**Which poles can occur.** (i) `Γ(δ₁₃)Γ(δ₂₄)` with `δ₂₄ = δ₁₃ − σ_α`, `σ_α = (Δ′₁ + Δ′₃ − Δ′₂ − Δ′₄)/2`: poles at `δ₁₃ = −j` and
`σ_α − j` — the (13)-channel double-trace families. (ii) After `δ₁₂ = −n`: `Γ(δ₁₄)Γ(δ₂₃)` at `δ₁₃ = Δ′₁ + n + j`, … — all to the right.
(iii) The polynomial `𝓜_α` (contact) or the heavy-exchange residues: no poles near the contour; a u-channel heavy pole sits at
`δ₁₃ ≈ 2 − Δ_χ/2 − m`, far left, giving `y^{Δ_χ/2 + …}`. (iv) The kernel
`K_α(δ₁₂, δ₁₃) = ∮_C∮_C t₁^w t₃^w f_k (polarisation factors)(U/V)^{−δ₁₂}V^{δ₁₃}`: its δ₁₃-singularities come only from the endpoint
`t → 0`, where the integrand is `t^{w−1−δ₁₃+a}` (`a ≥ 0` from polarisation factors, §2.4), i.e. at `Re δ₁₃ ≥ w` — to the right.
**Hence, as long as the sum over residues is finite (§2.2), the only singularities met when moving left are the (13)-channel
double-trace poles (i).** For an infinite residue sum (a non-dispersive term) the large-n tail is the joint Regge corner, its decay
exponent depends on `Re δ₁₃`, and the sum itself can be singular inside the strip (§2.6 (d)): this is the precise sense in which the
a = b question "is the same corner as LRC\*" (LEMMAS §4.1) — and why it disappears for dispersive terms.

### 2.2 The representation [Q] and per-term dispersiveness [D]
* **(R1) [Q].** A tree contact diagram of spin-2 fields is a finite sum of `monomial × [ε] × D_{Δ′}`: Paulos 1107.1504 writes the
  spin-1 and spin-2 bulk-to-boundary propagators as `D^{MA}_Δ = η^{MA} + Δ^{−1}P^A∂_M` (and its spin-2 analogue) acting on the scalar
  one; bulk covariant derivatives and metrics produce polynomials in `X`, `P_i·X`, `Z_i·X`, which reduce to `∂_{P}` acting on
  `∏(−2P_i·X)^{−a_i}`; and `∂/∂X_ij` raises indices, `I_{Δ+e_i+e_j} = −∂I_Δ/∂X_ij` (EIN §2.2). No denominators arise. Homogeneity
  fixes `Δ′_i = Δ_i + deg_{P_i}(monomial)`. Parity-odd vertices give one ε of five embedding vectors.
* **(R2) [Q].** `D_{Δ′} ∝ ∫[dδ]∏Γ(δ_ij)X_ij^{−δ_ij}`, `Σ_j δ_ij = Δ′_i` (Symanzik; Penedones 1011.1485).
* **Integrality [P].** Each monomial has `Σ_i deg Z_i = 8`, so the number of `Z·P` factors is even and the total P-degree is even;
  hence `σ_α` and the s-channel shift `σ′_α = (Δ′₁ + Δ′₂ − Δ′₃ − Δ′₄)/2` are integers (and half-integers only in the ε-terms through
  `Δφ`, §2.3).
* **Per-term dispersiveness [D].** At `δ₁₂ = −n` (and `σ′ − n`) the residue integrand contains `X₁₂^{n+c}X₃₄^{n+c′}` with integers
  `c, c′`, and every other factor is analytic in the closed disc bounded by `C` except at `t = 0`, where the cut of `V^{δ₁₃}` points
  outward. For `n ≥ −min(c, c′)` the integrand is analytic inside the disc and integrable at `t = 0`, so `∮_C = 0` by Cauchy. **Only
  finitely many residues are non-zero.** This is the og-gauge form of CLS-D's "double zeros on all but finitely many double twists"
  (2311.04271 §2.1 and its footnote on subtractions). [DV]: in the product family the residues vanish to `10⁻¹⁶–10⁻²⁰` beyond `n = 4`
  (class s, k = 2), `n = 5` (k = 4), at all `n ≥ 1` for 1324 and at all n for u (`cache/explore_k*.txt`); the classes with a single
  free-field (12)-link (`(2t₁−1)^{−1/2}`) are not dispersive (1243: residues `∝ n^{−2.2}`).

So for a contact term the light function is a **finite** sum of one-dimensional Mellin–Barnes integrals (2.1): no arc at infinity and no
`n`-uniformity question arises. (For the R⁴ term this is the same value as LSD's arc, by Cauchy.)

### 2.3 Lemma Y (counting) [P]
**Lemma Y.** For every term `monomial × [ε] × D_{Δ′}` of a contact diagram of four spin-2 operators of dimension 3, in the light-ray
configuration of §1,

    max(0, σ_α) ≤ m_α^{pol},

where `m_α^{pol}` = number of factors among `{Z₁·Z₃, Z₁·P₃, Z₃·P₁, P₁·P₃}` (each `O(y)`), plus `½` if the ε contains three of
`{P₁, Z₁, P₃, Z₃}`. Hence each term of the correlator, and each of its residues, is `O(y⁰ |log y|)` as `y → 0`.

*Proof.* Count contributions to `2σ = deg P₁ + deg P₃ − deg P₂ − deg P₄`. Let `a` (`a′`) be the number of `Z₂`/`Z₄` legs ending on
`P₁/P₃` (`P₂/P₄`), `b` the number of `Z₁/Z₃` legs on `P₂/P₄`, `e` the number of cross factors `Z₁·P₃, Z₃·P₁`, `g` of `Z₁·Z₃`, `c` of
`Z_{1,3}·Z_{2,4}`, `z` of `Z₂·Z₄`, `f` (`f′`) of `P₁·P₃` (`P₂·P₄`), and `k₁₃` (`k₂₄`) the number of ε-slots taken by `{P₁,Z₁,P₃,Z₃}`
(`{P₂,Z₂,P₄,Z₄}`), split as `e₁₃ + p₁₃` (Z, P) and `e₂₄ + p₂₄`. Leg counting (`J = 2` for each operator):
`a + a′ + c + 2z + e₂₄ = 4`, `b + c + 2g + e + e₁₃ = 4`. Then

    2σ = a − a′ − b + e + 2f − 2f′ + p₁₃ − p₂₄ = 2(g + e + f) − 2a′ − 2z − 2f′ + (e₁₃ + p₁₃) − (e₂₄ + p₂₄) − e₂₄ − … ≤ 2(g + e + f) + k₁₃ − k₂₄.

Without ε: `σ ≤ g + e + f = m^{pol}`. With one ε, `k₁₃ + k₂₄ = 5`: `k₁₃ ≤ 2` gives `σ − m ≤ −½`; `k₁₃ = 3` gives `σ ≤ m_dots + ½`, and
such an ε vanishes like `z(a) ∧ z(b) ∝ Δφ ∝ y^{1/2}` (the four vectors `P₁, Z₁, P₃, Z₃` span three dimensions that degenerate to two as
`z(b) → z(a)`; `P₁ − P₃ = −t₁Z₁ − t₃Z₃`); `k₁₃ = 4` gives `ε ≡ 0`. The leftmost `δ₁₃` pole is at `max(0, σ)` (§2.1 (i)), and the term
is `O(y^{m − max(0,σ)} log y)`. ∎

*Check* (`rh_lemma_count.py` → `cache/lemmaY.txt`): exhaustive enumeration of all Z-leg assignments, all `P·P` factors up to
multiplicity 2, parity-even and every one-ε placement: **997 836 structures, `max[max(0,σ) − m] = 0`**. The y-orders of the ε's were
measured numerically in the actual configuration: `k₁₃ = 3`: `0.500` in all 83 samples; `k₁₃ = 4`: exactly zero up to roundoff;
`k₁₃ ≤ 2`: O(1) (orders between −0.12 and +0.5, i.e. bounded; `{Z₁, Z₃}` pairs vanish like `y^{1/2}`).

This is the same count as the (13) light-cone OPE: an operator of twist τ in `T₁ × T₃` contributes `y^{τ/2 − 1}` at least (each Z not
contracted into the exchanged operator costs `y^{1/2}`, each Z contracted into it costs `y^{1/2}` through the `x₁₃²` power), and a
contact has only `[TT]` double traces (twist ≥ 2 in d = 3) in every channel (HPPS).

### 2.4 Endpoint convergence on the original contours [D]
At the single endpoint `t₁ → 0` the residue integrand is `t₁^{w − 1 − δ₁₃ + #(Z₃·P₁) + f}` (from `t₁³`, `f₄ ∝ t₁^{−1}`, `X₁₃^{−δ₁₃}`,
and the cross factors `Z₃·P₁ ∝ t₁`, `P₁·P₃ ∝ t₁`). On the term's own contour `Re δ₁₃ ↓ max(0,σ)`, the exponent is
`≥ 2 − #(Z₁·P₃) − g ≥ 0` by the leg count of `Z₁`, so the t-integrals converge on every term separately (t₃ likewise; the joint corner
factorises). Moving left only improves this. So (2.1) is legitimate term by term.

### 2.5 Proposition (L6′) [D]
> **Proposition.** Let `L_{R⁴}(a,b)` be the light function of the AdS₄ R⁴ contact ⟨TTTT⟩ (any of the three real 8-derivative
> couplings, parity even or odd) under `Ψ₄` (the j′ = 0 functionals of `Ĉ₄^{ds}`, `∂Ĉ₄^{ds}`). Then
>
>     |L_{R⁴}(a,b)| ≤ C (1 + |log Δφ|) Δφ^{−2}     uniformly in the mean angle,
>
> so `∫dφ₁dφ₃ |C_ν L_{R⁴}| ≤ C′ ∫ dΔφ Δφ^{1/2} log²(1/Δφ) < ∞` uniformly in ν: **p = 2 < 7/2, margin 3/2.** The same holds for every
> contact term (all derivative orders, so all of the EFT tower), for the double-trace part of every heavy exchange (the u-channel heavy
> pole adds `y^{Δ_h/2 − …}`), and for `Ĉ₅^{ds}` (no `y^{−1}`, p ≤ 0). The Einstein term (graviton = twist 1 in (13)) would give p = 3
> (margin ½), but its light function under `Ψ₄, Ψ₅` is zero anyway (Lemma A).

*Proof.* (2.1) with finitely many residues (§2.2), endpoint convergence (§2.4), and the pole census (§2.1): moving every `δ₁₃` contour
to `Re δ₁₃ = −½` crosses only the poles at `0` and `σ_α − j ≥ 0`, each giving `y^{m_α − δ*}` times at most `log y`; Lemma Y gives
`m_α^{pol} − δ* ≥ 0`; the explicit `f₄` gives the extra `y^{−1}`. The remaining line integrals are `O(y^{m_α + 1/2})` by absolute
convergence (`Γ(δ₁₃)Γ(δ₂₄)` decay; checked numerically in §2.6). The constants are continuous in the mean angle on the circle. ∎

**Grade.** **[D]**, given (R1), (R2) [Q]. Lemma Y [P]. The k = 4 kernel bound and the criterion are LEMMAS §4.1 [D].
The statement is free-standing: it does not use H_EFT(a), LRC\* or any light-cone hypothesis. **The [I] of L6′ is removed.**
v2's "subsumed by H_EFT(a)" is no longer needed.

### 2.6 Numerical checks [DV]
All on LIGHTSIDE_DIST's product family `T′ = σ ⊗ T` (spin 2, Δ = 5), P1's (2,4) polarisations, mean angle 0.1, symmetrised over the
1↔3 pairs of Wick classes (1↔3 maps `s ↔ t` and `1243 ↔ 1324`). Twist content of the (13) channel: one free-field (13)-link
(twist 1) in 1243/1324; none (double-trace-like) in s, t.

**(a) Δφ scan** (`rh_t1_collinear.py scan`, `cache/t1_report.txt`). Pieces: the top-Regge monomials `D^J/4`, `D = δ₁₂ − δ₁₄`, of the
first non-vanishing degree in each pair (`J = 4` for 1243+1324 at k = 2, the LSD borderline; `J = 8` for s+t).

| Δφ | y | 1243+1324, k = 2, `D⁴/4` (non-dispersive) | slope | s+t, k = 2, `D⁸/4` (dispersive) | slope | class s, k = 4, `D⁸/4` (dispersive) | slope |
|---|---|---|---|---|---|---|---|
| 1.2 | 0.638 | −8.770e-3 | — | +1.724e-1 | — | +6.171e-6 | — |
| 0.6 | 0.175 | −2.319e-1 | 4.73 | +9.935e-1 | 2.53 | +4.552e-5 | 2.88 |
| 0.3 | 0.0447 | −4.769 | 4.36 | +5.390 | 2.44 | −1.767e-4 | 1.96 |
| 0.15 | 1.12e-2 | −7.363e1 | 3.95 | +3.067e1 | 2.51 | −3.208e-3 | 4.18 |
| 0.075 | 2.81e-3 | −9.325e2 | 3.66 | +1.847e2 | 2.59 | −2.808e-2 | 3.13 |
| 0.0375 | 7.03e-4 | −1.042e4 | 3.48 | +1.069e3 | 2.53 | −1.934e-1 | 2.78 |

("slope" = local `−d log|V|/d log Δφ`.) The Δφ = 1.2 entry of the first column is LIGHTSIDE_DIST's P1 value (−8.7642e-3 + tail).
**Reading.** Along this physical scan the mean angle is fixed but `φ₁, φ₃` move, so the Mellin coefficient functions drift (by
`O(Δφ²) = O(y)`, with large coefficients), and the asymptotic regime of the Mellin transform sets in late (at Δφ = 0.3 the double-pole
term and the rest still cancel by factors 1.2–27, table (b)). The dispersive columns are consistent with `y^{−1}(A + B log y)`
(a sign change of `A + B log 8y` sits near y ≈ 0.05–0.2, which inflates the local slopes); the non-dispersive column grows faster
than its twist-rule value 3 (see (d)). The clean statement is the frozen-angle one:

**(a′) Frozen angles** (`rh_t1_frozen.py`, `cache/t1_frozen_k*_0.0375.txt`). At fixed `(φ₁, φ₃)` (Δφ₀ = 0.0375) the og gauge gives
`V(y) = (y₀/y)M(8y)` exactly, `M = (1/2πi)∫(8y)^{−δ}Γ(δ)²Φ(δ)`. Computed as `Res₀(y) + I_{−0.4}(y)`; check against the line at 2/3 at `y₀`:
`1.06862311e3` both (k = 2), `−1.9334683e−1` vs `−1.9334681e−1` (k = 4).

| y | 7.0e-4 | 1e-3 | 1e-4 | 1e-5 |
|---|---|---|---|---|
| `I_{−0.4}/Res₀`, k = 2 | 1.5e-2 | 2.3e-2 | 1.9e-3 | 1.8e-4 |
| `I_{−0.4}/Res₀`, k = 4 | 4.4e-2 | 6.7e-2 | 4.9e-3 | 4.4e-4 |

and `Res₀ = A + B log 8y` on `y = 10⁻³…10⁻⁹` to `3·10⁻¹³` (k = 2: `A = −467.7`, `B = −293.4`; k = 4: `A = 0.1996`, `B = 0.0743`).
So `V = O(y^{−1}|log y|)`, i.e. p = 2.

**(b) Contour shift** (`rh_t1_collinear.py shift`, Δφ = 0.3): the line integrals at `Re δ₁₃ = 2/3` and `0.4` agree, and equal the line
at `−0.4` plus the double-pole residue at `δ₁₃ = 0` (computed on a circle of radius 0.2): for the dispersive classes there is no other
singularity in `−0.4 ≤ Re δ₁₃ ≤ 2/3`; for the non-dispersive class there is one in `(−0.4, 0.4)`.

| pair / class, k | `|I(0.4)/I(2/3) − 1|` | `|(I(−0.4) + Res₀)/I(2/3) − 1|` | `Res₀/I(2/3)` |
|---|---|---|---|
| s+t, k = 2 (dispersive) | 1.6e-9 | 1.6e-9 | −1.18 |
| s, k = 4 (dispersive) | 2.3e-7 | 2.3e-7 | −26.7 |
| t, k = 4 (dispersive) | 1.2e-8 | 1.2e-8 | −26.7 |
| 1243, k = 2 (non-dispersive) | 1.1e-7 | **6.9e-4** (the `n`-sum on `Re δ₁₃ = −0.4` does not converge) | +0.47 |

(The individual t (k = 2) and 1324 entries are `10⁻¹²–10⁻¹⁵`, i.e. zero, and are not listed.)
**Laurent test** (`rh_t1_laurent.py`): `(1/2πi)∮_{|δ₁₃|=0.1, 0.2} Φ = 10⁻¹³–10⁻¹²` against `|Φ| ~ 10⁻³–10` for s (k = 2, 4) and 1243
(k = 2): Φ is analytic at 0, so `Γ(δ₁₃)²` gives exactly one log. Along the real axis `Φ_s(δ₁₃)` (k = 2) is smooth and monotone from
`δ₁₃ = −0.7` to `0.6` (8.17 … 17.25).

**(d) The non-dispersive tail** (`rh_t1_tail.py`, `cache/t1_tail_0.3.txt`): class 1243, k = 2, `D⁴/4`, residues `r_n(δ₁₃)` at real
`δ₁₃`, fit `|r_n| ∝ n^{−p}` on `n ∈ [40, 160]`: `p = 2.60, 2.18, 1.64, 1.25, 0.83, 0.47, 0.18` at `δ₁₃ = 1.2, 0.9, 2/3, 0.5, 0.35, 0.2,
0.05` (the last 40 residues are grid-noise limited). The sum over n stops converging near `Re δ₁₃ ≈ 0.4`: a singularity of Φ inside the
strip, which is what the extra growth of the 1243 column and the failure of its contour shift show.

**(c) Lemma Y** — §2.3. **(e) Dispersiveness** — §2.2.

### 2.7 What this does not cover
* The actual R⁴ ⟨TTTT⟩ is not evaluated; (R1) is used in the form "finite sum of monomial × D-function", not term by term.
* The exponent p = 2 is an upper bound; cancellations can only lower it.
* The bound **requires dispersiveness**. For non-holographic terms with infinitely many non-zero s-channel residues (the toy's 1243
  class) the collinear growth exceeds the twist rule and can reach the edge of the criterion (local slope 3.5 at Δφ = 0.04, §2.6). For
  ⟨TTTT⟩ dispersiveness of every tree term follows from integer Mellin shifts (§2.2) and is CLS-D's statement [Q].
* The physical Δφ-scan of the toy is not yet asymptotic at Δφ = 0.04 (local slopes 2.5–2.8 in the dispersive columns); the p = 2
  statement rests on the analytic structure, which is checked directly ((a′), (b)).

---

## 3. Target 2: H_EFT(a) for spinning externals

### 3.1 Statement
H_EFT(a) (v2 row L7(a); F_NORM route II step 5): for a heavy exchange `W_h` (mass m, spin `J_h`) in ⟨TTTT⟩, the term-wise value of
its double-trace part, `DT[W_h] := Σ_{double-trace residues} Ω_k`, equals the term-wise (P) value of its EFT expansion,
`TP_R := Σ_{d ≤ R} Ω_k[EFT_d]`, as an asymptotic series in `1/m²`; in particular `DT = Ω_k[EFT_{d*}](1 + o(1))` as `m → ∞`, where `d*` is the
lowest degree with non-zero term-wise value (`d* = 4`, the R⁴ term, for Ψ₄ and a scalar χ).

### 3.2 Derivation [D]
1. **Finitely many residues.** `W_h`'s s-channel Mellin poles: the heavy poles `δ₁₂ = δ_m ≈ 3 − Δ_h/2 − m` (non-integer for generic
   `Δ_h`: the heavy block, with CLS-D's `sin²` factor) and the double-trace poles of the measure at integers. By §2.2 only finitely
   many double-trace residues, `δ₁₂ ∈ {0, −1, …, −n₀}`, are non-zero — for `W_h` and for every EFT polynomial alike. (R3), (R4).
2. **Uniform convergence there.** At `δ₁₂ = −n`, `n ≤ n₀`, the exchange `Σ_m ρ_m(y)/(x − δ_m)` (x = δ₁₂, δ₁₄, δ₁₃ per channel) and
   its `δ₁₂`-derivative are given by the convergent EFT series `−Σ_j x^j ρ_m(y)/δ_m^{j+1}` whenever `|x|, |y| < |δ₀| ≈ Δ_h/2`, with
   geometric remainder `X^{(R+1)} = Σ_m ρ_m x^{R+1}/(δ_m^{R+1}(x − δ_m))`; organised by total Mellin degree (powers of `1/m²`).
3. **The δ₁₃ integral.** The weight `Γ(δ₁₃)Γ(δ₂₄)Γ(δ₁₄)Γ(δ₂₃)` decays like `e^{−2π|Im δ₁₃|}`; the region `|δ₁₃| ≳ Δ_h/2`, where the
   series diverges, contributes `O(e^{−c Δ_h})`. So `|TP_R − DT| ≤ C_R (const/m²)^{R+1−d*}|DT| + O(e^{−cΔ_h})`: an asymptotic series
   with optimum `R* ∝ Δ_h` (the behaviour FH §4.3 found for scalars, "R* ≈ π τ_h").
4. **Pointwise ⇒ paired.** Steps 1–3 hold at every fixed `(a, b)`; by §2.5 every term, and the remainder `X^{(R+1)}`, obeys the same
   `O(Δφ^{−2} log)` bound with constants proportional to the same sup-norms, so the bound survives the `C_ν` and `ν` pairings: the
   ν-paired form follows.
5. **Leading order.** The degree < 4 EFT terms have zero term-wise value (Lemma A for contacts, n < 4); the degree-4 term's value is
   `𝒜₄` of its top Regge coefficient ((P), LSD [D + DV]), which is `g₄^χ(1 + O(m⁻²))` times the R⁴ structure (F_NORM A3). Hence
   `Ĉ₄^{ds}[DT(W_χ)] = 𝒜₄[R⁴(g₄^χ, ĝ₄^χ)] + O(g₄^χ/m²)`, which is route II step 5.

**Why F_NORM's fixed-angle test (B3) failed.** In the product `σ ⊗ T_free` the connected classes 1243/1324 contain a single free-field
(12)-link `(2t₁−1)^{−1/2}`: their s-channel content (σ double traces ⊗ twist-1 free bilinears) sits at half-odd twist relative to the
double traces of `T′`, the double commutator does not vanish there, and **all** residues `δ₁₂ = −n` contribute, up to `n ~ |δ₀|`,
where the EFT series diverges. That is a property of this non-holographic family, not of the pointwise form: a holographic ⟨TTTT⟩ is
dispersive (CLS-D), and so are the s, t classes of the same family.

### 3.3 Checks [DV] (`rh_t2_heft.py`, `rh_t2_report.py`, `cache/t2d_pairs_P1_k*.txt`)
* **Exchanges.** Heavy σ-exchange in s, t, u with the exact d = 3 Mack polynomials and weights (exactmom, `Δ_σ = 2`), `J_h = 0, 2, 4`,
  `τ_h = 41, 61, 81` (`|δ₀| = 18.5, 28.5, 38.5`). The EFT remainders `X^{(R+1)}` are computed directly (no cancellation), organised by
  total Mellin degree. `DT` is read on the plateau `V[X^{(j)}]`, `j ≤ d*` (orders below `d*` are term-wise zero; plateau flat to
  `10⁻⁵–10⁻⁸`).
* **Dispersive pair s+t (T-factor `⟨T₁T₂⟩⟨T₃T₄⟩ + ⟨T₁T₄⟩⟨T₂T₃⟩`; spin-2 externals, genuine double-trace structure), k = 2, one angle
  pair (P1: Δφ = 1.2, random null (2,4) polarisations). Other angle pairs were not run; the derivation is angle-independent.**

| `τ_h` | `J_h` | `d*` | `DT/EFT_{d*} − 1` | `|TP_R − DT|/|DT|` at `R = 8, 12, 16, 24, 39` | Lemma A `(S_h + DT)/|DT|` |
|---|---|---|---|---|---|
| 41 | 0 | 8 | +0.123 | 1.1e-1, 5.9e-4, 6.4e-7, 1.2e-12, 4.1e-21 | 5.6e-5 |
| 41 | 2 | 8 | −0.031 | 3.2e-2, 3.2e-4, 1.2e-7, 2.2e-13, 5.2e-21 | 8.5e-5 |
| 41 | 4 | 8 | −0.241 | 3.2e-1, 2.4e-3, 4.6e-7, 4.0e-13, 1.4e-21 | 3.8e-5 |
| 61 | 0 | 8 | +0.0257 | 2.5e-2, 1.3e-5, 8.0e-10, 7.6e-18, 4.2e-30 | (m-sum not converged) |
| 61 | 2 | 8 | −0.041 | 4.3e-2, 9.4e-6, 1.0e-10, 2.1e-18, 6.3e-30 | (″) |
| 61 | 4 | 8 | −0.157 | 1.9e-1, 6.5e-5, 7.3e-10, 3.8e-18, 4.3e-30 | (″) |
| 81 | 0 | 8 | +0.0086 | 8.6e-3, 7.9e-7, 5.6e-12, 9.5e-22, 5.3e-37 | (″) |
| 81 | 2 | 8 | −0.029 | 3.0e-2, 6.4e-7, 5.9e-13, 3.3e-22, 9.5e-37 | (″) |
| 81 | 4 | 8 | −0.103 | 1.1e-1, 4.5e-6, 6.0e-12, 6.3e-22, 9.3e-37 | (″) |

  The single-trace sum `S_h` is truncated at `m ≤ 160`; it is converged only at `τ_h = 41`, where Lemma A per pair holds to
  `≤ 9·10⁻⁵`. In this sector the first non-zero EFT degree is 8 (the Regge counting of the s+t classes is lower than 1243's by 4), so
  `d* = 8` plays the role R⁴ plays for ⟨TTTT⟩ under Ψ₄.

* **Same pair, k = 4 (`Ψ₄`, the functional of `Ĉ₄^{ds}`):** the first non-zero EFT degree moves to `d* = 10` (two units above k = 2,
  as the Regge counting requires).

| `τ_h` | `J_h` | `DT/EFT_{d*} − 1` | `|TP_R − DT|/|DT|` at `R = 10, 14, 18, 26, 39` | Lemma A |
|---|---|---|---|---|
| 41 | 0 | +0.306 | 2.3e-1, 3.3e-3, 1.2e-5, 1.5e-10, 6.3e-18 | 3.8e-4 |
| 41 | 2 | +0.043 | 4.1e-2, 3.0e-4, 1.6e-6, 1.0e-12, 6.8e-19 | 5.2e-4 |
| 41 | 4 | −0.258 | 3.5e-1, 7.7e-3, 3.2e-6, 7.4e-12, 6.0e-18 | 8.2e-4 |
| 61 | 0 | +0.075 | 7.0e-2, 8.4e-5, 1.9e-8, 1.4e-15, 3.4e-26 | — |
| 61 | 2 | −0.031 | 3.2e-2, 1.9e-5, 8.8e-10, 1.2e-16, 4.3e-27 | — |
| 61 | 4 | −0.198 | 2.5e-1, 2.5e-4, 6.5e-9, 5.7e-17, 3.5e-26 | — |
| 81 | 0 | +0.031 | 3.0e-2, 5.4e-6, 1.4e-10, 2.1e-19, 1.5e-32 | — |
| 81 | 2 | −0.029 | 3.0e-2, 1.6e-6, 7.3e-13, 2.8e-20, 2.1e-33 | — |
| 81 | 4 | −0.138 | 1.6e-1, 1.8e-5, 5.8e-11, 5.9e-21, 1.6e-32 | — |

  For 1243+1324 at k = 4 the series again diverges (`7.6·10²` at R = 14, `7.6·10²²` at R = 39, `τ_h = 41, J_h = 0`).

* **Non-dispersive pair 1243+1324 (B3's dominant classes), k = 2:** `DT` plateau at `j ≤ 4`, then `|TP_R − DT|/|DT|` grows without
  bound: `1.6·10²` (R = 12), `1.4·10⁵` (16), `1.8·10²³` (39) at `τ_h = 41, J_h = 0`; the onset moves to larger R as `τ_h` grows
  (`R ≈ 12, 16, 24` for `τ_h = 41, 61, 81`), as expected when residues up to `n ~ |δ₀|` enter. This reproduces F_NORM B3.
* **Leading order.** For the scalar exchange (the only case route II uses) `DT/EFT_{d*} − 1 = 0.123 → 0.026 → 0.0086` (k = 2) and
  `0.31 → 0.075 → 0.031` (k = 4), local powers `≈ τ_h^{−3.6}` and `τ_h^{−3…−3.5}`. For `J_h = 4` it falls like `≈ 1/τ_h` (−0.24 → −0.10); for `J_h = 2` it stays at `−0.03 ± 0.01` up to `τ_h = 81`
  [NE, not understood; the asymptotic convergence of `TP_R` is unaffected].

### 3.4 Grade
**H_EFT(a), spinning: [D-sketch] → [D]** (pointwise in angle and ν-paired, for dispersive terms), given (R1)–(R4) [Q]; **[DV]** on the
dispersive spin-2 sector of the product family. Not exercised: the graviton's tensor structures, a spinning χ coupled through `W²`
(the check uses σ-exchanges dressed by spin-2 external factors), and `J_h ≥ 6`.

---

## 4. H_EFT(b): how the literature states it, and how our use compares

**CHMRS (2106.10274).** §4 ("Light contributions to holographic sum rules"): *"We assume that the light OPE data agrees with a tree-level
AdS effective field theory (EFT) with a derivative expansion"*, and, for a dispersive functional that vanishes on all double-trace
families with `τ > 2Δ_φ + 2p`, *"We assume that the OPE data of double-traces with τ ≤ 2Δ_φ+2p is well-approximated by a bulk EFT with
a truncated derivative expansion."* With this, ω|_light is computed from the Regge limit of `𝒢^EFT` (their eq. `eftformula`, proved in
their App. "Computing ω|_light"). Separately (§5), a "technical assumption" `O(Σ g/R_AdS²) = O(1/R_AdS²)` — that infinite towers of
higher-derivative couplings do not conspire to produce contributions enhanced at large `Δ_gap` — is stated ("It seems utterly
implausible to us that a counter-example … could be found") and argued, not proved, in §6.4 ("Comments on truncation errors"): lowering
the cutoff `M² → (1−ε)M²`, the flat EFT series converges for `|s|,|t|,|u| < M²`, and AdS corrections are expected to grow only like
powers `n^j`; "A rigorous proof … is left to the future."

**CLS-D (2311.04271), §2.2:** *"One can argue that ω|_light is determined by the Wilson coefficients of the low-energy EFT in AdS"*; the
spinning light side is not computed ("we leave exploration of light states to future work").

**Our use.** H_EFT(b) as used in (2′) is the conjunction of exactly these two CHMRS assumptions, applied to ⟨TTTT⟩ and the k = 4 (and
k = 5) functionals: (i) the finitely many double-trace families on which `Ψ₄` is non-zero (dispersiveness, §3.2 step 1) carry the OPE
data of the tree-level bulk EFT with a truncated derivative expansion; (ii) the heavy sum and the contact tower do not conspire
(CHMRS's technical assumption). **Nothing beyond CHMRS is assumed**, with one difference of scope: CHMRS state it for a scalar
⟨φφφφ⟩ (with gravity), we for ⟨TTTT⟩, where the EFT data include graviton exchange and the three 8-derivative couplings. What this note
adds is that the *mechanism* by which (i) is used — evaluating the light side of heavy-exchange double traces term-wise — is derived
(H_EFT(a), §3), so H_EFT(b) is the only remaining assumption between `g̃₄` and `g₄`.

---

## 5. Net effect on Result Q's grade chain

| link | v2 grade | now |
|---|---|---|
| L6′ (a = b, R⁴ arc) | **[I]** | **[D]** (given [Q] representation facts; Lemma Y [P]; mechanism [DV]) |
| L7(a) H_EFT(a) | [DV] scalar toy, [D-sketch] spinning | **[D]** spinning, pointwise and ν-paired, given (R1)–(R4) [Q]; [DV] on the dispersive spin-2 sector |
| L6 (F) | [D given Lemma A, H_EFT(a) ν-paired, heavy dictionary] | same, with H_EFT(a) now [D]: **(F) [D]** |
| L7(b) H_EFT(b) | [A] | [A] (unchanged; CHMRS-standard) |
| "(Pos) ⇒ (2), upper end" | **[I] in one link**; otherwise [A] | **[D-sketch]**: weakest links H_K through ⟦EESTEP⟧ (and LRC-a EE at k = 4), `g₅`-(P) [D-sketch] (δ = 0.2 survives either reading), κ > 0 [D-sketch + M] if `(t₄, α_T) ≠ 0`; **plus [A] H_EFT(b) for (2′)** |
| "(Pos) ⇒ (2), lower end" | same, no constant | same, no constant |

---

## 6. Proposed wording

*Letter (one sentence, replaces "one technical link open [I] (I4)"):*
> "In the gauge of the functional the collinear variable enters only through `X₁₃` and the polarisations of the light-ray pair, so the
> light function of any contact term — the R⁴ term included — is a finite sum of Mellin–Barnes integrals bounded by `Δφ^{−2}log(1/Δφ)`,
> which pairs absolutely with the celestial kernel (growth below `Δφ^{−7/2}` suffices)."

*Supplement (table rows):*
> **L6′** | the R⁴ (and every contact / heavy-exchange double-trace) light function obeys `|L| ≤ CΔφ^{−2}(1+|log Δφ|)`; integrable
> against `C_ν` with margin 3/2 | Q | **[D]** (given the polynomial×D-function form of spinning contact diagrams [Q]); counting lemma
> [P]; mechanism [DV]
>
> **L7** | (a) term-wise light value of a heavy exchange's double traces = actual value, as an asymptotic series in `1/m²`, pointwise in
> angle for dispersive terms (hence ν-paired); (b) light double-trace data = tree-level EFT, no conspiracy | Q | (a) **[D]**, [DV] scalar
> Mellin and dispersive spin-2 sector; (b) [A]

---

## 7. Suggested corrections to existing documents (not edited here)
1. **UNIFIED_THEOREM_v2.md §3 net-grade table and §4 rows L6′, L7, L6:** L6′ [I] → [D]; L7(a) → [D] (spinning, pointwise for dispersive
   terms); L6 (F) → [D] without qualification; "(Pos) ⇒ (2), upper end" → [D-sketch] (+ [A] for (2′)). §11 item 7 can be closed.
2. **F_NORM.md §0 table B3, §3.4, §4 row H_EFT(a); RESULT_Q_STATEMENT box at line ~132; UT's F_NORM box:** "H_EFT(a) is false pointwise
   in angle [NE]" should read "fails pointwise in the non-dispersive (1243/1324) classes of the toy, where the (12) channel carries
   twist-1 free-field content; holds pointwise in its dispersive classes (s+t) to `10⁻¹⁸–10⁻³⁷`; holographic ⟨TTTT⟩ terms are dispersive".
3. **LEMMAS_CLOSE.md §4.1 item 4** and **LRC_CLOSE.md §6.5 item 3:** the borderline arc's a = b behaviour is not a light-cone-edge
   question of the kernel: in the og gauge it is controlled by the (13)-channel twist of the term (`p = 4 − τ_min`), and for contacts
   `p = 2`. The heuristic "η^{5/2}, p = 5/2" is superseded.
4. **prl/supplement.tex** table rows L6′, L7 and the "Collinear set" paragraph (lines ~173–176, 235–239); **prl/letter.tex** I4/I5
   items (~lines 239, 266).
5. **LIGHTSIDE_DIST.md §6 item 4** ("integrability near a = b [I]"): closed by §2.5 here.

---

## 8. Reproduce

    cd ads4/code/r4heft && ./run_all.sh      # about 2 CPU-h, 2 processes at a time

| file | content |
|---|---|
| `rh_common.py` | Setup wrapper (lsd_core, production grid), residues for arbitrary Mellin amplitudes, symmetrised values, tails |
| `rh_explore.py` | which `δ₁₂ = −n` residues survive per Wick class (dispersiveness), k = 2, 4 → `cache/explore_k*.txt` |
| `rh_degree.py` | symmetrised values of `D^J/4`, J = 2…8, per class (first non-zero degree) → `cache/degree_k*.txt` |
| `rh_lemma_count.py` | Lemma Y: exhaustive enumeration and numerical y-orders of the ε's → `cache/lemmaY.txt` |
| `rh_t1_collinear.py` | Δφ scan and δ₁₃ contour-shift test → `cache/t1_report.txt` |
| `rh_t1_frozen.py` | frozen-angle Mellin asymptotics `Res₀ + I_{−0.4}` → `cache/t1_frozen_k*_0.0375.txt` |
| `rh_t1_laurent.py` | `∮Φ` around `δ₁₃ = 0` (no extra pole) → `cache/t1_laurent_*.txt` |
| `rh_t1_tail.py` | residue-tail exponent `p(δ₁₃)` of the non-dispersive class → `cache/t1_tail_0.3.txt` |
| `rh_t2_heft.py`, `rh_t2_report.py` | heavy exchanges with exact Mack residues, EFT remainders by total degree, single-trace sums; pair reports → `cache/t2d_pairs_P1_k*.txt` (`RH_MODE=j`: by propagator order, `cache/t2_pairs_P1_k2.txt`) |
