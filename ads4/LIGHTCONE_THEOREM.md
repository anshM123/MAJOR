# First-order large-spin correction as a theorem (conditional on a twist gap)

Status: written 15 Nov 2026; checked by one independent automated referee with a hostile brief (same model
family as the author; not a human expert). Verdict: "I could not break it. The proof of T is correct as
mathematics; established after the listed fixes" — all applied below (wording of hypotheses, two added
lines in T′, an untilted corollary, description). The referee verified block positivity and integer spacing
by exact-rational Casimir recursion in `d = 3, 4`, the kernel identities to `10^{−27}`, the Step-4 inequality
on 20,000 random positive measures, and that `−2b_*` is the Fitzpatrick–Kaplan–Poland–Simmons-Duffin
coefficient after converting the `z`-average to spin (ratio `1.0003` at `z = 10^{−6}`).
The strategy of §3 (a multi-point functional with double zeros at the higher towers, applied to the exact
monomial expansion) was proposed by an audit agent after reading Pal–Qiao–Rychkov (arXiv:2212.04893) and
van Rees (arXiv:2412.06907).

**How to describe it without overstatement (referee's sentence).** "If the lowest-twist exchanged operators
are finitely many and separated by a twist gap from everything else, unitarity and crossing prove that the
block-weighted total of `(τ − 2Δ)` over all operators near twist `2Δ`, smeared over spins up to `z^{−1/2}`, equals
`−2b_* z^{τ_*/2}(1+o(1))`: the leading large-spin coefficient in averaged form, as a theorem with no usable rate
and no statement about individual operators or about gravity." It is **not** "rigorous long-range attraction
in an emergent AdS4": no AdS or large-`c_T` input is used, the theorem applies equally to the 3d Ising model,
and `1/c_T` enters only through the Ward identity. **Physical content:** none beyond the lightcone bootstrap
(arXiv:1212.3616, 1212.4103) and the Lorentzian inversion formula (arXiv:1703.00278), which give more (each
spin `J > 1`, every trajectory) at the level of derivation. Novelty: no such theorem was found in the two
rigorous papers or their citing papers through Nov 2026 ("not found", not "does not exist"). False without
hypothesis (H): §5.

## 1. Setting and hypotheses

* A unitary CFT in `d ≥ 2`, a Hermitian scalar primary `φ` of dimension `Δ`, and
  `G(z, z̄)`, `0 < z, z̄ < 1`, the four-point function of `φ` in the usual cross-ratio normalisation.
* **(P1) Positive monomial expansion** [imported theorem: Pal–Qiao–Rychkov eqs. (2.2)–(2.3) and their
  references for the Euclidean axioms and reflection positivity; van Rees eq. (1) and eq. (33)]:
  `G(z, z̄) = Σ_{h, h̄} a_{h,h̄} z^h z̄^{h̄}` with `a ≥ 0`, convergent for `0 ≤ z, z̄ < 1`, `G(z,z̄) = G(z̄,z)`. A primary of
  twist `τ` and spin `ℓ` with coefficient `p` contributes only powers `z^{τ/2 + n}`, `n = 0, 1, 2, …`, and its
  `n = 0` term is exactly `p z^{τ/2} k_{τ+2ℓ}(z̄)`, `k_{2h}(x) := x^h ₂F₁(h, h; 2h; x)`. (Integer spacing and
  positivity re-checked by the referee by exact recursion in `d = 3, 4`. This does not rely on Hogervorst's
  dimensional-reduction formula.)
* **(P2) Crossing:** `((1−z)(1−z̄))^Δ G(z, z̄) = (z z̄)^Δ G(1−z, 1−z̄)`.
* **(H) Isolated minimal twist:** the identity is the only term with `h = 0` (`a_{0,0} = 1`; automatic in
  `d ≥ 3` by the unitarity bounds, **false in every local 2d CFT**, where the stress tensor and the vacuum
  descendants have twist zero — so for `d = 2` the theorem is empty for theories with a stress tensor); and
  there are `τ_* > 0`, `δ > 0` such that the non-identity primaries in `φ × φ` of twist `< τ_* + δ` are finitely
  many, at least one, all of twist exactly `τ_*`, and `τ_* + δ < 2Δ`. Put
  `b_* := Σ_{τ_O = τ_*} p_O Γ(τ_* + 2ℓ_O)/Γ(τ_*/2 + ℓ_O)² > 0`.
  (H) cannot at present be verified in any interacting CFT (for the 3d Ising model it is supported by
  numerics only), and no closed-form non-trivial example satisfying it is known.

No assumption on `d`, on the central charge, on Regge behaviour, on degeneracies or on mixing.

**Measures.** For `0 < z < 1` let `A_h(w) := Σ_{h̄} a_{h,h̄} w^{h̄}` and

`M_z := z^Δ Σ_h A_h(1−z) δ_{x = 2h − 2Δ}` (a positive measure on `R`; "monomial-twist measure"),
`μ_z := z^Δ Σ_{primaries} p_{τ,ℓ} k_{τ+2ℓ}(1−z) δ_{x = τ − 2Δ}` (the leading components of the primaries).

## 2. Statements

**Theorem T.** Assume (P1), (P2), (H). Fix `0 < η < 2`, `β ∈ (0, min(δ, 2)/(2Δ))`, an integer `M > τ_*/(2β)`,
and put `y := z^β`, `γ := min( min(1, δ/2) − βΔ, Mβ − τ_*/2 ) > 0`. Then as `z → 0`

`∫_{|x|<η} x · y^{x/2} dM_z(x) = −2 b_* z^{τ_*/2} + O(z^{τ_*/2 + γ})`,
`∫_{|x|<η} x² · y^{x/2} dM_z(x) = O(z^{τ_*/2 + γ})`,

for `z` small enough that `z^β e^{2M} ≤ 1/2`; the implied constants depend on the correlator, `η`, `β`, `M`.

**Untilted form (referee).** If moreover `γ > βη/2` — i.e. `β(Δ + η/2) < min(1, δ/2)` and `M > τ_*/(2β) + η/2` — then

`∫_{|x|<η} x dM_z(x) = −2 b_* z^{τ_*/2} (1 + o(1))`.

(Because `|x(1 − e^{−Sx})| ≤ S e^{Sη} x² e^{−Sx}` on `|x| < η` and `e^{Sη} = z^{−βη/2}`, the difference from the tilted
moment is `O(log(1/z) z^{−βη/2} z^{τ_*/2+γ})`. Without this step the tilted statement would only bound the
untilted mean from one side.)

**Corollary T′ (primaries).** If moreover `2Δ + η < τ_* + 2`, then `M_z` and `μ_z` coincide on `{|x| < η}`, so the
same two formulas hold with `μ_z`: for the primaries of twist `τ ∈ (2Δ − η, 2Δ + η)`,

`z^Δ Σ p_{τ,ℓ} k_{τ+2ℓ}(1−z) (τ − 2Δ) y^{(τ−2Δ)/2} = −2 b_* z^{τ_*/2} (1 + o(1))`.

In particular there are infinitely many primaries with `2Δ − η < τ < 2Δ`, of unbounded spin, and the
weighted second moment of `τ − 2Δ` is of smaller order than the first. (*Proof of the last claim.* The
weighted first moment is negative for all small `z`. A single primary in the window has weight at most
`C z^Δ log(1/z) · y^{−(2Δ−τ_*)/2}`, which is `o(z^{τ_*/2})` because `β < 1`; the same holds for any family of bounded
spin, since `Σ p` over it is finite. So no finite set and no bounded-spin set of primaries can carry the
moment.)

*Stress tensor in `d = 3`.* If the stress tensor is the only non-identity primary of twist `< 1 + δ` in
`φ × φ`, then `τ_* = 1`, `b_* = p_T · Γ(5)/Γ(5/2)² = 128 p_T/(3π)` with `p_T ∝ Δ²/c_T`, and T′ applies for
`(1+δ)/2 < Δ < (3 − η)/2` — inside `½ < Δ < 3/2`. For `Δ ≥ 3/2` only the monomial measure `M_z` (primaries together
with their SL(2) descendants) is controlled. Referee's normalisation: `p_T = (3/8) Δ² c_free/c_T`.

## 3. Proof

**Step 1 (an exact identity with a positive measure).** By (P1) with the roles of the variables exchanged,
`G(1−z, y) = Σ_h y^h A_h(1−z)`, so by (P2) with `z̄ = 1 − y`,

`∫ y^{x/2} dM_z(x) = z^Δ y^{−Δ} G(1−z, y) = ((1−z)/(1−y))^Δ G(z, 1−y)`.   (3.1)

**Step 2 (the direct channel to first order).** By (P1) and (H), `G(z, 1−y) = 1 + z^{τ_*/2} f(y) + R(z,y)` with
`f(y) := Σ_{τ_O = τ_*} p_O k_{τ_*+2ℓ_O}(1−y)` and `R = Σ_{h ≥ h_1} z^h A_h(1−y) ≥ 0`, `h_1 := τ_*/2 + min(1, δ/2)` (the
other powers of `z` are `τ_*/2 + n`, `n ≥ 1`, from the minimal-twist primaries and `≥ (τ_*+δ)/2` from the
rest). Fix `z_0 ∈ (0,1)`. For `z ≤ z_0`, by positivity, `R ≤ (z/z_0)^{h_1} G(z_0, 1−y)`, and by (P2)
`G(z_0, 1−y) = (z_0(1−y)/((1−z_0)y))^Δ G(1−z_0, y) ≤ C y^{−Δ}` for `y ≤ 1/2` (`G(1−z_0, ·)` is increasing). Hence

`0 ≤ R(z,y) ≤ C z^{h_1} y^{−Δ}`.   (3.2)

The function `k_{2h}(1−y) = (1−y)^h (Γ(2h)/Γ(h)²) Σ_{n≥0} ((h)_n/n!)² y^n [2ψ(n+1) − 2ψ(h+n) − log y]` is a
convergent series in `y^n` and `y^n log y`; so is `(1−y)^{−Δ}(1 + ε f(y))`, and with `y = e^{−2s}`, `ε := z^{τ_*/2}`:

`((1−z)/(1−y))^Δ G(z,1−y) = Σ_{m=0}^{M−1} e^{−2ms} [ W_m + ε(C_m + 2 B_m s) ] + r(z, s)`,   (3.3)

with constants `W_m, C_m, B_m` depending on `z` only through the factor `(1−z)^Δ`, `B_0 = (1−z)^Δ b_*`, and

`|r(z,s)| ≤ C ( y^M + z^{h_1} y^{−Δ} )` for `y ≤ 1/2`   (3.4)

(using `ε log(1/y) ≤ 1` for small `z`).

**Step 3 (two functionals).** Fix the step `1`, `S := (β/2) log(1/z)`, sample points `s_j := S + j`, `j = −M, …, M`,
and let `Q(x) := Π_{m=1}^{M−1} 4 sinh²((x − 2m)/2)`, `P := Q(0) > 0`,

`K_2(x) := e^{−Sx} · 4 sinh²(x/2) · Q(x)`, `K_1(x) := e^{−Sx} · 2 sinh(x) · Q(x)`.

Both are finite combinations `Σ_{j=−M}^{M} c_j e^{−s_j x}` with coefficients bounded in terms of `M` only. `K_2 ≥ 0`
has double zeros at `x = 0, 2, …, 2M−2`; `K_1` has a simple zero at `0`, double zeros at `2, …, 2M−2`, and
`K_1'(0) = 2P`. A combination `Λ = Σ c_j (\text{evaluation at } s_j)` with kernel `K(x) = Σ c_j e^{−s_j x}` maps
`e^{−2ms}(α + γ s)` to `α K(2m) − γ K'(2m)`. Therefore, applying `Λ_2`, `Λ_1` to (3.1) = (3.3):

`∫ K_2 dM_z = Λ_2[r]`, `∫ K_1 dM_z = −2P · 2εB_0 + Λ_1[r]`.   (3.5)

At the sample points `y = z^β e^{−2j}`, so by (3.4) `|Λ_i[r]| ≤ C ρ`, `ρ := z^{Mβ} + z^{h_1 − βΔ} = O(ε z^γ)`.

**Step 4 (conclusion).** On `|x| < η < 2`: `4 sinh²(x/2) ≥ x²` and `Q ≥ Q_min = (4 sinh²(1 − η/2))^{M−1} > 0`, so by (3.5)
`∫_{|x|<η} x² e^{−Sx} dM_z ≤ Cρ`, which is the second formula (`e^{−Sx} = y^{x/2}`). Next,
`|K_1(x) − 2P x e^{−Sx}| ≤ C x² e^{−Sx}` on `|x| < η`, and for `|x| ≥ η`, `|K_1/K_2| = |coth(x/2)| ≤ coth(η/2)`, so

`| ∫_{|x|<η} 2P x e^{−Sx} dM_z − ∫ K_1 dM_z | ≤ C ∫_{|x|<η} x² e^{−Sx} dM_z + coth(η/2) ∫ K_2 dM_z ≤ Cρ`.

With (3.5): `∫_{|x|<η} x e^{−Sx} dM_z = −2εB_0 + O(ρ) = −2 b_* z^{τ_*/2} + O(z^{τ_*/2+γ})`. ∎

*Remarks on the proof.* (i) The positivity of `K_2` on all of `R` is what allows the higher towers
(`x ≈ 2m`), the images of the crossed-channel identity (`x = −2Δ`) and everything else to be absorbed: they
enter only through `∫K_2 dM_z`, which is small. No lower limit on `y` beyond `y = z^β` is needed, which is why
the obstruction recorded in `ads4/TARGET.md` §4 (derived with a three-point functional) does not arise.
(ii) The identity image contributes `z^Δ y^{−Δ} = z^{Δ(1−β)}` to (3.1); it is part of `R`'s bound (3.2), and
is `o(ε)` because `β < δ/(2Δ)` and `τ_* + δ < 2Δ`.

**Proof of Corollary T′.** By (P1) a primary of twist `τ` contributes to `M_z` at `x = τ + 2n − 2Δ`, `n ≥ 0`, with
its `n = 0` weight equal to its weight in `μ_z`. A contribution with `n ≥ 1` in `{|x| < η}` needs
`τ < 2Δ − 2 + η < τ_*`, i.e. the identity, which contributes only at `h = 0`. ∎

## 4. What the theorem does and does not say

* It is an **average**: over all primaries with twist within `η` of `2Δ`, weighted by
  `p · k_{τ+2ℓ}(1−z) · y^{(τ−2Δ)/2}`, which selects spins `ℓ ≲ z^{−1/2}`. It does not say that every operator of
  the family has `τ < 2Δ`, does not resolve degenerate or mixed operators, and gives no statement at any
  fixed spin.
* It is **non-uniform in the gap**: the window `β < δ/(2Δ)` and the rate `γ` degenerate as `δ → 0`. In the
  three-dimensional Ising model with external `σ` the hypothesis holds only with `δ < 0.023` (the spin-4
  operator), so the asymptotic regime is reached extremely slowly. In large-`N` vector models `δ ~ 1/N`.
* If the minimal twist is not the stress tensor's (e.g. a scalar of dimension below `d − 2` in `φ × φ`), the
  theorem is about that operator; the stress tensor is then subleading and not isolated.
* It assumes identical external scalars. For non-identical or spinning operators the sign is not fixed by
  unitarity (Li–Meltzer–Poland, arXiv:1511.08025, relate it to the conformal-collider bounds).

## 5. Falsification record

* Free scalar (`Δ = ½` in `d = 3`): zero anomalous dimensions with a stress tensor present. Not a
  counterexample: (H) fails twice (infinitely many twist-1 currents; `τ_* = 2Δ`). Exact computation in
  `ads4/TARGET.md` §0. Same for `φ²` in the free theory and for any theory with infinitely many
  minimal-twist currents.
* Generalised free fields: no operator below `2Δ`; (H) holds vacuously with `b_* = 0` and the statement is
  `0 = 0`.
* Audit's search for a theory satisfying (H) with non-negative average: none found (free `φ²`, GFF,
  large-`N` `O(N)`, 3d Ising `σ` each fail (H) or have the predicted sign). This is a search, not a proof.
* Numerical test of Steps 3–4 on synthetic positive measures (`ads4/code/tauberian_test.py`): towers of 400
  atoms each at `x ≈ 2m` with planted `O(ε)` shifts, optional `o(ε)` mass at negative `x` and between towers.
  With all planted towers cancelled by the kernels (`M = 6`): the functional reproduces the directly computed
  windowed first moment with relative error `8·10^{−2}, 7·10^{−3}, 7·10^{−4}, 7·10^{−5}, 8·10^{−6}` at
  `ε = 10^{−2}, …, 10^{−6}` (towers only), and `72, 4.6, 0.30, 1.7·10^{−2}, 1.0·10^{−3}` at
  `ε = 10^{−2}, 10^{−4}, …, 10^{−10}` with the extra mass. With two towers left uncancelled (`M = 4`) the error
  still tends to zero but only reaches `0.65` at `ε = 10^{−14}`.
  **Reading:** the algebra of the functionals is right. The constants are large, of size `e^{O(M)}` (the first
  uncancelled tower is amplified by `sinh(2M)`), and no explicit constants are claimed.
* Referee's exact and numerical cases: free `φ`, free `φ²` (window contains only `x = 0`, (H) fails); GFF (`0 = 0`);
  2d Ising `σ` and the compact boson (no accumulation at `2Δ`, (H) necessarily fails). 3d Ising with published
  data: `b_T = 1.445`, stress-tensor coefficient `γ_0(T) = −0.00285` against `γ_0(ε) = −0.0926`; the stress tensor
  dominates only for spins above about `4.6·10³`, and with `δ = 0.0227` the theorem's regime needs `z ≲ 10^{−88}`.
  Large-`N` `O(N)`: `δ ≈ 0.15/N`, regime `log(1/z) ≫ 13N`. **No verification on a non-trivial crossing-symmetric
  example satisfying (H) exists; none is available.**
