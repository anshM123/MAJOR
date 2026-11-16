# Harmonic kernels for spinning Regge exchange in CFT3, and why scalar positivity does not carry over

Status: written 16 Nov 2026; mathematics verified numerically (`ads4/code/spinning_kernels.py`,
`ads4/code/spinning_moments.py`); **not yet refereed**. Part A is a theorem of harmonic analysis on the
hyperbolic plane. Part B (which CFT quantities these functions govern) is a symmetry argument, **not** a
derivation from spinning conformal blocks, and the weights with which the kernels enter a dispersive
functional are **not derived**. Nothing here is a locality bound.

## Why this is the first step

In the scalar case (arXiv:2106.10274) the action of a dispersive functional on a heavy block, in the regime
where the AdS curvature matters, is an integral over `ν` of the functional's weight against the harmonic
function `P_{−1/2+iν}(cosh β)` of the transverse hyperbolic plane `H²` (`d = 3`), and positivity reduces to the
signs of its `ν`-moments (`POSITIVITY_LEMMA.md`). For external operators with spin the harmonic function is
replaced by something else, channel by channel. Guessing it from the flat-space Bessel functions is not
safe: several different functions on `H²` have the same flat limit.

## Part A. The functions (theorem)

`G = SL(2,R)` acts on `H²` and on its boundary circle; `K = SO(2)` is the rotation group about a point; `a_β` is
the boost by geodesic distance `β`. Let `π_ν`, `ν ∈ R`, be the spherical principal series on `L²(S¹)`:
`(π_ν(g)f)(θ) = |d(g^{−1}θ)/dθ|^{1/2+iν} f(g^{−1}θ)`. Its `K`-types are `e^{imθ}`, `m ∈ Z`. Define the matrix coefficients

`Φ^ν_{m,n}(β) := ⟨e^{imθ}, π_ν(a_β) e^{inθ}⟩ = (1/2π) ∫_0^{2π} e^{−imθ} (cosh β − sinh β cos θ)^{−1/2−iν} e^{inθ′} dθ`,
`e^{iθ′} = (e^{iθ} cosh(β/2) − sinh(β/2)) / (cosh(β/2) − e^{iθ} sinh(β/2))`.

`Φ^ν_{0,0}(β) = P_{−1/2+iν}(cosh β)`.

**Theorem K.** For every integer `m` and `β > 0`:
1. (*Chebyshev–Mehler formula.*) `Φ^ν_{m,m}(β) = (√2/π) ∫_0^β cos(νt) · T_{2m}(cosh(t/2)/cosh(β/2)) · (cosh β − cosh t)^{−1/2} dt`,
   with `T_{2m}` the Chebyshev polynomial. Also
   `Φ^ν_{m,m}(β) = cosh^{2m}(β/2) ₂F₁(½+m+iν, ½+m−iν; 1; −sinh²(β/2))`; it is real, even in `ν`, equal to `1` at `β = 0`.
2. (*Flat limit.*) `Φ^ν_{m,m}(β) → J_0(νβ)` as `β → 0` with `νβ` fixed, for every `m`.
3. (*Moments*, defined by Abel regularisation as in Lemma Q.) With `u_0 := sech(β/2)`,
   `Q_0^{(m)}(β) := ∫_0^∞ Φ^ν_{m,m} dν = T_{2m}(u_0) / (2 sinh(β/2)) = cos(2m·gd(β/2)) / (2 sinh(β/2))`,
   `Q_2^{(m)}(β) := ∫_0^∞ ν² Φ^ν_{m,m} dν = −[ u_0 T′_{2m}(u_0)/(8 sinh(β/2)) + T_{2m}(u_0)/(8 sinh³(β/2)) ]`,
   `Q_1^{(m)}(β) := ∫_0^∞ ν Φ^ν_{m,m} dν = (√2/π) [ G_m(0)/β − ∫_0^β (G_m(t) − G_m(0)) t^{−2} dt ]`, `G_m(t) := T_{2m}(u(t))(cosh β − cosh t)^{−1/2}`,
   and `Q_1^{(m)}(0⁺) = 1/24 − m²/2`.
4. (*Signs.*) For `m = 0` these are the scalar moments: `Q_0 > 0`, `Q_1 > 0`, `Q_2 < 0` for all `β`. For `m ≠ 0`:
   `Q_0^{(m)}` changes sign exactly at `β_j = 2 arcsinh tan((2j+1)π/(4m))`, `j = 0, …, m−1`; `Q_1^{(m)}` is **negative** at
   small `β` and changes sign; `Q_2^{(m)}` changes sign.

*Proof.* (1) By the symmetry `θ → −θ`, `Φ^ν_{m,m} = (1/π) ∫_0^π (cosh β − sinh β cos θ)^{−1/2−iν} cos(mχ) dθ` with
`χ := θ′ − θ`, `e^{iχ} = z̄/z`, `z := cosh(β/2) − sinh(β/2)e^{iθ}`. Put `e^t := |z|² = cosh β − sinh β cos θ`, `t ∈ [−β, β]`. Then
`Re z = (e^t + 1)/(2cosh(β/2))`, so `cos(χ/2) = Re z/|z| = cosh(t/2)/cosh(β/2)`, which is even in `t`; and
`dθ = e^{t/2} dt / (2(cosh β − cosh t))^{1/2}`. Hence
`Φ^ν_{m,m} = (1/π)∫_{−β}^{β} e^{−iνt} cos(mχ(t)) (2(cosh β − cosh t))^{−1/2} dt`, and `cos(mχ) = T_{2m}(cos(χ/2))`. The
hypergeometric form is the classical expression of the diagonal matrix coefficients as Jacobi functions
(checked numerically to 12 digits). (2) Rescale `t = βx`; `T_{2m}(u) → 1`. (3) As in Lemma Q, Step 4′:
`∫_0^∞ cos(νt) dν = πδ(t)`, `∫_0^∞ ν² cos(νt) dν = −πδ″(t)`, `∫_0^∞ ν cos(νt) dν = −f.p. t^{−2}`, applied to the even function
`G_m`; `u′(0) = 0`, `u″(0) = u_0/4`, `g(0) = (√2 sinh(β/2))^{−1}`, `g″(0) = (4√2 sinh³(β/2))^{−1}`. For `Q_1(0⁺)`: as `β → 0`,
`G_m(t) = √2(β²−t²)^{−1/2} − (m²/√2)(β²−t²)^{1/2} + …`; the first term gives `1/24` (Lemma Q) and the second gives
`−(m²/π)[1 + ∫_0^1 (1−√(1−x²))x^{−2} dx] = −m²/2`. (4) `T_{2m}(cos ψ) = cos 2mψ` with `ψ = gd(β/2) ∈ (0, π/2)`. ∎

Numerical values for `m = 4` (`spinning_moments.py`):

| `β` | 0.1 | 0.3 | 0.4 | 0.6 | 1.0 | 1.5 | 2.0 | 3.0 | 5.0 |
|---|---|---|---|---|---|---|---|---|---|
| `Q_0^{(4)}` | 9.21 | 1.22 | −0.046 | −1.17 | −0.733 | 0.433 | 0.341 | −0.219 | 0.021 |
| `Q_1^{(4)}` | −7.79 | −6.55 | −5.59 | −3.41 | 0.0015 | 0.759 | 0.048 | −0.166 | 0.033 |
| `Q_2^{(4)}` | −1076 | −54.5 | −24.4 | −4.40 | 3.05 | 0.877 | −0.496 | −0.068 | 0.026 |

Zeros: `Q_0^{(4)}` at `β = 0.395, 1.253, 2.386, 4.636`; `Q_1^{(4)}` at `1.000, 2.043, 3.819`; `Q_2^{(4)}` at `0.718, 1.718, 3.196`.
For comparison `Q_1^{(0)} ∈ (0, 1/24]` and `Q_2^{(0)} < 0` throughout.

5. (*Orthogonality; no positive weight.*) For every integer `m ≥ 1`, every `j = 1, …, m` and every real `ν`,

   `∫_0^∞ Φ^ν_{m,m}(β) cosh^{−2j}(β/2) sinh β dβ = 0`,

   and the integral does not vanish for `j = m + 1`. Consequently, for any finite signed measure `F` on the
   real `ν`-axis, the function `A(β) := ∫ Φ^ν_{m,m}(β) F(dν)` satisfies the `m` linear constraints
   `∫_0^∞ A(β) cosh^{−2j}(β/2) sinh β dβ = 0`, `j = 1, …, m`. **In particular `A` cannot be non-negative on `(0, ∞)` unless
   it vanishes identically.** For `m = 0` there is no such constraint (`Φ^ν_{0,0}` is not orthogonal to any
   `cosh^{−2j}`), in agreement with the scalar positivity lemma.

*Proof of 5.* Insert the Chebyshev–Mehler formula and exchange the integrals (absolutely convergent for
`j ≥ 1`): the integral equals `(√2/π) ∫_0^∞ dt cos(νt) I_j(t)`, with
`I_j(t) = ∫_t^∞ T_{2m}(cosh(t/2)/cosh(β/2)) cosh^{−2j}(β/2) (cosh β − cosh t)^{−1/2} sinh β dβ`. Substituting
`x = cosh(t/2)/cosh(β/2) ∈ (0, 1]`: `cosh β − cosh t = 2cosh²(t/2)(x^{−2} − 1)` and `sinh β dβ = −4cosh²(t/2) x^{−3} dx`, so

`I_j(t) = 2√2 cosh^{1−2j}(t/2) ∫_0^1 T_{2m}(x) x^{2j−2} (1 − x²)^{−1/2} dx`.

The last integral is half of `∫_{−1}^{1} T_{2m}(x) x^{2(j−1)} (1−x²)^{−1/2} dx`, which vanishes for `2(j−1) < 2m` by the
orthogonality of Chebyshev polynomials to polynomials of lower degree, and is non-zero for `j = m+1`. ∎

*Interpretation.* The functions `cosh^{−2j}(β/2)`, `j ≤ m`, span the weight-`(m,m)` matrix coefficients of the discrete
series of `SL(2,R)` of lowest weights `1, …, m` (those containing the `K`-type `m`); they are square-integrable
eigenfunctions of the Casimir with eigenvalues different from those of the principal series, hence
orthogonal to it. A superposition over real `ν` has no component along them. The weight-zero channel has no
discrete series. In flat-space language: the partial waves of the helicity-violating channel start at spin
four; the "missing" low spins reappear in AdS as these constraints.

Numerical check (`ads4/code/discrete_series_nogo.py`): the integrals for `m = 1, 2, 4` and `ν = 0, 0.4, 1.3, 3.7` are
`< 10^{−17}` against absolute integrals of order one; for `m = 4` they vanish for `j = 1, …, 4` (`< 10^{−30}`) and equal
`0.0044`, `0.0108` for `j = 5, 6`; for `m = 0` the corresponding integral is `0.357`.

**Off-diagonal coefficients.** `Φ^ν_{0,μ}` is proportional (with a `ν`-dependent factor of modulus one times a
ratio of Gamma functions) to the associated conical function `P^{−μ}_{−1/2+iν}(cosh β)`; a cross-check that
derivatives of the scalar harmonic function produce these coefficients:
`−(F″ + F′/sinh β)/(ν² + ¼) = Φ^ν_{1,1}` for `F = P_{−1/2+iν}(cosh β)` (12 digits). Their moments are not treated here.

## Part B. What these functions govern (symmetry argument; to be confirmed from blocks)

In the Regge limit of a CFT3 four-point function the residual symmetry acting on the transverse data is
`SO(2,1)`, and the dependence on the impact parameter of the exchange of an operator of dimension
`3/2 + iν` is an eigenfunction of its Casimir with eigenvalue `−(ν² + ¼)`. For external scalars this fixes the
kernel as `Φ^ν_{0,0}`. For external operators with spin, each pair of external operators carries a weight under
the rotation group about its transverse position, `m_L, m_R ∈ {0, ±4}` for pairs of stress tensors (the
difference of the two helicities, in the convention where `h_1 = h_2` is a flip), and covariance fixes the
kernel in the channel `(m_L, m_R)` as `Φ^ν_{m_L,m_R}(β)` **up to a factor depending on `ν` and on the three-point
structures**. This is the continuation `J → −½ + iν`, `θ → iβ` of the Wigner functions `d^J_{λμ}(θ)` of flat-space
helicity partial waves, in the same way that `P_J(cos θ) → P_{−1/2+iν}(cosh β)`. In particular the kernel for
the helicity-violating diagonal channel (`+− → +−` in the notation of arXiv:2201.06602, partial wave
`d^J_{4,4}`) is `Φ^ν_{4,4}`, **not** `P_{−1/2+iν}`, although both tend to `J_0`.

**Not derived:** the `ν`-dependent factors; the measure and weights with which each channel enters a given
functional (the spinning Regge-moment formula); uniform remainders; the threshold region; the light side.
Until the first two are derived from the spinning blocks, the consequences below are conditional.

## Part C. Consequences for the locality programme (conditional on Part B)

0. **A normalisation-independent obstruction.** Suppose that, in the regime where a heavy operator's
   contribution is controlled by the Regge moments, the action of a dispersive functional on heavy states of
   the helicity-violating diagonal channel at impact parameter `β` is of the form
   `A(β) = ∫ F(ν) Φ^ν_{4,4}(β) dν` for some weight `F` on real `ν` — *whatever* `F` is: the functional's own weight
   multiplied by any `ν`-dependent factor coming from three-point structures, normalisations or measures.
   Then by Theorem K(5) `∫_0^∞ A(β) cosh^{−2j}(β/2) sinh β dβ = 0` for `j = 1, 2, 3, 4`, so `A` takes both signs (or vanishes).
   The positive weight `cosh^{−8}(β/2) sinh β` is concentrated at `β ≈ 0.5` in AdS units, so the sign change
   happens in the curved regime and cannot be pushed to the flat regime `β ≲ 1/M`.
   **Hence, within this class, no dispersive functional is positive on all heavy states of that channel,**
   and no choice of polarisation combination helps: a positive-semidefinite quadratic form in the couplings
   `(c^{++}, c^{+−})` needs a non-negative coefficient of `|c^{+−}|²`, which is `A`. (A negative off-diagonal entry
   would not by itself violate positivity; this is a diagonal entry.) The same holds for conserved currents
   (`m = 2`, two constraints).
   What would evade it: (a) the kernel of that channel is not `Φ^ν_{4,4}` (Part B is wrong); (b) the relevant
   contributions are not of Regge-moment form at impact parameters of order the AdS radius; (c) an
   independent bound on the *magnitude* of the heavy helicity-violating contributions in the curved regime,
   replacing positivity there; (d) functionals whose weight is not supported on real `ν`.
   A valid proof of an AdS4 bound with external gravitons must use one of these; the scalar strategy
   ("positive in flat space, and positive in the curved regime by the sign of the moments") cannot work.

1. **Only the helicity-preserving channel inherits the scalar lemma.** In the `(0,0)` channel the curved-regime
   positivity is Lemma Q. In the `(4,4)` channel all three moments change sign at impact parameters of order
   the AdS radius, and the first moment — whose positivity is, in the scalar case, "the crucial distinction
   between AdS4 and the momentum cutoff" — is negative already at small impact parameter:
   `Q_1^{(4)}(0⁺) = 1/24 − 8 = −191/24`.
2. **So positivity of a flat-space graviton functional does not uplift to AdS4 channel by channel.** If a
   functional acts on the `+−` heavy states with a weight whose small-`p` expansion starts as `a_1 p + a_2 p² + …`
   with `a_1 > 0 > a_2` (the pattern that works for scalars), its curved-regime action
   `a_1 Q_1^{(4)}/M² + a_2 Q_2^{(4)}/M³` is positive at `β ≪ 1` (where `Q_2 ≈ −β^{−3}` dominates) but becomes negative
   once `β ≳ (|a_2|/(8 a_1 M))^{1/3}`, i.e. just outside the flat regime. A proof of an AdS4 bound must
   therefore either use functionals whose leading weight in the helicity-violating channel vanishes to
   higher order, or combine channels, or control these states by other means.
3. This is an **obstruction and a corrected kernel**, not a constraint on CFT data. It becomes a statement
   about a specific functional only after the spinning Regge-moment formula is derived.
