# Harmonic kernels for spinning Regge exchange in CFT3, and why scalar positivity does not carry over

Status (16 Nov 2026, after hostile review). **Part A (harmonic analysis): established**, and largely
classical. **Part B as first written was wrong** in its assignment of kernels to physical channels, and
**the physical conclusions drawn from it (the former Part C items 0–2, "no positive functional in the
helicity-violating channel") are retracted.** The corrected dictionary below is the referee's and is itself
not yet derived from spinning blocks. Nothing here is a locality bound, and nothing here is an obstruction
to one.

Referee's three separate verdicts. *Mathematical correctness:* Theorem K items 1–3 and 5 established; item 4
after restricting to `|m| ≥ 2`. *Physical applicability:* the kernel family is the right one for spinning
Regge exchange, but `Φ_{4,4}` does not multiply a diagonal coupling; the obstruction claimed does not exist.
*Novelty:* items 1, 2, 5 are classical (Jacobi functions; Koornwinder's generalised Mehler–Dirichlet kernel,
*Ark. Mat.* 13 (1975) and the 1984 survey; orthogonality of continuous and discrete spectrum for
type-`m` spherical functions: Flensted-Jensen 1977, Shimeno 1994; equivalently Landau levels on `H²`). The
closed-form Abel moments and `Q_1(0⁺) = 1/24 − m²/2` were not found in the literature: "elementary corollaries,
apparently new but minor". No AdS4 application was found.

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
   `Q_0^{(m)}` changes sign exactly at `β_j = 2 arcsinh tan((2j+1)π/(4m))`, `j = 0, …, m−1`; `Q_1^{(m)}` is negative at
   small `β`; for `|m| ≥ 2`, `Q_1^{(m)}` and `Q_2^{(m)}` each change sign `m − 1` times. (For `|m| = 1` they do not:
   `Q_2^{(1)} = −(1 + 3sinh²(β/2))/(8 sinh³(β/2) cosh²(β/2)) < 0` and `Q_1^{(1)} < 0` throughout — referee's correction.)

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
   it vanishes identically.** (Measure class: finite signed measures, using the uniform bound
   `|Φ^ν_{m,m}(β)| ≤ P_{−1/2}(cosh β)`, which follows from `|T_{2m}| ≤ 1`. For Abel-regularised or polynomial weights `A` is a
   distribution with a singular part at `β = 0` and the pointwise statement fails: e.g.
   `∫_0^∞ Q_1^{(m)} cosh^{−2j}(β/2) sinh β dβ = −1`, the compensating `+1` being a delta function at the origin. Sharper
   form: `{cosh^{−2j}}` is a Chebyshev system, so `A` has at least `m` sign changes for finite measures. — Referee.) For `m = 0` there is no such constraint (`Φ^ν_{0,0}` is not orthogonal to any
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

6. (*Complex spectral parameter; a sharp threshold.* Added after the review; **not refereed**; presumably also
   classical.) The orthogonality of item 5 with a given `j ≤ m` holds for all complex `ν` in the strip
   `|Im ν| < j − ½` (the integral is analytic in `ν` there and vanishes on the real axis). Hence no superposition
   over `|Im ν| < m − ½` is non-negative. Conversely, for `ν = iκ` with `κ ≥ m − ½`, `Φ^{iκ}_{m,m}(β) > 0` for all `β`
   (Pfaff's transformation gives `cosh^{2m−2b}(β/2) ₂F₁(κ + ½ − m, b; 1; tanh²(β/2))`, `b = ½ + m + κ`, a series with
   non-negative terms); at `κ = m − ½` it equals `cosh^{−2m}(β/2)`. Numerically (`strip_threshold.py`): for `m = 4` the
   kernel changes sign for `κ = 0, …, 3.4` and is positive for `κ = 3.5, 3.6, 4.5, 6`; for `m = 2` the threshold is
   `1.5`; for `m = 0` every `κ ≥ 0` is positive.

**Off-diagonal coefficients and ladder relations (referee).** `Φ^ν_{0,μ} = (−1)^μ (½ − iν)_μ P^{−μ}_{−1/2+iν}(cosh β)` (the
factor is not of modulus one). With `s = ½ + iν`:
`dΦ_{m,n}/dβ = ((s+n)/2) Φ_{m,n+1} + ((s−n)/2) Φ_{m,n−1}`, and similarly on the left index. Hence four raisings on each
side of the scalar harmonic function give `Π_{k=0}^{3} (ν² + (k+½)²) · Φ^ν_{4,4}`; the zeros of this factor are the
discrete-series points `ν = ±i(k+½)`. A cross-check done here:
`−(F″ + F′/sinh β)/(ν² + ¼) = Φ^ν_{1,1}` for `F = P_{−1/2+iν}(cosh β)` (12 digits).

## Part B. What these functions govern (corrected after review; still to be derived from blocks)

**What is sound.** In the Regge limit of a CFT3 four-point function the residual symmetry acting on the
transverse data is `SO(2,1)`; components of weight `(m_L, m_R)` are `Φ^ν_{m_L,m_R}(β)` times polynomial or Gamma
factors in `ν`, regularity at `β = 0` giving uniqueness. This agrees with the differential-operator form of
spinning Regge theory (arXiv:1707.07689, 1705.02934).

**What was wrong.** The transverse group acts on the Regge pairs `(1,3)` and `(2,4)` ("operators 3 and 4 become
null from 1 and 2"). The weights are therefore the **helicity flips** `m_L = h_1 − h_3`, `m_R = h_2 − h_4`. They are
*not* the `s`-channel labels `λ = h_1 − h_2`, `μ = h_3 − h_4` of the Wigner functions `d^J_{λμ}`, which belong to the heavy
state's little group; the "continuation `J → −½ + iν`, `θ → iβ`" used in the first version swaps label and
argument (`J ~ b`, `θ ~ p`, whereas `ν ~ p`, `β ~ b`) and is not a valid dictionary. Consequently (notation of
arXiv:2201.06602, (2.14)–(2.16)):

| Process | Coupling | Regge-pair weights | Kernel |
|---|---|---|---|
| `++ → ++` | `|c^{++}|²` (diagonal) | `(0, 0)` | `P_{−1/2+iν}` |
| `+− → +−` (MHV, helicity-preserving, forward non-vanishing) | `|c^{+−}|²` (diagonal) | `(0, 0)` | `P_{−1/2+iν}` (plus a `(−1)^J Φ_{4,−4}` term from `+− → −+`) |
| `++ → −−` (all-plus) | `(c^{++})²` (**interference**, phase not fixed by unitarity) | `(4, 4)` | `Φ_{4,4}` |
| `++ → +−` | `c^{++} c̄^{+−}` (interference) | `(4, 0)` | `Φ_{4,0} ∝ P^{−4}` |

So `Φ_{4,4}` multiplies an off-diagonal, interference amplitude. Positivity of a quadratic form in the
couplings requires only `A_{00} ≥ |A_{44}|`, not `A_{44} ≥ 0`. A consistency check that should have stopped the
first version: the eikonal phase shift of general relativity is helicity-independent and positive, which
is impossible for any real-`ν` superposition of `Φ_{4,4}` by Theorem K(5); so the helicity-preserving diagonal
kernel cannot be `Φ_{4,4}`.

**Not derived (unchanged):** the spinning Regge-moment formula itself — arXiv:2311.04271 leaves the Regge
regime of spinning four-point functions "to future work", and no such formula was found in the literature —
with its `ν`-dependent factors, measure, remainders, threshold region and light side.

## Part C. Consequences (corrected)

**Retracted.** The first version claimed that no dispersive functional can be positive on heavy states of
the "helicity-violating diagonal channel", and that `Q_1^{(4)}(0⁺) = −191/24` shows the scalar positivity
mechanism fails for gravitons. Both claims rested on the wrong channel assignment and are withdrawn. In
addition: (i) the remark that the sign change "cannot be pushed to the flat regime" was invalid reasoning
(for weights `f(ν/M)` the region `β ≲ 1/M` contributes to the constraint at the same order); (ii) the value
`−191/24` is not robust — it sits at relative order `1/ν²`, exactly where the normalisation factors act (e.g.
`|Γ(9/2+iν)/Γ(1/2+iν)|² = ν⁸(1 + 21/ν² + …)`; and with the Plancherel measure `ν tanh πν` the first moment of the
scalar kernel is identically zero rather than `1/24`).

**What does follow, conditionally on the corrected dictionary and on a spinning Regge-moment formula.**
1. For sum rules whose heavy densities are diagonal (`|c^{++}|²`, `|c^{+−}|²`: the helicity-preserving, MHV-type
   components), the leading curved-regime kernel is the scalar one, so the positivity input is Lemma Q
   (`POSITIVITY_LEMMA.md`): `Q_1 > 0`, `Q_2 < 0`. The scalar lemma therefore **does** carry over to those
   components. This agrees with the earlier, independent inference in `EINSTEIN_OBSTRUCTIONS.md` §4.2.
2. For functionals that also use helicity-changing components, positivity is a matrix condition, e.g.
   `A_{00}(β) ≥ |A_{44}(β)|` for the `c^{++}` sector, where `A_{44}` is a superposition of `Φ_{4,4}`. Theorem K gives the
   explicit functions needed to test such a condition; by item 5, `A_{44}` necessarily oscillates in sign,
   which is harmless.
3. A generic position-space functional acting on weight-`m` data has components along the discrete series
   at `ν = ±i(k − ½)`, `k ≤ m` (the zeros of the Pochhammer factor above); restricting to real `ν` is natural only
   for `m_L m_R ≤ 0`.

**What would settle the dictionary:** the heavy-block limit of one spinning block at fixed `β` (currents or
stress tensors, `d = 3`) in the helicity basis; or the eikonal test.
