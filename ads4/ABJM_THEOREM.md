# A rigorous `N^{3/2}` law for the ABJM three-sphere partition function

Status: proof written 15 Nov 2026 and checked by one independent automated referee (hostile brief; same
model family as the author; not a human expert). Verdict: "the proof of (U), (L) and the fixed-k asymptotics
is correct; no invalid step", with errors found only in the remarks and framing, corrected below. The
referee recomputed the bounds independently and compared them with exact values and with the Airy formula.

**What is and is not new.** The `N^{3/2}` law and its coefficient are established physics
(Drukker–Mariño–Putrov, arXiv:1007.3837; Herzog–Klebanov–Pufu–Tesileanu, arXiv:1011.5487; Mariño–Putrov,
arXiv:1110.4066) and no novelty is claimed for them. What appears to be new, after searches of arXiv and
INSPIRE by two agents (not exhaustive): (a) a complete proof of the leading asymptotics as a stated theorem;
(b) explicit inequalities valid at every finite `N`. A competing route exists and should be credited: `ρ` is a
positive trace-class operator, so `Π_{i≤N} λ_i ≤ Z ≤ inf_z z^{−N} Π(1 + zλ_i)`, and any rigorous Weyl law for its
eigenvalues gives the leading term; Laptev–Schimmer–Takhtajan (arXiv:1510.00045) prove such laws for
neighbouring operators, and Kashaev–Mariño–Zakany (arXiv:1505.02243) relate the ABJM operator to them at a
physics level of rigour, but the conclusion for `Z(N,k)` is not drawn there. The lower-bound technique used
here is the standard positive-definiteness (Onsager-type) argument of rigorous statistical mechanics.

## Statement

For integers `N ≥ 1` and real `k > 0` let

`Z(N,k) := (1/N!) (4πk)^{−N} ∫_{R^N} Π_i dx_i/(2cosh(x_i/2)) · Π_{i<j} tanh²((x_i − x_j)/(2k))`.

(Relation to the physics. This is the Fermi-gas form of the localised `S³` partition function of
`U(N)_k × U(N)_{−k}` ABJM theory, Mariño–Putrov eqs. (3.2)–(3.5), after the Cauchy determinant identity
`det[1/(2cosh((x_i−x_j)/2k))] = 2^{−N} Π_{i<j} tanh²((x_i−x_j)/2k)`. It is positive and equals the standard
`Z_{S³}`, not merely its modulus (checked directly for `N = 1`, where it is `1/(4k)`). Two steps upstream are
**physics-level and outside the theorem**: supersymmetric localisation itself, and the passage from the
Kapustin–Willett–Yaakov integral — a `2N`-dimensional oscillatory integral that is not absolutely
convergent — to the Fermi-gas form by Fourier transforms and Fresnel integrals. Equivalently `Z = (1/N!) Σ_σ (−1)^σ ∫ Π ρ(x_i, x_{σ(i)})` with the Fermi-gas kernel
`ρ(x,y) = (2πk)^{−1} (2cosh(x/2))^{−1/2} (2cosh((x−y)/(2k)))^{−1} (2cosh(y/2))^{−1/2}`. The theorem is about the integral displayed above.)

**Theorem.** Let `F(N,k) := −log Z(N,k)` and `f_0 := π√(2k)/3`. For all `N ≥ 1`, `k > 0`:

(U) `F ≤ f_0 N^{3/2} + log N! + N log(8πk) − (N/2) log(2π²kN)`;

(L) for every `θ ∈ (0,1)` and `a > 0`, with `c_1 := π²k − 2∫_0^a v = 2∫_a^∞ v > 0` and `v(u) := −2 log tanh(|u|/(2k))`,
`F ≥ (2/3) √(θ c_1/2) N^{3/2} − N( 7ζ(3)k²θ/(2c_1) + v(a)/2 ) + log N! + N log(πk(1−θ))`.

Consequently, **for each fixed `k`**, `F(N,k) = (π√(2k)/3) N^{3/2} + O(N log N)` as `N → ∞` (take
`1 − θ = a = N^{−1/2}`). Explicitly (referee's evaluation):
`−(√2/(3π√k)) N (log N + 2 log 2k + 2) − (7ζ(3)/(2π²)) kN − O(√k N) ≤ F − f_0N^{3/2} ≤ (N/2) log(kN) + 0.733 N`.
The bounds (U), (L) themselves hold for all `N ≥ 1` and all real `k > 0`, but the asymptotic statement is
**not uniform in `k`**: the leading law follows only when `k/N → 0` and `k√N/log N → ∞`. In the 't Hooft
regime `k = N/λ` the error term in (L) is of the same order as the leading term (relative size `≈ 0.29 λ^{−1/2}`),
so the theorem neither captures nor contradicts the planar free energy.

With `L²/G_N` related to `N, k` as in the `AdS_4 × S⁷/Z_k` solution, `f_0 N^{3/2} = πL²/(2G_N)`, the
regularised on-shell action of Einstein gravity on Euclidean `AdS_4`. So the theorem says that the exact
free energy of the quantum theory reproduces the classical gravitational action at leading order, with
relative error `O(N^{−1/2} log N)`.

## Proof

Write `U(x) := log(2cosh(x/2))` and `W(x_1,…,x_N) := Σ_i U(x_i) + Σ_{i<j} v(x_i − x_j)`, so that
`Z = (1/N!)(4πk)^{−N} ∫ e^{−W}`. Elementary facts:
(i) `|x|/2 ≤ U(x) ≤ |x|/2 + log 2`;
(ii) `v ≥ 0` is even, and on `(0,∞)` it is decreasing and convex with `v(u) → 0` (`v' = −(2/k)/sinh(u/k)`,
which is negative and increasing);
(iii) `c := ∫_R v = π²k` and `m_1 := ∫_0^∞ u v(u) du = 7ζ(3)k²/2` (from `−log tanh t = 2 Σ_{n odd} e^{−2nt}/n`;
checked numerically for `k = 1, 2, 3`).

### Upper bound (U): Gibbs variational principle

For a probability density `ν` on `R^N` with finite entropy, Jensen's inequality gives
`log ∫ e^{−W} ≥ −∫ W ν − ∫ ν log ν`. Take `ν = Π_i f(x_i)`, `f = σ/N`,
`σ(x) := (μ − |x|/2)_+/c`, `μ := √(cN/2)`, so that `∫σ = 2μ²/c = N`. Then, by (i) and Young's inequality
`⟨f, v * f⟩ ≤ ‖v‖_1 ‖f‖_2²`,

`∫ W ν = N ∫ f U + (N(N−1)/2) ⟨f, v*f⟩ ≤ ∫ σ|x|/2 + N log 2 + (c/2) ∫ σ² = (2/3) μ N + N log 2`,

since `∫σ|x|/2 = (c/2)∫σ² = 2μ³/(3c)`. Also `f ≤ μ/(cN) = (2cN)^{−1/2}`, hence `−∫ν log ν = −N∫ f log f ≥ (N/2) log(2cN)`.
Therefore `log Z ≥ −log N! − N log(4πk) − (2/3)μN − N log 2 + (N/2) log(2cN)`, which is (U) because
`(2/3) μ N = (2/3)√(π²k/2) N^{3/2} = f_0 N^{3/2}`.

### Lower bound (L): positive definiteness and a dual certificate

*Step 1 (splitting the one-body term).* By (i), `U ≥ (1−θ)U + θ|x|/2`, so

`∫ e^{−W} ≤ (∫ e^{−(1−θ)U})^N · exp(−E_*)`, `E_* := inf_x [ θ Σ_i |x_i|/2 + Σ_{i<j} v(x_i − x_j) ]`,

and `∫ e^{−(1−θ)U} ≤ ∫ e^{−(1−θ)|x|/2} = 4/(1−θ)`.

*Step 2 (a bounded positive-definite minorant).* Put `v_1(u) := v(|u| + a)`. Since `v` is decreasing on
`(0,∞)`, `v_1 ≤ v`. As a function of `|u|`, `v_1` is continuous, convex, decreasing, with limit `0`; by
Pólya's criterion it is a positive-definite function. It is bounded by `v_1(0) = v(a)`, integrable, with
`∫ v_1 = c − 2∫_0^a v = c_1` and `∫_R v_1(u)|u| du ≤ 2 m_1`.

*Step 3 (mean-field lower bound).* Let `σ ≥ 0` be bounded with compact support and `x_1,…,x_N ∈ R`.
Positive definiteness of the bounded continuous function `v_1`, applied to the finite signed measure
`Σ_i δ_{x_i} − σ`, gives `⟨Σδ_{x_i} − σ, v_1 * (Σδ_{x_i} − σ)⟩ ≥ 0`, i.e.

`Σ_{i<j} v_1(x_i − x_j) ≥ Σ_i (v_1 * σ)(x_i) − ½⟨σ, v_1 * σ⟩ − (N/2) v_1(0)`.

*Step 4 (certificate).* Take `σ(x) := (μ − θ|x|/2)_+ / c_1` with `μ > 0` free. Since
`(μ − θ|x|/2)_+ ≥ μ − θ|x|/2` and `∫ v_1(u)|x − u| du ≤ c_1|x| + 2m_1`,

`(v_1 * σ)(x) ≥ μ − (θ/(2c_1)) (c_1|x| + 2m_1) = μ − θ|x|/2 − κ`, `κ := θ m_1/c_1`, for every `x`.

Hence `θ|x_i|/2 + (v_1*σ)(x_i) ≥ μ − κ` for each `i`, and with `½⟨σ, v_1*σ⟩ ≤ (c_1/2)‖σ‖_2² = 2μ³/(3θc_1)`:

`θ Σ|x_i|/2 + Σ_{i<j} v(x_i − x_j) ≥ Nμ − 2μ³/(3θc_1) − N(κ + v(a)/2)` for all configurations.

The right side is maximal at `μ = √(θc_1N/2)`, where `Nμ − 2μ³/(3θc_1) = (2/3)√(θc_1/2) N^{3/2}`.

*Step 5.* `Z ≤ (1/N!) (πk(1−θ))^{−N} exp(−(2/3)√(θc_1/2) N^{3/2} + N(κ + v(a)/2))`, which is (L).

*Asymptotics.* With `a = 1 − θ = N^{−1/2}`: `c_1 = c − O(N^{−1/2} log N)`, `v(a) = log N + O(1)`,
`log N! + N log(πk(1−θ)) = (N/2) log N + O(N)`, so (L) gives `F ≥ f_0 N^{3/2} − O(N log N)`; (U) gives
`F ≤ f_0 N^{3/2} + (N/2) log N + O(N)`. ∎

## Remarks

1. *Why the constant is exact.* The integrand is the Boltzmann weight of a classical one-dimensional gas
   in the linear trap `|x|/2` with a positive-definite, short-range pair repulsion of total strength
   `∫v = π²k`. For such a gas the mean-field (local density) energy `min_σ ∫σ|x|/2 + (π²k/2)∫σ²` is exact to
   leading order, and equals `f_0 N^{3/2}`. The Fermi-gas computation of Mariño–Putrov gives the same number
   as the phase-space volume `{|x| + |p| ≤ 2μ}/(4π²k)`; here it appears without any semiclassical expansion.
2. *What is not proved.* Localisation; the equivalence of the oscillatory integral with the Fermi-gas form;
   everything subleading — the `N^{1/2}` term fixed by `B = k/24 + 1/(3k)`, the logarithm (`+¼ log N` in `F`,
   i.e. `−¼ log N` in `log Z`), the Airy function, instanton corrections. At `N = 10⁴`, `k = 1` the true value of
   `F − f_0N^{3/2}` is `−80.3`, while the proved band is `[−3.6·10⁴, +5.3·10⁴]`: the theorem is consistent with
   the conjectured Airy form and far weaker than it. Nothing is proved about the 't Hooft regime, about any
   observable other than `Z_{S³}`, or about bulk dynamics. The best error the method gives is
   `+(½ + o(1)) N log N` above (a product trial state cannot do better) and, with `a ∝ N^{−1/2}/log N`,
   `−N log log N − O(N)` below. A sharper upper bound of order `N` via Araki's log-majorisation and coherent
   states is **only a sketch** (`ads4/TARGET.md` §5) and is not claimed.
3. *Consistency checks.* Author (`abjm_bounds.py`): `∫v = π²k`, `m_1 = 7ζ(3)k²/2` for `k = 1,2,3`; `Z(1,1) = 1/4`,
   `Z(2,1) = 1/(16π)` by quadrature; (L) and (U) bracket the exact `F(N,1)`, `N ≤ 4`; ratios of the bounds to
   `f_0N^{3/2}`: `(0.80, 1.21)` at `N = 10²`, `(0.976, 1.036)` at `10⁴`, `(0.997, 1.005)` at `10⁶`.
   Referee (independent scripts): exact `Z` by diagonalising `ρ` — `Z(1,k) = 1/(4k)` to `10^{−15}`; `Z(N,1)`,
   `N = 2,3,4`, equal to `1/(16π), (π−3)/(64π), (10−π²)/(1024π²)` to `≈10^{−14}`; Airy formula with
   `A(1) = log 2/4 − ζ(3)/(8π²)` agrees with exact `F` to `≈10^{−12}` for `N = 10–20`; (L) `≤ F ≤` (U) for
   `k = 1,2,3`, `N ≤ 14`, and against the Airy formula up to `N = 10⁶`; random scan of `(θ,a)` for six values of
   `k`, `N ≤ 30`: no violation; direct minimisation of `E_*` at `N = 5, 20, 60` respects the certificate bound.
   For `N ≲ 8` the best form of (L) is the trivial one from `tanh² ≤ 1`.
4. *Scope.* The proof uses only: `v ≥ 0` even, decreasing and convex on `(0,∞)` with limit `0`; `∫v < ∞`;
   `∫_0^∞ u v(u) du < ∞`; a singularity at `0` mild enough that `N v(a) = o(N^{3/2})` while `∫_0^a v → 0`; and a
   trap `U(x) = s|x| + O(1)`. For such an integral the same argument gives `F = (2/3)√(s ∫v) N^{3/2} (1 + o(1))`
   (ABJM: `s = ½`, `∫v = π²k`). Which other localised partition functions (mass deformations, flavours, other
   quivers) satisfy these hypotheses has not been examined.

## Extension: a one-parameter deformation (trial R-charges / imaginary mass) — **not yet refereed**

For `|λ| < 1/2` let `Z(N,k;λ)` be the same integral with the one-body weight `e^{λx_i}/(2cosh(x_i/2))`, and
`F(λ) := −log Z(N,k;λ)`. Put `s_± := ½ ∓ λ`, `S := 1/s_+ + 1/s_− = 4/(1 − 4λ²)`, `ℓ(x) := s_+ x_+ + s_− x_−`.

**Theorem B_λ.** For all `N ≥ 1`, `k > 0`, `|λ| < 1/2`, `θ ∈ (0,1)`, `a > 0`:

(U_λ) `F(λ) ≤ (2/3)√(2π²k/S) N^{3/2} + log N! + N log(8πk) − (N/2) log(π²kNS/2)`;
(L_λ) `F(λ) ≥ (2/3)√(2θc_1/S) N^{3/2} − N( 7ζ(3)k²θ·max(s_±)/c_1 + v(a)/2 ) + log N! − N log( S/(4πk(1−θ)) )`.

Hence, for fixed `k` and `λ`, `F(N,k;λ) = (π√(2k)/3) √(1 − 4λ²) N^{3/2} + O(N log N)`.

*Proof.* As for `λ = 0`, with `U_λ(x) := log(2cosh(x/2)) − λx`, which satisfies `ℓ ≤ U_λ ≤ ℓ + log 2`. Upper bound:
trial density `σ = (μ − ℓ)_+/c`, `∫σ = μ²S/(2c) = N`, `∫σℓ = (c/2)∫σ² = μ³S/(6c)`, so `∫Wν ≤ (2/3)μN + N log 2`, and
`f ≤ μ/(cN)`. Lower bound: certificate `σ = (μ − θℓ)_+/c_1`; since `|ℓ(x−u) − ℓ(x)| ≤ max(s_±)|u|`,
`(v_1 * ℓ)(x) ≤ c_1 ℓ(x) + 2 max(s_±) m_1`, so `θℓ + v_1*σ ≥ μ − κ` with `κ = 2θ max(s_±) m_1/c_1`;
`(c_1/2)‖σ‖² = μ³S/(6θc_1)`; the maximum over `μ` is at `μ² = 2θc_1N/S`; and `∫e^{−(1−θ)ℓ} = S/(1−θ)`. ∎

**Exact finite-`N` corollaries** (elementary). With `X := Σ_i x_i` and `E` the expectation in the `λ = 0` gas,
`Z(λ)/Z(0) = E[e^{λX}]`. Hence `λ ↦ −F(λ)` is convex (a cumulant generating function), `F` is even, and
`F(λ) ≤ F(0)`: **`F` is maximal at the undeformed point for every `N` and `k`** (Jensen, `E[X] = 0`).
`Z(1,k;λ) = 1/(4k cos πλ)`.

**Asymptotic corollaries.** (i) `N^{−3/2} F(N,k;λ) → φ(λ) := (π√(2k)/3)√(1 − 4λ²)` pointwise on `(−½, ½)`.
(ii) By convexity, the derivatives converge: `N^{−3/2} E_λ[X] → −φ′(λ)`. (iii) By the Gärtner–Ellis theorem
(`−φ` is differentiable and steep at `±½`), `X/N^{3/2}` satisfies a large-deviation principle at speed `N^{3/2}`
with rate function `I(y) = ½√(4f_0² + y²) − f_0`, `f_0 = π√(2k)/3`; near `y = 0`, `I(y) ≈ y²/(8f_0)`.

**Not proved:** `Var(X) = 4f_0 N^{3/2}(1 + o(1))`, i.e. the exchange of `N → ∞` with the second derivative at
`λ = 0`. Convexity gives it only on average over windows `|λ| ≲ h` with `h ≫ N^{−1/4}√(log N)`. A sharp lower bound
on `Var(X)` by the Helffer–Sjöstrand/Cramér–Rao inequality with the displacement field `u(x) = 2μ + |x|`, and a
sharp upper bound by the Brascamp–Lieb inequality (the Gibbs measure is log-concave on the ordered chamber;
its mean-field evaluation gives exactly `4f_0N^{3/2}`), are **sketched only** and need a law of large numbers
for the one- and two-point densities.

**Physical dictionary (checked against the sources by an independent agent, 16 Nov 2026; physics-level).**
`Z(N,k;λ)` is Nosaka's Fermi-gas integral (arXiv:1512.02862, (2.11), (2.16), `ξ = 2λ`, `η = 0`): ABJM with trial
R-charges `Δ(A_1) = ½−λ`, `Δ(B_1) = ½+λ`, `Δ(A_2) = Δ(B_2) = ½`, equivalently real mass `m_+ = 2iλ`, `m_− = 0` of
Binder–Chester–Pufu (arXiv:1906.07195, (3.1)); for real `m`, `Z(m)/Z(0) = E[e^{imX/2}]`. `φ(λ)` is
Jafferis–Klebanov–Pufu–Safdi (arXiv:1103.1181, (5.7)) on this line, the leading term of the Airy form with
`C = 2/(π²k(1+m_+²)(1+m_−²))`, and is reproduced by supergravity (Freedman–Pufu, arXiv:1302.7310).
`F(λ) ≤ F(0)` is an elementary instance of F-maximisation (arXiv:1205.4142, (4.9)).
**Stress tensor versus flavour current.** In a theory with only `N = 2` supersymmetry the second mass
derivative of `F` gives a flavour-current coefficient, `∂_m² Re F = (π²/2) τ_ff` (arXiv:1205.4142, (4.6)), **not**
`c_T`. ABJM has `N ≥ 6`, and `m_±` couple to currents in the stress-tensor multiplet; there
`c_T = (64/π²) ∂²_{m_±} F|_0` (Agmon–Chester–Pufu, arXiv:1711.07343, (26), for `N = 8`; Binder–Chester–Pufu (3.12)
for all `k`), in the normalisation where a free real scalar has `c_T = 1`. That Ward identity is physics-level
and is an input, not something proved here. Einstein gravity: `F = πL²/(2G_N)`, `c_T = 32L²/(πG_N)`, i.e.
`F = (π²/64) c_T`. The established large-`N` law on this line is known physics; only proofs and finite-`N`
inequalities are claimed as new.

Checks (`ads4/code/abjm_tilt.py`): the bounds bracket the exact `F(1,k;λ)` for three `(k,λ)` and `F(2,1;λ)`,
`λ = 0, 0.3` (quadrature; `Z(2,1;0.3) = 0.0671 > Z(2,1;0) = 0.0199`); ratios of the bounds to
`f_0√(1−4λ²)N^{3/2}` at `N = 10⁷`: within `0.2%` below and `1%` above for `λ = 0, 0.25, 0.4`.

## Stress-tensor normalisation: an exact formula and the Einstein relation `F/c_T → π²/64` — **not yet refereed**

(The identity and the convexity sandwich below were found by the agent that verified the dictionary; they
were then re-derived by hand and checked numerically, `ads4/code/abjm_ct.py`.)

**Definition (for the matrix integral).** `c_T(N,k) := (64/π²) ∂_m² F(N,k;m)|_{m=0}`, where
`Z(N,k;m)` is the integral with one-body weight `e^{imx/2}/(2cosh(x/2))`, `m` real. Equivalently
`c_T = (16/π²) Var(X)`, `X = Σ x_i`, in the undeformed gas. That this number is the stress-tensor two-point
coefficient of ABJM theory is the physics input quoted above.

**Proposition (dual form of the mass deformation).** With `E` the expectation in the undeformed gas and
`w_ij := sech²((x_i − x_j)/(2k))`, for real `m`

`Z(m)/Z(0) = sech^N(πm/2) · E[ Π_{i<j} (1 + sinh²(πm/2) w_ij)^{−1} ]`,

and the same with `sinh²(πm/2) → −sin²(πλ)` for the tilt `m = 2iλ`, `|λ| < ½`.

*Proof.* `Z = Tr_{∧^N}(V^{⊗N} B^{⊗N})` with `V = e^{imx/2}/(2cosh(x/2))`, `B = 1/(2cosh(p/2))`, `[x,p] = 2πik`. In the
momentum representation `V` has kernel `(2πk)^{−1}/(2cosh((p − p′ − πkm)/(2k)))`, so
`Z = (1/N!) ∫ Π dp_i B(p_i) (2πk)^{−N} det[1/(2cosh((p_i − p_j − πkm)/(2k)))]`. The matrix is a Cauchy matrix in
`α_i = q e^{p_i/k}`, `β_j = e^{p_j/k}`, `q = e^{−πm}`, and the Cauchy determinant gives
`det = (2cosh(πm/2))^{−N} Π_{i<j} 2 sinh² t_ij/(cosh πm + cosh 2t_ij)`, `t_ij = (p_i − p_j)/(2k)`. Dividing by the
`m = 0` value `Π tanh² t_ij · 2^{−N}`: `(1 + cosh 2t)/(cosh πm + cosh 2t) = (1 + sinh²(πm/2) sech² t)^{−1}`. At `m = 0` the
momentum-space integrand is the original gas. ∎

**Corollaries, exact for all `N`, `k`.**
1. `c_T = 16N + 32 E[Σ_{i<j} sech²((x_i − x_j)/(2k))]`; hence `16N ≤ c_T ≤ 16N²`.
2. `Z(m) > 0` for real `m`, and `F(λ) ≤ F(0) − N log sec(πλ)`.
3. The real-mass integral is again a gas of the same type: `Z(m) = sech^N(πm/2) × (gas with the same trap and
   pair potential v_m(u) = log(1 + cosh²(πm/2)/sinh²(u/2k)))`; `v_m` is positive, even, decreasing and convex on
   `(0,∞)`, and `∫ v_m = π²k(1 + m²)` (from `∫_0^∞ log(1 + cosh²a/sinh²t) dt = a² + π²/4`).

**Theorem G (the general gas).** Let `v` be even, positive, decreasing and convex on `(0,∞)`, `v → 0`,
`c := ∫v < ∞`, `m_1 := ∫_0^∞ u v < ∞`; let `ℓ(x) = s_+x_+ + s_−x_−`, `S = 1/s_+ + 1/s_−`, and `ℓ ≤ U ≤ ℓ + b_0`. For
`Z := (1/N!) ζ^{−N} ∫ exp(−Σ U(x_i) − Σ_{i<j} v(x_i − x_j))`, every `θ ∈ (0,1)`, `a > 0`, `c_1 := 2∫_a^∞ v`:

`(2/3)√(2θc_1/S) N^{3/2} − N(2θ max(s_±) m_1/c_1 + v(a)/2) − N log(S/(1−θ)) ≤ −log Z − log N! − N log ζ ≤ (2/3)√(2c/S) N^{3/2} + N b_0 − (N/2) log(cNS/2)`.

*Proof.* Verbatim the proof of Theorem B_λ, which used only these properties of `v` and `U`. ∎

**Theorem B″ (two-sided bounds on `c_T`; Einstein ratio).** Let `H(s) := log E[Π_{i<j}(1 − s w_ij)^{−1}]`, `s < 1`.
`H` is convex (each `−log(1 − s w)` is convex in `s`; Hölder), `H(0) = 0`, `H′(0) = E[Σ w_ij] = (c_T − 16N)/32`. Hence
for every `m > 0` and `λ ∈ (0, ½)`:

`[F(m) − F(0) − N log cosh(πm/2)] / sinh²(πm/2) ≤ (c_T − 16N)/32 ≤ [F(0) − F(λ) + N log cos(πλ)] / sin²(πλ)`.

Inserting the bounds of Theorem G for the three gases (undeformed; tilt `λ`; real mass `m`, for which
`F(m) − N log cosh(πm/2)` is the free energy of the `v_m`-gas) and letting `λ, m → 0` slowly (`λ, m ~ N^{−1/8}`):

`c_T(N,k) = (64√(2k)/(3π)) N^{3/2} (1 + o(1))`, and therefore **`F(N,k)/c_T(N,k) → π²/64`** as `N → ∞` at fixed `k`.

*Proof of the asymptotics.* Upper side: `F(0) − F(λ) ≤ f_0N^{3/2}(1 − √(1−4λ²)) + O(N log N)`, and
`(1 − √(1−4λ²))/sin²(πλ) → 2/π²`. Lower side: `F_gas(v_m) − F(0) ≥ f_0N^{3/2}(√(1+m²) − 1) − O_m(N log N)`, and
`(√(1+m²) − 1)/sinh²(πm/2) → 2/π²`. So `(c_T − 16N)/32 = (2f_0/π²) N^{3/2}(1 + o(1))`. ∎

No exchange of the large-`N` limit with a derivative, and no analytic continuation in the mass, is used:
convexity of `H` converts values of `F` at finite deformation into bounds on the derivative at zero.

**Checks (`ads4/code/abjm_ct.py`).** `N = 2`, `k = 1`: the identity at `m = 0.7` (`0.296081869` both sides);
`(16/π²)Var(X) = 37.3333331`, `16N + 32E[w] = 37.3333333` (literature value `112/3`). The agent additionally
reproduced the published values at `(N,k) = (2,2), (3,1)`. `∫v_m = π²k(1+m²)` to 12 digits. Evaluated two-sided
bounds divided by the leading term: `(0.38, 2.06)` at `N = 10⁴`, `(0.72, 1.31)` at `10⁶`, `(0.91, 1.10)` at `10⁸`,
`(0.93, 1.05)` at `10^{10}` — convergence is slow (rate `≈ N^{−1/4}`).

**What this is and is not.**
* It is a theorem about the matrix integral: the quantity `(64/π²)∂_m²F` has the asymptotics that Einstein
  gravity predicts for `c_T` given `F`, with explicit finite-`N` bounds and the exact inequality
  `16N ≤ c_T ≤ 16N²`. In bulk language: the Newton constant read from the on-shell action and the one read
  from the graviton two-point normalisation agree at leading order.
* The physics was not in doubt (the leading `c_T` and the Airy form are known). New, as far as searched: the
  proof, and possibly the finite-`N` inequalities; the dual representation of the mass deformation is in the
  literature (Nosaka; Honda et al.), and the formula `c_T = 16N + 32E[Σw]` may well be implicit there — **no
  search was made**.
* It gives **no new constraint on the gravitational coupling**: the error terms exceed every known
  correction. It says nothing about `⟨TTT⟩`, i.e. about whether the bulk theory is Einstein rather than
  higher-derivative beyond this one relation, and nothing about locality.
* Inputs not proved: localisation; the `N ≥ 6` Ward identity relating the mass derivative to `c_T`.
