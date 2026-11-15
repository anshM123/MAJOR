# A rigorous `N^{3/2}` law for the ABJM three-sphere partition function

Status: proof written 15 Nov 2026; numerical consistency checks in `ads4/code/abjm_bounds.py`; **not yet
refereed**. Novelty: the *result* is accepted physics (Drukker–Mariño–Putrov, arXiv:1007.3837;
Herzog–Klebanov–Pufu–Tesileanu, arXiv:1011.5487; Mariño–Putrov, arXiv:1110.4066). What is claimed here is
a complete proof with explicit error bounds. The literature audit found no such proof (it found the
statement "their derivation remains largely conjectural", arXiv:2602.10196), but the search of the
mathematical literature was not exhaustive.

## Statement

For integers `N ≥ 1` and real `k > 0` let

`Z(N,k) := (1/N!) (4πk)^{−N} ∫_{R^N} Π_i dx_i/(2cosh(x_i/2)) · Π_{i<j} tanh²((x_i − x_j)/(2k))`.

(This is the localised `S³` partition function of `U(N)_k × U(N)_{−k}` ABJM theory after the Cauchy
determinant identity; equivalently `Z = (1/N!) Σ_σ (−1)^σ ∫ Π ρ(x_i, x_{σ(i)})` with the Fermi-gas kernel
`ρ(x,y) = (2πk)^{−1} (2cosh(x/2))^{−1/2} (2cosh((x−y)/(2k)))^{−1} (2cosh(y/2))^{−1/2}`. The localisation
argument that leads to this integral is **not** part of the theorem.)

**Theorem.** Let `F(N,k) := −log Z(N,k)` and `f_0 := π√(2k)/3`. For all `N ≥ 1`, `k > 0`:

(U) `F ≤ f_0 N^{3/2} + log N! + N log(8πk) − (N/2) log(2π²kN)`;

(L) for every `θ ∈ (0,1)` and `a > 0` with `c_1 := π²k − 2∫_0^a v > 0`, where `v(u) := −2 log tanh(|u|/(2k))`,
`F ≥ (2/3) √(θ c_1/2) N^{3/2} − N( 7ζ(3)k²θ/(2c_1) + v(a)/2 ) + log N! + N log(πk(1−θ))`.

Consequently `F(N,k) = (π√(2k)/3) N^{3/2} + O_k(N log N)` as `N → ∞` (take `1 − θ = a = N^{−1/2}`).

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
2. *What is not proved.* The next terms (`N^{1/2}` with coefficient fixed by `B = k/24 + 1/(3k)`, the
   `−¼ log N`, the Airy function) — the error term here, `O(N log N)`, is larger than all of them. A sharper
   quantum upper bound `F ≤ f_0 N^{3/2} + (2 log 2 + √(2k)) N` follows from Araki's log-majorisation
   `e^{−(U+T)} ≺ e^{−U/2}e^{−T}e^{−U/2}` and a coherent-state trial density matrix; it is sketched in
   `ads4/TARGET.md` §5 and not written out.
3. *Consistency checks* (`abjm_bounds.py`): `Z(1,1) = 1/4` and `Z(2,1) = 1/(16π)` by quadrature; (L) and (U)
   bracket the exact values of `F(N,1)` for `N = 1,…,4` (using the published closed forms
   `1/4, 1/(16π), (π−3)/(64π), (10−π²)/(1024π²)`, the last two taken from the literature without
   re-derivation); the ratio of the two bounds to `f_0N^{3/2}` is `(0.80, 1.21)` at `N = 10²`,
   `(0.976, 1.036)` at `N = 10⁴`, `(0.997, 1.005)` at `N = 10⁶`.
4. *Scope.* The proof uses only: positivity, monotonicity and convexity of `v` on the half-line, `∫v`, and
   `U(x) = |x|/2 + O(1)`. It therefore applies verbatim to any matrix integral of this "ideal Fermi gas"
   form with a convex decreasing pair potential and a linear trap, e.g. with real `k`, and gives
   `F = (2/3)√(∫v · N³ · s/2)`-type laws for traps of slope `s`. Which other localised partition functions
   (mass deformations, `N_f` flavours, other quivers) satisfy the hypotheses has not been examined.
