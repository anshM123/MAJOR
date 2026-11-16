# Heavy-state positivity of the cubic-coupling functional `F_{g3}` of arXiv:2201.06602, and the bound on `|ĝ₃|²`

Status: written 16 Nov 2026 by an automated analyst; **not refereed**. Computer-assisted proof; code in
`ads4/code/flatbound/`, driver `bash run_all.sh` (about 3.5 min; prints one PASS/FAIL line per certificate and a final
verdict; current result: all PASS). Units `M = 1` throughout unless `M` is written.

## 0. Results

Notation (justified in §1–2). `ψ₂(p) = (1−p)³ p (65+155p+47p²)`, `ψ₃(p) = 2(1−p)³(81−1927p+4601p²)`. For a heavy state of
mass `m ≥ 1`, spin `J` and helicity channel `h ∈ {++, +−}` put `u = 1/m²`, `𝒥 = J(J+1)`, `σ₊₊ = +1`, `σ₊₋ = −1`,
`D₊₊ = P_J` (Legendre), `D₊₋ = d̃^J_{44} = P^{(0,8)}_{J−4}` (Jacobi), and

`F[m,J,h] = ∫₀¹ dp [ψ₂(p)(2m²−p²)/m⁶ + σ_h ψ₃(p)/m⁶] D_h(1−2p²/m²) + Φ_h(J)/m¹⁰`   (0.1)  (argument `cos θ = 1−2p²/m²`, `sin(θ/2) = p/m`)

* **literal functional** (eq. (3.3a) as printed, with `m_IR → 0`): `Φ₊₊ = 10𝒥 − 5`, `Φ₊₋ = 10𝒥 − 205`;
* **g5-free completion** `F_{g3} + 5 B₅^{(1)}|_{p=0}` (see §1.4 for why it is needed): `Φ₊₊ = 10𝒥`, `Φ₊₋ = 10𝒥 − 210`.

`F[m,J,h]` is the coefficient of `|c^h_{J,m²}|²` in the heavy-state action (a common positive factor
`16(2J+1)dm²/m²` dropped). Allowed spins: `J ∈ {0,2,4,…}` for `++`, `J ∈ {4,5,6,…}` for `+−`.

**Theorem 1 (heavy-state positivity; the paper's first claim below (3.3)).** For every `m ≥ M`, every allowed `J`,
both channels and both variants, `F[m,J,h] > 0`.

**Theorem 2 (exact low-energy identity).** Under hypotheses H1–H5 of §4, for every `0 < ε = m_IR/M < 1`, with the
g5-free functional `F^ε` (the `p`-integrals over `[ε,1]`):

`π G k(ε) |ĝ₃|² = 16 π G L(ε) − ⟨F^ε⟩_heavy − 𝓜`, `L(ε) = ∫_ε¹ ψ₂/p² dp = 65 log(1/ε) − 1171/15 + ℓ(ε)`,

`k(ε) = 20 − 2∫_ε¹ψ₂p⁶ + 2∫_ε¹ψ₃p⁴`, `k(0) = 330151/12012 = 27.48509…`, `ℓ(ε) = 40ε + 223ε²/2 − 259ε³/3 + 7ε⁴/2 + 47ε⁵/5`,
and `𝓜 ≥ 0` the light-matter term (§4.2). For the literal functional the right side contains in addition `−5 g₅`.

**Theorem 3 (bound).** Under H1–H5 and the infrared hypothesis H-IR (§4.3), for `0 < m_IR ≤ M/100`:

`|ĝ₃|² M⁸ ≤ 16 L(ε)/k(ε) ≤ A log(M/m_IR) + B`, **`A = 37.8387`, `B = −45.2059`** (rounded up; `A = 1040/k_lo`,
`B = 16(ℓ(1/100) − 1171/15)/k_hi`, exact rationals in `cert_low.py`); as `ε → 0` the constants tend to
`1040/k(0) = 37.8387` and `−16·1171/(15 k(0)) = −45.4452`. Light spin-0/spin-2 matter can only lower the right side.
This is exactly the paper's (3.4a) (`37.8 log(M/m_IR) − 45.4 − F_matter`). Without H-IR the statement holds with
`+ Δ_IR(ε)` on the right, `Δ_IR = ⟨[F^ε]₋⟩/(πG k(ε)) ≥ 0`, an explicit heavy-state quantity (§4.3).
For the **literal** functional one gets instead `|ĝ₃|²M⁸ ≤ A log + B − 5 g₅ M¹⁰/(πG k)` (`g₅` has either sign).

How positivity was established (Theorem 1), `b̂ := (2J+1)/m` (impact parameter):

| region | method | status |
|---|---|---|
| (P) `J < 400`, all `m ≥ 1` | `m⁴F` is an exact rational polynomial in `u` of degree `≤ J+1`; positivity on the closed interval `u∈[0,1]` by Descartes' rule of signs with dyadic bisection (exact integers) | certified, `cert_P.py` (596 polynomials) |
| (C) `J ≥ 400`, `1 ≤ m ≤ 2` | `Φ ≥ 10𝒥−210` beats crude bounds `|P_J| ≤ 1`, `|d̃^J_{44}| ≤ min(binom(J+4,8),(1−p²)^{−4})` | analytic + interval constants, `cert_C.py` |
| (M) `m ≥ 2`, `b̂ ≥ 120` (all `J`) | exterior Mehler formula, four integrations by parts: exact `1/b̂`, `1/b̂³` terms + explicit `O(b̂⁻⁴)` remainder; `+−` reduced to Legendre by a Clebsch–Gordan identity | analytic + interval constants, `cert_M.py` |
| (S) `J ≥ 400`, `m ≥ 2`, `b̂ ≤ 120` | continuous extension in `(b̂,u)`; exact double Taylor coefficients + Cauchy bounds on a polydisc from a Mehler–Dirichlet representation valid for complex arguments; 604 `b̂`-cells per channel, adaptive boxes | certified, `cert_S.py` |

**Findings.** (i) Positivity is true but **nearly marginal** in the `+−` channel: the exact value at `m = M`, `J = 4` is
`m¹⁰F = 19/315 ≈ 0.060` for the g5-free functional (`5.060` literal), against `≈ 12` for generic states; and at
`m ≈ 21.79 M`, `J = 245` (`b̂ ≈ 22.5`) one has `m⁴F ≈ 2.67·10⁻⁵`, i.e. `b̂³m⁴F ≈ 0.31` against `≈ 80` for the flat
asymptotics (floating point, not used). The second near-zero is invisible in the paper's `m¹⁰F` normalisation (value
≈ 2.9·10³ there). (ii) The functional as printed in (3.3a) does **not** by itself give (3.4a): its low-energy side contains
`−5g₅` (from `∂_{p²}B₄|₀`, eq. (2.32a)); (3.4a) follows only after adding `5B₅^{(1)}|₀` (or assuming `g₅ ≥ 0`). The completed
functional is also positive on all heavy states (Theorem 1) and gives exactly the paper's constants. (iii) The IR
cutoff destroys positivity at large impact parameter (as the paper says); a finite bound therefore needs a physics
hypothesis about that region (H-IR); it is not a theorem.

## 1. The sum rules (arXiv:2201.06602, read from the PDF)

### 1.1 Amplitudes, partial waves, unitarity
`M(1⁺2⁻3⁻4⁺) = ⟨23⟩⁴[14]⁴ f(s,u)` (2.1), `s+t+u = 0`, `u = −(p₁+p₃)²` (2.4), `f(s,u) = f(u,s)` (2.5). Low energy (2.7a):
`f_low = 8πG/(stu) + 2πG su|ĝ₃|²/t + g₄ + g₅t + g₆t² − g₆′su + … + f_matter`. Partial waves (2.12), unitarity (2.14)
`Im a_J^{+−+−} = |c^{+−}_{J,s}|²`, `Im a_J^{++−−} = |c^{++}_{J,s}|²`, and (2.16): at `s = m²`
`Im f(m²,−p²) = (16π/m⁸) Σ_{J≥4}(2J+1)|c^{+−}_{J,m²}|² d̃^J_{4,4}(1−2p²/m²)` (s-channel, `J ≥ 4`, all integers), at `t = m²`
`Im f(−p²,p²−m²) = (16π/m⁸) Σ_{J≥0 even}(2J+1)|c^{++}_{J,m²}|² d̃^J_{0,0}(1−2p²/m²)`. App. A (A.1)–(A.3): `d̃` is `d` without the
factors `((1+x)/2)^{(h+h′)/2}((1−x)/2)^{(h−h′)/2}`; so `d̃^J_{00} = P_J`, `d̃^J_{44}(x) = ₂F₁(4−J, J+5; 1; (1−x)/2) = P^{(0,8)}_{J−4}(x)`
(normalisation `N^J_{44} = 1`), and `d^J_{44} = ((1+x)/2)⁴ d̃^J_{44}`.

### 1.2 Fixed-`u` sum rules and their heavy side
(2.26a) `B_k^{(1)}(p²) = ∮ ds/(4πi) (2s−p²)/[s(s−p²)]^{(k−2)/2} [f(s,−p²) + f(p²−s,−p²)] = 0`, `k ≥ 2` even;
(2.26b) `B_k^{(1)} = ∮ ds/(4πi) [s(s−p²)]^{−(k−3)/2} [−f(s,−p²) + f(p²−s,−p²)] = 0`, `k ≥ 3` odd; valid after smearing
(2.21)–(2.25). Contour deformation (2.28)–(2.29): `−B|_low = B|_high`. Inserting (2.16): the first term has the
`s`-channel cut (`+−` states), the second, `f(p²−s,−p²) = f(−p², p²−s)`, the `t`-channel cut at `t = s` (`++` states).
Hence, with the average `⟨X⟩ = 16Σ_J(2J+1)∫_{M²}^∞ (dm²/m²) X` of (2.34), `x = 1−2p²/m²`:

* (2.33) `B_k|_high = ⟨(2m²−p²)(|c⁺⁺|²d̃₀₀(x) + |c⁺⁻|²d̃₄₄(x)) / (m^{k+4}(m²−p²)^{(k−2)/2})⟩`, `k` even;
* (2.35a) `B_k|_high = ⟨(|c⁺⁺|²d̃₀₀(x) − |c⁺⁻|²d̃₄₄(x)) / (m^{k+3}(m²−p²)^{(k−3)/2})⟩`, `k` odd.

I re-derived both from (2.26) and (2.16) (the `B₃` term has **opposite signs** in the two channels because of the
`−f(s,…)` in (2.26b)). Per state (dropping `16(2J+1)dm²/m²`): `C₂ = (2m²−p²)D_h/m⁶`, `C₃ = σ_h D_h/m⁶`,
`C₄ = (2m²−p²)D_h/(m⁸(m²−p²))`, `C₅ = σ_h D_h/(m⁸(m²−p²))`. With `D₊₊′(1) = 𝒥/2`, `D₊₋′(1) = (J−4)(J+5)/2`:
`∂_{p²}C₄|₀ = (1−2𝒥)/m¹⁰` (`++`), `(41−2𝒥)/m¹⁰` (`+−`); `C₅|₀ = σ_h/m¹⁰`. This gives (0.1).

### 1.3 Low-energy side (re-derived, `cert_low.py`)
`−B|_low = ½ Σ Res` of the integrand at the poles inside the arcs. From (2.7a):
`−B₂ = 16πG/p² + 2πG|ĝ₃|²p⁶`, `−B₃ = −2πG|ĝ₃|²p⁴` ((2.31a,b)), `−B₄ = 2g₄ + (4πG|ĝ₃|² + g₅)p² + (g₆+g₆′)p⁴` ((2.32a)),
`−B₅ = g₅ + (g₆−g₆′)p²` ((2.32b)) — all checked symbolically. Matter (B.1) gives `−B₂ = (p²−2m_ℓ²)(a₀m_ℓ⁴ + a₂(m_ℓ⁴ − 6m_ℓ²p² + 6p⁴))/m_ℓ⁴`,
`−B₃ = −(…)/m_ℓ⁴`, `−B₄ = 12a₂p²/m_ℓ⁴`, `−B₅ = 0` (`a_J = |g_J(m_ℓ)|²`), i.e. (3.8) up to the overall sign convention.

### 1.4 The functional; bracketing; the `g₅` term
(3.3a): `F_{g3} = ∫_{m_IR}^M dp (1−p)³[p(65+155p+47p²)B₂ + 2(81−1927p+4601p²)B₃] − 5∂_{p²}B₄|_{p=0}`. The PDF text does not
show the bracket glyph; the reading "`(1−p)³` multiplies both polynomials" is fixed **by the paper's own numbers**:
with it the low-energy side gives `A = 1040/k(0) = 37.8387`, `B = −45.4452` (paper (3.4a): 37.8, −45.4) and the matter
polynomials `2591/210 m² + 1/18` and `2591/210 m⁶ − 239/18 m⁴ − 64573/1540 m² + 330151/4004` of (3.9) exactly; the other
reading gives `A = 0.728` and `11711/9` instead of `1/18`.

The low-energy action of the printed functional is (exact, `cert_low.py`)
`16πG L(ε) − πG k(ε)|ĝ₃|² − 5g₅ + (matter)`. The `−5g₅` is absent from (3.4a). Since `g₅ = ⟨(|c⁺⁺|²−|c⁺⁻|²)/m¹⁰⟩`
((2.41)) has no sign, (3.4a) does not follow from (3.3a) as printed. Adding `5B₅^{(1)}|_{p=0}` removes it (its low
side is `5g₅`, heavy side `5σ_h/m¹⁰`), changes nothing else on the low side, and keeps heavy positivity (Theorem 1).
All certificates use `Φ₊₊ = 10𝒥−5`, `Φ₊₋ = 10𝒥−210`, the smaller of the two variants in each channel.

## 2. The heavy action in normal form

`𝔉_h(J,u) := m⁴F[m,J,h] = ∫₀¹ Ψ_h(p;u) D_h(1−2up²) dp + Φ_h u³`, `Ψ_h = ψ₂(2−up²) + σ_h uψ₃`.

**Lemma P (exact polynomial).** With `μ₂_j = ∫ψ₂p^j`, `μ₃_j = ∫ψ₃p^j` (rationals), `e₊₊(i) = i(i−1)`, `e₊₋(i) = (i+3)(i+4)`:
`𝔉_h = Σ_{k≥0} (−u)^k Π_{i≤k}(𝒥 − e_h(i))/(k!)² · [2μ₂_{2k} + u(−μ₂_{2k+2} + σ_hμ₃_{2k})] + Φ_h u³` (finite for integer `J`).
*Proof:* `P_J(1−2y) = ₂F₁(−J,J+1;1;y)`, `P^{(0,8)}_{J−4}(1−2y) = ₂F₁(4−J,J+5;1;y)`, Pochhammer products as stated. ∎ [proved; checked
against quadrature to 10⁻¹¹, `sanity_numeric.py`]. Example: `𝔉₊₋(4,u) = 2591/210 − 41u/18 − 10u³`, positive on `[0,1]`
with minimum `19/315` at `u = 1`.

**Continuous extension.** With `X = b̂²/4`, `β := 𝒥u = X − u/4`, `λ_i = ¼ + e_h(i)`:
`𝔉_h = G_h(b̂,u) + 10Xu² − (5/2 + c_h)u³`, `G_h = Σ_n (−1)^n/(n!)² Π_{i≤n}(X − λ_i u)·(T0_n + T1_n u)`, `c₊₊ = 5`, `c₊₋ = 210`.
`G_h` is entire in `b̂` and analytic in `|u| < 1`; at physical points it equals the integral. Flat slice:
`G_h(b̂,0) = 2∫ψ₂J₀(b̂p)dp = 2f̂(b̂)` in both channels (the scaling limit (3.5)).

**Lemma C (Clebsch–Gordan reduction of the `+−` kernel).** For all complex `J` with `c_j` finite and `|y| < 1`:
`(1−y)⁸ ₂F₁(4−J,J+5;1;y) = Σ_{j=−4}^{4} c_j(J) ₂F₁(−J−j, J+j+1; 1; y)`, with `c_j = ⟨J4;4−4|J+j,0⟩²` the explicit rational
functions in `common.py` (e.g. `c₀ = 35(J−3)(J−2)(J+3)(J+4)/(8(2J−3)(2J−1)(2J+3)(2J+5))`); `c_j ≥ 0` for `J ≥ 4`, `Σc_j = 1`.
*Proof:* for integer `J ≥ 4`, the CG product rule `d^J_{44}d⁴_{−4,−4} = Σ_L ⟨J4;4−4|L0⟩² d^L_{00}` with `d⁴_{−4,−4} = cos⁸(θ/2) = (1−y)⁴`
gives the identity with weights `CG²`. The coefficients of `y⁰…y⁸` determine the weights uniquely (the 9×9 matrix is
triangular-Vandermonde in `L(L+1)`; its determinant, computed symbolically, has no zero at integers `J ≥ 4`), and
`cert_cg.py` checks symbolically that the `c_j(J)` solve these 9 equations, so `CG² = c_j(J)` for all integer `J ≥ 4`. Each
coefficient of `y^k` of LHS−RHS is then a rational function of `J` vanishing at all integers `J ≥ 4`, hence `≡ 0`. ∎
[proved; `cert_cg.py` also checks `y⁹…y¹⁴` identically in `J`; numerically at complex `J`, `sanity_numeric.py`]
Consequently `d̃^J_{44}(1−2y) = (1−y)^{−8} Σ_j c_j P_{J+j}(1−2y)`.

**Classical inputs (quoted).** (Q1) Mehler–Dirichlet, DLMF 14.12.1 with `μ = 0`: `P_ν(cos θ) = (√2/π)∫₀^θ cos((ν+½)φ)(cos φ − cos θ)^{−1/2}dφ`,
`0 < θ < π`, all complex `ν`. (Q2) exterior form for integer `L`: `P_L(cos θ) = (√2/π)∫_θ^π sin((L+½)φ)(cos θ − cos φ)^{−1/2}dφ`
(Whittaker–Watson §15.231; checked numerically to 10⁻¹⁷). (Q3) `|d^J_{mm′}| ≤ 1`; `max_{[−1,1]}|P_n^{(0,8)}| = binom(n+8,n)`
(Szegő (7.32.2)). (Q4) CG product rule. Equation numbers from memory; formulas checked numerically.

## 3. Proof of Theorem 1

Coverage: a physical point `(J,u)` with `J < 400` is in (P). If `J ≥ 400`: `u ≥ 1/4` is (C); `u < 1/4`, `b̂ ≥ 120` is (M);
`u < 1/4`, `b̂ < 120` gives `u = (b̂/(2J+1))² ≤ (b̂/801)²`, which is (S). All regions are closed and overlap on boundaries.

**(P) [certified, `cert_P.py`].** For each `J < 400` (200 even `J` for `++`, 396 for `+−`) the polynomial of Lemma P
(cleared of denominators) is positive at `u = 0` and has no root in `[0,1]`: `[0,1]` is covered by dyadic intervals on
which the Möbius-transformed polynomial `(1+x)^d Q(1/(1+x))` has no sign variation (Descartes), endpoint values exactly
nonzero. Exact integer arithmetic (FLINT), 11 582 intervals in total, 59 s. This includes `m = M` and the near-zero at
`J = 245`.

**(C) [analytic + interval constants, `cert_C.py`].** `m¹⁰F = Φ + m⁶∫Ψ_hD_h`, `m⁶ ≤ 64`, `|Ψ_h| ≤ 2ψ₂ + |ψ₃|`. `++`: `|P_J| ≤ 1`,
`64(2∫ψ₂ + ∫|ψ₃|) ≤ 4402.1 < 10𝒥−5`. `+−`: `|d̃^J_{44}| ≤ min(binom(J+4,8), (1−p²)^{−4})` (Q3); splitting at
`(1−p*²)^{−4} = binom(J+4,8)` gives `∫ ≤ 62.90 + (377.88/4)·log binom(J+4,8) + 6044.5/(4(3/2)⁴) = 3888.6` at `J = 400`, and
`10𝒥 − 210 = 1.604·10⁶ > 64·3888.6`; the difference increases with `J`. Integrals/sups by Arb on 2000–4000 cells.

**(M) [analytic Lemma M + interval constants, `cert_M.py`].**
*Lemma M.* Let `L ≥ 0` be an integer, `N = L+½`, `m > 1`, `b = 2N/m`, and `w ∈ C⁴[0,1]` with `w(1)=w′(1)=w″(1)=0`. Then
`∫₀¹ w(p)P_L(1−2p²/m²)dp = w(0)/b − w″(0)/(2b³) + R`, `|R| ≤ C_R(w,u)/b⁴`,
`C_R ≤ (2/π)[(2/3)(|w‴(1)| + ‖w⁗‖₁) + 6u‖pw‴‖₁ + 4u‖w″‖₁] + u^{3/2}[(3π/4)sup|w″| + sup|w′|]`.
*Proof.* By (Q2) and Fubini, `∫wP_L = (m/π)∫₀^π sin(Nφ)G(φ)dφ`, `G(φ) = A(m sin(φ/2))`, `A(τ) = ∫₀^{min(1,τ)} w(p)(τ²−p²)^{−1/2}dp
= ∫₀^{t*}w(τ sin t)dt`, `t* = arcsin(min(1,1/τ))`. Because `w,w′,w″` vanish at 1, `A ∈ C³`, `A^{(k)} = ∫₀^{t*}sin^k t w^{(k)}(τ sin t)dt`
(`k ≤ 3`), `A⁗ = that + [τ>1] w‴(1)τ^{−3}t*′`, integrable. `G` is even about `φ = π` and `cos Nπ = 0`; four integrations by
parts give `∫₀^π sin(Nφ)G = G(0)/N − G″(0)/N³ + N^{−4}∫sin(Nφ)G⁗`, with `G(0) = πw(0)/2`, `G″(0) = πm²w″(0)/16`. Writing
`G⁗` by Faà di Bruno and changing to `τ` (`c = cos(φ/2) = (1−uτ²)^{1/2}`):
`(m/π)N^{−4}∫|G⁗| = (2/(πb⁴))∫₀^m |c³A⁗ − 6uτcA‴ + (3u²τ²/c − 4uc)A″ + (u²τ/c)A′|dτ`.
Swapping integrals (`p = τ sin t`): `∫|A⁗| ≤ (2/3)(|w‴(1)| + ‖w⁗‖₁)`, `∫τ|A‴| ≤ ‖pw‴‖₁`, `∫|A″| ≤ ‖w″‖₁`,
`τ²|A″| ≤ (π/4)sup|w″|`, `τ|A′| ≤ sup|w′|`, `∫₀^m dτ/c = π/(2√u)`. ∎ [proved; the expansion is checked numerically in
`sanity_numeric.py`: actual `b⁴R ∈ [−230, 310]`, bound 6852]
*Application.* `++`: `w₊ = Ψ₊`, `w₊(0) = 162u`, `w₊″(0) = 2(−80+21250u)`. `+−` (Lemma C, `m ≥ 2`): `𝔉₊₋ = Σ_j c_j H_{w₋}(J+j) + Φu³`,
`w₋ = Ψ₋(1−up²)^{−8}`, `w₋(0) = −162u`, `w₋″(0) = −2(80+21250u+1296u²)`, shifts `b̂_j = b̂ + 2j√u`. Arb (Taylor jets in `p`,
4000 cells, `u` as an interval) gives `C_R ≤ 6852` (`++`), `5867` (`+−`) for `u ≤ 1/100` and `≤ 5·10⁴`, `2.2·10⁵` for
`1/100 ≤ u ≤ 1/4`. With `10𝒥u³ = 2.5b̂²u² − 2.5u³`, `c_j ≥ 0`, `Σc_j = 1`, `b̂(1±δ)`, `δ = 8√u/b̂`, and minimising the
quadratic `2.5b̂²u² − 162u/b̂₋` over `u` (value `−2624.4/b̂₋⁴`), positivity reduces to `b̂⁴𝔉₊₋ ≥ 80b̂/(1+δ)³ − C_R/(1−δ)⁴
− 2624.4/((1−δ)²(1−η)) = 721 > 0` at `b̂ = 120` (increasing in `b̂`), `80b̂ − C_R = 2748 > 0` for `++`, and for `u ≥ 1/100`
the term `2.5b̂²u²` dominates. [certified, `cert_M.py`]

**(S) [certified, `cert_S.py`].** *Lemma S (bounds for complex arguments).* For `|u| ≤ ρ < 1` and complex `b̂`,
`P_ν(1−2up²) = (2/π)∫₀^{π/2} cos(b̂ p s(up²) sin σ)·E(θ_p,σ)^{−1/2} dσ`, `ν+½ = b̂/(2√u)`, `θ_p = 2arcsin(p√u)`,
`s(z) = arcsin(√z)/√z`, `E^{−1/2} = (a/sin a)^{1/2}(c/sin c)^{1/2}`, `a,c = θ_p(1±sin σ)/2` (from (Q1) with `φ = θ sin σ`, then
analytic continuation: both sides analytic in `ν` and in `θ`, `|θ| < π`). `s`, `arcsin` and `(z/sin z)^{1/2}` have nonnegative
Taylor coefficients, hence `|P_ν| ≤ cosh(Y)·θ̄/sin θ̄` with `Y = |Im b̂| p s(ρp²) + |b̂|p(s(ρp²)−1)`, `θ̄ = 2arcsin(p√ρ)`; for
`+−` (Lemma C) `|Im((J+j+½)θ_p)| ≤ Y + 4θ̄` and `Σ|c_j(J)| ≤ S_CG(Z_min)`, `|J+½| ≥ Z_min`. So `sup|G_h|` on a polydisc
`|b̂−b_c| ≤ r`, `|u| ≤ ρ` is bounded by an explicit `p`-integral `M`. [proved]
*Certificate.* `b̂ ∈ [1,120]`, cells of half-width `h = min(r/10, 0.1)`, `r = min(2, b_c/2)`, `u ∈ [0,u₀]`, `u₀ = ((b_c+h)/801)²`,
`ρ = 3u₀` (so `Z_min ≥ 110`). Exact rational coefficients `a_k(X)` of `G = Σ_k a_k u^k` (`k < 12`, `n ≤ 420`, rigorous geometric
tail), converted to Arb (768 bits); Taylor coefficients `g_{ik} = ∂^i_b̂ a_k(b_c)/i!` (`i < 9`). By Cauchy,
`𝔉 ≥ Σ_{i<9,k<12} g_{ik}t^iu^k + Φ-part − M[(h/r)⁹ + (u/ρ)¹²]/((1−h/r)(1−u/ρ))` on the cell; the polynomial is bounded below
on adaptively bisected boxes by a mean-value form with interval gradients. 604 cells per channel, at most 12 boxes per
cell, `M ≤ 418`. `b̂ ≤ 1` (`u ≤ 1/801²`): absolute-series bound `𝔉 ≥ 11.77`. Negative control: replacing `210` by
`10⁴` in `Φ₊₋` makes the certificate FAIL at `b̂ ≈ 22.6`. (`G` and its `u`-Taylor polynomial agree with the direct series
to 15 digits at test points.)

Theorem 1 follows. Since `𝔉(b̂,0) = 2f̂(b̂)`, (S) and (M) at `u = 0` re-prove `f̂ > 0` (`UNIFORM_POSITIVITY.md`, Theorem 2)
by an independent method.

## 4. The bound

### 4.1 Hypotheses (those of the paper, made explicit)
* **H1 (analyticity/crossing)** `f(s,u)` analytic as stated in §2.3 of the paper, crossing (2.5). *Physics hypothesis.*
* **H2 (unitarity)** (2.14): the heavy measure `⟨·⟩` is positive (`|c^{++}|², |c^{+−}|² ≥ 0`), spectrum `m ≥ M` for `J ≥ 3` and
  any state not in the EFT. *Physics hypothesis.*
* **H3 (Regge boundedness)** (2.23): `B₂^{(1)}, B₃^{(1)}` hold smeared against `ψ₂ 1_{[ε,1]}`, `ψ₃ 1_{[ε,1]}` (compact support, but not
  smooth at `ε`, as in the paper), and the forward limits `∂_{p²}B₄^{(1)}|₀`, `B₅^{(1)}|₀` hold unsmeared (the paper: "a technical
  modification, which we will ignore"); exchange of `∫dp` and `⟨·⟩`. *Physics hypothesis.*
* **H4 (tree-level EFT)** below `M` the amplitude is (2.7a) (+ (B.1) matter), loops neglected (`M ≪ M_Pl`). *Physics hypothesis.*
* **H5** spin selection: `c⁺⁻ = 0` for `J < 4`, `c⁺⁺ = 0` for odd `J` (footnote 8; angular momentum). *Kinematics.*
* **H-IR** see 4.3. *Physics hypothesis* (the paper's §3.1–3.2 prescription).

### 4.2 Derivation
Apply the g5-free functional with lower limit `ε` to (2.28). Low side (§1.3, exact): `16πG L(ε) − πG k(ε)|ĝ₃|² + C₀a₀ + C₂a₂`,
`C₀ = −(7773m_ℓ²+35)/630`, `C₂ = −(2223078m_ℓ⁶ − 2392390m_ℓ⁴ − 7555041m_ℓ² + 14856795)/(180180m_ℓ⁴)`, both `< 0` on `0 < m_ℓ ≤ 1`
(Sturm count, exact). High side: `⟨F^ε⟩`. This is Theorem 2 with `𝓜 = −Σ(C₀a₀+C₂a₂) ≥ 0`. If `⟨F^ε⟩ ≥ 0` then
`|ĝ₃|² ≤ 16L(ε)/k(ε)`; for `ε ≤ 1/100`: `ℓ` increasing, `ψ₃ > 0` on `[0,1/25]`, so `k(ε) ∈ [k_lo,k_hi]` with explicit rationals,
`16·65 log(1/ε)/k ≤ A log(1/ε)` and `16(ℓ(ε) − 1171/15)/k ≤ B` (negative numerator). [proved; exact in `cert_low.py`]

### 4.3 The infrared hypothesis
`F^ε[m,J,h] = F[m,J,h] − R_ε`, `m⁴|R_ε| ≤ ∫₀^ε|Ψ_h||D_h| ≤ (1−ε²)^{−4}(65ε² + 162uε + 17uε⁴)` (`ψ₂ ≤ 65p`, `0 < ψ₃ ≤ 162` on
`[0,0.04]`). By Theorem 1, `[F^ε]₋ ≠ 0` only for states with `0 < m⁴F < m⁴|R_ε|`, i.e. whose positive action is `O(ε²)`; in
region (M) (`++`, `u ≤ 1/100`) one has `m⁴F ≥ 22.9/b̂³ + 0.99·162u/b̂`, so these states have `b̂ ≳ 0.7 ε^{−2/3}` — the paper's
"negative only at `b ∼ m_IR^{−2/3}`". **H-IR:** the heavy-state contribution of this large-impact-parameter region is
non-negative, `⟨[F^ε]₋⟩ = 0` (equivalently `Δ_IR = 0`). This cannot be proved from H1–H5: in `D = 4` a functional positive at
all `b` must have `ψ₂(p) ∼ p` at `p → 0` and then the graviton pole makes the low side diverge (paper §3.1). It is the
statement that physics at `b ≳ m_IR^{−2/3}/M` is the known eikonal regime.

## 5. Status of every statement

| statement | status |
|---|---|
| sum rules (2.26), heavy sides (2.33)/(2.35a), sign of `B₃` per channel, forward-limit actions | [proved analytically] from (2.16), (2.26) |
| low-energy sides (2.31), (2.32), matter (3.8), functional low side, `k(0)`, `A`, `B`, matter signs | [certified, exact, `cert_low.py`] |
| bracketing of (3.3a) | [settled] by (3.4a) and (3.9) constants |
| `−5g₅` in the printed functional; g5-free completion | [proved, exact] |
| Lemma P (polynomial form) | [proved]; numerically cross-checked |
| Lemma C (CG reduction) | [proved] + [certified, `cert_cg.py`] |
| (Q1)–(Q4) classical formulas | quoted (numerically re-checked) |
| Lemma M (Mehler expansion with remainder) | [proved analytically] |
| Lemma S (complex Mehler bound) | [proved analytically] |
| region (P), `J < 400`, all `m ≥ M` | [certified, `cert_P.py`, exact integers] |
| region (C) | [proved] + [certified constants, `cert_C.py`] |
| region (M) | [proved] + [certified constants, `cert_M.py`] |
| region (S) | [certified, `cert_S.py`] |
| **Theorem 1** | [proved], given the above and correctness of FLINT/Arb |
| Theorem 2 | [proved] given H1–H5 |
| **Theorem 3** | [proved] given H1–H5 and H-IR |
| near-zeros (`b̂³m⁴F ≈ 0.31` at `m ≈ 21.8, J = 245`), scans, `b⁴R` sizes, localisation `b̂ ≳ ε^{−2/3}` outside (M) | [numerical only — not part of the theorem] |
| optimal bound (4.4) `24.9 log − 27.6` (SDPB, improved sum rules) | not addressed |

## 6. Reproduction
`bash ads4/code/flatbound/run_all.sh [--sanity]`. Files: `common.py` (exact data, `poly_u`, CG weights), `cert_P.py`,
`cert_C.py`, `cert_M.py`, `cert_S.py`, `cert_cg.py`, `cert_low.py`, `sanity_numeric.py` (floating point, not part of the proof),
`logs/`. Tools: python-flint 0.9 (FLINT/Arb ball arithmetic; no `mpmath.iv`), sympy. Parameters `JP = 400`, `B_A = 120`,
`u₁ = 1/100`.

## 7. What remains
* The hypotheses H1–H4 and H-IR are physics input; H-IR in particular is what turns positivity into a finite bound.
  Smearing of the forward-limit sum rules (H3) and the non-smooth IR cutoff are not justified here.
* Theorem 1 relies on FLINT/Arb being correct and on the quoted classical formulas (Q1)–(Q4).
* The optimal SDPB bound (4.4) and `F_{g4}` (3.3b) are not treated (the same code would apply to `F_{g4}` with new weights).
* Nothing here concerns AdS/CFT; the flat-space positivity is the `u`-family of `UNIFORM_POSITIVITY.md`'s Theorem 2, not
  its curved Theorem 1.
