# Heavy-state positivity of the cubic-coupling functional `F_{g3}` of arXiv:2201.06602, and what it implies for `|ĝ₃|²`

Status: written 16 Nov 2026 by an automated analyst; **revised twice the same day after the referee reports
in `ads4/REFEREE_FLAT_GRAVITON.md`** (round 1: IR hypothesis and published version, §1.5, §4; round 2: the eikonal
hypothesis is now imposed only on a certified large-`b̂`, small-`u` region, §4.3). Computer-assisted proof; code in `ads4/code/flatbound/`, driver
`bash run_all.sh` (≈ 4.5 min, ≈ 5.5 min with `--sanity`; one PASS/FAIL line per certificate and a final verdict; current
result: all PASS). Units `M = 1` unless `M` is written.

Versions of the paper: arXiv:2201.06602 **v1** (the only arXiv version; LaTeX source `gravity_4d.tex`) and the version of
record **JHEP05(2023)122** (full text via INSPIRE). They print *different* forward terms in (3.3a)/(3.3b); see §1.4–1.5.
Equation numbers are v1 unless marked "JHEP" (§2 numbers shift by −1 in JHEP).

## 0. Results

Notation (justified in §1–2). `ψ₂(p) = (1−p)³ p (65+155p+47p²)`, `ψ₃(p) = 2(1−p)³(81−1927p+4601p²)`. For a heavy state of
mass `m ≥ 1`, spin `J`, channel `h ∈ {++, +−}`: `u = 1/m²`, `𝒥 = J(J+1)`, `σ₊₊ = +1`, `σ₊₋ = −1`, `D₊₊ = P_J`,
`D₊₋ = d̃^J_{44} = P^{(0,8)}_{J−4}`, and

`F[m,J,h] = ∫₀¹ dp [ψ₂(p)(2m²−p²)/m⁶ + σ_h ψ₃(p)/m⁶] D_h(1−2p²/m²) + Φ_h(J)/m¹⁰`   (0.1)  (`cos θ = 1−2p²/m²`, `sin(θ/2) = p/m`)

`F` is the coefficient of `|c^h_{J,m²}|²` in the heavy action (common positive factor `16(2J+1)dm²/m²` dropped). Allowed
spins: `J ∈ {0,2,4,…}` for `++`, `J ≥ 4` for `+−`. Three versions of the forward term:

| version | forward term | `Φ₊₊` | `Φ₊₋` | low side contains |
|---|---|---|---|---|
| v1 as printed ("literal") | `−5∂_{p²}B₄^{(1)}|₀` | `10𝒥−5` | `10𝒥−205` | `−5g₅` |
| completed (this note) | `−5∂_{p²}B₄^{(1)}|₀ + 5B₅^{(1)}|₀` | `10𝒥` | `10𝒥−210` | no `g₅` |
| JHEP as printed | `−5∂_{p²}(B₄^{(1)}+B₄^{(2)})|₀` | `10𝒥−5` | `10𝒥−205`, plus an extra term ᵃ | `−5g₅` |

ᵃ `m¹⁰F^{JHEP}₊₋ = m¹⁰F^{lit}₊₋ + 5(−1)^J binom(J+4,8)(2𝒥−85)/(9m⁸)` (§1.5); not of the form `Φ/m¹⁰`.

**Theorem 1 (heavy-state positivity).** For every `m ≥ M`, every allowed `J` and both channels, `F[m,J,h] > 0` for the
v1-literal and for the completed functional. [certified; §3]

**Proposition 1′ (the JHEP-printed functional is not positive).** `m¹⁰F^{JHEP}(+−, J=4, m=M) = −6281/315 < 0`; more
precisely `m⁴F^{JHEP}₊₋(J=4) = 2591/210 − 41u/18 − 5u³ − 25u⁷`, negative for `M ≤ m < 1.090359 M`; also
`m¹⁰F^{JHEP}(+−, J=7, m=M) = −28385719/13860` and `(J=9) = −11495259817/171171`. Its low side still contains `−5g₅`
(`−B₄^{(2)}|_low ≡ 0`). [exact, `cert_jhep.py`]

**Theorem 2 (exact low-energy identity).** Under H1–H5 (§4.1), for every `0 < ε = m_IR/M < 1`, with the completed functional
`F^ε` (the `p`-integrals over `[ε,1]`):

`π G k(ε) |ĝ₃|² = 16 π G L(ε) − ⟨F^ε⟩_heavy − 𝓜`,  `L(ε) = ∫_ε¹ ψ₂/p² dp = 65 log(1/ε) − 1171/15 + ℓ(ε)`,

`k(ε) = 20 − 2∫_ε¹ψ₂p⁶ + 2∫_ε¹ψ₃p⁴`, `k(0) = 330151/12012 = 27.48509…`, `ℓ(ε) = 40ε + 223ε²/2 − 259ε³/3 + 7ε⁴/2 + 47ε⁵/5`,
`𝓜 ≥ 0` the light-matter term (§4.2). For the literal or JHEP functional the right side has in addition `−5g₅`.

**Theorem 3 (bound on `|ĝ₃|²`).** Let `0 < m_IR ≤ M/100`, `A = 37.8387`, `B = −45.2059` (rounded up; `A = 1040/k_lo`,
`B = 16(ℓ(1/100) − 1171/15)/k_hi`, exact rationals in `cert_low.py`; as `m_IR → 0` they tend to `37.8387` and `−45.4452`, the
constants of the paper's (3.4a)).
 (i) Under H1–H5: `|ĝ₃|²M⁸ ≤ 16L(ε)/k(ε) + Δ_IR ≤ A log(M/m_IR) + B + Δ_IR`, `Δ_IR := ⟨[F^ε]₋⟩/(πG k(ε)) ≥ 0`,
 where `[x]₋ = max(−x,0)`. `Δ_IR` is an explicit heavy-state quantity; by Theorem 1 it only receives contributions from
 states with `0 < m⁴F < m⁴|R_ε| ≤ (1−ε²)^{−4}(65ε² + 162uε + 17uε⁴)`, i.e. whose action is `O(ε² + uε)` — `O(ε)` for `u ∼ 1`, which
 at `ε = 1/100` includes ordinary `M`-scale states (referee, round 2). It is **not zero** if `G > 0` (§4.3).
 (ii) Let `0 < m_IR ≤ M/20000`. Then (certified, Lemma N) every state with `F^ε < 0` lies in the region
 `Ω_ε = {b̂ ≥ 120, u ≤ min(1/100, u_N(b̂;ε))}`, `u_N(b̂;ε) ≤ 65.4/b̂³ + 64.9ε/b̂² + 5.11ε/b̂` (so `s = m² ≥ 2.48·10⁴ M²` there).
 Under H1–H5 and the eikonal hypothesis `E(λ,η)` of §4.3 (on `Ω_ε` the `2/m⁴`-weighted heavy spectrum, as a function of
 `b̂`, is at most `(1+η)` times the tree-level graviton-exchange density `16πG b log⁺(1/(λ m_IR b))`):
 `Δ_IR ≤ (1+η)(16/k_lo)(65 I(λ) + 3.5ε/λ² + X(ε,λ))`, `I(λ) = ∫₀^{1/λ} J₁(x) log(1/(λx)) dx`, `X` explicit (finite-mass
 corrections, `X ≤ 0.118` at `λ = 1`). With `A = 37.8387`, `B = −45.4441` for `m_IR ≤ M/20000` and `λ = 1, η = 0`:
 **`|ĝ₃|²M⁸ ≤ 37.8387 log(M/m_IR) − 40.790`** (certified `−40.79081`). The coefficient of the logarithm is unchanged; the
 correction is O(1) in `B` (table in §4.3). Light matter can only lower the right side. For `M/20000 < m_IR ≤ M/100` only (i)
 is proved; the analogous numbers there (e.g. `−40.60` at `λ = 1`) are **estimates**, because the negative set then contains
 non-eikonal `M`-scale states.
For the literal or JHEP functional (i) and (ii) hold with an extra `−5g₅M¹⁰/(πG k)` (`g₅` has either sign); the same
omission (`−g₅`) is present in the low side of `F_{g4}` (3.3b)/(3.4b).

How positivity was established (Theorem 1), `b̂ := (2J+1)/m`:

| region | method | status |
|---|---|---|
| (P) `J < 400`, all `m ≥ 1` | `m⁴F` is an exact rational polynomial in `u` of degree `≤ J+1`; positivity on the closed interval `u∈[0,1]` by Descartes' rule of signs with dyadic bisection (exact integers) | certified, `cert_P.py` (596 polynomials) |
| (C) `J ≥ 400`, `1 ≤ m ≤ 2` | `Φ ≥ 10𝒥−210` beats crude bounds `|P_J| ≤ 1`, `|d̃^J_{44}| ≤ min(binom(J+4,8),(1−p²)^{−4})` | analytic + interval constants, `cert_C.py` |
| (M) `m ≥ 2`, `b̂ ≥ 120` (all `J`) | exterior Mehler formula, four integrations by parts: exact `1/b̂`, `1/b̂³` terms + explicit `O(b̂⁻⁴)` remainder; `+−` reduced to Legendre by a Clebsch–Gordan identity | analytic + interval constants, `cert_M.py` |
| (S) `J ≥ 400`, `m ≥ 2`, `b̂ ≤ 120` | continuous extension in `(b̂,u)`; exact double Taylor coefficients + Cauchy bounds on a polydisc from a Mehler–Dirichlet representation valid for complex arguments; 604 `b̂`-cells per channel, adaptive boxes | certified, `cert_S.py` |

**Findings.**
1. Positivity holds but is **nearly marginal** in the `+−` channel: the exact value at `m = M`, `J = 4` is `m¹⁰F = 19/315 ≈ 0.060`
   (completed; `1594/315 ≈ 5.06` literal), against `≈ 12` for generic states; and at `m ≈ 21.79 M`, `J = 245` (`b̂ ≈ 22.5`)
   `m⁴F ≈ 2.67·10⁻⁵`, `b̂³m⁴F ≈ 0.31` against `≈ 80` for the flat asymptotics (floating point, not used). This point lies
   outside the plotted range `m ≤ 16` of the paper's Fig. 3; in the plotted normalisation `m⁸F` it would be `≈ 6` (the figure's
   colour bar reads `m⁸F`, its caption says `m¹⁰`).
2. **The forward term printed in v1 (3.3a) is probably a misprint.** Its low side contains `−5g₅`, absent from (3.4a). The
   paper's own Fig. 3 of v1 (arXiv source `figures/func_g3_pm.pdf`, `func_g3_pp.pdf`, rendered and inspected here, as the referee
   did) matches the **completed** functional: in the `+−` panel the point `m = 1, b = 2J/m = 8` (`J = 4`) lies in the lowest
   band `< 0.1`, matching completed `19/315 ≈ 0.060`, not literal `5.06` (band 1–10); in the `++` panel `m = 1, b = 0`
   (`J = 0`) lies in the band 10–100, matching completed `3904/315 ≈ 12.39`, not literal `2329/315 ≈ 7.39` (band `< 10`). The
   `++` island `< 10` at `m ≈ 5–7, b ≈ 8–9` also matches (here `m⁸F ≈ 2.7` at `m = 6.15, J = 26`). So (3.4a) and the plots
   were almost certainly computed with a `g₅`-free combination; the defect is a typographical omission in v1.
3. **JHEP did not repair it** (Proposition 1′): the term `−5∂_{p²}B₄^{(2)}|₀` added in the published (3.3a) has vanishing low
   side, leaves `−5g₅` in place, and makes the functional negative on `+−` states with `m ≈ M` (`J = 4, 7, 9, …`), contradicting
   the claimed and plotted positivity.
4. The IR cutoff destroys positivity at large impact parameter (as the paper says, (3.6)); a finite bound needs a hypothesis
   about that region. The earlier version of this note assumed `⟨[F^ε]₋⟩ = 0`, which is incompatible with `G > 0` (§4.3) and
   has been withdrawn. Theorem 3 now states the unconditional form (i) and, for `m_IR ≤ M/20000`, a bound (ii) under an
   eikonal hypothesis imposed only on a certified large-`b̂`, large-`s` region that contains every negative state.

## 1. The sum rules (arXiv:2201.06602)

### 1.1 Amplitudes, partial waves, unitarity
`M(1⁺2⁻3⁻4⁺) = ⟨23⟩⁴[14]⁴ f(s,u)` (2.1), `s+t+u = 0`, `u = −(p₁+p₃)²` (2.4), `f(s,u) = f(u,s)` (2.5). Low energy (2.7a):
`f_low = 8πG/(stu) + 2πG su|ĝ₃|²/t + g₄ + g₅t + g₆t² − g₆′su + … + f_matter`. Partial waves (2.12), unitarity (2.14)
`Im a_J^{+−+−} = |c^{+−}_{J,s}|²`, `Im a_J^{++−−} = |c^{++}_{J,s}|²`, and (2.16): at `s = m²`
`Im f(m²,−p²) = (16π/m⁸) Σ_{J≥4}(2J+1)|c^{+−}_{J,m²}|² d̃^J_{4,4}(1−2p²/m²)` (s-channel, all integers `J ≥ 4`), at `t = m²`
`Im f(−p²,p²−m²) = (16π/m⁸) Σ_{J≥0 even}(2J+1)|c^{++}_{J,m²}|² d̃^J_{0,0}(1−2p²/m²)`. (A.1)–(A.3): `d̃^J_{00} = P_J`,
`d̃^J_{44}(x) = ₂F₁(4−J, J+5; 1; (1−x)/2) = P^{(0,8)}_{J−4}(x)` (`N^J_{44} = 1`), `d^J_{44} = ((1+x)/2)⁴ d̃^J_{44}`.

### 1.2 Fixed-`u` sum rules and their heavy side
(2.26a) `B_k^{(1)}(p²) = ∮ ds/(4πi) (2s−p²)/[s(s−p²)]^{(k−2)/2} [f(s,−p²) + f(p²−s,−p²)] = 0`, `k ≥ 2` even;
(2.26b) `B_k^{(1)} = ∮ ds/(4πi) [s(s−p²)]^{−(k−3)/2} [−f(s,−p²) + f(p²−s,−p²)] = 0`, `k ≥ 3` odd; valid after smearing
(2.21)–(2.25). Contour deformation (2.28)–(2.29): `−B|_low = B|_high`. Inserting (2.16): `f(s,−p²)` has the `s`-channel cut
(`+−` states), `f(p²−s,−p²) = f(−p², p²−s)` the `t`-channel cut at `t = s` (`++` states). With
`⟨X⟩ = 16Σ_J(2J+1)∫_{M²}^∞ (dm²/m²) X` (2.34), `x = 1−2p²/m²`:

* (2.33) `B_k|_high = ⟨(2m²−p²)(|c⁺⁺|²d̃₀₀(x) + |c⁺⁻|²d̃₄₄(x)) / (m^{k+4}(m²−p²)^{(k−2)/2})⟩`, `k` even;
* (2.35a) `B_k|_high = ⟨(|c⁺⁺|²d̃₀₀(x) − |c⁺⁻|²d̃₄₄(x)) / (m^{k+3}(m²−p²)^{(k−3)/2})⟩`, `k` odd.

Both re-derived from (2.26) and (2.16) (the `B₃` term has **opposite signs** in the two channels because of the `−f(s,…)`
in (2.26b)). Per state (dropping `16(2J+1)dm²/m²`): `C₂ = (2m²−p²)D_h/m⁶`, `C₃ = σ_h D_h/m⁶`,
`C₄ = (2m²−p²)D_h/(m⁸(m²−p²))`, `C₅ = σ_h D_h/(m⁸(m²−p²))`. With `D₊₊′(1) = 𝒥/2`, `D₊₋′(1) = (J−4)(J+5)/2`:
`∂_{p²}C₄|₀ = (1−2𝒥)/m¹⁰` (`++`), `(41−2𝒥)/m¹⁰` (`+−`); `C₅|₀ = σ_h/m¹⁰`. This gives (0.1).

### 1.3 Low-energy side (re-derived, `cert_low.py`)
`−B|_low = ½ Σ Res` at the poles inside the arcs. From (2.7a):
`−B₂ = 16πG/p² + 2πG|ĝ₃|²p⁶`, `−B₃ = −2πG|ĝ₃|²p⁴` ((2.31a,b)), `−B₄ = 2g₄ + (4πG|ĝ₃|² + g₅)p² + (g₆+g₆′)p⁴` ((2.32a)),
`−B₅ = g₅ + (g₆−g₆′)p²` ((2.32b)) — all checked symbolically. Matter (B.1), `a_J = |g_J(m_ℓ)|²`:
`−B₂ = (p²−2m_ℓ²)(a₀m_ℓ⁴ + a₂(m_ℓ⁴ − 6m_ℓ²p² + 6p⁴))/m_ℓ⁴`, `−B₃ = −(a₀m_ℓ⁴ + a₂(m_ℓ⁴ − 6m_ℓ²p² + 6p⁴))/m_ℓ⁴`,
`−B₄ = 12a₂p²/m_ℓ⁴`, `−B₅ = 0`. These agree **exactly** with (3.8) once the paper's `B|_matter` is read as `B|_low`, except
that (3.8a,b) print the factor as `(m_ℓ² − 6m_ℓp² + 6p⁴)`, which is dimensionally inconsistent: a misprint for
`m_ℓ⁴ − 6m_ℓ²p² + 6p⁴` (with the corrected factor, (3.9) is reproduced exactly).

### 1.4 The functional in v1; bracketing; the `g₅` term
v1 (3.3a): `F_{g3} = ∫_{m_IR}^M dp (1−p)³[p(65+155p+47p²)B₂ + 2(81−1927p+4601p²)B₃] − 5∂_{p²}B₄^{(1)}|_{p=0}`. The bracket is
explicit in the LaTeX source (`(1-p)^3 \big(p(65+…)B_2^{(1)}+2(81-…)B_3^{(1)}\big)`, referee); it is also forced by the paper's
numbers: with it `A = 1040/k(0) = 37.8387`, `B = −45.4452` ((3.4a): 37.8, −45.4) and (3.9) is reproduced exactly; the other
reading gives `A = 0.728` and `11711/9` instead of `1/18`. (A commented-out older version in the source uses
`−(5/2)∂_{p²}B₄` with a different normalisation; irrelevant.)
The low-energy action of the printed functional is (exact) `16πG L(ε) − πG k(ε)|ĝ₃|² − 5g₅ + (matter)`. Since
`g₅ = ⟨(|c⁺⁺|²−|c⁺⁻|²)/m¹⁰⟩` ((2.41)) has no sign, (3.4a) does not follow from (3.3a) as printed. Adding `5B₅^{(1)}|₀` (low side
`5g₅`, heavy side `5σ_h/m¹⁰`) removes it and changes nothing else; by Finding 2 this is almost certainly what the authors
computed. The same holds for `F_{g4}` (3.3b): its low side is `−14g₄ − g₅ − (70498/15015)πG|ĝ₃|² + 16πG∫ψ₂^{(g4)}/p²`, so (3.4b) needs
`+B₅^{(1)}|₀` (or `g₅ ≥ 0`); the check reproduces `99/7 = 14.14` and `35249/840840 = 0.042` of (3.4b) (`cert_low.py`).
All certificates use `Φ₊₊ = 10𝒥−5`, `Φ₊₋ = 10𝒥−210`, the smaller of the literal and completed values in each channel.

### 1.5 The published version JHEP05(2023)122
JHEP prints (3.3a) `F_{g3} = ∫ … − 5∂_{p²}(B₄^{(1)} + B₄^{(2)})|_{p=0}` and (3.3b) `F_{g4} = ∫ … − (7B₄^{(1)} + ∂_{p²}B₄^{(1)} + ∂_{p²}B₄^{(2)})|_{p=0}`,
while the text still says "These use only `B^{(1)}`". Here `B_k^{(2)}(p²) = ∮ ds/(4πi) (2s−p²) f(s,p²−s)/[s(s−p²)]^{(k+2)/2}`
(JHEP (2.26), fixed `t = −p²`). Exactly (`cert_jhep.py`):
* (a) **Low side**: `−B₄^{(2)}|_low ≡ 0` (EFT (2.7a) through `g₆`, including graviton and `ĝ₃` terms, and matter (B.1)). So the
  added term is a null constraint and the `−5g₅` (resp. `−g₅`) survives.
* (b) **Heavy action**: at fixed `t` only the `s`- and `u`-channel cuts appear, i.e. only `+−` states. From (2.16a), the
  residue of `f(s,u)` at `s = m²` is `∝ |c⁺⁻|² d̃^J_{44}(1+2u/m²)`; at fixed `t`, `u = p²−m²`, so the per-state density of `B₄^{(2)}`
  is `(2m²−p²) P^{(0,8)}_{J−4}(−1+2p²/m²)/(m¹²(m²−p²)³)` — the same as JHEP (2.34b) with `(−1)^J d̃^J_{4,−4}(x) = d̃^J_{44}(−x)`.
  At `p = 0` this is the backward value `(−1)^J binom(J+4,8)`. Hence
  `m¹⁰F^{JHEP}₊₋ = m¹⁰F^{lit}₊₋ + 5(−1)^J binom(J+4,8)(2𝒥−85)/(9m⁸)` (checked exactly for `J = 4…40`), `F^{JHEP}₊₊ = F^{lit}₊₊`.
* (c) **Not positive**: Proposition 1′. The added term is `−25/m⁸` at `J = 4` and `−2475/m⁸` at `J = 7` and grows like
  `J¹⁰/m⁸` with alternating sign, so the published functional is negative on infinitely many `+−` states near `m = M`.
  Like the plots (Finding 2), this indicates that the printed forward terms are not the ones used for (3.4a).

## 2. The heavy action in normal form

`𝔉_h(J,u) := m⁴F[m,J,h] = ∫₀¹ Ψ_h(p;u) D_h(1−2up²) dp + Φ_h u³`, `Ψ_h = ψ₂(2−up²) + σ_h uψ₃`.

**Lemma P (exact polynomial).** With `μ₂_j = ∫ψ₂p^j`, `μ₃_j = ∫ψ₃p^j` (rationals), `e₊₊(i) = i(i−1)`, `e₊₋(i) = (i+3)(i+4)`:
`𝔉_h = Σ_{k≥0} (−u)^k Π_{i≤k}(𝒥 − e_h(i))/(k!)² · [2μ₂_{2k} + u(−μ₂_{2k+2} + σ_hμ₃_{2k})] + Φ_h u³` (finite for integer `J`).
*Proof:* `P_J(1−2y) = ₂F₁(−J,J+1;1;y)`, `P^{(0,8)}_{J−4}(1−2y) = ₂F₁(4−J,J+5;1;y)`. ∎ [proved; checked against quadrature to
10⁻¹¹ and by the referee to 12 digits]. Example: `𝔉₊₋(4,u) = 2591/210 − 41u/18 − 10u³` (completed), minimum `19/315` at `u = 1`.

**Continuous extension.** With `X = b̂²/4`, `β := 𝒥u = X − u/4`, `λ_i = ¼ + e_h(i)`:
`𝔉_h = G_h(b̂,u) + 10Xu² − (5/2 + c_h)u³`, `G_h = Σ_n (−1)^n/(n!)² Π_{i≤n}(X − λ_i u)·(T0_n + T1_n u)`, `c₊₊ = 5`, `c₊₋ = 210`,
`T0_n = 2μ₂_{2n}`, `T1_n = −μ₂_{2n+2} + σ_hμ₃_{2n}`. `G_h` is entire in `b̂` and analytic in `|u| < 1`; at physical points it equals
the integral. Flat slice: `G_h(b̂,0) = 2∫ψ₂J₀(b̂p)dp = 2f̂(b̂)` in both channels (scaling limit (3.5)).

**Lemma C (Clebsch–Gordan reduction of the `+−` kernel).** For all complex `J` with `c_j` finite and `|y| < 1`:
`(1−y)⁸ ₂F₁(4−J,J+5;1;y) = Σ_{j=−4}^{4} c_j(J) ₂F₁(−J−j, J+j+1; 1; y)`, `c_j = ⟨J4;4−4|J+j,0⟩²` the explicit rational functions in
`common.py` (e.g. `c₀ = 35(J−3)(J−2)(J+3)(J+4)/(8(2J−3)(2J−1)(2J+3)(2J+5))`); `c_j ≥ 0` for `J ≥ 4`, `Σc_j = 1`.
*Proof:* for integer `J ≥ 4`, the CG product rule `d^J_{44}d⁴_{−4,−4} = Σ_L ⟨J4;4−4|L0⟩² d^L_{00}` with `d⁴_{−4,−4} = cos⁸(θ/2) = (1−y)⁴`
gives the identity with weights `CG²`. The coefficients of `y⁰…y⁸` determine the weights uniquely (the 9×9 determinant,
computed symbolically, has no zero at integers `J ≥ 4`), and `cert_cg.py` checks symbolically that the `c_j(J)` solve these
equations; so `CG² = c_j(J)` for integer `J ≥ 4`, and each coefficient of `y^k` of LHS−RHS is a rational function of `J`
vanishing at all integers `J ≥ 4`, hence `≡ 0`. ∎ [proved; also `y⁹…y¹⁴` checked identically in `J`]

**Classical inputs.** (Q1) DLMF 14.12.1 with `μ = 0` (checked on the DLMF site): `P_ν(cos θ) = (√2/π)∫₀^θ cos((ν+½)t)(cos t − cos θ)^{−1/2}dt`,
`0 < θ < π`, all complex `ν` (Ferrers function). (Q2) For integer `n`: `P_n(cos θ) = (√2/π)∫_θ^π sin((n+½)φ)(cos θ − cos φ)^{−1/2}dφ`.
*Proof:* apply DLMF 18.10.2 (the case `ν = n` of (Q1)) at `π−θ`, use `P_n(−x) = (−1)^nP_n(x)` and substitute `φ ↦ π−φ`
(`cos((n+½)(π−φ)) = (−1)^n sin((n+½)φ)`). Checked numerically to 10⁻¹⁷ (referee 10⁻²²). (Q3) `|d^J_{mm′}| ≤ 1` (unitarity);
`max_{[−1,1]}|P_n^{(0,8)}| = |P_n^{(0,8)}(−1)| = binom(n+8,n)` (DLMF 18.14.2, Szegő Thm 7.32.1; checked on the DLMF site). (Q4) the
CG product (Clebsch–Gordan series) rule for Wigner `d`-functions (standard; numerically checked in `sanity_numeric.py` §2).

## 3. Proof of Theorem 1

Coverage: a physical point `(J,u)` with `J < 400` is in (P). If `J ≥ 400`: `u ≥ 1/4` is (C); `u < 1/4`, `b̂ ≥ 120` is (M);
`u < 1/4`, `b̂ < 120` gives `u = (b̂/(2J+1))² ≤ (b̂/801)²`, which is (S) (with `b̂ ≤ 1` handled by `small_b`).

**(P) [certified, `cert_P.py`].** For each `J < 400` (200 even `J` for `++`, 396 for `+−`) the polynomial of Lemma P
(cleared of denominators) is positive at `u = 0` and has no root in `[0,1]`: `[0,1]` is covered by dyadic intervals on
which `(1+x)^d Q(1/(1+x))` has no sign variation (Descartes), endpoint values exactly nonzero. Exact integers (FLINT), 11 582
intervals, ≈ 60 s. This includes `m = M` and the near-zero at `J = 245`.

**(C) [analytic + interval constants, `cert_C.py`].** `m¹⁰F = Φ + m⁶∫Ψ_hD_h`, `m⁶ ≤ 64`, `|Ψ_h| ≤ 2ψ₂ + |ψ₃|`. `++`: `|P_J| ≤ 1`,
`64(2∫ψ₂ + ∫|ψ₃|) ≤ 4402.1 < 10𝒥−5`. `+−`: `|d̃^J_{44}| ≤ min(binom(J+4,8), (1−p²)^{−4})` (Q3); splitting at
`(1−p*²)^{−4} = binom(J+4,8)` gives `∫ ≤ 62.90 + (377.88/4)·log binom(J+4,8) + 6044.5/(4(3/2)⁴) = 3888.6` at `J = 400`, and
`10𝒥 − 210 = 1.604·10⁶ > 64·3888.6`; the difference increases with `J`.

**(M) [analytic Lemma M + interval constants, `cert_M.py`].**
*Lemma M.* Let `L ≥ 0` be an integer, `N = L+½`, `m > 1`, `b = 2N/m`, and `w ∈ C⁴[0,1]` with `w(1)=w′(1)=w″(1)=0`. Then
`∫₀¹ w(p)P_L(1−2p²/m²)dp = w(0)/b − w″(0)/(2b³) + R`, `|R| ≤ C_R(w,u)/b⁴`,
`C_R ≤ (2/π)[(2/3)(|w‴(1)| + ‖w⁗‖₁) + 6u‖pw‴‖₁ + 4u‖w″‖₁] + u^{3/2}[(3π/4)sup|w″| + sup|w′|]`.
*Proof.* By (Q2) and Fubini, `∫wP_L = (m/π)∫₀^π sin(Nφ)G(φ)dφ`, `G(φ) = A(m sin(φ/2))`, `A(τ) = ∫₀^{min(1,τ)} w(p)(τ²−p²)^{−1/2}dp
= ∫₀^{t*}w(τ sin t)dt`, `t* = arcsin(min(1,1/τ))`. Because `w,w′,w″` vanish at 1, `A ∈ C³`, `A^{(k)} = ∫₀^{t*}sin^k t w^{(k)}(τ sin t)dt`
(`k ≤ 3`), `A⁗ = that + [τ>1] w‴(1)τ^{−3}t*′`, integrable. `G` is even about `φ = π` and `cos Nπ = 0`; four integrations by parts
give `∫₀^π sin(Nφ)G = G(0)/N − G″(0)/N³ + N^{−4}∫sin(Nφ)G⁗`, with `G(0) = πw(0)/2`, `G″(0) = πm²w″(0)/16`. By Faà di Bruno and
`τ = m sin(φ/2)` (`c = cos(φ/2) = (1−uτ²)^{1/2}`): `(m/π)N^{−4}∫|G⁗| = (2/(πb⁴))∫₀^m |c³A⁗ − 6uτcA‴ + (3u²τ²/c − 4uc)A″ + (u²τ/c)A′|dτ`.
Swapping integrals (`p = τ sin t`): `∫|A⁗| ≤ (2/3)(|w‴(1)| + ‖w⁗‖₁)` (using `∫₁^∞ dτ/(τ⁴√(τ²−1)) = 2/3`, `∫₀^{π/2}sin³ = 2/3`),
`∫τ|A‴| ≤ ‖pw‴‖₁`, `∫|A″| ≤ ‖w″‖₁`, `∫₀^m dτ/c = π/(2√u)`, `τ|A′| ≤ sup|w′|·(τ − √(τ²−min(1,τ)²)) ≤ sup|w′|`, and
`τ²|A″| ≤ sup|w″|·∫₀^{min(1,τ)} s²(τ²−s²)^{−1/2}ds ≤ (π/4)sup|w″|` (for `τ ≤ 1` the integral is `πτ²/4`; for `τ ≥ 1` it decreases
in `τ` and equals `∫₀¹ s²/√(1−s²) = π/4` at `τ = 1`). ∎ [proved; referee re-derived; numerically `b⁴R ∈ [−230, 310]`, bound 6852]
*Application.* `++`: `w₊ = Ψ₊`, `w₊(0) = 162u`, `w₊″(0) = 2(−80+21250u)`. `+−` (Lemma C, `m ≥ 2`): `𝔉₊₋ = Σ_j c_j H_{w₋}(J+j) + Φu³`,
`w₋ = Ψ₋(1−up²)^{−8}`, `w₋(0) = −162u`, `w₋″(0) = −2(80+21250u+1296u²)`, shifts `b̂_j = b̂ + 2j√u ∈ [b̂(1−δ), b̂(1+δ)]`, `δ = 8√u/b̂`.
Arb (Taylor jets in `p`, 4000 cells, `u` an interval) gives `C_R ≤ 6852` (`++`), `5867` (`+−`) for `u ≤ 1/100`, and `≤ 5·10⁴`,
`2.2·10⁵` for `1/100 ≤ u ≤ 1/4`. With `10𝒥u³ = 2.5b̂²u² − 2.5u³`, `c_j ≥ 0`, `Σc_j = 1`, and
`min_u [2.5(1−η)b̂²u² − 162u/b̂₋] = −2624.4/((1−η)b̂²b̂₋²) = −2624.4/((1−η)(1−δ)²b̂⁴)`, positivity reduces to
`b̂⁴𝔉₊₋ ≥ 80b̂/(1+δ)³ − C_R/(1−δ)⁴ − 2624.4/((1−δ)²(1−η)) = 721 > 0` at `b̂ = 120` (increasing in `b̂`), `80b̂ − C_R = 2748 > 0` for
`++`, and for `u ≥ 1/100` the term `2.5b̂²u²` dominates. [certified, `cert_M.py`]

**(S) [certified, `cert_S.py`].** *Lemma S.* For `|u| ≤ ρ < 1` and complex `b̂`,
`P_ν(1−2up²) = (2/π)∫₀^{π/2} cos(b̂ p s(up²) sin σ)·E(θ_p,σ)^{−1/2} dσ`, `ν+½ = b̂/(2√u)`, `θ_p = 2arcsin(p√u)`, `s(z) = arcsin(√z)/√z`,
`E^{−1/2} = (a/sin a)^{1/2}(c/sin c)^{1/2}`, `a,c = θ_p(1±sin σ)/2` (from (Q1) with `t = θ sin σ`, then analytic continuation in
`ν` and `θ`, `|θ| < π`). `s`, `arcsin`, `(z/sin z)^{1/2}` have nonnegative Taylor coefficients, so `|P_ν| ≤ cosh(Y)·θ̄/sin θ̄`,
`Y = |Im b̂| p s(ρp²) + |b̂|p(s(ρp²)−1)`, `θ̄ = 2arcsin(p√ρ)`; for `+−` (Lemma C) `|Im((J+j+½)θ_p)| ≤ Y + 4θ̄` and
`Σ|c_j(J)| ≤ S_CG(Z_min)` for `|J+½| ≥ Z_min`. So `sup|G_h|` on a polydisc `|b̂−b_c| ≤ r`, `|u| ≤ ρ` is bounded by an explicit
`p`-integral `M`. [proved]
*Certificate.* `b̂ ∈ [1,120]`, cell half-width `h = min(r/10, 0.1)`, `r = min(2, b_c/2)`, `u ∈ [0,u₀]`, `u₀ = ((b_c+h)/801)²`,
`ρ = 3u₀` (`Z_min ≥ 110`). Exact rational `a_k(X)` of `G = Σ_k a_k u^k` (`k < 12`, `n ≤ 420`, ratio-checked geometric tail
using `|T0_n|, |T1_n| ≤ 2∫ψ₂, ∫ψ₂+∫|ψ₃| ≤ 100`, now asserted in code), in Arb at 768 bits; `g_{ik} = ∂^i_b̂ a_k(b_c)/i!` (`i < 9`). By
Cauchy, `𝔉 ≥ Σ_{i<9,k<12} g_{ik}t^iu^k + Φ-part − M[(h/r)⁹ + (u/ρ)¹²]/((1−h/r)(1−u/ρ))`; the polynomial is bounded below on adaptively
bisected boxes by a mean-value form with interval gradients. 604 cells per channel, ≤ 12 boxes per cell, `M ≤ 418`.
`b̂ ≤ 1`: absolute-series bound `𝔉 ≥ 11.77` (tail ratio `< 1/2` asserted). Mutation tests (referee and here) fail as they should;
the certificate is sharp to ≈ 5% (referee).

Theorem 1 follows. Since `𝔉(b̂,0) = 2f̂(b̂)`, (S) at `u = 0` re-proves `f̂ > 0` on `[0,120]` (`UNIFORM_POSITIVITY.md`, Theorem 2)
by an independent method.

## 4. The bound

### 4.1 Hypotheses (those of the paper, made explicit)
* **H1 (analyticity/crossing)** as in §2.3 of the paper, crossing (2.5). *Physics hypothesis.*
* **H2 (unitarity)** (2.14): `⟨·⟩` is a positive measure; states not in the EFT have `m ≥ M`. *Physics hypothesis.*
* **H3 (Regge boundedness)** (2.23): `B₂^{(1)}, B₃^{(1)}` hold smeared against `ψ₂1_{[ε,1]}`, `ψ₃1_{[ε,1]}` (compact support, not
  smooth at `ε`, as in the paper); the forward limits `∂_{p²}B₄^{(1)}|₀`, `B₅^{(1)}|₀` hold unsmeared (the paper: "a technical
  modification, which we will ignore"); `∫dp` and `⟨·⟩` commute. *Physics hypothesis.*
* **H4 (tree-level EFT)** below `M` the amplitude is (2.7a) (+ (B.1) matter), loops neglected (`M ≪ M_Pl`). *Physics hypothesis.*
* **H5** spin selection: `c⁺⁻ = 0` for `J < 4`, `c⁺⁺ = 0` for odd `J` (footnote 8). *Kinematics.*

### 4.2 Derivation of Theorems 2 and 3(i)
Apply the completed functional with lower limit `ε` to (2.28). Low side (§1.3, exact): `16πG L(ε) − πG k(ε)|ĝ₃|² + C₀a₀ + C₂a₂`,
`C₀ = −(7773m_ℓ²+35)/630`, `C₂ = −(2223078m_ℓ⁶ − 2392390m_ℓ⁴ − 7555041m_ℓ² + 14856795)/(180180m_ℓ⁴)`, both `< 0` on `0 < m_ℓ ≤ 1`
(Sturm count, exact; these are `−` the polynomials of (3.9)). High side: `⟨F^ε⟩ ≥ −⟨[F^ε]₋⟩`. This is Theorem 2, and
Theorem 3(i) follows with `𝓜 ≥ 0` dropped. For `ε ≤ 1/100`: `ℓ` increasing, `ψ₃ > 0` on `[0,1/25]`, so
`k(ε) ∈ [k_lo,k_hi]` (explicit rationals), `16·65 log(1/ε)/k ≤ A log(1/ε)` and `16(ℓ(ε) − 1171/15)/k ≤ B` (negative numerator).
[proved; exact in `cert_low.py`]

### 4.3 The infrared term `Δ_IR`
**Why `Δ_IR ≠ 0` (the earlier hypothesis H-IR is withdrawn).** `m⁴F^ε = 𝔉 − R_ε` with `|R_ε| ≤ (1−ε²)^{−4}(65ε² + 162uε + 17uε⁴)`.
Theorem 2 holds for every `ε`, so the heavy measure must reproduce the `16πG/p²` pole of `−B₂|_low` for all `p ∈ (0,1]`. In
the scaling regime (`m → ∞`, `b = 2J/m` fixed) put `dν(b) := ⟨(2/m⁴)(|c⁺⁺|²+|c⁺⁻|²)δ(b − 2J/m)⟩db`; then `B₂|_high(p) → ∫J₀(bp)dν(b)`
and `⟨F^ε⟩ → ∫f̂_ε(b)dν(b)`, `f̂_ε(b) = ∫_ε¹ψ₂J₀(bp)dp`. Reproducing `16πG/p²` requires the eikonal density
`dν ≈ 16πG b log(1/(μb)) db` at large `b` (`∫₀^{1/μ} b log(1/(μb))J₀(bp)db = (1−J₀(p/μ))/p²`), which populates
`ε^{−2/3} ≲ b ≲ ε^{−1}`, where `f̂_ε ≈ 40/b³ − 65εJ₁(εb)/b < 0` (paper (3.6)). Hence `⟨[F^ε]₋⟩ = O(G) > 0` when `G > 0`: assuming
`⟨[F^ε]₋⟩ = 0` contradicts H1–H4 (the referee's argument, accepted). The paper does not assume this either: it says that,
if the large-`b` region is "controlled by known low-energy-computable eikonal physics, positivity should not be crucial
there" (§2.5) and "accept[s] that this causes negativity at large impact parameters" (§3.1).

**Where is `F^ε` negative? (Lemma N, certified, `cert_eps.py`).** Let `ε* = 1/20000` and, for `0 < ε ≤ ε*`,
`r₀(ε) = 65ε²(1−ε²)^{−4}`, `r₁(ε) = (162ε + 17ε⁴)(1−ε²)^{−4}`; then `m⁴|R_ε| ≤ r₀(ε) + r₁(ε)u` (`ψ₂ ≤ 65p`, `0 < ψ₃ ≤ 162` on
`[0, 0.04]`, `|P_J| ≤ 1`, `|d̃₄₄| ≤ (1−y)^{−4}`). Since `[F^ε]₋ ≤ [R_ε]₊` (Theorem 1), a state can be in `N_ε := {F^ε < 0}` only if
`m⁴F < r₀(ε) + r₁(ε)u`. Let `u_N(b̂;ε)` be the largest root of `2.5(1−η)b̂²u² − (162/b̂₋ + r₁(ε))u − r₀(ε)` and
`Ω_ε := {b̂ ≥ 120, u ≤ min(1/100, u_N(b̂;ε))}` — all functions of `ε`.
*Lemma N.* For every `ε ≤ ε*`, `N_ε ⊂ Ω_ε`, and `u_N(b̂;ε) ≤ 65.4/b̂³ + 64.9ε/b̂² + 5.11ε/b̂`.
*Proof/certificate.* `r₀(ε)`, `r₁(ε)` are increasing in `ε`, so it suffices to certify `m⁴F − r₀(ε*) − r₁(ε*)u > 0` outside the
(`ε`-independent) set `{b̂ ≥ 120, u ≤ 1/100}`: (P) `J < 400`, all `u ∈ [0,1]`: exact Descartes on the shifted polynomials (596);
(C) `J ≥ 400`, `m ≤ 2`: `m⁴F ≥ 2.1·10⁴`; (S) `J ≥ 400`, `m ≥ 2`, `b̂ ≤ 120`: the cert_S cells with `r₀ + r₁u` subtracted (604 per
channel) and the small-`b̂` bound; (M) `b̂ ≥ 120`, `u ≥ 1/100`: the `Φ` term dominates. These margins then exceed
`r₀(ε) + r₁(ε)u` for every `ε ≤ ε*`. For `b̂ ≥ 120`, `u ≤ 1/100` the argument is made separately at each `ε`:
`m⁴F ≥ [80/b̂₊³ − C_R/b̂₋⁴] − 162u/b̂₋ + 2.5(1−η)b̂²u²` (both channels, §3), the bracket is `≥ 0` for `b̂ ≥ 120` (certified,
independent of `ε`), so `m⁴F < r₀(ε) + r₁(ε)u` forces `u < u_N(b̂;ε)`. Finally `u_N(b̂;ε) ≤ 162/(2.5(1−η)(1−δ)b̂³) + r₁(ε)/(2.5(1−η)b̂²)
+ (r₀(ε)/(2.5(1−η)))^{1/2}/b̂`, and `r₁(ε)/ε`, `r₀(ε)/ε²` are increasing, so the constants evaluated at `ε*` (`65.33`, `64.90`,
`5.103`) bound them for all `ε ≤ ε*`. ∎ The binding constraint is the `+−` near-zero at `m ≈ 21.8, J = 245`
(`m⁴F ≈ 2.67·10⁻⁵ > 162uε*`); `ε*` cannot be raised much above `≈ 7·10⁻⁵`. For larger `ε` the negative set contains ordinary
states (referee's float scan: at `ε = 1/100`, `++` states with `m ≈ 3.8–6.5M`, `J = 14–28`, and `+−` states near `b̂ ≈ 22.5`).

**Hypothesis `E(λ, η)` (eikonal model on `Ω_ε`).** Let `μ := λ m_IR` (`λ > 0`) be a *model* IR scale (not an observable: in
`D = 4` a positive eikonal measure needs an IR regulator). For every Borel `φ ≥ 0`,
`⟨(|c⁺⁺|² + |c⁺⁻|²)(2/m⁴) φ(b̂) 1_{Ω_ε}⟩ ≤ (1+η) ∫ φ dν_eik`, `dν_eik(b) = 16πG b log⁺(1/(μb)) db`.
`ν_eik` is the tree-level graviton-exchange density: `∫J₀(bp)dν_eik = 16πG(1 − J₀(p/μ))/p²`, i.e. the `16πG/p²` pole of `−B₂|_low`
for `p ≫ μ`; it is the transform of the 4D eikonal phase `δ_tree = −2Gs log(b/L_IR)`, valid for `b ≤ L_IR = 1/μ` (arXiv:2410.21499,
footnote 8; there the graviton pole at large `b` is generated by eikonal scattering, and the Tauberian analysis is for `d > 4`
with the IR part assumed calculable in the EFT). `E` constrains only states with `b̂ ≥ 120/M` and `s = m² ≥ 1/u_N ≥ 2.48·10⁴ M²`;
it involves no exact scaling-limit equality (finite-mass effects are bounded in Lemma IR), it does not force `N_ε` to be
empty, and `ν_eik` is what Theorem 2 requires at large `b`, so it is not vacuous. `η` is the tolerance for non-eikonal
(e.g. stringy) states in `Ω_ε`.

**Lemma IR.** For `ε ≤ ε*`, `1/λ < j₁,₁ ≈ 3.83`, under H1–H5 and `E(λ,η)`:
`Δ_IR ≤ (1+η)(16/k_lo)(65 I(λ) + 3.5ε/λ² + X(ε,λ))`, `I(λ) = ∫₀^{1/λ}J₁(x)log(1/(λx))dx = Σ_k (−1)^k λ^{−2k−2}/(2^{2k+1}k!(k+1)!(2k+2)²)`.
*Proof.* By Lemma N, `⟨[F^ε]₋⟩ = ⟨[F^ε]₋1_{Ω_ε}⟩` and per state `m⁴[F^ε]₋ ≤ [m⁴R_ε]₊`. Split
`m⁴R_ε = 2∫₀^εψ₂D_h + u∫₀^ε(σψ₃ − p²ψ₂)D_h`; the second part is `≤ r₁(ε)u`. For the first, compare `D_h(1−2up²)` with `J₀(b̂p)` at
`p ≤ ε`: by (Q1) in the form of Lemma S, `P_L(cos θ) = (2/π)∫cos(Nθ sin σ)E^{−1/2}dσ` with `1 ≤ E^{−1/2} ≤ θ/sin θ`, so
`|P_L(cos θ) − J₀(Nθ)| ≤ θ²/5`; `Nθ = b̂p s(up²)`, `|J₁| ≤ 0.582`, `s(z) − 1 ≤ z/5`, `θ² ≤ 4.1up²`; for `+−` (Lemma C) also
`|(1−y)^{−8} − 1| ≤ 8.1y` and the shifts `|jθ| ≤ 8.04p√u`. Hence `|D_h − J₀(b̂p)| ≤ 9.1up² + 0.12b̂up³ + 4.8p√u`, and with
`∫₀^εψ₂J₀(b̂p) ≤ 65εJ₁(εb̂)/b̂ + 14ε³` (`ψ₂ = pφ`, `0 ≤ 65 − φ ≤ 40p + 223p²`):
`m⁴[F^ε]₋ ≤ 2(65ε[J₁(εb̂)]₊/b̂ + 14ε³) + r₁(ε)u + 130(2.3uε⁴ + 0.024b̂uε⁵ + 1.6√u ε³)`, `u ≤ u_N(b̂;ε)`.
Apply `E` with `φ` equal to half this right side (a function of `b̂` once `u` is replaced by `u_N(b̂)`) and integrate:
`x = εb` gives `65 I(λ)`, `∫₀^{1/μ} b log(1/(μb))db = 1/(4μ²)` gives `3.5ε/λ²`, and the `u_N` terms give
`X(ε,λ) ≤ Λ·[ (r₁(ε)/2)∫b u_N + 65(2.3ε⁴∫b u_N + 0.024ε⁵∫b²u_N + 1.6ε³∫b√u_N) ]` over `b ∈ [120, L]`, `L = 1/μ = 1/(λε)`,
`Λ = log(L/120)`, with `u_N ≤ 65.5/b³ + 65.5ε/b² + 5.11ε/b`: `∫b u_N ≤ 65.5/120 + 65.5εΛ + 5.11/λ`, `∫b²u_N ≤ 65.5Λ + 65.5/λ +
2.56/(λ²ε)`, `∫b√u_N ≤ 2(65.5)^{1/2}(λε)^{−1/2} + (65.5ε)^{1/2}/(λε) + (2/3)(5.11ε)^{1/2}(λε)^{−3/2}`. Since `r₁(ε)/2 = ε·g(ε)` with `g`
positive and nondecreasing, every term of the expanded majorant is `c·g(ε)^{0 or 1}·ε^aΛ^k` with `k ∈ {1,2}`, `a ≥ k`
(the terms are `εΛ`, `ε²Λ²`, `ε⁴Λ`, `ε⁵Λ²`, `ε⁵Λ`, `ε^{5/2}Λ`, `ε²Λ`). As `d(ε^aΛ^k)/dε = ε^{a−1}Λ^{k−1}(aΛ − k) ≥ aε^{a−1}Λ^{k−1}(Λ − 1)`
and `Λ(ε) ≥ Λ(ε*) > 1` (checked for each `λ`), each term is increasing on `(0, ε*]`; hence `X(ε,λ) ≤ X(ε*,λ)`, and likewise
`3.5ε/λ²`. (With `r₀, r₁` frozen at `ε*` this would fail — the term `(r₁/2)(5.11/λ)Λ` would not vanish as `ε → 0`; referee
round 3.) ∎ [proved; `I(λ)`, `u_N`, `X` and the table certified in `cert_eps.py` (Arb)]

| `λ = μ/m_IR` | `I(λ)` | `X(ε*,λ)` | bound on `Δ_IR` (`η = 0`) | `B + Δ_IR` (`m_IR ≤ M/20000`) | model-dependent, in terms of `μ`: `B + A log λ + Δ_IR` |
|---|---|---|---|---|---|
| 2/7 | 1.06070 | 0.476 | 40.414 | −5.030 | −52.43 |
| 1/2 | 0.441919 | 0.254 | 16.870 | −28.574 | −54.80 |
| 1 | 0.121165 | 0.118 | 4.654 | **−40.790** | −40.79 |
| 2 | 0.031007 | 0.056 | 1.206 | −44.238 | −18.01 |
| 4 | 0.0077973 | 0.028 | 0.312 | −45.132 | +7.33 |

So under `E(1,0)`, for `m_IR ≤ M/20000`: `|ĝ₃|²M⁸ ≤ 37.8387 log(M/m_IR) − 40.790`; with `η > 0` add `4.654η`. At fixed `λ` the
coefficient of the logarithm is unchanged and the correction is O(1) in `B`. The last column re-expresses the bound through
the model scale `μ` (`m_IR = μ/λ`); its best entry, `37.84 log(M/μ) − 54.8` (`λ = 1/2`, `μ ≤ M/40000`), is obtained by optimising
over a model parameter, so it and any comparison with the paper's `−45.4` are model-dependent. The lemma is crude (it
discards `F > 0` inside `Ω_ε`), so these corrections are upper bounds.

## 5. Status of every statement

| statement | status |
|---|---|
| sum rules (2.26), heavy sides (2.33)/(2.35a), sign of `B₃` per channel, forward-limit actions | [proved analytically] from (2.16), (2.26) |
| low sides (2.31), (2.32), matter (3.8) (with the `m_ℓ⁴ − 6m_ℓ²p² + 6p⁴` correction), `k(0)`, `A`, `B`, matter signs, `F_{g4}` `g₅` term | [certified, exact, `cert_low.py`] |
| bracketing of (3.3a) | [settled] LaTeX source + (3.4a), (3.9) constants |
| `−5g₅` in v1 (3.3a); completion `+5B₅|₀`; misprint evidence from Fig. 3 | [proved, exact]; figures inspected (visual) |
| JHEP (3.3a): `−B₄^{(2)}|_low ≡ 0`, heavy action, negativity at `m = M`, `J = 4, 7, 9` | [certified, exact, `cert_jhep.py`] |
| Lemma P, Lemma C | [proved]; C also [certified, `cert_cg.py`] |
| (Q1)–(Q4) | (Q1),(Q3) DLMF 14.12.1, 18.14.2 (checked on the DLMF site); (Q2) proved from DLMF 18.10.2; (Q4) standard, numerically checked |
| Lemma M, Lemma S | [proved analytically] |
| regions (P), (S) | [certified, `cert_P.py`, `cert_S.py`] |
| regions (C), (M) | [proved] + [certified constants, `cert_C.py`, `cert_M.py`] |
| **Theorem 1** | [proved], given the above and correctness of FLINT/Arb |
| Proposition 1′ | [proved, exact] |
| Theorem 2, Theorem 3(i) | [proved] given H1–H5 |
| Lemma N (`N_ε ⊂ Ω_ε` for `ε ≤ 1/20000`) | [certified, `cert_eps.py`] |
| Lemma IR, Theorem 3(ii) (`m_IR ≤ M/20000`) | [proved] given H1–H5 and `E(λ,η)` on `Ω_ε` |
| analogous numbers for `M/20000 < m_IR ≤ M/100` | [estimate only — not proved] |
| former H-IR (`⟨[F^ε]₋⟩ = 0`) | **withdrawn** (incompatible with `G > 0`) |
| near-zeros, scans, `b⁴R` sizes | [numerical only — not part of any theorem] |
| optimal bound (4.4) `24.9 log − 27.6` (SDPB, improved sum rules) | not addressed |

## 6. Reproduction
`bash ads4/code/flatbound/run_all.sh [--sanity]`. Files: `common.py`, `cert_P.py`, `cert_C.py`, `cert_M.py`, `cert_S.py`,
`cert_cg.py`, `cert_low.py`, `cert_ir.py`, `cert_eps.py`, `cert_jhep.py` (its PASS means the stated findings, including the negativity of the
JHEP-printed functional, were verified), `sanity_numeric.py` (floating point, not part of the proof), `logs/`. Tools:
python-flint 0.9 (FLINT/Arb ball arithmetic; no `mpmath.iv`), sympy. Parameters `JP = 400`, `B_A = 120`, `u₁ = 1/100`.

## 7. What remains
* H1–H4 are physics input. A finite bound additionally needs a model of the large-`b` heavy spectrum (`E`), because in `D = 4` a
  functional positive at all `b` has a divergent graviton-pole term; `E` is a model hypothesis, not a theorem. A more
  rigorous route would replace the sharp IR cutoff by an AdS₄ regulator (cf. `UNIFORM_POSITIVITY.md`) or a Tauberian
  argument; neither is done here. Smearing of the forward-limit sum rules (H3) and the non-smooth IR cutoff are not justified.
* Theorem 1 relies on FLINT/Arb being correct and on (Q1)–(Q4).
* The optimal SDPB bound (4.4) and positivity of `F_{g4}` are not treated (the code applies with new weights).
* Which forward combination the authors actually used is inferred from their figures and constants, not known.
