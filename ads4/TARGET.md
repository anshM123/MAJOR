# AdS4/CFT3 programme: selected target, exact statement, missing theorem, first analysis

Status: research in progress. Nothing in this file is refereed. Sources and the full candidate list:
`ads4/LITERATURE.md`, `ads4/PROBLEM_CANDIDATES.md` (literature audit of 15 Nov 2026).

## 1. Ranking (five criteria, in the order given by the brief)

| Rank | Candidate | Outcome type | Verdict |
|---|---|---|---|
| 1 | P1: HPPS locality for the stress-tensor four-point function in d = 3 | condition for an Einstein dual | most significant; **no realistic route** with the available means |
| 2 | **P2: universal long-range attraction from stress-tensor exchange, as a theorem** | unproved aspect of emergent gravitational dynamics | **selected as primary target**; the simplest route is blocked (see §4) |
| 3 | P6: `N^{3/2}` law for the ABJM sphere partition function, as a theorem | nonperturbative datum | **fallback**; realistic (see §5) |
| 4–10 | Penrose inequality in AdS4; HRT cone in pure AdS4 on arbitrary slices; QES-cone conjecture; a positivity lemma; split property; strips in AdS4 (corollary-level); 2+1 with matter | — | lower on significance or tractability; the entropy-cone items are consistency constraints on states that already have a geometric dual, not derivations of gravity |

The audit's own caveat applies to everything tractable on the list: P2 and P6 "turn accepted physics into
theorems"; neither promises a new physical conclusion.

## 2. Exact statement of the primary target

Setting: a unitary CFT in `d = 3`; `G(z, z̄) = Σ p_{τ,ℓ} g_{τ,ℓ}(z, z̄)` the four-point function of identical
scalars of dimension `Δ`, `p ≥ 0`, discrete spectrum, crossing
`((1−z)(1−z̄))^Δ G(z,z̄) = (z z̄)^Δ G(1−z,1−z̄)`, blocks normalised by
`g_{τ,ℓ}(z,z̄) = z^{τ/2}[k_{τ+2ℓ}(z̄) + O(z)]`, `k_{2h}(x) = x^h ₂F₁(h,h;2h;x)`.
Hypothesis (H): besides the identity, the operators of twist `< τ_* + δ` are finitely many, all of twist `τ_*`;
`b_* := Σ_{τ_O = τ_*} p_O Γ(τ_*+2ℓ_O)/Γ(τ_*/2+ℓ_O)²`. Measure on twists: `μ_z := z^Δ Σ p_{τ,ℓ} k_{τ+2ℓ}(1−z) δ_τ`.

**Target statement (tilted form).** There is a window `y = y(z) → 0` such that, with `x = τ − 2Δ`,

`∫_{|x| < η} x · y^{x/2} dμ_z(τ) = −2 b_* z^{τ_*/2} (1 + o(1))`,   `z → 0`,

and `∫_{|x|<η} x² y^{x/2} dμ_z = o(z^{τ_*/2})`. For stress-tensor exchange `τ_* = 1`, `b_* ∝ Δ²/c_T`.
Meaning: in every such CFT3 the two-particle states of large spin are bound, on average, by an energy
`∝ −(1/c_T) ℓ^{−1}` — Newtonian attraction in the emergent `AdS_4`, with `G_N ∝ 1/c_T`, derived from
unitarity and crossing alone, at finite `N`.

Known at the level of derivation since arXiv:1212.3616, 1212.4103. Known as theorems: accumulation of twists
at `2Δ` (Pal–Qiao–Rychkov, arXiv:2212.04893) and convergence of spectrum and coefficients to the generalised
free values (van Rees, arXiv:2412.06907), **without the first correction**.

## 3. What is proved here: a Tauberian lemma (elementary, complete)

**Lemma.** Let `ν` be a positive measure on `[−η, η]`, `L(s) := ∫ e^{−sx} dν(x)`, and let `S > 0`, `0 < h ≤ 1`.
Suppose `L(s) = m + ε(c + 2bs) + r(s)` with `|r(s)| ≤ ρ` for `s ∈ {S−h, S, S+h}`. Then
(a) `∫ x² e^{−Sx} dν ≤ 4ρ/h²`;
(b) `| ∫ x e^{−Sx} dν + 2bε | ≤ ρ/h + (2/3) η cosh(hη) ρ`.

*Proof.* (a) `L(S+h) + L(S−h) − 2L(S) = ∫ e^{−Sx}(e^{hx} + e^{−hx} − 2) dν ≥ h² ∫ x² e^{−Sx} dν`, and the left side is
`r(S+h) + r(S−h) − 2r(S) ≤ 4ρ` because the affine part cancels. (b) `L(S−h) − L(S+h) = ∫ e^{−Sx}(e^{hx} − e^{−hx}) dν
= −4bεh + r(S−h) − r(S+h)`, and `|e^{hx} − e^{−hx} − 2hx| ≤ (h³|x|³/3) cosh(hη) ≤ (h³ η x²/3) cosh(hη)`; use (a). ∎

So the target follows once crossing is shown to give `L` in the stated form with `ρ/h = o(ε)`, `ε = z^{τ_*/2}`,
on three points of a window. Positivity does the rest; no assumption on the shape of the spectrum inside
the window is needed.

## 4. The missing theorem, and why the simplest route fails for gravity in d = 3

With `z̄ = 1 − y`, crossing reads

`z^Δ Σ p_{τ,ℓ} y^{τ/2−Δ}[k_{τ+2ℓ}(1−z) + O(y)] = ((1−z)/(1−y))^Δ G(z, 1−y)`,

and the direct channel gives `G(z,1−y) = 1 + z^{τ_*/2} Σ_* p_O k_{τ_*+2ℓ_O}(1−y) + R(z,y)`, with
`k_{2h}(1−y) = −(Γ(2h)/Γ(h)²)(log y + 2ψ(h) − 2ψ(1)) + O(y log y)`; the `log y` term is `2 b_* s`, `y = e^{−2s}`.

Two requirements pull in opposite directions.
* *Crossed channel.* The next double-twist tower sits at `x ≈ 2` with mass of order one and contributes
  `≈ y`; it is negligible against `ε` only if `y ≪ z^{τ_*/2}`.
* *Direct channel.* The remainder `R` (twists `≥ τ_* + δ`) must be `o(z^{τ_*/2})` uniformly in the window. The
  bound available from positivity alone is `R ≤ (z/z')^{(τ_*+δ)/2} G(z', 1−y) ≲ z^{(τ_*+δ)/2} y^{−Δ}`, which is
  `o(z^{τ_*/2})` only if `y ≫ z^{δ/(2Δ)}`.

They are compatible iff `δ > Δ τ_*`. Hypothesis (H) itself needs `δ < 2Δ − τ_*` (otherwise the double-twist
tower lies below `τ_* + δ`). Hence the crude route needs `Δ > τ_*/(2 − τ_*)` (the audit's condition) **and**
a twist gap `δ ∈ (Δτ_*, 2Δ − τ_*)`.

**For the stress tensor in `d = 3` this window is empty in any theory where the double-stress-tensor
operators `[TT]`, of twist `≥ 2`, appear in `φ × φ`:** then `δ ≤ 1`, while `δ > Δ` and `Δ > 1` are required.
(For exchange of a scalar of small twist the window is non-empty, but that is not gravity.)

**Missing theorem (precise).** *A second-order double-lightcone estimate:* for `y` in a window
`z^{a} ≤ y ≤ z^{b}` with `b > τ_*/2`,

`R(z,y) − [(z/y)^Δ-terms that are the images of the crossed-channel identity and of the crossed-channel
operators of twist τ_*] = o(z^{τ_*/2})` uniformly.

The crude bound fails because `R` is dominated by exactly those image terms, which are explicit; what is
needed is a bound on `R` *after* subtracting them. By the symmetry of the crossing equation under `z ↔ y`,
that is the same statement in the mirror channel with the roles of the variables exchanged, which suggests
an iterative ("bootstrap of estimates") argument: an a-priori bound of the form `C z^{α} y^{β}` on the
subtracted remainder in one channel, inserted in the other, improves the exponents. Whether the iteration
closes for `τ_* = 1`, `d = 3` is the open question. Next step: extract the exact form of the first-order
double-lightcone estimates of Pal–Qiao–Rychkov and van Rees and test whether their method tolerates the
subtraction.

## 5. Fallback (P6), with the route already visible

`Z(N,k) = e_N(λ_1, λ_2, …)`, the `N`-th elementary symmetric function of the eigenvalues of the positive
trace-class operator `ρ = A^{1/2} B A^{1/2}`, `A = 1/(2cosh(x/2))`, `B = 1/(2cosh(p/2))`, `[x,p] = 2πik`.
* Exact identity: `Tr A^s B^s = (2πħ)^{-1} (∫ (2cosh(x/2))^{-s} dx)²`, `ħ = 2πk`, which is `2C/s² (1 + o(1))` as
  `s → 0`, `C = 2/(π²k)`.
* Araki–Lieb–Thirring: `Tr ρ^s ≥ Tr A^s B^s` for `0 < s ≤ 1`: a rigorous lower bound with the conjectured
  constant.
* Missing: the matching upper bound on `Tr ρ^s` as `s → 0` (equivalently an upper bound on the eigenvalue
  counting function). By the Birman–Schwinger principle it is a sharp bound on the number of non-positive
  eigenvalues of `2cosh(p/2) − Λ/(2cosh(x/2))`, of the type proved by Laptev–Schimmer–Takhtajan
  (arXiv:1510.00045) for a neighbouring operator.
* Given the Weyl law `n(E) = CE² + O(E)`, the bounds `Π_{n≤N} λ_n ≤ e_N ≤ inf_μ e^{−μN} Π(1 + e^{μ}λ_n)` give
  `log Z(N,k) = −(π√(2k)/3) N^{3/2} + O(N)`.
