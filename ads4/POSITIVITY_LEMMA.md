# A positivity lemma for the AdS4 bulk-locality bound (d = 3)

Status: proved here 16 Nov 2026; checked by one independent automated referee with a hostile brief (same
model family as the author; not a human expert): "I could not break it. The mathematics is correct;
established after trivial fixes" — a typo, an equation number, the citation for Mehler's formula, and the
identification with the paper's (B.15), all applied below. The referee verified Mehler's formula against
the hypergeometric definition to `10^{−29}`, the finite-`κ` identity to `4·10^{−17}`, the lemma against its own
finite-part evaluation of (B.15) to `10^{−8}`–`10^{−13}`, `Q_1(0) = 1/24` to 16 digits, and the sharpness at
`n = 1`. Novelty: no proof found in the literature, with a limited search; the result is elementary.
Scope: a lemma inside the large-gap locality programme (Heemskerk–Penedones–Polchinski–Sully), for the
*scalar* four-point function with graviton exchange. It is **not** a result about the stress-tensor
four-point function and does not by itself make any locality bound a theorem.

## Context

Caron-Huot, Mazáč, Rastelli, Simmons-Duffin, "AdS bulk locality from sharp CFT bounds" (arXiv:2106.10274)
obtain the AdS4 bound (their eq. (5.64); the Regge limit of the functional is (5.63)) on the leading higher-derivative contact coupling of a scalar coupled
to gravity. Positivity of their functional in the regime where the AdS curvature matters (impact parameter
`β ~ 1`) is reduced to the sign of the "ν-moments"

`Q_n(β) := ∫_0^∞ dν ν^n P_{(2−d)/2 + iν}(cosh β)`   (their (B.14); defined by analytic continuation in `n`),

and is stated on numerical grounds: "By checking numerical examples, we find that `Q_n(β)` is everywhere
positive as a function of `β > 0` for `0 ≤ n ≤ 2` in `d ≥ 4`, and for `0 ≤ n ≤ 1` in `d = 3`"; for AdS4 they use
"Both are everywhere positive" for `Q_1` and `−Q_2`, and call the positivity of `Q_1` "the crucial distinction
between AdS4 and the momentum cutoff".

## Statement

Let `d = 3`, so that the harmonic function is the conical Legendre function `P_{−1/2+iν}(cosh β)`. For
`−1 < n < −1/2` the integral defining `Q_n(β)` converges absolutely; `Q_n(β)` denotes its analytic continuation
in `n`.

**Lemma Q.** For `β > 0` and `Re n > −1`,

`Q_n(β) = (√2/π) ∫_β^∞ K_n(t) (cosh t − cosh β)^{−1/2} dt`,

`K_n(t) := Γ(n+1) [ cos(πn/2) t^{−n−1} + 2 Σ_{m ≥ 1} sin((n+1) θ_m) r_m^{−n−1} ]`, `θ_m := arctan(t/(2πm))`,
`r_m := (4π²m² + t²)^{1/2}`.

Consequences, for every `β > 0`:
1. `Q_n(β) > 0` for all real `−1 < n ≤ 1`. In particular `Q_1 > 0`, with
   `K_1(t) = Σ_{m≥1} 8πm t/(4π²m² + t²)² = −Im ψ′(1 + it/(2π))/(2π²) > 0`.
2. `Q_0(β) = 1/(2 sinh(β/2))` and `Q_2(β) = −1/(8 sinh³(β/2)) < 0`.
3. More generally `K_{n+2} = −K_n″`, so `sign Q_{2k} = (−1)^k`.

So both sign conditions used for AdS4 in arXiv:2106.10274 (`Q_1 > 0`, `−Q_2 > 0`) hold for all `β > 0`, and
the range `0 ≤ n ≤ 1` found numerically there is sharp at the upper end (for `n` slightly above `1` the term
`cos(πn/2) t^{−n−1}` is negative and dominates at small `β`).

## Proof

*Step 1 (Mehler's formula).* For `ν > 0`, `β > 0`:
`P_{−1/2+iν}(cosh β) = (√2/π) coth(πν) ∫_β^∞ sin(νt) (cosh t − cosh β)^{−1/2} dt`.
Here `P_{−1/2+iν}(x) = ₂F₁(½−iν, ½+iν; 1; (1−x)/2)`, which is the paper's harmonic function (their (2.13)) at
`d = 3`, normalised to `1` at `x = 1`, and their measure is plain `dν`. A checkable derivation of the formula:
the connection formula DLMF 14.20.6 at `μ = 0`, `P_{−1/2+iν} = (i coth(πν)/π)(Q_{−1/2+iν} − Q_{−1/2−iν})`, combined
with the Laplace-type integral `Q_σ(cosh β) = ∫_β^∞ e^{−(σ+½)t}(2cosh t − 2cosh β)^{−1/2} dt`. (Equation numbers in
Gradshteyn–Ryzhik and Lebedev could not be confirmed and are not relied on; both ingredients were verified
numerically to `10^{−26}` by the referee.)

*Step 2 (Abel regularisation).* For `κ > 0` and real `n > −1` put `F_κ(n) := ∫_0^∞ ν^n e^{−κν} P_{−1/2+iν}(cosh β) dν`.
The double integral obtained by inserting Step 1 converges absolutely (`ν^n coth(πν)|sin νt| ≤ C ν^n` near `0`,
exponential decay at infinity, and `(cosh t − cosh β)^{−1/2}` is integrable on `(β,∞)`), so by Fubini

`F_κ(n) = (√2/π) ∫_β^∞ K_n^κ(t) (cosh t − cosh β)^{−1/2} dt`, `K_n^κ(t) := ∫_0^∞ ν^n e^{−κν} coth(πν) sin(νt) dν`.

With `coth(πν) = 1 + 2 Σ_{m≥1} e^{−2πmν}` and `∫_0^∞ ν^n e^{−aν} sin(tν) dν = Γ(n+1) Im (a − it)^{−n−1}
= Γ(n+1) sin((n+1) arctan(t/a)) (a² + t²)^{−(n+1)/2}` (term-by-term integration is justified because
`Σ_m ∫ ν^n e^{−2πmν}|sin tν| dν ≤ t Σ_m Γ(n+2)(2πm)^{−n−2} < ∞`):

`K_n^κ(t) = Γ(n+1) Σ_{m ≥ 0} c_m sin((n+1) arctan(t/(κ + 2πm))) ((κ+2πm)² + t²)^{−(n+1)/2}`, `c_0 = 1`, `c_{m≥1} = 2`.

*Step 3 (`κ → 0`).* For fixed `t > 0` each term converges, the `m = 0` term to
`sin((n+1)π/2) t^{−n−1} = cos(πn/2) t^{−n−1}`, and the series is dominated uniformly in `κ ∈ (0,1]` and in
`t ≥ β` by `C(1 + t)^{|n|+1}` (bound `|sin((n+1)θ)| ≤ min(1, (n+1) t/(2πm))`). Since
`(cosh t − cosh β)^{−1/2} ≤ C e^{−t/2}` for large `t`, dominated convergence gives
`lim_{κ→0} F_κ(n) = (√2/π) ∫_β^∞ K_n(t) (cosh t − cosh β)^{−1/2} dt =: R_n(β)`.

*Step 4 (identification with `Q_n`).* For `−1 < n < −1/2`, `|P_{−1/2+iν}(cosh β)| ≤ C ν^{−1/2}` at large `ν`, so
`∫ ν^n |P| dν < ∞` and `F_κ(n) → Q_n(β)` by dominated convergence; hence `Q_n = R_n` there. The same estimates
with complex `n` show that `R_n(β)` is analytic in `n` on `Re n > −1`. Two analytic functions agreeing on an
interval agree on the connected domain, so `R_n` is the analytic continuation of `Q_n` to `Re n > −1`. (In fact Steps 2–3 show directly that `R_n` is the
Abel limit of the moment for every real `n > −1`, so no continuation is needed to reach `n = 1, 2`.)

*Step 4′ (identification with the paper's (B.15)).* The cosine Mehler–Dirichlet formula
`P_{−1/2+iν}(cosh β) = (√2/π)∫_0^β cos(νt)(cosh β − cosh t)^{−1/2} dt` with the same Abel regulator gives, for
`−1 < n < 0`, `lim_{κ→0} F_κ(n) = −(√2/π) Γ(n+1) sin(πn/2) ∫_0^β t^{−n−1}(cosh β − cosh t)^{−1/2} dt`, which is their
(B.15) at `d = 3`. Both sides are analytic in `n` (the finite-part poles at `n = 0, 2, …` cancel against
`sin(πn/2)`), so (B.15), as continued by the authors, equals `R_n`.

*Step 5 (signs).* For real `−1 < n ≤ 1`: `cos(πn/2) ≥ 0`, and `0 < (n+1)θ_m < (n+1)π/2 ≤ π`, so every term of
`K_n(t)` is non-negative and those with `m ≥ 1` are positive. Hence `K_n > 0` and `Q_n > 0`. For `n = 1`,
`sin 2θ_m = 4πm t/r_m²`, giving the series for `K_1`; and `Σ_m m/(m² + a²)² = −Im ψ′(1+ia)/(2a)`.
The identity `Im (a − it)^{−n−3} Γ(n+3) = −∂_t² [Γ(n+1) Im (a − it)^{−n−1}]` gives `K_{n+2} = −K_n″`. With
`K_0(t) = 1/t + Σ_m 2t/(4π²m² + t²) = ½ coth(t/2)`: `K_2(t) = −cosh(t/2)/(4 sinh³(t/2)) < 0`.

*Step 6 (closed forms).* With `w = sinh(t/2)`, `w_0 = sinh(β/2)`, `cosh t − cosh β = 2(w² − w_0²)`:
`Q_0 = (√2/π) ∫_{w_0}^∞ dw/(√2 w √(w² − w_0²)) = 1/(2w_0)`, and
`Q_2 = −(√2/π) (1/(2√2)) ∫_{w_0}^∞ dw/(w³ √(w² − w_0²)) = −(1/(2π)) · π/(4w_0³) = −1/(8 w_0³)`. ∎

## Numerical confirmation (`ads4/code/qn_positivity.py`)

The paper's own definition (their (B.15) at `d = 3`: a finite-part integral over `t ∈ [0, β]`, an independent
representation from the one above) against the formula of Lemma Q:

| `n` | `β` | agreement |
|---|---|---|
| 1 | 0.05, 0.5, 2, 6, 15 | relative difference `10^{−4}` to `3·10^{−11}` (limited by the reference evaluation at small `β`) |
| 2 | 0.05 … 15 | the integral of Lemma Q equals the closed form `−1/(8 sinh³(β/2))` to `≤ 6·10^{−12}` (this row compares the integral with the closed form, not with (B.15); the referee compared with (B.15) at `n = 1.999, 2.001`) |
| 0 | 0.5, 3 | 12 digits, `= 1/(2 sinh(β/2))` |
| 0.5 | 1, 4 | `10^{−8}` |

`Q_1(0.01) = 0.0416662` (limit `1/24`, the value the literature audit had found independently). `K_1 > 0.0019` on
`0 < t ≤ 60`.

## Other dimensions (referee)

For `d ≥ 5` and `0 < n < 2`, positivity follows from (B.15) by one integration by parts: with
`g(u) = (cosh β − cosh uβ)^{(d−4)/2}`, `g(1) = 0` and `g′ ≤ 0`, so `f.p.∫_0^1 u^{−n−1} g = (1/n)∫_0^1 u^{−n} g′ < 0` and
`Q_n > 0` (checked in `d = 5` against direct Abel moments to 5–6 digits). In `d = 4` the functions are elementary,
`Q_n = sin(πn/2)Γ(n)/(β^n sinh β)`, with `Q_2 ≡ 0`. So `d = 3` is the only case in which (B.15) has no definite
sign term by term, and it is the case treated by Lemma Q.

## What this does and does not establish

* In (5.63) the functional's limit is `a_1 Q_1/M² + a_2 Q_2/M³` with `a_1 = 1310 > 0`, `a_2 = −1080 < 0`, so
  "`Q_1 > 0` and `Q_2 < 0` for all `β > 0`" is exactly the sign input used there (referee's reading of the paper).
* **Still unproved in the same AdS4 passage:** monotonicity of the ratio in their (5.35), stated "by
  inspection" (confirmed numerically for `d = 3`, `β ∈ [0.1, 30]`, by the referee); uniformity in `β` of the
  subleading and oscillatory terms (their footnote 27 and appendix F.1); positivity of the flat-space
  functional (5.55), which is numerical.
* A slip in the paper noticed by the referee: at large `β` the correction to the leading term of (5.44) is
  `1 − 2n ln 2/β`, not `O(1/β²)`.


* It replaces "by checking numerical examples" by a proof, for `d = 3`, of the positivity statements the
  AdS4 bound relies on in the curved regime.
* The rest of the derivation of that bound (the dispersive sum rules, the Regge and large-`M` asymptotics of
  the functional, the treatment of light states, the flat-space extremal functional found numerically)
  remains at the physics level of rigour.
* `d ≥ 4` (`0 ≤ n ≤ 2`) is not treated; in `d = 4` the functions are elementary and the statement is immediate.
* Nothing here concerns external stress tensors.
