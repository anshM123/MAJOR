# A positivity lemma for the AdS4 bulk-locality bound (d = 3)

Status: proved here 16 Nov 2026; numerical confirmation in `ads4/code/qn_positivity.py`; **not yet refereed**.
Scope: a lemma inside the large-gap locality programme (Heemskerk–Penedones–Polchinski–Sully), for the
*scalar* four-point function with graviton exchange. It is **not** a result about the stress-tensor
four-point function and does not by itself make any locality bound a theorem.

## Context

Caron-Huot, Mazáč, Rastelli, Simmons-Duffin, "AdS bulk locality from sharp CFT bounds" (arXiv:2106.10274)
obtain the AdS4 bound (their eq. (236)) on the leading higher-derivative contact coupling of a scalar coupled
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
`P_{−1/2+iν}(cosh β) = (√2/π) coth(πν) ∫_β^∞ sin(νt) (cosh t − cosh β)^{−1/2} dt`
(classical; e.g. Gradshteyn–Ryzhik 8.842; Lebedev, *Special Functions*, §7.4 — **cited, not re-derived**;
its consequences are checked numerically below).

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
interval agree on the connected domain, so `R_n` is the analytic continuation of `Q_n` to `Re n > −1`.

*Step 5 (signs).* For real `−1 < n ≤ 1`: `cos(πn/2) ≥ 0`, and `0 < (n+1)θ_m < (n+1)π/2 ≤ π`, so every term of
`K_n(t)` is non-negative and those with `m ≥ 1` are positive. Hence `K_n > 0` and `Q_n > 0`. For `n = 1`,
`sin 2θ_m = 4πm t/r_m²`, giving the series for `K_1`; and `Σ_m m/(m² + a²)² = −Im ψ′(1+ia)/(2a)`.
The identity `Im (a − it)^{−n−3} Γ(n+3) = −∂_t² [Γ(n+1) Im (a − it)^{−n−1}]` gives `K_{n+2} = −K_n″`. With
`K_0(t) = 1/t + Σ_m 4t/(4π²m² + t²) = ½ coth(t/2)`: `K_2(t) = −cosh(t/2)/(4 sinh³(t/2)) < 0`.

*Step 6 (closed forms).* With `w = sinh(t/2)`, `w_0 = sinh(β/2)`, `cosh t − cosh β = 2(w² − w_0²)`:
`Q_0 = (√2/π) ∫_{w_0}^∞ dw/(√2 w √(w² − w_0²)) = 1/(2w_0)`, and
`Q_2 = −(√2/π) (1/(2√2)) ∫_{w_0}^∞ dw/(w³ √(w² − w_0²)) = −(1/(2π)) · π/(4w_0³) = −1/(8 w_0³)`. ∎

## Numerical confirmation (`ads4/code/qn_positivity.py`)

The paper's own definition (their (B.15) at `d = 3`: a finite-part integral over `t ∈ [0, β]`, an independent
representation from the one above) against the formula of Lemma Q:

| `n` | `β` | agreement |
|---|---|---|
| 1 | 0.05, 0.5, 2, 6, 15 | relative difference `10^{−4}` to `3·10^{−11}` (limited by the reference evaluation at small `β`) |
| 2 | 0.05 … 15 | `≤ 6·10^{−12}`; both equal `−1/(8 sinh³(β/2))` |
| 0 | 0.5, 3 | 12 digits, `= 1/(2 sinh(β/2))` |
| 0.5 | 1, 4 | `10^{−8}` |

`Q_1(0.01) = 0.0416662` (limit `1/24`, the value the literature audit had found independently). `K_1 > 0.0019` on
`0 < t ≤ 60`.

## What this does and does not establish

* It replaces "by checking numerical examples" by a proof, for `d = 3`, of the positivity statements the
  AdS4 bound relies on in the curved regime.
* The rest of the derivation of that bound (the dispersive sum rules, the Regge and large-`M` asymptotics of
  the functional, the treatment of light states, the flat-space extremal functional found numerically)
  remains at the physics level of rigour.
* `d ≥ 4` (`0 ≤ n ≤ 2`) is not treated; in `d = 4` the functions are elementary and the statement is immediate.
* Nothing here concerns external stress tensors.
