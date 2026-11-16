# Uniform positivity of `F_M(β) = ∫_0^1 f(p) P_{−1/2+iMp}(cosh β) dp` for the weight `f(p) = (1−p)³ p (65 + 155p + 47p²)`

Status: written 16 Nov 2026. Computer-assisted proof (interval arithmetic, `ads4/code/unipos/`, driver
`run_all.sh`, about 10 s for the certificates). **Not yet refereed by anyone other than the author** (an
automated analyst; no human check). The statement is a theorem about an explicit integral; nothing here depends
on physics, and nothing here says that `F_M` is the complete action of any functional on heavy states (see
`EINSTEIN_PROGRAM.md` §3b for that caveat).

## 0. Results

Throughout, `f(p) = (1−p)³ p (65 + 155p + 47p²) = 65p − 40p² − 223p³ + 259p⁴ − 14p⁵ − 47p⁶`, `φ(p) := f(p)/p`,
`P_{−1/2+iν}(cosh β) = ₂F₁(½−iν, ½+iν; 1; −sinh²(β/2))`, and

`F_M(β) := ∫_0^1 f(p) P_{−1/2+iMp}(cosh β) dp`,  `f̂(b) := ∫_0^1 f(p) J_0(pb) dp`.

**Theorem 1 (curved).** `F_M(β) > 0` for **all** `M > 0` and all `β > 0`. (So the required statement holds with
any `M_0 > 0`; there is no regime restriction.)

**Theorem 2 (flat-space lemma).** `f̂(b) > 0` for all `b ≥ 0`.

**Theorem 1′ (second, independent proof of a weaker statement).** `F_M(β) > 0` for all `β > 0` and `M ≥ 45`, using
only Theorem 2-type one-dimensional certificates and the comparison (I1) of the task statement.

What the proofs rest on, apart from the arguments written out below:

* three classical integral representations, (R1)–(R3) in §1 (Mehler–Dirichlet, Mehler's formula with `coth`,
  Mehler–Sonine). They are quoted, not re-proved; each was re-verified numerically here to `≤ 2·10^{−31}`
  ((R1), (R2), 30-digit arithmetic) and `3·10^{−7}` ((R3), through its consequence (2.3), limited by a truncated
  tail) — `sanity_numeric.py`, part 1;
* five interval-arithmetic certificates C1–C5 (`mpmath.iv`, outward rounding, 120-bit precision; exact rational
  arithmetic for all coefficients). Correctness of `mpmath.iv`'s enclosures of `+ − × ÷ sin cos exp π` is assumed
  (spot-checked by `iv_selftest.py`, which is a sanity check and not a proof).

No statement of Theorems 1, 2, 1′ rests on a floating-point grid. Floating-point scans are reported in §6 and
are labelled as such.

**Structure of the proof of Theorem 1** (`b := Mβ`; the three regions overlap on their common boundaries):

| region | argument | certificates |
|---|---|---|
| (A) `Mβ ≥ 10.85`, every `M > 0` | exterior (Mehler) representation; the kernel `K_{f,M}(t) = S(Mt) + R_M(t)` is positive for `Mt ≥ 10.85` because `S > 0` there and `R_M > 0` always | C1 |
| (B) `0 < M ≤ 6`, every `β > 0` | the kernel `K_{f,M}(t)` is positive for **all** `t > 0`: for `M ≤ 1.15` by a monotonicity argument, for `1.15 ≤ M ≤ 6` because `S > 0` outside the window `Mt ∈ (9.6, 10.85)` and an explicit lower bound `K_low` is certified positive on the window | C1, C2, C3 |
| (C) `M ≥ 6`, `Mβ ≤ 10.85` (so `β ≤ 1.8084`) | exact convergent expansion `F_M(β) = Σ_n β^{2n} G_n(Mβ)` with explicit tail; positivity of the two-variable function certified on `{b ≤ 10.85, β ≤ b/6}`, which includes the flat slice `β = 0` | C4 |

The three-region scheme suggested in the task (flat / transition `β ~ M^{−1/3}` / curved, with explicit bounds on
`Q_3`, `Q_5`, `T_M`) turned out to be unnecessary: the sine transform `S` of `f` is negative only on one short
window `(9.6105, 10.8282)`, so beyond `Mβ = 10.85` positivity is pointwise in the Mehler kernel for every `M`, and
below it the dependence on the curvature is captured exactly by a power series in `β²` at fixed `Mβ`.

## 1. Classical inputs and exact data

**(R1) Mehler–Dirichlet.** For real `ν` and `β > 0`:
`P_{−1/2+iν}(cosh β) = (√2/π) ∫_0^β cos(νt) (cosh β − cosh t)^{−1/2} dt`.

**(R2) Mehler's formula.** For `ν > 0`, `β > 0`:
`P_{−1/2+iν}(cosh β) = (√2/π) coth(πν) ∫_β^∞ sin(νt) (cosh t − cosh β)^{−1/2} dt`.
(This is Step 1 of `POSITIVITY_LEMMA.md`, where a derivation from DLMF 14.20.6 and the Laplace integral for `Q_σ`
is given.)

**(R3) Mehler–Sonine.** For `x > 0`: `J_0(x) = (2/π) ∫_1^∞ sin(xu) (u² − 1)^{−1/2} du` (conditionally convergent);
and Bessel's integral `J_0(x) = (2/π) ∫_0^{π/2} cos(x sin θ) dθ`.

References: Erdélyi et al., *Higher Transcendental Functions* I, §3.7; DLMF §14.12, §14.20, §10.9. Equation numbers
were not checked against the sources (no access); the identities themselves were checked numerically, see §0.

**Exact data** (all verified symbolically/exactly in `common.py`, `cert_smallM.py`, and to `2·10^{−15}` against
quadrature in `sanity_numeric.py`):

* `f ≥ 0` on `[0,1]`; `∫_0^1 f = 2591/420 = 6.16905`; `∫_0^1 p f = 283/120 = 2.35833`.
* `f^{(j)}(0) = 0, 65, −80, −1338, 6216, −1680, −33840` and `f^{(j)}(1) = 0, 0, 0, −1602, −12384, −35520, −33840` for
  `j = 0,…,6`.
* For a polynomial `f` of degree 6 and `z ≠ 0`: `∫_0^1 f(p) e^{−zp} dp = Σ_{j=0}^{6} [f^{(j)}(0) − e^{−z} f^{(j)}(1)] z^{−j−1}`
  (integration by parts). With `z = −iτ` this gives the closed forms
  - `S(τ) := ∫_0^1 f(p) sin(pτ) dp = 80/τ³ + 6216/τ⁵ + 33840/τ⁷ + (1602/τ⁴ − 35520/τ⁶) sin τ + (12384/τ⁵ − 33840/τ⁷) cos τ`,
  - `f̃(τ) := ∫_0^1 f(p) cos(pτ) dp = −65/τ² − 1338/τ⁴ + 1680/τ⁶ + (1602/τ⁴ − 35520/τ⁶) cos τ + (−12384/τ⁵ + 33840/τ⁷) sin τ`.
* `φ′(p) = −(1−p)² (40 + 526p + 235p²)`, so `φ` is positive and strictly decreasing on `[0,1)`.

## 2. Lemmas and proofs

### Lemma 1 (sine transform of a decreasing function) — [proved analytically]

Let `h : (0,1) → (0,∞)` be bounded and strictly decreasing. Then `∫_0^1 h(p) sin(px) dp > 0` for every `x > 0`.

*Proof.* Extend `h` by `0` to `[1,∞)`; the extension `h̄` is non-increasing on `(0,∞)`. With `δ = π/x`,
`∫_0^∞ h̄(p) sin(px) dp = Σ_{k≥0} ∫_{2kδ}^{(2k+1)δ} [h̄(p) − h̄(p+δ)] sin(px) dp`. Every term is `≥ 0`, and the term `k = 0` is
`> 0` because for `0 < p < min(δ,1)` one has `h̄(p) > h̄(p+δ)` and `sin(px) > 0`. ∎

### Lemma 2 (exterior representation; the kernel is `S` plus a positive function) — [proved analytically, given (R2)]

For `M > 0`, `β > 0`:

`F_M(β) = (√2/π) ∫_β^∞ K_{f,M}(t) (cosh t − cosh β)^{−1/2} dt`, `K_{f,M}(t) := ∫_0^1 f(p) coth(πMp) sin(Mpt) dp`,   (2.1)

and `K_{f,M}(t) = S(Mt) + R_M(t)` with `R_M(t) := ∫_0^1 f(p) (coth(πMp) − 1) sin(Mpt) dp > 0` for all `t > 0`.   (2.2)

*Proof.* Insert (R2) with `ν = Mp` into the definition of `F_M`. `|f(p) coth(πMp) sin(Mpt)| ≤ f(p)(1 + 1/(πMp))` is
bounded on `(0,1]` (because `f(p) ≤ 65p`), and `(cosh t − cosh β)^{−1/2}` is integrable on `(β,∞)` (an inverse square
root singularity at `β`, exponential decay at infinity), so the double integral converges absolutely and Fubini
gives (2.1). For (2.2), `coth(πMp) − 1 = 2/(e^{2πMp} − 1)` and
`f(p)(coth(πMp) − 1) = φ(p) · 2p/(e^{2πMp} − 1)` is a product of two positive, strictly decreasing, bounded functions of
`p ∈ (0,1)` (`x/(e^x − 1)` is decreasing); apply Lemma 1 with `x = Mt`. ∎

The same computation with (R3) in place of (R2) gives the flat analogue: for `b > 0`

`f̂(b) = (2/π) ∫_b^∞ S(τ) (τ² − b²)^{−1/2} dτ`.   (2.3)

(Justification of the interchange: `|∫_1^U sin(xu)(u²−1)^{−1/2} du| ≤ ln 2 + 4` uniformly in `U > 1`, `x > 0` — write
`(u²−1)^{−1/2} = 1/u + ρ(u)` with `ρ ≥ 0`, `∫_1^∞ ρ = ln 2`, and use `|Si| ≤ 2` — so dominated convergence applies on
`p ∈ [0,1]`; then `S(bu) = O(u^{−3})` makes the outer integral absolutely convergent.)

### Lemma 3 (`S` is positive outside one window) — [certified by interval arithmetic, `cert_S.py` (C1), plus an analytic tail]

`S(τ) > 0` for `τ ∈ (0, 9.6] ∪ [10.85, ∞)`.

*Proof.* (i) `0 < τ ≤ π`: `f > 0` on `(0,1)` and `sin(pτ) > 0` for `0 < p < 1`. (ii) `τ ∈ [3.14, 9.6] ∪ [10.85, 25]`: certificate
C1. It uses `|S′(τ)| = |∫ p f cos(pτ) dp| ≤ 283/120 =: L`: if `s ≤ S(c)` is an interval lower bound at a point `c`, then
`S > 0` on `(c − s/L, c + s/L)`; the script covers the two intervals by 65 + 2855 such cells (cell inclusion checked in
exact rational arithmetic). (iii) `τ ≥ 25`: bounding `|sin|, |cos| ≤ 1` in the closed form,
`τ⁶ S(τ) ≥ 80τ³ − 1602τ² − 6168τ − 35520 =: c(τ)`, and `c(25) = 59030 > 0`, `c′(25) = 63732 > 0`, `c″(τ) = 480τ − 3204 > 0` for
`τ ≥ 25`. ∎

(Numerically `S < 0` exactly on `(9.61049, 10.82823)`, with minimum `−0.025061` at `τ = 10.150`. Not used.)

### Corollary A (region (A)) — [proved: Lemmas 2, 3]

If `Mβ ≥ 10.85` then `F_M(β) > 0`, for every `M > 0`. Also `f̂(b) > 0` for `b ≥ 10.85`.

*Proof.* In (2.1) the integration variable satisfies `Mt ≥ Mβ ≥ 10.85`, so `K_{f,M}(t) = S(Mt) + R_M(t) > 0` by Lemma 3
and (2.2). The flat statement follows in the same way from (2.3). ∎

### Lemma 4 (small `M`: the whole kernel weight is decreasing) — [certified by exact rational arithmetic, `cert_smallM.py` (C2)]

For `0 < M ≤ M_1 := 23/20`, the function `h_M(p) := f(p) coth(πMp) = φ(p) · p coth(πMp)` is strictly decreasing on `(0,1)`.
Consequently `K_{f,M}(t) > 0` for all `t > 0` (Lemma 1), and `F_M(β) > 0` for all `β > 0` (2.1).

*Proof.* `d/dp log h_M = φ′/φ + πM ℓ(πMp)` with `ℓ(x) := d/dx log(x coth x) = 1/x − 2/sinh 2x`. Claim: `ℓ(x) < 2x/3` for `x > 0`.
With `y = 2x` this reads `(6 − y²) sinh y < 6y`; for `y² ≥ 6` it is trivial, and for `y² < 6` it follows from
`sinh y = y Σ_k y^{2k}/(2k+1)! < y Σ_k (y²/6)^k` since `(2k+1)! ≥ 6^k` with strict inequality for `k ≥ 2`. Hence
`d/dp log h_M < −(40 + 526p + 235p²)/((1−p)(65 + 155p + 47p²)) + (2π²M²/3) p`, which is `≤ 0` on `(0,1)` as soon as

`D_c(p) := 40 + 526p + 235p² − c · p(1−p)(65 + 155p + 47p²) ≥ 0` on `[0,1]` for some `c ≥ 2π²M²/3`.

C2 checks: `2π²M_1²/3 = 8.7017… < c := 871/100` (interval arithmetic); `D_c` has no real root in `[0,1]` (Sturm
sequence, exact) and `D_c(0) = 40 > 0`; and, independently, `D_c ≥ 5.11` on `[0,1]` by interval evaluation on 4096
subintervals. (`D_c` decreases in `c`, so `c = 871/100` covers all `M ≤ M_1`.) ∎

### Lemma 5 (an explicit lower bound for the kernel) — [proved analytically]

For every integer `m ≥ 0`, `M > 0`, `τ > 0`:

`K_{f,M}(τ/M) > K_low^{(m)}(M, τ) := ∫_0^1 f(p) [1 + 2 Σ_{k=1}^m e^{−2πkMp}] sin(pτ) dp`
`= S(τ) + 2 Σ_{k=1}^m Im Σ_{j=1}^{6} [f^{(j)}(0) − e^{−z_k} f^{(j)}(1)] z_k^{−j−1}`, `z_k := 2πkM − iτ`.

*Proof.* `coth(πx) = 1 + 2 Σ_{k=1}^m e^{−2πkx} + 2 e^{−2π(m+1)x}/(1 − e^{−2πx})`. The remainder weight is
`f(p) · 2e^{−2π(m+1)Mp}/(1 − e^{−2πMp}) = φ(p) · [2p/(e^{2πMp} − 1)] · e^{−2πmMp}`, a product of positive, strictly decreasing,
bounded functions, so its sine transform is positive by Lemma 1. The closed form is the integration-by-parts
formula of §1 with `z = z_k` (note `Im e^{−z_k p} = e^{−2πkMp} sin(pτ)`). ∎

### Lemma 6 (kernel positivity on the window for `1.15 ≤ M ≤ 6`) — [certified by interval arithmetic, `cert_K.py` (C3)]

`K_low^{(4)}(M, τ) > 0` for all `(M, τ) ∈ [1.15, 6] × [9.6, 10.85]`.

C3 evaluates the closed form of Lemma 5 (natural interval extension, complex arithmetic written out in real
and imaginary parts) on 2425 boxes of size `0.05 × 0.05` covering the rectangle; all pass; the smallest certified
lower bound is `0.01454` (the true minimum of `K_{f,M}` on the rectangle is about `0.0244`, at `M = 6`, `τ ≈ 10.12`).

### Proposition B (region (B)) — [proved: Lemmas 2–6]

For `0 < M ≤ 6`: `K_{f,M}(t) > 0` for all `t > 0`, hence `F_M(β) > 0` for all `β > 0`.

*Proof.* `M ≤ 1.15`: Lemma 4. `1.15 ≤ M ≤ 6`: if `Mt ∉ (9.6, 10.85)` then `K_{f,M}(t) = S(Mt) + R_M(t) > 0` by Lemma 3 and
(2.2); if `Mt ∈ [9.6, 10.85]` then `K_{f,M}(t) > K_low^{(4)}(M, Mt) > 0` by Lemmas 5 and 6. Then use (2.1). ∎

*Remark (this cannot be pushed to all `M`).* Numerically `min_t K_{f,M}(t)` is `0.0074` at `M = 7`, `0.0017` at `M = 7.5`,
`−0.0002` at `M = 7.7`, and tends to `min S = −0.0251` as `M → ∞`: pointwise positivity of the Mehler kernel **fails**
for `M ≳ 7.6` (as already observed in `EINSTEIN_PROGRAM.md` §3b for `M ≥ 20`). For `M ≥ 6` and `Mβ < 10.85` an integrated
argument is needed; this is Lemma 7 and C4.

### Lemma 7 (exact curvature expansion at fixed `b = Mβ`) — [proved analytically, given (R1)]

Put `a = (β+t)/2`, `c = (β−t)/2`. Then `cosh β − cosh t = 2 sinh a sinh c = ((β² − t²)/2) · E`, `E := (sinh a/a)(sinh c/c) ≥ 1`,
and (R1) with `t = β sin θ` becomes

`P_{−1/2+iν}(cosh β) = (2/π) ∫_0^{π/2} cos(νβ sin θ) E(β, β sin θ)^{−1/2} dθ`   (`E ≡ 1` gives `J_0(νβ)`).   (2.4)

Let `(sinh a/a)^{−1/2} = Σ_{k≥0} d_k a^{2k}` (`|a| < π`). Then:

1. `(−1)^k d_k > 0`, and `Σ_k |d_k| y^{2k} = (y/sin y)^{1/2}` for `|y| < π`. (Indeed `(y/sin y)^{1/2} = Π_{j≥1} (1 − y²/(π²j²))^{−1/2}` has
   non-negative Taylor coefficients, and `(sinh a/a)^{−1/2}` is the same series at `y = ia`.)
2. For `0 ≤ s ≤ 1` and `0 ≤ β < π`: `E(β, βs)^{−1/2} = Σ_{n≥0} β^{2n} q_n(s)`,
   `q_n(s) := Σ_{k+l=n} d_k d_l ((1+s)/2)^{2k} ((1−s)/2)^{2l}` (an even polynomial of degree `2n` with rational coefficients;
   `q_0 = 1`, `q_1 = −(1+s²)/24`, `q_2 = (7 + 22s² + 7s⁴)/5760`), with
   `|q_n(s)| ≤ μ_n := Σ_{k+l=n} |d_k||d_l| = [y^{2n}] (y/sin y)`, `μ_0 = 1`, `μ_1 = 1/6`, `μ_2 = 7/360`, and `μ_n ≤ 2/π^{2n}` for `n ≥ 1`.
   (From `y/sin y = 1 + 2y² Σ_{j≥1} (−1)^{j+1}/(j²π² − y²)`: `μ_n = 2η(2n)/π^{2n}` with `0 < η(2n) = Σ_j (−1)^{j+1} j^{−2n} < 1`.)
   In particular `Σ_n μ_n β^{2n} = β/sin β`.
3. Consequently, for `M > 0` and `0 < β < π`, with `b = Mβ`:

   `F_M(β) = 𝔽(Mβ, β)`, `𝔽(b, β) := Σ_{n≥0} β^{2n} G_n(b)`, `G_n(b) := (2/π) ∫_0^{π/2} f̃(b sin θ) q_n(sin θ) dθ`,   (2.5)

   where `G_0 = f̂`, `|G_n(b)| ≤ (2591/420) μ_n`, `|G_n′(b)| ≤ (283/120) μ_n`, and `G_n` is entire:
   `G_n(b) = Σ_{k≥0} g_{n,k} b^{2k}`, `g_{n,k} = (−1)^k (m_{2k}/(2k)!) Σ_{j=0}^n q_{n,j} W_{k+j}`, `m_{2k} := ∫_0^1 f p^{2k}` , `W_m := binom(2m,m)/4^m`,
   `q_n(s) = Σ_j q_{n,j} s^{2j}`, with `|g_{n,k}| ≤ (2591/420) μ_n/(2k)!`.

*Proof.* 1 and 2: absolute convergence of the double series `Σ d_k d_l a^{2k} c^{2l}` for `|a|, |c| < π` and regrouping by
`k + l = n`; `((1±s)/2)^2 ≤ 1` gives the bound on `q_n`. 3: insert 2 into (2.4), integrate against `f` (uniform
convergence for fixed `β < π`), use `f̃(b sin θ) = ∫ f(p) cos(pb sin θ) dp`; the bounds follow from `|q_n| ≤ μ_n`,
`|cos|, |sin| ≤ 1`; the Taylor coefficients from `(2/π) ∫_0^{π/2} cos(x sin θ) sin^{2j} θ dθ = Σ_k (−1)^k x^{2k} W_{k+j}/(2k)!`, and
`|(2/π) ∫ sin^{2k}θ q_n(sin θ) dθ| ≤ μ_n`, `m_{2k} ≤ m_0`. ∎

(For orientation: `G_1(b) = [65 + ∫_0^1 f″(p) J_0(pb) dp]/(24 b²) − f̂(b)/24`; the first term is the origin of `65 Q_1(β)/M²` with
`Q_1(0) = 1/24`. The expansion (2.5) was checked against the hypergeometric definition: agreement `10^{−16}`–`10^{−30}` for
`β ≤ 1`, `1.3·10^{−9}` at `M = 6`, `β = 1.8` with 15 terms, consistent with the tail bound below.)

### Lemma 8 (region (C)) — [certified by interval arithmetic, `cert_expansion.py` (C4)]

`𝔽(b, β) > 0` on `R := {0 ≤ b ≤ 10.85, 0 ≤ β ≤ b/6}`. Hence `F_M(β) > 0` whenever `M ≥ 6` and `Mβ ≤ 10.85`
(then `β ≤ 10.85/6 = 1.8084 < π` and `(Mβ, β) ∈ R`), and `f̂(b) = 𝔽(b, 0) > 0` for `0 ≤ b ≤ 10.85`.

*Certificate.* With `N = 14`, `K = 36`, `x = β²`: for a `b`-cell `[b_c − h, b_c + h]` and `x ∈ [x_0, x_1]`, Lemma 7 gives for every `b` in the cell

`𝔽(b, β) ≥ P_c(x) − (Σ_{n≤N} μ_n x_1^n) · [(283/120) h + (2591/420) T_K(b_c + h)] − (2591/420) · 2 (x_1/π²)^{N+1}/(1 − x_1/π²)`,

`P_c(x) := Σ_{n≤N} x^n Σ_{k≤K} g_{n,k} b_c^{2k}`, `T_K(b) := (b^{2K+2}/(2K+2)!)/(1 − b²/((2K+3)(2K+4)))`.

The three subtracted terms bound, respectively, the variation of `G_n` over the cell, the truncation of the
`b`-series, and the truncation `n > N` of the `β²`-series (at most `1.3·10^{−6}` on `R`). The coefficients `g_{n,k}`, `μ_n` are
computed exactly (rationals, `expansion.py`); `P_c` is enclosed on `[x_0, x_1]` by interval Horner evaluation and by the
centred form `P_c(x_m) + P_c′([x_0,x_1])·([x_0,x_1] − x_m)`; `x`-intervals and `b`-cells are bisected adaptively. Result: 541
`b`-cells, 681 boxes, all PASS. (The smallest certified box bound, `1.3·10^{−5}`, is a property of the cell sizes; the
true minimum of `𝔽` on `R` is `f̂(9.755) ≈ 0.02316`, attained on the flat slice.)

### Proof of Theorem 1 — [proved: Corollary A, Proposition B, Lemma 8]

Let `M > 0`, `β > 0`. If `M ≤ 6`: Proposition B. If `M ≥ 6` and `Mβ ≥ 10.85`: Corollary A. If `M ≥ 6` and `Mβ ≤ 10.85`:
Lemma 8. ∎

### Proof of Theorem 2 — [proved: Corollary A (flat part), Lemma 8 or C5a]

`b ≥ 10.85`: Corollary A, i.e. (2.3) and Lemma 3. `0 ≤ b ≤ 10.85`: the slice `β = 0` of Lemma 8; independently,
certificate C5a (`cert_flat.py`): Taylor series of `f̂` with exact rational coefficients and rigorous tail,
`|f̂′| ≤ (2/π)(283/120)`, 148 cells, all PASS. ∎

(Numerically: `min_{b ≤ 10.85} f̂ = 0.023158` at `b = 9.755`; `b³ f̂(b) ≥ 20.9` on the scanned range `5 ≤ b ≤ 40` (minimum at `b ≈ 9.59`); the
asymptotics `f̂(b) = 40/b³ + 2331/b⁵ + 10575/b⁷ + O(b^{−9/2})` oscillatory, of amplitude `≤ 1602·√(2/π)·b^{−9/2}(1+o(1))`, is
confirmed, cf. (I5) below. None of this is used.)

### Lemma 9 and Theorem 1′ (second proof for `M ≥ 45`) — [Lemma 9 proved analytically given (R1); C5b certified by interval arithmetic, `cert_flat.py`]

**Lemma 9 (= ingredient (I1)).** For real `ν` and `β > 0`:
`|P_{−1/2+iν}(cosh β) − J_0(νβ)| ≤ 1 − P_{−1/2}(cosh β) ≤ sinh²(β/2)/4`; hence `|F_M(β) − f̂(Mβ)| ≤ (2591/420) sinh²(β/2)/4`.

*Proof.* By (2.4), `P − J_0 = (2/π) ∫ cos(νβ sin θ)(E^{−1/2} − 1) dθ` and `0 ≤ 1 − E^{−1/2}`, so the modulus is at most the value
at `ν = 0`, which is `1 − P_{−1/2}(cosh β)`. By Euler's integral, `P_{−1/2}(cosh β) = (2/π) ∫_0^{π/2} (1 + z sin²φ)^{−1/2} dφ ≥ 1 − z/4`,
`z = sinh²(β/2)`. ∎

**Theorem 1′.** For `M ≥ 45`: if `Mβ ≥ 10.85`, Corollary A; if `b = Mβ ≤ 10.85`, then
`F_M(β) ≥ f̂(b) − (2591/1680) sinh²(b/(2M)) ≥ f̂(b) − (2591/1680) sinh²(b/90) > 0` by certificate C5b (361 cells, all PASS). ∎

This proof uses only (R1), (R2), C1 and the one-dimensional certificate C5b; it does not use C2–C4. Without
(2.5) the method stops near `M ≈ 40` (`max_{b ≤ 10.85} [(∫f) b²/(16 f̂(b))]^{1/2} ≈ 40.1`, at `b ≈ 9.9`, numerically).

## 3. Certificates (reproduce with `bash ads4/code/unipos/run_all.sh`)

| id | script | statement certified | method | boxes | smallest certified bound |
|---|---|---|---|---|---|
| C1 | `cert_S.py` | `S > 0` on `[3.14, 9.6] ∪ [10.85, 25]`; `c(25), c′(25), c″(25) > 0` | point enclosures + Lipschitz constant `283/120`; exact integers | 65 + 2855 cells | `1.09·10^{−3}`, `9.4·10^{−4}` |
| C2 | `cert_smallM.py` | `D_{871/100} > 0` on `[0,1]`, `2π²(23/20)²/3 < 871/100`, `φ′` identity | Sturm (exact) + interval check | 4096 subintervals | `D ≥ 5.11` |
| C3 | `cert_K.py` | `K_low^{(4)} > 0` on `[1.15, 6] × [9.6, 10.85]` | interval evaluation of the closed form | 2425 | `0.01454` |
| C4 | `cert_expansion.py` | `𝔽 > 0` on `{b ≤ 10.85, β ≤ b/6}` | exact rational coefficients, interval Horner / centred form, explicit tails (`N = 14`, `K = 36`) | 541 `b`-cells, 681 boxes | `1.3·10^{−5}` |
| C5a | `cert_flat.py` | `f̂ > 0` on `[0, 10.85]` | Taylor series + Lipschitz constant `(2/π)(283/120)` | 148 cells | `2.6·10^{−5}` |
| C5b | `cert_flat.py` | `f̂(b) > (2591/1680) sinh²(b/90)` on `[0, 10.85]` | same | 361 cells | `9.5·10^{−6}` |

Each script prints one `PASS`/`FAIL` line per box or cell and a final `VERDICT` line; `run_all.sh` stores the full
output in `code/unipos/logs/` and prints the summaries and the overall verdict (all PASS; about 10 s on one
machine; `--no-sanity` skips the floating-point checks). Constants are in `common.py` (`TAU1 = 9.6`, `TAU2 = 10.85`,
`M1 = 23/20`, `MK = 6`). The margins are comfortable: the true minima are `S(9.6) = 1.07·10^{−3}`, `S(10.85) = 1.38·10^{−3}`,
`K ≥ 0.024`, `𝔽 ≥ 0.023`.

## 4. Status of every statement

| statement | status |
|---|---|
| (R1), (R2), (R3) classical representations | quoted from the literature; numerically re-verified (not re-proved) |
| closed forms of `S`, `f̃`; moments; derivatives of `f`; `φ′` | [proved analytically]; exact checks in `common.py`, `cert_smallM.py` |
| Lemma 1 (decreasing ⇒ positive sine transform) | [proved analytically] |
| Lemma 2 ((2.1), `R_M > 0`), (2.3) | [proved analytically] given (R2), (R3) |
| Lemma 3 (`S > 0` off `(9.6, 10.85)`) | [certified by interval arithmetic, `cert_S.py`] + analytic for `τ ≤ π`, `τ ≥ 25` |
| Corollary A (`Mβ ≥ 10.85`, all `M`; `f̂(b) > 0`, `b ≥ 10.85`) | [proved], depends on Lemma 3 |
| Lemma 4 (`M ≤ 1.15`) | [proved analytically] + [exact rational certificate, `cert_smallM.py`] |
| Lemma 5 (`K ≥ K_low^{(m)}`) | [proved analytically] |
| Lemma 6 (`K_low^{(4)} > 0` on the window, `1.15 ≤ M ≤ 6`) | [certified by interval arithmetic, `cert_K.py`] |
| Proposition B (`M ≤ 6`, all `β`) | [proved], depends on Lemmas 3–6 |
| Lemma 7 (expansion (2.5) and its bounds) | [proved analytically] given (R1) |
| Lemma 8 (`M ≥ 6`, `Mβ ≤ 10.85`) | [certified by interval arithmetic, `cert_expansion.py`] |
| **Theorem 1** (all `M > 0`, `β > 0`) | [proved], depends on all of the above |
| **Theorem 2** (`f̂ > 0` on `[0,∞)`) | [proved]: Corollary A + [certified, `cert_expansion.py` slice `β = 0`, and independently `cert_flat.py` C5a] |
| Lemma 9 (= (I1)) | [proved analytically] given (R1) |
| Theorem 1′ (`M ≥ 45`, second proof) | [proved]: Corollary A + [certified, `cert_flat.py` C5b] |
| failure of pointwise kernel positivity for `M ≳ 7.6`; location of the zeros of `S`; all tables in §6 | [numerical only — NOT part of the theorem] |
| correctness of `mpmath.iv` | assumed; `iv_selftest.py` is a sanity check only |

## 5. The ingredients (I1)–(I5) of the task, and what was tried and dropped

* (I1) correct; proved as Lemma 9; used only in Theorem 1′.
* (I2) correct, including the coefficients `−65/τ² − 1338/τ⁴ + 1680/τ⁶` and the oscillatory terms starting at `τ^{−4}`
  (closed form in §1). `∫_0^∞ f̃ = (π/2) f(0) = 0` (`Ψ(X) = ∫_0^X f̃ ~ 65/X`). Used in the form (2.4)–(2.5).
* (I3) `φ′ = −(1−p)²(40 + 526p + 235p²)` correct (exact check in C2); `Ψ > 0` follows from Lemma 1. The integration
  by parts of (I3) is not used. The same monotonicity idea, applied to the full weight `φ(p)·p coth(πMp)`, is
  Lemma 4, and applied to the remainders `φ(p)·2p e^{−2πmMp}/(e^{2πMp} − 1)` is Lemma 2/5.
* (I4) not used. (Given Lemma Q the signs are as stated: `a_2 Q_2`, `a_4 Q_4`, `a_6 Q_6 > 0`. In the present language
  `Σ a_n Q_n/M^{n+1} − T_M` is the exterior integral of `S(Mt) + R_M(t)`; the proof needs only `R_M > 0`, not its
  expansion.)
* (I5) correct: the non-oscillatory part of `S` is `80/τ³ + 6216/τ⁵ + 33840/τ⁷`, and (2.3) with
  `(2/π) ∫_1^∞ u^{−k}(u²−1)^{−1/2} du = 1/2, 3/8, 5/16` (`k = 3, 5, 7`) gives `40/b³ + 2331/b⁵ + 10575/b⁷`; `f‴(1) = −1602`.
  Confirmed numerically at `b = 40, 80, 160`. Not needed for Theorem 2, because `S > 0` on `[10.85, ∞)` is stronger.
* Dropped: a van der Corput/Fresnel bound for the oscillatory part of `∫_β^∞ S(Mt)(cosh t − cosh β)^{−1/2} dt` (with the
  constant `sup_x |∫_0^x e^{is} s^{−1/2} ds| ≈ 2.379`). It works, but only gives positivity for `Mβ ≳ 14–15`, and is
  superseded by Lemma 3. Also dropped: a quantitative lower bound `R_M ≥ c_M · 65 K_1(t)/M²` (not needed).
* Failed, as expected: positivity of the Mehler kernel `K_{f,M}` for all `M` (false for `M ≳ 7.6`, §2 remark).
* A counterexample search found none: on all grids `F_M > 0` (§6).
* Cross-check of the certificates themselves (`sanity_numeric.py`, part 6): on randomly chosen boxes of the logs
  of C1, C3, C4 the certified lower bounds lie below brute-force floating-point values, as they must.

## 6. Numerical illustration (floating point; NOT part of the proof) — `sanity_numeric.py`

Minimum over the grid `b = Mβ ∈ (0, 200]` (1150 points) of `F_M(β) e^{β/2} max(1, b)²`, and where it is attained:

| `M` | 0.02 | 0.1 | 0.3 | 1 | 1.15 | 2 | 4 | 6 | 10 | 15 | 20 | 30 | 45 | 80 | 200 | 1000 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| min | 9.57 | 6.80 | 6.38 | 6.23 | 6.22 | 6.20 | 6.18 | 6.18 | 6.18 | 4.25 | 3.44 | 2.83 | 2.52 | 1.82 | 1.01 | 0.335 |
| at `b` | 0.02 | 0.02 | 0.02 | 0.02 | 0.02 | 0.02 | 0.02 | 0.02 | 0.02 | 9.59 | 9.59 | 9.63 | 16.2 | 28.9 | 60.3 | 186 |

(For large `M` the scan covers only `β ≤ 200/M`; larger `β` is Corollary A. The decrease of the normalised minimum
with `M` reflects the normalisation — in the flat regime `F ≈ 40/b³` — not a loss of positivity.) So nothing special
happens at small `M`: for `M → 0` at fixed `b = Mβ`, `F_M(β) ≈ (√2/π)(sinh β)^{−1/2} Ψ(b)/M` with `Ψ(x) = ∫_0^1 φ(p) sin(px) dp > 0`
(checked: `473.0` against `474.0` at `M = 0.02`, `b = 5`, in units of `e^{−β/2}`), which is the limiting form of Lemma 4.

Pointwise kernel: `min_τ K_{f,M}(τ/M) = 0.140, 0.120, 0.0553, 0.0244, 0.0074, 0.0017, −0.0002, −0.0027, −0.0132, −0.0235, −0.0250`
for `M = 1.15, 3, 5, 6, 7, 7.5, 7.7, 8, 10, 20, 80` (the first two at the lower end of the `τ` grid, the others at `τ ≈ 10.13`).

Regions: for `M ≥ 6` and `b ≤ 10.85` the minimum of `F_M` is `0.071` (`M = 6`), `0.054` (8), `0.038` (12), `0.029` (20),
`0.024` (50), `0.0232` (`M → ∞`), near `b ≈ 9.8`; and `F_M(β) − f̂(Mβ)` is negative around `b ≈ 3.9` (`−0.066` at `M = 6`), so
`F_M ≥ f̂(Mβ)` is **not** true and could not have been used.

## 7. What is and is not established; caveats

* Established (subject to (R1)–(R3) and to the correctness of `mpmath.iv`): Theorems 1, 2, 1′ as stated. Nothing
  remains unproved for the statement posed; there is no excluded range of `M` or `β`.
* The proof is specific to this `f`. The method (sign of the sine transform beyond the last negative window,
  monotone weights, exact `β²`-expansion at fixed `Mβ`) is general, but every certificate would have to be rerun,
  and Lemma 3-type statements can fail for other weights (e.g. if `S` has negative windows at large `τ`).
* Not attempted: a quantitative lower bound for `F_M` (e.g. in terms of `f̂(Mβ)` and `65 Q_1(β)/M²`); the scalar
  weight of arXiv:2106.10274; any statement about the full functional, the `1/m` corrections, or light states.
* No independent referee has checked this document. The points a referee should attack first: the two
  regrouping steps of Lemma 7 (absolute convergence; `|q_n| ≤ μ_n`), the error budget of C4, and the reduction
  in Lemma 4 (`ℓ(x) < 2x/3`).

## 8. Files

`ads4/code/unipos/`: `common.py` (exact data, interval helpers, `S`, `K_low`), `expansion.py` (exact coefficients `d_k`,
`q_n`, `μ_n`, `g_{n,k}`), `cert_S.py` (C1), `cert_smallM.py` (C2), `cert_K.py` (C3), `cert_expansion.py` (C4), `cert_flat.py`
(C5a, C5b), `iv_selftest.py`, `sanity_numeric.py` (non-rigorous), `run_all.sh` (driver; writes `logs/`).
