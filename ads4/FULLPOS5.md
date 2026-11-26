# FULLPOS5: the high-twist certificate at large impact parameter, β > 8.4 (L ≳ 7)

*26 Nov 2026. Research agent. Unreviewed, not committed. No existing file was edited; nothing under `prl/` was touched.*
*New code: `ads4/code/fullpos3/fp3_ht_large.py` (certificate), `ads4/code/fullpos3/run_large_scan.py` (parallel launcher).
Outputs go to `fullpos3/cache/` (git-ignored). `fp3_ht.py`, `fp3_htdata.py` and `run_laptop_scan.py` are unchanged, so the laptop job
(β ≤ 8.4) is unaffected.*
*Compute here: ≤ 2 processes, `ulimit -v 4000000`, `OMP_NUM_THREADS=1`.*

Grades: **[certified]** Arb balls with every remainder proved; **[P]** proved on paper here; **[DV]** numerical.
Notation as in FULLPOS3 §4: c = tanh²(β/2), d = 1 − c = e^{−L} (so L = 2 log cosh(β/2) ≈ β − 1.386), h = 1/(N(N+s)), a = twist.

---

## 0. Verdict

1. **Large-β failure explained and fixed.** `stadium geometry failed` came from an absolute floor (10⁻¹²) on the arc-parameter width,
   which is larger than the whole stadium once d ≲ 10⁻¹¹. Four other pieces of `fp3_ht.py` also stop working at large β; none of them
   could produce a wrong PASS (each makes the run fail or loosens the bound), but each must be replaced (§1).
2. **Three sharper rigorous bounds** (§1.4–1.6): Lemma B′ (oscillation-aware bound on the conical function for the realised-weight
   residual), Lemma A′ (the coth correction integrated by parts once), and a frequency split of the residual "rest" term. Together they cut
   the error at β = 12, a₀ = 4500 from 0.27·F_cert to 0.013·F_cert.
3. **Probes [certified], cells of width 0.1, both channels:**

   | β-cell | L | a*(b = 0) | a*(b = 4) | err/F_cert at a* |
   |---|---|---|---|---|
   | [9, 9.1] | 7.6 | 3200 | 3200 | 0.31 |
   | [12, 12.1] | 10.6 | 3200 | 3200 | 0.45 |
   | [16, 16.1] | 14.6 | 3200 | 3200 | 0.62 |
   | [20, 20.1] | 18.6 | 3200 | 3200 | 0.79 |
   | [30, 30.1] | 28.6 | 3600 (3400 passes, b = 0 run: 0.31) | 3600 | 0.098 |
   | [45, 45.1] | 43.6 | 3600 | 3600 | 0.149 |
   | [63.9, 64] | 62.6 | 3600 (3400 passes, b = 0 run: 0.57) | 3600 | 0.185 |

   (The probes used the old twist ladder …, 2800, 3200, 3600, …; single runs show 3400 passing at β = 30 and 64.) **a* ≈ 3200–3600 on
   all of [9, 64], channel-independent to 7 digits.** a* grows slowly (§2): at fixed a the error grows like β relative to F_cert, while it
   falls steeply in a (≈ a⁻¹²). The new bounds also lower a* in the old range: the cell [8.24, 8.4] gives a* = 3400 (fp3_ht: 4000).
4. **β ≥ 64: not covered by this certificate, and it cannot be (§3).** Every term has an explicit β → ∞ behaviour. The low-frequency
   residual mass (EPSW cell ν ∈ [0, 0.31], ∫|r| ≈ 8.9·10⁻⁷) is known only as an absolute integral, and its bound grows like βe^{−β/2},
   against F_cert ∝ e^{−β/2}/β. So |G^r|/F_cert ≈ 8·10⁻⁶β² for every a, which reaches 1 near β ≈ 350. Cells can be pushed past 64
   (the launcher takes `--beyond B`, with LRINF's tail floor) up to roughly β ≈ 200, but not to ∞. The β → ∞ end needs the contour
   route (CONTOUR_BOUNDS, ℰ₂ < 0.01𝓜 for L ≥ 30.16, all a ≥ 998); its connection terms j ≥ 1 and k = 3 are still unchecked.
5. **Cell plan and launcher (§4).** 556 cells of width 0.1 on [8.4, 64]; `python3 run_large_scan.py` on the laptop. Estimated CPU:
   ≈ 90–110 CPU-hours, i.e. about 4–5 h on 22 workers. The part below the contour threshold (β < 31.6) is ≈ 25–30 CPU-h, about
   1.3 h, and runs first.

---

## 1. What was changed (`fp3_ht_large.py`) [P]

`fp3_ht_large.py` imports `fp3_ht` and replaces functions **inside its own process only** (`fp3_ht.circle_sups`, `jval`,
`norm_major_float`, `far_remainder`, `run`, `A0LIST`). Nothing on disk in `fp3_ht.py` / `fp3_htdata.py` changes, so the β ≤ 8.4
results and the running laptop job keep their meaning exactly. Every replacement is an upper bound in its own right; most are the same
formula evaluated more carefully.

### 1.1 Why the old code fails at large β (none of these can produce a wrong PASS)
| piece of `fp3_ht.py` | large-β behaviour | effect |
|---|---|---|
| `circle_sups`: stop subdividing when the arc-parameter width < 10⁻¹² (absolute) | the straight leg of the stadium has length D_l ≈ d(1 − e^{−Δβ}) ~ 10⁻¹¹ at β = 25 | `RuntimeError` (the reported failure) |
| `ctx.prec = 256`, z = c e^t, then 1 − z | 1 − z ~ d loses log₂(1/d) ≈ 1.44β bits | wide balls (FAIL, not wrong) |
| `jval`: J-tables on a grid up to β = 9.9 | for β > 9.9 it returns J(9.9) | valid but loses e^{−(β−9.9)/2} |
| `far_remainder`: candidate ρ ∈ [10⁻²⁰, 1] (absolute) | admissible ρ needs ρu_k < 1, u_k up to n² ~ 4a²/d² | "no admissible rho" for β ≳ 20 → run reports FAIL |
| `far_remainder`: iteration cap k ≤ 5·10⁸ | the sum runs to k ~ 10³/d | FAIL |
| `norm_major_float` (float, η choice only) | float c rounds to 1 for β > 37; log(1 − c^l) = −∞ | garbage η → huge/inf Arb bound → FAIL |
| G^r and residual components: Lemma B, P_{−½}(cosh β)·|ξ|/β | grows like βe^{−β/2}, against F_cert ∝ e^{−β/2}/β | valid; ratio grows like β² |
| coth correction: C_cor(η)·J_C | ∝ e^{−β/2}, no 1/β | valid; ratio grows like β |
| cells of fixed width Δβ | one R for the whole cell: the m-th order loses e^{3mΔβ} | width 0.5 fails at β = 30 (err/F = 4·10¹⁴ at a = 4500), 0.2 gives 0.048, 0.1 gives 0.020 |

Every other ingredient is scale-invariant in t/d: the radii (R ≤ 0.92 d), the Im ξ tolerance, h_max (≈ d²/(4a²)), the exact P_m, the
majorant, and the dressing operator 𝒟_b. Indeed hd·ϑ³ ~ h·d·d⁻³ ~ 1/(4a²) uniformly in β.

### 1.2 Geometry and precision
* ctx.prec = 256 + ⌈2β_b/ln 2⌉ bits (440 bits at β = 64).
* 1 − z = (1 − c_l) − c_l·expm1(t) (acb), X = (2 − (1 − z))/(1 − z), ξ = acosh X. Same quantities as before, without cancellation.
* The subdivision floor is relative: width < 10⁻¹²·(length of that piece of the stadium boundary).

### 1.3 J-tables [P]
For β ≥ 4: J_X(β) ≤ e^{−(β−g)/2}J_X(g), where g = the grid point (step 0.02) below β.
**Proof.** sinh(x + δ) ≥ e^δ sinh x for x, δ ≥ 0 (since cosh x ≥ sinh x), so the Mehler kernel K_β(u) = (2 sinh(β+u/2) sinh(u/2))^{−½}
≤ e^{−(β−g)/2}K_g(u). A_X(β + u) ≤ A_X(g + u) because A_X is a positive combination of negative powers. ∎
The grid values are `fp3_htdata.jtable` (unchanged). Below β = 4 the old tables are used.

### 1.4 Lemma B′: an oscillation-aware bound for the residual [P, DV-checked]
For ν > 0, Re ξ = β > 0, |Im ξ| = η < π/2 and any δ > 0:

  |P_{−½+iν}(cosh ξ)| ≤ (√2/π) coth(πν) γ_β e^{−β/2} e^{νη} Q(ν, δ; β),
  γ_β = (1 − e^{−2β})^{−½},
  Q = 2√2·arcsin√(1 − e^{−δ}) + ν⁻¹[2s(δ) + (√2/2)(coth β − 1)·arcsin(e^{−δ/2})],  s(δ) = √2e^{−δ/2}/√(1 − e^{−δ}).

**Proof.** Start from Lemma A (FULLPOS3 §4.2): P = (√2/π)coth(πν)∫₀^∞ sin(ν(ξ+u))k(u)du, k = (2 sinh(ξ+u/2) sinh(u/2))^{−½}.
1. |sinh w| ≥ sinh Re w gives |k(u)| ≤ K_β(u). Then 2 sinh(β+u/2) ≥ e^{β+u/2}(1 − e^{−2β}) and 2 sinh(u/2) = e^{u/2}(1 − e^{−u}) give
   K_β(u) ≤ γ_βe^{−β/2}s(u), with s(u) = √2e^{−u/2}(1 − e^{−u})^{−½}.
2. On [0, δ]: |sin| ≤ e^{νη}, and ∫₀^δ s = 2√2 arcsin√(1 − e^{−δ}) (substitute x = e^{−u}: ∫x^{−½}(1−x)^{−½} = 2 arcsin√x).
3. On [δ, ∞): integrate by parts once, with −cos(ν(ξ+u))/ν. The boundary term is ≤ e^{νη}|k(δ)|/ν. The integral is ≤ (e^{νη}/ν)∫_δ^∞|k′|,
   with k′ = −(k/4)[coth(ξ+u/2) + coth(u/2)] and |coth w| ≤ coth Re w ≤ coth β (since |cosh w|² ≤ cosh²Re w and |sinh w|² ≥ sinh²Re w).
4. Exact integrals: ∫_δ^∞ s = 2√2 arcsin(e^{−δ/2}), and ∫_δ^∞ s·coth(u/2) = √2[4e^{−δ/2}/√(1−e^{−δ}) − 2 arcsin(e^{−δ/2})]. The latter
   uses (1 + x) = 2 − (1 − x) and d/dx[2√x/√(1−x)] = x^{−½}(1−x)^{−3/2}. Collecting terms gives Q. ∎

* It is used per frequency cell [ν_lo, ν_hi] of the residual: coth(πν) ≤ coth(πν_lo), 1/ν ≤ 1/ν_lo, and the cell weight already
  carries e^{ην_hi}. The code takes **min(Lemma B, Lemma B′)** per cell and per stadium point. δ is minimised over 81 candidates.
  The factor decreases in β, so it is evaluated at β = 7 for every β ≥ 7.
* Size: the factor K(ν) = e^{β/2}·bound is 5.92 (ν = 0.112), 2.01 (ν = 1), 0.36 (ν = 48.75), 0.081 (ν = 998.75), for every β ≥ 7.
  Lemma B gives e^{β/2}P_{−½}(cosh β) = 6.2 (β = 8.4), 11.1 (16), 20.0 (30), 41.6 (64). So at ν ~ 10³ the gain is ×80–500.
* [DV] mpmath (`legenp`, type 3) at β ∈ {0.5, 2, 8.4, 12, 20, 30}, η ∈ {0, …, 1.5}, ν ∈ {0.06, …, 998.75}: the ratio of the true |P| to the
  bound is ≤ 0.98.

### 1.5 Lemma A′: the coth correction, integrated by parts [P, DV-checked]
g(p) = ψ(p)(coth(πΛp) − 1) = 2ψ(p)/(e^{2πΛp} − 1), with g(1) = 0 and g(0) = c₁/(πΛ) (ψ = Σc_np^n). For |Im t| ≤ η < π/2:

  |∫₀¹ g(p) sin(Λpt) dp| ≤ C₁(η)/(Λ|t|),  C₁ = |c₁|/(πΛ) + Σ_n 2|c_n|(2πΛ)^{−n}I_n(κ),  κ = η/(2π),
  I₁ ≤ 1 + (κ/2)[(1−κ)⁻² + (1−κ)⁻¹],  I_n ≤ 2(2(n−1)/e)^{n−1}cosh(κn) + (κ/2)n![(1−κ)^{−n−1} + (1−κ)^{−n}/n]  (n ≥ 2).

**Proof.**
1. One integration by parts gives g(0)/(Λt) + (Λt)⁻¹∫g′cos(Λpt), with |cos(Λpt)| ≤ cosh(Λpη).
2. In x = 2πΛp, take termwise |d/dx(x^n/(e^x−1))|. For n = 1, x/(e^x−1) decreases from 1. Integrate by parts against cosh(κx), and use
   ∫x sinh(κx)/(e^x−1) ≤ ½Σ_m(m−κ)⁻².
3. For n ≥ 2, x^n/(e^x−1) is unimodal with its peak at x* < n and peak value ≤ max x^{n−1}e^{−x/2} (since x/(e^x−1) ≤ e^{−x/2}).
   Split at x* and use ∫x^n e^{κx}/(e^x−1) = n!Σ(m−κ)^{−n−1}. ∎

With |t| ≥ Re ξ + u ≥ β, the coth part of the sup bound becomes min(C_cor(η)J_C, C₁(η)J_C/(Λβ)). C₁(0) = 0.120 (ψ₂) and 0.032 (ψ₃).
So C₁/(Λβ) beats C_cor = 1.9·10⁻⁵ for β > 6.3. At β = 64 the gain is ×10.
[DV] mpmath quadrature at Re t ∈ {0.3, 2, 10, 50}, η ∈ {0, 0.7, 1.5}: |integral|/bound ≤ 0.50.

### 1.6 The residual "rest" term, split by frequency support [P]
`fp3_htdata.rint_components` lumps the last three residual pieces into one 'rest' term (Lemma B only). They are now separate, with the
same formulas, so they sum to the old 'rest' exactly (checked):
* Gaussian tail 2.1·tail, supported on ν ≥ 1059;
* θ-part 1, 10¹⁰e^{−49}e^{η(Λ+30)}, on ν ≤ Λ + 30 (Lemma B only);
* θ-part 2, 2.1·3.2·10⁻¹⁰e^{1059η}, on Λ < ν ≤ 1059 (Lemma B′ with ν_lo = Λ).

Before the split, θ-part 2 under Lemma B was the binding local component at β = 12 (0.26·F_cert at a₀ = 4500).

### 1.7 Far remainder and majorant
* Candidate ρ grid: h_max·10^{i/4}, i = −24…56. This is the same 81-point span around h_max that the absolute grid gives at β ≈ 8; the
  lemma allows any admissible ρ.
* Iteration cap: 10⁴/d.
* The η-choice float majorant is now written in d: log((c^{−l} − 1)/l) = x + log(−expm1(−x)) − log l, with x = −l·log1p(−d). This form is
  stable for x → 0 (large β) and for x → ∞ (small β; a first version overflowed there, which only made low-β runs fail). The bound
  itself is still the Arb `fp3_ht.norm_major` at c₋, so the float only picks η. It agrees with fp3_ht's float to all digits at
  β ≈ 0.4 and 8.
* `run_large` is `fp3_ht.run` line for line except: the precision; G^r via Lemma B′; the two replaced functions; F_cert beyond 64; and a
  per-component budget `loc_comp` in the output. A `diff` is 15 lines.
* Result: far/F_cert ≤ 1.1·10⁻⁴ in every probe (the M-selection rule asks < 10⁻³). At β = 63.9 it is 5·10⁻²²⁸ with M = 450,
  n ≈ 10³¹, h_max ≈ 10⁻⁶², and the k-sum switches to its geometric tail at k ≈ 7.6·10³⁰.

### 1.8 Floor beyond β = 64 (only for `--beyond`)
`fcert_tail(β)` = (1 − 2.1·10⁻⁷)·A₁,min·(√2/π)K₀(β)ℓ/(βΛ²). Here A₁,min = 2ψ₂′(0) − ψ₃′(0) (exact), K₀ = (cosh β − 1)^{−½}, and
ℓ = 64L(64) is recomputed in Arb (= 0.8740…). This is LRINF_TIGHT §1.5 with the flat term dropped and r ≤ 2.03·10⁻⁷ taken as 2.1·10⁻⁷.
It is decreasing in β, so the cell value is taken at β_b. At β = 64 it gives 3.61·10⁻²⁰, against LRINF's cell value 3.60·10⁻²⁰.

### 1.9 Twist ladder
The bisection list gets intermediate values: 998, 1100, 1250, 1400, 1500, …, 2000, 2125, 2250, 2375, 2500, 2650, 2800, 2900, 3000, 3100,
3200, 3400, 3600, 3800, 4000, 4250, 4500, 4750, 5000, …. The probes in §2 used the old list (…, 2800, 3200, 3600, …).

---

## 2. Probe results [certified]

`fp3_ht_large.py scan` (bisection, both channels, cell width 0.1, exact P_m to m = 70, 81 radii). Results are in
`cache/large/probe{A,B}.jsonl`. They were moved out of the launcher's `cache/large_*.jsonl` pattern so that the production run recomputes
these cells with the final ladder. They came from a revision that differed only in the float η-choice and the ladder, and remain
valid certificates.

| β-cell | channel | a* | err/F | G^r/F | Σ local/F | far/F | M | top local components (/F) | cell time |
|---|---|---|---|---|---|---|---|---|---|
| [9, 9.1] | 0, 4 | 3200 | 0.307 | 0.0027 | 0.305 | 2.8·10⁻⁹ | 300 | θ₂ 0.134, ν∈[949,999] 0.066, [999,1049] 0.039 | 248 s |
| [12, 12.1] | 0, 4 | 3200 | 0.453 | 0.0043 | 0.448 | 2.1·10⁻⁷ | 300 | θ₂ 0.192, 0.094, 0.056 | 284 s |
| [16, 16.1] | 0, 4 | 3200 | 0.619 | 0.0063 | 0.613 | 1.1·10⁻⁴ | 300 | θ₂ 0.255, 0.125, θ₁ 0.080 | 334 s |
| [20, 20.1] | 0, 4 | 3200 | 0.788 | 0.0088 | 0.779 | 2·10⁻²⁰⁹ | 450 | θ₂ 0.315, 0.154, θ₁ 0.121 | 413 s |
| [30, 30.1] | 0, 4 | 3600 | 0.098 | 0.0186 | 0.079 | 8·10⁻²⁶ | 300 | θ₂ 0.025, θ₁ 0.018, 0.014 | 422 s |
| [30, 30.1] | 0 (single run) | **3400** | 0.314 | 0.0186 | — | 3·10⁻¹² | 300 | | 121 s |
| [45, 45.1] | 0, 4 | 3600 | 0.149 | 0.0315 | 0.118 | 1·10⁻¹⁵ | 300 | θ₁ 0.036, θ₂ 0.034, 0.018 | 689 s |
| [63.9, 64] | 0, 4 | 3600 | 0.185 | 0.0435 | 0.141 | 5·10⁻²²⁸ | 450 | θ₁ 0.055, θ₂ 0.036, 0.019 | 1012 s |
| [63.9, 64] | 0 (single run) | **3400** | 0.566 | 0.0435 | — | 7·10⁻²⁰⁵ | 450 | | 270 s |
| [8.24, 8.4] (old range, comparison) | 0, 4 | 3400 (fp3_ht: 4000) | 0.786 | | | | | | 214 s |

θ₂ = θ-part 2 of the residual (ν ∈ (Λ, 1059], Lemma B′); θ₁ = θ-part 1 (all ν, Lemma B only); the other entries are EPSW residual cells
(Lemma B′). The two channels agree to 7–8 digits everywhere (b only enters through s = 2b+1 in P_m and 𝒟_b, both subleading here). Cell
times are wall-clock on a shared 4-core machine (load ≈ 5) with the old ladder; a "single run" is one a₀ in one channel.

**How a* depends on β.** At a₀ = 3200 the error grows linearly in β, err/F ≈ 0.04β: 0.31, 0.45, 0.62, 0.79 at β = 9, 12, 16, 20.
* The binding components are the high-frequency residual pieces: θ-part 2 (ν ∈ (Λ, 1059]) and the EPSW cells at ν ≈ 900–1060, all under
  Lemma B′. Their kernel is β-independent times e^{−β/2}, while F_cert ∝ e^{−β/2}/β, hence the factor β.
* G^r/F_cert is 0.003–0.044 on [9, 64] (§3).
* The local terms fall steeply in a: at β = 12 the total error is 0.45F at a₀ = 3200 and 0.013F at a₀ = 4500, i.e. roughly ∝ a^{−12}.
  So a*(β) ∝ β^{1/12} in this regime: 3200 for β ≤ 20, 3400 at β = 30 and at β = 64 (err/F = 0.31 and 0.57).
* Above β ≈ 40 the Lemma-B-only θ-part 1 (r25 below) overtakes θ-part 2, and G^r rises to 0.044F at β = 64.

---

## 3. β → ∞, and the status of β ≥ 64

### 3.1 Behaviour of each term [P]
Scale t = dτ. Then 1 − z = d(1 − τ) + O(d²), ξ = β − log(1 − τ) + O(d), and h ≈ d²/(4a²). After multiplying by w³ = d³:

| term | β → ∞ behaviour (×e^{β/2}) | relative to F_cert ∝ 1/β |
|---|---|---|
| F_cert (LRINF; F_cert·βe^{β/2} = 1.3–1.4·10⁻⁴ over [8.4, 64], 1.83·10⁻⁴·(1 − 2·10⁻⁷) analytically for β ≥ 64) | ∝ 1/β | 1 |
| ψ p = 0 part A_s | ∝ (Λβ)⁻³ | → 0 |
| ψ p = 1 endpoint, e^{Λη}J_e | ∝ (Λβ)⁻⁶ | → 0 |
| coth correction: old C_cor·J_C / new Lemma A′ | const / ∝ 1/β | ∝ β / **O(1)** |
| residual, EPSW cells ν_lo ≥ 0.05 and θ-part 2, tail (Lemma B′) | const | **∝ β** |
| residual, cells [0, 0.3], [0.1, 0.312] and θ-part 1 (Lemma B only) | ∝ β | **∝ β²** |
| local orders h^m‖P_m‖_R | scale-invariant (O(d) corrections) | factor of the components above |
| far remainder | the crude envelope \|g_k\| ≤ B_k = O(e^{1573}k⁵) gives a d^{−3}-growing prefactor, beaten by (h/ρ)^M | small at every probe; needs M ∝ β eventually |

At fixed a the certified error therefore grows like β² relative to F_cert. The ν ≈ 0 residual mass is absolute and cannot be cured by
taking a larger. Numerically [certified evaluations of the m = 0 term, cell width 0.1]:

| β | 12 | 30 | 64 | 100 | 150 | 200 | 300 | 400 |
|---|---|---|---|---|---|---|---|---|
| G^r/F_cert, Lemma B only (old) | 0.0056 | 0.038 | 0.11 | 0.27 | 0.61 | 1.07 | 2.4 | 4.3 |
| G^r/F_cert, min(B, B′) (new) | 0.0043 | 0.019 | 0.043 | 0.096 | 0.20 | 0.35 | 0.76 | 1.32 |

### 3.2 What is and is not proved for β ≥ 64
* **Not proved here.** No uniform certificate on [64, ∞) exists with these ingredients: G^r/F_cert alone exceeds 1 at β ≈ 350 for every
  a. Beyond that point positivity of the m = 0 term itself (the leading Regge function with the realised weight) is not implied by
  absolute bounds on r. That needs the sign or oscillation of r near ν = 0, i.e. a contour treatment of the realised f.
* **Available by cells (not run here).** `run_large_scan.py --beyond B` adds width-0.1 cells on [64, B] with the LRINF tail floor (§1.8).
  From §3.1, a*(β) stays a few thousand up to β ≈ 150–200. Cost is ≈ 10 CPU-min per cell, i.e. ≈ 60 CPU-h per Δβ = 36.
* **What closes β → ∞.** The contour route. CONTOUR_BOUNDS certifies |ℰ₂| < 0.01𝓜(L) for all L ≥ 30.16 (β ≥ 31.5), every a ≥ 998:
  cells on [L₀, 400], then decaying terms beyond. ANGULAR_ROUTE §40's connection terms j ≥ 1 and k = 3 are still unchecked. Once those are
  done, they cover β ≥ 31.5 at every twist, which is more than this certificate could. So the high-twist cells are needed on
  [8.4, 31.6]. Running them to 64 gives a wide overlap with the contour threshold, in case that threshold moves up.

---

## 4. Cell plan, launcher, CPU estimate

* **Cells.** Width 0.1 on [8.4, 64]: 556 cells, aligned with every LRINF cell edge (0.05/0.1/0.5 grids). F_cert(cell) = min over the
  overlapping LRINF cells. Below β = 24 the LRINF cells are ≤ 0.1, so there is no loss. Above, a 0.1-cell inherits the 0.5-cell floor
  (≤ e^{0.25} loss).
* **Why 0.1.** Each order m loses e^{3mΔβ} through the single stadium radius (§1.1). Width 0.2 is about 2.4× worse in err/F at β = 30;
  width 0.5 fails outright.
* **Launcher.** `cd ads4/code/fullpos3; python3 run_large_scan.py [workers]`. It writes `cache/large_<ba>_<bb>.jsonl/.log` and skips
  cells already in `cache/large_*.jsonl` or `results/ht_large_done.jsonl`. `--summary` prints the table; `--merge` writes
  `results/ht_large_done.jsonl`; `--beyond B` appends cells on [64, B]. It reads none of the laptop job's files (`ht_scan_done.jsonl`,
  `laptop_*`, `ht_scan*`), and the laptop job reads none of its files.
* **CPU estimate.** Measured cell times here (wall, shared core, old ladder): 248 s (β = 9), 284 (12), 334 (16), 413 (20), 422 (30),
  689 (45), 1012 (64). Time grows with β through the precision (256 → 440 bits) and the far-remainder k-range (k ≈ 10³/d).
  * Interpolating over the 556 cells gives 90 h.
  * The finer twist ladder adds ≈ 20 %, for ≈ 90–110 CPU-hours in all, i.e. **≈ 4–5 h on 22 workers**.
  * [8.4, 31.6] alone (232 cells) is ≈ 25–30 CPU-h, ≈ 1.3 h.
  * Memory is ≈ 0.2 GB per worker.
  * If both launchers run at once, split the cores between them (e.g. `python3 run_large_scan.py 10`).

---

## 5. Blocks and open points
* **β → ∞:** see §3.2. The certificate's own limit is β ≈ 350 (m = 0 term); in practice ~200.
* **a* is set by the high-frequency residual pieces** (ν ≈ 900–1059), as in FULLPOS3 §4.6 at low β. Lowering a* below ~3000 needs the
  grouped / carrier-shifted local bounds (ANGULAR_ROUTE §64b, §66), not more cells.
* **The new lemmas also apply at β ≤ 8.4.** Lemma B′, Lemma A′ and the rest-split are valid for every β > 0. LOWBETA
  The laptop job was not changed. Re-running β ≤ 8.4 with `fp3_ht_large.py` would only sharpen a*.
* **Not done here:** the β ≤ 8.4 cells (laptop), the 556-cell production run, cells beyond 64.

---

## 6. Reproduce (`ads4/code/fullpos3/`; `ulimit -v 4000000; export OMP_NUM_THREADS=1`)
```
python3 fp3_ht_large.py 0 12 12.1 3200                    # one run (b, beta_a, beta_b, a0): JSON with the error budget per component
python3 fp3_ht_large.py scan "30,30.1;45,45.1" out.jsonl  # bisection in a0, both channels, -> cache/out.jsonl
python3 run_large_scan.py 22                               # production, laptop
python3 run_large_scan.py --summary
```
