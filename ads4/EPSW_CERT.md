# EPSW_CERT: ε_W for the k = 2 continuum weights in ball arithmetic, and the realised-weight T0 corollary at every β (M20)

*25 Nov 2026. Unreviewed, not committed. No existing file was edited, and nothing under `prl/` was touched.*
* *New code:* `ads4/code/epswcert/`:
  * `ew_common.py`, `ew_cert.py`, `ew_repair.py`, `ew_t0.py`, `ew_light.py`;
  * `ew_dv.py` ([DV] float check only);
  * `t3_flat.py`, a byte-identical copy of `t0cert/t3_flat.py`, re-run with a larger shift;
  * `run_all.sh`.
* *Outputs:* `epswcert/cache/`, whose `.gitignore` ignores everything except itself. The repaired weights are in there too (§3). Move them
  out if they are to be kept.
* *Read-only imports:*
  * `t0cert/t0_common.py` and, through it, `gaplarge/gl_common.py`, `lrtight/*`;
  * `kcert/kc_common.py`, `kc_dk0.py`;
  * `deshift` is not used; `robust/` is used through the `t3_flat` copy.
* *Data, read-only:* `prl/data/k2/k2_continuum_D_{300,1000}.json`, `prl/data/kernels.json`, `prl/data/nf02_psi.json`.
* *Compute:* about 15 CPU-min in total, ≤ 2 processes, `ulimit -v 3000000`, `OMP_NUM_THREADS=1`.

Grades are as in T0_CERT. **[certified]** means every inequality is proved here (or cited from a [P]/[certified] document), or is decided in
python-flint/Arb ball arithmetic or exact rationals with proved remainders. Nothing is sampled.

---

## 0. Verdict

**M20 is closed for a minimally repaired weight. It is only partially closed for the published weight, and that cannot be improved:
the published weight provably violates T0 at astronomically large impact parameter.**

| weight | Δ_g | c′ | β covered by the realised T0 corollary [certified] | beyond |
|---|---|---|---|---|
| published `k2_continuum_D_300` | 300 | **0.199** | **0 < β ≤ 24 658** (J + ½ ≤ m sinh(12 329)) | **H_realised < 0 on every `++` state with β ≥ 37 222** [certified] |
| published `k2_continuum_D_1000` | 1000 | **0.1999** (0.199: β ≤ 23 597) | **0 < β ≤ 23 585** | **H_realised < 0 on every `++` state with β ≥ 35 745** [certified] |
| **repaired** (`cache/k2_continuum_D_300_rep.json`) | 300 | **0.199** | **every β, every heavy state** (`++` and `+−`, every J, every m ≥ Λ) | — |
| **repaired** (`cache/k2_continuum_D_1000_rep.json`) | 1000 | **0.1999** | **every β, every heavy state** | — |

The c′ values are T0's own constants. The realised-weight terms cost no constant, only margin.

**Cause [certified].** The shifted density E(μ) = 2Re[χ(μ−i)R₋(μ−i)] of the published weights has a **negative endpoint**:

| | r(0) = E(0) − ψ(0) | r′(0⁺) = E′(0⁺) − ψ′(0) |
|---|---|---|
| Δ_g = 300 | **−1.25867479927·10⁻⁵** (± 4·10⁻¹⁷) | +1.89758165413·10⁻³ (ψ′(0) = 0.632819) |
| Δ_g = 1000 | **−3.92410426390·10⁻⁶** (± 5·10⁻¹⁸) | +5.99990561970·10⁻⁴ (ψ′(0) = 0.189180) |

This confirms the authors' enclosure −3.924105·10⁻⁶ < E(0) < −3.924103·10⁻⁶ at Δ_g = 1000. Their structural remarks also check out:
* R₋(−i) = 18/25 and Q(−i) = −16/9 (asserted with sympy);
* the log terms vanish at ν = −i.

Since ∫₀^∞(W − Φ)dt = (π/2)r(0), the realised image error is r(0)/(2 sinh(β/2)) + O(e^{−β/2}/β). It therefore overtakes the leading-Regge
floor A₁M₁/Λ² ∝ e^{−β/2}/β at β ~ 2–4·10⁴. No ε_W bound, however good, can remove this, because the defect is a sign.

**Repair [certified].** Change only `x_re[0]`, the coefficient of the c = 0 Gaussian, to the IEEE double
* 300: 135.99298935715314 → **135.99298894091532** (δx₀ = −4.162·10⁻⁷, 3.1·10⁻⁹ relative);
* 1000: 44.7154942909411 → **44.71549416117283** (δx₀ = −1.2977·10⁻⁷, matching the authors' −1.297683·10⁻⁷).

The new double is chosen so that r(0) is **certified positive**: +3.115·10⁻¹³ at 300 and +9.692·10⁻¹⁴ at 1000. r(0) = 0 exactly is not
representable with double coefficients, and the positive sign is the one that helps. The kink-matching step the authors proposed
(adjusting y₀) was tested and **rejected**: ε_W gets 7 600× worse and T0 Part K then fails (§3.2). It is not needed: r′(0⁺) > 0 only
raises the large-β floor, by 0.3 %.

**Light side.** The repair changes
* K = K₀ + dK₀ by < 10⁻¹⁵ relative;
* the realised Einstein value E = 16ΛV₃ by +6.9·10⁻⁷ (300) and +5.4·10⁻⁷ (1000);
* τ₀ by +3.5·10⁻⁷ and +2.7·10⁻⁷ relative.

**At 6 digits τ₀ is unchanged:** 3.13206·10⁻⁶ / 3.13326·10⁻⁶ (K = 1 / K_c) at 300, and **2.89126·10⁻⁸ / 2.89136·10⁻⁸** at 1000. Only the
7th digit of CERT_CLOSE's "R₀(K = 1+dk₀) = 3.133256·10⁻⁶ / 2.891361·10⁻⁸" moves, to 3.133257·10⁻⁶ / 2.891362·10⁻⁸ (§5). No headline
number changes at its quoted precision.

**ε_W certified (task 1)**, for all t ≥ 0, together with a decaying bound |W − Φ|(t) ≤ ε(t) := min(ε_W, B/t²):

| weight | ε_W [certified] | earlier [DV] value | B = \|r′(0⁺)\| + ∫\|r″\| [certified] | sup \|r\| |
|---|---|---|---|---|
| published, 300 | **≤ 3.5206·10⁻⁶** | 3.28·10⁻⁶ | ≤ 6.0047·10⁻³ | 1.2587·10⁻⁵ (at μ = 0) |
| published, 1000 | **≤ 1.1749·10⁻⁶** | 1.1·10⁻⁶ | ≤ 1.8985·10⁻³ | 3.9241·10⁻⁶ (at μ = 0) |
| repaired, 300 | ≤ 8.3425·10⁻⁶ | — | ≤ 6.0083·10⁻³ | 1.785·10⁻⁵ |
| repaired, 1000 | ≤ 2.6814·10⁻⁶ | — | ≤ 1.8998·10⁻³ | 5.556·10⁻⁶ |

**Δ_g = 2000, 5000 (task 3): not done.** It is not cheap, because no continuum weight exists there (§6).

**Proposed replacement for the T0 scope clause in the Supplement (2–3 sentences):**
> T0 is proved for the de-shifted (target-weight) leading-Regge action at every Δ_g ≥ 300. For the realised functional it holds at
> Δ_g = 300 and 1000 on every heavy state and at every impact parameter, with c = 0.199 and 0.1999 respectively [certified], provided
> the k = 2 continuum weight has a non-negative shifted endpoint E(0). This requires changing the published weights' first RBF
> coefficient by 3·10⁻⁹ relative; the published weights have E(0) = −1.3·10⁻⁵ / −3.9·10⁻⁶ and satisfy T0 only for β ≤ 2.4·10⁴
> (they fail it, H < 0, for β ≥ 3.6·10⁴), and no quoted number changes at its stated precision.

---

## 1. Set-up

* **Weight.** As published:
  `χ = Q·B + Σ_{j≤2} y_j ν^{2j} l(ν) e^{−ν²/2}`, `Q = ν²/(ν²+¼)²`, `l = (ν²+1)Log(ν²+1)`,
  `B = Σ_k x_k[g_k(ν−c_k) + g_k(ν+c_k)]`.
  All doubles of the JSON are taken **exactly**. Λ = √((Δ_g−2)(Δ_g−1)) exactly, k₀ = 27.48509817 exactly, and ψ is NF02's exact ψ₂.
* **Residual.** Only χ_re enters the leading-Regge image. The image of the real χ is real, and the functional takes Re.
  * Then `f(μ) = χ(μ−i)R₋(μ−i)`, `E = 2Re f`, `r = E − ψ` on μ ≥ 0.
  * For every real t, `W − Φ = ∫₀^∞ r(μ)cos(μt)dμ` (CERT_CLOSE §3.2, shifted form).
* **Exact algebra (sympy, asserted in code).**
  * `QR₋ = −32ν(2ν+3i)(2ν+5i)(2ν+11i)/((2ν−i)(2ν−7i)(2ν−3i)(2ν+7i))`, which removes the 0/0 at ν = −i/2.
  * `R₋(−i) = 18/25`.
  * Near μ = 0: `ν²+1 = μ(μ−2i)`, so `f = F₁ + F₂ log μ`, where
    * `F₁ = QR₋B + μ(μ−2i)Log(μ−2i)YR₋`;
    * `F₂ = μ(μ−2i)YR₋`;
    * `Y = Σy_jν^{2j}e^{−ν²/2}`.

    F₁ and F₂ are analytic in |μ| < 1, whose nearest singularity is R₋'s pole at ν = 0. `2Re F₂ = O(μ²)`, since F₂(0) = 0 and
    F₂′(0) = −2iY(−i)R₋(−i) is purely imaginary. So `r = r(0) + r′(0⁺)μ + O(μ² log μ)`. **r is C^∞ on (0, Λ) and (Λ, ∞), C⁴ at Λ
    (ψ has (1−p)⁵), and r″ ∈ L¹ with a log singularity at 0.**

## 2. The certificate for ε_W, B, r(0) (`ew_cert.py`) [certified]

### 2.1 Method
Template: HK_RESIDUAL §3.3, but with exact power-series arithmetic instead of contour trapezoids.

**[0, ⅕].**
* Taylor coefficients of F₁ and F₂ at 0 to order 64, by Arb series arithmetic. These are exact truncations.
* Cauchy truncation bounds use the maximum modulus of |F₁| and |F₂| on |μ| = 0.8, enclosed by a ring of 256 boxes:
  M₁ ≤ 1.35·10⁶ and M₂ ≤ 23 at 300.
* 131 cells: [0, 10⁻³·2⁻³¹], then geometric cells up to 1/500, then uniform cells of width 1/500.
* On every cell with u > 0, log μ = log c + log(1+s/c) is replaced by its degree-48 Taylor polynomial (|s/c| ≤ ⅓, explicit remainder).
  The cancellation between the analytic part and the μ² log μ part is therefore kept. A first version that bounded the two separately
  lost a factor 10³.
* ∫|r| and ∫|r″| come from polynomial Taylor shifts to the cell centres, plus all remainder terms, including
  `(T₂ log μ)″ = T₂″ log μ + 2T₂′/μ − T₂/μ²`.

**[⅕, ⌈Λ⌉+60].**
* 3602 panels at 300 and 10 602 at 1000, of width 1/40–1/10, split exactly at Λ.
* Taylor series of f to order 40 at each panel centre. The near RBF terms are those with |m − c_k| ≤ 15σ_k. The far terms are bounded on
  the disc by an explicit Gaussian majorant, ≤ 10⁻³⁴.
* Cauchy truncation on a disc of radius ρ = min(½, 0.6m), with h/ρ ≤ 0.2. The disc sup is bounded by
  * an explicit RBF bound times the ring-box enclosures of |QR₋|;
  * plus the ring-box enclosure of |Y l R₋|.

  The maximum-modulus principle applies because both are analytic on the closed disc (m − ρ > 0.05).
* 5 subcells per panel, by Taylor shift. The truncation-plus-far error is ≤ 1.6·10⁻²² per point.

**Tail μ > ⌈Λ⌉+60.** An erfc bound for the RBF terms (centres ≤ Λ+30), a Gaussian bound for the y-terms, and Cauchy for r″. Both tails
are ≤ 10⁻²⁴⁶.

### 2.2 Results
`cache/ew_cert_<Dg>_re[_rep|_repk].txt`.

| published | [0, ⅕] | (⅕, 1) | [1, 40) | [40, Λ) | (Λ, ∞) | **∫\|r\|** | **∫\|r″\|** |
|---|---|---|---|---|---|---|---|
| 300 | 5.87·10⁻⁷ | 9.50·10⁻⁷ | 1.635·10⁻⁶ | 3.45·10⁻⁷ | 3.5·10⁻⁹ | **≤ 3.5206·10⁻⁶** | **≤ 4.1071·10⁻³** |
| 1000 | 1.82·10⁻⁷ | 2.95·10⁻⁷ | 5.21·10⁻⁷ | 1.77·10⁻⁷ | 3·10⁻¹⁰ | **≤ 1.1749·10⁻⁶** | **≤ 1.2986·10⁻³** |

The certified ε_W is within 7 % of the [DV] quadratures of CERT_CLOSE §3.2, and above them, as it must be.

### 2.3 Consequences [P]
1. **Uniform bound.** `|W − Φ|(t) ≤ ε_W` for all t.
2. **Decaying bound.** Two integrations by parts give `W − Φ = −[r′(0⁺) + ∫₀^∞r″cos μt dμ]/t²`, since r and r′ vanish at ∞. Hence
   **|W − Φ|(t) ≤ B/t²** for t > 0.
3. **Mean of the error.** `∫₀^∞(W − Φ)dt = (π/2)r(0)`. This is Fourier inversion for the even extension r_e: r_e is continuous,
   r_e ∈ L¹, and W − Φ ∈ L¹ by item 2.

### 2.4 [DV] check of the mechanism (`ew_dv.py`, float)
(2/π)∫₀^T(W − Φ)dt = (2/π)∫r sin(μT)/μ at Δ_g = 300 is

| T | 10² | 10³ | 10⁴ | 10⁵ |
|---|---|---|---|---|
| value | −1.50·10⁻⁶ | −1.138·10⁻⁵ | −1.247·10⁻⁵ | −1.2575·10⁻⁵ |

It converges to r(0) = −1.2587·10⁻⁵, slowly: the deficit is ≈ 2r′(0⁺)/(πT), from the 0.007-wide feature of r near 0. This is why
direct β ≤ 60 quadratures of the image (also run) still see a positive error, +5·10⁻⁷/(2 sinh(β/2)). The negative constant only takes
over at β ≳ 10³. The Δ_g = 1000 run behaves the same way.

## 3. The repair (`ew_repair.py`)

### 3.1 Endpoint repair ('rep') [certified]
* **Exact derivative.** r(0) = 2Re[Q(−i)R₋(−i)B(−i)] (the l-terms vanish at ν = −i). Hence `∂r(0)/∂x₀ = −(64/25)e^{1/(2σ₀²)}` exactly,
  = −30.2393 for σ₀ = 0.45.
* **Construction.** x₀ is moved to the nearest double of the zero, then by ulps until r(0) > 0 is certified. One step suffices at both
  gaps.
* **Format.** Same as `prl/data/k2`. Only `x_re[0]` differs, and a `repair` field records the change. The old `eps_W_*`/`certificate`
  fields are replaced by a pointer to this certificate.
* **Effect on ε_W.** It grows to 8.34·10⁻⁶ / 2.68·10⁻⁶, because the correction δx₀·Q g₀ is a bump of width ~0.5 near μ = 0. It is still
  far inside every T0 margin (§4).
* **Possible refinement (not done).** A two-coefficient repair (x₀, x₁) that also holds V₃ fixed would leave the light side identical.

### 3.2 Kink repair ('repk'): rejected
* **What it was.** y₀ adjusted so that r′(0⁺) = 0 (∂r′(0⁺)/∂y₀ = −2π(18/25)e^{½}, checked in Arb). This makes r′(0⁺) ≈ 10⁻¹⁶.
* **Why it fails.** The log-basis function is not localised. ε_W becomes **6.36·10⁻²** (300) and **2.01·10⁻²** (1000), and
  B becomes 0.22 / 0.071. T0's Part K margin (4.58/Λ²) would then fail: 2(1+c)Λε_W = 46.
* **Not needed.** The realised/target floor ratio tends to 1 + r′(0⁺)/κ = 1.0030, and r′(0⁺) > 0 helps.

## 4. The realised-weight corollary (`ew_t0.py`) [certified]

### 4.1 Convention and image identity
* **Convention (T0_CERT §4.2).** The realised weight changes only the leading-Regge image carried by the hybrid's "2"-curvature:
  `X → X + 2δI/Λ`. For pattern g, Λ²·(curvature part) gains `(1−cs₂)·2Λ δI`.
* **Image identity.** `δI(β) = (√2/π)∫₀^β (W−Φ)(t)K(t)dt`, with K = (cosh β − cosh t)^{−½} and K₀ = K(0) = 1/(√2 sinh(β/2)).

### 4.2 Two image bounds [P + certified inputs]
* **(U) Uniform.** `|δI| ≤ ε_W P_{−½}(cosh β) ≤ ε_W min(1, β/(2 sinh(β/2)))`.
  This uses ∫₀^βK = √2∫₀^{π/2}dθ/√(1+S²sin²θ) ≤ (π/2)βK₀, from sin θ ≥ 2θ/π.
* **(A) Asymptotic.** `δI = r(0)K₀/√2 + R` with `|R| ≤ (√2/π)K₀ρ_R`, where for β ≥ 64
  `ρ_R = B/β + 8B ln(4cosh(β/4))/β² + 0.51(ε_W sinh 2 + 2B cosh(β/2)/β)/(cosh β − 1)`.
  The three pieces are:
  * K₀∫_β^∞|W−Φ| ≤ K₀B/β;
  * on [β/2, β], |W−Φ| ≤ 4B/β² and ∫_{β/2}^βK ≤ 2K₀ ln(4cosh(β/4));
  * on [0, β/2], K − K₀ ≤ 0.51K₀(cosh t − 1)/(cosh β − 1).

  βρ_R is decreasing for β ≥ 64.

### 4.3 Region by region (T0's certified margins, which are cited inputs)

| region | realised requirement | 300 / pub | 300 / rep | 1000 / pub | 1000 / rep |
|---|---|---|---|---|---|
| b̂ ≤ 16 (Parts F, T): Part-F shift − T_abs ≥ (1+c)2ε_W/Λ | slack 6.43·10⁻⁷ (300), 9.25·10⁻⁷ (1000, c = .199), **2.03·10⁻⁷** (1000, c = .1999, `t3_flat` copy re-run with shift 2.93·10⁻⁵: PASS) | need 2.8·10⁻⁸ ✓ | 6.7·10⁻⁸ ✓ | 2.8·10⁻⁹ ✓ | 6.4·10⁻⁹ ✓ |
| 16 ≤ b̂ ≤ 2Λ (Part K): min K_g ≥ (1+c)2Λε_W | K_g ≥ 4.585 / 4.934 / 4.923 | 2.5·10⁻³ ✓ | 6.0·10⁻³ ✓ | 2.8·10⁻³ ✓ | 6.4·10⁻³ ✓ |
| 1.76 < β ≤ 64 (Part L, 564 cells), with (U) | realised K/target K ≥ | 0.902 ✓ | 0.768 ✓ | 0.891 ✓ | 0.750 ✓ |
| β ≥ 64, with (A) and T0's tail | see below | ✓ to 24 658 | **✓ all β** | ✓ to 23 597 (c = .1999: 23 585) | **✓ all β** |

### 4.4 β ≥ 64
**Target floor (re-evaluated with the cited T0 code at Λ_a = Λ(Δ_g)).**
* Λ²curv_g ≥ (S_g − (1+c)tlA)M₁ − (1+c)FL. Here tlA is LR∞'s tail majorant (`gl_common.tail_majorant`), and FL is Lemma 5's flat remainder.
* M₁ ≥ (√2/π)K₀L(β) (LRINF_TIGHT region (c)), with βL(β) increasing.
* Dividing by (√2/π)K₀/β, every pattern g (16 combinations of s₂, s₃, σ, u ∈ {0, 1}; affine in u) needs
  `(S_g − (1+c)tlA)·βL(β) − (1+c)βFL·π sinh(β/2) − (1−cs₂)2Λ[(π/2)max(0, −r(0))β + βρ_R(β)] > 0`.

**Repaired weights.** r(0) > 0, so the r(0) term is dropped. Every other term is monotone in the safe direction, so **one evaluation at
β = 64 covers all β ≥ 64**. The minimum over patterns is ≥ 202.6 (300), 202.0 (1000, c = .199) and 201.7 (1000, c = .1999).

**Published weights.** The r(0) term grows linearly in β. Geometric β-cells (ratio 1.0005) pass up to **β₁ = 24 658.9 (300),
23 597.7 (1000, .199), 23 585.9 (1000, .1999)**.

### 4.5 Certified failure of the published weights
Take every `++` state (even J, every u ∈ (0, 1]). Bound H^R from above:
* **Flat parts.** Legendre's equation and two integrations by parts in x = 1 − 2up² give
  |B₂| + |B₃| ≤ (T₂+T₃)/(4uJ(J+1)) = (T₂+T₃)/(b̂² − u), with T₂ + T₃ ≤ 8.85·10⁸. The boundary terms vanish: at p = 0 through (1 − x²),
  at p = 1 through ψ's (1 − p)⁵. Only |P_J| ≤ 1 is used.
* **Curvature.** ≤ (2a₁₂ + a₁₃ + tlA)M₁ + FL, with M₁ ≤ (√2/π)K₀[1/β + (0.51/24)sinh(β/2)/(cosh β − 1)]. This uses:
  * c(t) ≤ c(0) = 1/24;
  * c < 0 on [β/2, β] (h(t) ≤ e^{−t/2} < t⁻² for t ≥ 32);
  * ∫₀^βc ≤ 1/β.
* **Realised term.** ≤ 2Λ[(π/2)r(0) + ρ_R], with r(0) < 0.

Every term except the r(0) term decreases in β. So at **β₂ = 37 222 (300), 35 763 (1000, .199), 35 745 (1000, .1999)** the bound is
< 0 (e.g. Λ²H^R/((√2/π)K₀) ≤ −4.1·10⁻⁶ at 300), and it stays < 0 for all β ≥ β₂. **H_realised < 0 ≤ cN_s there.**

## 5. Light side with the repaired weight (`ew_light.py`; kcert's `kc_dk0` imported read-only, χ rebound to the repaired file)

| Δ_g | quantity | published | repaired | relative change |
|---|---|---|---|---|
| 300 | dK₀ | −7.66632448297·10⁻⁴ | identical | < 10⁻¹⁶ |
| 300 | E = 16ΛV₃ | 10734.6054683 | 10734.6128999 | +6.92·10⁻⁷ |
| 300 | τ₀(K=1) / τ₀(K_c) | 3.13205501·10⁻⁶ / 3.13325626·10⁻⁶ | 3.13205609·10⁻⁶ / 3.13325735·10⁻⁶ | +3.46·10⁻⁷ |
| 1000 | dK₀ | −6.83014998391·10⁻⁵ | identical | < 10⁻¹⁶ |
| 1000 | E = 16ΛV₃ | 14339.5890757 | 14339.5968260 | +5.40·10⁻⁷ |
| 1000 | τ₀(K=1) / τ₀(K_c) | 2.89126258·10⁻⁸ / 2.89136133·10⁻⁸ | 2.89126336·10⁻⁸ / 2.89136211·10⁻⁸ | +2.70·10⁻⁷ |
| both | M₀ … M₃ (H17 sensitivities) | — | — | ≤ 10⁻⁹ (M₀), ≤ 2·10⁻¹⁴ (M₁₋₃) |

**τ₀ to 6 digits is unchanged:**

| | τ₀(K = 1) | τ₀(K_c) |
|---|---|---|
| Δ_g = 300 | 3.13206·10⁻⁶ | 3.13326·10⁻⁶ |
| Δ_g = 1000 | **2.89126·10⁻⁸** | **2.89136·10⁻⁸** |

* The 7-digit R₀(K = 1+dk₀) values quoted in CERT_CLOSE §3.2 move in their last digit: 3.133256 → 3.133257 (·10⁻⁶) and
  2.891361 → 2.891362 (·10⁻⁸).
* The realised Einstein offset V₃ − L_A/Λ grows from 7.4·10⁻¹⁰ to 4.9·10⁻⁷ (1000). It stays within the (4/9)ε_W bound.
* H17 was not re-run (kb2a4). Its inputs M_j move by ≤ 10⁻⁹, and τ₀ itself moves by the factor √(E_rep/E_pub) above. So every kb2a4/kcert
  τ₀ bound moves by ≤ 3.5·10⁻⁷ relative.
* K, 2L_AdS = 92.8135 (computed from ψ only) and every flat or quartic number are unaffected.

## 6. Δ_g = 2000, 5000 (task 3): not done
* **Why the grid weights cannot be used.** They are grids, and their de-shift identity fails beyond t ≈ 18 (README §8), so no all-t ε_W
  exists for them. The realised image is then uncontrolled at β ≳ 18.
* **What it would take.**
  1. A continuum fit with `code/cert/k2_cont.py` (≈ 2.5 CPU-h per gap). At 5000 the dense QR has ~1.8·10⁴ unknowns and needs > 3 GB, so
     it needs a blocked solver.
  2. Add the constraint E(0) ≥ 0, e.g. by the x₀-repair of §3.1.
  3. Run `ew_cert.py` (~2–5 min) and `ew_t0.py` with T0's Δ_g ≥ 1000 inputs. T0's bounds are valid for every Λ ≥ Λ(1000).
* **Expected outcome.** The margins scale favourably: the realised term is ∝ Λε_W, with ε_W ∝ Λ⁻¹ observed. So the expected outcome is
  "every β" at c = 0.1999.

## 7. Reproduce
    cd ads4/code/epswcert && bash run_all.sh        # ~ 15 CPU-min, <= 2 processes

| file | content | grade |
|---|---|---|
| `ew_common.py` | exact kernels and weight, Arb series for f, F₁, F₂, disc enclosures | — |
| `ew_cert.py DG re [pub\|rep\|repk]` | r(0), r′(0⁺), ε_W, ∫\|r″\|, B, sup\|r\| | [certified] |
| `ew_repair.py` | repaired JSONs (`cache/k2_continuum_D_<Dg>_{rep,repk}.json`) | [certified] sign of r(0) |
| `t3_flat.py` | byte copy of `t0cert/t3_flat.py`, run with shift 293/10⁷ (Δ_g = 1000, c = 1999/10000) | [certified] |
| `ew_t0.py` | §4: Parts F/K/L margins, β ≥ 64 tail, β₁, certified failure β₂ | [certified] |
| `ew_light.py` | §5, via kcert read-only | [certified] (ball) |
| `ew_dv.py` | §2.4 | [DV] |

## 8. Caveats
1. **Cited inputs.** These are T0's certified margins (Parts F, K, L, read from `t0cert/cache` and recomputed for Part L), LR∞'s tail
   majorant, the M₁ lower bound L(β) with βL increasing (LRINF_TIGHT §1.5), and Lemma 5. Their grades are as in T0_CERT and
   INDEP_CHECK_T0. The convention that the realised weight acts only on the "2"-curvature is T0_CERT §4.2's.
2. **What "every heavy state" means.** For the repaired weight it means `++` and `+−` at every J, every m ≥ Λ and every b̂, exactly T0's
   class. The OO caveat of T0 is unchanged. The certified failure of §4.5 is shown on `++` states only, which suffices.
3. **The repaired weight is a new definition of Ω^EE's k = 2 piece.** Adopting it means publishing the two changed coefficients. They are
   in `epswcert/cache/`, which is git-ignored, so they must be copied out to be kept. I did not touch `prl/`.
4. **χ_im** enters no leading-Regge or light value (Re), and it was not certified (`ew_cert.py DG im` exists but was not run).
