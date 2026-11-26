# FULLPOS4: the outgoing-split, frequency-integrated §55 contour at Δ_g = 1000

*26 Nov 2026. Research agent. Unreviewed, not committed. No existing file was edited; nothing under `prl/` was touched.*
*New code: `ads4/code/fullpos4/` (outputs in `fullpos4/cache/`, whose `.gitignore` ignores everything but itself).*
*Compute: ≤ 2 processes, `ulimit -v 4000000`, `OMP_NUM_THREADS=1`. The two FULLPOS3 scan processes were not touched.*

Grades: **[certified]** Arb balls or exact rationals with every remainder proved; **[P]** proved on paper here;
**[DV]** numerical; **[conj]** conjecture. Author/coordinator claims are marked as such and kept separate.

*(Running document; sections marked "in progress" are partial.)*

---

## 0. Verdict

**Not closed. No new certified coverage.** The sliver 998 ≤ a < a*(L) (a* = 3200–3600 after ANGULAR_ROUTE §73 and FULLPOS5),
J ≥ 1001, 0.04 < L ≲ 30 is still open. But the session ends with a concrete, numerically validated route (item 6), and quantified
obstructions for the others.

1. **The §55 outgoing contour works per frequency [DV].**
   * Author's recipe reproduced at the benchmark state; every piece matches to 8 digits.
   * 12-state sweep against the exact ladder sum: p = 3, 4; b = 0, 4; ν = 0.3–12; non-integer N. Rel. err 1e-14 to 2e-7.
   * The outgoing split also holds for ν ≤ 1 (rel. err 1e-12).
   * FULLPOS3's NaN state is fine (1.415451002331).
   * Trap: mpmath `fp.hyp2f1` is wrong by 1e-4 near e^{−2β} → 1.
2. **Frequency integration first works [DV].** 𝒫⁺_h = −(i/π)∫W_h(β+u)(4 sinh(β+u/2)sinh(u/2))^{−½}du is checked to 1e-7–1e-10.
   End-to-end, with the author's §74 diagonal extraction, against exact ladder sums: 6e-11, and 1.5e-7 at a ladder cancelling by 4e7.
3. **Contour norms [DV]** (units |G(c)| = |ℒ(G^ψ)| ≈ N_s/w³):
   * fully absolute: 2.2e10 (L = 0.2), 1.2e10 (L = 2.4), Lemma-A⁺-type majorant 3.1e10;
   * soft-signed: 2.6e6 and 1.3e4;
   * after the §74 diagonal extraction, falling like ξ₀^{−2}: 61 at ξ₀ = 0.35 (L = 0.2), 97 → ~16 at L = 2.4.
   * **Brute-force contour certification does not fit:** ~1e7 amplitude evaluations per state, ≥ 1e3–1e5 cells.
4. **§78 quadrature-free reference: exact pointwise** [ball arithmetic + DV], to 1e-77 … 0 at three states including the turning point
   ν = 993. Integrated with ψ₂ at (1500, 2000): H_ref = H_direct to all double digits.
   * Naive interval cells in a fail unless δa ≲ 1e-550 (ball wrapping in the M-step recurrence).
5. **b = 4:** 𝒞_{s=9} = a/U_{0,4}, so the helicity shift does not change the integrand (checked).
6. **Lemma U (author §80) twist test: PASSES at every L tested, both channels [DV, ideal NF02 weights, M = 1000].**
   At fixed c, a single degree-1001 polynomial in y covers all a ≥ 998:
   * margin (P̂₀ − Σ_{ℓ≥1}|P̂_ℓ|)/P̂₀ = 0.972 / 0.971 (L = 0.1), 0.985 / 0.985 (L = 0.5), 0.9976 / 0.9976 (L = 2), 0.99999 / 0.99998 (L = 7),
     for b = 0 / 4.
   * The h-series at a = 998 has a 1e11–1e18 transient in A_m (m ≈ 40–100), then converges: partial sums are constant from M′ = 300.
   * Complex impact parameter (Lemma C test): max_y|H_M(β_c + iη, y)|/F_cert = 1.26 (L = 2, η = 0.02) and 1.27–1.40 (L = 0.5, η = 0.02).
     Low L is in progress (§8.1).
7. **What is missing for a certificate (honest):**
   1. rigorous jet enclosures (GL error on an ellipse, plus ≥ 6000–8000 bits against recurrence wrapping);
   2. the realised EPSW-repaired weight instead of the ideal NF02 one;
   3. an Arb reproduction of the §78 truncation error;
   4. the c-direction: ~1.4e3 c-points × 15–30 min, i.e. several hundred core-hours DV and more certified (§8.2). The strip constant
      K = sup|H_M|/F_cert is ~1.3–1.4 at L ≥ 0.5, so cells there can be ~η wide; at low L the author measured a transfer > 463.

**b = 4.** The coordinator's identity ℒ₄(F) = ℒ₀(z⁴F)/U_{0,4} is correct, but 𝒞_{s=9} = a/U_{0,4} exactly (checked to 1e-12 at 3 states).
So the shifted s = 1 integrand equals the s = 9 integrand pointwise and the norm is unchanged (author §72 agrees). The cancellation is the
I_U/I_T corner (|ξ|⁸/U_{0,4} ~ c^{−4} away from |ξ|² ≈ c). The soft-signed evaluation pairs the corner exactly; the b = 4 norm after the
§74 extraction was not measured (the runs were stopped to free cores).

---

## 1. The §55 contour, per frequency: reproduced [DV]

### 1.1 Conventions (checked)
Notation of FULLPOS3. With x = e^{iα}, y = e^{iγ}:
* ξ = (y−x)/(xy−1) = sin((γ−α)/2)/sin((γ+α)/2), z = ξ², β = Log((1+ξ)/(1−ξ)) = log tan(γ/2) − log tan(α/2);
* u/v = xy and (1+u)/(1+v) = y/x, so the §53 phase becomes e^{iSα − iaγ}.

The §59 radial-soft form, rewritten in (α, γ) with γ = α − iv, is
  𝒦_p = (𝒞/2π)∫_{α=−π/2−iε}^{π/2−iε}dα∫_{γ=α}^{α−i∞}dγ W,  W = e^{iSα−iaγ}ξ^{s−1}F_{p,ν}(ξ²)/sin α.
* With the complete conical function it reproduces the exact ladder sum: (p, ν, J, b, Δ) = (3, 1.2, 10, 0, 32.5) gives 1.1592763187111, rel. err 1.3e-14 (`fp4_pf.py … base full`).
* **Outgoing split on the same path.** P⁺_ν(β) = Γ(iν)e^{(−½+iν)β}₂F₁(½, ½−iν; 1−iν; e^{−2β})/(√πΓ(½+iν)), and 𝒦_p = 2Re I⁺ with P → P⁺ in W.
  The rel. err is 1.2e-14 at the same state.
* Identity used below: **P⁺_ν(β) = −(i/π)coth(πν)Q_{−½−iν}(cosh β)** (Γ(iν)Γ(1−iν)/(πΓ(½+iν)Γ(½−iν)) = −i coth(πν)/π). It is checked
  against mpmath to 2e-14 at (ν, β) = (3.3, 0.4+0.2i) and (7, 0.05+0.01i).

**Numerical trap found.** mpmath's double-precision `fp.hyp2f1` is wrong by 1e-4 relative near e^{−2β} → 1 (the zero-balanced log case).
That alone produced 1e-5 errors in the first tests. `fp4_pf._h2f1_zb` uses the direct series for |w| ≤ 0.6, the ψ-function log expansion about w = 1
for |1−w| ≤ 0.75, and 30-digit mpmath otherwise.

### 1.2 The author's §55 recipe (coordinator's message): reproduced
I_U: α = A−iε, γ = C−iε, 0 ≤ A ≤ C ≤ π/2; I_T: γ = π/2 − it, t ≥ ε; I_L: α = A−iε, γ = A − i(ε+v), −π/2 ≤ A ≤ 0.
𝒦_p = (𝒞/π)Re(I_U − iI_T − iI_L) for ν > 1 (P⁺), and (𝒞/2π)(I_U − iI_T − iI_L) with the complete conical function for ν ≤ 1.
`fp4_s55.py` (Gauss–Legendre panels graded at A = 0; t⁴ maps at the diagonal; ξ rational, β principal log).

**Benchmark state:** n = 2, N = 7.5, b = 0, p = 3.
| ν | I_U | I_T | I_L | contour | exact (ladder) | rel. err |
|---|---|---|---|---|---|---|
| 1.2 | 0.50721927 − 1.81183203i | −1.69119565 + 0.02922216i | 0.04932043 + 0.10796520i | 1.1281655251235 | 1.1281655251250 | 1.4e-12 |
| 0.8 | (complete conical) | | | 1.21738448201959 | 1.21738448201953 | 5e-14 |

All three pieces agree with the author's values to the 8 printed digits.

**Sweep** (`fp4_s55_sweep.sh` → `cache/fp4_s55_sweep.txt`). Exact = terminating ladder sum Σ_{k≤n}U_kg_k (mpmath, 30 digits).
| n | N | b | p | ν | rel. err | ∫\|W\|/\|𝒦\| |
|---|---|---|---|---|---|---|
| 2 | 5.5 | 0 | 3 | 0.8 | 1.1e-14 | 16 |
| 10 | 30.5 | 0 | 3 | 1.2 | 1.4e-14 | 6.2e2 |
| 7 | 23.25 | 0 | 4 | 2.7 | 2.5e-12 | 9.2e2 |
| 5 | 26.7 | 4 | 3 | 0.6 | 5.7e-8 | 2.6e7 |
| 6 | 34.3 | 4 | 4 | 4.1 | 3.4e-8 | 1.2e8 |
| 6 | 17.9 | 0 | 3 | 5.5 | 2.9e-12 | 20 |
| 4 | 19.75 | 4 | 4 | 0.3 | 2.4e-7 | 8.5e7 |
| 10 | 43.5 | 0 | 3 | 12.0 | 1.3e-11 | 1.1e2 |
| 8 | 40.2 | 4 | 3 | 7.7 | 2.1e-10 | 8.8e5 |

The state (J, b, Δ, ν, p) = (2, 0, 7.5, 0.8, 3) is the one where FULLPOS3 §5.3 got NaN. It now gives 1.415451002331 (exact 1.415451002331).
The cause of the NaN was the literal diagonal/branch reading, not the contour.

**Observation (b = 4).** I_U ≈ i·I_T to 5–7 digits, and 𝒦 is their small difference: ∫|W|/|𝒦| = 1e5–1e8.
* The corner γ = π/2 − iε carries large ξ ≈ 1 values ((1−ξ²)^{−p}) times 𝒞 ≈ a(N/n)^{2b}.
* The b = 4 integrand is concentrated where the two pieces cancel. A better tail position (γ₁ < π/2), or the ξ^{2b} zero exploited
  through the diagonal, should be used for b = 4.

### 1.3 My own deformed contour (independent check)
For A > 0, γ runs α → A + 3ε − iε (horizontal) → A + 3ε (up) → π/2 (real axis) → π/2 − i∞; for A < 0 the §59 ray.
* On the first three pieces Re β ≥ 0 and Im β ≥ 0. This is proved on paper: |tan((x−iy)/2)|² = (cosh y − cos x)/(cosh y + cos x) and
  arg tan((x−iy)/2) = −atan(sinh y/sin x).
* It gives 1.415450861 (exact 1.415451002) and 1.15927627 (exact 1.15927632), i.e. 1e-7 and 5e-8.
* ∫|W| = 571 vs 2.6 on the base path at (3, 1.2, 10, 0, 32.5): staying on Im γ ≥ −ε costs (1−ξ²)^{−p} ~ (tan(γ/2)/|tan(α/2)|)^{p−½}
  where e^{−iaγ} gives no damping.
* **A pitfall:** going straight up to the real axis at A ≪ ε makes Re β ≪ 0 (|tan(γ/2)| < |tan(α/2)|), with |F⁺| ~ (ε/A)^p
  (rel. err 0.19 in a first version).

---

## 2. Frequency-integrated outgoing amplitude 𝒫⁺_h [DV]
`fp4_int.py`. With h(ν)dν = 2ψ₂(p)dp (so G^ψ = (1−z)^{−3}(𝒫⁺_h + 𝒫⁻_h)):
  𝒫⁺_h(β) = −(i/π)∫₀^∞W_h(β+u)(4 sinh(β+u/2)sinh(u/2))^{−½}du,  W_h(t) = 2∫₀¹ψ₂coth(πΛp)e^{iΛpt}dp.
This equals the author's §58/§71 form (−(i/√2π)∫K W_h).
* For |t| ≥ 0.04, W_h = W_S + W_E + W_C exactly: the p = 0 endpoint, the p = 1 endpoint (carrier e^{iΛt}, from j = 5), and the coth correction.
* W_E is integrated on the ray u = re^{iπ/4}, where e^{iΛt} decays.
* For |β| < 0.04, u ≤ 0.06 uses W_h directly (Gauss–Legendre in p).

**Checks.**
* 2Re𝒫⁺_h(β) against the Mehler value 2·OmLR(β) (largespin/ls_core) at β = 0.01 … 8: agreement 1e-10 (β ≤ 0.1) to 7e-7 relative (β = 8).
* 𝒫⁺_h at complex β against direct p-quadrature of mpmath P⁺_{Λp} (1280 nodes): 0.5+0.3i and 0.2+1.0i agree to 1e-8.

**Size of the outgoing amplitude (real β):**
| β | 0.01 | 0.03 | 0.1 | 0.3 | 0.6 | 1 | 2 | 4 | 8 |
|---|---|---|---|---|---|---|---|---|---|
| \|𝒫⁺_h\| | 1.10 | 0.134 | 1.2e-2 | 1.3e-3 | 2.9e-4 | 8.9e-5 | 1.3e-5 | 2.2e-6 | 2.2e-7 |
| 2Re𝒫⁺_h (= G-floor) | 0.89 | 3.9e-2 | 1.0e-3 | 5.3e-5 | 2.0e-5 | 1.5e-5 | 1.1e-5 | 4.3e-6 | 4.4e-7 |

* The outgoing amplitude is **not** floor-sized for β ≲ 1. Its imaginary part ≈ ψ₂′(0)/(Λβ)²-type (the odd-j p = 0 terms, which cancel in
  P = P⁺ + P⁻ but not in P⁺ alone) is 25× the floor at β = 0.3 and 7× at β = 1.
* At β ≲ 1/Λ, |𝒫⁺_h| → O(∫|ψ₂|) (log-singular).

### 2.1 The outgoing split also holds for ν ≤ 1 on the author's contour [DV]
With P⁺ forced at ν ≤ 1 on the §55 contour, (𝒞/π)Re(I_U − iI_T − iI_L) reproduces the exact kernel:
* (n, N, b, p, ν) = (2, 7.5, 0, 3, 0.8): rel. err 1.3e-12;
* (2, 5.5, 0, 3, 0.5): 6e-13;
* (9, 31.1, 0, 4, 0.95): 1.1e-11.

The complete-conical treatment below ν = 1 is therefore optional. The integrated amplitude can use P⁺ at every frequency, since h(ν) ∝ ν cancels
the 1/ν of coth(πν).

### 2.2 End-to-end check of the integrated pipeline [DV]
Test (`fp4_toy.py`): the same ψ₂ weight at reduced bandwidth Λ_toy = 8, i.e. G = 2∫ψ₂F_{3,8p}dp, at (n, N, b) = (10, 30.5, 0).
The exact value is ℒ(G) = Σ_{k≤n}U_kM_k with M_k = 2∫ψ₂g_k(8p)dp.
* First runs: 8.665890 against the exact 8.665246 (rel. err 7e-5). Cause: a bug in `fp4_int.Weight`. The coth-correction window was
  min(0.015, 6/λ); at λ = 8 this truncates 2/(e^{2πλp}−1) at e^{−0.75}. It is now min(1, max(0.015, 6/λ)). This matters only for toy
  bandwidths; at Λ = 998.5 the window was already complete.
* After the fix, with the §74 diagonal extraction (§4 below): D + R = **8.665245954894 vs exact 8.665245954376, rel. err 6e-11.**
So frequency integration first, then the §55 contour with the outgoing amplitude (plus the diagonal block), reproduces the action.

---

## 3. The norm of the frequency-integrated contour: structure [P-level analysis] and measurements [DV]

### 3.1 Two regions where the integrated outgoing amplitude is not floor-sized [P-level]
Write the contour integrand as W = (𝒞/π)ξ^{s−1}(1−ξ²)^{−p}𝒫⁺_h(β)e^{iSα−iaγ}/sin α.

**(i) The diagonal γ = α (β = 0, z = 0) is a forced boundary of every §53/§55/§59 contour.** It is the image of the corner u = v = 0.
* Near it, dγ = sin γ dβ ≈ sin α dβ. So per unit A the diagonal strip contributes (𝒞/π)e^{1−aε}∫|𝒫⁺_h(β)||ξ|^{s−1}dβ, independent of α.
* For b = 0 (ξ⁰), measured (`fp4_int`): ∫₀^{0.04}|𝒫⁺_h|dβ ≈ 0.075, of which 0.044 comes from β ∈ [10⁻³, 10⁻²]. The β ≤ 10⁻⁴ part is 0.003, and
  ∫_{0.04}^{∞} ≈ 0.003.
* Reason: at β ≲ 1/ν, P⁺_ν = −(i/π)coth(πν)Q_{−½−iν} ~ (i/π)log(νβ) carries the full weight. At 1/Λ ≲ β ≲ 1 the odd-j p = 0 endpoint terms
  ψ^{(j)}(0)/(iΛt)^{j+1} survive: they cancel only in P = P⁺ + P⁻ (Lemma A), not in P⁺ alone.
* **Consequence.** Diagonal norm ≈ (𝒞/π)·(retained A-range)·0.075·e^{1−aε}.
  * On the author's contour the A-range is π, giving ≈ 0.2𝒞. That is ~400 for b = 0, a = 2000, against |ℒ(G^ψ)| ~ 4·10⁻⁵.
  * With the §59 hard localisation, |A| ≤ A₀ = κΛ/S, it is ≈ 0.13κΛτ (τ = a/S ≈ (1−c)/4).
    The α-contour can be lifted (Im α = Y, damping e^{−SY}) only where the outgoing growth e^{νY/|A|} is beaten, i.e. |A| ≳ ν/S, for every
    frequency present (ν up to Λ).
  * Splitting h into frequency bands with their own A₀ ~ ν_I/S does not help: each band's diagonal mass is ~∫|h_I|/ν_I, so the total is
    ~(4κτ/π)∫|h| with ∫|h| = 2∫|ψ₂| ≥ 2|g₀| = 41.6.
* **The diagonal norm is therefore ≳ τ∫|h| ~ (1−c)·O(10) in any variant, against H/w³ ~ F_cert(1−c)^{−3}.** The ratio scales as
  (1−c)⁴/F_cert. It is ~10⁵–10⁶ at L ≲ 0.5, and drops to O(10²) only around L ≈ 3. (Measured values in §3.2.)
* For b = 4 the factor |ξ|⁸ removes the diagonal. But 𝒞 ≈ a(N/n)⁸ ≈ a c^{−4} then amplifies every contour region with |ξ|² ≫ c by (|ξ|²/c)⁴.
  This is the I_U ≈ iI_T cancellation seen at small states (§1.2): it is intrinsic at small c, not an artefact of small n.

**(ii) Large β (z → 1) on I_U.** On I_U, Im γ = −ε, so the carrier gives no damping. β = log(tan(C/2)/tan(A/2)) reaches log(2S) at A ~ ε,
and (1−ξ²)^{−p} = cosh^{2p}(β/2) ~ e^{pβ}.
* Per frequency, P⁺_ν decays like e^{−ν Im β}, with Im β ≈ atan(ε/A) − atan(ε/C) > 0. This kills A ≲ νε, so the amplitude stays at the
  physical scale (S/ν)^{p−½}… only for ν ≳ a.
* For low frequencies the per-frequency amplification on the author's contour (fixed tail at γ₁ = π/2) is ~(πa/(2ν))^{p−½}. Small states:
  ∫|W|/|𝒦| = 20–900 at ν = 1–12, n ≤ 10 (§1.2).
* The integrated amplitude's low-frequency content does not decay in Im β: it is the p = 0 endpoint, power-law ψ′(0)/(Λβ)². So
  |𝒫⁺_h(β ≈ 8)| ≈ 2·10⁻⁷ multiplies e^{3β}/64 ≈ 10⁹.
* A per-frequency-optimal tail sits at the soft saddle, γ₁(ν) = arcsin(ν/a). No single contour serves all ν; only band-dependent contours can.

### 3.2 Measurements at open-region states [DV]
* `fp4_abs.py`: raw ∫|W| on the author's contour (I_U, I_T, I_L), with the actual |𝒫⁺_h| and with a Lemma-A⁺-type majorant
  (1/π)∫|K|(|W_S| + |W_C| + e^{−Λ Im t}|W_E e^{−iΛt}|)du. The majorant is 1–3× the actual amplitude.
* `fp4_soft.py`: **soft-signed / hard-absolute** norm. For each α the whole γ path (arc + tail, so the corner is paired automatically, or
  the lower ray) is integrated SIGNED with half-wavelength panels; then
  N_soft = (𝒞/π)∫|e^{iSα}Ψ(α)/sin α|dA.
  This is the norm an α-cell certificate pays if the soft integral is evaluated exactly per α (the author's IBP/corner-pairing).
* Weight: ψ₂ (the main part G_A^ψ of H/w³), p = 3.
* Reference: |G(c)| = |ℒ(G^ψ)| to ~1e-6 (FULLPOS3 §3.1). Note w³|G(c)| ≈ N_s at these states (ratio 0.99 at (2000, 1500)).

**(a, J, b) = (2000, 1500, 0); S = 5001, c = 0.184, L = 0.20, |G(c)| = 2.88e-5, N_s = 1.59e-5.**
| \|A\|S band | [0,1) | [1,10) | [10,100) | [100,1000) | ≥1000 | total | total/\|G(c)\| |
|---|---|---|---|---|---|---|---|
| raw ∫\|W\| (I_U + I_T + I_L) | 5.0e5 | 1.4e5 | 1.1e3 | 26 | 137 | 6.45e5 | **2.2e10** |
| soft-signed N_soft | 0.025 | 0.21 | 2.2 | 18.7 | 54.2 | 75.4 | **2.6e6** |

(The first raw run had an overflow in I_T, e^{Λ|Im β|} formed before the carrier damping. It was fixed by carrying e^{iΛβ} in the exponent; the
I_U and I_L values were unaffected: I_U = 6.44e5, I_L = 25.4. In the band split, 139 of the raw norm sits at |β| < 0.01, 20 at
0.01–0.1, and 1.5 at 0.1–0.3.)

Reading:
1. **The raw norm is useless** (2e10). It is dominated by |A|S < 10: I_U at large β, plus the I_U/I_T corner. The soft-signed norm shows that this
   is an exact cancellation along each soft path (0.23 left of 6.4e5). The author's IBP/corner pairing is therefore essential, not optional.
2. **What survives is the diagonal-endpoint region |A|S ≳ 100.**
   * For |A| ≫ 1/a the γ path never reaches the physical region (γ ≥ A ≫ 1/a). Ψ(α) is then the diagonal boundary expansion: log-singular,
     amplitude ~ M_h log/(ia), with M_h ~ ∫|h|.
   * Hard localisation (lifting α where the outgoing growth e^{νY/|A|} < e^{SY}, |A| ≳ Λ/S) removes the ≥ 1000 band, leaving ≈ 21, i.e.
     **≈ 7e5 × |G(c)|**.
   * Separating the diagonal boundary terms analytically and lifting them in α (they are frequency-free) would remove more. But their total is
     the small-z (z → 0) information of G, the same information whose ladder sum Σ U_kg_k cancels from |U_kg_k| ≫ 1 down to 3e-5. So at least
     one of the two pieces has norm ≳ O(G(0)) = O(20). The effective ratio stays ≳ 1e5–1e6 at L ≈ 0.2.
   * **Superseded in part by the author's §74 extraction (§4).** The diagonal block is a REAL one-dimensional integral against the closed-form
     𝔥_B. Its norm (O(10²–10³)·N_s) costs nothing in a certificate, because it is evaluated to ~1e-20 by Arb quadrature. The decisive
     number is then the norm of the remaining contour (§4).

**(a, J, b) = (998, 20000, 0); S = 40999, c = 0.907, L = 2.38, |G(c)| = 6.21e-3** (N_s not computed: `hyb_fast` needs (2J)×J arrays).
| \|A\|S band | [0,1) | [1,10) | [10,100) | [100,1000) | ≥1000 | total | total/\|G(c)\| |
|---|---|---|---|---|---|---|---|
| raw ∫\|W\| | 5.9e7 | 1.7e7 | 1.3e5 | 451 | 125 | 7.6e7 | 1.2e10 |
| soft-signed | 0.011 | 0.017 | 0.20 | 1.86 | 80.0 | 82.1 | 1.3e4 |

* Hard localisation at |A| ≤ Λ/S removes the ≥1000 band: **2.1, i.e. 3.4e2 × |G(c)|.**
* The diagonal/hard structure scales as predicted: the |A|S < 1000 soft norm falls from 21 (L = 0.20) to 2.1 (L = 2.38), while |G(c)| grows
  by 200.

---

## 4. The author's §74 diagonal extraction: implemented and measured [DV]
`fp4_diag.py`. Soft variable ξ: Y = (X−ξ)/(1−ξX), X = e^{ε+iA}, dγ/sin α = 2dξ/((1−ξX)(1−ξ/X)), and the carrier is X^B((1−ξX)/(1−ξ/X))^a. With
K = X^B(1−ξX)^{a−1}(1−ξ/X)^{−a−1}:
* ℒ(G) = (𝒞/π)∫dA∫Kξ^{s−1}F dξ (full F), equivalently 2Re(𝒞/π)∫dA∫Kξ^{s−1}F⁺dξ (outgoing).
* Each soft path is cut at its LAST crossing of |ξ| = ξ₀. The start becomes the real segment [0, ξ₀] plus an arc on |ξ| = ξ₀.
* ℒ = D + R:
  * D = (𝒞/π)∫₀^{ξ₀}ξ^{s−1}G(ξ²)𝔥_B(ξ)dξ, using the REAL G;
  * R = 2Re(𝒞/π)∫dA[∫_arc Kξ^{s−1}F⁺dξ + ∫_{rest}W⁺dγ/2·…], i.e. the §55 pieces with the diagonal prefix replaced by the arc.

**Checks.**
* 𝔥_B by direct A-quadrature is real to 3e-12, and agrees with the author's 41-term closed form to 2.8e-12 (ξ ∈ {0.002, 0.01, 0.02}
  at (2000, 1500); 5e-15 at the toy state).
* Toy (λ = 8, n = 10, N = 30.5): D + R = 8.665245954894 vs exact 8.665245954376 (**6e-11**).
* (2000, 1500, 0): D = −3.05e-6 = −0.106·|G(c)| = −0.104·N_s/w³. The author's diagnostic is −1.66e-6 ≈ −0.10·N_s at his state.

**Norms (soft-signed per A, as in §3.2; N_s/w³ ≈ |G(c)| at (2000, 1500)).**
| state | L | \|G(c)\| | D | D-norm (real 1D) | N_rem (soft) | N_rem/\|G(c)\| | bands of N_rem, \|A\|S: <10 / [10,100) / [100,1000) / ≥1000 |
|---|---|---|---|---|---|---|---|
| (2000, 1500, 0), ξ₀ = 0.02 | 0.20 | 2.88e-5 | −3.05e-6 | 2.46e-2 | 0.457 | 1.6e4 | 0.013 / 0.075 / 0.171 / 0.197 |
| (2000, 1500, 0), ξ₀ = 0.05 | 0.20 | 2.88e-5 | −1.92e-7 | 2.46e-2 | 0.0846 | 2.9e3 | 0.005 / 0.019 / 0.029 / 0.031 |
| (2000, 1500, 0), ξ₀ = 0.1 | 0.20 | 2.88e-5 | +2.29e-7 | 2.46e-2 | 0.0228 | 7.9e2 | 0.0023 / 0.0059 / 0.0072 / 0.0074 |
| (2000, 1500, 0), ξ₀ = 0.2 | 0.20 | 2.88e-5 | +1.70e-6 | 2.47e-2 | 0.00571 | **2.0e2** | 0.0010 / 0.0016 / 0.0017 / 0.0015 |
| (2000, 1500, 0), ξ₀ = 0.35 | 0.20 | 2.88e-5 | +1.51e-6 | 2.46e-2 | 0.00175 | **61** | 0.00044 / 0.00048 / 0.00047 / 0.00035 |
| (998, 20000, 0), ξ₀ = 0.02 | 2.38 | 6.21e-3 | +2.8e-7 | 1.3e-3 | 0.599 | 97 | 0.015 / 0.012 / 0.077 / 0.496 |

* **The remaining norm falls like ξ₀^{−2}** (the outgoing amplitude at |β| ≥ 2ξ₀ is ~ψ′(0)/(Λβ)²). Meanwhile the diagonal block's
  1D norm stays ≈ 855|G(c)|, dominated by small ξ where G ≈ G(0). D itself moves with ξ₀ (it is the real-segment share), and only D + R is
  invariant; the toy checks of that invariance at larger ξ₀ are in §4.1.
* Hard localisation (|A| ≤ Λ/S) would remove the ≥ 1000 band: at ξ₀ = 0.02 that leaves ≈ 0.26 (9e3|G(c)|) at L = 0.2 and ≈ 0.10
  (**16|G(c)|**) at L = 2.38.
* The diagonal block's norm (850 × |G(c)|) costs only precision: it is a 1D real integral of explicit functions.
* The **raw** norm of the remaining contour is unchanged (6.45e5): the soft-signed evaluation stays essential.

### 4.1 Invariance of D + R under ξ₀, and a strong-cancellation check [DV]
* Toy λ = 8, (n, N) = (10, 30.5): D + R = 8.665245954894 at ξ₀ = 0.02 and 8.665245954894 at ξ₀ = 0.3 (exact 8.665245954376; rel. err 6e-11).
  * D moves from 0.381 to 0.628; N_rem/|ℒ| from 7.4 to 2.2.
  * The author's 41-term closed form for 𝔥_B is an asymptotic series in 1/B: it fails at B = 21, ξ = 0.3 (0.8 off), so the direct
    A-quadrature is used.
* Toy λ = 30, (n, N) = (40, 85.5), a = 45.5: the exact ladder sum is 0.0990087008 with max|U_kM_k| = 3.7e6 (cancellation 4e7).
  * At ξ₀ = 0.3, D + R = 0.0990086857, **rel. err 1.5e-7** on the coarse α-grid (`fp4_toy2.py`).
  * N_rem/|ℒ| = 19, D-norm/|ℒ| = 7.2.
* The signed total at the physical states was not evaluated. It needs ~3000 α-nodes × ~4a soft nodes ≈ 1–2e7 amplitude evaluations,
  i.e. 4–7 h per state on one core.

---

## 5. The author's §78 quadrature-free reference: benchmark and cell test

### 5.1 Pointwise benchmark [ball arithmetic; DV comparison]
`fp4_ref.py` (python-flint). 𝒦_ref = Σ_{j≤3M}W_j[t^j]Y, with the R_{m,j} recurrence and the ODE jet exactly as in §78. The exact comparison is
w³Σ_{k≤n}U_kg_k, with midpoint propagation in the g_k recurrence; ball wrapping there grows like (1+√2)^k. p = 3, b = 0, M = 1000:
| (J, a) | bits | max\|R_{m,j}\| | ν = 20 | 300 | 700 | 950 | 993 |
|---|---|---|---|---|---|---|---|
| (1001, 998) | 8000 | 1e571 | 0 (1600 b) | 0 | 0 | 1.4e-271 | 6.6e-237 |
| (1500, 2000) | 1600 | 1e187 | 0 | 0 | 2e-185 | 4e-77 | — |
| (20000, 998) | 8000 | ~1 | 0 | 0 | — | 2.2e-276 | — |

(Entries are |w³𝒦_ref − w³𝒦_exact|; 0 means agreement to all printed and internal digits.)
* The reference is exact for practical purposes at all three states, including the turning point ν = 993 at a = 998.
* Cost: W_j takes 5 s at 1600 bits and 40 s at 8000 bits per state. Each frequency then costs 0.04–0.4 s.
* At 1600 bits the rigorous radii blow up (1e65 at ν = 300, 1e686 at ν = 950) while the midpoints are right; ~8000 bits gives tight balls.
* **Integrated (`fp4_ref_int.py`)** at (1500, 2000, 0), 1600 GL nodes in p, 4000 bits:
  * H_ref = 2∫ψ₂w³𝒦_ref = 1.566917662696e-5 = H_direct (same nodes) to all double digits;
  * w³G(c) = 1.566918295e-5 (rel. 4e-7);
  * 2∫|ψ₂||w³𝒦| = 0.66 (the 4e4 ν-cancellation).

### 5.2 Uniformity over the continuum: naive interval cells fail [ball arithmetic]
`fp4_ref_cell.py`: a ∈ [a₀, a₀+δ] at fixed integer J, all of c, d, h, w as balls, ν = 300, M = 1000, 8000 bits.
| J | δ = 1e-1000 | 1e-600 | 1e-300 | 1e-100 | 1e-30 | 1e-3 |
|---|---|---|---|---|---|---|
| 1001 | tight | rad 4e-58 | **rad 4e242** | | | |
| 20000 | | | rad 1e379 | 1e579 | 1e664 | 1e718 |

* **The usable twist-cell width is δ ~ 1e-550.** At J = 20000 the true cancellation is mild (max|R| ≈ 1), yet the wrapping in the M-step
  R recurrence (≈ 2^M) still kills every cell.
* Taylor models are not implemented. Their remainder terms still pay Σ_j|W_jy_j| (up to 1e187–1e571) × (δσ)^{K+1}. At K = 100 and
  sensitivity σ ~ 1 per unit a that allows δ ≲ 1e-6, i.e. ~1e9 a-cells for each J.
* J cannot be enumerated: n ≈ 2ae^L gives J up to ~1e7 at L = 7 and ~1e12 at L ≈ 20. It must be a continuous cell parameter too.
* **Verdict on §78 [DV + ball tests]:** an excellent pointwise reference, but no route to a continuum certificate without an a priori
  uniform derivative bound on H_ref in (a, 1/J). Such a bound is the per-order-majorant problem in another form.

---

## 6. Obstructions recorded
1. **Raw ∫|integrand| of the frequency-integrated §55 contour is ~1e10 × N_s at both test states.** The integrated outgoing amplitude has
   power-law low-frequency content: the odd-j p = 0 terms ψ^{(j)}(0)/(iΛt)^{j+1} survive in P⁺ alone. That content has no e^{−ν Im β} decay,
   so I_U at large β (z → 1, (1−ξ²)^{−p} ~ e^{pβ}) and the I_U/I_T corner carry huge mass. It cancels exactly along each soft path.
   **Rule: never bound the soft integral absolutely; evaluate it signed per α.**
2. **Diagonal γ = α (z = 0):** a forced boundary of every §53–§59 contour (the image of u = v = 0). There |𝒫⁺_h| ~ ∫|h| ~ 20 (= G(0)/2 +
   (i/π)M_h log 2β) for β ≲ 1/Λ, and ~ψ′(0)/(Λβ)² up to β ~ 1.
   * Soft-signed norm before extraction: 2.6e6 at L = 0.2. It is removed by the author's §74 real-segment extraction.
   * The leading log expansion (i/π)[M_h log 2β + D_h] (M_h = 20.77) holds only for β ≲ 1e-4. The block must be the full real-G integral,
     not a few constants.
3. **Per-frequency vs integrated geometry.** The author's fixed tail γ₁ = π/2 is optimal only for ν ≈ a; low frequencies would want
   γ₁(ν) = arcsin(ν/a). This is irrelevant once the soft integral is signed.
4. **Ball wrapping** in long three-term recurrences (g_k ladder: (1+√2)^k; R_{m,j}: ~2^M) inflates radii by 1e300–1e3000 regardless of the
   true cancellation. Point evaluations need ~8000 bits; parameter cells are infeasible without Taylor models.
5. **Cost / uniformity (the decisive one).** See §0 item 5.

## 7. Files and reproduction (`ads4/code/fullpos4/`; `ulimit -v 4000000; export OMP_NUM_THREADS=1`)
| file | content | grade |
|---|---|---|
| `fp4_kern.py` | exact ladder kernel; §59 form with mpmath (slow) | DV |
| `fp4_pf.py` | (α, γ) contour, base §59 path and own outgoing deformed path; accurate P⁺ (`_h2f1_zb`) | DV |
| `fp4_s55.py`, `fp4_s55_sweep.sh` | author's §55 recipe; 12-state sweep → `cache/fp4_s55_sweep.txt` | DV |
| `fp4_int.py` | frequency-integrated 𝒫⁺_h (Euler/Q form, W_S/W_E/W_C split, rotated endpoint ray, factored e^{iΛβ}), Lemma-A⁺-type majorant | DV |
| `fp4_abs.py` | raw ∫\|W\| (actual and majorant) on the author's contour | DV |
| `fp4_soft.py` | soft-signed / hard-absolute norm | DV |
| `fp4_diag.py` | §74 diagonal extraction + remaining contour (D, R, norms); `fp4_toy.py`, `fp4_toy2.py` exact end-to-end checks | DV |
| `fp4_ref.py`, `fp4_ref_cell.py`, `fp4_ref_int.py` | §78 reference (python-flint balls), cell enclosure test, integrated comparison | DV + balls |
| `fp4_lemU.py` | Lemma U twist polynomial at fixed real c, both channels, Chebyshev test, M-convergence → `cache/fp4_lemU.jsonl` | DV |
| `fp4_lemU_c.py` | the same at complex u (Lemma C strip test) → `cache/fp4_lemU_c.jsonl` | DV |

Results: `cache/fp4_abs*.jsonl`, `fp4_soft.jsonl`, `fp4_diag.jsonl`, `fp4_toy2.jsonl`, `fp4_ref*.jsonl`, `*.log`.
```
python3 fp4_s55.py 2 7.5 0 3 1.2 32 40          # §55 benchmark (I_U, I_T, I_L)
python3 fp4_int.py                              # amplitude self-tests
python3 fp4_diag.py 2000 1500 0 4 - 0.35         # diagonal extraction + remaining norm (~15 min)
python3 fp4_toy2.py 30 40 85.5 0 0.3 16          # strong-cancellation end-to-end check (~5 min)
python3 fp4_ref.py 1001 998 0 300,950,993 1000 8000
python3 fp4_ref_cell.py 1001 998 0 300 1000 8000 1e-600,1e-300
python3 fp4_ref_int.py 1500 2000 0 1600 1000 4000
python3 fp4_lemU.py 2.0 1000 8000 1500           # twist test at L = 2, both channels (~11 min)
python3 fp4_lemU_c.py 3.314908908306154 0.02 500 8000 1500   # complex-u strip test
```

---

## 8. Lemma U (author, ANGULAR_ROUTE §80): the twist direction as one polynomial per c [DV]
`fp4_lemU.py`. At fixed c (u = β_c), every state with a ≥ 998 has x = h/d² = y·x_* with y ∈ (0, 1], and R_{m,j}(c, x) = y^mR_{m,j}(c, x_*).
So H_M = P(c, y)/D⁴ with
* P = ΣA_my^m + d²x_*yΣ(q_bA_m + B_m)y^m and D = 1 + q_bd²x_*y;
* A_m = Σ_jR_{m,j}(c, x_*)Y^A_j, Y_j = d³(cd)^jG^{(j)}(c)/j!.

Implementation:
* Jets Y^A (ψ₂, p = 3) and Y^D (Λ²ψ₃, p = 4) of the integrated G: per-frequency ODE jets, Gauss–Legendre in p.
  N_p must resolve e^{iΛpβ_c} and degree-2j polynomials: 8000 at L = 2, 15000 at L = 7.
* Y^B = 𝒟_bY^A + σ_bY^D by the exact jet contraction
  [t^j]𝒟_bg = [(j+2)(j+1)g_{j+2} + ((d−c)j + (7−4b)c + 1)(j+1)g_{j+1} + (−cdj(j−1) + (7−4b)cdj)g_j]/(2cd) + (9/2 − b(b−1))g_j.
* Chebyshev conversion on y ∈ [0, 1]: P̂_ℓ = (2−δ_{ℓ0})Σ_{m≥ℓ}C(2m, m−ℓ)4^{−m}p_m. Test: P̂₀ − Σ_{ℓ≥1}|P̂_ℓ| > 0.
* Midpoint arb arithmetic (1500–3000 bits), M = 1000. IDEAL NF02 weights (no EPSW residual).

**Sanity.** A₀ = d³G(c) = 2·OmLR(β_c) (Mehler) to 1e-9 at L = 2 (6.067386e-6).

**Convergence in M (L = 2, b = 0).**
* |A_m| = 3e-5 (m = 10), 9e11 (40), 4e17 (100), 1e-3 (200), 3e-38 (300), 7e-78 (400), 2e-253 (1000).
* The partial sums P_{M′}(1) are wild up to M′ ≈ 200 and agree to all digits from M′ = 300 on. The truncated h-series at a = 998 really does
  converge, after a 1e18 transient. (The truncated per-order bounds of FULLPOS3 saw only the transient.)

**Twist-positivity test, all a ≥ 998 at fixed c:**
| L (c) | b | A₀ = H(a = ∞) | H_M(a = 998) | \|P̂₁\| | \|P̂₂\| | noise floor of \|P̂_ℓ\| | margin (P̂₀ − Σ\|P̂_ℓ\|)/P̂₀ |
|---|---|---|---|---|---|---|---|
| 0.1 (0.0952) | 0 | 1.8895e-5 | 1.9968e-5 | 5.4e-7 | 2.3e-13 | ~1e-13; < 1e-33 (ℓ ≥ 100) | **0.97240 PASS** |
| 0.1 (0.0952) | 4 | 1.8895e-5 | 1.7822e-5 | 5.4e-7 | 2.8e-13 | ~1e-13 | **0.97078 PASS** |
| 0.5 (0.3935) | 0 | 1.2919e-5 | 1.3309e-5 | 1.9e-7 | 1.6e-13 | ~3e-16; < 2e-30 (ℓ ≥ 100) | **0.98515 PASS** |
| 0.5 (0.3935) | 4 | 1.2919e-5 | 1.2529e-5 | 1.9e-7 | 4.4e-13 | ~5e-16 | **0.98469 PASS** |
| 2 (0.8647) | 0 | 6.0674e-6 | 6.0969e-6 | 1.5e-8 | 2.3e-14 | ~1e-17 (ℓ = 5–50); < 4e-32 (ℓ ≥ 100) | **0.99757 PASS** |
| 2 (0.8647) | 4 | 6.0674e-6 | 6.0378e-6 | 1.5e-8 | 3.0e-14 | ~1e-17 | **0.99756 PASS** |
| 7 (0.99909) | 0 | 3.4905e-7 | 3.4906e-7 | 4.8e-12 | 1.0e-17 | ~1e-18; < 5e-35 (ℓ ≥ 100) | **0.999986 PASS** |
| 7 (0.99909) | 4 | 3.4905e-7 | 3.4904e-7 | 5.9e-12 | 8.3e-18 | ~1e-18 | **0.999983 PASS** |

(F_cert(β_c = 0.64) = 1.54e-5; F_cert(1.47) = 1.12e-5; F_cert(3.31) = 4.86e-6; F_cert(8.39) = 2.66e-7. At L = 0.5: max|R_{m,j}| = 1e213, |A_m| peaks at 2e18,
constant partial sums from M′ = 300, 3000 bits, N_p = 6000, 17 min.
At L = 0.1: max|R| = 1e1204, A_m transient ≤ 1e11, 3000 bits, N_p = 5000, 10 min. The P̂ noise floor (1e-13 ≈ 5e-9 relative) is
higher there, consistent with 3000 bits against |R| ~ 1e1204; lower L needs ≥ 6000 bits.)
* At L = 7, |A_m| peaks at 1e14 (m ≈ 100), and the partial sums are constant from M′ = 300 on.
* Cost: 11 min per c at L = 2 (N_p = 8000) and 28 min at L = 7 (N_p = 15000), on one core; both channels share the jets.

### 8.1 Complex impact parameter (Lemma C test) [DV]
`fp4_lemU_c.py`: the same construction in complex arithmetic at u = β_c + iη. Here c = tanh²(u/2), x_* is analytically continued, and the seeds are
P_ν(cosh u), ∂_βP_ν. Reported:
* max_{y∈[0,1]}|H_M(u, y)|/F_cert(Re u);
* max_y|H_M(u, y) − H_M(u, 0)|/F_cert (401 y-points, M = 500).

(The first L = 2 runs ran at 256 bits, because `import fp3_htdata` resets ctx.prec. Midpoint arithmetic at 256 bits still covers the 1e18
transient; a 1500-bit recheck is queued.)
| L | u | η | b | \|H_M(u, 0)\|/F_cert | max_y\|H_M(u, y)\|/F_cert | max_y\|H_M(u, y) − H_M(u, 0)\|/F_cert |
|---|---|---|---|---|---|---|
| 0.5 | 1.4738 | 0.02 | 0 | 1.216 | **1.275** | 0.154 |
| 0.5 | 1.4738 | 0.02 | 4 | 1.216 | **1.397** | 0.455 |
| 2 | 3.3149 | 0.02 | 0 | 1.249 | **1.258** | 0.0090 |
| 2 | 3.3149 | 0.02 | 4 | 1.249 | **1.249** | 0.0083 |
| 2 | 3.3149 | 0.01 | 0 | 1.249 | 1.255 | 0.0061 |
| 2 | 3.3149 | 0.01 | 4 | 1.249 | 1.249 | 0.0061 |

The author's own measurement (§87, coordinator's message) is at low L: c = 0.037 + 0.0037i (u ≈ 0.39 + 0.02i), y = 1/20, b = 0, with
|H_M(u, y) − H_M(u, 0)| > 463F_cert. At L = 2 the strip transfer is ~1e-2 F_cert. So the strip bound is strongly L-dependent; low-L points
are queued (`cache/lemUc_L0.04.log`).

### 8.2 Cost of the c-direction (estimate)
* **Per c (twist test, DV):** 11 min (L = 2) to 28 min (L = 7). It is dominated by the ν-quadrature for the jets, with N_p ≈ 3000 + 1.4Λβ_c
  (up to ~4e4 at L = 30). A certified version needs rigorous jet enclosures (GL error bound on a Bernstein ellipse in p, plus ≥ 8000 bits
  against the (1+√2)^{3M} wrapping of the jet recurrence): perhaps 3–10× the cost.
* **Cells in u from Cauchy on the strip |Im u| ≤ η.** With K := sup_strip max_y|H_M|/F_cert, a second-order Taylor step using the exact
  first derivative (derivative contractions, ~2× cost) gives cells of half-width δ ≈ η·√(2·margin/K).
  * At L ≈ 2 (K ≈ 1.3, H/F ≈ 1.25): δ ~ η, i.e. ~25 points per unit u.
  * At low L (K ≳ 463 by the author's point): δ ~ 0.05η, i.e. ~1000 points per unit u.
  * Rough total: (31 units × 25) + (0.6 units × 1000) ≈ 1.4e3 c-points × 15–30 min ≈ 350–700 core-h (DV). The certified version is several
    times more.
* **The author's design:** one degree-6000 jet per Δβ = 0.05 cell (~600 cells), with the grid points inside a cell from the Taylor shift.
  That moves the jet cost to ~600 × ~1 h, plus cheap R recursions (~1 min each at M = 1000). Also several hundred core-hours.
* **Verdict:** not an overnight job on 22 cores (≈ 176 core-h), but within a few nights, or one night on ~100 cores, if
  1. the jets are certified;
  2. the realised (EPSW-repaired) weight replaces the ideal NF02 one (a residual at ν ≈ 1000–1130 > a needs the same jets with more nodes);
  3. the §78 truncation error E_ref is reproduced in Arb.
