# FULLPOS4: the outgoing-split, frequency-integrated §55 contour at Δ_g = 1000

*26 Nov 2026. Research agent. Unreviewed, not committed. No existing file was edited; nothing under `prl/` was touched.*
*New code: `ads4/code/fullpos4/` (outputs in `fullpos4/cache/`, whose `.gitignore` ignores everything but itself).*
*Compute: ≤ 2 processes, `ulimit -v 4000000`, `OMP_NUM_THREADS=1`. The two FULLPOS3 scan processes were not touched.*

Grades: **[certified]** Arb balls or exact rationals with every remainder proved; **[P]** proved on paper here;
**[DV]** numerical; **[conj]** conjecture. Author/coordinator claims are marked as such and kept separate.

*(Running document; sections marked "in progress" are partial.)*

---

## 0. Verdict (in progress)

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
