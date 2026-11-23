# LARGESPIN_L: Route L (large impact parameter / lightcone) to the large-spin region of Ω^EE at Δ_g = 1000 — what is proved, a conditional theorem with explicit J\*, and where the route stalls

*23 Nov 2026. Unreviewed, not committed by this agent. No existing file was edited; `LARGESPIN_W.md` and `code/lsw/` were not touched
(read once, after the results below were fixed, for the cross-check in §7.4). Note: Route W's commits `64bbf24`/`c1ae5ac` swept early
snapshots of `code/lsl/lsl_core.py` and `lsl_frame.py` into git; the working tree holds the final versions.*
* *New code:* `ads4/code/lsl/` — `lsl_core.py` (NF02 weights from `prl/data/nf02_psi.json`, closed-form `Φ_k`, `Ω_k(β)`, `M₁`, `P_{−½+iν}`),
  `lsl_frame.py` (gauge identity and the three classical integrals, §2), `lsl_lrinf.py` + `lsl_lrinf_grid.py` (Lemma LR∞, §3),
  `lsl_hybcmp.py` (hybrid vs leading Regge, J\*(m)), `lsl_smear.py` (Bernstein smearing, lightcone limit of CHMRS (3.20), §2.6, §4.1),
  `lsl_window.py` (window tolerance map J = 26 … 10⁶, de-shift vs raw, §4.2), `lsl_ringtol.py` (large-β tolerance profile, §4.3),
  `run_all.sh`. Caches in `code/lsl/cache/` (its `.gitignore` ignores everything except itself).
* *Read-only imports:* `n2sub/n2_common.py` (X grid, R∓, Pnu), `largespin/ls_core.py` (OmLR, hyb_fast). No bytecode was written.
* *Literature:* CHMRS arXiv:2106.10274, LaTeX source (`sec_heavy.tex`: eqs. `functional action`, `neat B integral`, `heavy limit I`,
  `Fourier identity`, `moment maintext`, `regge subleading` = (3.20), and the two bullet properties after it); classical integrals [Q]:
  Laplace's second integral for `P_J`, Euler's Γ integral, the Riesz/Gindikin Laplace transform on the Lorentz cone, the Paley–Wiener
  theorem for the H₂ spherical (Mehler–Fock) transform (Koornwinder; Helgason).
* *Compute:* ≈ 0.5 CPU-h, ≤ 2 processes, `ulimit -v 3000000`, `OMP_NUM_THREADS=1`.

Grades as UNIFIED_THEOREM_v3 §0. `β`: H₂ impact parameter, `sinh(β/2) = (J+½)/m`; `b̂ = 2Λ sinh(β/2) = (2J+1)Λ/m`; `u = Λ²/m²`;
`τ = Δ − J ≥ Δ_g`; `F_LR(β, u, σ) := 2Ω₂(β) + uσΩ₃(β)`, `Ω_k(β) := ∫₀¹ψ_k(p)P_{−½+iΛp}(cosh β)dp` (the leading-Regge per-state action;
σ = +1 / −1 for `++` / `+−`). "De-shift form" = N2's R1: a raw-row correction `𝒜[f] = g₊f(ν+i) + g₋f(ν−i)` with `f` strip-analytic, whose
k = 2 cost is exactly the ψ₂ cost of `f`. "Relative" = relative to the local amplitude of the leading-Regge image.

---

> **Update (PW_IMAGES.md, unreviewed): R1 is replaced by an exact identity plus one bounded monopole term. PW-R1 and PW-Loc are false as literally stated, but neither is needed any more. Theorem W keeps J\* = 50 (100 at ×2) WITHOUT P3. Theorem L keeps J\* = 857 with P_L reduced to magnitudes.**
> * **Exact de-shift identity for the EE kernel.** Every b-shifted scalar block in EE = Σ_b satisfies CHMRS's magic-block identity with b-independent coefficients:
>   2η(1+r²)Ĝ¹_b(ν) = (1+i/2ν)Ĝ⁰_b(ν+i) + (1−i/2ν)Ĝ⁰_b(ν−i).
>   It is [P] to r⁸, [DV] to r⁴⁰ at 50 digits, [D] at all orders. So every plane-wave image, in any gauge and with any insertion, is exactly 𝒜_g[F] (CHMRS's scalar de-shift), with F even and analytic for |Im μ| ≤ 3/2. End to end against ee_quad: 2·10⁻⁷–1.1·10⁻⁶. **Route W's stripped (R2) worst case cannot occur.**
> * **PW-R1 is refuted** with respect to the EE operator 𝒜 = 𝒜_g∘c. It needs F(±i/2) = 0, which O(1/m²) radial insertions violate [P]: |y|² gives a ratio of 0.914, r gives 0.500. Whether the full ⟨TTTT⟩ sum cancels this is [I].
> * **The failure is exactly a monopole term, and it is bounded [DV].** The exact cost is ∫ψ₂F/c + D·F(i/2), with D = −9.98·10⁻⁴ (1000) / −1.10·10⁻² (300), ≈ −990/Λ². D is constant in impact parameter, and the equivalent cost kernel stays below 10⁻³ (no e^{t/2} growth).
> * **PW-Loc is false at O(1/m²):** index-shifted harmonics carry density 0.271 per unit coefficient at t < 1 (closed form [P]). It is no longer needed, since every finite-order structure now has an exact cost.
> * **Theorem W without P3 [P given P0, P1, P2, P2′, P4] + [DV]:** the scalar-shaped worst case is tolerated at ×24–29 for J ≥ 150 (stripped form gave ×1.8–7). J\* = 50 at (3, 10, 10); J\* = 100 at ×2. P2′ (the radial-shift family exhausts the O(1/m²) multipliers) is [D-sketch].
> * **Theorem L [D given P_L′]:** J\* = 857 at the gap. The monopole tolerance gives a tolerated Λ_amp of 1240 (J = 857), 555 (10⁴), 166 (10⁵) and 43 (10⁶).
> * **Still needed by both:** the spinning O(1/m²) magnitudes (Λ_amp ≈ 2–3 [NE] for J ≤ 100, [M] above) and the beyond-all-orders premise. Untested: polarisation-tilt O(1/m²) terms.

> **Update (INDEP_CHECK2.md, independent reproduction with new code): all five targets reproduced, nothing load-bearing refuted. Two substantive qualifications: Theorem W's J\* = 50 → 64 once the k = 3 piece is included, and Lemma LR∞'s r(β) table uses a cruder constant than the stated formula.**
> * **Magic-block identity:** exact with l and b symbolic to O(r⁶), so it is an identity in b. It is exact to O(r⁸) for every b = −4…4, and holds to 4·10⁻⁴⁸ at 50 digits to r⁴⁰ (including b = 0.37, 7, 10 and complex ν); controls fail at O(1). **The monopole constant D is reproduced to 8 digits, and its t-independence is now derived:** 1/c = 1/32 + (15/2)/(4μ²+1) − (9/2)/(4μ²+9), so D = (15/2)·defect[1/(4ν²+1)].
> * **Lemma LR∞:** reproduced; the true r(β) is ≤ 8·10⁻⁷ on β ∈ [1, 32]. **Discrepancy:** the stated formula (with ∫|ψ|) gives r(1) = 0.016, while the table (0.577) used Σ|c_n|/(n+1). Both are rigorous; the text must name the constant used. With the tighter constant the bound holds down to β = 0.5 (r = 0.28), which would give J\* = 324 at the gap instead of 857 (needs a certified ∫|ψ₃| and a denser grid).
> * **Order of limits (route N):** sound, and the quotes are verbatim or faithful. Minor fixes:
>   * the CHMRS ε quote is in §1, not §1.1.1;
>   * the "technical assumption" is eq. (5.3), not §5.3;
>   * CLS-D restates the logic in §2.2.2 but never discusses the order of limits, so do not say "the logic of CLS-D";
>   * f⁻ = −1.
> * **Theorem W:** the 2×2 criterion d ≤ H is correct and sharp (200 000 random tests). All J\* table entries are reproduced. The Gershgorin factor is ≤ 1.5 only if f₊₊ ≈ f₊₋ (otherwise up to 3; 1.50–1.53 on the samples). **Main qualification: d omits the k = 3 piece of the corrections. Including it moves J\* 50 → 64 at (3, 10, 10), 30 → 33, and 64 → 80**; the J = 50, 51 states fail at d/H ≈ 1.24. This must be stated as a named premise or used as the J\*. The cache stores the 4a/4b shapes with 1/m² instead of 1/m⁴ (documentation slip only).
> * **⟨TTO⟩ counts:** verified against DKKPS's TTTT.tex table and by an independent helicity count (ℓ ≤ 11).

## 0. Summary

### 0.1 Table

| # | claim | grade | explicit constants / J\* | residual |
|---|---|---|---|---|
| **L1** | **Lemma LR∞ (leading-Regge action at large impact parameter).** For k = 2, 3 and every β > 0: `Ω_k(β) = ψ_k′(0)M₁(β)/Λ² − ψ_k″(0)/(2b̂³) + E_k(β)` with an explicit bound on `|E_k|` (§3). Hence for **all β ≥ 1 (b̂ ≥ 1.0421Λ), all m ≥ Λ (u ≤ 1), both helicities**: `F_LR ≥ (1 − r(β))·A₁M₁(β)/Λ² − A₂/(2b̂³) > 0`, `A₁ ≥ 327.2`, `−A₂ = 2047.1`, with `r ≤ 0.577 (β ≥ 1), 0.070 (β ≥ 1.5), 0.015 (β ≥ 2), 1.6·10⁻³ (β ≥ 3), 1.0·10⁻⁴ (β ≥ 5)`, → 0 as β → ∞. | **[P + DV]** (inequalities proved; constants evaluated in floating point on a β grid, step 0.01 on [0.9, 3], 0.05 to 16, 0.5 to 64; analytic tail [D]) | J\*(m) = ⌈(1.0421·m − 1)/2⌉ (m in AdS units, J\* ≈ 0.52m): **857** at the gap, 987 (1.1m_gap), 1275 (1.3m_gap), 780 (1.5Λ), 2602 (5Λ) | none for the leading-Regge term |
| L1′ | Coverage of the **leading-Regge / hybrid floor**: with Lemma C′ (b̂ ≤ 119, m ≥ 1.3Λ) and ROUTEA ([91.2, 2Λ], m ≥ 1.5Λ) it is now rigorous on **every b̂ for m ≥ 1.5Λ**, and on b̂ ≥ 1.04Λ for every m ≥ Λ. Hybrid − F_LR ≤ 2·10⁻⁶ relative at J ≥ 1000 (sampled). | [P + DV] / [certified] / [NE] (hybrid vs F_LR) | — | m ∈ [Λ, 1.5Λ), b̂ ∈ (119, 1.04Λ): [NE] (scan); m ∈ [Λ, 1.3Λ), b̂ < 119: [NE] |
| **L2** | **Spin-Gamma representation (exact).** In the gauge `x₅ = 0, x₁ = −e, x₂ = ∞` the crystal's heavy factor is `B^{Δ/2}P_J(t) = (1/2π)∮dφ (−w²)^{−(τ−1)/2}(−n_φ·w)^{−(J+1)}`, `w = e + x + y`, and `(−n·w)^{−(J+1)} = ∫₀^∞ s^J e^{s n·w}ds/J!`. So **at fixed twist the per-state image is a Gamma(J+1) average over s of one J-independent, CFT-data-free function**; with the Riesz transform of `(−w²)^{−λ}` every heavy block is a **positive** superposition of plane waves `e^{p·w}`, p ∈ V₊ (Bernstein). | **[P]** for the scalar heavy factor (gauge identity [DV] to 10⁻¹⁴; integrals [Q], checked [DV]); **[D-sketch]** for ⟨TTO⟩⟨OTT⟩ (finite set of J-independent polynomial insertions, index shifts J → J − r, r ≤ 4) | concentration: sd β(p) = 0.048–0.059, ⟨\|p\|⟩/m = 0.990–0.994 (τ = 1000, J = 50–3000) | — |
| L3 | **R1 reduction.** Every correction generated by the heavy factor (the Bernstein smearing, at every order) is of de-shift form automatically. R1 can fail only through the plane-wave images of the EE kernel itself, a property of the functional, not of the CFT. | [D given L2] | — | **R1 for the EE plane-wave images [I]** |
| **N-L1** | **Negative: the exact positivity of L2 gives no sign-free control.** The Bernstein smearing alone shifts the image by 1–50 % of the local amplitude at ν = 5–30 (J = 50–3000); the full scalar correction (CHMRS (3.20)) is ≤ 2·10⁻³. The kernel-side corrections cancel the smearing to that level. The "positive light-ray" route does not bound the non-ring part by the ring part. | [DV] given CHMRS (3.20) [Q] | table §2.6 | — |
| **L4** | **Lightcone (J → ∞) uniformity.** At fixed τ the relative beyond-LR correction tends to a **J-independent lightcone amplitude** `F±(ν, τ) − 1`, not to O(τ/J) or O(1/J). For the scalar it is explicit: `[(±iν − ½)(ν² + c₁) + (ν² − a(a+1))c₂ + c₃]/(12(τ−2)²)`; its size is ≤ 1.2·10⁻⁴ on ν ≤ 10 at τ = 1000, and CHMRS (3.20) is uniform in J (`η/m² ≤ 1.5/(τ−2)²`). The J-dependent part is O(e^{−2β}) ~ (τ/J)² plus Gamma smoothing. | [Q + DV] (scalar); [D-sketch] (spinning) | F(0, 1000) − 1 = 7.75·10⁻⁶; table §4.1 (J = 26 … 10⁶) | the spinning lightcone amplitude [I] |
| **L5** | **Large-J window map (relative tolerances, sign-free, ν < ν_W, κ = ½):** de-shift form: ε_ds(ν<10) has its **minimum 6.7·10⁻⁴ at J ≈ 350–500** (gap) and a plateau 1.3–1.6·10⁻³ for J ≥ 3000; ε_ds(ν<20) min 3.0·10⁻⁴, plateau 5–6·10⁻⁴. **Raw (non-de-shift) form: ε_raw(ν<10) ≈ 2.7·10⁻⁴ at J = 150–500, then ≈ 0.7/J** (7·10⁻⁵ at 10⁴, 6.8·10⁻⁷ at 10⁶). | [DV] (sampled J; F_LR from L1) | table §4.2 | R1 is decisive from J ≈ 150 on, with required accuracy growing ∝ J |
| L6 | **Large-β tolerance profile and an H₂ Paley–Wiener locality lemma.** At β ≥ 1 the Lemma-I tolerance flattens: density at t ≥ 1 is tolerated at 0.04–0.9 of the ring weight (J = 857 … 10⁶), ring-located corrections (t ≥ 0.97β) at 0.67–0.85; only t ≲ 0.1 (b̂ ≲ 100) is penalised (×10²–10⁷). Lemma W: a correction `G(ν)·P_{−½+iν}(cosh β)` with G the spherical transform of a radial measure supported in B(T) costs ≤ ‖k‖·max_{\|t−β\|≤T}\|Ω_LR(t)\|. | Lemma W **[P given PW theorem [Q]]**; profile [DV] | table §4.3 | — |
| **TL** | **Theorem L (conditional).** For every heavy primary with **b̂ ≥ 1.0421Λ** (J ≥ J\*(m), **J\* = 857 at the gap**) the per-state matrix `M_Ω ⪰ 0` **for every coupling vector** (2×2 ++/+−, 1×1 odd J, OO), given premise **P_L**: per entry, (R1) de-shift form, and (Loc) the de-shifted beyond-LR correction has Lemma-I density with total variation ≤ ε_out at t ≥ 1 and ≤ ε_in at t < 1, sufficient: ε_out ≤ tol(1)/4 (≥ 0.009 for J ≤ 10⁶, ≥ 0.11 for J ≤ 10⁴) and ε_in ≤ tol(0⁺)/4 (≥ 4·10⁻⁹ for J ≤ 10⁶). Equivalent ν-space form: window relative size ≤ ε_ds/2 (§4.2) + AG bulk. | **[D given P_L]** | J\* above | **P_L for the graviton**: structure [D-sketch] (L2–L4, N2, W2), magnitude [I]; scalar analogue satisfies P_L by ×10³–10⁴ |

### 0.2 Bottom line

* **No unconditional theorem.** Route L does not deliver "Ω^EE ≥ 0 for all J ≥ J\*" without an input on the spinning beyond-leading-Regge
  correction. The precise stall is P_L (§5): R1 for the EE plane-wave images, plus the absence of sub-AdS (t < 1) density in the
  correction. Both are now **J-uniform, CFT-data-free** statements about the functional's own Laplace transform (L2, L3), with ⟨TTO⟩
  ratios entering only through finitely many coefficients. They are no longer per-J window data.
* **What is rigorous and new.**
  1. **Lemma LR∞** makes the leading-Regge action positive, with an explicit remainder, for every β ≥ 1 and every mass. Combined with
     Lemma C′ and ROUTEA, the leading-Regge hybrid floor is rigorous on **all** b̂ for m ≥ 1.5Λ. The old "b̂ > 2Λ: hybrid [NE]" entry is
     closed. This also supplies Route W's premise P0 beyond J ≈ 857.
  2. **The spin-Gamma / plane-wave representation** (L2): large spin is an exact positive smoothing of a J-independent object.
     Consequences: R1 for the heavy-factor part is automatic (L3); the lightcone limit is a J-independent amplitude (L4).
  3. **Lemma W**: at large β, "ring-located" can be relaxed to "AdS-local" (within T of the ring) at bounded cost.
* **What is negative (reported prominently).**
  * **N-L1:** Bernstein positivity of the heavy factor is real but useless for sign-free bounds (smearing 10⁻²–0.5 vs correction ≤ 10⁻³).
  * **N-L2:** large J is *not* easier in the window. The de-shift window tolerance is minimal at J ≈ 350–500 (6.7·10⁻⁴) and then plateaus
    at ≈ 1.5·10⁻³ (ν < 10). It does not grow, because the decisive ν-region ~1/β and the window both sit at small ν.
  * **N-L3:** in raw (non-de-shift) form the window tolerance collapses like ≈ 0.7/J. A J-independent raw-form relative error of the
    scalar's own size (≈ 8·10⁻⁶) would exceed it at J ≈ 10⁵. **R1 is the load-bearing large-spin premise, and increasingly so.**
  * **N-L4:** the corrections do not decay like O(τ/J) or O(1/J). They tend to a J-independent lightcone amplitude of size
    O(a³/τ²)·poly(ν).
* **Explicit J\*:** 857 at the gap (b̂ ≥ 1.0421Λ), set by Lemma LR∞'s threshold β = 1. Below it, between Route W's sample range
  (J ≤ 344) and J\*, the leading-Regge floor is [NE] near the gap and the window is at its tightest. **Near-gap 150 ≤ J < 857 remains
  numeric/model.**

---

## 1. Setting

Per heavy state the Ω^EE action is a Hermitian form in the ⟨TTO⟩ couplings (UT3 §2.5). At leading Regge, EE is helicity-flip free
(off-diagonal ≡ 0, EE_FUNCTIONAL §4 [P]). The diagonal is `F_LR(β, u, σ_h)`, and beyond-LR corrections enter every entry. Lemma I
(LARGESPIN §3.1 [P]): if the de-shifted image is `∫P_{−½+iν}(cosh t)σ(dt)`, its k = 2 action is `∫Ω₂(t)σ(dt)`. Near the gap at large J,
`β ≈ log(2J/τ)` grows without bound. At J ≫ τ the state is in the lightcone regime (`e^β ≈ (Δ+J)/(Δ−J)`). The task asks whether this
can be turned into a proof for J ≥ J\*.

---

## 2. Exact large-spin structure (`lsl_frame.py`, `lsl_smear.py`)

### 2.1 The gauge identity [P + DV]
The crystal (`l1_crystal.build`) uses the symmetric gauge `x₃ = −x₁ = y, x₄ = −x₂ = x/(−x²), x₅ = e`. The heavy state enters only through
`B^{Δ/2}P_J(t)/β_{Δ,J}` (HLJ_ANALYTIC §1.2 [P]), with `B = m₃₅m₄₅m₁₂/(m₁₅m₂₅m₃₄)` a five-point conformal invariant and
`t = Z·W/√(Z²W²)` the cosine at x₅ between two covariant vectors.

The conformal integral can be re-gauged. In CHMRS's first frame (`x₅ = 0, x₁ = −e, x₂ = ∞, x₃ = (e+y)/(e+y)², x₄ = x/(−x²)`, x, y ∈ V₊):

    B = 1/(−w²),   |t| = w⁰/|w| = cosh χ_w,   w := e + x + y.                                                        (2.1)

Checks:
* F1: 2000 random x, y; B to 2·10⁻⁷ (the L = 10⁷ stand-in for ∞) and |t| to 10⁻¹⁵.
* F2: invariance of B and |t| under 500 random conformal maps, to 10⁻¹⁴. The sign of t flips only when a point crosses infinity,
  which is the (±1)^J patch convention.

The kernel (everything else) is state-independent.

### 2.2 Proposition G (spin-Gamma representation) [P for the scalar heavy factor]
Laplace's second integral (d = 3) and Euler's integral give:

    B^{Δ/2}P_J(t) = (1/2π)∮dφ (−w²)^{−(τ−1)/2} (−n_φ·w)^{−(J+1)},     n_φ = (1, cos φ, sin φ),                       (2.2)
    (−n·w)^{−(J+1)} = ∫₀^∞ ds s^J e^{s n·w}/J!     (n·w < 0 on V₊).                                                   (2.3)

Checks: F3 to 7.5·10⁻¹⁴ up to (τ, J) = (1000, 344); F4 to 3·10⁻¹⁰ up to J = 3000.

Hence, with 𝒦_ν the state-independent kernel in this gauge (Tonelli: every heavy-side factor is positive):

    Image_{τ,J}(ν) = (1/β_{Δ,J}) ∫₀^∞ ds [s^J e^{−s}/J!] 𝒢_{τ,ν}(s),
    𝒢_{τ,ν}(s) := (1/2π)∮dφ ∫∫dx dy 𝒦_ν(x,y)(−w²)^{−(τ−1)/2} e^{s n_φ·(x+y)}.                                        (2.4)

**At fixed twist, J enters only through a Gamma(J+1) average (relative width J^{−1/2}) of one J-independent, CFT-data-free function of s**,
and through the explicit normalisation β_{Δ,J}.

### 2.3 Corollary B (Bernstein / plane waves) [P]
For λ > ½, `(−w²)^{−λ} = [(2λ−1)/(2πΓ(2λ))]∫_{V₊}d³p (−p²)^{λ−3/2}e^{p·w}`. This is the Riesz transform; F5 checks it to 5·10⁻⁷. So every
scalar heavy factor is `∫dμ_{Δ,J}(p) e^{p·w}`, with **μ ≥ 0 explicit**:

* `p = s n + q`;
* s ~ Gamma(J+1);
* q ~ `(−q²)^{λ−3/2}e^{−q⁰}`, λ = (τ−1)/2;
* the image is a positive average over p of **plane-wave images** `𝒜_ν(p) = ∫∫𝒦_ν e^{p·(x+y)}`.

Concentration (Monte Carlo, CHMRS exponents, τ = 1000):

| J | β_cl | ⟨β(p)⟩ − β_cl | sd β(p) | ⟨\|p\|⟩/m |
|---|---|---|---|---|
| 50 | 0.096 | +0.023 | 0.054 | 0.990 |
| 100 | 0.183 | +0.012 | 0.059 | 0.991 |
| 344 | 0.525 | +0.008 | 0.057 | 0.992 |
| 3000 | 1.948 | +0.010 | 0.048 | 0.994 |

The width is J-independent, ≈ 1.5/√τ. Tails are Gamma/Pareto, so the mass at β(p) < β_cl − δ is ≤ e^{−O(τδ²)}.

### 2.4 Corollary (R1 reduction) [D given G]
The de-shift operator 𝒜 acts on ν only, and the heavy-side measure is ν-independent. So:
* if every plane-wave raw image is of de-shift form, `𝒜_ν^{raw}(p) = 𝒜[f_p](ν)` with f_p strip-analytic, then so is every state's raw
  image, with `f = ∫dμ f_p`;
* at leading order the plane-wave raw image is the leading-Regge raw image `𝒜[P_{β(p)}]`, which is in de-shift form by construction of X;
* **hence R1 holds automatically for the whole heavy-factor (smearing) part of the correction, at every order**, and can fail only through
  the EE kernel's own plane-wave corrections.

This is consistent with N2-e: the measured D_h excess is an a_h normalisation, i.e. exactly of de-shift form.

### 2.5 Spinning extension [D-sketch]
In this gauge the left structure `⟨T(−e)T(∞)O(0,n)⟩` is a J-independent polynomial in n: the factor (n·e)^{J−q} is a constant. The right
structure `⟨O^S(0,n)T(x₃)T(x₄)⟩` depends on x, y only through w and the inverted polarisations. Structure b gives
`(−n·w)^{−(J+1)+r_b}(−w²)^{−λ_b}Q_b(n, w, ζ₃, ζ₄)`, with r_b ≤ 4 and Q_b of bounded degree.

(2.3) then applies with `J → J − r_b`. So **each matrix entry is a Gamma(J+1−r_b) average of a J-independent function**, and the coupling
ratios enter only as the finitely many coefficients λ_aλ̄_b. This matches HLJ_ANALYTIC §1 ("Jacobi polynomials of bounded index in the same
angle t, times bounded-degree polarisation factors" [D]). It has not been checked on the spin-J engine.

### 2.6 Negative N-L1: the positive smearing is not small (`cache/smear.txt`) [DV given CHMRS (3.20) [Q]]
Relative to the local LR amplitude, the "smearing-only" image is `S(ν) = E_μ[w(p)P(cosh β(p))]/E_μ[w(p)]`, with w = \|p\|^{−2(k+2Δφ−1)}
(CHMRS `Fourier identity`), for the scalar Δφ = 3, k = 2 at τ = 1000. It is compared with the full CHMRS (3.20) correction:

| J | | ν = 0 | 2 | 5 | 10 | 20 | 30 |
|---|---|---|---|---|---|---|---|
| 100 | smearing only | −3.6e-4 | −8.5e-3 | −4.6e-2 | −9.6e-2 | +0.51 | +0.02 |
| 100 | CHMRS full | +8.9e-6 | +1.8e-5 | +4.5e-5 | +5.6e-5 | −6.6e-4 | +6.9e-4 |
| 3000 | smearing only | −1.8e-3 | +7.8e-3 | +1.0e-2 | −0.10 | −0.39 | −0.38 |
| 3000 | CHMRS full | +6.6e-6 | −1.1e-5 | −2.9e-5 | +1.1e-4 | +4.3e-5 | −1.6e-3 |

The heavy-factor smearing is an O(ν²/τ) effect. The plane-wave images' kernel corrections (the kernel in this gauge is not the exact Regge
kernel) cancel it to O(η ν²/m²). So the Bernstein decomposition is **exact and positive, but not a decomposition into "positive + small"**.

It does not bound the non-ring part by the ring part. The symmetric gauge has the exact Regge kernel (CHMRS `moment maintext`) but a heavy
factor that is not completely monotone. Neither gauge has both properties, so the subleading Regge correction is exactly the mismatch.

---

## 3. Lemma LR∞ (`lsl_lrinf.py`, `lsl_lrinf_grid.py`; `cache/lrinf*.txt`) [P + DV]

**Statement.**
* Notation:
  * ψ = ψ_k (k = 2, 3) is the exact NF02 polynomial: degree 14, ψ(0) = 0, (1−p)⁵ | ψ;
  * a₁ = ψ′(0), a₂ = ψ″(0) (ψ₂: 188.896, −1023.546; ψ₃: 50.64, 7.3·10⁻⁵);
  * `D = cosh β − 1`, `K₀ = D^{−½}`, `X = Λβ`;
  * `M₁(β) = ∫₀^∞2νP_{−½+iν}(cosh β)dν/(e^{2πν}+1)` (> 0 for all β, [certified] h13_cert_m1).
* Then:

      Ω_k(β) = a₁M₁(β)/Λ² − a₂/(2b̂³) + E_k(β),
      |E_k| ≤ (√2/π)[T_a + T_b + T_c],

  with:
  * `T_a = (K₀/Λ)d₁(X)`;
  * `T_b = K₀/(4DΛ³)·∫_X^∞x²R_∞`;
  * `T_c = ∫₀^{β/2}R(Λt)ρ_b(t)dt + R_∞(Λβ/2)(π/√2)P_{−½}(cosh β)`;
  * `R_∞(x) = Σ_{j≥3 odd}|ψ^{(j)}(0)|/x^{j+1} + Σ_{j≥5}|ψ^{(j)}(1)|/x^{j+1}`;
  * `R = min(R_∞, a₁/x² + ∫|ψ|)`;
  * `d₁(X) = Σ_{j≥2 even}|g^{(j)}(0)|/X^{j+1} + Σ_{j≥5}|g^{(j)}(1)|/X^{j+1}` (g = ψ/p);
  * `ρ_b(t) = K₀[t⁴cosh t/(48D) + (3/8)2^{5/2}((cosh t − 1)/D)²]`.

**Proof.**
1. Mehler–Dirichlet: `Ω = (√2/π)∫₀^βΦ(Λt)(cosh β − cosh t)^{−½}dt`, with `Φ(x) = ∫₀¹ψ cos(px)dp` (Fubini).
2. Exact closed form, by integrating a polynomial by parts: `Φ(x) = −a₁/x² + r(x)`, with
   `r(x) = Σ_{j≥3 odd}(−1)^{(j+1)/2}ψ^{(j)}(0)/x^{j+1} + Re Σ_{j≥5}(−1)^jψ^{(j)}(1)e^{ix}/(ix)^{j+1}`. Hence |r| ≤ R_∞ for every x > 0.
3. The Hadamard finite part of the −a₁/x² term is `(a₁/Λ²)(π/√2)M₁(β)`, using the identity
   `(√2/π)[K₀/β − ∫₀^β(K − K₀)t^{−2}dt] = M₁(β)`. It follows from Mehler–Dirichlet, the distributional `∫₀^∞ν cos(νt)dν = −t^{−2}` (f.p.)
   and `∫₀^∞ν tanh(πν)P_{−½+iν}(cosh β)dν = 0` for β > 0 (Mehler–Fock) [D]; ROUTEA §4.1 has the equivalent c(t) form. It is checked to
   12 digits at β = 0.1, 1, 4, 10. Since ΦK is integrable, the two finite parts (of −a₁/x² and of r) recombine exactly.
4. Split `K = K₀ + K₀t²/(4D) + ρ`, with 0 ≤ ρ ≤ ρ_b on t ≤ β/2 (Taylor of (1−ε)^{−½} at ε ≤ ½, and of cosh).
5. Assemble:
   * the `K₀` term gives `(K₀/Λ)[∫₀^XΦ − a₁/X]`, which is T_a, using `∫₀^XΦ = ∫₀¹(ψ/p)sin(pX)dp` and its closed form;
   * the `K₀t²/(4D)` term gives `(K₀/(4DΛ³))C(X)`, with `C(X) = ∫₀^X(x²Φ + a₁)` → `C_∞ = −(π/2)a₂`. This is checked: 1596.5 / 1602.1 /
     1605.0 at X = 200 / 400 / 800 → 1607.78. It produces the flat tail −a₂/(2b̂³) exactly, and its remainder is T_b;
   * the ρ term is T_c. On [β/2, β], ρ ≤ K and `∫K = (π/√2)P_{−½}(cosh β)`.

∎ [P] (the identity of step 3: [D + DV]). The bound terms are explicit one-dimensional integrals of elementary positive functions, evaluated in floating point [DV]. For β ≥ 12,
every ratio below decreases, using h13's region-(c) lower bound on M₁ [D].

**Hybrid-level remainder.** Define `r(β) := (2|E₂| + |E₃|)/(A₁,min M₁/Λ²)`, with `A₁,min = 2a₁⁽²⁾ − a₁⁽³⁾ = 327.15` (u ≤ 1, worse σ). The
fine grid has step 0.01 on [0.9, 3], 0.05 to 16 and 0.5 to 64:

| β ≥ | 1.0 | 1.2 | 1.5 | 2.0 | 3.0 | 5.0 | 64 |
|---|---|---|---|---|---|---|---|
| sup r | 0.577 | 0.225 | 0.070 | 0.015 | 1.6e-3 | 1.0e-4 | 1.6e-7 |

Hence **`F_LR(β, u, σ) ≥ (1 − r(β))A₁M₁(β)/Λ² + 2047.1/(2b̂³) > 0` for all β ≥ 1, u ≤ 1, σ = ±1.**

**Direct check** [DV]: the Mehler-φ quadrature with closed-form Φ agrees with the two-term form to ≤ 4·10⁻⁶ (k = 2) on
β ∈ [0.3, 4], and (k = 3) to ≤ 2·10⁻⁵ on [1, 4] (8·10⁻⁴ at β = 0.3, where the k = 3 form's omitted terms matter), and with `ls_core.OmLR` to ≤ 6·10⁻⁷ at β ≤ 2. At β ≳ 10 the direct quadrature itself fails (cancellation to e^{−β/2}/Λ²);
there the lemma is the only reliable evaluation.

The bound is loose at β < 1, because the oscillatory p = 1 terms are counted in absolute value. One integration by parts in t would
lower the threshold. That was not done.

### 3.3 Hybrid vs F_LR, and J\*(m) (`cache/hybcmp.txt`)
The project's hybrid (flat NF02 partial wave + curvature) differs from F_LR by `(hybrid − F_LR)/F_LR`:
* `+1·10⁻⁶ … +3·10⁻⁷` at J = 1000–3000 (gap and 1.3m_gap, both helicities);
* `+1.7·10⁻³` at (J = 100, 1.3m_gap).

At large J the floor is F_LR.

| m | m_gap | 1.1m_gap | 1.3m_gap | 1.5Λ | 2Λ | 5Λ |
|---|---|---|---|---|---|---|
| J\* (β ≥ 1) | 857 | 987 | 1275 | 780 | 1041 | 2602 |
| J at β ≥ 1.5 (r ≤ 7 %) | 1737 | 2034 | 2702 | 1232 | 1642 | 4105 |
| J at β ≥ 2 (r ≤ 1.5 %) | 3188 | 3776 | 5113 | 1760 | 2347 | 5867 |

**Caveat (weight realisation).** The lemma concerns the target weights ψ_k. The published realisations satisfy the de-shift identity only
up to `ε_W P_{−½}(cosh β)` (HK_RESIDUAL (ε-DS)). For the continuum weight (ε_W ≈ 1.1·10⁻⁹ per unit ∫|ψ|) this is ≲ 3·10⁻⁶β² of F_LR. That is
< 1 % for β ≲ 60 (J ≲ 10²⁸). The grid weight `d_abel_X_1000.txt` satisfies the Abel identity only to t ≈ 22–25 (J ≲ 10¹²; UT3 §2.2).
"All J" therefore refers to the exact min-norm functional.

---

## 4. The window at large spin

### 4.1 Lightcone uniformity (L4) [Q + DV scalar; D-sketch spinning]
Large J at fixed τ is not a 1/J expansion. CHMRS state that the coefficient of `P^{(n)}/m^{2j}` grows at most like `η^{j+n}`, and
`η/m² = 1/m² + (2J+1)²/(2m⁴) ≤ 1.5/(τ−2)²` for all J. In the J → ∞ limit `(η²−1)P′ → η(±iν − ½)P` per travelling wave, so CHMRS (3.20)
becomes a J-independent relative factor:

    F±(ν, τ) − 1 = [(±iν − ½)(ν² + c₁) + (ν² − a(a+1))c₂ + c₃] / (12(τ−2)²),
    a = 2Δφ + k − 3/2,   c₁ = −3a² + 6a(2Δφ−3) − 2,   c₂ = 3a − 6Δφ + 19/2,   c₃ = −2(Δφ−4)₃ − 2(a+½)(1 − a² − a/2).

Scalar Δφ = 3, k = 2: a = 6.5, c₁ = −11.75, c₂ = 11, c₃ = 623.

CHMRS (3.20) correction / local LR amplitude, τ = 1000 (`cache/smear.txt`):

| J | β | ν = 0 | 0.5 | 1 | 2 | 5 | 10 |
|---|---|---|---|---|---|---|---|
| 26 | 0.052 | 9.8e-6 | 1.2e-5 | 1.4e-5 | 2.0e-5 | 5.5e-5 | 1.8e-4 |
| 100 | 0.183 | 8.9e-6 | 1.1e-5 | 1.3e-5 | 1.8e-5 | 4.5e-5 | 5.6e-5 |
| 1000 | 1.100 | 6.6e-6 | 8.1e-6 | 9.4e-6 | 3.9e-6 | 4.3e-6 | −2.9e-5 |
| 10⁴ | 3.05 | 6.9e-6 | 7.8e-6 | −4.1e-6 | 4.7e-6 | −1.4e-5 | 5.6e-5 |
| 10⁶ | 7.60 | 7.5e-6 | −7.4e-6 | 8.4e-6 | −1.2e-6 | 2.9e-5 | 1.1e-4 |
| ∞ (\|F₊ − 1\|) | — | 7.8e-6 | 8.0e-6 | 8.7e-6 | 1.1e-5 | 3.0e-5 | 1.2e-4 |

The sign pattern at ν ≥ 1 reflects the travelling-wave phase (`F₊ ≠ F₋`: the iν-odd term is a ring displacement, class C).

**Answer to the task's question.** The Regge-window correction **cannot** be bounded by an O(τ/J) or O(1/J) relative term. It tends to a
J-independent lightcone amplitude, uniformly O(1/τ²)·poly(ν, a). Only the J-dependent remainder is O(τ/J)-type (e^{−2β} Harish-Chandra
terms, Gamma smoothing).

For the graviton:
* the structural statement (ring classes at every order, J-uniform coefficients) is the N2/W2 [D-sketch] together with L2;
* Route W's W1a/W1c remove the O(1/m) pieces;
* the size of the spinning F − 1 is [I].

### 4.2 Window tolerance map, J = 26 … 10⁶ (`lsl_window.py`, `cache/window_map.txt`) [DV]
Floor `F_LR` (worse helicity; Lemma LR∞ at β ≥ 2). The tolerated relative size ε, against the local LR amplitude, of a sign-free correction
on ν < ν_W, with κ = ½ of the floor:
* **de-shift form**: `ε_ds = κF_LR/∫_{ν<ν_W}(2ψ₂ + u|ψ₃|)E dν/Λ`;
* **raw form**: `ε_raw = κF_LR/∫(2|χ|E_raw + u|ψ₃|E)dν/Λ`, where χ = X/(R₋+R₊) and E_raw is the amplitude of `R₋P(ν+i) + R₊P(ν−i)`;
* `Γ_tol` (c₂) is the HLJ-style absolute sup-norm tolerance.

Gap states:

| J | 26 | 100 | 150 | 200 | **344** | **500** | 1000 | 2000 | 5000 | 10⁴ | 10⁵ | 10⁶ |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| β | 0.052 | 0.18 | 0.26 | 0.34 | 0.52 | 0.69 | 1.10 | 1.61 | 2.40 | 3.05 | 5.31 | 7.60 |
| ε_ds(ν<10) | 0.18 | 4.3e-3 | 1.7e-3 | 1.0e-3 | **6.7e-4** | **7.0e-4** | 9.0e-4 | 1.2e-3 | 1.4e-3 | 1.6e-3 | 1.5e-3 | 1.3e-3 |
| ε_ds(ν<20) | 4.5e-2 | 1.1e-3 | 5.5e-4 | 4.0e-4 | **3.0e-4** | 3.1e-4 | 3.8e-4 | 4.7e-4 | 5.6e-4 | 6.0e-4 | 5.8e-4 | 4.9e-4 |
| ε_raw(ν<10) | 1.6e-3 | 3.2e-4 | 2.7e-4 | 2.8e-4 | 2.7e-4 | 2.6e-4 | 2.4e-4 | 2.0e-4 | 1.2e-4 | **7.1e-5** | **7.8e-6** | **6.8e-7** |
| Γ_tol(ν<10), c₂ | 1.8e5 | 5190 | 2200 | 1380 | 879 | 849 | 1080 | 1570 | 2650 | 3840 | 1.1e4 | 3.0e4 |
| ε_ds(ν<ν_B) | 5.6e-3 | 1.6e-3 | 1.3e-3 | 1.3e-3 | 1.7e-3 | 2.4e-3 | 5.7e-3 | 1.3e-2 | 3.4e-2 | 6.2e-2 | 0.41 | 2.9 |

* **1.3m_gap** has the same pattern: ε_ds(10) min 7.1·10⁻⁴ (J = 500), plateau 1.5–1.6·10⁻³; ε_raw(10) = 1.3·10⁻⁵ (10⁵), 1.2·10⁻⁶ (10⁶).
* **Calibration against HLJ's T_B(c = 3)** (different floor, 0.9λ₀N_s, plus contour terms): Γ_tol(ν < ν_B) = 1970 / 2740 c₂ at
  J = 100 / 344, against T_B = 1402 / 1095. Same order of magnitude.
* **Reading:**
  * the HLJ window ν < ν_B = 3m/(2J+1) shrinks to 0 at large J and becomes irrelevant. But the region not covered by AG checks or
    direct data, ν ≲ 10–20 (fold region; n1 fails below ν ≈ 20), does not shrink;
  * the action at large β is decided at ν ≲ 1/β, inside it;
  * hence the de-shift tolerance plateaus rather than grows. The **hardest band is J ≈ 350–500** near the gap;
  * in raw form the tolerance falls like ≈ 0.7/J (≈ e^{−β}). The raw image grows like e^{β/2} while the de-shifted floor decays like
    e^{−β/2}/β; X removes the growth only for de-shift-form images.

### 4.3 Large-β tolerance profile and Lemma W (`lsl_ringtol.py`, `cache/ringtol.txt`)
`tol(t) = F_LR(β)/(2|Ω₂(t)| + u|Ω₃(t)|)` is the weight tolerated at radius t, per unit ring weight, at the gap, worse helicity:

| J (β) | t = 0.003 | 0.03 | 0.1 | 0.3 | 1 | β/2 | β |
|---|---|---|---|---|---|---|---|
| 857 (1.00) | 8.0e-7 | 3.7e-4 | 1.4e-2 | 0.27 | 0.91 | 0.60 | 0.91 |
| 3000 (1.95) | 6.7e-7 | 2.8e-4 | 1.0e-2 | 0.20 | 0.70 | 0.69 | 0.96 |
| 10⁴ (3.05) | 4.5e-7 | 1.8e-4 | 6.6e-3 | 0.13 | 0.45 | 0.54 | 0.99 |
| 10⁵ (5.31) | 1.4e-7 | 5.4e-5 | 2.0e-3 | 0.039 | 0.14 | 0.25 | 1.00 |
| 10⁶ (7.60) | 3.7e-8 | 1.4e-5 | 5.3e-4 | 0.010 | 0.036 | 0.12 | 1.00 |

Ring-located corrections (support t ≥ (1−ε_r)β) are tolerated at w_tol:
* 0.67–0.85 for ε_r = 0.03;
* 0.49–0.83 for ε_r = 0.1;
* 0.20–0.75 for ε_r = 0.3;

over β = 1.2–7.6. M₁ is decreasing on [0.05, 12] (sampled).

At large β the profile is flat on t ≳ 0.3 (b̂ ≳ 0.3Λ). **Only sub-AdS impact parameters t ≲ 0.1 (b̂ ≲ 100, the flat/bulk-point zone) carry
the 10³–10⁷ penalty.** Contrast H_LJ at J ≤ 344, where tol(β/3) ≈ 6·10⁻³ (LARGESPIN §3.2).

**Lemma W (H₂ Paley–Wiener locality) [P given PW [Q] and Lemma I].** Let k be a finite complex radial measure on H₂ supported in the ball
B(T), with spherical transform G(ν) = k̂(ν); equivalently, G is even, entire, and of exponential type ≤ T. Then the de-shifted correction
`f(ν) = G(ν)P_{−½+iν}(cosh β)` is the spherical transform of `k ∗ δ_β` (bi-K-invariant convolution is multiplicative), which is supported
in t ∈ [β − T, β + T] with total variation ≤ ‖k‖. By Lemma I, its action is ≤ `‖k‖·max_{|t−β|≤T}(2|Ω₂(t)| + u|Ω₃(t)|)`.

For T ≤ β − 1 and β ≥ 1 this is ≤ ‖k‖·F_LR(β)/tol(β − T), with the tol(t ≥ 1) values above (Lemma LR∞ supplies the monotonicity). Polynomial
G (CHMRS-type ν^{2n}) are the case k = p(−Δ_H)δ_o. Their exact cost is `Λ^{2n}Ω_{p^{2n}ψ}(β)` (LARGESPIN identity B). At large β it is
suppressed relative to F_LR, because p^{2n}ψ has zero first moment.

---

## 5. Theorem L (conditional)

**Theorem L.** Assume A1′, A2 and Δ_g = 1000. Let O be a heavy primary in T × T (τ ≥ Δ_g) with **b̂ ≥ 1.0421Λ**, i.e. **J ≥ J\*(m) =
⌈(1.0421 m − 1)/2⌉** (m in AdS units; J\* = 857 at the twist edge). Let `M_Ω(O)` be its per-state matrix (UT3 §2.5): 2×2 (++, +−) for even J ≥ 4, 1×1
(+−) for odd J, OO. Write `M_Ω = diag(F_LR(β,u,σ_h)) + δM`. Assume premise **P_L** for every entry δM_{ab}:
* **(R1)** the raw correction is of de-shift form, so the entry's k = 2 cost is the ψ₂ cost of a de-shifted image f_{ab} (the k = 3 piece
  has no de-shift);
* **(Loc)** the H₂ density σ_{ab} of f_{ab} (Lemma I; any admissible choice), in units of the leading-Regge ring weight, satisfies
  `|σ_{ab}|([0,1)) ≤ ε_in` and `|σ_{ab}|([1,∞)) ≤ ε_out`, with

      ε_in·sup_{t<1}(2|Ω₂| + u|Ω₃|)(t) + ε_out·sup_{t≥1}(2|Ω₂| + u|Ω₃|)(t) ≤ ½·min_h F_LR(β, u, σ_h).                      (5.1)

Then **`M_Ω(O) ⪰ 0` for every coupling vector λ.** Sufficient numbers (each half of the budget), using the tolerance profile at the gap (Ω_k decreasing in t: [NE] on [0.002, 4] (LARGESPIN), the two-term form beyond):
* `ε_out ≤ tol(1)/4`: 0.23 (J = 857), 0.11 (J = 10⁴), 0.009 (J = 10⁶);
* `ε_in ≤ tol(0⁺)/4`: ≈ 10⁻⁷ (J = 857) … 4·10⁻⁹ (J = 10⁶).

*Proof.*
1. Lemma LR∞ gives F_LR(β, u, σ) > 0 for both helicities at β ≥ 1, u ≤ 1.
2. Lemma I and (5.1) give `|δM_{ab}| ≤ ½ min_h F_LR` for each entry.
3. A Hermitian 2×2 with diagonal ≥ F/2 and off-diagonal ≤ F/2 is ⪰ 0; the 1×1 cases follow from step 2 directly. ∎

**[D given P_L]**; the [P + DV] floor enters through Lemma LR∞.

**ν-space form (the one the data test).** Under (R1), (Loc) follows from two pieces. The first is a window bound `|f_{ab}(ν)| ≤ ε·E(ν)` on
ν < ν_W, with ε ≤ ½ε_ds(ν_W) (§4.2): ≥ 4·10⁻⁴ for ν_W = 10 and J ≥ J\*(m) (min ε_ds(ν<10) ≈ 8.5·10⁻⁴ there, sampled). The second is the AG bulk on ν ≥ ν_W (travelling waves of frequency
≥ 0.97β), whose cost the ring tolerance (≥ 0.67 at ε_r = 0.03) absorbs. The one-line algebra is the HLJ contour lemma with κ = ½.

**Status of P_L.**

| ingredient | scalar analogue | graviton |
|---|---|---|
| R1 | exact (CLS-D magic-block identity [Q], N2) | heavy-factor part [D given L2] (§2.4); a_h part [D] (N2); **EE plane-wave kernel part [I]**. Route W W3: R1 is load-bearing at J ≳ 200 [DV], consistent with N-L3 |
| Loc, finite orders | every CHMRS term is ring-supported, i.e. ε_in = 0 at each order; costs ≤ 6·10⁻⁵ of the floor at J ≥ 150 (N2) | ring classes at every order [D-sketch] (N2, W2, L2); sizes [I] |
| Loc, beyond all orders (density at t < 1) | far region `O(Δ^N e^{−cΔ})`, c ≥ 0.336 (RD_CLOSE [D]). This is ≤ e^{−336} against ε_in ≈ 10⁻⁸ | RD for EE [D given P_cyc] (RD_CLOSE); not checked at J ≥ 150 |
| window size | ≤ 1.2·10⁻⁴ on ν ≤ 10 for all J (L4) vs ε_ds/2 ≥ 4·10⁻⁴ (J ≥ J\*) | measured non-ring remainder 2–5·10⁻⁴ at J = 29–100 (§7.2); J-independence [D-sketch] |

**What would make TL a theorem.** P_L reduces, via L2 and L3, to two CFT-data-free statements about the EE functional's plane-wave images
`𝒜_ν(p)`:
* **(PW-R1)** each is of de-shift form;
* **(PW-Loc)** its H₂ density at t < β(p) − β(p)/2 is ≤ e^{−c|p|}.

Both are statements about the Laplace transform of an explicit kernel (with ≤ 4 polarisation insertions). They are open [I]. The obstacle
is the same as Route W's item 1: the EE-kernel polarisation contraction beyond leading order.

---

## 6. Route-L ideas assessed

| idea (task) | verdict | where |
|---|---|---|
| β(J, m) grows ⇒ localisation | yes: β ≈ log(2J/τ); at J ≫ τ the action is an endpoint (ν ≲ 1/β) quantity, A₁M₁(β)/Λ² — Lemma LR∞ [P + DV] | §3 |
| Lemma I + tolerance profiles; Lemma C′ to b̂ ≤ 119 | extended: rigorous floor for β ≥ 1 at all m; profile flattens at large β (sub-AdS only penalised) | §3, §4.3 |
| lightcone factorisation ⇒ O(τ/J) or O(1/J) relative window corrections | **no**: the limit is a J-independent amplitude F(ν,τ) − 1 = O(a³/τ²)·poly(ν); only the remainder is O(τ/J) | §4.1 |
| Phragmén–Lindelöf / maximum modulus between lightcone and Regge regimes | not usable as such: at fixed ν the LR image has zeros in J (∝ cos νβ(J)), so relative corrections have poles in J. Growth bounds in Im J are not available, and the Gamma representation (2.4) gives analyticity in J but no bound. The usable analytic statement is in ν: Paley–Wiener on H₂ (Lemma W) | §4.3 |
| positivity of the Lorentzian-inversion / light-ray integrand | the per-state positive structure is the heavy factor's Bernstein measure (exact, §2.3). It is cancelled by kernel corrections at O(ν²/τ), so it gives no sign-free bound (N-L1). ANEC/EEC positivity holds for the sum over states, not per state | §2.6 |

---

## 7. Checks against existing data

### 7.1 Scalar (exact skeleton / CHMRS)
* CHMRS (3.20) is reproduced on the exact skeleton at J = 24–100 (N2 §3). L4's table shows its relative size is J-uniform from J = 26 to 10⁶
  and converges to the lightcone F (7.5·10⁻⁶ vs 7.75·10⁻⁶ at ν = 0).
* In c₂: 8.9·10⁻⁶ × m²(J = 100) = 10.6 c₂, consistent with the "exact scalar 6.5–10 c₂, J-independent" of HLJ/N2.

### 7.2 Graviton window data in relative units (quoted; m² from the state)

| state | quantity (source) | c₂ | relative | ε_ds(ν<10) at that J | L = ε_ds/relative if J-independent: hardest band (6.7e-4) / J ≥ J\* (≥ 8.5e-4) |
|---|---|---|---|---|---|
| (1029, 29) | \|Γ\| on ν = 25–47 (SMALLNU) | ≤ 310 | ≤ 2.9e-4 | ≈ 0.1 | 2.3 / 2.9 |
| (1040, 32) | \|Γ\| on ν = 25–40 (SMALLNU) | 450–550 | 4.2–5.1e-4 | ≈ 0.08 | 1.3–1.6 / 1.7–2.0 (values at ν = 25–40, outside ν < 10) |
| (1064, 64) | C2 \|Γ − κ\| on ν = 13–33 (HLJ) | ≤ 350 | ≤ 3.1e-4 | ≈ 1e-2 | 2.2 / 2.7 |
| (1100, 100) | D_h-class excess (de-shift form, harmless by congruence/N2) | 2154 | 1.8e-3 | 4.3e-3 | not counted (ring class, de-shift form) |
| (1100, 100) | non-ring remainder after A+C (N2 F1) / A+B+C (DIRCV §8: \|r\| ≲ 100 ± 120) | 340–380 / ≲ 280 | 2.8–3.2e-4 / ≲ 2.3e-4 | 4.3e-3 | 2.1–2.4 / ≥ 2.9 (central) |

**Reading.**
* If the graviton's non-ring relative window correction is J-independent at fixed τ (as L4 shows for the scalar), its measured value
  2–5·10⁻⁴ sits at:
  * L ≈ 1.3–3 against the sign-free de-shift tolerance in the hardest band (J ≈ 350–500);
  * L ≈ 3–7 on the large-J plateau (ν < 10; ε_ds ≈ 1.5·10⁻³);
  * L ≈ 1–2 with ν_W = 20.
* It is 40× the scalar's.
* In raw form it would exceed ε_raw(10) already at J ≈ 150 and fail by ×10³ at J = 10⁶. This is the quantitative content of N-L3.

### 7.3 Internal
* The window map's ε_ds at J = 100 equals Γ_tol/m² (4.3·10⁻³ both), as it must at small β where E ≈ 1.
* The floor matches the hybrid to ≤ 2·10⁻⁶ at J ≥ 1000 (§3.3) and LARGESPIN's `barrier_1000` values to ≤ 6·10⁻⁷.

### 7.4 Cross-check with Route W (`LARGESPIN_W.md`, read after §§2–6 were fixed)
* **Agreement:**
  * no unconditional theorem;
  * the residual is the spinning beyond-LR structure at O(1/m²);
  * **R1 is load-bearing from J ≈ 150–200 on.** W: stripped ν² classes are tolerated only at 1.4–7.5 c₂ for J = 250–344; here ε_raw ≈ 0.7/J
    beyond 10³;
  * every finite order is ring class (W2; §4.1).
* **Complementary:**
  * W's theorem rests on P0 (the leading-Regge hybrid floor) and is tabulated on the H_LJ samples (J ≤ 344). Lemma LR∞ proves P0's
    leading-Regge part for β ≥ 1, J ≥ 857 at the gap;
  * §4 extends the tolerances to J = 10⁶ and shows they do not degrade in de-shift form (plateau) but collapse in raw form;
  * W's W1a/W1c (Weyl parity, congruence) remove the O(1/m) normalisation pieces. L2–L3 remove the heavy-factor part of R1. The two routes
    cover different parts of the same residual.
* **Different J\*:** W gets J\* = 50 given scalar-sized O(1/m²) coefficients (P2), R1 and an O(1/m⁴) bound. L's J\* = 857 is set by
  where the leading-Regge floor is *proved*, not by the size premise. At J < 857 L relies on the [NE] hybrid floor, as W does.

---

## 8. Revised large-spin inventory (Result D, Δ_g = 1000)

| region | before (LARGESPIN, N2, DIRCV, SMALLNU, HIGHSPIN) | now | basis |
|---|---|---|---|
| leading-Regge floor, m ≥ 1.5Λ, all b̂ | [P + DV] b̂ ≤ 119; [certified] [91.2, 2Λ]; [NE] b̂ > 2Λ | **rigorous on all b̂** ([P + DV] / [certified] / [P + DV]) | Lemma C′, ROUTEA, **LR∞** |
| leading-Regge floor, m ∈ [Λ, 1.5Λ) | [NE] for b̂ > 119 (1.3–1.5Λ) and for all b̂ near the gap | **[P + DV] for b̂ ≥ 1.0421Λ**; [NE] below | LR∞ |
| near gap, J = 12–32 | [NE±] / [NE + M] (window now measured) | unchanged | SMALLNU, HIGHSPIN |
| near gap, even J ≥ 50, odd J ≥ 33, J ≤ 150 (H_LJ, data to J = 100) | [D-sketch + M], leverage 1.6–4.5 central at (1100, 100) | unchanged | DIRCV §8 |
| **near gap, 150 ≤ J < 857** | [M] ("nothing computed at J ≥ 150") | **[D-sketch + M]**, explicit sign-free tolerances (ε_ds(ν<10) ≥ 6.7·10⁻⁴, min at J ≈ 350–500; ε_raw ≈ 2.7·10⁻⁴); floor [NE] (hybrid scan); **the hardest band; remains numeric/model** | §4.2; W (J ≤ 344) |
| **J ≥ J\*(m), b̂ ≥ 1.0421Λ (J\* = 857 at the gap)** | [M] | **[D given P_L]** (Theorem L): floor [P + DV]; P_L is J-uniform and CFT-data-free (L2–L3); scalar analogue satisfies it with margin ×3 (window, sign-free) to ×10⁴ (Lemma I). P_L for the graviton: structure [D-sketch], magnitude and R1 for the EE plane-wave kernel [I] | §5 |
| J beyond the published weight's validity (grid X: J ≳ 10¹²; continuum ε_W: J ≳ 10²⁸) | — | the theorem refers to the exact min-norm functional; the realisation error is explicit (§3.3) | HK_RESIDUAL |

**Letter wording (suggested).** "At large spin, a heavy state's action is a positive average of leading-Regge actions at nearby impact
parameters. That leading action is positive with an explicit error bound at every impact parameter beyond the AdS radius, which covers all
states with J ≳ 0.52m. Beyond leading Regge, the corrections relevant at large spin are a single spin-independent function, the
lightcone amplitude. Positivity follows if this amplitude preserves the de-shift structure of the functional and is AdS-local. This holds
with large margin for the scalar analogue, and is supported but not derived for gravitons." Do not write "proved for J ≥ J\*".

---

## 9. Corrections and remarks for existing documents (not applied)

1. **UT3 §5.2 / LARGESPIN §7 row "b̂ > 2Λ: hybrid [NE]"** → "leading-Regge action [P + DV] for b̂ ≥ 1.0421Λ, every m ≥ Λ (LARGESPIN_L
   Lemma LR∞)". The hybrid differs from it by ≤ 2·10⁻⁶ there.
2. **EE_DG1000_LARGEJ §0.3** ("beyond-Regge weight relative to N_s grows like J³"). This holds only in the flat-tail regime b̂ ≲ 4Λ^{2/3}.
   In relative units the window tolerance has a minimum at J ≈ 350–500 and plateaus beyond (§4.2). It grows ∝ J in c₂ units.
3. **LARGESPIN §6(d)** ("the light-cone regime occurs only at b̂ > 2Λ, m > 2Λ"). Near the gap the lightcone regime β ≳ 1 starts at
   J ≈ 860. It is reached by every near-gap trajectory at large J, and δ there reduces to a J-independent amplitude (§4.1).
4. **SMALLNU_WINDOW §0 item 5** ("for J ≳ 32 the window is again a live input"). In relative units it stays live up to J ≈ 10³. Its
   tolerance then plateaus (de-shift form) or collapses ∝ 1/J (raw form). The latter is the strongest argument for settling R1.

---

## 10. Reproduce, caveats

    cd ads4/code/lsl && bash run_all.sh        # ≈ 0.5 CPU-h, one process at a time

Caveats:
1. Δ_g = 1000 only.
2. Lemma LR∞: the inequalities are proved, but the bound's constants are floating-point evaluations of explicit 1D integrals on a β grid
   (step 0.01–0.5), so [P + DV], like Lemma C′, not [certified]. For β > 64 the analytic tail is [D]. M₁ > 0 and its β ≥ 12 lower bound are
   [certified] (h13).
3. The window map samples 16 spins per mass at 2 masses. ε_ds uses the local-amplitude envelope (running max over ±π/β). The raw map uses
   χ = X/(R₋+R₊) interpolated across the real zeros of R₋+R₊ (±0.03). The k = 3 piece is counted without de-shift, with |ψ₃|.
4. §2.6 uses CHMRS's scalar exponents and Fourier-identity weight. It is a Monte Carlo with 4·10⁵ samples (statistical error ~10⁻³ on the
   smearing row); the negative conclusion is O(10⁻¹) vs O(10⁻⁴).
5. Proposition G is proved for the scalar heavy factor. The spinning version (§2.5) is a sketch and was not checked on the spin-J engine.
6. Theorem L's (Loc) uses "any admissible" Lemma-I density. For corrections known only on ν ∈ [0, Λ] the density is not unique. The premise
   is that some admissible choice satisfies the bounds.
