# Result Q, Einstein instance of Lemma A + D_LR: what can be evaluated, a Witten-diagram spinning check with a J₀ = 2 spin-2 exchange, and what remains

*22 Nov 2026, Einstein-check agent. Unrefereed, not committed. No existing document or code file was edited; no `spinfunc` / `rs2_*` file
was touched. New code: `ads4/code/resultq/ein_dbar.py` (D-function engine with Lorentzian continuation), `ein_run.py` (production).
Outputs: `resultq/cache/ein_*.txt`, `ein_run_*.log`. Read-only imports: `og_free` (free-scalar `⟨TTTT⟩`, gauge, `f_k`, contour), `og_run`
(control helpers). Compute: ≈ 4 h CPU, ≤ 2 workers, 1 thread each, < 0.6 GB per job.*

Grades: **[P]** proved; **[D]** derived; **[D-sketch]**; **[DV]** derived and numerically verified; **[NE]** numerical evidence;
**[M]** model; **[A]** assumption; **[I]** open.

---

## 0. Summary

1. **The Einstein `⟨TTTT⟩` itself was still not evaluated [I].** A literature check (§2.1) finds no representation that the
   position-space, Lorentzian, complex-`t` light-ray pipeline can use:
   * the explicit tree-level AdS₄ graviton `⟨TTTT⟩` exists only in momentum space (Raju 1201.6452; Albayrak–Kharel 1908.01835);
   * Caron-Huot–Li 2102.08160 (helicity basis) gives mean-field `[TT]` data and tree-level anomalous dimensions for *currents*, not
     the Einstein `[TT]` data in every structure;
   * the closed-form graviton exchange (Giecold 1205.6444, D'Hoker et al. hep-th/9903196) is for *scalar* externals. In AdS₄ it does
     not truncate to integer-index D-functions (it needs half-integer `Δ`, i.e. elliptic `D̄`-functions).

   So option (a) is out of reach in one session, and the Einstein instance stays unchecked in the strict sense.
2. **Instead: a direct check on correlators that carry every *qualitative* Einstein feature the free-field checks lacked (option (b))
   [DV].** Test operator `T′ = σ ⊗ T_free` (Δ = 5, J = 2) in the tensor product of the free scalar (og_free) and an AdS₄ scalar theory
   whose `O(λ)` four-point function `G_σ` is a sum of **Witten diagrams evaluated exactly**:
   * contact diagrams of Mellin degree 0, 2, 3 (`J₀ = 0, 2, 2`) and 4 (`J₀ = 4`, the R⁴-type control);
   * a **massive spin-2 exchange** (Δ_e = 4, τ = 2) in the s-, t- and u-channel, with Regge intercept **J₀ = 2** and a single-trace
     Mellin pole.

   These bring in what free fields and GFF lacked:
   * Li₂/log transcendentality, with logs (anomalous-dimension structure) at **fixed** external dimension;
   * genuine spin-2 Regge growth;
   * a single-trace spin-2 pole in the (13) channel.
3. **Lemma A holds pointwise for all of them, at the 10⁻¹³–10⁻¹¹ level [DV].** The test is swapped form, D_LR, CLS-D T-shift, generic null (2,4)
   polarisations, k = 4 (Ψ₄) and k = 5 (Ψ₅), at three generic angle pairs:
   * Lemma-A residual ≤ 1.9·10⁻¹² of the largest family (n = 60), 2.6·10⁻¹⁴ – 9·10⁻¹⁴ at n = 80;
   * **every 1↔3 pair identity separately** (σ-exchange s ↔ t, u ↔ u, each contact, times each free Wick class) ≤ 5.5·10⁻¹³;
   * EE structure: zero at every angle (Δφ = 0.05 – 3), ≤ 2·10⁻¹¹ of the largest family. The s- and t-exchange families are
     individually large and cancel;
   * collinear approach: the zero persists to Δφ = 0.025 (≤ 2·10⁻¹⁰), with families growing like `Δφ^{−3}`, as in the free check.
4. **Controls fail at O(1)**, with one important refinement (§3.4):
   * non-single-valued σ-part (not a microcausal Wightman function): 1.3 – 2.5;
   * non-primary free factor: 1.3 – 8.5;
   * light-transform weight `t^{5+q}`, q = +1, +2: 3.3 – 5.2;
   * wrong symmetrisation: 0.33 – 5.9.

   **The weight control is one-sided.** q = −1 still gives 10⁻¹⁰, and the same holds for the free T of the OG check (q = −1, −2, −3:
   10⁻⁶ – 10⁻¹⁵, §3.4). So the zero tests "weight ≤ Δ + J − 2" (decay of the commutator light-ray integral), not the exact
   light-transform weight. RESULT_Q_OG_CHECK §3.4's reading of that control should be qualified.
5. **Where J₀ enters, and a finding about J₀ ≥ k (§3.5–3.6) [DV in the test family].**
   * The swapped double integral has three endpoint regimes. The single endpoint is a lightcone limit and J₀-blind. The same-side joint
     endpoint `t₁, t₃ → 0` is the Regge limit, with exponent **`k − 1 − J₀`** in |t| (all seven pieces, k = 2 and 4). Here J₀ is the
     Mellin growth at fixed (13)-channel variable: 1 for the s- and t-exchange, 2 for the u-exchange and for c2, c3, 4 for c4.
   * The joint corner is integrable iff that exponent is > −2, i.e. **J₀ ≤ k**. So in this family the pointwise zero holds even at
     J₀ = k:
     * the J₀ = 4 contact under Ψ₄: residual ≤ 3·10⁻¹² at generic angles, ≤ 2·10⁻¹⁰ down to Δφ = 0.025;
     * the J₀ = 2 exchange under k = 2 and k = 3: ≤ 6·10⁻¹³.
   * **Consequence [DV in the test family; D-sketch for Einstein/R⁴].** If the same holds for the actual terms, a light-side value of
     a J₀ = k term cannot sit in the fixed-angle integrand. The calibrated exponent rule (§3.5) says it should hold. The terms in question
     are the Einstein k = 2 pole of Result D and the g₄ term of Ψ₄.
     * Their value must then come from the collinear set a = b, or from the ν/celestial integration. D_LR drops a = b by definition.
     * For the Einstein instance of *Ψ₄, Ψ₅* (J₀ = 2 < 4) this is harmless.
     * It bears on how the light side is *defined* for J₀ ≥ 4 terms (§5).
6. **Option (c) gives no reduction [D] (§4.2).** Einstein minus the t₄ = 0 free mixture still contains, with opposite sign, the free
   theory's twist-1 higher-spin currents J_s. Each J_s exchange alone has J₀ = s ≥ 4, and the tower is not a contact. So the difference
   does not split into pieces that Lemma A or (F)+(P) treat separately.
7. **Regrading (§5).**
   * **Lemma A + D_LR stays [D]**, and the Einstein instance is still not evaluated.
   * The mechanism is now **[DV] on interacting (Witten-diagram) spinning correlators**, including a J₀ = 2 spin-2 exchange in all
     channels. Exercised: transcendental structure, logs at fixed Δ, spin-2 Regge growth, a single-trace spin-2 pole in the (13) channel,
     and the (13)/corner part beyond the free scalar.
   * **Not exercised:**
     * the graviton's own twist (τ = 1; here τ = 2);
     * the tensor structures generated by graviton exchange (here the spin of T′ comes from `T_free`);
     * the exact Einstein coefficients and external Δ = 3 (here Δ = 5).
   * Because the proof is structural and linear, what remains is a residual risk of a coefficient-level failure [D], not an untested
     mechanism.
   * Recommended line: "Lemma A + D_LR [D]; mechanism [DV] at O(1/C_T) on free fields and on spinning Witten-diagram correlators with a
     J₀ = 2 spin-2 exchange (RESULT_Q_EINSTEIN_CHECK); Einstein `⟨TTTT⟩` itself not evaluated [I]."

---

## 1. What was and was not evaluated

| item | evaluated? | how / why not |
|---|---|---|
| Einstein tree `⟨TTTT⟩` (AdS₄ graviton, Δ = 3 externals), any representation | **no** | no position-space form exists (§2.1); momentum space → light-ray complex-t not practical |
| T-exchange block + Einstein `[TT]` data in every structure | **no** | the d = 3 spinning `[TT]_{n,J}` O(G) data are not published; inversion in all helicity structures is a multi-day project |
| scalar graviton exchange (τ = 1) in AdS₄, position space | **no** | needs half-integer-index `D̄` (elliptic); integer-index truncation forces even τ (§2.2) |
| **spinning Witten-diagram correlators**: `T′ = σ ⊗ T_free`, σ-part = contact (Mellin degree 0, 2, 3, 4) and **spin-2 exchange** (τ = 2, J₀ = 2) in s, t, u | **yes** | `ein_dbar.py` (exact D-functions, Lorentzian continuation, 32–60-digit jets), `ein_run.py` |
| swapped form of `Ψ⁺_k`, k = 4, 5 (and 2, 3 as J₀ ≥ k probes), D_LR domain, T-shift | **yes** | og_free contour and `f_k`, ellipse-family contours (§2.5) |
| generic null (2,4) polarisations, 3 angle pairs; every σ-piece × every free Wick class; 1↔3 pair identities | **yes** | §3.1 |
| EE structure, Δφ = 0.05 – 3 | **yes** | §3.2 |
| controls: non-single-valued σ, non-primary T_free, weight, symmetrisation | **yes** | §3.4 |
| "no T-shift" control on the Witten part | **no** | a D-function is not a product of links; the og phase control has no analogue. Microcausality is controlled instead by the non-single-valued σ (§3.4) |
| endpoint (Regge) exponents vs J₀; collinear-corner exponents | **yes** | §3.5, §3.6 |
| scalar-external Ψ₄ on the scalar Witten diagrams (Δ = 2) | evaluated, **outside the convergence domain** | single endpoint `t⁻¹·log`; MFT s + t still cancels (validates the scalar weight) (§2.6) |
| per-ν pairing, de-shifted `Ĉ₄^{ds}` | implied | pointwise zero ⇒ zero for every ν-weight, as in OG §3.3 |
| R⁴ normalisation `−Ĉ₄^{ds} = 2g₄` | **no** | the J₀ = 4 *test* contact is annihilated pointwise (§3.5). If the same holds for the R⁴ `⟨TTTT⟩` (expected, D-sketch), its value is not in the fixed-angle integrand |

---

## 2. Representation, sources, validation

### 2.1 Literature (looked up, not recalled)
* **Raju, arXiv:1201.6452**, "Four Point Functions of the Stress Tensor and Conserved Currents in AdS4/CFT3". Tree-level `⟨TTTT⟩` of
  pure AdS₄ gravity, explicit, via BCFW-type recursion, **in momentum space**. The abstract stresses the simple momentum-space
  structure; it does not mention a position-space form (abstract checked).
* **Albayrak–Kharel, arXiv:1908.01835** (JHEP 12 (2019) 135), "Towards the higher point holographic momentum space amplitudes II:
  gravitons". AdS₄ graviton Witten diagrams **in momentum space**.
* **Caron-Huot–Li, arXiv:2102.08160**, "Helicity basis for three-dimensional CFT". It has MFT OPE data for currents and T, and tree-level
  anomalous dimensions for conserved *currents* by Lorentzian inversion. Its abstract does not cover the Einstein `[TT]` data.
* **Giecold, arXiv:1205.6444**: a universal Mellin formula for graviton exchange between *scalars* (11 hypergeometrics).
* **D'Hoker–Freedman–Mathur–Matusis–Rastelli, hep-th/9903196**: graviton exchange between scalars. The search confirmed the paper; its
  formulas were not re-derived here.
* No position-space Einstein `⟨TTTT⟩` in d = 3 was found.
* A Fourier transform of Raju's helicity amplitudes to complex light-ray kinematics at 10⁻¹⁰ is not practical (OG §5.1 item 3, agreed).

### 2.2 Why a product correlator, and why τ = 2 rather than the graviton's τ = 1
* **Lemma A uses no Einstein-specific identity** (V4FINAL §2.2). It needs:
  1. a Wightman function analytic on the two-variable contour region (light-ray analyticity + microcausality);
  2. conformal covariance with the matching light-transform weight;
  3. endpoint integrability.
* **Any finite sum of `X`-monomials × D-functions satisfies (1).** A D-function is a contact Witten diagram. Multiplying by integer
  powers of `x_ij²` keeps that: positive powers are polynomials, and `(x_ij²)^{−n}` is the Wightman 2-pt function of a decoupled
  Δ = n field.
  * So any rational, truncating Mellin amplitude gives a legitimate test correlator.
  * The product with the free `⟨TTTT⟩` is the `O(λ/N)` term of a tensor-product theory. It is microcausal and conformal, and `T′` is
    a primary.
* **Mellin dictionary (exact, used to build the pieces).**
  * `θ_ij := −X_ij∂/∂X_ij` multiplies the Mellin integrand by `δ_ij`, so a polynomial amplitude `P(δ)` is `P(θ)I_Δ`.
  * A pole `1/(δ₁₂ − 1)` is `X₁₂⁻¹I_{Δ−e₁−e₂}`.
  * Raising: `I_{Δ+e_i+e_j} = −∂I_Δ/∂X_ij`, exactly.
* **Truncation needs integer indices.** A pole at `s = τ + 2m` is absorbed by `Γ(δ₁₂)` into integer-index D-functions when 2Δ − τ is a
  positive even integer. This is D'Hoker–Freedman–Rastelli-type truncation, used here through
  `Γ(δ)/(δ − N) = (δ−1)…(δ−N+1)Γ(δ − N)`.
  * With integer Δ this forces even τ. The graviton (τ = 1 in d = 3) would need half-integer Δ, i.e. half-integer `D̄`-indices.
    Those are expected to be elliptic; not attempted.
  * So the exchange is a **massive spin-2 of Δ_e = 4 (τ = 2) between Δ = 2 scalars**. It has a single pole at s = 2 and **Regge
    intercept J₀ = 2**.
* **Residue.** The residue polynomial is the spin-2 continuous-Hahn polynomial at the leading pole, `Q = (δ₁₃ − δ₁₄)²/4 + 1/12` [D].
  * It follows from orthogonality to 1 with the measure `Γ(δ₁₃)²Γ(δ₁₄)²` at δ₁₂ = 1: `⟨y²⟩ = 1/12` for `sech²(πy)`.
  * Crossing images: `t: (δ₁₃ − δ₁₂)²/4 + 1/12` over `(δ₁₄ − 1)`; `u: (δ₁₂ − δ₁₄)²/4 + 1/12` over `(δ₁₃ − 1)`.
  * The test does not depend on this choice. Any polynomial residue is microcausal.
* **Contacts (on `I₂₂₂₂`):**

  | piece | Mellin amplitude | Regge growth |
  |---|---|---|
  | c0 | M = 1 | — |
  | c2 | σ₂ = δ₁₂² + δ₁₃² + δ₁₄² | s² |
  | c3 | σ₃ = δ₁₂δ₁₃δ₁₄ | s² at fixed t |
  | c4 | σ₂² | s⁴ (the R⁴-type control) |

### 2.3 The D-function engine (`ein_dbar.py`)
* **Base.** `I₁₁₁₁ = Φ(z, z̄)/(2X₁₃X₂₄)`, with `Φ = [2Li₂(z) − 2Li₂(z̄) + log(zz̄)log((1−z)/(1−z̄))]/(z − z̄)`, U = zz̄, V = (1−z)(1−z̄).
* **Pieces.** Every piece is a finite sum `c·ΠX_ij^{e_ij}·D_U^mD_V^nΦ`, with `D_U = U∂_U`, `D_V = V∂_V`, m + n ≤ 6.
* **Gauge (og_free).** The effective squared distances are:
  * `X₁₂ = 2t₁ − 1`, `X₃₄ = 2t₃ − 1`, `X₁₃ = 2yt₁t₃` (exact), `X₁₄ = 2t₁ + 1`, `X₂₃ = 2t₃ + 1`, `X₂₄ = 4`;
  * hence `U = (2−a)(2−b)/(8y)`, `V = (2+a)(2+b)/(8y)`, with a = 1/t₁, b = 1/t₃. These are og_free's u′, v′.
  * Only Φ needs continuation; all monomial powers are integers.

### 2.4 Validation (`cache/ein_val.txt`, `ein_run_C.log`)

| check | result |
|---|---|
| `Φ(1,1)` vs the known `D̄₁₁₁₁(1,1)` | 2.343907238689 (12 digits) |
| `I_Δ` from the operator algebra vs a direct Feynman-parameter integral, (1111), (1122), (2222), (3142), 2 Euclidean configurations | 1·10⁻¹⁴ – 3·10⁻¹⁴ (1111: 2·10⁻⁹, quadrature-limited reference) |
| Mellin consistency θ₁₂ = θ₃₄ on `I₂₂₂₂` | 0 |
| crossing: `xs(x₃,x₂,x₁,x₄) = xt(x)`, `xu`, c3, c4 invariant | 6·10⁻¹⁵ – 3·10⁻¹⁰ (c4: order-6 polydisc at the Euclidean check point) |
| scalar MFT (Δ = 2, weight `t^{Δ−2}`), k = 4: s + t = 0 (Lemma A for MFT [P]) | 1.4·10⁻⁹ / 4.5·10⁻⁹ (n = 40/60); q = +1: 2.5 |
| analytic continuation: path independence of Φ, 3 paths, 4 nodes | ≤ 5·10⁻¹⁶ |
| continuation order (t₁ first vs t₃ first along C), 49 node pairs | ≤ 1.2·10⁻¹¹ |
| contour independence (analyticity between hairpin and circle), all families | κ = 0.75 vs κ = 1: ≤ 1.2·10⁻⁹ S; κ = 0.5 (n = 100): 1.4·10⁻⁶ S (quadrature-limited) |
| n-convergence of families | n = 60 vs 80: 1.1·10⁻⁸ S |
| high-precision jets vs independent Cauchy polydisc (benign node) | agree to 10⁻¹² (order ≤ 6) |

### 2.5 Continuation route and two numerical pitfalls (recorded so they are not repeated)
* **Physical route.** The continuation runs from the Euclidean point t₁ = t₃ = ∞ (z̄ = z*):
  1. along the totally spacelike real axis to t = 1;
  2. then t₁ along its contour from t = 1;
  3. then t₃.

  On this route `log(1−z)`, `log(1−z̄)` and `Li₂` never change branch. Only `log z` and `log z̄` wind, oppositely: the s-channel
  monodromy (branch census over all 3600 nodes at Δφ = 1). So the integrand stays on the s-channel sheet, and the bulk-point (z = z̄)
  sheet is never reached. Minimum |z − z̄| along the routes ≥ 0.2 at the angles checked.
* **Pitfall 1: other routes.** Routes through t₃ = ∞ (log spirals) agree with the physical route at every contour node. A route that
  moves both variables jointly through the disk interior first, (t₁, t₃) = (0.8 ∓ 0.4i) and then outward, does **not** agree at
  endpoint nodes. So the integrand is not single-valued on all of `(ℂ∖(−∞, ½])²`.
  * What the swapped form needs is analyticity in the region swept between the double-commutator segment and the contour. That is
    confirmed by the κ-family (ellipses `t = (1 − cos θ)/2 − i(κ/2)sin θ`, κ = 1 → 0.5) and by order independence.
  * The free theory never showed this, because there everything is a product of single-pair powers.
  * Lemma C's phrase "one analytic function h(t₁, t₃)" should be read as **"analytic on the region between the hairpin and C"**. That is
    verified here for interacting correlators [DV].
* **Pitfall 2: derivatives.** Derivatives `D_U^mD_V^nΦ` are ill-conditioned in float64.
  * The log-polydisc fails where U ≈ V ≫ 1 (the light-ray endpoint).
  * Float (z, z̄)-jets lose all digits at order ≥ 3 already at |z| ≈ 5.
  * Production therefore evaluates **every node in mpmath** (32 + 4·log₁₀ max(|z|, |z̄|) digits, jet arithmetic).
  * A first production run with a float/mp split showed spurious Lemma-A residuals of 10⁻⁷ (c3) and 5·10⁻⁴ (c4) at one angle pair.
    They disappeared, down to 10⁻¹³, with all-mp jets.

### 2.6 Scalar externals are outside the convergence domain at Δ = 2
* With the scalar weight `t^{Δ−2} = t⁰`, the single endpoint behaves as `|f₄G| ~ t⁻¹ × logs`, and the u-channel identity/exchange is
  worse (`t⁻²`).
* The k = 4 values of c0, c2, … do not converge in n (c0: 8.3 → 9.9 from n = 40 to 60; xu: 4·10⁵ → 4·10⁶). So no scalar Lemma-A
  statement is made.
* The MFT s- and t-terms converge and cancel (§2.4). This validates the scalar weight and gauge.
* The spinning product inherits the light-ray polarisation suppression of `T_free` and converges with margin (§3.5).

---

## 3. Results

### 3.1 Generic polarisations, k = 4 (Ψ₄) and k = 5 (Ψ₅) (`ein_main_end.txt`)
* Residual = `|E(a,b) + (−1)^k E(b,a)|/S`, summed over all J₀ < 4 σ-pieces × all six free Wick classes (1234, 1243, 1324, MFT s, t, u).
  S = the largest single `|E[piece × class]|`.
* "Pair" = max over the 1↔3 pair identities `E[p × c](a,b) + (−1)^k E[π(p) × π(c)](b,a)`, where π swaps xs ↔ xt, 1243 ↔ 1324, s ↔ t.

| (φ₁, φ₃), Δφ | n | k | S | residual | max pair | wrong symm. | J₀ = 4 contact c4 |
|---|---|---|---|---|---|---|---|
| (−2.699, −2.326), 0.373 | 60 | 4 | 2.66·10⁻² | 2.6·10⁻¹³ | 7.8·10⁻¹⁴ | 1.17 | 2.9·10⁻¹² |
| | 60 | 5 | 2.90·10⁻³ | 7.5·10⁻¹⁴ | 4.6·10⁻¹⁴ | 0.33 | 5.2·10⁻¹⁴ |
| | 80 | 4 | 2.66·10⁻² | 2.6·10⁻¹⁴ | 1.8·10⁻¹⁴ | 1.17 | 1.6·10⁻¹³ |
| | 80 | 5 | 2.90·10⁻³ | 9.2·10⁻¹⁴ | 1.3·10⁻¹⁴ | 0.33 | 2.3·10⁻¹⁴ |
| (−2.275, 1.810), 4.085 | 60 | 4 | 7.26·10⁻⁵ | 4.1·10⁻¹³ | 1.1·10⁻¹³ | 5.51 | 1.2·10⁻¹³ |
| | 60 | 5 | 2.30·10⁻⁵ | 1.9·10⁻¹² | 5.5·10⁻¹³ | 3.37 | 2.2·10⁻¹² |
| (0.338, −0.103), 0.441 | 60 | 4 | 1.09·10⁻² | 6.7·10⁻¹³ | 1.9·10⁻¹³ | 2.73 | 2.8·10⁻¹² |
| | 60 | 5 | 7.77·10⁻³ | 4.8·10⁻¹⁴ | 2.5·10⁻¹⁴ | 1.00 | 2.7·10⁻¹⁵ |

* **Per σ-piece** (each piece with its 1↔3 image), the residual is ≤ 2·10⁻¹¹, mostly ≤ 3·10⁻¹³:
  * contacts c0, c2, c3;
  * the u-exchange xu alone;
  * xs + xt.
* **Spin-2 exchange in all three channels with J₀ = 2: annihilated pointwise.**
  * The s- and t-exchange families are individually O(S); they cancel as `xs(a,b) = ∓xt(b,a)`.
  * xu (the channel of the (13) single-trace pole) cancels against itself across free classes.
* **Largest single entries at pair 0, k = 4.**
  * xu × 1234 = −2.65·10⁻², c0 × 1234 = 1.89·10⁻², c2 × 1234 = 1.07·10⁻², c3 × 1234 = 9.8·10⁻³, xt × 1234 = 4.0·10⁻³.
  * So the cancellation is between O(S) numbers. The 1243/1324 classes (the (13) link of `T_free`) enter at O(10⁻¹ S).

### 3.2 EE structure, k = 4 (`ein_ctrl_ee.txt`)
* E(a,b) = E(b,a) for EE, so Lemma A predicts a pointwise zero.
* Values below are summed over free classes; S = the largest |piece × class|.

| Δφ | S | c0 | c2 | c3 | xu | xs | xt | \|total J₀<4\|/S | c4 / max\|c4 × class\| |
|---|---|---|---|---|---|---|---|---|---|
| 0.05 | 30.2 | 6·10⁻¹⁴ | −6·10⁻¹⁴ | 6·10⁻¹⁵ | −1·10⁻¹³ | +0.3758513 | −0.3758513 | 1.4·10⁻¹³ | 3.4·10⁻¹² |
| 0.3 | 3.27·10⁻³ | 5·10⁻¹⁶ | 7·10⁻¹⁶ | 1·10⁻¹⁶ | 8·10⁻¹⁶ | +5.231308·10⁻⁴ | −5.231308·10⁻⁴ | 7.6·10⁻¹³ | 1.8·10⁻¹² |
| 1.0 | 4.97·10⁻⁶ | −1·10⁻¹⁷ | −1·10⁻¹⁷ | −4·10⁻¹⁸ | −3·10⁻¹⁸ | +7.776357·10⁻⁷ | −7.776357·10⁻⁷ | 8.9·10⁻¹² | 9.4·10⁻¹³ |
| 2.0 | 5.09·10⁻⁷ | −3·10⁻¹⁹ | 3·10⁻¹⁹ | −4·10⁻¹⁹ | −6·10⁻¹⁹ | +4.395375·10⁻⁷ | −4.395375·10⁻⁷ | 1.4·10⁻¹¹ | 2.7·10⁻¹² |
| 3.0 | 7.49·10⁻⁷ | 2·10⁻¹⁸ | 4·10⁻¹⁸ | 4·10⁻¹⁹ | 1·10⁻¹⁸ | +7.186240·10⁻⁷ | −7.186240·10⁻⁷ | 1.9·10⁻¹¹ | 1.4·10⁻¹¹ |

* In EE, every σ-piece except the s/t exchange pair vanishes on its own (summed over free classes). The s- and t-channel spin-2
  exchanges are each between S/80 and S and cancel to 7 digits as printed.
* The J₀ = 4 contact also vanishes in EE.

### 3.3 Spin-2 exchange at k = 2, 3 (J₀ ≥ k): still annihilated pointwise (`ein_ctrl_ee.txt` (4d))

| k | n = 40 | n = 60 | n = 80 | scale |
|---|---|---|---|---|
| 2 (J₀ = k) | 2.7·10⁻¹⁰ | 1.6·10⁻¹³ | 3.9·10⁻¹³ | 3.8·10⁻⁵ |
| 3 | 1.4·10⁻⁹ | 2.8·10⁻¹³ | 5.6·10⁻¹³ | 4.0·10⁻⁵ |

The J₀ = 0 contact is ≤ 5·10⁻¹⁴ for both. §3.5 explains why the zero survives at J₀ = k.

### 3.4 Controls (pair (0.9, −1.4), n = 60; residual of the J₀ < 4 total relative to S)

| control | what it breaks | k = 4 | k = 5 |
|---|---|---|---|
| none | — | 1.4·10⁻¹¹ | 2.6·10⁻¹¹ |
| wrong symmetrisation | — | 5.90 | 5.68 |
| σ-part **not single-valued** (`Φ` without its log term) | microcausality (the Euclidean continuation has monodromy) | 1.93 (pieces 1.8 – 2.5) | 1.33 (1.9 – 2.2) |
| non-primary free factor ξ = 0 / ¼ | conformal primarity of T′ | 7.1 / 8.5 | 1.4 / 1.3 |
| weight `t^{5+q}`, q = +1 / +2 | light-ray decay (see below) | 3.3 / 4.5 | 4.1 / 5.2 |
| weight `t^{5+q}`, **q = −1** | — | **1.0·10⁻¹⁰** | **2.6·10⁻¹¹** |

**The weight control is one-sided.**
* Re-run on the free T alone (og_run.E_custom, the OG check's control pair):

  | q | k = 4 | k = 5 |
  |---|---|---|
  | −3 | 2·10⁻⁶ | 1·10⁻⁷ |
  | −2 | 7·10⁻⁷ | 2·10⁻¹⁰ |
  | −1 | 1·10⁻⁹ | 1·10⁻¹⁴ |
  | 0 | 1·10⁻¹⁴ | 9·10⁻¹⁵ |
  | +1 | 0.39 | 9·10⁻³ |

  The q < 0 values are small but not at roundoff. That is probably quadrature at the stronger endpoint singularity; it was not
  checked in n.
* So the swapped-form zero holds for every weight ≤ Δ + J − 2 and fails above it. The weight control probes the decay of the
  commutator light-ray integral (convergence at the far end), not light-transform covariance itself.
* The genuine sensitivity to conformal structure is carried by the non-primary control. The sensitivity to microcausality is carried
  by the non-single-valued control.

### 3.5 Endpoint exponents: the Regge limit of the swapped form and the role of J₀ (`ein_main_end.txt`)
`|t₁⁵t₃⁵ f_k G_σ[piece] G_free[1243]|` as t → 0 along C (Δφ = 1), exponent in |t|:

| limit | k | c0 (J₀ 0) | c2 (2) | c3 (2) | c4 (4) | xs (1*) | xt (1*) | xu (2) |
|---|---|---|---|---|---|---|---|---|
| single (t₃ fixed) | 2, 4 | +1.34 | +1.31 | +1.48 | +1.29 | +1.30 | +1.19 | +0.49 |
| joint, same side (θ₁ = θ₃ → 0) | 4 | **+3** | **+1** | **+1** | **−1** | **+2** | **+2** | **+1** |
| | 2 | +1 | −1 | −1 | −3 | 0 | 0 | −1 |
| joint, opposite sides | 2, 4 | +4.9 – +5.0 / +8.8 – +9.0 | | | | | | |

\* xs, xt: at fixed u (the (13)-channel Mellin variable) the s- and t-exchange grow like s¹.

* **The single endpoint is a lightcone limit.** Leading twist dominates and the exponent is J₀-blind; the scalar case of §2.6 diverges
  there.
* **The same-side joint endpoint is the Regge limit.** The exponent is exactly `k − 1 − J₀` for every piece (with J₀ the growth at fixed
  u) [DV].
* **Integrability.** The joint corner is a 2-real-dimensional neighbourhood. If the diagonal law holds in all same-side directions,
  the corner is integrable iff the exponent > −2, i.e. iff **J₀ ≤ k**. The opposite-side directions are much softer. Hence:
  * the J₀ = 4 contact under Ψ₄ and the J₀ = 2 exchange under k = 2 converge absolutely and are annihilated pointwise (§3.1, §3.3);
  * c4 under k = 2 (exponent −3) is the first case that diverges.
* **Calibration on the free T.** OG §3.4(d) gives exponent `k − J₀` (f₄: +3 at J₀ = 1). The product here sits one unit lower because of
  the extra σ weights.
* **Einstein (J₀ = 2, Δ = 3).** On either calibration the Ψ₄ margin is ≥ 3 units above the −2 threshold [D, by the calibrated rule;
  the Einstein exponent itself is not computed].

### 3.6 Collinear corner (D_LR) (`ein_corner.txt`)
Generic polarisations, n = 60, Δφ = 0.2 → 0.025. Fitted p in `max_class |E[piece × class]| ~ Δφ^{−p}` (log-modified power laws, so
p is a local slope).

| k | c0 | c2 | c3 | c4 (J₀ = 4) | xs | xt | xu | Lemma-A residual (J₀ < 4 total), Δφ = 0.2 … 0.025 |
|---|---|---|---|---|---|---|---|---|
| 4 | 3.20 | 3.43 | 2.76 | 3.27 | 2.66 | 2.64 | 3.29 | 4·10⁻¹³ … 1.8·10⁻¹¹ (c4: 6·10⁻¹² … 1.7·10⁻¹⁰) |
| 2 | 1.40 | 1.49 | 1.08 | divergent* | 2.66 | 2.62 | 3.06 | 8·10⁻¹³ … 2.1·10⁻¹⁰ |

\* c4 at k = 2 has joint-endpoint exponent −3 (§3.5). Its n = 60 values (3·10⁶ → 5·10⁹) are quadrature artefacts of a divergent
integral and are not used.

* **Pointwise zero down to Δφ = 0.025.** The pointwise zero of every J₀ < 4 total persists into the corner, at both k. The families grow
  by up to three decades over this range.
* **Family growth at k = 4.** Every family, the J₀ = 4 contact included, grows like `Δφ^{−(2.6…3.4)}`. That is the free-scalar twist-1 law
  `Δφ⁻³` of RESULT_Q_OG_CHECK §3.6, here log-modified. **No piece shows an enhanced corner singularity.**
* **What this means for D_LR [DV].** The single families of the J₀ ≤ 2 pieces are no more singular at a = b than in the free check.
  The total is zero at every Δφ ≥ 0.025. So D_LR's "the collinear set carries nothing" is supported on the Witten family exactly as on
  the free scalar.
* **A caveat on integrability.** Whether `Δφ^{−3}` times T′'s own celestial kernel is integrable was not computed: the j′ = 0 kernel of
  OG §3.3 is for Δ = 3, not for T′ (Δ = 5). Only the comparison with the free law is claimed.
* **At k = 2 the exchanges stand out [NE].** The exchange families (J₀ = 2 = k) grow faster than the contacts (p ≈ 2.6–3.1 vs
  1.1–1.5). This is weak evidence that J₀ = k content concentrates towards a = b.
* **The J₀ = 4 contact at k = 4 does not stand out [I].** Where its (F)+(P) value would sit (an a = b distribution, or the ν
  integration) is not determined by these data.

---

## 4. The reduction argument (option (b)) and the difference with the free mixture (option (c))

### 4.1 What Lemma A uses, and where each input is now tested

| input of V4FINAL §2.2 | Einstein `⟨TTTT⟩` (tree) | tested on |
|---|---|---|
| (i) Wightman function analytic between the hairpin and C, with the four orderings as boundary values (Lemma C) | bulk locality of tree Witten diagrams; the s-channel sheet only | free scalar, fermion, GFF (trivially); **Witten diagrams: order independence 1·10⁻¹¹, contour (κ) independence 1.2·10⁻⁹ (§2.4–2.5)** |
| (ii) microcausality at fixed angles | [D] | **non-single-valued control fails at 1.3 – 2.5 (§3.4)** |
| (iii) T primary, light-transform weight Δ + J − 2 | Δ = 3, J = 2 | non-primary controls (OG, FH, §3.4); the weight control is one-sided (§3.4) |
| (iv) endpoint convergence | J₀ = 2 [Q: Costa–Hansen–Penedones] | **Regge exponent `k − 1 − J₀` measured (§3.5); margin ≥ 3 units for Ψ₄** |
| (v) relabelling 1↔3 | identical T | pair identities per channel and class (§3.1) |
| logs at fixed external Δ (O(G) anomalous dimensions) | yes | GFF ∂/∂Δφ (Δ moving); **D-functions (Δ fixed) [DV]** |
| single-trace spin-2 pole in all channels incl. (13) | τ = 1 graviton | free T (τ = 1, J₀ = 1); **massive spin-2 (τ = 2, J₀ = 2) [DV]** |
| (13) link / collinear corner | twist-1 (13) channel | free scalar; **Witten product (classes 1243, 1324 and xu) [DV]**; corner §3.6 |

What is **not** covered by any computed correlator:
* the exact Einstein coefficients (graviton τ = 1 with Einstein `⟨TTT⟩` couplings, EH contact, external Δ = 3);
* the spinning tensor structures generated by graviton exchange itself (here the spin of T′ comes from `T_free`).

Because Lemma A is linear and its proof is structural, this is a residual-risk statement, not a gap in the mechanism:
* Einstein would have to violate (i)–(v);
* the tests show the zero survives every *kind* of term that Einstein contains.

### 4.2 Option (c): Einstein − (t₄ = 0 free mixture) [D]
* Both have the same `⟨TT⟩` and `⟨TTT⟩`, so single-T exchange cancels in every channel.
* But the free mixture contains, at the same order 1/C_T, the twist-1 higher-spin currents J_s (s = 4, 6, …) in every channel.
  * Each J_s exchange alone has Regge spin s ≥ 4. Only their infinite sum has J₀ = 1.
  * The J_s exchanges are not contacts, so (F)+(P) does not apply to them.
* The O(G) double-trace data differ in every (n, J) as well.
* Hence `Ψ[𝒢_E] = Ψ[𝒢_E − 𝒢_free]` exactly, and the difference is no easier than 𝒢_E. It does not decompose into pieces that Lemma A
  (J₀ < 4 each) or (F)+(P) (contacts) handle separately.
* **(c) gives no reduction.** The difference is itself a J₀ ≤ 2 microcausal object, so Lemma A covers it only as a whole, which is the
  original claim.

---

## 5. Regrading

| claim | before | after |
|---|---|---|
| Lemma A + D_LR on J₀ < 4 tree terms (the Result Q input) | [D]; mechanism [DV] on free fields; Einstein unchecked | **[D]**; mechanism **[DV] also on interacting spinning Witten-diagram correlators** (contacts J₀ = 0, 2; spin-2 exchange J₀ = 2 in s, t, u; logs at fixed Δ; Li₂) |
| Einstein instance | [I], leverage 10¹⁰–10²² exposed | **[I] for the exact function**. The qualitative ingredients are now [DV], except the graviton's τ = 1 and its own tensor structures. The residual risk is a coefficient-level failure, which the structural proof excludes [D] |
| Lemma C ("one analytic h") | [D] | **[DV] in the form "analytic between hairpin and C"** for Witten diagrams; not analytic on all of `(ℂ∖cut)²` (§2.5) |
| twist-1 corner / D_LR at O(1/C_T) | [DV] (free scalar) | **[DV] also on the Witten family**: families `Δφ^{−(2.6–3.4)}` like the free twist-1 law; total zero to Δφ = 0.025 (§3.6) |
| weight control as a test of light-transform covariance | stated in OG §3.4 | **qualified**: one-sided (weights below Δ + J − 2 also give 0) |
| pointwise annihilation needs k > J₀ | stated (CLS-D) | in the test family the swapped-form pointwise zero holds for **J₀ ≤ k**, because the joint corner is integrable [DV]. CLS-D's k > J₀ is sufficient, not sharp, for the pointwise form |

**Implications to record in the unified statement (not about the Einstein Ψ₄, Ψ₅ instance):**
* **[D-sketch]** Light-side values of J₀ = k terms live at the collinear set or in the ν integration. This would cover the Einstein
  k = 2 pole of Result D and the R⁴ value `−Ĉ₄ = 2g₄` of Q.
  * Basis: the J₀ = k test pieces are annihilated pointwise (§3.1, §3.3), and the exponent rule of §3.5 extends this to the actual
    terms.
  * The pointwise vanishing has not been computed on the actual Einstein or R⁴ `⟨TTTT⟩`.
* D_LR drops a = b by definition. So **D_LR cannot be applied uniformly**:
  * for J₀ < k terms (Einstein in Ψ₄, Ψ₅) it is consistent;
  * the J₀ = k values must be supplied by (F)+(P) (Q) or by the H7 pole (D), as the documents already do.
* The unified theorem should say explicitly that the light side is **D_LR on J₀ < 4 terms plus term-wise (F)+(P) on J₀ ≥ 4 terms**, and
  not "the D_LR functional applied to the full correlator". Read the second way, the D_LR value of the R⁴ term would be 0 (§3.1 c4
  column).
* This is a statement about the definition of the light side. It does not change the Einstein zero.

---

## 6. What remains
1. **The Einstein function itself** (cost: days). Two routes:
   * **Scalar-external graviton exchange in AdS₄ with half-integer Δ.** It truncates (Δ = 5/2: poles at s = 1, 3). Needed:
     `D̄_{Δ_i}` with half-integer indices (elliptic `K`, `E`) continued along the physical route. The engine here is ready apart from the
     base function. This would test the τ = 1 graviton pole, but with scalar externals (outside the convergence domain unless Δ ≥ 3,
     cf. §2.6), or inside a σ ⊗ T_free product as here.
   * **Spinning:** weight-shifting / spinning differential operators acting on such seeds to build the actual graviton
     `⟨TTTT⟩` structures, or Lorentzian inversion of the crossed-channel T block in all d = 3 helicity structures (Caron-Huot–Li basis)
     plus a convergent `[TT]` resummation on the contour.
2. **The R⁴ normalisation in position space.** It is now clear that it is not in the fixed-angle integrand (§3.5). It must be
   extracted at the collinear set (a = b) or in the ν integration: one more reason (F)+(P) is a separate hypothesis.
3. **Lemma C's precise domain.** The integrand of interacting correlators is analytic between the hairpin and C but not on all of
   `(ℂ∖cut)²`. A short derivation from the s-channel OPE (z, z̄ ∉ [1, ∞) on the s-channel sheet) would upgrade §2.5 to [D].
4. **The weight control** should be restated in RESULT_Q_OG_CHECK / RESULT_Q_FERMION_HEFT as one-sided (§3.4).

---

## 7. Reproduce

```bash
cd ads4/code/resultq
(ulimit -v 3000000; OMP_NUM_THREADS=1 python3 ein_run.py val)            # §2.4 validation, ~1 min  -> cache/ein_val.txt
(ulimit -v 3000000; OMP_NUM_THREADS=1 python3 ein_run.py main end)       # §3.1, §3.5, ~30 min     -> cache/ein_main_end.txt
(ulimit -v 3000000; OMP_NUM_THREADS=1 python3 ein_run.py ctrl ee)        # §3.2–3.4, ~45 min       -> cache/ein_ctrl_ee.txt
(ulimit -v 3000000; OMP_NUM_THREADS=1 python3 ein_run.py cauchy)         # §2.4 contour/order, ~30 min -> cache/ein_cauchy.txt
(ulimit -v 3000000; OMP_NUM_THREADS=1 python3 ein_run.py corner)         # §3.6, ~45 min          -> cache/ein_corner.txt
```

| file | content |
|---|---|
| `ein_dbar.py` | `Φ` and branch tracking (`State`, `step`); term algebra for `ΠX^e·D_U^mD_V^nΦ` (`I1111`, `raise_pair`, `theta`, `apply_poly`); test pieces (`pieces`: c0, c2, c3, c4, xs, xt, xu); ellipse contours and the physical continuation route (`contour_ellipse`, `continue_ellipse`); derivative tables: `phi_derivs` (log-polydisc, Euclidean use only), `phi_derivs_jet` (float jets, diagnostic), `phi_derivs_mixed` / `_mp_jet_derivs` (production, mpmath jets); `grid_D_ellipse_hp` |
| `ein_run.py` | sections `val` (§2.4), `main` (§3.1), `ee` (§3.2), `ctrl` (§3.3–3.4), `end` (§3.5), `cauchy` (§2.4–2.5), `corner` (§3.6) |
| `cache/ein_*.txt`, `ein_run_*.log` | recorded output quoted here |

**Conventions.**
* `T′ = σ ⊗ T_free`, weight `t⁵`, `f_k` and the gauge of og_free.
* "Residual" is relative to the largest single `|E[piece × class]|` over the J₀ < 4 pieces.
* The angle pair (a, b) is `x₁ = −t₁z(a)`, `x₃ = t₃z(b)`.
* The routes `continue_ab` / `grid_D` (log spirals + log-polydisc) are kept for the diagnostics of §2.5 only. Do not use them for
  production.
