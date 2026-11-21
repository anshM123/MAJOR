# Result Q, V₄-final: H_LT from the definitions, the ν-order circularity, H_νout, and the final hypothesis list of the upper bound

*21 Nov 2026, V₄-final agent. Unrefereed, not committed. No existing document or code file was edited; no `bj_*` file was touched or run.
New code: `ads4/code/resultq/v4e_lt.py`, `v4e_swap.py`; outputs `resultq/cache/v4e_lt.txt`, `v4e_swap.txt` (git-ignored). Read-only
sources: CLS-D 2311.04271 `draft.tex` and KS-D 1805.00098 `draft.tex` (LaTeX sources already in the session scratchpad from earlier
stages; treated as data). Read-only inputs: `functionals_q4_deshift.json` (the `R_j`), `functionals_q7_xprec.json` (XP weights),
`v4d_core.py`, `v4d_og.py`. Compute: about 2 CPU-minutes, at most 2 threads, under 1 GB.*

Grades: **[P]** proved (including exact algebra); **[D]** derived; **[D-sketch]**; **[DV]** derived and numerically verified; **[NE]**
numerical evidence; **[Q]** quoted from the literature; **[M]** model; **[I]** open.

---

## 0. Answer in brief

1. **H_LT is not an extra hypothesis. It is what the CLS-D and KS-D definitions say [D]; given those definitions, its consequence is [P].**
   * **What the definitions fix.** CLS-D place operators 1 and 3 on null rays from a common point x₀, through `δ(x₁₀²)δ(x₃₀²)θθ`, with
     their polarisations equal to the ray directions `x₁⁺₀`, `x₃₀`. They state that 1 and 3 are spacelike "in the entire integration
     range", which excludes the collinear set. KS-D define the light transform with the transformed operator on its real null contour.
     The only complex continuation KS-D use is of the ray parameter along its own null direction (`ζ₀ = Im(v)z`).
   * **Consequence on the complexified light-transform manifold [P, sympy].** The manifold is
     `{x₀ + a z₁, x₀ + b z₃}`, with complex x₀, a, b and complex null z. On it the numerator of ⟨T₃T₁⟩ vanishes identically, and
     `x₁₃² = −2ab z₁·z₃`. So the u-channel identity times f₄ is **identically zero** under every regularisation that stays on the
     manifold:
     * the cone-preserving shift `t → t + iε`;
     * KS-D's continuation of the ray parameter;
     * complex celestial angles.
   * **Why the tube is not a regularisation of the CLS-D integral.** The timelike tube `x₁ → x₁ + iεe` fails in three separate ways:
     * it leaves the manifold: `x₁₀² = 2iεt₁ + ε² ≠ 0`, so it violates CLS-D's `δ(x₁₀²)`, and z₁ is no longer the ray direction;
     * its divergence comes from giving f₄'s `1/x₁₃²` pole an iε. That is the only `(13)` singularity of f₄, and Lemma C assigns no
       Wightman iε to it;
     * under it, CLS-D's own vanishing theorem fails on MFT, even for k = 2 in commutator form (V₄′ table).
   * **Prescription independence of everything else [DV].** Every other corner piece is absolutely integrable. The O(G) twist-1
     structures give the same value under the tube limit and under the light-ray (pointwise) value, to 1–4·10⁻⁵ (§1.4).
   * **What the bound would lose under the tube reading.** Nothing, after subtracting one universal, theory-independent term: the
     `C_T²` u-channel identity contact. Without that subtraction the functional is not defined at all, so there would be no bound
     (not a weaker one).
2. **The ν-inside / ν-outside circularity is resolved.** The light side never needs either order.
   * **Pointwise-in-angle Lemma A [D].** For every tree term with J₀ < 4 (Einstein, E×W³, MFT), take the integrand, integrate it over
     t₁, t₃ (this step is ν-independent) and symmetrise it under 1↔3. The result vanishes **identically for φ₁ ≠ φ₃**. The ingredients
     are microcausality, CLS-D's endpoint convergence (J₀ < k) and Lemmas B and C. None of the corners is used: neither the large-η
     corner nor the collinear one.
   * At φ₁ = φ₃ there is nothing, by H_LT.
   * For MFT this is V₄′'s exact `R_s(a,b) = −R_t(b,a)` [P].
   * **So the J₀ < 4 light side is the zero function.** It pairs to 0 with any ν-weight in any order of the ν-, angle- and t-integrals.
   * **V₄‴'s per-ν statement is superseded.** V₄‴ evaluated the light side family by family, and V₄″ used the tube. Those were the only
     sources of the obstructions (i) and (ii) in V₄″. Here the light side is never evaluated family by family.
   * **Answer to (a): yes.** ν outside on both sides is a consistent definition. It is literally what Stages 3–7 computed: heavy
     image = `∫χ·(per-ν image)` per state, absolutely convergent. It was re-verified here to 10⁻⁸ (§2.3).
   * **Answer to (b).** The **only** interchange left is the heavy OPE sum against the (celestial-corner, ν) integration.
     * **Fubini over (ν, states) needs the raw corner.** Per state, `∫|χ|Σ|R_j||P| ≍ e^{5β/2}`, while the ν-outside image is
       `≍ e^{−β/2}β` (ratio `5·10⁶` at β = 6). So Fubini over (ν, states) holds iff the raw per-ν heavy sum converges: J₀ < ½ under K1′.
       Call the sufficient hypothesis **H_ip** (impact-parameter decay of the heavy tower).
     * **The kernel-level ν-inside form is not a position-space object.** The position-normalised weight grows like `e^{κν}`, with
       κ = 1.60–1.73 > π/2 (§2.4). As an analytic functional it would need celestial-angle rotations beyond the (13) light-cone cut.
     * **So under K1′ alone the interchange is a Proposition K for the heavy-normalised kernel.** This is the same class as Result D's
       I14(b). Call it **H_K**, grade **[D-sketch]**.
3. **H_νout is not independent.**
   * **J₀ < 4 terms.** The summed per-ν light side is identically 0, so it is trivially bounded [D].
   * **No family-by-family bound exists [DV].** A single family grows like `e^{(π−1.26)ν}` in heavy normalisation:
     * `|I_t|N₀/N = 10⁻², 7, 2.9·10⁸, 1.2·10¹⁵` at ν = 0.5, 8, 20, 30;
     * its ν-outside de-shift diverges with the cutoff: `8.6·10⁻¹⁰ → 1.1·10²` (ν_c = 5 → 40).
     * Only the exact cancellation survives.
   * **J₀ ≥ 4 EFT terms.** H_νout is just (P). The absolute moments `∫|χ|Σ_j|R_j||ν−ij|^{2n}` are finite (2.2, 83, 1.4·10⁴, 3.6·10⁶
     for n = 0–3, value weight). So any polynomial light side has an absolutely convergent ν-outside de-shift [P given (P)].
   * **H_νout is deleted.**
4. **O(G) Einstein ⟨TTT⟩ u-channel T block with the spinning kernel: not computed [I], and no longer load-bearing.**
   * The O(G) light-family cancellation is the pointwise Lemma A (microcausality). An explicit check would only be a cross-check.
   * The corner integrability of every twist-1 structure is [D]. Its prescription independence is now [DV].
5. **The bound changes nothing numerically.** Its hypothesis list shrinks:
   * V₄‴, H_LT and H_νout are replaced by the pointwise Lemma A [D] together with H_LT as a definition [D].
   * The one remaining load-bearing Regge input is H_K [D-sketch]. Its strong sufficient form is H_ip, which turns it into [D given H_ip].

---

## 1. H_LT from the definitions (`v4e_lt.py` → `cache/v4e_lt.txt`)

### 1.1 What the sources fix [Q]

**CLS-D (2311.04271).** The four places that matter:
* **Light transform** (eq. `light_transform_definition`): `L[O](x,z) = ∫dα (−α)^{−Δ−J} O(x − z/α, z)`. The operator sits on the null line
  from x in the direction of its own polarisation; in embedding space it is `O(Z − αX, −X)`.
* **Defining identity** (eq. `superconvergence_vdirection`): `∫dv₁dv₃ f ⟨φ₄[φ₃(0,v₃,y₃), φ₁(0,v₁,y₁)]φ₂⟩ = 0`, "since
  `x₁₃² = y⃗₁₃² > 0`, the two operators are spacelike separated in the entire integration range".
* **Spacetime form** (eqs. `scalarfunctional_plus_deltafunctions`, `scalarfunctional_plus_final`, and the spinning analogue below
  `spinningfunctional_definition`): `δ(x₁⁺₀²)δ(x₃₀²)θθ`, and "the polarizations of O₁, O₃ are fixed to be `x₃₀`, `x₁⁺₀` (null vectors
  due to the delta functions)".
* **Subtracted functional** (§3.4.2, the paragraph introducing it): `Ψ_{k,ν,f}[𝒢] = Ψ⁺[f𝒢] − (1↔3)`, again because "points 1 and 3 are always spacelike in
  the range of integration".

**KS-D (1805.00098, Lemma `lemma:light`, eq. `lemmaeplsilonprescription`).**
* `L[O]|Ω⟩ = 0` is proved with the other operators shifted by `−ikεê₀`.
* The light-transformed operator is at a **real** point of its contour.
* The analytic continuation is `v → v + i Im v` **along z**: "`ζ₀ = Im(v)z` will remain in the future of `ζ₁`".
* The commutator forms (`commutatorexamples`) follow from wrapping the same v-contour.

### 1.2 Why a definition is needed at all [P]

* **Where it is needed.** On the integration manifold the map `(x₀, α₁, z₁, α₃, z₃) ↦ (x₁, x₃)` is a submersion except at `z₁ = z₃`. There
  the image of its differential in the `x₁₃` direction is `span{z, ∂_φz}`. Its annihilator is the covector `(k, −k)` with `k ∝ z`.
* **Why the restriction is undefined there.** That covector lies in the wave-front set of the Wightman function ⟨T₃T₁⟩ at null
  separation `x₁₃ ∝ z`. By Hörmander's criterion the restriction of ⟨T₃T₁⟩ is therefore not defined at the collinear set.
* **The point is not removable.** This persists after integrating over x₀, x₂ and x₄. (V₄″ made the same observation with x₀ fixed.)
* **What CLS-D do instead.** They do not need the restriction: their domain is `x₁₃² > 0`, which excludes `z₁ = z₃`.

**Definition D_LR.** `Ψ⁺_{k,ν}[X]` is the Lebesgue integral of the pointwise integrand over that domain. It is defined whenever the
integrand is absolutely integrable at the collinear set.

### 1.3 Proposition LT [P, exact algebra, `v4e_lt.py` (1)]

On the complexified light-transform manifold, with light-ray polarisations, the following hold identically.

| quantity | value on the manifold |
|---|---|
| null vectors | `z = (1, (1−m²)/(1+m²), 2m/(1+m²))`, m ∈ ℂ |
| `x₁₃²` | `4ab(m₁−m₃)²/((1+m₁²)(1+m₃²)) = −2ab z₁·z₃` |
| `N = (z₁·z₃)x₁₃² − 2(z₁·x₁₃)(z₃·x₁₃)` | **0** |
| cone-preserving `t → t + iε` | `N = 0`, `x₁₀² = 0` |
| tube `x₁ → x₁ + iεe` | `N = ε²(1 + cos Δφ)`; `x₁₀² = 2iεt₁ + ε² ≠ 0` |
| f₄ as a function of `Y = x₁₃²` | `−A(A − Y s₂₄)(A + Y s₂₄)/(4Y·B)`: a simple pole at Y = 0 only |

**Consequences.**
* `W_u·f₄` is identically zero on the manifold away from `z₁·z₃ = 0`, and nothing else in f₄ is singular in `(13)`.
* So **D_LR** and every regularisation that stays on the manifold assign it exactly 0:
  * KS-D's ray-parameter continuation;
  * cone-preserving shifts;
  * complex x₀;
  * complex Lorentz rotations of z₁, z₃.

### 1.4 The tube is not a regularisation of the CLS-D integral, and nothing else depends on the prescription

**Three independent reasons against the tube.**
1. **It leaves the support.**
   * It violates `δ(x₁₀²)` (table above), and the polarisation is no longer `x₁₀`.
   * Operatorially, `T(x + iεe) = e^{−εH}T(x)e^{εH}`. So the tube is the light transform of T₁ based at the different point
     `x₀ + iεe`, not CLS-D's common base point.
   * With `f = 1` it still gives 0: `L[T₁]|Ω⟩ = 0` survives conjugation.
2. **It invents an iε.** The divergence `ε^{−1/4}` exists only through f₄'s `1/x₁₃²`. The control without that pole goes like
   `ε^{+3/4}` (V₄′ §3). Lemma C gives the `(12)`, `(14)`, `(23)` and `(34)` poles the Wightman iε of their ordering. The `(13)` pole has
   none, because `x₁₃² > 0` on the whole real domain. The tube supplies one by fiat.
3. **It contradicts CLS-D on MFT.**
   * In the tube the commutator form of raw k = 4, and even of raw k = 2, diverges on MFT at interior large-η corner points (V₄′ §3
     table: `ε^{−0.30}`).
   * MFT is unitary with J₀ ≤ 1. So under the tube reading CLS-D's statement "Ψ_{k,ν,f} vanishes on physical four-point functions"
     would be false for the simplest such correlator.
   * Under D_LR the commutator integrand is pointwise 0 and the statement holds.

**Prescription independence of the integrable pieces [DV].** The O(G) u-channel twist-1 structures `U_mn` of `v4d_og` are
absolutely integrable at the corner, with margin `r^{1/2}`. Compare:
* the tube limit `ε → 0` (fit `G₀ + Aε^{1/4+iν/2} + Bε^{1/4−iν/2}`, ε = 10⁻⁴ … 3·10⁻⁷);
* the light-ray pointwise value (soft cutoff `(x²/(x²+d²))²`, d → 0).

| structure | ν = 2: tube `G₀` | ν = 2: light-ray `G₀` | rel. diff | ν = 6: rel. diff |
|---|---|---|---|---|
| (1,1) | `−1.88613 + 1.89315i` | `−1.88610 + 1.89316i` | 1.1·10⁻⁵ | 3.7·10⁻⁶ |
| (2,2) | `16.9757 − 17.0378i` | `16.9749 − 17.0384i` | 4.2·10⁻⁵ | 1.3·10⁻⁵ |
| (0,0) | `1.88607 − 1.89320i` | `1.88610 − 1.89316i` | 2.1·10⁻⁵ | 4.1·10⁻⁶ |

* The differences are at the level of the fit residuals.
* By dominated convergence the same holds for every absolutely integrable piece: the O(1) regular families (`|Δφ|^{9/2}`), all
  twist ≥ 1 u-channel terms (homogeneity lemma, V₄″ §2) and heavy blocks.
* **The prescription dependence is confined to the u-channel identity.** That term is `⟨T₄T₂⟩⟨T₃T₁⟩` exactly. It is fixed by C_T and is
  the same in every theory.
* The branch of the half-integer powers is fixed by the operator ordering. The comparison uses the same branch (`+i0`) on both sides.

### 1.5 Verdict on H_LT

| claim | grade |
|---|---|
| the collinear set needs a definition (WF ∩ conormal ≠ ∅) | [P] |
| CLS-D/KS-D define the functional on the domain `x₁₃² > 0` with operators on their null contours (D_LR, light-ray reading) | [D] (reading of explicit statements; neither paper discusses `z₁ = z₃`) |
| `N ≡ 0` and `W_u f₄ ≡ 0` on the complexified light-transform manifold | [P] |
| the tube violates `δ(x₁₀²)`, gives f's (13) pole an iε not in the definition, and contradicts CLS-D's theorem on MFT | [P] + [DV] (V₄′) |
| every absolutely integrable corner piece is prescription-independent; twist-1 structures agree to ≤ 4·10⁻⁵ | [D] + [DV] |
| loss under the tube reading: none after subtracting the universal `C_T²` identity contact; no functional without it | [D] |

**H_LT → definition D_LR, grade [D].** It is no longer a hypothesis about the theory.

---

## 2. The ν-order circularity (`v4e_swap.py` → `cache/v4e_swap.txt`)

### 2.1 The three integrations and where each is needed

The de-shifted functional is `∫dν χ(ν) ∫dφ₁dφ₃ C_ν(φ₁,φ₃) ∫dt₁dt₃ [t³t³ f₄ 𝒢⁺]`, divided by `N_com(ν)` (heavy normalisation).
* The t-integrals are **ν-independent**. All of Lemmas B and C, the vacuum-annihilation contour arguments and the edge residues live
  there, at fixed `(φ₁, φ₃)`.
* ν enters only through the celestial structure C_ν, which is 1↔3 symmetric for j′ = 0 (antisymmetric for j′ = 1, Ψ₅).
* There are two corners, and both are `φ₁ → φ₃`:
  * the collinear corner at the light-cone edge, `t = ½`;
  * the large-η Regge corner at `t < ½`, fixed ρ = 2t (V₄′ §0.1).

### 2.2 Pointwise-in-angle Lemma A [D]

Let T be a tree term with J₀ < 4 (Einstein, E×W³, MFT). Define `E⁺_T(a,b) := ∫dt₁dt₃ t³t³ f₄ ⟨[T₄,T₃][T₁,T₂]⟩_T`, with 1 at celestial angle
a and 3 at b. Then for every a ≠ b:

    E⁺_T(a,b) + E⁺_T(b,a) = 0        (j′ = 0; for j′ = 1, Ψ₅: E⁺(a,b) − E⁺(b,a) = 0).

**Proof.**
* **Step 1: the commutator vanishes at fixed angles.** At fixed a ≠ b, points 1 and 3 are spacelike for all t₁, t₃. So `f⟨4[3,1]2⟩ = 0`
  pointwise by microcausality, which is [D] for tree correlators.
* **Step 2: the t-integrals converge.** At fixed angles the only boundaries are the ray endpoints and the edges t = ½. The endpoint
  convergence is CLS-D's condition `k > J₀` [Q]. The endpoint rows and the edges are Lemmas B and C [P].
* **Step 3: from the commutator to the double commutator.** The t-contour arguments of Lemma C turn `∫dt f⟨4312⟩` into `E⁺(a,b)`. The same
  arguments turn `∫dt f⟨4132⟩` into the 1↔3 image.
* **Step 4: relabelling.** For identical operators, with f antisymmetric under `u′ ↔ v′`, the 1↔3 image equals `−E⁺(b,a)`.

∎ Neither corner is used.

**Check [P].** For MFT, `E⁺ = R_s + R_t` and `R_s(a,b) = −R_t(b,a)` exactly (V₄′ §2.1, all 81 tensor components). For EE,
`R_s + R_t ≡ 0` already pointwise, so `|I_s + I_t| = 0` at every ν (§2.5).

**Consequences.**
* **At a = b nothing is added (D_LR).** So for J₀ < 4 terms the symmetrised light-side integrand is **the zero function on the whole
  celestial torus**.
* Its pairing with C_ν, with χ, or with any kernel-level total kernel is 0 in every order of integration. It needs no corner
  convergence and no ν-decay.
* The per-ν, per-family evaluation of V₄‴ (and the J-resummation of `[TT]₀` at the corner) is never needed.
* This supersedes Stage 12's per-ν Lemma A. That lemma needed V₄, i.e. convergence at the corners. The pointwise form needs only the
  fixed-angle convergence that CLS-D already require.

### 2.3 (a) ν outside on both sides is a consistent definition [P/D]

Define `Ĉ_out[X] := ∫dν χ(ν) Ψ̂_ν[X]` on pieces X for which `Ψ̂_ν[X]` converges absolutely for every ν and `∫|χΨ̂_ν[X]| < ∞`.

| piece | in the domain? | value | grade |
|---|---|---|---|
| heavy block (leading Regge) | yes, margin `e^{−(π/2)ν}` | image `I(β) = Σ_j∫χR_jP_{−½+i(ν−ij)}(cosh β)` = `P_{−½}(cosh β)` (value), `−D(β)` (deriv) | [P] given `R_j`, χ; re-verified to 10⁻⁸ at β = 0.5–6 |
| J₀ < 4 tree terms | yes (zero integrand) | 0 | [D] (§2.2) |
| J₀ ≥ 4 EFT terms (g₄, `|ĝ₃|²`, g₅, …) | the position integral diverges (J₀ ≥ k); the per-ν value is given by (F)+(P) | `Σ_j∫χR_j p((ν−ij)²)` = `p(0)` (value), `p′(0)` (deriv) | [P given (P)]; absolute convergence §3 |

This is exactly what Stages 3–7 computed (per-state ν-outside images) and what the light side needs. **No piece requires ν inside.**

### 2.4 (b) The remaining interchange: heavy OPE sum against (corner, ν)

The sum rule needed is `Σ_O Ĉ_out[G_O] = −Ĉ_out[EFT]`. It comes from the pointwise identity
`Σ_O E⁺_{O,sym}(a,b) = −E⁺_{EFT,sym}(a,b)` for a ≠ b (§2.2 removes the J₀ < 4 terms). To get there, the state sum must be interchanged
with the celestial-corner and ν integrations.

**Route 1: Fubini.** It needs `Σ_O ∫dν|χ||Ψ̂_ν[G_O]| < ∞`. Per state (`v4e_swap` (1), XP value weight):

| β | 0.5 | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|---|
| `I(β)` (= `P_{−½}`) | 0.9846 | 0.9409 | 0.7957 | 0.6234 | 0.4641 | 0.3337 | 0.2341 |
| `A(β) = ∫|χ|Σ|R_j||P|` | 2.30 | 5.59 | 58.9 | 703 | 8.5·10³ | 1.04·10⁵ | 1.27·10⁶ |
| `A/|I|` | 2.3 | 5.9 | 74 | 1.1·10³ | 1.8·10⁴ | 3.1·10⁵ | 5.4·10⁶ |

* `A ∝ e^{2.50β}`: this is the raw `|P_{5/2+iν}|`. The derivative weight gives the same rate.
* **So Fubini over (ν, states) ⇔ the raw per-ν heavy sum converges at the large-η corner.** Under K1′ that needs **J₀ < ½**
  (Stage 3 §1.5), which K1′'s J₀ ≤ 1 does not give.
* Sufficient hypothesis **H_ip**: the heavy tower's leading-Regge impact-parameter density, weighted by `e^{5β/2}`, is summable.
  * This is physically natural at tree level with Δ_g ≥ 300: each heavy exchange has profile `Π_{Δ−1}(β) ~ e^{−(Δ−1)β}`, and the large-b
    tail is graviton-dominated, i.e. light.
  * It is **not** implied by uniform Regge bounds.
* Given H_ip, the per-ν sum rules hold and Fubini gives the de-shifted rule: **[D given H_ip]**.

**Route 2: kernel level (Stage 4, J₀ < 7/2).** The total kernel `∫w(ν)C_ν dν` is not a position-space function.
* **The weight grows.** `w = χN_com(0)/N_com(ν)` grows at rate κ (`v4e_swap` (2)):

  | window | [10,20] | [20,30] | [30,40] | [50,60] | [70,80] | [90,100] | [100,110] |
  |---|---|---|---|---|---|---|---|
  | κ | 2.10 | 1.73 | 1.67 | 1.62 | 1.61 | 1.600 | 1.596 |

  κ tends to π/2 = 1.571 from above, while `|C_ν| ~ ν^{−1/2}`.
* **The analytic-functional reading crosses the cut.** As an analytic functional, the pairing needs the celestial angle rotated by
  θ > κ > π/2. That moves `x₁₃² ∝ Δφ²` past `arg x₁₃² = π`, the (13) light-cone cut. Beyond it the orderings 4312 and 4132 differ, so
  microcausality is no longer available. The weight sits exactly at the edge of the microcausal sector.
* **What Stage 4's "J₀ < 7/2" actually establishes.** It is a statement about `Σ_O` of per-state ν-outside images (the series converges).
  Equating that series to the EFT side is a **Proposition K for the heavy-normalised kernel**.
  * This is the class of Result D's I14(b): D's de-shifted `Ĉ′_EE` has the same `e^{π|ν|}` caveat, Stage 4 §2.3.
  * It must be formulated in the impact-parameter (heavy-normalised) representation, where the de-shifted kernel exists pointwise
    (`P_{−½}(cosh L)` on H₂). It cannot be formulated as a position kernel.

**Precise statement that makes the bound consistent.**

> Result Q's upper bound uses the ν-outside functional Ĉ_out throughout. Its light side is 0 (J₀ < 4 terms, pointwise Lemma A + D_LR)
> plus (F)+(P) for the EFT terms, both [D]. Its heavy side is the convergent series of per-state ν-outside images [P/D]. The single
> load-bearing interchange is **H_K: `Σ_O Ĉ_out[G_O] = −Ĉ_out[EFT]` under K1′ with J₀ < 7/2** [D-sketch]. H_K follows by Fubini from
> **H_ip** (raw per-ν corner convergence of the heavy tower) [D given H_ip].

**Grades.**
* The circularity (ν inside on the heavy side vs ν outside on the light side) is resolved: [D].
* What is left is a pure heavy-side Regge interchange: H_K [D-sketch], or H_ip.

### 2.5 Why V₄″/V₄‴ saw obstructions

Both of V₄″'s obstructions, (i) the large-ν growth and (ii) the ν = 0 junction, are properties of **individual** light families or of
the tube-regularised `W_u`. Example (`v4e_swap` (4)): the MFT `[TT]₀` family I_t alone.

| ν | 0.5 | 2 | 5 | 8 | 12 | 16 | 20 | 25 | 30 |
|---|---|---|---|---|---|---|---|---|---|
| `|I_t|` (position norm.) | 1.0·10⁻² | 1.4·10⁻³ | 1.1·10⁻³ | 2.1·10⁻⁴ | 6.3·10⁻⁶ | 9.5·10⁻⁸ | 9.4·10⁻¹⁰ | 2.0·10⁻¹² | 3.2·10⁻¹⁵ |
| `|I_t|N₀/N` (heavy norm.) | 1.1·10⁻² | 3.2·10⁻³ | 0.12 | 7.2 | 2.0·10³ | 7.2·10⁵ | 2.9·10⁸ | 5.6·10¹¹ | 1.2·10¹⁵ |

* Decay rate ≈ 1.26 < κ. The truncated ν-outside de-shift of I_t alone runs `8.6·10⁻¹⁰, 2.3·10⁻⁹, 4.1·10⁻⁵, 7.5·10⁻³, 1.1·10²` for
  ν_c = 5, 10, 20, 30, 40.
* `I_s + I_t = 0` exactly at every ν: the integrand cancels pointwise.
* The bound never uses a single family.

---

## 3. H_νout

* **Reduction [D].** In the original wording H_νout is "the per-ν light side, summed over families, is polynomially bounded in heavy
  normalisation". For J₀ < 4 terms that sum is identically 0 (§2.2), not merely bounded. For J₀ ≥ 4 EFT terms the per-ν light side is
  `Σ_j R_j(ν)p((ν−ij)²)` by (P), which is polynomial in heavy normalisation since `R_j ~ ν⁶`.
* **Explicit bound [P given (P)].** `|L̂_ν| ≤ Σ_terms |c_term| Σ_j |R_j(ν)||p_term(ν−ij)|`. Its ν-outside de-shift converges absolutely.

  | | n = 0 | 1 | 2 | 3 |
  |---|---|---|---|---|
  | exact moments `Σ_j∫χR_j(ν−ij)^{2n}`, value | 1 | 0 | 0 | 0 |
  | exact moments, deriv | 0 | 1 | 0 | 0 |
  | absolute moments `∫|χ|Σ_j|R_j||ν−ij|^{2n}`, value | 2.15 | 82.6 | 1.36·10⁴ | 3.55·10⁶ |
  | absolute moments, deriv | 0.173 | 4.04 | 533 | 1.37·10⁵ |

  * So for the R⁴ entry, `|−Ĉ^{ds}|_{R⁴}| ≤ 2|g₄|·2.15`, while the exact value is `2g₄`.
  * Higher contacts enter with `Λ^{−2j}` times moments of degree `2j`. The finite absolute moments show that each term is absolutely
    convergent. The tower sum is hypothesis 1 ("no conspiracy").
* **Obstruction to anything weaker [DV].** No family-by-family or partial-sum bound exists (§2.5). H_νout holds only through the exact
  pointwise cancellation.
* **Verdict:** H_νout is **deleted** as an independent hypothesis. It is pointwise Lemma A [D] plus (P) [D-sketch, spinning].

---

## 4. O(G): the Einstein ⟨TTT⟩ u-channel T block with spin [I, not load-bearing]

* **Not computed.** That would need the d = 3 Einstein ⟨TTT⟩ (fixed combination of the two parity-even structures) in light-ray
  polarisation, the u-channel T block at the corner, and the `[TT]₀` `γ_J`, `δp_J` in this gauge.
* **Why it is no longer load-bearing.** The O(G) light-family cancellation is the pointwise Lemma A, which is microcausality at fixed
  angles (§2.2). It is not an independent Einstein-specific identity.
* **What is covered.**
  * Every twist-1 (13) structure is corner-integrable with margin `r^{1/2}` (homogeneity lemma [D], V₄″ §2).
  * It is prescription-independent (§1.4) [DV].
  * So no O(G) contact term can appear at φ₁ = φ₃.
* **What an explicit evaluation would test.** The spinning O(G) analogue of V₄′'s exact `R_s(a,b) = −R_t(b,a)`. That is a consistency
  check of Lemmas B and C at O(G).

---

## 5. Final hypothesis list for Result Q's upper bound (ĝ₃ = 0, tree level, Δ_g = 300/1000)

| # | hypothesis | role | grade |
|---|---|---|---|
| 1 | A1 (twist gap for all J), A2 (tree level); local EFT light sector, no conspiracy in the contact tower | setup | assumption |
| 2 | KKSDZ + K1′/K2/K3 for the de-shifted functionals (J₀ < 7/2) **and H_K** (heavy OPE sum ↔ (corner, ν) interchange for the heavy-normalised `Ĉ₄^{ds}`, `∂Ĉ₄^{ds}`, `Ĉ₅^{ds}`) | the only load-bearing Regge input | **[D-sketch]**, Result D I14(b) class; **[D given H_ip]** |
| 2′ | (optional strong form) **H_ip**: raw per-ν corner convergence of the heavy tower (impact-parameter decay beyond K1′) | makes 2 a Fubini argument | [D-sketch] physically (tree level, Δ_g ≥ 300); not implied by K1′ |
| 3 | **Pointwise-in-angle Lemma A** for J₀ < 4 tree terms (Einstein, E×W³, MFT): microcausality (D), endpoint convergence J₀ < k (Q), Lemmas B, C (P), `J₀(E×W³) ≤ 3` (D) | light side = 0 in every ν-order; replaces V₄, V₄′, V₄″, V₄‴ | **[D]** (MFT instance [P]) |
| 3a | **H_LT → definition D_LR** (CLS-D's domain `x₁₃² > 0`, operators on their null contours) | the collinear set carries nothing | **[D]** (consequence [P]); tube reading loses nothing after the universal `C_T²` subtraction |
| 3b | O(G) corner integrability, prescription independence | no O(G) contact | [D] + [DV] |
| — | ~~H_νout~~ | deleted: = 3 + (P) | — |
| 4 | Spinning (F) + (P) | EFT light-side values (`2g₄`, …) and H_νout for J₀ ≥ 4 terms | [D-sketch] |
| 5 | Continuum de-shift exists: `W[χ] = Φ` for all t (C4 value/deriv; C5 3-shift, spinning j′ = 1 coefficients) | heavy images `P_{−½}`, `−D(β)` beyond t = 30 | [NE] to t = 30 (Stage 7); C5 spinning [I] as before. Strip analyticity of `χR_j` [NE] is consistency only, no longer load-bearing on the light side |
| 6 | Inherited from Result D: H7, EE admissibility, c₁, t₄ bound (or carry c₁t₄) | | as in D |
| 7 | H_smear^Q | smeared parts | [NE + M] |
| 8 | H_fwd^Q (unmeasured forward channels; real-section ≡ crystal for EE) | | [NE] / [I] |
| 9 | H_OO^Q (parity-odd J = 0, 2 channels equal ++) | | [I] |
| 10 | HYB accuracy; β > 30 | | [NE] |

**Changes relative to RESULT_Q_STATEMENT / V₄‴.**
* Items 3 (V₄ → V₄‴), H_LT and H_νout collapse into the pointwise Lemma A [D] plus the definition D_LR [D].
* The circularity is resolved.
* The residual Regge content is isolated as H_K (or H_ip).
* Hypothesis 5 no longer needs strip analyticity for the light side.
* **Numerically nothing changes:** `0 ≤ g̃₄Λ⁶/(8πG) ≤ 17.24·log(2Λ/3) − 19.22/−19.27`.
* **Overall grade of the upper bound: unchanged at "conditional".** The weakest load-bearing links are H_K [D-sketch], (F)+(P)
  [D-sketch], H_fwd^Q/H_OO^Q [I] and the smeared parts [M].

---

## 6. Files (`ads4/code/resultq/`)

| file | content | output | time |
|---|---|---|---|
| `v4e_lt.py` | (1) sympy: N ≡ 0 on the complexified light-transform manifold; tube leaves the cone; f₄'s only (13) pole. (2) tube-limit vs light-ray value for the twist-1 structures (1,1), (2,2), (0,0) at ν = 2, 6 | `cache/v4e_lt.txt` | 46 s |
| `v4e_swap.py` | (1) per-state ν-outside image vs absolute integral, β = 0.5–6 (Fubini cost). (2) growth rate of the position-normalised weight. (3) exact and absolute moments (H_νout for (P) terms). (4) single family I_t vs summed light side | `cache/v4e_swap.txt` | 64 s |

Run: `OMP_NUM_THREADS=2 python3 v4e_lt.py; python3 v4e_swap.py`. Both import `v4d_core` (and `v4e_lt` imports `v4d_og.U_mn`). Neither
modifies any existing file.
