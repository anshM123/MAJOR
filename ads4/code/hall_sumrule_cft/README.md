# CFT3 thermal Hall-channel sum rule: derivation from the T × T → T OPE

Written 27 Nov 2026 by a research agent. **Not refereed.** Nothing outside this folder was edited, and nothing was
committed. (The scripts import `../sumrule_freefield/cs_fermion_probe.py` and `../scope_sumrule/hall_uv_sign.py`
read-only, with bytecode writing switched off.)

Grades:
* **[D]** derived here, complete given the cited inputs;
* **[C]** computed by a script in this folder (output `out_*.txt`);
* **[Q]** quoted from a source whose arXiv LaTeX was read (for Jain–John 2107.00695, read for this note);
* **[H]** hypothesis;
* **[I]** inferred, not derived.

| script | content | runtime |
|---|---|---|
| `lemmas.py` | L1 the energy matrix in the helicity basis; L2 the rank/parity lemma for contact terms; L3 Fourier transforms of the T and Δ = 2 tails; L4 assembly | 10 s |
| `cs_fermion_exact.py` | **exact** large-frequency expansion of the Ghosh–Mazumdar N = ∞ CS-fermion thermal correlator (sympy, using the gap equation) and a numerical check against direct evaluation | 20 s |
| `bulk_hall_exact.py` | **exact** large-Ω constant for the bulk `Q = C C C̃` coupling, plus a real-frequency check of the convergent (spectral) form | 3 s |

---

## 0. Verdict

**The claim holds with the stated factor:** `|S_H| = |α_T|/32`, and in complex form

    S_E − 3/8 + i σ S_H = (t₄ + i α_T)/32 ,        |S_E − 3/8 + i σ S_H| = τ/32 ≤ 1/8,

Here σ = ±1 is a sign that is the same for every CFT (it depends only on conventions; §5). The equivalent convergent form

    S_H = J_H/ε = (1/πε) ∫_{−∞}^{∞} η_H(ω) dω ,     η_H(ω) = −Im δh(ω)/ω,

holds under the hypotheses of §6. Overall grade: **[D]**.

There are two independent exact checks [C]:
* **Large-N CS-fermion, all λ.** The Hall constant is exactly `ε sin(πλ)/8 = ε α_T/32`, not just to 1e-5. The full complex helicity constant is `−ε/8 − ε(t₄ − iα_T)/32` exactly.
* **AdS₄ bulk, `CCC̃`.** `S_H = −45 λ_Q` exactly. The real-frequency integral gives −45.0002.

The collaborator's sketch is right in substance. Four points need adding (§2.4–2.6, §4):
1. The rotation gives a factor ±i only because the homogeneous and odd structures are both supported on all-equal helicities.
2. The relative normalization of α_T and t₄ must come from the collider matrix, written in the helicity basis.
3. The even UV constant itself depends on the contact-term scheme. So the factor 1/32 must be read off the physical t₄-slope of `S_E`, not off "the even OPE constant".
4. Three side statements in the earlier notes need correcting (§7):
   * temperature-dependent parity-odd local terms do exist, but only through Δ = 2 scalars;
   * the CS-fermion `b₁ q₃` term is ψ̄ψ exchange, not an arbitrary contact term;
   * the 1e-5 mismatch was a fit artefact.

---

## 1. Setup and conventions

* Euclidean coordinates `x = (τ, x, y)`. Zero spatial momentum, `p = (Ω, 0, 0)`.
* Shear operators `T₁ = (T_xx − T_yy)/2`, `T₂ = T_xy`.
* The plane-integrated correlator is `G_AB(τ) = ∫d²x ⟨T_A(τ,x) T_B(0)⟩_β`, and `G_AB(ω_n)` is its Fourier transform.
* The retarded function satisfies `G^R(iΩ) = G_E(Ω)` for Ω > 0.
* Decomposition: `G_AB = g δ_AB − h ε_AB` (with `ε₁₂ = +1`), i.e. `G^R = [[g, −h],[h, g]]`. Helicity channels: `G_± = g ± i h`. `δ` means thermal minus vacuum.
* Reality and Bose symmetry (`G_AB(τ) = G_BA(−τ)`) give:
  * `g(τ)` is even in τ, and `h(τ)` is odd;
  * `h_E(ω_n)` is real and odd;
  * `h^R(−ω) = h^R(ω)*`.
* Definitions: `J_E = lim δg(iΩ)`, `J_H = lim_{Ω→+∞} δh(iΩ)`, `S_H = J_H/ε`, and `S_E = (δg(0) − J_E)/ε`.
* Thermal one-point function. Lorentzian `⟨T_μν⟩ = diag(ε, P, P)` with `P = ε/2`. In Euclidean form, `⟨T_αβ⟩_β = −(3ε/2)(u_α u_β − δ_αβ/3)`, with `u = τ̂`.
* CDP's d = 3 result, with `t₄` in the HM normalization: `S_E = 3/8 + t₄/32` (`../scope_sumrule` §1) [Q + C].
* The energy matrix of CDP 1707.03007 is `(E/4π)[[1 − t₄/4, α_T/4],[α_T/4, 1 + t₄/4]]`.

## 2. Derivation

### 2.1 Thermal OPE at large Euclidean frequency [standard; Q]
The behaviour of `δG(iΩ)` at large Ω is set by the small-τ singularities of the plane-integrated `δG(τ)`:
* The region `|x| ≳ β` gives a part that is smooth and periodic in τ, so its Fourier transform decays faster than any power.
* The region `|x| ~ |τ| ≪ β` is governed by the OPE.

So

    δG_AB(iΩ) ≃ Σ_O ⟨O⟩_β Ĉ^O_AB(Ω̂) |Ω|^{3−Δ_O}  (+ descendants, which have no vev).

This is the CDP / Caron-Huot (0903.3958) / Iliesiu et al. (1802.10266) framework.

### 2.2 Which operators can give a constant [D]
* **Odd spin vanishes.** The rotation by π in the (τ, x) plane maps `S¹_β × R²` to itself, so `⟨O_{τ…τ}⟩_β = (−1)^ℓ ⟨O_{τ…τ}⟩_β`. Odd-spin one-point functions therefore vanish at zero chemical potential, with no symmetry assumption. Only even ℓ contributes.
* **Unitarity in d = 3** requires `Δ ≥ ℓ + 1` for ℓ ≥ 1. Hence:
  * ℓ = 2 gives a constant (Δ = 3) only for a conserved spin-2 operator, i.e. T (by H1: a unique stress tensor);
  * ℓ ≥ 4 has Δ ≥ 5, so its terms decay;
  * scalars give a constant only at Δ = 3 (marginal). Scalars with Δ < 3 give growing terms (§4.2).

### 2.3 The most general T term at zero spatial momentum, and its Fourier transform [D; C: L3]
* **Symmetry.** The term is linear in `⟨T⟩_β ∝ ε(u u − 1/3)`. At `p ∥ u` it can depend only on `Ĉ^{αβ}_AB(p̂) p̂_α p̂_β`. Rotations about τ̂ make the 2×2 matrix of the form `a δ_AB + b ε_AB`. Bose symmetry together with oddness under `p → −p` forces `b ∝ sgn Ω`. Hence the most general T contribution is

      δG_AB(iΩ)|_T = J_E δ_AB − sgn(Ω) J_H ε_AB ,      J_E, J_H ∝ ε  (P = ε/2 enters only through tracelessness).

* **Position space** (`lemmas.py` L3). The T term is `|x|^{−3} f(x̂)`. Integrating it over the plane at fixed τ gives `F(sgn τ)/|τ|`, where `F(±)` is the integral of f over the upper or lower hemisphere.
  * **Even channel:** `∫ e^{iωτ}|τ|^{−1}` has a pole (it gives `log|ω|`). Since no theory has such a log, the T term in `δg(τ)` has no `1/|τ|` tail: the angular average f₀ vanishes. `J_E` is then the coefficient of a `δ(τ)` distribution, and that is exactly where contact terms `δ³(x)·T` live. **So `J_E` depends on the scheme; only `S_E` is physical.**
  * **Odd channel:** `PV∫ e^{iωτ}/τ dτ = iπ sgn ω`, which is finite. The Hall constant is a separated-point tail, `δh(τ) ∝ J_H/(πτ)` for `|τ| ≪ β`, multiplied by the factor i that is usual for Euclidean parity-odd correlators. **So `J_H` is a non-local, scheme-independent quantity.** By L2 the odd angular function has no l = 0 part, so no log appears.

### 2.4 Extracting the OPE coefficient from the vacuum ⟨TTT⟩ [D]
`Ĉ_AB(p̂)` is the coefficient of the non-analytic `|q|³` term in the soft limit `q → 0` of

    ⟨T_A(p) T_B(−p−q) T_γδ(q)⟩ = [analytic in q] + Ĉ_AB^{αβ}(p) ⟨T_αβ(q) T_γδ(−q)⟩ + O(|q|⁴).

The third leg is soft, and `⟨TT⟩(q) = c_T |q|³ Π⁽²⁾(q̂)`.
* Descendants of T give `|q|⁴` and higher.
* Semi-local terms with `x₃ → x₁` or `x₃ → x₂` are analytic in q.
* Semi-local terms with `x₁ → x₂` are polynomial in p. They are the `δ³(x)·T` contact terms of the OPE, which by L2 are parity-even only.
* The ranges of `Π⁽²⁾(q̂)` over all q̂ span the traceless symmetric tensors. So this functional determines `Ĉ·⟨T⟩_β` completely.
* **The extraction is a linear functional acting on the soft leg only.**

### 2.5 Helicity locking [Q: Jain–John 2107.00695; D for the consequence]
Jain–John write `⟨TTT⟩ = c_nh⟨⟩_nh + c_h⟨⟩_h + c_odd⟨⟩_odd`, with `⟨⟩_h = (FB − FF)/2` and `⟨⟩_nh = (FB + FF)/2` (their eqs. `hnh12`, `gnspin1`). The parity-odd part of a three-point function of exactly conserved currents is purely homogeneous. Two of their results are used [Q]:
* **eq. `s1s2s3sh`:** `⟨T⁻T⁻T⁻⟩_h = (c_even + i c_odd) × (one fixed spinor structure)`, and `+++` is its complex conjugate. Every other helicity component is pure nh.
* **eq. `3ptmapaf`:** `odd = R_k(even_h)` and `R_k(odd) = −even_h`, where `R_k = k_k⁻¹ ε^{μ α k_k}` acts on one index of any one leg. So `R_k² = −1`, and R has eigenvalue `±i` on helicity ±. Their footnote to `3ptmapaf` gives the collinear version: `sgn(k₃)·even_h = odd` for momenta along the third axis. That is exactly the factor `sgn Ω` found in §2.3.

Applying R to a **hard** leg and the extraction of §2.4 to the **soft** leg act on different legs, so they commute. The hard leg has momentum `Ω τ̂`, so `R = sgn(Ω)·(±ε_AB)` on `(T₁, T₂)`. The h-part of Ĉ is ∝ δ_AB by parity, so the odd part is `± sgn(Ω) ε_AB` times the same number. Equivalently, in helicity channels:

    J_∓ − J_nh = K_th (c_even ± i c_odd) ,   K_th real  (real because c_odd = 0 must give a parity-even, real J).

### 2.6 Normalization of α_T relative to t₄ [D; C: L1]
In the helicity basis `v_± = (e_a ± i e_b)/√2`, the CDP energy matrix is `(E/4π)·[[1, −(t₄ + iα_T)/4],[−(t₄ − iα_T)/4, 1]]` (L1).
* The diagonal is independent of t₄ and α_T. It comes from nh: the Einstein value t₄ = α_T = 0 is pure nh.
* The helicity-flip element is a linear functional on the equal-helicity sector, which by §2.5 is ∝ `(c_even + i c_odd) K_coll`.
* K_coll is real, because a parity-even theory must give real t₄.

Hence `t₄ : α_T = c_even : c_odd`, with a universal real constant. **The collider definition of α_T (disc `t₄² + α_T² ≤ 16`) is exactly the one that makes the thermal combination `t₄ ± iα_T`:**

    J_∓ = J_nh + κ (t₄ ± i σ α_T) ε ,      κ real, σ = ±1 universal.

### 2.7 Fixing κ from the physical even slope [Q + C; D for the scheme argument]
`S_E = (δg(0) − J_E)/ε` is free of contact terms: by L2 every T-dependent even local term is a constant, so it cancels in the difference.
* In the covariant scheme, `δg(0) = P` (CDP, hydrodynamics), so `S_E = 1/2 − J_nh/ε − κ t₄`. Then `S_E = 3/8 + t₄/32` gives **κ = −1/32** and `J_nh = ε/8` (CDP's `J/ε = 1/8 − t₄/32`).
* For κ to be the slope of the ⟨TTT⟩_h soft functional, `J_nh` must not depend on t₄. This holds because:
  * the parity-even ⟨TTT⟩ in d = 3 has no ultralocal terms (degree-3 polynomial, rank 9, which is odd, so by L2 it would need ε);
  * its semi-local terms are fixed by Ward identities in terms of `c_T`. This is Jain–John's nh = (FB + FF)/2, which contains all Ward-identity terms.
* Independent confirmations of the slope −1/32 [C]:
  * free scalar and free Dirac (`../scope_sumrule`, `../sumrule_freefield`: S_E = 1/2 and 1/4);
  * Einstein (3/8);
  * the **exact** CS-fermion family (§3.1), where the λ-dependent part of `J_E` is exactly `−t₄ε/32` and the λ-independent offset is a scheme constant.

### 2.8 Result for the UV constant [D]
From §2.6 and §2.7, `J_H = −σκα_Tε` up to the convention sign:

    J_H = σ α_T ε / 32 ,      S_H = σ α_T / 32 .

### 2.9 IR side and the convergent form [D]
* **`δh^R(0) = 0`.** There are two arguments:
  1. The static response matrix is the Hessian of the free energy (the Kubo–Mori inner product is real and symmetric on Hermitian operators), so its antisymmetric part vanishes.
  2. Directly, `h_E` is odd in `ω_n`.

  The bulk confirms this: `Π(Ω)/Ω → 108`.
* **Dispersion.** Assume `δh − J_H` is analytic in the upper half-plane and `O(|ω|^{3−Δ*})` there, real axis included (H5). Kramers–Kronig at ω = 0 then gives `Re δh(0) − J_H = −(1/π)∫η_H dω`, so `S_H = (1/πε)∫_ℝ η_H dω`.
* **Convergence.**
  * At large ω, `η_H = O(ω^{2−Δ*})`, which is integrable because Δ* > 3. In holography Δ* = 6 (the [TT] double trace), and the bulk tail is `η_H w⁴ → 1.55×10³ ε`.
  * At small ω, `η_H(ω) → η_H`. At finite c_T there can be a logarithmic long-time tail, which is integrable.
* **No positivity.** η_H(ω) is not sign-definite. In the bulk, η_H(0) = +54ε, yet the integral is −45·π·ε. So the sum rule alone gives no bound.

## 3. Checks

### 3.1 Large-N CS-fermion, exact in λ [C: `cs_fermion_exact.py`]
The inputs are the Ghosh–Mazumdar correlator (eq. `strf1a`, as coded in the earlier probe) and three closed forms, with `z = −e^{−μ+iπλ}`:
* `∫ψ = −Im Li₂(z)`;
* `∫ψE² = −Im[μ²Li₂ + 2μLi₃ + 2Li₄](z)`;
* `li = 2i Im[μLi₃ + Li₄](z)`.

All three are checked against quadrature to 1e-32. The key step is the gap equation, `Im Li₂(z) = −(π/2)(1 − λ)μ_F`. With it, sympy proves:

    G(q₃) = ⟨J₋₋J₊₊⟩ = 2 G_+ = A₃q₃³ + A₁q₃ + A₀ + O(1/q₃),   A₂ = 0,
    A₃ = sin(πλ/2) e^{iπλ/2}/(64πλ)              (Re: c_T ∝ sin πλ/(πλ); Im: vacuum odd contact)
    A₁ = μ_F² sin(πλ/2) e^{iπλ/2}/(8πλ)          (ψ̄ψ, Δ = 2)
    A₀ = (e^{iπλ} − 1) ε/4  ⇒  A₀/(2ε) = −1/8 − (t₄ − iα_T)/32   (t₄ = −4cos πλ, α_T = 4 sin πλ).

* **Hall constant:** `S_H = sin(πλ)/8 = α_T/32` **exactly**, for all λ. The earlier ratio of 1.00001 was fit bias. The quantity `|G − A₃q³ − A₁q − A₀|·q₃` stays bounded up to n = 400, which confirms the closed form numerically.
* **Shear constant:** `−1/8 − t₄/32` exactly. The slope matches CDP. The offset −1/8 (against CDP's +1/8) is the λ-independent contact convention, as for free Dirac.
* **Locking for a second operator.** The ψ̄ψ term carries the phase `e^{iπλ/2}`, so `b₁/a₁ = tan(πλ/2)` exactly. This is the same mechanism as §2.5, applied to `⟨TTO⟩ = c_e h_even + c_o h_odd` (JJ eq. `sso`). The phase is half that of the T term, consistent with MZ's θ = πλ/2 [I]. Hence the Hall `b₁q₃` term is a Δ = 2 exchange, analytic because `sgn(Ω)|Ω|^{3−Δ} = Ω` at Δ = 2. It is not a free contact term.
* **Sign in these conventions:** `S_H = +α_T/32` in GM's `J₋₋` orientation combined with CDP's `α_T = +4 sin πλ`. The two papers' orientation conventions were not cross-checked [I].

### 3.2 Bulk `L ⊃ λ_Q C C C̃` on the AdS₄ brane [C: `bulk_hall_exact.py`]
* **UV constant.** At large Ω the response is localized at `z ~ 1/Ω`, and the constant is the O(mass) term, i.e. T exchange.
  * The pure-AdS regular mode is `Φ₀ = 1 + Ωz`.
  * The O(m) source is `S₁ = −144Ω³z⁴`; everything else is O(Ω⁻⁴).
  * The unique polynomially bounded solution is `X = −30z³ − 30Ωz⁴ − 18Ω²z⁵`.
  * This gives `Π(∞) = −90`, so **`S_H = −45 λ_Q` exactly**, for every Ω. The earlier fit gave −44.99.
* **Map to α_T.** With `λ_Q = −α̃₃/6` and `|α_T| = 240|α̃₃|` (`DICTIONARY_T4.md` §4, graded [D] there), `|S_H| = 7.5|α̃₃| = |α_T|/32`. The bulk mechanism behind it (the duality rotation `C³ → C_±³` on all-plus/all-minus vertices) is the image of §2.5.
* **Spectral form on the real axis.**
  * `(1/πε)∫_ℝ η_H dω = −45.0002`, against the exact −45.
  * `η_H(0)/ε = 54.0000`.
  * The tail is `∝ ω⁻⁴`, as expected for Δ* = 6.

  This confirms the hypotheses H3 and H5 in a holographic example.

## 4. Contact terms, other operators, convergence, sign (task 3)

### 4.1 Contact terms [D; C: L2]
At finite T, a local term in `⟨T_A T_B⟩_β` must be a local functional of the background. On flat `S¹ × R²` such a functional is T-independent. So T-dependence enters local terms only as `c·p^k·⟨O⟩_β`: the seagulls `δT/δg` and the `δ(x)O` OPE contact terms, with k = 3 − Δ_O.

The term has rank `4 + ℓ_O + k`. An SO(3)-invariant tensor of even rank contains no ε, and one of odd rank always contains ε (L2: ranks 4–7 checked: 3, 6, 15, 36 invariants with parity +, −, +, −). Hence:
* **O = T (k = 0, rank 6): no parity-odd local term exists.** `J_H` cannot be shifted by any scheme choice or seagull.
* The only T-dependent parity-odd local terms come from **scalars with Δ = 2** (k = 1, linear in ω). Such a term does not touch the constant, but it does change `η_H(ω)` by a constant. This is CS-fermion's `b₁` term.
* The gravitational Chern–Simons contact term `κ_g ω³` is T-independent, so it drops out of δ.
* This corrects `../scope_sumrule` §4.1, which said "no temperature-dependent contact term exists in the odd channel". That is true for constants, but not for terms linear in ω.

### 4.2 Operators other than T [D]
In the Hall channel a scalar O contributes `c_o ⟨O⟩_β sgn(Ω)|Ω|^{3−Δ}`, where c_o is its ε-structure (helicity-locked) coefficient.

| Δ | effect on J_H | remedy |
|---|---|---|
| 3 (marginal) | shifts the constant; the odd FT has no log | excluded by H4 |
| 2 < Δ < 3 | J_H limit does not exist | subtract `c_o⟨O⟩(−iω)^{3−Δ}`-type term |
| = 2 | analytic `Ω`, i.e. constant shift of η_H (CS-fermion at N = ∞) | subtract constant from η_H |
| 1 < Δ < 2 | growing `Ω|Ω|^{2−Δ}` | subtract |
| > 3 | decays, makes ∫η_H converge | none |

* Spin ≥ 4 has Δ ≥ 5, so its contributions decay. In CS-matter at N = ∞ this includes the slightly broken higher-spin currents.
* Double traces [TT] have Δ ≥ 6.
* So for CS-fermion the integral form needs the ψ̄ψ subtraction (and the N = ∞ zero-frequency structure). Its UV constant is clean (§3.1).

### 4.3 Convergence [D, given H5]
See §2.9. The integrand decays as `ω^{2−Δ*}` once every Δ ≤ 3 scalar with `c_o⟨O⟩ ≠ 0` is absent or subtracted.

### 4.4 Relative sign [D for universality; I for the value]
σ is fixed by the orientation (`ε_{τxy}`), the helicity labels, and CDP's sign of α_T (equivalently of p_T). Every step is linear with fixed kinematic factors, so **σ is the same for every CFT**. Fixing it once, in any example, fixes it everywhere.
* In GM + CDP conventions, σ = +1, provided their orientations agree, which was not checked.
* In the bulk, `S_H = −45λ_Q`, but `α_T = ∓240α̃₃` has an unfixed sign (`DICTIONARY_T4` §4).

The physically meaningful sign is `S_H` relative to `η_H` in the same h. That ratio is not fixed by the OPE: `S_H = −(5/12)·4πη_H/s` holds only for the pure `CCC̃` brane (`../scope_sumrule` §4.3), and is not universal.

## 5. Comparison with the collaborator's argument
* **"Extraction commutes with R_k on a hard leg": correct** (§2.5). The two operations act on different legs.
* **It needs two facts:**
  1. h and odd are supported only on (−−−) and (+++), so R produces a pure phase ±i and no mixing with nh. This is JJ `s1s2s3sh`.
  2. α_T is defined by the collider, whose helicity-flip element is ∝ `t₄ + iα_T` (L1). Without fact 2, "same magnitude" would hold only up to the unknown ratio of JJ's `c_odd` to the paper's α_T.
* **The diagonal t₄ coefficient is not itself an observable** (§2.3, §2.7). The statement that survives is about the t₄-slope of the physical `S_E`, and that slope is what transfers to `S_H`.

## 6. Theorem

**Theorem (Hall-channel thermal sum rule, CFT3).** Assume the following.

* **H1** — a unitary CFT₃ with a unique conserved spin-2 operator T (no decoupled sectors).
* **H2** — the KMS state on R² at temperature T > 0, with zero chemical potentials. Parity and time reversal may be broken.
* **H3** — the thermal OPE governs the large-Euclidean-frequency expansion of the plane-integrated `⟨T_A T_B⟩_β`, and the large-|x| region contributes only terms that decay faster than any power.
* **H4** — no scalar O in T × T with `Δ_O ≤ 3` has both `⟨O⟩_β ≠ 0` and a non-zero parity-odd (`ε`-structure) ⟨TTO⟩ coefficient. For `Δ_O < 3` it suffices that such terms are subtracted.

Then:

1. `δh(iΩ) → J_H` as Ω → +∞, with **`J_H = σ α_T ε/32`**. Here σ = ±1 is universal, α_T is normalized by the CDP energy matrix (`t₄² + α_T² ≤ 16`), and ε is the energy density.
2. `J_H` has no contact-term ambiguity.
3. `δh(iΩ) − J_H = O(Ω^{3−Δ*})`, where Δ* > 3 is the smallest dimension among the remaining contributing operators.
4. If in addition H5 holds (`δh^R − J_H` is analytic in the upper half-plane and `O(|ω|^{3−Δ*})` uniformly up to the real axis), then `S_H := J_H/ε = (1/πε)∫_ℝ η_H(ω)dω` with `η_H = −Im δh^R/ω`, and the integral converges absolutely.
5. Together with CDP (whose hypotheses are H1–H4 with "parity-even" in place of "odd" in H4, plus `δg(0) = P`): `S_E − 3/8 + iσS_H = (t₄ + iα_T)/32`, so `τ = 32 |S_E − 3/8 + iσS_H|`.

| step | content | grade |
|---|---|---|
| 1 | thermal OPE asymptotics (H3) | [Q] standard (CDP; Caron-Huot 0903.3958) |
| 2 | only T and Δ ≤ 3 scalars give non-decaying terms; odd spins vanish by the π-rotation | [D] |
| 3 | most general T term `J_E δ − sgnΩ J_H ε`; even part is δ(τ), odd part is a 1/τ tail | [D] + [C] L3 |
| 4 | no parity-odd local term ∝ ⟨T⟩; odd ones only from Δ = 2 scalars | [D] + [C] L2 |
| 5 | OPE coefficient = soft `|q|³` functional on the soft leg | [D] |
| 6 | helicity locking; R_k on a hard leg commutes with step 5 | [Q] JJ 2107.00695 + [D] |
| 7 | collider normalization: same real constant for t₄ and α_T | [D] + [C] L1, from the quoted CDP matrix |
| 8 | κ = −1/32 from the physical slope of `S_E`; `J_nh` independent of t₄ | [Q] CDP + [C] free/Einstein/CS-exact + [D] |
| 9 | `δh(0) = 0`; dispersion; convergence | [D] (H5 assumed) |
| 10 | value of σ | universal [D]; value conventional, not fixed [I] |
| check | CS-fermion, all λ, N = ∞: `S_H = α_T/32` exactly; complex constant `∝ t₄ − iα_T` exactly | [C] |
| check | bulk `CCC̃`: `S_H = −45λ_Q` exactly; ∫η_H gives −45.0002; `|S_H| = |α_T|/32` given the [D] dictionary factor 240 | [C] |

**Weakest links:**
* H5 on the real axis is a standard assumption, verified here only holographically.
* The CS-fermion check involves the CDP 1707.03007 values of `(t₄, α_T)(λ)` [Q via `../sumrule_freefield`].
* Steps 6–7 rest on JJ's spinor-helicity theorem, which was read but not re-derived.

## 7. Corrections to earlier notes
* `../scope_sumrule/README.md` §4.1: "no temperature-dependent contact term exists in the odd channel at p = 0" is true for **constants**. Parity-odd T-dependent terms **linear in ω** exist exactly when a Δ = 2 scalar has a vev (§4.1).
* The same note, §4.2 and §4.4: `|S_H| = |α_T|/32` "verified only holographically / CFT derivation [estimated]" is now [D], with exact CS and bulk checks.
* `../sumrule_freefield/README.md` §3:
  * the `b₁q₃` Hall term is the ψ̄ψ (Δ = 2) exchange, phase-locked to the shear `a₁|q₃|` term (`b₁/a₁ = tan(πλ/2)` exactly), not an unresolved contact term;
  * the ratio 1.00001 was fit bias, and the exact value is 1.
