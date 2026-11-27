# Free-field and Chern–Simons-matter checks of the thermal shear/Hall sum rules (CFT3)

Written 27 Nov 2026 by a research agent. **Not refereed.** Units `T = β = 1`. Marks: **[C]** computed here (script
and output in this folder), **[Q]** quoted from a paper whose arXiv LaTeX source was read, **[I]** inferred.

Definitions are the ones in the task: `T1 = (Txx − Tyy)/2`, `T2 = Txy`, `G^R = [[g, −h], [h, g]]`, `ρ± = Im(g ± ih)`,
`δ` = thermal − vacuum, `S_E = (1/πε) ∫₀^∞ (δρ₊ + δρ₋)/ω dω` (claimed `= 3/8 + t₄/32`), and `S_H` = the
large-frequency (UV) constant of the Hall correlator, normalized by `ε` (as in `../scope_sumrule/hall_uv_sign.py`),
claimed `|S_H| = |α_T|/32`. Convention for `t₄`: `DICTIONARY_T4.md` §1, which is also CDP's. Free real scalar
`t₄ = +4`, free Majorana or Dirac `t₄ = −4`, collider disc `t₄² + α_T² ≤ 16`.

| file | what | runtime |
|---|---|---|
| `euclid_tau.py` → `out_euclid_tau.txt` | Task 1: `S_E` for a free scalar and a free Dirac fermion, in imaginary time | 5 s |
| `light_ops.py` → `out_light_ops.txt` | Task 1 caveats: φ², ψ̄ψ, improvement term at `p ≠ 0`, large-`ω_n` behaviour | ~15 min |
| `cs_fermion_probe.py` → `out_cs_fermion_probe.txt` | Task 3: large-N CS-fermion UV constants vs `t₄/32` and `α_T/32` | ~20 min |

## 1. Free fields: exact match [C]

**Method.** This is independent of `../scope_sumrule/free_fields_sumrule.py`, which uses Hamiltonian matrix
elements. Here everything comes from the Euclidean path integral in the mixed `(τ, k)` representation. The scalar
and the Dirac fermion (γ = Pauli, `S(τ,k) = −γ₀∂_τΔ_F − iγ·kΔ_F`, checked to satisfy `γ·∂S = δ`) are Wick
contracted, with the gamma traces done in sympy. The output is `G_AB(τ) = ∫d²x⟨T_A(τ,x)T_B(0)⟩_β`. At `p = 0`
every term has the form `a(k)e^{−2kτ} + b(k)e^{2kτ} + c(k)`. KMS holds (`b = a e^{−2k}`). The τ-independent part
`c(k)` is the zero-frequency (Drude) weight `ρ/ω = πβC δ(ω)`. The overall sign is fixed by reflection positivity.
Normalizations are checked against the literature, not assumed:

* the vacuum `G₀(τ)` reproduces the Osborn–Petkou `C_T`: scalar `3/(32π²)`, Dirac `3/(16π²) = 2 × 3/(32π²)`;
* `ε` comes from the free energy: scalar `ζ(3)/π`, Dirac `3ζ(3)/(2π)`.

**Spectral functions (exact).** Per real scalar and per 2-component Dirac fermion:

| | `ρ₀(ω)` | `ρ_T/ρ₀` (pair continuum) | Drude `βC/ε` | pair `(2/π)∫δρ/ω / ε` | **`S_E`** | `3/8 + t₄/32` |
|---|---|---|---|---|---|---|
| scalar | `ω³/512` | `coth(ω/4)` | `3/8` | `+1/8` | **`1/2`** | `1/2` (t₄ = +4) |
| Dirac | `ω³/256` | `tanh(ω/4)` | `3/8` | `−1/8` | **`1/4`** | `1/4` (t₄ = −4) |

The ratio is `1.000…` to 20 digits. Closed forms:
* pair term: scalar `ζ(3)/(8π) = ε/8`, Dirac `−3ζ(3)/(16π) = −ε/8`;
* Drude term: scalar `3ζ(3)/(8π)`, Dirac `9ζ(3)/(16π)`, both `= 3ε/8` (integrate `∫k³(−n′)` by parts).

The Euclidean zero mode `∫₀^β G dτ − ∫_ℝ G₀ dτ` reproduces `S_E` exactly.
* `T1` and `T2` give identical correlators (isotropy).
* The cross correlator `⟨T1 T2⟩` vanishes identically, so `h = 0` and `S_H = 0` for both theories, as it must.
* **Conclusion: `S_E = 3/8 + t₄/32` holds exactly at both ends of the allowed `t₄` range, with no discrepancy
  factor.**

**The Drude weight is essential.** In the free theory, at strictly `p = 0`, the intraband term is a `δ(ω)` in
`ρ/ω`. Without it the integral is `+1/8` (scalar) and `−1/8` (Dirac), which is wrong. In imaginary time it is simply
the τ-independent part of `G(τ)`. It is the free limit of the hydrodynamic peak `η/(1 + ω²τ_R²)`, and it is
included automatically in the Matsubara zero mode, i.e. the static susceptibility.

**Improvement term (scalar) [C].** At `p = 0` the `ξ(∂_i∂_j − δ_ij∂²)φ²` term drops out of `T1` and `T2`, since it
contains only spatial derivatives. At `p ≠ 0` (`light_ops.py` §4, `p` along the diagonal, `τ = 1/2`) the result is:
* the canonical piece converges to the `p = 0` value 0.286970;
* the cross term scales as `p²`;
* the `ξ²` piece is `p⁴⟨φ²φ²⟩`, which is IR-divergent at fixed `p` (`∝ p² log(p/m_IR)`) but vanishes as `p → 0`.

So the improvement term matters only through the order of limits. It does not change `S_E`.

**The light scalar φ² (Δ = 1) [C].** In 2+1 D, `⟨φ²⟩_β` of the massless scalar is IR divergent:
`(1/2π) log(1/m_IR)`, checked numerically. The `p = 0` correlators of `T1` and `T2` contain only IR-finite
`k³`-weighted integrals, so they cannot depend on it. In the TT OPE, φ² enters only through descendants
`∂_i∂_j φ²`, which vanish at zero spatial momentum. In any case its `ω^{3−Δ} = ω²` would be analytic, i.e. a contact
term. Directly: `δρ` decays like `e^{−ω/2}`, and `ω_n² δG_E(ω_n) → (2/π)∫ωδρ` (scalar 1.98039, Dirac −3.71322).
The leading thermal OPE datum therefore has Δ = 5. There is no `ω_n²` and no `ω_n⁰` term in the Wick-contracted
correlator. **The extra φ² term in the sum rule is exactly zero.**

**Fermion ψ̄ψ (Δ = 2, parity odd) [C].** In the free theory `⟨ψ̄ψ⟩_β = 0` and `⟨T_A(τ) ψ̄ψ(0)⟩ = 0` at `p = 0`
(trace of a single γ). It does not contribute.

**Contact terms [C]/[I].**
* In the Wick-contracted correlators `δG_E(ω_n → ∞) = 0`, so `δG_E(0) = εS_E`.
* For the scalar this equals `P = ε/2`, the hydrodynamic value used by CDP.
* For the fermion it is `ε/4`. The covariant correlator must then carry an extra temperature-dependent contact
  `ε/4`, a seagull from the vielbein dependence of `T`. This shifts both ends equally, so it cannot affect the
  spectral integral. It equals CDP's OPE constant `ε(1/8 − t₄/32)` at `t₄ = −4`.

## 2. Universal region [Q] + [I]

* **Collider disc** (CDP, arXiv:1707.03007, eqs. `disc` and `causalityconstraints`) [Q]. The energy matrix for
  stress-tensor states is `[[E/4π(1 − t₄/4), α_T E/16π], [α_T E/16π, E/4π(1 + t₄/4)]]`. Positivity gives
  **`t₄² + α_T² ≤ 16`**. CDP note the constraint was "also obtained in Appendix C of" Maldacena–Zhiboedov
  arXiv:1112.4078, without the map to `p_T`. Their `t₄` is the paper's: `t₄ = −4(n_f − n_s)/(n_f + n_s)`.
* **Induced region.** With `S_E − 3/8 = t₄/32` and `|S_H| = |α_T|/32`:
  **`(S_E − 3/8)² + S_H² ≤ (1/8)²`**, a disc of radius 1/8 centred at `(3/8, 0)`. The sign convention for `S_H`
  is irrelevant because the region is a disc.
  * Free scalar `(1/2, 0)` and free fermion `(1/4, 0)` lie on its boundary.
  * Einstein gravity is the centre.
  * Mixed free theories fill the diameter.
  * Large-N CS + fundamental fermion/boson (`t₄ = −4cos πλ`, `α_T = 4 sin πλ`; 1707.03007) trace the whole
    circle, at `(3/8 − cos πλ/8, ± sin πλ/8)`. But see §3: these theories violate the light-scalar condition.
* **Novelty.**
  * The parity-even projection **is already published**: CDP arXiv:1612.00609 §5, eq. `3dcaus`, give
    `P/2 ≤ δG_R(0) ≤ P` in d = 3, i.e. `1/4 ≤ S_E ≤ 1/2`, saturated by free fermions and bosons.
  * The same paper evaluates `S_E = (3 − cos θ)/8` for CS-fermions from ⟨TTT⟩ data, and states that "the parity
    odd term does not play any role in the shear channel".
  * 1707.03007 has the `(t₄, α_T)` disc but no Hall/thermal statement; its introduction only cites the sum-rule
    paper.
  * Bradlyn–Goldstein–Read arXiv:1207.7021 discuss a viscosity sum rule in which the Hall viscosity cancels. That
    is a different (non-relativistic) object.
  * **I found no paper stating a thermal Hall-channel sum rule, or the `(S_E, S_H)` disc.** The new content is
    the parity-odd half: the Hall UV constant `= α_T ε/32` together with the disc. Searches were limited (arXiv
    sources of CDP ×2, MZ cited, two web searches), so this is "not found", not "proven absent".

## 3. Chern–Simons–fermion: data exist, and a UV-side check was done [C]

**Literature [Q].** Ghosh–Mazumdar, arXiv:1912.06589 (v2, 2022), eq. `strf1a`, give the renormalized thermal
correlator `⟨J_{--}(−q)J_{++}(q)⟩_β` of the regular CS-fermion theory, exact in λ at N = ∞. It is computed at
`q^± = 0`, i.e. **zero spatial momentum**, for all Matsubara `q₃`. That is exactly the needed object:
`J_{--} = T1 + iT2`, so `Re G = 2g_E` and `Im G = 2h_E`. The analogous current correlator (σ_xx, σ_xy) is in
Gur-Ari–Hartnoll–Mahajan, arXiv:1605.01122, including an N = ∞ Drude delta tied to the higher-spin charges. The
energy density is from Aharony–Giombi–Gur-Ari–Maldacena–Yacoby, arXiv:1211.4843, eq. `F-nmass-nc-fer`. My code
reproduces their `O(λ²)` coefficient.

**What `cs_fermion_probe.py` does.** It evaluates `strf1a` with the exact `μ_F(λ)`, fits the large-`q₃` expansion
(relative residual ~1e-16, stable under changes of the fit window), and compares UV constants. The normalization
check is `a₃ → (2/256) sin πλ/(πλ)`, the free Dirac value times `c_T(λ)`.

| λ | `S_H = lim h_E/ε` | `sin πλ/8 = |α_T|/32` | ratio |
|---|---|---|---|
| 0.01 | 0.0039263869 | 0.0039263449 | 1.0000107 |
| 0.3 | 0.1011281 | 0.1011271 | 1.0000100 |
| 0.5 | 0.1250011 | 0.125 | 1.0000091 |
| 0.9 | 0.0386274 | 0.0386271 | 1.0000074 |

* **Hall UV constant:** `|S_H| = |α_T|/32` holds to 1e-5 (limited by the fit) for all λ tested.
* **Shear UV constant:** `a₀/(2ε) = (cos πλ − 1)/8 = −1/8 − t₄/32` holds to ~1e-6. So the λ-dependence (the
  `t₄/32` term) is confirmed. The λ-independent offset `−1/4` relative to CDP's `1/8 − t₄/32` is the same contact
  convention as in the free Dirac case.

**Caveat: CS-matter does not satisfy the sum-rule hypotheses [C].** For λ ≠ 0 there are two extra thermal terms
linear in `q₃`:
* **Shear channel:** a non-analytic `a₁|q₃|` term, with `a₁ → μ_F²/16 ∝ λ²` and `a₁/ε` up to 0.09 at λ = 0.9.
  This is the parity-odd scalar ψ̄ψ (Δ = 2 + O(1/N), `⟨ψ̄ψ⟩_β ∝ μ_F ≠ 0`) entering TT: the `⟨O⟩ω^{3−Δ}` term. It
  gives `δρ ∝ ω` at large ω, so **the spectral integral `S_E` diverges** for CS-fermion at λ ≠ 0, unless that term
  is subtracted.
* **Hall channel:** an analytic `b₁ q₃` term (`∝ λ³`), a temperature-dependent contact-type term of the kind
  Ghosh–Mazumdar flag as unresolved. It is linear in ω in `G^R`, so it shifts what one would read off as `η_H` and
  must be handled before an integral form of the Hall sum rule can be tested.

So CS-matter confirms the **OPE/UV-constant side** of both claims exactly. It cannot test the **spectral-integral
form** without:
1. subtracting the ψ̄ψ terms;
2. the N = ∞ Drude weight for T (a Matsubara `n = 0` or time-averaged computation, as in 1605.01122);
3. fixing the temperature-dependent contact terms of their cutoff scheme.

That is a multi-day project.
