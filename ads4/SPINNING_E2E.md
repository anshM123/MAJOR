# Spinning Regge-moment formula in d = 3: end-to-end test on exact `⟨JφJφ⟩` and `⟨TφTφ⟩` blocks, with the normalisation derived

Status: written 16 Nov 2026 by a derivation/numerics agent; **not refereed**. Code and recorded outputs: `ads4/code/e2e/`
(index in §12; `sh run_all.sh` regenerates every number; `sh run_all.sh fast` reuses the cached block tensors in `data_modes/`).
Sources: LaTeX of arXiv:2106.10274 ("CHMRS"), 2311.04271 ("CLS-D"), 2102.08160 ("CHL"), 2201.06602 ("CLPS") in the session
scratchpad `…/scratchpad/reggemoment/src/`; equation numbers are those of the arXiv versions (for CLS-D the LaTeX labels, since
numbers were not compiled). Costa–Penedones–Poland–Rychkov (arXiv:1107.3554, 1109.6321) were used from memory only for the *form*
of the differential operators; every property used is re-derived or checked here.

Marks: **[DVE]** derived and verified end-to-end; **[DVP]** derived, verified piecewise; **[D]** derived, not verified;
**[I]** inferred; **[ND]** not done; **[Q]** quoted. "Numerical evidence" = a statement read off from the numbers, with no derivation.

## 0. Summary

1. **Correlator.** `⟨O_1 φ_2 O_3 φ_4⟩`, `O = J` (`l = 1`, `Δ = 2`) or `O = T` (`l = 2`, `Δ = 3`), `Δ_φ = 1.75`, `s`-channel `(12) → (34)`: the conserved
   operators form the **Regge pair** `(1,3)`, which is where "helicity preserved / flipped at the spinning pair" is defined; exchange
   of a parity-even `(Δ, J)` with the unique conserved parity-even `⟨OφO_{Δ,J}⟩` structure. This is the one-pair half of `⟨TTTT⟩`. The
   channel `TT → φφ` was not done (§11).
2. **Exact blocks [validated].** Scalar blocks with general external dimensions by a radial Casimir recursion (C; double and
   `__float128`); spinning partial waves from the differential operators `D_11, D_12` on 4 (currents) or 9 (stress tensors) scalar
   seeds with analytic derivatives to fourth order. Validation (§2): conservation `10⁻¹⁴`, OPE limit (linear approach, 4 cases),
   mean-field decomposition of the seeds `5·10⁻¹⁶`, exact consistency of the shadow coefficients with conservation, and a direct
   quadrature of the Lorentzian shadow representation (`0.1–0.3 %`, quadrature-limited).
3. **Functional.** The CHMRS Regge moment applied to the Fourier modes `(q_1,q_3)` of the stripped tensor (`Π̂^{(q_1,q_3)}_{k,ν}`, §3). These
   span all leading-Regge `ν`-harmonic tensor kernels; the CLS-D spinning functional (`spinningfunctional_definition`) reduces in `d = 3`
   to an explicit combination of them (§3.3, §9).
4. **Derived leading formula with absolute normalisation [DVE]** (§4.5, no constant fitted anywhere):
   `(−1)^{J+l} Π̂^{(q_1,q_3)}_{k,ν}[s_13^{l+1}s_24^{Δ_φ}W_diff] → (π/4)·4^{l+1+Δ_φ}/b^{(1−Δ_φ)}_{Δ,J} · γ_{a_x}(ν)γ_{a_y}(ν) m^{−a_x−a_y} · |T_l| · [M_pres P_{−½+iν}(cosh β) + M_flip F^ν_{2l}(β)]`,
   and `⟨OφOφ⟩ ⊃ λ²_{OφO}(−½)^J W_diff/L_0²` (§1.4). Measured/predicted at `m ≈ 100`, `k = 6` (§5): `l = 1`: `R_pres = 0.991, 0.988` and
   `R_flip = 0.985, 0.982` (`β = 0.4, 1`, `ν = 0`); `l = 2`: `R_pres = 0.985, 0.987`, `R_flip = 0.968, 0.975`; with `m²(R−1)` converging, i.e.
   `R = 1 + c_2/m² + …` for both helicity components, `l = 1, 2`, `k = 4, 6`, `ν = 0, 1, 3`, `β = 0.4, 1, 2`, `m = 10…100`.
5. **`O(1/m)` [D + numerical evidence]**: absent for both `l`, both components, in the Weyl-odd variables `(m, β)` of CHMRS with the exact
   `1/b` — proved by a reflection argument (§6.1), verified pointwise, and by fits: `|c_1| ≤ 0.008` in 11 of 12 preserved cases (`k = 6`, `β ≤ 1`; all `≤ 0.027`), against `|c_2|/m ≈ 1`.
6. **`O(1/m²)` [numerical evidence]**: measured coefficients (`−25 … −740`, §6.2). Spin-dependent part of the preserved one: the
   exact normalisation polynomials (derived) plus the flat-space helicity term in the curved form `−2l²P′/P`: fitted coefficient
   `−2.01` (`l = 1`, predicted `−2`) and `−8.00` (`l = 2`, predicted `−8`) (§6.3); flip corrections measured but not understood (§6.4).
7. **Uniformity [numerical evidence + D]** (§7): fails at fixed `J` by the factor `Π[(J+½)²−(i−½)²]/(J+½)^{2l}` (derived; after restoring it the
   formula is uniform down to `J = 1`, `β ≈ 0.04`); not uniform in `ν/m` (`O(ν²/m²)`), but the spin dependence is: spinning = scalar
   skeleton × spin factor to `10⁻³` for `ν/m ≤ 0.2`.
8. **CLS-D** (§9): no discrepancy; their polarisation prescription is helicity preserving only up to `O(ν^{−2l})`.
9. **`⟨TTTT⟩` not reached** (§8); **MFT normalisation (`2P/m^{2k}` form) not done** (§0.1).

### 0.1 Answers to the coordinator's questions (16 Nov, after the restart)

**(1) Is `C^±_{k,ν} = (2P/m^{2k})(|c^{++}|² ± |c^{+−}|²)(1+…)` confirmed end-to-end with no fitted constants?** *Partly.* What is confirmed
end-to-end, with no fitted constant, is the absolute heavy-block value of the leading-Regge tensor moments of the exact
`⟨OφOφ⟩` partial wave, for both helicity components at the spinning pair, and through §1.4 its value per unit `λ²_{OφO}`. The
normalisation was fixed as follows, every step derived and separately checked:
  * blocks: `g ∼ z^{(Δ−J)/2}z̄^{(Δ+J)/2}` (CHMRS (2.33)), fixed by the scalar **mean-field decomposition** (`2P_FK` exactly, §2.1);
  * spinning structures: the **conserved** parity-even `⟨OφO′⟩` structure (conservation imposed symbolically; unique; checked on
    the four-point function to `10⁻¹⁴`), normalised by its `V_1^lV_5^J` coefficient; the map to the differential basis (`L_0`) is exact;
  * the relation "partial wave ↔ `λ²`" by the **OPE limit** (§2.3, `λ²(−½)^J W_diff/L_0²`);
  * the dispersive side by the **Lorentzian shadow coefficient** `b^{(1−Δ_φ)}`, fixed by conservation (exact identity, §2.3) and by direct
    quadrature;
  * `T` itself: the **Ward identity** fixes `λ_{φφT}` (Osborn–Petkou, quoted) and `⟨TT⟩ = C_T`; it enters only as the conversion
    `λ²_{TφO} → λ²_{TφO}/C_T` for a unit-normalised `T` and is **not tested numerically** (nothing in a block-by-block test depends on it).
  What is **not** confirmed: (a) the MFT-normalised form "`2P/m^{2k}`": it requires the large-`Δ` mean-field density of `[Oφ]_{n,J}` in the
  conserved basis (both parities), which was not computed; the formula above is the complete replacement of it, expressed through
  `b^{(1−Δ_φ)}` instead of through `q^{MFT}`. For scalars the two are equivalent (CHMRS (2.59), verified in `REGGE_MOMENT_TTTT.md` §1.3);
  (b) the combinations `|c^{++}|² ± |c^{+−}|²` themselves: they need both pairs spinning; only the one-pair factors are verified (`R`
  end-to-end; `S`, which produces the `−` sign, only by direct `H²` integration, since parity-even `⟨OφOφ⟩` blocks do not probe it).
  So: **one-pair normalisation [DVE]; `2P/m^{2k}` MFT form [ND]; `±` structure for `⟨TTTT⟩` [D, factors verified piecewise].**

**(2) Approach to 1.** `1 − R ∼ c_2/m²` at fixed `ν`, both components, `k = 4, 6`, `l = 1, 2`: `m²(R−1)` converges (e.g. `l = 2`, `β = 0.4`, `ν = 0`,
`k = 6`: `−96.3, −130.2, −146.4, −151.2, −152.5` at `m = 12…102`). The `**` flags at `k = 2` (and at `k = 4`, `m ≳ 70`) are numerical, not
physical: the moment integrand behaves like `u^{k−1/2}` at `u = √(ρ̄/ρ) → 0`, where the radial series needs `N ≳ 70/(4su)` levels; the region
that cannot be reached with `N ≤ 10⁴` is dropped and its contribution is estimated (`tail`); `**` means `tail > 10⁻³`, and for `k = 2` it reaches
`10⁻¹`. The negative or `> 1` entries at `β = 2`, `ν = 3` (and `β = 1`, `ν = 3`, `m ≈ 10`) are the `1/m²` term exceeding the leading one:
`c_2 ≈ −730` (preserved, `k = 6`), so `|c_2|/m² > 1` for `m < 27`; in addition `P_{−½+3i}(cosh 2)` is near a zero (`zeroP ≈ 0.35` at `ν ≈ 2.9`).
The same entries approach 1 at larger `m` (`0.61 → 0.86` at `m = 40 → 70`, `l = 1`), and `m²(R−1)` stabilises (`−645, −707`, `l = 2`).

**(3)** `O(1/m)`: cancels (derived + verified, §6.1). `O(1/m²)`: §6.2–6.4. On "opposite-sign helicity terms": in `⟨OφOφ⟩` the term
`∝ (|c^{++}|² − |c^{+−}|²)` of `REGGE_MOMENT_TTTT.md` (the `g_1²S_{13}S_{24}` / `B_3` term) cannot appear — it needs a parity-odd structure at
each pair; what can be measured is the helicity-*flip* component (`F^ν_{2l}`, a different harmonic function) and the helicity
dependence of the `1/m²` term of the preserved one, whose flat-space part `+2λ²(ν/β)J_1(νβ)` (sign and coefficient `2l²`) is confirmed
in the curved form `−2l²P′/P` (§6.3).

**(4)–(6)**: §7, §9, §10–11.


## 1. Conventions

### 1.1 Kinematics, embedding space, real section

* `d = 3`, signature `(−,+,+)`. Embedding: `P = (P⁺,P⁻,P^μ) = (1, x², x^μ)`, `Z = (0, 2x·z, z^μ)`, `P_i·P_j = −x_ij²/2`,
  `Z_i·P_j = z_i·(x_j − x_i)`, `Z_i·Z_j = z_i·z_j`. In the code `s_ij := 2P_i·P_j = −x_ij²`.
* Scalar blocks: `a = −Δ_12/2`, `b = Δ_34/2`; `⟨φ_1φ_2φ_3φ_4⟩ ⊃ K g^{a,b}_{Δ,J}(u,v)`,
  `K = (x_24²/x_14²)^{Δ_12/2}(x_14²/x_13²)^{Δ_34/2}/(x_12^{Δ_1+Δ_2}x_34^{Δ_3+Δ_4})`, `u = x_12²x_34²/(x_13²x_24²)`, `v = x_14²x_23²/(x_13²x_24²)`,
  `g ∼ z^{(Δ−J)/2} z̄^{(Δ+J)/2}` for `z ≪ z̄ ≪ 1` (the normalisation of CHMRS (2.33); `g = (4r)^Δ Σ B_{n,j} r^n P_j(η)`, `B_{0,J} = J!/(½)_J`).
* CHMRS variables (2.28)–(2.30): `α = Δ−J−2`, `ᾱ = Δ+J−1`, `m² = αᾱ`, `cosh β = η_AdS = (α²+ᾱ²)/(2αᾱ)`. Weyl-odd combinations:
  `A_+ := Δ − 3/2 = m cosh(β/2)`, `A_− := J + ½ = m sinh(β/2)`.
* **Frame** (CHMRS (C.2)): `(x_1,x_2,x_3,x_4) = (−y, x/x², y, −x/x²)`, `x, y` future timelike; `r = |x||y|`, `η = −x·y/(|x||y|) = cosh ξ` are the
  `u`-channel radial variables (`ρ, ρ̄ = r e^{±ξ}`, `w = −(1−ρ)²/(4ρ) < 0`). All `x_ij` are timelike, `u, v > 0`, `z, z̄ < 0`.
* **Real section.** A partial wave on this section is defined with every non-integer power `(x_ij²)^p → |x_ij²|^p` and `g` the
  real-analytic block at `z, z̄ < 0` (CHMRS's `G^s` "with absolute values"). The double discontinuity of the physical correlator is
  `2 sin²(π(τ − Δ_O − Δ_φ − l)/2)` times this object (`τ = Δ − J`; both commutators give the same sine because the two vertices have the
  same exponents; for `l = 0` this is CHMRS's factor). All "Regge moments" below are those of the real-section object, i.e. the
  quantity called `Π̂[G^s]/(2 sin²)` in CHMRS (3.16).
* Functions: `γ_a(ν) = |Γ((a − ½ + iν)/2)|²` (CHMRS (2.55), `d = 3`); `F^ν_μ(ξ) := (1/2π)∫_0^{2π}dθ (cosh ξ − sinh ξ cos θ)^{−½+iν} cos μθ`
  `= Γ(½+iν)/Γ(½+iν+|μ|) · P^{|μ|}_{−½+iν}(cosh ξ)` (type-3 Legendre function), `F^ν_0 = P_{−½+iν}(cosh ξ)`; `F^ν_μ → (−i)^μ J_μ(νξ)` for `ξ → 0`, `νξ` fixed.
  In the notation of `SPINNING_KERNELS.md`, `F^ν_μ = Φ^{−ν}_{0,μ}`.

### 1.2 Three-point structures

For `⟨O_1(P_1,Z_1) φ_2(P_2) O_5(P_5,Z_5)⟩`, `O_1` of spin `l` and dimension `Δ_1`, `O_5 = (Δ, J)`:

`V_1 = [(Z_1·P_2)(P_1·P_5) − (Z_1·P_5)(P_1·P_2)]/(P_2·P_5)`, `V_5 = [(Z_5·P_1)(P_2·P_5) − (Z_5·P_2)(P_1·P_5)]/(P_1·P_2)`,
`H_15 = (Z_1·Z_5)(P_1·P_5) − (Z_1·P_5)(Z_5·P_1)`,

`S_n := V_1^{l−n} H_15^n V_5^{J−n} / (s_12^{κ_12} s_15^{κ_15} s_25^{κ_25})`, `n = 0…l`, `κ_12 = (Δ_1+l+Δ_φ−Δ−J)/2`, `κ_15 = (Δ_1+l+Δ+J−Δ_φ)/2`, `κ_25 = (Δ_φ+Δ+J−Δ_1−l)/2`

("standard basis"; the same combinations as the `V, H` of Costa et al. up to powers of 2. For `l = 0`: `S_0 = 2^J S3` with `S3` the seed
defined below, and the usual position-space structure `(z·x_13/x_13² − z·x_23/x_23²)^J/(…)` equals `(−2)^J S3 = (−1)^J S_0`, see §1.4).

**Differential basis.** `D_11 = (P_1·P_2)(Z_1·∂_{P_2}) − (Z_1·P_2)(P_1·∂_{P_2})`, `D_12 = (P_1·P_2)(Z_1·∂_{P_1}) − (Z_1·P_2)(P_1·∂_{P_1}) + (Z_1·P_2)(Z_1·∂_{Z_1})`
(the operators of arXiv:1109.6321 for a scalar at point 2; they are tangent to the cone, preserve transversality, commute with each
other — checked symbolically — and lower `Δ_1`, resp. `Δ_2`, by one while raising the spin at 1). With the scalar seed
`S3(Δ_1′,Δ_2′) = v_5^J/(s_12^{α}s_15^{β}s_25^{γ})`, `v_5 = (Z_5·P_1)(P_2·P_5) − (Z_5·P_2)(P_1·P_5)`, `(α,β,γ) = ((Δ_1′+Δ_2′−Δ+J)/2, (Δ+J+Δ_12′)/2, (Δ+J−Δ_12′)/2)`:
`W_w := D_w S3(Δ_1 + #a, Δ_φ + #b)`, `w` a word in `a ≡ D_11`, `b ≡ D_12`. One finds (`structures.py`, exact)

`D_11 S3 = S3·[(Δ+J−1−Δ_12) V_1/2 − J H_15/V_5]`, `D_12 S3 = S3·(s_25/s_15)[(Δ−J−1+Δ_12) V_1/2 + J H_15/V_5]` (dimensions of the result).

**Conserved structures [D; conservation verified on the four-point function, §2.3].** The combination `Σ_w c_w W_w` is conserved at point 1
(physical-space divergence computed symbolically at a rational point, unique null vector) for

* `l = 1` (`Δ_1 = 2`): `c_a = −(Δ−Δ_φ−J)(Δ−Δ_φ+J+1)`, `c_b = Δ² − 3Δ − Δ_φ² + 3Δ_φ + J² + J`;
* `l = 2` (`Δ_1 = 3`): `c_aa = (Δ−Δ_φ−J)(Δ−Δ_φ−J+2)(Δ−Δ_φ+J+1)(Δ−Δ_φ+J+3)`,
  `c_ab = −2(Δ−Δ_φ−J)(Δ−Δ_φ+J+1)(Δ² − 3Δ − Δ_φ² + 5Δ_φ + J² + J − 6)`,
  `c_bb = Δ⁴ − 6Δ³ − 2Δ²Δ_φ² + 6Δ²Δ_φ + 7Δ² + 6ΔΔ_φ² − 18ΔΔ_φ + 6Δ + Δ_φ⁴ − 6Δ_φ³ − 2Δ_φ²J² − 2Δ_φ²J + 11Δ_φ² + 6Δ_φJ² + 6Δ_φJ − 6Δ_φ + J⁴ + 2J³ − 5J² − 6J`.

In the standard basis `Σ_w c_w W_w = 2^{−J} L_0 · 𝔖^{cons}`, `𝔖^{cons} := S_0 + ℓ_1 S_1 (+ ℓ_2 S_2)`:

* `l = 1`: `L_0 = −J(Δ−1)(Δ−Δ_φ−J−1)`, `ℓ_1 = −2(Δ−Δ_φ)/(Δ−Δ_φ−J−1)`;
* `l = 2`: `L_0 = ½ΔJ(Δ−1)(J−1)(Δ−Δ_φ−J−3)(Δ−Δ_φ−J−1)`, `ℓ_1 = −4(2Δ−2Δ_φ+J)/(Δ−Δ_φ−J−3)`,
  `ℓ_2 = 4(2(Δ−Δ_φ)² − J² − J)/((Δ−Δ_φ−J−3)(Δ−Δ_φ−J−1))`.

(The normalisation "coefficient of `S_0` equal to 1" is singular at `Δ − J = Δ_φ + 1` and `Δ_φ + 3`, i.e. exactly at double-twist values;
it is used only to state the result in a familiar basis. `c_{b^l}` is invariant under the Weyl reflection `(Δ,J) → (3−Δ, −1−J)`.)

### 1.3 Partial waves

`W_diff(P_i;Z_1,Z_3) := Σ_{w,w′} c_w c_{w′} D^{(12)}_w D^{(34)}_{w′} [K^{(w,w′)} g^{(w,w′)}_{Δ,J}]`, with `D^{(34)}` the same operators for `(1,2) → (3,4)` and the seed
dimensions `(Δ_1+#a_w, Δ_φ+#b_w, Δ_1+#a_{w′}, Δ_φ+#b_{w′})`. Stripped frame tensor: `𝒢(z_1,z_3; r,η) := s_13^{l+1} s_24^{Δ_φ} W_diff` (Cartesian
components; invariant under the residual dilatation and Lorentz-covariant).

### 1.4 Normalisation in terms of OPE coefficients [DVP]

* Scalars (derived in §2.3 from the leading OPE term and checked numerically): if `⟨φ_1φ_2O(x_3,z)⟩ = λ_{12O}(z·x_13/x_13² − z·x_23/x_23²)^J/(x_12^{Δ_1+Δ_2−Δ+J}x_13^{Δ+Δ_12−J}x_23^{Δ−Δ_12−J})`
  and `⟨O(x,z)O(0,z′)⟩ = (z·I(x)·z′)^J/x^{2Δ}`, then `⟨φ_1φ_2φ_3φ_4⟩ ⊃ (−½)^J λ_{12O}λ_{34O} K g`. In embedding terms the usual
  structure is `(−2)^J S3`, hence the partial wave of `(S3, S3)` is `(−1/8)^J K g`.
* Spinning: if `⟨O_1(x_1,z_1)φ(x_2)O(x_5,z_5)⟩ = λ_{OφO} · 𝔖^{cons}` (unit-normalised exchanged operator), then
  **`⟨O_1φ_2O_3φ_4⟩ ⊃ λ_{OφO}² (−½)^J W_diff / L_0²`** (Euclidean configurations; on the real section the same with the phases
  stripped). Checked by the OPE limit in §2.3: `W_diff = 8^J ⟨C_L, R⟩` on the real section, where the factor `(−1)^J` relative to the
  Euclidean `(−8)^J` is the usual sign between `P_J(η)` and `P_J(−η)`.
* Normalisation of `T`: for a canonically normalised stress tensor (Ward identity
  `⟨φ(x_1)φ(x_2)T(x_3,z)⟩ = −(dΔ_φ/((d−1)S_d)) (z·x_13/x_13² − z·x_23/x_23²)²/(x_12^{2Δ_φ−d+2}x_13^{d−2}x_23^{d−2})`, `S_d = 2π^{d/2}/Γ(d/2)` [Q, Osborn–Petkou])
  one has `⟨T(x,z)T(0,z′)⟩ = C_T (z·I·z′)²/x^{2d}`, `C_T ∝ c_T`. `λ_{TφO}` above refers to this `T`; in terms of the unit-normalised
  `T̂ = T/√C_T` the coefficient entering a normalised correlator is `λ_{TφO}²/C_T`. Nothing in the tests depends on `C_T`.
* Sign: `L_0(Δ,J) L_0(3−Δ,−1−J)` has the sign `(−1)^l`, and the positive quantity on the real section is `(−1)^{J+l} Π̂[𝒢]`
  (for `l = 0` this is CHMRS's `(−1)^J`; the extra `(−1)^l` is the sign of `⟨Oφ O′⟩ → ⟨φ O O′⟩`).


## 2. Task 1 — exact blocks: implementation and validation

### 2.1 Scalar blocks with general external dimensions (`casimir_terms.py`, `cblock.c`, `blocks.py`)

The Dolan–Osborn Casimir operator with `a, b ≠ 0`, `D_z = z²(1−z)∂_z² − (1+a+b)z²∂_z − abz`, is transformed symbolically to radial
variables and multiplied by `(1−r²)(1+2rη+r²)²(1−2rη+r²)`; acting on `r^E P_j(η)` it becomes 22 terms `r^k η^i × {1, (1−η²)∂_η}` with
coefficients linear in `E², E, j(j+1), C` (`terms.json`), i.e. a recursion of depth 8 in `n` and band `j → j ± 3`. The C code performs the
recursion and evaluates `Σ B_{n,j} (±1)^n (Δ+n)^k r^n P_j^{(l)}(cosh χ)` for `k, l ≤ 4` (all derivatives analytic, positive Legendre
recurrences), in `double` and `__float128`.

Checks (`out_t1_scalar_blocks.txt`, `out_t1b_mft_scalar.txt`):

| check | result |
|---|---|
| `a = b = 0` against the independent recursion of `reggemoment/blocks3d.py` | `2·10⁻¹⁶ … 10⁻¹²` (growing with `|ρ|`, round-off of the old code) |
| general `(a,b)`: double vs quad; quad vs independent mpmath implementation | `≤ 10⁻¹⁵`; `≤ 1.2·10⁻¹⁵` |
| Casimir equation in `(z,z̄)` (numerical differentiation of the quad block), four cases incl. `z,z̄ < 0` | residual `2.6·10⁻⁹ … 1.2·10⁻¹¹` |
| leading term `g/(z^{(Δ−J)/2}k_{Δ+J}(z̄)) − 1`, `k_β = x^{β/2}₂F₁(β/2+a, β/2+b; β; x)` | `∝ z`: `(2.0–2.9)·10⁻⁴, ·10⁻⁶, ·10⁻⁸` at `z = 10⁻⁴, 10⁻⁶, 10⁻⁸` |
| derivative sums vs 4th-order Taylor remainder | remainder `∝ δ⁵` (ratios 31.3, 31.7) |
| **mean-field decomposition**, channel `(12)(34)` of `⟨φ_1φ_2φ_2φ_1⟩`: `Σ_{n,J} P_{n,J} g^{a,a}_{Δ_1+Δ_2+2n+J,J} = u^{(Δ_1+Δ_2)/2}v^{−Δ_2}`, four pairs `(Δ_1,Δ_2)`, two points | `5·10⁻¹⁶ … 7·10⁻¹⁶` |

The mean-field coefficients are those of Fitzpatrick–Kaplan (arXiv:1111.6972) in the form
`P_{n,J} = (Δ_1−h+1)_n(Δ_2−h+1)_n(Δ_1)_{J+n}(Δ_2)_{J+n} / [n! J! (J+h)_n (Δ_1+Δ_2+n−2h+1)_n (Δ_1+Δ_2+2n+J−1)_J (Δ_1+Δ_2+n+J−h)_n]`, `h = 3/2`;
at `Δ_1 = Δ_2` the sum over even `J` of `2P` is CHMRS (2.33) exactly (ratio 2.000000000000 for `n ≤ 2`, `J ≤ 6`), which fixes the block
normalisation; the test at `Δ_1 ≠ Δ_2` fixes the conventions for `(a, b)`.

Two numerical points that matter for dispersive functionals (both were found the hard way):

* On the real section the block is needed at `ρ, ρ̄ < 0`. We use `g^{a,b}(u,v) = (−1)^J v^{−b} g^{−a,b}(u/v,1/v)` so that the series is
  evaluated at positive radial coordinates; its coefficients are **positive iff the two parameters of the reflected block have
  non-negative sum**. For negative sum (the seed `(b^l,b^l)`) the terms alternate and cancel by many orders of magnitude near the
  Regge limit (e.g. the scalar test of §4.1 with `Δ_1 < Δ_2` failed at `m ≈ 80`); there we use in addition the Euler transformation
  `g^{c_1,c_2} = v^{−c_1−c_2} g^{−c_1,−c_2}`.
* Truncation: with `X := −N log|ρ|_max` the error of a derivative sum of order `p ≤ 4` truncated at level `N` behaves like `e^{−X}X^{p+1}`.
  Unconverged derivative sums are *not small* (the conserved combination no longer cancels), so points with `X − 10 log X < 25` are
  dropped from the integration grid. The dropped region is `u ≡ √(ρ̄/ρ) ≲ 30α/(Nt)` (§5.1).

### 2.2 Spinning partial waves (`embed.py`, `structures.py`, `spinblock.py`)

`W_diff` (§1.3) is obtained by applying the words `D_w ⊗ D_{w′}` symbolically to `K̃ g̃(P,Q)`, where the seed is written as
`K g^{a,b} = (−1)^J K̃ g̃`, `K̃ = K v^{−b}`, `g̃ = g^{−a,b}(u/v,1/v)`, and `P = ½ log(uv)`, `Q = ½ log(v/u)`. The result is a polynomial in the seven
invariants `Z_1·P_j, Z_3·P_j, Z_1·Z_3`, in `1/s_ij`, and in the jets `∂_P^i∂_Q^j g̃`, `i+j ≤ 2l` (5317 terms for the largest word pair of `l = 2`;
generated code `gen_phi_l*.py`). The jets are obtained from the radial derivative sums by exact Taylor composition with
`cosh L = e^{(Q−P)/2} + e^{Q}`, `|η| = e^{Q} − e^{(Q−P)/2}`. The variables `(P,Q)` are essential: in the Regge region `∂_P = O(1)` while `∂_Q = O(1/r)`,
and with `(log u, log v)` jets the conserved combination is a difference of terms `(1/r)^{2l}` times larger.

### 2.3 Validation of the spinning partial waves (`out_t2_spinning_blocks.txt`, `out_t2b…`, `out_t5…`)

| check | cases | result |
|---|---|---|
| conservation at `x_1` of the four-point `W_diff`: `|∂^μW_{μ…}| / max|∂_μW_{ν…}|` (quad precision, 4th-order differences) | `l = 1`: `(Δ,J,Δ_φ) = (9.75, 3, 1.75)`, `(6.5, 4, 2.5)`, two configurations each | `2.5·10⁻¹⁴, 1.3·10⁻¹⁴, 1.4·10⁻¹⁷, 3.4·10⁻¹⁵` |
| | `l = 2`: `(9.75, 3, 1.75)`, `(8.5, 4, 2.5)` | `2.8·10⁻¹⁴, 4.9·10⁻¹⁴, 2.3·10⁻¹⁴, 3.8·10⁻¹⁴` |
| control: same with `c_a ↔ c_b` | `l = 1` | `0.70`, `1.9` |
| OPE limit `x_12 = δ·v → 0`: `W_diff/(8^J⟨C_L,R⟩) − 1` at `δ = 10⁻³, 10⁻⁴, 10⁻⁵` (must vanish linearly) | `l = 1`, `(9.75, 3, 1.75)` | `0.1407, 0.01398, 0.001397` |
| | `l = 1`, `(6.5, 4, 2.5)` | `4.08·10⁻³, 4.07·10⁻⁴, 4.07·10⁻⁵` |
| | `l = 2`, `(9.75, 3, 1.75)` | `7.23·10⁻³, 7.20·10⁻⁴, 7.20·10⁻⁵` |
| | `l = 2`, `(8.5, 4, 2.5)` | `7.40·10⁻³, 7.38·10⁻⁴, 7.38·10⁻⁵` |
| double vs quad evaluation of the stripped tensor in the dispersive region (`r = 0.002`, `η = 1.4`, 9 polarisation pairs) | `(Δ,J) = (30.25, 6)` | `1.1·10⁻¹²` (`l = 1`), `1.0·10⁻⁸` (`l = 2`) |
| shadow coefficients vs conservation: `c_w(Δ,J)/(b_w c_w(3−Δ,−1−J)) · b^{(1−Δ_φ)}` for every word `w` | 4 cases per `l`, `Δ` up to 100.5 | `1.0000000` (all digits) for all words |
| control: the same with `b^{(0)}` for all words | | `1.159, 1.037, 0.929, 1.012` (`l = 1`, word `a`); `4.78, 1.44, 1.16, 1.13` (`l = 2`, word `aa`) |
| direct 4-d quadrature of the Lorentzian shadow representation / `W_diff` | `l = 1`: `(7.3, 2)`, `(8.6, 3)`, two polarisation pairs | `0.9987, 1.0009, 0.9993, 1.0008` |
| | `l = 1`, `(7.3, 2)` towards the Regge region (`√r = 0.15, 0.1`; `ξ = 0.5, 1.5`) | `0.9989, 0.9983, 1.0012` |
| | `l = 2`: `(10.3, 2)`, `(11.6, 3)` | `1.0028, 1.0007, 1.0029, 1.0010` |
| | scalar seeds, five sets of external dimensions incl. `Δ_12 ≠ Δ_34`, `J = 2, 3` | `1.0007 … 1.0032` |

(The deviations of the last three lines are the error of the four-dimensional Gauss quadrature with endpoint singularities, not of
the blocks; they are of both signs and do not scale with `Δ` or `l` in a systematic way.)

* *Casimir equation*: automatic, since the operators commute with the Casimir of `(12)`; it is not an independent check and is not
  listed. *Conservation* is: it tests the coefficients `c_w`, the 4-point implementation of the operators and the jets.
* *OPE limit*: `C_L(n) := |x_25²|^Δ L(x_1,x_2,x_5→∞; z_1, I(x_25)n)` is the leading OPE tensor of the left structure, `R(n)` the right
  structure with the exchanged operator at `x_2`; the spin-`J` contraction on the null circle is
  `⟨A,B⟩ = Σ_m a_{−m}b_m 2^J(−1)^{J−m}/binom(2J,J−m)` (`a_m, b_m` Fourier coefficients of `A(n(φ))`, `B(n(φ))`). This is independent of the
  four-point differential operators.
* *Shadow coefficients vs conservation*: the right structure of the shadow representation is `Σ_w c_w(Δ,J) b_w^{−1} D_w S3^S`; since the
  shadow transform preserves conservation and the conserved structure is unique, `c_w(Δ,J)/(b_w c_w(3−Δ,−1−J))` must not depend on `w`.
  With `b^{(Δ_34)}` of §4.1 the three (two) values agree to all 20 digits printed, and equal `1/b^{(1−Δ_φ)}`; with `b^{(0)}` they differ
  by `O(1)`. This is an exact test of the `Δ_34`-dependence of the Lorentzian shadow coefficient.
* *Direct shadow integral*: `W_diff` equals `(−1)^J/(4b^{(1−Δ_φ)}) ∫_{x_3<x_5<x_4} d³x_5 ∫D n L(z_1)F^S(z_3)` (§4.2) evaluated by 4-dimensional
  Gauss quadrature, at moderate cross-ratios and towards the Regge region.
* *Mean-field decomposition of `⟨OφOφ⟩` itself*: **not done** (it needs the parity-odd partial waves as well, or a lightcone
  analysis of the leading-twist tower; note that `L_0 = 0` on that tower). The conventions are fixed instead by the OPE limit and the
  direct shadow integral.


## 3. Task 2 — the functional: tensor Regge moments

### 3.1 Definition

Scalars [Q, CHMRS (2.39)–(2.41), (2.50), (C.1)]: `Π̂_{k,ν}[𝒢] = ∫_1^∞ [dη] P_{−½+iν}(η) ∫_0^{r_max(η)} dr r^{k−2} 𝒢(r,η)`, `[dη] = 2dη` (`d = 3`);
`C_{k,ν} = (a_{Δ_φ}/γ²) Π̂_{k,ν} + O(Π̂_{k+2})` [(2.58), (C.16)].

Tensor case. At the frame point `(r,η)` let `ŷ = y/|y|`, `e_∥` the unit tangent at `ŷ ∈ H²` pointing to `x̂`, `e_⊥` the third unit vector, and
`z(θ) = ŷ + cos θ e_∥ + sin θ e_⊥` (the null polarisations = the celestial circle of the Regge pair). The stripped tensor is a trigonometric
polynomial, `𝒢(z(θ_1), z(θ_3); r, η) = Σ_{|q_1|,|q_3| ≤ l} 𝒢_{q_1q_3}(r,η) e^{i(q_1θ_1+q_3θ_3)}` (`q = ∓l`: the transverse tensors `ε_±^l/2^{l/2}`; `q = 0, ±1`:
components along `ŷ`). Define

**`Π̂^{(q_1,q_3)}_{k,ν}[𝒢] := ∫_1^∞ [dη] F^ν_{q_1+q_3}(ξ) ∫_0^{r_max} dr r^{k−2} 𝒢_{q_1q_3}(r,η)`.**

For scalars this is `Π̂_{k,ν}`. These `(2l+1)²` moments are exactly the leading-Regge, `ν`-harmonic tensor functionals: the kernel is
`K^ν_{q_1q_3}(x̂,ŷ) = (1/2π)∫D ℓ (−x̂·ℓ)^{−½+iν}(−ŷ·ℓ)^{−½−iν} k_{q_1q_3}(ŷ,ℓ)` (Fourier–Helgason representation: `ℓ` a boundary null vector,
`k_{q_1q_3}(ŷ,ℓ)` the mode functional in the frame `(ŷ, n_ℓ = ℓ/(−ŷ·ℓ) − ŷ, n_ℓ^⊥)`), and rotating the frame by `θ` multiplies the mode by
`e^{i(q_1+q_3)θ}`. Every `SO(2,1)`-covariant kernel of spectral parameter `ν` built on the pair `(1,3)` is a combination of them
(e.g. the structures of `REGGE_MOMENT_TTTT.md` §4: `D_1 P = 2[(1,−1)+(−1,1)]` for currents, `D_1² P ≅ 4[(2,−2)+(−2,2)] + (2/9)(0,0)`,
`D_2² P ≅ (4/9)(0,0)` on traceless tensors). The exact moments are symmetric under `(q_1,q_3) → (q_3,q_1)` and `→ (−q_1,−q_3)`.

The exact identity behind everything that follows (CHMRS (C.4)–(C.6), whose derivation is purely kinematical and holds verbatim for a
Lorentz-covariant tensor kernel, because the residual gauge group of the frame — Lorentz transformations and the dilatation
`x → λx, y → y/λ` — acts on the Cartesian components by the global Lorentz matrix only):

`∫[dη]∫dr r^{k−2} K^I(x̂,ŷ)𝒢_I = 2 ∫_{x,y>0} d³x d³y (|x||y|)^{k−4} K^I(x̂,ŷ) 𝒯_I(x,y)`,   (★)

where `𝒢_I = ∫_{x_3<x_5<x_4} d³x_5 ∫D n (…)` is the Lorentzian shadow representation of the stripped partial wave and `𝒯_I(x,y)` its integrand
in the gauge `x_5 = e = (1,0,0)`, `n = (1,1,0)`. For scalars, (★) with the explicit `𝒯` is CHMRS (C.6).

### 3.2 Numerical evaluation (`moments.py`)

Variables `(u,t)`: `√ρ = s = t/(2(α+ᾱu))`, `√ρ̄ = su` (so `e^ξ = 1/u`); a heavy block decays like `e^{−t}` uniformly in `u`
(`e^{−4m√r cosh((ξ−β)/2)} = e^{−2s(α+ᾱu)}`). Gauss–Legendre panels in `u` (72 nodes) and `t` (62 nodes); 25 null polarisation pairs per point
and a discrete Fourier transform give the modes. The harmonic functions are computed by the periodic trapezoid rule.
The block tensors on the grid are cached in `data_modes/`.

**Small-`u` region.** At fixed `t` the integrand behaves like `u^{k−1/2}` for `u → 0` (each tensor index costs `e^{ξ/2}`; for scalars it is
`u^{2Δ_φ+k−3/2}`), and the block must be evaluated at `1−|ρ̄_s| ≈ 4su → 0`, where the radial series needs `N ≳ 70/(4su)` levels. With the
`N` affordable here (`N = max(4000, 400Δ)`, capped at `10⁴` for `l = 1` and `8000` for `l = 2`; the cost grows like `N³`) the region
`u ≲ 30α/(Nt)` is dropped. Consequences (§5.1): `k = 6` moments are accurate to `10⁻⁴…10⁻⁶`, `k = 4` to `10⁻³…10⁻⁴` at the largest `m`, and `k = 2`
only for `m ≲ 20`. **This is a
limitation of the present numerics, not of the formula; it is the reason why the tables are for `k = 4, 6`.** (For gravitons the
lowest CLS-D sum rules have `k = 2, 3`.)

### 3.3 Relation to the CLS-D spinning functionals [Q + D]

[Q] CLS-D define (their eq. `spinningfunctional_definition`, §3.5; scalar case `scalarfunctional_definition`, §3.1)
`Ψ^{(a),(c)}_{k,ν,λ′}[𝒢] = ∫ d^dx_2 d^dx_4 d^dx_0 D^{d−2}z_0/vol · (⟨0|O_4(x_4^+) L[O](x_0,z_0) O_2(x_2)|0⟩^{(a)})^{−1} ∫D^{d−2}z_1 D^{d−2}z_3 ⟨P̃†_{δ_1,λ_1}(z_1)P̃†_{δ_3,λ_3}(z_3)P_{δ,λ′}(z_0)⟩^{(c)} ⟨Ω|O_4[L[O_3](x_0,z_3), L[O_1](x_0,z_1)]O_2|Ω⟩`,
`O = (d/2+iν, J_1+J_3−1, λ′)`, `δ = (d−2)/2+iν`, `δ_i = Δ_i−1`; the light transform (`light_transform_definition`) contracts `O_1, O_3` with
their own null directions, `z_1 = x_{1^+0}`, `z_3 = x_{30}`; subtraction factors `f^±_{k;J_e}(u′,v′)` (`subtraction_factor_evenjprime/oddjprime`);
the functional has signature `(−1)^{j′}` under `1 ↔ 3`. They compute only the bulk-point limit `ν, m ≫ 1` (`spinningfunctional_bulkpoint_withrtensor`)
and state: "We leave an exploration of the Regge regime for four-point functions of spinning operators to future work." In the scalar
case `Ψ_{k,ν}` "is equivalent to the `C̃_{k,ν}` functional defined in the appendix of" CHMRS.

[D] In `d = 3`: `SO(d−2)` is trivial, so `λ′` is only a parity label `j′ ∈ {0,1}`, there is one dual structure for scalars `O_2, O_4`, and the
polarisation dependence is completely explicit. With `x_1 = −y`, `x_3 = y`: `x_0 = |y| n̂`, `n̂ ⟂ y`, so `z_1 ∝ ŷ + n̂`, `z_3 ∝ ŷ − n̂`: **antipodal points
`θ_1 = φ`, `θ_3 = φ+π` of the celestial circle**. The celestial structure is `∝ |sin(θ_0−φ)|^{−δ}` (times `sgn sin` for `j′ = 1`) and the dual
structure `∝ (−z_0·x̂)^{−½+iν}` at leading order in the Regge limit. Expanding `|sin ψ|^{−δ} = Σ_n c_n e^{2inψ}`, `c_n/c_0 = (δ/2)_{|n|}/(1−δ/2)_{|n|}`
(checked numerically):

`Ψ_{k,ν,j′=0} ∝ Σ_{q_1+q_3 even} (−1)^{q_3} [c_{(q_1+q_3)/2}/c_0] Π̂^{(q_1,q_3)}_{k_r,ν}`, and the analogous sum over odd `q_1+q_3` for `j′ = 1`   (`clsd_kernel.py`).

Not derived: the `(k,ν)`-dependent normalisation, the relation between their `k` and the index `k_r` of the moment, and the effect of the
subtraction factor beyond the leading Regge moment (`f^+ ∝ (v′−u′) ∼ η/(2r)` multiplies the kernel by `η`, which shifts `ν → ν ± i` by the
Legendre recurrence, as in CHMRS (C.13)–(C.16) for `C̃_{k,ν}` versus `Π̂_{k,ν}`). So **what is applied to the exact blocks in §5 is the
complete basis `Π̂^{(q_1,q_3)}_{k,ν}`**, from which the leading-Regge action of any such functional follows; the CLS-D combination is
discussed in §9.


## 4. Task 3 — derivation of the leading heavy formula, including the normalisation

The chain has five links; each is derived and each is tested separately, and the product is tested end-to-end in §5.

### 4.1 Scalar skeleton with unequal dimensions [DVE]

Shadow representation [Q for equal dimensions: CHMRS (3.4)–(3.6); extended here]:
`K g^{a,b}_{Δ,J} = ((−1)^J/(4 b^{(Δ_34)}_{Δ,J})) ∫_{x_3<x_5<x_4} d³x_5 ∫D n S3(1,2;5,n) S3^S(5,n;3,4)`, where `S3^S` is `S3` with `(Δ,J) → (3−Δ, −1−J)`
(Lorentzian shadow), all factors taken positive, `∫D n = ∫_0^{2π}dφ_n`, `n = (1, cos φ_n, sin φ_n)`, and (`4^J S3` is CHMRS's `|⟨φφO(x_5,n)⟩|`)

**`b^{(Δ_34)}_{Δ,J} = π Γ(Δ−1) Γ(J+½) Γ((J+Δ+Δ_34)/2) Γ((J+Δ−Δ_34)/2) Γ((Δ−1−J+Δ_34)/2) Γ((Δ−1−J−Δ_34)/2) / [2 Γ(Δ−½) Γ(J+1) Γ(J+Δ) Γ(Δ−1−J)]`**,

i.e. CHMRS (3.6) with the two squared Gamma functions split by `±Δ_34/2`. Evidence for the extension: (i) the consistency with
conservation of §2.3 (exact, fixes the dependence on `Δ_34` up to a function of period 2); (ii) direct quadrature of the
representation for five sets of dimensions, including seeds with `Δ_12 ≠ Δ_34`, `J = 2, 3`: ratio `1.0004…1.0009` (quadrature error;
`out_t5_shadow_direct.txt`); (iii) the end-to-end test below, where `b^{(0)}` instead of `b^{(Δ_34)}` leaves an `O(1/m)` error.

Frame integrand. With `N(v) := −(e+v)² = 1 − 2e·v − v²` (so `N(±x) = 1 ± 2x⁰ + |x|²`), `|x|² := −x²`, one finds exactly, for dimensions
`(Δ_1,Δ_2,Δ_1,Δ_2)`,

`ℱ := s_13^{Δ_1}s_24^{Δ_2} S3_L S3^S_R = 2·4^{Δ_1+Δ_2} |x|^{2Δ_2}|y|^{2Δ_1} (1+2x·y+x²y²)^{2−Δ_1−Δ_2} 𝖵^J 𝖵′^{−1−J} N(y)^{−β}N(−y)^{−β′}N(x)^{−γ}N(−x)^{−γ′}`,
`𝖵 = N(x) − (n·y)N(x) + (n·x)N(y) − |x|²N(y)`, `𝖵′(x,y) = 𝖵(−x,−y)`, `β = (Δ+J+Δ_12)/2`, `β′ = (2−Δ−J+Δ_12)/2`, `γ = (Δ+J−Δ_12)/2`, `γ′ = (2−Δ−J−Δ_12)/2`

(for `Δ_1 = Δ_2` this is `2^{4Δ_φ+1}(|x||y|)^{2d}` times CHMRS's `𝒯`, (C.6)). For `x = X/m`, `y = Y/m`: `ℱ → 2·4^{Δ_1+Δ_2}|x|^{2Δ_2}|y|^{2Δ_1}e^{2p·x+2p̄·y}`,
`p = (Δ, J, 0)`, `p̄ = (Δ, −J, 0)` in coordinates `(t, z, ⊥)`: operators 2, 4 see `p`, the Regge pair 1, 3 sees `p̄`. With (★) and CHMRS (3.18),

`(−1)^J Π̂_{k,ν}[u^{−(Δ_1+Δ_2)/2} g^{a,b}] → (π/4)·4^{Δ_1+Δ_2}/b^{(Δ_12)}_{Δ,J} · γ_{2Δ_2+k−1}(ν) γ_{2Δ_1+k−1}(ν) · P_{−½+iν}(η_AdS)/m^{2(Δ_1+Δ_2+k−1)}`.

End-to-end test with exact blocks (`t3_scalar_mixed.py`; `R` = exact/predicted, `R_0` = the same with `b^{(0)}`), `(Δ_1,Δ_2) = (3, 1.75)`, `k = 2`:

| `(Δ, J)` | `m` | `β` | `ν` | exact `(−1)^JΠ̂` | `R` | `m²(R−1)` | `R_0` | `m(R_0−1)` |
|---|---|---|---|---|---|---|---|---|
| (14.166, 2) | 12.42 | 0.4 | 0 | 0.253173328613 | 1.00139630 | 0.215 | 0.880826 | −1.480 |
| (24.299, 4) | 22.35 | 0.4 | 0 | 1265.6065592 | 0.99871183 | −0.644 | 0.929984 | −1.565 |
| (44.565, 8) | 42.22 | 0.4 | 0 | 4.9479799891e+12 | 0.99945340 | −0.974 | 0.962426 | −1.586 |
| (85.097, 16) | 81.95 | 0.4 | 0 | 2.33620579857e+34 | 0.99983961 | −1.077 | 0.980582 | −1.591 |
| (125.629, 24) | 121.69 | 0.4 | 0 | 1.38979519581e+57 | 0.99992569 | −1.100 | 0.986914 | −1.592 |
| (14.166, 2) | 12.42 | 0.4 | 3 | 0.0241466800196 | 1.12168586 | 18.762 | 0.986633 | −0.166 |
| (44.565, 8) | 42.22 | 0.4 | 3 | 425852591933.0 | 1.01022869 | 18.231 | 0.972802 | −1.148 |
| (125.629, 24) | 121.69 | 0.4 | 3 | 1.18492524191e+56 | 1.00122956 | 18.207 | 0.988201 | −1.436 |
| (13.402, 5) | 10.55 | 1.0 | 0 | 0.674165533893 | 1.02102403 | 2.342 | 0.864239 | −1.433 |
| (45.861, 20) | 39.34 | 1.0 | 0 | 9.411709197e+13 | 0.99951058 | −0.758 | 0.955738 | −1.741 |
| (136.747, 62) | 119.94 | 1.0 | 0 | 1.2686593392e+64 | 0.99992472 | −1.083 | 0.985343 | −1.758 |
| (136.747, 62) | 119.94 | 1.0 | 3 | −4.03193563946e+62 | 1.00085588 | 12.312 | 0.986261 | −1.648 |

`m²(R−1)` converges (pure `1/m²` corrections), `m(R_0−1)` converges to a non-zero constant: the formula and `b^{(Δ_34)}` are verified,
and the corrections start at `1/m²`. Same for `(2, 1.75)` and `(1.75, 3)` (output file).

### 4.2 Shadow representation of the spinning partial wave [D; DVP by direct quadrature]

Apply `D_w ⊗ D_{w′}` to the scalar representation (the operators act on `x_1,x_2`, resp. `x_3,x_4`, and commute with the `x_5`-integral; the
boundary of the diamond does not contribute because the integrand vanishes there like a power `∼ (Δ+J)/2`):

`W_diff = ((−1)^J/4) ∫ d³x_5 D n [Σ_w c_w D_w S3_L] · [Σ_{w′} c_{w′} b_{w′}^{−1} D_{w′} S3^S_R]`, and by §2.3 the second bracket equals `(1/b^{(1−Δ_φ)}) F^S`,
`F^S := Σ_w c_w(3−Δ,−1−J) D_w S3^S` — the *same analytic family* of conserved structures at the Lorentzian-shadow point. So

`𝒯(z_1,z_3;x,y) = ((−1)^J/(4b^{(1−Δ_φ)})) · s_13^{l+1}s_24^{Δ_φ} · 𝓛(z_1) F^S(z_3)`, `𝓛 = Σ_w c_w(Δ,J)D_wS3_L = 2^{−J}L_0(Δ,J)𝔖^{cons}_{Δ,J}`, `F^S = 2^{1+J}L_0(3−Δ,−1−J)𝔖^{cons}_{3−Δ,−1−J}`.

The spinning "shadow coefficient" is therefore not an independent input: `κ = 1/b^{(1−Δ_φ)}` in the differential basis. (In the standard
basis it is `2 L_0(Δ,J)L_0(3−Δ,−1−J)/b^{(1−Δ_φ)}`, with
`L_0L_0^S = −J(J+1)(Δ−1)(Δ−2)(α−Δ_φ+1)(α+Δ_φ−1)` for `l = 1` and `¼[Δ(Δ−1)(Δ−2)(Δ−3)][(J−1)J(J+1)(J+2)][(α−Δ_φ)²−1][(α+Δ_φ)²−1]` for `l = 2`.)

### 4.3 Heavy limit of the integrand [D; DVP pointwise, `vertex.py`]

The operators act on `S3 ∼ e^{S}`, `S = O(m)`, by their symbols. In the frame, with `e = (1,0,0)`, `n = (1,1,0)`:
`V_1/2 → −(z_1·e)/2`, `H_15/V_5 → (z_1·n)/2 − z_1·e`, hence `Q_a ≡ D_11S3/S3 → −½ z_1·p` and `Q_b ≡ D_12S3/S3 → −(s_25/2) z_1·p̄`: the two operators
insert the two momenta. With the leading coefficients `c_a → −m²`, `c_b → Δ²+J²` (currents) and the identity `−m²p + (Δ²+J²)p̄ = −2ΔJm n_∥`,
`n_∥ := (−J, Δ, 0)/m` (unit, orthogonal to `p̄`, pointing from `p̄̂` to `p̂` on `H²`):

`𝓛 → S3^{(b)} s_25 · ΔJm (z_1·n_∥)` (`l = 1`),
`𝓛 → S3^{(bb)} s_25² · ½Δ²J²m² [(z_1·n_∥)² − (z_1·e_⊥)²]` (`l = 2`; from `m⁴(z·p)² − 2m²(Δ²+J²)(z·p)(z·p̄) + (Δ⁴+J⁴)(z·p̄)² = 2Δ²J²m²[2(z·n_∥)² − (z·p̄̂)²]`),

and the same for `F^S` with `z_3` (sign `(−1)^l`). So the leading tensor is `(z_1·ε_+)^l + (z_1·ε_−)^l` at `p̄̂ ∈ H²`, `ε_± = (n_∥ ± ie_⊥)/√2`, for both
operators of the Regge pair (transverse, helicity `±l`; this is Lemma V of `REGGE_MOMENT_TTTT.md`, now with its coefficient). Altogether

`𝒯 → ((−1)^{J+l}/(4b^{(1−Δ_φ)})) · 2·4^{l+1+Δ_φ} |x|^{2Δ_φ}|y|^{2(l+1)} e^{2p·x+2p̄·y} · 𝖭_l² · τ̂_l(z_1)τ̂_l(z_3)`,
`τ̂_1 = z·n_∥ = (z·ε_+ + z·ε_−)/√2`, `𝖭_1 = A_+A_−m = ½m³ sinh β`;  `τ̂_2 = (z·ε_+)² + (z·ε_−)²`, `𝖭_2 = ½A_+²A_−²m² = m⁶sinh²β/8`.

(At this order `Δ, J` may be replaced by `A_± = Δ − 3/2, J + ½`; the form with `A_±` is the one without `O(1/m)` corrections, §6.)
Pointwise check against the exact integrand, `Δ/J = 2.5`, `(X,Y)` fixed, three polarisation pairs:

| `l` | `m = 44.2` | `456.7` | `4581` | `45824` |
|---|---|---|---|---|
| 1 | 1.4901, 0.9824, 0.9169 | 1.0427, 0.99853, 0.99241 | 1.00421, 0.999856, 0.999248 | 1.000421, 0.9999856, 0.9999249 |
| 2 | 1.1535, 1.1102, 2.9784 | 1.0143, 1.0106, 1.1044 | 1.00142, 1.00106, 1.00994 | 1.000142, 1.000106, 1.000989 |

(exact/leading for three polarisation pairs; the approach is `O(1/m)` pointwise for a fixed ordered pair `(z_1,z_3)` — see §6.1 for why
this does not survive in the moments.)

### 4.4 Fourier–Laplace transform of the tensor kernels [D; DVP by direct integration, `laplace.py`, `t6_laplace_check.py`]

Master formula (from CHMRS (3.45)): for `q` future timelike and `ℓ` null,
`∫_{y>0} d³y e^{q·y}|y|^{b−3}(−ŷ·ℓ)^{−σ} = c(b,σ)|q|^{−b}(−q̂·ℓ)^{−σ}`, `c(b,σ) = 2^{b−1}√π Γ((b−σ)/2)Γ((b+σ−1)/2)`; `c(a, ½+iν) = 2^{a−1}√π γ_a(ν)`.
(Proof: radial integral, then `∫_{H²}(−q̂·ŷ)^{−b}(−ŷ·ℓ)^{−σ} = 2π∫_1^∞ dη η^{−b}P_{σ−1}(η)`.) In the Fourier–Helgason representation of §3.1 the `x`-integral
gives `c(a_x, ½−iν)|2p|^{−a_x}(−p̂·ℓ)^{−½+iν}`, and the `y`-integral of `|y|^{a_y−3}(−ŷ·ℓ)^{−½−iν}k_{q_1q_3}(ŷ,ℓ)[t]` for a constant tensor `t` is a finite sum
of derivatives of the master formula (`y_μ → ∂/∂q^μ`; `n_ℓ·ε = |y|(ℓ·ε)/(−y·ℓ) − (y·ε)/|y|`, `n_ℓ^⊥·ε = det[ε,y,ℓ]/(−y·ℓ)`). For
`t = (z_1·ε_{h_1})^l(z_3·ε_{h_3})^l`:

`∫_{y>0}d³y e^{q·y}|y|^{a−3}(−ŷ·ℓ)^{−½−iν} k_{q_1q_3}(ŷ,ℓ)[t] = c(a,½+iν)|q|^{−a}(−q̂·ℓ)^{−½−iν} · M^{q_1q_3}_{h_1h_3}(ν;a) e^{iμθ}`,  `ℓ/(−q̂·ℓ) = q̂ + cos θ n_∥ + sin θ e_⊥`,

`M = R + g_1 S`, `R, S` rational in `(ν, a)`, `g_1 = γ_{a−1}(ν)/γ_a(ν)`, `μ = 0` for helicity preserved (`h_3 = −h_1` in this convention, i.e. `ε_{h_1}^l ⊗ \bar ε_{h_1}^l`)
and `μ = ±2l` for helicity flipped. `S` is odd under the reflection `(q_i, h_i) → (−q_i, −h_i)` and cancels for parity-even blocks; the
remaining `ℓ`-integral gives `F^ν_μ(β)`. For the parity-even tensor `τ̂_l(z_1)τ̂_l(z_3)` the result is `M_pres P_{−½+iν}(cosh β) + M_flip F^ν_{2l}(β)` with
(`𝖽_b := ν² + (a−b)²`, `a = a_y`; `out_mforms.txt`; the sums over `h` are included, and for `l = 1` the factor ½ of `τ̂_1τ̂_1` is in `T_1`):

| `l` | mode | `M_pres` | `M_flip` |
|---|---|---|---|
| 1 | `(1,−1)` | `(ν² + 2(a−2)² + ¼)/(4𝖽_{5/2})` | `(½+iν)(3/2+iν)/(4𝖽_{5/2})` |
| 1 | `(1,1)` | `(½−iν)(3/2−iν)/(4𝖽_{5/2})` | `(ν² + 2(a−3)² + ¼)/(4𝖽_{5/2})` |
| 1 | `(0,0)` | `(2a − ν² − 17/4)/𝖽_{5/2}` | `(½+iν)(3/2+iν)/𝖽_{5/2}` |
| 2 | `(2,−2)` | `N_a/(512 𝖽_{5/2}𝖽_{9/2})` | `(½+iν)_4/(32 𝖽_{5/2}𝖽_{9/2})` |
| 2 | `(2,2)` | `(½−iν)_4/(32 𝖽_{5/2}𝖽_{9/2})` | `N_b/(512 𝖽_{5/2}𝖽_{9/2})` |
| 2 | `(0,0)` | `9(128a² − 128aν² − 800a + 16ν⁴ + 360ν² + 1113)/(128 𝖽_{5/2}𝖽_{9/2})` | `9(½+iν)_4/(8 𝖽_{5/2}𝖽_{9/2})` |

`N_a = 128a⁴ − 1536a³ + 128a²ν² + 6432a² − 512aν² − 10880a + 16ν⁴ + 616ν² + 6297`,
`N_b = 128a⁴ − 2560a³ + 128a²ν² + 19488a² − 1280aν² − 66880a + 16ν⁴ + 3304ν² + 87225`.
(`(½+iν)_{2l}F^ν_{2l} = P^{2l}_{−½+iν}`, the function called `P^{(2l)}` in `REGGE_MOMENT_TTTT.md`; its table for `D_1`, `D_1²`, `D_2²` is reproduced
exactly: e.g. `M⁰_{D_1²} = 4·2R^{(2,−2)}_{+−} + (2/9)M^{(0,0)}`.) Direct `H²` integration of eight entries, including `S ≠ 0` and all stress-tensor
rows, agrees to 7–10 digits (`out_t6_laplace_check.txt`).

### 4.5 Result [DVE]

Combining (★), §4.2–4.4 (`2·(1/4)·2 = 1`, `c(a_x,·)c(a_y,·)|2p|^{−a_x}|2p̄|^{−a_y} = (π/4)γ_{a_x}γ_{a_y}m^{−a_x−a_y}`):

**`(−1)^{J+l} Π̂^{(q_1,q_3)}_{k,ν}[s_13^{l+1}s_24^{Δ_φ}W_diff] = (π/4)·4^{l+1+Δ_φ}/b^{(1−Δ_φ)}_{Δ,J} · γ_{a_x}(ν)γ_{a_y}(ν)/m^{a_x+a_y} · |T_l| · [M^{q_1q_3}_pres(ν;a_y) P_{−½+iν}(cosh β) + M^{q_1q_3}_flip(ν;a_y) F^ν_{2l}(β)] · (1 + O(1/m²))`**,

`a_x = 2Δ_φ + k − 1`, `a_y = 2l + k + 1`, `|T_1| = ½(½m³ sinh β)²`, `|T_2| = (m⁶ sinh²β/8)²`, `m, β` from `(α, ᾱ)`.

Equivalently, for the partial wave of the standard structures `𝒲[𝔖^{cons},𝔖^{cons}] = W_diff/(2^J L_0²)` on the real section: replace `|T_l|` by
`|T_l|/(2^J L_0²) → 2^{−J} e^{lβ}·{½, 1}` (`l = 1, 2`; the factor `e^{lβ} = (ᾱ/α)^l` is the price of normalising by the `V_1^lV_5^J` coefficient).
The contribution of an operator `O` to the double discontinuity is this times `λ_{OφO}² · 2 sin²(π(τ−Δ_O−Δ_φ−l)/2)`.

Structure: (i) helicity preserved at the Regge pair ↔ `P_{−½+iν}(cosh β)`, the scalar harmonic function, with no other `β`-dependence;
(ii) helicity flipped ↔ `F^ν_{2l}(β) ∝ P^{2l}_{−½+iν}(cosh β)`; (iii) which combination a functional sees is fixed by its kernel through the
rational functions `M`; the transform is not local on `H²` (a purely transverse kernel at `ŷ` feeds both components);
(iv) all `ν`-dependence other than `γ_{a_x}γ_{a_y}` is rational; `M_pres^{(l,−l)} > 0` for real `ν`.


## 5. Task 3 — end-to-end comparison: exact blocks versus the derived formula

Setting: `Δ_φ = 1.75`; `(Δ, J)` chosen at fixed `β ≈ 0.4, 1, 2` with `m ≈ 10, 20, 40, 70, 100` (`β ≈ 2`: up to 70); `k = 2, 4, 6`; `ν = 0, 1, 3`;
all `(2l+1)²` modes. "Measured" = the tensor Regge moment of the exact partial wave `s_13^{l+1}s_24^{Δ_φ}W_diff` (no large-`Δ`
approximation anywhere on this side); "predicted" = the formula of §4.5 with the exact `b^{(1−Δ_φ)}_{Δ,J}`. **Nothing is fitted.**
Files: `out_t4_spinning_e2e_l{1,2}.txt` (all modes), `out_t9_tables.txt` (absolute values), `out_t9b_summary.txt`, `out_t7_corrections_l{1,2}.txt`.

### 5.1 Numerical accuracy of the measured side

* *Dropped small-`u` region* (§3.2). For every entry the script computes `tail` = the contribution of the region added between
  the truncations `N/1.5` and `N`, extrapolated with the power law of the integrand. For the largest `m` (≈ 100): `k = 6`: `10⁻⁵…10⁻⁶`
  (`l = 1`), `10⁻⁴…10⁻⁵` (`l = 2`); `k = 4`: `2·10⁻⁴…2·10⁻³`; `k = 2`: `10⁻²…10⁻¹` (useless). Entries with `tail > 10⁻⁴` are marked `*`
  (`> 10⁻³`: `**`) below and are not used in the fits of §6.
* *Truncation and grid study* (`out_t4b_convergence.txt`): at `β ≈ 1`, `m ≈ 10, 20`, comparing `N = 200Δ, 400Δ` with `800Δ` (`m ≈ 10`) or `200Δ` with `400Δ` (`m ≈ 20`), and the base quadrature
  grid with a grid of 1.4× the nodes in each direction (modes `(l,−l)`, `(l,l)`, `ν = 0, 3`). Truncation: `k = 6`: `≤ 3·10⁻⁷` already at
  `N = 200Δ`; `k = 4`: `≤ 2·10⁻⁴` at `200Δ`, `≤ 10⁻⁵` at `400Δ`; `k = 2`: `2·10⁻² … 2·10⁻¹` at `200Δ`, `≤ 3·10⁻²` at `400Δ` (the `k = 2` moments converge
  slowly in `N`, as expected from the `u^{3/2}` integrand). Grid: `k = 4, 6`: `≤ 10⁻⁶` except `l = 2`, `m ≈ 10` (`4·10⁻⁴`); at `m ≈ 20`: `≤ 5·10⁻⁷`.
  (The production runs use `N = max(4000, 400Δ)` capped at `10⁴`/`8000`, so the relative truncation is largest at the largest `m`;
  this is what the per-entry `tail` estimate measures.)
* *Round-off*: double versus quad evaluation of the tensor at a point of the dispersive region: `10⁻¹²` (`l = 1`), `10⁻⁸` (`l = 2`) (§2.3).
* *Internal consistency*: the exact moments satisfy `(q_1,q_3) ↔ (q_3,q_1)` and `↔ (−q_1,−q_3)` to the accuracy above, and imaginary
  parts are `< 10⁻⁵` of the real parts.

### 5.2 Measured and predicted values

Sample of absolute values (mode `(l,−l)`, `k = 6`; the complete list for `k = 4, 6`, both transverse modes, is `out_t9_tables.txt`):

| `l` | `(Δ, J)` | `m` | `β` | `ν` | measured `(−1)^{J+l}Π̂^{(l,−l)}_{6,ν}` | predicted | measured/predicted |
|---|---|---|---|---|---|---|---|
| 1 | (14.1662, 2) | 12.42 | 0.4 | 0 | 0.186776789 | 0.319940322 | 0.583786 |
| 1 | (44.5652, 8) | 42.22 | 0.4 | 0 | 5.6042332e+12 | 5.91197816e+12 | 0.947946 |
| 1 | (105.3630, 20) | 101.82 | 0.4 | 0 | 5.53160531e+45 | 5.58341276e+45 | 0.990721 |
| 1 | (105.3630, 20) | 101.82 | 0.4 | 3 | 1.03611212e+45 | 1.0388919e+45 | 0.997324 |
| 1 | (115.1076, 52) | 100.75 | 1.0 | 0 | 5.92675959e+52 | 6.00119806e+52 | 0.987596 |
| 1 | (115.1076, 52) | 100.75 | 1.0 | 3 | −3.24873762e+51 | −3.32548355e+51 | 0.976922 |
| 1 | (109.8254, 82) | 70.20 | 2.0 | 0 | 1.63940048e+52 | 1.73790894e+52 | 0.943318 |
| 2 | (14.1662, 2) | 12.42 | 0.4 | 0 | 2369.04106 | 6306.76749 | 0.375635 |
| 2 | (44.5652, 8) | 42.22 | 0.4 | 0 | 1.42943091e+19 | 1.55735037e+19 | 0.917861 |
| 2 | (105.3630, 20) | 101.82 | 0.4 | 0 | 4.90293003e+53 | 4.97614453e+53 | 0.985287 |
| 2 | (105.3630, 20) | 101.82 | 0.4 | 3 | 9.82712214e+52 | 9.87523294e+52 | 0.995128 |
| 2 | (115.1076, 52) | 100.75 | 1.0 | 0 | 4.13894716e+61 | 4.19477108e+61 | 0.986692 |
| 2 | (115.1076, 52) | 100.75 | 1.0 | 3 | −3.01901383e+60 | −3.09013508e+60 | 0.976984 |
| 2 | (109.8254, 82) | 70.20 | 2.0 | 0 | 2.56895821e+61 | 2.72389384e+61 | 0.943120 |

(The moments range over 60 orders of magnitude; signs and the sign change of `P_{−½+3i}(cosh β)` between `β = 0.4` and `1` are reproduced.)

**Helicity components.** At leading order each mode is `pref·[M_pres X + M_flip Y]` with `X = P_{−½+iν}(cosh β)` (helicity preserved at the
spinning pair) and `Y = F^ν_{2l}(β)` (flipped). The measured `X, Y` are obtained from the two transverse modes `(l,−l)`, `(l,l)` by solving the
`2×2` system with the *predicted* `M`'s; `R_pres = X_meas/X`, `R_flip = Y_meas/Y`. (The mode `(l,−l)` is dominated by `X`, `(l,l)` by `Y`.)

**l = 1, k = 6:  R_pres / R_flip  (measured/predicted)**

| β ≈ 0.4, (Δ,J) → | (14.17, 2) | (24.30, 4) | (44.57, 8) | (74.96, 14) | (105.36, 20) |
|---|---|---|---|---|---|
| m | 12.4 | 22.4 | 42.2 | 72.0 | 101.8 |
| ν = 0 | 0.5838 / 0.4291 | 0.8325 / 0.7437 | 0.9479 / 0.9162 | 0.9816 / 0.9699 | 0.9907 / 0.9847 |
| ν = 1 | 0.6036 / 0.4390 | 0.8433 / 0.7505 | 0.9517 / 0.9188 | 0.9830 / 0.9709 | 0.9914 / 0.9853 |
| ν = 3 | 0.8073 / 0.5282 | 0.9436 / 0.8073 | 0.9852 / 0.9399 | 0.9950 / 0.9788 | 0.9975 / 0.9893 |

| β ≈ 1.0, (Δ,J) → | (13.40, 5) | (24.22, 10) | (45.86, 20) | (80.48, 36) | (115.11, 52) |
|---|---|---|---|---|---|
| m | 10.6 | 20.1 | 39.3 | 70.0 | 100.7 |
| ν = 0 | 0.4355 / 0.3105 | 0.7568 / 0.6692 | 0.9237 / 0.8906 | 0.9747 / 0.9631 | 0.9876 / 0.9818 |
| ν = 1 | 0.4912 / 0.3282 | 0.7946 / 0.6847 | 0.9378 / 0.8971 | 0.9796 / 0.9654 | 0.9900 / 0.9830 |
| ν = 3 | 0.0211 / 0.5408 | 0.5810 / 0.8422 | 0.8737 / 0.9572 | 0.9588 / 0.9865 | 0.9799 / 0.9935 |

| β ≈ 2.0, (Δ,J) → | (16.60, 11) | (32.36, 23) | (63.87, 47) | (109.83, 82) |
|---|---|---|---|---|
| m | 9.8 | 20.0 | 40.4 | 70.2 |
| ν = 0 | 0.2077 / 0.1523 | 0.5685 / 0.5041 | 0.8468 / 0.8159 | 0.9433 / 0.9308 |
| ν = 1 | 0.4195 / 0.1849 | 0.8026 / 0.5516 | 0.9558 / 0.8413 | 0.9869 / 0.9415 |
| ν = 3 | -0.5634 / -0.0532 | -0.0998 / 0.4555 | 0.6085 / 0.8359 | 0.8582 / 0.9438 |

**l = 1, k = 4:  R_pres / R_flip  (measured/predicted)**

| β ≈ 0.4, (Δ,J) → | (14.17, 2) | (24.30, 4) | (44.57, 8) | (74.96, 14) | (105.36, 20) |
|---|---|---|---|---|---|
| m | 12.4 | 22.4 | 42.2 | 72.0 | 101.8 |
| ν = 0 | 0.7733 / 0.6400 | 0.9185 / 0.8587 | 0.9758 / 0.9563 | 0.9915 / 0.9830* | 0.9956 / 0.9864** |
| ν = 1 | 0.7942 / 0.6520 | 0.9274 / 0.8647 | 0.9786 / 0.9585 | 0.9926 / 0.9851* | 0.9963 / 0.9908* |
| ν = 3 | 1.0016 / 0.7573 | 1.0074 / 0.9140 | 1.0029 / 0.9747 | 1.0011 / 0.9911 | 1.0007 / 0.9949* |

| β ≈ 1.0, (Δ,J) → | (13.40, 5) | (24.22, 10) | (45.86, 20) | (80.48, 36) | (115.11, 52) |
|---|---|---|---|---|---|
| m | 10.6 | 20.1 | 39.3 | 70.0 | 100.7 |
| ν = 0 | 0.6871 / 0.5645 | 0.8864 / 0.8270 | 0.9669 / 0.9476 | 0.9892 / 0.9819* | 0.9943 / 0.9881** |
| ν = 1 | 0.7492 / 0.5889 | 0.9149 / 0.8403 | 0.9760 / 0.9523 | 0.9923 / 0.9843* | 0.9965 / 0.9917* |
| ν = 3 | 0.3673 / 0.8590 | 0.7989 / 0.9675 | 0.9458 / 0.9929 | 0.9829 / 0.9980* | 0.9907 / 0.9987* |

| β ≈ 2.0, (Δ,J) → | (16.60, 11) | (32.36, 23) | (63.87, 47) | (109.83, 82) |
|---|---|---|---|---|
| m | 9.8 | 20.0 | 40.4 | 70.2 |
| ν = 0 | 0.4803 / 0.4044 | 0.7803 / 0.7291 | 0.9323 / 0.9130* | 0.9746 / 0.9646** |
| ν = 1 | 0.8069 / 0.4667 | 0.9757 / 0.7742 | 1.0011 / 0.9311 | 1.0034 / 0.9759** |
| ν = 3 | -0.6265 / 0.2606 | 0.3252 / 0.7665 | 0.8122 / 0.9452* | 0.9281 / 0.9771** |

**l = 2, k = 6:  R_pres / R_flip  (measured/predicted)**

| β ≈ 0.4, (Δ,J) → | (14.17, 2) | (24.30, 4) | (44.57, 8) | (74.96, 14) | (105.36, 20) |
|---|---|---|---|---|---|
| m | 12.4 | 22.4 | 42.2 | 72.0 | 101.8 |
| ν = 0 | 0.3756 / 0.1373 | 0.7394 / 0.5335 | 0.9179 / 0.8345 | 0.9708 / 0.9391 | 0.9853 / 0.9685* |
| ν = 1 | 0.3937 / 0.1400 | 0.7529 / 0.5382 | 0.9230 / 0.8369 | 0.9728 / 0.9402 | 0.9863 / 0.9693 |
| ν = 3 | 0.5946 / 0.1650 | 0.8853 / 0.5785 | 0.9705 / 0.8566 | 0.9902 / 0.9480 | 0.9951 / 0.9736 |

| β ≈ 1.0, (Δ,J) → | (13.40, 5) | (24.22, 10) | (45.86, 20) | (80.48, 36) | (115.11, 52) |
|---|---|---|---|---|---|
| m | 10.6 | 20.1 | 39.3 | 70.0 | 100.7 |
| ν = 0 | 0.4122 / 0.2066 | 0.7426 / 0.5772 | 0.9184 / 0.8519 | 0.9728 / 0.9491 | 0.9867 / 0.9747 |
| ν = 1 | 0.4717 / 0.2159 | 0.7840 / 0.5883 | 0.9340 / 0.8571 | 0.9783 / 0.9511 | 0.9894 / 0.9758 |
| ν = 3 | -0.0962 / 0.3122 | 0.5248 / 0.6882 | 0.8558 / 0.9011 | 0.9529 / 0.9671 | 0.9770 / 0.9839 |

| β ≈ 2.0, (Δ,J) → | (16.60, 11) | (32.36, 23) | (63.87, 47) | (109.83, 82) |
|---|---|---|---|---|
| m | 9.8 | 20.0 | 40.4 | 70.2 |
| ν = 0 | 0.2213 / 0.1178 | 0.5740 / 0.4603 | 0.8470 / 0.7904 | 0.9431 / 0.9198 |
| ν = 1 | 0.4532 / 0.1360 | 0.8118 / 0.4905 | 0.9572 / 0.8080 | 0.9872 / 0.9275 |
| ν = 3 | -0.6320 / 1.0284 | -0.1052 / 1.4408 | 0.6054 / 1.2297 | 0.8565 / 1.0922 |

**l = 2, k = 4:  R_pres / R_flip  (measured/predicted)**

| β ≈ 0.4, (Δ,J) → | (14.17, 2) | (24.30, 4) | (44.57, 8) | (74.96, 14) | (105.36, 20) |
|---|---|---|---|---|---|
| m | 12.4 | 22.4 | 42.2 | 72.0 | 101.8 |
| ν = 0 | 0.4941 / 0.2359 | 0.8149 / 0.6443 | 0.9448 / 0.8785* | 0.9807 / 0.9219** | 0.9902 / 0.8527** |
| ν = 1 | 0.5154 / 0.2379 | 0.8273 / 0.6486 | 0.9490 / 0.8843* | 0.9822 / 0.9454** | 0.9911 / 0.9178** |
| ν = 3 | 0.7484 / 0.2608 | 0.9474 / 0.6843 | 0.9882 / 0.9011 | 0.9962 / 0.9634 | 0.9982 / 0.9809** |

| β ≈ 1.0, (Δ,J) → | (13.40, 5) | (24.22, 10) | (45.86, 20) | (80.48, 36) | (115.11, 52) |
|---|---|---|---|---|---|
| m | 10.6 | 20.1 | 39.3 | 70.0 | 100.7 |
| ν = 0 | 0.6417 / 0.4248 | 0.8679 / 0.7484 | 0.9612 / 0.9202* | 0.9873 / 0.9701** | 0.9935 / 0.9752** |
| ν = 1 | 0.7121 / 0.4390 | 0.9012 / 0.7589 | 0.9719 / 0.9245 | 0.9910 / 0.9740* | 0.9958 / 0.9829** |
| ν = 3 | 0.1666 / 0.5863 | 0.7314 / 0.8515 | 0.9270 / 0.9574 | 0.9768 / 0.9860* | 0.9881 / 0.9924* |

| β ≈ 2.0, (Δ,J) → | (16.60, 11) | (32.36, 23) | (63.87, 47) | (109.83, 82) |
|---|---|---|---|---|
| m | 9.8 | 20.0 | 40.4 | 70.2 |
| ν = 0 | 0.5013 / 0.3439 | 0.7833 / 0.6923 | 0.9319 / 0.8962* | 0.9745 / 0.9568** |
| ν = 1 | 0.8560 / 0.3864 | 0.9823 / 0.7231 | 1.0021 / 0.9097 | 1.0034 / 0.9665** |
| ν = 3 | -0.7051 / 2.0282 | 0.3207 / 1.5250 | 0.8077 / 1.1737* | 0.9310 / 1.0627** |

(`*`: estimated numerical error `> 10⁻⁴`; `**`: `> 10⁻³`. The `k = 2` tables are in `out_t9b_summary.txt`: there `R_pres` behaves like the others for `m ≲ 40` while `R_flip` is dominated by numerical error.)

Reading: for both currents and stress tensors, both helicity components, all three `β` (including `β = 1, 2`, far from flat space),
`ν = 0, 1, 3`, and `m` from 10 to 100, the ratio tends to 1; at `m ≈ 100`, `k = 6`: `|R_pres − 1| ≤ 2.3 %`, `|R_flip − 1| ≤ 3.2 %`, and the
deviation decreases like `1/m²` (`m²(R−1)` converges, §6; e.g. `l = 2`, `β ≈ 0.4`, `ν = 0`: `m²(R_pres−1) = −96.3, −130.2, −146.4, −151.2, −152.5`).
The slow start is due to the size of the `1/m²` coefficient (`∼ −(k+2l)²`-ish, `−100…−700`), not to a wrong leading term. Where the
numerics allow (`k = 2`: `m ≲ 40` for `R_pres`), `k = 2` shows the same approach with a much smaller coefficient
(`R_pres = 0.900, 0.966, 0.988` at `m = 12, 22, 42`, `l = 1`, `β ≈ 0.4`, `ν = 0`).

**Overdetermination.** The other `(2l+1)² − 2` modes (longitudinal and mixed, e.g. `(0,0)`, `(1,0)`, `(2,−1)`), computed from the measured
`X, Y` and the predicted `M`'s, agree with the exact moments up to `O(1/m²)` relative deviations (column "other modes" of
`out_t7_corrections_*.txt`: `1.3·10⁻², 2.6·10⁻³, 2.8·10⁻⁴, 3.5·10⁻⁵, 5.5·10⁻⁶` for `l = 1`, `β ≈ 0.4`, `k = 6`, `ν = 0`, `m = 12…102` —
decreasing faster than `1/m²` because the leading `1/m²` correction is mode-independent to a good approximation; for `l = 2`:
`1.9·10⁻², 5.1·10⁻³, 5.8·10⁻⁴, 7.4·10⁻⁵`). E.g. all eight independent modes of `l = 2`, `(Δ,J) = (74.96, 14)`, `k = 6`, `ν = 0` have
measured/predicted `= 0.97085, 0.95776 [(2,2)], 0.97092, 0.97089, 0.97065, 0.97077, 0.97086, 0.97091`. This tests the whole table of `M`
(including the non-transverse rows) and the statement that only two functions of `β` appear.

**Sensitivity controls** (what the ratio would be with one ingredient altered; `k = 6`, `ν = 0`, mode `(l,−l)`):

| l | (Δ,J) | m | derived | b^(0) instead of b^(1−Δφ) | (Δ,J) instead of (Δ−3/2, J+1/2) in T_l | M_pres only (flip term dropped) |
|---|---|---|---|---|---|---|
| 1 | (74.96, 14) | 72.0 | 0.98160 | 0.97381 | 1.01124 | 0.98171 |
| 1 | (105.36, 20) | 101.8 | 0.99072 | 0.98515 | 1.01145 | 0.99083 |
| 1 | (80.48, 36) | 70.0 | 0.97467 | 0.96588 | 0.96493 | 0.97529 |
| 1 | (115.11, 52) | 100.7 | 0.98760 | 0.98140 | 0.98061 | 0.98822 |
| 1 | (109.83, 82) | 70.2 | 0.94332 | 0.93173 | 0.92895 | 0.94512 |
| 2 | (74.96, 14) | 72.0 | 0.97085 | 0.96314 | 1.03038 | 0.97085 |
| 2 | (105.36, 20) | 101.8 | 0.98529 | 0.97975 | 1.02695 | 0.98529 |
| 2 | (80.48, 36) | 70.0 | 0.97284 | 0.96407 | 0.95350 | 0.97285 |
| 2 | (115.11, 52) | 100.7 | 0.98669 | 0.98050 | 0.97279 | 0.98669 |
| 2 | (109.83, 82) | 70.2 | 0.94312 | 0.93153 | 0.91461 | 0.94314 |

Each alteration is an `O(1/m)` effect, i.e. it would show up as a non-zero limit of `m(R−1)`; §6.1 bounds that limit by `0.008` for the
helicity-preserving component, against `m(R_{altered}−R) = −0.57` (the `b^{(0)}` column, `l = 1`, `m = 102`) or `+2.1` (`l = 1`) / `+4.2` (`l = 2`) (the `(Δ,J)` column, `β ≈ 0.4`). Factors of 2, signs
and the powers of `m` and `sinh β` in `|T_l|` are tested at `O(1)`.


## 6. Task 4 — first corrections

Notation: `R − 1 = c_1/m + c_2/m² + …` for `R_pres` and `R_flip`, with `m² = αᾱ`, `cosh β = (α²+ᾱ²)/(2αᾱ)` and the exact `b^{(1−Δ_φ)}`.

### 6.1 `O(1/m)`: absent, for currents and for stress tensors [D; verified pointwise and end-to-end]

*Argument.* By §4.2, `b·𝒯(z_1,z_3;x,y) ∝ s_13^{l+1}s_24^{Δ_φ} 𝓛_{Δ,J}(z_1;−y, x/x²) 𝓛_{3−Δ,−1−J}(z_3; y, −x/x²)`: the two factors are the *same* analytic
function, evaluated at Weyl-reflected labels and at reflected points. Hence `b𝒯` is invariant (up to the overall sign `(−1)^l`) under
`(x, y; A_+, A_−; z_1, z_3) → (−x, −y; −A_+, −A_−; z_3, z_1)`, `A_+ = Δ−3/2`, `A_− = J+½`. In the heavy expansion at fixed `X = mx`, `Y = my`, `β`, this
reflection is `m → −m`. Therefore `b𝒯 = 𝒯_0 [1 + t_1/m + t_2/m² + …]` with `t_n(z_1,z_3) = (−1)^n t_n(z_3,z_1)`: **the `O(1/m)` term is antisymmetric
under the exchange of the two polarisations**, and it cancels in every kernel symmetric under `1 ↔ 3`. Since the exact moments obey
`Π̂^{(q_1,q_3)} = Π̂^{(q_3,q_1)}`, it cancels in all of them — *provided* the expansion variables are Weyl-odd (`A_±`, or `m, β` built from
`α = A_+−A_−`, `ᾱ = A_++A_−`) and `1/b` is kept exact. In terms of `Δ, J` themselves there are `O(1/m)` terms (control table of §5.2). The
same argument applies verbatim to any number of spinning legs, hence to `⟨TTTT⟩`, for kernels symmetric under the simultaneous
exchange `1 ↔ 3`, `2 ↔ 4`.

*Pointwise test* (`out_vertex.txt`): `m·(exact/leading − 1)` at fixed `(X,Y)`, `J/Δ = 0.4`, for `(z_1,z_3)`, for `(z_3,z_1)`, and for the symmetrised tensor:

| `l` | polarisations | `m` | `(z_1,z_3)` | `(z_3,z_1)` | symmetrised | `m²·(sym − 1)` |
|---|---|---|---|---|---|---|
| 1 | (0.3, 1.9) | 457 / 4581 / 45824 | 1.278 / 1.296 / 1.298 | −1.248 / −1.293 / −1.298 | 0.0149 / 0.00149 / 0.000149 | 6.82 / 6.83 / 6.83 |
| 1 | (4.0, 5.1) | 457 / 4581 / 45824 | −4.855 / −4.856 / −4.856 | 4.905 / 4.861 / 4.857 | 0.0249 / 0.00248 / 0.000247 | 11.39 / 11.34 / 11.33 |
| 2 | (0.3, 1.9) | 457 / 4581 / 45824 | 1.599 / 1.609 / 1.611 | −1.628 / −1.612 / −1.611 | −0.0145 / −0.00148 / −0.000148 | −6.61 / −6.76 / −6.78 |
| 2 | (4.0, 5.1) | 457 / 4581 / 45824 | 37.97 / 36.57 / 36.43 | −37.65 / −36.53 / −36.43 | 0.158 / 0.0152 / 0.00151 | 72.3 / 69.6 / 69.3 |

*End-to-end test* (`out_t7_corrections_l*.txt`). Three-parameter fits to the three largest `m` (`k = 6`; ansatz (A): `c_1/m + c_2/m² + c_4/m⁴`,
ansatz (B): `c_2/m² + c_4/m⁴ + c_6/m⁶`):

| `l` | `β` | `ν` | `c_1` pres (A) | `c_1` flip (A) | `c_2` pres (B) | `c_2` flip (B) |
|---|---|---|---|---|---|---|
| 1 | 0.4 | 0 / 1 / 3 | −0.0016 / −0.0009 / +0.0008 | −0.016 / −0.006 / −0.002 | −96.94 / −89.53 / −25.47 | −160.3 / −154.5 / −111.9 |
| 1 | 1.0 | 0 / 1 / 3 | −0.0068 / −0.0013 / −0.0019 | −0.026 / −0.011 / +0.003 | −127.43 / −102.01 / −205.37 | −187.8 / −174.9 / −66.0 |
| 2 | 0.4 | 0 / 1 / 3 | −0.027 / −0.0017 / +0.0017 | −0.20 / −0.12 / −0.018 | −153.85 / −143.50 / −50.01 | −327.5 / −327.3 / −277.9 |
| 2 | 1.0 | 0 / 1 / 3 | −0.0076 / −0.0021 / −0.0038 | −0.074 / −0.039 / −0.009 | −136.87 / −108.60 / −235.23 | −263.1 / −250.9 / −165.5 |

For the helicity-preserving component `|c_1| ≤ 0.008` in 11 of 12 cases (to be compared with `|c_2|/m ≈ 1…2` at `m = 100`, and with
`c_1 ≈ −0.57`, `+2.1`, `+4.2` for the altered formulas of §5.2); the fitted `c_1` is of the size of the neglected `c_6/m⁵` term and of the
numerical error (`10⁻⁵ m`). For the flip component of `l = 2` at `β ≈ 0.4` the fit is contaminated by higher orders (the flip
expansion parameter there is `1/J²`, `J ≤ 20`) and by the larger numerical error of the `(2,2)` mode. At `β ≈ 2` only four values
`m ≤ 70` exist and `c_2 ≈ −300…−740`, so the fits are not asymptotic (`|c_1| ≲ 1`, inconclusive there).
**Conclusion: no `O(1/m)` term for stress tensors (nor for currents), in the variables `(m, β)` of CHMRS built from `(α, ᾱ)`.**

### 6.2 `O(1/m²)` coefficients [measured; analytic value not derived]

`c_2` from the three-parameter fit (B) where available (otherwise the two-point even fit), `k = 6` and `k = 4`:

| `β` | `k` | `ν` | `l = 1`: `c_2^pres` | `c_2^flip` | scalar skeleton (2, 1.75) | `l = 2`: `c_2^pres` | `c_2^flip` | scalar skeleton (3, 1.75) |
|---|---|---|---|---|---|---|---|---|
| 0.4 | 6 | 0 | −96.9 | −160.3 | −90.95 | −153.9 | −327.5 | −72.28 |
| 0.4 | 6 | 1 | −89.5 | −154.5 | −84.56 | −143.5 | −327.3 | −66.05 |
| 0.4 | 6 | 3 | −25.5 | −111.9 | −30.32 | −50.0 | −277.9 | −13.47 |
| 1.0 | 6 | 0 | −127.4 | −187.8 | −126.66 | −136.9 | −263.1 | −100.23 |
| 1.0 | 6 | 1 | −102.0 | −174.9 | −102.38 | −108.6 | −250.9 | −77.27 |
| 1.0 | 6 | 3 | −205.4 | −66.0 | −197.89 | −235.2 | −165.5 | −162.09 |
| 2.0 | 6 | 0 | −295 | −363 | −293.6 | −297 | −425 | −230.4 |
| 2.0 | 6 | 1 | −60 | −304 | −63.0 | −59 | −381 | −14.9 |
| 2.0 | 6 | 3 | −730 | −281 | −723.3 | −741 | +498 | −618.1 |
| 0.4 | 4 | 0 | −44.0 | −80.7 | −38.12 | — | — | −23.77 |
| 0.4 | 4 | 3 | +6.0 | −46.5 | +0.95 | −18.6 | −198 | +13.48 |
| 1.0 | 4 | 1 | −38.1 | −77.1 | −38.55 | −44.6 | −124 | −19.61 |
| 1.0 | 4 | 3 | −84.7 | −10.3 | −76.92 | −114.4 | −68 | −47.29 |

("Scalar skeleton": the `1/m²` coefficient of the exact scalar block with external dimensions `(l+1, Δ_φ, l+1, Δ_φ)` relative to the scalar
formula of §4.1, same `k, ν, β`; `out_t7b_scalar_skeleton.txt`. Estimated errors: `≲ 0.3` for `β ≤ 1`, `k = 6`; a few units for `β = 2` and for
`k = 4`; the `l = 2`, `k = 4` entries rely on `m ≤ 42`.) The coefficients are large, grow with `k` and `β`, oscillate with `ν` through
`P′/P`, and most of the helicity-preserving coefficient is that of the scalar skeleton, i.e. it is *not* specific to spin.

### 6.3 Spin-dependent part, helicity preserved; comparison with the flat-space helicity term

Define `D_pres := c_2^pres − c_2^{skeleton}`. Two contributions can be predicted:

* **(pol) [D]** The exact partial wave vanishes identically for `J < l` (and at the reflected points), because the conserved structure
  does (`L_0 ∝ J(J−1)…`, §1.2); the exact normalisation is `Π_{i=1}^{l}[(J+½)² − (i−½)²]·Π_{i=1}^{l}[(Δ−3/2)² − (i−½)²]` where the leading formula has
  `(A_−A_+)^{2l}` (§4.2). This gives `pol = −κ_l [1/sinh²(β/2) + 1/cosh²(β/2)]`, `κ_1 = ¼`, `κ_2 = 5/2`.
* **(hel) [I, from flat space]** The helicity-dependent `1/m²` term of flat partial waves, `REGGE_MOMENT_TTTT.md` §6.2,
  `(1/m²)·2λ²(ν/β)J_1(νβ)` relative to `J_0(νβ)` with helicity `λ = l`. Its natural curved form is `hel = −2l² P′_{−½+iν}(η)/P_{−½+iν}(η)`, `η = cosh β`
  (derivative with respect to `η`), which is the first-order effect of `η_AdS → η_AdS − 2l²/m²`, i.e. of `A_±² → A_±² − l²` in the argument of the
  harmonic function.

Test (`out_t7c_spin_dependence.txt`): least-squares fit of the *form* `D_pres − pol = a_0 + a_1 cosh β + a_2 P′/P + a_3 sinh²β P′/P` to the six
points `β ∈ {0.4, 1}`, `ν ∈ {0, 1, 3}` at `k = 6` (prediction of (hel): `a_2 = −2l²`; nothing is known a priori about `a_0, a_1, a_3`):

| `l` | `a_2` predicted | `a_2` measured | `a_0` | `a_1` | `a_3` | residuals | coefficient of `pol` if left free (9 points; prediction 1) |
|---|---|---|---|---|---|---|---|
| 1 | −2 | **−2.01** | 0.07 | 0.10 | −0.00 | `≤ 0.09` | 1.04 |
| 2 | −8 | **−8.00** | 0.15 | −17.26 | −2.02 | `≤ 0.07` | 1.002 |

(with all nine points including `β = 2`: `a_2 = −1.90`, `−7.94`, residuals `≤ 1.1`; at `k = 4`, `l = 1`: `a_2 = −2.01`.) So:

**`c_2^pres(l=1) = c_2^{skeleton(2,Δ_φ)} − ¼[1/sinh²(β/2) + 1/cosh²(β/2)] − 2 P′/P + 0.2(1)`**,
**`c_2^pres(l=2) = c_2^{skeleton(3,Δ_φ)} − (5/2)[1/sinh²(β/2) + 1/cosh²(β/2)] − 8 P′/P − 2.0(1)·sinh²β P′/P − 17.3(3)·cosh β + 0.2(2)`** (`k = 6`).

* The flat-space helicity term is present with the predicted coefficient `2l²` for both currents and stress tensors, in the
  curved form `−2l²P′/P` valid for all `β` tested (not only `β → 0`): at `β = 1`, `ν = 3` it contributes `−6.6` (`l = 1`) and `−26.3` (`l = 2`),
  with the sign of `P′/P` opposite to that at `β = 0.4`.
* The normalisation polynomials are confirmed (coefficient `1.002` for `l = 2`, where `pol = −64.1` at `β = 0.4`).
* For `l = 2` there is an additional, `ν`-independent-in-form piece `−2 sinh²β P′/P − 17.3 cosh β` relative to the scalar skeleton, absent
  for `l = 1`. It is regular for `β → 0` (a constant, plus `O(β²)`), so it is an AdS-curvature effect and does not affect the flat-space
  helicity structure. It has the form of a shift `α² → α² + 2`, `ᾱ² → ᾱ² + 2` in `m² = αᾱ` and `η_AdS` (which gives `−2 sinh²β P′/P − c·cosh β`),
  but the coefficient of `cosh β` was not derived and its `k`-dependence could not be measured reliably (`k = 4`: too few points). [I]
* `a_0 ≈ 0.1…0.2` is compatible with a small constant for both `l` (and with zero within the systematic error of the fits).

### 6.4 Helicity-flip component

`D_flip := c_2^flip − c_2^{skeleton}` (`k = 6` / `k = 4`):

| `β` | `ν` | `l = 1` | `l = 2` |
|---|---|---|---|
| 0.4 | 0 | −69.4 / −42.6 | −255 / — |
| 0.4 | 1 | −70.0 / −42.6 | −261 / — |
| 0.4 | 3 | −81.6 / −47.4 | −264 / −212 |
| 1.0 | 0 | −61.1 / −32.3 | −163 / — |
| 1.0 | 1 | −72.5 / −38.6 | −174 / −104 |
| 1.0 | 3 | +131.9 / +66.6 | −3 / −21 |
| 2.0 | 0 | −70 / — | −194 / — |
| 2.0 | 3 | +442 / — | +1116 / — |

Unlike the preserved component, the flip correction (i) depends on `k`, (ii) is much larger, (iii) grows like `1/β²` at small `β`
(`F^ν_{2l} ∝ β^{2l}`, so any shift of the argument costs `2l/β²`; the natural parameter is `1/(J+½)²`), and (iv) changes sign where
`F^ν_{2l}(β)` approaches a zero. A fit with the functions of §6.3 (with `F′/F` in place of `P′/P`) does **not** describe it (residuals of
order 100). In flat space the corresponding statement is that the corrections to `d^J_{l,−l}` are controlled by `(J+l)!/(J−l)!` and
by the kernel; no analytic prediction was derived here. [measured only]

### 6.5 Summary of task 4

| statement | status |
|---|---|
| no `O(1/m)` for `⟨JφJφ⟩`, `⟨TφTφ⟩`, both helicity components, in Weyl-odd variables with exact `1/b` | derived (reflection argument), verified pointwise and by fits (`|c_1| ≤ 0.008`, preserved) |
| `O(1/m²)` coefficients as functions of `(ν, β, k)`, preserved and flip | measured (§6.2, §6.4); not derived |
| spin-dependent part of the preserved coefficient: `pol` and the flat-space helicity term `−2l²P′/P` | predicted (`pol`: derived; `hel`: inferred from flat space), confirmed: `a_2 = −2.01`, `−8.00` |
| extra `l = 2` piece `−2 sinh²β P′/P − 17.3 cosh β` (`k = 6`) | measured, form fitted, not derived |
| flip correction | measured, not understood analytically |


## 7. Task 5 — uniformity

Files: `out_t8_uniformity_l{1,2}.txt` (`t8_uniformity.py`). Two directions: `ν = x·m` at fixed `β`, and fixed `J` (small impact parameter,
`β ≈ (2J+1)/m`, down to `β ≈ 0.04`).

### 7.1 `ν` growing with `m` (`ν = x m`, `x = 0.05…0.8`)

Columns: `R_pres` (spinning, preserved component), `R_scal` (the exact *scalar* block with the skeleton dimensions `(l+1, Δ_φ)` against the
scalar formula of §4.1), and `R_pres/(R_scal·corr)` where `corr` is the spin-dependent improvement of §6.3 (exact normalisation
polynomials and `P(η_AdS − 2l²/m²)/P(η_AdS)`; no fitted parameter). Sample (`l = 1`, `k = 6`):

| `β` | `m` | `ν` | `x` | `R_pres` | `R_scal` | `R_pres/R_scal` | `R_pres/(R_scal·corr)` | `zeroP` |
|---|---|---|---|---|---|---|---|---|
| 0.4 | 12.4 | 0.62 / 1.24 / 2.48 / 3.73 | 0.05 / 0.1 / 0.2 / 0.3 | 0.591 / 0.615 / 0.724 / 0.992 | 0.614 / 0.635 / 0.731 / 0.954 | 0.963 / 0.968 / 0.991 / 1.040 | 1.0010 / 0.9985 / 0.987 / 0.960 | 0.99 / 0.96 / 0.87 / 0.68 |
| 0.4 | 42.2 | 2.11 / 4.22 / 8.44 | 0.05 / 0.1 / 0.2 | 0.965 / 1.035 / 1.081 | 0.966 / 1.022 / 1.100 | 0.9992 / 1.0122 / 0.9829 | 1.000004 / 0.99930 / 0.99875 | 0.90 / 0.58 / 0.89 |
| 0.4 | 101.8 | 5.09 / 10.18 | 0.05 / 0.1 | 1.0206 / 1.0349 | 1.0145 / 1.0341 | 1.0060 / 1.0007 | 0.99993 / 1.00020 | 0.34 / 0.97 |
| 1.0 | 39.3 | 1.97 / 3.93 / 7.87 | 0.05 / 0.1 / 0.2 | 1.023 / 0.980 / 1.293 | 1.019 / 0.981 / 1.284 | 1.0039 / 0.9994 / 1.0071 | 0.99957 / 0.99957 / 0.99858 | 0.40 / 0.99 / 0.68 |
| 1.0 | 70.0 | 3.50 / 7.00 | 0.05 / 0.1 | 0.980 / 1.032 | 0.981 / 1.032 | 0.99937 / 0.99971 | 0.999998 / 0.99993 | 0.95 / 1.00 |

Findings [numerical]:
* The leading formula is **not uniform in `ν/m`**: `R_pres − 1` grows like `x²·O(1…10)` (column `(R_scal−1)/x²` of the output) — at `x = 0.2`
  the deviation is 5–30 %, at `x ≥ 0.3` the leading formula is useless (and near a zero of `P_{−½+iν}(cosh β)`, small `zeroP`, the ratio
  is meaningless). This is the same non-uniformity as for scalars (CHMRS: corrections `ν²/m²`, `ν⁴/m⁴`, …).
* **What replaces it**: the exact scalar skeleton times the spin-dependent factor. `R_pres/(R_scal·corr)` stays within `10⁻³` of 1 up to
  `x = 0.2` for `m ≥ 40` (and within `2·10⁻⁴` at `m ≥ 70`, `x ≤ 0.1`), whereas `R_pres` itself is off by up to 30 %. I.e. the `ν`-non-uniformity
  sits entirely in the scalar skeleton; the spin-dependent dressing (helicity-preserving component) is uniform in `ν/m ≲ 0.2` at least.
  At `x ≥ 0.3` the double-precision moments lose accuracy (`supp < 10⁻⁴`, the suppression `e^{−πν}` of the moment against the
  integrand) and near zeros of `P` the comparison is ill-conditioned; there is no evidence either way.
* The flip component has no such simple completion (§6.4).

### 7.2 Small impact parameter (fixed `J`, `β ≈ (2J+1)/m`), including `β ≲ 1/ν`

`m ≈ 20, 40, 80`, `J = 1, 2, 4` (currents) / `J = 2, 3, 5` (stress tensors), `ν = 0, 1, 3, 6`, so `νβ` runs from 0 to 2.7 and `β` from 0.04 to 0.55.
"lead" = leading formula; "impr" = with the exact normalisation polynomials `Π[(J+½)²−(i−½)²]Π[(Δ−3/2)²−(i−½)²]` (derived, §4.2, §6.3);
"impr2" = in addition `η_AdS → η_AdS − 2l²/m²` in the preserved harmonic function (the flat-space helicity term, §6.3); "impr3" = in
addition the factor `Π[(J+½)²−(i−½)²]/(J+½)^{2l}` on the flip harmonic function (flat-space `d^J_{l,−l}`). `l = 1`, `k = 6`, mode `(1,−1)`:

| `m` | `J` | `β` | `ν` | `νβ` | lead | impr | impr2 | impr3 |
|---|---|---|---|---|---|---|---|---|
| 20 | 1 | 0.150 | 0 / 3 | 0 / 0.45 | 0.728 / 0.805 | 0.820 / 0.906 | 0.819 / 0.885 | 0.819 / 0.885 |
| 40 | 1 | 0.075 | 0 / 3 / 6 | 0 / 0.22 / 0.45 | 0.844 / 0.866 / 0.937 | 0.949 / 0.975 / 1.054 | 0.949 / 0.969 / 1.030 | 0.949 / 0.969 / 1.031 |
| 80 | 1 | 0.0375 | 0 / 3 / 6 | 0 / 0.11 / 0.22 | 0.877 / 0.883 / 0.901 | 0.987 / 0.994 / 1.013 | 0.987 / 0.992 / 1.008 | 0.987 / 0.992 / 1.008 |
| 80 | 4 | 0.112 | 0 / 3 / 6 | 0 / 0.34 / 0.67 | 0.975 / 0.981 / 1.002 | 0.987 / 0.994 / 1.014 | 0.987 / 0.992 / 1.008 | 0.987 / 0.992 / 1.009 |

Findings:
* **The leading formula fails at fixed `J` by an `O(1)` factor that does not go away as `m → ∞`**: for `J = 1` the ratio tends to
  `0.877 → J(J+1)/(J+½)² = 8/9 = 0.889` [`l = 1`; for `l = 2`, `J = 2`: `(J−1)J(J+1)(J+2)/(J+½)⁴ = 24/39.06 = 0.614`]. This is the statement that the
  heavy-limit vertex `(A_+A_−)^l` must be replaced by the exact normalisation of the conserved structure, which vanishes at `J < l`. With
  that replacement ("impr") the ratio is `1 − O(1/m²)` *uniformly in `J`*: at `m = 80`, `k = 6`, `ν = 0` it is `0.98697, 0.98696, 0.98691` for
  `J = 1, 2, 4` and `0.99131` at `J = 20` (`m = 102`, `β = 0.4`), i.e. the remaining correction is the `J`-independent `c_2/m²` with
  `c_2 ≈ −80…−90` (`k = 6`), the same as the scalar skeleton at small `β`.
* In the regime `β ≲ 1/ν` (`νβ ≲ 1`) nothing special happens once "impr" is used; the `ν`-dependence of the residual is then removed further
  by "impr2" (e.g. `m = 40`, `J = 1`, `ν = 6`: `1.054 → 1.030`; `m = 80`, `J = 4`, `ν = 6`: `1.014 → 1.008`); the remaining `ν`-dependence is
  `O(ν²/m²)` of the scalar skeleton.
* For `β ≳ 1/ν` with `ν = 6` and `m = 20` (`νβ ≳ 1.5`) the leading formula and all improvements fail (e.g. `J = 4`: lead `−333` near a zero of
  `P`): this is `ν/m = 0.3`, the non-uniformity of §7.1, not a small-`β` effect.
* Flip mode `(1,1)`: "impr3" moves the ratio in the right direction (`0.9855 → 0.9868` at `m = 80`, `J = 1`, `ν = 0`), and at `m = 80` all
  `J` agree to `2·10⁻³`, consistent with `1 − O(1/m²)` once the flat-space factor `(J+l)!/((J−l)!(J+½)^{2l})` is restored.

**Summary of task 5**: the formula of §4.5 is uniform in `β` down to `β ≈ 1/m` *after* replacing `(J+½)^{2l}(Δ−3/2)^{2l}` by the exact
normalisation polynomials (derived); it is not uniform in `ν/m` (`O(ν²/m²)` corrections of order one at `ν/m ≈ 0.3`), but the
spin-dependent part is: exact spinning moment `≈` (exact scalar skeleton moment) × (spin factor of §6.3) to `10⁻³` for `ν/m ≤ 0.2`.
**Stress tensors (`l = 2`, `out_t8_uniformity_l2.txt`, `k = 6`).** (a) `ν = xm`, `β ≈ 0.4`: `R_pres = 0.942, 1.048, 0.998` at `m = 42`, `x = 0.05, 0.1, 0.2`
(and `1.035, 1.031` at `m = 102`, `x = 0.05, 0.1`), while `R_pres/(R_scal·corr) = 0.9898, 0.9892, 0.9820` and `0.99895, 0.99730`. The residual is
now the `ν`-independent `l = 2` piece of §6.3 not included in `corr` (`≈ −17.3 cosh β/m² ≈ −0.010` at `m = 42`); its `ν`-dependence over
`x ≤ 0.2` is `0.8 %` against `10 %` for `R_pres` itself. Same conclusion as for currents. (b) Fixed `J`: "lead" tends to the predicted
`(J−1)J(J+1)(J+2)/(J+½)⁴ = 0.614` (`J = 2`; measured `0.606` at `m = 80`, i.e. `0.614 × 0.987`), and "impr" at `m = 80`, `ν = 0` is
`0.98709, 0.98707, 0.98700` for `J = 2, 3, 5`: uniform in `J` down to `J = l = 2`, `β = 0.06`. "impr2" removes most of the `ν`-dependence at
`νβ ≲ 1` (`m = 80`, `J = 2`, `ν = 6`: `1.0313 → 1.0079`), and for the flip mode "impr3" brings `(2,2)` in line with `(2,−2)`
(`m = 20`, `J = 3`, `ν = 0`: `0.729 → 0.773` vs `0.815`; at `m = 80` all within `4·10⁻⁴`).


## 8. Task 6 — `⟨TTTT⟩`: not reached; what the results imply [D, not verified]

**Not reached.** Tasks 1–5 are complete only with the numerical limitation on `k` stated in §3.2 and without an analytic `1/m²`
term, and the doubly spinning partial wave was not built. What is missing, concretely: (i) the differential basis for two spinning
operators at each vertex (`D_11, D_12, D_22, D_21, H_12`; 14 parity-even structures for `⟨TTO⟩`, 2 conserved, plus parity-odd ones),
i.e. words of order up to 4 at *each* vertex and jets of order 8 of the seeds; (ii) the precision problem of §2.1–2.2 is already
marginal in double precision for one spinning pair at order 4 (errors `10⁻⁷…10⁻¹⁰` pointwise at `Δ ≈ 30…100`), so order 8 requires
the whole post-processing in extended precision (only the radial series exists in `__float128` here); (iii) the heavy limit of the
two conserved `⟨TTO⟩` structures in the differential basis, i.e. the 2×2 map from `λ^{(1,2)}_{TTO}` to `(c^{++}, c^{+−})`; (iv) mean-field data
in the helicity basis for the MFT-normalised form.

**Strongest statement supported by the `⟨TφTφ⟩` results.** The derivation of §4 is local in the two Regge pairs: (a) the shadow
representation of a doubly spinning partial wave follows from the scalar one by acting with the differential operators, with the
same scalar coefficient `b^{(Δ_34)}` (this mechanism is now tested, §2.3); (b) the heavy vertex factorises into a plane wave times
`τ_L = Σ_{h_1,h_2} c^{h_1h_2} ε_{h_1}²(p̄̂) ⊗ ε_{h_2}²(p̂)` (Lemma V of `REGGE_MOMENT_TTTT.md`, verified there numerically for the left vertex only);
(c) in the Fourier–Helgason representation the `x`- and `y`-transforms are performed at fixed boundary point `ℓ` and are independent,
each giving the factor `M^{q q′}_{h h′}(ν;a) e^{iμθ}` *verified here* (`a = k+5` for both pairs). Hence for a kernel
`K^{(q_1q_3)}_{13} ⊗ K^{(q_2q_4)}_{24}`:

`Π̂ → 𝒩 γ_{k+5}(ν)² m^{−2(k+5)} Σ_{h_i} c^{h_1h_2} \bar c^{h_3h_4} M^{q_1q_3}_{h_1h_3}(ν) M^{q_2q_4}_{h_2h_4}(ν) Φ^ν_{μ_L,μ_R}(β)`, `μ_L = l(h_1+h_3)`, `μ_R = l(h_2+h_4)` in the convention of §4.4,

with `Φ^ν_{0,0} = P_{−½+iν}(cosh β)` (the functions of `SPINNING_KERNELS.md`). For the helicity-preserving components `μ_L = μ_R = 0`, and since
`M_{±∓} = R ± g_1 S`:

`Π̂|_{pres} → 𝒩 γ² m^{−2(k+5)} P_{−½+iν}(cosh β) { R_{13}R_{24} (|c^{++}|² + |c^{+−}|²) + g_1² S_{13}S_{24} (|c^{++}|² − |c^{+−}|²) }` (each summed over the two parity images).

So: **yes, for parity-even kernels (`S`-free, i.e. symmetric under `q → −q`) the components `++ → ++` and `+− → +−` have the same harmonic
function `P_{−½+iν}(η_AdS)` and the same kernel factor `R_{13}R_{24}` at leading order**; kernels odd at both pairs see the difference, with the factor
`g_1²S_{13}S_{24}`. This is the structure of `REGGE_MOMENT_TTTT.md` §5 with `κ_+M_𝔖² ↔ R_{13}R_{24}`, `κ_−|M_𝔄|² ↔ g_1²S_{13}S_{24}`, now with the one-pair
factors `R` verified end-to-end (`S`: verified by direct integration only, since parity-even `⟨TφTφ⟩` blocks do not probe it).
The first helicity-dependent correction is `O(1/m²)` (no `O(1/m)`, §6.1); its form was not derived; in the flat corner it must
reduce to the helicity-dependent term of the flat partial waves recorded in `REGGE_MOMENT_TTTT.md` §6.2
(`(1/m²)[2λ²(ν/β)J_1(νβ) + …]`, `λ = h_1 − h_2 ∈ {0, 4}`), and §6.3 shows that a term of exactly this form, `−2l²P′/P` with `λ = l`, is present in
`⟨OφOφ⟩` for `l = 1` and `l = 2` with the predicted coefficient, for `β` up to 1 at least. The couplings: `c^{h_1h_2}` are the coefficients of the helicity tensors in the heavy limit of the conserved `⟨TTO⟩` structures
(normalised as in §4.3, i.e. relative to the differential-basis structure), linear in `λ^{(a)}_{TTO}`; with `T` canonically normalised
`|c|² ∝ λ²` and the sum-rule density is `|c|²/C_T²`; the explicit 2×2 matrix and the MFT normalisation ("`2P/m^{2k}`") are **not done**.


## 9. Task 7 — cross-check against Chang–Landau–Simmons-Duffin

What can and cannot be compared: CLS-D compute the *bulk-point* limit (`ν, m ≫ 1`, `ν < m`) of functionals for four spinning operators
by a saddle point; they do not treat the Regge regime, mixed correlators, or `d = 3` explicitly. The present results are for one
spinning pair, fixed `ν`, any `β`. The overlap is structural.

1. **Scalars, normalisation and factors of 2** [Q + DVE]. CLS-D (`superconvergence_bulkpoint_review`, `scalar_sumrule_CFTtoflatspace`):
   `Ψ_{k,ν;f_k}[Δ,J] → ((2m²−ν²)/(m²−ν²)) 𝒫_J(1−2ν²/m²)/(m²(m²−ν²))^{k/2}`, heavy action defined with the same `q_{Δ,J}` as CHMRS (`heavyaction_scalar_definition`).
   For `ν ≪ m` this is `2𝒫_J(1−2ν²/m²)/m^{2k} → 2J_0(νβ)/m^{2k}`, the `β → 0` limit of CHMRS (2.59) `2P_{−½+iν}(η_AdS)/m^{2k}`, which is the formula verified
   with exact blocks in `REGGE_MOMENT_TTTT.md` §1.3 and, for unequal dimensions, in §4.1 above (same factor 2, no factor `i`). No discrepancy.
2. **Block side.** CLS-D write the block through the Lorentzian shadow representation with contracted indices and an `r`-tensor
   (`block_shadowrep0`, `shadowcoeff_definition`, `rtensor_naturalformula`); the spinning shadow coefficient "will involve mixing of different
   tensor structures" and is replaced by a large-`Δ` saddle. Here the continuous-spin (`∫Dn`) representation of CHMRS is used and
   the spinning coefficient is obtained exactly: for the unique conserved parity-even `⟨OφO′⟩` structure it is the scalar `b^{(1−Δ_φ)}`
   in the differential basis (§4.2). Their large-`Δ` conservation statement (§4.1, `conservation_largeDelta`: `v^μ_{1,23}⟨O_{1;μ…}…⟩ = 0`, fixed `J`) is the
   `J/Δ → 0` case of §4.3 (tensor transverse to `p̄`, helicity `±l`), extended here to `J/Δ` fixed with the coefficient.
3. **Polarisations and helicity labels** [D]. In CLS-D the light-transformed operators carry the polarisations `x_{1^+0}`, `x_{30}`, which
   map in flat space to the *same* null vector `s̄` for both shock gravitons (`sumrule_bulk_polarizationmapall`: `e_1 = e_3 ∝ s̄`). In `D = 4` a null
   transverse polarisation is a helicity eigenstate, so the CLS-D functionals probe the **helicity-preserving** component of the shock
   pair at large `ν`. The `d = 3` reduction of §3.3 reproduces this: with the weights `(−1)^{q_3}c_{(q_1+q_3)/2}/c_0`, the heavy action is
   `A_pres(ν)P_{−½+iν}(cosh β) + A_flip(ν)F^ν_{2l}(β)` with (`clsd_kernel.py`)
   `A_flip/A_pres = −3/(2ν²) + O(ν⁻³)` (`l = 1, a = 5`), `−15/(2ν²)` (`l = 1, a = 7`), `−15/(8ν⁴)` (`l = 2, a = 7`), `+105/(8ν⁴)` (`l = 2, a = 9`):
   the flip admixture is `O(ν^{−2l})` and disappears in the bulk-point regime, but **it is present at `ν = O(1)`**, i.e. in the curved
   regime the CLS-D-type functional is not purely helicity preserving. The odd structure (`j′ = 1`) gives zero on parity-even
   `⟨OφOφ⟩` partial waves, consistently with its signature `(−1)^{j′}` under `1 ↔ 3`.
4. **Flat limit of the harmonic functions.** `P_{−½+iν}(cosh β) → J_0(νβ)` and `F^ν_{2l}(β) → (−i)^{2l}J_{2l}(νβ) = (−1)^l J_{2l}(νβ)`, with `β ≈ (2J+1)/m`:
   these are the small-angle limits of the Wigner functions `d^J_{l,l}` and `d^J_{l,−l}` of the flat `Oφ → Oφ` partial waves (helicity preserved /
   flipped), as in CLPS (2.12)–(2.19) for gravitons. The phase `(−1)^l` (no factor `i` for even `2l`) and the absence of any other
   `β`-dependence at leading order agree with the flat-space structure; the `i`'s in `M_flip` (`(½+iν)_{2l}`) combine with
   `Γ(½+iν)/Γ(½+iν+2l)` in `F^ν_{2l}` into the real function `P^{2l}_{−½+iν}`.
5. **Not compared / open.** (a) The overall `(k,ν)`-dependent constants of CLS-D's spinning functionals (they fix none: "we do not include any
   coefficient in front of the integral"), their matrices `M^±_{4γ}`, `M^±_{4g}` (general `D`, four spinning operators; several entries vanish
   or degenerate at `D = 4`, and the kernels `j′ ≥ 2` do not exist in `d = 3`) and the `e^{iπ/4}` of their saddle (a property of their
   unnormalised `Ψ^+`, absorbed in `A_{k,ν}` for scalars). (b) The identification of their `k` with the moment index `k_r` used here
   (expected: equal, since both count powers `1/m^{2k}` of the heavy density), and the role of the subtraction factors beyond leading
   Regge order. (c) Whether the spinning analogue of CHMRS (2.60) holds for all `ν < m`: CLS-D "expect the same to be true in the spinning
   case"; §7 shows that the *Regge-moment* formula has `ν²/m²` corrections of order one already at `ν/m ≈ 0.3`, as for scalars, so this
   needs the functional `C`-type combination, not a single moment.
   **Discrepancies found: none**; one refinement (item 3) and the caveats (a)–(c).


## 10. Status list

| # | statement | mark |
|---|---|---|
| 1 | Radial Casimir recursion for general `(a,b)`; scalar blocks; mean-field decomposition of mixed scalars | [DVE] |
| 2 | Conserved parity-even `⟨OφO′⟩` structures (`l = 1, 2`) in differential and standard bases; `L_0`, `ℓ_i` | [D], conservation verified on the 4-pt function |
| 3 | Spinning partial waves `W_diff` from exact blocks; OPE normalisation `λ²(−½)^J W_diff/L_0²` | [DVP] (OPE limit, conservation) |
| 4 | Lorentzian shadow coefficient for unequal dimensions `b^{(Δ_34)}` | [DVE] (exact conservation identity, direct quadrature, scalar end-to-end with `R_0` control) |
| 5 | Spinning shadow coefficient = `1/b^{(1−Δ_φ)}` in the differential basis | [DVP] (exact identity + direct 4-d quadrature, 0.1–0.3 %) |
| 6 | Exact identity (★) between tensor Regge moments and the shadow integrand | [D] (kinematic, as CHMRS (C.4)–(C.6)); verified through 7 |
| 7 | Leading heavy formula §4.5 incl. absolute normalisation, both helicity components, `l = 1, 2`, `k = 4, 6` (and `k = 2`, `m ≲ 40`), `ν = 0,1,3`, `β = 0.4,1,2` | **[DVE]**, ratio → 1 as `1 + c_2/m²` |
| 8 | Table of `M^{q_1q_3}` (all `(2l+1)²` modes; overdetermination) | [DVE] for the parity-even (`R`) part; `S` part [DVP] (direct `H²` integration) |
| 9 | No `O(1/m)` correction (Weyl-odd variables, exact `1/b`) | [D] + verified pointwise and end-to-end (`|c_1| ≤ 0.008` in 11/12 preserved cases, all `≤ 0.027`; `≲ 0.2`, flip `l = 2`) |
| 10 | `O(1/m²)` coefficients, preserved and flip, vs `(ν, β, k)` | numerical evidence only |
| 11 | Spin-dependent `1/m²` term, preserved: `pol` + `−2l²P′/P` | `pol` [D, confirmed]; `−2l²P′/P` [I from flat space, confirmed: `−2.01`, `−8.00`] |
| 12 | Extra `l = 2` piece `−2 sinh²β P′/P − 17.3 cosh β` (`k = 6`) | [I] (form fitted; not derived) |
| 13 | Flip `1/m²` correction | measured, not understood |
| 14 | Uniformity at small `β` after restoring the normalisation polynomials | [D] + numerical evidence (`J ≥ 1`, `β ≥ 0.04`) |
| 15 | Non-uniformity in `ν/m` is that of the scalar skeleton (spin factor uniform to `ν/m ≤ 0.2`) | numerical evidence |
| 16 | `d = 3` reduction of the CLS-D polarisation prescription; flip admixture `O(ν^{−2l})` | [D] (`clsd_kernel.py`), not compared with a CLS-D number |
| 17 | MFT-normalised form "`2P/m^{2k}`" for spinning operators | [ND] |
| 18 | Mean-field decomposition of `⟨OφOφ⟩` | [ND] |
| 19 | Channel `TT → φφ` | [ND] |
| 20 | Parity-odd structures; `S`-terms end-to-end; `(|c^{++}|² − |c^{+−}|²)` terms | [ND] (`S` only [DVP]) |
| 21 | `⟨TTTT⟩` helicity-preserving components | [D] at leading order by factorisation (§8), [ND] end-to-end |
| 22 | `T` normalisation via Ward identity / `C_T` | [Q], enters only as `λ²/C_T`; not tested |


## 11. What remains; exact blocking steps

1. **`k = 2, 3` numerically.** The radial series is a poor representation in the corner `ρ̄/ρ → 0` (`z̄ → −∞` at fixed `z`), which weighs
   `u^{k−1/2}` in tensor moments. *Blocking step:* an expansion of the `d = 3` block around `z̄ = ∞` at fixed `z` (the crossed lightcone
   expansion with its logarithms), or `N ∼ 10⁵` radial levels (the cost grows like `N³`), or an extrapolation of the level sums in `n`.
2. **Analytic `1/m²` term.** Mechanical with the present tools: expand the exact integrand `𝓛·F^S·ℱ` of §4.2 to second order at fixed
   `(X,Y)` in the variables `A_±` (the first order is antisymmetric under `1 ↔ 3`, §6.1) and apply `laplace.M_coeff` with the extra
   polynomial factors of `x, y` (derivatives in `p`, `p̄`). Not done; §6.3 gives the measured coefficients and an empirical formula for
   the spin-dependent part of the helicity-preserving component.
3. **Parity-odd structures.** The conserved parity-odd `⟨OφO′⟩` structure needs the second-order operator
   `D̃_1 = ε(Z_1, P_1, ∂_{P_1}, P_2, ∂_{P_2})` (jets of order `2l+2`). It is needed for: the `S`-parts of `M` end-to-end (now [DVP]), the
   mean-field decomposition of `⟨OφOφ⟩`, and parity-violating theories.
4. **Mean-field normalisation** (`2P/m^{2k}` form): needs `p^{MFT}` of `[Oφ]_{n,J}` in the conserved basis (both parities), e.g. from CHL
   §3.2, and its large-`Δ` limit; equivalently the unitarity normalisation of `𝔖^{cons}` (the factor that replaces `e^{lβ}`).
5. **`TT → φφ` channel** (`λ_{TTO}λ_{φφO}`, Regge pairs `(T,φ)`): needs the functional for non-identical Regge pairs (different phases in
   the double discontinuity) and the operators at both points of one vertex; the kernels are then `Φ^ν_{l,±l}` rather than `P`.
6. **`⟨TTTT⟩`**: §8, items (i)–(iv).
7. **CLS-D functionals themselves**: their normalisation, the relation `k ↔ k_r`, and the Regge-moment expansion beyond the leading
   term (the analogue of CHMRS (C.13)–(C.16)) in `d = 3`; §3.3 and §9 give only the leading-Regge polarisation structure.
8. **Uniformity** in `ν ∼ m` and in the corner `β ≲ 1/m`: §7 identifies what replaces the leading formula (the exact scalar skeleton,
   resp. the exact normalisation polynomials and `η_eff`), on numerical evidence only; a derivation would follow from item 2 and from
   the bulk-point saddle of CLS-D applied to one spinning pair.
9. Positivity statements (signs of `ν`-moments after smearing, the light side, convergence/swapping for spinning functionals) were
   not addressed here.

## 12. Code index (`ads4/code/e2e/`; `sh run_all.sh` regenerates all outputs)

| file | content | output |
|---|---|---|
| `casimir_terms.py`, `terms.json` | radial Casimir recursion for general `(a,b)` (symbolic derivation) | `out_casimir_terms.txt` |
| `cblock.c`, `blocks.py` | recursion and derivative sums, `double`/`__float128`; mpmath reference | — |
| `t1_scalar_blocks.py`, `t1b_mft_scalar.py` | scalar block validation; mean-field decomposition | `out_t1_*.txt`, `out_t1b_*.txt` |
| `embed.py`, `structures.py` | embedding algebra, differential operators, conserved structures, standard basis | `out_structures.txt` |
| `spinblock.py`, `gen_phi_l*.py` | spinning partial waves (generated code), `(P,Q)` jets | — |
| `t2_spinning_blocks.py`, `t2b_shadow_consistency.py`, `t5_shadow_direct.py` | conservation, OPE limit, precision; shadow coefficients; direct shadow integral | `out_t2_*.txt`, `out_t2b_*.txt`, `out_t5_*.txt` |
| `vertex.py` | exact shadow integrand, heavy limit, pointwise `O(1/m)` test | `out_vertex.txt` |
| `laplace.py`, `mforms.py`, `t6_laplace_check.py` | tensor Fourier–Laplace transform (`M` tables), closed forms, direct integration | `out_laplace_table.txt`, `out_mforms.txt`, `out_t6_*.txt` |
| `formulas.py`, `predict.py`, `cases.py` | `b^{(Δ_34)}`, `γ_a`, `F^ν_μ`; the leading formula; `(Δ,J)` at fixed `β` | — |
| `moments.py` | grids, scalar and tensor Regge moments, truncation mask, tail estimate | — |
| `t3_scalar_mixed.py` | scalar skeleton, unequal dimensions, end-to-end | `out_t3_scalar_mixed.txt` |
| `t4_spinning_e2e.py`, `t4b_convergence.py` | spinning end-to-end runs (cache in `data_modes/`), numerical-error study | `out_t4_*.txt`, `out_t4b_*.txt` |
| `t7_corrections.py`, `t7b_scalar_skeleton.py`, `t7c_spin_dependence.py` | helicity components, `1/m`, `1/m²`; scalar reference; spin dependence | `out_t7*.txt` |
| `t8_uniformity.py` | `ν = xm`, small `β` | `out_t8_uniformity_l*.txt` |
| `t9_tables.py`, `t9b_summary.py` | measured and predicted values; ratio tables and controls (§5) | `out_t9_tables.txt`, `out_t9b_summary.txt` |
| `clsd_kernel.py` | `d = 3` Regge reduction of the CLS-D polarisation prescription | `out_clsd_kernel.txt` |

Run time on 4 cores: block evaluations (`t4`, all four runs) about 2.5 h (cached in `data_modes/`, git-ignored, as are the `*.pkl` caches);
`t4b` about 30 min; everything else about 1 h. All caches are regenerated by `run_all.sh` if absent.
