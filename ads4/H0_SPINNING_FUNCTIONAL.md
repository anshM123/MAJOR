# H0: the spinning (`⟨TTTT⟩`) functional — route decision and first milestone

Status: 18 Nov 2026, derivation/numerics agent. **Unrefereed. Not committed.** Targets hypothesis H0 of `RESULT_D.md` (with H2, H4, H14,
H9b) and finding F1 of `REFEREE_FOUNDATION.md`. Code: `ads4/code/spinfunc/` (`run_all.sh`, single core, ≈ 10 min; outputs in
`spinfunc/cache/`, git-ignored). It imports `omegau/`, `exactskel/` (e1/e7 caches) and `objid/` **read-only**. Nothing else was modified.

Source read for this round: Chang–Landau–Simmons-Duffin, arXiv:2311.04271 (**CLS-D**), LaTeX source (`draft.tex`, §2–§4.6, App. B),
downloaded to the session scratchpad. Equation labels below are CLS-D's LaTeX labels.

Marks: **[Q]** quoted; **[D]** derived; **[DV]** derived and checked (exact algebra, sympy); **[D-sketch]**; **[NE]** numerical evidence;
**[M]** model; **[I]** assumed/inferred.

## 0. Summary

1. **Route decision: A, in a concrete form ("A′").** The spinning C-type functional the project has been calling "not constructed" is,
   up to one step, already in the literature:
   * CLS-D define light-ray-commutator sum rules `Ψ^{(a),(c)}_{k,ν,j′}` for external spinning operators (stress tensors included), with
     **meromorphic** subtraction factors `f^±_{k;J_e}(u′,v′)` (`subtraction_factor_evenjprime`, `…oddjprime`). Their scalar version
     *is* CHMRS's `C̃_{k,ν}` (CLS-D footnote to `scalarfunctional_nokernelexpression`) [Q].
   * They prove (saddle point, general `d`, `ν, m ≫ 1` at fixed `ν/m < 1`) that the heavy action of these functionals equals the flat
     graviton sum rule (shock amplitudes; the lowest-subtracted fixed-`u` sum rules of CLPS) [Q].
   * **The project's `Ψ_k` are exactly these CLS-D functionals with the subtraction factor replaced by the weights `A_k`**
     (`SUPERCONVERGENT_UPLIFT.md` item 1); `A₃`, `A₇` coincide with CLS-D's, `A₄, A₅, A₆` do not. The replacement was made to get pure
     harmonics in the Regge limit (CLS-D's factors carry `η^j`, `j = 1,0,3,2,1,0` for `k = 2…7`). It is what breaks meromorphy (H9b) and
     the flat limit (F1).
   * **Route A′:** `Ĉ_k :=` CLS-D's `Ψ_{k,ν,j′; f^±_{k;2}}` (same light-ray operators, same dual/celestial structures as the project's
     `Ψ_k`), followed by the CHMRS `η`-inversion in `ν` (Legendre recursion), for **all** `k = 2…7`, including `Ĉ₂` (CLS-D's unsubtracted
     `f⁺_{2;2} = v′−u′`). Route B is dominated (§1).
2. **Milestone M1 executed: the graviton bulk-point dictionary of both objects, exactly in `x = ν/m`** (§3). Using CLS-D's factorisation
   (the subtraction factor enters the bulk-point action only through its value at the "scattering crystal"):
   * the crystal is the CHMRS Regge-frame point `η = 1`, `r = r* = (1−s)/(1+s)`, `s = √(1−x²)`: all 9 printed CLS-D saddle values
     `f^±_{k;J_e}|crystal` are reproduced **[DV]**;
   * cross-check: the same method reproduces `OBJECT_IDENTITY.md`'s independently derived (and refereed) scalar density ratio
     `(1−r*²)^{k+2}/(1+r*²)` for `k = 2, 4, 6` **[DV]**;
   * **graviton result [D given CLS-D]:** relative to B′ (= CLS-D's flat sum rule = route A′ at leading order), the bulk-point density of
     the project's `Ψ_k` is `R_k(x)`:

     | `k` | `R_k` (in `r = r*`) | in `s = √(1−x²)` | `x = 0.5 / 0.7 / 0.9 / 0.99` |
     |---|---|---|---|
     | 3, 7 | 1 | 1 | 1 |
     | 4 | `(1−r²)²/(1+r²)` | `8s²/((1+s)²(1+s²))` | 0.985 / 0.920 / 0.620 / 0.120 |
     | 5 | `(1+6r²+r⁴)/(1−r²)²` | `(1+s⁴)/(2s²)` | 1.042 / 1.235 / 2.73 / 25.1 |
     | 6 | `(1−r²)⁴/(1+r²)` | `128s⁴/((1+s)⁶(1+s²))` | 0.975 / 0.869 / 0.443 / 0.023 |

     `R_k − 1 = O(x⁴)` (`R₄ = 1 − 3x⁴/16 − …`, `R₅ = 1 + x⁴/2 + …`, `R₆ = 1 − 5x⁴/16 − …`).
   * **The scalar-analogue model of `OBJECT_IDENTITY.md` was wrong for gravitons**, qualitatively at `k = 5` (it has `Ψ₅`'s density → 0 as
     `x → 1`; the graviton `Ψ₅` density diverges like `1/(2(1−x²))`), and by four powers of `(1−r²)` at `k = 4, 6`.
3. **Consequences [NE + M]** (§4):
   * **F1 sharpened: the flat limit of the functional as currently defined is not positive.** With the graviton densities, the D-safe flat
     action is negative for `+−`, `m/M < 1.0003`, `J = 16–24` (min `F/N_s = −0.0162` at `m = M`, `J = 20`; quadrature-independent, 600–6000
     nodes). The cause is the `k = 5` term. Under the twist gap (`m_min(J)/Λ = √(1+2J/(Δ_g−1))`) these states exist only for
     `Δ_g ≳ 6·10⁴`, so the `Δ_g = 100–1000` numbers are not affected at leading order, but the current functional is **not** the uplift of
     a positive flat functional, and the NF02 certificate does not apply to it. **Route A′'s flat limit is NF02 exactly**, so the
     certificate applies [Q + D].
   * **Near gap, `Δ_g = 100 / 300 / 1000`** (bulk-point-factor model on the exact scalar skeleton, all 1745 / 1223 / 1042 records):
     `n = 1` for both objects. Min `F/N_s`: route A′ `0.2267 / 0.2057 / 0.2022` (= `RESULT_D.md`, by construction); current `Ψ_k` with
     graviton densities `0.2281 / 0.2077 / 0.2042`. Largest object shift `0.035 / 0.052 / 0.121 N_s`.
   * **Where the near-gap action lives** (route A′; worst-sign and coherent tolerances, §4.3): the region `x > 0.7` carries a median
     `≤ 0.5 %` of `Σ|density|`; a coherent relative error of the graviton leading density on `x > 0.7` of `±22 %` (`Δ_g = 1000`, worst
     of three sign patterns) to `±180 %` (`Δ_g = 100`) does not move `n = 1`. The load-bearing window is **`x ∈ [0.3, 0.7]`**: there a
     coherent density error of `≈ 4.5–8 %` (`x > 0.3`) or `≈ 8–30 %` (`x > 0.5`) would break `n = 1`. The graviton-vs-scalar remainder (H4) may differ coherently by a
     factor `≈ (−1.2, 2.6)` on `x > 0.3` at `Δ_g = 100`, and by factors `> 10` at `Δ_g ≥ 300`.
   * **`η`-inversion cost** (large-`ν` symbol estimate [D-sketch + NE]): `≤ 0.087 / 0.0088 / 0.0012 N_s` (`∝ 1/Λ²`); `n = 1` unchanged.
4. **Impact on the hypotheses** (§6): H0 goes from "object not constructed" [I] to "object identified in the literature (general `d`);
   `d = 3` specialisation, `η`-inversion analyticity, convergence and light side to be shown" [D-sketch]. H2 at leading order in `1/m` and
   fixed `x < 1` becomes [Q] for route A′ (conditional on the `d = 3` specialisation, M2). For the current `Ψ_k` the B′ assumption is
   **false** at `x ~ 1` for `k = 4, 5, 6` (densities `R_k`·B′). H4 is unchanged ([I]) but localised: `x ∈ [0.3, 0.7]`, and at
   `O(1/m²)` it is a computable fluctuation correction around the CLS-D crystal (M6).
5. **The check suggested for this round ("verify `Ĉ₄` against the existing `Ψ₄`/`Π̂` data at `x ≲ 0.1`") cannot discriminate.** All
   candidate objects agree to `O(x⁴)`: `|R_k − 1| ≤ 5·10⁻⁵` at `x = 0.1`, against `10⁻³` data precision. The leading-Regge identity of `Ĉ_k`
   with `Ψ_k` holds by construction; it was not tested numerically here.

## 1. Route decision

| | Route A (spinning `C`-type `Ĉ_k`) | Route B (keep `Ψ_k/N_k`) |
|---|---|---|
| object | CLS-D `Ψ_{k,ν,j′}` with CLS-D's meromorphic `f^±_{k;2}`, then `η`-inversion | CLS-D `Ψ_{k,ν,j′}` with the project's `A_k` (non-meromorphic `r^{k−3}` for even `k`) |
| literature | construction, validity (commutator vanishing; Regge convergence `k > J₀` for even `j′`, `k − 1 > J₀` for odd `j′`) and bulk-point dictionary to flat graviton sum rules: CLS-D [Q]; `η`-inversion: CHMRS App. superconvergence (scalar) [Q] | same functional class; nothing on non-meromorphic weights |
| flat limit | **NF02 exactly** at leading order (CLS-D) → certified (`NOFORWARD_FUNCTIONAL.md`) | B′ × `R_k` (§3): **not** NF02, **not positive** (§4.1); needs new weights and a new certificate. And for even `k` the flat image is not a meromorphic dispersive sum rule (`R₄`, `R₆` contain `√(1−x²)`), so it is not a standard flat functional at all |
| light side (`k(0)`, `A`, `B_AdS`, graviton pole, contacts, `[TT]`) | CLS-D/CHMRS-type: meromorphic weight ⇒ `dDisc` form, finitely many subtracted double-trace families, polynomial contacts ⇒ **(P)/H9b analysable** [D-sketch]; flat-limit light side by the same dictionary (CLS-D: "it follows that the contributions of light states must match") | even `k`: (P) **unanalysable** (no `dDisc` rewriting; H9b stays [I] forever); the flat image of the light side is not a standard EFT sum rule |
| exact near-gap value (H2/H4 at `x ~ 1`, finite `m`) | not available; but: (i) leading order in `1/m` is [Q]; (ii) `O(1/m²)` is a fluctuation integral around the crystal with algebraic slow factors (M6); (iii) Mellin route possible in principle (meromorphic weights) | not available even for scalars (non-meromorphic weight: no Mellin form, `OBJECT_IDENTITY.md` §2.4) |
| cost to a statement at the level of `RESULT_D.md` | M2–M5 below: ≈ 5–8 weeks (§7) | re-optimise flat weights, certify, recompute light side with non-standard flat sum rules; (P) still open; ≈ 4–6 weeks for a weaker endpoint |

**Decision: route A′.**
* **Why A is feasible now.** The two hard components were thought to be missing: the spinning functional and its flat limit. Both are in
  CLS-D. What is left is the `d = 3` specialisation (CLS-D work in general `d`, with `D ≥ 8` for their amplitude basis), the
  `η`-inversion (an explicit difference equation in `ν`), convergence at the single-ray endpoints, and the light side.
* **Why B is dominated.** Its only advantage, a pure-harmonic Regge kernel without inversion, was bought by changing the subtraction
  factor. That change is what destroys:
  * meromorphy (H9b);
  * the identity of the flat limit with NF02 (F1);
  * and, as found here, flat positivity itself.

  Route B would have to repair all three. Route A′ keeps the CLS-D factor and moves the Regge-side harmonic purity into the
  `ν`-inversion, at an `O(1/Λ²)` cost (§4.4).
* **What neither route provides cheaply:** exact finite-`m` graviton moments at `x ∈ [0.3, 0.7]` (H4). Position space needs 50–90 digits
  there (`SPINNING_MELLIN.md` §1.1). The realistic next step is the `1/m` expansion around the CLS-D crystal (M6), not exact Mellin.

## 2. Route A′: definition and properties

**Definition.** For `k = 2…7`, `j′ ≡ k (mod 2)`, with `(a)` the dual structure at `(2,4)` realising the project's flip-free kernel (`K*` or
`𝔖_ff` for even `k`, `𝔄` for odd `k`) and `(c)` the celestial structure at `(1,3)` of the project's `Ψ_k`:

`Ψ^{CLS}_{k,ν} := Ψ^{(a),(c)}_{k,ν,j′}[f^±_{k;2}(u′,v′) G]`  (CLS-D `spinningfunctional_definition` with `subtraction_factor_…jprime`, `J_e = 2`),

with `f⁺_{2;2} = v′−u′`, `f⁻_{3;2} = −1`, `f⁺_{4;2} = ((v′−u′)³ − (v′−u′))/(4u′v′)`, `f⁻_{5;2} = −((u′−v′)²−1)/(4u′v′)`,
`f⁺_{6;2} = (v′−u′)/(u′v′)`, `f⁻_{7;2} = −1/(u′v′)`.

* **Leading Regge kernel [DV, `f1` item (3)].** `f_{k;2} = c_k η^j r^{k−3}(1 + O(r²))`, with `(c_k, j) = (½,1), (−1,0), (8,3), (−16,2),
  (128,1), (−256,0)`. The project's `A_k` are `r^{k−3}(1 + O(r²))` with no `η`.
* **`η`-inversion [D-sketch].** In `d = 3`, `ηP_ν = ½(1 + 1/(2iν))P_{ν−i} + ½(1 − 1/(2iν))P_{ν+i}` (`P_ν ≡ P_{−½+iν}(η)`). For a
  `ν`-smeared functional `∫dν φ(ν) X_ν` (`X`: pure-harmonic leading kernel, the object the D-safe weights `φ_k` refer to), the same leading
  kernel is obtained from `∫dν φ̃(ν) Ψ^{CLS}_{k,ν}` with `φ̃` solving
  `φ(ν) = ½φ̃(ν+i)(1 + 1/(2i(ν+i))) + ½φ̃(ν−i)(1 − 1/(2i(ν−i)))` (iterated `j` times). Equivalently, `φ̃` is the Mehler–Fock transform of
  `K_φ(η)/η^j`, which exists because `η ≥ 1`. At large `ν`, `φ̃ = φ + j[φ″/2 − φ′/(2ν) + φ/(2ν²)] + O(ν⁻⁴)`. This needs analyticity of
  `Ψ^{CLS}_{k,ν}` in `|Im ν| ≤ j` (as for CHMRS `C₂`) [I].
* **Bulk-point limit.** At the crystal `η = 1`, so the inversion is the identity up to `O(1/ν²)`. The bulk-point action is CLS-D's flat
  sum rule (B′) [Q], with an inversion correction estimated in §4.4.
* **Heavy Regge action.** Same leading kernel as the project's `Ψ_k` (normalisation `N̂_k(ν) = c_k N_k(ν)`), hence the same
  `(2P/m^{2k})(|c_s|² ± |c_o|²)` form and the same [DVE] scope of `SPINNING_TTTT_EXT.md` [D-sketch]. Subleading-in-`r` structure differs
  (affects `Π̂_{k+2}`-level terms only).
* **Validity / convergence [Q, conditional].** Commutator vanishing for any meromorphic `f` (CLS-D §3.3). Regge convergence for
  `k > J₀` (even `j′`) and `k − 1 > J₀` (odd `j′`), CLS-D §3.5.2. Endpoint convergence for `3 ≤ J₀ < k` (the `|ĝ₃|²` question, H15) is not
  addressed by CLS-D [I] (M4).
* **Light side.** `Ĉ₂` uses the unsubtracted `f⁺_{2;2}`: zero on all double traces; the graviton pole arises from the marginal `J₀ = 2`, as
  for CHMRS `C₂` [D-sketch]. For `k ≥ 3` the weights are meromorphic, so the CHMRS/CLS-D double-trace pattern applies and (P) is
  analysable [D-sketch] (M5).

## 3. Milestone M1: the graviton bulk-point dictionary (`f1_bulkpoint.py`)

**The factorisation used [Q].** CLS-D §3.3, "Saddle point analysis", and §3.5.1: in the bulk-point limit the quickly varying part of the
integrand is identical for scalar and spinning operators. The integral localises on the scattering crystal (`scattering_crystal`), and
"since the subtraction factors do not have any `ν` and `m` dependence", `f(u′,v′)` enters only through its crystal value. **Hence two
functionals with the same `(a)`, `(c)` and weights `f₁`, `f₂` have bulk-point actions in the ratio `f₁/f₂|crystal`** [D from Q]. This
holds for any weight that is smooth and nonzero near the crystal. `sgn(v′−u′)r^{k−3}` qualifies: `sgn = +1` on the moment domain,
`SUPERCONVERGENT_UPLIFT.md` `s3`.

**Normalisation.** Both `Ψ_k/N_k` and `Ĉ_k/N̂_k` are normalised by their leading Regge kernel. The Regge and bulk-point regimes overlap at
`x → 0` (CHMRS, scalar case; CLS-D expect the same for spin). So the density of `Ψ_k` relative to `Ĉ_k` is
`R_k(x) = [A_k/f_{k;2}](η=1, r*(x)) / [A_k/f_{k;2}](η=1, r → 0)`.

**Checks [DV] (`cache/f1_bulkpoint.txt`):**
1. **Crystal ↔ frame.** With the CHMRS frame `u′ = (1−2ηr+r²)²/(16r²)`, `v′ = (1+2ηr+r²)²/(16r²)` at `η = 1`, `r = r*(x)` and CLS-D's
   orientation (`u′ ↔ v′`), every printed crystal value `f^±_{k;J_e}|crystal` is reproduced to `≤ 4·10⁻¹⁴` relative: `J_e = 0`, `k = 2, 4, 6`;
   `J_e = 2`, `k = 2…7`. E.g. `u′v′ = s⁴/x⁸`, `|v′−u′| = (2−x²)/x²`, `(u′−v′)² − 1 = 4s²/x⁴`.
2. **Scalar cross-check.** `A = r^{k+1}` against the scalar `f^{scalar}_k` (= CHMRS `C̃_k`) gives exactly `(1−r*²)^{k+2}/(1+r*²)`
   (`k = 2, 4, 6`). This is `OBJECT_IDENTITY.md` §2.5's `Ψ/C` ratio, derived there from CHMRS's large-`ν` inversion kernels and confirmed by
   the referee (`rf4`). Two independent derivations agree.
3. **Graviton `R_k`:** the table of §0 item 2. The difference from the scalar analogue is a factor `(u′v′)^{±1}` between the `J_e = 0` and
   `J_e = 2` reference weights. At the crystal that factor is `∝ (1−r²)⁴`, which removes four powers of `(1−r²)` at even `k`. At `k = 5` the
   project's `A₅ = (u′+v′)/(32u′v′)` (pure `r²` at leading Regge order) replaces CLS-D's `η²`-carrying `f⁻_{5;2}`. At the crystal
   `(u′−v′)² − 1 → 0` as `x → 1`, hence the divergence.

**Scope [I].** CLS-D's derivation is for general `d` at `ν, m ≫ 1`, `ν < m`. The `d = 3` case (`S^{d−3} = S⁰` locus, `SO(d−2)` trivial,
`j′ ∈ {0,1}`) is not separately checked here (M2). Continuity to finite `ν` is CLS-D's stated expectation, proven for scalars by CHMRS.

## 4. Impact

### 4.1 Flat limit (`f2_flat.py`, `f2b_check.py`) [NE]

Min flat `F/N_s` over `J ≤ 80` (NF02 weights; object `A` = NF02 itself, B′ identity `≤ 10⁻¹³`):

| object | `++` | `+−` |
|---|---|---|
| A′ (`R = 1`) | +0.2000 (`m = M`, `J = 4`) | +0.2000 (`m = M`, `J = 5`) |
| current `Ψ_k`, graviton densities (`PsiG`; `C⁺₂` taken as B′) | +0.2540 | **−0.0162** (`m = M`, `J = 20`) |
| same, `C⁺₂` as `Π̂₂` (`PsiG2`) | +0.2303 | −0.0158 |
| `OBJECT_IDENTITY` scalar analogue (`PsiS`) | +0.2715 | +0.1782 |

* **Negative states.** `+−`, `J = 16–24`, `m/M < 1.000033 … 1.000284` (bisection per `J`). Identical at 600, 2000 and 6000 quadrature
  nodes and for `J ≤ 120`.
* **Mechanism.** At `m = M`, `J = 20` the `k = 5` part goes from `−0.020 N_s` (B′) to `−0.288 N_s`.
* **AdS reach.** Under the twist gap a spin-`J` state has `m ≥ m_min(J) = √((Δ_g−2)(Δ_g+2J−1))`. The window is therefore reached only for
  `Δ_g > 6.3·10⁴` (`J = 18`), `7.9·10⁴` (`J = 20`), `1.4·10⁵` (`J = 22`), `4.7·10⁵` (`J = 16, 24`).

### 4.2 Near gap (`f3_ads.py`) [NE + M]

* **Model** (as `objid/o2`): object moment = `R_X(k,x)` × the exact-skeleton model moment, applied to density, curvature and exact
  remainder. `A′` reproduces e7 (node sum vs `F_ex`: `≤ 2·10⁻⁹ N_s`).
* **Scope:** all e1 states (both skeletons, `J ≤ 10`, high-`J` samples).

| `Δ_g` | object | `n` (`δ = 0` / `0.03`) | min `F/N_s` [state] | max `|ΔF|/N_s` vs A′ |
|---|---|---|---|---|
| 100 | A′ / PsiG / PsiG2 | 1 / 1 / 1 | 0.2267 / 0.2281 / 0.2381 [`+−`, 1.05Λ, `J = 5`] | — / 0.035 / 0.058 |
| 300 | A′ / PsiG / PsiG2 | 1 / 1 / 1 | 0.2057 / 0.2077 / 0.2171 [`+−`, 1.02Λ, `J = 5`] | — / 0.052 / 0.067 |
| 1000 | A′ / PsiG / PsiG2 | 1 / 1 / 1 | 0.2022 / 0.2042 / 0.2077 [`+−`, 1.01Λ, `J = 5, 4`] | — / 0.121 / 0.123 |

`√(t₄²+α_T²)` is unchanged (`n = 1`) in this model for both objects. But only for A′ is the flat limit the certified NF02 and the light side
(`A`, `B_AdS`) the one computed (§1).

### 4.3 Where H2/H4 are load-bearing (`f3_ads.py`, `f5_coherent.py`; route A′) [NE + M]

* **Share in `x`.** The share of `Σ_{k,i}|density + curvature|` from `x > 0.5 / 0.7 / 0.9` is (median over states)
  `0.004 / 0.000 / 0.000` (`Δ_g = 100`) and `0.124 / 0.005 / 0.000` (`Δ_g = 1000`). The exact scalar remainder is `≤ 0.020 / 0.006 / 0.003 N_s`.
* **Worst-sign sufficient tolerances** (any sign per node; minimum over states). A relative error `δ` of the graviton leading density,
  or of the graviton remainder relative to the scalar one, on `x > x₀`, does not move `n = 1` if:

  | `x₀` | density `δ` (100 / 300 / 1000) | remainder `δ` (100 / 300 / 1000) |
  |---|---|---|
  | 0 | 0.003 / 0.003 / 0.002 | 0.044 / 0.27 / 1.35 |
  | 0.5 | 0.049 / 0.013 / 0.008 | 0.15 / 0.40 / 1.7 |
  | 0.7 | 0.70 / 0.086 / 0.032 | 5.3 / 1.4 / 3.1 |
  | 0.9 | 182 / 27 / 2.2 | 265 / 82 / 22 |

  The small-`x₀` density numbers are dominated by `k`-cancellations at high `J` (`J = 20–28`) and are very conservative.
* **Coherent deviations** (`F → F + ε Σ_{k, x>x₀} c_k T_k`, with `c_k = 1`, `(−1)^k` or `[k≥4]`). Interval of `ε` keeping `n = 1`:

  | `T`, `x₀` | `Δ_g = 100` | 300 | 1000 |
  |---|---|---|---|
  | density, `x > 0.3` (`c = 1`) | (−0.079, +0.065) | (−0.049, +0.079) | (−0.048, +0.076) |
  | density, `x > 0.5` | (−0.32, +0.22) | (−0.15, +0.10) | (−0.12, +0.077) |
  | density, `x > 0.7` | (−1.8, +2.3) | (−0.60, +0.63) | (−0.54, +0.62) |
  | remainder, `x > 0.3` | (−2.3, +1.7) | (−16, +49) | (−144, +340) |
  | remainder, `x > 0.5` | (−1.2, +1.9) | (−10, +6.7) | (−83, +52) |

* **Reading.**
  * The statement "98 % of the near-gap action is the flat B′ density **at `x → 1`**" should read "**at `x ≲ 0.7`**".
  * Beyond `x ≈ 0.7` the graviton density may be off coherently by `≈ 20–50 %` (`Δ_g = 300–1000`, depending on the sign pattern; the
    `k ≥ 4`-only pattern is the tightest, `(−0.26, +0.22)` at `Δ_g = 1000`), or by a factor `≈ 2` (`Δ_g = 100`).
  * The load-bearing window for both H2 (`≈ 5–30 %` accuracy needed) and H4 (`O(1)` factor at `Δ_g = 100`, loose at `Δ_g ≥ 300`) is
    `x ∈ [0.3, 0.7]`.
  * The current position-space reach is `x ≲ 0.1–0.17`.

### 4.4 Cost of the `η`-inversion (`f4_etainv.py`) [D-sketch + NE]

* **Change of weights.** With the large-`ν` symbol of §2, `|φ̃/φ − 1|` is `≤ 10⁻²` on most of `p ∈ [0.1, 0.9]` but reaches `0.4–0.7` at
  `p = 0.9`, `Δ_g = 100`. These are the `(1−p)⁵` zeros, `φ″/φ ∼ 20/(Λ(1−p))²`.
* **Action shift.** `max|ΔF|/N_s = 0.087 / 0.0088 / 0.0012` (`Δ_g = 100 / 300 / 1000`); nodes with `ν ≥ 3` only: `0.021 / 0.0058 / 0.0011`.
  `n = 1` unchanged; min `F/N_s` `0.2256 / 0.2056 / 0.2022`.
* **Reliability.** The symbol is unreliable at `ν = O(1)` and in the `p → 1` boundary layer. The exact difference-equation solution is M3.
* **Alternative.** Keep `φ_k` on `Ψ^{CLS}` with no inversion. Then the bulk-point action is NF02 exactly, but the Regge-side term
  `2ε_kσ^k P_{−½+iν}(cosh β)` becomes a combination of `P_{−½+i(ν±im)}`. Which is cheaper is to be decided in M3.

## 5. What the `x ≲ 0.1` data can and cannot test

* **All objects agree there.** `Ψ_k`, `Ĉ_k`, `Π̂` and `C` agree at `x ≲ 0.1` up to `O(x⁴)` in the bulk-point density
  (`|R_k − 1| ≤ 5·10⁻⁵` at `x = 0.1`; `|Π̂/C − 1| ≤ 5·10⁻⁵`) and up to common `O(1/m²)` terms. The `⟨TTTT⟩` position-space moments have
  `10⁻³` precision. So **no existing data can distinguish route A′ from the current objects**; the comparison is by construction
  (same leading Regge kernel).
* **What the data do test.** They test the common leading Regge formula (H1) and reading B (H2) at small `x`. The latter is now implied
  at leading order by CLS-D's block ↔ partial-wave map at the crystal (`Blockvspartialwave`).
* **What would discriminate.** An `x ≈ 0.5` test of the `k = 5` density (`R₅ − 1 = 4 %`) or of `k = 4, 6` at `x ≈ 0.7` (`8–13 %`).

## 6. Grades and effect on `RESULT_D.md`

| statement | grade |
|---|---|
| CLS-D functionals with `f^±_{k;J_e}` are valid spinning sum rules; scalar version = CHMRS `C̃`; bulk-point action = flat graviton (shock) sum rules, general `d` | [Q] |
| the project's `Ψ_k` = CLS-D functionals with `f → A_k`; `A₃`, `A₇` = CLS-D's; `A₄, A₅, A₆` ≠ | [D] |
| crystal = CHMRS frame `(η=1, r*)`; 9 CLS-D crystal values reproduced | [DV] |
| weight-ratio rule `f₁/f₂|crystal` for the bulk-point action | [D from Q] |
| scalar cross-check reproduces `OBJECT_IDENTITY` (`k = 2, 4, 6`) | [DV] |
| graviton densities `R_k` of the current `Ψ_k` | [D given CLS-D, `d = 3` specialisation [I]] |
| flat limit of the current `Ω^{NF,SC}` negative at `+−`, `m/M < 1.0003`, `J = 16–24`; AdS-reachable only for `Δ_g ≳ 6·10⁴` | [NE] (given `R_k`) |
| near-gap `n = 1` for both objects at `Δ_g = 100/300/1000`; margins of §4.2 | [NE + M] (bulk-point-factor model) |
| tolerances / load-bearing window `x ∈ [0.3, 0.7]` | [NE + M] |
| route A′ leading Regge kernel = project's (after `η`-inversion); `η`-inversion cost `≤ 0.09 N_s` (`Δ_g = 100`) | [D-sketch] / [D-sketch + NE] |

**Hypotheses (proposed regrading; `RESULT_D.md` not edited):**
* **H0** [I] → **[D-sketch]**: the object exists in the literature (CLS-D) up to the `η`-inversion. Open: `d = 3` specialisation (M2);
  inversion analyticity and normalisation (M3); endpoint convergence (M4); light side (M5).
* **H14 / F1:**
  * For the current `Ψ_k` the flat limit is now known exactly (§3). It is neither NF02 nor positive (§4.1), so F1 is **confirmed and
    sharpened**.
  * Under route A′ the flat limit is NF02 [Q + D], so the H14′ part of H0 is reduced to M2.
* **H2:**
  * For route A′, at leading order in `1/m` and fixed `x < 1`: [Q] (CLS-D), pending M2.
  * For the current `Ψ_k`: the B′ assumption is **false** for `k = 4, 5, 6` at `x ≳ 0.5`. In the model this does not change `n`.
* **H4:** still [I]. It is load-bearing only on `x ∈ [0.3, 0.7]`, with the tolerances of §4.3. The `x → 1` region is irrelevant to the
  threshold at the stated tolerances.
* **H9b:** analysable under route A′ (meromorphic weights) [D-sketch]; for the current even `Ψ_k` it stays [I].
* **H15** (`|ĝ₃|²`-free `Ψ_{5,6,7}`): under A′, CLS-D's Regge counting gives convergence only for `J₀ < k` (even `j′`), `J₀ < k − 1` (odd
  `j′`); single-ray endpoints still [I] (M4).

## 7. Plan for the remaining milestones

| # | milestone | content | validation | estimate |
|---|---|---|---|---|
| M2 | `d = 3` specialisation of CLS-D | `SO(d−2)` trivial, `j′ ∈ {0,1}`, `S⁰` crystal locus; 4D helicity form of the 7 + 3 lowest-subtracted sum rules; identify the `(a)` combinations that give `K*`/`𝔖_ff`/`𝔄` and the MHV `B_k^{(1)}[f]` | reproduce the project's leading Regge pair responses (`laplace.py` M-tables, `N_k`) from CLS-D's large-`ν` dual structures; reproduce reading B (`d^J_{44}`) at the crystal | 3–5 days |
| M3 | exact `η`-inversion | solve the difference equation for `j = 1, 2, 3` (or Mehler–Fock of `K_φ/η^j`); analyticity strip of `Ψ^{CLS}_ν`; support of `φ̃` beyond `p = 1`; `N̂_k(ν)` zeros; recompute `k(0)`, `A`, `B_AdS`, H17 | the large-`ν` symbol of §4.4; scalar `C₂` = CHMRS identity | ≈ 1 week |
| M4 | validity | endpoint convergence of `Ψ^{CLS}` with `f^±_{k;2}` for `J₀ ≤ 2` (Einstein-free for `k ≥ 4`) and `J₀ ≤ 4` (H15); I14 swappability | KKSDZ / CLS-D counting; scalar CHMRS case | 1–2 weeks (to [D-sketch]) |
| M5 | light side | graviton pole in `Ĉ₂` (H6/H7) via the unsubtracted `f⁺_{2;2}`; contacts; `[TT]` double-trace content (H9b) | flat dictionary of contacts (`B̂_{k,T}`-type, polynomial in `T`) | 2–3 weeks |
| M6 | finite-`m` correction at fixed `x` (H4 at `O(1/m²)`) | expand the CLS-D crystal integral to next order; the fast part is spin-independent, so the spin dependence is in derivatives of algebraic slow factors | the scalar version must reproduce the CHMRS `1/m²` series and exactskel at `x ≤ 0.6`; `⟨TTTT⟩` position-space data at `x ≤ 0.1` | 2–4 weeks |
| M7 | exact Mellin evaluation of `Ĉ_k` on heavy `⟨TTTT⟩` blocks | spinning Mellin in a structure basis; only needed if M6's accuracy at `x ∈ [0.3, 0.7]` misses the §4.3 targets | — | months |

**Shortest useful sequence:**
* **M2 → M3.** This redefines D-safe on route A′ with the same `n = 1` numbers (now for a defined functional whose flat limit is the
  certified NF02).
* **M6.** This gives H4 on the load-bearing window.
* **M4, M5.** These close H0's validity and light-side clauses.

## 8. Files (`ads4/code/spinfunc/`)

| file | content | output (`cache/`) |
|---|---|---|
| `f1_bulkpoint.py` | sympy: CLS-D crystal ↔ CHMRS frame; leading Regge forms; scalar cross-check; graviton `R_k` | `f1_bulkpoint.txt` |
| `spf_common.py` | `R_X(k,x)` for objects A′ / PsiG / PsiG2 / PsiS; near-gap integrand pieces for any `Δ_g`; state loader (e1 + e7) | — |
| `f2_flat.py`, `f2b_check.py` | flat limit with graviton densities; robustness; twist-gap reach | `f2_flat.txt`, `f2b_check.txt` |
| `f3_ads.py` | near-gap model at `Δ_g = 100/300/1000`; worst-sign tolerances; shares in `x` | `f3_ads_<Δ_g>.txt`, `f3_<Δ_g>_recs.pkl` |
| `f4_etainv.py` | `η`-inversion cost (large-`ν` symbol) | `f4_etainv.txt` |
| `f5_coherent.py` | coherent deviation tolerances | `f5_coherent.txt` |
| `run_all.sh` | runs all of the above | — |
