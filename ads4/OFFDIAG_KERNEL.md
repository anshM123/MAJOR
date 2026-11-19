# The exact leading-Regge helicity off-diagonal of `Ĉ₂′` and `Ĉ₃` (route A′), and what it does to the large-`J` obstruction

*19 Nov 2026, off-diagonal-kernel agent. Unrefereed, not committed. Nothing in `RESULT_D.md`, `STATUS.md`, `H0_SPINNING_FUNCTIONAL.md`,
`H7_I14_LIGHTSIDE.md`, `LARGEJ_OFFDIAG.md`, `S1_S4_FIXES.md` was edited. New code: `ads4/code/spinfunc/s1_offdiag_kernel.py`,
`s1_validate.py`, `s2_lp.py`, `s3_structures.py`; caches `spinfunc/cache/s1_*`, `s2_*`, `s3_*`. `v3_regge_kernel`, `q1_c3_regge_kernel`,
`v4_chi`, `r3_redesign`, `r5_universal`, `q3_analysis`, `robust/nf_lp` are imported read-only.*

Grades: **[P]** proved, **[D]** derived, **[DV]** derived + validated numerically, **[NE]** numerical evidence, **[M]** model, **[I]** open.
Conventions as `LARGEJ_OFFDIAG.md`: `Λ = √((Δ_g−2)(Δ_g−1))`, `b̂ = (2J+1)Λ/m`, `β = arccosh(1 + (2J+1)²/(2m²))`, `t = νβ`, `Γ_off` in l9 units
(`F₊₋` density per unit `ψ₂`-weight `= eo2·Γ_off/m`).

## 0. Summary

1. **The tension is only apparent** [D]. H7 §18 (and S1_S4 §1) computed only helicity-*preserving* tensors (the two diagonal entries);
   "`++ = +−`" and "`h₀ = 0`" say nothing about the `++ ↔ +−` element.
2. **Exact kernel** [D + DV]. The split representation extends to independent polarisations without new symbolic input. At leading Regge
   order: **`Ĉ₃` has no helicity off-diagonal at all** (`o_H` is `2↔4`-antisymmetric and vanishes on the flipped dual tensor) [P]; for `Ĉ₂′`
   only the flip of the dual line (2,4) survives, so in the parallel-transported basis the heavy `2×2` is exactly `[[P_l, Re Y],[Re Y, P_l]]`.
   Scaling form `Y = β²[R(t) + iνF(t)]`, `R = (47/4)F`, `∫tF = 1`, **`G₁ = ∫tR dt = 47/4`** [NE to 10⁻³].
3. **Validation** [NE]. The kernel reproduces the CFT chain's *imaginary* off-diagonal with no free parameter (2.5–4 % on the `Δ = 1600`
   rows). The CFT *real* part equals `Re Y` plus a per-state basis shear `s·P_l`, `s/β² = −1.07 ± 0.13` (residual 0.02–0.1 `β²`); the
   direct RQMC point `(672,36)` is reproduced to `0.008β²`. r4's universal `g ≈ −R + 1.00·J₀`: its `J₀` part is the shear (a change of OPE
   basis, verdict-neutral), and its `t ≤ 12` truncation inflated `|G₁|` from 11.75 to 12.5.
4. **Large `β`** [NE]. The exact first moment `M_off(β)` is 11.75 → 10.5 (`β = 0.1`) → 7.2 (0.25) → 4.5 (0.5) → 4.5 (1.1) → 0.5 (1.7);
   the weight-independent ratio `M_off/M₁` stays `≥ 110` for `β ≤ 1.1`. The extrapolated small-`β` `g` overestimates the off-diagonal moment by
   1.7–2.8× at `β = 0.25–0.5`, but not by enough to matter.
5. **The obstruction is real** for `Ω^{A′}` as defined [M]: with the exact kernel, NF02 fails from `b̂ ≈ 32–37` (300), and the NF02-family
   LP is infeasible on any window reaching `β ≈ 0.4` (300), `0.2` (1000), `0.12` (3000) — the twist-allowed band reaches `β ≈ 1.8`.
6. **Other `k = 2` structures** [DV + NE + M]. Regge-side check of H0 §12: the unsuppressed flat flip of `P², PH, H²` is `1 : 2 : 4`, so
   exactly the crystal flip-free family cancels it. In that family, `κ* = −0.4332` (`S = q² + κ*·Hh·q`) kills the scaling-limit off-diagonal
   **identically**, but the finite-`β` remainder (`M_off ≈ 1.4–1.8` at `β = 0.25–0.5`) still makes the band infeasible. Feasible windows
   double; at the nominal windows the bounds are `3.06·10⁻⁶` (300, `β ≤ 0.3`) / `2.92·10⁻⁸` (1000, `β ≤ 0.2`) [M].
   **`EE = V₀⁴(z₂′·z₄′)²` has identically zero leading-Regge off-diagonal** and the same Regge-overlap normalisation as `q²`; with it,
   NF02 is positive on the whole band and the bounds stay `3.12·10⁻⁶ / 2.89·10⁻⁸` [M]. But its forward normalisation has real zeros
   (`ν ≈ 1.17, 6.64`: singular own de-shift) and its crystal flip content is unknown [I]. These two items decide whether route A′ can be
   repaired at leading Regge order.

## 1. The kernel

### 1.1 The tension with H7 §18 (answer to the preliminary question) [D]
`v3_regge_kernel.heavy_image` (H7 §18.2) and `q1_c3_regge_kernel.heavy_image3` (S1_S4 §1.2) only ever use **preserving** tensors
(`τ₁₃ = ε_{h₁}²⊗ε_{−h₁}²`, `τ₂₄ = ε_{h₂}²⊗ε_{−h₂}²`). The two helicity classes `(h₁,h₂) = (+,−)` [`++`] and `(+,+)` [`+−`] are the two
**diagonal** entries. "`++ = +−`" and "`h₀ = 0`" are statements about the diagonal; the equivariance proof of §18.3 uses weight-0 tensors
explicitly. **The helicity-changing matrix element was never computed**, so there is no conflict between §18 and the large-`J`
off-diagonal: §18 is silent about it.

### 1.2 Extension to helicity-changing polarisations [D]
Same split representation (§18.1 / S1_S4 §1.1, same assumptions: `Ψ⁺` reading, ν-only constants dropped, Lemma-V referral), with the four
polarisations independent: `z₁ = ε_{a₁}(p̄̂)`, `z₃ = ε_{a₃}(p̄̂)`, `z₂′ = ε_{a₂}(p̂)`, `z₄′ = ε_{a₄}(p̂)` in the parallel-transported frame
(`v3.frames`). The transform engine and monomial expansions of v3/q1 are already general in the polarisations (they use only
`ε² = 0`), so no new symbolic input is needed (`s1_offdiag_kernel.img2/img3/img_vec`). The parity-even heavy `2×2` (even `J`) is
`F_ab = ½Σ_{(h₁,h₂)∈a}Σ_{(k₁,k₂)∈b} L(a₁=h₁, a₃=−k₁, a₂=h₂, a₄=−k₂)`, so the off-diagonal contains only single flips (one line preserved,
one flipped) and the diagonal contains preserved terms plus double flips.

**Which configurations survive** (`s1_validate.py (a)`, `cache/s1_validate.txt`):

| configuration | `Ĉ₂′` (`q²`, `j′=0`, η-kernel) | `Ĉ₃` (`o_H q`, `j′=1`) |
|---|---|---|
| flip on the dual line (2,4) (`z₂′ = z₄′`) | **non-zero**, `O(β)` (relative to `P_l`) | `0` identically |
| flip on the light-ray line (1,3) (`z₁ = z₃`) | `≤ 10⁻¹⁴` | `≤ 10⁻¹⁷` |
| double flip | `≤ 10⁻¹⁵` | `≤ 10⁻³³` |

* **`Ĉ₃` has no leading-Regge helicity off-diagonal at all [P given §1.1 of S1_S4]:** `o_H(x₀=0) = (z₂′·z₀)ε₃(z₄′,x̂,z₀) −
  (z₄′·z₀)ε₃(z₂′,x̂,z₀)` is antisymmetric in `z₂′ ↔ z₄′` and vanishes identically on the flipped dual tensor `z₂′ = z₄′`; line-(1,3) flips
  vanish numerically (`≤ 10⁻¹⁷`) [DV]. The `g₃` that `LARGEJ` fitted is therefore not a leading-Regge object; the C3H K₂ ladders do
  not show leading-Regge scaling (noisy, §1.4).
* **`Ĉ₂′`: the whole off-diagonal is the dual-line flip** `F₊₋ = L(+,−; −,−)` (both parity partners equal), and the diagonal is
  exactly helicity-blind with no double-flip admixture [DV]. In the geometric (parallel-transported) helicity basis the leading-Regge
  heavy `2×2` of the de-shifted functional is therefore exactly
  **`[[P_l(cosh β), Y],[Y, P_l(cosh β)]]`**, `Y(ν,β) := L_flip(ν,β)/L_pres(ν,0)` (forward-normalised), with `Re` taken at the end.

### 1.3 Scaling limit and two exact-looking moments (`s1_validate.py (b)`, `cache/s1_scaling.pkl`) [NE]
At `ν → ∞`, `t = νβ` fixed (`ν = 500` and `2000` agree to `≤ 10⁻³`):
**`Y(ν, β) = β²[R(t) + iνF(t)] (1 + O(1/ν))`, with `R(t) = (47/4)·F(t)` pointwise** (ratio `11.750 ± 0.001` where `F ≠ 0`, `ν = 2000`), and
**`∫₀^∞ tF(t)dt = 1.0000 ± 0.005`, `∫₀^∞ tR(t)dt = 11.750 ± 0.06`** (Abel mean of the oscillating tail, `t ≤ 80`). So the imaginary
part is `O(β)` relative to the diagonal and the real part `O(β²)`; the real part is the `47/(4ν)`-rotated image of the imaginary one.

| `t` | 0.5 | 1 | 1.5 | 2 | 2.5 | 3 | 3.7 | 4.5 | 5 | 6 | 7 | 8 | 10 | 12 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `R(t) = Re Y/β²` (geometric basis) | 0.030 | 0.116 | 0.217 | 0.399 | 0.525 | 0.689 | 0.797 | 0.817 | 0.735 | 0.467 | 0.151 | −0.077 | −0.103 | 0.060 |

`R ~ t²` at small `t`, single positive lobe up to `t ≈ 7.5`, then a decaying oscillating tail.

### 1.4 Validation against the CFT chain (`s1_validate.py (c)–(e)`, `cache/s1_validate.txt`) [NE]
**Imaginary part — direct test of normalisation and helicity combination.** The K₂ rows of q3 (r4's 25 rows) carry `Im Γ_off` (dropped
by `Re`, kept "for the record"). The kernel reproduces it with **no free parameter** (only the arbitrary row sign): `Im_CFT/Im_ker` median
`0.99–1.003` and rms deviation `2.5–4 %` on all `Δ = 1600` rows with `J ≥ 10`, `5–8 %` at `Δ = 800`, `10–30 %` at `Δ = 336–400` (finite-`m`
/ K₂ truncation, growing with `x = ν/m`). This fixes the normalisation (forward-normalised by the `++` reference, as l9) and confirms that the
CFT off-diagonal is the dual-line flip.

**Real part — a per-state change of OPE basis.** `Re_CFT ≠ Re Y`: e.g. at `t ≈ 1` the CFT has `g ≈ 0.8` where `R = 0.12`. The two-parameter
fit `Re Γ_off^{CFT}/m = Re Y + (c/ν)Im Y + s·Re D` (`D` = forward-normalised preserving image at the state, i.e. `P_l` + shift terms;
`s` constant per state = a real shear of the OPE basis, `c` = a ν-only phase) gives on `Δ ≥ 800`: `c ≈ 0` (`|c| ≤ 1`, i.e. no phase), and
**`s/β² = −1.07 ± 0.13`** (shear alone; range `[−1.29, −0.86]`), residual rms `0.02–0.12 β²` against `0.6–0.8 β²` for the raw CFT
values. **One constant reproduces r4's universal function:** `g_r4(t) ≈ −R(t) + σJ₀(t)` with `σ = 1.00`, rms `0.116` (r4's own spread 0.092):

| `t` | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 10 | 12 |
|---|---|---|---|---|---|---|---|---|---|---|
| `g_r4` (median) | +0.80 | −0.30 | −1.10 | −1.23 | −0.85 | −0.25 | +0.11 | +0.08 | −0.11 | +0.09 |
| `−R + 1.00·J₀` | +0.65 | −0.18 | −0.95 | −1.22 | −0.91 | −0.32 | +0.15 | +0.25 | −0.14 | −0.01 |
| `−R` alone (geometric basis) | −0.12 | −0.40 | −0.69 | −0.83 | −0.74 | −0.47 | −0.15 | +0.08 | +0.10 | −0.06 |

**Direct (RQMC) ladder points** (LARGEJ §2.3, `Δ_g = 300`, `ν = 25 … 66.8`, C2): kernel + one shear per state reproduces them with
`s/β² = −1.22 / −1.21` and rms residual **`0.069 β²` at `(336,18)` and `0.008 β²` at `(672,36)`** (the more Regge-like point):
`(672,36)` direct `0.65, 1.18, 1.12, 0.62, 0.05, −0.27` vs kernel + shear `0.64, 1.18, 1.13, 0.63, 0.06, −0.28`.

**Interpretation [D for the algebra, I for the origin of `s`].** The CFT chain's helicity basis is the rank-1 basis of the crystal
order-0 matrix on `x ≥ 0.15` (`l8_analysis.basis`), fixed per state. A real per-state change of basis `e₊₊ → e₊₊ + a e₊₋`, `e₊₋ → e₊₋ +
b e₊₊` maps `[[d, Y],[Y, d]]` to `[[d + 2aY + a²d, Y + (a+b)d], [·, d + 2bY + b²d]]`: it adds `s·d` (`s = a+b ≈ −β²`) to the off-diagonal
and only `O(β⁴)` to the diagonal, and leaves `det` and the PSD verdict invariant. The fitted `s` is per-state constant, universal in
`s/β²`, and the residual is at the K₂ noise level, so **the `J₀`-like part of r4's `g` is a basis effect, and the physical
(basis-invariant) real off-diagonal is `R`** [NE]. Why the crystal basis is sheared by `≈ −β²` relative to the parallel-transported one
is not derived [I]; it is irrelevant for the verdict *provided the diagonal is taken in the same basis*. LARGEJ's model took the diagonal in
the geometric basis (`flat + curvature`) and the off-diagonal in the sheared one; after `ψ`-integration the shear term is `s·(diagonal)`,
`|s| ≤ 0.1`, harmless — except that r4/r5 integrated the `J₀` piece only on `t ∈ [1, 12]` (held below 1), which does not cancel and adds
`≈ −2.7σ` to `G₁`.

**`Ĉ₃` rows.** The kernel's `Ĉ₃` off-diagonal is zero (§1.2). The CFT `Re Γ₃,off` scales by `5.6 / 1.04 / 1.9 / 2.65` under `m → 2m`
at fixed `(ν, β)` on four K₂ ladders (leading Regge would be 2.00 each; C2 gives 1.75–2.5), and is equal at `(336,18)` and `(672,36)`
(LARGEJ §2.3) — not a leading-Regge object; noisy [NE]. It enters with `uψ₃`; in the exact leading-Regge model it is set to 0.

### 1.5 Answer to task 1 [DV + NE]
* **Exact leading-Regge off-diagonal of `Ĉ₂′`:** `Γ_off/m = Re Y(ν, β)` in the geometric basis, `Y = L_flip/L_pres(β=0)`, computed
  exactly (no fit; φ-trapezoid with `nphi = 192`: at `nphi = 128` the error vs `512` is `≤ 3·10⁻⁵` for `t ≤ 100`, `≤ 10⁻⁸` for `t ≤ 20`) on `β ∈ [0.005, 2.1]`, `ν ≤ 3000`
  (`cache/s1_table_C2a.pkl`; s1 mode `table`). Scaling form `β²R(νβ)` with `G₁ ≡ ∫tR dt = 47/4` [NE to 10⁻³].
* **`Ĉ₃`: identically zero** at leading Regge [P for the dual-line flip; DV for the light-ray-line flip].
* Finite-`β` values: §2.1.

## 2. The obstruction with the exact kernel (`s2_lp.py`; `cache/s2_*`)

### 2.1 Finite `β` and the first moment (`s2_lp.py moment`, `cache/s2_moment.txt`) [NE]
Exact `M_off(β) := ∫₀^∞ν Re Y(ν,β) dν` (no de-shift on the off-diagonal; Abel mean of the tail) against the diagonal's curvature moment
`M₁(β) = ∫2νP_l/(e^{2πν}+1)` (ROUTEA_LARGE_B; `M₁(0) = 1/24`). In the curvature regime both enter as `2ψ₂′(0)·M/Λ²`, so `M_off/M₁` is the
weight-independent asymptotic off-diagonal/diagonal ratio (LARGEJ's `24|G₁|`):

| `β` | 0.005 | 0.01 | 0.02 | 0.05 | 0.1 | 0.16 | 0.25 | 0.35 | 0.5 | 0.75 | 1.1 | 1.7 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `M_off` | 11.6 ± 0.7 | 11.75 | 11.68 | 11.38 | 10.51 | 9.19 | 7.19 | 5.44 | 4.53 | 7.43 | 4.52 | 0.46 |
| `M_off/M₁` | 278 | 282 | 280 | 273 | 253 | 221 | 174 | 132 | 112 | 189 | 123 | 15 |

* `M_off(0) = 47/4` (§1.3), against LARGEJ's `|G₁| ≈ 12.5` from the median `g` (whose `J₀` basis part, truncated at `t = 12`, contributes
  `≈ 2.7`; the rest of the difference is the hold `g(t<1) = g(1)`). **The extrapolated small-`β` form (`G₁ = 12.5`, `β`-independent) is right to
  ≤ 20 % for `β ≤ 0.1` and overestimates the off-diagonal moment by a factor 1.7–2.8 at `β = 0.25–0.5`**; the ratio stays `≥ 110` up to `β ≈ 1.1` (`b̂ ≈ 1.3Λ`) and only
  collapses near `β ≈ 1.7`. (`β ≥ 0.5`: the un-de-shifted `Re Y` is sensitive to the near-zeros of `h₋+h₊` at `ν ≈ 0.73, 1.9, 6.4`, which
  explains the non-monotonic 0.75 entry; the de-shifted (χ) smearing used below removes that sensitivity: `X = χ(h₋+h₊)` dips where
  `h₋+h₊` does.)

### 2.2 NF02 per state (`s2_lp.py check`, `cache/s2_check_300_b2.2_{psi,chi}.txt`) [NE + M]
Model M′: diagonal = r5's (stripped flat `B₂+B₃` + leading AdS curvature `2∫ψ₂(P_l−J₀) + uσ∫ψ₃(P_l−J₀)`; = the exact leading-Regge
diagonal in the geometric basis), off-diagonal = exact kernel (`OFF=chi`: `Re∫X Y dν/Λ`, `X` = v4_chi min-norm solution for each weight;
`OFF=psi`: `∫ψ₂ Re Y`), `Ĉ₃` off-diagonal 0, Gauss nodes `96 + 3b̂` (so large `b̂` is resolved — r5 used 64). States: r5's grid
(`m/Λ ∈ [1.02, 5]`, `b̂ ∈ [10, 2Λ]`, twist-allowed, 319 states at `Δ_g = 300`).
* `Δ_g = 300`: first negative `b̂ = 37 / 34 / 35 / 33 / 32.5 / 32 / 46` at `m/Λ = 1.1 / 1.2 / 1.3 / 1.5 / 2 / 3 / 5` (LARGEJ: 30–35, 46);
  `λ/N_s ≈ −0.5` at `b̂ ≈ 37`, `−4` at 56, `−17` at 100, `−30` at 139, `−60 … −190` at `b̂ = 200–480`. `psi` and `chi` agree to
  `≤ 10 %` for `β ≤ 0.5`.

### 2.3 The LP (NF02's family, r3.solve, flat margin `δ = 0.2`, AdS margin `d_w = 0`) (`s2_lp.py lp`) [M]

| `Δ_g` | AdS rows | `d` | result | bound 300 / 1000 / 3000 |
|---|---|---|---|---|
| 300 | `β ≤ 0.12` | 8 | `A = 3049.7`, `B_AdS = −7012.7` | `3.08·10⁻⁶` / `2.86·10⁻⁸` / — |
| 300 | `β ≤ 0.3` (`b̂ ≤ 90`) | 8 | `A = 10794`, `B_AdS = −41845` | `4.24·10⁻⁶` / — |
| 300 | `β ≤ 0.4` | 8, 12 | **infeasible** | |
| 300 | `β ≤ 0.5 / 1.0 / 2.2` | 8 (0.5: also 12) | **infeasible** | |
| 1000 | `β ≤ 0.12` (`b̂ ≤ 120`) | 8 | `A = 4452.8`, `B_AdS = −14003` | — / `3.13·10⁻⁸` |
| 1000 | `β ≤ 0.2 / 0.3 (d=12) / 0.5 / 1.0 / 2.2` | 8 | **infeasible** | |
| 3000 | `β ≤ 0.12 / 0.2 / 0.3` | 8 | **infeasible** | |

(LARGEJ with the fitted `g`: 300, `β ≤ 0.3` infeasible at `d = 8`, `9.05·10⁻⁶` at `d = 12`; 1000, `β ≤ 0.12` → `3.25·10⁻⁸`, `β ≤ 0.2`
infeasible.) The twist-allowed band reaches `β ≈ 1.8` (`b̂ = 2Λ`) on this grid.

### 2.4 Answer to task 2
* **Exact `g` at large `β` vs the extrapolation:** right in the scaling regime (`G₁ = 47/4` vs `12.5`), smaller by 1.4–2.5 at
  `β = 0.25–0.5` (first moment), non-zero and same sign through `β ≈ 1.1`.
* **The obstruction is real for `Ω^{A′}` as defined** (`Ĉ₂′` with `q²`) [M + NE]: with the exact leading-Regge off-diagonal the
  NF02-family LP is feasible only on windows `β ≤ 0.3` (300) / `0.12` (1000), and infeasible on any window reaching `β ≈ 0.4` (300) /
  `0.2` (1000) / `0.12` (3000), far inside the twist-allowed band. The weight-independent mechanism (LARGEJ §3, `M_off/M₁ ≫ 1` for
  `ψ₂′(0) ≠ 0`) survives with `M_off/M₁ ≥ 110` on `β ≤ 1.1`.
* What is still model: graviton diagonal corrections beyond leading Regge and the scalar remainder (dropped, as in LARGEJ), the `O(1/m)`
  `Ĉ₃` off-diagonal (dropped), the χ null-space gauge (min-norm) at `ν ≲ 18`.

## 3. Other k = 2 (2,4) structures and a cancelling combination (`s3_structures.py scan`, `cache/s3_scan.txt`)

### 3.1 Structures [D]
At `x₀ = 0` (H7 §18.1 reduction, inversion-referred `z′`) the `j′ = 0` building blocks are `H₂₀ = z₂′·z₀`, `H₄₀ = z₄′·z₀`,
`V₀ = −z₀·x̂`, `V₂ = z₂′·x̂`, `V₄ = z₄′·x̂` (`V_{2,04} = 2z₂·x₂ → −2|x|V₂` by the referral) and `Ĥ₂₄ = z₂′·I(x̂)·z₄′ = e₂₄ + 2V₂V₄`.
With `P := H₂₀H₄₀`, `Hh := V₀²Ĥ₂₄`, v3's `q = 2P − Hh`. Computed (all with the `16/V₀⁴` dual normalisation of v3):
`q²`, `PP = P²`, `PH = P·Hh`, `HH = Hh²`, **`Nq := Hh·q = 2PH − HH`**, and the `2↔4`-asymmetric `A1 = H₂₀²V₀²V₄²` (+ mirror),
`A2 = P V₀²V₂V₄`, `A3 = H₂₀²H₄₀V₀V₄` (+ mirror), `VV = V₀⁴V₂²V₄²` (`s1_offdiag_kernel.STRUCTS`).

### 3.2 Results [DV + NE]
* **Every structure:** line-(1,3) flips `≤ 6·10⁻¹⁶`; preserving image exactly 3-term with `h₀ = 0` (`|h₀|/(|h₋|+|h₊|) ≤ 10⁻¹³`), as
  H7 §18.3 proves for any weight-0 tensors. So each structure has its own spinning de-shift and satisfies (HN) by construction.
* **Forward normalisation** `n_S(ν)/n_{q²}(ν)` at `ν = 1 … 3000`: `PP` and all asymmetric ones `→ 1/4`; **`PH`, `HH`, `Nq` → 0**
  (`Nq`: `0.70, 0.72, 0.21, 0.031, 3.0·10⁻³, 3.3·10⁻⁴, 3.0·10⁻⁵, 3.3·10⁻⁶` at `ν = 1, 3, 10, 30, 100, 300, 1000, 3000`, i.e. `≈ 30/ν²`).
* **Flat helicity flip in the Regge overlap.** `Re o_S/n_{q²}` for `PP : PH : HH` is `1 : 2.0000 : 4.000` times an `O(ν²β²) = O(t²)`
  function (e.g. `2476.55 : 4952.91 : 9905.56` at `t = 1`, `ν = 2000`), i.e. an unsuppressed helicity flip that cancels **exactly in
  the combinations with `a + 2b + 4c = 0`** — precisely H0 §12's crystal flip-free condition for `aP² + bPH + cH²` (with `P ↔ 𝒱₂𝒱₄`,
  `Hh ↔ H`). `q²` and `Nq` are the two flip-free directions; their residual off-diagonals are `O(β²)`. This is an independent
  Regge-side confirmation of the crystal statement [DV]. All asymmetric structures carry the same unsuppressed flip as `PP` (`A − PP`
  cancels it; their crystal flip content is not computed [I]).
* **Scaling-limit off-diagonals of the two flip-free directions** (`ν = 400, 2000`, agree to `10⁻³`): `Re o_{Nq}/n_{q²} = β²R_N(t)` with
  **`R_N(t) = 2.308·R(t)` pointwise** (`0.269/0.116, 0.922/0.399, 1.840/0.797, 1.078/0.467` at `t = 1, 2, 3.7, 6`), `∫tR_N dt = 27.12`
  (`≈ 217/8`); imaginary parts `∫t Im/ν dt = 1.25` (`Nq`) vs `1.00` (`q²`) — different shape weight, irrelevant after `Re`.
* **Cancellation.** For `S(κ) = q² + κ·Nq` (flip-free at the bulk point, bulk-point density identical to `q²`'s since `n_N/n_q → 0`),
  `G₁(κ) = 11.749 + 27.12κ` (scan: `−1.81` at `κ = −0.5`, `+4.97` at `−0.25`), so **`κ* = −0.4332` (≈ `−94/217`) makes the
  scaling-limit real off-diagonal vanish identically in `t`**, not just its first moment. `S(κ*) = q·(2P − 1.433 Hh)`.

### 3.3 Finite `β` and the LP with `S(κ*)` (`s3_structures.py table2` → `cache/s3_table_Nqa.pkl`, `s3_ratio_Nq.pkl`; `S2_KAPPA=-0.4332 s2_lp.py`) [NE + M]
`Y_S = (Y_q + κ W_N)/(1 + κ r)`, `W_N = o_N/n_q`, `r = n_N/n_q` (exact tables, same nodes as §2; `OFF=psi`, i.e. no de-shift on the
off-diagonal — the `S` de-shift `u_S` was not built).
* First moment `M_off^S(β)`: `0.0013, 0.031, 0.175, 0.55, 1.00, 1.39, 1.46, 1.76` at `β = 0.005, 0.02, 0.05, 0.1, 0.16, 0.25, 0.35, 0.5`
  (`q²`: 11.6 … 4.5), i.e. `≈ 55β²` at small `β`, ratio `M_off/M₁ = 13 / 24 / 34 / 35 / 44` at `β = 0.1 / 0.16 / 0.25 / 0.35 / 0.5`. The
  scaling-limit cancellation is exact, but the `O(β)`-relative remainder still exceeds the curvature diagonal by 1–2 orders in the curvature
  regime. (`β ≥ 0.75` values are not reliable for `S`: `|1 + κ*r|` has a near-zero `0.049` at `ν = 1.89`, which `OFF=psi` does not tame.)
* NF02 per state, `Δ_g = 300`: first negative `b̂ ≈ 64–70` (all `m/Λ ≥ 1.2`; `q²`: 32–37); `λ/N_s ≈ −3` at `b̂ ≈ 100`, `−5 … −9` at 120–140.
* LP (NF02 family, `d = 8`, `δ = 0.2`, `d_w = 0`):

| `Δ_g` | AdS rows | `q²` (§2.3) | `S(κ*)` |
|---|---|---|---|
| 300 | `β ≤ 0.3` | `A = 10794`, `4.24·10⁻⁶` | **`A = 3435.6`, `B_AdS = −9359`, bound `3.06·10⁻⁶`** |
| 300 | `β ≤ 0.4` | infeasible | `A = 4905`, `3.30·10⁻⁶` |
| 300 | `β ≤ 0.5` | infeasible | `A = 16417`, `4.99·10⁻⁶` |
| 300 | whole band (`β ≤ 1.8`) | infeasible | **infeasible** |
| 1000 | `β ≤ 0.12` | `3.13·10⁻⁸` | `A = 3022.3` (= NF02), `2.89·10⁻⁸` |
| 1000 | `β ≤ 0.2` | infeasible | **`A = 3548.8`, `B_AdS = −9902`, `2.92·10⁻⁸`** |
| 1000 | `β ≤ 0.3` / `0.5` / band | infeasible | no convergence / infeasible / infeasible |
| 3000 | `β ≤ 0.12` / `0.2` | infeasible | `3.92·10⁻¹⁰` / `6.0·10⁻¹⁰` (`A = 3045` / `9833`) |

**Verdict for the crystal-verified family [M]:** the best fixed combination is `κ* = −0.4332` (scaling-limit off-diagonal ≡ 0); it roughly
doubles the feasible `b̂`-window but **does not make the band feasible** at `Δ_g = 300, 1000`.

### 3.4 A structure with identically zero leading-Regge off-diagonal: `EE = V₀⁴(z₂′·z₄′)²` [P + NE + I]
* Any structure with a factor `e₂₄ = z₂′·z₄′` vanishes on the flipped dual tensor (`ε·ε = 0`); line-(1,3) flips vanish for every
  structure (§3.2). So **`Ĉ_EE` has no leading-Regge helicity off-diagonal at any `(ν, β)`** [P for the dual line, DV `≤ 10⁻³⁰` numerically].
  `e₂₄ = Ĥ₂₄ − 2V₂V₄` is a covariant structure (`H₂₄ − 2V_{2,04}V_{4,02}` up to normalisation), presumably the covariant form of the old
  Regge-frame `(D₁+D₂)²` ("`z′_in·z′_out`", H0 §12), which H0 conjectured but did not show to map to the crystal `Q`.
* Preserving image 3-term with `h₀ = 0` [DV]; **forward normalisation `n_EE/n_q → 1`** (`0.99996 − 0.0073i` at `ν = 3000`), i.e. the
  same Regge-overlap density as `q²` [NE].
* **With zero off-diagonal, NF02 is positive on every state of the band** in model M′ (`Δ_g = 300`: no negative state, 319 states;
  `Δ_g = 1000`: min `λ/N_s = +0.2056`, 374 states, the LP returns NF02: `A = 3022.34`, `B_AdS = −6552.66`, bounds **`3.12·10⁻⁶` / `2.89·10⁻⁸`**,
  unchanged) [M].
* **Two defects block it [NE + I]:** (i) its forward normalisation `h₋ + h₊` has **real zeros** at `ν ≈ 1.172` and `≈ 6.643` (`n_EE/n_q`
  passes through 0 linearly: e.g. `−0.0011−0.0004i → +0.0060+0.0025i` between `ν = 1.17` and `1.18`), the same kind of defect as the old
  scalar-`g` de-shift at `ν₀ = 5.979`; adding `μ·EH` (`EH = V₀⁴e₂₄Ĥ₂₄`, also Regge-flip-free, `n_EH/n_q → 0`) does not remove the zero near
  `ν ≈ 6` for any `μ ∈ [−3, 3]` (`min|n/n_q| ≤ 6.6·10⁻³`). So (HN) by its own de-shift is singular at two points [NE]. (ii) Its **crystal
  (bulk-point) helicity content is not computed**: flip-freeness with the stripped density shape at finite `x` is exactly what route A′'s
  flat certificate needs, and the Regge overlap only tests `x → 0` [I].
* Interpolations `(1−λ)q² + λEE` scale the large-`ν` off-diagonal by `≈ (1−λ)`; covering the band needs `1−λ ≲ 10⁻²`, where the
  normalisation nearly vanishes at `ν ≈ 1.17, 6.64` — no clean intermediate.

### 3.5 Answer to task 3
* Within the `k = 2`, `j′ = 0` structures whose flat limit is verified (H0 §12's flip-free family = span{`q²`, `Nq`}), **no fixed
  combination makes the LP feasible on the whole band**; the best (`κ* = −0.4332`, cancelling the scaling-limit off-diagonal exactly)
  gives at the nominal windows **`3.06·10⁻⁶` (300, `β ≤ 0.3`) and `2.92·10⁻⁸` (1000, `β ≤ 0.2`)**, i.e. essentially the nominal
  `3.13·10⁻⁶ / 2.89·10⁻⁸`, but only on those windows [M].
* All requirements other than the flat limit are structure-independent once each structure carries its own spinning de-shift: (HN) needs
  `h₀ = 0` (true for every structure) and no real zero of `h₋ + h₊` (true for `q²`, `S(κ*)` [min `|1+κ*r| = 0.049` at `ν = 1.89`],
  false for `EE`); the pole `16πG/(ν²+9/4)` and `A` follow from D3 + the weights (H7 §17, unchanged); `Ĉ₃|_E = 0` is untouched.
* A complete escape at leading Regge exists in principle (`EE`-type, §3.4), conditional on two open items. The asymmetric structures all
  carry an unsuppressed flat flip individually; their flip-free differences (`A − PP`) are Regge-side candidates whose crystal content is
  also open [I].

## 4. Graded summary

| statement | grade |
|---|---|
| H7 §18 / S1_S4 §1 computed only diagonal (preserving) images; no conflict with the large-`J` off-diagonal | [D] (code read) |
| `Ĉ₃` leading-Regge off-diagonal ≡ 0 (dual-line flip: `o_H` antisymmetric; light-ray-line flip numerically `≤ 10⁻¹⁷`) | [P] + [DV] |
| `Ĉ₂′` off-diagonal = dual-line flip only; diagonal helicity-blind, no double flips; `2×2 = [[P_l, Re Y],[Re Y, P_l]]` | [DV] (given H7 §18.1 [D]) |
| Scaling `Y = β²(R + iνF)`, `R = (47/4)F`, `∫tF = 1`, `G₁ = 47/4` | [NE] (`10⁻³`) |
| Kernel = CFT imaginary off-diagonal, no free parameter (2.5–4 % at `Δ = 1600`) | [NE] |
| CFT real off-diagonal = kernel + per-state OPE-basis shear `s ≈ −1.07β²`; direct `(672,36)` to `0.008β²` | [NE]; origin of `s` [I] |
| r4's `g = −R + J₀`; the `J₀` part is basis, verdict-neutral; r4/r5 truncation inflated `|G₁|` 11.75 → 12.5 | [NE] + [D] |
| Exact `M_off(β)`, `M_off/M₁ ≥ 110` on `β ≤ 1.1` | [NE] |
| Obstruction real: NF02 fails from `b̂ ≈ 32–37` (300); LP infeasible on windows reaching `β ≈ 0.4 / 0.2 / 0.12` (300/1000/3000) | [M] (diag = flat + leading curvature; graviton diag corrections, scalar remainder, `O(1/m)` `Ĉ₃` off-diagonal dropped; χ min-norm gauge) |
| Regge confirmation of H0 §12: flat flip `PP : PH : HH = 1 : 2 : 4` | [DV] |
| `S(κ*)`, `κ* = −0.4332`: scaling off-diagonal ≡ 0; band still infeasible; windows 300 `β ≤ 0.5`, 1000 `β ≤ 0.2`, 3000 `β ≤ 0.2` | [NE] + [M] |
| `EE`: off-diagonal ≡ 0 at leading Regge; NF02 positive on the band; bounds unchanged | [P] + [M] |
| `EE`: real zeros of `h₋+h₊` at `ν ≈ 1.17, 6.64`; crystal flip content | [NE]; **[I]** |

## 5. Files and hand-over

**Code** (`ads4/code/spinfunc/`, new; all imports of project modules read-only):
* `s1_offdiag_kernel.py`: general-polarisation heavy images (`img2`, `img3`, vectorised `img_vec`), structure library `STRUCTS` +
  `A_monos_struct`, mode `table <tag> <β-indices>` → `cache/s1_table_<tag>.pkl` (`{β: (ν, Y)}`; `C2a` = even indices of `BETAS`, `C2b`
  = 0.01, 0.03; `C2c` was started for the odd indices and stopped before completing any `β`).
* `s1_validate.py` → `cache/s1_validate.txt`, `s1_scaling.pkl`, `s1_cftrows.pkl`, `s1_validate_fits.pkl` (~2 min + q3 loading).
* `s2_lp.py` (modes `moment`, `check`, `lp`; env `S2_DG`, `S2_BMAX`, `S2_OFF=psi|chi|none`, `S2_KAPPA`, `S2_TABLE`, `S2_TAG`, `R3_D`):
  `cache/s2_moment*.txt`, `s2_check_*`, `s2_lp_*`, `s2_chi_<Δg>_q2_d<d>.pkl` (v4_chi solutions per basis function). Each LP < 2 min.
* `s3_structures.py` (`scan` → `cache/s3_scan.txt/.pkl`; `table2 <tag> <idx>` → `cache/s3_table_Nqa.pkl`; `ratio` →
  `cache/s3_ratio_Nq.pkl`); `cache/s3_ratio_EE.pkl`, `s3_ratio_EH.pkl` were produced inline (same loop as `ratio_table`, `nphi = 96`).

**Not done / next steps (value order):**
1. **Crystal test of `EE` (and `EH`)**: map `V_{2,04}V_{4,02}` to the D = 4 crystal (g1/g2 machinery) and check flip content and the
   stripped `(1−x²)⁻⁴` ratio. If flip-free with the right density, route A′ with `Ĉ_EE` has no leading-Regge off-diagonal (§3.4).
2. **The real zeros of `h₋ + h₊` for `EE`**: exact `h∓^{EE}` by the v4_hclosed method (rational × Γ), then decide whether the Abel-form
   χ is well posed with a real pole of `u` (index count), or find a Regge-flip-free combination (structures with a factor `e₂₄`) without
   real zeros.
3. The de-shift of the off-diagonal for `S(κ*)` (`u_S`) and the null-space gauge dependence of `Re∫X Y` at `ν ≲ 18` (not quantified).
4. The origin of the crystal-basis shear `s ≈ −β²` (§1.4) — irrelevant for the verdict, but needed to compare the CFT chain's diagonal
   corrections with the leading-Regge ones in one basis.
