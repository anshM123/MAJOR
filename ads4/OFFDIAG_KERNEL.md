# The exact leading-Regge helicity off-diagonal of `Ĉ₂′` and `Ĉ₃` (route A′), and what it does to the large-`J` obstruction

*19 Nov 2026, off-diagonal-kernel agent. Unrefereed, not committed. Nothing in `RESULT_D.md`, `STATUS.md`, `H0_SPINNING_FUNCTIONAL.md`,
`H7_I14_LIGHTSIDE.md`, `LARGEJ_OFFDIAG.md`, `S1_S4_FIXES.md` was edited. New code: `ads4/code/spinfunc/s1_offdiag_kernel.py`,
`s1_validate.py`, `s2_lp.py`, `s3_structures.py`; caches `spinfunc/cache/s1_*`, `s2_*`, `s3_*`. `v3_regge_kernel`, `q1_c3_regge_kernel`,
`v4_chi`, `r3_redesign`, `r5_universal`, `q3_analysis`, `robust/nf_lp` are imported read-only.*

Grades: **[P]** proved, **[D]** derived, **[DV]** derived + validated numerically, **[NE]** numerical evidence, **[M]** model, **[I]** open.
Conventions as `LARGEJ_OFFDIAG.md`: `Λ = √((Δ_g−2)(Δ_g−1))`, `b̂ = (2J+1)Λ/m`, `β = arccosh(1 + (2J+1)²/(2m²))`, `t = νβ`, `Γ_off` in l9 units
(`F₊₋` density per unit `ψ₂`-weight `= eo2·Γ_off/m`).

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
