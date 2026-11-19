# Route A′ with the `(2,4)` structure `EE = V₀⁴(z₂′·z₄′)²`: crystal (flat) limit, the real zeros of its forward normalisation, and the bounds

*19 Nov 2026, EE-functional agent. Unrefereed, not committed. Nothing in `RESULT_D.md`, `STATUS.md`, `H0_SPINNING_FUNCTIONAL.md`,
`H7_I14_LIGHTSIDE.md`, `OFFDIAG_KERNEL.md` or any other existing document was edited. New code: `ads4/code/spinfunc/u1_crystal_EE.py`,
`u2_hclosed_struct.py`, `u3_chi_EE.py`, `u4_combo_lp.py` (+ `u4_queue.sh`), `u5_ee_points.py`, `u6_ee_neargap.py`; caches
`spinfunc/cache/u*`. Project modules (`g1`, `g4`, `v3_regge_kernel`, `v4_hclosed`, `v4_chi`, `s1`, `s2_lp`, `l9_points`, `q3_analysis`, …)
are imported read-only (where a cache name was hard-wired it is redirected by monkey-patching inside the process).*

Grades: **[P]** proved (here: exact symbolic computation), **[P given X]**, **[D]** derived, **[DV]** derived + checked numerically,
**[NE]** numerical evidence, **[M]** model, **[I]** open. Conventions as `OFFDIAG_KERNEL.md`.

## 0. Summary

1. **Blocker (ii) — flat limit: resolved, positively** [P for the crystal identity + DV]. The covariant form of `EE` is
   `EE = [V_{0,42}²(H₂₄ − 2V_{2,04}V_{4,02})/(X₂₀X₄₀)]²` (reduces to `V₀⁴(z₂′·z₄′)²` at `x₀ = 0`). At the `D = 4` scattering crystal,
   `V_{2,04} ∝ H₂₀` and `V_{4,02} ∝ H₄₀` on the whole polarisation plane with `V_{2,04}V_{4,02}X₂₀X₄₀/(X₂₄H₂₀H₄₀) = −4` and `λ = 1` for **all**
   `x` (sympy, `x` symbolic). Hence the photon factor of `EE` **equals `q` identically at the crystal**, so **`EE ≡ q²` at the bulk point**:
   helicity-flip-free, opposite/same channel ratio exactly `(1−x²)⁻⁴`, same normalisation (ratio `1.00000000`), both `S⁰` points. NF02's flat
   certificate applies unchanged. No admixture is needed. (`EH = V₀⁴e₂₄Ĥ₂₄` is identically **null** at the crystal.)
2. **Blocker (i) — real zeros of `h₋ + h₊`: not fatal** [P given H7 §18.1 + DV + NE]. Exact `h∓^{EE}` (v4_hclosed method): `h₋`, `h₊`
   separately have no real zero; only their sum vanishes (`ν = 1.171513`, `6.642858`). Moreover **`h₋^{EE} = g·c(ν+i)`, `h₊^{EE} = g̃·c(ν−i)`
   with `c(ν) = 32(4ν²+1)(4ν²+9)/((4ν²+25)(4ν²+81))`** (even, positive on ℝ): `EE`'s leading-Regge kernel is CHMRS's scalar kernel times
   a positive ν-only normalisation, `Ĉ′_EE = c(ν)⁻¹𝒜_g⁻¹[Ψ_EE]`. The zeros are an artefact of forward normalisation. The Abel equation
   written in the `χ`-normalisation (`X = χ(h₋+h₊)` is *not* the unknown; `Y = χ(|h₋|+|h₊|)` is) has a bounded operator; `X` vanishes
   automatically (linearly) at the two zeros. Solvable, residual `≤ 3·10⁻¹⁵`, **null space of dimension 2** (index `w = 5`; SVD: two
   singular values `≤ 6·10⁻¹⁵`, next `0.61`), graviton-pole identity to `≤ 1.1·10⁻¹⁴`, grid-converged. Option (b) (`EE + λq²`) is worse,
   not better: no `λ` with `|λ| ≥ 0.001` (both signs, tested to `±0.001`) gives a feasible band LP in model M′ (§4.2).
3. **Best candidate: `Ω^{EE} := ∫ψ₂ Ĉ′_EE + ∫ψ₃Ĉ₃` (NF02 weights, `Ĉ₃` unchanged)** [M for the heavy side]. Leading-Regge helicity
   off-diagonal `≡ 0` [P]; NF02 per-state positive on the **whole twist-allowed band** (`β ≤ 2.2`, `b̂ ≤ 2Λ`): min `λ/N_s = 0.2263 / 0.2056 /
   0.2056` at `Δ_g = 100 / 300 / 1000` (253 / 319 / 374 states); the LP returns NF02 itself: **`A = 3022.34`, `B_AdS = −6552.66`**, bounds
   **`2.21·10⁻⁴ / 3.13·10⁻⁶ / 2.89·10⁻⁸`** (nominal convention; `r3.bound`: `2.18·10⁻⁴ / 3.12·10⁻⁶ / 2.89·10⁻⁸`), times the `δk(0)` factor
   `(1+δk)^{−½} = 1.0036 / 1.0004 / 1.00003` (`δk(0)_EE = −7.245·10⁻³ / −7.666·10⁻⁴ / −6.830·10⁻⁵`, h1 convention; `q²`: `−7.200·10⁻³ /
   −7.661·10⁻⁴ / −6.830·10⁻⁵`), i.e. `2.22·10⁻⁴` at 100 in h1's convention, inside the H17 budget. `Ĉ₃|_E = 0` untouched.
4. **Near gap beyond leading Regge (§5)** [NE + M]: the `K₂` crystal rows of the S1_S4 chain were recomputed with `EE` for the two
   states that break `n = 1` for `q²` at `Δ_g = 300` (`J = 18, 20`; rows `Δ = 336/400/800/1600`, 156 + 153 crystal points). Same chain,
   same helicity basis: **`J = 18`: `q²` −0.848 → `EE` +0.450; `J = 20`: `q²` −1.335 → `EE` +0.481** (`m ≈ 1.06–1.07Λ`); the off-diagonal
   correction drops 33× / 92× (`|dF₊₋|/√(NN)` 1.62 → 0.049, 2.15 → 0.024), the diagonal correction grows (0.33 → 0.50, 0.21 → 0.48 `N_s`).
   At `Δ_g = 1000`: +0.683 → +0.732, +0.688 → +0.937. So the S1 failure at 300 disappears for the computed states; `n = 1` at 300 is
   restored **for `J = 18, 20`**, with margin `t*_J ≈ +1.8 / +2.0` (the diagonal now limits). `J ≤ 16`, odd `J` and `J > 20` were not
   recomputed with `EE` [I].
5. **Answer.** Flat verdict: `EE` is flip-free with exactly `q²`'s stripped density (no admixture needed). Blocker (i) is resolved by the
   Abel form in the χ-normalisation (option (a)); option (b) has no feasible `λ ≠ 0`. Best unconditional candidate (within the project's
   standing hypotheses): `Ω^{EE}`, bounds `2.21·10⁻⁴ / 3.13·10⁻⁶ / 2.89·10⁻⁸` (unchanged), now with a model-level per-state positivity on the
   whole twist-allowed band at leading Regge and at the two previously failing near-gap states.

## 1. Blocker (ii): the crystal (flat bulk-point, `D = 4`) image of `EE` (`u1_crystal_EE.py` → `cache/u1_crystal_EE.txt`)

### 1.1 Covariant form [D + DV]
With CLS-D's building blocks (definitions as `g4_images.py`: `X_ij = x_ij²`, `H_ij/X_ij = z_i·I(x_ij)·z_j`,
`V_{i,jk} = [(z_i·x_ik)x_ij² − (z_i·x_ij)x_ik²]/x_jk²`), at `x₀ = 0`, `x₄ = −x₂`:
`V_{2,04}V_{4,02}/X₂₄ = −(z₂·x₂)(z₄·x₂)/x₂²` and `H₂₄/X₂₄ = z₂·z₄ − 2(z₂·x₂)(z₄·x₂)/x₂²`, so

**`e₂₄ := (H₂₄ − 2V_{2,04}V_{4,02})/X₂₄ = z₂·z₄`** (signature-independent; `= z₂′·z₄′` after the inversion referral, since `I² = 1`),

and the graviton structure is **`EE = [V_{0,42}²X₂₄e₂₄/(X₂₀X₄₀)]² = [V_{0,42}²(H₂₄ − 2V_{2,04}V_{4,02})/(X₂₀X₄₀)]²`**, whose `x₀ = 0` reduction
is `V₀⁴(z₂′·z₄′)²` with the same `V₀`-normalisation as `q = (2H₂₀H₄₀ + V_{0,42}²H₂₄)/(X₂₀X₄₀)` (v3/s1). Check (0) of u1: `e₂₄ = z₂·z₄` to
`≤ 4·10⁻¹⁵` at random complex configurations. `H₂₄` and `V_{2,04}V_{4,02}` have identical homogeneity (degree 1 in `X₂`, `X₄`, `Z₂`, `Z₄`) and
are both transverse (`Z_i → Z_i + αX_i`), so `EE` is a legitimate fixed `⟨TT𝒪⟩` (dual) structure [D]. Caveat: H0 §9.1 quotes CLS-D's
`j′ = 0` list as `H₂₀^{2−n}H₄₀^{2−n}H₂₄ⁿ` only; `EE` contains `V_{2,04}V_{4,02}`, which is outside that list away from the crystal (at the
crystal it is not, §1.2). That CLS-D's construction (commutator vanishing, Regge convergence `J₀ < 2`, double-trace zeros, Lemma R)
applies to `EE` as to `q²` is [D-sketch]: none of those arguments uses the specific `(2,4)` structure (OFFDIAG §3.5 makes the same point).

### 1.2 The crystal identity [P (symbolic) + DV]
At the scattering crystal of g4/g5 (`x₄ = −x₂ = x/(−x²)`, `x = (−ia, ia, 0)`, `x₀ = (0, 0, ∓a)`, `z₀ = (1,0,0)`, `a = (m − √(m²−ν²))/ν`,
CLS-D's polarisation plane `z = (0, −u, w)`), sympy with `x = ν/m` symbolic (u1 item (5)):
* `λ := 4V_{0,42}²X₂₄/(X₂₀X₄₀) = 1`;
* `V_{2,04}/H₂₀` and `V_{4,02}/H₄₀` are **independent of the polarisation** `(u, w)`;
* `V_{2,04}V_{4,02}X₂₀X₄₀/(X₂₄H₂₀H₄₀) = −4` for all `x`.

Hence the photon factor `E := (λ/4)(Ĥ − 2V̂V̂) = Ĥ/4 + 2P̂ = q` **identically** (`Ĥ = H₂₄/X₂₄`, `P̂ = H₂₀H₄₀/(X₂₀X₄₀)`), and `EE = E² = q²` at the
crystal. Numerically (u1 (1)–(4), `x = 0.05 … 0.95`, both `S⁰` points): on the flat bilinears `{𝒱𝒱, 𝒲𝒲, 𝒲𝒱, 𝒱𝒲}`,
`P̂ → (1,0,0,0)`, `Ĥ → (−4,−4,0,0)`, `V̂V̂ := V_{2,04}V_{4,02}/X₂₄ → (−4,0,0,0)`, `E = q → (1,−1,0,0)` (×`1/(4ν²(m²−ν²)²)`), residuals
`≤ 2·10⁻¹⁵`. Graviton helicity tables (g1 helicity vectors, CLS-D conjugate pairing):

| `x` | `EE` flips / preserved | `EE` opp/same | `(1−x²)⁻⁴` | `EE`/`q²` (same channel) | `|EE − q²|/|q²|` (all channels) | `|EH|/|q²|` |
|---|---|---|---|---|---|---|
| 0.05 | `10⁻³²` | 1.01006281 | 1.01006281 | 1.00000000 | `7·10⁻¹⁶` | `6·10⁻¹⁶` |
| 0.3 | `2·10⁻³¹` | 1.45825821 | 1.45825821 | 1.00000000 | `7·10⁻¹⁶` | `4·10⁻¹⁶` |
| 0.5 | `2·10⁻³¹` | 3.16049383 | 3.16049383 | 1.00000000 | `4·10⁻¹⁶` | `3·10⁻¹⁶` |
| 0.8 | `3·10⁻³²` | 59.53741808 | 59.53741808 | 1.00000000 | `3·10⁻¹⁶` | `8·10⁻¹⁷` |
| 0.95 | `6·10⁻³⁴` | 11065.767 | 11065.767 | 1.00000000 | `2·10⁻¹⁵` | `8·10⁻¹⁷` |

The `S⁰` point `+1` is the `x₃`-reflection of `−1` for `E` exactly (as for `q`, g5).

### 1.2b Check in the crystal-integral engine (full null-circle contraction) (`u5_ee_points.py diag` → `cache/u5_diag_scaling.txt`) [NE]
u1 evaluates the structures on CLS-D's physical polarisation plane. The l6/l9 engine contracts the `(2,4)` structure over the whole
spin-2 null circle, and at the engine's saddle `V_{2,04}` is **not** proportional to `H₂₀` off that plane (e.g. `J = 18`, `Δ = 400`,
`x = 0.1`: `V_{2,04}/H₂₀ = −0.0998i, −1.6·10⁴i, −0.1995i` on the three lightcone basis vectors). So the identity is a statement about the
bulk-point limit, not about finite `m`. The engine's leading Laplace order `M₀` (and `K₁ = M₀ + M₂`) of `Ĉ_EE` vs `Ĉ₂[q²]` on parity-even
blocks, same `(Δ, ν)`:

| `J` | `m ≈ 400` (`x = 0.15 / 0.3`) | `m ≈ 1600` | `m ≈ 6400` |
|---|---|---|---|
| 4 | `K₁`: `0.90−0.32i` / `0.95−0.15i` | `0.99−0.09i` / `1.00−0.05i` | `1.000−0.023i` / `1.000−0.012i` |
| 8 | `0.84−0.21i` / `1.01−0.16i` | `0.98−0.05i` / `1.00−0.03i` | `0.999−0.018i` / `1.000−0.007i` |
| 18 | `0.39−0.15i` (near a zero of the density) / `1.02−0.23i` | `0.89−0.19i` / `1.00−0.05i` | `0.994−0.045i` / `1.000−0.011i` |

(entry `(1,1)` of the `2×2`; the four entries agree to a few % except near density zeros). **The ratio → 1 as `m → ∞` at fixed `J`
(deviation `∝ 1/m`, mostly a phase, i.e. the imaginary `O(1/m)` term removed by `Re`, H0 §22.2)**: the flat-limit identity `EE ≡ q²` holds
in the engine too [NE]. At finite `m` the two structures differ by `O(J²/m)` in the real part — a genuine near-gap difference (§5).

### 1.3 Verdict for task 1
* **`EE` is helicity-flip-free at the bulk point and its stripped density is exactly `q²`'s** (same `(1−x²)⁻⁴` channel ratio, same
  normalisation) [P for the crystal identity, given the crystal configuration and polarisation map of H0 §17; DV]. The bulk-point image of
  `Ĉ′_EE` is therefore the stripped `B⁽¹⁾₂` with the same normalisation as `Ĉ′₂` (`q²`), and the flat limit of `∫ψ₂Ĉ′_EE + ∫ψ₃Ĉ₃` is NF02:
  **the NF02 flat certificate applies** (same grade as for `q²`: H0 §R.7 [Q + DV], "products of structures evaluate to products at the
  crystal").
* No admixture is needed. For the record: `EH` is null at the crystal (`E` vanishes on flips, `Ĥ` on preserved), so every combination
  `a·EE + b·q² + c·EH + d·Nq` is flip-free with density `(a + b)×` stripped `B₂` [P + DV].
* Interpretation [D]: at the crystal the dual-structure space degenerates (`V_{2,04} ∝ H₂₀`), which is why CLS-D's three-structure `j′ = 0`
  basis suffices there; away from it (Regge, `x₀ = 0`) `EE` and `q²` differ, and only `EE` has `e₂₄ = z₂′·z₄′`, which vanishes on the
  flipped dual tensor (OFFDIAG §3.4).

## 2. Exact `h∓` of `EE` and why the zeros of `h₋ + h₊` are harmless (`u2_hclosed_struct.py` → `cache/u2_hclosed.txt`, `.pkl`)

**Method** [P given H7 §18.1]: v4_hclosed's exact computation (equivariant gradient at one celestial point, rising-factorial ratios of the
`K(c,λ)`) with the `(2,4)` monomials of `STRUCTS['EE']`. All structures share `K^ref_A K^ref_B`, so the rational `R∓` are directly
comparable. Checks: `q²` reproduces v4's `R^T∓` exactly (difference `0`); the forward ratio `(R₋+R₊)_S/(R₋+R₊)_{q²}` equals s1's `img_vec`
(β = 0) to 8 digits for `EE`, `EH`, `Nq` at `ν = 1, 1.17, 3, 6.64, 20` [DV]. (`h∓ = N_com R∓`, `N_com` even, positive, zero/pole-free on
`|Im ν| < 11/2`.)

* **`R₋^{EE} = −2(2ν−i)(2ν+i)²(2ν+3i)(2ν+5i)(2ν+11i) / [ν(2ν−7i)(2ν−3i)(2ν+7i)]`, `R₊^{EE}(ν) = R₋^{EE}(−ν)`.**
* **Scalar structure [P given §18.1]:** `R₋^{EE}/R₋^{sc} = c(ν+i)`, `R₊^{EE}/R₊^{sc} = c(ν−i)` exactly, with
  **`c(ν) = 32(4ν²+1)(4ν²+9)/((4ν²+25)(4ν²+81))`** (even; positive on ℝ; poles `±5i/2, ±9i/2` outside the strip; zeros `±i/2, ±3i/2`).
  So `h₋^{EE} = g·c(ν+i)`, `h₊^{EE} = g̃·c(ν−i)` (CHMRS `g`), i.e. **`Ψ_EE`'s leading-Regge kernel is `𝒜_g[c(ν)P_l]`**, (HN) holds with the
  scalar de-shift and the normalisation `N̂ = c(ν)`, and `𝒜_h = 𝒜_g ∘ c`: **`Ĉ′_EE = c(ν)⁻¹ 𝒜_g⁻¹[Ψ_EE]`**. (`q²` does not have this form:
  `a₋(μ−i)/a₊(μ+i) ≠ 1`, H7 §19.2. `EH` has it too, with `c_EH = 3840(4ν²+1)(4ν²+9)/((4ν²+25)(4ν²+81)²)`.)
* Zeros/poles (`|Im ν| ≤ 2`): `h₋^{EE}`: zeros `±i/2` (`−i/2` double), `−3i/2`; pole `0`, `3i/2`. `h₊^{EE}`: mirror. **No real zero of `h₋`
  or `h₊` on `ν > 0`.** `h₋ + h₊` (= `c(ν+i)g + c(ν−i)g̃`): zeros `±i/2` and the **real** `±1.171513`, `±6.642858` — the zero of the scalar
  `g + g̃` at `5.979` reshuffled by the positive factor `c(ν±i)`. `EH`: real zeros `1.10717`, `5.68104`.
* `ρ = h₊/h₋`: **`ρ(ν)ρ(−ν) = 1`** (unlike `q²`), winding `w = Z₊ − P₊ = 5` (numeric `5.0000`) → predicted null dimension `(w−1)/2 = 2`
  (`q²`: `w = 7`, 3). `u_EE = i s/ν`: `s(0) = 165/646`, `s(∞) = 23/2` (same large-ν symbol as `q²`).

## 3. The non-local weight of `Ĉ′_EE` (`u3_chi_EE.py` → `cache/u3_chi_EE_w.txt`, `u3_null_EE.txt`, `u3_chi_<S>_<Δg>.pkl`)

**Formulation [D].** The Abel equation of H7 §23.2, `cosh t·C[X] − i sinh t·S[uX] = C[ψ₂]`, involves only `X = χ(h₋+h₊)` and
`uX = χ(h₋−h₊)`; neither requires dividing by `h₋+h₊`. With `n := (|h₋|+|h₊|)·e^{i arg(h₋+h₊)_{q²}}` (positive modulus, smooth zero-free phase),
`Y := χn`, `a := (h₋+h₊)/n`, `b := (h₋−h₊)/n` (`|a|, |b| ≤ 1`):

**`A[Y] := C[aY] − i tanh t·S[bY] = C[ψ₂]/cosh t`, `X = aY`.**

The operator is bounded; `X` vanishes wherever `h₋+h₊` does, with `Y` (equivalently `χ`) regular. This is exactly "require `X` to vanish
there": it is automatic, not a constraint. The functional pairs `Re∫XΨ̂ = Re∫χΨ_EE` and never sees the forward normalisation. Lemma N
(H7 §24.1) holds verbatim (it is formulated in `f = χh₋`, `f(−ν) = χh₊`).

**Results [NE]** (min-norm correction `Y − ψ₂`, i.e. the χ-normalisation gauge of H7 §24.4; same code run on `q²` as control):

| | `Δ_g = 100` | 300 | 1000 |
|---|---|---|---|
| residual (V2), `EE` | `2.5·10⁻¹⁵` | `1.6·10⁻¹⁵` | `2.9·10⁻¹⁶` |
| **V3 graviton pole** `∫Y[h₋D(ν+i)+h₊D(ν−i)]/n` vs `∫ψ₂D`, `EE` (no inversion) | **`3.6·10⁻¹⁵`** (0.97) | **`1.1·10⁻¹⁴`** (0.82) | **`3.5·10⁻¹⁵`** (0.67) |
| `δk(0) = −2∫(Re X − ψ₂)p⁶dp`, `EE` | `−7.2452·10⁻³` | `−7.6663·10⁻⁴` | `−6.8302·10⁻⁵` |
| same code, `q²` (v4_chi: `−7.2007·10⁻³ / −7.661·10⁻⁴ / −6.830·10⁻⁵`) | `−7.2005·10⁻³` | `−7.6611·10⁻⁴` | `−6.8297·10⁻⁵` |
| `max|X−ψ₂|/max ψ₂` on `[0,Λ]` (on `ν ≥ 20`), `EE` | 0.56 (0.031) | 0.22 (0.0074) | 0.068 (0.0018) |
| refinement `dν/2`: `max|ΔX|/max ψ₂` (`ν ≥ 20`); `δk(0)` | `1.9·10⁻⁵` (`6·10⁻⁹`); same 6 digits | `2.6·10⁻⁵` (`5·10⁻⁹`); same | `3.2·10⁻⁵` (`2·10⁻⁹`); same |
| interior `(X−ψ₂)/ψ₂` at `p = 0.4` vs first-order symbol `½ψ″ − (σ̃ψ)′` | `+1.346·10⁻²` / `+1.319·10⁻²` | `+1.369·10⁻³` / `+1.366·10⁻³` | `+1.214·10⁻⁴` / `+1.214·10⁻⁴` |

* **Null space [NE + index].** Dense SVD of `A` (depends on `h∓` only): `dν = 0.05, N = 2048`: smallest `8·10⁻¹⁵, 7·10⁻¹⁴`, then `0.86`;
  `dν = 0.025, N = 4096`: `4.5·10⁻¹⁶, 5.7·10⁻¹⁵`, then `0.61`. **Exactly two null directions**, as `(w−1)/2 = 2` predicts; peaks at
  `ν = 0.91` and `6.56`, supported on `ν ≲ 12.4` and `≲ 17.9` (`|Y_h(20)| ≤ 1.3·10⁻⁴ max`). Control (same code, `q²`): three small
  singular values (`6·10⁻¹⁵, 5·10⁻¹⁴, 1.1·10⁻⁷`), as v4. By Lemma N the two null directions are invisible to the leading-Regge action, the
  graviton pole, `|ĝ₃|²` and the flat limit [P given NE-observed regularity]; the gauge is a convention exactly as for `q²` (one fewer
  direction). Since `ρ(ν)ρ(−ν) = 1` for `EE`, the CHMRS class (`χ` even, analytic on the closed strip) is well defined here (it is not for
  `q²`, H7 §24.4); it imposes `F(−i/2) = F′(−i/2) = 0` (double zero of `h₋` at `−i/2`), 2 complex conditions for 2 parameters, so a
  canonical gauge may exist — **not tested [I]**.
* **Near the zeros** (`Δ_g = 300`, units `max ψ₂`): `|a| = 2.5·10⁻³`, `X = +1.0·10⁻⁴` at `ν = 1.170`; `|a| = 1.4·10⁻³`, `X = −3.4·10⁻⁴` at
  `ν = 6.650`; `X` changes sign linearly through each zero while `Y` stays smooth (`−0.033`, `0.236`). The weight is finite and the
  functional is well defined there.

**Answer to blocker (i), option (a) [D + NE]: the real zeros are not fatal.** The Abel-form weight is well defined (bounded operator,
`X` vanishes at the zeros automatically), solvable (consistent RHS, residual `10⁻¹⁵`), has a 2-dimensional null space with the same Lemma-N
status as `q²`'s 3-dimensional one, and satisfies the graviton-pole identity to `10⁻¹⁴`. What is singular at the zeros is only the
*forward-normalised local* de-shift `(1 ± u_EE)/2` (`|u_EE| → ∞`); any chain that uses it must stay away from `ν ≈ 1.17, 6.64` (the
`Δ_g = 300, 1000` chains do: all nodes have `ν ≥ 15`; at `Δ_g = 100` the `x = 0.05` node sits at `ν ≈ 4.9–6.5`, where the local form was
already invalid for `q²`, H7 §20.3).

## 4. Per-state positivity at leading Regge on the whole band, and the LP (`u4_combo_lp.py`, `u4_queue.sh`, `u4_queue2.sh`; `cache/u4_*`)

Model M′ of OFFDIAG §2.2 (diagonal = stripped flat `B₂ + B₃` + leading AdS curvature `2∫ψ₂(P_l − J₀) + uσ∫ψ₃(P_l − J₀)`, `Ĉ₃` off-diagonal 0;
r5's state grid `m/Λ ∈ [1.02, 5]`, `b̂ ∈ [10, 2Λ]`, twist-allowed). For `EE` the diagonal is the same as for `q²` (same Regge diagonal `P_l`
by construction of the de-shift, same flat density by §1) and the off-diagonal is `0` at every `(ν, β)` [P].

### 4.1 `EE` alone [M]

| `Δ_g` | states (`β ≤ 2.2`) | NF02 min `λ/N_s` [state] | LP (NF02 family, `δ = 0.2`, `d_w = 0`) |
|---|---|---|---|
| 100 | 253 | **+0.2263** (`m = 1.05Λ`, `b̂ ≈ 10`) | NF02: `A = 3022.34`, `B_AdS = −6552.66` |
| 300 | 319 | **+0.2056** (`1.02Λ`) | NF02: same |
| 1000 | 374 | **+0.2056** (`1.02Λ`) | NF02: same |

No state of the band is negative; the large-`b̂` states sit at `λ/N_s ≈ 0.93–1.00` (curvature-dominated diagonal, nothing off-diagonal).
The binding states are the near-gap flat corner, as in NF02's certificate.

### 4.2 Option (b): `S(λ) = EE + λq²` [NE + M]
* **Real zeros** (`u4_combo_lp.py zeros`, `cache/u4_zeros.txt`) [P + NE]: `h₋+h₊` of `S(λ)` vanishes on `ν > 0` iff `r(ν) := n_EE/n_q = −λ`;
  `Im r = 0` only at `ν = 0.926, 1.172, 2.344, 6.643, 6.684` with `Re r = −0.285, 0, −0.678, 0, −0.023`; `r(0⁺) = 0.270`, `r(∞) = 1`.
  So **`S(λ)` has no real zero for `λ ∈ (−1, 0)` (except `λ = −0.270`, a zero at the endpoint `ν = 0`) and for `λ ∈ (0, 0.0228)`**; for
  `λ > 0` zeros reappear at `λ = 0.0228, 0.285, 0.678`. But it is badly conditioned: `min_ν |n_EE + λn_q|/(|n_EE| + |λ n_q|) = 0.048 / 0.081 /
  0.135` at `λ = −0.02 / −0.05 / −0.1` (at `ν ≈ 6.5–6.6`), `0.003 / 0.030 / 0.084` at `λ = +0.02 / 0.05 / 0.1`.
* **The off-diagonal is not `λ`-suppressed where `n_EE` vanishes.** `o_EE ≡ 0`, so `Y_S = Y_q·λ/(r + λ)`: `≈ Y_q·λ/(1+λ)` at large `ν`, but
  `= O(Y_q)` in a window of width `O(λ)` around `ν = 1.17, 6.64`, where the forward normalisation of `S` is carried by `λq²` alone.
* **LP on the whole band** (OFF=psi, the s3 protocol): **infeasible for every `λ ∈ {±0.001, ±0.002, ±0.005, ±0.01, −0.02}` at both
  `Δ_g = 300` and `1000`** (`cache/u4_queue*.log`). NF02 per state fails at `β ≈ 0.9–1.2` (`b̂ ≈ 280–400` at 300; `≈ 330–620` at 1000),
  erratically in `b̂` (e.g. `λ = 0.001`, 300: `−0.04 … −3.1`; `λ = −0.001`: `−0.9 … −15`): the narrow `O(1)` spike of `Y_S` near `ν ≈ 6.6`
  against the small large-`b̂` curvature diagonal. In the exact functional this spike is multiplied by the non-local weight `X_S`, which
  (like `X_EE`) may dip there; that was not computed. **No feasible `λ ≠ 0` was found; `λ = 0` (EE itself) is the only clean member** and
  needs neither a zero-free normalisation nor the spike [NE + M].

## 5. Near-gap corrections with `EE`: the `J = 18, 20` crystal rows (`u5_ee_points.py`, `u6_ee_neargap.py` → `cache/u5_pts_EE_J*_k2*.pkl`, `cache/u6_ee_neargap.txt`)

**What was computed [NE].** The l9 `K₂` crystal engine (l6 integrand, `L1_DEG = 6`) with the `(2,4)` photon factor `q` replaced by
`E = V_{0,42}²(H₂₄ − 2V_{2,04}V_{4,02})/(X₂₀X₄₀)` (`D = E⊗E`), at exactly the `(Δ, ν)` keys of the q² rows of S1_S4 (`q2_pts_C2_J*_k2*`):
`J = 18`: rows `Δ = 336, 400` (19 `x` each), `800, 1600` (16 each); `J = 20`: `400` (19), `800, 1600` (16); each with its `Δ_ref = 2·10⁴`
reference at the same `ν` (~21–38 s per point, one core). `Ĉ₃` rows (`C3H`, `o_H q`) unchanged. Then the S1_S4 chain (`q3_analysis`:
`K₂` linear chain, aligned off-diagonals, `c_off = mΓ_off`, combined near-gap states, `Q3_C3SHIFT = none`, rows `all`) run twice, with the
`q²` rows and with the `EE` rows; local `k = 2` de-shift `(1 ± u)/2` with `u_T` resp. `u_EE` (all nodes at `ν ≥ 15`, away from the zeros
`1.17, 6.64`); helicity basis of each row from the `q²` records in both runs (the rank-1 basis of `M₀` at `x ≥ 0.15`; `M₀^{EE} → M₀^{q²}` as
`m → ∞`, §1.2b).

| `Δ_g` | `J` | min `λ/N_s` (uncorrected) | `q²`: min `λ/N_s` / `max|dF₊₋|/√(NN)` / `max|dF_hh|/N_h` | **`EE`**: same | `EE` `t*_J` |
|---|---|---|---|---|---|
| 300 | 18 | +0.892 | **−0.848** / 1.620 / 0.334 (`1.06Λ`) | **+0.450** / 0.049 / 0.497 (`1.06Λ`) | `(−∞, +1.8)` |
| 300 | 20 | +0.924 | **−1.335** / 2.152 / 0.206 (`1.07Λ`) | **+0.481** / 0.024 / 0.481 (`1.07Λ`) | `(−∞, +2.0)` |
| 1000 | 18 | +0.751 | +0.683 / 0.142 / 0.016 | +0.732 / 0.0018 / 0.020 | `(−∞, +38)` |
| 1000 | 20 | +0.840 | +0.688 / 0.208 / 0.039 | +0.937 / 0.0027 / 0.130 | `(−6.5, +98)` |

(`q²` column = S1_S4's chain, reproduced: their `J = 18, 20` range `−0.4 … −1.3`.)

**Reading.**
1. **The off-diagonal correction collapses** (33× / 92× at 300, ~80× at 1000), as the leading-Regge kernel predicts (`o_EE ≡ 0`): what is
   left (`≤ 0.05 N_s`) is beyond-Regge / the per-row basis shear. This confirms OFFDIAG §1's diagnosis that the S1 failure is the
   leading-Regge `q²` off-diagonal [NE].
2. **The diagonal correction is larger for `EE`** (`≈ 0.5 N_s` at `J = 18, 20`, `Δ_g = 300`, vs `0.2–0.33` for `q²`): `EE` and `q²` agree at the
   bulk point but differ at `O(J²/m)` (§1.2b), and that difference shows up in the diagonal. It is negative-definite-ish here (`t*_J` has a
   finite positive end, `+1.8`/`+2.0`), so the margin is a factor ≈ 2 on the computed correction. Its growth with `J` (the `J²/m` scaling)
   is the obvious next risk for `J ≳ 24` near the gap at `Δ_g = 300` [I].
3. Not recomputed with `EE` [I]: even `J ≤ 16` (the `q²` corrections there were harmless, `J = 12, 14` at 0.46 / 0.23 and `J = 16`
   marginal), odd `J` (parity-odd sector, `p3_odd_action` also uses `q`), `J > 20`, the direct (RQMC) rows, and `Δ_g = 100`. The binding
   states (`+−`, `J = 5`, 0.2054 / 0.2021) are parity-odd and not recomputed; their `q²` correction was `≤ 2·10⁻⁴`.

## 6. Light side and bounds for `Ω^{EE} := ∫ψ₂Ĉ′_EE + ∫ψ₃Ĉ₃` (`Ĉ′_EE := 𝒜_h⁻¹[Ψ_{2,ν}[EE]]`, `Ĉ₃` unchanged)

| item | value | grade |
|---|---|---|
| graviton pole `−Ĉ′_EE|_E = 16πG/(ν²+9/4)` | by its own de-shift ((HN) exact by construction: `h₀ = 0`, helicity-blind, H7 §18.3 holds for every weight-0 tensor), as for `q²` (H7 §20.1); smeared form V3 to `≤ 1.1·10⁻¹⁴` at 100/300/1000 | [P given Lemma R, Q3, Q4, D2, §18.1 kernel (D)] + [DV] |
| `Ĉ₃|_E = 0` | unchanged (`Ĉ₃` not touched) | as before |
| `A`, `B_AdS` | **`3022.34`, `−6552.66`** (NF02; `L_A` uses `ψ₂`) | [P given the pole] |
| flat limit | NF02 (crystal image of `EE` ≡ `q²`, §1), certificate applies | [P (crystal identity) + Q + DV] |
| `k(0)` | `1` (flat) `+ δk(0)`, `δk(0)_EE = −7.245·10⁻³ / −7.666·10⁻⁴ / −6.830·10⁻⁵` (h1 convention; = `q²`'s to `0.6 % / 0.07 % / 0.01 %`) | [NE]; beyond that H17 as for `q²` |
| bounds `√(t₄²+α_T²)` | **`2.21·10⁻⁴ / 3.13·10⁻⁶ / 2.89·10⁻⁸`** at `Δ_g = 100/300/1000` (nominal), × `1.0036 / 1.0004 / 1.00003` in h1's convention (`2.22·10⁻⁴` at 100) | [M] for the heavy side (§4–§5) |
| null space of the weight | 2-dim, invisible to all `Â_h`-images of strip-analytic densities (Lemma N) | [P given NE regularity] |

## 7. Graded summary

| statement | grade |
|---|---|
| covariant form `EE = [V_{0,42}²(H₂₄ − 2V_{2,04}V_{4,02})/(X₂₀X₄₀)]²` (reduces to `V₀⁴(z₂′·z₄′)²`) | [D] + [DV `4·10⁻¹⁵`] |
| at the `D = 4` crystal `V_{2,04} ∝ H₂₀`, `V_{4,02} ∝ H₄₀` on the polarisation plane, ratio `−4`, `λ = 1` ⇒ `E = q`, `EE = q²` for all `x` | [P] (sympy, `x` symbolic) + [DV] |
| `EE` flip-free at the bulk point, stripped `(1−x²)⁻⁴`, same normalisation as `q²` ⇒ flat limit NF02, certificate applies | [P given H0 §17's crystal + polarisation map] |
| the same in the crystal-integral engine (full null-circle contraction) as `m → ∞` at fixed `J`; finite-`m` difference `O(1/m)` (imaginary) + `O(J²/m)` | [NE] |
| `EE` is an admissible CLS-D dual structure (convergence, double-trace zeros, Lemma R as for `q²`) | [D-sketch] |
| exact `h∓^{EE}`; `h₋ = g c(ν+i)`, `h₊ = g̃ c(ν−i)`, `c = 32(4ν²+1)(4ν²+9)/((4ν²+25)(4ν²+81))`: (HN) with the scalar de-shift and `N̂ = c` | [P given H7 §18.1 (D)] |
| `h₋`, `h₊` zero-free on `ν > 0`; `h₋+h₊` real zeros `1.171513`, `6.642858`; `w = 5` | [P] |
| Abel weight in the χ-normalisation: bounded, solvable, `X` vanishes at the zeros, 2-dim null space, V3 to `10⁻¹⁴`, grid-converged | [D] + [NE] |
| null directions invisible to leading Regge / graviton pole / `|ĝ₃|²` / flat limit (Lemma N) | [P given NE regularity] |
| canonical (CHMRS-class) gauge for `EE` | [I] (class exists since `ρ(ν)ρ(−ν) = 1`; not tested) |
| `EE + λq²`: no real zero for `λ ∈ (−1,0)∖{−0.27}` or `(0, 0.0228)`; band LP infeasible for all tested `λ ≠ 0` (`|λ| ≥ 0.001`) | [P] (zeros) + [M] (LP, OFF=psi) |
| `Ω^{EE}`: leading-Regge off-diagonal `≡ 0`; NF02 positive on the whole twist-allowed band at 100/300/1000 (min `0.2263/0.2056/0.2056`) | [P] + [M] (model M′) |
| `A = 3022.34`, `B_AdS = −6552.66`, bounds `2.21·10⁻⁴ / 3.13·10⁻⁶ / 2.89·10⁻⁸` (× `1.0036/1.0004/1.00003` from `δk(0)`) | [P given the pole] for `A`, `B_AdS`; [M] for the bound (heavy side) |
| near gap, `K₂` chain with `EE` rows: `J = 18, 20` at `Δ_g = 300` positive (+0.450 / +0.481; `q²`: −0.848 / −1.335); off-diagonal correction ÷33 / ÷92; diagonal correction ≈ 0.5 `N_s` (margin ×1.8–2.0) | [NE + M] (S1_S4 chain model M; `K₂` truncation, no direct rows) |
| `EE` near-gap corrections for `J ≤ 16`, odd `J`, `J > 20`, `Δ_g = 100` | [I] (not computed) |

## 8. Hand-over

**Files** (`ads4/code/spinfunc/`, new; project modules imported read-only, cache names redirected by in-process monkey-patching):
* `u1_crystal_EE.py` → `cache/u1_crystal_EE.txt`: covariant reduction check, crystal decomposition on `{𝒱𝒱, 𝒲𝒲, 𝒲𝒱, 𝒱𝒲}`, graviton helicity
  tables of `q²`, `EE`, `EH`, `S⁰` reflection, sympy identity (seconds).
* `u2_hclosed_struct.py` → `cache/u2_hclosed.txt`, `u2_hclosed.pkl` (srepr of `R∓` for `q2`, `EE`, `EH`, `Nq`; `numfuncs()` gives numpy callables)
  (~10 s).
* `u3_chi_EE.py [Δg …]` → `cache/u3_chi_EE_w.txt` (weights; run with `U3_SKIPN=1`), `u3_null_EE.txt`/`.pkl` (SVD), `u3_chi_{EE,q2}_<Δg>.pkl`
  (`nu`, `X`, `Y`, `phi`, `dk0`, `a`, `b`). `SolverAB` is the χ-normalised Abel solver for any `(R₋, R₊)` (`Rfun(name, lam)`).
* `u4_combo_lp.py zeros|check|lp` (env `U4_DG`, `U4_LAM`, `U4_BMAX`) + `u4_queue.sh`, `u4_queue2.sh` → `cache/u4_zeros.txt`,
  `u4_{check,lp}_<Δg>_b2.2_<λs>.txt/.pkl`, `u4_queue*.log`.
* `u5_ee_points.py J "D1,D2"` (env `U5_TAG`) → `cache/u5_pts_EE_J{18,20}_k2*.pkl` (EE crystal records, same keys as `q2_pts_C2_*`);
  `u5_ee_points.py diag` → `cache/u5_diag_scaling.txt`. (`u5_pts_EE_J18_k2.pkl` is one timing record of the same engine.)
* `u6_ee_neargap.py` (env `U6_JEE`, `Q3_JEVEN`, `Q3_JODD`) → `cache/u6_ee_neargap.txt` (S1_S4 chain, `q²` vs `EE` rows side by side).

**Next steps (value order).**
1. **`EE` rows for the rest of the near-gap set at `Δ_g = 300`**: even `J = 4 … 16` (rows as `q2_pts_C2_J*` / l8 records; ~10–12 min per
   row per core), odd `J = 5 … 13` (patch `p3_odd_action`'s `q` the same way as `u5.build_spin_EE`), then `J = 24 … 32` near the gap with
   the hybrid base (S1_S4 §2.6). The diagonal correction (`≈ 0.5 N_s` at `J = 18, 20`, growing like `J²/m`) is the quantity to watch.
2. **Direct (RQMC) cross-check** of one `EE` point at the failing `q²` state (`q4_direct.py` with the `E` structure), as S1_S4 §2.3 did.
3. **`Δ_g = 100`**: the `x = 0.05` node sits at `ν ≈ 4.9–6.5`, next to the real zero `6.64` of the `EE` forward normalisation: the local
   de-shift is unusable there (it was already invalid for `q²`); needs `Re∫XΨ̂` with the non-local `X` of u3 and contour data.
4. **Canonical gauge for `EE`**: test the CHMRS-class conditions `F(−i/2) = F′(−i/2) = 0` on the two null vectors (v5_null machinery);
   if they fix the gauge, `Ω^{EE}` has a canonical weight (an advantage over `q²`).
5. For whoever edits `RESULT_D.md`/H0 §R (not edited here): route A′ could be redefined with `Ĉ′_EE` in place of `Ĉ′₂`; everything on the
   light side and the flat side carries over unchanged (§6), the leading-Regge large-`J` obstruction of OFFDIAG disappears (§4), and the
   S1 near-gap failure at `Δ_g = 300` is removed for the computed `J = 18, 20` (§5). The covariant `EE` lies outside CLS-D's printed
   three-structure `j′ = 0` list (§1.1 caveat); its admissibility as a dual structure should be stated as an explicit [D-sketch] input.
