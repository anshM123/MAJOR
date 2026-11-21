# Closing the parity-odd items: the origin of `Γ_OO ≡ Γ_++`, and the `+−` J = 5 state of Result Q

*21 Nov 2026, parity-odd-close agent. Unrefereed, not committed. No existing document or code file was edited; nothing of `bj_*` was
touched or run. New code `ads4/code/spinfunc/po_*.py` (+ `po_queue.sh`); caches `spinfunc/cache/po_*` (git-ignored). Read-only
imports: `p1_tto_odd`, `p2_shadow_odd`, `p3_odd_action` (via `p8_odd_pw`'s point-wise copy), `w0_ee_common`, `l5_spin`, `l8_analysis`,
`l7_dc2_analysis`, `bg_ser`, `bg_geom`, `bg_ee`, `bg_eeq`, `bg_qan`, `bf_qan`, `bf_core`, `bh_ffscan` (copied), `resultq/stage6_*`.
≤ 2 worker processes, ≤ 0.6 GB each.*

Grades: **[P]** proved; **[D]** derived (an argument with a stated input); **[DV]** derived and checked numerically to round-off;
**[NE]** numerical evidence; **[M]** model; **[I]** open.

## 0. Summary

| item | result | grade |
|---|---|---|
| **Goal 1: origin of `Γ_OO ≡ Γ_++`** | **Helicity lemma** (§1.2): if `T₁`'s null polarisation satisfies `V_{1,25} = 0` (a helicity eigen-polarisation in the 3-point frame), the symmetric parity-odd `⟨TTO⟩` equals `±i ×` a fixed combination of the even structures, as a function of `(z₂, z₅)`; the combination is the `++` (same-helicity) structure and the sign is the chirality. **The functional's crystal saddle has exactly `V_{1,25}(x₁₀) = 0`** (`x₀ = (0, 0, n₀a)`; `n₀` = chirality). Hence at order 0 the OO entry is `κ ×` the `T0`-basis `++` diagonal exactly, and the per-row constants are `κ(Δ)/κ(Δ_ref)`, `κ = (w·t₊₊)(v·t₊₊)`: computed **`2.2036 / 2.0917 / 1.7439 / 1.7755 / 1.1408`** = HYPOTHESES_CLOSE2's `2.20 / 2.09 / 1.74 / 1.77 / 1.14` (a basis-normalisation constant, not dynamics). Beyond the saddle the lemma fails point-wise (`V₁ ≠ 0`), but the integrated statement — no `++ ↔ −−` helicity flip, `X = 0` — holds to `10⁻⁸` (`Ĉ'_EE`) / `4·10⁻⁵` (`Ĉ₃`) relative, `≤ 0.05 / 0.12 c₂`, on all 11 K₂ rows (`J = 4–22`), with a residual falling like `Δ⁻⁴`: **`Γ_OO = Γ_++` exactly at order 0 and up to `O(m⁻²) c₂` beyond** (not an exact identity). Odd `J`: same lemma with the antisymmetric even pair; the parity-odd `+−` block is a `+−` helicity diagonal. | lemma **[D]** (+ **[DV]** to `10⁻⁸–10⁻¹²`); saddle `V₁ = 0` **[P]**; order-0 identity **[P given lemma]**; per-row constants **[DV]**; `O(m⁻²) c₂` residual **[NE]**; its mechanism **[I]** |
| **Goal 2: Result Q, `+−` J = 5 (parity-odd)** | Parity-odd prefactor built on the analytic-radial engine (`po_ee_odd.py`), point-wise = exact odd integrand to `3.6·10⁻² / 1.9·10⁻³ / 5.4·10⁻⁵` (`a⁰/a¹/a²`), integrated k = 2 image = direct n3 points within 1.5σ. **`s_d` (H0) = `−0.539 ± 0.033` (gap, Δ_g = 300), `−0.512 ± 0.004` (gap, Δ_g = 1000), `−0.531 ± 0.026` (1.19Λ, Δ_g = 300)**; gauges H2/H3/H4/SB/XP within `0.016`. **Favourable sign**; same per-harmonic profile as the even `+−` J = 4 (`σ_{+−} = −1` universal correction). FF/N at the gap: `+0.144 / +0.153` (δ = 0.05), `+0.284 / +0.291` (δ = 0.2) (Δ_g = 300 / 1000); channel minimum (1.19–1.21Λ) `+0.065` (δ = 0.05), `+0.218` (δ = 0.2). | **[NE + M]** |
| **Result Q full-set minimum** | **unchanged: `+0.0061 ± 0.0013`** (δ = 0.05; SB `+0.0048`) and **`+0.1427 ± 0.0017`** (δ = 0.2; SB `+0.1410`), `++` J = 0 at the gap, both Δ_g, all gauges; the last unmeasured near-binding channel (`+−` J = 5) is now measured and moves up. | **[NE + M]** |
| **Item 3** | `J = 6` with direct small-`ν` nodes: EE block `+0.3075 ± 0.0001`, OO `+0.3078 ± 0.0001` (convention A `+0.307`). `J = 15` (odd, rows 336/400, row injection): near-gap min `+0.829 ± 0.016` (convention A `+0.648`). Even `J = 8, 10` and odd `J = 7–13` not done. | **[NE + M]** |


---

## 1. Why `Γ_OO ≡ Γ_++`: a helicity lemma for `⟨TTO⟩` and the crystal geometry

### 1.1 The precise statement that holds

Notation (HYPOTHESES_CLOSE2 §1): at even `J` the functional acts on the 3×3 structure space `[E₁, E₂, O]` (`l5_spin.conserved`
even pair, `p1` odd structure = CHL's `X^O_same`); `M^{n₀}` is the action at one point `n₀ = ±1` of `S⁰` (orders `M0, M2, M4` after the
shadow inverse); the physical action is `Σ_{n₀} M^{n₀}`, in which the `EO` entries cancel (PARITY_ODD_ROUTEA §3). `T0(Δ)` is the l8
order-0 helicity basis of the even 2×2 block, `t₊₊*` the dual (`T0⁻¹`) row of the `++` channel.

| statement | grade |
|---|---|
| **(A) Helicity lemma.** For any three points `x₁, x₂, x₅` and a null `z₁` with `V_{1,25}(z₁) = 0` (one of the two null eigen-directions of the little group `SO(2)_ℂ` of the three points, i.e. a helicity eigen-polarisation of `T₁`), the symmetric conserved parity-odd structure restricted to `z₁` is a fixed imaginary combination of the even ones, **as a function of `(z₂, z₅)`**: `⟨T₁T₂O⟩_odd = ±i (w₁⟨…⟩_{E₁} + w₂⟨…⟩_{E₂})`. `w` depends only on `(Δ, J)` and the basis normalisation; the sign is the chirality of the root. | **[D]** (orbit + covariance part [P]; existence of the relation from the helicity decomposition of conserved structures [D, input below]); **[DV]**: 4 random complex configurations × 2 roots, `J = 4` (`Δ = 400`, residual `≤ 2·10⁻⁸`) and `J = 5` (`Δ = 400.37`, antisymmetric even set, `≤ 2·10⁻¹²`); control (generic null `z₁`): `10⁻³ … 0.6` |
| **(B) The crystal saddle has `V_{1,25}(x₁₀) = 0`.** `x₁ = (−ia, ia, 0)`, `x₂ = (i/a, −i/a, 0)`, `e = (1, 1, 0)`, `x₀ = (0, 0, n₀a)`: `z₁ = x₁₀` is null and `z₁·(x₁₂/x₁₂² − x₁₅/x₁₅²) = −a/(a + b) + a²/(1 + a²) = 0` for `b = 1/a`. The `S⁰` label `n₀` is the chirality. Likewise on the 34 side. | **[P]** (2 lines; numerically `0.0` at the saddle) |
| **(C) Order 0.** By (A)+(B) the odd row and column of the order-0 (saddle) 3×3 action are `w`- and `v`-combinations of the even ones, `w, v` imaginary, `v = −w` up to the 34-side basis; with the order-0 helicity decomposition this forces `w ∝ t₊₊*` (lemma L2 of `po_sympy.py`), so `M0_OO = κ (T0⁻¹ M0_EE T0⁻ᵀ)₊₊` with `κ = (w·t₊₊)(v·t₊₊)`. | **[P given (A)]**; numerically `≤ 10⁻¹⁴` on every K₂ record |
| **(D) All orders, integrated.** `M_OO = κ (T0⁻¹ M_EE T0⁻ᵀ)₊₊` holds for the **summed Laplace series** (`M0 + M2 + M4`) with the order-0 `κ`, at every node of 7 rows (`J = 12–22`): `κ(ν)` constant to `≤ 1.1·10⁻⁸` relative on the states (`Δ = 336, 400`) and `≤ 10⁻¹³` on the references for `Ĉ'_EE`; `≤ 3.5·10⁻⁵` / `8·10⁻¹⁰` for `Ĉ₃` (the larger values at nodes where the `Ĉ₃` `++` element is small). Absolute: `|Γ_OO − Γ_++| ≤ 0.05 c₂` (`Ĉ'_EE`) / `0.12 c₂` (`Ĉ₃`) (HYPOTHESES_CLOSE2 §1.4); e.g. `J = 12`, `Ĉ'_EE`: `≤ 2·10⁻⁴ c₂` except one node (`ν = 298.7`, `M2/M0 = 0.8`): `0.013 c₂`. Equivalently (L1 of `po_sympy.py`): the **`++ ↔ −−` helicity-flip element of the S⁰-summed functional vanishes**, `X = W_{++,−−} + W_{−−,++} = 0` to `10⁻⁸ D` (`Ĉ'_EE`). On the previous agent's rows `J = 4–10` (`hc_p3`, `po_identity_lowJ.py`): `κ` spread `2·10⁻¹⁰ … 4·10⁻⁹` (Δ = 400) → `2·10⁻¹³ … 5·10⁻¹²` (Δ = 1600) for `Ĉ'_EE`, `1.5·10⁻⁶ … 1.3·10⁻⁵` → `6·10⁻⁹ … 6·10⁻⁸` for `Ĉ₃`: the residual **falls like `Δ⁻⁴…⁻⁵`**, as the point-wise residual (E) does. So `X/D = O(Δ⁻⁴)`, i.e. **`Γ_OO − Γ_++ = O(m⁻²) c₂`** (a relative `O(1/m⁴)` effect, beyond the order at which the corrections are defined). | **[NE]** (exact vanishing: not supported; `O(Δ⁻⁴)` scaling **[NE]**; a proof of the scaling **[I]**) |
| **(E) It is not a point-wise identity.** Off the saddle `V_{1,25}(x₁₀) ≠ 0` and (A) fails point-wise: F-level residual `2·10⁻⁷ / 10⁻⁵ / 4·10⁻³` at `|V₁| ≈ 0.01 / 0.04 / 0.12`; prefactor level `OO − wᵀ pref_EE v` relative `2·10⁻¹² … 3·10⁻¹⁰` (`Δ = 400`, `|ξ| = 0.1 … 1.6`), `10⁻⁹ … 10⁻⁷` (`Δ = 100`), falling `∝ Δ⁻⁴`; the single-`n₀` 3×3 action is rank 2 only at order 0 (its `EO` entries leave the rank-2 pattern at `10⁻⁴ … 10⁻³` relative in `M2, M4`), while the S⁰-summed `OO`/`++` relation (D) survives to `10⁻⁸`. | **[NE]** |
| **(F) Consequence: `Γ_OO ≡ Γ_++` and the per-row constant.** With (D), `G_OO,st/G_OO,ref = [κ(Δ)/κ(Δ_ref)] G_++,st/G_++,ref` at every `ν`; the per-row fit constant absorbs the bracket (L3 of `po_sympy.py`), so `Γ_OO = Γ_++` node by node and **`a_O/a_++ = κ(Δ)/κ(Δ_ref)`**: computed `2.20362 / 2.09170 / 1.74393 / 1.77552 / 1.14080` (rows `(336,18) / (400,18) / (336,22) / (400,12) / (400,14)`; also `1.18413` at `(400,16)`, `2.16511` at `(400,20)`), identical for `Ĉ'_EE` and `Ĉ₃` to `10⁻⁶`: exactly HYPOTHESES_CLOSE2's `2.20 / 2.09 / 1.74 / 1.77 / 1.14`. `κ` is the squared projection of the `p1` odd structure (unit-max coefficient vector) on the `T0` `++` column (unit norm, SVD-defined per `Δ`), over the shadow `β` ratio — a **basis-normalisation constant, not dynamics**; it is not `tan²(πτ/2)` (that vanishes at the even-integer rows) but carries the same information once CHL's normalisation is used (PARITY_ODD_ROUTEA §12.4, `N_full = 1`). | (F) **[P given (D)]**; numbers **[DV]** |
| **(G) Odd `J`.** (A) holds with the antisymmetric even pair (`X^E_same`, `X^E_opp` antisymmetric at odd `J`): `w = ±(4.834278i, 2.184051i)` at `(400.37, 5)`, the same `w` as the K₂ rows `p3_pts_C2_J5_anti`; on those rows (`J = 5, 7`; `Δ = 400 … 3200.37`) the odd entry equals `κ ×` one helicity diagonal of the antisymmetric even block (partner channel `≤ 4·10⁻⁴` of `w`, limited by the 6-node rank-1 fit; `κ` ν-constant to `10⁻⁴ … 2·10⁻⁷`). So the parity-odd `+−` exchange at odd `J` is, node by node, the `+−` diagonal of the (unphysical for identical `T`) antisymmetric even block — PARITY_ODD_ROUTEA §4's "`N` is a normalisation of the whole block", now explained. It gives no independent check of §2 (that block is not computed on the engine). | (A) at odd J **[DV]**; block statement **[NE]** |

**Verdict on goal 1.** The identity is **derived at order 0** (= the crystal saddle, where the functional's `T₁`/`T₃` polarisations are
helicity eigen-polarisations; [P given (A)], (A) itself [D]/[DV]), and the per-row constants are **explained exactly** ([DV], 6 digits).
Beyond order 0 it is **not** a point-wise consequence of (A): it holds as the integrated statement "`X ≈ 0`" (no `++ ↔ −−` flip)
to `10⁻⁸` (`Ĉ'_EE`) / `4·10⁻⁵` (`Ĉ₃`) relative — `≤ 0.05 / 0.12 c₂` — on every computed row, with a residual falling like `Δ⁻⁴`
(`Γ_OO − Γ_++ = O(m⁻²) c₂`) [NE]. So the precise statement is **`Γ_OO = Γ_++` exactly at order 0 and up to `O(m⁻²) c₂`**, not an exact
identity; why the first corrections to the lemma cancel in `X` is open [I]. For Result D nothing changes: the
residual is `10³–10⁵×` below every tolerance (`H_OO` margins `+0.84 … +1.01`, failure thresholds `≳ 10⁴ c₂`), so "`H_OO` positivity ⇐
EE-block positivity" (HYPOTHESES_CLOSE2 §1.5) now rests on (D) [NE, `10⁻⁸`] plus (C) [P given (A)] instead of a bare numerical identity.

### 1.2 The lemma (A): argument
1. **Orbits [P].** Three (complex) points are conformally equivalent to `(0, e′, ∞)`, whose stabiliser is `SO(2)_ℂ` (rotations about
   `e′`). A null `z₁` at `x₁ = 0` with `V_{1,25} ∝ z₁·e′ = 0` is one of the two null eigenvectors `η_±` of that rotation (eigenvalues
   `e^{±iθ}`); scale is irrelevant (structures are homogeneous of degree 2 in `z₁`). So `{(x₁, x₂, x₅, z₁): z₁² = 0, V_{1,25} = 0}` is
   two orbits of the connected complex conformal group, exchanged by parity.
2. **Covariance [P].** Under the connected group the odd structure transforms like the even ones (`det = +1`). A linear relation
   `O| = Σ w_a E_a|` valid at one point of an orbit therefore holds on the whole orbit with the same `w`; parity maps `O → −O`, `E → E`
   and one orbit to the other, so `w → −w`. This is why `w` is configuration-independent and why its sign is `n₀` (§1.1 (B)).
3. **Existence [D].** In the frame `(0, e′, ∞)`, restricting `T₁` to `η₊` keeps only its helicity-`+2` component along `e′`. Conserved
   `⟨TTO⟩` structures are fixed by their components with both `T`'s at helicity `±2` (Caron-Huot–Li's helicity structures `|h₁h₂⟩`;
   input, not re-proved here), so the restricted space is spanned by `|++⟩|` and `|+−⟩|`: two functions for three structures. With
   `E_same ∝ |++⟩ + |−−⟩`, `O_same ∝ i(|++⟩ − |−−⟩)` and `E_opp ∝ |+−⟩ + |−+⟩`, `O| = i c E_same|`. That is (A) with `Σ w_a E_a ∝ E_same`.
   The tests (`po_helicity_lemma.py`) confirm it, and confirm `w` against the K₂ rows' `w` (`J = 4`: `±(4.289479i, −0.965631i)` at
   `Δ = 400` from both the random configurations and the functional's saddle).

### 1.3 Symbolic checks (`po_sympy.py`, `cache/po_sympy.txt`) [P]
* **L1** (helicity parametrisation, any complex `W` and its parity image): `E_sE_s = 2(D + X)`, `OO = 2|c|²(D − X)`; so `OO ∝ E_sE_s`
  at every `ν` ⟺ `X = 0`.
* **L2** (order 0): if `M0_EE = λ₁u₁u₁′ᵀ + λ₂u₂u₂′ᵀ` and the `O` row is `∝ u₁′ᵀ`, the row null vector satisfies `w·u₂ = 0`, `w·u₁ = γ`.
* **L3**: `wᵀAv = αβ (T⁻¹AT⁻ᵀ)₁₁` for `w = αt₁*`, `v = βt₁*`; and `Γ_OO − Γ_++ = 0` with `a_O = a κ_s/κ_r`.
* Instantiation on the rows: `po_identity.py` → `cache/po_identity.txt` (all 7 `oo_p3` rows, both kinds); odd `J`:
  `po_identity_odd.py` → `cache/po_identity_odd.txt`; point-wise: `po_pwtest.py`, `po_pwscan.py` → `cache/po_pwscan.txt`.

## 2. Result Q: the parity-odd `+−` J = 5 state on the analytic-radial engine

### 2.1 Construction and validation (`po_ee_odd.py`) [DV]
* **Prefactor.** `p3_odd_action.tto_odd_samples` / `eps_cof` re-built term by term on `bg_ser.ASer` (embedding vectors
  `(1, x², x)`, `(0, 2x·z, z)`, `ε` by cofactors, `ε(Z_a, Z_b, X_a, X_b, X₅)X_ab × even monomials`, real coefficient vector
  `p1_tto_odd.conserved_odd` at `Δ` and `3 − Δ`), contracted exactly as `bg_ee` / `w0_ee_common` (`E ⊗ E`, null-circle Fourier
  contraction). The odd-structure factors `(m₁₂m₁₅m₂₅)^{−1/2}`, `(m₃₄m₃₅m₄₅)^{−1/2}` (p3's `ex12`, `ex34`) go into the log list with
  exponent `−½` each; everything else (EE slow exponent, `v′ − u′`, Lemma-B `R` for `k = 4`, `B`-factor, closed-form integrals) is
  `bg_eeq`'s. At odd `J` the symmetric conserved set is the odd structure alone, so the image is 1×1 and the shadow coefficient
  (ν-independent) drops out with the per-state constant `c`.
* **Point-wise check against the project's exact odd integrand** (`p8_odd_pw.P8` + `w0.make_build_odd_EE`, kind `EE`, on the real
  section of `bf_core`, 8 points, `Δ = 300`, `J = 5`; `cache/po_ee_odd_check.txt`): ratio series/exact `→ 1` order by order,
  max deviation **`3.6·10⁻² / 1.9·10⁻³ / 5.4·10⁻⁵`** at `a⁰ / a¹ / a²` (the even EE check of BELOW_FOLD2 §1.3: `0.11 / 1.1·10⁻² /
  6.3·10⁻⁴`); the ratio is real-positive (no branch or phase slip of the half-integer powers) [DV].
* **Integrated check against direct non-perturbative points** (`po_k2cmp.py`, `cache/po_k2cmp.txt`): engine EE `k = 2` image (K = 2)
  at `(400, 5)` vs the existing n3 direct points `hc_pts_EE_J5_D400` (same odd integrand, CRN state/reference),
  `m²[D(ν) − D(39.85)]`: `ν = 19.92`: engine `+12 539` vs direct `+13 300 ± 520`; `ν = 9.96`: `+15 654` vs `+16 594 ± 3 823`
  (1.5σ / 0.25σ; the engine's own K = 2 truncation at `ν = 20`, `Δ = 400` is not estimated here) [NE].
* **Cost.** Depth 5 (K = 4): 29 min, 0.9 GB per state; reference (`Δ = 2·10⁴`, K = 2, depth 3) 4 min. Same order structure as the
  even `J = 4` images (`a⁻⁸` leading, one order lost to a leading cancellation, as bg_ee).

### 2.2 Chain (`po_qan5.py`) [D, as bg_qan / bh_qan4]
`Γ_exact = m²(c T_ref I_st/I_ref − T_st) + γ(m/m_ref)² Σ_j R_j(ν − ij)²`, `Γ_K = Γ_exact − γ Σ_j R_j(ν − ij)² − const·N̂₄`
(`γ = γ_{+−} = 20.75`, HYB; REFCORR = 1 as bh_qan4), `s_d = 0.04 Σ χ_d Re Γ_K`, gauges H0/H2/H3/H4 (Stage 6), SB (canonical), XP
(Stage 7 x-precision); FF/N from the Stage-6 HYB model at the state, `FF → FF − a·2u³ s_d − c(2u³/Λ²) s_v`. Truncation error =
size of the a⁴ term.

### 2.3 Result [NE + M]

| state | Δ_g | `s_d` H0 | H2 | H3 | H4 | SB | XP | FF/N δ = 0.05 (K=0 → meas.) | δ = 0.2 |
|---|---|---|---|---|---|---|---|---|---|
| **`+−` J = 5, gap `1.0167Λ`** (Δ = 305.03) | 300 | **−0.539 ± 0.033** | −0.540 | −0.539 | −0.538 | −0.554 | −0.545 | +0.125 → **+0.144 ± 0.001** (SB +0.145) | +0.261 → **+0.284 ± 0.001** (SB +0.284) |
| **`+−` J = 5, gap `1.0050Λ`** (Δ = 1005.01) | 1000 | **−0.512 ± 0.004** | −0.513 | −0.512 | −0.511 | −0.527 | −0.518 | +0.134 → **+0.153 ± 0.0001** (SB +0.153) | +0.270 → **+0.291 ± 0.0002** (SB +0.292) |
| `+−` J = 5, `1.1935Λ` (Δ = 357.8; binding region of the J = 5 channel, 1.187–1.206Λ) | 300 | −0.531 ± 0.026 | −0.532 | −0.531 | −0.530 | −0.546 | −0.537 | +0.0505 → **+0.0654 ± 0.0007** (SB +0.0659) | +0.2003 → **+0.2182 ± 0.0009** (SB +0.2187) |

* **Sign: favourable.** `s_d < 0` raises FF in every gauge at both `Δ_g`; the gauge spread is `≤ 0.016`. Per-order `s_d` (H0, Δ_g = 300):
  `+19.25 / −22.36 / +2.30 / +0.24 / +0.033` (k = 0…4; Δ_g = 1000: `+19.24 / −23.67 / +3.81 / +0.12 / +0.004`): large, cancelling
  low orders (as for `+−` J = 4: `+3.83 / −9.91 / +5.02 / +0.51 / +0.04`), then geometric convergence (×10 per order).
* **Cross-check with the helicity-odd universal correction (BELOW_FOLD3 §1) [NE].** Per harmonic, `Γ_K/Σ|R_j|` at `ν = 0.5 … 10`:
  J = 5 (odd) `+1.05 +4.58 +8.41 +10.77 +10.38 +7.71 +0.77 −3.21 −3.15 +1.44 +2.41` vs `+−` J = 4 (even) `+1.03 +4.56 +8.34 +10.61
  +10.23 +7.77 +1.66 −1.86 −2.23 +0.48 +1.40` (Δ_g = 300); at Δ_g = 1000 the two agree to `≤ 0.1` per harmonic. So the parity-odd
  `+−` exchange carries the **same** `σ_{+−} U(ν)` (`σ_{+−} = −1`) as the parity-even `+−` exchange: `s_d(J=5) = −0.539 / −0.512` vs
  `s_d(J=4) = −0.516 / −0.510`. This is what §1 (G) predicts structurally (the odd `+−` block is a `+−` helicity diagonal) and what the
  universal-`σ_h` extrapolation of BELOW_FOLD3 assumed with the unfavourable sign bracketed; the bracket is now resolved favourably.
* Covariant part (bg_covfit rule, ν < 12): `δγ = −0.10` (−0.5 % of `γ_{+−}`) at 300, `−0.02` at 1000; the rest is the universal
  non-covariant part.
* Caveats: the engine's EE beyond-LO equivalence with the crystal cycle (BELOW_FOLD2 §5.4) is assumed here as there; `s_d` is applied
  m-independently within (`+−`, J = 5) [M] (tested at one more mass, row above); the HYB model elsewhere [M].

### 2.4 Result Q full-set minimum (`po_ffscan.py` = bh_ffscan + the J = 5 pickles → `cache/po_ffscan.txt`) [NE + M]
Measured `s_d` now on `++` J = 0, 2, 4, `+−` J = 4 and **`+−` J = 5** (Δ_g = 300: mean of the gap and 1.19Λ measurements, `−0.535`,
error = max(truncation, half-spread)); applied m-independently within each (h, J) [M]; other states at K = 0.

| δ | Δ_g | **min FF/N** (H0 / SB / XP) | at | `+−` J = 5 channel (K = 0 → measured) | `+−` J = 4 | `++` J = 2 |
|---|---|---|---|---|---|---|
| 0.05 | 300 | **+0.0061 ± 0.0013** / +0.0048 / +0.0056 | `++` J = 0, 1.000Λ | +0.0499 → **+0.0648 ± 0.0009** (1.206Λ) | +0.117 | +0.020 |
| 0.05 | 1000 | **+0.0061 ± 0.0013** / +0.0048 / +0.0056 | `++` J = 0 | +0.0500 → **+0.0643 ± 0.0001** (1.206Λ) | +0.099 | +0.020 |
| 0.2 | 300 | **+0.1427 ± 0.0017** / +0.1410 / +0.1421 | `++` J = 0 | +0.2000 → **+0.2182 ± 0.0011** (1.187Λ) | +0.277 | +0.167 |
| 0.2 | 1000 | **+0.1427 ± 0.0017** / +0.1410 / +0.1421 | `++` J = 0 | +0.2000 → **+0.2175 ± 0.0001** (1.187Λ) | +0.259 | +0.167 |

* **The Result Q full-set minimum is unchanged**: `+0.0061 ± 0.0013` (δ = 0.05; SB `+0.0048`) and `+0.1427 ± 0.0017` (δ = 0.2; SB
  `+0.1410`), binding state `++` J = 0 at the gap, both Δ_g, every gauge. The last unmeasured near-binding channel of BELOW_FOLD3 (`+−`
  J = 5, which would have failed δ = 0.2 only for `s_d ≳ +5`) is measured with the **favourable** sign (`s_d ≈ −0.53`): its minimum
  rises from +0.050 / +0.200 to **+0.065 / +0.218** (δ = 0.05 / 0.2) [NE + M].
* δ = 0.2: closed on every measured binding channel with margin ≥ +0.141 [NE + M]. δ = 0.05: still +0.005 ± 0.001 at `++` J = 0 (the
  `σ_h`-odd universal-U double-counting question of BELOW_FOLD3 §4.2 is untouched here) [M].
* Remaining odd-J `+−` (J = 7, 9, …): not measured; the J = 5 measurement shows the parity-odd `+−` exchange carries the even `+−`
  `σ_{+−} = −1` correction, so the favourable-sign universal extrapolation applies to them [M] (the worst-sign bracket of bh_ffscan,
  +0.036 / +0.183, is no longer needed for J = 5 and is above the minimum for the others anyway).

## 3. Item 3 (partial): direct small-`ν` nodes at `Δ_g = 300` for `J = 6` (even) and `J = 15` (odd)

Method unchanged (HYPOTHESES_CLOSE2 §1.1): n3 hybrid contour, `2¹⁰·2·8` points, state and reference `2·10⁴` on common random numbers.
* **`J = 6`** (`oo_points.py 6 400 x`, `x = 0.021 … 0.05`, `ν = 8.37 … 19.92`; 3×3, i.e. EE block and OO block on the same samples;
  RQMC relative errors `8 % / 1.7 % / 0.6 % / 0.5 % / 0.3 % / 0.2 %` for EE, scalar `4.4 % … 0.2 %`). Chain: `po_t6.sh` = oo_t.sh
  variants for `J = 6` (`OO_EESRC=oo`, `OO_Z3EE=6:400`, twist edge `1.0199Λ` added), summary `oo_sum.py po6`:

  | `J = 6` | best | dir / cv / zero / lin | convention A | jackknife |
  |---|---|---|---|---|
  | EE block (run0), near-gap min | **`+0.3075`** [1.25Λ] | `+0.307` (all) | `+0.307` | `± 0.0001` |
  | OO block, near-gap min | **`+0.3078`** [1.25Λ] (LO `+0.3091`) | `+0.308` (all) | `+0.307` | `± 0.0001` |

  The direct nodes (`Γ_c = −29 −43i, +14, +32, +42, +51, +49 c₂` OO; EE `++` `−106 −43i, −62, −44, −33, −22, −21`) do not move the
  `J = 6` minimum, which sits at `1.25Λ` where the small-`ν` window has no weight [NE + M]. (The OO and `++` rows differ here by a
  smooth `≈ 75 c₂` offset with equal imaginary parts: the `a`-fit offset `∝ D_h` of HYPOTHESES_CLOSE2 §1.4 — the `J = 6` EE rows carry
  a different node set from the OO rows — not a violation of §1 (D), which holds at the level of `M`.)
* **`J = 15`** (odd; `hc_points.py` n3 path, rows 400 (`x = 0.021 … 0.05`) and 336 (`x = 0.025 … 0.0625`, `ν = 8.35 … 20.9`), 1×1
  odd EE / `Ĉ₃` integrand; RQMC errors `7 % / 1.4 % / 0.6 % / 0.5 % / 0.3 % / 0.2 %`). z4's non-local small-`ν` term is wired for even
  `J` only (`HC_Z3` / `Z4_JS = 15` leave odd `J` unchanged, checked: `+0.6481` with and without), so the nodes enter by **row injection**
  (`HC_INJ = 15:336,400`, `hc_chain.py z4 300`, `po_t15.sh`): near-gap min **`+0.829 ± 0.016`** (jackknife over the 8 replicas)
  at `1.0586Λ`, against convention A `+0.648` (same run settings; `+0.630` in HYPOTHESES_CLOSE's sweep); without the Hampel rule on
  injected nodes `+0.829`; all-`m` minimum `+0.829`. Per state (convention A → direct): `1.05Λ +0.707 → +0.928`, `1.0586Λ +0.648 →
  +0.829`, `1.1Λ +0.747 → +0.913`, `1.2Λ +0.761 → +0.854`, `1.3Λ +0.900 → +0.977`. Convention A was pessimistic again [NE + M]
  (below the lowest node, `ν < 8.35`, the row is held as in convention A [M]).
* **Not done:** even `J = 8, 10`; odd `J = 7, 9, 11, 13` (each ≈ 30 min of nodes per row + chain; same recipe). The binding state at
  `Δ_g = 300` is unchanged (`J = 5`, `+0.2038`, hybrid twist edge; HYPOTHESES_CLOSE2).

## 4. Remaining open items [I]
1. **Why `X = O(Δ⁻⁴)`** (no `++ ↔ −−` flip of the S⁰-summed functional up to `O(m⁻²) c₂`) beyond the crystal saddle; the lemma gives
   order 0 only, and the integrated statement is numerical (`10⁻⁸` relative on the K₂ rows, `Δ⁻⁴…⁻⁵` scaling; within RQMC noise on the
   direct nodes). A likely
   route: the off-saddle violations of the lemma on the 12 and 34 sides cancel in the OO combination (the single-`n₀` `EO` entries
   leave the rank-2 pattern at `10⁻⁴`, the OO/`++` relation does not) — not shown.
2. The existence step of the lemma uses the input "conserved `⟨TTO⟩` structures are fixed by their helicity-`±2` components"
   (Caron-Huot–Li); verified numerically, not re-derived.
3. Result Q: the engine's EE beyond-LO equivalence with the crystal cycle (BELOW_FOLD2 §5.4) is assumed for the odd structure as for
   the even; `s_d` m-independence within `+−` J = 5 is checked at two masses (300) only [M]; odd `J ≥ 7` `+−` unmeasured [M]; HYB
   double counting of the universal `σ_h U` (BELOW_FOLD3 §4.2) untouched — it decides the δ = 0.05 sign, not δ = 0.2.
4. Item 3: even `J = 8, 10` and odd `J = 7–13` not computed; odd-`J` direct nodes enter the chain only through `HC_INJ` (row
   injection), not through z4's non-local term.

## 5. Files (`ads4/code/spinfunc/`; caches in `cache/`)
| file | content | run |
|---|---|---|
| `po_pwtest.py` | point-wise row/column/OO relations of the 3×3 prefactor at random points (`EE`, `C3H`) | 1 min |
| `po_helicity_lemma.py` | lemma (A): random configurations with `V_{1,25} = 0`, controls, the saddle; `off` mode: off-saddle `V₁` → `cache/po_helicity_lemma.txt` | 1 min |
| `po_identity.py` | null vectors, `w·t₊₋ = 0`, `κ(ν)` constancy, per-row constants on the 7 `oo_p3` rows → `cache/po_identity.txt` | 10 s |
| `po_identity_lowJ.py` | the same on the J = 4–10 `hc_p3` rows → `cache/po_identity_lowJ.txt` | 10 s |
| `po_identity_odd.py` | the same at odd `J` on `p3_pts_C2_J{5,7}_anti` → `cache/po_identity_odd.txt` | 5 s |
| `po_pwscan.py` | point-wise OO-relation residual vs `|ξ|`, `Δ` → `cache/po_pwscan.txt` | 2 min |
| `po_sympy.py` | symbolic lemmas L1–L3 → `cache/po_sympy.txt` | 5 s |
| `po_ee_odd.py` | parity-odd EE prefactor on the a-series engine; `check` (point-wise) and `run` (images) → `cache/po_eeq_D*_J5_K*.pkl`, `po_ee_odd_check.txt` | 29 min (depth 5), 4 min (ref) |
| `po_queue.sh` | `A`: reference + Δ_g = 1000 gap; `B`: Δ_g = 300 gap; `C`: 1.19Λ | — |
| `po_qan5.py` | Result Q chain for the 1×1 odd image (s_d, s_v, FF/N, J = 4 cross-check) → `cache/po_qan5.txt`, `po_qan5_D*_J5_Dg*.pkl` | 1 min |
| `po_ffscan.py` | bh_ffscan + J = 5 → `cache/po_ffscan.txt` | 2 min |
| `po_k2cmp.py` | engine EE k = 2 vs direct n3 points, `(400, 5)` → `cache/po_k2cmp.txt` | 7 min |
| `po_item3.sh`, `po_t6.sh` | `J = 6` direct nodes (`oo_points.py`, unchanged) and chain (`oo_chain.py`, tag `po6_*`) → `cache/oo_pts_*_J6_D400.pkl`, `oo_chain_po6_*` | 27 min + 15 min |
| `po_item3b.sh`, `po_t15.sh` | `J = 15` direct nodes (`hc_points.py`, n3, rows 400/336) and chain (`hc_chain.py z4`, tags `po15_*`) → `cache/hc_pts_*_J15_D*.pkl`, `hc_chain_z4_po15_*` | 30 min per row |

Reproduce goal 2: `./po_queue.sh A; ./po_queue.sh B; ./po_queue.sh C` (or `python3 po_ee_odd.py run <Δ> 5 4 5`, reference
`python3 po_ee_odd.py run 20000 5 2 3`), then `python3 po_qan5.py 305.0297506031432 5 4 300`,
`python3 po_qan5.py 1005.0089775159122 5 4 1000`, `python3 po_qan5.py 357.8 5 4 300`, `python3 po_ffscan.py`.
Nothing of mine is left running at hand-over.
