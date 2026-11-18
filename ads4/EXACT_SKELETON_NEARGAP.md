# Exact scalar-skeleton heavy moments on the near-gap states (Δ_g = 100): Ω^U threshold without the tolerance convention

Status: 18 Nov 2026, derivation/numerics agent. **Unrefereed.** Second pass: §6 (coverage to `m/Λ = 2.8`, the bound, `Δ_g = 300`, NF02). Some runs are still in progress; see §6.
Code: `ads4/code/exactskel/` (`run_all.sh`; caches in `exactskel/cache/`, git-ignored). It imports, read-only, `omegau/`
(common, adsaction, remk, u4_scan), `robust/` and `uniform/skel.py`. Nothing outside `exactskel/` was modified, except one line in `.gitignore`.

Marks: [D] derived, [DV] derived and checked numerically, [NE] floating-point numerics, [M] model or cut, [I] inferred or assumed.

## 1. Method [D + DV]

**Why direct position-space moments fail.** The real-section moments `Π̂_{k,ν}` lose a factor of about `e^{−πν}` to cancellation
(with `ν` up to ≈ 98). That needs roughly 100–150 digits, so the radial-block route of `e2e/` and `lowk/` is out of reach for these nodes.

**What is used instead: the CHMRS Mellin representation** (arXiv:2106.10274, eqs. CfromBMellin and Bkt action). This is what CHMRS
themselves used for "exact" actions.
* `C_{k,ν} = 2^{1−2k}π^{3/2} a_{Δφ} γ_a(ν)^{−2} ∫dT/(4πi) Γ(T/2)² γ_{a−T}(ν) B̂_{k,T}`.
* `B̂` is a sum over descendants of Mack-polynomial residues.
* Mack polynomials (`mack.py`) come from the explicit formula of Costa–Gonçalves–Penedones, built exactly in rational arithmetic.
  The convention map is `s = −T`, `t = τ + 2n`, with normalisation from arXiv:2008.04931 (qCal). Checks:
  * `J = 1` agrees with 2008.04931;
  * the `y`-parity symmetry holds to 1e−31.
* **Contour shift.** Moving the contour from `Re T = c₀` to `c₁ ≈ min(T*(ν), τ − ½)` picks up the residues at `T = A ± iν + 2j`.
  This removes the `e^{−πν/2}` cancellation. The remaining cancellation is ≤ 10⁹ (≤ 10¹³ only at nodes where the moment itself nearly
  vanishes).
* **`n`-sum.** It is summed directly for `n < 400`; the tail uses midpoint Euler–Maclaurin with Gauss quadrature.
* **Odd `k`.** Pochhammers are continued with Gamma ratios.

## 2. Validation [DV]

**Against the CHMRS order-1/m⁸ series where it is controlled** (`t_conv`, `t_check2`): `Δ = 120–160`, `J = 0, 2, 4, 5, 10`,
`Δφ = 2, 2.5, 3`, `k = 2…7`, `ν ≤ 15`.

| `k` | relative agreement |
|---|---|
| 2 | 1e−12 |
| 4 | 1e−11 |
| 6, 7 | 1e−10 |

* This includes odd `J`.
* It includes odd `k`. **So the Gamma-ratio continuation in `k` is, to 1e−10, the same function as OMEGA_U's "k-continued series".**
* The forward terms `R₄` and `R₅` agree with the series to about 1e−4 (relative).

**Stability** (`t_largenu`, `e2`): independent of the contour position, the line step, `N₀` and the number of Euler–Maclaurin nodes.
* At `m = Λ`, `x = 0.998` (worst case) the variation is ≤ 2e−6 relative.
* Elsewhere it is ≤ 1e−8.
* Rounding estimate on `ρ_k`: ≤ 2e−7 (`Δφ = 3`) and ≤ 1.2e−5 (`Δφ = 2`).

**What these numbers mean for the skeleton part.** The moments themselves are [DV]. The skeleton part of the action carries numerical
error ≲ 1e−4 `N_s`. The p-quadrature check (`e4_quad.py`, 64 vs 96 vs 128 nodes) was written but **not run**.

## 3. The truncated series against exact (`e5`, `cache/e5_100.txt`)

**For `x < 0.2`:** `|Δρ| ≲ 2e−4`. This is the `1/m⁶` truncation of the fixed-`J` series.

**For `x ∈ [0.4, 0.6)`:** `|Δρ|` reaches `1e−2` to `0.1`.

**For `x > 0.6`:** the series is useless. It gives `ρ` of about 0.5 %, while the exact `|ρ|` reaches about 0.9, and exact `ρ → −1` as
`x → 1` at `m ≈ Λ`. There the exact moment stays finite while `E_k` diverges.

**The series' own error estimate fails on "controlled" nodes.**
* Several thousand nodes are affected for `x ≤ 0.6`, with up to 5× excess.
* It also fails at `x → 1` on nodes that the absolute `1e−4·E` switch declared controlled.
* So the earlier `E`/`U` bookkeeping was **not conservative**.

## 4. Ω^U action with the exact skeleton (`e3`, `cache/e3_100.txt`) [NE + I]

**Model.** It is the OMEGA_U model, with `ρ_k` exact on every node and exact `R₄` and `R₅`. There is no control cut and no `κ`.
* `F_exact = F_hyb + S_exact`.
* `+−` uses the full `d^J_{44}` (reading B) with the scalar `ρ` transferred.
* Odd `k` uses the B′ envelope.
* **Reading E**, as an I3 test, takes odd moments to be `ε_k × (the exact continued scalar)`.

**Scope.** `Δ_g = 100`, `m/Λ ∈ [1, 1.55]` on the u4 grid, `J ≤ 10`, `Δφ = 3` and `2`. This is about 190 `++` and 200 `+−` states per
skeleton.

**M02 results:**
* **Non-positive states.** They exist only at `m ≤ 1.07Λ`:
  * `+−, J = 4`: down to `−1.39 N_s` at `1.04Λ`;
  * `++, J = 0` at `m = Λ`: `−0.10 N_s` (`Δφ = 3`) and `−0.15 N_s` (`Δφ = 2`).
* **For `m ≥ 1.08Λ` every state is positive.** The minimum is `F/N_s = +0.25` (`++, J = 2`, `1.2–1.4Λ`). The exact remainder there is
  ≤ 0.012 `N_s`.
* **True threshold `n = 1.08`.** It is the same for both skeletons and both odd-`k` readings. With `δ = 0.03` it becomes 1.08–1.09.
* **Comparison.** OMEGA_U's series criterion on the same states gives `n(κ_req)` = 1.15 / 1.27 / 1.40 for `κ_req` = 0.1 / 0.3 / 1.
  The tolerance-convention `n` was therefore conservative, even though its error bars were not.
* **What the tolerance actually had to cover.** On the uncontrolled nodes, `κ_eff = S_exact(unc)/U` is between −0.1 (`+−, J = 4`,
  1.04Λ) and −0.001.
* **Odd-`k` sensitivity (H3).** With exact `ρ`, the I3 prefactor ambiguity is suppressed near `x → 1`. The reason is that `D + ρ → 0`
  there, so the shift is ≤ 0.06 `N_s` for `m ≥ 1.08Λ`.

**M01:** `n = 1.10`. **F_g3:** `n = 1.08–1.09`. Its `++` island at `5.7Λ` was not scanned.

**Consequence for the bound (not recomputed).** The `n = 1.08` threshold sits below every `n(κ_req)` row of OMEGA_U's §0 item 7 table
(1.15–1.40). By the OMEGA_U formula, the bound would therefore come out below the `κ_req = 0.1` entry, about `5·10⁻⁵`, at `Δ_g = 100`.
This is still conditional on H1–H4 and H6–H12. **The skeleton part no longer needs H5.**

## 5. Not done / caveats

* **The scan stopped at `m/Λ = 1.55`** (runs were in progress at hand-back).
  * `[1.55, 2]` is covered only by OMEGA_U, with `κ ≥ 4`.
  * `Δφ = 2.5` was not run.
  * The quadrature check `e4` and the full `e2` were not run. `run_all.sh` runs all of them.
* **`J > 10` (from `m ≈ 1.105Λ`)** keeps OMEGA_U's `κ ≥ 3.5` [M].
* **`Δ_g = 300`** needs the contour to cross the `B̂` poles. The cancellation is about 1e26 at `c₁ < τ`. Not done.
* **NF02 / H13** (large `b̂`, `J` up to about 100). Not done.
* **Still [I]:** H4 (skeleton to graviton transfer, and the `+−` transfer), H3 (odd-`k` two-pair moments), and H1/H2.

## 6. Second pass (coverage, bound, `Δ_g = 300`, NF02)

### Coverage at `Δ_g = 100` [NE + I]

**Scan.** `J ≤ 10`, `m/Λ ∈ [1, 2.8]`, `Δφ = 3` and `2`: about 700 helicity-states per skeleton.
* Unchanged: M02 is non-positive only at `m ≤ 1.07Λ`.
* `n_exact` = 1.08 (`δ = 0`); 1.09 at the worst of `δ = 0.03` and odd-`k` reading E.
* On `m ∈ [1.75, 2.8]` the minimum `F/N_s` is `+0.21` (`++, J = 4`).

**High-J support (`mack.py`).**
* `mack_P_fast` reorganises the exact rational construction. Its cost is still roughly `J⁶`: J = 20 takes 10 s, J = 24 takes 30 s.
* `Qpoly_cheb` evaluates `Q^n` in a scaled Chebyshev basis. This removes a 1e7 loss of precision in the monomial basis at large `n`.
* J = 20 and 24 agree with the CHMRS series to 1e−10 – 1e−12.
* J = 60 is out of reach by this construction.

**Running when this was written:**
* the high-J sample (`J = 11…28`, 15 values of `m/Λ` in `[1.11, 3]`, `Δφ = 3`; `cache/e1_100_Df3_hiJ.pkl`);
* `Δφ = 2.5` to `m/Λ = 3`;
* the `e2`/`e4` checks.

`e3_action.py` merges the high-J file automatically. Until it is in, `J > 10` rests on OMEGA_U's `κ ≥ 3.5` [M].

### Bound (`e6_bound.py`, same formula as `omegau/u7`) [D given all inputs]

| `n` | `\|ĝ₃\|²Λ⁸ ≤` | `√(t₄²+α_T²) ≤` |
|---|---|---|
| 1.08 | 276 | 4.2·10⁻⁵ |
| 1.09 | 296 | 4.4·10⁻⁵ |

* The same code reproduces OMEGA_U's values: `5.4·10⁻⁵` (1.15), `8.1·10⁻⁵` (1.28) and `1.14·10⁻⁴` (1.40).
* So the new bound is about 2× below OMEGA_U's central value.
* It is conditional on H1–H4 and H6–H12. **H5 is no longer needed for the skeleton.**

### `Δ_g = 300` and `1000`: not achieved

**Box contour (`skelexact2.py`).** Crossing the real `B̂` poles is done by circle residues; the contour avoids the Γ-pole rows. It is
validated against the `Δ_g = 100` results to ≤ 1e−6 (`t_cross`, `t_box`). At `Δ_g = 300` the small-`ν` values match the series to 1e−12.

**Why the large-`x` nodes fail.** At `x ≳ 0.6` there is an irreducible cancellation:
* the complex-residue sums behave like `Σ (iν/4)^j/j!`, a cancellation of about `e^{ν/4}`;
* the box's horizontal segment near the pole rows is about 1e19 larger than the result;
* a vertical line at `c ≈ ν` leaves about `10^{0.15ν}`, i.e. 1e43 at `ν = 297`.

**What would be needed.** A true steepest-descent contour with the `1/sin` expansion of the Γ-pole rows, or 50–60-digit arithmetic.
mpmath is estimated at about 1 h per state, and installing python-flint was not permitted. **So `n(Δ_g)` beyond 100 is not determined.**

### NF02 / Ω^{NF,SC} (H13) (`e7_nf.py`, `cache/e7_100.txt`) [NE + I]

**On every scanned state the exact-skeleton NF02 action is positive, down to `m = Λ`.**
* States: `J ≤ 10`, `m ≤ 2.8Λ`, both skeletons, `b̂` up to about 21.
* Margins: `F/N_s ≥ +0.22`.
* `F_exact ≈ F_hyb`: the exact remainder is ≤ 0.02 `N_s`.
* The truncated series had given `κ < 0` there; that was an artefact of the series.

**What this covers, and what it does not.**
* For the scanned states, `n_exact = 1` for NF02.
* The open region of H13 (`b̂ > 20`, `J` up to hundreds) is not covered, apart from the high-J sample (`b̂ ≤ 37` at `m ≈ 1.5Λ`) once it is in.
* **So H13 is reduced but not closed.** The obstruction is the Mack-polynomial construction at large `J`.
