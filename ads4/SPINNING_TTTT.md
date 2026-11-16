# `⟨TTTT⟩`: helicity map, mean-field normalisation and the `O(1/m)` question for stress tensors; status of the doubly-spinning blocks

Status: written 16 Nov 2026 by a derivation/numerics agent (revised the same day: §3.5–3.6, §4.4–4.5); **not refereed**; not committed. Code and outputs: `ads4/code/fourpt/tttt/`
(caches in `tttt/cache/`, git-ignored with `fourpt/cache/`). Imports, does not modify, `ads4/code/e2e/`, `ads4/code/mft/`, `ads4/code/fourpt/`.
Builds on `SPINNING_4PT.md` (revised after `REFEREE_SPINNING_4PT.md`), `SPINNING_E2E.md`, `SPINNING_MFT.md`. CHL = Caron-Huot–Li,
arXiv:2102.08160.

Marks: **[D]** derived; **[DV]** derived and verified numerically; **[DVE]** derived and verified end to end on exact blocks; **[NE]**
numerical evidence (exact rationals where stated); **[ID]** closed form identified from exact data and checked at further points; **[Q]**
quoted; **[I]** inferred; **[ND]** not done.

## 0. Summary

1. **`⟨TTO⟩` helicity map, exact at rational points [NE, exact].** Two conserved parity-even structures (plus one `d = 3` relation among the 14
   differential-basis words); exact opposite- and same-helicity lines (momentum-space projection of CHL, calibrated on currents); CHL's
   `T^{E,h}` identified exactly in the conserved basis; the word-wise shadow map sends each helicity line exactly to the same helicity line at
   `(3−Δ,−1−J)` (shadow commutes with helicity for `l = 2`). Two points, `(97/10, 13/3)` and `(53/4, 11/2)`. §3.1–3.2.
2. **MFT normalisation for stress tensors is Weyl-even [D; exact Fock data, no input from CHL's MFT formulas].** (Revised: this used to rest
   on CHL's quoted `l = 2` MFT factors.) CHL's helicity structures are now built **in position space with exact rational coefficients**
   (`chl_xspace.py`, term-wise inverse Fourier transform), which works at the double-twist points; the two-particle Fock computation for
   `[TT]_{N,J}` (`N ≤ 6`, 39 entries) gives exactly **`p^MFT[T^{E,h}] = A_h P_FK(1,1;Δ,J)/2304`** with `A_opp = (−J)_4/(J+1)_4`,
   `A_same = (Δ−1)_4/(2−Δ)_4` [ID], helicity-diagonal Gram matrices (`n ≥ 1`); for currents the same construction gives `A_h P_FK(1,1)/16`
   (23 entries), equivalent to the earlier closed forms. The shadow eigenvalue of CHL's structures obeys **`κ_h A_h = 2^{3J+3/2}`**: proved
   symbolically for `l = 1` (identity of rational functions in `(Δ,J)`, `t10`), exact (rational arithmetic) at four points for `l = 2`
   (`t8`); hence `p^MFT[e_h] = 8^J × (scalar factor) × (Weyl-even)` for the Weyl-covariant helicity families. **The MFT normalisation of
   the graviton-channel helicity couplings adds no `O(1/m)` term;** with the reflection argument: no `O(1/m)` in the leading heavy
   `⟨TTTT⟩` Regge moments for parity-even exchanges (conditional on Assumption T and on the `l = 2` identity off the checked points; not end
   to end). §3.3, §3.5.
3. **Currents (goal 4 for `l = 1`) [D, exact Fock data; DV at block level].** Exact two-particle Fock MFT data in our normalisation, helicity-
   diagonal; absolute constant `1/4` per helicity quadruple, `C^± = (2P/m^{2k})(|c_same|² ± |c_opp|²)`, all factors Weyl-even; leading-twist block
   decomposition of `⟨J_1J_3⟩⟨J_2J_4⟩` reproduced with no free constant (`10⁻⁶`). Physical sign of the opposite-sign term `+(|c_same|²−|c_opp|²)` =
   flat `σ_{++} = +1`, `σ_{+−} = −1` (`SPINNING_4PT.md` §4.5–4.6). §2.
4. **`l = 2` heavy vertices and absolute constant [NE + exact inputs].** CHL's structures are helicity-pure at the heavy vertex. The leading
   two-pair coefficient is `T^{hLhR} = 8√(λ_{hL}λ_{hR})/m⁸`. The constant extrapolates to `E = 7.9996` along `β = 0.4`, with
   `m = 22 … 179` and `O(1/m²)` corrections. Combined with the exact Fock data, this gives **`1/4` per helicity quadruple for stress tensors**,
   with the factor `A_h` cancelling. Hence `C^± = (2P/m^{2k})(|c_same|² ± |c_opp|²)` for `⟨TTTT⟩` [D given `E = 8`]. §3.6.
5. **`⟨TTTT⟩` end to end [DVE at the percent level, `β = 1`, two masses].**
   * The partial waves are computed by **numerical Taylor-mode application of the `D`-words** (`taylor4.py`), validated to `5·10⁻¹⁵`.
   * The exact tensor Regge moments (`k = 6`, `ν = 0`, 625 modes) of `W[X_{hL},X_{hR}]` were computed at `β = 1`, `m = 39.3`
     (`Δ = 45.86`, `J = 20`) and `m = 70.0` (`Δ = 80.48`, `J = 36`).
   * Against the factorised prediction, with no free constant, at `m = 40 → 70`:
     - same–same `1.083 → 1.024`; opp–opp `0.957 → 1.002`; mixed `0.940 → 0.982`;
     - `Σ_+` `0.973 → 1.011`; the opposite-sign term **`Σ_−` `+0.92 → +0.97`**, so its sign, as for currents, is confirmed;
     - fit residuals `2–71 % → 0.7–7 %`.
   * Prediction-free checks: factorisation `0.85 → 0.94`; `T^{oo}/T^{ss}` against `λ_o/λ_s`: `2.80/3.15 → 1.859/1.900`; `T^{so} = T^{os}`
     exactly, so the Weyl-odd `O(1/m)` term of the integrand cancels in the moments.
   * The deviations shrink by `3.3–3.5` between the two masses (`m²` ratio `3.17`), consistent with `1 + O(1/m²)` [NE, two masses only].
     §4.4–4.5.
6. **Not done:** closed forms of the `l = 2` helicity lines in the differential basis (only the position-space form of `chl_xspace`);
   an exact proof of `E = 8`; parity-odd exchanges; `k = 2, 3`; more masses and the `1/m` fits.

## 1. Tools (all new, `tttt/`)

* **Structures at infinity** (`infframe.py`): `⟨O_1(x)O_2(0)X(∞,ζ)⟩` for every differential-basis word, from the embedding point at infinity
  `P_∞ = (0,1,0)`, `Z_∞ = (0,0,ζ)` (`s_{i5} = −1`, `Z_i·P_∞ = 0`, `Z_∞·P_i = ζ·x_i`), first in the invariants `a·x, b·x, ζ·x, x², a·b, a·ζ, b·ζ`.
  Checked: the finite-point words converge to these values as `x_5 → ∞` (`O(1/L)`, 150-digit test).
* **Conservation at infinity** (`infcons.py` symbolic in `(Δ,J)`; `infcons2.py` exact at rational points by the chain rule in the invariants,
  harmonic projection for spin 2, evaluation on exact 3d configurations with Gaussian-rational null vectors). Reproduces `v_1, v_2` of
  `SPINNING_4PT.md` exactly for `l = 1`. (A bug — one random scale per *component* of the null vectors, which destroys nullness — was found
  and fixed; results before the fix are discarded.)
* **Momentum-space helicity projection** (`ft_helicity.py`): term-wise Fourier transform `FT[(ζ·x)^e|x|^{2t}] = (−i)^e π^{3/2}2^{2t+3+e}
  Γ(t+e+3/2)/Γ(−t) (ζ·p)^e|p|^{−2t−3−2e}` (harmonic polynomial, Bochner), `(a·x) → i a·∂_p`; evaluated at `p = ẑ` with transverse null `a, b`;
  exact version with the Gamma functions reduced to Pochhammer symbols. Calibration on `l = 1`: the lines `T(+,−) = 0` and `T(+,+) = 0`
  are exactly `ĥ_opp` and `ĥ_same` of `SPINNING_4PT.md` §2.3 (exact rationals at two points).
* **CHL structures in momentum space** (`chl_momentum.py`): eq. `Tpm` with the E normalisation of eq. `even odd`, compared with the Fourier
  transform of our structures (frames: `x_CHL = −x`, so `FT_CHL(p) = FT_ours(−p)`).
* **Two-particle Fock space** (`fock2.py`): exact MFT data of `[O O]_{N,J}`, the method of `SPINNING_MFT.md` §2 with two spinning particles
  (`mft/fock.py` one-particle Gram forms).

## 2. Currents: exact MFT data, absolute constant, block-level test (coordinator goal 4 for `l = 1`)

See `SPINNING_4PT.md` §4.6 (written there as part of the referee revision). In short [D from exact Fock data]:
`p^MFT[ĥ_opp] = 8^J P_FK(1,1;Δ,J)/(8Y²(J−1)J(J+1)(J+2))`, `p^MFT[ĥ_same] = 8^J P_FK(1,1;Δ,J)/(8Y²Δ(Δ−1)(Δ−2)(Δ−3))` (all 25 primaries with `n ≥ 1`, `N ≤ 6`,
helicity-diagonal Gram matrix; `n = 0`: one primary carrying the sum); absolute constant `1/4` per helicity quadruple, `C^± =
(2P/m^{2k})(|c_same|² ± |c_opp|²)`, `Q_h = 1 + O(1/m²)` (`out_t2_mft_heavy_l1.txt`); block-level leading-twist test of `⟨J_1J_3⟩⟨J_2J_4⟩` with the exact
doubly-spinning waves: ratio `1 + O(ε)` (`10⁻⁶` at `ε = 10⁻⁵`), no free constant, control `12.6` (`out_t3_lightcone_l1.txt`) [DV].
Physical sign of the opposite-sign term: `+(|c_same|² − |c_opp|²)` for the covariant odd⊗odd kernel (`SPINNING_4PT.md` §4.5, `out_t4_sign_l1.txt`).

## 3. Stress tensors: exact helicity map and the Weyl parity of the MFT normalisation (coordinator goal 2)

**3.1 Conserved space and exact helicity lines [NE, exact] (`t5_helicity_l2.py`, `out_t5_helicity_l2_pointA.txt`, `out_t5_helicity_l2_53_4_11_2.txt`).**
At `(Δ,J) = (97/10, 13/3)` (and `(53/4, 11/2)`, §3.3): conservation at `x_1` and `x_2` in the frame at infinity leaves a 3-dimensional null space
of the 14 words, containing the single `d = 3` relation among the words (found exactly at infinity as well): **two conserved parity-even
`⟨TTO⟩` structures**. The momentum-space projection (calibrated on currents, §1) gives exactly one *opposite*-helicity line (`T(+,−) = 0`) and
one *same*-helicity line (`T(+,+) = 0`) modulo the relation, with exact rational coefficients.

**3.2 Weyl covariance and CHL identification [NE, exact / 30 digits].** For each helicity line `e_h(Δ)`, the word-wise shadow image `T e_h(Δ)`
lies exactly on the same helicity line at `(3−Δ,−1−J)` (the wrong components vanish identically): the shadow map commutes with helicity for
stress tensors too, so `e_h(Δ̃) := T e_h(Δ)` is a Weyl-covariant helicity family. CHL's `T^{E,same}`, `T^{E,opp}` (momentum space, eqs. `Tpm`,
`even odd`) are proportional to `e_same`, `e_opp` with the same ratio `N_h` on both non-zero components — the exact identification of CHL's
`⟨TTO⟩` structures in the differential basis at this point.

**3.3 The Weyl parity of the MFT normalisation [D from Q + exact checks].** With `T^{E,h} = N_h e_h` and the Weyl-covariant `e_h`, the shadow acts on
CHL's structure as `T^{E,h}(Δ) → κ_h T^{E,h}(Δ̃)`, `κ_h = N_h(Δ)/N_h(Δ̃)`. CHL's MFT data are `p[T^{E,h}] = c(1,1)^s(Δ,J)·A_h·(const)` [Q, eq. `free data`],
`A_opp = (−J)_{2l}/(J+1)_{2l}`, `A_same = (Δ−1)_{2l}/(2−Δ)_{2l}`. Then `p[e_h] = N_h² p[T^{E,h}]` and

  `p[e_h](Δ,J)/p[e_h](3−Δ,−1−J) = κ_h² A_h² × [c(1,1)^s(Δ)/c(1,1)^s(Δ̃)]`,

so the normalisation of the Weyl-covariant family is Weyl-even (up to the helicity-independent scalar factor, which is that of the scalar
problem) iff `κ_h A_h` is helicity-independent. Result (`out_t6_kappa_l1.txt`, `out_t5_helicity_l2_*.txt`):

| `l` | `(Δ, J)` | `κ_opp A_opp` | `κ_same A_same` | `2^{3J+3/2}` |
|---|---|---|---|---|
| 1 | (97/10, 13/3) | 23170.4750059 | 23170.4750059 | 23170.4750059 |
| 1 | (53/4, 11/2) | 262144 | 262144 | 262144 |
| 2 | (97/10, 13/3) | 23170.4750059 | 23170.4750059 | 23170.4750059 |
| 2 | (53/4, 11/2) | 262144 | 262144 | 262144 |

(ratio `κ_opp A_opp/κ_same A_same = 1` to `10⁻³⁶` at the `l = 2` point A, to `10⁻¹⁶` at point B). So **`κ_h = 2^{3J+3/2}/A_h` exactly, for currents and for stress tensors**:
the helicity dependence of CHL's shadow factor is exactly the inverse of their MFT helicity factor, and the remainder `2^{3J+3/2} =
2^{3J/2}/2^{3J̃/2}` is the reflection ratio of the convention factor `8^{J/2}` — the same `8^J` that appears in the exact `l = 1` Fock data
(`SPINNING_4PT.md` §4.6), where it is cancelled by the block normalisation `(−1/8)^J`. Hence `p^MFT[e_h] = 8^J × (scalar MFT factor) × (Weyl-even)`
for `l = 2`, as proved for `l = 1` by the exact Fock data. **Consequence: the MFT (unitarity) normalisation of the stress-tensor helicity
couplings introduces no `O(1/m)` term**, and with the reflection argument of `SPINNING_4PT.md` §6.1 (which holds for any `l` once the families
are Weyl-covariant and helicity-pure — now shown) the leading heavy `⟨TTTT⟩` Regge moments with MFT-normalised helicity couplings have
corrections `1 + O(1/m²)` for parity-even exchanges.
Status (revised, §3.5): (a) the MFT helicity factors `A_h` are now derived from our own exact Fock data for `l = 1` and `l = 2` [ID]; (b)
the identity `κ_h A_h = 2^{3J+3/2}` is proved symbolically for `l = 1` and holds exactly (rational arithmetic) at four points for `l = 2`
(`out_t8_kappa_exact.txt`; a symbolic `l = 2` run did not finish, §3.5); (c) Assumption T. Not verified end to end (no `⟨TTTT⟩` moments, §4).

**3.4 What is not done for `l = 2`.** Closed forms of the helicity lines in the differential-basis words as functions of `(Δ,J)`; only the
position-space closed form of §3.5 exists. Also: parity-odd exchanges. (The heavy vertices and the absolute constant are now done: §3.6.)

**3.5 CHL's structures in position space; direct Fock computation for `[TT]`; the `κ_h A_h` identity (`chl_xspace.py`, `t7_fock_l2.py`,
`t7b_closed_forms.py`, `t8_kappa_exact.py`, `t10_kappa_symbolic.py`).**

*Construction [D].* CHL's `T^{h_1h_2}(p) ∝ (e_1Π_{−h_1p}e_3)^l (e_2Π_{h_2p}e_3)^l (p·e_3)^{J−2l}|p|^{β−3}` is expanded in `p·e_i`, `|p|`
(products of two `ε(e_i,p,e_3)` are Gram determinants), and each term is inverse-Fourier transformed with
`∫d³p/(2π)³ e^{ipx}(e_3·p)^n|p|^{2s} = i^n(2π)^{−3}π^{3/2}2^{2s+3+n}Γ(s+n+3/2)/Γ(−s)(e_3·x)^n|x|^{−2s−3−2n}`, `p·e_i → −i e_i·∂_x`.
The Gamma functions of all terms are reduced to Pochhammer ratios relative to one common factor; collecting every `(Δ,J)`-dependent
prefactor (CHL's `pref`, `G_E`, the transform), all `π`'s cancel and

  `T^{E,h} = C_l(Δ,J) X_h`,  `C_l = (−1)^l 2^{3l−1/2} 2^{J/2} ((2−Δ+J)/2)_{2l} ((Δ+J+1−4l)/2)_{2l}`,

with `X_h` a structure with coefficients rational in `(Δ,J)` (powers `(e_3·x)^{J+k}|x|^{…}`). Checks: for `l = 1`, `X_h` is exactly proportional
to our `ĥ_h` (symbolically: `ĥ_h = X_h·Δ(Δ−1)` resp. `J(J−1)` `×(Δ−J−4)(Δ−J−2)(Δ+J−3)(Δ+J−1)Y/2^{J+2}` in the CHL frame) and reproduces the
momentum-space `κ_h` of `t6` (all 12 printed digits). At double-twist points individual Pochhammer ratios have zeros (contact terms, dropped
correctly) or poles; the limit `Δ → 2l+2+N` is taken as the leading term of `X_h` in `Δ − Δ_N` (order `ε^{−1}` for `n < 2l`, compensated by the
simple zero of `((2−Δ+J)/2)_{2l}` in `C_l`).

*Fock data [D + ID].* The limit structures are polynomial (for `J ≥ 2`; the `l = 2` same-helicity structure is not polynomial for `J = 0, 1`,
where CHL's formula does not apply) and pass the two-particle range check (conservation) of `fock2.py`. Results (`out_t7_fock_l1.txt`,
`out_t7_fock_l2.txt`, `out_t7b_closed_forms.txt`), `n = (2l+N−J)/2`:
* Gram matrix **diagonal in helicity** for every `n ≥ 1` (`l = 2`: 9 points with both helicities, `N ≤ 6`); at `n = 0` the two structures are
  proportional (one primary), as for currents.
* **`p^MFT[T^{E,h}] = A_h P_FK(1,1;Δ,J)/K_l`, `K_1 = 16`, `K_2 = 2304 = 48²`** — exact at all 23 (`l = 1`) and 39 (`l = 2`) entries
  (`n = 0`: twice the formula, the two helicities carrying one primary). In particular `p[T^{E,opp}]/p[T^{E,same}] = A_opp/A_same` exactly:
  **CHL's `l = 2` MFT helicity factors are confirmed by an independent computation** and are no longer an input. For `l = 1` the closed form is
  equivalent to `p[ĥ_h] = 8^J P_FK(1,1)/(8Y²Π_h)` of `SPINNING_4PT.md` §4.6 (via the ratio above).
* (The constants `K_l` refer to CHL's normalisation conventions; only `A_h` and the `(Δ,J)` dependence enter the Weyl-parity argument.)

*The identity `κ_h A_h = 2^{3J+3/2}`.* `κ_h` is the eigenvalue of the word-wise shadow map on CHL's structure, `T T^{E,h}(Δ,J) = κ_h T^{E,h}(3−Δ,−1−J)`.
With `T^{E,h} = C_l X_h`: decompose `X_h` into the differential-basis words (for `l = 2` modulo the `d = 3` relation, which at infinity is the Gram
determinant of `(x, a, b, ζ)` times `x⁴/(ζ·x)²`), apply `T`, read off `λ_h` from `T X_h(Δ) = λ_h X_h(Δ̃)` (the residual must vanish
identically), `κ_h = λ_h C_l(Δ)/C_l(Δ̃)`. Since `C_l(Δ)/C_l(Δ̃) = 2^{J+1/2} ×` (rational) and the words carry `2^{−J}`, the identity reduces to an
identity between rational functions of `(Δ,J)`.
* `l = 1` [D, symbolic]: `λ_same/2^{2J+1} = (Δ−3)(Δ−2)(Δ−J)(Δ+J+1)/(Δ(Δ−1)(Δ−J−4)(Δ+J−3))`, `λ_opp/2^{2J+1} = (Δ−J)(J+1)(J+2)(Δ+J+1)/(J(J−1)(Δ−J−4)(Δ+J−3))`,
  and `κ_h A_h/2^{3J+3/2} = 1` identically (`out_t10_kappa_symbolic_l1.txt`).
* `l = 2` [NE, exact]: `κ_h A_h/2^{3J+3/2} = 1` exactly (rational arithmetic, algebraic `2^{J}` tracked) for both helicities at `(97/10,13/3)`,
  `(53/4,11/2)`, `(61/7,29/5)`, `(127/10,47/7)` (`out_t8_kappa_exact.txt`); `T X_h(Δ) ∝ X_h(Δ̃)` holds exactly at each point (shadow commutes
  with helicity). Symbolic `l = 2`: not completed. The word decomposition over `Q(Δ,J)` (`t10_kappa_symbolic.py 2`) did not finish its first step within the session's limits (about 2 h, or killed for memory). So the `l = 2` identity rests on the four exact points [NE], not on a proof.

**3.6 Heavy vertices of the `l = 2` helicity structures and the absolute MFT constant (`t12_vertex_l2.py`, `t14_vertex_fit_l2.py`,
`t15_mft_constant_l2.py`; outputs `out_t12_vertex_l2_b0.4.txt`, `out_t14_vertex_fit_l2.txt`, `out_t15_mft_constant_l2.txt`).**

*Method.* This is the same as for `l = 1` (`fourpt/vertex2`, `t3_vertex_l1`).
* The exact word coefficients of `X_h(Δ,J)` are used on the left vertex. The right vertex uses `T X_h(Δ) = λ_h X_h(Δ̃)` (exact, §3.5).
* The shadow integrand is evaluated at fixed `(X,Y) = m(x,y)` (`fourpt/t7_vertex_l2.helicity`; 30 digits), on the ray
  `A_± = s A_±⁰` with `A_+⁰ = 22.79`, `A_-⁰ = 4.5` (`β = 0.40`) and `s = 1, 2, 4, 8` (`m = 22.3 … 178.7`, `J = 4 … 35.5`).
* One point on a ray with `β = 0.9` (`m = 44.2`) was also computed.
* Control: the same pipeline for `l = 1` reproduces the closed-form vertices of `predict4` (ratios `1.003`–`1.01` on the diagonal waves at
  `m = 22, 44`, with the known antisymmetric `O(1/m)` term on the mixed ones).

*Results.*
* **Helicity purity** [NE]: `X_same` couples only `h_2 = −h_1`, and `X_opp` only `h_2 = h_1`. Wrong-helicity components are `10⁻⁵…10⁻⁹` of the
  leading ones at `m = 22`.
* **Factorisation**: `T^{hLhR} = V_{hL}V'_{hR}` holds exactly.
* **Leading vertex [ID, numerical]**: `T^{hh} = λ_h E/m⁸`, with `E` helicity-independent and `E → 8`:
  - `m⁸T^{hh}/λ_h = 10.24, 8.495, 8.120, 8.030` (same) and `10.33, 8.514, 8.125, 8.031` (opp) at `m = 22.3, 44.7, 89.4, 178.7`;
  - `m²(E/8 − 1) = 139.6, 123.4, 119.9, 119.0`, i.e. `O(1/m²)` corrections;
  - Richardson extrapolation gives `E(∞) = 7.9996`.
* **Asymmetry of the mixed waves.** `T^{so}/T^{os} = 75.7, 6.60, 2.52, 1.58` along the ray, approaching 1. It tracks `λ_o/λ_s = 141.8, 8.94,
  2.92, 1.70`. This is the antisymmetric (Weyl-odd) `O(1/m)` term, as for `l = 1`. For stress tensors it is large at small `J`, because
  `A_opp = (−J)_4/(J+1)_4` is `0.014` at `J = 4`. The geometric mean is symmetric: `m⁸√(T^{so}T^{os}/(λ_sλ_o)) = 10.46 … 8.032`.
* Hence the leading two-pair coefficient is

  **`T^{hLhR}_h = 8 √(λ_{hL}λ_{hR})/m⁸`** on the helicity classes of `(hL, hR)`, with `λ_h = 2^{2J+1}P(Δ̃)/(A_hP(Δ))` and
  `P(Δ) = ((2−Δ+J)/2)_4((Δ+J−7)/2)_4`.

  Equivalently, the Weyl-covariant family `e_h ∝ X_h/√λ_h` has the Weyl-even vertex `V_e² = 8/m⁸`, the `l = 2` analogue of
  `V_ĥ ∝ A_±²m²(A_+²+A_−²)/2`.

*Absolute constant.* With the exact Fock data `p[X_h] = |C_2|²A_hP_FK(1,1)/2304` (§3.5):
* `p[X_h]8^{−J}T^{hh} = 2^{12}P(Δ)P(Δ̃)P_FK(1,1)E/(2304m⁸)`. **`A_h` cancels**, so the normalisation is helicity-independent and Weyl-even.
* The algebra of `t2_mft_heavy_l1` with `l = 2` (`a_{3,3}`, `4^6`, `m^{−10}`) gives `Q_h = (E/8)(1 + O(1/m²))`. At `E = 8` the result is
  `Q_h = 0.98837, 0.99988, 0.9999988` at `m = 102, 10³, 10⁴` (`β = 0.4`), and similarly at `β = 1, 2` (`out_t15`).
* So **the constant is `1/4` per helicity quadruple for stress tensors too** [NE: the exact Fock data and the exact structures enter;
  `E = 8` is a numerical limit to `5·10⁻⁵`, not a proof].
* Consequently **`C^± = (2P_{−½+iν}(cosh β)/m^{2k})(|c_same|² ± |c_opp|²)`** for `⟨TTTT⟩` with MFT-normalised helicity couplings [D given
  `E = 8`]. The sign of the opposite-sign term follows from the same covariant kernel argument as for currents (`SPINNING_4PT.md` §4.5) [I for
  `l = 2`; tested end to end in §4.5].

## 4. `⟨TTTT⟩` blocks and moments (coordinator goals 1, 3)

**4.1 Symbolic word pairs: infeasible [ND].** Measured cost of the symbolic word pairs at jet order 8 (`feas_l2.py`, `cache/feas_l2.log`):
`HH × HH` 9 terms (0 s), `AACC × HH` 22 800 (45 s), `ACH × ACH` 73 095 (171 s), `AACC × AACC` not finished after 67 min. The 196 pairs are mostly of
the last type.

**4.2 Taylor-mode application of the operators [D; validated] (`taylor4.py`, `t9_taylor_validate.py`, `out_t9_taylor_validate.txt`).** The
operators act on functions of the 24 dot products (`s_ij`, `Z_i·P_j`, `Z_i·Z_j`; treated as independent, as in `embed2`). Every `D` is a first-
order operator `Σ (∂σ/∂W·V) ∂/∂σ` with polynomial coefficients. The state is represented as a polynomial in the 18 `z`'s (exact exponents) whose
coefficients are truncated Taylor series in the 6 deviations `δs_ij`. The order needed is the number of `D`'s in the pair (`≤ 8`). Each `D`
consumes one order, so the series are truncated progressively (the `Z`-derivative terms are truncated with the others, since at most the
remaining number of `D`'s still act). The seed `K̃ g̃(P,Q)` is built from its `(P,Q)`-jet by series composition (`log`, `exp`, powers). The jets
of the 25 seed blocks to order 8 come from `cblock12.c`, a copy of `e2e/cblock.c` with the jet arrays enlarged (`kmax, lmax ≤ 11`), and from
`pq_jets_M`, the any-order version of `e2e/spinblock.pq_jets` (agrees with it to `3·10⁻¹⁵` at order 4).
Validation with a synthetic seed `g̃ = e^{aP+bQ}(1+cPQ²)` and random values of all 24 symbols:
* against the symbolic `word4` expressions: all 25 `l = 1` pairs (worst `5.3·10⁻¹⁵`); `l = 2` pairs `HH×HH`, `HH×ACH`, `BEH×HH`, `AACC×HH`,
  `HH×ABCE` (`≤ 1.0·10⁻¹⁵`);
* order-8 pairs, which have no symbolic reference: `D^{(12)}D^{(34)} = D^{(34)}D^{(12)}` (the two orderings exercise the truncation
  differently): `AACC×AACC` `1.1·10⁻¹⁵`, `ABCE×BBEE` `4.1·10⁻¹⁵`, `ACH×AAEE` `5.9·10⁻¹⁶` (4539 `z`-monomials at order 8).

**4.3 Feasibility at a heavy point [NE] (`t11_feas_taylor.py`, `out_t11_feas_taylor.txt`).** Case `β = 1`, `m ≈ 40` (`Δ = 45.861`, `J = 20`), two
points of the moment grid of `e2e/moments.py`, `N = 8000` radial terms. The 196 pairs are contracted with CHL's **exact** helicity structures on
both vertices: word coefficients of `X_h` from the exact decomposition of §3.5, 88 s per `(Δ,J)`. They are evaluated on the `5⁴` polarisation grid
of `moments4.pol_grid4`.

| | point 0 (`U = 0.054`, `S = 0.21`) | point 1 (`U = 0.58`, `S = 0.69`) |
|---|---|---|
| order-8 jets, 25 seeds (mp composition) | 4.1 s | 3.6 s |
| 196 pairs, float64, one point | 27 s (slowest `AACC×AACC` 1.9 s) | 22 s |
| 196 pairs, long double, one point | 31 s | 32 s |
| cancellation `Σ|terms|/|W|` (same-same, same-opp, opp-opp) | `1.5·10⁵`, `1.3·10⁷`, `2.6·10⁶` | `1.3·10⁵`, `7.7·10⁵`, `3.2·10⁶` |
| float64 vs long double, `W[s,s]`, `W[s,o]`, `W[o,o]` | `9·10⁻⁸`, `3·10⁻⁵`, `2.5·10⁻⁵` | `1·10⁻¹⁰`, `4·10⁻⁹`, `2·10⁻⁸` |

Vectorised over 200 points (`AACC×AACC`, order 8): seed series 13 s (float64) / 41 s (long double); operator application 0.6 s / 1.4 s.
Hence the per-point cost is about 1.7 s for the 25 seeds plus about 0.6 s for the 196 pairs (float64), and about 5 s + 3 s in long double. On the
moment grid (≤ 4464 points) that is about 3 h (float64) or 10 h (long double) per `(Δ,J)` on one core. The seed composition dominates and could be
made much cheaper (e.g. series in `log s`). Precision: the conserved helicity combinations cancel by up to `10⁷` already at `m ≈ 40`, so float64
gives `10⁻⁵`-level accuracy on the mixed-helicity waves. That is marginal for `O(1/m)` fits at `m ~ 100`; long double gives `~10⁻¹²` [I].

**4.4 Production code (`tmoments_l2.py`).** The production code adds the following:
* a vectorised order-8 jet composition (`pq_jets_vec`, agrees with the mp version to `10⁻¹⁵`);
* the 25 seed blocks from `cblock12` (double radial sums, `N = 8000`);
* one seed series per seed block, truncated per word pair;
* the Taylor-mode application of all 196 pairs, contracted on the fly with CHL's exact helicity word vectors;
* the `5⁴` polarisation grid and the Fourier modes of `fourpt/tmoments4`.

The moments and kernels are those of `fourpt/tmoments4` (`t13_analysis_l2.Modes`).

A bug was found and fixed during this work: the derivatives of `arccosh` had been lambdified in expanded form, `P_k(c)/(c²−1)^{k−½}` with
`(c²−1)^k` multiplied out. That cancels catastrophically near `c = 1` (small `S`). It is now evaluated as `P_k(c)/sh^{2k−1}` with
`sh = (1/r−r)/2`, and agrees to `10⁻¹²` at `r = 0.9995`. The `t11` feasibility numbers are unaffected (their points have `r ≤ 0.64`).

Cost on the shared machine (load average ≈ 17 on 4 cores): `β = 1`, `m ≈ 40` (`Δ = 45.861`, `J = 20`, 3336 grid points) took 8150 s wall
(about 1.5–2 s per point when not starved). Peak memory was 3.8 GB. A concurrent second run was killed by the out-of-memory handler, so the
runs were done sequentially.

**4.5 End-to-end `⟨TTTT⟩` check (`t13_analysis_l2.py`, `out_t13_analysis_l2_*.txt`) [DVE: 5–8 % at `m = 39`, 0.2–2.6 % at `m = 70`].**

*Setup.*
* `k = 6`, `ν = 0`, all 625 modes `q = (q_1,q_3,q_2,q_4)`.
* Waves `W[X_{hL},X_{hR}]` for the four helicity pairs.
* A mixed wave `W[e,e]` with `e = X_same + wX_opp`, where `w` is chosen so that the predicted `T_oo = 2T_ss`. This makes the opposite-sign
  combination `Σ_− = T_oo − T_ss` non-zero and positive, so its sign is tested.
* Prediction: `Pref Σ_h T_h M M 𝒦` with `T^{hLhR} = 8√(λ_{hL}λ_{hR})/m⁸` (§3.6). It contains no fitted constant.
* Fits: one-scale fits (uniform weights) and 8-class fits (uniform weights). The doubly flipped classes have tiny design columns and are not
  determined.

| `β = 1`, `m = 39.3` (`Δ = 45.861`, `J = 20`) | one-scale fit (resid.) | dominant class fit/pred | 8 largest predicted modes, meas/pred |
|---|---|---|---|
| `W[X_s,X_s]` | `1.083` (`2.2 %`) | `ss`: `1.082` | `1.080 … 1.085` |
| `W[X_o,X_o]` | `0.957` (`20 %`) | `oo`: `0.962` | `0.954 … 1.000` |
| `W[X_s,X_o]` = `W[X_o,X_s]` (identical, as required by the pair exchange) | `0.940` (`71 %`) | `pf/fp`: `0.940` | `0.992 … 1.024` |
| mixed `W[e,e]` | `0.969` (`18 %`) | `Σ_+`: `0.973`, **`Σ_−`: `0.922`** | `0.954 … 1.000` |

Prediction-free checks:
* `T^{so}T^{os}/(T^{ss}T^{oo}) = 0.85` (fitted dominant classes).
* `T^{oo}/T^{ss} = 2.80`, against `λ_o/λ_s = 3.15`. This tests the helicity dependence `A_h` of the MFT/shadow normalisation.
* `T^{so} = T^{os}` exactly. The integrand is not symmetric at this mass (`T^{so}/T^{os} ≈ 2.4` pointwise, §3.6), so its Weyl-odd
  `O(1/m)` part cancels in the moments, as the reflection argument requires.

Reading:
* The diagonal waves and the dominant modes of all waves agree with the factorised prediction to 1–8 %, with no free constant. This
  includes the absolute normalisation (`E = 8`, i.e. the `l = 2` heavy vertex and, through `p[X_h]`, the `1/4` per quadruple) and the sign
  of the opposite-sign term (`Σ_−` fit/pred `= +0.92`; a sign error would give `−1`).
* The 20–70 % residuals come from modes that are predicted small, in particular preserved–preserved components of the mixed-helicity waves,
  which are predicted zero. They are consistent with the large Weyl-odd `O(1/m)` terms of `l = 2` at `J = 20`, where
  `λ_o/λ_s − 1 ≈ 2`, but at one mass this is [I].
* Precision: float64 is accurate to `≤ 3·10⁻⁵` pointwise at this mass (`t11`), negligible here.

| `β = 1`, `m = 70.0` (`Δ = 80.484`, `J = 36`) | one-scale fit (resid.) | dominant class fit/pred | 8 largest predicted modes, meas/pred |
|---|---|---|---|
| `W[X_s,X_s]` | `1.024` (`0.7 %`) | `ss`: `1.024` | `1.024` (all 8) |
| `W[X_o,X_o]` | `1.002` (`1.5 %`) | `oo`: `1.0015` | `0.999 … 1.002` |
| `W[X_s,X_o]` = `W[X_o,X_s]` | `0.982` (`7.4 %`) | `pf/fp`: `0.982` | `0.973 … 0.984` |
| mixed `W[e,e]` | `1.008` (`1.6 %`) | `Σ_+`: `1.011`, **`Σ_−`: `0.974`** | `0.999 … 1.003` |

Prediction-free checks at `m = 70`: factorisation `T^{so}T^{os}/(T^{ss}T^{oo}) = 0.94` (it was `0.85` at `m = 40`); `T^{oo}/T^{ss} = 1.859`, against
`λ_o/λ_s = 1.900`; `T^{so} = T^{os}`.

Cost and precision at `m = 70`: 3475 points in 7520 s wall (chunks of 100 points; about 2 s per point). In float64 against long double at one
grid point (`t11`, `out_t11_feas_taylor_m70.txt`), the cancellation is `2.4·10⁵ … 1.2·10⁷` and the deviation is `≤ 1.3·10⁻⁵`.

**Between the two masses:**
* The deviations from the prediction shrink as follows: same–same `+8.3 % → +2.4 %` (factor 3.5); mixed `−6.0 % → −1.8 %` (factor 3.3);
  opp–opp `−4.3 % → +0.2 %`; `Σ_−` `−7.8 % → −2.6 %`.
* The residuals drop from 20–71 % to 0.7–7 %. The `m²` ratio is `3.17`, consistent with `1 + O(1/m²)` and no `O(1/m)` [NE, two masses; no fit].

**Marks.**
* [DVE] at the percent level (`m = 70`), with no free constant:
  - the factorised `⟨TTTT⟩` moment formula with the `l = 2` vertices `8√(λλ')/m⁸`;
  - the helicity dependence `λ_h ∝ 1/A_h`;
  - the absolute normalisation;
  - the sign of the opposite-sign term `Σ_−`.
* [I]: the `1/m²` approach, from two masses only, and the attribution of the `m = 40` residuals to Weyl-odd terms.
* Not tested: the MFT-normalised `C^±` itself, i.e. MFT sums over operators. This rests on the exact Fock data plus `E = 8` (§3.6).

## 5. Status

| # | statement | mark |
|---|---|---|
| 1 | Structures at infinity, conservation there (symbolic `l = 1`, exact at rational points `l = 2`), momentum-space helicity projection | [D], checks [DV] (limit `x_5 → ∞`; `l = 1` reproduces `v_1, v_2`, `ĥ_h`) |
| 2 | `⟨JJ⟩` MFT data in the helicity basis, closed forms, absolute constant `1/4`, Weyl-even | [D from exact Fock data] |
| 3 | Two-pair block normalisation `p(−1/8)^J W` (leading twist) | [DV] |
| 4 | Physical sign `+(|c_same|²−|c_opp|²)` of the opposite-sign term (`l = 1`) | [D] + DVE; functional ↔ flat sum rule [I] |
| 5 | `⟨TTO⟩`: 2 conserved structures + 1 relation; exact helicity lines; CHL `T^{E,h}` identified; shadow commutes with helicity | [NE, exact] at two points |
| 6 | `κ_h = 2^{3J+3/2}/A_h` (`l = 1, 2`) | `l = 1`: [D] symbolic; `l = 2`: [NE, exact] at four points |
| 7 | MFT normalisation of `l = 2` helicity couplings Weyl-even ⇒ no `O(1/m)` from it | [D] from our Fock data + row 6 |
| 7a | `p^MFT[T^{E,h}] = A_h P_FK(1,1)/K_l`, helicity-diagonal (`l = 1, 2`) | [ID], exact Fock data, 23 + 39 entries |
| 8 | No `O(1/m)` in leading heavy `⟨TTTT⟩` moments (parity-even exchanges, MFT-normalised) | [D] conditional on 7 and Assumption T; not end to end |
| 9 | `⟨TTTT⟩` exact partial waves by Taylor mode (operator application) | [D], validated to `10⁻¹⁵`; feasibility at `m ≈ 40` [NE] |
| 9a | `⟨TTTT⟩` end-to-end moments vs factorised prediction, both `±` (`k = 6`, `ν = 0`) | [DVE] at 5–10 % (dominant modes), two masses `m = 39, 70`: deviations `8 % → 2.4 %`, `O(1/m²)`-like; `1/m` fits [ND] |
| 9b | `l = 2` heavy vertex `T^{hLhR} = 8√(λλ')/m⁸`; `1/4` per quadruple; `C^±` for `⟨TTTT⟩` | [NE] (`E = 8` extrapolated to `5·10⁻⁵`), exact inputs [D/ID] |
| 10 | parity-odd exchanges; `k = 2, 3`; proof of `E = 8` | [ND] |

## 6. Code index (`ads4/code/fourpt/tttt/`)

| file | content | output |
|---|---|---|
| `infframe.py`, `infcons.py`, `infcons2.py` | structures at infinity; conservation (symbolic / exact at rational points) | — |
| `ft_helicity.py`, `chl_momentum.py` | momentum-space helicity projection; CHL `T^{E,h}` | — |
| `fock2.py`, `mft2.py`, `t1_mft_l1.py` | two-particle Fock space; exact MFT data of `[JJ]` in the helicity basis | `out_t1_mft_l1.txt` |
| `t2_mft_heavy_l1.py` | absolute two-pair heavy constant, `l = 1` | `out_t2_mft_heavy_l1.txt` |
| `t3_lightcone_l1.py` | block-level test of the two-pair MFT normalisation | `out_t3_lightcone_l1.txt` |
| `t4_sign_l1.py` | physical sign of the opposite-sign term | `out_t4_sign_l1.txt` |
| `t5_helicity_l2.py` | `⟨TTO⟩` helicity map, Weyl parity of the MFT normalisation | `out_t5_helicity_l2.txt` |
| `t6_kappa_l1.py` | calibration `κ_h A_h` for currents | `out_t6_kappa_l1.txt` |
| `feas_l2.py` | feasibility of the symbolic `⟨TTTT⟩` word pairs | `cache/feas_l2.log` |
| `taylor4.py`, `cblock12.c` | Taylor-mode application of the `D`-words; any-order `(P,Q)`-jets; block library with order-11 jets | — |
| `t9_taylor_validate.py` | Taylor mode vs symbolic word pairs; vertex commutation at order 8 | `out_t9_taylor_validate.txt` |
| `t11_feas_taylor.py` | `⟨TTTT⟩` feasibility at a heavy point: timing, cancellation, float64 vs long double | `out_t11_feas_taylor.txt` |
| `t12_vertex_l2.py`, `t14_vertex_fit_l2.py`, `t15_mft_constant_l2.py` | `l = 2` heavy vertices of `X_h`; `E → 8`; absolute MFT constant | `out_t12_vertex_l2_b0.4.txt`, `out_t14_vertex_fit_l2.txt`, `out_t15_mft_constant_l2.txt` |
| `tmoments_l2.py`, `t13_analysis_l2.py` | Taylor-mode `⟨TTTT⟩` waves on the moment grid; end-to-end comparison | `cache/prod_*.log`, `out_t13_analysis_l2_*.txt` |
| `chl_xspace.py` | CHL `T^{E,h}` in position space, exact; double-twist limit | — |
| `t7_fock_l2.py`, `t7b_closed_forms.py` | Fock data of `[JJ]`, `[TT]` in CHL's helicity basis; closed forms | `out_t7_fock_l1.txt`, `out_t7_fock_l2.txt`, `out_t7b_closed_forms.txt` |
| `t8_kappa_exact.py`, `t10_kappa_symbolic.py` | `κ_h A_h = 2^{3J+3/2}`: exact at points / symbolic | `out_t8_kappa_exact.txt`, `out_t10_kappa_symbolic_l*.txt` |
