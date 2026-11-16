# Two spinning Regge pairs: `⟨JJJJ⟩` end-to-end, the `±` helicity structure, the shadow matrix of conserved structures, and the `O(1/m)` question

Status: written 16 Nov 2026 by a derivation/numerics agent; **revised after the referee report `ads4/REFEREE_SPINNING_4PT.md`** (no fatal
issues; changes listed in §11). New computations for the revision are in `ads4/code/fourpt/tttt/` (`t1_mft_l1.py`, `t2_mft_heavy_l1.py`,
`t3_lightcone_l1.py`, `t4_sign_l1.py`). Code and recorded outputs: `ads4/code/fourpt/` (index §10;
`sh run_all.sh`; caches in `ads4/code/fourpt/cache/`, git-ignored). Imports, does not modify, `ads4/code/e2e/`. Conventions, frame,
functionals and the one-pair results are those of `SPINNING_E2E.md` (§1, §3, §4) and `SPINNING_MFT.md`; CHL = Caron-Huot–Li,
arXiv:2102.08160 (LaTeX in the session scratchpad; labels quoted).

Marks: **[D]** derived; **[DV]** derived and verified numerically; **[DVE]** derived and verified end-to-end on exact blocks;
**[NE]** numerical evidence (read off numbers); **[ID]** closed form identified from exact data (no continuous parameter fitted) and
then checked exactly at further points — not derived; **[Q]** quoted; **[I]** inferred; **[ND]** not done.

## 0. Summary

1. **`⟨JJJJ⟩` exact partial waves [DV].** Doubly-spinning conformal partial waves for a parity-even `(Δ,J)` exchange, both vertices the
   general conserved parity-even `⟨JJO⟩` (two structures each), built with the CPPR differential operators at both pairs on the `e2e`
   scalar blocks (25 word pairs). Validated: conservation at `x_1, x_2, x_3` (`≤ 2·10⁻⁸`, control `O(1)`), pair exchange (`10⁻¹⁴`), reduction of the
   `(H,H)` word pair to `H_12H_34 ×` scalar exchange (`10⁻¹⁵`), exact two-pair mode symmetries (`10⁻¹⁰`). §1–3.
2. **Shadow map and helicity basis, `l = 1` [D, exact].** The Lorentzian shadow map acts word by word with `b^{(Δ_34)}` of the seeds (exact
   consistency with conservation); in the polynomial conserved basis `(v_1,v_2)` it reads `(1/b^{(0)}) ×` Weyl reflection. That *some* basis has this
   form is automatic (Hilbert 90, §8) and carries no physics by itself (referee M1); the content is (i) the word-wise mechanism and
   `T(Δ̃)T(Δ) = 𝟙`, (ii) CHL's exact helicity structures are `v_1 + f_h v_2` with Weyl-even `f_h` [ID, exact], and (iii) the physical (MFT)
   normalisation of those structures, item 4. §2.
3. **Factorised formula end to end [DVE].** The tensor Regge moments of the exact waves (81 modes, `k = 4, 6`, ν = 0, 1, 3, β ≈ 0.4, 1, `m` up to
   ≈ 100) match the product of two one-pair vertex factors, `(π/4)4⁴/b^{(0)} γ_a² m^{−2a} Σ_h V V M(ν) M(−ν) 𝒦_{μ_L,μ_R}(β)`, with no fitted constant:
   ratios `0.96 … 1.009` at `m ≈ 100`, approaching 1 as `1 + c_2/m²`. §4–5.
4. **`±` structure, with the physical sign [DVE].** For kernels defined covariantly and identically at both Regge pairs (`D_1`, resp.
   `D_5 = ε(z_in, z_out, x̂_in−out)` with the actual polarisations), the odd⊗odd kernel acts as **`+(|c_same|² − |c_opp|²)`** relative to the even one
   (`K_odd/K_even = +0.961` on the exact same-helicity wave, `−0.961` on the opposite-helicity wave, β ≈ 0.4, 1, ν = 0; prediction and exact
   moments agree to 5 digits). "same" = CHL `T^{++}` = flat all-incoming `++` (`J_z = 0`, all `J`), "opposite" = flat `+−` (`J ≥ 2`): this is the sign
   of the flat-space sum rules with `σ_{++} = +1`, `σ_{+−} = −1` (`FLAT_GRAVITON_BOUND.md`) and of the `SPINNING_E2E.md` §8 conjecture. The
   frame-level formula `R_LR_R(Σ_s+Σ_o) − g_1²S_LS_R(Σ_s−Σ_o)` of §4.3 has `S_R` defined in the inversion-referred frame (`det I = −1`), which is why it
   carries a minus sign (referee M2; §4.5). §4.3–4.5, §5.
   **MFT normalisation [D from exact Fock data, own conventions]:** `p^MFT[ĥ_opp] = 8^J P_FK(1,1;Δ,J)/(8Y²(J−1)J(J+1)(J+2))`,
   `p^MFT[ĥ_same] = 8^J P_FK(1,1;Δ,J)/(8Y²Δ(Δ−1)(Δ−2)(Δ−3))` (two-particle Fock space, exact on all 25 primaries with `n ≥ 1`, `N ≤ 6`; Gram matrix
   exactly diagonal in the helicity basis). Absolute heavy constant: **`1/4` per helicity quadruple**, i.e. `C^±_{k,ν} = (2P_{−½+iν}(cosh β)/m^{2k})
   (|c_same|² ± |c_opp|²)` in the kernel normalisation of §4.4, with **all** factors Weyl-even and corrections `1 + O(1/m²)` (referee M3 resolved
   without CHL's conventions). §4.4.
5. **`O(1/m)` verdict: absent for parity-even exchanges in `⟨JJJJ⟩`** — [D] by the reflection argument in the Weyl-covariant helicity basis
   + [NE] fits (`|c_1| ≤ 0.2` over ansätze and components, raw modes included; referee: up to `0.19` with a `c_3/m³` term) against mutations giving
   `|c_1| = 2.7 … 10` + [D] the exact MFT normalisation is Weyl-even (item 4). Parity-odd exchanges were not built (referee M4). §6.
6. **General `l` [D].** For any spin a Weyl-covariant helicity basis exists (Hilbert 90 + helicity commuting with the shadow map); this is
   bookkeeping and holds for any involutive `T` — every helicity-dependent "shadow eigenvalue" can be renormalised away. The physical content
   is the Weyl parity of the MFT (unitarity) normalisation of the helicity structures: shown for `l = 1` (item 4), open for `l = 2` (§7, §8).
7. **`⟨TTTT⟩` [NE, ingredients only].** `⟨TTO⟩`: 14 words, one 3d relation, two conserved structures; the word-wise shadow map with
   `b^{(0,±2,±4)}` maps conserved to conserved to `10⁻⁹` (three points). Heavy shadow integrand: helicity-pure leading tensor; a weak Richardson
   probe (ratio of right-vertex parity partners) behaves as expected. The `±` algebra and the absence of `O(1/m)` for stress tensors remain
   **conjectural** (helicity-pure Weyl-covariant `l = 2` families and their MFT normalisation not yet shown). End-to-end `⟨TTTT⟩`: **not done**.

## 1. Set-up

* **Correlator.** `⟨J_1 J_2 J_3 J_4⟩`, all four conserved currents (`l = 1`, `Δ_i = 2`), `s`-channel `(12) → (34)`, exchange of a parity-even
  `(Δ, J)`. Regge pairs `(1,3)` and `(2,4)` in the CHMRS frame `(x_1,x_2,x_3,x_4) = (−y, x/x², y, −x/x²)`. No Bose symmetry is imposed
  (the currents may be distinct), so both parity-even conserved `⟨J_1J_2O⟩` structures appear independently at each vertex.
* **Differential basis** [D, CPPR arXiv:1109.6321 form, every property re-checked]: at the vertex `(1,2)` the operators
  `A = D_11`, `B = D_12`, `C = D_22`, `E = D_21` (`D_ij`: spin at `i` raised, `Δ_j` lowered; `D_11, D_22` with their `Z`-derivative terms) and
  multiplication by `H_12`; for spin 1 at both points the five words `AC, AE, BC, BE, H` act on scalar seeds `S3(2+s_1, 2+s_2)` with
  `(s_1,s_2) = (1,1), (2,0), (0,2), (1,1), (1,1)`. They are linearly independent (rank 5) and exactly transverse in `Z_1` and `Z_2`
  (`out_t1_structures_l1.txt`). Partial waves: `W[c_L,c_R] = Σ_{w,w′} c_L[w] c_R[w′] D^{(12)}_w D^{(34)}_{w′}[K g^{(w,w′)}_{Δ,J}]` with the
  scalar blocks of `e2e/cblock.c` (seeds with `a, b ∈ {0, ±1}`); 25 word pairs, ≤ 2100 terms each, jets of order 4.
* **Two-pair moments.** Polarisation circles: pair `(1,3)` at `ŷ`, pair `(2,4)` at `x̂` referred by the inversion tensor
  (`z_{2,4} = I(x̂) z′(θ)`), frames with `e_∥` pointing to the other point and `e_⊥·v = det[v, â, e_∥]`. Modes `𝒢_q`, `q = (q_1,q_3,q_2,q_4)`,
  and the `ν`-harmonic kernel (Fourier–Helgason, generalising `SPINNING_E2E.md` §3.1)
  `𝒦^ν_{μ_L,μ_R}(ξ) = (1/2π)∫dθ (cosh ξ − sinh ξ cos θ)^{−½+iν} e^{iμ_Lθ} e^{iμ_Rφ_R(θ)}`, `e^{iφ_R} = (sinh ξ − cosh ξ cos θ − i sin θ)/(cosh ξ − sinh ξ cos θ)`,
  `Π̂^q_{k,ν} = ∫[dη] 𝒦^ν_{q_1+q_3, q_2+q_4}(ξ) ∫dr r^{k−2} 𝒢_q(r,η)`. Check: `𝒦_{μ,0} = F^ν_μ` of the one-pair study to 12 digits
  (`out_t3_vertex_l1.txt`). The identity (★) of `SPINNING_E2E.md` §3.1 holds verbatim (its proof only uses Lorentz covariance of the
  kernel; `I(x̂)` is Lorentz covariant).

## 2. Three-point structures `⟨J J O⟩` and the shadow matrix (`t1_structures_l1.py`, `t1b_helicity_l1.py`, `chl_map.py`)

**2.1 Conserved space [D, exact].** Conservation at `x_1` and `x_2` (symbolic in `Δ, J`) leaves a two-dimensional space with polynomial basis
(`Y := Δ² − 3Δ + J² + J = A_+² + A_−² − 5/2`, `m² = Δ² − 3Δ − J² − J + 2 = A_+² − A_−²`):

  `v_1 = (m⁴ − 4)·AC − (Δ+J)(Δ−J−1)Y·(AE + BC) + Y²·BE`,   `v_2 = 2·AC + Y·H`.

**2.2 Shadow map [D, exact].** The Lorentzian shadow representation acts word by word (the operators act on points 1, 2 and commute with the
`x_5` integral; for scalar seeds the coefficient is `b^{(Δ_34)}`, [DVE] in `SPINNING_E2E.md` §4.1): the right structure of `W[c_L, c_R]` is
`Σ_w c_R[w]/b_w D_w S3^S`, `b_w = b^{(Δ_34 of the seed)}`, here `b^{(0)}` for `AC, BE, H` and `b^{(±2)}` for `AE, BC`, with
`b^{(2)}/b^{(0)} = (Δ+J)(Δ−J−1)/((Δ+J−2)(Δ−J−3))`. **For both basis vectors, `v_a(Δ,J)·b^{(0)}/b_w = v_a(3−Δ, −1−J)` exactly** (all five
components; `out_t1_structures_l1.txt` (c)). So in the basis `(v_1, v_2)` the shadow matrix is `(1/b^{(0)}) × (Weyl reflection) × 𝟙`.
Referee M1: by §8 such a basis exists for any involutive word-wise map, so "no mixing / no relative eigenvalue" is a statement about a basis
choice, not physics. What is checked here and not automatic: the word-wise mechanism with the `Δ_34`-dependent `b` (an exact identity on all
five components) and `T(Δ̃)T(Δ) = 𝟙`.

**2.3 Helicity basis [D + ID, exact].** CHL's parity-even helicity structures `T^{E,opp}, T^{E,same}` (eq. `VVstcor`, given there in the
frame `(0, x, ∞)` on five monomials `B_V`) were located inside the conserved space by evaluating `v_1, v_2` in the same frame
(embedding point at infinity `P_∞ = (0,1,0)`, `Z_∞ = (0,0,ζ)`), at exact rational `(Δ,J)` (Gaussian-rational null polarisations, exact
linear algebra; `chl_map.py`). Closed forms [ID: fitted on six exact points with a six-parameter polynomial ansatz (an exploratory script, not kept), then checked
exactly at 6 further points, 2 of them at Weyl-reflected labels `(3−Δ, −1−J)`, `out_t1b_helicity_l1.txt`; the referee added `(29/3, 5/2)` and two
more reflected points]:

  `ĥ_opp := v_1 + J(J+1)(m² + 2)·v_2 = v_1 + (A_−² − ¼)(m² + 2)·v_2`,   `ĥ_same := v_1 − (Δ−1)(Δ−2)(m² − 2)·v_2 = v_1 − (A_+² − ¼)(m² − 2)·v_2`,

  `T^{E,opp} ∝ ĥ_opp/(J(J−1)αY)`,  `T^{E,same} ∝ ĥ_same/(Δ(Δ−1)αY)` (common factor), `α = Δ − J − 2`.

Both coefficient ratios are **Weyl-even**: the exact helicity structures are Weyl-covariant families, `T ĥ_h(Δ,J) = ĥ_h(3−Δ,−1−J)`.
This is forced [D]: the helicity operators `h_1, h_2` (conformally invariant in `d = 3`, CHL §2.2 [Q]) act on points 1, 2 and commute with
any transform at the exchanged point, so the shadow map sends helicity lines to helicity lines; with the map equal to the identity on
`(v_1, v_2)`, a line `v_1 + f v_2` goes to `v_1 + f(Δ,J) v_2` at the reflected labels, which is the same line iff `f(Δ,J) = f(3−Δ,−1−J)`.
CHL's helicity-dependent shadow factors `(2−Δ)_{|h_1+h_2|}/(Δ−1)_{|h_1+h_2|}` (eq. `shadow T` [Q]; in the E/O basis eq. `hard shadow`) are
therefore of the form `f(Δ̃)/f(Δ)` — as any such factor must be (§8); whether they matter is decided by the physical normalisation, §4.4. Consistency check [I]: with
`T^{E,same} ∝ ĥ_same/(Δ(Δ−1)αY)` and `T^{E,opp} ∝ ĥ_opp/(J(J−1)αY)`, the *Euclidean* shadow (`Δ → 3−Δ`, `J` fixed) acting covariantly on `ĥ_h` produces
relative eigenvalues `(3−Δ)(2−Δ)/(Δ(Δ−1)) = (2−Δ)_2/(Δ−1)_2` (same) and `1` (opposite) — exactly CHL's `hard shadow` matrix.

## 3. Exact partial waves: validation (`t2_validate_l1.py`, `out_t2_validate_l1.txt`)

| check | result |
|---|---|
| conservation at `x_1`, `x_2`, `x_3` of `W[v_a,v_b]` (`a,b = 1,2`), `(Δ,J) = (9.75,3), (14.5,4), (21.25,6)`: `|∂^μW_μ|/max|∂_μW_ν|` (4th-order differences, double) | `10⁻¹³ … 2·10⁻⁸` (difference-limited) |
| control: word `AC` alone at the left vertex (not conserved) | `0.7 … 1.6` at `x_1, x_2`; `10⁻⁸` at `x_3` (right vertex conserved) |
| pair exchange `W[v_1,v_2](1234) = W[v_2,v_1](3412)` | `4·10⁻¹⁵ … 2·10⁻¹⁴` |
| reduction to scalar exchange: word pair `(H,H)` vs `H_12H_34 · K g_{(3,3,3,3)}` (independent scalar evaluation) | `4·10⁻¹⁶ … 4·10⁻¹⁵` |
| exact mode symmetries (two-pair Lemma S): `Π̂^q[v_1,v_2] = Π̂^{q′}[v_2,v_1]`, `q′ = (q_3,q_1,q_4,q_2)`, and `= Π̂^{−q}` | `≤ 7·10⁻¹⁰` and `≤ 2.4·10⁻¹¹` relative to the dominant mode, all 10 cases of §5 (`k = 6`, ν = 1) |

Block-level normalisation: the leading-twist (`n = 0`) decomposition of `⟨JJJJ⟩_MFT` with the exact Fock coefficients is verified (§4.6). Not
done: the full MFT decomposition (needs the parity-odd waves as well).

## 4. The factorised leading formula (`vertex2.py`, `t3_vertex_l1.py`, `predict4.py`)

**4.1 Heavy vertex [ID + DV pointwise].** In the frame, `x = X/m`, `y = Y/m`, the exact shadow integrand
`𝒯^{ab} = s_13²s_24²((−1)^J/(4b^{(0)})) L[v_a](Δ,J) · L[v_b](3−Δ,−1−J)` (right factor by §2.2) divided by
`((−1)^J/(4b^{(0)}))·2·4⁴|x|⁴|y|⁴e^{2P·x+2P̄·y}`, `P = (A_+, A_−, 0)`, `P̄ = (A_+,−A_−,0)`, tends to `Σ_h V_a(h_1,h_2)V_b(h_3,h_4) Π(z·ε_h)`
with ε's in the frames at `p̄̂` (pair 13) and `p̂` (pair 24), and

  `V_1(h_1,h_2) = (A_+A_−)² m²/2` (both `h_2 = −h_1`, "same", and `h_2 = h_1`, "opposite"),   `V_2 = −A_+²/2` (same), `+A_−²/2` (opposite).

(Identified from the exact integrand at `m` up to `5·10⁵`; deviation `∝ 1/m` for all 16 helicity components and all `(a,b)`, three `β`;
longitudinal leakage `∝ 1/m`; the Weyl-symmetrised deviation `∝ 1/m²`: `m²·max|sym − 1| = 16.2, 63.6, 14.0` converged at `β = 0.84, 0.40, 1.59`;
`out_t3_vertex_l1.txt`.) Consequently `ĥ_opp → A_−²m²(A_+²+A_−²)/2 · (opposite only)`, `ĥ_same → A_+²m²(A_+²+A_−²)/2 · (same only)` up to
`O(m⁻²)` relative leakage: the CHL labels are confirmed in the heavy limit ("same" = charge 0 under the common transverse rotation, exists
for all `J`; "opposite" = charge ±2, `J ≥ 2`).

**4.2 Prediction [D].** With (★), the two Fourier–Laplace transforms at fixed boundary point `ℓ` (`SPINNING_E2E.md` §4.4; the `x`-pair with `ν → −ν`),
and the remaining `ℓ`-integral = `𝒦^ν_{μ_L,μ_R}(β)`:

**`(−1)^J Π̂^q_{k,ν}[s_13²s_24² W[c_L,c_R]] → (π/4)·4⁴/b^{(0)}_{Δ,J} · γ_a(ν)² m^{−2a} · Σ_h V_{c_L}(h_1,h_2) V_{c_R}(h_3,h_4) M^{q_1q_3}_{h_1h_3}(ν;a) M^{q_2q_4}_{h_2h_4}(−ν;a) 𝒦^ν_{μ_L(h),μ_R(h)}(β)`**,

`a = k + 3`, `V_c = Σ c_a V_a`, `M = R + g_1S` the one-pair table (`μ = 0` for `h_3 = −h_1`, ±2 otherwise). Nothing is fitted.

**4.3 The `±` structure [D].** In the component preserved at both pairs (`h_3 = −h_1`, `h_4 = −h_2`; harmonic `𝒦_{0,0} = P_{−½+iν}(cosh β)`),
`M_{h,−h} = R + h g_1 S` and parity `V(h_1,h_2) = V(−h_1,−h_2)` give

  `Π̂|_{pres,pres} ∝ P_{−½+iν}(cosh β) [ R_LR_R (Σ_s + Σ_o) − g_1² S_LS_R (Σ_s − Σ_o) ]`,  `Σ_s = 2V(+,−)²` (same), `Σ_o = 2V(+,+)²` (opposite),

with `S_L`, `S_R` the odd parts of the `M` tables in the frames of §1 (pair 24 referred by the inversion tensor). Which sign this is
physically is derived in §4.5 (referee M2): for covariantly defined kernels the opposite-sign term enters as `+(|c_same|² − |c_opp|²)`, the
conjectured sign. The opposite-sign term is a parity-even effect of the doubly-spinning block (it needed a parity-odd structure in `⟨OφOφ⟩`).
For `W[v_1,v_1]` (`V_1` same = opposite) the opposite-sign term vanishes identically; for `W[v_2,v_2]` it is `∝ A_−⁴ − A_+⁴ ≠ 0`.

**4.5 The physical sign of the opposite-sign term [D + DVE] (`tttt/t4_sign_l1.py`, `tttt/out_t4_sign_l1.txt`).** Define at each Regge pair
`(in, out) = (1,3), (2,4)` — by the same rule, with the operators' actual polarisations and `n̂ ∥ x_in − x_out` — the even kernel
`D_1 = z_in·z_out + (z_in·n̂)(z_out·n̂)` and the odd kernel `D_5 = ε(z_in, z_out, n̂)`. In the frame, `x_13 = −2y` (`n̂ = −ŷ`) and `x_24 = −2x/|x|²` (`n̂ = −x̂`),
and `z_{2,4} = I(x̂)z′` with `det I = −1`, `I x̂ = −x̂`; the two orientation flips cancel, so `D_5^{13}D_5^{24} = ε(z_1,z_3,ŷ)ε(z′_2,z′_4,x̂)`, i.e.
`sin(θ_3−θ_1)sin(θ_4−θ_2)` on the circles (and `D_1 = cos(θ_out−θ_in)`). On the leading tensor `ε(ε_h, ε_{−h}, â) = −ih` in the det-rule frames, so pointwise
`D_5⊗D_5` weighs `−h_1h_2` (frame labels) `= +1` on "same" (`h_2 = −h_1` in the inversion-referred frame) and `−1` on "opposite". After the
Laplace transforms (end to end, exact moments, `k = 6`):

| β, `m` | ν | `K_odd/K_even` on `W[ĥ_same,ĥ_same]`: pred / exact | on `W[ĥ_opp,ĥ_opp]`: pred / exact |
|---|---|---|---|
| 0.4, 101.8 | 0 / 1 / 3 | +0.96122 / +0.96122;  +0.96731 / +0.96733;  +0.99755 / +0.99766 | −0.96127 / −0.96127;  −0.96765 / −0.96765;  −1.00475 / −1.00477 |
| 1.0, 100.7 | 0 / 1 / 3 | +0.96125 / +0.96126;  +0.96759 / +0.96760;  +0.99115 / +0.99123 | −0.96127 / −0.96127;  −0.96765 / −0.96764;  −1.01049 / −1.01059 |

(the same at `m ≈ 70`). So the covariant odd⊗odd kernel acts as `+(|c_same|² − |c_opp|²)` and the even one as `(|c_same|² + |c_opp|²)`, with
`|K_odd| ≈ K_even` — `K_even ± K_odd ≥ 0` per helicity, as unitarity requires. **Flat-space dictionary [D, modulo the identification of
functionals].** The labels are physical: "same" (CHL `T^{++}`, `J_z = h_1 − h_2 = 0`, all `J`) is the flat all-incoming `++` state of
arXiv:2201.06602 (`D_{++} = P_J`), "opposite" (`J ≥ 2l`) is flat `+−` (`D_{+−} = d̃^J_{2l,2l}`). In flat space `ε(ε_in, ε̄_out, n̂)` with helicity preserved
is proportional to the incoming helicity, with the same proportionality at both pairs (the rule is symmetric under `(1,3) ↔ (2,4)`), so
`D_5⊗D_5` weighs `h_1h_2`: `+1` for `++`, `−1` for `+−` — the `σ_{++} = +1`, `σ_{+−} = −1` of the flat graviton functional (`FLAT_GRAVITON_BOUND.md` (0.1)).
Thus the AdS opposite-sign term has **the same sign as the flat-space one and as the `SPINNING_E2E.md` §8 conjecture** read with physical
labels; the minus sign of §4.3 is the frame convention for `S_R`. What is *not* derived: which AdS functional uplifts which flat sum rule
(e.g. that `B^{(1)}_3` corresponds to a `D_5⊗D_5`-type kernel with positive coefficient) — `REGGE_MOMENT_TTTT.md` §5 [I].

**4.4 Mean-field normalisation — superseded by §4.6 for the absolute part [D from Q].** CHL's MFT data are diagonal in `T^{E,h}` with `H_opp = (−J)_2/(J+1)_2`, `H_same = (Δ−1)_2/(2−Δ)_2`
times the scalar `(Δ_1,Δ_2) = (1,1)` data [Q, eq. `free data`]. With §2.3 and §4.1, `p^MFT[ĥ_h]·V_h² = (h-independent factor) × Q_h`,

  `Q_opp = 16A_−⁴m⁴(A_+²+A_−²)²/((4A_−²−1)(4A_−²−9)(2Y)²)`,  `Q_same = 16A_+⁴m⁴(A_+²+A_−²)²/((4A_+²−1)(4A_+²−9)(2Y)²)` (`out_t8_mft_l1.txt`),

both **Weyl-even** and both `→ m⁴/4`: `Q_opp/Q_same = 1 + 10/(m² sinh²β) + O(m⁻⁴)`. The scalar factor is `P_FK(1,1)/b^{(0)} = (8/π²)A_+A_−`, also
Weyl-even (`SPINNING_MFT.md` §4.1 with `S = 2`). Hence, MFT-normalised (`|c^h|² = p^h/p^h_MFT`), the two helicities enter the preserved
component with **equal weight** — `C^± ∝ (|c^{++}|² ± |c^{+−}|²)` in the sense of §4.3 — and with **no `O(1/m)` term from the normalisation**.
(The `n = α/2` prefactor of CHL's `VVstcor` cancels the non-Weyl-even `1/α` of `a_h`; `α²` drops out — referee m7.) Referee M3: the
`h`-independent factor (`K_c² 2^{J+1}` and CHL's scalar data vs `P_FK(1,1)/b^{(0)}`) was asserted, not mapped, and a `J`-dependent mismatch would
not be harmless. **This is now replaced by an exact computation in our own conventions, §4.6.**

**4.6 Exact MFT data and the absolute two-pair constant [D, exact Fock data] (`tttt/fock2.py`, `tttt/t1_mft_l1.py`, `tttt/t2_mft_heavy_l1.py`).**
The method of `SPINNING_MFT.md` §2 extended to two spinning particles: a primary `[J_1J_2]_{N,J}` has wavefunction `⟨X(∞,ζ)J(x_1,z_1)J(x_2,z_2)⟩`,
a polynomial conserved at both points; the norm is the tensor product of the one-particle conserved-current Gram forms (`mft/fock.py`, unchanged).
The structures `ĥ_h` are evaluated at infinity directly from the differential basis (`P_∞ = (0,1,0)`, `Z_∞ = (0,0,ζ)`; the conserved space found there
reproduces `v_1, v_2` exactly, `tttt/infcons.py`). Results (`out_t1_mft_l1.txt`, all `N ≤ 6`, exact rationals):
* `n = (τ−2)/2 ≥ 1`: the 2×2 Gram matrix of `(ĥ_opp, ĥ_same)` is **exactly diagonal** (MFT is helicity-diagonal, confirming CHL and the identification
  of §2.3); `ĥ_opp` is absent for `J = 0, 1`;
* `n = 0`: `ĥ_opp = ĥ_same` (one primary; helicity undefined, as CHL state);
* closed forms, exact on all 25 primaries with `n ≥ 1` [NE, exact; no parameter]:
  **`p^MFT[ĥ_opp] = 8^J P_FK(1,1;Δ,J)/(8Y²(J−1)J(J+1)(J+2))`**, **`p^MFT[ĥ_same] = 8^J P_FK(1,1;Δ,J)/(8Y²Δ(Δ−1)(Δ−2)(Δ−3))`**;
  at `n = 0` the single primary carries the sum of the two (exact on the 7 values).
`8^J` is the convention factor of the scalar seed (`p[S3] = 8^J P_FK`), which the block normalisation `(−1/8)^J` cancels. With `P_FK(1,1)/b^{(0)} =
(8/π²)A_+A_−` and the verified heavy formula §4.2, the MFT-normalised moment of one helicity-`h` operator, in the normalisation of `SPINNING_MFT.md`
(`a_{Δ_O,Δ_O}/γ_a²`, density `q = 2n_J(2Δ−3)/(πm²p^MFT)`), is

  **`C^{(q)}_{k,ν}[Δ,J; h] = (1/(4m^{2k})) Σ_{h∈class h} M^{q_1q_3}_{h_1h_3}(ν) M^{q_2q_4}_{h_2h_4}(−ν) 𝒦^ν_{μ_L,μ_R}(β) · Q_h`**,  `Q_h = 1 + q_2/m² + …`

(`out_t2_mft_heavy_l1.txt`: `m²(Q_opp−1) → 66.30, 12.45, 3.139` and `m²(Q_same−1) → 7.028, 5.206, 2.379` at β = 0.4, 1, 2, converged from `m = 100`
to `10⁴`; `m(Q−1) → 0`). So the constant is `1/4` per helicity quadruple (the same `1/4` as per helicity pair in the one-pair formula of
`SPINNING_MFT.md`), every factor is Weyl-even, and the normalisation adds no `O(1/m)` term. In the preserved–preserved component a helicity-`h`
operator contributes `½R_LR_R P ∓ ½g_1²S_LS_R P`; normalising the kernels by `R_LR_R/4` resp. `g_1²S_LS_R/4` gives
**`C^± = (2P_{−½+iν}(η_AdS)/m^{2k})(|c_same|² ± |c_opp|²)`** (sign of the odd kernel as in §4.5), the form conjectured in `REGGE_MOMENT_TTTT.md` §5, now
derived for currents with its constant. Remaining assumption: the doubly-spinning version of the block normalisation "contribution
`λ_Lλ_R(−1/8)^J W[c_L,c_R]`" (`SPINNING_E2E.md` §1.4; tested at block level for one pair in `SPINNING_MFT.md` §5.1); its two-pair block test is
**done** (`tttt/t3_lightcone_l1.py`, `tttt/out_t3_lightcone_l1.txt`) [DV]: on the real section with `x_12 = T(1,(1−ε)n̂)` nearly
null and `z_1 = z_2 = n_12`, the MFT correlator `⟨J_1J_3⟩⟨J_2J_4⟩` equals `Σ_{J=2}^{30} p(n=0,J)(−1/8)^J W[ĥ,ĥ]^{(J+2,J)}` (exact blocks of `fourpt/spinblock2`,
exact Fock `p`) with ratio `1 + O(ε)` and no free constant or sign: `ratio − 1 = 1.6·10⁻⁶, −4.7·10⁻⁶, −1.6·10⁻⁴` (point A) and `1.7·10⁻⁶, 2.3·10⁻⁶,
−1.6·10⁻⁶` (point B) at `ε = 10⁻⁵` for three `(z_3, z_4)`, linear in `ε`; control `(+1/8)^J`: `12.6`, `9.4`. This tests, for the whole leading-twist
family, the chain Fock `p` → structure `ĥ` → contribution `p(−1/8)^J W` for two spinning pairs. Families with `n ≥ 1` (where the helicities
separate) enter only through the exact closed forms.

## 5. End-to-end comparison (`t4_compute_l1.py`, `t5_analysis_l1.py`, `out_t5_analysis_l1.txt`)

**Cases.** `(Δ,J)` at fixed `β ≈ 0.4` (`m = 12.4, 22.4, 42.2, 72.0, 101.8`; `J = 2…20`) and `β ≈ 1` (`m = 10.6, 20.1, 39.3, 70.0, 100.7`; `J = 5…52`),
as in `SPINNING_E2E.md` §5; `k = 4, 6`; `ν = 0, 1, 3`; all 81 modes; the four basis waves `W[v_a,v_b]` and, by bilinearity, the exact-helicity
waves `W[ĥ_opp,ĥ_opp]` ("oo"), `W[ĥ_same,ĥ_same]` ("ss"), `W[ĥ_opp,ĥ_same]` ("os"). Radial truncation `N = max(4000, 400Δ)` capped at `10⁴`;
cost ≈ 30 min per case on the shared machine.

**Extraction.** The 8 parity classes of helicity components `T_h` are fitted to the 81 measured modes with the *predicted* weights
`M M 𝒦` (weighted least squares); `R_c = T_c(fit)/T_c(pred)`. Classes: `oo` = `(+,+,−,−)` and `ss` = `(+,−,−,+)` (preserved at both pairs; opposite /
same helicity at the vertex), `ff_oo`, `ff_ss` (flipped at both), `pf`, `fp` (one flipped); `Σ_± = T_oo ± T_ss` (§4.3: `Σ_+` multiplies `R_LR_R`, `Σ_−` the
opposite-sign `g_1²S_LS_R`). The quoted LSQ residuals are dominated by the `1/|pred|` weighting and are not a goodness-of-fit measure (referee m3): with uniform
weights the 8-class residual is `0.9…1.6·10⁻³`, and a **one-scale fit** (the predicted pattern of all 81 modes times one free number) has residual
`1.9·10⁻³` with scale `0.9932` (`W[v_2,v_2]`, β ≈ 0.4, `m = 101.8`, ν = 0; referee's `ref_analysis.py`). Raw modes (no class fit) give the same ratios.

**Numerical accuracy.** Tail estimates (dropped small-`u` region) at the largest `m`: `k = 6`: `10⁻⁶ … 1.4·10⁻⁵`; `k = 4`: `2·10⁻⁵ … 3·10⁻³` (β ≈ 1).
Exact symmetries hold to `10⁻¹⁰` (§3).

**Results at the largest mass (`k = 6`; complete tables for `k = 4, 6`, all `m`, in `out_t5_analysis_l1.txt`):**

| β, `m` | ν | `W[v_1,v_1]`: oo / ss | `W[v_2,v_2]`: ss / Σ₊ / Σ₋ | `W[v_1,v_2]`: ss / Σ₊ / Σ₋ | `W[ĥ_opp,ĥ_opp]`: oo | `W[ĥ_same,ĥ_same]`: ss | `W[ĥ_opp,ĥ_same]`: pf | LSQ resid. |
|---|---|---|---|---|---|---|---|---|
| 0.4, 101.8 | 0 | 0.9814 / 0.9942 | 0.9933 / 0.9932 / 0.9933 | 0.9927 / 0.9931 / 0.9924 | 0.9608 | 0.9932 | 0.9846 | 0.009–0.07 |
| 0.4, 101.8 | 1 | 0.9822 / 0.9951 | 0.9939 / 0.9939 / 0.9939 | 0.9933 / 0.9937 / 0.9929 | 0.9801 | 0.9939 | 0.9805 | 0.005–0.015 |
| 0.4, 101.8 | 3 | 0.9917 / 0.9996 | 0.9982 / 0.9982 / 0.9983 | 0.9965 / 0.9966 / 0.9964 | 0.9939 | 0.9981 | 0.9841 | 0.002–0.006 |
| 1.0, 100.7 | 0 | 0.9705 / 1.0034 | 0.9929 / 0.9923 / 0.9935 | 0.9960 / 0.9982 / 0.9946 | 0.9616 | 0.9922 | 0.9843 | 0.02–0.07 |
| 1.0, 100.7 | 1 | 0.9702 / 1.0094 | 0.9938 / 0.9934 / 0.9943 | 0.9829 / 0.9831 / 0.9828 | 0.9834 | 0.9935 | 0.9821 | 0.006–0.05 |
| 1.0, 100.7 | 3 | 0.9752 / 0.9828 | 0.9825 / 0.9819 / 0.9832 | 0.9832 / 0.9841 / 0.9827 | 0.9746 | 0.9849 | 0.9922 | 0.002–0.008 |

(For `W[v_1,v_1]` the predicted `Σ_−` vanishes identically, §4.3; the measured `Σ_−` is `O(1/m²)` of `Σ_+` there. In `W[v_2,v_2]`,
`T_oo/T_ss ≈ (A_−/A_+)⁴ ≈ 1.6·10⁻³` at β ≈ 0.4, so `R(Σ_−) ≡ R(Σ_+) ≡ R_ss` there and is not separate evidence (referee m4); the independent
information on the `S` term is `R_oo`, the exact-helicity waves and §4.5.) Approach to 1, e.g. `W[ĥ_same,ĥ_same]`
`ss`, β ≈ 0.4, ν = 0: `0.704, 0.880, 0.9625, 0.9867, 0.9932` at `m = 12…102`, `m²(R−1) = −46, −60, −67, −69, −70`; `W[ĥ_opp,ĥ_opp]` `oo`, same
β, ν: `−0.12, 0.39, 0.79, 0.923, 0.961`, `m²(R−1) = −173, −306, −378, −400, −406`. Every leading component approaches 1 as `1 + c_2/m²`, with `c_2`
between `−20` and `−420` (fits B, `out_t5_analysis_l1.txt`); the opposite-helicity component converges more slowly at small β (its leading
tensor is suppressed by `tanh²(β/2)` relative to the same-helicity one, so the effective expansion parameter is `1/(m sinh(β/2))² ≈ 1/J²`, as for the
flip component of the one-pair study). Subdominant classes (e.g. `ff_oo` in `W[v_2,v_2]`, predicted `∝ A_−⁴` against `A_+⁴`) are dominated by the
`1/m²` corrections of the dominant ones and are not meaningful as ratios.

**The `±` structure, end to end [DVE at `k = 4, 6`, β ≈ 0.4, 1].** In `W[v_2,v_2]` and `W[v_1,v_2]` the predicted opposite-sign combination `Σ_−` is
non-zero and the measured one agrees with it at the same level as `Σ_+` (table); in the exact helicity basis the same-same and opposite-opposite
components separately approach 1, which is the statement `Π̂|_{pres,pres} ∝ R_LR_R(|c^{++}|²+|c^{+−}|²) − g_1²S_LS_R(|c^{++}|²−|c^{+−}|²)` with the
derived `M`'s (referee's mutations: `S → 0.5S, 2S` gives `R(Σ_−) = 3.98, 0.248`; `S = 0` raises the residual to 0.59). The oo and ss columns of the fit differ only through `S` (`M_{+−} − M_{−+} = 2g_1S`), so the separate agreement of `R_oo` and
`R_ss` (equivalently of `Σ_−`) with the prediction is a direct test of the `g_1²S_LS_R` term — the first end-to-end test of the `S` factors, which
`SPINNING_E2E.md` could only check by `H²` quadrature.

## 6. `O(1/m)`

**6.1 Argument [D, given Assumption T of `SPINNING_E2E.md` §6.1].** By §2.2, `b^{(0)}𝒯^{ab}(z_1,z_2,z_3,z_4; x,y; A_±) = L[v_a](A_±)(z_1,z_2; −y, x/x²) ·
L[v_b](−A_±)(z_3,z_4; y, −x/x²)`, so the reflection `(x,y,A_±) → (−x,−y,−A_±)` maps `b𝒯^{ab}(z_1,z_2,z_3,z_4)` to `b𝒯^{ba}(z_3,z_4,z_1,z_2)` (the inversion
tensor `I(x̂)` is invariant). At fixed `(X,Y)` this is `m → −m`, so the `O(1/m)` term `t^{ab}_1` obeys `t^{ab}_1(z_1,z_2,z_3,z_4) = −t^{ba}_1(z_3,z_4,z_1,z_2)`.
The exact two-pair mode symmetry (pair exchange `(12) ↔ (34)` of the `s`-channel wave combined with `v → −v` of the frame, as in Lemma S;
verified to `10⁻¹⁵` pointwise in §3 and on the moments) gives `Π̂^q[ab] = Π̂^{q′}[ba]`. Together: no `O(1/m)` term in any mode of any
`W[v_a,v_b]`, in the Weyl-odd variables with the exact `1/b^{(0)}`. Any basis with Weyl-even rational coefficients (in particular the exact
helicity basis, §2.3) inherits this, since Weyl-even rational functions of `A_± = m(cosh, sinh)(β/2)` expand in `1/m²`.
Pointwise check: §4.1 (symmetrised deviation `∝ 1/m²`).

**6.2 Fits and mutations [NE].** `R − 1 = c_1/m + c_2/m² + …` (`m² = αᾱ`); ansatz A: `c_1, c_2, c_4` from the three largest `m`; C: `c_1…c_4` from
the four largest; `k = 6`. Leading components (full list in `out_t5_analysis_l1.txt`):

| β | ν | `W[v_1,v_1]` ss / oo | `W[v_2,v_2]` ss / Σ₋ | `W[v_1,v_2]` ss / Σ₋ | `W[ĥ_opp,ĥ_opp]` oo | `W[ĥ_same,ĥ_same]` ss | `W[ĥ_opp,ĥ_same]` pf |
|---|---|---|---|---|---|---|---|
| 0.4 | 0 | −0.024 / −0.029 | −0.033 / −0.033 | −0.044 / −0.043 | −0.040 | −0.034 | **−2.6** |
| 0.4 | 1 | −0.005 / −0.008 | −0.007 / −0.007 | −0.010 / −0.011 | −0.012 | −0.007 | +0.015 |
| 0.4 | 3 | +0.005 / +0.005 | +0.002 / +0.002 | +0.005 / +0.005 | +0.002 | +0.002 | −0.002 |
| 1.0 | 0 | −0.16 / −0.19 | −0.099 / −0.099 | −0.21 / −0.18 | −0.16 | −0.15 | +0.023 |
| 1.0 | 1 | −0.031 / −0.064 | −0.014 / −0.015 | −0.094 / −0.083 | −0.038 | −0.030 | −0.016 |
| 1.0 | 3 | −0.030 / −0.031 | −0.030 / −0.022 | −0.008 / −0.004 | −0.022 | −0.031 | −0.000 |

(`c_1` of ansatz A; ansatz C differs by up to ±0.1 at β ≈ 1, ν = 0, and with `(c_1, c_2, c_3)` on four masses the raw modes give `c_1 = +0.05 … +0.19`
(referee) — the bound is ansatz-dependent, `|c_1| ≤ 0.2` over all ansätze tried.) **Mutations** (same fits, the prediction
altered at `O(1/m)`): `b^{(0)} → b^{(2)}` (the shadow coefficient of the words `AE, BC`): `c_1 = +4.06 … +4.60` in every component; normalisation
`A_+² → Δ²` in the same-helicity part of `V_2` (not Weyl-even): `c_1 = −5.4 … −5.9` in `W[v_2,v_2]` ss, `−2.7 … −3.0` in `W[v_1,v_2]` ss, and
unchanged where `V_2`'s same part does not enter (a built-in control). The measured `|c_1|` are at least 20 (typically 100–1000) times smaller than the mutation values,
of the size of the neglected `c_3/m²`, `c_6/m⁵` terms, and decrease with ν as `c_2` does.

Weak spots: (i) `W[ĥ_opp,ĥ_same]` (one pair preserved, one flipped) at β ≈ 0.4, ν = 0: the class fit gives `c_1 = −2.6`, but this is an
artefact of the 8-class fit — the raw modes of this wave (which at leading order contain only the pf/fp classes) give `c_1(A) = −0.009` at ν = 0 and
`−0.006` at ν = 1 (referee; the earlier "leakage" explanation was wrong and is withdrawn);
(ii) β ≈ 1, ν = 0: `|c_1| ≤ 0.21` with `c_2` of either sign and small (`R_ss(v_1,v_1)` crosses 1: `0.96, 1.04, 1.020, 1.007, 1.003`), so the 3-mass fit
cannot separate `c_1` from higher terms — still a factor ≥ 20 below the mutations (referee's further mutations: a CHL-type non-Weyl-even
factor `Δ(Δ−1)/((3−Δ)(2−Δ))` on one same-helicity vertex gives `c_1 = −3.97`; `A_−² → J²` in `V_2(opp)` gives `+10.2`); (iii) only `m ≤ 102`, `k = 6` (at `k = 4` the tails reach `3·10⁻³`).

**6.3 Verdict.** For `⟨JJJJ⟩` there is **no `O(1/m)` term**: [D] for every `W[v_a,v_b]` and hence for the exact helicity waves (§6.1, given
Assumption T); [NE] confirmed by the fits at `k = 6`, β ≈ 0.4, 1 against mutations that are detected at `O(1)` in `c_1` (with the weak spots
above); and [D from Q] it survives the MFT normalisation of the physical helicity couplings (§4.4: all factors Weyl-even). The referee's mechanism
(helicity-dependent shadow eigenvalues `κ_h` with non-Weyl-even ratios, SPINNING_E2E §11 item 10) does not occur: the eigenvalues are
normalisation artefacts `f_h(Δ̃)/f_h(Δ)`, and the physical normalisation compensates them (§8). With the exact MFT normalisation of §4.6 (all factors Weyl-even) the
conjecture of `SPINNING_E2E.md` §6.1 is **settled for parity-even exchanges in `⟨JJJJ⟩`** (derivation + end-to-end numerics + exact MFT data) and
**reduced for stress tensors** to the Weyl parity of the `l = 2` MFT normalisation in the helicity basis (§7, §8). **Parity-odd exchanges
(referee M4)** — operators at odd `τ − 2`, with the two parity-odd `⟨JJO⟩` structures of CHL (`T^{O,h}`) — were not built (they need the
`ε`-operator `D̃` at both points); the end-to-end match, the `±` structure and the `O(1/m)` verdict are untested there. By §8 a Weyl-covariant odd
helicity basis exists; the open point is again the Weyl parity of its MFT normalisation (CHL: odd data `½c(1,2)^s` × the same helicity matrix [Q]).

## 7. `⟨TTTT⟩` (`build_l2.py`, `t6_structures_l2.py`, `t7_vertex_l2.py`)

**7.1 `⟨TTO⟩` structures [NE] (`out_t6_structures_l2.txt`).** 14 words `{AA,AB,BB}×{CC,CE,EE}`, `{A,B}×{C,E}×H`, `HH` (built exactly by `build_l2.py`,
≤ 5300 terms). On physical 3d kinematics they satisfy **one** linear relation (singular value `6·10⁻⁵⁴` at 60 digits; it involves
`AACC, ABCE, BBEE, ACH, AEH, BCH, BEH, HH`). Conservation at `x_1` and `x_2` (spin-2 divergence of the traceless part, 6th-order differences,
double precision) leaves a 3-dim null space = **2 conserved structures + the relation** (smallest singular values `10⁻¹²…2·10⁻¹¹` against
`≥ 10⁻³`), at `(Δ,J) = (9.7,6), (13.25,5), (11.5,8)` and at the reflected labels; checked at a second, independent configuration. This agrees
with the two parity-even conserved `⟨TTO⟩` structures of arXiv:1708.05718 / `REGGE_MOMENT_TTTT.md` §3.2.

**7.2 Shadow map [NE, 10⁻⁹].** The word-wise map `c_w ↦ c_w b^{(0)}/b_w` (seeds with `Δ_34 = 0, ±2, ±4`) sends the conserved space at `(Δ,J)` into the
conserved space at `(3−Δ,−1−J)`: residual `1·10⁻¹⁰ … 1.7·10⁻⁹` for all three null vectors at the three points, against `0.16 … 0.95` for the control with
`b^{(0)}` for every word; the relation itself maps to the relation with a common factor (all components to 10 digits). (A first run reported a
failure; the cause was a bug in the column scaling of the null vectors, found by checking them at an independent configuration.) So the
hypotheses of §8 hold for stress tensors, and Weyl-covariant families exist.

**7.3 Heavy shadow integrand with Weyl-covariant families [NE, partial] (`out_t7_vertex_l2.txt`).** Two families `v_1, v_2` built numerically by
the Hilbert-90 recipe `v = c + T(Δ̃)c(Δ̃)` (two different rational normalisation rules for `c`), `T v(Δ) = v(Δ̃)` exactly by construction. At exactly
fixed `β = 0.843` (`A_± = s·(48.75, 20.5)`, `s = 1, 2, 4, 8`, `m = 44 … 354`; `J` non-integer is allowed at the level of three-point functions):
* the leading tensor has pure helicity content `(z_1·ε_{±}(p̄̂))²(z_2′·ε_{±}(p̂))²`: non-helicity modes `/` helicity modes `= 0.046, 0.0056, 0.0029` at `m = 44,
  88, 177` (0.007 at `m = 354`, where the double-precision null space degrades: gap `sv_3/sv_4 = 6.7·10⁻³`);
* parity `T_{−h} = \bar T_h` exactly; factorisation `T^{ab}_h = V_a(h_1,h_2)V_b(h_3,h_4)` is exact pointwise (the integrand is a product of two three-point
  functions) and not an independent test;
* `1/m` structure, without closed forms: Richardson ratios `(g(2s)−g(s))/(g(4s)−g(2s))` (→ 2 for a `1/m` correction, → 4 for `1/m²`) of the dominant
  helicity component relative to a reference — the ratio `T(+−+−)/T(+−−+)` of right-vertex parity partners whose leading values are equal by
  construction, hence a **weak probe** (referee m8): unsymmetrised `2.13, 1.99` (the pointwise `1/m` term is there), **Weyl-symmetrised `4.13, 3.94`** — the
  `O(1/m)` term is antisymmetric under the exchange of the two vertices, exactly as in §6.1. For the subdominant components (`≲ 0.1` of the
  dominant one) the ratios are erratic (e.g. `13.0, 0.89`; `2.2, 10.5`), consistent with the precision loss of the double-precision `⟨TTO⟩` null
  vectors (gap `10⁻⁴` at `m = 177`) — **inconclusive** there.

**7.4 What this does and does not establish for `⟨TTTT⟩`.** Checked at three `(Δ,J)`: conserved space, word-wise shadow map. Pointwise: a
helicity-pure leading tensor for numerically built Weyl-covariant families. **Conjectural** (referee (h)): that these families can be chosen
helicity-pure (exact `l = 2` helicity structures), hence the `±` algebra for `l = 2`, and the absence of `O(1/m)`. Not established: the end-to-end
test on exact `⟨TTTT⟩` blocks (order-8 jets, 196 word pairs — not attempted); the closed forms of the `l = 2` heavy vertex in the helicity basis
and the exact `l = 2` helicity structures; the `l = 2` MFT normalisation in that basis, hence the absence of `O(1/m)` for physical `⟨TTTT⟩`
couplings, which remains **conditional** on the Weyl-evenness of `p^{MFT}[ĥ_h]·V_h²` for `l = 2` (true for `l = 1`, §4.4; CHL's `l = 2` MFT data have the same
structure `H_h = (−J)_{J_1+J_2}(Δ−1)_{|J_1−J_2|}/…` [Q], which makes it likely [I]).

## 8. The shadow matrix of conserved `⟨O O O′⟩` structures in general (stage 4)

Setting: `V(Δ,J) ⊂ ℂ^{n_w}` the coefficient vectors (on the differential-basis words, modulo the `d = 3` relations `𝓡(Δ,J)`) of the
parity-even structures `⟨O_1O_2O_{Δ,J}⟩` conserved at points 1 and 2 (`O_i` conserved of spin `l`); `T(Δ,J): c_w ↦ c_w b^{(0)}/b_w` the word-wise
Lorentzian shadow map, `b_w/b^{(0)} = r_{k_w}(Δ,J)`, `Δ_34 = 2k_w`.

**Proposition [D].** (i) `r_k(Δ,J)` is rational and `r_k(Δ,J) r_k(3−Δ,−1−J) = 1` (from `Γ(u+k)Γ(u−k)/Γ(u)²` with `u ↦ 1−u` under the Weyl
reflection and the reflection formula); hence `T(Δ̃)T(Δ) = 𝟙`. (ii) `T(Δ,J)` maps `V(Δ,J)` onto `V(Δ̃,J̃)` and `𝓡(Δ,J)` onto `𝓡(Δ̃,J̃)` (the Lorentzian
shadow transform at the exchanged point preserves conservation at 1, 2 and maps the zero function to zero; given the scalar `b^{(Δ_34)}`).
(iii) **There is a basis of rational families `v_a(Δ,J)` of `V` with `T v_a(Δ,J) = v_a(3−Δ,−1−J)`**, i.e. the shadow matrix is
`(1/b^{(0)}) × (Weyl reflection) × 𝟙`; it is unique up to `GL_n` matrices with **Weyl-even** rational entries. Proof: `Φ[c](Δ,J) := T(Δ̃,J̃) c(Δ̃,J̃)` is a
semilinear involution of the space of rational families over the field `ℂ(Δ,J)` with respect to the Weyl reflection `σ`; by Speiser's lemma
(Hilbert 90 for `GL_n`) its fixed points span, and `c + Φ[c]` is fixed for any family `c`. (iv) **The exact helicity structures can be chosen
inside such a basis**: the helicity operators at points 1, 2 commute with the shadow transform, so `T` maps each helicity line `𝓛_h(Δ)` to
`𝓛_h(Δ̃)`; on a line, `T e_h = κ e_h(Δ̃)` with `κ(Δ)κ(Δ̃) = 1`, so (Hilbert 90 for `GL_1`) `κ = f(Δ̃)/f(Δ)` and `f e_h` is Weyl-covariant. Any other
normalisation of the helicity structures has a "helicity-dependent shadow coefficient" `f_h(Δ̃)/f_h(Δ)` — this is what CHL's Euclidean
`(2−Δ)_{|h_1+h_2|}/(Δ−1)_{|h_1+h_2|}` is. In a Weyl-covariant helicity basis the shadow transform acts diagonally with the common
coefficient `1/b^{(0)}`, and helicity-dependent shadow factors are of the form `f_h(Δ̃)/f_h(Δ)`. **This is bookkeeping** (referee M1): it holds
for any involutive word-wise map and contains no physics; the physics is in the normalisation, (v). (v) Consequently (§6.1, given Assumption T and the two-pair Lemma S) the heavy Regge moments of
`W[v_a, v_b]` have no `O(1/m)` term; the question of `O(1/m)` for physical couplings reduces to whether the **physical normalisation** of the
helicity structures is Weyl-even at order `1/m`. For `l = 1` it is (§4.6, exact Fock data in our own conventions). Caveat: (iv) needs one-dimensional helicity lines in the parity-even
conserved space (true for generic `J ≥ 2l`; low `J` and the parity-odd sector are not discussed).

Verification: (i)–(iv) explicitly and exactly for `l = 1` (§2, all closed forms); for `l = 2` (ii) numerically to `10⁻⁹` (§7) and (iii) by the
numerical Hilbert-90 construction (§7). Not done: (iv) for `l = 2` explicitly (CHL give `x`-space helicity structures only for two currents),
and the `l = 2` MFT normalisation (CHL's `free data` covers it [Q], but the map `T^{E,h} ↔ v_a` for stress tensors was not identified).

## 9. Status list and what remains

| # | statement | mark |
|---|---|---|
| 1 | Differential basis with two spinning operators per vertex; `⟨JJO⟩`: 5 independent transverse words, 2-dim conserved space `(v_1,v_2)` | [D], exact |
| 2 | Shadow map on conserved `⟨JJO⟩` structures = `(1/b^{(0)})` × Weyl reflection in `(v_1,v_2)`, exactly | [D] (exact symbolic) |
| 3 | Exact CHL helicity structures `= v_1 + f_h v_2`, `f_opp = J(J+1)(m²+2)`, `f_same = −(Δ−1)(Δ−2)(m²−2)`, Weyl-even | [ID] exact at 6 checked points (2 reflected; referee +3); independently reproduced by the momentum-space helicity projection (`tttt/ft_helicity.py`) |
| 4 | Exact `⟨JJJJ⟩` partial waves (25 word pairs): conservation at three points, pair exchange, scalar reduction, two-pair Lemma S | [DV] |
| 5 | Heavy vertex closed forms `V_1, V_2`; helicity labels | [ID], pointwise `∝ 1/m` to `m = 5·10⁵` |
| 6 | Factorised leading formula (§4.2), absolute per unit coupling in the `v`-basis, all 81 modes, `k = 4, 6`, ν = 0, 1, 3, β ≈ 0.4, 1 | **[DVE]**: ratios `0.96 … 1.009` at `m ≈ 100`, `1 + c_2/m²` |
| 7 | `±` structure incl. the opposite-sign `g_1²S_LS_R(|c^{++}|²−|c^{+−}|²)` term | **[DVE]** (same range as 6) |
| 8 | Exact MFT data in the helicity basis (diagonal; closed forms); absolute constant `1/4` per helicity quadruple; `C^± = (2P/m^{2k})(|c_same|² ± |c_opp|²)`; normalisation Weyl-even | [D from exact Fock data, own conventions] (§4.6); block-level two-pair link: §4.6 |
| 8b | Physical sign of the opposite-sign term: `+(|c_same|²−|c_opp|²)` for the covariant odd⊗odd kernel = flat `σ_{++} = +1, σ_{+−} = −1` | [D] + DVE (§4.5); functional ↔ flat sum rule identification [I] |
| 8c | Parity-odd exchanges in `⟨JJJJ⟩` | [ND] (referee M4) |
| 9 | No `O(1/m)` for parity-even exchanges in `⟨JJJJ⟩` | [D] (reflection + Lemma S + Assumption T) + [NE] fits vs mutations (§6.2; weak spots stated) |
| 10 | General: Weyl-covariant helicity bases exist for any `l` (bookkeeping, no physics by itself; needs 1-dim helicity lines, i.e. `J ≥ 2l` generic) | [D] (Hilbert 90 + commutation with helicity; uses CHL's conformal invariance of helicity [Q]) |
| 11 | `⟨TTO⟩`: 14 words, one 3d relation, 2 conserved structures; word-wise `b^{(Δ_34)}` shadow map consistent (`Δ_34 = 0, ±2, ±4`) | [NE] at `10⁻⁹` (three `(Δ,J)`) |
| 12 | `⟨TTTT⟩` heavy integrand with Weyl-covariant families: helicity content, factorisation, Weyl symmetrisation | see §7 |
| 13 | `⟨TTTT⟩` end-to-end on exact blocks | **[ND]** (jets of order 8) |
| 14 | `l = 2` exact helicity structures and MFT normalisation in the Weyl-covariant basis | [ND] |
| 15 | MFT decomposition of `⟨JJJJ⟩_MFT` with positive coefficients; absolute MFT constant (`2P/m^{2k}`-type) for two pairs | [ND] |
| 16 | `k = 2, 3`; `ν ∼ m`; other Δ of the external operators (only conserved currents here) | [ND] |

**What remains, in order of use for the graviton bound.** (1) `⟨TTTT⟩` end-to-end: the doubly-spinning words need jets of order 8
(`e2e/pq_jets` and the radial derivative sums are hard-wired to order 4) and, per `SPINNING_E2E.md` §8, extended precision; 196 word pairs
(14 × 14, one 3d relation per vertex). (2) The `l = 2` helicity map `T^{E,h} ↔` Weyl-covariant basis (CHL do not give `x`-space `⟨TTO⟩`
helicity structures; a momentum-space evaluation of their eq. `Tpm` would do) and with it the `l = 2` MFT normalisation — the only missing
step for "no `O(1/m)`" in `⟨TTTT⟩` with physical couplings. (3) The absolute MFT constant for two pairs (expected `1/16` per helicity
quadruple × the identical-particle factor), best by an exact Fock-space computation as in `SPINNING_MFT.md` §2, which would also replace
CHL's "tentative" doubling. (4) `k = 2, 3` (the graviton sum rules), for which the radial representation is inadequate
(`SPINNING_E2E.md` §11 item 1).

## 10. Code index (`ads4/code/fourpt/`; `sh run_all.sh`)

| file | content | output |
|---|---|---|
| `embed2.py` | embedding algebra with two spinning operators per vertex (`D_11, D_12, D_22, D_21, H_12`), 3- and 4-point words | — |
| `structures2.py`, `t1_structures_l1.py` | words, conservation (symbolic), conserved basis `v_1, v_2`, exact shadow map | `out_t1_structures_l1.txt` |
| `chl_map.py`, `t1b_helicity_l1.py` | CHL helicity structures inside the conserved space, exact rational identification and closed forms | `out_t1b_helicity_l1.txt` |
| `spinblock2.py`, `gen_l1.py` | doubly-spinning partial waves (25 word pairs; generated code `cache/gen4_l1.py`) | — |
| `t2_validate_l1.py` | conservation, pair exchange, scalar `(H,H)` reduction | `out_t2_validate_l1.txt` |
| `vertex2.py`, `t3_vertex_l1.py` | exact shadow integrand, heavy vertex closed forms, kernel conventions | `out_t3_vertex_l1.txt` |
| `moments4.py`, `tmoments4.py`, `t4_compute_l1.py` | two-pair kernels, modes on the moment grid (cache `cache/modes4/`) | — |
| `predict4.py`, `t5_analysis_l1.py` | factorised prediction; component fits, `±`, `1/m` fits, mutations, Lemma S | `out_t5_analysis_l1.txt` |
| `t8_mft_l1.py` | MFT normalisation of the helicity components | `out_t8_mft_l1.txt` |
| `build_l2.py`, `t6_structures_l2.py`, `t7_vertex_l2.py` | `⟨TTO⟩` words, 3d relation, conserved space, shadow-map consistency; `⟨TTTT⟩` heavy integrand | `out_t6_*.txt`, `out_t7_*.txt` |

## 11. Revisions after the referee report (`ads4/REFEREE_SPINNING_4PT.md`)

* M1: "no mixing / CHL eigenvalues are artefacts" reframed as Hilbert-90 bookkeeping valid for any basis (§0 items 2, 6; §2.2, §2.3, §8); the
  physical content moved to the MFT normalisation.
* M2: the physical sign derived with covariantly defined pair kernels (§4.5): `+(|c_same|² − |c_opp|²)` for `D_5⊗D_5`, matching flat
  `σ_{++} = +1`, `σ_{+−} = −1` and the `SPINNING_E2E.md` §8 conjecture; the minus sign of §4.3 is the inversion-referred frame of `S_R`
  (`det I = −1`). Verified on the exact moments (`tttt/out_t4_sign_l1.txt`). Not derived: which AdS functional uplifts which flat sum rule.
* M3: the MFT normalisation computed exactly in our conventions by a two-particle Fock space (§4.6): helicity-diagonal, closed forms on all 25
  `n ≥ 1` primaries, absolute constant `1/4` per helicity quadruple, `C^± = (2P/m^{2k})(|c_same|² ± |c_opp|²)`, all factors Weyl-even.
* M4: parity-odd exchanges stated as untested (§0 item 5, §6.3, §10 row 8c).
* m1 label `shadow T`; m2 "6 checked points"; m3 one-scale and uniform-weight residuals quoted; m4 `Σ_−` in `W[v_2,v_2]` not independent; m5 `|c_1|`
  range over ansätze (`≤ 0.2`); m6 the `c_1 = −2.6` is a class-fit artefact (raw modes `−0.009`); m7 `α²` cancellation stated; m8 Richardson probe
  named and `⟨TTTT⟩` claims weakened (§0 item 7, §7.4).
