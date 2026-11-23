# H_K residual checks: the k = 4 EE heavy image, spinning K1′, ε-DS in the Tonelli argument, and (P) at J = 2

*23 Nov 2026. Unreviewed, not committed. No existing file was edited. Answers the remaining items of EESTEP_CLOSE §7 (2, 3, 5) and
DESHIFT_CLOSE §7.4, as far as they bear on H_K.*
* *New code:* `ads4/code/hkres/`: `hk1_heavy.py`, `hk2_k1spin.py`, `hk3_epscert.py`, `hk4_p_j2.py`.
* *Outputs:* `ads4/code/hkres/cache/` (git-ignored by its own `.gitignore`).
* *Read-only imports:* `eestep/{ee6_lrprofile, ee_quad, ee_blocks, ee3_deshift}`, `h7i14/v3_regge_kernel`, `spinfunc/{s1_offdiag_kernel,
  bg_eeq}` (bg_eeq's cache directory redirected), `resultq/{stage4_hj, og_free}`, `repro_engine/sdchain`, `deshift/kern`.
  Data: `resultq/cache/stage4_hj.pkl`, `resultq/functionals_q4_deshift.json`, `prl/data/q_deshift_continuum_{c4v,c4d}.json`.
* *Compute:* ≈ [[CPU]] CPU-minutes, at most 2 processes, one thread each, < 3 GB each.

Grades: **[P]** proved; **[D]** derived; **[DV]** derived and checked numerically; **[NE]** numerical evidence; **[M]** model;
**[I]** open; **[certified]** interval-rigorous (ball arithmetic with explicit truncation and aliasing bounds). A chain carries the grade of
its weakest link.

---

## 0. Verdict

| item | before | after |
|---|---|---|
| 1. Heavy image of the k = 4 EE position kernel; EE3 ν-normalisation; f_pos(±3i/2) | LRC-a EE k = 4 [D]; normalisation "rests on the same code path, not rechecked" | **Confirmed. LRC-a EE k = 4 → [D + DV].** EESTEP's b-shifted leading-Regge kernel equals the split kernel that Stage 4 transformed, with ratio **1 to 1.8·10⁻⁷**, ν-, β- and polarisation-independent (4 pol. sets × β = 1, 3, 6 × ν = 0.5–6). Its rank-3 plane-wave transform, done numerically and independently of stage4_hj's symbolic route, reproduces the stored R_j^{EE3} (and R_j^{sc3}) **to ≤ 1.1·10⁻¹¹**. So the position weight is χ·const/N_com, with the scalar's const, and **f_pos(±3i/2) = 0 stands** (§1) |
| 2. Spinning K1′, numerical | scalar examples only | **[D + DV]** on the exact free-scalar ⟨TTTT⟩ (N = 1, 10; c = 0.5, 0.9; 63 polarisation sets; box down to 10⁻⁶): Cauchy–Schwarz holds everywhere (saturated only where it must be), B(c) is finite and sits at the corner (c, c), \|𝒢_Regge\| ≤ 0.53–0.999 B(c), \|dDisc 𝒢\| ≤ 0.007–0.25 · 2B(c), and the edge row \|M − 1\| ≤ (B − 1)(min(ρ, ρ_b)/c)^{τ_min/2} holds over the whole box with τ_min = 1 (§2) |
| 3. DS for Q in the Tonelli argument | "DS [NE]", the only non-derived input of H_K for Q | **For the explicit continuum weights the swap needs no de-shift identity.** The analyticity that HK derived from DS (Lemmas S2/S4) holds by construction [P], and the criterion is unchanged (J₀ < 5/2) [D]. ε-DS fixes only the values (images and moments = targets ± ε). ε enters a Tonelli bound only additively, through a finite C_E. The bounds move by ≤ ε × the factors of §3.2 (E). **ε_W is now [certified] for C4:** [[EPS]] (§3.3) |
| 4. (P): C5 piece and J = 2 | J = 0 derivative piece only | [[P4]] (§4) |

**H_K.**
* **Result D: [D]** (unchanged; DS for k = 2 is [D + NE] through Lemma S2 and the cutoff).
* **Result Q: [D]** for Ω^Q_up *defined* by the published continuum weights (DESHIFT §2.4; equal to the Stage-7 XP weights to 10⁻¹⁶).
  * The swap for the C4 pieces uses the weights' explicit analytic structure [P], not DS.
  * The C5 piece is raw Fubini (J₀ < 3/2).
  * Q's Ĉ′_EE (k = 2, ψ₂) is in the same position as D's: DS [D + NE] via Lemma S2 and the cutoff.
  * **ε-DS is an input of the bound values, not of H_K.** It is [D + certified] for C4 (§3.3) and [DV] for C5 (2.4·10⁻⁷, per-state
    numbers only).
  * If Ω^Q_up is instead kept as the Stage-7 *grid* weights, H_K is [D given DS] as in EESTEP, with DS [NE].
* Criteria unchanged: J₀ < 3/2 (k = 2, 5, Ĉ₃), **J₀ < 5/2** (k = 4, Q's weights). Under J₀ ≤ 1 the margins are ½ and 3/2.

---

## 1. The heavy image of the k = 4 EE position kernel (`hk1_heavy.py` → `cache/hk1_heavy.txt`)

### 1.1 What has to be checked
* EESTEP's f_pos argument: the position weight of `Ĉ₄^{ds}` is `χ·const/N_com`, so a pole of the weight at 3i/2 is a pole of χ.
  This presupposes that the heavy image of the EESTEP position kernel `K^EE_ν = N(ν) r² E_ν` is `const·N_com Σ_j R_j^{EE3} P_{l+j}` with
  a **ν-independent** const, equal to the scalar one. A hidden factor g(ν) with a zero at ±3i/2 would flip the conclusion.
* The heavy (plane-wave) map is the transform of the split representation `K₀ = ∫dz₀ A(x̂, z₀) B(ŷ, z₀)` (v3 §18.1). Stage 4 applied it,
  at rank 3, to the monomials of `A_EE = 16(z₂′·z₄′)² A_sc` and of `B_EE` = v3's θ-average of τ₁₃. EESTEP never compared its kernel with
  that split kernel.

### 1.2 Link (a): position space [DV]
Prediction from the split form, in EESTEP's frame and scalar normalisation (`E = K/(N r²)`, scalar `L = P_l`):

    L_ν(β; pols) / P_l(cosh β)  =  16 (z₂·z₄)² · ∫dφ (−z₀·x̂)^{δ−1} B_EE(φ) / ∫dφ (−z₀·x̂)^{δ−1},    δ = ½ + iν,

with `B_EE` evaluated from `v3_regge_kernel.B_monomials(spin=True)` at `ŷ`, with EESTEP's polarisations w₁, w₃.
`L_ν` is EESTEP's own object (`ee6_lrprofile.L`: b-shifted blocks, r⁰ Cauchy coefficient).

| | β = 1 | β = 3 | β = 6 |
|---|---|---|---|
| max \|ratio − 1\|, 4 polarisation sets × ν ∈ {0.5, 1.5, 3, 6} | 3·10⁻⁸ | 5·10⁻¹⁰ | 1.8·10⁻⁷ |

* The ratio is **1**, not merely constant: the ν-only factors that v3 dropped are exactly the scalar's, and EESTEP's dressing
  `D = [(H₂₄ − 2V₂V₄)/(X₂₀X₄₀)]²` reduces at x₀ = 0 to `16(z₂·z₄)²` (and z′₂·z′₄ = z₂·z₄, the inversion being an isometry).
* The β = 6 entries need 2¹⁷ φ-nodes: `(−z₀·x̂)^{δ−1}` has width ≈ 2e^{−β} (a first run with 512 nodes gave O(1) errors at β = 6 only).

### 1.3 Link (b): the heavy image itself [DV]
* **Method.** `L₄(β) = −16 ∫dz₀ ⟨∂³Â, ∂³B̂⟩` (stage4 §1.1), but evaluated numerically: a truncated 5-variable Taylor series (orders 2, 2, 3
  in e·∂, f·∂, ∂_μ) of `S^{−λ}R^κ` on a 192-point z₀ grid, the full η-contraction, and a least-squares fit
  `L₄ = Σ_{j=±1,±3} h_j P_{l+j}(cosh β)` over 14 β ∈ [0.15, 3].
* **Independence.** It uses none of stage4_hj's equivariant decomposition, Gram contraction or `to_harmonics`. It shares only the
  monomials and the transform constants K(c, λ).

| ν | fit residual (sc / EE) | max_j \|(h_j/h_{−1}^{sc}) / (R_j/R_{−1}^{sc3}) − 1\|, scalar and EE |
|---|---|---|
| 0.37 | 1.9·10⁻¹⁵ / 3.3·10⁻¹⁵ | 1.1·10⁻¹¹ |
| 1.3 | 1.7·10⁻¹⁵ / 2.3·10⁻¹⁵ | 2.1·10⁻¹² |
| 2.9 | 2.1·10⁻¹⁵ / 2.7·10⁻¹⁵ | 8.8·10⁻¹³ |
| 5.5 | 2.8·10⁻¹⁵ / 2.7·10⁻¹⁵ | 7.7·10⁻¹³ |

The per-j ratios `h_j^{EE}/h_j^{sc}` match the exact rationals to all printed digits (e.g. ν = 1.3, j = −3: −1.7023 + 5.08097i).

**(a) + (b):** the heavy image of the EESTEP k = 4 EE position kernel is `κ₀ N_com(ν) Σ_j R_j^{EE3}(ν) P_{l+j}(cosh β)`, with κ₀ the
ν-independent scalar constant. **LRC-a for EE at k = 4: [D] → [D + DV].** The ν-normalisation identification of EESTEP §3: **[DV]**.

### 1.4 The f_pos conclusion (exact algebra, `hk1` (c))
* `R₋₃^{EE3}(−3i/2) = −524288`, `R₋₁^{EE3}(−3i/2) = 6193152/25` (EESTEP's values, reproduced). `R₋₃^{sc3}(−3i/2) = 0`.
* The scalar normalisation `n(ν) := R_j^{sc3}/(a_j(l)E₉(ν − ij))` is j-independent (checked at real and complex ν) and equals
  `1024/(Γ(11/4 − iν/2)²Γ(11/4 + iν/2)²) = 2²¹π²/N_com`. So the scalar position coefficient is a pure constant, and
  `N(ν)/N_com ∝ 1/N_com` is regular and non-zero at ±3i/2 (`n(±3i/2) = 92.714…`).
* Hence the position weight of Q's C4 pieces is `χ·const/N_com` exactly, and f_pos(±3i/2) = 0 ⇔ χ regular at −3i/2.
  * **Stage-7 XP weights:** [D given DS] (HK Lemma S4) + [NE] (EESTEP §3), as before; now with the normalisation [DV].
  * **Continuum weights (DESHIFT §2):** **[P]**. `χ = Q·B` with B entire and Q's poles at ±3i, ±i/2, ±5i/2, ±i(3 + 2√2)/2 only; there is
    no pole at ±3i/2 by construction.
* **k = 4 criterion: J₀ < 5/2** stands.
* *Side remark.* `R₊₁^{EE3}(−3i/2) = 0` and `R₊₃^{EE3}` has a pole at −3i/2. Neither matters: S4 concerns `χR₋₃` in −2 < Im ν < 0.

---

## 2. Spinning K1′ (`hk2_k1spin.py` → `cache/hk2_k1spin.txt`)

### 2.1 Set-up
* **Correlator.** The exact free-scalar ⟨TTTT⟩ in d = 3 (conformal coupling, N copies, disconnected O(N²) + connected O(N)), Wick engine
  of `og_free` (link derivatives, vertex, matchings). A unitary CFT with spin-2 external operators. (Its J₀ is not < 1; K1′ is the
  Rindler bound and does not need that.)
* **Frame (as lem2).** Light-cone coordinates in the (t, x¹) plane: `x₁ = (1, 1) = −x₂`, `x₄ = (ρ, ρ̄) = −x₃`. Regge sheet:
  `0 < ρ < 1 < ρ̄`, `ρ_b = 1/ρ̄`, box `ρ, ρ_b ≤ c`. Only (14) and (23) are timelike.
* **Wightman orderings and Rindler positivity.** With Θ the CRT reflection `(t, x¹, x²) → (−t, −x¹, x²)`, `z → R z*`,
  `R = diag(−1, −1, 1)` (so x₂ = Θx₁, x₃ = Θx₄), the four orderings are ⟨Θ[P]Q⟩ with P, Q ∈ {O₁O₄, O₄O₁}:

  | ordering | i0 of (14), (23) | sheet |
  |---|---|---|
  | ⟨O₃O₂O₁O₄⟩, ⟨O₂O₃O₄O₁⟩ | (+, +), (−, −) | Regge (ρ̄ ∓ i0) |
  | ⟨O₂O₃O₁O₄⟩, ⟨O₃O₂O₄O₁⟩ | (+, −), (−, +) | first |

* **Bound.** `|⟨Θ[P]Q⟩|² ≤ ⟨Θ[P]P⟩⟨Θ[Q]Q⟩`. For polarisations (z₁, z₂, z₃, z₄), Q carries (z₁, z₄) and P carries (Rz₂, Rz₃) at (x₁, x₄).
* **Normalisation.** Every value is divided by `√(N_P N_Q)`, where N is the identity ((12)(34)) part of the positive comparison
  correlator. Then `|𝒢| ≤ √(𝒢_P𝒢_Q) ≤ M := ½(𝒢_P + 𝒢_Q)`, with `𝒢_{P,Q} ≥ 1`. dDisc is `½(first-sheet pair) − ½(Regge pair)`.
* **Polarisations.** z₁, z₂ = (1, cos a, sin a); z₃, z₄ boosted with their pair (the frame of x₄), as K1′/EESTEP prescribe. Lab-frame
  z₃, z₄ are run as a control.

### 2.2 Results [DV]
* **Scalar control.** O = φ² through the same engine. The first-sheet orderings equal the all-spacelike correlator at (ρ, ρ_b) to
  1.2·10⁻¹¹, i.e. the comparison correlators are first-sheet correlators at the same cross-ratios, as K1′ states.
  `|𝒢_Regge| ≤ 𝒢_E(ρ, ρ_b)` pointwise.
* **Spin 2** (63 polarisation sets, 61 × 61 log grid, ρ, ρ_b ∈ [10⁻⁶, c]):

  | c, N | max \|𝒢\|²/(𝒢_P𝒢_Q) | B(c) = sup M (at) | sup \|𝒢_Regge\| / B | sup \|dDisc\| / 2B | max (M − 1)/((B − 1)(min(ρ,ρ_b)/c)^{½}) | edge slope |
  |---|---|---|---|---|---|---|
  | 0.5, 1 | 1.000000000000 | 3.51·10⁵ (0.5, 0.5) | 0.53 | 0.25 | 1.0000 | 0.50 |
  | 0.5, 10 | 1.000000000000 | 2.71·10⁵ (0.5, 0.5) | 0.94 | 0.033 | 1.0000 | 0.52 |
  | 0.9, 1 | 1.000000000000 | 2.19·10¹⁵ (0.9, 0.9) | 0.987 | 0.0066 | 1.0000 | 0.52 |
  | 0.9, 10 | 1.000000000000 | 2.18·10¹⁵ (0.9, 0.9) | 0.999 | 0.0007 | 1.0000 | 0.66 |

  * Cauchy–Schwarz is never violated. Equality occurs only where it must: for the polarisation sets with P ∝ Q.
  * No comparison value falls below 1, and the identity parts N_P, N_Q are positive everywhere.
  * **B(c) is finite and attained at the corner (c, c).** Its size is the t-channel identity `(u/v)^Δ` at z = z̄ = 4c/(1+c)², as for
    GFF with Δ = 3 (lem2). K1′ needs finiteness only.
  * **Edge row of K1′.** `|M − 1| ≤ (B(c) − 1)(min(ρ, ρ_b)/c)^{τ_min/2}` holds on the whole box with τ_min = 1 (φ² and the twist-1
    currents), with equality only at (c, c). The local slope at the edge is 0.50.
  * **Uniform in η.** At fixed ρ = c, M and |𝒢_R| decrease monotonically to 1 as ρ_b → 0 (η up to 4.7·10³).
  * **Lab-frame polarisations** (control): identical B(c), edge ratio and slope. With the self-normalisation √(N_PN_Q) the bound is
    insensitive to the polarisation frame. The frame only changes the normaliser.
* **Grade: K1′ [D] → [D + DV]** (spinning check on a free CFT). Not exercised: an interacting spinning correlator. The σ⊗T_free
  product of EIN would add one; its σ part is a Witten-diagram sum with its own continuation, and it was not run.

---

## 3. DS for Q, stated cleanly: ε-DS in the Tonelli argument

### 3.1 Statement
Define Ω^Q_up's k = 4 de-shifted pieces by the explicit continuum weights `u_p` (p = C4v, C4d; `prl/data`). These are
`u = Q·B` with Q rational (prescribed poles at ±3i, ±i/2, ±5i/2, ±i(3+2√2)/2 and zeros at 0, ±i, ±2i), B an entire Gaussian-RBF sum,
super-Gaussian decay beyond ν ≈ 60. By DESHIFT §2.1 (shift identity, [D]), for all t ≥ 0,

    |W[u_p](t) − Φ_p(t)| ≤ ε_p,       Φ_{C4v} = 1,  Φ_{C4d} = −t²/2.                                  (ε-DS)

Then:
* **(i) The swap (H_K proper) does not use the de-shift identity at all.** Absolute convergence and the interchange follow from K1′,
  K2ν, K3, LRC-a and LRC\* together with the analytic structure of u. That structure was *derived* from DS for the grid weights (HK
  Lemmas S2, S4); for the explicit weights it holds **by construction [P]**. The criterion is unchanged: J₀ < 5/2 (k = 4).
* **(ii) ε-DS fixes the values.** It identifies the per-state images and the light-side moments with their targets, up to ε_p
  (images) and `M_{n,p} = ∫|E_p,even|μ^{2n}` (moments). If one insists on writing the kernel as "target + remainder", ε enters the
  absolute-convergence bound only additively, through a finite constant C_E (§3.2 (C)).
* **(iii)** The bounds move by at most ε_p times the factors of §3.2 (E).

**Grade: H_K for Q's C4 pieces [D given K1′, K2ν, K3, LRC\*, LRC-a], with no DS input. The values need ε-DS: [D + certified] for C4
(§3.3), [DV] for C5.**

### 3.2 Derivation [D]
Write `E_p(μ) = Σ_j [u_p(μ + ij)R_j(μ + ij)] − S_j` (DESHIFT §2.1), so that `W − Φ = ½∫E_p cos μt` and `ε_p ≥ ∫₀^∞|E_p,even|`.

**(A) Kernel.** `w = u·const/N_com` (§1: the const is the scalar's, ν-independent). We have `|u| ≤` (rational) × (super-Gaussian) and
`N_com⁻¹ ~ e^{πν}ν^{−9}`. So by K2ν the position kernel `K_W = ∫w C_ν dν` converges absolutely at every point of the CLS-D domain. The HK
cutoff is not even needed: it changes u by < e^{−300} (DESHIFT §2.4).

**(B) Leading Regge: growth from analyticity, value from ε-DS.**
* In each j-term of the LR image (heavy, `Σ_j∫χR_jP_{l+j}`; position, `∫w ηⁿ Σ_k p_k P_{l+k}`), shift the contour by ij.
  * u is meromorphic with the poles listed. Those at ±i/2, ±5i/2, ±i(3+2√2)/2 sit at zeros of every crossed R_j (heavy) and at zeros of
    c_l or decaying edge terms (position: ±i/2, ±5i/2 → l = −1, −3 where c_l = 0; 2.914i → η^{−0.41}, EESTEP §3).
  * The shifted lines carry no growth. The target pole ±3i gives the target profile.
  * The resulting bound `|I_W| ≤ (|target| + Σ_j∫|f_j|)P_{−½}(cosh β)` holds without any smallness of E. This is what the corner count
    uses; Lemma PW's raw regrowth is excluded because u is analytic, not because the identity holds.
* For the EE k = 4 kernel the same holds per tensor component. Its polarisation dependence sits in ν-independent C_b(x), and every edge
  coefficient is c_l × (B^{(b)}, F_b), analytic in the strip apart from 3i/2 (+2p) where c_l vanishes (EESTEP §1.4). §1 makes w's
  analytic structure that of χ.
* **Value.** Mehler–Dirichlet, `I_W(β) − I_Φ(β) = (√2/π)∫₀^β (W − Φ)(t)(cosh β − cosh t)^{−½}dt`, gives for **all β**:

      |I_W(β) − I_Φ(β)| ≤ ε_p P_{−½}(cosh β) ≤ ε_p .

  * For C4d the residue part is exact only to the certified ball sizes of §3.3. These add `d₀ + d₂t² + ½|a₂|(√(2π) + πt)` to
    `|W − Φ|`, i.e. at most `(d₀ + d₂β² + |a₂|(1.3 + 1.6β))P_{−½}` to the image. All are ≤ 5·10⁻³⁵.

**(C) If the kernel is split as "target + remainder".** After the shift, the remainder of the LR position profile is a real-line integral at
harmonic label μ. Its coefficient is the heavy residual E_p(μ) divided by the heavy map at μ: `E_p(μ)/(κ₀E₉(μ))` for the scalar-type map,
and the same fixed linear map at label μ for EE k = 4 (§1). So

    |K_W^{LR} − K_target^{LR}| ≤ C_E r^a ηⁿ P_{−½}(η)·(de-shifted power),      C_E := ∫|E_p(μ)|/|κ₀E₉(μ)| dμ < ∞.

C_E is finite, since E_p decays super-Gaussianly. It is **not small**, since `1/E₉ ~ e^{πμ}μ^{−14}`, and it is not a bound on any computed
number. This is the only way ε enters a Tonelli bound, and it does so additively.

**(D) Beyond leading Regge, and the swap.**
* LRC\* uses only the analyticity of `w c_l` in the strips, apart from target singularities, and regularity at 3i/2. Both are explicit
  here. So the (exact − LR) remainder has the exponents of EESTEP §2: k = 4 shell `η^{½}`.
* The corner count therefore gives
  `Σ_O ∫|K_W||dDisc G_O| ≤ 2B(c)∫_box |K_W| + (K2 region) < ∞` for J₀ < 5/2.
* By K3 and Tonelli, `Σ_O Tr[P_O Ω_W[G_O]] = −Ω_W[light]`, with absolute convergence. **This holds exactly for Ω_W. The de-shift
  identity is never used in it.**

**(E) From Ω_W to the bounds: where ε_p acts.** Write `a_p` for the LP coefficients of the pieces, and N_O for the per-state
normalisation (forward image = 1).
* **Heavy side, per state:** `|Ω_W[G_O] − Ω_target[G_O]|_{LR} ≤ Σ_p |a_p| ε_p P_{−½}(cosh β_O) ≤ Σ_p |a_p| ε_p`.
  * Per-state positivity margins m_O (binding FF/N ≥ 0.1427; LR minima ≥ 0.9866) survive iff `Σ_p|a_p|ε_p < m_O N_O`.
  * DESHIFT §2.5 evaluated the worst-sign shift over all 8 sign patterns on the full LR/HYB scan: **≤ 6·10⁻¹¹** (dominated by C5's
    2.4·10⁻⁷; the C4 pieces alone give ≤ 10⁻¹²).
* **Light side:** `|m_{n,p} − m_{n,p}^{target}| ≤ M_{n,p} := ∫₀^∞|E_p,even|μ^{2n}dμ` (M_{0,p} = ε_p). So
  `|ΔL| ≤ Σ_p |a_p| Σ_n |λ_n| M_{n,p}`, with λ_n the EFT coefficients of degree 2n (Λ^{−2n}-suppressed). For a coupling g fixed by
  `L(g) = 0`, `|δg| ≤ |ΔL|/|∂L/∂g|`.
  * **R⁴:** `2g₄W(0)`, `|W(0) − 1| ≤ ε_{C4v}`, so the relative shift of the R⁴ coefficient is ≤ ε_{C4v}.
  * **g₅ and |ĝ₃|² (n ≤ 1):** `M_{1,C4d} = 2.1·10⁻¹⁰`, `M_{1,C5} = 8.3·10⁻⁶`, in Λ⁻² units [DV, DESHIFT; M_{n≥1} not certified].
* **Result.** Every bound moves by **≤ Σ_p|a_p|ε_p / min_O(m_O N_O)** through the heavy side and by
  **≤ max(ε_{C4v}, M_{1,C4d}/Λ², M_{1,C5}/Λ²)** (relative) through the light side. Numerically: ≤ 6·10⁻¹¹ on the margins, ≤ 10⁻¹² (C4) and
  < 10⁻⁶ (C5, the g₅ cancellation) on the light side. These are DESHIFT §1.2 and §2.5's numbers, now with the C4 ε certified. **No
  computed number moves.**

### 3.3 ε_W made interval-rigorous for C4 (`hk3_epscert.py` → `cache/hk3_epscert_{c4v,c4d}.txt`) [certified]
DESHIFT's ε_W were 384-bit **midpoint** sums of |E_even| on a 0.01 grid. Certification turned out to be mechanical, but not trivial, for
two reasons.
* **Cancellation.** The RBF sum B cancels by about 10³⁰ at |Im z| = 3, and E is a further cancellation of O(1) terms down to 10⁻¹³. At
  192 bits the balls blow up beyond μ ≈ 45. The run uses 384 bits, as DESHIFT did.
* **Termwise enclosures.** Interval evaluation of E over a region cannot see the cancellation: it gives M ≈ 10⁵–10¹³. It is therefore used
  only inside Cauchy-type remainders, where it is multiplied by `(r/ρ)^K ≤ 10⁻²⁸/M`.

**Method.**
* **Definition.** The functional is the published one: `χ = Q·B`, the 60-digit decimal coefficients, h and σ the IEEE doubles that
  ab_solve used. With exact decimals 0.15, 0.45 instead, the residue part moves by ~10⁻¹⁶.
* **Exact cancellation.** `4z² + (3 + 2√2)²` in Q against `4z² ∓ 12iz − 1 = (2z ∓ ia)(2z ∓ ib)` in `R_{∓3}` is done exactly (sympy).
  After it, E_even is analytic in |Im μ| < ½ (the only other near singularity, μ = 0, is removed by S_j).
* **(A) Laurent data at μ = 0.** Trapezoid rule on |μ| = 0.2 with 160 points, plus the aliasing bound from an enclosure of
  `max|μᵖf_j|` on |μ| = 0.4. This gives ball-valued a_k, residue part and a₂.
* **(B) Panels of width 0.15 on [0, 110].**
  * Taylor coefficients at the panel centre: K-point trapezoid (K = 54–80) on |z − m| = 0.09, with aliasing bound.
  * Truncation by Cauchy: `|c_n| ≤ M/ρⁿ`, ρ = 0.3. M is the maximum-modulus enclosure on the boundary of 2-unit chunk rectangles.
  * `∫_panel|E| ≤ Σ_n (|ĉ_n| + e_n)·2(w/2)^{n+1}/(n+1) + w M (w/2ρ)^{N+1}/(1 − w/2ρ)`.
* **(C) Tail μ > 110.** Explicit Gaussian/coefficient bound.

| piece | DESHIFT midpoint ε_W | **certified ε_W** | residue part: \|const − Φ(0)\|, \|t²-coeff − Φ₂\|, \|a₂ sum\| | tail > 110 |
|---|---|---|---|---|
| C4 value | 1.205·10⁻¹³ | [[EV]] | [[RV]] | [[TV]] |
| C4 derivative | 4.171·10⁻¹¹ | [[ED]] | [[RD]] | [[TD]] |

**Certified statement.** For all t ≥ 0, `|W(t) − Φ(t)| ≤ ε_W^{cert} + (residue terms)`, and `|image − target| ≤ ε_W^{cert}P_{−½}(cosh β)`
(+ the ≤ 10⁻³⁴ residue terms of §3.2 (B)).
* The analytic step, the shift identity of DESHIFT §2.1, is [D].
* The numerical inequality is [certified].
* **ε-DS for C4: [DV] → [D + certified]. C5 not certified.** Its pole at a quintic root needs a root-isolation certificate, and it is not
  load-bearing for H_K.

---

## 4. (P) completion (`hk4_p_j2.py` → `cache/hk4_p_j2.txt`)

[[SEC4]]

---

## 5. Proposed table rows (UT §3.1; not edited here)

| # | before | proposed |
|---|---|---|
| LRC-a, EE k = 4 | [D] | **[D + DV]**: EESTEP kernel = Stage-4 split kernel (ratio 1 ± 1.8·10⁻⁷); numerical rank-3 heavy image = R_j^{EE3} (≤ 1.1·10⁻¹¹) |
| f_pos(±3i/2) = 0 (Q's C4) | [D given DS + NE] | normalisation [DV]; **[P] for the continuum weights** (no pole at ±3i/2 by construction) |
| K1′ | [D] | **[D + DV]** (free-scalar ⟨TTTT⟩, generic polarisations, whole box, edge row with τ_min = 1) |
| DS for Q (H_K input) | [NE] | **not an H_K input** for the explicit continuum weights (analyticity by construction); **ε-DS is the input of the values: C4 [D + certified]** (ε_W ≤ [[EV]] / [[ED]]), C5 [DV] |
| H_K, Result Q | [D given DS] | **[D]** (Ω^Q_up defined by the continuum weights); values: ε-DS [certified] (C4) |
| H_K, Result D | [D] | [D] (unchanged) |
| (P) for g₅ | [D-sketch + NE] (J = 0, derivative piece) | [[P5]] |

## 6. Reproduce

    cd ads4/code/hkres; ulimit -v 3000000; export OMP_NUM_THREADS=1
    python3 hk1_heavy.py                     # ~1 min  -> cache/hk1_heavy.txt
    python3 hk2_k1spin.py                    # ~1 min  -> cache/hk2_k1spin.txt
    python3 hk3_epscert.py c4v & python3 hk3_epscert.py c4d; wait      # [[T3]] each -> cache/hk3_epscert_{c4v,c4d}.txt
    [[RUN4]]
