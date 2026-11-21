# Below the fold, part 3: the `+−` J = 4 gap state of Result Q, the "512 c₂" K₂ disagreement, and a closed-form spin-J contraction

*21 Nov 2026, below-fold-3 agent. Unrefereed, not committed. No existing document or code file was edited; nothing of the `oo_*` agent
was touched or run. New code `ads4/code/spinfunc/bh_*.py`; caches `spinfunc/cache/bh_*` (git-ignored). Read-only imports: `bg_eeq`,
`bg_ee`, `bg_geom`, `bg_ser`, `bg_qan`, `bg_scalar`, `bf_qan`, `bf_exact`, `n3_exact`, `l1_crystal`, `l6_spin_action`, `l7_dc2_analysis`,
`l8_analysis`, `l5_spin` (via bg_ee), `resultq/stage6_model`, `stage6_stress`, `stage6_scan`, `common_q`. Compute ≈ 3.5 CPU-h on ≤ 2 cores,
≤ 1.4 GB per process (`ulimit -v 3 GB`).*

Grades: **[P]** proved; **[D]** derived; **[DV]** derived and checked numerically; **[NE]** numerical evidence; **[M]** model; **[I]** open.
`c₂` = `m² ×` relative correction; `s_d = 0.04 Σ_ν χ_d(ν) Re Γ_K(ν)` (Stage 7 §2, BELOW_FOLD2 §4). "Engine" = the analytic-radial 1/Δ
engine of BELOW_FOLD2 (real section D), whose EE beyond-LO equivalence with the crystal cycle is **assumed** there (§5.4) and here.

---

## 0. Answer in brief

1. **`+−` J = 4 gap state measured [NE + M].** Helicity basis = the project's order-0 Laplace basis `T0(Δ)` (l8 rule) at the same Δ; it
   diagonalises the engine's image (reference: off-diagonal `3·10⁻⁹`, `+−/++` ratio ν-independent to `10⁻⁶` at a⁰ — the a⁰
   helicity-blindness, made explicit) [DV]. Result (`bh_qan4.py`, all four entries of the 2×2 prefactor matrix × inverse shadow matrix):

   | state | Δ_g | `s_d` H0 | H2 | H3 | H4 | SB | XP | FF/N δ = 0.05 (K=0 → meas.) | δ = 0.2 |
   |---|---|---|---|---|---|---|---|---|---|
   | **`+−` J = 4, gap `1.0134Λ`** (Δ = 304.03) | 300 | **−0.516 ± 0.039** | −0.517 | −0.516 | −0.515 | −0.532 | −0.522 | +0.079 → **+0.118 ± 0.003** (SB +0.119) | +0.229 → **+0.278 ± 0.004** (SB +0.279) |
   | **`+−` J = 4, gap `1.0040Λ`** (Δ = 1004.01) | 1000 | **−0.510 ± 0.004** | −0.511 | −0.510 | −0.509 | −0.525 | −0.516 | +0.059 → **+0.099 ± 0.0003** (SB +0.100) | +0.209 → **+0.259 ± 0.0004** (SB +0.261) |
   | `++` J = 4, same Δ (by-product) | 300 | +0.517 ± 0.017 | +0.517 | +0.517 | +0.516 | +0.532 | +0.522 | +0.553 → +0.533 | +0.743 → +0.721 |
   | `++` J = 4, same Δ | 1000 | +0.510 ± 0.001 | +0.511 | +0.510 | +0.509 | +0.525 | +0.516 | +0.555 → +0.536 | +0.746 → +0.725 |

   errors = size of the a⁴ term; omitted a⁵ ≈ 1/13 of it.
   * **The correction is helicity-odd and universal:** `Γ_K^{h}(ν) ≈ σ_h U(ν)` with `σ_{++} = +1, σ_{+−} = −1` and `U` the **same function**
     as the `++` J = 0 image of BELOW_FOLD2 (per harmonic −1.0 / −4.6 / −8.3 / −10.4 / −9.9 at ν = 0.5 / 1 / 1.5 / 2 / 2.5 for J = 0,
     J = 4 `++`, and minus J = 4 `+−`, agreeing to ≤ 0.6). The helicity-even part projects to `(s_d^{++} + s_d^{+−})/2 = 0.0001` (H0).
     So **`+−` gains** (`s_d < 0` raises FF), `++` loses, by the same amount in every gauge and at both Δ_g [NE].
   * Helicity mixing: `m²κ ≈ +1.3–1.4` (ν-flat), `|Γ_off|/Σ|R_j| ≤ 1.4` on ν < 18, `s_off = −0.003 … −0.009` (300), `+0.001 … −0.005`
     (1000): negligible in the 2×2 condition (`FF_off/N ≲ 10⁻³` against diagonal FF/N ≥ 0.1) [NE + M].
2. **Full-set minimum and δ verdict [NE + M].** With measured `s_d` on `++` J = 0, 2, 4 and `+−` J = 4 (others at K = 0):
   **δ = 0.05: +0.0061 ± 0.0013 (SB +0.0048); δ = 0.2: +0.1427 ± 0.0017 (SB +0.1410)**, both Δ_g, binding state unchanged (`++` J = 0 at the
   gap). Applying the universal `σ_h × 0.51` to every state (odd-J `+−` bracketed with the worst sign) changes nothing: odd-J `+−`
   worst +0.036 / +0.183. **δ = 0.2: closed on every measured binding state with margin ≥ +0.14 [NE + M]**; the only unmeasured near-binding
   channel, `+−` J = 5 (parity-odd: the other agent's sector, **not attempted**), fails δ = 0.2 only if `s_d ≳ +5` (10× any measured
   value). δ = 0.05 stays positive but at +0.005 ± 0.001 (10 % of δ) [M].
3. **The 512 c₂ is a comparison artefact; K₂ is right [DV].** `bg_k2check` built the Laplace side as `log(Σ M_k)`; the rows' `M2/M0` is a
   Δ-constant (+0.052) and `M4/M0` a ν-only term (−0.076 at ν = 20, the same in the reference), so the plain sum carries a spurious
   `−(A2 − R2)A4 ≈ 600 c₂` cross term. In the form Result D actually uses (`l9_analysis.pred`, K = 2), K₂ = engine to **2.0 c₂** (J = 0,
   Δ = 400, ν 20 → 40; engine a⁴ step 2.1), **15 c₂** at Δ = 1600, ν 80 → 160 (engine a⁴ step 16), and −22 c₂ at J = 2 (engine a⁴ step
   107); the scalar analogue (−107 c₂ quoted in BELOW_FOLD2) is likewise +2.5 c₂ vs the exact skeleton. **Impact on Result D: none** (§2.2).
4. **Spin-J contraction [DV] / engine integration [I].** `⟨(A·z)^J,(B·z)^J⟩_κ = N_J (A²B²)^{J/2} P_J(Â·B̂)`, `N_J = J!/(2J−1)!!`, and its
   polarised form (derivatives `P_J^{(n)}`, n ≤ 8) reproduce the 60-digit sampled pairing to 4·10⁻¹² for J ≤ 32 (§3). Not yet wired into
   the ASer engine, so the Result D point (336, 28), ν < 8.3, Δ_g = 300 was **not computed [I]**. `HYPOTHESES_CLOSE2.md` appeared at 12:16 (the `oo_*` agent,
   whose scope includes exactly this window, by direct thimble/n3 nodes) but is a header-only stub (`@@SUMMARY@@`) at 13:45: it settles nothing yet.

---

## 1. `+−` J = 4 (`bh_eeq4.py`, `bh_helbasis.py`, `bh_qan4.py`, `bh_ffscan.py`)

### 1.1 Images [DV for the algebra]
`bg_eeq.images` (read-only) at Δ = 304.0297506 (Δ_g = 300 gap, `m = 1.01329Λ`) and 1004.0089775 (Δ_g = 1000, `1.00400Λ`), K = 4,
depth 5, kind EE4 (`k = 4`: EE2 × Lemma-B R), **all four (a, c) entries** of `pref[a][c]` (a: `⟨T₁T₂O⟩` conserved structures at Δ, c:
`⟨T₃T₄Õ⟩` at 3 − Δ); reference Δ = 2·10⁴, K = 2, depth 3. Every entry has the same leading a-power (`d = −4`), so the entries combine
linearly with no dropped factor. Cost: 77 / 73 min (pieces), 0.9–1.1 GB; reference 9 min.
`M_{ab}(ν) = Σ_c I_{ac}(ν)(b⁻¹)_{bc}` with the shadow matrix `b = l6_spin_action.shadow(Δ, 4)` (fit residuals `1.3·10⁻⁶`; constant in ν, so
it only fixes the structure mixing), exactly as `l6_spin_action.action` forms `M = Ser·b⁻ᵀ`.

### 1.2 Helicity basis [DV]
`T0(Δ)` = l8_analysis.basis rule (rank-1 decomposition over x of the order-0 crystal matrix, `x = 0.15 … 0.8`, nine fresh l6 Laplace
points per Δ — above the fold; `++` = flattest `d/W_h`), computed at the **same** Δ as the image (the conserved-structure basis is
SVD-defined per Δ, so bases are not transportable): rank-2 residual `sv₃/sv₁ = 1.7·10⁻⁶ / 1.5·10⁻⁷ / 4·10⁻¹⁰`, order-0 off-diagonal
`1.5·10⁻³ / 1.5·10⁻⁴ / 4·10⁻⁷` (Δ = 304 / 1004 / 2·10⁴), flatness 0.034 vs 2.33. Checks on the engine side:
* reference image in `T0(2·10⁴)`: off-diagonal `3.2·10⁻⁹` and `I_{+−}/I_{++} = 0.732625` constant to `10⁻⁶` at every ν (a⁰ helicity
  blind), with the a² drift at ν = 40 equal to `20 x²` (`W_{+−}/W_{++} = 1 + 20x²`, J = 4);
* state: `Re Γ_exact/Re Σ_j R_j(ν − ij)² ≈ 19–24` (`+−`) and 0.7–0.9 (`++`) for ν ≥ 3, i.e. HYB's `γ_{+−} = 20.75`, `γ_{++} = 0.75`: the
  channels are identified correctly.
A basis change can only move diagonal elements by a constant (absorbed by `c_h`) or at second order in the (≤ 1.5·10⁻³) mixing.

### 1.3 Chain [D]
As bg_qan per helicity: `Γ_exact,h = m²(c_h T_ref I_h,st/I_h,ref − T_st) + γ_h (m/m_ref)² Σ_j R_j(ν − ij)²` (the last term restores the
reference's own finite-J term; it moves `s_d` by 0.005 and is dropped in bg_qan), `Γ_K,h = Γ_exact,h − γ_h Σ_j R_j(ν − ij)² − const·N̂₄`,
`Γ_K(0.02) = 0`; projections on the whole grid (ν < 41; ν < 18: H0 −0.519, SB −0.477 — SB's 1.3 % weight above 18 meets the ν²-growing
x⁴ remainder there). FF/N = Stage-6 HYB at the state, `FF → FF − a·2u³ s_d − c(2u³/Λ²)s_v` (`s_v = +13.6 ± 0.4`).

### 1.4 Per-harmonic profile and order structure (Δ_g = 300, `+−`, `Γ_K/Σ|R_j|`, [a⁴ term])
ν = 0.5 / 1.02 / 1.5 / 2.02 / 2.5 / 3.02 / 4.02 / 6.02 / 10 / 16: +1.03 [−0.20] / +4.56 [−0.36] / +8.34 / +10.61 / +10.23 / +7.77 / +1.66 /
−2.23 / +1.41 / −2.35. This exceeds Stage 7's uniform tolerance `G_U ≈ 5–6.5 c₂` per harmonic at ν ≈ 1.5–3, **with the favourable sign**.
Per-order `s_d` (H0): +3.83 / −9.91 / +5.02 / +0.51 / +0.039 (k = 0…4; Δ_g = 1000: +3.83 / −10.52 / +6.00 / +0.18 / +0.004). Covariant
part (bg_covfit rule, ν < 12): `δγ = −0.047` (−0.2 % of 20.75); the rest (−0.470) is the universal non-covariant `−U`.

### 1.5 Full-set scan (`bh_ffscan.py` → `cache/bh_ffscan.txt`) [NE + M]
| δ | Δ_g | min FF/N (H0 / SB / XP) | at | `+−` J = 4 channel | `++` J = 2 | `++` J = 4 |
|---|---|---|---|---|---|---|
| 0.05 | 300 | +0.0061 ± 0.0013 / +0.0048 / +0.0056 | `++` J = 0, 1.000Λ | +0.117 (gap) | +0.020 (1.089Λ) | +0.263 (1.57Λ) |
| 0.05 | 1000 | +0.0061 / +0.0048 / +0.0056 | `++` J = 0 | +0.099 (gap) | +0.020 | +0.263 |
| 0.2 | 300 | +0.1427 ± 0.0017 / +0.1410 / +0.1421 | `++` J = 0 | +0.277 | +0.167 (1.118Λ) | +0.380 (1.63Λ) |
| 0.2 | 1000 | +0.1427 / +0.1410 / +0.1421 | `++` J = 0 | +0.259 | +0.167 | +0.380 |
`s_d` is applied m-independently within each (h, J) [M] (as bg_ffscan); unmeasured channels stay at K = 0 (Stage 6 minima ≥ +0.0497 /
+0.1970). Universal extrapolation (`σ_h × s_0` on all states, odd-J `+−` worst sign): same minima; odd-J `+−` worst +0.036 (δ = 0.05),
+0.183 (δ = 0.2).

### 1.6 Reading
* **δ = 0.2 holds on every binding state that the parity-even engine can reach** (`++` J = 0, 2, 4; `+−` J = 4), every gauge, both
  Δ_g, margin ≥ +0.141 [NE + M]. `+−` J = 5 would need `s_d ≳ +5` (δ = 0.2) or `≳ +1.8` (δ = 0.05) to fail.
* δ = 0.05: not broken by `+−` J = 4 (it improves to +0.10–0.12); still +0.005 ± 0.001 at `++` J = 0 — positive but thin [M].
* The `σ_h`-odd, J- and Δ-independent shape of `Γ_K` suggests a **structural** O(1/m²) term of the exact k = 4 EE image (of the same
  helicity parity as the `B₅ = σu³P₀` piece) rather than a state-specific correction; whether HYB already carries it elsewhere
  (double counting) is **[I]** and would change both the `++` loss and the `+−` gain by the same amount.

## 2. The 512 c₂ "K₂ vs engine" disagreement: a resummation-form artefact of the comparison, not of K₂

### 2.1 Diagnosis [DV]
BELOW_FOLD2 §3.3 compared `m²[D(ν) − D(19.92)]`, `D = log(state/ref)`, with the Laplace side formed as the **plain sum**
`log(M0 + M2 + M4 + M6) − (ref)` (`bg_k2check.py`; `M6 = 0` in these rows: `L1_DEG = 6`). The rows' relative orders are not small in
the way that form needs: `A2 = M2/M0 = +0.0524` is a Δ-dependent constant (≈ 8 200/m²), while `A4 = M4/M0 = −0.0757 / −0.0175` at
ν = 19.9 / 39.8 is dominated by a **ν-only** piece that is the same in the reference (`R4 = −0.0776 / −0.0194`). The plain sum contains
the formally `O(ε⁶)` cross term `−A2·A4 + R2·R4 ≈ −(A2 − R2) A4`, i.e. `m² × 0.0515 × 0.0757 ≈ 620 c₂` at ν = 20 and ≈ 140 at ν = 40:
the 512. Result D never uses this form: its rows go through `l9_analysis.pred(·, K = 2)` (state series × `exp(−κ_ref)`, truncated at
`ε⁴`; `q3_analysis` → `A.grid(J, 2, …)`), equivalently the cumulant form of l8.

`bh_k2forms.py` (→ `cache/bh_k2forms.txt`), `m²[D(ν) − D(ν₀)]`, Laplace minus trusted value:

| check | ν | trusted | plain (bg_k2check) | cumulant | **chain = Result D's `pred`, K = 2** | K₁ only |
|---|---|---|---|---|---|---|
| scalar Δ = 400, vs exact Mellin, ν₀ = 39.85 | 19.92 | −1800.3 | **+106.9** | +2.5 | **+2.5** | +2.1 |
| | 30.0 | −1041.1 | +27.0 | +0.6 | +0.6 | +0.4 |
| | 60.0 | +3073.8 | −19.5 | −0.4 | −0.3 | +0.3 |
| | 79.7 | +7347.2 | −19.4 | +6.1 | +6.3 | +7.8 |
| EE k = 2, J = 0, Δ = 400, vs engine K = 4 (K = 3), ν₀ = 39.85 | 19.92 | +598.9 (+596.8) | **+512.4** | −2.0 (+0.1) | **−2.0 (+0.1)** | −2.0 |
| | 79.7 | −2425 (−2406) | −105 | +14 (−5) | +14 (−5) | +15 |
| | 119.6 | −6590 (−6499) | −63 | +78 (−14) | +78 (−14) | +80 |
| EE k = 2, J = 0, Δ = 1600, vs engine K = 4 (K = 3), ν₀ = 159.85 | 79.92 | +9565.7 (+9549.6) | +96.8 | −14.9 (+1.2) | **−14.9 (+1.2)** | −15.1 |
| EE k = 2, J = 2, Δ = 400, vs engine K = 4 (K = 3), ν₀ = 39.85 | 19.92 | +7942 (+8049) | **+503** | −22 (−129) | **−22 (−129)** | −20 |
| | 79.7 | −34406 (−36966) | +43 | +166 (+2725) | +166 (+2725) | +157 |

* **Scalar:** the same "−107 c₂" BELOW_FOLD2 quoted for the scalar Laplace is the plain-form artefact; in Result D's form K₂ misses
  the exact skeleton by ≤ 2.5 c₂ on ν = 20–60 (6 c₂ at 80) [DV]. (The exact skeleton record itself fails at ν ≥ 100 at Δ = 400: all
  three methods disagree with it identically there; not used.)
* **EE J = 0:** K₂ (chain form) = engine to **2.0 c₂** on ν = 20 → 40, i.e. within the engine's own a⁴ step (2.1) [DV]. At ν ≥ 80 the
  engine is the unconverged side (its a⁴ step is 19 / 91 c₂ at ν = 80 / 120; x = 0.2–0.3 is outside the fixed-ν expansion's comfort).
* **EE J = 0, Δ = 1600** (x = 0.05 → 0.1, the regime of Result D's Δ_g = 1000 rows): chain −14.9 c₂ against an engine a⁴ step of
  16.1 (K = 3: +1.2) — agreement within the engine's truncation [NE]. (ν = 320 overflows the engine.)
* **EE J = 2:** chain −22 c₂ against an engine whose a⁴ step is itself 107 c₂ on 20 → 40 (the J(J+1) terms slow the fixed-order series);
  agreement within the engine's truncation; the plain form is again off by ≈ 500 [NE].
* **Verdict: K₂, as Result D uses it, is right at the tested points (J = 0, 2; Δ = 400, ν = 20–40; J = 0, Δ = 1600, ν = 80–160) to the engine's accuracy
  (2 c₂ at J = 0); the 512 c₂ was produced by the comparison script's resummation, not by K₂** [DV for J = 0, NE for J = 2].

### 2.2 Impact on Result D [D]
None from this item: Result D's near-gap and high-J tables use the chain form, which agrees with the engine. For scale only: a
uniform 500 c₂ on ν = 20–40 would have been sub-threshold for the computed states (failure thresholds for uniform windows are
`≳ 4·10³ c₂`: HYPOTHESES_CLOSE §4–5; edge `≈ 5.6·10³`), but well above the `H_LJ` tolerances (1–30 c₂, uniform) — that worry is removed.
**Not tested here:** K₂ at large `J²/m` (HYPOTHESES_CLOSE §5 already shows it fails by ×5 at (336, 28) from x = 0.2) and K₂ at J ≥ 4
against the engine (the engine needs the closed-form spin-J contraction, §3, to reach the Result D spins). The disagreement BELOW_FOLD2
attributed to "Laplace truncation" is withdrawn.

## 3. Closed-form spin-J contraction for the EE polarisation (`bh_spinJ.py`)

### 3.1 The identity [DV]
The O-index contraction of the engine (and of l5_spin / w0) is the null-circle pairing
`⟨f, g⟩ = Σ_{n=−J}^{J} f_n g_{−n}/κ_n` of two degree-J polynomials in the null polarisation z₅ (Fourier coefficients on
`null_circle(J)`, `κ_n = (−1)^{J+n} 2^{−J} C(2J, J+n)`). For linear forms:

  `⟨(A·z)^J, (B·z)^J⟩ = N_J (A²B²)^{J/2} P_J(A·B/√(A²B²))`,  `N_J = J!/(2J−1)!!`

(checked for J = 1 … 12 to 10⁻¹⁵, complex vectors, the light-cone metric `−½(p₀q₁ + p₁q₀) + p₂q₂`). By polarisation, for the
structures actually present — `f = (V·z)^{J−p} Π_{i≤p}(A_i·z)` with `V = V_{5,12}`, `A_i ∈ {H₁₅, H₂₅}` (powers ≤ 2 each, so `p ≤ 4`),
and the same on the `⟨Õ T₃T₄⟩` side —

  `⟨f, g⟩ = N_J (J−p)!(J−q)!/(J!)² · Π_i(A_i·∂_a) Π_j(B_j·∂_b) H(a, b)|_{a=V, b=W}`,  `H = (a²b²)^{J/2} P_J(a·b/|a||b|)`,

a finite sum of `P_J^{(n)}(t)`, `n ≤ p + q ≤ 8`, `t = V·W/√(V²W²)` (the scalar spin angle), times polynomials in the dot products of
`{V, A_i}` with `{W, B_j}` and powers `(V²W²)^{J/2 − k}`. Checked against the **60-digit** sampled pairing for J = 4 … 32,
`(p, q)` up to `(4, 4)`, near-forward and random kinematics: **max relative deviation 4·10⁻¹²** [DV] (`cache/bh_spinJ.txt`,
exact multilinear nilpotent evaluation of the derivatives, `P_J^{(n)} = (2n−1)!! C^{(n+½)}_{J−n}`).

Side finding [NE, random kinematics only]: the **double-precision** sampled pairing (the engine's and l5_spin's method) loses all
digits for generic complex vectors at J ≳ 24 when `|⟨f,g⟩| ≪ |V|^J|W|^J` (relative errors 10⁻⁷ … 10⁷ at J = 24–32 in the table);
the project's large-J records carry `fw_frame.patch_S` for exactly this floor, so this is a confirmation, not a new defect, but the
closed form removes the issue structurally.

### 3.2 What it buys and what is left [D / I]
* In the engine, `V₅^J` powers on 2J+1 samples (J = 4 takes 75 min at depth 5 vs 12 min for J = 2; the sampled V₅^J products grow steeply with J) are replaced by ≤ 25 terms
  per structure pair, each `P_J^{(n)}(1 + δ)` times ASer dot products: **J-independent cost** apart from the `P_J` series.
* `P_J^{(n)}(1 + δ) = Σ_m c_{m+n}(J) (m+n)!/m! δ^m`, `δ = t − 1 = a²τ₂ + a⁴τ₄ + …` — the same object bg_scalarJ already resums at fixed
  J/Δ (N = 12, R = 4 suffices for the scalar at (336, 28/32), §2.2 of BELOW_FOLD2); the EE version needs the same resummation for
  `n ≤ 8` derivatives and for `(V²W²)^{J/2}` against the 3-pt exponents `−J/2 log m₁₅ …` (which partly cancel it).
* **Not implemented in the engine [I].** The ASer implementation (vectors `V, H₁₅, H₂₅` as 3-vectors of ASer, already built in
  `bg_ee.tto` as `V5b, H15b, H25b`), its J = 4 check against the sampled `pref`, and the J/Δ resummation are the next step;
  estimated 1–2 days of work. Consequently **the Result D test at (336, 28), ν < 8.3, Δ_g = 300 was not computed [I]**.
  `HYPOTHESES_CLOSE2.md` did not exist at 12:00; the `oo_*` agent created it at 12:16 with that window (J = 28, 32, direct nodes) in
  its title, but at 13:45 it is a header-only stub — not settled there yet either; check it before computing.

## 4. Remaining open items [I]
1. `+−` J = 5 (and odd J generally): parity-odd `⟨TTO⟩` structure — the other agent's sector; not attempted. `+−` J = 6, 8: not run
   (≈ 3 h each at depth 5 with sampled contraction); the universal-`σ_h` extrapolation says favourable [M].
2. Origin of the `σ_h`-odd universal `U(ν)/m²` and whether HYB double counts it (§1.6). If it is structural, both δ verdicts move by the
   same ±0.04 (δ = 0.05) / ±0.057 (δ = 0.2) at `++` J = 0 — the δ = 0.05 sign hinges on it, δ = 0.2 does not.
3. EE beyond LO on the real section ≡ crystal cycle (BELOW_FOLD2 §5.4) — still assumed; the n3 EE J = 4 points at Δ = 100
   (`y4_n3_pts_EE_J4`, ν = 9.8, 14.8) carry 0.2–1 % errors, i.e. ±20–90 c₂: too coarse to test it; not used.
4. Engine at J ≥ 6: the closed-form contraction (§3) in ASer plus the J/Δ resummation of `P_J^{(n)}`; then K₂ vs engine at the Result D
   spins, and the (336, 28) ν < 8.3 window.
5. Inherited: HYB elsewhere [M], K at unmeasured states, `b₀, b₅, c₁`, C5 beyond t ≈ 18, etc. (Stage 6/7 lists).

## 5. Graded claims
| claim | grade |
|---|---|
| a⁰ EE image is helicity-blind: in `T0(Δ)` the reference image is diagonal to 3·10⁻⁹ with ν-independent `+−/++` ratio (10⁻⁶) | [DV] |
| `T0(Δ)` (order-0 Laplace) diagonalises the engine's J = 4 images; channels identified by `γ_h` (Γ_exact/cov ≈ 20 vs 0.8) | [DV] |
| `+−` J = 4 gap: `s_d = −0.516 ± 0.039` (300), `−0.510 ± 0.004` (1000); gauge spread ≤ 0.017 | [NE] (given D ≡ crystal for EE beyond LO) |
| `++` J = 4 at the same Δ: `s_d = +0.517 ± 0.017 / +0.510 ± 0.001` | [NE] |
| `Γ_K^h ≈ σ_h U(ν)`, U universal (J = 0, 2, 4; Δ_g = 300, 1000); helicity-even projection ≈ 0 | [NE] |
| `+−` J = 4 FF/N: +0.118 / +0.099 (δ = 0.05), +0.278 / +0.259 (δ = 0.2) (Δ_g = 300 / 1000) | [NE + M] |
| full-set min with measured states: +0.0061 ± 0.0013 (δ = 0.05; SB +0.0048), +0.1427 ± 0.0017 (δ = 0.2; SB +0.1410) | [NE + M] |
| δ = 0.2 holds on all parity-even binding states; `+−` J = 5 unmeasured (fails only for `s_d ≳ +5`) | [NE + M] / [I] |
| helicity mixing negligible in the 2×2 condition | [NE + M] |
| the 512 c₂ was the plain-sum resummation of the comparison; K₂ in Result D's `pred` form = engine to 2 c₂ (J = 0, Δ = 400) | [DV] |
| same at Δ = 1600 (15 c₂ vs engine step 16) and J = 2 (−22 vs engine step 107); scalar K₂ vs exact ≤ 2.5 c₂ (ν 20–60) | [NE] / [DV] |
| no impact of this item on Result D's minima | [D] |
| closed-form spin-J pairing (Legendre identity + polarisation) = sampled pairing, J ≤ 32, to 4·10⁻¹² | [DV] |
| double-precision sampled pairing loses all digits at J ≳ 24 for generic kinematics | [NE] |
| engine with closed-form contraction; Result D (336, 28) ν < 8.3 point | [I] |

## 6. Files (`ads4/code/spinfunc/`; caches in `cache/`)
| file | content | run |
|---|---|---|
| `bh_eeq4.py` | J ≥ 4 EE k = 4 images, all (a, c) entries (wraps `bg_eeq.images`) → `cache/bh_eeq_D{Δ}_J4_K{K}.pkl` | 9 min (ref) / 75 min (state) |
| `bh_helbasis.py` | `T0(Δ)` from fresh l6 Laplace points → `cache/bh_helbasis.pkl` | 1 min per Δ |
| `bh_qan4.py` | helicity-resolved Result Q chain, covariant split, mixing → `cache/bh_qan4.txt`, `bh_qan4_D*_J4_Dg*.pkl`, `bh_shadow.pkl` | 1 min |
| `bh_ffscan.py` | full-set min FF/N with measured `s_d` (`++` J = 0, 2, 4; `+−` J = 4) + universal-`σ_h` extrapolation → `cache/bh_ffscan.txt` | 2 min |
| `bh_k2scal.py` | scalar Δ = 400: exact vs engine vs Laplace K0–K2 → `cache/bh_k2scal.txt` | 5 min |
| `bh_k2forms.py` | plain / cumulant / chain (`pred`) Laplace forms vs exact (scalar) and engine (EE J = 0, 2; Δ = 400, 1600) → `cache/bh_k2forms.txt` | 3–15 min |
| `bh_spinJ.py` | closed-form spin-J contraction and its check (mpmath) → `cache/bh_spinJ.txt` | 2 min |
Order: `bh_eeq4.py 20000 4 2 3`; `bh_eeq4.py 304.0297506031432 4 4 5`; `bh_eeq4.py 1004.0089775159122 4 4 5`;
`L1_DEG=6 bh_helbasis.py 4 304.0297506031432 1004.0089775159122 20000`; `bh_qan4.py <Δ> 4 4 <Δ_g>`; `bh_ffscan.py`;
`bh_k2scal.py`; `bh_k2forms.py scalar`; `bh_k2forms.py ee 0 [1600 330]`; `bh_k2forms.py ee 2`; `bh_spinJ.py`.
