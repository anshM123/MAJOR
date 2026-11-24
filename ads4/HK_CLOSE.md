# H_K: the heavy-sum ↔ (corner, ν) interchange under the KKSDZ bound `J₀ ≤ 1`

*23 Nov 2026. Unreviewed, not committed. This note answers REFEREE_REPORT §M8. No existing file was edited.*
* *New code:* `ads4/code/hk/hk_toy.py`, `ads4/code/hk/hk_real.py`.
* *Outputs:* `ads4/code/hk/cache/` (git-ignored): `hk_toy.txt`, `hk_real_q_d_flat.txt`, `hk_real_qcut.txt`.
* *Read-only inputs:*
  * `code/resultq/functionals_q7_xprec.json` (XP C4 weights, 64 digits);
  * `code/resultq/functionals_q4_deshift.json` (the `R_j`);
  * `code/spinfunc/cache/u3_chi_EE_{300,1000}.pkl` (Result D's Abel weights).
* *Compute:* about 7 CPU-minutes, at most 2 processes, under 1 GB.

Grades:
* **[P]** proved;
* **[D]** derived;
* **[D-sketch]** derivation sketched;
* **[DV]** derived and checked numerically;
* **[NE]** numerical evidence;
* **[M]** model;
* **[A]** assumption;
* **[Q]** quoted;
* **[I]** open.

A chain carries the grade of its weakest link.

---

> **Update (REGGE_TREE.md, unreviewed): referee-3 M2 is a presentation defect (the order of limits was never stated), not a derivation gap. Done as CHMRS do it, the Result never needs H_K at tree level; all grades stay.**
> * **H_K's J₀** is the exponent of the uniform bound on dDisc of the **exact** correlator (KKSDZ (4.70)). **Route N [D]:**
>   1. Apply the sum rule at finite c_T. H_K is [D] under J₀ ≤ 1, with margin ½ at k = 2, 3, 5 and 3/2 at k = 4.
>   2. (Pos) makes the heavy side ≥ 0 at every c_T, multi-trace states included.
>   3. Expand only the light side in 1/c_T.
>
>   The tree-level heavy sum is never expanded; CHMRS §4.4 show it need not be analytic in 1/c_T. The order-of-limits step is [D-sketch] given A2. A2 must be restated as tree-level EFT light data up to o(1/c_T), uniformly in spin.
> * **Term-by-term reading (not made):** with a uniform planar J₀ = 2, the corner count fails at k = 2, 3 **and 5**, and for generic k = 4 weights (Q's k = 4 weights pass with margin ½). "J₀ = 2 (Einstein)" is the Wightman/elastic intercept; the heavy tree dDisc is empty for pure Einstein, and for strings its fixed-η intercept is 2 − 9/(2Λ²) with an η⁻² k = 2 moment (CHMRS (4.45)).
> * **Literature [Q]:**
>   * CHMRS (§1.1.1, (1.12)–(1.14), §4.1, §4.4) and CLS-D (§3.4, §3.5.2) use J₀ ≤ 1 non-perturbatively and assume no tree-level J₀ < 2.
>   * KKSDZ §4.5 gives J₀^planar ≤ 2 (4.80); their planar sum rules need J₀ < 3.
> * **Model check:** with the uniform J₀ = 1 bound, the corner integrals converge for every k; with the uniform J₀ = 2 bound they diverge at k = 2, 3, 5; with a stringy Regge-pole model they converge for every k.
> * **Wording:** proposed replacements for letter.tex:109 (A2), :133–135 and :183–185, plus supplement paragraphs "Order of limits" and "Term-wise use (not made)", are in REGGE_TREE §5. Say "every state with Δ−J ≥ Δ_g", not "heavy single traces".
> * **Open:** the size of the light-side remainder O(Δ_g^p/c_T²) is not computed, so keep "at leading order in 1/c_T".

> **Update (NU_TAIL_CHECK.md, unreviewed): the exponential growth above a state's mass affects NO computed number, but the published grid-weight DEFINITION at Δ_g = 1000 is defective, and the fix is free.**
> * **No effect on computed numbers [P, by inspection]:** every computed per-state number (crystal/D′, K₂ rows, direct n1/n3, GK edge nodes, windows, EDGE-R, hybrid floors) uses ν ≤ Λ ≤ m. The light-side constants are unaffected, and so is everything built from the continuum weights or the exact min-norm X: X and X^Q decay like e^{−(π/2)(ν−Λ)} beyond Λ, and tails beyond m are ≤ 4.4·10⁻⁷ (D) and ≤ 1.6·10⁻⁴ (Q) of a state's action.
> * **Problem (c):** the stored grids `prl/data/d_abel_X_1000.txt` and `prl/data/k2/xq_abel_X_1000.txt` run to Λ+312 and carry a non-decaying double-precision floor (|X| ≈ 10⁻¹⁵–10⁻¹⁶) beyond Λ. Read literally with θ(Λ+100, 10), that floor times the exact image growth (Γ ≈ 10¹⁶ at Λ+100) gives tails up to 650× (D) and 1300× (Q) a near-gap state's action. The sign is set by round-off, so (Pos) for that grid-defined object is untested at m − Λ ≲ 60–100. HK_CLOSE §5 and §6 claims are false for that object, though true for the computed values. Δ_g = 300 grids are unaffected (they end at Λ+29).
> * **Fix (changes no number):** define the functional by the continuum weights; or read the grids as zero beyond Λ+10; or use the cutoff (Λ+60, 6) for the k = 2 pieces (same leakage, tail 5·10⁻⁹ / 4.6·10⁻⁷).
> * **Growth numbers [DV]:** exact scalar Arb series. Γ(m) = 1.8·10³ and Γ(Λ+100) = 1.6·10¹⁶ at Δ_g = 1000, J = 0; the rate per unit ν is 0.13–0.19 at x = 1, rising to π/2 only as x → ∞.

## 0. Verdict

1. **As currently written, H_K cannot be obtained from `J₀ ≤ 1` by any absolute-convergence argument. The obstruction is the functional's definition, not Regge physics.**
   * **Compact ν-support [P].** Take any de-shift weight supported on a compact ν-set; every stored grid weight is one. By Paley–Wiener, its
     leading-Regge heavy image grows at large `β` at the **raw** rate:
     * `e^{β/2}` for `Ĉ′_EE`;
     * `e^{5β/2}` for `Ĉ₄^{ds}` and `∂Ĉ₄^{ds}`.

     Tonelli under the uniform KKSDZ bound then needs `J₀ < ½` (§3). This is the old `H_ip`, and it cannot be improved for such weights.
   * **The stored weights fit this pattern [NE].**
     * Q's XP weights satisfy the de-shift identity up to a residual `1.7·10⁻⁵⁴ e^{3t}`. Their image equals the ideal `P_{−½}` only
       for `β ≲ 42`.
     * D's double-precision Abel weights do so only for `β ≲ 25`.
   * **The uncut exact de-shift weight is no better.**
     * It is not a position-space kernel: its weight grows like `e^{κν}`, `κ → π/2⁺` (V4FINAL §2.4).
     * In the flat model, its tail beyond the light-cone edge `x = ν/m > 1` makes per-state images diverge for near-maximal-spin states
       with `m ≳ 2.2–2.7Λ` (§6) [M].
2. **The obstruction is removable at no numerical cost [D, conditional on the inputs below].**
   * **Definition.** Define each de-shifted piece with an **analytic cutoff**: weight `χ_c = χ·c`, where
     `c(ν) = ½[erf((ν+ν_c)/a) − erf((ν−ν_c)/a)]` is even, entire, bounded on strips and has super-Gaussian decay.
   * **Proposition K_c (§4).**
     * The cutoff functional `Ω_c` is a genuine position-space functional.
     * Its leading-Regge image differs from the ideal one by `δ_c`. The bound on `δ_c` is:

       | piece | bound on `δ_c` | grade |
       |---|---|---|
       | `Ĉ′_EE` | `ε_c·P_{−½}(cosh β)` | [D] |
       | `Ĉ₄^{ds}`, `∂Ĉ₄^{ds}` | `ε_c·C·e^{(½+0)β}` | [D] |
       | `Ĉ₄^{ds}`, `∂Ĉ₄^{ds}`, given the measured analytic structure | `ε_c(P_{−½} + |D|)` | [D + NE] |

     * Hence, by Tonelli, the heavy sum converges absolutely and H_K holds for `Ω_c` whenever:
       * `J₀ < 3/2` for `k = 2`;
       * `J₀ < 5/2` for `k = 4` (`7/2` given the NE structure).

       In particular it holds under KKSDZ's `J₀ ≤ 1`, with margin **½**.
   * **Size of `ε_c` [NE/D].**
     * Q (`ν_c = 150`, `a = 6`): `ε_c ≲ e^{−180}`.
     * D (`ν_c = Λ + 100`, `a = 10`): `ε_c ≲ 10⁻⁴¹`.

     Every computed number is unchanged at that level, because the stored weights agree with `χ_c` to that accuracy on their grid.
   * **Two pieces need nothing new.**
     * `Ĉ₅^{ds}`: its raw per-ν kernel already passes (`J₀ < 3/2`), so H_K is plain Fubini.
     * `Ĉ₃`: it is not de-shifted.
3. **Status of H_K.** It moves from **[D-sketch] ("[D given H_ip, J₀ < ½]")** to
   **[D given K1′, K2ν, LRC, K3, DS]**, i.e. derived from inputs that are already in the hypothesis list.
   * K2ν is the only new item: a mild ν-uniformity of the K2 kernel bounds, [D-sketch].
   * The requirement `J₀ < ½` is gone. No independent Regge hypothesis remains.
   * Net grade: **[D-sketch]**, by the weakest link (K1′ and K2ν are [D-sketch]). It is no longer an open problem of a different kind.
4. **The theorem's functionals must be redefined with the cutoff.**
   * This is a definitional change. No computed number moves.
   * Without the change, H_K for the stored functionals is false-or-unprovable under `J₀ ≤ 1` (item 1).
   * Without the change, the literal "exact de-shift" is ill-defined on part of the heavy spectrum (§6, [M]).

---

## 1. Precise statement of H_K

**Notation.**
* **Per-ν images.** A de-shifted piece is `Ĉ = ∫₀^∞ χ(ν) Ψ̂_ν dν`, where `Ψ̂_ν = Ψ_{k,ν}/N_com(ν)` is CLS-D's per-ν light-ray functional
  in heavy normalisation.
* **Leading-Regge (LR) image of a heavy state `O`** (`cosh β = 1 + (2J+1)²/2m²`):
  * `a_O(ν) = Σ_j R_j(ν)P_{−½+iν+j}(cosh β)`;
  * the shifts are `j = ±1` for `k = 2` and `j = ±1, ±3` for `k = 4`;
  * `R_j` are as in RESULT_D_STATEMENT §2 (`h_∓`) and RESULT_Q_STAGE4 §1.3.
* **Mehler–Dirichlet form.**
  * `I_χ(β) = ∫χ a_O dν = (√2/π)∫₀^β W_χ(t)(cosh β − cosh t)^{−½}dt`, with `W_χ(t) = Σ_j∫₀^∞χR_j cos((ν−ij)t)dν`.
  * The de-shift is the identity `W_χ = Φ` on `t ≥ 0`, where `Φ = C[ψ]` is the target. The targets are:

    | piece | `Φ` | image |
    |---|---|---|
    | `Ĉ′_EE` | `C[ψ₂]` | `𝒦_{ψ₂}` |
    | `Ĉ₄^{ds}` | `1` | `P_{−½}` |
    | `∂Ĉ₄^{ds}` | `−t²/2` | `−D(β)` |

**H_K (precise).** For each `Ĉ ∈ {Ĉ′_EE (D's and Q's ψ₂), Ĉ₄^{ds}, ∂Ĉ₄^{ds}, Ĉ₅^{ds}}` and every theory satisfying A1 and A2:

    Σ_O Tr[P_O · Ĉ[G_O]] = −Ĉ[light],    with  Σ_O |Tr[P_O · Ĉ[G_O]]| < ∞,

* `Ĉ[G_O] = ∫χ Ψ̂_ν[G_O] dν` is the per-state ν-outside image;
* `P_O` is the PSD OPE matrix;
* `Ĉ[light] = ∫χ Ψ̂_ν[light] dν`. Its value is 0 on the `J₀ < 4` terms (Lemma A + D_LR) plus (F)+(P) on the EFT terms.

**Where `J₀` enters: the corner count.**
* **The count.** For a position kernel `|K| ≲ r^a η^c` in the Regge corner `{ρ = re^β ≤ c}`, the uniform bound
  K1′ `|dDisc 𝒢| ≤ B r^{1−J₀}` gives absolute convergence iff `c < a + 1 − J₀` (Stage 3 §1.5; H7 §14). The Rindler/Cauchy–Schwarz
  bound (Lemma S, K3) then gives `Σ_O ∫|K||dDisc G_O| ≤ ∫|K| dDisc 𝒢 < ∞`.
* **Heavy-measure form.** Write `dμ_k` for the heavy measure (`Σ_O Tr P_O`·normalisation at `β_O`). Then K1′ ⇒ `∫e^{γβ}dμ_k < ∞`
  for every `γ < a + 1 − J₀`, and this is sharp when K1′ is saturated.
* **Thresholds:** `a = 0` for `k = 2, 3`; `a = 2` for `k = 4, 5`.

| piece | `c` (growth of the kernel the swap actually uses) | condition | under `J₀ ≤ 1` |
|---|---|---|---|
| raw `Ĉ₂` per ν (`P_{l+1}`) | ½ | `J₀ < ½` | fails |
| `Ĉ′_EE`, ideal (uncut) de-shift | −½ (`𝒦_w ~ η^{−½}/log η`) | `J₀ < 3/2` | ok, but not a position kernel |
| `Ĉ′_EE`, compact-support (stored) weight | ½ (§3) | `J₀ < ½` | **fails** |
| **`Ĉ′_EE`, analytic cutoff** | −½ (`ε_c P_{−½}`, §4) | **`J₀ < 3/2`** | **ok, margin ½** |
| `Ĉ₃` (no de-shift) | −½ | `J₀ < 3/2` | ok |
| raw `Ĉ₄`, `∂Ĉ₄` per ν | 5/2 | `J₀ < ½` | fails |
| `Ĉ₄^{ds}`, compact-support (stored) weight | 5/2 (§3) | `J₀ < ½` | **fails** |
| **`Ĉ₄^{ds}`, `∂Ĉ₄^{ds}`, analytic cutoff** | ½ guaranteed; −½ given NE (§4) | **`J₀ < 5/2`** (`7/2`) | **ok, margin 3/2** |
| raw `Ĉ₅` per ν | 3/2 | `J₀ < 3/2` | ok, margin ½ |

The old "`J₀ < ½`" (`H_ip`) is exactly the raw / compact-support rows.

---

## 2. What the four routes give

| route | outcome |
|---|---|
| **1. Proposition K for the heavy-normalised kernel** | **Works, after a redefinition.** The de-shift does improve the large-`β` behaviour (exact image `P_{−½}` / `𝒦_w`). It can only be used through a position-space kernel, and that requires the analytic cutoff (§4). Without the cutoff, no Proposition K exists, for two reasons: the uncut kernel is not a distribution on the CLS-D domain (V4FINAL), and every compactly supported approximation regrows at the raw rate (§3). |
| **2. Tauberian / positivity (Fatou, monotone convergence)** | **Not needed, and alone not sufficient.** For the inequality the bound needs (`Σ p_O I ≤ −L`), Fatou needs only a one-sided summable minorant of the approximants' images. The natural approximants (truncations) have images `±ε e^{5β/2}` that oscillate in sign, so no such minorant exists under `J₀ ≤ 1`. Positivity alone therefore does not rescue the compact-support functionals. |
| **3. KKSDZ / CLS-D style** | **This is the route used.** Per-block `dDisc` positivity (Rindler / Cauchy–Schwarz, K3) + the uniform Regge bound (K1′) + Tonelli in position space. Our cutoff kernels fall in its class (genuine kernels with corner growth below threshold); the uncut ones do not. KKSDZ (1904.05905) and CLS-D (2311.04271) supply the bound and the per-ν functionals; neither treats ν-smeared, de-shifted weights, and the swap for them is our Proposition K_c. |
| **4. Sharpest criterion if nothing works** | Superseded by route 1. For the **stored** (compact-support) functionals, the sharp criterion is the old `H_ip`: `∫e^{(j_max−½)β}dμ_k < ∞`, which holds if `J₀ < ½`, or if `J₀ ≤ 1` with transverse decay `|dDisc 𝒢| ≤ B r^{1−J₀}η^{−γ}`, `γ > J₀ − ½` [D]. That criterion is stronger than anything known non-perturbatively. |

---

## 3. The no-go for compactly supported weights [P]

**Lemma PW.**
* *Hypotheses.* Let `χ ∈ L¹` be supported in `[0, ν_c]` and not identically zero. Let `R_j` be rational, with `R_{j_max}` zero-free
  on the real axis (true for `EE`: no `R_j` has a real zero on `ν > 0`, Stage 4 §1.4).
* *Conclusion.* `limsup_{β→∞} β⁻¹ log|I_χ(β) − I_target(β)| = j_max − ½` for every target with `I_target = O(e^{−β/2}β^p)`.
  * For `k = 2`: the image regrows like `e^{β/2}`.
  * For `k = 4`: it regrows like `e^{5β/2}`.
  * The regrowth comes with at most a sub-exponential decay factor, set by the smoothness of `χ`.

*Proof.*
1. The `e^{j_max t}` coefficient of `W_χ(t)` is `F_{−j_max}(t) = ∫_ℝ χR_{−j_max}e^{−iνt}dν` (χ extended evenly).
2. `W_χ − Φ = O(e^{(j_max−ε)t})` would require `F_{−j_max}(t) = O(e^{−εt})` together with the lower terms:
   `e^{3t}F_{−3} + e^{t}F_{−1}` must be `O(e^{(3−ε)t})`, and `F_{−1}` is bounded.
3. By Paley–Wiener, that would make `χR_{−j_max}` analytic in a strip. A compactly supported function that is analytic in a strip is
   identically zero.
4. The Abel transform maps `e^{c't}` growth to `e^{(c'−½)β}` growth. ∎

**Corollary.**
* No compactly supported weight passes the K1′ Tonelli test at any `J₀ > ½`. With a `C^∞` weight it passes only up to `J₀ = ½`.
* So the claim "Proposition K for the heavy-normalised kernel" cannot be realised by a weight on the physical window `ν ≤ Λ`.

**The stored weights [NE] (`hk_real_q_d_flat.txt`; Q also by Stage 7's off-grid check).**

| weight | residual of the de-shift identity | rate | identity fails (`|W − Φ| ~ |Φ|`) at | image correct for |
|---|---|---|---|---|
| Q `C4_value` (XP, 64 digits) | `1.7·10⁻⁵⁴ e^{3t}` (`3.4·10⁻¹⁵` at `t = 30`, `0.024` at 40) | 3 = `j_max` | `t ≈ 41` | `β ≲ 42` |
| Q `C4_deriv` | `9.9·10⁻⁵⁵ e^{3t}` | 3 | `t ≈ 45` | `β ≲ 45` |
| D `X` (300, 1000; double precision) | `≈ 10⁻¹⁶e^t‖X‖₁` floor; `|E−Φ| = 2.5·10⁻⁴` vs `Φ = 1.3·10⁻³` at `t = 22` (300) | 1 = `j_max` | `t ≈ 23` (300), `≈ 22` (1000) | `β ≲ 25` |

* These functionals are exactly what HYB's "`β > 30` [M]" region is about.
* As position-space objects (finite ν-sums of `C_ν`) they are swappable **only** if `J₀ < ½`. The states that cause trouble have
  `J/m ≳ e^{12}–e^{21}`. They are irrelevant numerically, but they are decisive for a theorem.

---

## 4. Proposition K_c (analytic cutoff)

### 4.1 Inputs
* **K1′** (uniform KKSDZ bound in the corner, with KKSDZ's loophole) [Q + D-sketch].
* **K3 / Lemma S:** per-block Rindler positivity ⇒ Cauchy–Schwarz domination [D].
* **LRC** (the K2 row "η → ∞").
  * In the corner the per-ν kernel is `r^a N_com(ν)Σ_jR_jP_{−½+iν+j}(η)`.
  * Subleading terms carry an extra factor `r` and at most the raw growth.
  * This is the input of H7 §14's Fubini corollary: [P given K1′] for scalar-type kernels, [D] for spinning components.
* **K2ν** (new, mild).
  * Away from the corner, `|C_ν(x)| ≤ A(x)e^{γν}` for some finite `γ`, with `A` integrable against `|dDisc 𝒢|` (endpoints and edges as in
    K2).
  * V4FINAL quotes `|C_ν| ~ ν^{−½}`, so any `γ > 0` suffices. Grade [D-sketch].
* **DS** (de-shift existence).
  * An exact `χ` with `χR_j ∈ L¹(ℝ)`, exponentially decaying, even, `χ = O(ν²)` at 0, and `W_χ = Φ` for all `t ≥ 0`.
  * D: the Abel operator is bounded and solvable [D + NE].
  * Q C4: [NE] (identity to `4·10⁻¹⁵` at `t ≤ 30`; the residual is a pure `1.7·10⁻⁵⁴e^{3t}` artefact, §3).
  * Data: `χ(ν)/ν² → 3.24·10⁻⁷`, so `χR_j` is regular at 0 [NE].

### 4.2 Lemma S2: strip analyticity forced by the k = 2 identity [D] + [NE]
**Definitions.**
* `f = χh₋` on `ν > 0` and `f(−ν) = χh₊` (EE_FOUNDATIONS §3.1).
* `f̂(t) = ∫_ℝ f e^{−iνt}dν = C[X] − iS[uX]`, so that `E(t) = ½[e^t f̂(t) + e^{−t}f̂(−t)]`.

**Statement.** `E = Φ` on `t ≥ 0` gives `g(t) := e^t f̂(t) = 2Φ(t) − e^{−t}f̂(−t)`, with:
* `|g| ≤ 2‖Φ‖_∞ + ‖f‖₁` on `t ≥ 0`;
* `|g| ≤ e^t‖f‖₁` on `t < 0`;
* `g = O(t^{−2})` at `+∞` (Prop K step 2).

**Consequences.**
* `f` continues analytically to `−1 < Im ν < 0`.
* Its boundary function `f(·−i) = 𝓕⁻¹[g]` exists, and `‖f(·−is)‖₁` is uniformly bounded for `s ∈ [0,1]`.

**Data [NE]** (D, 300 and 1000): `|g − 2Φ|` falls from `8·10⁻⁵` (`t = 8`) to `2.4·10⁻⁷` (`t = 14`), against `|Φ| ~ 3·10⁻³`. The data
then reach the double-precision floor.

### 4.3 Lemma S4: the k = 4 identity [D] + [NE]
**Definitions.**
* `F_j(t) = ∫_ℝ χR_je^{−iνt}dν`.
* `W = ½Σ_je^{−jt}F_j`.
* `G(t) = e^{3t}F_{−3} + e^{t}F_{−1}`, so that `W = ½[G(t) + G(−t)]`.

**Statement [D].** `W = Φ` on `t ≥ 0` gives `G(t) = 2Φ(t) − G(−t)` and
`|F_{−3}(t)| ≤ e^{−3t}(2|Φ| + ‖χR_3‖₁ + ‖χR_1‖₁) + e^{−2t}‖χR_{−1}‖₁`.
So `χR_{−3}` is analytic in `−2 < Im ν < 0`, with uniformly `L¹` lines. This uses only the identity.

**Data [NE]** (XP, `hk_real_q_d_flat.txt`): the actual structure is stronger.
* **`|F_{−3}(t)|e^{3t} = 2.0000` for `t ∈ [8, 35]` (value) and `= t² = 2|Φ|` (derivative).** So `χR_{−3}` is analytic down to `Im ν = −3`,
  apart from the single singularity at `−3i` that carries the target exactly.
* `F_{−1}` decays like `e^{−2.87t}` (value) and `e^{−2.74t}` (derivative), between `t = 10` and `30`.
* `|G − 2Φ|` falls to `2·10⁻²⁷` at `t = 18`, then reaches the floor.

### 4.4 Lemma C: the cutoff residual [D]
**Setting.** Let `c` be even and entire, with `sup_{|Im ν|≤j_max}|c| ≤ K_c`, and let `δ = χ(c − 1)`.

**k = 2.**
* For `t ≥ 0`, shifting the contour `ℝ → ℝ − i` (Lemma S2) gives `e^tδ̂(t) = ∫f(ν−i)(c(ν−i)−1)e^{−iνt}dν`.
* Hence `|E_c(t) − Φ(t)| ≤ ε_c := ½[‖f(·−i)(c(·−i)−1)‖₁ + ‖f(c−1)‖₁]` for all `t ≥ 0`.
* Since `Abel[1] = P_{−½}`, this gives **`|I_c − I| ≤ ε_c P_{−½}(cosh β)`**.

**k = 4.**
* By Lemma S4, `e^{3t}δF_{−3}(t) = O(ε_c e^{(1+η)t})` for every `η > 0`.
* Also `|e^tδF_{−1}| ≤ ε_c e^t`.
* Hence `|W_c − Φ| ≤ C_η ε_c e^{(1+η)t}` and `|I_c − I| ≤ C_ηε_c e^{(½+η)β}` [D].
* With the NE structure of S4, the `−3i` residue contributes `(c(−3i)−1)·2Φ` and the rest decays. So `|W_c − Φ| ≤ ε_c(1 + |Φ|)` and
  `|I_c − I| ≤ ε_c(P_{−½} + |D|)` [D + NE].

**Size of `ε_c`.**
* It is dominated by `∫_{ν≳ν_c}|χR|` and by the leakage `erfc((ν_c − ν_bulk)/a)`.
* Q (`ν_c = 150`, `a = 6`): `χ ~ e^{−1.46ν}` gives `ε_c ≲ e^{−180}`, with leakage `≲ e^{−400}`.
* D (`ν_c = Λ+100`, `a = 10`): the leakage `erfc(10)‖X‖ ≈ 10⁻⁴¹` dominates, since the exact `X`-tail at `ν_c` is `≲ e^{−170}`.

### 4.5 Proposition K_c [D given K1′, K2ν, LRC, K3, DS]
**Statement.** For each cutoff piece `Ĉ_c = ∫χ_cΨ̂_ν dν`:
1. **(kernel)** `K_c = ∫w_cC_ν dν` with `w_c = χ_c/N_com` converges absolutely at every point of the CLS-D domain.
   * The reason: `N_com⁻¹ ~ e^{πν}ν^{−9}` and `|χ| ~ e^{−(π/2)ν}`, while `c` decays super-Gaussianly (K2ν).
   * So the e^{κν} obstruction of V4FINAL is removed.
2. **(corner)** By LRC, `|K_c| ≤ C r^a(|I_c(η)| + r·η^{c_raw}∫|χ_c|poly)`.
   * Lemma C bounds `I_c`.
   * The subleading term counts as `∫η^{J₀−5/2}`.
   * So the corner integral against K1′ converges for `J₀ < 3/2` (`k = 2`) and `J₀ < 5/2` (`k = 4`; `7/2` given NE).
3. **(swap)** K3 + Tonelli give `Σ_O Tr[P_OĈ_c[G_O]] = Ĉ_c[𝒢_heavy] = −Ĉ_c[light]`, with absolute convergence.
   * Per state, `Ĉ_c[G_O] = ∫χ_cΨ̂_ν[G_O]dν` (Fubini for one block).
   * **This is H_K for `Ĉ_c`.**
4. **(light side)**
   * The `J₀ < 4` terms are still 0 pointwise (Lemma A + D_LR, ν-independent).
   * The (P) moments change by `∫|χ(c−1)R_j||ν−ij|^{2n} ≤ ε_c ν_c^{2n}`.
   * The graviton pole and the `|ĝ₃|²` term change by `O(ε_c)`.

### 4.6 `Ĉ₅^{ds}` and `Ĉ₃` [D given K1′, LRC]
* The raw per-ν `Ψ_{5,ν}` has `c = 3/2 < 2`, so each per-ν sum rule is swappable under K1′.
* `|a_O(ν)| ≤ Σ_j|R_j(ν)|P_{j−½}(cosh β)` (since `|P_{μ+iν}| ≤ P_μ` on `z ≥ 1`). Hence
  `Σ_O∫|χ₅||a_O| ≤ ∫|χ₅|poly · ∫βe^{3β/2}dμ₅ < ∞`, and Fubini over (ν, states) gives H_K for `Ĉ₅^{ds}`.
* This holds for any weight with `∫|χ₅|poly < ∞`, whatever the quality of the C5 de-shift. Stage 7's C5 residual (`8·10⁻⁴` at `t = 20`)
  matters for per-state positivity at `β ≳ 15`, not for H_K.
* `Ĉ₃` is not de-shifted: K2 alone suffices.

### 4.7 Checks
**Toys [DV]** (`hk_toy.txt`). These are exactly solvable analogues with the real weights' decay `e^{−(π/2)ν}`:
* `k = 2`: `χ = sech(πν/2)`;
* `k = 4`: `χ = ½(1+ν²)sech(πν/2)`, with symbol `cosh³t`.

| toy, `ν_c`, `a` | truncation: growth rate of `δW` on `[20, 40]` | cutoff: `δW` for `t ≥ 10` | `|δI_C|/P_{−½}` at `β = 10–40` |
|---|---|---|---|
| k=2, 20, 2 | 0.97 (raw 1) | constant `2.2·10⁻⁴⁵` `= (1 − c(±i))` | `≤ 4·10⁻¹⁵` |
| k=2, 30, 2.5 | 0.97 | `1.6·10⁻⁶⁴` | `≤ 3·10⁻²¹` |
| k=4, 20, 2 | **2.97** (raw 3); `δI_T = 1.7·10⁷` at `β = 20` against `P = 6·10⁻⁴` | `5.8·10⁻⁴⁵` | `≤ 3·10⁻¹²` |
| k=4, 30, 2.5 | 2.97; `δI_T = 4·10²²` at `β = 40` | `4.6·10⁻⁶⁴` | `≤ 2·10⁻¹⁸` |

Lemma C's prediction is borne out: the asymptotic cutoff residual is `(c(±i·shift) − 1)·Φ`, which is bounded.

**Real Q weights, cut inside their own grid [NE]** (`hk_real_qcut.txt`). XP `C4_value`, `ν_c = 50`, `a = 3`, with
`∫_{ν>50}|χ| = 2.3·10⁻²⁹`:

| `t` | 0 | 5 | 8 | 14 | 20 | 30 |
|---|---|---|---|---|---|---|
| truncation `|W − Φ|` | `1.8·10⁻¹⁶` | `1.3·10⁻¹¹` | `6·10⁻⁹` | `4.1` | `1.6·10⁶` | `6.7·10²⁰` |
| cutoff `|W − Φ|` | `4.4·10⁻¹⁵` | `8·10⁻³⁴` | `1.5·10⁻⁵⁴` | `3·10⁻⁴⁷` | `4·10⁻³⁹` | `6·10⁻²⁶` |

* The cutoff residual falls by 39 orders of magnitude. Its only `e^{3t}` component has coefficient `≈ 3·10⁻⁶⁵`: the 64-digit data floor,
  below even the uncut stored weights' own `1.7·10⁻⁵⁴`.
* Truncation grows at the raw rate.
* The derivative weight is identical in pattern.

---

## 5. From Ω_c back to the bounds

* **Per-state images.**
  * On the grid where per-state images are computed (`ν ≲ 25–30`), `c = 1 − O(e^{−400})`. Every computed per-state value is therefore
    unchanged to `≲ e^{−180}` (Q) and `≲ 10⁻⁴¹` (D) relative.
  * At leading Regge, `|I_c − I| ≤ ε_cP_{−½}` (D) and `≤ ε_c(P_{−½} + |D|)` (Q, given NE).
  * So the per-state positivity hypotheses for `Ω_c` are those already stated, up to relative `O(ε_cβ²)`. The only exception is at
    `β ≳ ε_c^{−1/2}` (`J/m ≳ exp(e^{40})`).
* **Remainder.** If one keeps the ideal functional `Ω`, the sum rule holds up to `|R| ≤ ε_c Σ_O Tr P_O·P_{−½}(cosh β_O) + |L − L_c|`.
  * The sum is finite by K1′.
  * For `Ĉ′_EE` it equals the light side of the forward (`Φ = 1`) de-shifted functional, i.e. `O(16πG/(9/4))`.
  * So `R/(16πG) = O(ε_c)` [D].
* **Recommendation.** State the theorem with the cutoff functionals:
  * `χ_c = χ·c`;
  * Q: `ν_c = 150`, `a = 6`;
  * D: `ν_c = Λ+100`, `a = 10`.

  Then H_K holds **exactly** for them [D given the inputs], and no `R` appears.

---

## 6. Side finding: the uncut weight beyond the light-cone edge [M]

**The model.** In the flat model, a heavy state's image at `x = ν/m > 1` (momentum transfer above `m`, outside the physical region)
grows like `|P_J(1−2x²)| ~ e^{2J arccosh x}`. The twist gap gives `J ≤ (m²/(Δ_g−2) − Δ_g + 1)/2`. The table gives the worst exponent of
`log w(ν) + 2J arccosh(ν/m)` (`hk_real_q_d_flat.txt`):

| weight tail | `Δ_g = 300` | `Δ_g = 1000` |
|---|---|---|
| uncut, `κ = π/2` | `m_crit = 2.74Λ` (exponent up to `+2.4·10⁴`) | `2.74Λ` |
| uncut, `κ = 1.46` (measured) | `2.59Λ` | `2.59Λ` |
| D's `X`-tail, uncut | `2.19Λ` | `2.17Λ` |
| Q cutoff (150, 6) | none; worst `−1086` | none; `−2.2·10⁴` |
| D cutoff (`Λ+100`, 10) | none; worst `−19.6` (`= log 3·10⁻⁹`, at `J = 0`) | none |

**Consequence.** Taken literally to `ν → ∞`, the exact de-shift is not even per-state finite on near-maximal-spin states with
`m ≳ 2.2–2.7Λ` (flat model). The computed numbers never saw this, because every stored weight is effectively cut at
`ν ≤ 123` (Q) or `≈ Λ + 20` (D). This is a second, independent reason to put the cutoff in the definition.

---

## 7. Grades

| item | grade |
|---|---|
| Lemma PW (compact support ⇒ raw regrowth; Tonelli needs `J₀ < ½`) | **[P]** |
| Stored weights regrow at the raw rate (residuals, breakdown `β ≈ 25` D / `42` Q) | [NE] |
| Lemma S2 (`k = 2` strip analyticity from the identity) | [D]; data [NE] to `t ≈ 14` |
| Lemma S4 (`χR_{−3}` analytic to depth 2 from the identity) | [D] |
| Stronger Q structure (`e^{3t}F_{−3} → 2Φ`, `F_{−1} ~ e^{−2.8t}`) | [NE] |
| Lemma C (cutoff residual bounds) | [D]; toys [DV]; real Q weights [NE] |
| Proposition K_c (H_K for the cutoff functionals, `J₀ < 3/2`) | **[D given K1′, K2ν, LRC, K3, DS]** |
| H_K for `Ĉ₅^{ds}` (plain Fubini) | [D given K1′, LRC] |
| Transfer to the computed per-state values (unchanged to `e^{−180}` / `10⁻⁴¹`) | [D + NE] |
| K2ν | [D-sketch] (new; mild) |
| §6, flat-model blow-up of the uncut tail | [M] |
| **Net H_K** | **[D-sketch]** by weakest link (K1′, K2ν). It is **no longer an independent hypothesis** and **no longer needs `J₀ < ½`**. |

---

## 8. Wording for the letter (replaces the H_K sentence; addresses M8)

> *Interchange (H_K).* The de-shifted pieces are defined with an analytic ν-cutoff (`χ → χ·c`, `c` an entire erf-plateau at
> `ν_c ≫ 1`). This changes no computed number (relative effect `< 10⁻⁴⁰`) and makes each piece a genuine position-space functional.
> For these functionals, the interchange of the heavy OPE sum with the (Regge-corner, ν) integrations follows from Tonelli's theorem:
> * the inputs are the KKSDZ Regge bound `J₀ ≤ 1` (made uniform in the transverse angle, K1′), per-block Rindler positivity, and the
>   leading-Regge form of the kernels;
> * the de-shift identity forces strip analyticity of the weights, so the cutoff spoils the cancellation only by `O(ε_c)` relative and
>   never at the raw rate;
> * the heavy sum converges absolutely for any `J₀ < 3/2`, which leaves a margin of ½ over the known bound.
>
> Without the cutoff no such argument exists. A weight of compact ν-support (in particular any grid weight) regrows at the raw
> rate `e^{5β/2}` by Paley–Wiener, so it would require `J₀ < ½`. We grade the interchange [D-sketch]. It inherits K1′ and the
> uniformity of the kernel bounds, which are written out in the Supplement.

---

## 9. What remains

1. **K1′** (KKSDZ's bound made uniform in the angle), and KKSDZ's loophole. Unchanged; still [D-sketch] / [Q].
2. **K2ν and LRC in position space.**
   * Write out the corner expansion of the spinning per-ν kernels, with subleading-in-`r` terms bounded uniformly in ν up to `e^{γν}`.
   * LRC is H7 §14's existing input; K2ν is new.
3. **DS for Q, analytically.** Existence of an exact C4 solution. The [NE] is strong: the `e^{3t}F_{−3} = 2Φ` structure holds to
   5 digits over `t ∈ [8, 35]`, and the residual is a pure precision artefact.
4. **Spinning components** of Lemma S/K3 for `EE` at `k = 4`. Inherited; [D].
5. **Not edited here.** RESULT_D_STATEMENT, RESULT_Q_STATEMENT and UNIFIED_THEOREM should add the cutoff `c` to the definitions of
   `Ω^EE` and `Ω^Q` (fixed numbers), and regrade H_K as in §7.
6. **Beyond H_K.**
   * The C5 de-shift is poor at `t ≳ 15`. This affects per-state positivity at `β ≳ 15` (HYB), not H_K.
   * §6 [M] says that per-state images beyond the light-cone edge must never be sampled by the weights. The cutoff ensures this.

## 10. Reproduce

    cd ads4/code/hk
    ulimit -v 3000000
    OMP_NUM_THREADS=1 python3 hk_toy.py            # ~6 min -> cache/hk_toy.txt
    OMP_NUM_THREADS=1 python3 hk_real.py q d flat  # ~30 s  -> cache/hk_real_q_d_flat.txt
    OMP_NUM_THREADS=1 python3 hk_real.py qcut      # ~30 s  -> cache/hk_real_qcut.txt
