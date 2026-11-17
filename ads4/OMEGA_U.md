# `Ω^U`: the corrected AdS uplift (B′ basis), heavy side from the truncated CHMRS scalar-skeleton series, light side, conditional bound

Status: written 17 Nov 2026 by a derivation/numerics agent. **Refereed once (`REFEREE_OMEGA_U.md`); revised (list in §5). Not a theorem about AdS4.** This note replaces the AdS analyses of
STATUS links 9–11, which were done for the old functional `Ω`; `FLAT_UNSTRIPPED.md` showed that `Ω` has the wrong flat limit (`F_U`).
Code: `ads4/code/omegau/` (`bash run_all.sh`, ≈ 30 min on 2 cores; `bash run_all.sh fast` ≈ 1 min once the caches exist). Outputs are in
`omegau/cache/`, which is git-ignored. The code imports, read-only:
* `robust/{rcommon,check,lp}.py` and `robust/functionals_certified.json`;
* `uniform/chmrs_series.py`, plus `uniform/cache/{chmrs_series,remainder_coeffs}.pkl` (the CHMRS order-1/m⁸ ancillary series).

`remk.py` is an adapted copy of `uniform/remainder.py` and `uniform/skelfast.py`, extended to `k = 2…7`. Nothing outside `omegau/` was
modified, apart from one line added to `.gitignore`.

Marks:
* **[D]** derived, here exact (sympy);
* **[DV]** derived and checked numerically;
* **[DVE]** checked end to end on exact blocks (quoted);
* **[NE]** floating-point numerics;
* **[M]** a statement about a stated model or cut;
* **[I]** inferred or assumed;
* **[Q]** quoted.

Nothing here is fitted.

Notation follows `ADS4_FUNCTIONAL.md`, `UNIFORM_REMAINDER.md` and `FLAT_UNSTRIPPED.md`:
* `Λ = √((Δ_g−2)(Δ_g−1))`, `ν = Λp`, `u = Λ²/m²`, `x = ν/m = p√u`;
* normal form = action per `|c^h|²` times `m⁴/Λ⁴`;
* `D₊₊ = P_J`, and `D₊₋ = d^J_{44} = (1−x²)⁴P^{(0,8)}_{J−4}` (the full Wigner function).

## 0. Results

**The functional.**

`Ω^U = ∫₀¹dp Σ_{k=2}^{7} φ_k(p) ε_kΛ^{2k} C^{(−)^k}_{k,Λp} − aΛ¹⁰∂_{ν²}C⁺_{4,ν}|₀ + (5a/2)Λ¹⁰C⁻_{5,0}`

* `ε_k = 1` for even `k` and `½` for odd `k`.
* The weights are those of `FLAT_UNSTRIPPED.md` §0 item 5, built from the M02, M01 or F_g3 weights `ψ₂, ψ₃, a`.
* The forward coefficients are `−5, 25/2` for F_g3.

1. **Flat limit [D, exact; DV].**
   * Under the B′ dictionary `ε_k m^{2k}C^h_k → e_k(x²)D_h(1−2x²)`, the integrand of `Ω^U` equals the stripped integrand
     `ψ₂(2−x²)D̃ + σuψ₃D̃` **pointwise in `p`**, for both helicities. Here:
     * even `k`: `e_k = (2−x²)/(1−x²)^{1+k/2}`, which is the CHMRS bulk-point density;
     * odd `k`: `e_k = σ/(1−x²)^{(k+1)/2}`.
   * The forward terms equal the stripped ones: `2aJJ` for `++` and `a(2JJ−42)` for `+−` (`u1`, sympy).
   * Numerically, the sum over `k` of the B′ parts minus the stripped action is `≤ 3·10⁻¹⁰ N_s` for `m ≥ 1.01Λ` (M02, M01, F_g3, `J ≤ 60`).
   * **Cancellation amplification [NE].** The individual B′ terms are much larger than their sum near the gap, in `++` at small `J`.
     The ratio `N_U/N_s = Σ_k|part_k|/(stripped component norm)` is:

     | `m/Λ` | 1.0005 | 1.01 | 1.05 | 1.2 | 1.5 | 2 | ≥ 3 |
     |---|---|---|---|---|---|---|---|
     | `N_U/N_s` (M02, `++, J = 0`) | 182 | 62 | 20 | 4.7 | 1.75 | 1.19 | ≤ 1.08 |

     * At `m = Λ` (`J = 0`) the `k = 6, 7` terms diverge logarithmically, while their sum stays finite.
     * In `+−` the ratio is `≤ 1.7`, because `d₄₄` carries `(1−x²)⁴`.
2. **Odd-`k` prefactor [I, quantified].**
   * The B′ odd density `σ(1−x²)^{−(k+1)/2}` and the CHMRS scalar formula continued to odd `k`, `ε_kE_k`, differ at `O(x⁴)`.
   * If the AdS odd moments tended to the second, the flat action of M02 would shift by `≤ 1.6N_s` at `m = 1.05Λ` (`++`), by `0.09N_s` at
     `1.2Λ`, by `0.007N_s` at `1.5Λ` and by `7·10⁻⁴N_s` at `2Λ`.
   * Above `1.5Λ` the question is therefore immaterial for positivity; below `1.2Λ` it is `O(1)` (`u1`).
   * Only the B′ form is a legitimate fixed-`u` sum rule: half-integer powers of `Z` are not.
3. **AdS heavy action from the scalar skeleton (CHMRS order-1/m⁸ truncated asymptotic series, not exact moments) [Q + D + I; `adsaction.py`].** For each `k` the model is
   `ε_km^{2k}C^h_k = e_kD_h + 2ε_kσ^k[P_{−½+iν}(cosh β) − J₀(b̂p)] + σ^k env_k ρ_k`.
   * Here `ρ_k` is the scalar-skeleton remainder from the truncated CHMRS order-1/m⁸ series (fixed-`J` to `1/m⁴`, fixed-η to `1/m⁸`). Its coefficients are curvature-subtracted for
     `k = 2…7` (`u2`), with `Δφ = 3` (leading-order `⟨TTTT⟩` skeleton) and `Δφ = 2` (MFT choice).
   * For even `k` and `++`, this is the scalar moment up to the series truncation.
   * For `+−`, the full `d^J_{44}` is used (reading B) and `ρ_k` is transferred from the scalar [I].
   * For odd `k`, the scalar series is continued in `k`, and `ρ_k` is transferred to the B′ odd envelope [I].
   * Forward terms: `R₄`, `R₅` at `ν = 0`, with the B′ weight `(5a/2)` on `C⁻_{5,0}`, i.e. 5× the weight in `Ω`.
4. **Sanity: flat limit of the AdS action [NE, `u3`].** At fixed `(m/Λ, J, h)`, the skeleton action (all `k`, forward terms included)
   approaches the flat action of M02 or F_g3:
   * `|F_skel − F_flat|` (absolute, normal form; `F_flat = 0.13–1.1` on the tested states) is `2·10⁻⁴ … 2·10⁻²` at `Δ_g = 100`;
   * it is `≲ 2·10⁻⁶` at `Δ_g = 10⁴`, and `≲ 3·10⁻⁷` at `10⁵`, which is the floating-point floor.

   The decrease is not a clean `1/Λ²`: at `x ≳ 0.6` the truncated `x²`-series is not convergent at fixed `m`. In `+−` at `Δ_g = 1000` it
   is not even monotone (`3.4·10⁻³` at `1.2Λ, J = 4`; `7.8·10⁻⁴` at `1.5Λ, J = 6`, larger than at `Δ_g = 100`), because the truncated series is
   used on nodes the cut calls uncontrolled. The check is also close to tautological: the model is flat plus terms that vanish by construction.
   **It is a consistency check of the code, not evidence for the dictionary I1.**
5. **Positivity scan [NE + M + I; `u4`, `u7`].**
   * **Scope.**
     * `Δ_g = 20, 50, 100, 300, 1000`;
     * states with `m/Λ ∈ [1, 10]`, `b̂ ≤ 40`, `J ≤ 400`: 1454 `++` and 2279 `+−` states at `Δ_g = 100`;
     * both helicities, both skeletons;
     * 64 Gauss nodes in `p`.
   * **Controlled region (heuristic cut, unchanged from `UNIFORM_REMAINDER`).**
     * Fixed-`J` series: the ratio test passes, `≤ 0.8`.
     * Otherwise the fixed-η series at `x ≤ 0.6`.
     * Every state below `1.5Λ` has uncontrolled nodes. They start at `x ≈ 0.43–0.45` for `k = 7` and at `x ≈ 0.55–0.6` for `k = 2`; the
       uncontrolled weight is dominated by the `k = 3–6` terms (`cache/u4_*`).
   * **Criteria.**
     * "abs": the per-`k`, per-node remainder bound `E`, with no cross-`k` cancellation, so the cancellation of item 1 is *not* used.
     * "coh": the signed skeleton remainder.
     * `κ = (F_hyb − E)/U` on states with uncontrolled weight `U` (B′ envelopes).
   * **For `m ≥ 1.5Λ` (M02):**
     * fully controlled margin `≥ 0.18N_s` (`++`) and `≥ 0.35N_s` (`+−`);
     * `κ ≥ 2.6` for `Δ_g ≥ 50` (up to 15–19 at `Δ_g = 1000`);
     * at `Δ_g = 20`, `κ = 0.74` (`Δφ = 3`) or `2.0` (`Δφ = 2`).
   * **Below `1.5Λ` (M02, `κ` "abs", `Δφ = 3`).** `κ` per band, `++ / +−`:

     | `Δ_g` | `[1,1.1)` | `[1.1,1.2)` | `[1.2,1.3)` | `[1.3,1.5)` |
     |---|---|---|---|---|
     | 20 | 0.0065 / – | 0.055 / – | 0.088 / 0.11 | 0.16 / 0.45 |
     | 50 | 0.011 / 0.028 | 0.091 / 0.039 | 0.16 / 0.19 | 0.32 / 0.58 |
     | 100 | 0.015 / 0.013 | 0.099 / 0.046 | 0.20 / 0.20 | 0.42 / 0.57 |
     | 300 | 0.024 / 0.007 | 0.10 / 0.049 | 0.20 / 0.20 | 0.40 / 0.57 |
     | 1000 | 0.44 / 0.013 | 0.66 / 0.054 | 0.98 / 0.21 | 1.6 / 0.68 |

     * The `Δφ = 2` values are similar (`cache/u7_summary.txt`).
     * "coh" and "abs" differ by `≤ 10 %`.
     * The `C_s = 40` spin envelope changes nothing except at `m ≈ Λ`, where `κ < 0`.
     * The `Δ_g = 1000`, `++` row is an artefact of the cut: the absolute switch "last term `≤ 10⁻⁴E_k` counts as controlled" in
       `remk.residual_fixedJ` declares more nodes controlled at large `m`. It is not a physical improvement.
     * **Negative `κ`** (`F_hyb − E < 0` on some state) occurs for M01 at `Δ_g = 20` (down to `−0.11N_s`, `Δφ = 2`) and for F_g3 in `+−`, `J = 4`,
       `m ≈ 1.02–1.08Λ`, at `Δ_g = 50–1000` (down to `−0.06N_s`), and at `Δ_g = 20` (down to `−0.27N_s`). It also occurs at the F_g3 `++`
       island (`−0.013N_s`, `Δφ = 2`, `Δ_g = 50`). M02 has none.
   * **Threshold as a function of the tolerance** (revised after referee M1).
     * Definition: `n(κ_req, δ)` = the smallest grid `m/Λ` such that every scanned state with `m ≥ nΛ` has
       `F_hyb − E − κ_req U − δN_s > 0`.
       * `κ_req` = the tolerated relative AdS error on the uncontrolled nodes, per `k` and per node;
       * `δ` = an allowance for what `E` and `U` omit;
       * worst case over both helicities and both skeletons (`u7`, `cache/u7_summary.txt`).
     * M02 values (`δ = 0`; with `δ = 0.03`, `n` rises by `≤ 0.05`):

       | `Δ_g` | `κ_req = 0.05` | 0.1 | 0.3 | 0.5 | 1 | 2 |
       |---|---|---|---|---|---|---|
       | 20 | 1.22 | 1.26 | 1.40 | 1.50 | 1.60 | 2.15 |
       | 50 | 1.12 | 1.17 | 1.30 | 1.35 | 1.45 | 1.50 |
       | 100 | 1.11 | 1.15 | 1.27 | 1.35 | 1.40 | 1.45 |
       | 300 | 1.11 | 1.15 | 1.27 | 1.35 | 1.40 | 1.45 |
       | 1000 | 1.10 | 1.15 | 1.24 | 1.28 | 1.35 | 1.40 |

     * At `κ_req = 1`: M01 gives 1.65 / 1.50 / 1.45 / 1.45 / 1.40. F_g3 gives 1.60 / – / 1.35 / 1.35 / 1.35, where "–" means its `++` island
       (`m ≈ 5.7Λ`) fails at `Δ_g = 50` with `Δφ = 2`. **With `δ = 0.03`, F_g3 has no threshold at any `Δ_g`**: its island margin is only
       `0.004–0.011N_s`.
     * **The tolerance is a convention, so it needs calibrating [M].**
       * On the uncontrolled nodes of the threshold states, the truncated series itself gives `|ρ_k| ≈ 0.4–2 %` (`≤ 3.5 %` at `1.01Λ`).
         "Uncontrolled" is a verdict of the `x²`-ratio test, not an indicated `O(1)` error.
       * What `E` and `U` omit:
         * the graviton-versus-scalar transfer (I4, I5), which `SKELETON_TEST.md` sizes at a few `/m²` relative for same helicity;
         * the odd-`k` prefactor (I3), ≈ 1–3 % of `N_s` at `1.3–1.4Λ` (item 2).
       * A tolerance of ≈ 10 times the indicated `|ρ|` gives `κ_req = 0.3`; with `δ = 0.03` to cover I3, this is the **central (calibrated, not
         derived) choice**.
       * `κ_req = 1` is the conservative choice.
       * `κ_req ≤ 0.1` trusts the truncated series on nodes where its convergence is not established.
     * **Result for M02:**
       * `n ≈ 1.25–1.4` (central, `κ_req = 0.3`, `δ = 0.03`): 1.40 / 1.35 / 1.28 / 1.28 / 1.24 for `Δ_g = 20 / 50 / 100 / 300 / 1000`;
       * `n ≈ 1.35–1.6` (conservative, `κ_req = 1`);
       * `n ≈ 1.15` if `κ_req = 0.1` is accepted.
     * The binding state at the conservative threshold (`Δ_g = 100`, `1.40Λ`, `++`, `J = 2`) has `F/N_s = 0.283`, `E/N_s = 0.010` and
       `U/N_s = 0.229`. Its spare margin `F − E − U ≈ 0.044N_s` is comparable to the omitted I3/I4 effects.
     * **`n = 1` is not supported at any `κ_req ≥ 0.05`.** On `[Λ, 1.3Λ)`, `κ` is `0.004–0.2` for `Δ_g ≤ 300`. Positivity there needs the
       uncomputed relative AdS correction on the uncontrolled nodes (`x ≳ 0.43–0.6`, including `ν → m`) to be below 0.4–20 %, *per `k`*.
6. **Light side [D given (F), (P), (R), counting rule; `u5`].**
   * **Flat.**
     * `−B′_k|_low`: Einstein (`16πG/q`) appears only in `B′₂`.
     * Contacts `g₄…g₆′` appear in `B′₂…B′₆`.
     * `|ĝ₃|²` appears in `B′₂, B′₃, B′₄` (`2πGq³`, `2πGq²`, `4πGq`).
     * In `Ω^U` all contacts cancel exactly.
     * `k(0)`, `A`, `B` reproduce `FLAT_UNSTRIPPED.md`: M02 `A = 51.7366`, `B_AdS = −86.2981`; F_g3 `37.8387`, `−60.7875`.
   * **Graviton pole.** The pole sits only in `C⁺₂`, with weight `φ₂ = ψ₂` and residue `16πG/(ν²+9/4)` [D given (F), (R)]. So
     `L_AdS(Λ) = ∫ψ₂/(p²+9/(4Λ²))` is the same as for `Ω`.
   * **`k ≥ 3` [I].** By the counting rule (`LIGHT_INPUTS.md` §0 item 1), an Einstein term in a dispersive functional with leading moment `k`
     vanishes iff `J₀(G^E) < k`. Applying it to the project's moments needs:
     * that each `C_k` (`k = 3…7`) is the leading moment of a spin-`k` dispersive functional. For `C⁻₃` this is [I]: it is not the CLS-D bare
       commutator (`LIGHT_INPUTS.md` §0 items 3–4). For `C⁻₅`, `C⁺₆` and `C⁻₇`, never even defined as two-pair moments, it is weaker still.
     * `J₀ ≤ 2` (RB₂) [I, KKSDZ].

     If these hold, the `φ_k` weights multiply zero, and `A`, `B_AdS` are unchanged.
     * Leverage of a test residue `16πG(ν²+9/4)^{−k}` against the leading term, at `Δ_g = 100` (M02): `6.7·10³, 1.9·10³, 6.2·10², 1.4·10², 36`
       for `k = 3…7`; `10¹³` for the forward terms.
     * Only exact zeros are usable, so any failure of these inputs is fatal to the bound.
     * By the same rule, `C⁻₅`, `C⁺₆`, `C⁻₇` carry no `|ĝ₃|²` term (its Regge spin is 4).
   * **Contacts in AdS: the cancellation is not exact.**
     * By (F) and (P), each `−ε_kΛ^{2k}C_k|_g` is a polynomial in `ν²` whose top coefficient is the flat one.
     * The lower coefficients `a_{k,g,j}` (spinning contact Witten diagrams, not computed) do **not** cancel in `Σφ_k`. They leave
       `Σ_k lev_{k,g}a_{k,g}/Λ² + O(Λ⁻⁴)`, now also for `g₄` and `g₅`. In `Ω`, `g₄` dropped out and `g₅` cancelled exactly.
     * The leverages for M02 (normalised, `k(0) = 1`) are in the table below.
     * For the size of `a`, the scalar analogue gives `|a| ≈ 36` at `Δφ = 3` (`−C₂[c₃STU] = ν² − (2Δφ+½)(2Δφ−½)`). For a degree-`d`
       polynomial the next-to-top coefficient is taken as `d·36`.
     * At smearing `Λ/n` the relative residues grow by `n²`. Write `g_f ≤ c_g G_f log(Λ/n)` with `c_g` unknown; the flat two-sided bounds
       carry `O(10–100)` constants [A].
     * At `n = 1.4`, the shift of the bound relative to `16L_AdS/k` is:

       | `Δ_g` | `c_g = 1` | `c_g = 30` |
       |---|---|---|
       | 50 | 0.6 % | 17 % |
       | 100 | 0.1 % | 4 % |
       | 300 | 0.01 % | 0.4 % |
       | 1000 | 0.001 % | 0.03 % |

       (`u5`.) So **`Δ_g = 50` is not clearly negligible; `Δ_g ≥ 100` is**, under EFT power counting.

     | coupling | `k = 2` | `k = 3` | `k = 4` | `k = 5` | forward `k = 4` | forward `k = 5` | `Σ|lev|` |
     |---|---|---|---|---|---|---|---|
     | `g₄` | 0.26 | −0.38 | — | — | — | — | 0.64 |
     | `g₅` | 0.11 | 0.71 | −1.39 | — | — | — | 2.2 |
     | `g₆` | 0.06 | 0.57 | −1.46 | 1.00 | −1.48 | 2.46 | 7.0 |
     | `g₆′` | — | — | — | — | — | — | 2.3 |
     | `|ĝ₃|²` | 0.043 | 0.471 | −1.12 | — | — | — | 1.63 |

     The `|ĝ₃|²` row (corrected by a factor `π` after referee m1; F_g3's `k = 2` entry × `k(0)` = 9677/13860, as in `ADS4_LIGHT.md`) gives
     `k_AdS − k(0)` of about 17 % (`Δ_g = 20`), 0.6 % (100) and 0.006 % (1000) at `n = 1`, |a| = 36. At `n = 1.4` with `d·36` it is
     about 7 % (50), 1.6 % (100) and 0.02 % (1000).
7. **Conditional AdS bound [D given all inputs of §3].** With heavy positivity for `m ≥ nΛ` (smearing `Λ/n`):
   `|ĝ₃|²Λ⁸ ≤ n⁸·16L_AdS(Λ/n)/k(0)`, and `√(t₄²+α_T²) ≤ 240(…)^{1/2}/Λ⁴`.

   M02, `√(t₄²+α_T²) ≤`, at the thresholds of item 5 (`δ = 0`; `u7`):

   | `Δ_g` | `κ_req = 0.1` | central `κ_req = 0.3` (`δ = 0.03`) | conservative `κ_req = 1` | `n = 1.5` | `n = 1` (not supported) |
   |---|---|---|---|---|---|
   | 20 | `0.041` (`n = 1.26`) | `0.060` (1.40) | `0.098` (1.60) | `0.078` | `0.018` |
   | 50 | `8.5·10⁻⁴` (1.17) | `1.5·10⁻³` (1.35) | `1.9·10⁻³` (1.45) | `2.2·10⁻³` | `4.7·10⁻⁴` |
   | 100 | `5.4·10⁻⁵` (1.15) | **`8.1·10⁻⁵`** (1.28) | `1.1·10⁻⁴` (1.40) | `1.5·10⁻⁴` | `3.2·10⁻⁵` |
   | 300 | `7.5·10⁻⁷` (1.15) | `1.1·10⁻⁶` (1.28) | `1.6·10⁻⁶` (1.40) | `2.1·10⁻⁶` | `4.4·10⁻⁷` |
   | 1000 | `6.9·10⁻⁹` (1.15) | `9.2·10⁻⁹` (1.24) | `1.3·10⁻⁸` (1.35) | `1.9·10⁻⁸` | `4·10⁻⁹` |

   * **The honest statement for `Δ_g = 100` is the range `5.4·10⁻⁵ – 1.5·10⁻⁴`**, with central value `≈ 8·10⁻⁵`, depending on the tolerance
     convention [M]. All of it is conditional on the inputs of §3.
   * F_g3 at `n = 1` gives `0.015 / 4.1·10⁻⁴ / 2.7·10⁻⁵ / 3.8·10⁻⁷ / 3.4·10⁻⁹`; it has no supported threshold (item 5).
   * `Δ_g = 20` is meaningless: `k_AdS − k(0) ≈ 17 %` already at `n = 1`, and `Δ_g = 50` is marginal (item 6).
   * A threshold quoted to `±0.05` changes the bound by `±25 %`.
8. **Comparison with the superseded `Ω` numbers.**
   * The light-side constants `A`, `B_AdS`, `L_AdS` are unchanged: the pole weight is `ψ₂` in both functionals.
   * The M02 threshold is comparable to that for `Ω` (`n ≈ 1.5`, `1.5·10⁻⁴` at `Δ_g = 100`). It is now given as `n(κ_req)`:
     1.27 (central) – 1.40 (conservative) at `Δ_g = 100`. `Ω`'s `n ≈ 1.5` came from the same kind of `κ`-convention.
   * The near-gap `κ` drop. For `Ω` on `[1, 1.5)`, `κ` was `0.36–0.84 (++) / 0.29–0.66 (+−)`. For `Ω^U` the minimum over `[1, 1.1)` is
     `0.004–0.03` (the `0.44` at `++`, `Δ_g = 1000` is a cut artefact, item 5). The reason is that `Ω^U`'s B′ envelopes `E_k ∝ (1−x²)^{−1−k/2}` are much larger near
     `x → 1` (item 1), and the "abs" bound uses no cross-`k` cancellation.
   * Regge regime (`u6`, leading Regge action):
     * `Ω^U`(M02) has no negative grid point in `++`, and in `+−` its last negative point is `1.16Λ` (outside the validity of the Regge kernel);
       with `C₁ = 1` the last negative point is `1.81 / 1.44Λ`, against `2.14 / 1.18Λ` for `Ω`.
     * For F_g3 the leading-order `++` island disappears (`0 < 7.38Λ`), but with `C₁ = 1` it returns (`7.2Λ`).
   * `Ω` additionally had a wrong flat limit (`F_U`, which is not heavy-positive and keeps its contacts). So all of `Ω`'s AdS numbers are
     superseded, even where they coincide.

## 1. Heavy side: details

* **Grid and quadrature.**
  * The grids are those of `uniform/t3_scan.py`.
  * The flat parts come from the stripped functional, which is exact Gauss for polynomials.
  * The per-`k` B′ parts (used only in `N_U`) come from 400-node Gauss. At `m = Λ`, `J = 0` they cancel only to `3·10⁻³N_s`, the
    log-divergent pieces.
* **Hybrid.** The Regge curvature per `k` is `2ε_kσ^k u^{k−2}∫φ_k[P_{−½+iΛp}(cosh β) − J₀(b̂p)]`, by exact Mehler–Dirichlet transforms of
  the polynomial `φ_k` (`common.Transforms`). The forward curvature is `a u³(2m²W − 2JJ − ½) + 5aσu³(P₀ − 1)`.
* **Remainder coefficients** (`u2`, `cache/u2.log`). `r̃₁(x = 0, J = 0)` for `k = 2…7` is:
  * `Δφ = 3`: `+12.9, +1.1, −20.6, −54.4, −102, −166`;
  * `Δφ = 2`: `−19.1 … −238`.

  The `1/m⁴` coefficients reach `7·10⁴` at `k = 7`. The high-`k` moments are therefore worse-behaved asymptotic series, which matters only
  near the gap, where `u^{k−2} ≈ 1`.
* **Forward step.** `∂_{ν²}` uses a step `h = max(0.05, 0.1/β)`, which avoids round-off of `m²∂_{ν²}` at large `m`.
* **What is not covered.**
  * `m > 10Λ` and `b̂ > 40`. These are covered as before by the certified flat margin plus the Regge statements, [I].
  * The uncontrolled nodes themselves.
  * Parity-odd heavy exchanges.

## 2. Light side: details (`u5`, `cache/u5.log`)

* **Low sides.** The `−B′_k|_low` for `k = 2…8` are recomputed independently of `flatunstripped/sumrules.py`. Their values: `B′₇` and `B′₈`
  are 0 through `g₆`; `B′₆ = 2g₆`.
* **Leverage.** `lev_{k,g} = (top flat coefficient) × ∫φ_k p^{2d−2}/k(0)`.
  * The forward pieces enter only where the AdS polynomial has a `ν²`-coefficient at `O(Λ⁻²)`, i.e. top degree 2 in `B′₄`, or a constant
    in `B′₅`.
  * `g₅` in the forward terms is exact: `∂_{ν²}` of a degree-1 polynomial, and the constant in `C⁻₅`.
* **`|ĝ₃|²`.** The `k(0)` coefficient is exact. The AdS corrections `b, f, e` (three unknowns at `O(Λ⁻²)`) are not computed.
* **Linear `t₄`, (R), double traces.** These are as in `ADS4_LIGHT.md` §2.3 and §3. The `k = 3` term could again be replaced by the CLS-D
  bare commutator (`OMEGA_PRIME.md`), which removes the RB₂ item for `k = 3` only.

## 3. Inputs

**Updated since `Ω`.** `SPINNING_TTTT.md` §0 item 5 (unrefereed) checks `C^± = (2P/m^{2k})(|c_same|² ± |c_opp|²)` for `⟨TTTT⟩` end to end
[DVE], including the absolute constant ¼ and the sign of the opposite-helicity term (`Σ_−`: `0.92 → 0.97`). The check is at the percent
level at `m = 70`, and covers only `k = 6`, `ν = 0`, `β = 1`, with two masses. So the "heavy formula for `⟨TTTT⟩`" input is now [DVE] at
that point. It remains [I] for:
* `k = 2, 3, 4, 5, 7`;
* the odd-`k` two-pair `C⁻` (`ε_odd = ½`);
* `ν ≠ 0`, and above all `ν ~ m`.

| # | input | status | used for |
|---|---|---|---|
| I1 | B′ dictionary `ε_kΛ^{2k}C_k → B′_k`, i.e. the spinning bulk-point limit with full `d^J_{44}` (reading B for gravitons) | [I]; photons [NE] (`SKELETON_TEST` §8); gravitons not tested | flat limit = F_g3 / M02 (item 1) |
| I2 | two-pair heavy formula for `⟨TTTT⟩`, leading order | [DVE] at `k = 6`, `ν = 0`, `β = 1` (percent level, unrefereed); [I] at `k = 2…5, 7` | Regge part, normalisation |
| I3 | odd-`k` two-pair `C⁻` moments (`k = 3, 5, 7`): never computed; their normalisation `ε_odd = ½`; the odd B′ prefactor `σ(1−x²)^{−(k+1)/2}` rather than the `k`-continued scalar `E_k` | [I]; sensitivity `O(1)` below `1.2Λ`, `< 1 %` above `1.5Λ` (item 2) | `k = 3, 5, 7` terms |
| I4 | the scalar skeleton's relative remainder transfers to graviton moments, per `k`, both helicities; skeleton `Δφ ∈ [2, 3]` (currents: `Δφ* ≈ 2.4–2.55`) | [I] | `E`, `S` |
| I5 | odd-`k` remainder = CHMRS series continued in `k` | [I] | `ρ_{3,5,7}` |
| I6 | heuristic control cut (ratio `≤ 0.8`, fallback `x ≤ 0.6`, absolute `10⁻⁴E_k` switch); tolerance `κ_req` (0.3 central, 1 conservative) and allowance `δ` | [M] | `n` (item 5 table) |
| I7 | spin-specific `1/m²` terms `≤ 40 cosh β/m²` per `k` (or `k`-independent, then harmless) | [I] | `C_s` column |
| I8 | (F) flat limit of the spinning light side; (P) polynomial light side of contacts and `W³⊗W³` | [I] (Q for scalars) | items 6–7 |
| I9 | (R) `−C⁺₂|_E = 16πG/(ν²+9/4)` exactly | [D-sketch given IP] | `B_AdS` |
| I10 | counting rule applies to the project's `C_k`, `k = 3…7` (each the leading moment of a spin-`k` dispersive functional: [I] for `C⁻₃`, weaker for `C⁻₅, C⁺₆, C⁻₇`) + RB₂ (`J₀(G^E) ≤ 2`) ⇒ no Einstein term in `C_{3…7}`, no `|ĝ₃|²` in `C_{5,6,7}` | **[I]** (the rule itself is [D from Q] for CLS-D-type functionals; RB₂ [I], KKSDZ) | leverage `36 – 6.7·10³` (smeared, `Δ_g = 100`), `10¹³` (forward): load-bearing |
| I11 | EFT power counting of contacts, `g_f ≤ c_g G_f log Λ` (`c_g` unknown, `O(1–30)`); AdS next-to-top coefficients `≈ d·36` | [A] | contact residues; `Δ_g = 50` marginal |
| I12 | no other light single traces; tree level `c_T ≫ Δ_g²`; double traces at double zeros (unshown for spinning kernels) | [A] / [I] | light side |
| I13 | parity-odd heavy exchanges | not treated | — |
| I14 | convergence/swappability of `C⁺₂` (graviton-marginal) and of the higher `C_k` | [I] | all |

## 4. Code (`ads4/code/omegau/`)

| file | content | log |
|---|---|---|
| `common.py` | `φ_k`, B′ densities, flat per-`k` parts, AdS kinematics, Regge/flat transforms | — |
| `u1_flat.py` | exact pointwise identity, forward map, numeric identity, `N_U/N_s`, odd-prefactor sensitivity | `cache/u1.log` |
| `u2_coeffs.py` | curvature-subtracted remainder coefficients, `k = 2…7`, `Δφ = 3, 2, 2.5` | `cache/u2.log` |
| `remk.py` | per-`k` remainder, control cut, forward remainders (adapted from `uniform/`) | — |
| `adsaction.py` | per-state AdS action of `Ω^U` (hybrid, `S`, `E`, `U`) | — |
| `u3_flatlimit.py` | AdS skeleton action → flat action | `cache/u3_Df3.log` |
| `u4_scan.py` | positivity scan; `OMEGAU_DF = 3` or `2` | `cache/u4_<Δ_g>_Df*.txt` |
| `u5_light.py` | low sides, contact cancellation, AdS leverages, Einstein leverages, bound vs `n` | `cache/u5.log` |
| `u6_regge.py` | leading Regge action, `Ω^U` vs `Ω` | `cache/u6.log` |
| `u7_summary.py` | band `κ`; `n(κ_req, δ)` and the bound at it; negative-`κ` states | `cache/u7_summary.txt` |
| `run_all.sh` | driver, ≤ 2 cores | — |

## 5. Revisions after `REFEREE_OMEGA_U.md`

* **M1.** The single "supported threshold" is replaced by `n(κ_req, δ)` for all `Δ_g` (item 5). The tolerance is calibrated against the
  truncated-series `|ρ|` (0.4–2 %) and the omitted I3/I4 effects (`δ = 0.03`). The bound is now a range: `5.4·10⁻⁵ – 1.5·10⁻⁴` at `Δ_g = 100`,
  central `8·10⁻⁵` (item 7). The binding state's spare margin, `0.044N_s`, is stated. Code: `u7_summary.py`, second part.
* **M2.** I10 is relabelled [I]. The completion requirement is stated for `C⁻₃` and the higher `C_k`, with the leverages (item 6, §3).
* **m1.** The factor `π` in the `|ĝ₃|²` leverage row is fixed (`u5_light.py`). `k_AdS − k(0)`: 17 % / 0.6 % / 0.006 % at `n = 1`.
* **m2.** The residues are re-estimated at `n = 1.3–1.6`, with `n²` scaling, `d·36` coefficients and `c_g = 1, 30`. `Δ_g = 50` is called
  marginal.
* **m3.** The near-gap `κ` statements are corrected. The `Δ_g = 1000` `++` row is flagged as a cut artefact, and the negative `κ` are
  reported.
* **m4.** `u3` is declared a code-consistency check, not evidence for I1, and its non-monotonicity is stated.
* **m5.** "Exact skeleton moments" → truncated order-1/m⁸ series.
* **m6.** The onset of uncontrolled nodes is `x ≈ 0.43` (`k = 7`) to `0.6` (`k = 2`).

