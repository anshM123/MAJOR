# `Ω^{NF,SC}`: forward-free superconvergent uplift of NF02 — heavy scan, light side, conditional bound

Status: written 18 Nov 2026 by the `OMEGA_U` agent. **Refereed once (`REFEREE_OMEGA_NF.md`); revised (list in §3). Not a theorem about AdS4.** It applies the machinery of `OMEGA_U.md`
(revised) to the forward-free functional NF02 of `NOFORWARD_FUNCTIONAL.md`, with the superconvergent replacements of
`SUPERCONVERGENT_UPLIFT.md`.

Code: `ads4/code/omegau/nf_*.py`, `nf1`, `nf4`–`nf7` (run via the tail of `run_all.sh`; single core, ≈ 3 min given the `u4` skeleton
caches). Logs: `omegau/cache/nf{1,5,6}.log`, `cache/nf7_summary.txt` (git-ignored). Nothing outside `ads4/code/omegau/` and this file
was modified.

Marks as in `OMEGA_U.md`: **[D]**, **[DV]**, **[NE]**, **[M]**, **[I]**, **[A]**. Nothing is fitted.

**Definition.** `Ω^{NF,SC} = ∫₀¹dp Σ_{k=2}^{7} φ_k(p) ε_kΛ^{2k} Ĉ_{k,Λp}`. There are no forward terms (`a = 0`). Here:
* `φ_k` are built from NF02's `ψ₂`, `ψ₃` exactly as in `FLAT_UNSTRIPPED.md` item 5;
* `Ĉ₂` is `C⁺₂` with the flip-free kernel `𝔖_ff⊗𝔖_ff`;
* `Ĉ_{3,5,7} = Ψ_k/N_k` and `Ĉ_{4,6} = Ψ_k/N^{ff}_k` (`SUPERCONVERGENT_UPLIFT.md`, end of §0).

At leading heavy order the normalised replacements equal the project moments (ibid. items 2, 3, 6). So the heavy model is that of
`adsaction.py` with `a = 0`, plus a residual one-pair flip in `k = 4, 6`, bounded in item 2.

**Numerics.** NF02's expanded coefficients reach `6·10¹⁰`, with strong cancellation. All weights are therefore evaluated in factored form
`(1−p)⁵pX(p)`, and the Mehler–Dirichlet and `J₀` transforms use exact rational moments and endpoint derivatives (`nf_common.ExactCosT`).

## 0. Results

1. **Flat limit [D + DV; `nf1`].**
   * The B′ parts sum to the stripped NF02 action to `≤ 10⁻¹³N_s`. This is the pointwise identity of `OMEGA_U` item 1, for generic `ψ`.
   * The flat minimum is `F/N_s = 0.2000` at `m = M` (`J = 4` for `++`, `J = 5` for `+−`), as in `NOFORWARD_FUNCTIONAL.md`.
   * **No cancellation amplification.** Because NF02's weights vanish like `(1−p)⁵`, `N_U/N_s ≤ 1.53` (`++`) and `≤ 2.2` (`+−`) on every
     scanned state, including `m = Λ`. For `Ω^U` (M02) the ratio was 62–519 near `m = Λ`.
   * The odd-`k` prefactor sensitivity (I3) is `≤ 0.12N_s` (`+−` at `m = Λ`), `≤ 0.009N_s` (`++`, `m ≥ 1.05Λ`), and `≤ 0.003N_s` at
     `≥ 1.4Λ`. For M02 it was up to `1.6N_s` at `1.05Λ`.
2. **Heavy scan [NE + M + I; `nf4`, `nf7`].**
   * **Scope.** The states, nodes, skeletons (`Δφ = 3, 2`) and per-`k` skeleton remainders are those of `OMEGA_U` (re-using its caches).
   * **Near the gap and at small impact parameter, much better than `Ω^U`.**
     * With `b̂ = (2J+1)Λ/m ≤ 5`, every state is positive with `κ ≥ 1` down to **`m = Λ`**, for all `Δ_g = 20 … 1000`.
     * With `b̂ ≤ 10`, the threshold is `n(κ_req = 1) = 1.35–1.45` (`Δ_g ≥ 50`; 2.15 at `Δ_g = 20`), and `n(0.3, δ = 0.03) = 1.13–1.14`.
     * The controlled margins there are `≥ 0.33N_s`.
   * **At large impact parameter the skeleton method fails.**
     * The failing states have `b̂ ≥ 10`: 603 / 10636 states in `[10, 20)`, 2853 / 8428 in `[20, 30)`, 3376 / 5998 in `[30, 40]`. None fails
       below `b̂ = 10`.
     * There NF02's action is the forward-free tail `F ≈ −[2ψ₂″(0)+σuψ₃″(0)]/(2b̂³)`.
       * At `m = 1.5Λ`, `Δ_g = 100`, the leading Regge action is `0.078 / 0.022 / 0.005` (`++`) and `0.22 / 0.033 / 0.006` (`+−`) of the
         pointwise coherent norm `∫|Σ_k u^{k−2}2ε_kσ^kφ_k||P|` at `b̂ ≈ 10 / 20 / 40` (`nf6`). Against the per-`k` pointwise norm it is
         similar.
       * This sits below the per-node error of the truncated series: the unresummed `(𝒥x²)ⁿ/m²` cross terms, and the fixed-η fallback's
         bulk-point truncation.
       * That these are expansion artefacts rather than AdS physics is plausible but **not shown**: no resummed or independent computation
         was done. If they are artefacts, there is no evidence for the size of the true correction there at all.
     * Even the signed ("coh") criterion fails there, because `S` carries the same truncation artefacts.
     * For M02 this region was dominated by the forward term `2aJJu³`, so the issue never arose.
     * With `b̂ ≤ 40`, no `n ≤ 3` is supported at `κ_req ≥ 0.3`, for any `Δ_g`. With `b̂ ≤ 20`, `n(1) = 2.75–2.9` and `n(0.3, 0.03) = 2.0–2.1`.
   * **Large `b̂`: leading Regge action [NE] and envelopes (revised after referee M1; `nf6`).**
     * `FR = Σ_k u^{k−2}2ε_kσ^k∫φ_kP_{−½+iΛp}(cosh β)` is positive on every physical grid point (`β ≤ 8`, `m ≤ 10³Λ`), for `Δ_g = 20, 50, 100, 1000`.
     * With an envelope `FR − C₁uN_x > 0`, the conclusion depends entirely on the norm `N_x`:
       * **Component norm** `N_comp = Σ_k|∫…P|`. This shares NF02's `1/b̂³` cancellation, so the envelope then assumes that the correction
         density cancels like the leading action. Nothing argues this. On `b̂ > 20` it tolerates `C₁ ≤ 1.9 / 4.1 / 8.7` at `n = 1.5 / 2.1 / 3`.
       * **Pointwise norms** (coherent or per `k`). The envelope gives no uniform `n`: the last negative point is `34Λ` (`Δ_g = 100`, `C₁ = 1`)
         and `213Λ` (1000), and grows with `Λ`. On `b̂ > 20` it tolerates a relative pointwise correction of only `≤ 4·10⁻³` (`Δ_g = 100`) or
         `≤ 9·10⁻⁵` (1000) at `n = 2.1`.
     * **The skeleton's own data imply `C₁ ≫ 1` where it fails** (referee):
       * on `b̂ ∈ (15, 20]` its threshold of 2.0–2.1 corresponds to `C₁ ≈ 4`, and on `(20, 30]` its 3.0 to `C₁ > 5`;
       * the `C₁_eff` measured in `REFEREE_UNIFORM_REMAINDER` (17 / 3.2 / 0.9 at `Δ_g = 20 / 50 / 100`, for M02, whose `N` was
         forward-dominated) would give `n ≈ 4.4 / 2.2 / < 1.5`.

       So `C₁ = 1` is **not** calibrated, and the earlier "central `n ≈ 1.5`" is **withdrawn**.
   * **What is supported (skeleton only, by `Δ_g`; `nf7`, `b̂ ≤ 20`).**

     | `Δ_g` | central (`κ_req = 0.3`, `δ = 0.03`) | conservative (`κ_req = 1`) |
     |---|---|---|
     | 20 | 3.0 (`abs`) / 2.6 (`coh`) | 3.0 |
     | 50 | 2.05 | 2.75 |
     | 100 | 2.00 | 2.80 |
     | 300 | 2.10 | 2.90 |
     | 1000 | 2.10 | 2.90 |

     * So the supported threshold is **`n ≳ 2.0–2.1` (central) to `2.75–2.9` (conservative), for `Δ_g ≥ 50`, for the states with `b̂ ≤ 20` only.**
     * `δ = 0.03` covers the I3 odd-prefactor shift for `m ≥ 1.05Λ` (`++`) and `m ≥ 1.1Λ` (`+−`). It reaches `0.12N_s` only at `+−`, `m = Λ`,
       far below these thresholds.
     * **`b̂ > 20` is open [I15, I].** Positivity there would need one of the following:
       * a relative AdS correction, in component norm, of `C₁ ≤ 4` (`n = 2.1`) or `≤ 9` (`n = 3`), i.e. a correction density with the same
         `1/b̂³` cancellation as the leading action;
       * or a pointwise relative correction `≲ 4·10⁻³` (`Δ_g = 100`);
       * or, best, a fixed-`b̂` (eikonal-type) resummation of the `(𝒥x²)ⁿ` series of the skeleton, so that the remainder is computed rather
         than bounded. **None of these exists.**
     * `n = 1` holds only for `b̂ ≤ 5`. Those near-gap states are `++` with `J ≤ 2`; every `+−` state with `m < 1.8Λ` has `b̂ > 5`.
   * **Binding states.**
     * Below `b̂ = 10`: `m ≈ 1.0–1.1Λ`, with `κ = 0.02–0.09`.
     * For the `b̂ ≤ 20` thresholds: states with `b̂ ∈ (15, 20]` at `m ≈ 2–2.9Λ`.
   * **Residual flip** of `Ψ₄/N^{ff}₄`, `Ψ₆/N^{ff}₆` with NF02 weights (the `superconv/c2` model): `≤ 1.7·10⁻⁶N_s` (`Δ_g = 100`) and
     `≤ 5·10⁻¹⁰N_s` (`Δ_g = 1000`). There is no forward flip piece (`nf5` (d)).
3. **Light side [D given (F), (P), (R)/IP, RB₂; `nf5`].**
   * **Flat.** All contacts cancel in `Σφ_kB′_k`, with no forward terms. Verified (normalised, `k(0) = 1`; unnormalised
     `k(0) = 88865643959/3233230000`):

     | constant | value |
     |---|---|
     | `A` | **3022.3377** |
     | `B` | −5327.2096 |
     | `B_AdS` | **−6552.6620** |

   * **Graviton pole.** It sits only in `Ĉ₂`, weight `ψ₂`, residue `16πG/(ν²+9/4)` with `𝔖_ff` [D given IP, (F)].
   * **`k ≥ 3`.** The Einstein term vanishes: `Ψ_k[G^E] = 0` given RB₂, microcausality of `G^E` and I14; for `k = 3`, with no RB₂ at all.
   * **Contacts.** In AdS the cancellation is again not exact. The leading `Λ⁻²` leverages `Σ|lev|` are:

     | coupling | NF02 (no forward) | `Ω^U`, M02 |
     |---|---|---|
     | `g₄` | 19.6 | 0.64 |
     | `g₅` | 11.9 | 2.2 |
     | `g₆` | 6.9 | 7.0 (3.9 of it forward) |
     | `g₆′` | 2.1 | 2.3 |
     | `|ĝ₃|²` | 8.4 | 1.6 |

     * The smeared leverages are larger: NF02's weights are larger relative to `k(0)`.
     * But relative to NF02's much larger `16L_AdS/k` (`5.5·10³` at `Δ_g = 100`, `n = 2.1`), the contact shift is `2·10⁻⁴` (`c_g = 1`) or
       `7·10⁻³` (`c_g = 30`) at `n = 2.1`, and `5·10⁻⁴ / 1.4·10⁻²` at `n = 2.9` (`Δ_g = 100`). It is `≤ 10⁻⁴` at `Δ_g = 1000`.
     * `k_AdS − k(0)` (`|a| = d·36`) is `19 %` (`n = 2.1`) to `36 %` (`n = 2.9`) at `Δ_g = 100`, and `0.2–0.35 %` at 1000. **At
       `Δ_g = 100` this unknown `O(Λ⁻²)` correction to the `|ĝ₃|²` coefficient is no longer small.**
   * **Leverage reduction.** Take a test residue `16πG(ν²+9/4)^{−k}/N_k(ν)` in each replaced term — what a failure of (P) for the
     replacements would cost — relative to each functional's own leading `16πGΛ²L_AdS`:

     | `Δ_g` | `Ω^{NF,SC}`, `k = 3…7` | `Ω^U` (M02, same normalisation), `k = 3…7` | `Ω^U` forward `C₄`, `C₅` |
     |---|---|---|---|
     | 100 | 48, 200, 5, 16, 0.4 | `2.1·10⁴, 3.2·10⁴, 540, 1.2·10³, 23` | `4.6·10¹⁴`, `1.8·10¹³` |
     | 1000 | `2.6·10³, 1.1·10⁴, 270, 870, 20` | `1.4·10⁷, 2.3·10⁷, 4·10⁵, 10⁶, 2·10⁴` | `2.9·10²²`, `1.1·10²¹` |

     * The forward pieces, `10¹³–10²²`, are gone.
     * **Attribution (referee m1).** The comparison is per each functional's own bound numerator.
       * Most of the smeared reduction is the factor `L_AdS(NF02)/L_AdS(M02) = 43 / 49 / 53` (`Δ_g = 20 / 100 / 1000`). That is the same large
         `A` that costs the factor `≈ 7` in the bound.
       * The weight shape adds `9 / 3.2 / 2.2 / 1.4 / 1.3` for `k = 3…7` at `Δ_g = 100`, and `100 / 40 / 30 / 22 / 20` at `Δ_g = 1000`.
       * The `k = 3` shape gain is due to `ψ₃(0) = 0`. For `k = 4` it does not apply: `φ₄(0) = 0` for both functionals.
     * The smeared leverages are `48–200` only at `Δ_g = 100`; at `Δ_g = 1000` they are `2.6·10³–1.1·10⁴`.
     * The largest remaining leverage is `k = 4`. Its weight is non-meromorphic, so (P) is unanalysed there.
4. **Conditional bound** `√(t₄²+α_T²) ≤ 240(n⁸·16L_AdS(Λ/n)/k)^{1/2}/Λ⁴`, at the `b̂ ≤ 20`-supported `n` (item 2). It is valid **only if**
   `b̂ > 20` is positive (I15 [I]), together with all of §1:

   | `Δ_g` | `Ω^{NF,SC}`, central `n` | `Ω^{NF,SC}`, conservative `n` | `Ω^U` (M02), central / conservative (`OMEGA_U.md` item 7) |
   |---|---|---|---|
   | 50 | `0.047` (`n = 2.05`) | `0.14` (2.75) | `1.5·10⁻³` / `1.9·10⁻³` |
   | 100 | `3.1·10⁻³` (2.00) | `1.1·10⁻²` (2.80) | `8.1·10⁻⁵` / `1.1·10⁻⁴` |
   | 300 | `5.4·10⁻⁵` (2.10) | `1.9·10⁻⁴` (2.90) | `1.1·10⁻⁶` / `1.6·10⁻⁶` |
   | 1000 | `5.2·10⁻⁷` (2.10) | `1.8·10⁻⁶` (2.90) | `9.2·10⁻⁹` / `1.3·10⁻⁸` |

   * `Δ_g = 20` (`n ≥ 2.6`) gives `√ ≥ 3.5`, i.e. nothing.
   * At `Δ_g = 100`, `k_AdS − k(0) ≈ 20–36 %` is an additional unquantified uncertainty (item 3).
   * The cost of dropping the forward terms is now a factor `≈ 40–140` in `√(t₄²+α_T²)` at matched conventions: `≈ 7` from `A`, the rest
     from `n ≈ 2–2.9` instead of `1.3–1.4`. Against this, the `Ω^U` thresholds rest on the forward term at large `b̂`, with the light-side
     forward leverage `10¹³⁺`.
   * The withdrawn version quoted `1.0·10⁻³` (`n = 1.5`), `3.7·10⁻³` (`n = 2.1`) and `2.2·10⁻⁴` (`n = 1`) at `Δ_g = 100`. The `n = 1` row
     was not a bound and is dropped.
   * What `Ω^{NF,SC}` buys is the removal of the `10¹³⁺` forward leverage, and with it the most dangerous light-side inputs (§1).

## 1. Inputs: changes relative to `OMEGA_U.md` §3

| `OMEGA_U` input | status for `Ω^{NF,SC}` |
|---|---|
| I1 B′ dictionary / reading B | unchanged [I] |
| I2 `⟨TTTT⟩` heavy formula | unchanged ([DVE] only at `k = 6`, `ν = 0`); the forward `k = 4, 5` at `ν = 0` are no longer needed |
| I3 odd-`k` prefactor | **weakened**: odd moments are the explicit `Ψ_k/N_k`, whose leading heavy action is [DV given P1(j′=1)–P3]; the prefactor sensitivity is `≤ 0.12N_s` (`≤ 0.003N_s` at `≥ 1.4Λ`) instead of `O(N_s)` |
| I4, I5 skeleton transfer, odd `k` by continuation | unchanged [I] |
| I6 control cut, `κ_req`, `δ` | unchanged [M]; **plus the new I15** |
| I7 spin-specific `1/m²` | unchanged [I] |
| I8 (F), (P) | (F) unchanged; (P) is now needed only for smeared moments (no `ν = 0` pieces); still [I] and unanalysed for even `k` (non-meromorphic `A_k`), leverage ≤ `2·10²` (`Δ_g = 100`) instead of `10¹³` |
| I9 (R) | unchanged for `C⁺₂`, re-derived for `𝔖_ff` [D given IP, (F)] |
| **I10** | **Einstein part reduced** to RB₂ + microcausality of `G^E` + I14 for `k = 4…7`, and to nothing beyond the T-block argument for `k = 3`. The "completion of the project `C_k`" assumption is removed (`Ψ_k` are explicit sum rules). **`|ĝ₃|²` absence in `k = 5, 6, 7` remains [I]** (needs reach `J₀ ≤ 4`). The forward leverage is removed |
| I11 contacts, power counting | unchanged [A]; the shift is negligible at `Δ_g ≥ 100` |
| I12 light spectrum, tree level, double traces | unchanged; the double traces of even `Ψ_k` are unanalysed [I] |
| I13 parity-odd heavy | not treated |
| I14 swappability | unchanged [I], and now also used in the Einstein argument |
| **new** I15 | **[I]**: positivity of the states with `b̂ > 20` (`b̂ ∈ (10, 20]` is skeleton-controlled at `n ≳ 2`). The leading Regge action is positive there, but no calibrated envelope exists: it needs `C₁ ≤ 4` (component norm, `n = 2.1`) or a pointwise relative correction `≲ 4·10⁻³` (`Δ_g = 100`), against skeleton-implied `C₁ ≳ 4–5` and `C₁_eff` up to 17 |
| **new** I16 | flat NF02 positivity is not fully certified: `+−`, `J ≥ 400`, `M ≤ m < 1.633M` is float-only (`NOFORWARD_FUNCTIONAL.md`) |
| **new** I17 | P1(j′=0) and the sign property (`s3`) for even `Ψ_k`; light-ray Jacobian `j_k(ν)` (possibly `k`-dependent) untracked [I] |
| removed | the `D₁⊗D₁` flip contamination of `OMEGA_U` (`SUPERCONVERGENT_UPLIFT.md` item 6; up to `0.4N_s` there) — `𝔖_ff` is flip-free; the residual flip is `≤ 2·10⁻⁶N_s` |

## 2. Code (`ads4/code/omegau/`)

| file | content | log |
|---|---|---|
| `nf_common.py` | NF02 (exact, normalised), factored `ψ`, `φ_k`, `ExactCosT`, flat parts, per-state action (`a = 0`) | — |
| `nf1_flat.py` | factored vs exact weights, flat identity, `N_U/N_s`, odd-prefactor sensitivity | `cache/nf1.log` |
| `nf4_scan.py` | heavy scan (reuses `skel_u4_*` caches); `OMEGAU_DF = 3, 2` | `cache/recs_nf4_*.pkl` |
| `nf7_summary.py` | bands, `n(κ_req, δ)` (`abs` / `coh`), binding states, failures by `b̂`, `b̂`-restricted thresholds | `cache/nf7_summary.txt` |
| `nf6_regge.py` | leading Regge action; envelopes with component, coherent-pointwise and per-`k` pointwise norms; tolerated `C₁` on `b̂ > 20` | `cache/nf6.log` |
| `nf5_light.py` | low sides, constants, contact leverages, normalised test residues vs `Ω^U`, residual flip | `cache/nf5.log` |

## 3. Revisions after `REFEREE_OMEGA_NF.md`

* **M1.** The `C₁ = 1` large-`b̂` envelope is withdrawn as uncalibrated. `nf6` now uses the component, coherent-pointwise and per-`k`
  pointwise norms, and reports the tolerated `C₁`. The supported statement is now `n ≳ 2.0–2.1` (central) / `2.75–2.9` (conservative), for
  `b̂ ≤ 20`. `b̂ > 20` is open. I15 is relabelled [I], and the bound table is recomputed at the supported `n`.
* **M2.** All `n` are given per `Δ_g`. `Δ_g = 20` gives nothing.
* **m1.** The leverage comparison is per own bound. The `L_AdS` ratio (43–53) is separated from the shape factor, and the `Δ_g = 1000`
  values are stated.
* **m2.** The `n = 1` row is dropped. The `b̂ ≤ 5` states are `++`, `J ≤ 2` only.
* **m3.** The I3 shift is covered by `δ = 0.03` at all supported thresholds.
* **m4.** The norm behind the "fraction of integrand norm" statement is specified (coherent pointwise; per-`k` similar).
* **m5.** Contact shifts are given at `n = 2.1` and `2.9`. They remain small, but `k_AdS − k(0)` reaches 20–36 % at `Δ_g = 100`.

