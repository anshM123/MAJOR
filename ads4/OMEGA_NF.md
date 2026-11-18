# `Ω^{NF,SC}`: forward-free superconvergent uplift of NF02 — heavy scan, light side, conditional bound

Status: written 18 Nov 2026 by the `OMEGA_U` agent. **Unrefereed. Not a theorem about AdS4.** It applies the machinery of `OMEGA_U.md`
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
     * There NF02's action is the forward-free tail `F ≈ −[2ψ₂″(0)+σuψ₃″(0)]/(2b̂³)`. This is `10⁻²–10⁻⁵` of the integrand norm, and it
       sits below the per-node error of the truncated series: the unresummed `(𝒥x²)ⁿ/m²` cross terms, and the fixed-η fallback's
       bulk-point truncation. Both are artefacts of the expansion, not AdS physics.
     * Even the signed ("coh") criterion fails there, because `S` carries the same truncation artefacts.
     * For M02 this region was dominated by the forward term `2aJJu³`, so the issue never arose.
     * With `b̂ ≤ 40`, no `n ≤ 3` is supported at `κ_req ≥ 0.3`, for any `Δ_g`. With `b̂ ≤ 20`, `n(1) = 2.75–2.9` and `n(0.3, 0.03) = 2.0–2.1`.
   * **Large `b̂` instead, leading Regge action [NE + M; `nf6`].**
     * `FR = Σ_k u^{k−2}2ε_kσ^k∫φ_kP_{−½+iΛp}(cosh β)` is positive on every physical grid point (`β ≤ 8`, `m ≤ 10³Λ`), `Δ_g = 20, 100, 1000`.
     * Its relative margin at `b̂ ≥ 10` is `FR/N ≥ 0.15` (`++`) and `≥ 0.52` (`+−`).
     * With the envelope `C₁uN` the last negative point is `1.50Λ / 1.29Λ` for `C₁ = 1`, and `2.12Λ / 1.92Λ` for `C₁ = 3`.
     * This is the old "statement R" type of argument [M]. It does not control `O(ν²/m²)` terms beyond the envelope.
   * **Combined threshold [M].**
     * Central: skeleton `κ_req = 0.3`, `δ = 0.03` on `b̂ ≤ 10`, plus the `C₁ = 1` Regge envelope above. This gives **`n ≈ 1.5`**, set by the
       large-`b̂` `++` envelope.
     * Conservative: `κ_req = 1`, plus `C₁ = 3`. This gives **`n ≈ 2.1`**.
     * `n = 1` holds only for `b̂ ≤ 5`.
   * **Binding states.**
     * Below `b̂ = 10`: `m ≈ 1.0–1.1Λ`, with `κ = 0.02–0.09`.
     * Above `b̂ = 10`: the Regge-envelope threshold, `++`, `b̂ ≈ 10–12`, `m ≈ 1.5Λ`.
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
     * But relative to NF02's much larger `16L_AdS/k` (`7.5·10³` at `Δ_g = 100`) the contact shift is `≤ 10⁻⁴` (`c_g = 1`) or `3·10⁻³`
       (`c_g = 30`) at `Δ_g ≥ 100`, `n ≤ 1.5`.
     * `k_AdS − k(0)` is `4–10 %` (`Δ_g = 100`) and `0.04–0.09 %` (1000).
   * **Leverage reduction.** Take a test residue `16πG(ν²+9/4)^{−k}/N_k(ν)` in each replaced term — what a failure of (P) for the
     replacements would cost — relative to the leading `16πGΛ²L_AdS`:

     | `Δ_g` | `Ω^{NF,SC}`, `k = 3…7` | `Ω^U` (M02, same normalisation), `k = 3…7` | `Ω^U` forward `C₄`, `C₅` |
     |---|---|---|---|
     | 100 | 48, 200, 5, 16, 0.4 | `2.1·10⁴, 3.2·10⁴, 540, 1.2·10³, 23` | `4.6·10¹⁴`, `1.8·10¹³` |
     | 1000 | `2.6·10³, 1.1·10⁴, 270, 870, 20` | `1.4·10⁷, 2.3·10⁷, 4·10⁵, 10⁶, 2·10⁴` | `2.9·10²²`, `1.1·10²¹` |

     * The forward pieces, `10¹³–10²²`, are gone.
     * The smeared `k = 3, 4` leverage falls by `≈ 150–450` (`Δ_g = 100`) and `≈ 2000–5000` (`Δ_g = 1000`), because `ψ₃(0) = 0` makes `φ₃(0) = 0`.
     * The largest remaining leverage is `k = 4`: `2·10²` at `Δ_g = 100`. Its weight is non-meromorphic, so (P) is unanalysed there.
4. **Conditional bound** `√(t₄²+α_T²) ≤ 240(n⁸·16L_AdS(Λ/n)/k)^{1/2}/Λ⁴`:

   | | `Δ_g = 100` | `Δ_g = 1000` |
   |---|---|---|
   | `Ω^{NF,SC}`, central `n ≈ 1.5` | `1.0·10⁻³` | `1.4·10⁻⁷` |
   | `Ω^{NF,SC}`, conservative `n ≈ 2.1` | `3.7·10⁻³` | `5.2·10⁻⁷` |
   | `Ω^{NF,SC}`, `n = 1` (holds only for `b̂ ≤ 5`) | `2.2·10⁻⁴` | `2.9·10⁻⁸` |
   | `Ω^U` (M02), central `n = 1.28` / conservative `1.40` (`OMEGA_U.md` item 7) | `8.1·10⁻⁵` / `1.1·10⁻⁴` | `9.2·10⁻⁹` / `1.3·10⁻⁸` |

   * The cost of dropping the forward terms is a factor `≈ 12–40` in `√(t₄²+α_T²)` at matched conventions: `≈ 7` from `A` (as in
     `NOFORWARD_FUNCTIONAL.md`), the rest from the larger `n`.
   * What it buys is the removal of the `10¹³⁺` forward leverage, and with it the most dangerous light-side inputs (§1).

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
| **new** I15 | large-`b̂` (`b̂ ≳ 10`) states: positivity rests on the leading Regge action plus a `C₁uN` envelope [M], not on the skeleton remainder |
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
| `nf6_regge.py` | leading Regge action and envelopes | `cache/nf6.log` |
| `nf5_light.py` | low sides, constants, contact leverages, normalised test residues vs `Ω^U`, residual flip | `cache/nf5.log` |
