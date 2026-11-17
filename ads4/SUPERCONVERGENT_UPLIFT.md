# `Ω^{SC}`: superconvergent replacements for the moments `C_k`, `k = 3…7`, of `Ω^U`

Status: written 17 Nov 2026 by the `LIGHT_INPUTS`/`OMEGA_PRIME` agent. **Not refereed.** Goal: replace input I10 of `OMEGA_U.md` §3. I10 assumes that
"each project `C_k` is the leading moment of a spin-`k` dispersive functional". The replacement is an explicit construction. This document generalises
`OMEGA_PRIME.md` (`k = 3`) to all `k`. Code: `ads4/code/superconv/` (`bash run_all.sh`, single core, ≈ 10 min; logs in `superconv/cache/`, git-ignored). It imports
`ads4/code/e2e/laplace.py` read-only.

Marks:
* **[D]** derived; **[D-sketch]** outlined;
* **[DV]** derived and checked exactly (sympy) on the `laplace.py` `M`-tables;
* **[Q]** quoted: KKSDZ = 1904.05905, CLS-D = 2311.04271;
* **[I]** inferred; **[ND]** not done.

Premises P1–P3 of `OMEGA_PRIME.md` apply to every [DV] about heavy actions:
* P1: the light-ray mode map;
* P2: the leading moment index;
* P3: the two-pair heavy formula, conditional for stress tensors.

## 0. Results

1. **Construction [D-sketch].** `Ψ_{k,ν} := ∫ A_k · K^{(13)}_{j′,ν} ⊗ K^{(24)} · ⟨T(x₄,z₄)[L[T](x₀,z₃), L[T](x₀,z₁)]T(x₂,z₂)⟩`, with:
   * `K^{(13)}_{j′,ν}` the covariant celestial structure, with `j′ ≡ k (mod 2)`;
   * `K^{(24)}` a free covariant `(2,4)` kernel (`OMEGA_PRIME.md` §7), taken to be `𝔄 = D₅D₁` for odd `k` and the flip-free `𝔖_ff = D₁² + D₁D₂` for even `k`;
   * a conformally invariant weight `A_k = r^{k−3}` (odd `k`, symmetric under `1↔3`) or `A_k = sgn(v′−u′) r^{k−3}` (even `k`, antisymmetric). Here `r` is the CHMRS
     `u`-channel radial cross-ratio.

   Properties:
   * **Valid sum rule.** The commutator of two operators on a common null plane vanishes pointwise. So `Ψ_k[G] = 0` for every planar correlator for which the integral
     converges, whatever the kernel.
   * **Convergence.** Near the Regge limit the kernel is `O(r^{k−2})`, so the integral converges iff `J₀ < k` (`LIGHT_INPUTS.md` §2.3). Elsewhere `A_k` is bounded (`r ≤ 1` on the
     domain, [I]), so convergence there is that of the bare commutator.
   * **Leading Regge kernel.** It is `r^{k−2} × (antipodal weights w^{(j′)}_q(ν)) × F^ν_M × K^{(24)}`, a **pure `ν`-harmonic**.
   * **Why not CLS-D's own factors.** CLS-D's meromorphic subtraction factors `f^±_{k;2}` give leading kernels `∝ η^{j}`, with `j = 1, 0, 3, 2, 1, 0` for `k = 2…7`.
     Only `k = 3, 7` are pure. The factor `η^j` mixes `ν` with `ν ± i, …`, so the heavy action is not `∝ P_{−½+iν}(cosh β)`.
     * Any `1↔3`-antisymmetric meromorphic factor is odd in `η` [D]: `η = √v′ − √u′`, and symmetric functions are functions of `(r, η²)`.
     * Pure harmonics for even `k` therefore need the non-meromorphic `sgn(v′−u′)`; for odd `k` one may use `(u′+v′)^b/(u′v′)^c`, e.g. `(u′+v′)/(u′v′) ≈ 32r²` at `k = 5` [D].
     * Non-analytic weights are harmless for the Einstein statement (item 4). They may change which `[TT]` families enter, which affects (P), not (E).
2. **Odd `k = 3, 5, 7` (`C⁻₃, C⁻₅, C⁻₇`; forward `C⁻₅|₀`) — complete [DV given P1–P3].**
   * The leading heavy action is `N_k(ν) × (2P_{−½+iν}(cosh β)/m^{2k})(|c_same|² − |c_opp|²)`, exactly the project form, with
     * `N₃ = −(32/7)(2ν+i)(4ν²+9)/(4ν²+201)`;
     * `N₅ = −(32/9)(2ν+i)(4ν²+41)/(4ν²+393)`;
     * `N₇ = −(32/11)(2ν+i)(4ν²+89)/(4ν²+649)`.
   * **No real zeros**: `|N_k| ≥ 0.20, 0.37, 0.40` on ℝ. `N₅(0) ≠ 0`, so the forward piece `C⁻_{5,0}` is defined.
   * **No flip interference.** This is structural, for all `ν`, `a`, `l`: the antipodal odd weights at `(1,3)` and `𝔄`'s `z₂↔z₄` antisymmetry at `(2,4)` both kill it (referee
     argument, `REFEREE_OMEGA_PRIME.md`).
   * **Real and even in `ν`.** The shadow factor is `δ/(1−δ)`, which `N_k` shares, and `N_k(−ν) = −N̄_k(ν)` (`OMEGA_PRIME.md` §8; `k`-independent).
   * So the flat limit and the heavy analysis of `OMEGA_U.md` carry over **unchanged** for the odd terms.
3. **Even `k = 4, 6` (`C⁺₄, C⁺₆`; forward `∂_{ν²}C⁺₄|₀`) — preserved part matched, flip interference not removable [DV].**
   * **Preserved–preserved part.** It is `N_k(ν)(2P/m^{2k})(|c_s|² + |c_o|²)`, with
     * `N₄ = 16(4ν²+81)(4ν²+393)(16ν⁴+216ν²+263)/((4ν²+137)²(4ν²+561)²)`;
     * `N₆ = 16(4ν²+169)(4ν²+649)(16ν⁴+536ν²+2023)/(16ν⁴+4712ν²+251161)²`.

     Both are real, even and **positive on ℝ** (minima `0.023`, `0.056` at `ν = 0`). They are defined relative to `D₁²⊗D₁²`; with `𝔖_ff` normalisation they differ by a positive
     factor.
   * **Flip interference.** The covariant `j′ = 0` celestial structure at `(1,3)` has a nonzero helicity-flip response. Its ratio to the preserved response is fixed by
     covariance (F-normalised harmonics):

     | `ν` | 0 | 0.5 | 1 | 2 | 3 | 5 | 10 | 20 |
     |---|---|---|---|---|---|---|---|---|
     | `k = 4` | 0.19 | 0.21 | 0.20 | 0.099 | 0.046 | 0.012 | 1.1·10⁻³ | 7.9·10⁻⁵ |
     | `k = 6` | 0.22 | 0.28 | 0.34 | 0.28 | 0.18 | 0.066 | 8.4·10⁻³ | 6.7·10⁻⁴ |

   * **Consequence.** The leading heavy action contains a sign-indefinite term `∝ ρ_k(ν) c_same c̄_opp 𝒦^ν_{±4,0}(β)`, bounded by `½|ρ_k|(|c_s|²+|c_o|²)|𝒦_{4,0}|`. It
     vanishes in the flat limit (`ρ_k ~ ν⁻⁴`), so the flat limit is unchanged. But it is **not** in the `OMEGA_U` heavy analysis. In particular the forward term
     `∂_{ν²}C⁺₄|₀` sits at `ν = 0`, where `ρ₄ = 0.19`, but with the harmonic `F^0_4(β) = O(β⁴)` at small `β`. Its size in heavy normal form is **[ND]**.
   * **Why it cannot be removed within this family [D].**
     * The `(1,3)` structure is fixed by covariance (`j′ = 0` is unique).
     * The `(flip at 13) × (pres at 24)` term is `∝ A₁₃^flip A₂₄^pres`, and `A₂₄^pres ≠ 0` is needed.
     * `j′ = 1` at `(1,3)` would give the `σ_h` sign, i.e. `C⁻` rather than `C⁺`.
   * **Possible remedy [ND].** Use `z`-derivatives of the light-ray operators, which give transverse components at `(1,3)`. Each derivative costs one unit of Regge
     convergence, which `r^{k−3}` repays for `k ≥ 4`. Not worked out.
4. **Light side [D given RB₂, G^E microcausal, swappability].**
   * For every `k ≥ 3`, `Ψ_k[G^E] = 0`, because the tree Einstein correlator is microcausal [I] and has `J₀(G^E) = 2 < k` (RB₂, [I]/KKSDZ [Q]). So the Einstein part of
     `Ψ_k|_light` vanishes exactly, **with no assumption about the completion of the project's `C_k`**. I10 reduces to RB₂ (plus I14 swappability).
   * By the same rule, `|ĝ₃|²` terms vanish in `k = 5, 6, 7` if `J₀(G^{|ĝ₃|²}) ≤ 4` [I].
   * `k = 3` is stronger. The bare commutator's light side is the single `T` block, which is Einstein-free pointwise (`OMEGA_PRIME.md`), so not even RB₂ is needed.
   * For `k ≥ 4` the `T`-block argument fails, because `A_k` sits inside the null integrals (weighted null integrals are not `E(n)`).
5. **`C⁺₂` stays** (graviton pole; the convergence `J₀ < 2` is marginal). Reference only: a `j′ = 0` replacement would have `N₂` with a real zero at `ν ≈ 0.84`.
6. **Side finding about `OMEGA_U` (important) [DV].** The project's even kernel `D₁⊗D₁` (`D₁²` per pair for `T`) has a helicity-flip response whose ratio to the preserved one
   tends to **1** as `ν → ∞` (F-normalised): `0.34`, `0.70`, `0.94` at `ν = 10, 20, 50` for `k = 4`.
   * The full leading heavy action of `C⁺_k` with this kernel therefore contains `O(1)` `c_same c̄_opp` interference in the flat regime `ν = Λp`. Neither the flat dictionary I1
     nor the `OMEGA_U` heavy scan includes it. The end-to-end check `SPINNING_TTTT.md` (`ν = 0`) cannot see it, since the ratio is `4·10⁻⁴` there.
   * Remedy [D]: use the flip-free even kernel `𝔖_ff = D₁² + D₁D₂ = D₁(z₁·z₃)` at both pairs. Its flip responses cancel exactly (Lemma L table), and its preserved response is
     `(N₁−N₂)/(64𝖽𝖽) → ½`. This is within the project's stated freedom ("any combination of `D₁², D₁D₂, D₂²`", `REGGE_MOMENT_TTTT.md` §4).
   * `Ω^{SC}` uses `𝔖_ff` at `(2,4)`. At `(1,3)` the residual `ρ_k` of item 3 decays like `ν⁻⁴` instead.

**`Ω^{SC}`** = `Ω^U` with:
* `C⁻_{3,5,7} → Re(Ψ_k/N_k)` (= `Ψ_k/N_k`);
* `C⁺_{4,6} → Ψ_k/N_k` (with the residual flip term of item 3);
* the forward pieces `∂_{ν²}C⁺_{4}|₀ → ∂_{ν²}(Ψ₄/N₄)|₀` and `C⁻_{5,0} → (Ψ₅/N₅)|₀`;
* `C⁺₂` unchanged.

Status of its light side: no Einstein term in `k = 3…7` [D given RB₂ (k ≥ 4); [D] for `k = 3`]. Its heavy side equals `Ω^U`'s for odd `k` [DV given P1–P3], and up to the
`ρ_k` interference for even `k` ([ND] in normal form).

## 1. Method (`s0_check_D1sq.py`, `s1_responses.py`, `s2_even.py`)

* **Pair responses.** For a pair, `A(h₁,h₃) = Σ_q w_q M^q_{h₁h₃}(ν; a)`, with `a = 2l+k+1 = k+5` and exact `M = R + g₁S` from `laplace.py`. Preserved means `h₃ = −h₁` (harmonic `P`);
  flip means `h₃ = h₁` (harmonic `F^ν_4`).
* **Antipodal weights** (P1):
  * `j′ = 0`: `w_q = (−1)^{q₃}(δ/2)_n/(1−δ/2)_n`, with `n = |q₁+q₃|/2`;
  * `j′ = 1`: as in `OMEGA_PRIME.md` §6.
* **Project kernels in contraction weights.**
  * `D₁² ≅ 4[(2,−2)+(−2,2)] + (2/9)(0,0)`. `s0` checks it against the Lemma L table: the preserved response is exactly `N₁/(64𝖽𝖽)`; the flip agrees up to the `P^{(4)}` vs `F₄`
    normalisation.
  * `𝔄 ≅ (2,−2) − (−2,2)`, which matches the tabulated shape.
* **Reality and evenness [DV].**
  * `j′ = 0`: the shadow factor is 1 and `w_q(−ν) = w̄_{−q}(ν)`; with real `A_k`, `𝔖_ff` and `N_k`, `Ψ` is real and even.
  * `j′ = 1`: as in `OMEGA_PRIME.md` §8.
* **Normalisation.** `N_k` is the ratio of the preserved pair responses. For odd `k` that is the `(1,3)` ratio, since `(2,4)` is the project's own `𝔄`. For even `k`,
  `s1` prints the `(1,3)`-only ratio (j′=0 vs `D₁²`), and `s2` the full two-pair `N_k` of item 3 (including `𝔖_ff` vs `D₁²` at `(2,4)`). In both cases it is so that the action equals `(2P/m^{2k})(|c_s|² ± |c_o|²)`. The light-ray Jacobian `j(ν)` remains untracked
  (`OMEGA_PRIME.md` item 5, [I]).

## 2. What remains

| item | status |
|---|---|
| `Ψ_k` valid spin-`k` sum rules (`A_k` weights; domain `r ≤ 1`) | [D-sketch]; domain [I] |
| odd `k`: heavy action `= N_k·`project, flip-free, real/even, `N_k ≠ 0` on ℝ | [DV given P1–P3] |
| even `k`: preserved part `= N_k·`project, `N_k > 0` on ℝ; residual flip `ρ_k` (≤ 0.34, `~ν⁻⁴`) | [DV]; effect in heavy normal form [ND]; remedy via `z`-derivatives [ND] |
| Einstein-free light side, `k = 3…7` | [D given RB₂, microcausality of `G^E`, swappability] |
| project `D₁⊗D₁` kernels carry `O(1)` flip interference at large `ν` | [DV] (`s2`); affects `OMEGA_U` I1/I2 |
| `j(ν)`; `O(1/m²)` heavy terms; contact `(P)` for non-meromorphic `A_k` | [I] / [ND] |
