# Referee report on `ads4/SUPERCONVERGENT_UPLIFT.md` and `ads4/code/superconv/`

Referee: independent, hostile. Date 17 Nov 2026. No file under review was modified.
* Scratch code: `ads4/code/referee_superconv/` (`r1_direct.py`, `r2_nu20.py`, `r3_forward.py`). Logs are in `cache/`, which is git-ignored.
* `r1` and `r2` do not use `laplace.py`, `s1` or `s2`. They are a direct numerical Fourier–Laplace transform on `H²`, the method of `referee_omegaprime/r1`, re-implemented here.
  * The celestial weights come from quadrature of `|sin ψ|^{−δ}` (`j′ = 0`) and of `sgn(sin ψ)|sin ψ|^{−δ}` (`j′ = 1`).
* The `superconv` logs agree with the text.

## Independent checks

| check | result |
|---|---|
| `N₅(ν = 1.3)`, `a = 10`, `j′ = 1` vs `𝔄` | direct `−1.10444934 − 0.42478821i`, claim identical (rel. `2·10⁻¹¹`); flips `< 10⁻¹⁸` |
| `k = 4`, `a = 9`, `j′ = 0` at `ν = 1` | pres `0.06732404` (claim same); flip `−0.0092442 + 0.0092856i` (claim same to `4·10⁻⁶`) |
| `D₁²` flip/pres, `a = 9` | `ν = 10`: `0.34281` (claim `0.34281`); `ν = 20`: `0.69978` at the finest grid, converging to the claimed `0.69993` |
| `d_k` pattern | `d_k = 8(k+2)² + 1 = 8(a−3)² + 1`. Numerator `4ν² + 2k² − 9`, prefactor `−32/(a−1)`. Both reproduce `OMEGA_PRIME`'s `N₃` exactly, so the pattern is consistent. |
| `j′ = 0` weights | Fourier coefficients of `|sin ψ|^{−δ}` give `(δ/2)_n/(1−δ/2)_n` analytically, as used in the code |
| `𝔖_ff` flip-freedom | structural: the factor `z₁·z₃` contracts `ε_h·ε_h = 0` |
| `c2` forward estimate (`r3`) | the one-sided `ν²` difference underestimates `∂_{ν²}` of the (even) flip term by about 9 % |

## Verdicts

* **Correctness:** correct, independently reproduced. The closed forms, the `|N_k|` minima and the flip ratios are right. The threshold `J₀ < k` is **not** off by
  one: the weight `r^{k−3}` reproduces CLS-D's Regge shift `n = k − 2J_e + 1`.
* **Soundness:** sound for the Einstein claim at leading Regge order, given `r ≤ 1` [I] and I14. There are gaps M2 and M3.
* **Labels:** overstated in places (M3, M4, and "I10 reduces to RB₂", M1).
* **Applicability:** not established for the even and forward pieces (M1), nor for the kernel switch in `C⁺₂` (M4).
* **Novelty:** low to moderate.
  * Weighting the light-ray commutator to raise the superconvergence spin is the KKSDZ/CLS-D mechanism.
  * What is new: the pure-harmonic (non-CLS-D) weights, the closed forms `N_k`, and the flip diagnosis of `D₁⊗D₁` together with its remedy `𝔖_ff`.

## Issues

**Fatal:** none found.

**M1 (applicability).** "I10 reduces to RB₂" overstates the result. It removes the Einstein term only.
* For `k ≥ 4` the weight sits inside the null integrals. The doc concedes that `[TT]` and contact terms (P) are then unanalysed.
* This matters most for the forward pieces, which carry I10's leverage of about `10¹³`. `∂_{ν²}(Ψ₄/N₄)|₀` also divides by `N₄(0) = 0.023`.
* Even the odd-`k` weights `r²` and `r⁴` are algebraic in `(u′, v′)`, not meromorphic. `(u′+v′)/(u′v′) ≈ 32r²` holds only at leading order.

**M2 (soundness, claim 1).** "Kills anything with `J₀ < k`" is shown only for the `r → 0` region.
* Elsewhere the doc defers to the bare commutator. A single `L[T]` converges only for `J₀ < 3`.
* For `3 ≤ J₀ < k` the single-ray endpoints are not shown to be suppressed by `r^{k−3}`. This covers the `|ĝ₃|²` claim with `J₀ ≤ 4` in `k = 5, 6, 7`.
* `J₀ = 2`, i.e. Einstein, is unaffected.

**M3 (soundness/label, even `k`).** The code never implements `sgn(v′−u′)`.
* The code assumes that at leading Regge order the sign is constant on the `y > 0` domain of the `M`-tables, i.e. that it labels the sheet. That is plausible, but it is not derived.
* If the sign varied with the celestial angle, `j′ = 0 × sgn` would be the `j′ = 1` (`C⁻`) structure.
* P1 was stated and refereed for `j′ = 1` only.
* So the even-`k` [DV] labels should be conditional on P1(`j′ = 0`) and on the sheet-label property.

**M4 (claim 4).**
* `⟨JJJJ⟩` is `l = 1`, so it cannot test the `l = 2` kernel `𝔖_ff`.
* The `⟨TTTT⟩` test is `k = 6`, `ν = 0`, `β = 1`, at percent level. There `D₁²` flip/pres is `4·10⁻⁴`. So "the tests cover `𝔖_ff`" holds only at `ν = 0`, not where the kernel choice matters (`ν ≳ 5`).
* "The scan is not valid for `D₁⊗D₁`" rests on a model bound with unfixed `O(1)` factors. It should read "not shown valid".
* Switching the unreplaced `C⁺₂` to `𝔖_ff` requires re-deriving I9 (the residue `16πG/(ν²+9/4)`) and the light-side constants for the new kernel. This was not done.

**m1.** The forward bound is `≤ 6.5·10⁻⁴ N_s`, not `6·10⁻⁴` (`r3`). This is harmless.

**m2.** Say which normalisation `Ω^{SC}` uses. `N₄` is quoted relative to `D₁²⊗D₁²` (min `0.023`); `s1` prints the `(1,3)`-only ratio (min `0.027`).

**m3.** The Jacobian `j(ν)` may now depend on `k`. It is untracked, as stated.
