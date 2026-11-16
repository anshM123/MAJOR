# Referee report on `ads4/ADS4_LIGHT.md` and `ads4/code/adslight/`

Referee: independent, hostile. Date 16 Nov 2026. Scratch code: `ads4/code/referee_adslight/` (`r1_scalar_grav.py`,
`r2_constants.py`; logs in `cache/`, git-ignored). No file under review was modified.

## Verdicts

* **Correctness (arithmetic, numerics):** correct. Everything I recomputed agrees. That covers the scalar `8πG/(ν²+9/4)`, the dictionary, the 9/4
  pole, `k(0)`, `A`, `B_AdS`, `C₁`, `L_AdS(Λ)`, and the `Δ_g = 100` number.
* **Soundness:** the scalar part and the pole/residue analysis are sound. The `⟨TTTT⟩` light side is sound only conditionally, and one
  load-bearing input is mislabelled and under-discussed (M1 below).
* **Labels:** mostly honest. Two exceptions: §0 item 3 and item 5 claim more than [D given (F),(R)]/(P) (M1), and the `C₁` dependence on (R) is omitted (m2).
* **Applicability:** the headline number `√(t₄²+α_T²) ≤ 2.7·10⁻⁵` at `Δ_g = 100` assumes heavy positivity down to `m ≈ Λ ≈ Δ_g`.
  The companion `ADS4_FUNCTIONAL.md` does not support that (M2).
* **Novelty:** the scalar result is CHMRS's own (BknuGravity), here reproduced. The `H_{d−1}` pole `ν²+d²/4` is standard conformal Regge theory
  (Cornalba et al.; KPZ 1705.02934; CHMRS). The `t₄`-independence of the helicity-diagonal residue in `d = 3` is the Hofman–Maldacena/CDP
  statement (cf. Afkhami-Jeddi–Kundu–Tajdini 1805.07393 for AdS4 shocks; MP 1712.04861 for the AdS phase shift). New here: the bookkeeping
  (degree lemma, exact `g₅` cancellation, `L_AdS` closed form, `B_AdS`). Novelty: low/incremental.

## Independent checks

1. **Scalar graviton [reproduced].** `r1_scalar_grav.py` sums the `T` single-trace Mellin poles. The residue *shape* comes from the
   Costa–Gonçalves–Penedones leading-Regge Mack coefficient, `Γ(a+m)²/(m!Γ(5/2+m))`, `a = 3/2−Δφ`. That shape matches CHMRS (MGravityRegge)
   after reflection. The **normalisation is fixed independently of CHMRS's `R_m` constant**, by the flat limit:
   `Σ_m R_m = −8πG c₁` (Penedones, degree-1 amplitude `8πG s²/(−u)`). The analytic Gauss sum shows that CHMRS's `R_m` satisfies this
   exactly. Contacts are normalised by `s²+t²+u² ↔ −C₂ = 2` (Barnes lemma); `c₁/c₂ = 1/(4Δφ−1)`. The method is a trapezoid rule in `Im z` with
   a vectorised pole sum to `m = 2·10⁵`. Result: `−C₂/(8πG) = 1/(ν²+9/4)` to `≤ 2·10⁻¹³` for `ν ≤ 7` and `3·10⁻¹⁰` at `ν = 12`, at `Δφ = 1.2, 1.75, 2.6`.
   Only the kernel `Γ(z)²Γ(c−z)Γ(c̄−z)` is taken from CHMRS. Its `ν`-structure is consistent with (C2C4Light): by Barnes,
   `⟨z⟩ = |c|²/(c+c̄)` gives `ν² − (4Δφ²−¼)`.
2. **Dictionary.** `κ = dΓ(d/2)³/(4π^{d/2}Γ(d+2)) = 1/256` at `d = 3`. Hence `c_T = 32/(πG)`, which matches `c_T = 32L²/(πG)`. Check: the standard
   `C_T = 3L²/(π³G)` (Osborn–Petkou) divided by `C_T^{free scalar} = 3/(32π²)` gives `32L²/(πG)`, with `c_T` normalised so that a free scalar has `c_T = 1`, consistent
   with `DICTIONARY_T4.md` (`c_T = n_B + n_F`). Correct.
3. **9/4 from first principles.** `P_{−½+iν}(cosh β)` on `H₂` has `−∇² = ¼ + ν²`. A field of mass `m² = Δ₀(Δ₀−1)` (boundary dimension 1)
   has propagator `1/(ν² + (Δ₀−½)²)`. Spin-`J` exchange of dimension `Δ` gives `Δ₀ = Δ−1`, and for `T`, `Δ₀ = 2` gives `ν² + 9/4`. In general `d`:
   `(d−1 − (d−2)/2)² = d²/4`, consistent with CHMRS's `8πG/(ν²+d²/4)`. Correct.
4. **Residue `t₄`/`α_T`-free.** For `ε = e₊e₊`, `|n·ε·n|² = ¼` for every `n`. Rotating the CDP matrix to the helicity basis gives diagonal `1` and off-diagonal
   `(−t₄+iα_T)/4`. Correct. This is consistent with KPZ/MP: their polarisation-dependent AdS phase-shift terms are (polynomial in `ν`)×pole.
   On the `d = 3` helicity diagonal, the `t₄` part vanishes at the pole, so what remains is a contact in `H₂` (the analogue of `∂_z²∂_z̄² log b ∝ δ²(b)` in flat `R²`).
   The flat-limit counting then permits only the constant `c₁`, because flat MHV gets no linear `R³` correction.
5. **Item-5 constants (`r2_constants.py`, plain quadrature + fit, no sympy closed forms).** The fit gives `k(0) = 27.4850982…`, `A = 37.838686`,
   `B_AdS = −60.787499`, `C₁ = 54.8648`. At `Δ_g = 100` (`Λ = 98.499`), the asymptotic form gives `113.451` and the exact `16L_AdS/k(0)` gives `113.545`. The bound is
   `√(t₄²+α_T²) ≤ 2.716·10⁻⁵`. The conversion `240√(·)/Λ⁴` is correct (`57600 = 240²`, AdS units `L = 1`, `x² = |ĝ₃|²Λ⁸`).

## Issues

**Fatal:** none found at the level of what is actually claimed (everything is conditional).

**M1 (major: label and missing sensitivity).** The statement "no Einstein term in `C⁻₃`, `C⁺₄`, `C⁻₅`" is presented as [D given (F),(R)]. It is
not a consequence of (F). (F) controls only the fixed-`p` scaling. It allows any `O(G)` non-polynomial residue `G h(ν)` that decays fast enough:
`h = O(ν⁻⁶)` in `C⁻₃` and `O(ν⁻⁸)` in `C⁺₄`, for example `G/(ν²+9/4)³`. Such residues are not small in `Ω`:
* `ψ₃(0) = 162 ≠ 0`, so a `C⁻₃` residue contributes `≈ 81Λ⁵G∫h`, which is `~Λ³` times the leading `16πGΛ²L_AdS`. For `h = (ν²+9/4)⁻³` (`r2`), the ratio is 1.6·10³ at `Δ_g = 100`.
* The forward terms sit at `ν = 0`, where (F) gives no control at all. For `h = (ν²+9/4)⁻⁴` in `C⁺₄`, the ratio is `~10¹²`.

  What the claim actually needs is the spinning carry-over of CHMRS (eftformula): the light side is a Regge Laurent coefficient, and Regge spin 2 < k gives
  exactly zero. That carry-over is [I] (§2.1 says so). In flat space the per-channel `B₃` pieces `±8πG/p⁴` are non-zero and cancel only via the large-`s`
  decay, so the AdS analogue is not automatic. *Fix:* relabel §0 items 3 and 5 and the §4 table row as "[I]: spinning (eftformula) carry-over". Add
  a separate input row, and state the `Λ³`/`Λ⁸` sensitivity explicitly. The same applies to (P), the exactness of the `|ĝ₃|²` forward terms: any
  non-polynomial piece there has the same leverage.

**M2 (major: applicability).** The table and the `Δ_g = 100` headline use `Λ = √((Δ_g−2)(Δ_g−1))`. That presumes `Ω|_heavy ≥ 0` for heavy states down to `m ≈ Λ`.
`ADS4_FUNCTIONAL.md` §0 item 6 reports something different. The leading Regge action is negative near `m ≈ 6.3Λ` (`Δ_g ≥ 50`). With an
envelope, positivity starts at ≈8.1Λ, and at ≈22Λ for `Δ_g ≳ 1000`. Positivity elsewhere rests on an [M] model, and the `+−` near-zero needs
`Λ ≳ 444` if the `β²` sign is not trusted. *Fix:* attach these caveats to §0 item 7 and §4. Also state the scaling: if one must take
`Λ_eff = Δ_g/n`, the bound on `√(t₄²+α_T²)` grows by about `n⁴`. For example, at `Δ_g = 100` and `n = 8`, the bound is O(0.1), not 2.7·10⁻⁵.

**m1 (minor).** §2.3: the linear term shifts `x` by `431|c₁|/Λ²`, not 860 (860 is the `x²` shift divided by `x₀`). Numerically, at `Δ_g = 100`, the `x²`-shift is
`0.95|c₁|`, which exceeds the retained `C₁/Λ = 0.557`. So quoting `C₁/Λ` as the next term is meaningful only if `|c₁| ≲ 0.5`, or asymptotically. Say so.

**m2 (minor).** If (R) fails, `r(ν) = O(ν⁻⁴)` shifts `C₁` as well, by `∝ ∫ν² r dν` (the `ψ₂ ∝ p²` term). The §4 table names only `B_AdS`.

**m3 (minor).** The "`O(G t₄² ν²)`" in §0 item 3 should read "polynomial of degree ≤ 2 in `ν²`", with the top coefficient fixed or absent by flat matching. These terms feed `b₂, b₁, b₀`.

**m4 (minor, [I] flagged but worth a sentence).** For `⟨TTTT⟩`, the `[TT]` double-twist families have twists `2 + n` for all integers `n ≥ 0` (component twists). A double zero
`sin²(π(τ−τ₀)/2)` covers only one parity class. It must be shown that each helicity structure's kernel has double zeros at *its* families. Otherwise
some MFT families enter at `O(1/c_T)` or `O(1)`. The double-counting worry (item 4 of the brief) is otherwise unfounded. The heavy side is
single-trace only. The heavy-induced `O(1/c_T)` double-trace data are the EFT contacts, and the EFT Laurent coefficient already contains them. `c_T ≫ Δ_g²` (up to the
`10·log` factor) is the right tree-level condition.

**m5 (minor, reproducibility).** `adslight/cache/t2.log` is truncated: section (ii) is missing, although `t2` runs fine when invoked by hand.
`run_all.sh` says t1 takes 10–20 min, while §5 of the document says ≈5 min. `set -e` does not catch failures through `| tee` (add `set -o pipefail`).

**m6 (minor).** §1.2(d) checks the `H₂` propagator only to `10⁻²`–`10⁻³`. That suffices for pole ↔ dimension. The "CHMRS printed `π`" remark should be stated as unresolved, not as a CHMRS error.
