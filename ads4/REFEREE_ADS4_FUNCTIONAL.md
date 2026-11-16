# Referee report: `ads4/ADS4_FUNCTIONAL.md` and `ads4/code/adsfunc/`

Referee: independent, hostile; 16 Nov 2026. Nothing under review was edited. Scratch code: `ads4/code/referee_adsfunc/`
(`ref_check.py`: mpmath, exact hypergeometric `P_{−1/2+iν}`, no import from `adsfunc/`; `ref_island_phys.py`, `ref_kappa.py`: use the
validated `adsfunc/common.py` on finer grids). Logs in `referee_adsfunc/cache/` (git-ignored). Sources read: CHMRS arXiv:2106.10274 and
arXiv:2201.06602 v1 LaTeX (fetched).

## Reproduction
* Re-ran `t1, t4, t5, t6` and `t2` for `Δ_g = 50, 100`: identical to the stored logs.
* Independent (mpmath): exact `2∫ψ₂ = 2591/210`, `∫p²ψ₂ + ∫ψ₃ = 41/18`; `𝔉₊₋(4,u_max) = 2.31009 (Δ_g=100), 0.31486 (1000)` ✓.
  `++` island, `Δ_g = 100`, `m = 10^{0.8}Λ`, `Λβ = 8.70`: `𝔉_R = −4.81·10⁻³`, relative `−1.009·10⁻²` (doc `−4.90·10⁻³`, `−1.015·10⁻²` at `Λβ = 8.67`) ✓.
  Near-zero `+−`, `m_fl = 21.7`, `J = 244`: flat `2.67519·10⁻⁵`, hybrid shift `5.2738/Λ²` (`Δ_g = 1000`), `5.4314·10⁻⁴` (`Δ_g = 100`) ✓.
  `W = z − 5z²/8 + O(z³)` ✓; `m²sinh²(β/2) = (J+½)²` ✓ (exact; = CHMRS `cosh β = 1 + (2𝒥²+(d−2)²/2)/m²` at `d = 3`).
  `L_AdS(100) = 196.0139` ✓; graviton mass `9/4 = d²/4` agrees with CHMRS `eq:BknuGravity` (`8πG/(ν²+(d/2)²)`).
* Island on physical states (even integer `J`): min relative `𝔉_R = −0.72 %, −0.99 %, −1.08 %` at `J = 26`, `m ≈ 6.13Λ` (`Δ_g = 50, 100, 1000`) — the
  negativity is real, not a continuous-`J` artefact.

## Verdicts
* **Correctness:** dictionary, `ε_k`, signs, forward-term translation and `𝔉_R` are right (checked term by term: `½ψ₃Λ⁶·2σP/m⁶ → σu∫ψ₃P`;
  `−5Λ¹⁰·2(−W)/m⁸ → 10WΛ²u²`; `(5/2)Λ¹⁰·2σP₀/m¹⁰ → 5σu³P₀`). One formula error in the hybrid model (Issue 2) and a self-inconsistent line
  in §1.2 (Issue 5). Numbers reproduce.
* **Soundness:** statement R is overstated (Issue 1); the `7.4Λ` threshold is a grid artefact at the 1 % level (Issue 3).
* **Labels:** mostly honest (abstract: "not a theorem", §0 item 7 "nothing new proved", STATUS row 9 "models and scans, nothing
  proved" — consistent). Exception: "[D|H] given the established heavy formula" overstates what is established (Issue 4).
* **Applicability:** everything at `m ≲ 10Λ` rests on the underived spinning bulk-point formula and the [M] hybrid; at the binding
  hybrid states (`m ≈ 1.1–1.24Λ`, `u ≈ 0.7–0.8`) the Regge kernel is used at `ν ≈ m`, outside its validity, and cross terms `O(uβ²)` are
  not suppressed. The hybrid is an unjustified interpolation there (the doc says so: "a guess").
* **Novelty:** low/incremental. The construction is CHMRS §(holographic sum rules) applied verbatim (`u → −ν²`, Regge/bulk-point overlap,
  forward derivatives `∂_{ν²}C_{k,ν}|₀`, positivity of `(−∂_{ν²})ⁿP` = their Fig. `nu_derivatives`); `ε_k` is bookkeeping for the chosen
  `C⁻` normalisation. arXiv:2201.06602 has **no §4 "AdS translation"** (its §4 is "Results"); it has only a remark (§1, §5) that the
  bounds lift to AdS₄ with `log(M R_AdS)`. The doc's `L_AdS = 65 log(2Λ/3) − 1171/15` makes that remark quantitative, conditionally.
* **Uniqueness claim:** true in the weak sense that the flat limit `Σc_kΛ^{2k}C_k → Σc_kε_k⁻¹B_k` is injective (CHMRS's "uniqueness
  properties"); it fixes nothing at `O(1/Λ²)` and nothing about the light side. Fine as labelled [D|H].
* **Flat limit:** the g₅-free completion is reproduced (leading order, term by term). The literal v1 variant is `Ω − (5/2)Λ¹⁰C⁻_{5,0}`; the doc
  never states it (Issue 6). The flat forward constants (`−5`, `−205/−210`) are absent from `𝔉_R` and are as large as the leading `10𝒥` at
  `J = 4` (`+−`: Regge `202.5u³` vs flat `−10u³`), so `𝔉_R` is meaningless at small `J` — acknowledged.

## Issues
1. **Major (soundness/labels).** Statement R / §0 item 6(i) / STATUS row 9: "`Ω_{g3} > 0` for all heavy states with `m ≥ 7.4Λ`, every
   `Δ_g ≥ 50`, every `C₀ ≤ 1000`". H_env bounds only `O(cosh β/m²)` corrections, but CHMRS's Regge formula carries `O(ν²/m²)` with
   `ν ≤ Λ`, i.e. relative `Λ²/m²`, which exceeds `C₀/m²` as soon as `Λ² > C₀` (all `Δ_g ≳ 35`). The doc itself attributes the island to
   these terms. Measured: relative margin of `𝔉_R` (++) at `m = 8Λ` is `1.3–1.5 %` vs `Λ²/m² = 1.6 %`; at `9Λ`, `5 %` vs `1.2 %`. So
   "for any `C₀`" is vacuous and `7.4Λ` is not a threshold for `Ω`. *Fix:* state R for the leading action only, add `O(ν²/m²)` to H_env
   (e.g. `C₁Λ²/m²`) and quote `M₀(C₀, C₁)`; or say "`m ≳ 10Λ` if the `O(ν²/m²)` terms are ≤ the flat ones".
2. **Minor-to-major (formula error).** Hybrid forward term (§2.4; `common.hybrid_shift`): `5(1−σ_h)(1−P₀)` should be `5σ_h(P₀−1)` (Regge
   `5σu³P₀` minus its flat limit `5σu³`). The coded term exceeds the correct one by `+5(1−P₀)u³ > 0` in both channels (optimistic), and
   violates the doc's own property "equals the leading Regge action for `u → 0` at fixed `β`" at the order of the `B₅` term. Impact small:
   `J = 4` shift `−1.640 → −1.658` (`Δ_g = 20`), asymptote `−2.29·10³ → −2.32·10³/Λ²` (×Λ² at `Δ_g = 1000`); near-zero unchanged (`10⁻¹⁰`).
   *Fix:* correct formula, re-run `t3, t4`.
3. **Minor (numerics).** The `7.4Λ` threshold comes from `t2`'s coarse `β` grid. On finer `Λβ` grids (step 0.05, `ref_island_phys.py`; 0.01, `ref_kappa.py`): `𝔉_R < 0` at `m = 7.36Λ` for
   `Δ_g = 100` (`−4·10⁻⁴` rel.) and `1000` (`−1.9·10⁻³`); last negative `m/Λ` ≈ `7.38` (`Δ_g = 100`), `7.46` (`Δ_g = 1000`). "every
   `Δ_g ≥ 50`, `m ≥ 7.4Λ`" is false at `Δ_g = 1000`. *Fix:* `≈ 7.5Λ`, refine grid, mark "grid".
4. **Major (labels/applicability).** The "established" heavy formula (⟨JJJJ⟩, STATUS row 6) was verified at `k = 4, 6` (`ν = 0, 1, 3`). `Ω` uses
   `k = 2, 3, 5`; `C⁻` at **odd** `k` with two spinning pairs is nowhere established (k = 3 only one pair, unrefereed, `SPINNING_LOWK.md`;
   `k = 5` never computed). So only the `B₄` term is "[D|H]" in the stated sense; `ε_odd = ½` rests on an assumed odd-`k` normalisation.
   *Fix:* relabel `[D|H*]` "assuming the k=4,6 formula extends to k=2,3,5 and to odd k"; add `k = 5` and odd-`k` `C⁻` to gap 3.
5. **Minor (typo in key derivation).** §1.2: "`C_k → (2/ε_k)σ^kJ₀/m^{2k}`" and "`2/m^{2k}` vs `(2/ε_k)/m^{2k}`, so `ε_k = 1, ½`" is
   self-inconsistent (`2/ε₃ = 4`, but flat `C₃ → σ/m⁶`). Should read `2ε_k`. Final `ε_k` correct.
6. **Minor.** State the literal-variant uplift and that its light side carries `−5g₅` (cf. flat Theorem 3).
7. **Minor.** §2.3 cites Lemma 9 (`|P − J₀(νβ)|`) for the curvature at fixed `b̂`; the hybrid uses `J₀(b̂p)`, `b̂ = 2Λsinh(β/2) ≠ Λβ`
   (extra `O(b̂β²)`); `t3`'s 4 % `J`-steps and 3 % `m`-steps miss the near-zero (min at `(22.11, 247)`, flat `6.2·10⁻⁵`, not `2.68·10⁻⁵`;
   `t4` refines it); "`0.06 + 258/(Δ_g+7)`" is a linearisation (2.47 vs exact 2.31 at `Δ_g = 100`). `t2` also skips `Λβ < 0.3`.

## What is proof / model / scan
Proofs (imported): UNIFORM_POSITIVITY Thm 1 (B₂ term of `𝔉_R`), flat Thm 1, kinematics §1.4 [D]. Exact: `J = 4` lift formula (given the
[I] flat-mass identification). Model [M]: hybrid, curvature shift signs (`−2.36·10³/Λ²`, `+5.27/Λ²`), all `Δ_g ≳ 10` statements.
Floating-point grid scans [NE]: `7.4Λ`, island `−1 %`, t3 tables. No AdS4 positivity statement is proved; the document says so.
