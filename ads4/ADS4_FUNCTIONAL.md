# The AdS4/CFT3 uplift of the flat graviton functional `F_{g3}`: definition, heavy-side action, positivity, remaining gaps

Status: written 16 Nov 2026 by a derivation/numerics agent. **Unrefereed. Not a theorem about AdS4; no AdS4 bound on `ĝ₃` follows
from it.** Code: `ads4/code/adsfunc/` (`bash run_all.sh`; outputs in `adsfunc/cache/`, git-ignored). All numerics are floating
point (no certificates). It imports nothing from, and modifies nothing in, the other code directories.

Marks: **[D]** derived here (or elementary); **[D|H]** derived given the established heavy formula of `SPINNING_4PT.md` (currents) —
for `⟨TTTT⟩` that formula is itself conditional (`SPINNING_TTTT.md`), so every [D|H] statement about stress tensors inherits that condition;
**[DV]** derived and checked numerically; **[NE]** numerical evidence (floating point, read off from numbers); **[M]** statement about the
*model* of §2.4 (a stated interpolation, not derived); **[Q]** quoted (CHMRS = arXiv:2106.10274, LaTeX source read); **[I]** inferred.

Notation (AdS units `R = 1`, `d = 3`): heavy state `(τ, J)`, twist gap `τ ≥ Δ_g` (CHMRS's `Δ_gap` is a *twist* gap), `α = τ − 2`,
`ᾱ = τ + 2J − 1`, `m² = αᾱ`, `cosh β = (α² + ᾱ²)/(2αᾱ)` [Q CHMRS (2.28)–(2.30)]. Flat functional, weights, channels `h ∈ {++, +−}`,
`σ₊₊ = +1`, `σ₊₋ = −1`, `D_h`, `Φ_h` (completed, `g₅`-free: `Φ₊₊ = 10𝒥`, `Φ₊₋ = 10𝒥 − 210`), normal form `𝔉_h(J,u) = m⁴F`:
`FLAT_GRAVITON_BOUND.md` (0.1), §2. Heavy AdS formula: `C^±_{k,ν}[Δ,J] = (2P_{−1/2+iν}(cosh β)/m^{2k})(|c_same|² ± |c_opp|²)(1 + O(1/m²))`,
same ↔ flat `++`, opposite ↔ flat `+−` (`SPINNING_4PT.md` §4.5–4.6).

## 0. Results

1. **Definition [D|H, §1].** With the smearing scale `Λ := √((Δ_g−2)(Δ_g−1)) = Δ_g − 3/2 + O(1/Δ_g)` (the smallest `m` of any heavy
   state, so `ν ≤ Λ ≤ m` for every heavy state, which is the validity condition `ν < m` of the CHMRS bulk-point saddle):

   **`Ω_{g3} := ∫₀¹ dx [ψ₂(x) Λ⁴ C⁺_{2,Λx} + ½ ψ₃(x) Λ⁶ C⁻_{3,Λx}] − 5Λ¹⁰ ∂_{ν²}C⁺_{4,ν}|_{ν=0} + (5/2)Λ¹⁰ C⁻_{5,0}`.**

   Dictionary: `B_k^{(1)}(p²) ↔ ε_k Λ^{2k} C^{(−)^k}_{k, ν=Λp}`, `ε_k = 1` (k even), `½` (k odd); `∂_{p²} ↔ Λ²∂_{ν²}`; `m_flat ↔ m/Λ`;
   `b̂ = (2J+1)/m_flat ↔ 2Λ sinh(β/2)` (exact, since `m² sinh²(β/2) = (J+½)²` [D]). `C⁺` = preserved–preserved tensor Regge moment with
   the even kernel `D₁⊗D₁`, `C⁻` = the same with the odd kernel `D₅⊗D₅`, each normalised so that its leading heavy action is `2P/m^{2k}` per
   unit helicity weight. The forward terms `∂_{p²}B₄^{(1)}|₀`, `B₅^{(1)}|₀` are replaced by the ν-derivative at `ν = 0` of the `k = 4` even moment
   and the `ν = 0` value of the `k = 5` odd moment — not by higher-k moments, and with no smearing (in AdS `ν = 0` is regular; CHMRS use
   `C_{k,0}`, `∂_{ν²}C_{2,ν}|₀` as functionals [Q]). This resolves the [I] "which AdS functional uplifts which flat sum rule" of
   `SPINNING_4PT.md` §4.5 **at the level of leading heavy densities**: the matching of helicity label, power `m^{−2k}`, sign `σ_h` and
   normalisation in the overlap regime `1/Λ ≪ β ≪ 1` fixes `(k, ±, ε_k)` uniquely [D|H]. That the full functional (light side; bulk-point
   regime beyond leading order in `Λ²/m²`) is the uplift is **not** derived [I] (§1.3).
2. **Heavy action, Regge regime [D|H, §2.1].** For `m → ∞` at fixed `β` (and fixed `Λ`):
   `Ω_{g3}[Δ,J;h] = |c^h|² (Λ/m)⁴ 𝔉_R`, with
   **`𝔉_R = 2∫₀¹ψ₂ P_{−1/2+iΛp}(cosh β)dp + σ_h u ∫₀¹ψ₃ P_{−1/2+iΛp}(cosh β)dp + 10Λ²u² W(β) + 5σ_h u³ P_{−1/2}(cosh β)`**,
   `u = Λ²/m²`, `W(β) = −∂_{ν²}P_{−1/2+iν}(cosh β)|₀ = sinh²(β/2) − (5/8)sinh⁴(β/2) + …`. Each term carries its own relative correction
   `O(1/m², ν²/m²)`. For `β → 0` at fixed `b̂` it is the flat `𝔉_h` **to leading order in `u` term by term** (`J₀(b̂p)` ↔ `P_{−1/2+iΛp}`, `𝒥 + ¼` ↔
   `m²W(β)`) [D]; it is *not* the flat `𝔉_h` at `O(u)` inside the `B₂` term (`−up²ψ₂`, helicity terms of `D_h`), which are `O(ν²/m²)`.
3. **Heavy action, bulk-point regime [I, §2.2].** If the spinning analogue of CHMRS's bulk-point formula (§3.4, label `eq:bulkpointbknu`) holds (not derived), then
   `Ω_{g3}[Δ,J;h] = |c^h|² F[m/Λ, J, h] (1 + O(β²) + O(1/m²))` with `F` the flat functional (0.1): positive factor `|c^h|²` (times the positive
   heavy-density normalisation), flat mass `m/Λ = √(αᾱ)/Λ`. Supporting evidence: (i) the scalar theorem [Q]; (ii) overlap consistency with item 2
   [D]; (iii) the measured spin-dependent `1/m²` term `−2l²P′/P` = flat helicity term (`SPINNING_E2E.md` §6.3; `l = 1` confirmed, `l = 2` only as
   `β → 0`) [NE]; (iv) no `O(1/m)` term, so `m_flat = √(αᾱ)/Λ` is fixed to relative `O(1/m²)` [D for currents, conditional for T].
4. **Kinematic lift of the near-marginal `+−, J = 4` point [D given 3].** Heavy states have `m² ≥ (Δ_g−2)(Δ_g+2J−1)`, so the flat region
   actually sampled is `u ≤ u_max(J) = (Δ_g−1)/(Δ_g+2J−1)`. For `+−` (`J ≥ 4`): `u ≤ 1 − 8/(Δ_g+7)`, and the flat margin at threshold is
   `𝔉₊₋(4,u_max) = 2591/210 − 41u/18 − 10u³ ≈ 0.0603 + 32.28·8/(Δ_g+7)`: **2.31 at `Δ_g = 100`, 0.31 at `Δ_g = 1000`** instead of `19/315`
   (`t4`). The model's curvature shift there is `≈ −2.36·10³/Λ²` (dominated by the forward-term curvature `10m²(W − sinh²(β/2))`), so the
   `J = 4` threshold stays positive for `Δ_g ≳ 10` (model) and `Δ_g ≳ 0.25·C₀` with a conservative `1/m²` envelope of size `C₀ cosh β/m²` [M].
   Without the lift (states exactly at `m = Λ`, `J = 4`) one would need `Λ ≳ 200` (curvature) and `Λ ≳ 330 √(C₀/100)` (envelope).
5. **The near-zero `+−`, `m_flat ≈ 21.7`, `J ≈ 244` [M + NE, `t4`].** Flat value `𝔉 = 2.675·10⁻⁵` (`N_cons = 2.45`). The curvature shift is
   **positive**, `+5.27/Λ²` (it is the `65 Q₁`-type AdS term, `UNIFORM_POSITIVITY.md` Lemma 7, `G₁ ≈ 65/(24b²)`), and dominates the envelope
   `C₀ cosh β/m² · N_cons ≈ 5.2·10⁻³ C₀/Λ²` for `C₀ ≲ 10³`; so in the model the near-zero is *lifted* for every `Δ_g ≥ 10`
   (local minimum `3.2·10⁻⁵` at `Δ_g = 1000`, `5.7·10⁻⁴` at `Δ_g = 100`). If the sign of the `O(β²)` term is not trusted, one needs
   `5.27/Λ² < 2.675·10⁻⁵`, i.e. **`Λ ≳ 444`**; with the envelope in addition, `Λ ≳ 444·(1 + 10⁻³C₀)^{1/2}`.
6. **Scans [NE/M, §3].** (i) Leading Regge action `𝔉_R` on physical points (`t2`): positive everywhere for `Δ_g ≤ 20`; for `Δ_g ≥ 50`
   positive for all `m ≥ 7.4Λ` (and any `C₀ ≤ 1000` envelope), negative (relative `−1 %`) only in the `++` island `m ≈ 6.3Λ`, `Λβ ≈ 8.7`, where
   the dropped `O(ν²/m²)` terms are larger than the margin. (ii) Hybrid model on ≈ 10⁴–1.8·10⁴ physical states per `(Δ_g, h)` (`t3`): **positive on
   every scanned state for `Δ_g = 20, 50, 100, 300, 1000`**; with the `1/m²` envelope `C₀ cosh β/m²`: all pass for `C₀ ≤ 100` (all `Δ_g`),
   `C₀ = 300` from `Δ_g = 50`, `C₀ = 1000` from `Δ_g = 100`; failures sit at `m ≈ 1.1–1.24Λ`, `J = 2–5`. Smallest relative margins:
   `++` island `5.3·10⁻⁴`, `+−` near-zero `2.7·10⁻⁵` (`Δ_g = 1000`). Sign-agnostic criterion (B, `β ≤ 0.3`): passes for all `C₀ ≤ 1000` only at
   `Δ_g = 1000` (at `Δ_g = 300` the `+−` near-zero fails by `≤ 2·10⁻⁵` relative, at `Δ_g = 100` by `2·10⁻⁴`).
7. **What is proved, precisely (§3).** Nothing new at the level of a theorem. Rigorous ingredients: `UNIFORM_POSITIVITY.md` Theorem 1 (the
   `B₂` term of `𝔉_R` is `> 0` for all `Λ`, `β`); flat Theorem 1 (`F > 0` for `m ≥ M`). Everything combining them with error terms rests on
   unmeasured or unbounded `1/m²`, `β²` coefficients for `⟨TTTT⟩` and on the underived spinning bulk-point formula.
8. **Light side [I, §4].** If graviton exchange in `⟨TTTT⟩` gives `−C⁺_{2,ν}|_grav = 16πG/(ν² + 9/4)` (the scalar result of CHMRS with the flat
   `16πG/p²` normalisation; not computed), the AdS curvature cuts off the flat logarithm: `L_AdS(Λ) = ∫ψ₂/(p² + 9/(4Λ²)) = 65 log(2Λ/3) −
   1171/15 + O(1/Λ)` [D], and the would-be bound is `|ĝ₃|²Λ⁸ ≲ 37.8387 log Λ − 60.79 + O(1/Λ)` — **conditional on every gap of §4**, not a result.

## 1. The functional

### 1.1 Ingredients [Q + D|H]
* Scalar case [Q CHMRS §2.6, (2.56)–(2.59)]: `C_{k,ν}` = harmonic (`H₂`) transform of the Regge moment `Π_{k,η}` divided by `γ²`, heavy
  action `2P_{−1/2+iν}(η_AdS)/m^{2k}(1 + O(ν²/m², 1/m²))` for `ν ≪ m`; bulk-point formula (§3.4, `eq:bulkpointbknu`) `C_{k,ν} → (2m²−ν²)/(m²−ν²) · P_J(1−2ν²/m²)/(m²(m²−ν²))^{k/2}`,
  valid for all `ν < m` up to `O(J²/m²)`, which is the flat `𝒞_{k,u}` density under `u = −ν²`; light side `= ` flat light side `(1 + O(1/ν²))`.
  CHMRS (§5, item (⋆); graviton light term §4.2 `eq:BknuGravity`): heavy-positive flat sum rules uplift via `𝒞_{k,u} → C_{k,ν}`, `u → −ν²`, with relative corrections in `1/Δ_gap`.
* Spinning case [D|H]: the tensor Regge moments `Π̂^q_{k,ν}` of `SPINNING_E2E.md` §3.1 / `SPINNING_4PT.md` §1; in the component preserved at
  both Regge pairs (harmonic `𝒦_{0,0} = P_{−1/2+iν}`) the covariant even (`D₁⊗D₁`) and odd (`D₅⊗D₅`) kernels give, after normalisation by
  `R_LR_R/4` resp. `g₁²S_LS_R/4`, `C^±_{k,ν} = (2P/m^{2k})(|c_same|² ± |c_opp|²)` (currents: established; stress tensors: conditional).
  `k` even for `C⁺`, odd `≥ 3` for `C⁻` (subtraction factors `f^±_{k;J_e}` [Q CLS-D §3.5.2, as cited in `REGGE_MOMENT_TTTT.md` §5]).

### 1.2 The dictionary flat sum rule ↔ AdS functional [D|H]
Flat heavy densities per state (`FLAT_GRAVITON_BOUND.md` §1.2, units `M = 1`): `C₂ = (2m²−p²)D_h/m⁶`, `C₃ = σ_hD_h/m⁶`,
`C₄ = (2m²−p²)D_h/(m⁸(m²−p²))`, `C₅ = σ_hD_h/(m⁸(m²−p²))`; at fixed `b̂` and `u → 0`: `C_k → (2/ε_k) σ_h^k J₀(b̂p)/m^{2k}`.
AdS: `C^±_{k,ν} → 2(±1 for opp)P/m^{2k}`, and in the overlap `1/Λ ≪ β ≪ 1`, `P_{−1/2+iν}(cosh β) → J₀(νβ)`. With `p = ν/Λ`, `m_flat = m/Λ`,
`b̂ = 2Λ sinh(β/2)`, the densities coincide iff the flat `B_k^{(1)}` is represented by `ε_kΛ^{2k}C^{(−1)^k}_{k,Λp}`:
* label: AdS "same" = flat `++` (`D₊₊ = P_J`), "opposite" = `+−` (`D₊₋ = d̃^J_{44}`; `J ≥ 4` both sides) [D + DVE, `SPINNING_4PT.md` §4.5];
* sign: the odd kernel weighs `+1` (same), `−1` (opposite) = `σ_h` of the odd-k `B_k^{(1)}` (the `−f(s,…)` of (2.26b)) [D + DVE for currents];
* power and normalisation: `2/m^{2k}` vs `(2/ε_k)/m^{2k}`, so `ε_k = 1` (even), `½` (odd) [D];
* `k ↔ k`: forced by the power `m^{−2k}`; parity of `k` ↔ parity of the kernel: forced by crossing (even-k `B^{(1)}` are `s ↔ t` symmetric, odd-k
  antisymmetric; even/odd `C^±`) [D].
No other combination of the established `C^±` has the flat `B_k^{(1)}` density as its leading term, so at this level the uplift is unique.
**Not derived [I]:** that the uplifted functional reproduces the flat density beyond leading order in `Λ²/m²` (i.e. the spinning bulk-point
formula, which must produce `(2m²−ν²)`, `(m²−ν²)^{−(k−2)/2}` and `D_h(1−2ν²/m²)` including the helicity dependence of `d̃^J_{44}`), and that its
light side is the flat light side (no spinning Witten-diagram computation).

### 1.3 The forward terms
Flat `−5∂_{p²}B₄^{(1)}|₀ + 5B₅^{(1)}|₀` ↦ `−5Λ¹⁰∂_{ν²}C⁺_{4,ν}|₀ + (5/2)Λ¹⁰C⁻_{5,0}` [D|H at the level of 1.2]. Heavy action at leading Regge order
(per `|c^h|²`, normal form): `10Λ²u²W(β) + 5σ_hu³P_{−1/2}(cosh β)` [D|H]; `W(β) > 0` and `P_{−1/2}(cosh β) > 0` for all `β` (CHMRS figure `nu_derivatives`, §5.2 [Q]; `t1`
numerics). Exact identity `m²sinh²(β/2) = (J+½)²` [D] and `W = z − (5/8)z² + O(z³)`, `z = sinh²(β/2)` [D, `t1`] give
`10m²W(β) = 10𝒥 + 5/2 − (25/4)(J+½)⁴/m² + …`: the flat `10𝒥` is reproduced, with a curvature correction `O(J⁴/m²)`. The flat constants
(`1` and `41` in `∂_{p²}C₄|₀`, `σ_h` in `C₅|₀`, giving `Φ₊₋ = 10𝒥 − 210`) sit at relative `O(1/m²)` of the AdS forward term and are reproduced only
through the (underived) `ν²/m²` structure of `C⁺_{4,ν}` — the forward term is where the `1/m²` coefficients matter most at small `J`. In
AdS the forward limits need no smearing (the `H₂` spectrum is gapped, `ν ∈ ℝ`, and the graviton pole becomes `1/(ν² + 9/4)`), which removes the
flat hypothesis H3's "unsmeared forward limit" [Q for scalars, I for `⟨TTTT⟩`]. Higher-k moments at `ν → 0` are not needed.

### 1.4 Choice of `Λ` [D]
`m² = (τ−2)(τ+2J−1)` is increasing in `τ` and `J`, so `min m = √((Δ_g−2)(Δ_g−1)) =: Λ` (at `J = 0`). With this `Λ`, every heavy state has
`m ≥ Λ ≥ ν`. (CHMRS's `M = Δ_gap` would leave `J = 0` states with `Λ/m` up to `1 + 3/(2Δ_g)`, where part of the `ν`-range lies beyond the saddle
collision `ν = m` and the bulk-point formula changes; harmless for `++, J = 0` (flat margin ≈ 12.4) but avoidable.) On heavy states
`cosh β/m² = (1/α² + 1/ᾱ²)/2 ≤ 1/(Δ_g − 2)²` [D], so a `1/m²` envelope `C₀ cosh β/m²` is uniformly `≤ C₀/(Δ_g−2)²`; and
`β ≤ 2 log(m/(Δ_g − 2))` (from `α = m e^{−β/2}`).

## 2. Heavy-side action

### 2.1 Regge regime [D|H; `t2`]
Insert `C^±_{k,ν}` into `Ω_{g3}`; with `ν = Λp`: item 2 of §0. `∫ψ₂P_{−1/2+iΛp}(cosh β)dp = F_Λ(β)` of `UNIFORM_POSITIVITY.md`, so the `B₂` term is
`> 0` for all `Λ, β` (Theorem 1 there). Relative corrections: `O(cosh β/m²)` from the heavy formula (coefficients measured, §2.3), and
`O(ν²/m²)` — which, for `β ≫ 1/Λ`, act at the endpoint `ν = O(max(1, 1/β))` because the weights are polynomial at `p = 0` (CHMRS §5 mechanism:
`Σ a_n Q_n/Λ^{n+1}`), so they behave as `O(1/m²)` there [I]; the oscillatory endpoint `ν = Λ` is suppressed by `(1−p)³`.
**Numerics (`t2`, physical points `τ ≥ Δ_g`, `β ≤ 8`, `1 ≤ m/Λ ≤ 10⁴`):** `𝔉_R > 0` everywhere for `Δ_g ≤ 20`; for `Δ_g ≥ 50` `𝔉_R` is negative in a
small `++` region around `m/Λ ≈ 6.3`, `Λβ ≈ 8.7` (min `𝔉_R/Σ|parts| = −0.007 … −0.011`) and positive for all `m ≥ 7.4Λ`. That region is the flat
"`++` island" (`m⁸F ≈ 2.7` at `m ≈ 6`, `b ≈ 8–9`, `FLAT_GRAVITON_BOUND.md` Finding 2), where the dropped `O(u)` terms of the `B₂` density (relative
size `1/m_flat² ≈ 2.5 %`) exceed the margin (≈ 1 %): the leading Regge action alone is not sufficient there and the bulk-point/flat completion
is needed (the hybrid model of §2.4 is positive there). `+−`: `𝔉_R > 0` on all physical points, smallest relative margin `1.4·10⁻²` (`Δ_g = 1000`,
near `m/Λ ≈ 21.5`, i.e. the flat near-zero region).

### 2.2 Bulk-point regime [I]
The spinning analogue of CHMRS `eq:bulkpointbknu` — saddle point of the doubly-spinning shadow integral at `ν ~ m` — is not derived. Conditional on it,
`Ω_{g3}[Δ,J;h] = |c^h|²·F[m/Λ, J, h]·(1 + O(β²) + O(1/m²))`. The flat functional's positivity (Theorem 1 of `FLAT_GRAVITON_BOUND.md`) then gives
positivity at small `β` for `Λ` large enough, with the margins of §0 items 4–5.

### 2.3 Corrections: what is known
| correction | where | status | size |
|---|---|---|---|
| curvature `O(β²)` at fixed `b̂` (`P` vs `J₀`) | all components | exact at leading Regge order (`|P − J₀(νβ)| ≤ sinh²(β/2)/4`, `UNIFORM_POSITIVITY.md` Lemma 9) [D]; in the bulk-point regime [I] | `B₂`: `β²G₁(b)`, `G₁ ≈ 65/(24b²)` at large `b` (positive) |
| forward-term curvature `10m²(W − z)` | `B₄` | [D] at leading Regge order | `−(25/4)(J+½)⁴/m²` |
| `O(1/m²)` of `C^±` (no `ν`) | all | measured, not derived: `l = 1` one pair `k = 2`: `c₂ ≈ −17`, `k = 3`: `−28` (`SPINNING_LOWK.md`, in progress), `k = 4, 6`: `−40 … −300`, growing with `β` (`SPINNING_E2E.md` §6.2); `l = 2` one pair: extra `−17.3 cosh β − 2 sinh²β P′/P` (fitted form); two pairs (`⟨JJJJ⟩`): `−20 … −420` (`SPINNING_4PT.md` §5); `⟨TTTT⟩`: **not measured** | modelled by `C₀ cosh β/m²`, `C₀ ∈ {100, 300, 1000}` |
| `O(ν²/m²)` | `B₂…B₅` | flat part = the flat `O(u)` terms [I, bulk point]; AdS part not derived | — |
| helicity `1/m²` term `−2λ²P′/P` | `+−` | `l = 1` [NE], `l = 2` only `β → 0` [NE]; = flat `d̃_{λλ} − d̃₀₀` [D in flat space] | `λ = 4` for `+−` |
| `O(1/m)` | — | absent (currents [D + NE]; `T` conditional) | 0 |
| flip / helicity-changing components (`Φ_{4,±4}`, `P^{−4}`) | not in `Ω_{g3}` (pres–pres projection) | the functional is built only from pres–pres moments [D]; leakage at `O(1/m²)` [NE for currents] | — |

### 2.4 The hybrid model [M]
To interpolate between §2.1 and §2.2 the scans use `𝔉_hyb = 𝔉_flat(J,u;h) + [𝔉_R − (its β → 0 limit at fixed b̂)]` term by term, i.e. kernel
`D_h(1−2up²) + P_{−1/2+iΛp}(cosh β) − J₀(b̂p)` and forward term `Φ_h − 10𝒥 − 5/2 + 10m²W(β) + 5(1−σ_h)(1 − P_{−1/2}(cosh β))`. It equals the flat action
for `β → 0` and the leading Regge action for `u → 0` at fixed `β`; at intermediate points it is a guess whose error is of the order of the
curvature shift itself. Two criteria are evaluated: (A) `𝔉_hyb − C₀(cosh β/m²)N_comp > 0` (model at face value, `N_comp = Σ|parts|`);
(B) `𝔉_flat − |shift| − C₀(cosh β/m²)N_cons > 0` on `β ≤ 0.3` (sign of the curvature shift not trusted; `N_cons = ∫|Ψ_h||D_h| + |Φ|u³`).

## 3. Positivity

### 3(a) Heavy states with `m ≥ M₀`
**What is proved.** Only the ingredients: (P1) `UNIFORM_POSITIVITY.md` Theorem 1 — the `B₂` term of `𝔉_R` is positive for every `Λ > 0`, `β > 0`;
(P2) flat Theorem 1 — `F[m,J,h] > 0` for `m ≥ M`; (P3) [D] the kinematic facts of §1.4 and §0 item 4. No statement of the form "`Ω_{g3} > 0`
for `m ≥ M₀(Δ_g)`" is proved, because (a) the `1/m²` coefficients of `C^±_{k,ν}` for `⟨TTTT⟩` at `k = 2…5` are neither derived nor measured, and
(b) the `B₃` and forward terms of `𝔉_R` have no positivity theorem.

**Conditional statement R (Regge side) [NE, given H_env].** Hypothesis H_env(C₀): each term of `Ω_{g3}` equals its leading Regge expression
times `1 + r`, `|r| ≤ C₀ cosh β/m²` (on heavy states `≤ C₀/(Δ_g−2)²`, §1.4). Measured magnitudes (`|c₂| ≈ 17–420` for currents, larger with
`β`; §2.3) suggest `C₀ ~ 10²–10³`; for `⟨TTTT⟩` nothing is measured. Then on the `t2` grid (`β ≤ 8`, `1 ≤ m/Λ ≤ 10⁴`, physical points only):
`Ω_{g3} > 0` for all heavy states with **`m ≥ M₀ = 7.4Λ ≈ 7.4(Δ_g − 3/2)`, for every `Δ_g ≥ 50` and every `C₀ ≤ 1000`**; for `Δ_g = 20`:
`m ≥ 6.8Λ` (`C₀ = 300`), `74Λ` (`C₀ = 1000`); `Δ_g = 10`: `5.8Λ` (`C₀ = 100`), `86Λ` (`C₀ = 300`); `Δ_g = 6`: `93Λ` (`C₀ = 100`) — for small
`Δ_g` the envelope, not the functional, sets `M₀`. For `Δ_g ≥ 50` this `M₀` is set by the `++`
island, not by the envelope: below `7.4Λ` the leading Regge action itself is not positive for `Δ_g ≥ 50`, i.e. positivity there is decided by
the `O(ν²/m²)` terms, which only the bulk-point formula supplies. `β > 8` (requires `m > 55Λ`) is not scanned; all terms decay like `e^{−β/2}`
and `Λ_*(β) → 0` (`t6`), so nothing new is expected [I].

**Conditional statement B (bulk-point side) [M + NE].** If the spinning bulk-point formula holds with corrections `O(β²) + O(1/m²)` and the hybrid
model of §2.4 represents the `O(β²)` term, then (`t3`): positive on every scanned physical state (`m ≤ 60Λ`, `J ≤ min(3000, 200m/Λ)`) for
`Δ_g ∈ {20, 50, 100, 300, 1000}`; with H_env(C₀) in addition: `Δ_g ≥ 20` (`C₀ ≤ 100`), `Δ_g ≥ 50` (`C₀ = 300`), `Δ_g ≥ 100` (`C₀ = 1000`). The binding
states are near threshold, `m ≈ 1.1–1.24Λ`, `J = 2–5`, where the forward-term curvature `10m²(W − sinh²(β/2))u³` is large and negative (e.g.
`Δ_g = 20`, `+−`, `J = 4`, `m = 1.2Λ`: flat `7.41`, shift `−1.56`). If the sign of the `O(β²)` term is not trusted (criterion B, `β ≤ 0.3`),
the `++` island needs `Δ_g ≳ 100–300` and the `+−` near-zero `Λ ≳ 444` (`t4`; `Δ_g = 300` fails by `2·10⁻⁵` relative, `Δ_g = 1000` passes).
States with `J > 200 m/Λ` (`b̂ > 400`) are not scanned; there the forward term `≈ 2.5b̂²u²` dominates as in flat region (C) [I].

**The two near-marginal flat points**, quantified (§0 items 4–5, `t4`):
| point | flat margin | AdS effect | needed for AdS corrections not to destroy positivity |
|---|---|---|---|
| `+−, J = 4, m = M` | `19/315 = 0.060` (`N_cons = 65`) | kinematic lift `u ≤ 1 − 8/(Δ_g+7)`: margin `≈ 0.06 + 258/(Δ_g+7)` [D given §2.2]; curvature `≈ −2.36·10³/Λ²` [M] | `Δ_g ≳ 10` (curvature), `Δ_g ≳ 0.25 C₀` (conservative envelope); without the lift `Λ ≳ 200` and `Λ ≳ 330(C₀/100)^{1/2}` |
| `+−, m ≈ 21.7M, J ≈ 244` | `2.675·10⁻⁵` (`N_cons = 2.45`) | curvature shift `+5.27/Λ²` [M], i.e. favourable | any `Δ_g ≥ 10` if the sign is right (`C₀ ≲ 10³`); `Λ ≳ 444 (1 + 10⁻³C₀)^{1/2}` sign-agnostic |

### 3(b) The region `m ~ Λ` (and `ν/m` not small)
* `UNIFORM_POSITIVITY.md` Theorem 1 is a statement about the **Regge** kernel `P_{−1/2+iν}(cosh β)` (valid for `ν ≪ m`) and covers only the `B₂`
  term. At `m ~ Λ` the relevant kernel is the bulk-point one (`P_J(1−2ν²/m²)` with `O(J²/m²)` corrections [Q scalars, I spinning]), so the
  theorem is not the relevant statement there.
* The `B₃` term: no AdS statement. In the Regge regime it is the term that limits positivity (relative size `Λ³/m² × 81Q₀(β)/(65Q₁(β))` against
  the curvature-dominated `B₂` term at `β ~ 1`), and the forward term rescues it, exactly as in flat Lemma M: minimising the `+−` sum
  `−162Q₀u/Λ + 10WΛ²u²` over `u` gives `−(162Q₀)²/(40WΛ⁴)` against `130Q₁/Λ²`, i.e. positivity of the leading terms for `Λ > Λ_*(β) =
  162Q₀/√(5200WQ₁)` [D|H, asymptotic for `Λβ ≫ 1`; `t6`]: `Λ_* = 246, 89, 23, 10.8, 6.5, 3.4` at `β = 0.3, 0.5, 1, 1.5, 2, 3` (at small `β` the flat
  `80/(Λβ)³` term, flat Lemma M, takes over). `t2` confirms `𝔉_R > 0` at all physical points with `m ≥ 7.4Λ` for every `Δ_g` tried.
* The forward terms: `W(β) > 0`, `P_{−1/2} > 0` at leading Regge order; at `m ~ Λ` their AdS form is unknown beyond leading order in `J²/m²`.
* `m ~ Λ` with `β ~ 1` (`J ~ Λ`): in flat variables `b̂ ~ Λ`, where the forward term `≈ 2.5b̂²u² ~ Λ²` dominates every other term by a factor `~ Λ²/65`
  (flat region (C)); positivity there needs only that the AdS kernels stay `O(1)` [I]. The delicate states are `m ≈ Λ` with small `J` (§0 item 4) and
  `m ≈ 5–25 Λ`, `b̂ ≈ 8–25` (`++` island, `+−` near-zero) — both in the overlap/bulk-point region, where the margins are 1 % and `10⁻⁵`
  relative.
* No uniform remainder estimate exists (STATUS row 9).

## 4. Remaining gaps to an AdS4 bound on `ĝ₃` (exact list)
1. **Heavy formula for `⟨TTTT⟩`** (end to end; `l = 2` absolute constant and relative `C⁺`/`C⁻` kernel normalisation; parity-odd heavy exchanges,
   where `λ^{±±} ∝ λ^E ± λ^O` must again give `|c^h|²`), `SPINNING_TTTT.md`.
2. **Spinning bulk-point formula** (§2.2): the doubly-spinning saddle at `ν ~ m`, producing the flat graviton densities with `(2m²−ν²)`,
   `(m²−ν²)^{−(k−2)/2}`, `d̃^J_{44}(1−2ν²/m²)`, and its `O(J²/m²)`, `O(1/m²)` corrections. Without it nothing controls `m ≲ 10Λ`.
3. **`k = 2, 3` for stress tensors**: convergence/swappability of `C⁺_{2,ν}`, `C⁻_{3,ν}` for `⟨TTTT⟩` (Regge intercept; the graviton makes `k = 2`
   borderline exactly as for scalar `C_{2,ν}`), and the heavy formula at `k = 2, 3` (`SPINNING_LOWK.md` does currents; `l = 2` in progress).
4. **Uniform remainder** (STATUS row 9): rigorous bounds on the `1/m²` and `β²` terms over all heavy states, at least in the delicate regions above;
   the margins require `|c₂| ≲ 10²–10³` (Regge) and the sign or a bound of the `O(β²)` term at the near-zero (`Λ ≳ 444` sign-agnostic).
5. **Light side, graviton**: the `⟨TTTT⟩` graviton-exchange Witten diagram in the `C⁺_{2,ν}` moment (expected `16πG/(ν²+9/4)`: universal transverse
   propagator on `H₂` [I]) and its `ĝ₃`-dependent parts in `C⁺_2`, `C⁻_3`, `∂_{ν²}C⁺_4` (expected `πG k(0)|ĝ₃|²(1 + O(1/Λ²))` [I]); matching to
   `16πG/p²` and to the flat `−B|_low` of `FLAT_GRAVITON_BOUND.md` §1.3.
6. **Light side, double traces `[TT]_{n,ℓ}` (τ < Δ_g)** and contact terms: the spinning analogue of CHMRS §4.2 `eq:BknuContactsLeading`; the uplift of the
   `g₅`-free combination measures infinitely many contact couplings with `1/Δ_g²`-suppressed coefficients (CHMRS §5, label `technical assumption`);
   light matter (`𝓜 ≥ 0` in flat space) must be re-signed in AdS.
7. **Parity-odd**: heavy parity-odd exchanges (item 1) and the parity-odd cubic coupling on the light side (`|ĝ₃|²` contains `α_T`, `DICTIONARY_T4.md`).
8. **Large impact parameter**: in AdS the IR is regulated (no eikonal hypothesis E needed) — but this rests on positivity of the Regge-regime action
   for `β ≳ 1` including corrections, which is [NE] only (`t2`), and on the graviton light side of item 5.

## 5. Code (`ads4/code/adsfunc/`, `bash run_all.sh`)
| file | content | output |
|---|---|---|
| `common.py` | weights; closed-form cosine transforms; `∫w P_{−1/2+iΛp}(cosh β)` via Mehler–Dirichlet; `W`, `P_{−1/2}`; flat partial waves; AdS kinematics; `𝔉_R`; hybrid shift | — |
| `t1_validate.py` | kernels vs hypergeometric quadrature (`≤ 2·10⁻¹¹`); `W`, `P₀` vs mpmath; `f̂(9.755) = 0.0231575`; flat `19/315`, `3904/315`, near-zero `2.67·10⁻⁵` | `cache/t1.log` |
| `t2_regge_leading.py` | `𝔉_R` on physical `(m, β)` grids, `Δ_g = 6 … 1000`, envelopes | `cache/t2.log` |
| `t3_hybrid_scan.py` | hybrid model on physical states, criteria (A), (B) | `cache/t3*.log` |
| `t4_marginal.py` | `J = 4` threshold and `J ≈ 244` near-zero vs `Δ_g` | `cache/t4.log` |
| `t5_light.py` | `L_AdS(Λ)` and the would-be constants | `cache/t5.log` |
| `t6_regge_estimate.py` | `Q₀`, `Q₁`, `W` and the Regge-regime threshold `Λ_*(β)` | `cache/t6.log` |
