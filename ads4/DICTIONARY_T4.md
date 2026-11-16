# Dictionary: 4D cubic graviton couplings ↔ CFT3 ⟨TTT⟩ data, the conditional CFT3 form of the flat bound, and the literature

Status: written 16 Nov 2026 by a research agent; **not refereed**. Scripts and recorded output in
`ads4/code/dictionary/` (§10). Primary sources were read from their arXiv LaTeX sources (equation numbers are
arXiv-version numbers; where harvmac/letter numbering could not be fixed reliably the LaTeX label is quoted).

Marks: **[Q]** quoted (reference, equation); **[DV]** derived here and verified by an independent script;
**[D]** derived here, not independently verified; **[I]** inferred (argument given, no computation);
**[H]** heuristic/conditional — not a theorem.

---

## 0. Results in one page

**Bulk action (AdS4, radius `L`)**, normalised as in CLPS (arXiv:2201.06602 (2.9)–(2.11)) [Q]:

`S = (1/16πG) ∫ d⁴x √−g [ R + 6/L² − (α₃/6) W⁽³⁾ − (α̃₃/6) W̃⁽³⁾ ] + …,   ĝ₃ = α₃ + i α̃₃`,

`W⁽³⁾ = W_{μν}{}^{ρσ} W_{ρσ}{}^{αβ} W_{αβ}{}^{μν}`, `W̃⁽³⁾ = W_{μν}{}^{ρσ} W_{ρσ}{}^{αβ} W̃_{αβ}{}^{μν}`, `W̃_{μνρσ} = ½ ε_{μν}{}^{αβ} W_{αβρσ}`.
(CLPS write `R` for `W`; in flat space they coincide on the three-point vertex, and in AdS4 the difference is a
Gauss–Bonnet term, Ricci terms and a cosmological-constant shift, none of which affects `t₄`; §3.4.)

**The dictionary** (leading order in `α₃/L⁴`, tree level = leading order in `1/c_T`):

| quantity | formula | status |
|---|---|---|
| parity-even non-Einstein parameter | **`t₄ = −240 α₃/L⁴ = −240 Re ĝ₃ /L⁴`** | [DV] (two independent literature checks, §3.3) |
| parity-odd parameter (Chowdhury–David–Prakash normalisation) | `α_T = ∓240 α̃₃/L⁴` (sign convention-dependent) | [D]/[I] (§4) |
| convention-free combination | **`t₄² + α_T² = 57600 |ĝ₃|²/L⁸`**, i.e. `|c| = 60 |ĝ₃|/L⁴` with `c = (t₄+iα_T)/4` | [D] |
| general `d` (parity-even) | `t₄ = −2(d−1)d(d+1)(d+2) α₃/L⁴` (`−720 α₃/L⁴` in d = 4) | [DV] at d = 3, 4 |
| two-point normalisation | `C_T = 3L²/(π³G)`; `c_T ≡ C_T/C_T^{free real scalar} = 32 L²/(πG)` | [Q]+[DV] (ABJM check) |
| free-field basis | `c_T = n_B + n_F`, `t₄ = 4(n_B − n_F)/(n_B + n_F)` (`n_F` = Majorana count) ⇒ `n_B − n_F = c_T t₄/4 = −1920 α₃/(π G L²)` | [Q]+[D] |
| cutoff | `M L = Δ_gap` (CHMRS convention; heavy states `m²L² = (Δ−J−2)(Δ+J−1)`) | [Q] convention |
| IR cutoff | `m_IR → 3/(2L)` (AdS4 graviton propagator `1/(ν²+9/4)`, CHMRS §6.3), i.e. `log(M/m_IR) → log(2Δ_gap/3)` | [H] |

Hofman–Maldacena normalisation of `t₄` in `d = 3` (used throughout; §1):
`⟨E(n)⟩ = (E/2π) [ 1 + t₄ ( |ε_ij nⁱnʲ|²/(ε*_ij ε_ij) − 1/4 ) ]`, i.e. for a linearly polarised state
`⟨E(θ)⟩ = (E/2π)[1 + |c| cos(4θ − φ)]`, `|c| = √(t₄² + α_T²)/4`; free scalar `t₄ = +4`, free Majorana `t₄ = −4`;
Einstein gravity `t₄ = α_T = 0`; collider disc `t₄² + α_T² ≤ 16`.

**Conditional CFT3 form of the flat bound** [H] (§6). If Theorem 3(i) of `FLAT_GRAVITON_BOUND.md`
(`|ĝ₃|²M⁸ ≤ 37.84 log(M/m_IR) − 45.21 + Δ_IR`) lifts to AdS4 with controlled corrections, then for holographic CFT3s

`t₄² + α_T² ≤ (57600/Δ_gap⁸) [ 37.84 log Δ_gap + κ ] (1 + o(1))`,  i.e.  `|t₄ + iα_T| ≤ 1476 √(log Δ_gap + κ/37.84) / Δ_gap⁴`,

`57600 × 37.8387 = 2.18 × 10⁶`; `κ = O(1)` is **not** determined (with the CHMRS IR prescription and no AdS
correction, `κ = −45.21 − 37.84 log(3/2) = −60.55`, i.e. constant term `−3.49 × 10⁶`). The leading coefficient is
the robust part; the constant is not.

**Literature verdict** (§7): **no rigorous (or explicit-constant) `Δ_gap`-dependent bound on `t₄` or `α_T` in CFT3
was found.** Closest prior: Afkhami-Jeddi–Kundu–Tajdini arXiv:1805.07393 §7 (d = 3, heuristic, parametric,
`|t₄| ≲ ln Δ_gap/Δ_gap⁴`, no constants); CLPS 2201.06602 (flat D = 4 bound + a one-sentence remark that it
"translates" to holographic CFT3 with `log(M R_AdS)`, no formula, no proof); the explicit-constant analogue in
d = 4 is CLPS-higher-D arXiv:2205.01495, `|a−c|/c ≤ 23/Δ_gap² + O(Δ_gap⁻⁴)` (uplift asserted, not proven for
spinning correlators).

---

## 1. CFT3 conventions

1. **Hofman–Maldacena, general `d`** [Q] (Hofman–Maldacena 0803.1467 §2 for d = 4; general d in Buchel–Escobedo–Myers–
   Paulos–Sinha–Smolkin 0911.4257, and Li–Lü–Ma 2108.05384, equation preceding their (t2-t4-ABC)):
   `⟨E(n)⟩ = (E/Ω_{d−2})[1 + t₂(ε*_ij ε_il nʲnˡ/ε*ε − 1/(d−1)) + t₄(|ε_ij nⁱnʲ|²/ε*ε − 2/(d²−1))]`.
   In `d = 3` the `t₂` bracket vanishes identically (for any traceless symmetric 2×2 `ε`, `ε*_ijε_il nʲnˡ/ε*ε = ½`),
   so only `t₄` survives, with constant `2/(d²−1) = 1/4`. **[DV]** (`t4_dictionary.py`, `check_HM_d3`: angular mean of
   `|ε_ij nⁱnʲ|²/ε*ε` is `1/4`, range `[0, ½]` for linear and constant `¼` for circular polarisation.)
   *Misprint found:* Bueno–Cano–Ruipérez 1802.00018 eq. (ees) prints `−1/8`; with `−1/8` the flux would not integrate
   to `E` and positivity would give `−8/3 ≤ t₄ ≤ 8`, not the `|t₄| ≤ 4` they (correctly) state. Their `t₄` is the
   standard one with `−1/4`.
2. **Free fields** [Q] (Buchel et al. 0911.4257; BCR 1802.00018 footnote in §5; Chowdhury–David–Prakash 1707.03007,
   "`t₄ = −4(n_f−n_s)/(n_f+n_s)`"): real scalar `t₄ = +4`, Majorana `t₄ = −4`; `C_T^{scalar} = C_T^{Majorana} = 3/(32π²)`
   (Osborn–Petkou normalisation `⟨TT⟩ = C_T I/x^{2d}`), so `c_T = n_B + n_F` and `t₄ = 4(n_B − n_F)/(n_B + n_F)`.
3. **Parity-odd** [Q] (CDP 1707.03007: the energy matrix displayed in their introduction and re-derived in their §3): in the basis of the two linear
   polarisations, at a detector in direction `n`, the flux matrix is
   `[[ (E/4π)(1 − t₄/4), α_T E/(16π) ], [ α_T E/(16π), (E/4π)(1 + t₄/4) ]]`, `α_T = 8π⁴ p_T/(3(n_f + n_s))`; positivity
   gives the disc `t₄² + α_T² ≤ 16`. Equivalently, in the helicity basis `h = ±2`, off-diagonal element `c e^{∓4iθ}`
   with `c = (t₄ + iα_T)/4` up to a phase convention (`EINSTEIN_PROGRAM.md` §3, (E5)). Only `|c|` is convention-free.
4. **Other normalisations met in the literature**: Li–Meltzer–Poland 1511.08025 footnote: `λ_TTT = 3C_T(60 − t₄)/(2¹⁰π)`
   [Q, not checked]; AKT 1805.07393 §7 use `ñ_s, ñ_f, ñ_odd` with Einstein ⇔ `ñ_s = 4ñ_f`, `ñ_odd = 0` [Q] (their
   normalisation of `ñ_f` relative to `n_F` was not checked; only the statement "`t₄ ∝ (ñ_s − 4ñ_f)/c_T`" is used).

## 2. Flat side: ĝ₃ ↔ Lagrangian ↔ CEMZ's α₄ (D = 4)

CLPS [Q]: `M(1⁺2⁺3⁻) = √(8πG)[12]⁶/([13]²[23]²)`, `M(1⁺2⁺3⁺) = (ĝ₃/2)√(8πG)([12][13][23])²` (2.6);
`f_low ⊃ 8πG/(stu) + 2πG su|ĝ₃|²/t` (2.7a); action (2.9) with `−(1/3!)(α₃R⁽³⁾ + α̃₃R̃⁽³⁾)`, `g = η + √(32πG) h`,
and `ĝ₃ = α₃ + iα̃₃` (2.11). CEMZ (1407.5597, §2 eq. \totalr and App. B eq. \threescalarandvector) [Q]: `A_ggg = √(32πG)[A₁² + α₂A₁A₂ + α₄A₂²]`,
`A₁ = ε₁·ε₂ ε₃·p₁ + cyc`, `A₂ = ε₁·p₂ ε₂·p₃ ε₃·p₁`; `α₂`-structure absent in D = 4. CLPS-higher-D (2205.01495, after
(eq:action)) [Q]: `(1/16πG)(α₄/12)C³`, `C³ = 3I₁ − 4I₂′`, "chosen so that the on-shell vertex agrees with CEMZ".

Verified numerically (`three_point_vertices.py`, random complex three-particle kinematics, explicit 4×4 tensors;
output `out_three_point_vertices.txt`) **[DV]**:
* linearised Riemann `R_i = (κ/2)F_i⊗F_i`, `F_i = p_iε_i − ε_ip_i`; on-shell `I₁ = R₁^{ab}{}_{cd}R₂^{cd}{}_{ef}R₃^{ef}{}_{ab} = 8A₂²`,
  `I₂ = R₁{}_a{}^c{}_b{}^d R₂{}_c{}^e{}_d{}^f R₃{}_e{}^a{}_f{}^b = 4A₂²` (the ECG/QTG contraction), `I₂′ = 4A₂²` (the 2205.01495
  contraction), equal for all 3! assignments; hence **`I₂ = I₁/2` on the three-point vertex** (and `C³ = 8A₂²`).
  (`I₂ = I₁/2` also holds as an *algebraic identity for Weyl tensors in dimension 4 and 5*, fails in 6: `t4_dictionary.py`,
  `check_weyl_identity`.)
* `|A₂(+++)|² = |[12][23][31]|²/8`; `|A₁(−−+)| = |⟨12⟩³/(⟨23⟩⟨31⟩)|/√2` (so `√(32πG)A₁²` is CLPS's Einstein vertex);
  `I₁(+,+,−) = 0`; reference-spinor independence.
* duality: `*F = ∓iF` on helicity `±` (here `−i` for `+` with `ε^{0123} = +1`), and `[R R R̃](+++) = −i I₁(+++)`.

Consequences (exact): `L = (1/16πG)(−α₃/6)I₁` ⇒ `M(+++) = −(α₃/6)(1/16πG)·3!·(κ/2)³·8A₂² = −2α₃√(32πG) A₂²`, so
`|M(+++)| = (|α₃|/2)√(8πG)|[12][13][23]|²` — **CLPS's `|ĝ₃| = |α₃|` confirmed** — and

**`α₄^{CEMZ} = −2α₃`** (= the `α₄` of 2205.01495). With the parity-odd term, the all-plus vertex is `∝ α₃ ∓ iα̃₃`
(sign = orientation convention), i.e. `|ĝ₃|² = α₃² + α̃₃²` is the only convention-free combination.

## 3. AdS side: α₄ → t₄

### 3.1 The formula used [Q]
CEMZ §3.3 (AdS discussion), eq. \timedelay: for the linearised equation
`δR_{μν} + α₂ Ř_{(μ}{}^{ραβ}δR_{ν)ραβ} + (α₄/2)[∇_{(μ}∇_{ν)}Ř^{αβρσ}]δR_{αβρσ} = 0` (`α_i` in units of `L`) on the AdS shock
`h = zϖ(ρ)δ(u)/(1−ρ²)` (eq. \exactsolution, `ϖ₄ ∝ −8(1−ρ²+(1+ρ²)log ρ)`),
`t₄(ρ) = −Dα₄ (1−ρ²)² (ϖ′−ρϖ″)/(4ρϖ) · [D(1+ρ²)² + 2(1+ρ⁴)]/(4ρ²)`, with `ρ → 0` = flat space (`ρ = b/2`) and `ρ → 1`
= the Hofman–Maldacena energy-correlator configuration (CEMZ text after \timedelay; App. F).

### 3.2 Limits **[DV]** (`t4_dictionary.py`, `check_varpi`, `check_limits`; sympy, exact)
* `ϖ_D` from the hypergeometric formula equals the printed closed forms for D = 4, 5 (constant ratio).
* `ρ → 1`: **`t₄ → D(D²−1)(D−2) α₄`** for D = 4, 5, 6, 7 (120, 360, 840, 1680); `t₂ → (D−1)(D−2)(α₂ − (4D−2)α₄)`.
* `ρ → 0`: `16ρ⁴ t₄(ρ) → (D−4)(D−2)D(D+2)α₄`, which is CEMZ's flat formula (App. B, \energycorr) with `b = 2ρ` (D = 5, 6, 7);
  in D = 4 the helicity-flip term `t₄(ρ)ϖ(ρ)ρ⁴ → 24α₄` stays finite (the leading term is logarithmic, hence the
  vanishing `(D−4)` of the normalised flat formula). So CEMZ's AdS `α₄` is the amplitude `α₄`.

With `α₄ = −2α₃`, `D = d+1`: **`t₄ = −2(d−1)d(d+1)(d+2) α₃/L⁴`; d = 3: `t₄ = −240 α₃/L⁴`.**

### 3.3 Independent checks against published nonperturbative results **[DV]** (`check_GB`, `check_QTG_ECG`)
* **Gauss–Bonnet, all d** (tests the `ρ → 1` normalisation of the formula): CEMZ action `(1/16πG)(R + λ_C X₄)`, flat
  `t₂ = 4λ_C(D−4)(D−2)/b²` (CEMZ §2, eq. \totaimeGB) vs `t₂ = (D−2)(D−4)α₂/b²` ⇒ `α₂ = 4λ_C`; then `ρ → 1` gives
  `t₂ = 4(D−1)(D−2)λ_C`, which with `λ_C = λL²/((D−3)(D−4))` equals Buchel et al. 0911.4257 `t₂ = 4λ(D−1)(D−2)/((D−3)(D−4))`
  at leading order, for D = 5…8. (CEMZ's text at \sumscalar/below \resultforamplitudesums says `α₂ = ±λ_GB/4`; that
  is inconsistent with their own two formulas and is not used.)
* **Quasi-topological gravity, d = 4**: Myers–Paulos–Sinha 1004.2055 (action) `+(7/8)μL⁴ Z′₅`, `Z′₅ = I₁ + Ricci terms`;
  ⇒ `α₃ = −(21/4)μL⁴`, `t₄ = −720α₃ = 3780μ` — **equals MPS (t2t4) `t₄ = 3780μ f∞²/(1−2λf∞−3μf∞²)` at leading order.**
  (Same with Myers–Robinson's `(7/4)μL⁴Z₅`, `Z₅ = I₂ + …`, via `I₂ = I₁/2`.)
* **Einsteinian cubic gravity, d = 3**: BCR 1802.00018 (ECG) `−(μL⁴/8)P`, `P = 12I₂ + I₁ + Ricci terms ≅ 7I₁`
  ⇒ `α₃ = (21/4)μL⁴`, `t₄ = −240α₃ = −1260μ` — **equals BCR (t4ecg) `t₄ = −1260μf∞²/(1−3μf∞²)` at leading order**,
  including the sign. BCR obtained it by a different route (twist-operator dimensions via Chu–Miao's `h_q`
  formula), so this is an independent check of exactly the d = 3 coefficient `−240`.
* Caveat on `t₂` (irrelevant in d = 3): the QTG `t₂ = −2088μ` of MPS is *not* reproduced by naively combining the
  `α₄` term of CEMZ's `t₂(ρ)` with the Gauss–Bonnet-like pieces that `Z′₅` generates on the AdS background
  (`−2340μ` by my accounting). The two-derivative structure evidently receives an additional AdS-curvature
  contribution from the cubic term that CEMZ's parametrisation absorbs into `α₂`. `t₄` is immune: only the
  four-derivative structure produces the `(ε·n)⁴` angular dependence. [D]

### 3.4 W³ versus R³ on AdS4 [D]
`R_{ab}{}^{cd} = W_{ab}{}^{cd} − (2/L²)𝟙` on bivectors (on-shell), so `R⁽³⁾ = W⁽³⁾ − (6/L²)W⁽²⁾ + const`; `W⁽²⁾` = Gauss–Bonnet
+ Ricci terms (topological + field-redefinable in 4D), the constant shifts `L` (the `f∞` factors). Likewise
`RRR̃ = WWW̃ + (topological WW̃)`. Hence at leading order `t₄`, `α_T` depend only on `α₃`, `α̃₃`; `C_T` shifts at relative
`O(α₃/L⁴)`.

Structure counting: AdS4 has exactly three cubic graviton vertices (Einstein, `W³`, `WWW̃`) and CFT3 exactly three
`⟨TTT⟩` structures (`n_B`, `n_F`, `p_T`; Giombi–Prakash–Yin 1104.4317). The dictionary is therefore exact at tree level
(no `1/Δ_gap` corrections to the *identification*); corrections are `O(α₃²/L⁸)` from nonlinear terms (cf. the `f∞`).

## 4. Parity-odd coupling
In 4D the only parity-odd cubic vertex is `WWW̃` (`R̃⁽³⁾` of CLPS); `RR̃` (Pontryagin) is topological (contributes only a
contact term to `⟨TT⟩`; AKT 1805.07393 §7.1 footnote). On-shell the `RRR̃` vertex vanishes for mixed helicities and equals
`∓i I₁` on all-plus/all-minus [DV], so the all-plus vertex is `∝ α₃ ∓ iα̃₃`. In the energy-flux matrix the
`α₃`-term is the helicity-flip off-diagonal element; the `α̃₃`-term rotates its phase by `π/2` with the same modulus.
Therefore `c ∝ ĝ₃` or `ĝ₃*`, `|c| = 60|ĝ₃|/L⁴`, and `α_T = ∓240α̃₃/L⁴` in CDP's normalisation, the sign depending on the
orientation conventions (`ε^{0123}`, helicity labels, CDP's `p_T` sign) — **not fixed here**. [D] for `|c|`, [I] for the
sign. In CDP's `p_T`: `p_T = 3c_Tα_T/(8π⁴) = ∓2880 α̃₃/(π⁵GL²)` [D, sign unverified].

## 5. G_N, c_T, M, m_IR
* `C_T = (d+1)Γ(d+1)/((d−1)π^{d/2}Γ(d/2)) · L^{d−1}/(8πG)` [Q, standard] `= 3L²/(π³G)` in d = 3; `c_T = 32L²/(πG)`. Check:
  ABJM `F = πL²/(2G) = π√(2k)N^{3/2}/3` gives `c_T = 64√(2k)N^{3/2}/(3π)`, the known value [DV] (`check_CT`).
  `ĝ₃`, `t₄` and the bound are `G`-independent at leading order (`|ĝ₃|²` is normalised to the Einstein vertex).
* `M L = Δ_gap` [Q] (CHMRS 2106.10274 §1: "`Δ_gap = M R_AdS`"); heavy-state mass `m²L² = (Δ−J−d+1)(Δ+J−1)` (CHMRS (2.28)).
  For the lightest higher-spin single trace (`J = 4`, d = 3): `m L = Δ_gap(1 − 3/(2Δ_gap) + …)`; the flat `M` is
  "the mass above which all states start", so the precise identification is `M L = min_{J>2 single-trace} m L` [D]; this
  shifts `κ` only.
* `m_IR`: in AdS4 the graviton `t`-channel term in CHMRS's scalar functionals becomes `8πG/(ν² + (d/2)²)` (CHMRS §6.3,
  eq. \ads4 graviton) — an effective IR mass `3/(2L)`. For the graviton functional `F_{g3}` the IR-divergent piece is
  `∫ψ₂(p)/p²`, `ψ₂ ≈ 65p`; replacing the sharp cutoff `p > ε` by `1/(p²+a²)` changes it by `O(a)` (`check_regulator`:
  `0.296, 0.024, 0.0023` at `a = 10⁻², 10⁻³, 10⁻⁴`). So, **if** the spinning AdS4 functional has the same regulator,
  `log(M/m_IR) → log(2Δ_gap/3) + O(1/Δ_gap)` [H]. CHMRS's own scalar AdS4 bound has `log(0.124 Δ_gap)` (§6.3, eq. \eqref{eq:adsfourbound}): the
  `O(1)` inside the log is functional-dependent, so only "`log Δ_gap + O(1)`" should be claimed.

## 6. What the flat bound would say about CFT3 [H — conditional, not a theorem]

**Statement (conditional).** Let a unitary CFT3 have `c_T → ∞` with planar factorisation, and all single-trace
operators of spin `J > 2` (other than `T`) at `Δ ≥ Δ_gap ≫ 1`. Assume the AdS4 lift (L1)–(L5) below. Then, with
`A = 37.8387` (certified, `FLAT_GRAVITON_BOUND.md` Theorem 3(i)):

`t₄² + α_T² ≤ (240/Δ_gap⁴)² [A log Δ_gap + κ] (1 + o(1)),  κ = O(1) undetermined`
  `⟺ 16|c|² ≤ 2.18×10⁶ (log Δ_gap + κ/A)/Δ_gap⁸ ⟺ |t₄ + iα_T| ≤ 1476 √(log Δ_gap + κ/A)/Δ_gap⁴.`

In bulk language: `|α₃ + iα̃₃| ≤ √(A log(ML) + κ) L⁴/(ML)⁴`. With CLPS's own uncertified optimum (4.4)
(`24.9 log − 27.6`) the prefactor becomes `240√24.9 = 1198` instead of `1476`.

Naive numbers (`translation()`), `κ = −60.55` (no AdS correction, CHMRS IR prescription), Theorem 3(i) needs
`m_IR ≤ M/100`, i.e. `Δ_gap ≥ 150` with `m_IR = 3/(2L)`: `Δ_gap⁴√(t₄²+α_T²) ≤ 2.7×10³ (Δ_gap = 150)`, `3.4×10³ (10³)`,
`4.4×10³ (3×10⁴)`. For comparison, the collider disc is `√(t₄²+α_T²) ≤ 4`; formally the right side drops below 16 already at
`Δ_gap ≈ 5`, but there it is meaningless (it turns negative below `Δ_gap ≈ 4.95`, which shows that the `O(1)` constant is
uncontrolled at small `Δ_gap`).

**Hypotheses for the lift** (cross-referenced to `STATUS.md`):
* (L1) dispersive CFT sum rules for `⟨TTTT⟩` in d = 3 whose flat limit is CLPS's `B_k^{(1)}` (fixed-`u`, MHV) — CLS-D
  2311.04271 give the general-d bulk-point dictionary but leave the Regge regime and the non-shock sum rules open
  (chain link 6: normalisation not established);
* (L2) heavy-state positivity of the AdS4 version of `F_{g3}` for all `(Δ, J)` above the gap, including the curved region
  `β ≳ 1` (links 5, 9);
* (L3) uniform control of the `1/m²` corrections, including the opposite-sign `B₃` term (link 7);
* (L4) light side: graviton exchange gives `16πG L(ε)` with `ε → 3/(2Δ_gap)`; double traces `O(1/c_T)`; light matter
  of spin ≤ 2 (including KK modes) enters with the helpful sign (`𝓜 ≥ 0`) (link 10);
* (L5) the IR remainder `Δ_IR` is absent in AdS (as CHMRS found for scalars: `Q₁(β) > 0`), or bounded.
None of (L1)–(L5) is proved here.

**Shape of the result versus the folklore.** Parametric estimates in the literature are linear in the coupling:
`|t₄| ≲ ln Δ_gap/Δ_gap⁴` (AKT). The dispersive bound is quadratic, `|ĝ₃|² ≲ log/M⁸`, so it gives
`|t₄| ≲ √(ln Δ_gap)/Δ_gap⁴` — a factor `√ln Δ_gap` stronger, with an explicit constant — and treats `t₄` and `α_T`
together, matching the rotation-invariant collider disc.

## 7. Literature on gap-dependent bounds on non-Einstein ⟨TTT⟩ data

Search method: INSPIRE citation lists of 2201.06602 (207 records), 2205.01495, 1407.5597 (since 2023), 1805.07393,
1610.09378, 1709.03597, 2106.10274, 2311.04271, 1707.03007, 1904.05905; INSPIRE/arXiv abstract searches ("higher-spin
gap", "t_4", "stress tensor"+"Regge"+"holographic", AdS4 graviton bounds) to 16 Nov 2026. A negative search result,
not a proof of absence.

| ref. | d | object | rigour | Δ_gap scaling | constants |
|---|---|---|---|---|---|
| Hofman–Maldacena 0803.1467; Buchel et al. 0911.4257 | 4; general d | collider bounds (d = 3: `|t₄| ≤ 4`) | derivation from energy positivity; ANEC later derived (1605.08072, 1610.05308) | none | exact |
| Li–Meltzer–Poland 1511.08025 | 3 | collider bounds from lightcone bootstrap | physics-level derivation | none | exact |
| Chowdhury–David–Prakash 1707.03007 (and 1812.07774) | 3 | `t₄² + α_T² ≤ 16` | derivation (collider) | none | exact |
| CEMZ 1407.5597 | bulk D ≥ 4 (incl. AdS4) | `α₄^{1/4} ≲ 1/Δ_gap` | heuristic ("numerical coefficient that we cannot fix") | `|t₄| ≲ Δ_gap⁻⁴` | none |
| Afkhami-Jeddi–Hartman–Kundu–Tajdini 1610.09378 | **4** | `|t₄/c| ≲ Δ_gap⁻⁴`, `|t₂/c| ≲ Δ_gap⁻²`, `|a−c|/c ≲ Δ_gap⁻²` | infinite gap: derivation under assumptions; finite gap: "an estimate rather than a strict bound" | `Δ⁻⁴` | none |
| AJHKT 1709.03597 | general d | shockwave ↔ OPE dictionary | derivation | no TTT gap bound (no `t₄` in text) | — |
| **Afkhami-Jeddi–Kundu–Tajdini 1805.07393 §7** | **3** | `|ñ_s−4ñ_f|/c_T, |ñ_odd|/c_T ≲ ln Δ_gap/Δ_gap⁴`; bulk `W³`, `W̃WW` couplings `∼ ln Δ_gap/Δ_gap⁴` | heuristic (assumes heavy OPE coefficients cancel "large numerical factors") | `ln Δ/Δ⁴` (log violation noted) | none |
| Kulaxizi–Parnachev–Zhiboedov 1705.02934 | general d > 2 | Einstein `⟨TTT⟩`-type constraints at infinite gap; finite gap `|β₄| ≲ log Δ_gap/Δ_gap²` for `⟨JTJ⟩` | derivation (infinite gap) + parametric | `log Δ/Δ²` (currents) | none |
| Meltzer–Perlmutter 1712.04861 | general d (d = 3 remarks) | `⟨TTO⟩`, `α_n ≲ Δ_gap⁻ⁿ`, `|a−c|/c ≲ Δ_gap⁻²` | heuristic/parametric (infinite gap: derivation) | `Δ⁻ⁿ` | none |
| Belin–Hofman–Mathys 1904.05892 | 4 | `a = c` at infinite gap; `|α₁| ≲ Δ_gap⁻²` | derivation (infinite gap) + estimate | `Δ⁻²` | none |
| Kologlu–Kravchuk–Simmons-Duffin–Zhiboedov 1904.05905 | general d; explicit d = 4 | superconvergence sum rules for `(t₄+2t₂)²`; planar `t₂ = t₄ = 0` without single traces | derivation with stated loopholes | none explicit (`a_{J=3} ≲ 1/(Δ_gap⁴ c_T)` estimate) | none |
| Belin–Hofman–Mathys–Walters 2011.13862 | 4 | light-ray algebra | derivation | none | — |
| Caron-Huot–Mazáč–Rastelli–Simmons-Duffin 2106.10274 | AdS_{d+1} incl. AdS4 | **scalars only**: e.g. `g₂ ≥ −4516 c_T^{free}/(Δ_gap² c_T) log(0.124Δ_gap)` | physics-level rigorous, numerically checked positivity | `log Δ/Δ²` | explicit |
| **Caron-Huot–Li–Parra-Martinez–Simmons-Duffin 2201.06602** | flat D = 4 | `|ĝ₃|²M⁸ ≤ 37.8 log(M/m_IR) − 45.4` (3.4a), `24.9 log − 27.6` (4.4) | derivation with IR assumption (this repo: certified positivity + eikonal hypothesis) | `log M/M⁸` | explicit (flat) |
| same, §5 remark | 3 | "(4.4) translates to a bound on three-point OPE coefficients of stress tensors in holographic 3d CFTs … `log(MR_AdS)`" | assertion, no formula | — | — |
| **Caron-Huot–Li–Parra-Martinez–Simmons-Duffin 2205.01495** | flat D ≥ 5; d = 4 | `|a−c|/c ≤ 23/Δ_gap² + O(Δ_gap⁻⁴)` | flat: derivation; uplift asserted ("rigorous … up to corrections suppressed by `1/Δ_gap`"), not proved for spinning correlators | `Δ⁻²` | explicit |
| Chiang–Huang–Li–Rodina–Weng 2201.07177 | flat D = 4 | `R³`, `D²ⁿR⁴` bounds assuming perturbative UV completion | derivation under assumption | — | explicit (flat; not checked) |
| Chang–Landau–Simmons-Duffin 2311.04271 | general d | spinning dispersive sum rules, bulk-point dictionary for `T` | derivation | no bound | — |
| This repo (`EINSTEIN_PROGRAM.md` §3a, `KERNELS_D3.md`) | 3 | superconvergence LP: `|c|² ≲ 5×10³/Δ_gap⁴`; cannot beat `2700/Δ_gap⁴` | numerics, conditional on (E) | `Δ⁻⁴` in `|c|²` | numerical |
| 2023–2026: Belin–Borak–Henriksson–Piron–Zhiboedov 2608.26246; Meçaj–Moult–Walters–Xin 2512.09986; Meçaj–Moult–Stoffels–Walters–Xin 2608.19322; Chen–Lü 2311.05307; Perlmutter 2610.24496; Häring–Zhiboedov 2410.21499; Dong 2610.19149; Bhat–Zahed 2304.02003 | various | energy-correlator positivity (`C₄^T = t₄²/16` in d = 3, no gap); holographic `t₂, t₄` for d ≥ 4; `Δ_gap ≲ (δc)^{1/d}`; graviton pole; IR-finite half-space positivity (scalars); scalar AdS locality | various | **none bounds CFT3 `⟨TTT⟩` by `Δ_gap`** | — |

**Plain statement.** A rigorous `Δ_gap`-dependent bound on `t₄` (or `α_T`) in CFT3 does **not** exist in the
literature found. The d = 3 statements are parametric estimates (AKT 1805.07393, CEMZ) or an unproved one-line uplift
remark (CLPS 2201.06602 §5). Even in d = 4 the only explicit-constant statement (2205.01495, for `a − c`) rests on an
asserted uplift.

*Citation corrections to the task brief:* arXiv:1203.2964 is a charmonium paper; the Chowdhury–Raju–Sachdev–Singh–Strack
paper is **1210.5247**, and it concerns `⟨TJJ⟩` (Weyl–Maxwell coupling `γ`, `|γ| ≤ 1/12`, free fermion/boson `γ = ±1/12`),
not `⟨TTT⟩`. arXiv:2201.07177 is Chiang–Huang–Li–Rodina–Weng; the higher-D CLPS paper is **2205.01495**. 1709.03597 is
AJHKT "Shockwaves from the OPE" (general d; no d = 3 `⟨TTT⟩` gap bound). Belin–Hofman–Mathys–Walters (2011.13862) is d = 4.

## 8. Assessment

What a rigorous (modulo stated hypotheses) AdS4/CFT3 bound would add:
1. The first `Δ_gap`-dependent bound on non-Einstein `⟨TTT⟩` data in **any** dimension with explicit constants and a
   proof beyond flat space; in d = 3 it would replace AKT's parametric `ln Δ_gap/Δ_gap⁴` by
   `|t₄ + iα_T| ≤ 1476√(log Δ_gap + O(1))/Δ_gap⁴` — sharper by `√ln Δ_gap` and with the constant.
2. Joint treatment of the parity-even and parity-odd structures (`t₄² + α_T²`), which is natural only in d = 3 (complex
   `ĝ₃`), matching the collider disc.
3. A first stress-tensor (spinning) instance of the CHMRS "AdS bulk locality from sharp CFT bounds" programme;
   CHMRS cover scalars only, CLS-D build the dictionary without bounds.

Closest prior result: **AKT 1805.07393 §7** (same object, same dimension, same `Δ_gap⁻⁴` and log, heuristic, no
constants); closest *quantitative* analogue: **2205.01495** (d = 4, `a − c`, explicit constant, uplift asserted).
The flat-space input (CLPS 2201.06602 and `FLAT_GRAVITON_BOUND.md`) is not new; the dictionary in §0 is assembled from
known pieces (its d = 3 coefficient `−240` is, as far as I found, not printed anywhere, but follows from CEMZ + CLPS and
is confirmed by the ECG result of BCR). Novelty therefore rests entirely on the AdS4 lift (chain links 6, 7, 9, 10),
which does not exist.

## 9. Conventions not verified
* The sign of `t₄` relative to `α₃` is fixed by matching to BCR and MPS (both agree with `t₄ = −240α₃/L⁴`), not by an
  independent sign computation of the shock vertex.
* The sign/phase of `α_T` relative to `α̃₃` (orientation, `ε^{0123}`, helicity labelling, CDP's `p_T` sign).
* AKT's `ñ_s, ñ_f, ñ_odd` normalisation; LMP's `λ_TTT` relation (quoted, not rederived).
* The QTG `t₂` coefficient (§3.3 caveat): not reproduced; irrelevant for d = 3.
* `κ` (the `O(1)` in the log) and the precise `M ↔ Δ_gap` identification at `O(1/Δ_gap)`.
* That the spinning AdS4 functional regulates the graviton pole as `1/(p² + 9/(4L²))` (true for CHMRS's scalar functionals).

## 10. Files
`ads4/code/dictionary/three_point_vertices.py` (flat vertices, `I₁ = 8A₂²`, `I₂ = 4A₂²`, duality, normalisations;
output `out_three_point_vertices.txt`); `ads4/code/dictionary/t4_dictionary.py` (CEMZ limits, GB/QTG/ECG checks, Weyl
identity, HM d = 3 constant, `C_T`, IR regulator, translation numbers; output `out_t4_dictionary.txt`; ≈ 5 s).
