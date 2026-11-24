# DSKETCH_UP: the two [D-sketch] links of the chain: `c₁` universality / ν-independence, and the action of light single traces (A1′)

*24 Nov 2026, research agent. Unreviewed, not committed. No existing file was edited; `AG_PROOF.md` and `code/agproof/` were not
touched. New code: `ads4/code/dsketch/` (`asd_alg.py`, `t1a_einstein2.py`, `t1a_solve2.py`, `t1a_lemmaM.py`, `t1a_boundary.py`,
`t1b_barnes.py`, `t2_block_growth.py`, `t2_tolerance.py`); outputs in `ads4/code/dsketch/cache/` (git-ignored by its own
`.gitignore`). Read-only imports: `exactskel/exactmom.py`, `uniform/skel.py`; data read: `prl/data/nf02_psi.json`,
`prl/data/q3_full_delta0.2.json`. Literature: CHMRS 2106.10274 (LaTeX source: eqs. `BMellinDefinition`, `CfromBMellin`, `Bkt action`,
`eftformula`, `moregeneralformula`, §`sec_light`); abstracts checked for GPT 1410.4185, Sleight–Taronna 1804.09334, CGP 1209.4355.
Compute ≈ 0.6 CPU-h, ≤ 2 processes, 1 thread, `ulimit -v 3000000`.*

Grades as in UNIFIED_THEOREM_v3 §0: [P] > [certified] > [DV] > [D] > [D-sketch] > [NE±] > [NE] > [NE + M] > [M] > [A] > [I]; [Q] quoted.

---

## 0. Summary

| link (UT3 / supplement row) | old grade | **new grade** | what was shown | residual |
|---|---|---|---|---|
| **L3d input: a regular second-order anti-self-dual (ASD) solution** (C1_VALUE §3.2 step 3) | [D-sketch] | **[DV]** (explicit, exact) + **[D]** (all momenta) | For two ASD legs, the `O(ab)` Einstein equation on Euclidean Poincaré AdS₄ has the **explicit, elementary, regular solution** `Γ_ij = P_ij(z)·e^{iK·x − (|k₃|+|k₄|)z}/z²`, `P_ij` **quadratic** polynomials, unique in that class, and its second-order self-dual Weyl tensor **vanishes identically** (exact sympy, two kinematic points, §1.2–1.3). General momenta: mode completeness of free SD spin-2 fields at fixed `K` (2-dimensional, `λ = ±|K|`, exact) + Bianchi ⇒ a regular ASD completion always exists (§1.4). | none for this input |
| near-boundary premise of the ASD lemma (`I_±` needs no counterterm: orthonormal Weyl `= O(z³)`) | [D-sketch] (inside C1V §3.1) | **[DV]** on the explicit configuration; [D] (FG) | conformal-frame Weyl `W̃ = O(z)` (sharp) at orders `a`, `b`, `ab` (§1.5) | — |
| **L3d ASD closure / `c₁` universal (fixed by `⟨TT⟩`, `⟨TTT⟩` alone)** | [D-sketch] (weakest link) | **[D]** (inputs [DV], [Q + DV]) | steps 1–4 of C1V §3.2 now have no sketched input; in addition the remnant is **exactly factorised**, `Φ₃₄ = −b₊(k₃,k₄)·W₊[regular label-(+) mode at K]`, with `b₊` the label-(+) TT boundary value of the explicit ASD solution, **non-zero** (`b₊ = 7/300` at `k₃=(0,0,3)`, `k₄=(4,0,0)`) | value of `c₁` still [I] (not attempted) |
| **L3a: `J₀(E×W³) = 2` in the preserved sector** | [D-sketch] | **[D]** (inputs [Q], [P]) | written out (§2.2): Mellin degree ≤ 3; top degree = flat amplitude, whose preserved part is 0 [P]; no `u`-channel pole in the preserved sector [P/DV]; so every surviving term is a polynomial in `u` of degree ≤ 2 ⇒ growth ≤ 2 | the spinning-external Mellin analytic structure is **quoted** (GPT, Sleight–Taronna, CGP), not re-derived for `⟨TTTT⟩` |
| **L3a: polynomiality / ν-independence ("no ν^{2n}, n ≥ 2")** | [D given (E)] | **[D]**; Mellin→ν map **[P + DV]** | (E) is not needed for this term: polynomial top-Regge coefficient in Mellin `T` ⇒ polynomial in ν² of the same degree, **exactly** (Barnes' first lemma; `T ↦ (A²+ν²)/2A`; checked to 10⁻³⁰, §2.3); degree counting then leaves exactly a constant (§2.4) | as above |
| **LS4: at k = 2, 3 a light J ≤ 2 state contributes its block action** | [D-sketch] | **[D-sketch]**, resisting step now one named input (E′) | Exact **form** of the light value derived (§3.3): a polynomial in ν² of degree `n − k` whose top coefficient is the flat (CHLPS) one **exactly** and whose lower coefficients carry `O(1)` AdS-unit shifts; for J = 0 the whole AdS effect is one number `δ(Δ_L)`: `m² → μ² = m²Λ² + δ`. | (E′): CHMRS's light-side theorem (`moregeneralformula`) for the spinning `k = 2, 3` functionals on tree terms of Regge spin `n = 3, 4 > k` (equivalently: no double-trace subtraction poles in the per-ν spinning kernel) |
| **NG2 / UT3 §5.5 "LR form = true AdS image at Δ_L = O(1)"** | [D-sketch] | **REFUTED as stated** [DV, scalar CHMRS representation] | **Negative finding (§3.2):** the exact per-ν action of one block grows like `e^{cν}`, `c → π/2`, for ν above its mass (10¹⁰–10²⁸ at ν = 30–60 for Δ = 4.5–12), contour- and step-independent to 10⁻¹². The leading-Regge/hybrid "image" (flat + `P − J₀`) is not the AdS image of a light block at any Δ_L ≪ Λ. What survives is the Regge-polynomial light value of §3.3 (≠ LR form; its non-flat part is a constant, not `P − J₀`). | as LS4 |
| **LS6/LS7: `Ω^Q_up + 0.040 Ω^EE` positive on light scalars** | [D-sketch] | **[D-sketch]**, reduced to one inequality | positive **iff** `δ(Δ_L) > −692.5 − m²Λ²` (Δ_g = 300) / `−7748.7 − m²Λ²` (1000), AdS units (§3.4); the Mellin→ν map alone gives `|δ| ≈ A²/2 ≈ 21` (k = 2 coefficient) [M] ⇒ margin ×33 (300) / ×370 (1000) | `δ(Δ_L)` of the spinning `φW²` exchange not computed |
| ε* = 0.038034 repair | [I] | **[I], now with a definite adverse estimate** | positive iff `μ² ≥ 0`; the Mellin→ν shift alone is **negative** (`−A²/2`), so ε* is expected to **fail** on light scalars with `Δ_L ≲ 5–6` [M]; ε* must not be used (0.040 already replaces it) | — |
| Ω^EE on light J = 0, 2 (evidence for (Pos), not a link) | [D-sketch] (LS6) | [D given `|δ| < 4.9·10⁴ / 5.4·10⁵` (J = 0), `K < 3.2·10³ / 3.5·10⁴` (J = 2)] | tolerances ≥ 10² × the natural shift (§3.4) | as LS4 |

**Net effect on the theorem.**
* **Item (1) without light single traces** (A1, or A1′ with no light J ≤ 2 coupled to TT): the only [D-sketch] link (L3) is removed. **"(Pos) ⇒ (1)" becomes [D]** (H17 band stays [D-sketch + M] and only moves the constant by ±0.24 % / ±0.021 %).
* **Under A1′ with light single traces present:** item (1) and item (2) keep a [D-sketch] link, now precisely (E′) (LS4), and (2b) additionally the single inequality on `δ(Δ_L)` (LS6).
* **Negative findings, prominently:** (i) the leading-Regge light-state images of NEARGAP2 §2.2/§2.4 and A1PRIME §4.2 ("flat × (1 + O(Δ_L²/Λ²))") are **not** the exact block action in any CHMRS-type exact representation tested; their numbers remain meaningful only as the *flat* values, which happen to coincide with the top coefficients of the exact Regge polynomial; (ii) ε* is expected to fail; (iii) a caution outside scope: per-ν block actions grow exponentially once ν exceeds the block's mass, so any weight that reaches ν > m for a near-gap heavy state (e.g. de-shift tails up to Λ + 100 for states with `Λ ≤ m < Λ + 100`) should be re-checked against this growth (not examined here).

---

## 1. Target 1(a): the regular second-order ASD solution, constructed

### 1.1 What C1_VALUE §3.2 needs
Legs 3, 4 are linearised ASD gravitons (regular label-(−) modes, C1V `v2`). Step 3 assumed a regular second-order solution `δh^{ASD}₃₄`
of the `O(h₃h₄)` Einstein equation whose second-order SD Weyl tensor vanishes. With it, `χ = δh^W − δh^{ASD}` is a regular linear
solution and `Φ₃₄ = W₊[χ]`, which fixes the remnant through the boundary data `b₊` of `χ`.

### 1.2 Construction (`t1a_einstein2.py`, `t1a_solve2.py`)
* Metric `g = z⁻²(δ + a γ₃ + b γ₄ + ab Γ)`, `γ_{3,4} = ε^{(−)}ε^{(−)}(1+|k|z)e^{ik·x−|k|z}` (TT, radial gauge), `Γ` radial gauge
  (`Γ_{zμ} = 0`) with x-dependence `e^{iK·x}`, `K = k₃ + k₄`. A truncated two-parameter algebra (`asd_alg.py`) computes Christoffels,
  Riemann, Ricci, Weyl to `O(ab)` exactly; `E_μν = R_μν + 3g_μν`.
* Orders 1, `a`, `b` vanish identically (checks the legs).
* The `O(ab)` sources are polynomials in `z` times `e^{−(|k₃|+|k₄|)z}` — no `1/z` poles. The ansatz `Γ_ij = P_ij(z)e^{−(|k₃|+|k₄|)z}`
  (deg ≤ 5, complex rational coefficients) turns all ten equations into a linear system; **it has exactly one solution, and it is
  quadratic**:

  `k₃ = (0,0,3)`, `k₄ = (4,0,0)`, `K = (4,0,3)` (index order x₁x₂x₃):
  ```
  Γ₁₁ = −33(3234z² + 3493z + 499)/9604      Γ₁₂ =  i(4410z² + 4403z + 629)/196      Γ₁₃ = 2(23373z² + 22708z + 3244)/2401
  Γ₂₂ = (21462z² + 21427z + 3061)/588        Γ₂₃ = −2i(1911z² + 1988z + 284)/147     Γ₃₃ = −40(2940z² + 3367z + 481)/7203
  (all × e^{−7z}; h = Γ/z²)
  ```
  The second point `k₃ = (1,2,2)`, `k₄ = (2,−1,2)`, `K = (3,1,4)`, `|K| = √26` (non-orthogonal legs, equal energies) gives again a
  unique quadratic `P_ij·e^{−6z}` (`cache/t1a_solve2_1_2_2_2_-1_2.txt`).
* **Interpretation.** In the conformal frame `g̃ = δ + …` the solution is the AdS image of the flat ASD perturbiner: the product of two
  complex-null plane waves `e^{i(k·x) − |k|z}` has "momentum" `(K, i(|k₃|+|k₄|))`, exactly as the second-order ASD (self-dual
  gravity) perturbiner in flat space. The polynomial prefactor is the curvature correction.

### 1.3 Verification [DV] (`cache/t1a_solve2_*.txt`)

| check | (0,0,3)+(4,0,0) | (1,2,2)+(2,−1,2) |
|---|---|---|
| `O(ab)` Einstein equations (all 10 components), exact | ✓ | ✓ |
| `W₊` at orders `a`, `b` (ASD legs) | 0 | 0 |
| `W₋` at orders `a`, `b` | 64 non-zero comps | 144 |
| **`W₊` at order `ab`, projected with the full metric's Hodge star** | **0 (all 256)** | see `t1a_solve2_1_2_2_2_-1_2.txt` |
| `W₋` at order `ab` | 144 non-zero | — |
| control: `Γ + (regular label-(+) mode at K)` (still a solution) | `W₊^{(ab)}` 144 non-zero | — |

The Hodge star of the **full** metric is used (`√g`, `g^{-1}` to `O(ab)`): with ASD legs, `P₊[g]W₊ = W₊` at `O(ab)` forces the
second-order SD part to be background-SD, so "W₊^{(ab)} = 0" is the frame-independent statement `Ψ^{(2)} = 0`.

### 1.4 All momenta [D + DV] (`t1a_lemmaM.py`, `cache/t1a_lemmaM.txt`)
1. **Lemma M (mode completeness) [DV].** A free SD spin-2 field on EAdS₄ with x-dependence `e^{iK·x}` is, after the conformal rescaling
   `ψ_flat = z⁻¹ψ_AdS` (Penrose–Rindler I §5.7), a flat SD field `Σ φ_ab S^a ⊗ S^b` (`φ` symmetric traceless, 5 components) obeying
   `∂^μΦ_μνρσ = 0`, a constant-coefficient system `M(∂_z)φ = 0` (48 equations). Its solution space has dimension
   `deg gcd(5×5 minors of M(λ)) = deg(λ² − |K|²) = 2` (generic rank 5; rank 4 at `λ = ±|K|`), checked exactly for `K = (4,0,3)` and
   `(3,1,4)`; by rotation and scaling covariance it holds for every `K ≠ 0`.
2. The two solutions are `W₊` of the **regular label-(+)** mode (`z⁻¹W₊ = e^{−|K|z}·const`) and of the **irregular label-(−)** mode
   (`e^{+|K|z}`), checked exactly (both `K`), each traceless and in `ker M(∓|K|)`.
3. **Construction for any momenta.** Let `δh_p` be any regular (decaying) particular solution of the `O(ab)` Einstein equation (it exists:
   the sources decay like `e^{−(|k₃|+|k₄|)z}` and `|K| < |k₃|+|k₄|` for non-collinear legs; variation of parameters with the homogeneous
   modes `(1 ± |K|z)e^{∓|K|z}`). By the second-order Bianchi identity (C1V §3.2 step 2, [D]) `Φ_p := W₊^{(ab)}[δh_p]` is a free SD
   field; it decays, so by Lemma M `Φ_p = c·W₊[reg. (+) at K]`; then **`δh^{ASD} := δh_p − c·h^{reg,+}_K` is regular, solves the
   `O(ab)` Einstein equation and has `Ψ^{(2)} = 0`.** It is unique up to regular ASD linear solutions (label-(−) mode at K, gauge).
   The explicit cases of §1.2 have `c = 0` already.
4. **Explicit remnant.** For the normalisable solution `δh^W = δh^{ASD} + χ`, `χ` is the regular linear solution whose boundary value
   cancels that of `δh^{ASD}` (TT part plus boundary diffeomorphism/Weyl parts, which carry no Weyl tensor). Hence, exactly,
   `Φ₃₄ = −b₊·W₊[h^{reg,+}_K]`, `b₊ :=` label-(+) TT projection of `Γ(z = 0)`; gauge-invariant under boundary diffeomorphisms and Weyl
   rescalings. At the first point **`b₊ = 7/300 ≠ 0`** (`b₋ = 126/25`), confirming C1V's JJ21-based "the remnant is generically
   non-zero" by an independent route.

### 1.5 The near-boundary premise [DV] (`t1a_boundary.py`)
C1V §3.1 used that `I_± = ∫√g TrW_±³` needs no boundary counterterm because the orthonormal-frame Weyl tensor is `O(z³)` [D-sketch].
Orthonormal components are `z²W̃_μνρσ[g̃]`. On the explicit configuration every component of `W̃` vanishes at `z = 0` at orders `a`, `b`,
`ab`, and has a non-zero `O(z)` term (72/72/178 components; 64/64/144 with `O(z)`): `W = O(z³)`, sharp. (General statement: Fefferman–
Graham, `W ~ z³(Cotton[g₀], T)`, [D].)

### 1.6 Consequence for L3d / universality
All four steps of C1V §3.2 now have written-out inputs: step 1 (partition) [D]; step 2 (Bianchi) [D] (and verified trivially here: every
`O(ab)` solution is `δh^{ASD}` + a linear solution); **step 3 [DV]/[D]** (this section); step 4 (boundary relation) [D] with the JJ21 check
[Q + DV]. Hence **"`c₁` is fixed by `⟨TT⟩`, `⟨TTT⟩` data alone" is [D]**. The value of `c₁` was not attempted; F₂ stays [A].

---

## 2. Target 1(b): `J₀(E×W³) = 2` and polynomiality, written out

### 2.1 Inputs (quoted, not re-derived)
* **(M1) Analytic structure** [Q]: a tree-level `O(1/c_T)` Mellin amplitude is meromorphic, with poles only at the twists (and satellites)
  of exchanged single traces, residues polynomial in the other Mellin variable, plus an entire polynomial part (contacts); its large-Mellin
  asymptotics is polynomially bounded with degree fixed by the derivative count (Penedones 1011.1485; GPT 1410.4185 for one spinning
  operator; Sleight–Taronna 1804.09334 for spinning external legs, where the external-spin dependence factorises).
* **(M2) Regge limit** [Q]: the `u`-Regge limit of the correlator is `s → ∞` at fixed Mellin `u`, and the growth exponent of the correlator
  equals that of the Mellin amplitude (CGP 1209.4355).
* **(M3) Flat limit** [Q]: the top homogeneous part of the Mellin amplitude is the flat amplitude (Penedones; GPT §flat limit for spin).
* **(M4) The per-ν map** (CHMRS `CfromBMellin`): `C_{k,ν} ∝ γ_a(ν)⁻²∫dT Γ(T/2)²γ_{a−T}(ν)B̂_{k,T}`; for the spinning functional the same
  map acts on each helicity-preserving plane-wave pair (D3 [P given H7 §18.1]; heavy and light share `γ_a²M²`, H7 Q5).

### 2.2 Derivation of `J₀ ≤ 2` [D]
Count the full helicity amplitude (helicity factor included, as in A1PRIME §1: `n = 4 + growth of f`). Einstein `⟨TTTT⟩` has Mandelstam
degree 1 (2 derivatives); the `O(ĝ₃)` tree has degree 3 (2 + 6 − 2 derivatives).
1. **Top degree (3) = flat amplitude (M3).** At `O(ĝ₃)` the flat four-graviton amplitudes are `M(−+++)`, `M(++++)` (both Regge spin 2)
   and **`M(−−++) = 0`** [P]: every factorisation channel of `(−−++)` at `O(ĝ₃)` needs a `W³` vertex (all-same helicity) on one side and
   an Einstein vertex (mixed) on the other, which is impossible with two `−` and two `+` external legs, so `M` is a polynomial with
   little-group weight `⟨12⟩⁴[34]⁴` (dimension 8) but total dimension 6: zero. So the preserved (MHV) top-degree part vanishes.
2. **No `u`-channel pole in the preserved sector** [P + DV, C1 §3.1, C1V §2.2]: the `u`-channel `T` exchange at `O(t₄)` is pure flip on
   one line, on-shell (HM one-point matrix) and off-shell (Weyl factorisation; CHP preserved elements `¼□(□+ν²+¼)`, `½(□+ν²+¼)`
   vanish on every H₂ eigenfunction).
3. **Lower degrees (≤ 2).** By (M1) and step 2, every remaining term is a sum of `s`- and `t`-channel pole terms and polynomials; at large
   `s` and fixed `u` (`t = Σ − s − u`) each has a Laurent expansion whose coefficients are **polynomials in `u`**, and growth ≤ its degree
   ≤ 2. Hence **`J₀ ≤ 2`** in the preserved sector (and `= 2` in the flip sectors, which `EE` annihilates).
   The ASD closure (§1) is consistent: the preserved remnant is an energy-channel (`s`, `t`) exchange with no total-energy pole, i.e. no
   top-degree part.

### 2.3 Polynomiality [P + DV] (`t1b_barnes.py`, `cache/t1b_barnes.txt`)
With `z = T/2`, `c, d = (A ± iν)/2`, Barnes' first lemma gives **exactly**
`C[(T/2)_j]/C[1] = (c)_j(d)_j/(A)_j = Π_{i<j}((A/2+i)² + ν²/4)/(A)_j`,
so `T⁰ ↦ 1`, `T ↦ (A²+ν²)/(2A)`, `T² ↦ (A²+ν²)²/(4A(A+1))`, `T³ ↦ (A²+ν²)²(A²+ν²−4)/(8A(A+1)(A+2))`: a polynomial top-Regge coefficient
of degree `n` in Mellin `T` maps to a polynomial in ν² of degree `n`, **with no non-polynomial remainder**. Checked by direct numerical
contour integration (mpmath, 30 digits) at `A = 4.5, 6.5, 7.5`, `ν = 0.7, 3, 11`, `j ≤ 3`: agreement ≤ 4·10⁻³¹ (one point 2·10⁻²³);
symbolic images exact. **Input (E) is therefore not needed for the `O(t₄)` term**: its only non-polynomial source would be a `u`-pole,
excluded by step 2.

### 2.4 ν-independence [D]
At `O(t₄)` (`t₄ ∝ ĝ₃`, AdS units `L = 1`) the preserved top-Regge (`s²`) coefficient has Mandelstam degree `3 − 2 = 1` at the flat order —
which vanishes (step 1) — and degree `0` at `O(L⁻²)`: a **constant**. Degree `−1` (`O(L⁻⁴)`) would need a `u⁻¹`, i.e. a `u`-pole: absent.
By §2.3 the per-ν image is a constant: **`−Ĉ′_EE ⊃ 16πGc₁t₄`, ν- and Δ_g-independent** (equivalently C1 §3.5's flat counting: `ν^{2n}`,
`n ≥ 2`, diverges in the flat limit; `n = 1` is the vanishing flat MHV term). Grade **[D]**, inputs (M1)–(M4) [Q]/[P given H7].

### 2.5 What is quoted rather than derived
(M1) for four spinning external `T`'s (degree bound and polynomial residues) is taken from GPT/Sleight–Taronna/CGP; it was not
re-derived for `⟨TTTT⟩` here. If a referee asks for more, the explicit `O(ĝ₃)` preserved remnant of §1.4 (`b₊ × ⟨TTT⟩_{W₊}`) is the
object whose Regge limit would have to be taken directly.

---

## 3. Target 2: the action of light single traces

### 3.1 Claims under test
* **LS4** (A1P §3): at k = 2, 3 a light J ≤ 2 state contributes its block action `Tr[P_L M_Ω(O_L)]`, like a heavy state.
* **LS6 / NG2 §2.2, §2.4**: that action is the flat value × `(1 + O(Δ_L²/Λ²))`; evaluated with the heavy hybrid "flat + leading-Regge
  curvature `P_{−½+iν}(cosh β) − J₀(νB)`", giving `F/N ≥ 0.788 / 0.638` (EE) and the ε = 0.040 margin `0.022N`.

### 3.2 Negative finding: the per-ν action of one block above its mass [DV, scalar CHMRS representation] (`t2_block_growth.py`)
The only exact per-ν representation of a block action available in the project is CHMRS's Mellin formula (scalar externals,
`exactskel/exactmom.py`, validated there to 10⁻¹² against the CHMRS series for heavy states). Evaluated at ν **above** the block's
mass (Δφ = 3; heavy density `m^{2k}C_{k,ν}`):

| Δ, J | m | k | ν = 10 | ν = 20 | ν = 30 | ν = 60 | `d ln|C|/dν` (ν = 40→60) | flat `C_BP` at ν = 30 |
|---|---|---|---|---|---|---|---|---|
| 40, 0 | 38.5 | 2 | 2.222 (flat 2.222) | 3.163 (3.246) | 7.36 (9.03) | 7.3·10⁶ | 0.70 (onset) | 9.03 |
| 20, 0 | 18.5 | 2 | 3.13 | 24.1 | 6.3·10³ | 5.1·10¹⁶ | 1.06–1.12 | −0.24 |
| 12, 2 | 10.2 | 2 | 1.44 | 1.7·10⁴ | 3.3·10⁸ | 3.5·10²⁴ | 1.24–1.30 | −45 |
| 10, 0 | 8.5 | 2 | 11.4 | 5.6·10³ | 9.8·10⁷ | 9.2·10²³ | 1.24–1.30 | −0.08 |
| 4.5, 0 | 2.96 | 2 | 7.63 | 1.9·10⁵ | 3.2·10¹⁰ | 1.4·10²⁸ | 1.36–1.40 | −0.01 |
| 4.5, 0 | 2.96 | 3 | 1.36 | 2.6·10⁴ | 3.3·10⁹ | 8.0·10²⁶ | 1.34–1.38 | — |

* Below the mass the exact action approaches the flat bulk-point form (Δ = 40: 2.2223 vs 2.2223 at ν = 10). **Above it, it grows like
  `e^{cν}` with `c → π/2`**: the `γ_a(ν)⁻²` normalisation (`~ e^{πν}`) is no longer compensated, because a single block's Regge moment
  decays only like `e^{−πν/2}` beyond its mass.
* **Not an artefact:** five contour positions (`Re T = 1.5 … 9.5`, with and without the residues at `A ± iν + 2j`) and two line steps agree
  to ≤ 2.5·10⁻¹² (`cache/t2_block_growth.txt`).
* **Consequence.** A light state samples ν ∈ (m_L, Λ], i.e. this regime. Neither "flat × (1 + O(Δ_L²/Λ²))" nor the leading-Regge form
  `P_{−½+iν}(cosh β)` (bounded, oscillating) describes the exact block action there. **The LR form is not the AdS image of a light block,
  at Δ_L = O(1) or at any Δ_L ≪ Λ** — in the only exact representation tested. NG2's light-state numbers are evidence about the flat
  values only.

### 3.3 What the sum rule actually needs, and its exact form [D]
* **Mechanism (CHMRS representation) [D-sketch].** `B̂_{k,T}[M] = (−1)^{k/2}∮_∞ ds′ M(s′,T′)/[(2Δφ−s′)(Δφ−s′/2)_{k/2}(Δφ−T′/2)_{k/2}]`.
  Cauchy: (arc at ∞) = (residues at the poles of M — the single-trace block action, CHMRS `Bkt action`) + (residues at the kernel's
  subtraction poles, i.e. at double-trace points — the double-trace contributions that the subtracted `k = 2` kernel does **not**
  annihilate: at `[φφ]₀` the subtraction factors cancel the dDisc double zero, so it is not a double zero). For a light exchange the arc is polynomial in `T` (Regge data) ⇒ polynomial in ν² (§2.3);
  the block's `e^{πν/2}` growth is carried by the `T`-poles at `T = τ + 2n + 2j` that the subtraction factors place left of the saddle
  `T* ≈ ν²/2A`, and is cancelled by the double-trace terms.
* **The light value that enters the sum rule** is the term-wise sum over the light exchange's OPE content (single trace + its double-trace
  dressing). By CHMRS's theorem (app. `omega_light`, eq. `moregeneralformula`, [Q] for scalar externals) it equals the arc formula:
  `−ω|_light = ¼Σ_{k′=k}^{n}∫dη∮dr r^{k′−2}a_{k′}(η)𝒢₊^{X_L}`. For the project's **unsubtracted** k = 2, 3 functionals with `[TT]` at
  double zeros (L1), this light value **is** the block action (LS4), *provided* the spinning carry-over of that theorem holds for
  `n_X = 3, 4 > k` — call it **(E′)**. (E′) is the step that resists: the spinning per-ν kernel in Mellin form is not available
  (SPINNING_MELLIN §5), and LIGHTSIDE's arc identity was established only for `n_X ≤ k`.
* **Exact form given (E′) [D].** In the MHV sector a light J ≤ 2 state has poles only in the `t = (14)` channel (no `+−` structure for
  J < 4, A1P §2.1), so at fixed `u` all Laurent coefficients are polynomials in `u` (M1), hence in ν² (§2.3). Degree counting as in §2.4:
  * **J = 0** (`n = 3`): `k = 3` value = top coefficient = **exactly the flat `B₃ = |g₀|²`**; `k = 2` value `= −|g₀|²(ν² − 2μ²)/Λ²`
    with the ν² coefficient exactly flat and **one AdS number** `μ² = m²Λ² + δ(Δ_L)` (AdS units). So the full AdS effect on a light scalar
    is `m² → μ²/Λ²` in A1P's exact polynomials.
  * **J = 2** (`n = 4`): `k = 2`: degree 2 in ν², `k = 3`: degree 1; top coefficients exactly flat (CHLPS `−6p⁶`, `6p⁴`), lower ones shifted
    by `O(1)` AdS units.
  * The NG2 "curvature" `P − J₀` is **not** part of this light value.

### 3.4 Positivity with explicit margins (`t2_tolerance.py`, `cache/t2_tolerance.txt`)

**J = 0** (exact rationals from `nf02_psi.json`, `q3_full_delta0.2.json`; Λ = 298.4996 / 998.4999):

| functional | AdS light value per `|g₀|²` | positive iff | Δ_g = 300 | Δ_g = 1000 |
|---|---|---|---|---|
| `Ω^EE` | `11.317027 + 20.771893 μ²/Λ²` | `μ² > −0.5448Λ²` | `μ² > −48545` | `μ² > −543191` |
| `Ω^Q_up` | `−0.430430 + 2.032073 μ²/Λ²` | `μ² > +0.2118Λ²` | `μ² > 18873` | `μ² > 211183` |
| **`Ω^Q_up + 0.040 Ω^EE`** | **`0.0222508 + 2.862949 μ²/Λ²`** | **`μ² > −0.007772Λ²`** | **`μ² > −692.5`** | **`μ² > −7748.7`** |
| `Ω^Q_up + ε*Ω^EE` | `0.0000000 + 2.822108 μ²/Λ²` | `μ² ≥ 0` | `μ² ≥ 0` | `μ² ≥ 0` |

with `μ² = (Δ_L − 2)(Δ_L − 1) + δ(Δ_L)` (heavy-state mass convention plus the AdS shift). Over the whole unitary range `Δ_L ≥ ½`,
`(Δ_L−2)(Δ_L−1) ≥ −¼`.

**J = 2** (shifts `|s_j|, |t_j| ≤ κ(m² + κ)^{d_j−1}`, `κ = K/Λ²`, on every lower coefficient, worst sign):

| functional | min flat `F²` on `0 < m ≤ Λ` | survives every shift with `K ≤` (300 / 1000, AdS units) |
|---|---|---|
| `Ω^EE` | 0.4481 (m = 0.675Λ) | 3.16·10³ / 3.54·10⁴ |
| `Ω^Q_up` | 1.9884 (m = 0.976Λ)ᵃ | 5.13·10⁴ / 5.74·10⁵ |
| `Ω^Q_up + 0.040 Ω^EE` | 2.1916 | 4.28·10⁴ / 4.79·10⁵ |

ᵃ A1PRIME prints "min 1.97"; its printed polynomial `2.637798 + 0.480850m² − 3.156162m⁴ + 2.032073m⁶` has minimum 1.988 at
`m = 0.976Λ`. Immaterial.

**Natural size of `δ` [M].** The Mellin→ν map alone (§2.3) turns a Mellin-linear coefficient into `ν² + A²`, `A = a − ½ = 6.5` (k = 2),
`7.5` (k = 3) for `Δ_T = 3`: a shift of `μ²` by `−A²/2 ≈ −21` (k = 2). The exchange's own satellite-pole average adds `O(Δ_L)` terms,
not computed. Against this:
* **`Ω^EE` (J = 0, 2): positive for every Δ_L with margins ≥ ×2000 (J = 0) / ×75 (J = 2, 300) in `δ`.**
* **ε = 0.040 (J = 0): positive iff `μ² = (Δ_L−2)(Δ_L−1) + δ > −692.5` (300) / `−7748.7` (1000); margin ×33 / ×370 over `A²/2`. Not refuted; not proved: the
  spinning `δ(Δ_L)` is the one missing number.**
* **ε\*: expected to fail.** It needs `μ² ≥ 0`; the Mellin→ν shift alone is negative, so light scalars with `(Δ_L−2)(Δ_L−1) < A²/2`
  (`Δ_L ≲ 6`) would make `Ω^Q_up + ε*Ω^EE` negative [M]. NG2's switch to 0.040 is therefore necessary, not cosmetic.

### 3.5 Verdicts for Target 2
* **LS4**: stays **[D-sketch]**; the gap is now exactly (E′) (spinning carry-over of CHMRS's light-side theorem to Regge spins
  `n = 3, 4 > k`, equivalently "no double-trace subtraction poles in the spinning per-ν kernel"). The **form** of the light value is [D]
  given (E′).
* **"LR form = AdS image" (UT3 §5.5, NG2 §2.2)**: **refuted** in the exact CHMRS representation [DV]; replace by "the exact light value is
  the Regge polynomial of §3.3; its flat top coefficients are the A1P polynomials; the AdS part is a constant shift `δ`".
* **LS6/LS7 (ε = 0.040)**: [D-sketch], reduced to `δ(Δ_L) > −692.5 / −7748.7` (one number per Δ_L; natural size ~21, [M]).
* Ω^EE light positivity (part of (Pos), evidence): robust (tolerances 10²–10⁴ × natural).
* **What happens to NG2's numbers.** Under (E′) they stand as *flat* values: the exact light value is A1P's flat polynomial with
  `m² → μ²/Λ²`. For a natural `δ ≈ −21` the shift moves `F_D⁰` by `20.77·δ/Λ² = −4.9·10⁻³` (300; 4·10⁻⁴ of `F`) and the ε = 0.040
  combination by `−6.8·10⁻⁴` (300) / `−6.0·10⁻⁵` (1000) against its margin `0.0223` (3 % / 0.3 %). The NG2 "LR curvature" column
  (≤ 1.8·10⁻⁴N) is not part of the exact value and should be dropped. Without (E′) — e.g. in any representation where the light
  block's kernel has double-trace subtraction poles — the bare block action is exponentially far from these numbers (§3.2), and only the
  total light value (block + its double-trace dressing) is the Regge polynomial.

---

## 4. Proposed replacement text

**Supplement Table (inputs), row L3:**
> L3 & $O(t_4)$ light term $=16\pi Gc_1t_4$, $c_1$ universal (fixed by $\langle TT\rangle$, $\langle TTT\rangle$); shock ($u$-channel) and
> $k=3$ parts $=0$; $c_1=0$ given premise F$_2$ & D, Q & form [D] (inputs: Mellin analytic structure [Q]; Barnes map [P\,+\,DV]; no
> preserved $u$-pole [P\,+\,DV]); universality [D] (ASD closure; the regular second-order ASD solution is explicit and checked exactly
> [DV]); shock part [P]; $c_1=0$ [D-sketch given F$_2$], F$_2$ [A]; value otherwise [I], carried exactly\\

**Row LS4:**
> LS4 & at $k=2,3$ the light value of a light $J\le2$ exchange is a polynomial in $\nu^2$ of degree $n-k$ (Regge data): top
> coefficients flat, lower ones shifted by $O(1)$ in AdS units ($J=0$: $m^2\to\mu^2/\Lambda^2$); it equals the block action given
> $[TT]$ double zeros & D, Q & form [D given (E$'$)]; (E$'$) (CHMRS's light-side formula for spinning $n=3,4>k$) [D-sketch]\\

**Row LS6 (and the "AdS images" sentence of §light):**
> LS6 & $\Omega^Q_{\rm up}+0.040\,\Omega^{EE}$ is positive on every light scalar iff $\mu^2>-692.5$ ($300$) / $-7748.7$ ($1000$) (AdS
> units); $\Omega^{EE}$ iff $\mu^2>-0.545\Lambda^2$ & Q & [D-sketch]: the AdS shift of $\mu^2$ is not computed (natural size $\approx21$)\\

> *Replace* "AdS images of light states: flat value $\times(1+O(\Delta_L^2/\Lambda^2))$ for $\Delta_L\ll\Lambda$ [D-sketch]" *by*
> "The exact per-$\nu$ action of a single block grows like $e^{\pi\nu/2}$ above its mass in CHMRS's Mellin representation, so light
> states are not treated as heavy-type images: their light value is the Regge polynomial of their exchange, whose top coefficients are
> the flat values above and whose only other effect is an $O(1)$ shift $m^2\to\mu^2/\Lambda^2$ [D given (E$'$)]."

**§sec:c1, sentence (c):** replace "(one input, a regular second-order solution, [D-sketch])" by "(the regular second-order ASD solution
exists and is elementary, $\Gamma_{ij}=P_{ij}(z)e^{iK\cdot x-(|k_3|+|k_4|)z}/z^2$ with quadratic $P_{ij}$; checked exactly [DV])", and
"[D-sketch, by the weakest-link rule]" by "[D]"; and "[D given $J_0(E\times W^3)=2$, D-sketch]" by "[D]".

**Net-grade sentences:** "Implication (1): [D] if no light single trace couples to $TT$; otherwise [D-sketch] through (E$'$) (LS4).
(2b) adds the single inequality LS6."

**UT3:** §0 `c₁` row: "[D-sketch] (weakest link …)" → "[D]"; §4.3 L3a/L3d grades → [D]; §4.7 LS4 → "[D-sketch]: (E′)"; LS6 → as above;
§5.5 row "LR form = true AdS image at Δ_L = O(1)": "[D-sketch]" → "**refuted** [DV, CHMRS repr.]; replaced by the Regge polynomial (§3.3
here)".

---

## 5. Reproduce

All on one core, `ulimit -v 3000000; export OMP_NUM_THREADS=1`; outputs in `ads4/code/dsketch/cache/`.
```
cd ads4/code/dsketch
python3 -B t1a_einstein2.py 0 0 3 4 0 0      # O(ab) Einstein equations (seconds)
python3 -B t1a_solve2.py 0_0_3_4_0_0          # explicit ASD solution + W+ checks (~1 min)  -> t1a_solve2_0_0_3_4_0_0.txt
python3 -B t1a_einstein2.py 1 2 2 2 -1 2; python3 -B t1a_solve2.py 1_2_2_2_-1_2     # second point (~20-30 min, sqrt(26) basis)
python3 -B t1a_lemmaM.py                      # Lemma M, mode identification, b_+ (~5 min)        -> t1a_lemmaM.txt
python3 -B t1a_boundary.py 0_0_3_4_0_0        # Weyl = O(z^3) orthonormal (~2 min)                -> t1a_boundary_*.txt
python3 -B t1b_barnes.py                      # Barnes map, numerical + symbolic (~3 min)          -> t1b_barnes.txt
python3 -B t2_block_growth.py                 # exact per-nu block action above the mass (~2 min) -> t2_block_growth.txt
python3 -B t2_tolerance.py                    # positivity tolerances (seconds)                    -> t2_tolerance.txt
```
Every script prints `ALL CHECKS PASS: True` where it has checks; all did.

## 6. Grade table

| statement | grade |
|---|---|
| explicit regular second-order ASD solution (two kinematic points), exact Einstein + `W₊^{(ab)} = 0` | **[DV]** |
| existence for all non-collinear momenta (Lemma M + Bianchi + regular particular solution) | **[D]**; Lemma M **[DV]** |
| remnant `Φ₃₄ = −b₊W₊[reg (+) at K]`, `b₊ ≠ 0` | [D]; `b₊ = 7/300` [DV] |
| orthonormal Weyl `O(z³)` (no counterterm) | [DV] on the configuration; [D] |
| ASD closure; `c₁` universal | **[D]** |
| `J₀(E×W³) ≤ 2` in the preserved sector | **[D]** (inputs [Q], [P]) |
| polynomial Mellin coefficient ↦ polynomial in ν² (Barnes) | **[P + DV]** |
| `O(t₄)` light term `= 16πGc₁t₄`, ν- and Δ_g-independent | **[D]** |
| per-ν action of one block grows like `e^{(π/2)ν}` above its mass (CHMRS kernels) | **[DV]** (5 states, contour/step independent to 10⁻¹²) |
| LR / hybrid form = exact light-block image | **refuted** in the CHMRS representation |
| light value of a light J ≤ 2 exchange = Regge polynomial (form, flat top coefficients, constant shifts) | [D given (E′)] |
| (E′) (spinning carry-over of CHMRS's light formula for `n = 3, 4 > k`) | [D-sketch] — the residual of Target 2 |
| `Ω^Q_up + 0.040Ω^EE` > 0 on light scalars | [D-sketch]; iff `δ > −692.5 / −7748.7`; natural `|δ| ≈ 21` [M] |
| ε\* repair | expected to fail [M] |
