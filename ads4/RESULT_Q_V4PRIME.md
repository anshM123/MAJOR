# Result Q, V₄′: the swapped-form light side of `Ψ₄`, `Ψ₅` on MFT `⟨TTTT⟩` with the spinning kernel, and the corner `u′ → 0`, `r → 0`, `ρ → 1`

*21 Nov 2026, V₄′ agent. Unrefereed, not committed. No existing document or code file was edited; no `bj_*` file was touched or run. New code:
`ads4/code/resultq/v4c_core.py`, `v4c_mft.py`, `v4c_wu.py`; outputs `resultq/cache/v4c_mft.txt`, `v4c_wu.txt` (git-ignored). Scratch
sympy scripts in the session scratchpad (`v4c/sym*.py`; their content is reproduced in `v4c_mft.py` part (1)). Read-only source: CLS-D
`draft.tex` (eqs. `scalarfunctional_plus_final`, `spinningfunctional_definition`, `subtraction_factor_*`). Compute: ≈ 15 CPU-minutes, ≤ 2
cores, < 1 GB.*

Grades: **[P]** proved (incl. exact computer algebra); **[D]** derived; **[D-sketch]**; **[DV]** derived and numerically verified; **[NE]**
numerical evidence; **[M]** model; **[I]** open.

---

## 0. Answer in brief

1. **A gauge in which the corner is explicit.** Complete gauge fixing (EE_FOUNDATIONS §1.4): `x₀ = 0`, `x₂ = −e`, `x₄ = +e`, `z₀` fixed;
   `x₃ = t₃z₃` (future ray), `x₁ = −t₁z₁` (patch coordinate of `x₁⁺`), `z = (1, cos φ, sin φ)`. Then exactly
   `x₁₂² = 2t₁ − 1`, `x₃₄² = 2t₃ − 1`, `u′ = (1−2t₁)(1−2t₃)/(8t₁t₃(1−c))`, `v′ = (1+2t₁)(1+2t₃)/(8t₁t₃(1−c))`, `c = cos(φ₁ − φ₃)`.
   * The light-cone edge `u′ = 0` is `t₁ = t₃ = ½` for every celestial pair; there `v′ = 2/(1−c)`.
   * **The corner `r → 0`, `ρ → 1` is the collinear limit `φ₁ → φ₃`** (`Δφ ≈ 4r/ρ`). At `t < ½` the same limit is the large-η Regge
     corner at fixed `ρ = 2t`.
   * `x₁₂²` and `x₃₄²` are independent and linear along the rays, so `∫dt₁dt₃ Δ₁₂Δ₃₄[g] = (−2πi)² Res_{t₁=½}Res_{t₃=½} g` exactly. All
     MFT edge terms are exact rational residues. The kernel is the actual spinning one: light-ray polarisations, light-transform weight
     `t³`, `j′ = 0` celestial structure with `δ₁ = δ₃ = 2`, `δ = ½ + iν`, `f⁺₄ = f⁺₂R`, and any constant (2,4) tensor (EE ≡ `(z₂·z₄)²` here).
2. **MFT, regular light families: Lemma A holds in swapped form, pointwise in the celestial angles [P].**
   * k = 4: the s-channel identity and the **J-resummed `[TT]₀` family** (`⟨T₄T₁⟩⟨T₃T₂⟩`) satisfy exactly
     **`R_s(φ₁,φ₃) = −R_t(φ₃,φ₁)`** for every (2,4) tensor. The 1↔3-even part, the only one the `j′ = 0` structure sees, vanishes identically.
     For EE, in closed form at all angles: `R_s = −R_t = −(c−3)(c−1)(c+1)/2¹⁹`.
   * Ψ₅ (k = 5, `j′ = 1`): `R_s(φ₁,φ₃) = +R_t(φ₃,φ₁)`, so the 1↔3-odd part vanishes identically. EE: `R_s = R_t = (c−3)(c−1)²(c+1)/2²⁰`.
   * Calibration: the k = 2 and k = 3 identity edge terms are exactly 0, and the `[TT]_{n≥1}` families are exactly 0.
   * Corner: each term separately goes like `|Δφ|^{9/2}·log` (k = 4) and `|Δφ|^{13/2}·log` (k = 5). That is absolutely convergent with
     margin `11/2` in powers of `r`.
   * **Residual vs cutoff:** `|I_s + I_t|/|I_t| = 1.1·10⁻⁹ / 1.3·10⁻⁹ / 1.9·10⁻⁹` (ν = 0 / 1 / 3) for every cutoff from `1` down to
     `10⁻³`. This is the precision floor of the numerical residues; the exact value is 0.
3. **MFT, the u-channel identity: a corner contact term, and V₄′ fails for the RAW functional at O(1) [D + NE].**
   * `W_u = ⟨T₄T₂⟩⟨T₃T₁⟩` is the alternating-J part of the `[TT]₀` tail, i.e. the u-channel light-cone singularity the review worried
     about. For light-ray polarisations it vanishes identically for `φ₁ ≠ φ₃`: the numerator is `≡ 0`.
   * At `φ₁ = φ₃` it is a contact distribution of homogeneity −2 in `Δφ`. The edge residue of `f₄` carries a `1/x₁₃²` pole exactly there
     (`Res f₄ ∝ v′²x₁₃²`).
   * Power count: `C ~ |Δφ|^{5/2}`, `Res f₄ ~ Δφ⁻²`, contact `~ Δφ⁻²`, measure `dΔφ`. The total is `−½`: **not absolutely convergent**.
   * Timelike-tube regularisation (`x₁ → x₁ + iεe`): the corner term grows as **`ε^{−1/4}·log ↔ r_c^{−1/2}`** (fitted `−0.31` at ν = 0,
     `−0.37` at ν = 1). At `ε = 10⁻⁵` it is already 8× each regular term.
   * Cone-preserving prescription (`t → t + iε`): it is exactly 0, but that prescription is not an interior tube direction at the corner.
   * The same `ε^{−1/4}` appears in the **unswapped commutator form** at interior large-η corner points (`ρ = 0.6`), for raw k = 4 and raw
     k = 2 alike. Controls without the `f` pole (`f = 1`; k = 5) give `ε^{+3/4} → 0`.
4. **Q's de-shifted `Ĉ₄^{ds}` cures it [D-sketch].**
   * The 4-shift weight moves the ν-contour by 3 units on `Q^{iν/2}`, with `Q = A₁₃/(A₁₀A₃₀) ~ Δφ²`. This is the same shift that turns the
     raw `η^{5/2}` into `η^{−½}log η`.
   * The contact term then goes like `ε^{+5/4}` (fitted `1.19`), so it tends to 0 in any regularisation.
   * **So the MFT light side of `Ĉ₄^{ds}` is exactly 0, provided the ν-integral is taken inside the corner integral (kernel level).**
5. **O(G) [D-sketch, power counting only].**
   * A twist-1 (graviton/T) u-channel exchange is suppressed by `(x₁₃²)^{1/2} ~ r` relative to the identity.
   * The raw corner then converges, `ε^{+1/4}` (fitted `0.19`, margin `r^{1/2}`), and the de-shifted one has margin `r^{7/2}` (fitted
     `ε^{1.69}`). No corner remnant at O(G).
   * An exact spinning O(G) swapped-form evaluation was **not** done [I].
6. **Verdict on V₄′.**
   * **As stated (the raw CLS-D integral `Ψ_{4,ν}` converges absolutely at the corner): false, already at O(1)**, with an explicit MFT
     counterexample. The size is a regularisation-dependent, divergent `O(1)·C_T²` term, i.e. `O(1/G)` relative to the O(G) sum-rule
     content: infinite leverage for any per-ν use.
   * **Restated as V₄″ for Q's actual functional** (de-shifted total kernel, ν-integral inside): it holds at O(1) [P for the regular
     families + D-sketch for the contact term], and at O(G) by power counting [D-sketch].
   * The upper bound survives with V₄″ in place of V₄′. **Two clarifications enter Q's hypothesis list:**
     * the de-shift (hypothesis 5) is load-bearing on the light side too;
     * the per-ν Lemma A of Stage 12 §2.1 is not available at the large-η corner for MFT.

---

## 1. Setup (`v4c_core.py`)

**Configuration.** CLS-D's final form (`scalarfunctional_plus_final`, spinning version):
* causal order `4 > 3 > 0 > 1 > 2`;
* `T₁`, `T₃` polarised along their light rays;
* light-transform weight `t^{Δ+J−2} = t³`;
* celestial `j′ = 0` structure `C = A₁₃^{(2+δ)/2}(A₁₀A₃₀)^{−δ/2}`, with `A_{ij} = −2z_i·z_j`, `δ = ½ + iν`, `δ₁ = δ₃ = Δ_T − 1 = 2`.

**Fixed gauge.** With `x₀, x₂, x₄, z₀` fixed (10 conditions; finite residual stabiliser), the dual structure is a constant (2,4) tensor `S`.
* For EE, `e₂₄ = (I(x₀₂)z₂)·(I(x₀₄)z₄) = z₂·z₄`, because `I(±e)` is the time reflection. So EE is the full contraction `T₄^{μν}T₂_{μν}`.
* Validity of the sum rule cannot depend on `S` (EE_FOUNDATIONS §1.4), so all 81 components are tested.

**The `⁺` images.** They only produce signs `(−1)^{Δ+J} = −1` per 2-point factor linking a `⁺` point to a non-`⁺` point. Every Wick term
has an even number of such factors, so the naive patch formulas enter with sign +1 [D].

**Double commutator.** `Δ₁₂Δ₃₄` flips the `(12)`, `(34)` iε of the 2-point functions and of `f` alike (Lemma C). Since `x₁₂²` and `x₃₄²` are
linear in `t`, the `t`-integrals are exact residues at `t = ½`.

**MFT decomposition.** With `T(x)T(y,z) = v_μv_ν/(X²)³`, `v = I(X)z`:

| term | content | double-commutator support |
|---|---|---|
| `W_s` | s-channel identity | 6th-order poles (k = 4); 5th (k = 2, 3) |
| `W_t` | `[TT]₀`, all J | `f`'s simple poles |
| `W_u` | u-channel identity | 0 off `Δφ = 0` (§3) |
| `[TT]_{n≥1}` | — | double zeros, exactly 0 |

No J-sum is performed: `W_t` *is* the resummed `[TT]₀` family on the edge.

---

## 2. MFT: regular light families (`v4c_mft.py` → `cache/v4c_mft.txt`)

### 2.1 Exact results [P]

* **Calibration (sympy, rational celestial points `z = (1, (1−m²)/(1+m²), 2m/(1+m²))`).**
  * k = 2 (`f⁺₂ = v′ − u′`) and k = 3 (`f⁻₃ = −1`): `Res_s = 0` exactly, for EE and for random (2,4) tensors. The identity has no k = 2, 3
    light side, as it must.
* **k = 4:**
  * `Res_s(P) = −Res_t(Q)` and `Res_s(Q) = −Res_t(P)` exactly (`P = (1/3, 2/5)`, `Q` = swap), for EE and random tensors.
  * Symbolic in the angle (EE): **`R_s = −(c−3)(c−1)(c+1)/2¹⁹`, `R_t = +(c−3)(c−1)(c+1)/2¹⁹`, `R_s + R_t ≡ 0`.**
  * These are the edge residues; `Ψ⁺` carries a further `(−2πi)²`.
* **k = 5:**
  * `Res_s(P) = +Res_t(Q)` exactly.
  * EE: `R_s = R_t = (c−3)(c−1)²(c+1)/2²⁰`, even under `φ₁ ↔ φ₃`. The `j′ = 1` celestial structure is 1↔3-odd, so it annihilates it.
* **Reading.**
  * The s-channel identity edge term at `(φ₁,φ₃)` is the 1↔3 image of the `[TT]₀` edge term, with sign `∓` for `f^±`.
  * This is the position-space, spinning counterpart of RV §2.3's Mellin statement ("double-twist residues at `S′ = 2Δ_φ, 2Δ_φ+2` cancel
    those at `S′ = T, T−2`). Here it holds **pointwise in the celestial angles**, and for every dual structure.

### 2.2 Numerics [DV]

The Cauchy–FFT residues (radius 0.2, 32×32 nodes) on 300 random celestial pairs, over all 81 tensor components:
* `max |R_s(a,b) + R_t(b,a)|/scale = 1.2·10⁻⁵`, median `7·10⁻¹⁰`;
* without the swap, the median is `0.64`;
* k = 2: `|R_s| ≤ 3·10⁻⁸` relative.

### 2.3 Integrated EE light side vs corner cutoff (`|Δφ| > d_c`; `K_ν` = the `φ̄`-integral of `(A₁₀A₃₀)^{−δ/2}`, mpmath)

| ν | `I_t` (`d_c → 0`) | `|I_s + I_t|/|I_t|` at `d_c = 1, 0.3, 0.1, 0.03, 0.01, 10⁻³` | corner exponent of each integrand |
|---|---|---|---|
| 0 | `−1.516347·10⁻²` | `1.1·10⁻⁹` (all) | `4.41` (analytic `9/2`, × log) |
| 1 | `−(4.3518 + 3.9975i)·10⁻³` | `1.3·10⁻⁹` (all) | `4.40` |
| 3 | `+(4.1340 + 5.5170i)·10⁻⁴` | `1.9·10⁻⁹` (all) | `4.52` |

* Each term converges by `d_c ≈ 0.1`; the tail goes like `d_c^{11/2}`.
* The residual is flat at the residue precision. **Lemma A in swapped form, identity + `[TT]₀`, holds with an absolutely convergent corner.**
* The s-term exponent at ν = 1 (`3.73`) is precision-limited (`R_s` is a 6th-order residue).

---

## 3. MFT: the u-channel identity at the corner (`v4c_wu.py` → `cache/v4c_wu.txt`)

**Pointwise.** `N = (z₁·z₃)x₁₃² − 2(z₁·x₁₃)(z₃·x₁₃) ≡ 0` on the cones (unregularised: `10⁻¹⁸ … 10⁻²²`), so `W_u = 0` for `φ₁ ≠ φ₃`.

**Why a prescription is needed.** At `φ₁ = φ₃` both `N` and `x₁₃²` vanish. The Wightman wave-front covector `(k, −k)`, `k ∝ z`, is conormal
to the product of the two cones, so the restriction of `⟨T₃T₁⟩` there is not defined without a prescription.

**Edge term.** `U = (−2πi)²Res_t[t³t³f₄T₄₂N²/(X₁₃²)⁵]`, with `T₄₂ = −1/64` for EE.

**Pairing.** `G(ε) = ∫dΔφ A₁₃^{(2+δ)/2}K_ν(Δφ)U_ε(Δφ)`:

| case | fitted `|G| ~ ε^p` | power count | reading |
|---|---|---|---|
| k = 4 edge, timelike tube `x₁ → x₁ + iεe`, ν = 0 | **`−0.307`** | `−¼` (× log) | **divergent**: `|G| = 0.0143 → 0.1246` (`ε = 10⁻² → 10⁻⁵`) vs `|I_t| = 0.0152` |
| same, ν = 1 | **`−0.371`** | `−¼` | divergent, oscillating |
| cone-preserving `t → t + iε` | `G ≡ 0` (exact algebra) | — | `Im x₁₃ = ε(z₁+z₃)` is null at the corner: not an interior tube direction |
| control: no `f` pole (restriction only) | `+0.69` (ν = 0), `+0.63` (ν = 1) | `+¾` | → 0 |
| k = 5 edge (`Res f⁻₅ ∝ v′x₁₃² = O(1)`: no `x₁₃²` pole) | `+0.69` | `+¾` | → 0; besides, `U` is even and the `j′ = 1` structure odd: pairing ≡ 0 |
| **commutator form**, interior corner point `t = 0.3` (`ρ = 0.6`), raw `f₄` | **`−0.297`** | `−¼` | divergent (per-ν Lemma A ill-defined there) |
| commutator form, raw `f₂ = v′ − u′` (k = 2) | **`−0.297`** | `−¼` | divergent: the known raw-`Ĉ₂` large-η failure (Stage 3 §1.5, Result D) |
| commutator form, `f = 1` | `+0.69` | `+¾` | → 0 |
| k = 4 edge, contour-shifted kernel (`δ → δ + 3`, the de-shift) | **`+1.19`** | `+5/4` | → 0 |

**Power count** (`ε ~ Δφ² ~ r²` at `ρ = 1`). Collect the factors:
* celestial `|Δφ|^{5/2}`;
* `f₄` edge residue `Δφ⁻²` (the `x₁₃²` pole);
* the contact distribution, degree `−2`;
* the measure, `+1`.

The total is **`−½`**: a deficit of half a power of `r`. Without `f`'s pole it is `+3/2`, and with the de-shift (`+3`) it is `+5/2`.

**What it means.**
* The review's mechanism is real, but it sits at **O(1), not O(G)**. The `[TT]₀` tail with alternating sign builds the u-channel
  light-cone singularity on top of `f₄`'s `x₁₃²` pole at the corner.
* With light-ray polarisations its pointwise value is 0, but its contact term is not integrable against the raw kernel.
* **The raw per-ν `Ψ_{4,ν}` (and, in commutator form, raw `Ψ_{2,ν}`) is not absolutely convergent at the corner on MFT.** Its value depends
  on the prescription (0 or divergent).
* This is the light-side face of Stage 3 §1.5's caveat (raw `Ĉ₄`, raw `Ĉ₂` fail the large-η corner). Here it is realised by a physical,
  unitary O(1) correlator, not only under the uniform Regge bound.

**The de-shift [D-sketch].**
* Q's `Ĉ₄^{ds} = ∫dν χ(ν)Ψ_{4,ν}` (Stage 4) depends on ν near the corner only through `Q^{±iν/2}` (two branches, `η^{±iν}`) times
  `ν`-only factors (the gauge-fixed dual constants are `const^{iν}` and phases `e^{πν·b}`, entire).
* The weight that turns raw `η^{5/2}` into `η^{−½}log η` shifts each branch by 3 units, a gain of `Δφ³` up to a log.
* So the contact term goes like `ε^{5/4} → 0` in every interior tube regularisation, and the MFT light side of `Ĉ₄^{ds}` is exactly 0.
* **Inputs:**
  * hypothesis 5 (continuum de-shift, analytic in the strip `|Im ν| ≤ 3` up to strip-edge poles);
  * the ν-integral is taken at the kernel level, inside the corner integral;
  * the test here shifts only the `A₁₃` power, which is the factor that sets the corner exponent.

---

## 4. O(G) [D-sketch; exact spinning O(G) evaluation not done, I]

**What an O(G) remnant needs.** At O(G) the swapped-form light families are the s-channel T block and `[TT]₀` with `(γ_J, δp_J)`. The
`[TT]_{n≥1}` families enter only at O(G²), since `sin²` is then `O(γ²)`. An O(G) remnant needs a corner term that fails to vanish.

**s-channel T block.** It is regular at the corner: its collinear block grows only like `log` as `v′ → ∞`. Per block, it is absolutely
convergent with margin `≥ 3/2`.

**The u-channel singularity of the `[TT]₀` tail.** At O(G) it is the u-channel T exchange.
* In the (13) light-cone limit, a twist-τ exchange carries `(x₁₃²)^{τ/2}` relative to the identity. With `τ_T = 1`, the deficit `−½`
  becomes `+½`.
* Natural-size model (identity numerator structure, `(X₁₃²)^{5−τ/2}`): raw `ε^{0.19}` (power count `+¼`); de-shifted `ε^{1.69}`
  (`+7/4`).
* **No corner remnant at O(G)** in either form, with margins `r^{1/2}` (raw) and `r^{7/2}` (de-shifted).
* The light-ray-polarised `⟨T₃T₁T⟩` structures are not identically zero, unlike the identity's. Only their size is bounded here.

**Cancellation structure.**
* The exact O(1) identity `R_s(a,b) = ∓R_t(b,a)` is the spinning, position-space form of the scalar-Mellin cancellation (RV §2.3).
* At O(G) the same pairing (s-channel light families ↔ 1↔3-crossed t-channel ones) is what Lemma A requires. Its spinning O(G) check
  (Einstein `⟨TTTT⟩` in swapped form) remains **[I]**.

---

## 5. Verdict and regradings

| claim | grade |
|---|---|
| Complete-gauge reduction: corner = collinear limit `φ₁ → φ₃` at `t = ½`; exact residue form of all MFT edge terms | [D] |
| MFT, k = 4: s-identity + J-resummed `[TT]₀` cancel pointwise, `R_s(a,b) = −R_t(b,a)` for every (2,4) structure; EE closed form | **[P]** (exact algebra) + [DV] (81 components, 300 pairs) |
| MFT, k = 5 (`j′ = 1`): `R_s(a,b) = +R_t(b,a)`, odd part ≡ 0 | **[P]** |
| k = 2, 3 identity edge term = 0 | [P] |
| corner power counting of the regular terms: `|Δφ|^{9/2}log` (k = 4), margin `11/2` | [P] + [NE] |
| u-channel identity contact term: raw `Ψ_{4,ν}` not absolutely convergent at the corner (deficit `r^{1/2}`); timelike-tube value diverges as `ε^{−1/4}log` | **[D + NE]** |
| same for the commutator form of raw `Ψ_{4,ν}`, `Ψ_{2,ν}` along the large-η corner | [D + NE] |
| Ψ₅ raw: no contact term (no `x₁₃²` pole; parity) | [D + NE] |
| de-shifted `Ĉ₄^{ds}`: contact term → 0; MFT light side exactly 0 | **[D-sketch]** (hypothesis 5 + kernel-level ν-integral) |
| O(G) corner: no remnant (twist-1 suppression), margins `r^{1/2}` raw, `r^{7/2}` de-shifted | [D-sketch] (power counting; natural-size model [M]) |
| exact spinning O(G) (Einstein) swapped-form cancellation | [I] |

**V₄′ → V₄″.**
* **Statement.** On tree-level `⟨TTTT⟩` (`J₀ < 4`), the **total de-shifted** kernels `Ĉ₄^{ds}`, `∂Ĉ₄^{ds}` (and the k = 5 analogue), with the
  ν-integral inside, pair absolutely convergently with `𝒢_rest` at the corner `u′ → 0`, `r → 0`, `ρ → 1`.
* **The swapped form equals the commutator form there.**
  * O(1): identity + `[TT]₀` cancel exactly [P]; the u-channel contact term vanishes [D-sketch].
  * O(G): power counting [D-sketch].
* **Grade:** [D-sketch], with its O(1) regular part now [P].
* **The literal V₄′ for the raw per-ν `Ψ_{4,ν}` is false** (§3, MFT counterexample). Stage 12 §2.1's per-ν Lemma A must be read with the
  de-shifted kernel at the large-η corner.

**Effect on the bound.**
* **No O(G) remnant is found.** The leverage `10¹⁰–10²²` is not triggered: the O(1) failure lives only in the raw functional, which Q
  does not use.
* Hypothesis 5 (de-shift existence) is now also load-bearing for the **light** side. That adds weight to the open spinning `j′ = 1` 3-shift
  de-shift of `Ĉ₅` (Stage 4 [I]).
  * For `Ĉ₅` the light side is safe even raw at O(1) (§3).
  * But its heavy-side validity still needs the de-shift.

**Next (cheapest).**
1. Repeat §3 with the stored `χ(ν)` of `functionals_q4_deshift.json`, with the gauge-fixed ν-normalisation matched to `N̂₄(ν)`. This
   would make the de-shift cure numerical rather than a contour argument.
2. The O(G) swapped-form check with the u-channel T exchange in its collinear (light-ray OPE) form, in the same gauge.

---

## 6. Files (`ads4/code/resultq/`)

| file | content | output | time |
|---|---|---|---|
| `v4c_core.py` | complete-gauge geometry, spinning 2-pt tensors, `f_k`, Wick terms, Cauchy–FFT edge residues (EE and full tensor), celestial structure | — | — |
| `v4c_mft.py` | (1) sympy exact residues k = 2–5, EE closed forms; (2) 300-pair × 81-component numerics; (3) EE light side vs corner cutoff, ν = 0, 1, 3; (4) corner exponents | `cache/v4c_mft.txt` | 35 s |
| `v4c_wu.py` | u-channel identity contact term: timelike-tube vs cone-preserving; k = 4, 5; control; commutator form (`f₄`, `f₂`, `1`); de-shift contour shift; O(G) twist-1 model | `cache/v4c_wu.txt` | ≈ 10 min |

Run order: `python3 v4c_mft.py` (imports `Kfun` from `v4c_wu`), then `python3 v4c_wu.py`. Thread variables ≤ 2.
