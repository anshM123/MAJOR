# LARGESPIN_ANALYTIC: analytic large-spin / large-impact-parameter control of Ω^EE at Δ_g = 1000 — what can be proved, what cannot, and why

*23 Nov 2026. Unreviewed, not committed. No existing file was edited. HIGHSPIN.md and `code/highspin/` were read, not touched.*
* *New code:* `ads4/code/largespin/` — `ls_core.py` (NF02 weights, `Φ_k`, the leading-Regge action `Ω_LR,k(t)` by Mehler–Dirichlet, a
  mpmath-free hybrid `hyb_fast`), `ls_curvbound.py` (Lemma C′), `ls_barrier.py`, `ls_struct.py`, `ls_xterm.py`, `ls_curvfloor.py`
  (impact-parameter lemma, tolerance profile, structures), `ls_agscan.py` (AG saddle test on H_R′ states, thin wrapper of
  `spinfunc/lj_scan.py`), `ls_strip.py` (H_bulk), checks `t_core.py`, `t_hyb.py`, `t_env.py`, `t_xsigned.py`, driver `run_all.sh`.
* *Caches:* `code/largespin/cache/` (its `.gitignore` ignores everything except itself).
* *Read-only imports/data:* `spinfunc/ee1000_common.py` (→ `w2_chain`, `q3_analysis`, `nf_common`), `spinfunc/lj_scan.py`, `lj_saddle.py`,
  `poserr/cache/large_m_1000*.json`, `prl/data/d_abel_X_1000.txt`, `uniform/cache/residual_coeffs_Df3.pkl`.
* *Compute:* ≈ 1.3 CPU-h, ≤ 2 processes, `ulimit -v 3000000`, `OMP_NUM_THREADS=1`.

Grades as UNIFIED_THEOREM_v3 §0. `c₂ = m² ×` relative correction; `b̂ = (2J+1)Λ/m`; `β` the H₂ impact parameter, `sinh(β/2) = (J+½)/m`.
"Floor" = the certified flat floor `N_flat/5` (CERT_CLOSE) plus `min(0, curvature)`.

---

## 0. Summary

### 0.1 Table

| hypothesis | old grade | new grade | what is proved / derived here (explicit constants) | residual (model) |
|---|---|---|---|---|
| **H_LJ(1000)** (near gap, even J ≥ 50, odd J ≥ 33, m ≤ 1.3Λ) | [D-sketch + M] | **[D-sketch + M] — unchanged grade; restated with a derived criterion** | (i) **Lemma I [P]**: action of any image `I = ∫P_{−½+iν}(cosh t)σ(dt)` is `∫Ω_LR(t)σ(dt)` exactly; `Ω_LR,2 > 0` and decreasing on `t ∈ [0.002, 4]` [NE]. (ii) Tolerance profile (ψ₂-weighted, all 42 H_LJ sample states): weight tolerated **at the ring t = β: ≥ 0.062** of the ring weight; at β/3: ≥ 5.8·10⁻³; at β/10: ≥ 1.6·10⁻⁴; at t ≤ 1/Λ: ≥ 6.9·10⁻⁸. (iii) Structures, tolerated coefficient (ψ₂ weight): ring amplitude `(a/m)P`: **a ≥ 80**; `(a/m²)ν²P`: a ≥ 37; ring displacement `(a/m)β∂_βP`: a ≥ 67; index-shifted `(a/m)Re P_{½+iν}`: a ≥ 10.8; resolvent-smeared `(a/m)P/(4ν²+1)`: a ≥ 12.8; **window bump `(a/m²)e^{−(νβ/3)²}`: a ≥ 379 c₂** (192 with the sign-free de-shift term); **full-range smooth `a/m²`: a ≥ 0.11 c₂**. (iv) AG to all orders excludes the full-range smooth class on `ν ≥ ν_B` (HLJ, [D]). | the **Regge-window structure**: the non-ring part of the window correction must be ≲ 190–380 c₂ (J = 344) … 630–815 c₂ (J = 100); measured ≲ 550 at J = 100 (L ≈ 1.1–1.5), ≈ 10 c₂ for the exact scalar; unmeasured J ≥ 150. Plus AG beyond all orders and **the de-shift caveat (N2)** |
| **H_R′** (m ≥ 1.3Λ, b̂ ≥ 5) | [M] (b̂ < b̂_A: [NE + M]) | **b̂ ≥ b̂_A(m): [M] → [D-sketch + M]** (same structure and weakest link as H_LJ); 5 ≤ b̂ < b̂_A: [NE + M] unchanged | (i) **Lemma C′ [P + DV]**: hybrid floor `λ_hyb ≥ (N_flat/5 − 2B₂(β) − uB₃(β))/N_s > 0` on every sampled state with **b̂ ≤ 119**, all m ∈ [1.3, 20]Λ (rig′ ≥ 0.049 at b̂ ≈ 119; ≥ 0.135 at b̂ ≤ 91). With ROUTEA's certified window `[91.2, 2Λ]` (m ≥ 1.5Λ) the **leading-Regge hybrid is now rigorous on all b̂ ≤ 2Λ for m ≥ 1.5Λ**. (ii) AG checked on 5 new H_R′ states (m = 1.5–5Λ, b̂ = 50–300): non-degenerate travelling-wave saddles on ν ∈ [3, Λ], `min Re κ̃/β = 1.0007–1.0024`, vertical decay ≥ **1.0030β** [D + NE]. (iii) Tolerances over 89 H_R′ sample states: window bump ≥ 324 c₂ (150 with the de-shift term), ring amplitude ≥ 6.5 (flat-only floor) — ×3–460 larger with the computed curvature in the floor at b̂ ≥ 300 (§3.5). | as H_LJ; in addition `m ∈ [1.3, 1.5)Λ, b̂ ∈ [119, b̂_max]`: hybrid floor only [NE]; `b̂ > 2Λ`: hybrid [NE] |
| **H_bulk** (ν ≥ 20, J > 0, envelope `3(450+100J)` c₂) | [M], ×76 on S_D | **[M] — unchanged; restated in relative form, strip family** | On S_D the bulk correction may reach **ε_bulk ≥ 3.9 % of the flat envelope E_k(x)** (9.8 % at 1.3Λ) before any state fails; strips with a stated bulk factor f: f ≥ 10 for b̂ < 6.9 (1.3Λ) … 13.8 (20Λ); f ≥ 3 for b̂ < 8.5 … 15.3. Beyond b̂ ≈ 15 the relative tolerance drops to 10⁻³ … 3·10⁻⁵ (no sign-free argument). | the bulk-point (fixed x = ν/m) 1/m expansion of the **spinning** heavy image is not available; no analytic route (§5.3) |

### 0.2 Bottom line

* **No hypothesis is removed, and none is replaced by a proof.** The decisive missing input is the same for H_LJ(1000) and H_R′: the
  structure of the beyond-leading-Regge correction **inside the Regge window** `ν ≲ 3/β` (equivalently, its density in transverse-H₂
  impact-parameter space at radii `t < β`). Deriving it needs the spinning heavy-block Regge formula at subleading order (the ⟨TTTT⟩
  analogue of CHMRS (3.17)–(3.20)), which REGGE_FROM_BLOCKS §8 lists as not derived. Nothing in the project's toolkit (the 1/Δ engine,
  the saddle analysis, the exact scalar skeleton) reaches it at `J²/Δ ≳ 2` (§6).
* **What is new and rigorous:**
  1. **Lemma I** (impact-parameter representation of the per-state action) and the **tolerance profile** `tol(t)`: the correct norm for
     corrections at large b̂. It explains exactly why sign-free per-ν bounds fail (a correction at radius `t ≲ 1/Λ` is penalised ×10⁵–10⁶
     relative to one on the ring) and why the leading-Regge-shaped (`D_h`) excess at (1100, 100) is not a violation *in the ψ₂ convention*
     (ring amplitude a ≈ 2.7–5.9 vs a_tol ≈ 200).
  2. **Lemma C′**: the hybrid floor is rigorous up to b̂ ≈ 119 (Lemma C stopped at 10), closing the gap to ROUTEA's certified window.
  3. The AG saddle property is checked on H_R′ states (it was only checked near the gap), so H_R′ at large b̂ has exactly the
     decomposition of H_LJ.
* **What is negative (reported prominently):**
  * **N1.** The analytic replacement asked for is **not achieved**; §6 says precisely why each route fails.
  * **N2 (new, possibly affects HLJ_ANALYTIC / HWIN_THIMBLE).** Ω^EE's k = 2 part acts on *raw* images with the min-norm Abel weight `X`,
    not ψ₂; `|X − ψ₂|/ψ₂ = 0.3–0.9` at ν ≤ 10 and the Regge window at J ≥ 150 lies entirely in ν < 25 where `X ≠ ψ₂`. The ψ₂-convention
    tolerances of ring-type corrections hold **only if the correction is a ν-independent relative factor of the raw leading-Regge image on
    ν ≲ 25** (then the de-shift identity makes it exact). Counted sign-free, the de-shift mismatch collapses the ring tolerance: a_tol[A]
    80 → **0.31** (J = 344), 200 → **3.2** at (1100, 100), **below the measured D_h-shaped excess a ≈ 2.7–5.9**. HLJ's window tolerance
    `T_B` used ψ₂ inside the window; with |X| the window-bump tolerance halves (379 → 192 c₂ at J = 344; 815 → 628 at (1100, 100), vs the
    measured remainder ≲ 550: L ≈ 1.1). Which convention the HWIN/fw `Γ` are in must be checked before the refined H_win is quoted (§3.6).
  * **N3.** The index-shifted Regge functions `P_{l+n}(cosh β)` on which HLJ §5.2's [D-sketch] rests have H₂ densities **inside** the disc
    `t ≤ β`, not on the ring (Paley–Wiener); they are tolerated only at a ≤ 10.8 (J = 344, ψ₂). The [D-sketch] therefore does not by itself
    make window corrections harmless (§6(c)).
  * **N4 (doc correction).** ROUTEA's certified hybrid window `[91.2, 2Λ]` holds for **m ≥ 1.5Λ** (its `|w₃| ≤ 4/9`); UT3 §5.2 and PE §2.4
    quote it on `m ≥ 1.3Λ`. Lemma C′ covers `b̂ ≤ 119` at 1.3Λ; `m ∈ [1.3, 1.5)Λ, b̂ ∈ [119, 530]` is hybrid-[NE] only.

---

## 1. Setting

For a heavy state `(J, m)`, helicity `h`, the diagonal action is `diag_h = F_h(flat) + curv_h + δ_h` (PE §2.1): certified flat NF02 part
(`F_h ≥ N_flat/5`, CERT_CLOSE), exact leading-Regge curvature `curv = 2∫ψ₂[P_{−½+iΛp}(cosh β) − J₀(b̂p)] + uσ_h∫ψ₃[P − J₀]`, and the
beyond-leading-Regge correction δ, the only model content. The hybrid of `spinfunc/ee1000_common.hyb` is reproduced by `ls_core.hyb_fast`
to `10⁻¹¹` relative (`t_hyb.py`; curvature by Mehler–Dirichlet instead of mpmath `₂F₁`, which fails at large β).

**Leading-Regge action at impact parameter t.** `Ω_LR,k(t) := ∫₀¹ψ_k(p)P_{−½+iΛp}(cosh t)dp = (2/π)∫₀^{π/2}Φ_k(Λs)/cosh(s/2) dφ`,
`s = 2 asinh(sinh(t/2) sin φ)`, `Φ_k(x) = ∫₀¹ψ_k cos(px)dp` (Mehler–Dirichlet in the Lemma-C variables). Checked against direct
`₂F₁` quadrature to 2·10⁻⁷ (`t_core.py`). Table (`cache/barrier_1000.txt`):

| t | 0.003 | 0.01 | 0.05 | 0.1 | 0.2 | 0.3 | 0.5 | 1.0 | 2.0 | 4.0 |
|---|---|---|---|---|---|---|---|---|---|---|
| b̂ = 2Λ sinh(t/2) | 3.0 | 10 | 50 | 101 | 202 | 298 | 501 | 1035 | 2340 | 7243 |
| `Ω_LR,2(t)` | 7.46 | 0.441 | 4.07e-3 | 5.11e-4 | 7.01e-5 | 2.71e-5 | 1.18e-5 | 7.58e-6 | 5.39e-6 | 2.14e-6 |
| `b̂³Ω_LR,2` | 199 | 442 | 514 | 520 | 576 | 719 | 1478 | 8415 | 6.9e4 | 8.1e5 |

`Ω_LR,2 > 0` and non-increasing on the whole grid (300 points, `t ∈ [0.002, 4]`) [NE]; the flat `b̂⁻³` regime (`b̂³Ω ≈ 512 = −A₂/4`) crosses over
to the AdS `A₁M₁(β)/Λ²` regime at b̂ ≈ 300–500 (ROUTEA §4.1). `Ω_LR,3` changes sign near t ≈ 0.01–0.02.

---

## 2. Lemma C′: the hybrid floor up to b̂ ≈ 119 [P + DV]

**Statement.** With `σ(φ) = 2 sinh(β/2) sin φ` (so `Λσ = b̂ sin φ`), `s = 2 asinh(σ/2)`, `w₂ = 2`, `w₃ = uσ_h`:

    curv_k = w_k (2/π) ∫₀^{π/2} [Φ_k(Λs)/cosh(s/2) − Φ_k(Λσ)] dφ ,
    |…| ≤ E₀(Λs)·s²/8 + Λσ³/24 · E₁(Λs),       E_j(x) ≥ sup_{y ≥ x} |Φ_k^{(j)}(y)|,

so `|curv_k| ≤ |w_k| B_k(β)`, `B_k` the cell-wise upper sum (E_j non-increasing, `s, σ` increasing in φ), and

    λ_hyb ≥ rig′ := (N_flat/5 − 2B₂(β) − uB₃(β))/N_s .

**Proof.** Mehler–Dirichlet [Q] in the substitution of Lemma C (PE §2.1); `1 − sech(s/2) ≤ s²/8`; `0 ≤ σ − s ≤ σ³/24` (`asinh y ≥ y − y³/6`);
mean value theorem. Decay of the envelopes, by integration by parts with `ψ_k(0) = 0` and the factor `(1−p)⁵` (so all boundary terms at
p = 1 vanish): `|Φ(x)| ≤ (|ψ′(0)| + ∫|ψ″|)/x²`; `|Φ′(x)| ≤ ∫|(pψ)″|/x²` and `≤ (|(pψ)″(0)| + ∫|(pψ)‴|)/x³`; `|Φ″(x)| ≤ ∫|(p²ψ)″|/x²`. ∎
Constants (k = 2 / 3, padded 10⁻⁶): `∫|ψ| = 10.386 / 12.836`, `∫p|ψ| = 3.533 / 5.108`, `C₂ = 489.6 / 451.6`, `C₃ = 113.9 / 208.1`,
`C₃′ = 1481.7 / 2073.1`. `E₀, E₁` are the running suprema of `|Φ|, |Φ′|` on a grid of step 0.01 up to x = 6000 (float evaluation of `Φ`,
Lipschitz pads from the analytic decay bounds), and the analytic bounds beyond. Grade **[P]** (inequalities) **+ [DV]** (grid values of `Φ`
in floating point; the purely analytic-constant version is also printed and gives rig′ > 0 to b̂ ≈ 80).

**Result** (`cache/curvbound_1000.txt`, worse helicity, J = ⌊(b̂m/Λ − 1)/2⌋):

| m/Λ \ b̂ | 10 | 30 | 50 | 70 | 91 | 100 | 119 | 150 |
|---|---|---|---|---|---|---|---|---|
| 1.3 | 0.1999 | 0.1977 | 0.1891 | 0.1698 | 0.1350 | 0.1127 | **0.0490** | −0.098 |
| 2 | 0.2000 | 0.1980 | 0.1901 | 0.1722 | 0.1379 | 0.1175 | **0.0570** | −0.078 |
| 20 | 0.1999 | 0.1980 | 0.1904 | 0.1730 | 0.1396 | 0.1202 | **0.0619** | −0.067 |

**rig′ > 0 on every sampled state with b̂ ≤ 119 at every m/Λ ∈ {1.3, 1.5, 2, 3, 5, 10, 20}** (sampling in b̂ only; rig′ is monotone
decreasing in b̂ along every sampled mass, so the statement between samples is [NE]). The bound `B ≈ 8·10⁻⁵` is ≈ 5× the true curvature
(1.5·10⁻⁵): the crossover to the curvature-dominated regime is what stops it, and ROUTEA's positivity certificate (`curv > 0`, b̂ ∈
[91.2, 2Λ], m ≥ 1.5Λ) takes over with overlap [91.2, 119].

**Consequence.** The leading-Regge hybrid (flat + curvature) is now rigorous (**[certified] flat + [P + DV] curvature**) for
b̂ ≤ 119 at all m ≥ 1.3Λ, and for all b̂ ≤ 2Λ when m ≥ 1.5Λ. This removes the curvature from "the AdS correction" everywhere except
`m ∈ [1.3, 1.5)Λ, b̂ > 119` and `b̂ > 2Λ`. It does nothing for δ.

---

## 3. The impact-parameter representation and the right norm for δ

### 3.1 Lemma I [P]
Let σ be a finite signed Borel measure on `[0, ∞)` and `I(ν) := ∫P_{−½+iν}(cosh t)σ(dt)` for `ν ∈ [0, Λ]`. Then

    ∫₀¹ψ_k(p) I(Λp) dp = ∫ Ω_LR,k(t) σ(dt).

*Proof.* `|P_{−½+iν}(cosh t)| ≤ P_{−½}(cosh t) ≤ 1` (Mehler–Dirichlet), `ψ_k ∈ L¹`; Fubini. ∎

σ is the image's density in transverse-H₂ impact-parameter space (the spherical transform). It need only reproduce I on `[0, Λ]` (ψ vanishes
beyond), so it is not unique; any admissible choice gives the same action. Examples: the leading-Regge image `P(cosh β)` is **exactly the
ring** `σ = δ_β` [P]; a constant-in-ν correction is a density at `t ≲ 1/Λ`; a window-smooth one (smooth on scale ν ~ 1/β) a density spread
over `t ≲ β/3`.

**Corollary (tolerance profile).** A correction to the k = 2 image with density σ_δ changes the action by at most `2∫|Ω_LR,2||σ_δ|`.
The weight tolerated at radius t, in units of the ring weight 2, is `tol(t) = floor/(2|Ω_LR,2(t)|)`; since `Ω_LR,2` is decreasing, a density
carried by `t ≥ (1−ε)β` of total weight `w` is harmless if `w < τ(ε) := floor/(2Ω_LR,2((1−ε)β))`.

### 3.2 Tolerance profile over H_LJ(1000) and H_R′ (`ls_struct.py`, `cache/struct_1000.txt`)
42 H_LJ states (J = 26–344 incl. odd, `m = m_gap, 1.1m_gap, 1.2m_gap ≤ 1.3Λ`) and 89 H_R′ states (`m/Λ ∈ {1.3, …, 20}`, b̂ = 5 … 2000,
J ≤ 3000). Flat-only floor (`N_flat/5`):

| | ring `t = β` | `0.8β` | `β/3` | `β/10` | `t = 1/Λ` | `τ(0)` | `τ(0.2)` | `τ(0.5)` |
|---|---|---|---|---|---|---|---|---|
| H_LJ, min | **0.062** | 0.047 | 5.8e-3 | 1.6e-4 | 6.9e-8 | 0.062 | 0.047 | 0.019 |
| H_LJ, max | 0.208 | 0.107 | 1.0e-2 | 2.7e-4 | 1.3e-4 | | | |
| H_R′, min | 2.2e-3 (b̂ = 2000) | 1.9e-3 | 1.3e-3 | 1.3e-4 | 1.3e-9 | 2.2e-3 | 1.9e-3 | 1.6e-3 |

(The H_R′ minima are at b̂ ≥ 800, where the flat-only floor is a poor measure; §3.5.) `tol(β/3) ≈ 0.2/27` for every state is the flat
`t⁻³` law. **This is the quantitative content of "smooth corrections are amplified by b̂³…b̂⁵" (SUBLEADING_REGGE §3, PE §2.2):** the same
weight costs ×10⁵–10⁶ more at the origin than on the ring.

### 3.3 Candidate correction structures (k = 2 image; tolerated |coefficient|, worst state)

| | structure `δI` | H₂ density | H_LJ min (state) | H_R′ min | natural size / evidence |
|---|---|---|---|---|---|
| A | `(a/m)P(cosh β)` | ring | **80.3** (J = 344, gap) | 6.5 (b̂ = 2000) | measured O(1/m) term is ν-flat: `a = +21.0` (J = 0), J-independent 20.9 (J = 28–100) [NE]; **positive**, i.e. a uniform rescaling, harmless at any size |
| B | `(a/m²)ν²P` = `−(a/m²)(Δ_β+¼)P` | ring derivative | 37 (J = 26) | 2.2 (b̂ = 6.9) | O(1/m²) amplitude variation |
| C | `(a/m)β∂_βP` | ring displacement | 67 | 10.6 | κ-model phase: a ≈ mβ²/24 ≲ 12 at J = 344 |
| D | `(a/m)Re P_{½+iν}(cosh β)` | **inside the disc** t ≤ β | **10.8** (J = 344) | 0.46 (b̂ = 2000) | index-shifted Regge terms (HLJ §5.2) |
| E | `(a/m)P·¼/(ν²+¼)` | ring ∗ H₂ resolvent | 12.8 | 1.9 | rational ν-dependence of Regge residues |
| F | `a/m²` (a in c₂) | `t ≲ 1/Λ` | **0.11 c₂** (J = 344) | 0.011 c₂ | "uniform smooth" class; excluded to all orders on ν ≥ ν_B by AG |
| G | `(a/m²)e^{−(νβ/3)²}` | `t ≲ β/3` | **379 c₂** | 324 c₂ | window-smooth bump; = HLJ's `T_B` class |

Identity B (exact, from the H₂ Legendre equation `(Δ_t + ν² + ¼)P = 0`) checked to ≤ 8·10⁻⁴ by direct quadrature at b̂ ≤ 530.

### 3.4 The (1100, 100) "violation" in this language
HWIN_THIMBLE §4.2: C2 `+−` at ν = 16.3 is `−2154 ± 43` c₂ against the sign-free `T_B = 1409`; the excess is `D_h`-shaped with
`c₀ ≈ 3–6.5·10³`, remainder `≲ 550` (≤ 800 at 2σ). In units of §3.3 (`Γ = c₀D_h` ⇔ `δI = (c₀/m²)·env·D_h ≈ (c₀/m²)·2P` in the window):
**ring amplitude a = c₀/m ≈ 2.7–5.9 against a_tol[A] = 200** (ψ₂ convention): cost ≤ 3 % of the floor [P given the c₀D_h decomposition,
which is a fit, NE]. The remainder 550 c₂ as a worst-case window bump: a_tol[G] = 815 → L ≈ 1.5. *So in the ψ₂ convention the sign-free
"violation" is a mis-chosen norm, not a failure.* But see N2/§3.6.

### 3.5 Curvature in the floor at large b̂ (`cache/curvfloor_1000.txt`)
At b̂ ≥ 300 the hybrid is curvature-dominated; with floor `N_flat/5 + curv` (computed curvature [NE], certified positive on [91.2, 2Λ]
for m ≥ 1.5Λ) every tolerance of §3.3 scales up by **×3.1 (b̂ = 300) … ×13 (530) … ×40 (800) … ×460 (b̂ = 2000)**. E.g. J = 344 gap:
a_tol[F] 0.11 → 1.45 c₂, a_tol[D] 10.8 → 142, a_tol[A] 80 → 1048; over all sampled states with b̂ ≥ 300 (both regions) the minima become a_tol[A] ≥ 549, a_tol[D] ≥ 125, a_tol[F] ≥ 1.2 c₂.

### 3.6 The de-shift caveat (N2) (`ls_xterm.py`, `cache/xterm_1000.txt`, `t_xsigned.py`)
The k = 2 part of Ω^EE is `Re Σ_k dν X_k Ψ̂_{EE,ν_k}/(h₋+h₊)` with the min-norm Abel weight `X` acting on the **raw** image; the hybrid's ψ₂
form is exact for the leading-Regge raw image (de-shift identity, LRC\*) and, by linearity, for any correction that is a ν-independent
multiple of it. `X` does **not** act on the stripped ring: `∫(X − ψ₂)P(cosh β)dν/Λ ≈ −2·10⁻³ … +3·10⁻³` for β = 0.05–0.5, ≫ `Ω_LR,2(β)`
(`t_xsigned.py`) — so per-ν "stripped" corrections are not well defined below ν ≈ 25 without the raw form. Counting the mismatch sign-free
(`+2Σ_{ν<30}|X − ψ₂||S|dν/Λ`, as HLJ §2 did only on `[ν_B, 25]`):

| structure | H_LJ min, ψ₂ → with X-term | at (1100, 100) | H_R′ min |
|---|---|---|---|
| A (ring amplitude) | 80 → **0.31** (J = 344); 206 → 14.5 (J = 50) | 200 → **3.2** | 6.5 → 0.031 |
| D | 10.8 → 0.27 | — | 0.46 → 0.011 |
| E | 12.8 → 3.7 | — | 1.9 → 0.31 |
| G (window bump) | 379 → 192 c₂ | 815 → **628** c₂ | 324 → 150 c₂ |
| F | unchanged (0.11 c₂) | 1.9 c₂ | unchanged |

**Reading.** If the `D_h`-shaped excess at (1100, 100) is a ν-independent relative factor of the raw image over ν ≲ 25 (where `X ≠ ψ₂`;
above 25, `|X − ψ₂| ≤ 3·10⁻³ψ₂`), the de-shift identity makes it exact and the ψ₂ column applies (harmless, ×30–70). If its ν-variation
inside ν < 25 is arbitrary, the sign-free bound (3.2) is **below** the measured `a ≈ 2.7–5.9`, i.e. the refined H_win is not established at
(1100, 100) by this bound. HLJ's `T_B` and the "remainder inside T_B" statement (HWIN §4) were computed with ψ₂ in the window; with |X|
the margin is ×1.1 instead of ×1.5–2.5. **Action:** check in `fw_gamma`/`pw_gamma` whether `Γ` is defined relative to the raw image weighted
by X, and recompute HLJ's `W_B` with |X| on `ν < ν_B`. Not done here (out of scope: the HIGHSPIN/HWIN pipelines were not touched).

---

## 4. What is derived about the correction's impact-parameter density

| piece | statement | grade |
|---|---|---|
| leading Regge | density = ring δ_β exactly | [P] |
| flat partial wave | `P_J(1 − 2x²) = ∫₀^∞J₀(xτ)J_{2J+1}(τ)dτ` (0 ≤ x < 1; Weber–Schafheitlin, DLMF 10.22.56; checked numerically at n = 3, x = 0.4), so in `t = τ/m` the density is `m J_{2J+1}(mt)`, ring at `t = (2J+1)/m ≈ β`; weight inside `(1−ε)` of the ring ≤ `(2J+1)(1−ε)K(1−ε)^{2J+1}`, `K(z) = z e^{√(1−z²)}/(1+√(1−z²))` (Kapteyn, DLMF 10.14.5): **ε = 0.2: 6.6·10⁻³ (J = 50), 1.2·10⁻⁶ (100), 1.6·10⁻¹⁰ (150), 7·10⁻²⁶ (344)**. The "centrifugal barrier" of the flat partial wave is exponential but has an Airy width ≈ `(2J+1)^{−2/3}`: ε = 0.1 is not small at J ≤ 100 | [P] |
| bulk `ν ≥ ν_B`, all orders of the saddle expansion | every term is a travelling wave of local frequency ≥ 0.97β ⇒ density at `t ≥ 0.97β`; no F-class content | [D] (HLJ §3.6); checked near the gap (HLJ) and **now on H_R′**: (2Λ, J = 50 / 150), (3Λ, 150), (5Λ, 250), (1.5Λ, 225): Newton converges on ν ∈ [3, Λ], base Hessian ≥ 1.4·10³, `min Re κ̃/β = 1.0018 / 1.0018 / 1.0012 / 1.0007 / 1.0024`, vertical decay ≥ 1.0030β at y ≤ 30/β (`cache/agscan_1000.txt`) [NE] |
| beyond all orders | no other thimble with low frequency | [NE + M] (HLJ §3.4; not extended) |
| O(1/m) term in the window | ν-flat (structure A, a = +21) at J = 0 and J = 28–100 (`J²/Δ ≤ 1`, ν < 12); no O(m)-growing ν-dependent c₂ content between Δ = 1100 and 2·10⁴ ⇒ D/E-type O(1/m) coefficients ≲ 0.05 there, vs a_tol ≥ 10.8 | [NE] (PE §2.3, HLJ_CLOSE §1) |
| O(1/m²) term in the window | ν-dependent, ≤ 450 + 100J c₂ on ν < 12 (`J²/Δ ≤ 1`); its density (ring vs disc vs window-smooth) is **not determined**; the engine series is asymptotic at `J²/Δ ≳ 2` | [I] in form; magnitude [NE] only for `J²/Δ ≤ 1` |
| exact scalar skeleton (laboratory) | window correction 6.5–10 c₂, J-independent (HLJ §4.1): ×38 below a_tol[G] at J = 344 even as a worst-case bump; bulk corrections integrate to ≤ 10⁻³ c₂-equivalent (no inside-barrier weight) | [NE, exact] |

---

## 5. Consequences for the three hypotheses

### 5.1 H_LJ(1000)
**H_LJ(1000) ⇐ H_AG (all orders [D]; beyond all orders [NE + M]) ∧ R_win ∧ H_edge-amp (unchanged, HLJ §5.1)**, with
**R_win**: *in the Regge window, the correction decomposes into a ring part (structures A–C, tolerated at a ≥ 37–80 in the ψ₂ convention)
and a non-ring part whose cost is below the floor; sufficient: non-ring part ≤ a_tol[G] (379 c₂ at J = 344, 815 at J = 100; ψ₂) as a bump,
or D/E-type coefficients ≤ 10.8 / 12.8.*
* Compared with HLJ's `H_win` (`|Γ − κ| ≤ T_B`, violated at (1100, 100)) R_win is the *derived* form of the "refined" H_win: the ring class is
  not removed by hand; its tolerance is computed [P given the class].
* Evidence: exact scalar ×38; graviton J = 64: inside; J = 100: remainder ≲ 550 vs 815 (ψ₂) / 628 (|X|), L ≈ 1.1–1.5; J ≥ 150: none.
* **Grade unchanged, [D-sketch + M]**: the [D-sketch] (window correction is ring-shaped) is not upgraded (N3: the HLJ §5.2 argument yields
  index-shifted functions whose density lies inside the disc), and the magnitude at J ≥ 150 is [M]. With the de-shift caveat N2 open, the
  (1100, 100) data support R_win only at L ≈ 1.1–1.5, not 2.7–4.8 as quoted in PE §1.1 D1 for the refined window.

### 5.2 H_R′ (m ≥ 1.3Λ, b̂ ≥ 5)
* **Hybrid floor**: [P + DV] (Lemma C′) for b̂ ≤ 119; [certified] (ROUTEA) for b̂ ∈ [91.2, 2Λ], m ≥ 1.5Λ; [NE] for m ∈ [1.3, 1.5)Λ,
  b̂ ∈ (119, b̂_max] and for b̂ > 2Λ.
* **δ at b̂ ≥ b̂_A(m)** (8.5–15): the H_LJ decomposition applies verbatim (AG now checked on H_R′ states). Window tolerances are no smaller
  than at the gap: a_tol[G] ≥ 324 c₂ (150 with |X|) over the 89 states, and grow ∝ m² at fixed b̂ (e.g. 2.6·10⁴ c₂ at 5Λ, b̂ = 120).
  **Grade [M] → [D-sketch + M]**, the same as H_LJ(1000), with the same residual R_win (plus N2).
* **5 ≤ b̂ < b̂_A(m)**: unchanged [NE + M] (bulk-model dominated; §5.3 gives the explicit bulk factor per band).

### 5.3 H_bulk (`ls_strip.py`, `cache/strip_1000.txt`)
Minimum tolerable bulk factor `f_bulk` on the ν ≥ 20 envelope `3(450+100J)` c₂ (engine-covered states `J²/Δ ≤ 1`; Lemma-C floor for
b̂ < 10), and the equivalent **relative** tolerance `ε_bulk = f·g(J)/m²` (fraction of the flat envelope `E_k(x)`):

| m/Λ \ b̂ | [0,5) | [5,7.5) | [7.5,10) | [10,15) | [15,20) | [20,30) |
|---|---|---|---|---|---|---|
| 1.3 | 84 / **9.8 %** | 10.4 / 1.6 % | 4.6 / 0.77 % | 0.53 / 0.13 % | 0.13 | 0.024 |
| 2 | **76** / 4.9 % | 23 / 1.8 % | 4.8 / 0.49 % | 1.2 / 0.17 % | 0.27 | 0.055 |
| 5 | 239 / 4.5 % | 44 / 1.2 % | 12.8 / 0.44 % | 3.9 / 0.19 % | 0.86 | 0.27 |
| 20 | 980 / 3.9 % | 282 / 1.6 % | 58 / 0.45 % | 24 / 0.26 % | — | — |

* On S_D the bulk correction may be **≥ 3.9 % of the flat envelope** at every mass: the only measured O(1/m) correction (the ν-flat 21/m,
  1.6 % at 1.3Λ) is carried separately (×1) and the J = 0 engine's bulk content is ≤ 40 c₂ (≤ 2·10⁻⁵ relative). Strip family: `f ≥ 10` for
  b̂ < 6.9, 6.0, 7.5, 8.3, 10.6, 13.7, 13.8 (m/Λ = 1.3 … 20); `f ≥ 3` for b̂ < 8.5, 8.7, 11.5, 13.0, 15.0, 15.3, 13.8.
* **Why no analytic route** (§6(f)): at b̂ < 5 the state has β ≤ 0.005, so the impact-parameter separation of §3 is empty (ring and origin
  coincide on the scale of Ω_LR); the bound needed is a relative bound on the spinning image in the bulk-point regime `ν ~ xm`, i.e. the
  fixed-x 1/m expansion of the graviton heavy image, which does not exist in the project (the 1/Δ engine is a fixed-ν Regge expansion and
  its helicity series is unconverged at ν ≥ 20). The scalar-skeleton analogue (CHMRS series, `residual_coeffs`) has `r₁res ~ 10²–10³` on
  S_D (`𝒥x² ≤ 6.25`), ×10²–10³ below the tolerance, but the transfer to gravitons is H2/H4 [M].
* **Grade [M], unchanged.**

---

## 6. Approaches tried or examined, and exactly why they fail

| route | what it would give | why it fails / stops |
|---|---|---|
| (a) sign-free "flat floor − |AdS correction|" with a per-ν bound on δ | a theorem for all b̂ ≥ b\* | Structural: §3.2 — at t ≲ 1/Λ the tolerated weight is 10⁻⁷–10⁻⁹ of the ring; the relative bulk tolerance falls from ~5 % (b̂ < 5) to 10⁻³ (b̂ ≈ 15) and 3·10⁻⁵ (b̂ ≈ 30) (§5.3). Works for the **curvature** (Lemma C′, b̂ ≤ 119) because its density is known; cannot work for δ. |
| (b) Regge-jet / endpoint expansion at large b̂ (`∫ψ δI` determined by the ν → 0 jet of δI) | δΩ = O(1/m)·Ω for any analytic relative correction | True only for corrections that are AG (ring-supported). The window content of δ is exactly what is unknown; the endpoint expansion of a window-smooth term is the G class (tolerance 192–379 c₂), not O(1/m)·Ω. |
| (c) "subleading Regge terms are `P_{l+n}(cosh β)` of the same β" (HLJ §5.2) ⇒ ring-shaped | R_win [D] | **Fails as stated**: `P_{l+n}(cosh β)` is entire of exponential type β (Paley–Wiener ⇒ density in `t ≤ β`) but its density is **not** on the ring (`P_{½+iν} = P_{−½+i(ν−i)}`; the c-function `c(ν−i)` has a pole in the upper half-plane, so its travelling-wave split is not AG). Cost per unit coefficient: structure D, a_tol ≥ 10.8 (H_LJ, ψ₂). Harmless only if the coefficients are ≲ 1–10; that needs the spinning Regge kernel at O(1/m²). |
| (d) light-cone / large-spin asymptotics at fixed twist | J → ∞ control | At Δ_g = 1000 the H_LJ region (J ≤ 344, m ≤ 1.3Λ) has β ≤ 0.53: it is the **Regge**, not the light-cone, regime. The light-cone regime (β ≫ 1) occurs only at b̂ > 2Λ, m > 2Λ, where `Ω_LR,2 ∝ M₁(β) > 0` [certified] but δ is equally unknown. |
| (e) 1/Δ engine (`J²/Δ ≤ 1`, [D given RD]) | window corrections with an error bar | Covers ν < 12–20 only for `J ≤ √Δ`, i.e. `b̂ ≲ 2Λ/√m`; never the window at large b̂ (PE §2.4) and never the bulk at J > 0. |
| (f) bulk-point (fixed x) 1/m expansion for H_bulk | relative bound on the ν ≥ 20 correction | Not available for the spinning image (§5.3). |
| (g) exact H₂ "barrier" for the AdS image (density exponentially small inside (1−ε)β) | R_win [D] | Holds for the flat partial wave [P] and to all orders of the saddle expansion on ν ≥ ν_B [D]; in the window it is equivalent to R_win itself; a non-perturbative proof needs the thimble decomposition of the crystal cycle *in the window*, where HWIN found the per-wave method biased (10²–10³ c₂) and the fold blocks ν < 10. |

---

## 7. Revised Δ_g = 1000 positivity inventory (Result D)

| region | before (UT3 §5.2, PE, HIGHSPIN §6) | now | basis |
|---|---|---|---|
| flat limit | [certified] | [certified] | CERT_CLOSE |
| near gap J = 0, 2, 4, 5, 6, 8, 12, 16 | [NE±] | unchanged | NEARGAP2, HIGHSPIN |
| near gap J ≤ 32 others, J = 29 | [NE + M] | unchanged | HIGHSPIN |
| near gap J = 36, 40, 48; even J = 26–46 without rows | [NE + M] / [M] | unchanged | — |
| **near gap even J ≥ 50, odd J ≥ 33 (H_LJ)** | [D-sketch + M], L 2.7–4.8 | **[D-sketch + M]**; residual = R_win (§5.1); L ≈ 1.1–1.5 at (1100, 100) once the remainder is counted as a window bump with |X|, ≥ 30 for the ring part in the ψ₂ convention; **N2 open** | §3–§5.1 |
| `m ≥ 1.3Λ`, `b̂ < 5` (S_D) | [NE±], L ≥ 8.7, bulk ×76 | [NE±]; bulk tolerance restated: ε_bulk ≥ 3.9 % of E_k(x) | PE §2.4; §5.3 |
| `m ≥ 1.3Λ`, `5 ≤ b̂ < b̂_A(m)` | [NE + M] | [NE + M]; explicit bulk factor per band (≥ 3 up to b̂ ≈ 8.5–15) | §5.3 |
| **`m ≥ 1.3Λ`, `b̂ ≥ b̂_A(m)` — hybrid floor** | [P] b̂ < 10; [certified] [91.2, 2Λ]; [NE] elsewhere | **[P + DV] for b̂ ≤ 119 (all m)**; [certified] [91.2, 2Λ] (m ≥ 1.5Λ only); [NE]: m ∈ [1.3,1.5)Λ & b̂ > 119, and b̂ > 2Λ | §2; ROUTEA |
| **`m ≥ 1.3Λ`, `b̂ ≥ b̂_A(m)` — δ (H_R′)** | [M] | **[D-sketch + M]** (AG [D] + checked on 5 H_R′ states; residual R_win, tolerances ≥ 324 / 150 c₂) | §4, §5.2 |
| OO block J ≥ 12 | [D given identity] | unchanged | PE D5 |

**Letter wording (suggested).** "At large impact parameter the per-state action is an average of the leading-Regge action over the state's
transverse impact-parameter density (Lemma I). The leading-Regge part is bounded rigorously up to b̂ ≈ 120 and certified beyond; corrections
supported near the state's own impact parameter are tolerated at the 6–20 % level, while corrections supported near zero impact parameter
would be amplified by up to 10⁷. Whether the AdS corrections inside the Regge window are of the first kind is the one model input at large
spin; it is supported by the exact scalar skeleton (×38) and by graviton data at J ≤ 100 (×1.1–1.5), and is not derived."
Do not write "proved for J ≥ J\*" or "for b̂ ≥ b\*": no such statement is established for the full action.

---

## 8. Corrections to existing documents (not applied)

1. **UT3 §5.2 row "m ≥ 1.3Λ, b̂ ≥ b̂_A(m)"** and **PE §2.1/§2.4/§2.6**: "hybrid floor [certified] on b̂ ∈ [91.2, 2Λ]" holds for **m ≥ 1.5Λ**
   (ROUTEA §4.2: `|w₃| ≤ 4/9`). Add: "[P + DV] on b̂ ≤ 119 for all m ≥ 1.3Λ (LARGESPIN §2)"; PE's "Lemma C: b̂ < 10" can be replaced by Lemma C′.
2. **PE §1.1 D1, UT3 §5.6 (H_LJ leverage 2.7–4.8)**: the leverage is quoted for the refined window in the ψ₂ convention; with the de-shift
   weight counted sign-free in the window it is ≈ 1.1 at (1100, 100) (§3.6). Until the Γ convention is checked, quote "1.1–4.8".
3. **HLJ_ANALYTIC §2 / §5.2**: (i) `W_X` should also be applied on `ν < ν_B` (the window), not only on `[ν_B, 25]`; (ii) the §5.2 [D-sketch]
   "P_{l+n}(cosh β) of the same β ⇒ leading-Regge-shaped" is not sufficient (§6(c)); the c-class tolerance there is the ring class (A–C) only.
4. **HWIN_THIMBLE §0 item 5 / RESULT_D update box**: "after removing the D_h class the remainder is inside T_B" — inside the ψ₂-window `T_B`;
   with |X| the margin at (1100, 100) is ×1.1.
5. **PE §2.2 / §0 "no sign-free argument can work"** — correct for per-ν bounds on δ; Lemma C′ shows it is false for the curvature (sign-free
   up to b̂ ≈ 119).
6. **UT3 §5.6 H_R′**: "[M]" → "[D-sketch + M] at b̂ ≥ b̂_A (AG checked on H_R′ states; residual as H_LJ)".

---

## 9. Reproduce, caveats

    cd ads4/code/largespin && bash run_all.sh      # ≈ 1.3 CPU-h; ≤ 1 process at a time as written

| step | output |
|---|---|
| `ls_barrier.py` | `Ω_LR,2/3(t)` table, sign and monotonicity → `cache/barrier_1000.txt`, `omlr_1000.npz` |
| `ls_struct.py` | tolerance profile, structures A–G, identity-B check → `cache/struct_1000.{txt,json}` |
| `ls_xterm.py` / `t_xsigned.py` | de-shift term (§3.6) → `cache/xterm_1000.txt` |
| `ls_curvfloor.py` | floor with computed curvature (§3.5) → `cache/curvfloor_1000.txt` |
| `ls_curvbound.py` | Lemma C′ (§2) → `cache/curvbound_1000.txt` (`curvbound_analyticC_partial.log`: analytic-constant version) |
| `ls_strip.py` | §5.3 → `cache/strip_1000.txt` |
| `ls_agscan.py 1000 2:50,2:150,3:150,5:250,1.5:225` | AG on H_R′ → `cache/agscan_1000.txt` |

Caveats.
1. Only Δ_g = 1000 was treated. Only the k = 2 image enters the structure tables (the k = 3 piece has weight u ≤ 0.59 and `Ω_LR,3` of
   comparable size; its tolerances are of the same order, not tabulated). Off-diagonal (`++`/`+−`) corrections are not included in §3.
2. States are sampled (42 + 89 for §3; 7 masses × 26 b̂ for §2); nothing is certified between samples.
3. The envelopes `E₀, E₁` of Lemma C′ use floating-point values of `Φ` on the grid (with rigorous Lipschitz pads): [DV], not [certified].
4. §3.6 uses the published grid weight `d_abel_X_1000.txt`; the continuum weight differs by ≤ 4·10⁻⁷ of max ψ.
5. The AG checks reuse HLJ's scalar saddle code; graviton = same exponent [D] (HLJ §1.4), not separately checked.
