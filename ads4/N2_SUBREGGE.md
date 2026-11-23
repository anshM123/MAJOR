# N2_SUBREGGE: which weight multiplies the window data (N2), and the first subleading Regge term projected by Ω^EE

*23 Nov 2026. Unreviewed, not committed. No existing file was edited. `HIGHSPIN.md` and `code/highspin/` were not touched.*
* *New code:* `ads4/code/n2sub/`, with caches in `cache/` (its `.gitignore` ignores everything except itself). `run_all.sh` reproduces
  everything in about 3 CPU-minutes.
  * `n2_common.py`: the published X, Y, ψ₂ and R∓; the chains' local de-shift operator 𝒜; P_{−½+iν}(cosh β) for complex ν; D_h; env; κ.
  * Checks: `t_conv.py`, `t_deshift.py`.
  * N2: `n2_tb.py`, `n2_table.py`, `n2_refit.py`, `n2_refit_b.py`.
  * Task 2: `s_chmrs.py`, `s_cost.py`.
* *Read-only inputs:*
  * `prl/data/d_abel_X_1000.txt`, `kernels.json`;
  * `spinfunc/cache/pw_gamma.pkl`, `sm_exact_{300,1000}.pkl`;
  * `largespin/cache/struct_1000.json`;
  * read-only imports of `spinfunc/ee1000_common.py`, `u2_hclosed_struct.py` and `largespin/ls_core.py`.
* *Literature:* CHMRS arXiv:2106.10274, LaTeX source (`sections/sec_heavy.tex`, §3.3–3.4, eq. `regge subleading`, which is (3.20) in
  the printed numbering). Downloaded to the session scratchpad.
* *Compute:* < 0.1 CPU-h, one process at a time, `ulimit -v 3000000`, `OMP_NUM_THREADS=1`.

Grades as in UNIFIED_THEOREM_v3 §0. Units: `c₂` = `m² ×` relative correction, as in the rows. Floors are named every time they are
used. "HLJ floor" is `0.9 λ₀ N_s` (HLJ_ANALYTIC's `T_B`). "LS floor" is `N_flat/5` (LARGESPIN §3).

---

## 0. Summary

| item | result | grade | effect on H_LJ(1000) |
|---|---|---|---|
| **N2-a: what the fw/pw/l9 Γ are** | The window data are **raw images**, divided by the forward normalisation `R₋+R₊`. Only the leading-order flat image is de-shifted (by the local operator 𝒜), and the κ image is subtracted **stripped**. The exact k = 2 action of a row correction is `Re∫X·env·Γ_c/m²`; the ψ₂ form is a local approximation. The published X de-shifts exactly the chains' 𝒜: `Re Σ X 𝒜[P_β] = Σ ψ₂ P_β` to ≤ 2·10⁻¹² relative, for β = 0–1 (§1). | [P] (definitions and code) + [DV] | **N2's premise is correct**: X, not ψ₂, multiplies the data |
| N2-b: Γ is singular at the real zeros of `R₋+R₊` (ν = 1.1714, 6.6428) | Explains "local Γ singular at every J" (HLJ_CLOSE:295). A sign-free bound on `sup|Γ − κ|` over a window that contains these zeros is infinite unless the raw correction vanishes there. The regular object is `aΓ_c`, with weight `|Y|` (`XΓ_c = Y·aΓ_c`). | [D] | H_win must be restated (§2.2) |
| **N2-c: T_B with the exact weight** | c = 3, worse helicity: (1100, 100) **1409 → 1439** (`|X|`) / **1152** (`|Y|`); J = 344: 1095 (published, with W_X) → 937 / 787. Over J = 48–344 the change is −28 % … +3 % against the published values. There is **no collapse.** | [NE] | leverage changes by ≤ ×0.7 |
| N2-d: LARGESPIN's "1.1" and "3.2" | **Artefacts of the bound used.** (i) The window bump was bounded by `ψ₂ + |X−ψ₂|`, which overcounts relative to `|X|` (`|X| < ψ₂` on ν = 1–1.5, 5–12). With `|X|`: G at (1100, 100) is 815 → **861**, not 628. (ii) "1.1 vs 2.7" also changes the floor (LS vs HLJ, ×4.9) and the structure (a bump extending past the window, ×2.8). Only a factor 1.06 comes from N2. (iii) The ring amplitude "3.2" is a sign-free bound on a *known* shape. Signed, a stripped ring `c₀D_h` costs `R2/R1 = −6.5` × the ring action at J = 100, so a_tol = **31** (LS floor), against the measured a = 4.5–6.6. | [DV] | the (1100, 100) ring excess is harmless in **either** reading at J = 100 |
| **N2-e: the D_h-class excess is an a_h normalisation of the raw record, i.e. exactly 𝒜-form** | (1) A per-row a_h error multiplies the whole raw record, so it produces `ε m² 𝒜[env D_h]/env` exactly [D]. (2) Test: C3H (no de-shift, common a_h) must carry `σ_h ε m²`. Results: (364, 64): +4230 ± 239 / −4151 ± 269 against +4168 / −4175 predicted. (1100, 100): dir −2019 ± 346 vs −2129; fw-only +6133 ± 2744 vs +5344; cv −3362 ± 1357 vs −5080. (3) Refitting in the raw form (𝒜-columns, κ de-shifted) collapses the (364, 64) residual from 87–93 to **6–12 c₂**, with r = **0 ± 3 / 1 ± 5**. | [D] mechanism + [NE] | ring part exact by the de-shift identity: cost = ψ₂ ring cost (a_tol 200 at J = 100) |
| N2-f: (𝒜−1)κ artefact in the published `y = Γ − κ` | κ is subtracted stripped, but it enters the raw data as `𝒜[env κ]/env`. The difference is **±120–340 c₂** at the (1100, 100) nodes and **±90–460** at (364, 64). It is part of HWIN's "remainder". | [DV] | remainders were overstated |
| **Corrected leverage at (1100, 100)** (CV, 6 nodes, worse helicity) | Remainder after the ring class in raw form: A + C class 339–383 c₂. A + B + C (CHMRS class): r = +38 ± 83 / +26 ± 85, rms at the error level (rms/err = 1.0); `|r| + 2σ + rms` ≤ 280 c₂. **L = T_B/remainder = 3.0–4.2 (A + C, central, `|Y|`/`|X|`; 1.8–2.2 at 2σ); ≥ 4.1 at 2σ with the full class (A + B + C, `|Y|`)**. Published: 2.6 (1.8 at 2σ); LARGESPIN: 1.1. | [NE] | **H_LJ(1000) leverage restored to ≥ 3 (CV)** |
| caveat (not N2): estimator | On the same 4 nodes the **dir** estimator gives r = −1310 ± 257 / −810 ± 131 (A + C, 1 dof), against cv −160 / +328. With dir, L ≈ 0.9–1.4. This is the pre-existing dir/CV systematic (HWIN §0.2). It is now the dominant uncertainty at (1100, 100). | [NE] | **reported prominently; unresolved** |
| N2 at J ≥ 150 | A ring correction in **stripped** raw form (R2) costs `R2/R1 = +27 / +108 / +228 / +219` × the ring action at J = 150 / 200 / 300 / 344. Tolerated a ≤ 6.6 / 1.4 / 0.42 / 0.37 (LS floor), against 178–80 if it is 𝒜-form. **At J ≥ 200 the ring class is harmless only in de-shift (𝒜) form.** | [DV] | new explicit hypothesis R1 (§2.5) |
| **Task 2: scalar subleading Regge term** | CHMRS (3.20) is reproduced on the project's exact scalar skeleton (Δ_φ = 3, k = 2; 9 states, J = 24–100, Δ_g = 300/1000) up to a **ring-amplitude constant −8.0…−9.1·P** (normalisation) and an O(ν⁴/m²) rest: ≤ 0.04 c₂ (1000) / ≤ 0.5 c₂ (300) on ν < 16. Every term is **ring-supported** (P, ν²P, ∂_βP at the state's β: distributions on t = β, Lemma I). The **non-ring window part is 0 at O(1/m²)**. The project's "exact scalar ≈ 6.5–10 c₂" **is** this term (ν → 0: 86.75/6 − 8.0 = 6.5). | [Q] + [DV] | scalar laboratory: R_win's non-ring part vanishes at first subleading order |
| Task 2: its action at the H_LJ states (incl. J ≥ 150) | `|δF|/floor` ≤ 1.3·10⁻² (HLJ floor) / 5.9·10⁻² (LS floor) at J = 26; 6.5·10⁻⁴ / 3.2·10⁻³ at J = 100; **≤ 6·10⁻⁵ / 3·10⁻⁴ at J = 150–344**. Direct quadrature = Lemma I (derivatives of Ω_LR at β) to 4 digits. In LS units: A a ≈ 0.006–0.008 (tol ≥ 80), B a ≈ 1.8 (tol ≥ 37), C a ≲ 0.005 (tol ≥ 67). | [DV] given [Q] | margin ×17–10⁴ for the scalar analogue |
| Task 2: graviton window data in the CHMRS class | Ring class (A, B, C) in raw form fits (1100, 100)-cv and (364, 64) to the error level with **r ≈ 0** (above). Fitted B coefficient: 3.7–5.0 c₂/ν² (1000, J = 100; tol 939); 0.11 (300). | [NE] | supports R_win at J = 64 (300) and J = 100 (1000, cv) |
| **Task 2: spinning subleading formula** | **Not derived.** The EE = Σ_b b-shifted-scalar identity (EESTEP §1.2) decomposes the light-ray **kernel**, with position-dependent `c_b(x)`. The heavy **block** of ⟨TTTT⟩ is not a sum of scalar (a, b)-blocks. CHMRS's mechanism needs the spinning three-point product `𝒯(x, y)` (⟨TTO_{Δ,J}⟩ structures, with state-dependent ratios beyond the helicity normalisation) in the symmetric frame. Structural statement: at each fixed order in 1/m², the image is a finite sum of rational(ν) × (index-shifted) harmonic functions and β-derivatives at the state's β, so there is no window-smooth (F/G) class at any finite order. | [I] formula; [D-sketch] structure | J ≥ 150 graviton magnitude stays [M] |

**Verdict on N2: (a) for the quantitative claim, with one real consequence.**
* The convention point is correct: the rows are raw, and X is the weight.
* The claimed collapse (leverage 1.1, ring tolerance 3.2) is not. It came from a triangle-inequality overbound, a floor change and a
  structure change.
* The measured D_h excess is an a_h normalisation of the raw record, i.e. exactly of de-shift form, so its exact cost is the ψ₂ ring cost.
* Recomputed with the exact weight and the raw-form fit, the (1100, 100) leverage is **3.0–4.2** (cv, central). It is ≥ 4.1 at 2σ once the
  CHMRS ring-derivative class is included. The table is in §2.4.
* What N2 does add is a hypothesis that must now be explicit. At J ≥ 200, ring-type window corrections are harmless only if they come in
  the de-shift pattern (R1). That is true for the a_h mechanism and exact for the scalar (magic-block identity), but not derived for EE
  beyond leading Regge.

---

## 1. The N2 trace: which weight multiplies which data

### 1.1 Definitions (file:line)

| object | definition | evidence |
|---|---|---|
| functional, k = 2 | `Re Σ_k dν X_k Ψ̂_{EE,ν_k}/(h₋+h₊)(ν_k)`, `X = χ(h₋+h₊)` (min-norm Abel weight) | `prl/data/d_abel_X_1000.txt:1`; `prl/data/README.md:60`; EE_FUNCTIONAL.md:143 ("pairs `Re∫XΨ̂ = Re∫χΨ_EE`") |
| Abel equation | `C[aY] − i tanh t S[bY] = C[ψ₂]/cosh t`, `a, b = (R₋ ± R₊)/n`. Equivalent to `Re∫X·[R₋P(ν+i) + R₊P(ν−i)]/(R₋+R₊) = ∫ψ₂P` for every β (Mehler: `cos((ν±i)t) = cos νt cosh t ∓ i sin νt sinh t`) | `kernels.json` `D_k2_EE.abel` |
| local de-shift of the chains | `A_h[f](x) = g₊ f(x + i/m) + g₋ f(x − i/m)`, `g± = (1 ± u)/2`, `u = (R₋ − R₊)/(R₋ + R₊)` (EE mode); shift term `𝒜[f] − f` with `f = env·D_h` (flat LO) | `l9_analysis.py:102–111`; `q3_analysis.py:108–119`; `w2_chain.py:8, 293`; `fw_gamma.py:3–5` |
| the rows' Γ | `Γ = m²[(G − a·sh)/(a env) − D_h]` ⇒ `Γ env/m² = G/a − 𝒜[env D_h]`, with `G/a` = state/reference record (reference ≈ `(R₋+R₊)·C_ref`) = raw image ÷ (R₋+R₊) | `l9_analysis.py:16, 181, 186` |
| κ image | **stripped** `κ = m²{A D(x_e) − D + 2[P − A J₀(νβ)]/env₂}` subtracted from Γ | `sm_curvfix.py:7, 66–87` (reproduced exactly: `cache/t_conv.txt`) |
| chain action | `dF = Σ p₂ env Γ/m²`, `p₂` = ψ₂ GL weights (local form) | `ee1000_common.py:102`; `w2_chain.py:163` ("the chain's local form uses psi2 node weights") |
| HLJ's T_B | `W_B = Σ_{ν<ν_B} p₂ env/(m²N_s)` (ψ₂ on the window); `W_X = Σ_{[ν_B,25]}|X−ψ₂|` only beyond the window | `lj_contour.py:10, 52, 71, 75`; HLJ_ANALYTIC.md:145–151 |
| the exact non-local density (already in the project) | `Re[(X/φ) Γ_c]` with `Γ_c` the complex forward-normalised correction, φ = ψ₂ | `z4_chain.py:9–16` |
| LARGESPIN's X-term | `ψ₂`-cost + `2Σ_{ν<30}|X − ψ₂||S|` (triangle inequality) | `largespin/ls_xterm.py:43` |

### 1.2 Numerical confirmation [DV] (`t_deshift.py`, `t_conv.py`)
* With the **published** grid X (0.04 grid) and the published `R∓`, `Re Σ X 𝒜[P_β] dν/Λ = Σ ψ₂ P_β dν/Λ`. Relative agreement:
  * `0` at β = 0;
  * `2.4·10⁻¹³` at β = 0.1836, i.e. (1100, 100);
  * `2·10⁻¹²` at β = 1.
* `kernels.json`'s `R∓` equal the chain's `u_EE` source (`u2_hclosed_struct`) to `10⁻¹⁶`.
* So X de-shifts exactly the operator the rows use. The rows are therefore raw, forward-normalised images.
* Raw read as stripped (`Re Σ(X − ψ₂)P_β`), relative to the ring action `Σψ₂P_β`:

  | β | 0.05 | 0.10 | 0.18 | 0.30 | 0.50 | 1.0 |
  |---|---|---|---|---|---|---|
  | ratio | −0.48 | −3.5 | −7.4 | +61 | +227 | −127 |

  It is concentrated at ν < 25 and vanishes exactly at β = 0.
* Real zeros of `R₋+R₊`: ν = **1.1714, 6.6428**. The weight X vanishes there; `Γ_c` has poles there for any raw correction that does not
  vanish there.

### 1.3 Consequence for the action of the data [P]
The κ point drops out of the action, so only the window weight is at stake. Write the raw row content as
`Γ_c = 𝒜[env κ]/env + δ_raw`; the de-shift identity gives `Re∫X𝒜[envκ] = ∫ψ₂ envκ`. Then:

    δF_exact = Re∫X env (Γ_c − 𝒜κ)/m² = ∫ψ₂ env (Γ − κ)/m² + Re∫(X − ψ₂) env Γ_c/m².

* The chains (`W2_GAUGE=none`) and HLJ's T_B keep only the first term on the window.
* The second term vanishes for every raw correction of the form `𝒜[env f]/env` with f strip-analytic, i.e. with an H₂ density σ satisfying
  `∫e^t|σ| < ∞` (de-shift identity + Lemma I). It is O(1) relative for any other form (§1.2).
* EE_NEARGAP §4 tested this at small J (`W2_GAUGE=min`) and found agreement to 10⁻⁴ at 1000. That test does not cover large J, where
  `N_s ∝ b̂⁻³` makes the ν < 25 window matter.

---

## 2. Recomputed window tolerances and leverage

### 2.1 T_B with the exact weight (`n2_tb.py` → `cache/n2_tb.txt`) [NE]
`T_B = 0.9λ₀/(W_win + E_bd + W_X + E_top + E_edge)`, exactly as `lj_contour.py`. Only `W_win` changes. E_bd and W_X stay as they are,
because the contour on ν ≥ ν_B deforms the analytic ψ₂ and the X − ψ₂ difference there is W_X. Minimum over the two helicities, gap
states, c = 3 (c = 5):

| J | published ψ₂ (no W_X) | published ψ₂ + W_X | **`|X|`** (tol. on `sup|Γ_c|`) | **`|Y|`** (tol. on `sup|aΓ_c|`) |
|---|---|---|---|---|
| 48 | 2456 (1314) | 2463 | 2480 (1324) | 2257 (1243) |
| 64 | 1912 (1019) | 1919 | 1943 (1025) | 1695 (933) |
| **100** | **1408 (734)** | 1393 | **1439 (746)** | **1152 (631)** |
| 150 | 1210 (624) | 1174 | 1265 (645) | 946 (496) |
| 200 | 1231 (641) | 1111 | 1170 (670) | 881 (478) |
| 250 | 1359 (712) | 1063 | 1020 (764) | 775 (519) |
| 300 | 1582 (824) | 1057 | 942 (887) | 748 (592) |
| 344 | 1847 (968) | 1096 | 937 (967) | 787 (655) |

* Window-integrated, `|X| ≈ ψ₂`. Per band, `|X|/ψ₂` is: 6.1 (ν < 0.5), 0.14 (1–1.5), 1.9–2.2 (1.5–5), 0.05–0.6 (5–7), 0.3 (7–10),
  1.03 (12–15). ψ₂ ∝ p makes the small-ν excess cheap.
* `|Y| ≥ |X|` because `|a| < 1`: `|a| = 0.67 / 0.77 / 0.86 / 0.94` at ν = 13.3 / 16.3 / 20.8 / 33.4.
* **Net: T_B moves by −28 % … +3 % against the published values. N2 does not collapse the window tolerance.**

### 2.2 Restating H_win
* `sup_{ν<ν_B}|Γ − κ|` is infinite for generic raw corrections, because of the poles at 1.17 and 6.64 (HLJ_CLOSE excluded these by
  hand).
* **Proposed form:** `sup_{ν<ν_B} |a(Γ_c − 𝒜κ − ring class in 𝒜-form)| ≤ T_B^{|Y|}`. This is regular, and it is the object the exact
  weight sees.

### 2.3 Structure tolerances, both floors (`n2_table.py` → `cache/n2_table.txt`) [DV]
Columns:
* **G bump**, `a e^{−(νβ/3)²}/m²`, LS floor: `ψ₂` / LARGESPIN (`ψ₂ + |X−ψ₂|`) / **`|X|` (correct sign-free)**.
* **Ring A**, `(a/m)P`, LS floor: **R1** = raw is 𝒜-form (exact, equal to LARGESPIN's ψ₂ value) / **R2** = raw is the stripped shape
  (signed exact `Re Σ X env P`).

| J | G: ψ₂ / LS / **`|X|`** | A: **R1** / **R2** (R2/R1) |
|---|---|---|
| 26 | 2842 / 2755 / 2860 | 207 / 284 (+0.73) |
| 50 | 1480 / 1347 / 1510 | 207 / 98 (−2.1) |
| 64 | 1183 / 1029 / 1220 | 207 / 43 (−4.8) |
| **100** | **815 / 628 / 861** | **200 / 30.9 (−6.5)** |
| 150 | 604 / 401 / 642 | 178 / 6.6 (+27) |
| 200 | 502 / 299 / 508 | 148 / 1.38 (+108) |
| 250 | 442 / 244 / 409 | 120 / 0.64 (+189) |
| 300 | 403 / 211 / 337 | 97 / 0.42 (+228) |
| 344 | 379 / 192 / 290 | 80 / 0.37 (+219) |

* LARGESPIN's sign-free A values: 3.2 (J = 100), 0.31 (J = 344).
* The D_h class itself, per unit `c₀`, at (1100, 100) `+−` (`n2_tb.txt`, HLJ floor): R1 `c₀_tol = 1.2·10⁶`; R2 `1.6·10⁵`. Measured
  `c₀` = 4.9–7.3·10³.

### 2.4 Corrected leverage at (1100, 100) (`n2_refit.py`, `n2_refit_b.py`) [NE]
Fits on ν ∈ [12, 34] of the HWIN data (`pw_gamma.pkl`; 6 nodes for cv `_w`, 5 for fw-only `_none`). The fits are labelled as follows:
* **F0** is fw_gamma's `cs` fit (stored convention);
* **F1** is the same fit with 𝒜-columns (`𝒜[env φ]/env`) and κ de-shifted;
* **F1+B** adds `𝒜[env ν²D]/env`.

| data (worse helicity) | F0: \|r\| + rms | F1: \|r\| + rms | F1+B: r ± σ, rms (rms/err) | T_B (ψ₂ / `|X|` / `|Y|`) | **L** (F0 → F1 → F1+B, 2σ) |
|---|---|---|---|---|---|
| (1100, 100) cv `++` | 545 | 383 | +38 ± 83, 76 (1.00) | 1435 / 1441 / 1158 | 2.6 → 3.0–3.8 (1.8 at 2σ) → ≥ 4.1 |
| (1100, 100) cv `+−` | 268 | 339 | +26 ± 85, 78 (1.02) | 1409 / 1439 / 1152 | 5.3 → 3.4–4.2 (1.8 at 2σ) → ≥ 4.2 |
| (1100, 100) fw-only `++ / +−` | 479 / 248 | 372 / 303 | −56 ± 83 / −61 ± 87 | as above | → 3.1–3.9 → ≥ 3.9 |
| (364, 64) cv `++ / +−` (300; T_B(+X) = 255) | 136 / 141 | **7 / 13** | **+1 ± 1 / +1 ± 5**, 3 / 10 | 255 | 1.9 → 20–36 → **≥ 12** |
| (1100, 100) **dir**, 4 nodes (1 dof) `++ / +−` | 1731 / 907 | 1426 / 869 | — (too few nodes) | — | **0.8–1.3** |

* For the A + C leverages, `L = T_B^{|Y|}/(|r|+rms)` … `T_B^{|X|}/(|r|+rms)`. The 2σ entries use `|r| + 2σ_r + rms` with `T_B^{|Y|}`.
* **Corrected H_LJ(1000) leverage at (1100, 100): 3.0–4.2 (cv, A + C class, central; 1.8 at 2σ), ≥ 4.1 at 2σ with the full ring class.**
  This replaces 2.7 (PE D1) and 1.1 (LARGESPIN §3.6).
* **Negative, prominently:**
  * the dir estimator (same 4 nodes: r = −1310 ± 257 / −810 ± 131 against cv −160 ± 10 / +328 ± 108) gives L ≈ 0.8–1.4 (A + C);
  * the dir − cv difference of ≈ 10³ c₂ in `r` is the dominant uncertainty at (1100, 100);
  * it is unrelated to N2;
  * the dir C3H `+−` check (−2019 ± 346 against the predicted −2129) is the best of all (1100, 100) C3H checks.

  Deciding between dir and cv needs more window nodes with both estimators (≥ 6 per helicity).

### 2.5 What N2 adds (hypothesis R1)
* **R1:** "ring-type window corrections enter the raw rows in the de-shift pattern `𝒜[env·(ring)]/env`."
* **Status:**
  * a_h normalisation: [D];
  * the C3H σ_h test and the raw-form fits: [NE] at J = 64 (300) and J = 100 (1000);
  * scalar: exact (CLS-D = 𝒜_g[C_CHMRS], magic-block identity, H0 §20.2) [Q];
  * EE beyond leading Regge: [I].
* **Why it is load-bearing at J ≥ 200:** a stripped-form ring term is tolerated only at a ≤ 0.4–1.4 there (§2.3).

---

## 3. Task 2: the first subleading Regge term, as projected by Ω^EE's k = 2 part

### 3.1 The formula [Q]
CHMRS `sec_heavy.tex`, eq. `regge subleading`, with `a = 2Δ_φ + k − d/2` and `m² = (Δ−J−d+1)(Δ+J−1)` (= the project's m²):

    (m^{2k}/2) C_{k,ν} = P(η) + (1/6m²){ P′(η)[(η²−1)(ν² − 3a² + 6a(2Δ_φ−d) − 2) − (d−4)(d−2)]
                                        + ηP(η)[(ν² − a(a+d−2))(7d/2 + 3a − 1 − 6Δ_φ) − 2(Δ_φ−d−1)₃ − 2(a+(d−2)/2)(1 − a² − a(d−2)/2)] } + O(m⁻⁴).

For d = 3, Δ_φ = 3, k = 2 (a = 6.5) this becomes

    δ = (1/6)[(86.75 + 11ν²) cosh β · P + (1/sinh β + sinh β (ν² − 11.75)) ∂_βP].

The CHMRS derivation (symmetrical frame, rescaling `x, y → x/m, y/m`, Fourier–Laplace identity) makes every order a polynomial in
η and P^{(n)}(η) with ν-entire coefficients.

### 3.2 Check against the exact scalar skeleton (`s_chmrs.py` → `cache/s_chmrs.txt`) [DV]
`m²(ex₂/2 − P)` (exact, `sm_exact`) against δ, on nodes ν < 40:
* **Residual = c·P** with `c = −8.01, −8.03, −8.05, −8.13` (Δ_g = 1000, J = 40, 48, 64, 100) and `−8.07, −8.18, −8.50, −9.06`
  (Δ_g = 300, J = 24, 40, 64, 100). This is a ring-amplitude 1/m² term, a normalisation-convention difference (e.g. of m² or of the MFT
  density). It is not identified.
* After removing `cP`, the rest is:
  * ≤ 0.016–0.045 c₂ on ν < 16 and ≤ 2.3 on ν < 40 at 1000;
  * ≤ 0.13–0.47 on ν < 16 and ≈ 18–21 on ν < 40 at 300.

  This is the O(ν⁴/m²) next order.
* The "exact-scalar window correction 6.5–10 c₂" (HLJ_ANALYTIC §4.1, LARGESPIN §4) is this term: at ν → 0, `86.75/6 − 8.0 = 6.5`.

### 3.3 Impact-parameter content and action at the H_LJ states (`s_cost.py` → `cache/s_cost.txt`) [DV given Q]
* Every term is a distribution on the ring `t = β`. For P this is `δ_β`; for ν²P it is `−(Δ_t+¼)δ_β`; for ∂_βP it is `−δ′_β`.
* By Lemma I (extended to compactly supported distributions; Ω_LR is smooth for t > 0) the action is a fixed combination of
  `Ω_LR,2^{(n)}(β)`, n ≤ 3.
* Direct quadrature over ν ∈ [0, Λ] and the Lemma-I evaluation agree to 4 digits (J ≤ 150) and within a few % at J ≥ 250, where both
  are ~10⁻⁶ N_s.

| J (gap) | 26 | 50 | 64 | 100 | 150 | 200 | 250 | 300 | 344 |
|---|---|---|---|---|---|---|---|---|---|
| \|δF\|/HLJ floor | 1.3e-2 | 3.3e-3 | 2.1e-3 | 6.5e-4 | 5.6e-5 | 6.0e-6 | 2.4e-6 | 8.5e-6 | 1.0e-5 |
| \|δF\|/LS floor | 5.9e-2 | 1.5e-2 | 9.5e-3 | 3.2e-3 | 3.2e-4 | 4.1e-5 | 2.1e-5 | 9.5e-5 | 1.4e-4 |

* In LARGESPIN's units:
  * A: a = (86.75 cosh β/6 − 8)/m ≈ 0.006–0.008 (tol ≥ 80);
  * B: a = 11/6 ≈ 1.8 (tol ≥ 37);
  * C: a ≈ 1/(6mβ sinh β) ≲ 0.005 (tol ≥ 67).
* **For the scalar analogue, R_win holds at O(1/m²) at every H_LJ state including J ≥ 150, with margin ×17 (J = 26, LS floor) to ×10⁴.**
* Its non-ring (F/G) part is identically zero at this order.
* This replaces LARGESPIN's "exact scalar ×38 as a worst-case bump" by an exact statement at first subleading order.

### 3.4 Graviton window data in the same class (§2.4) [NE]
* The ring class (A: D, B: ν²D, C: xD′) in raw 𝒜-form, plus a constant r, fits (364, 64) with residual 3–10 c₂ and (1100, 100) cv at
  rms/err ≈ 1.0.
* **r is consistent with 0 in both** (+1 ± 5; +26–38 ± 85).
* The fitted B coefficients are 0.11 (300, J = 64) and 3.7–5.0 c₂/ν² (1000, J = 100), against the scalar 1.8 and a tolerance 939 at
  J = 100.
* **So the "check against the graviton data at J = 100 (≈ 550 c₂)":** the 550 was ring class plus convention artefacts. The non-ring
  part is `|r| + 2σ` ≤ 200 c₂ (cv), but see the dir caveat.

### 3.5 The spinning formula: precise obstacle [I] and what is derived [D-sketch]
* **Why the EE b-shift identity does not reduce it to scalar terms.**
  * EESTEP §1.2 writes the **kernel** `K^EE_ν(x) = N(ν) r² Σ_{b=−4}^{4} c_b(x) v^{−b/2} Ĝ_{−b}(x; ν)`, where `c_b` depends on position
    `(r, β)` and on the polarisations.
  * CHMRS's subleading formula is for a scalar functional acting on a **scalar** heavy block.
  * Here each b-term pairs `Ψ^{(b)}_ν` with `c_b(x)·G^{TTTT}_{Δ,J}(x; pols)`. This is not an (a, b)-shifted scalar block: it is the
    spinning ⟨TTO⟩⟨OTT⟩ block, which carries several three-point structures.
  * Reproducing (3.20) needs one of the following:
    * CHMRS's `𝒯(x, y)` for ⟨TTO_{Δ,J}⟩ in the symmetrical frame (the product of spinning three-point structures, with its own 1/m
      tower);
    * or a weight-shifting decomposition of the ⟨TTTT⟩ block into scalar blocks of shifted (Δ, J, a, b).

    Neither exists in the project.
  * In addition, beyond leading order the coefficients depend on **ratios of ⟨TTO⟩ structures** that the per-helicity a_h do not fix.
    So the spinning subleading term is not universal across heavy states. This is an extra unknown, not only an unperformed computation.
* **What survives [D-sketch].**
  * In CHMRS's mechanism every order of the 1/m expansion is (polynomial in x, y, polarisations) × `e^{x·p+y·p̄}`.
  * The kernel side is (rational in ν) × `P_{l+s}(η_{xy})` × powers of `|x||y|`. The index shifts `|s| ≤ n+1` come from `Ĝ_b`'s
    r-expansion; `B^{(b)}[n, m]` is polynomial in b and rational in ν (EESTEP §1.3).
  * The Fourier–Laplace identity then maps each term to (rational in ν) × p-derivatives of `|p|^{−a} P_{l+s}(η_p)`.
  * Hence **at each fixed order the k = 2 raw image is a finite sum of index-shifted harmonic functions and their β-derivatives at the
    state's own β**: classes A–E of LARGESPIN, never F/G (window-smooth or origin-supported).
  * Its exact X-action is a smooth function of β.
  * Whether it is harmless is decided by two things:
    * (i) the coefficients;
    * (ii) whether the index-shifted pieces come in the de-shift pattern (R1). A mismatched pattern costs up to ~200× ring at J ≥ 250
      (§2.3).
* **Status for J ≥ 150 (graviton).**
  * The magnitude is still [M].
  * The structure is [D-sketch]: there is no window-smooth class at any finite order. The finite-order coefficients therefore enter
    only through classes with computable per-unit costs (A/B/C ≥ 37–80 in 𝒜-form).
  * The residual risks are beyond-all-orders terms (already H_AG [NE + M]) and R1.

---

## 4. Corrections to existing documents (not applied)

1. **LARGESPIN_ANALYTIC §0.1 (N2 bullet), §3.4, §3.6 table, §5.1, §7 row H_LJ, §8 items 2 and 4.**
   * Replace "with |X| the leverage at (1100, 100) drops to ≈ 1.1; the ring tolerance 3.2 is below the measured excess" by: "the X-term
     of §3.6 is a triangle-inequality overbound. With the exact sign-free `|X|`, G is 815 → 861 (J = 100) and 379 → 290 (J = 344), and
     T_B is 1409 → 1439 (`|X|`) / 1152 (`|Y|`). The D_h excess is an a_h normalisation of the raw record (de-shift form, exact ψ₂ ring
     cost). A stripped-form ring costs a_tol = 31 (signed) at J = 100, against the measured a = 4.5–6.6 (N2_SUBREGGE §2)."
   * Keep the R2 caveat for J ≥ 200 (§2.5).
2. **POS_ERRORS.md update box (line 35) and HLJ_CLOSE.md:22.** "Until then quote 1.1–4.8" → "**3.0–4.2 (cv, raw-form refit, central), ≥ 4.1
   at 2σ with the full ring class; 0.8–1.4 with the dir estimator (dir/cv systematic, not N2)**".
3. **HWIN_THIMBLE §4.2–4.3, §5.**
   * The D_h class must be fitted as `𝒜[env D_h]/env`, and κ must be de-shifted (`(𝒜−1)κ` is ±120–460 c₂ at the nodes).
   * The (364, 64) remainder becomes 7–13 c₂ (was 136–141, "≈ T_B"); r = 0 ± 3–5.
   * The (1100, 100) cv remainder becomes 339–383 (A + C) and 26–38 ± 85 (A + B + C).
   * "Inside T_B" also holds for the regular `|Y|` form.
4. **HLJ_ANALYTIC §2 / H_win statement.**
   * `sup|Γ − κ|` over a window containing ν = 1.17, 6.64 is infinite for generic raw corrections.
   * State H_win on `a(Γ_c − 𝒜κ − ring class)` with `W_win = Σ|Y| env`. The §0 T_B table changes by −28 % … +3 % (§2.1).
5. **LARGESPIN §4 row "exact scalar skeleton"** and **SUBLEADING_REGGE §4.**
   * The exact-scalar window correction is CHMRS (3.20) plus a −8.0…−9.1·P normalisation term, all ring class [Q + DV].
   * Its cost is ≤ 6.5·10⁻⁴ of the HLJ floor at J = 100 and ≤ 6·10⁻⁵ at J ≥ 150 (§3.3), not "a bump ×38 below a_tol[G]".
6. **REGGE_FROM_BLOCKS §8 item 1.** Add: the EE b-shift identity decomposes the kernel, not the spinning block; the missing input is the
   ⟨TTO⟩ `𝒯(x, y)` (or a weight-shifting decomposition), and it depends on structure ratios beyond a_h (§3.5).

## 5. Caveats
1. The fits have 2 dof (F1+B, 1000) or 1 dof (dir), on 4–6 nodes, so they are not a substitute for more nodes. (364, 64) has 5 dof.
2. `|Y|` depends on the phase/modulus convention of `n` (README §4). `|X|` does not.
3. The scalar comparison uses Δ_φ = 3; the −8 normalisation term is fitted, not derived.
4. The R1/R2 dichotomy brackets the stripped/raw readings. A raw correction can also be partly 𝒜-form; its exact cost lies between the
   two.
5. Everything is at Δ_g = 1000 except the (364, 64) and 300 scalar checks. Nothing here is certified.
