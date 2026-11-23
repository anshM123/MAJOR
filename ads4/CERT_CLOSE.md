# CERT_CLOSE: the flat certificates (NF02 closed at 1/5; Q's flat statement corrected) and the k = 2 continuum weights

*23 Nov 2026. Unreviewed, not committed. No existing file was edited; `EESTEP_CLOSE.md` and `code/eestep/` were not touched.*

* *Code (new):* `ads4/code/cert/` — `nf_corner.py`, `nf_lemmaB_check.py`, `nf_exact_min.py`, `nf_Pext.py` (NF02);
  `q_cert.py`, `q_scan_bhat.py`, `q_scan_refine.py`, `q_backward.py`, `q_scan_uncov.py` (Q); `k2_cont.py`, `k2_xq_check.py`, `k2_export.py` (k = 2);
  `proto_corner.py` (float prototype used only to choose the cutoff; not a certificate); `cert_run_all.sh` (reproduction).
* *Caches/logs:* `ads4/code/cert/cache/` (git-ignored by its own `.gitignore`).
* *New data:* `ads4/prl/data/k2/` (README, `SHA256SUMS`).
* *Read-only imports:* `code/robust/{nf_cert,rcert,rcommon}.py` (the NOFORWARD / ROBUST certificate machinery, unchanged; `nf_cert`
  patches `rcert` to centred interval evaluation at import, as in its own runs), `code/deshift/{kern,d_abel}.py`,
  `code/resultq/common_q.py`; data `robust/functionals_noforward.json`, `prl/data/{q3_full_delta0.2.json,kernels.json,d_abel_X_*.txt}`,
  `resultq/cache/rv_k2eng_Y_{300,1000}.pkl`.
* *Compute:* ≈ 2.5 CPU-hours, ≤ 2 processes, `ulimit -v 3000000`, `OMP_NUM_THREADS=1`.

Grades as in UNIFIED_THEOREM_v2 §0. **[certified]** is used only for statements whose every inequality is checked in exact rational
or FLINT/Arb ball arithmetic, on top of analytic lemmas that are written out here or in the cited (refereed) documents.

---

## 0. Verdict

| item | result | grade |
|---|---|---|
| **NF02 flat floor** | **`𝔉 ≥ (1/5)N` on every flat heavy state** (both channels, all `J`, all `m ≥ M`). The uncertified corner `+−`, `J ≥ 400`, `M ≤ m < 1.633M` is closed by a new backward-endpoint lemma (Lemma B, §1.3). The floor is sharp: at `+−`, J = 5, `m = M` the ratio is exactly `22603106861933/113015528461279 = 0.2000000103` | **[certified]** (new corner code unreviewed) |
| what changes for D | "`𝔉 ≥ (3/20)N` except one corner, 0.200 float" → "`𝔉 ≥ N/5` everywhere, attained to `10⁻⁸`". **No bound value changes.** UT/RD-S's "min 0.200 [certified]" is now **true** (it was not when written) | — |
| **Q flat (REPRO D9)** | REPRO reproduced: `0.196585` (`+−`, J = 26, `u = 0.02535`), `0.19990` (`++`, J = 2). **But the true flat minimum is lower:** `0.1835` at `++`, `J ≈ 3600`, `m ≈ 100M` (`b̂ ≈ 72`), a region the earlier dense grids missed (coarse `u`-grid, `J ≤ 3000`) | [NE] |
| **Q flat, a definitional defect** | The stored `Q3-full-delta0.2` weights (expanded IEEE doubles, coefficients up to `7.8·10⁵`) carry their `(1−p)³` factor only to rounding: `ψ₂(1) = 1.8·10⁻¹¹`, `ψ₃(1) = 2.4·10⁻¹⁰`. At `+−`, `m = M` the backward peak of `d̃⁴⁴ ~ J⁸` amplifies this: **the flat action of the stored functional is negative for J ≥ 1990 at `m = M`**. With the `(1−p)³` factor made exact (`Q3*`, a `10⁻⁹` change) the minimum there is `0.2000`. | [NE] (finding); fix = definitional |
| **Q flat certificate (Q3\*)** | `FF ≥ (9/50)N` **certified** on: all `J < 400` (every mass); `J ≥ 400` with `b̂ ≤ 120` (this region contains the global minimum, `b̂ ≈ 72`); `++`, `J ≥ 400`, `m ≤ 2M`; `J ≥ 400`, `b̂ ≥ 120`, `2M ≤ m ≤ 10M`. **Not certified:** `+−`, `J ≥ 400`, `m ≤ 2M`, and `J ≥ 400`, `b̂ ≥ 120`, `m ≥ 10M` (the existing lemmas are too crude for Q's weights there, §2.3); float minima there are `0.304` and `0.997` | [certified] on the listed regions; [NE] elsewhere |
| **`X^Q` published** | `prl/data/k2/xq_abel_X_{300,1000}.txt`; independently re-solved (DESHIFT's solver): `X` to `≤ 8·10⁻⁷` of `max ψ₂^Q` (ν < 2), `≤ 10⁻¹¹` (ν > 25); V3 exact to `10⁻¹⁴`; `dk₀^Q = −1.630147·10⁻⁴ / −1.453690·10⁻⁵` | [DV] |
| **continuum k = 2 weights (D and Q)** | built in DESHIFT's shifted form, k = 2 version: **`ε_W = 3.3·10⁻⁶ / 1.1·10⁻⁶` (D, 300/1000; `1.1·10⁻⁹ / 1.1·10⁻¹⁰` of `∫|ψ|`) and `1.9·10⁻⁷ / 7.4·10⁻⁸` (Q)** for all `t ≥ 0`. Pole structure: double poles at `±i/2` (the double zeros of `R_∓`), branch points at `±i` (the `|μ|` kink of `ψ₂`). Equal to the stored grid weights to `≤ 4·10⁻⁷` of `max ψ`; V3 and `dk₀` reproduced (`dk₀` to `≤ 4.5·10⁻⁷` relative) | [DV] |
| cubic bound | **unchanged**: `R₀ = 3.13206·10⁻⁶ / 2.89126·10⁻⁸` (K = 1), and `(1+dk₀)^{−½} = 1.00038354 / 1.00003415` for both the grid and the continuum weight | [DV] |

---

## 1. NF02: `𝔉 ≥ N/5` on all flat states [certified]

Notation as NOFORWARD_FUNCTIONAL (NF): `u = M²/m²`, `N = |B₂| + |B₃|`, four sign patterns `g = (s₂, s₃)` with
`W₀ = 2(1−δ′s₂)ψ₂`, `W₁ = −(1−δ′s₂)p²ψ₂ + σ(1−δ′s₃)ψ₃`; `𝔉 − δ′N ≥ 0 ⇔ FF_g := ∫(W₀ + uW₁)D_h ≥ 0` for all four `g`.

### 1.1 Why 1/5 is the right constant
* `nf_exact_min.py` (exact Lemma-P polynomials): at `m = M` the ratio is `0.2000000103` (`+−`, J = 5), `0.2000000121` (`+−`, J = 4),
  `0.2000000302` (`++`, J = 4), and `0.20000001–0.20000006` for `+−` J = 6, 7, 8, 14, 19. The LP pinned many `m = M` states at its
  margin 0.2; the exact rounding left all of them **above** 1/5 by `≥ 10⁻⁸`. Away from `m = M` the ratio rises
  (`0.20004` at `u = 0.9995`). So `δ′ = 1/5` is the largest round constant that can be certified, and it was.

### 1.2 Coverage at `δ′ = 1/5`

| region | method (code) | result |
|---|---|---|
| P: `J < 400`, all `u ∈ [0,1]`, both channels | exact Lemma-P polynomials, Descartes (`robust/nf_cert.py … P`) | PASS, 2384 polynomials |
| M: `J ≥ 400`, `u ≤ 1/4`, `b̂ ≥ 120` | Lemma M, no forward term (`nf_cert.py … M`) | PASS (tightest: `+−`, `s₂ = +1, s₃ = −1`, `u ∈ [1/10, 1/4]`: `b̂⁴𝔉 ≥ 938` vs `C_R ≤ 2.14·10⁶`) |
| S: `J ≥ 400`, `u ≤ 1/4`, `b̂ ≤ 120` | Cauchy cells (`nf_cert.py … S±±`) | PASS, 2416 cells per channel |
| C`++`: `J ≥ 400`, `u ≥ 1/4` | Lemma M (Legendre) (`nf_cert.py … Cpp`) | PASS |
| C`+−`: `J ≥ 400`, `1/4 ≤ u ≤ 3/8` | Lemma M + Lemma C (`nf_cert.py … Cpm`) | PASS |
| **corner: `+−`, `J ≥ 400`, `3/8 ≤ u ≤ 1`** | **Lemma C + M on `w_A`, Lemma B on `w_B` (`cert/nf_corner.py 1/5 400`)** | **PASS**, all 4 patterns × 41 `u`-cells; tightest `val/main = 0.150` (`s₂ = +1, s₃ = −1`, `u ∈ [0.99, 1]`) |
| cross-check: `+−`, `400 ≤ J < 800`, all `u` | exact Lemma-P + Descartes (`cert/nf_Pext.py`) | PASS, 1600 polynomials (redundant with the corner) |
| low side (`k, A, B`, matter signs) | exact (NF §4, unchanged, `δ′`-independent) | PASS (NF logs) |

The union is all flat states (NF §3's coverage plus the corner). The regions P–Cpm are NOFORWARD's certificates (refereed in
REFEREE_NOFORWARD at `3/20`), re-run unchanged at `1/5`; all logs print the NF02 fingerprint `sha256 b8bd1237…20db4a`.

### 1.3 The corner: Lemma B and the certificate

*Split.* `w = W₀ + uW₁ = w_A + w_B`, `w_A = wχ`, `χ = 1 − I_t(7,7)`, `t = p/c`, `c = 3/5` (C⁶ at `p = 0` and `p = c`, `χ ≡ 0` on
`[c, 1]`; gluing checked exactly). `I_A = ∫w_A D`, `I_B = ∫w_B D`, `D = P^{(0,8)}_{n}(1−2up²)`, `n = J − 4`.

*Part A [certified].* `ŵ_A := w_A(1−up²)^{−8}` vanishes on `[c, 1]`, so it is uniformly `C⁴` up to `u = 1` (this is what failed for
the whole `w`, NF §4.1). Lemma C (FLAT_GRAVITON §2) and Lemma M (§3 there) give
`b̂⁴I_A ≥ γ(u)b̂/(1+ds)³ − C_R/(1−ds)⁴`, `γ = −w″(0)/2`, `ds = 8/(2J+1)`, with `C_R` from Arb Taylor jets on 1200 `p`-cells, sup over
six coarse `u`-cells. `C_R ≤ 3.6·10⁶ … 5.5·10⁶`.

*Lemma B [P; identity [DV] to `10⁻³⁰`].* Let `y(p) = P^{(0,8)}_n(1−2up²)`, `λ = n(n+9)`, `K = p(1−up²)⁹`, `ρ = p(1−up²)⁸`. Then
1. `(Ky′)′ = −4uλρy` on `(0,1)` (the Jacobi equation, DLMF 18.8.1, in the variable `p`; checked symbolically for general `n, u`).
2. For `f ∈ C²[0,1]` with `f(0) = f′(0) = 0` and `0 < u ≤ 1`:
   `∫₀¹ f y dp = −(4uλ)⁻¹[β_f + ∫₀¹ (Lf) y dp]`, with
   `Lf := d/dp[(f′ − f/p)(1−up²) + 16upf]` (so `Lp^m = (m−1)²p^{m−2} + (17−m)(m+1)up^m`) and
   `β_f = f(1)(1−u)y′(1) − [(f′(1) − f(1))(1−u) + 16uf(1)]y(1)`.
   *Proof.* Write `f = gρ`; integrate `∫g(Ky′)′` by parts twice on `[η, 1−η]`; `gK = f(1−up²)` and
   `g′K = (f′ − f/p)(1−up²) + 16upf` are polynomial expressions (checked symbolically), continuous on `[0,1]`, and vanish at `p = 0`
   under the hypothesis; let `η → 0`. ∎
3. `4uλ = b̂² − 81u`.

*Application.* `f₀ = w_B` (`= p⁸·poly` near 0, `≡ w` on `[c,1]`). Exact checks: `f₀(1) = f₀′(1) = 0`, `f₁ = Lf₀` has `f₁(1) = f₁′(1) = 0`,
`f₂ = L²f₀` has `f₂(1) = 0`, `f₂′(1) = (1−u)²r(u)` (`r` linear, exact). Three applications:
`|I_B| ≤ (b̂² − 81u)⁻³[|f₂′(1)|(1−u)|y(1)| + ∫|L³w_B||y|]`, with `|y| ≤ min(B_n, (1−up²)⁻⁴)`, `B_n = binom(J+4, 8)` (Szegő 7.32.1,
and `|d^J_{44}| ≤ 1`; both as used in FLAT_GRAVITON (Q3)). Checked numerically against 60-digit quadrature at `J = 60, 201`,
`u = 0.5, 0.9, 1`: agreement `≤ 3·10⁻³⁰`.
* Away from `(p, u) = (1, 1)`: Arb, 2-D centred evaluation of the exact `L³w_B` (exact rational Taylor shift in `p` and `u` to each cell
  centre), `u`-cells geometric in `1−u` down to `1/100`, `p`-cells uniform then geometric to `1 − 10⁻⁹`.
* Corner `u ≥ 0.99`, `p ≥ 0.99`: `L³w_B = p^{−6}G(δ, ε)`, `δ = 1−p`, `ε = 1−u`, `G` an exact polynomial with **no terms of degree < 2**
  (checked). With `q = 1−up² ≥ ε + κδ`, `κ = 1.97`: the quadratic part integrates to `B_n^{1/4}·sup_e G(e)`,
  `G(e) = ∫₀^∞|c₂₀x² + c₁₁xe + c₀₂e²|min(1, (e+κx)⁻⁴)dx` (Arb, closed-form tails); the rest is `≤ K₃(δ+ε)³`, giving
  `K₃(¼ + log(0.03·B_n^{1/4}))/κ`. The boundary term is `≤ sup|r|·B_n^{1/4}`.
* *Uniformity in J.* Every term is checked at `J = 400`: the main term increases with `J`, `C_R/(1−ds)⁴` decreases, and
  `b̂⁴(b̂²−81u)⁻³ × {1, B_n^{1/4}, log B_n}` decrease (`B_n^{1/4} ≤ (J+4)²/8!^{1/4}`), so the inequality holds for all `J ≥ 400`.
* Numbers (`J = 400`, worst cell `u ∈ [0.99, 1]`, pattern `s₂ = +1, s₃ = −1`): main `1.741·10⁷`, `C_R` term `4.69·10⁶`, far part
  `1.35·10⁶`, corner part `8.57·10⁶`; margin `0.15` of the main term. The patterns with `s₂ = s₃` have margin `0.75`: the corner is hard only
  for `s₂ ≠ s₃` because the backward peaks of `B₂` and `B₃` cancel in `𝔉` but not in `N`.

**What made the corner tractable.** NF §4.1 needed `(1−p)^{12}` weights because Lemma C multiplies the weight by `(1−up²)^{−8}`. Lemma B
works with the Jacobi operator itself, whose adjoint `L` maps polynomials to polynomials, so the `(1−up²)^{−8}` never appears; the
`(1−p)⁵` weights suffice. No successor functional is needed: **NF02 itself is certified.**

---

## 2. Q's flat statement

### 2.1 REPRO D9 checked, and two corrections [NE]
* `q_cert.py dense` (1201 log-spaced `u` plus dense windows, `J ≤ min(3000, (400/√u − 1)/2)`): `+−` min `0.196585` at J = 26,
  `u = 0.02535` (`m = 6.28M`) — **REPRO's number reproduced**; `++` J = 2: `0.199895` at `u = 0.788` — reproduced.
* **Lower minimum at large J.** `q_scan_bhat.py` (`u ∈ [10⁻⁶, 10⁻²]`, `J ≤ 2·10⁴`, `b̂ ≤ 200`) and `q_scan_refine.py`: `++` min
  `0.18349` at J = 3608, `u = 1.0·10⁻⁴` (`m = 100M`, `b̂ = 72.2`). Mechanism: at small `u` and `b̂ ~ 70` the three terms
  `γ₀/b̂³` (B₂), `σuψ₃(0)/b̂` (B₃, negative for `++` since `ψ₃(0) = −5.70`) and `ab̂²u²/2` (forward) are comparable, and `N` counts B₃ with
  the opposite sign. REPRO's and Stage 6's dense grids stopped at `J = 3000` / coarser `u`. `+−` in that region: `0.1998`.
* **Rounding defect at `m = M`.** `q_backward.py`: splitting the stored doubles exactly as `ψ = (Q3\* part) + rem` (below), `rem` of
  degree ≤ 3 with coefficients `≤ 1.7·10⁻⁸`, and evaluating both parts accurately:

  | `+−`, `u` | stored (binary) Q3 | Q3\* (`rem` dropped) |
  |---|---|---|
  | 1 (`m = M`) | **`𝔉 < 0` for J ≥ 1990**; `𝔉/N = −0.857` at J = 20000 | min `+0.20000` (J = 10) |
  | 0.9999 … 0.3 | min `+0.20002 … +0.417` | identical |

  The cause is the backward peak `d̃^J_{44}(−1) = binom(J+4, 8)` acting on `w(1) ≈ 10⁻⁹`. Previous dense checks cut `J` at
  `(400/√u − 1)/2 ≈ 200` at `u = 1`, so they never saw it. **Fix (definitional, changes nothing else above `10⁻⁹`):** state Q3's weights
  in factored form, `ψ₂ = p(1−p)³X(p)`, `ψ₃ = (1−p)³Y(p)`, with `X, Y` the exact quotients of the stored doubles (`q_cert.exact_factor`;
  the dropped remainders have coefficients `≤ 1.7·10⁻⁸`; `ψ₂(0) = 0` and `A = 2ψ₂′(0) = 17.240713160205082` are unchanged). This is
  **Q3\*** (exact-coefficient fingerprint `sha256 80e0834a…33f2`), the object certified below. In AdS the same defect would enter only through the backward edge of the smeared pieces at
  `m ≈ Λ`, `J ≳ 2000`; it is not a computed state, and Q3\* removes it.

### 2.2 The certified floor for Q
**`9/50 = 0.18`**, certified on the regions of §2.3, which contain every state where any scan found `FF/N < 0.30`; the global flat
minimum is `0.1835` [NE] (§2.1), so 0.18 is within 2 % of the truth. The nominal δ = 0.2 is **not** met (by 8 %). A constant closer to
0.1835 would need the S cells to be sharper than their ≈ 5 % (FLAT_GRAVITON §3); not attempted.

### 2.3 Region-by-region status of the Q3\* certificate at `δ′ = 9/50` (`q_cert.py`)
`q_cert.py` is ROBUST_FUNCTIONAL's certificate (`rcert.py`) with the one extra term `−2c(1+δ′)u²` of Q's action inserted in each region
(P: polynomial; C: `Φ → Φ − 4c_q`, since `m² ≤ 4`; M: `η → η + c_q/(a_φb̂²/2)`; S: the exact `u²` Taylor coefficient lowered by `c_q`).

| region | code | result at `δ′ = 9/50` |
|---|---|---|
| P: `J < 400`, all `u ∈ [0,1]`, both channels | exact Lemma-P polynomial − `c_q u²`, Descartes | **PASS** (2384 polynomials; includes `+−` J = 26, `u = 0.0253`, ratio 0.1966) |
| S: `J ≥ 400`, `u ≤ 1/4`, `b̂ ≤ 120` | Cauchy cells | **PASS** (2416 cells per channel; contains the global minimum 0.1835, `b̂ ≈ 72`) |
| C`++`: `J ≥ 400`, `u ≥ 1/4` | `|P_J| ≤ 1`, forward term | **PASS** (`2aJJ − c − 4c_q = 44558 > 64∫(|W₀|+|W₁|) ≤ 254`) |
| C`+−`: `J ≥ 400`, `u ≥ 1/4` | `|d̃⁴⁴| ≤ min(binom, (1−p²)⁻⁴)` | **FAIL** (`44552 < 64U = 2.5–3.6·10⁵`): Q's forward coefficient `a = 0.169` is too small for this crude bound. Lemma B (§1.3) would close it; not done |
| M, `1/100 ≤ u ≤ 1/4`, `b̂ ≥ 120` | Lemma M with the forward term | **PASS** (all 8 patterns) |
| M, `u ≤ 1/100`, `b̂ ≥ 120` | Lemma M | **FAIL**: `γ₀b̂ = 3.6–5.2 × 120 ≈ 430–620` against `C_R ≈ 1.0–1.5·10⁴`; Q's `ψ₂″(0)` is small in these units, so Lemma M bites only at `b̂ ≳ 3·10³`, beyond the S cells' validity at `J_P = 400` |

Float evidence in the two uncertified regions (`q_scan_uncov.py`) [NE]: `J ≥ 400`, `u ≤ 1/100`, `b̂ ∈ [120, 600]`, `J ≤ 2·10⁴`: min
`0.304` (`++`, `b̂ = 123`) and `0.998` (`+−`); `+−`, `J ≥ 400`, `u ≥ 1/4` (Q3\*, `J ≤ 6000`): min `0.997`. Both far above 9/50.

---

## 3. The k = 2 weights

### 3.1 `X^Q` (missing from `prl/data/`) — now published [DV]
`prl/data/k2/xq_abel_X_{300,1000}.txt`: the stored min-norm Abel weight of Q's `ψ₂` (`resultq/cache/rv_k2eng_Y_*.pkl`, RESULT_Q_REVISION
§1.3), same columns as `d_abel_X_*.txt`. Independent re-solve (`k2_xq_check.py`, DESHIFT §4.1's solver with Q's `ψ₂`):

| | Δ_g = 300 | Δ_g = 1000 |
|---|---|---|
| `max|X_new − X_stored|/max ψ₂^Q`, ν < 2 / 2–10 / 10–25 / 25–Λ | `2.1·10⁻⁷ / 5.2·10⁻⁸ / 2.4·10⁻¹⁰ / 3.3·10⁻¹²` | `8.3·10⁻⁷ / 6.0·10⁻⁸ / 2.9·10⁻¹⁰ / 1.3·10⁻¹¹` |
| V3 (graviton pole), stored / new | `1.3·10⁻¹⁴ / 2.6·10⁻¹⁵` | `2.1·10⁻¹⁴ / 6.0·10⁻¹⁵` |
| `dk₀^Q`, stored (new agrees to) | `−1.630146574·10⁻⁴` (`3·10⁻¹³`) | `−1.453689613·10⁻⁵` (`8·10⁻¹¹`) |
| off-grid `|E − Φ|` at t = 14 / 18 / 22 (new / stored) | `1.4·10⁻⁸ / 1.2·10⁻⁶ / 6.5·10⁻⁵` / `3.9·10⁻⁸ / 1.6·10⁻⁷ / 1.0·10⁻⁴` | `2.7·10⁻⁸ / 9.3·10⁻⁶ / 7.3·10⁻⁴` / `1.1·10⁻⁷ / 8.5·10⁻⁶ / 8.3·10⁻⁵` |

Like D's, the grid weight holds the identity only to `t ≈ 18–22`; §3.2 removes that limitation.

### 3.2 Continuum k = 2 weights for D and Q [DV]

**Shifted form, k = 2 [D].** In the χ-normalisation (`X = χ(R₋+R₊)`, `prl/data/README` §4) the identity is
`W[χ](t) := ∫₀^∞χ[R₋cos((ν+i)t) + R₊cos((ν−i)t)]dν = Φ(t) := ∫₀^∞ψcos(νt)dν`. Exact facts: `R₊(ν) = R₋(−ν)`, `R₊ = conj R₋` on ℝ
(checked symbolically). For real even χ the integrand is even, so `W = ½∫_ℝ(…)`. Shifting the `R₋` term by `ν = μ − i` and the `R₊` term
by `ν = μ + i` gives, for **every** real `t`,
`W(t) = ½∫_ℝ E(μ)cos(μt)dμ`, `E(μ) = f(μ) + f(−μ) = 2Re f(μ)`, `f(μ) = χ(μ−i)R₋(μ−i)`,
provided `χR₋` is analytic in `−1 < Im ν < 0`, continuous up to the boundary, and decaying. No residue terms arise (unlike C4, the
target is a function, not a pole). With `Φ = ½∫_ℝψ(|μ|)cos(μt)`: `|W(t) − Φ(t)| ≤ ε_W := ∫₀^∞|E − ψ|` for all `t ≥ 0`.

**Pole structure from the common zeros [D + DV].** In the lower strip `R₋` has a single pole, at `ν = 0` (cancelled by `χ(0) = 0`), and a
**double zero at `−i/2`** (`R₋ ∝ (2ν+i)²…`); `R₊` mirrors this at `+i/2`. So χ may carry poles of order ≤ 2 at `±i/2`. The data decide:
with the prefactor `ν²/(ν²+¼)^m`, `m = 0, 1, 2`, the least-squares fit reproduces the stored weight only for `m = 2`
(`max|χ − χ_stored|/max χ = 1.3, 0.09–1.3, 5·10⁻⁶`). **Kink.** `ψ(|μ|)` has `ψ′(0) = κ ≠ 0` (`κ = 0.6328 / 0.1892` for D,
`0.02888 / 0.008633` for Q); no χ analytic on the closed strip can produce it. The exact solution is singular at `ν = ±i`; the
functions `ν^{2j}(ν²+1)log(ν²+1)e^{−ν²/2}` (analytic in `|Im ν| < 1`) have boundary values `−2π b₀|μ| + O(μ²log|μ|)` in `E`, and carry
the kink exactly.

**Construction.** `χ = Q(ν)·B(ν) + Σ_{j≤2}y_jν^{2j}l(ν)e^{−ν²/2}`, `Q = ν²/(ν²+¼)²`, `B` = even Gaussian RBFs (spacing/width
`0.15/0.45` to ν = 40, `0.3/0.9` to `Λ+30`; 1233 / 3567 unknowns at 300 / 1000). Least squares (blocked Householder QR, double
precision) on `E − ψ` (≈ 5·10³–1.3·10⁴ μ-rows to `Λ+60`) plus weak rows pulling χ to the stored weight's real part (and to 0 beyond the
stored grid). The imaginary part of the stored χ (8 % of `max Re χ` at small ν; a null vector, `W[Im χ] = 0`) is fitted the same way
with target 0. `ε_W` is then a 4.2·10⁵-point quadrature of `|E − ψ|` to `Λ + 60` (Gaussian tail beyond `< 10⁻¹⁰⁰`).

| weight | `ε_W` (Re) | `/∫|ψ|` | `ε_W` (Im, null) | `max|X − X_stored|/max ψ` (ν < 2 / 25–Λ) | `dk₀` continuum (stored) | V3 rel. |
|---|---|---|---|---|---|---|
| D, 300 | `3.3·10⁻⁶` | `1.1·10⁻⁹` | `4.1·10⁻⁷` | `3.8·10⁻⁷ / 1.7·10⁻¹⁰` | `−7.666324482·10⁻⁴` (`−7.666324168·10⁻⁴`) | `6·10⁻¹¹` |
| D, 1000 | `1.1·10⁻⁶` | `1.1·10⁻¹⁰` | `1.3·10⁻⁷` | `4.3·10⁻⁷ / 5.6·10⁻¹¹` | `−6.830149984·10⁻⁵` (`−6.830146911·10⁻⁵`) | `8·10⁻¹⁰` |
| Q, 300 | `1.9·10⁻⁷` | `6.2·10⁻¹⁰` | `1.8·10⁻⁸` | `2.0·10⁻⁷ / 1.4·10⁻¹⁰` | `−1.630146593·10⁻⁴` (`−1.630146574·10⁻⁴`) | `8·10⁻⁹` |
| Q, 1000 | `7.4·10⁻⁸` | `7.3·10⁻¹¹` | `5.6·10⁻⁹` | `2.2·10⁻⁷ / 2.1·10⁻¹¹` | `−1.453689637·10⁻⁵` (`−1.453689613·10⁻⁵`) | `3.5·10⁻⁷` |

(V3: the graviton-pole identity `∫χ[R₋D(ν+i) + R₊D(ν−i)] = ∫ψD`, `D = 1/(z²+9/4)`; its error is `≤ (4/9)ε_W` by construction, observed
smaller. Double vs 50-digit evaluation of χ: `≤ 1.1·10⁻¹²` of `max|χ|`.)

**Consequences.**
* The continuum weights are the stored functionals to the grid's own discretisation accuracy (`≤ 4·10⁻⁷` of `max ψ` vs the
  `1.5·10⁻⁶` spread between DESHIFT's two grid solvers), now with `|W − Φ| ≤ ε_W` **for all t** (the grids fail at `t ≈ 22`). Hence
  `|image − target| ≤ ε_W P_{−½}(cosh β)` for all β, as for C4/C5.
* **Cubic bound unchanged:** `R₀` depends on the weight only through `K = k(0)(1 + δk(0))`; `(1+dk₀)^{−½} = 1.00038354 / 1.00003415`
  for both weights; `R₀(K = 1) = 3.13206·10⁻⁶ / 2.89126·10⁻⁸`, `R₀(K = 1+dk₀) = 3.133256·10⁻⁶ / 2.891361·10⁻⁸` (either weight).
* **DS for D** can be upgraded from "[D + NE] (exact min-norm X exists)" to "ε-DS for the explicit continuum weight [DV]", exactly as
  DESHIFT did for Q; exact existence is then not needed. The same holds for Q's smeared k = 2 piece.
* **Cutoff.** Use D's cutoff (`ν_c = Λ+100`, `a = 10`) for both k = 2 weights: `X^Q` extends to `ν ≈ Λ`, so Q's `ν_c = 150` (meant for
  C4/C5) would truncate it. UNIFIED_THEOREM_v2 §2.1 should say so.

---

## 4. Wording for the letter

* **D, flat.** "In the bulk-point limit `Ω^EE` reduces to the forward-free flat functional NF02, for which we prove by a computer-assisted
  (interval-arithmetic) certificate that `𝔉 ≥ N/5` on every heavy state; the bound is attained to `10⁻⁸` at `m = M`. **[certified]**"
  Delete "except one corner" and the "float-checked 0.765" caveat.
* **Q, flat.** "The flat limit of `Ω^Q_up` is the δ = 0.2 functional Q3; its flat margin `FF/N` is ≥ 0.18 on every state we
  computed (minimum 0.1835 at `m ≈ 100M`, `J ≈ 3600` [NE]) and is certified `≥ 9/50` on all states with `J < 400` or `b̂ ≤ 120` and on most of the rest [certified]; the uncertified remainder (`+−`, `m ≤ 2M`, `J ≥ 400`; `m ≥ 10M`, `b̂ ≥ 120`) is float-checked [NE]." Replace "δ = 0.2 (dense min 0.1966;
  nominal at the 2 % level)" by "nominal δ = 0.2; flat minimum 0.1835 [NE]". **State Q3's weights in factored form, `ψ₂ = p(1−p)³X(p)`, `ψ₃ = (1−p)³Y(p)` (Q3\*, §2.1).**
* **Weights.** "All weights defining the functionals are published, including the k = 2 de-shift weights as explicit continuum functions
  satisfying the de-shift identity to `ε_W/∫|ψ| ≤ 1.1·10⁻⁹` uniformly in `t` [DV]."

## 5. Suggested changes to UNIFIED_THEOREM_v2 (not made here)
1. §5 "Not certified" bullet and §11 item 4: replace by "NF02: `𝔉 ≥ N/5` on all flat states [certified] (CERT_CLOSE §1)". Table row
   "Flat" in §5.2: certificate grade [certified].
2. §2.3 / §12 "Q3 δ nominal: dense min 0.1966": → "flat min 0.1835 (`++`, `m ≈ 100M`) [NE]; `FF ≥ 0.18N` certified on `J < 400` and on
   `b̂ ≤ 120` (CERT_CLOSE §2.3)"; add the factored-weight definitional fix (§2.1 here) to §2.3's definition of Q3.
3. §2.2, §11 item 18, row DS: D now has an explicit continuum k = 2 weight with ε-DS [DV] (`prl/data/k2/`); §11 item 17 is closed
   (`X^Q` published, re-solved).
4. §2.1: the k = 2 pieces of Q use D's cutoff parameters.

## 6. Reproduce
    cd ads4/code/cert; bash cert_run_all.sh nf   # NF02 at 1/5 (~45 min; NF_PEXT=1 adds the 30-min P-extension cross-check)
    bash cert_run_all.sh q                        # Q scans + Q3* regions at 9/50 (~45 min)
    bash cert_run_all.sh k2                       # k = 2 continuum weights, X^Q re-solve, export + checksum check (~15 min)
