# CONTOUR_BOUNDS: certified contour bounds for the signed estimate of ℰ₂ (ANGULAR_ROUTE §§36–38, §40)

*26 Nov 2026. Research agent, about 1.5 h. Unreviewed and not committed. No existing file was edited; nothing under `prl/` or `code/fullpos/` was touched.*
* *Code:* `ads4/code/contour/`, with outputs in `cache/`. The cache's `.gitignore` ignores everything but itself.
* *Read-only inputs:* `prl/data/k2/k2_continuum_D_1000_rep.json`, `prl/data/kernels.json`, `prl/data/nf02_psi.json`, and `code/epswcert/ew_common.py` (imported; its cache path is redirected).
* *Compute:* about 15 CPU-min, ≤ 2 processes, `ulimit -v 3000000`.

Grades:
* **[certified]**: Arb ball arithmetic, with every remainder bounded.
* **[conditional]**: certified arithmetic, applied to an author input that has not been re-checked.
* **[DV]**: sampled.

---

## 0. Verdict

**New large-impact threshold for the leading connection term ℰ₂.** The factors 100/a², 300/a² and |D| ≤ 5·10¹² are kept exactly as in §38, with a ≥ 998.

| target | §38 as published (its own numbers, re-evaluated) | **certified numbers, §38 factorisation** | **certified numbers, \|fB\| bounded jointly** |
|---|---|---|---|
| \|ℰ₂\| < 0.1·𝓜(L) | L ≥ 37.4 | L ≥ 29.26 | **L ≥ 29.12** |
| \|ℰ₂\| < 0.01·𝓜(L) | L ≥ 38.5 | L ≥ 30.30 | **L ≥ 30.16** |
| \|ℰ₂\| < 0.001·𝓜(L) | L ≥ 39.8 | L ≥ 31.38 | L ≥ 31.24 |
| \|ℰ₂\| < 1.01·10⁻⁴·𝓜(L) (§38's "64") | L ≥ 48.4 | L ≥ 34.67 | L ≥ 34.53 |

How these thresholds were obtained:
* All thresholds use the horizontal leg at Y = 9/4, which is the best of the six heights tried (§3).
* Each is rigorous on a cell cover of [L₀, 400] (§6). Beyond L = 400 every term decays.
* The §38-as-published column uses W ≥ 0.03 on [0, 0.01], |B| ≤ 20 and the 10⁷ / 2·10⁹e^{17(y−3/4)} / 4·10³¹ constants. Its "L ≥ 64" was already conservative: its own numbers give 1.01·10⁻⁴ at L ≥ 48.4.

**Which piece dominates.** The horizontal leg is the whole story.
* At the threshold, the horizontal leg contributes 98–99 % of the bound ρ(L). The other contributions are:
  * short detour: ~10⁻¹²;
  * vertical detour: ~10⁻¹²;
  * the floor 100/a² = 1.004·10⁻⁴.
* The bottleneck is now the transform bound |D| ≤ 5·10¹², not f:
  * ∫|fB|(t+9i/4) dt is concentrated (> 99.99 %) at t < 45;
  * there, §38's own profile 200(1+t)^{5/4} is ≤ 2.4·10⁴, not 5·10¹².

**With §38's D-profile [conditional].** This uses |D(t+9i/4)| ≤ 200(1+t)^{5/4} for t ≤ 998 (below the turning point, g = 0) and 5·10¹² beyond.
* ρ < 0.1 for L ≥ **19.22**.
* ρ < 0.01 for L ≥ **20.27**.
* ρ < 0.001 for L ≥ 21.36.
* At L = 28, ρ ≤ 1.0040·10⁻⁴, which is the floor.

**Comparison with §40's L ≥ 28.** The authors claim a k = 2 budget < 0.0035·𝓜 at L = 28.
* *With the uniform 5·10¹²:* my certified ℰ₂ bound at L = 28 is ρ ≤ 1.19, which does **not** support the claim. ρ(30) ≤ 0.014.
* *With the D-profile:* ρ(28) ≤ 1.004·10⁻⁴, which supports it with a wide margin.
* *With the authors' own ∫|fBD| < 10²¹:* ρ(28) ≤ 5.5·10⁻⁴, which also supports it.
* **So §40's L ≥ 28 rests on the D-profile, not on the band sup 5·10¹².** The profile is §38's own claim, marked "not re-checked", and is the one input to verify.
* This covers the ℰ₂ part only. The connection terms j ≥ 1 and k = 3 of §40 are not checked here.

### Certified numbers versus the values §38 used

| # | quantity | §38 | **certified** |
|---|---|---|---|
| 1 | sup\|f(t+3i/4)\|, 0 ≤ t ≤ ¼ | 10⁷ | **≤ 3.924** (sampled max 3.725, at t = 0) |
| 2 | sup over 3/4 ≤ y ≤ 9/4 of \|f(¼+iy)\|e^{−17(y−3/4)} | 2·10⁹ | **≤ 227** (the profile itself is the table in §2) |
| 3 | ∫_{1/4}^∞ \|f(t+9i/4)\| dt | 3·10¹⁷ | **≤ 9.0084·10¹³** |
| 3′ | ∫_{1/4}^∞ \|fB\|(t+9i/4) dt | (3·10¹⁷·20) | **≤ 5.3018·10¹¹** |
| 4 | W ≥ 0 on [0, ¾] | claimed | **certified.** W = (½−y)²V(y) with V > 0, a double zero at y = ½ |
| 4a | min W on [0, 0.01] | > 0.03 | **≥ 0.05865** (W(0) = 0.0604088…) |
| 4b | max W on [0, ⅛] | ≤ 0.25 | **≤ 0.06049** |
| 4c | max W on [⅛, ¾] | ≤ 64 | **≤ 0.5893** (at y = ¾) |
| 5 | \|B\| on the detour (short and vertical pieces) | ≤ 20 | **≤ 0.1034** |
| 5′ | \|B\| on the horizontal leg Y = 9/4 | ≤ 20 | **≤ 0.0171** (≤ 0.0080 for t ≤ 1054) |

Notes:
* *Row 2:* the sup is attained at y = 9/4: 2.69·10¹³·e^{−25.5} ≈ 227. The growth of |f(¼+iy)| is closer to e^{24y} than to e^{17y}, so integrating the table of §2 directly beats any e^{17y} envelope.
* *Row 3:* the integral comes from t < 45 (the σ = 0.45 family); the σ = 0.9 family gives only 1.1·10⁴.
* *The lower bound on 𝓜:* I_W(L) = ∫₀^{3/4}W e^{−Ly} dy ≥ 0.0545/L at L = 28 (0.0523/L at L = 20). §38's 0.03·(1−e^{−0.01L})/L is 0.0073/L at L = 28.

---

## 1. Definitions as implemented (`cb_common.py`)

**The weight.**
* χ = χ_re of the repaired JSON (fields `x_re`, `y_re`).
* This is the realised, real-analytic weight, the one used by `epswcert` (E = 2Re[χR₋](μ−i)) and by `fullpos`.
* χ_im (a null vector) is not included. f is real on the real axis and reproduces ψ there: f(5) = 0.93311395 against ψ(5) = 0.93311396.

**The effective density.**
* θ(z) = ½[erf((z+ν_c)/a) − erf((z−ν_c)/a)], with ν_c = Λ+100 and a = 10.
* G(ν) = θ(ν′)χ(ν′)R₋(ν′), with ν′ = ν − i.
* **f(ν) = G(ν) + conj G(conj ν).** This equals χ_θR₋(ν−i) + χ_θR₊(ν+i), because R₊(z) = conj R₋(conj z).

**Branch.**
* Log(ν′²+1) is replaced by Log ν + Log(ν−2i). On Re ν > 0 the two agree; the right-hand side is continuous up to Re ν = 0.
* This gives the boundary values from Re ν > 0 required by §38, and it is safe for ball evaluation of boxes that straddle Re ν = 0.
* Checked: −2B(iy)Im f(iy) from this evaluation equals the closed form of §4 to 14 digits at 7 points.

**Singularities.**
* f has poles at ±3i/2, ±5i/2, … (from QR₋ and R₋); ν = ±i is removable; there are branch points at 0 and ±2i.
* B(ν) = Γ(2iν)/Γ(½+iν)² has poles at ν = im (m = 0, 1, 2, …) and zeros at i/2, 3i/2, ….
* Boxes are evaluated through the duplication formula with a shift by 5. Arb's Γ returns nan on wide boxes near its poles.

**Method** (`line_panel`):
* Taylor polynomial of order N = 32–40 at each panel centre, computed by power-series ball arithmetic. The large cancellations among the RBF terms (e^{(Y+1)²/2σ²} ≈ e^{26} at Y = 9/4, σ = 0.45) therefore happen inside the coefficients.
* The truncation is bounded by Cauchy on a disc of radius ρ:
  * the disc sup is |θQR₋| (ring of boxes) × Σ_near|x_k| sup|g_k| (term by term), plus |θYlR₋| (ring of boxes);
  * RBF terms with |Re ν′ − c_k| > 20σ_k are bounded on the disc.
* Each panel is cut into subcells by exact polynomial shifts. On a subcell, sup|p| ≤ Σ|c_n|h^n and ∫|p| ≤ Σ|c_n|2h^{n+1}/(n+1).
* f, B and fB are each enclosed this way.
* The per-point truncation and far errors are ≤ 2·10⁻¹⁸ for f and ≤ 10⁻²⁹ for B.
* `check_disc` asserts that every disc avoids the poles and cuts. The distance is ≥ ¼ everywhere on the contour, and ρ = 0.2 near the axis.

## 2. Vertical piece ν = ¼ + iy, ¾ ≤ y ≤ 9/4 (`cb_lines.py vert`) [certified]

Setup: 120 panels of width 1/80, aligned with the ⅛-grid. Each panel has 200 subcells, giving 24 000 enclosures, stored in `cache/cb_vert.json`. Every row of that file gives (y₀, y₁, sup|f|, sup|B|, sup|fB|, ∫|f|, ∫|fB|), so the authors can integrate against e^{−Ly} exactly.

Per ⅛-interval, compared with the authors' §40 table:

| y-interval | **sup\|f\| certified ≤** | sampled max ≥ | authors (§40) | agreement | sup\|B\| ≤ | sup\|fB\| ≤ |
|---|---|---|---|---|---|---|
| [0.750, 0.875] | 12.610 | 12.610 | 13 | ✓ | 0.103 | 1.303 |
| [0.875, 1.000] | 84.020 | 84.020 | 85 | ✓ | 0.103 | 7.937 |
| [1.000, 1.125] | 565.42 | 565.41 | 570 | ✓ | 0.0945 | 37.63 |
| [1.125, 1.250] | 4345.6 | 4345.5 | 4.4e3 | ✓ | 0.0665 | 174.7 |
| [1.250, 1.375] | 4.3010e4 | 4.3009e4 | 4.4e4 | ✓ | 0.0402 | 1047 |
| [1.375, 1.500] | 5.5841e5 | 5.5840e5 | 5.7e5 | ✓ | 0.0244 | 9668 |
| [1.500, 1.625] | 9.4741e6 | 9.4740e6 | 9.6e6 | ✓ | 0.0173 | 1.523e5 |
| [1.625, 1.750] | 1.8410e8 | 1.8410e8 | 1.85e8 | ✓ | 0.0175 | 3.222e6 |
| [1.750, 1.875] | 3.6432e9 | 3.6432e9 | 3.65e9 | ✓ | 0.0191 | 6.954e7 |
| [1.875, 2.000] | 7.1322e10 | 7.1322e10 | 7.14e10 | ✓ | 0.0192 | 1.271e9 |
| [2.000, 2.125] | 1.3873e12 | 1.3873e12 | 1.39e12 | ✓ | 0.0178 | 1.774e10 |
| [2.125, 2.250] | 2.6920e13 | 2.6920e13 | 2.7e13 | ✓ | 0.0128 | 2.112e11 |

**All twelve of the authors' numbers are confirmed as upper bounds.** They are tight, 0.1–3 % above the true sup.
* The sup sits at the right end of each interval: |f(¼+iy)| increases monotonically, by a factor of about 10–20 per ⅛.
* There is no visible feature at the pole y = 3/2. Its residue is small next to the Gaussian growth.

Coarse profile, from the same file: sup|f| over the cell containing y.

| y | 0.75 | 0.85 | 1.0 | 1.1 | 1.25 | 1.4 | 1.5 | 1.6 | 1.75 | 1.9 | 2.0 | 2.1 | 2.25 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| sup\|f\| | 2.0 | 9.2 | 87 | 4.0e2 | 4.7e3 | 7.3e4 | 5.8e5 | 5.5e6 | 2.0e8 | 6.9e9 | 7.5e10 | 8.0e11 | 2.7e13 |

## 3. Horizontal legs ν = t + iY, t ≥ ¼ (`cb_lines.py horiz Y`) [certified]

Panels:

| range | panel width | ρ | subcells per panel | order |
|---|---|---|---|---|
| [¼, 45] | 0.04 | 0.2 | 4 | 40 |
| [45, 1054] | 0.2 | 0.5 | 2 | 40 |

That is 14 568 cells per height.

**Tail t > T_END = c_max + 25 = 1054.** It is bounded analytically:
* |θ(x+iy)| ≤ e^{y²/a²}, from θ(z) = (a√π)⁻¹∫_{−ν_c}^{ν_c}e^{−((z+s)/a)²}ds;
* |QR₋| ≤ its value at |ν′| = 1054, since each factor ratio decreases;
* the Gaussian tails are summed exactly with erfc, over every centre;
* the Y-term is ≤ e^{−T²/4};
* |B| ≤ 2^{−2Y−1}π^{−1}·B(3−Y, ½)·(shift ratio), from Γ(z)/Γ(z+½) = (shift)·B(z+3, ½)/√π with Re(z+3) = 3−Y > 0.
* The tail is ≤ 4·10⁻¹⁷⁸ in every case.

| Y | ∫\|f\| ≤ | from t < 45 | from t > 45 | **∫\|fB\| ≤** | from t < 45 | from t > 45 | sup\|B\| (t ≤ 1054 / tail) | sup\|f\| |
|---|---|---|---|---|---|---|---|---|
| 1.00 | 1.0799e4 | 604 | 1.020e4 | 65.80 | 20.5 | 45.3 | 0.0947 / 0.054 | 110 |
| 1.25 | 8.7959e4 | 7.776e4 | 1.020e4 | 2039.9 | 2008 | 32.0 | 0.0410 / 0.041 | 1.63e4 |
| 1.50 | 1.3699e7 | 1.369e7 | 1.020e4 | 2.3791e5 | 2.379e5 | 22.7 | 0.0252 / 0.032 | 2.73e6 |
| 1.75 | 2.5016e9 | 2.502e9 | 1.020e4 | 3.0205e7 | 3.021e7 | 16.0 | 0.0176 / 0.025 | 4.79e8 |
| 2.00 | 4.7013e11 | 4.701e11 | 1.024e4 | 4.0213e9 | 4.021e9 | 11.4 | 0.0179 / 0.020 | 8.81e10 |
| **2.25** | **9.0084e13** | 9.008e13 | 1.132e4 | **5.3018e11** | 5.302e11 | 9.4 | 0.0080 / 0.017 | 3.01e13 |

**Trade-off.**
* The horizontal cost 5·10¹²·∫|fB|·e^{−YL} needs L ≈ ln(5·10¹²∫|fB|/(0.0027·target))/Y.
* This decreases with Y, because ∫|fB| grows by about e^{4.8} per ¼ in Y while e^{−YL} gains e^{7.5} at L ≈ 30.
* So **Y = 9/4 is the best of the six heights.** The ρ < 0.01 thresholds are:

| Y | 1 | 1.25 | 1.5 | 1.75 | 2 | 2.25 |
|---|---|---|---|---|---|---|
| L ≥ (joint \|fB\|) | 45.4 | 39.0 | 35.6 | 33.3 | 31.5 | 30.2 |

* The 5·10¹² bound is §38's bound for Y = 9/4. For Y < 9/4 it is assumed unchanged, which does not matter since those heights lose anyway.
* Moving the leg closer to the pole and zero pair at 5i/2, for example Y ≈ 2.4, would plausibly gain about 0.5 in L. This was not tried.

**The authors' horizontal number.**
* With the uniform |D| ≤ 5·10¹², the certified bound on ∫|fBD|(t+9i/4) is 5·10¹²·5.30·10¹¹ = **2.65·10²⁴, which is above the authors' < 10²¹.**
* With §38's profile |D| ≤ 200(1+t)^{5/4} (t ≤ 998), it is **≤ 8.14·10¹⁴ [conditional], consistent with < 10²¹.**
* So the authors' 10²¹ must use the profile, not the band sup. Both numbers rest on the profile claim in §38, which has not been re-checked.
* The two families separate cleanly. The σ = 0.45 family (t ≲ 45) carries the whole of ∫|fB| = 5.30·10¹¹. The σ = 0.9 family carries 9.4.

## 4. The local density W on [0, ¾] (`cb_w.py`) [certified]

**Closed form for Im f(iy).** On ν = iy, the Gaussian part and QR₋ are real. Only the Log in χ(ν+i) contributes an imaginary part: the point lies on the cut, and the boundary value is taken from Re ν > 0. So

  Im f(iy) = −π y(2+y)·θ(is)·Y(is)·R₊(is),  with s = 1+y and Y(is) = (y₀ − y₁s² + y₂s⁴)e^{s²/2}.

**Factorised W.** Using B(iy) = Γ(−2y)/Γ(½−y)² and the duplication formula:

  W(y) = (½−y)²·V(y),
  V(y) = 4√π(2+y)·2^{−2y}·[Γ(1−y)/Γ(3/2−y)]·θ(is)·Y(is)·(2s−11)(2s−5)(2s−1)²(2s+1)/[s(2s−7)(2s+3)(2s+7)].

**Certificate.**
* The pole of B at y = 0, B's zero at y = ½ and R₊'s zero at y = ½ cancel exactly. This is the "sign reversals cancel" of §35.
* V is enclosed on 7500 cells. **V > 0 on every cell**, so W ≥ 0 on [0, ¾], with a double zero at y = ½.
* The closed form matches −2B(iy)Im f(iy), evaluated directly, at y = 0.05, 0.2, 0.37, 0.49, 0.51, 0.6 and 0.74 (agreement to ≥ 12 digits).
* Only y_re enters W: the RBF part of χ drops out of the segment entirely.

**Values.**
* W(0) = −2Y(i)R₊(i) = 0.0604088150…
* W(0.1) = 0.0445, W(0.2) = 0.0321, W(0.3) = 0.0217, W(0.4) = 0.0095, W(0.5) = 0.
* W(0.6) = 0.0340, W(0.7) = 0.265, W(0.75) = 0.5888.

**Bounds.** These are the four certified W bounds in the table of §0:
* min W on [0, 0.01] ≥ 0.058658;
* max W on [0, ⅛] ≤ 0.060488;
* max W on [⅛, ¾] ≤ 0.589212.

**I_W(L) lower bounds.** The cell lower bounds give I_W(L) = ∫₀^{3/4} W e^{−Ly} dy:

| L | 10 | 20 | 28 | 40 | 64 |
|---|---|---|---|---|---|
| I_W ≥ | 4.606e-3 | 2.618e-3 | 1.945e-3 | 1.403e-3 | 9.008e-4 |

## 5. |B| on the contour pieces [certified]

| piece | sup\|B\| |
|---|---|
| short, t + 3i/4 | 0.1006 |
| vertical, ¼ + iy | 0.1034 (max near y = 0.89) |
| horizontal, Y = 9/4 | 0.0080 for t ≤ 1054; tail bound 0.0171 |

§38 used 20.

## 6. The final inequality (`cb_final.py`)

**The bound.** From §38, with the A ≥ 0 part of 𝓜 dropped and 𝓜 ≥ √w(2/Λ)I_W(L):

  |ℰ₂|/𝓜 ≤ ρ(L) := 100/a² + 2X(L)/I_W(L),

  X(L) = (300/a²)[e^{−3L/4}∫₀^{1/4}|fB|(t+3i/4)dt + ∫_{3/4}^{Y}|fB|(¼+iy)e^{−Ly}dy] + 5·10¹²·e^{−YL}∫_{1/4}^∞|fB|(t+iY)dt.

Here a = 998, and the vertical integral uses the cell sups against the exact ∫e^{−Ly}. Two ways of bounding fB:
* *"§38 factorisation":* sup|f|·sup|B| piece by piece.
* *"joint":* |fB| enclosed directly.

**Rigour in L.** X and I_W both decrease in L. So on a cell L_k ≤ L ≤ L_{k+1}, ρ(L) ≤ 100/a² + 2X(L_k)/I_W(L_{k+1}).
* The thresholds come from a cover of [L₀, 400]: steps of 0.5, then 0.01 near L₀.
* Beyond L = 400, X decays exponentially while I_W ≥ 0.058(1−e^{−0.004L})/L.

**Results.**
* The threshold table is in §0.
* Values of ρ(L) with the joint bound, Y = 9/4:

| L | 20 | 25 | 28 | 30 | 40 |
|---|---|---|---|---|---|
| ρ(L) ≤ | 5.8e7 | 917 | 1.19 | 0.0141 | 1.004e-4 (floor) |

* With the D-profile [conditional]: ρ(20) ≤ 0.0179, ρ(25) ≤ 1.007e-4.

## 7. Caveats

1. **The weight is χ = χ_re.** If §38's χ were meant to include iχ_im, every number in §§1–3 would change, but W would not (§4 uses only y; y_im is 10⁻² of y_re).
2. **Author inputs taken as given.** These are 100/a², 300/a², |D| ≤ 5·10¹² (and its profile, in the conditional variant), a ≥ 998 and A ≥ 0. For Y < 9/4 the 5·10¹² is an assumption.
3. **Scope.** Only the leading connection term ℰ₂ is covered. The j ≥ 1 connection terms and k = 3 of §40 are not examined.
4. **Reproduce:**

   ```
   cd code/contour
   python3 cb_w.py
   python3 cb_lines.py short
   python3 cb_lines.py vert
   for Y in 1.0 1.25 1.5 1.75 2.0 2.25; do python3 cb_lines.py horiz $Y; done
   python3 cb_final.py
   ```

   Run times: about 2 min per height, under 1 min for everything else.
