# Functional data for the Supplement (ads4/prl/data)

Written by `ads4/code/deshift/export_data.py` (23 Nov 2026) from the stored, fixed functionals; **updated 23 Nov 2026 (v4)**: errata
for stale values (§5), the Q3\* file (§1), pointers to `k2/` and to the figure script (§6). Nothing here is new physics: it is the
data that defines `Ω^EE` (Result D) and `Ω^Q_up` (Result Q), in plain text / JSON, plus the continuum (cut-off) de-shift weights of
`DESHIFT_CLOSE.md`. `MANIFEST.txt` lists sizes and SHA-256 prefixes. Existing data files were **not** modified in v4 (their hashes are
quoted in the Supplement); where a field inside a file is out of date, §5 says so.

Conventions: `Λ = √((Δ_g−2)(Δ_g−1))` (298.4996 / 998.4999); `p = ν/Λ`; `P(μ) := P_{−1/2+iμ}(cosh β)`; heavy images are leading-Regge
(LR); `C[f](t) = ∫₀^∞ f cos(νt) dν`, `S[f](t) = ∫₀^∞ f sin(νt) dν`.

## 1. Flat weights

| file | content |
|---|---|
| `nf02_psi.json` | NF02 (Result D): `ψ_k(p) = (1−p)⁵ p X_k(p)/k₀`, `k = 2, 3`, as exact rational coefficient lists `psi2_exact`, `psi3_exact` (lowest power first; divide by `k₀ = 27.48509817`). `Ω^EE = ∫₀¹dp ψ₂ Λ⁴Ĉ′_{EE,Λp} + ∫₀¹dp ψ₃ ½Λ⁶Ĉ_{3,Λp}`; normalisation `k(0) = 1`. |
| `q3_full_delta0.2.json` | `Q3-full-delta0.2` (Result Q, δ = 0.2): `psi2`, `psi3` (doubles as exact `repr` strings), `a`, `c` (called `c_Q` in the Supplement), `A = 2ψ₂′(0) = 17.2407`, `B`, `κ`, and the flat-action definition. These stored doubles carry the factor `(1−p)³` only to rounding; the **definition** of the flat target is Q3\* (next row). |
| `q3star_factored.json` | **Q3\*** (new in v4; `CERT_CLOSE.md` §2.1): `ψ₂ = p(1−p)³X(p)`, `ψ₃ = (1−p)³Y(p)` with `X_exact`, `Y_exact` the exact quotients (rational strings, lowest power first) of the stored doubles, read as exact binary rationals; also the expanded exact `ψ₂`, `ψ₃`, float copies, the dropped remainders (`≤ 9.4·10⁻⁹`), the largest coefficient change (`1.7·10⁻⁸`) and the canonical fingerprint `sha256 80e0834a…33f2` that `ads4/code/cert/q_cert.py` certified. `A = 2ψ₂′(0)` changes by `4·10⁻¹¹` relative. Produced by `make_q3star.py` (standard library only; asserts the fingerprint). |

## 2. Leading-Regge kernels (`kernels.json`)
* **Q, k = 4 (EE3):** `R_j`, `j = ±1, ±3`. Image of a weight χ: `Σ_j ∫₀^∞ χ R_j P(ν − ij) dν`; de-shift transform
  `W[χ](t) = ∫₀^∞ χ Σ_j R_j cos((ν − ij)t) dν`. Value piece: `W = 1` (image `P_{−1/2}`); derivative piece: `W = −t²/2` (image `−D(β)`).
* **Q, k = 5 (C5, `++`; `+−` = −`++`):** `R_j`, `j = 0, ±2`; `N̂₅ = Σ_j R_j`. Functional `∫₀^∞ Re[c Ψ_{5,ν}/N̂₅] dν`; `W[c](t) =
  Re∫₀^∞ (c/N̂₅) Σ_j R_j cos((ν − ij)t) dν = 1`.
* **D, k = 2 (EE):** `R_∓` and the `q²` kernels that fix the phase convention of the χ-normalisation (§4).

## 3. Result Q de-shift weights
**As fixed on 22 Nov (Stage 7, `functionals_q7_xprec.json`):** `q_deshift_xp_C4_value.txt`, `q_deshift_xp_C4_deriv.txt` (χ_k, 64
digits) and `q_deshift_xp_C5.txt` (`c_k`), on `ν_k = (k + ½)/25`; the functional is the midpoint sum `(1/25)Σ_k …`. The C4 weights
satisfy `W = Φ` to `4·10⁻¹⁵` for `t ≤ 30` and fail at `t ≈ 41`; the C5 weight only to `t ≈ 18` (`DESHIFT_CLOSE.md` §2).

**Continuum (cut-off) weights (`DESHIFT_CLOSE.md` §2; recommended definition):** `q_deshift_continuum_{c4v,c4d,c5}.json`.
* **Form.** C4: `χ(ν) = Q(ν)B(ν)`, with `Q = ν²(ν²+1)(ν²+4)/[(ν²+9)^p(4ν²+1)(4ν²+25)(4ν²+(3+2√2)²)]`, `p = 1` (value) / `3`
  (derivative). C5: `D := c/N̂₅ = Q B`, with `Q = iν(ν²+1)/(ν²+4)·Π_a i/(ν − ia)`, `a ∈ {a₅, 3/2, √73 − 15/2}`, where `a₅ = −1.19592…`
  is the root of `32ν⁵+1072iν⁴−10000ν³−29240iν²+20922ν−1809i` on the negative imaginary axis (125 digits in the file).
* **Basis.** `B(ν) = x₀g(ν) + Σ_{k=1}^{K}[x_k^{re}(g(ν−kh)+g(ν+kh)) + x_k^{im} i(g(ν−kh)−g(ν+kh))]`, `g(z) = exp(−z²/2σ²)`;
  `h = 0.15`, `σ = 0.45`, `K = 400` (C4 has only `x^{re}`). The coefficients are given to 60 digits; **evaluate in ≥ 300-bit arithmetic**
  (the basis is ill-conditioned by design).
* **ε_W.** For all `t ≥ 0`, `|W − Φ| ≤ ε_W`, hence `|image − target| ≤ ε_W P_{−1/2}(cosh β)` for all β. **Use the certified values
  `ε_W ≤ 2.673·10⁻¹³` (C4 value) and `≤ 9.058·10⁻¹¹` (C4 derivative)** (ball arithmetic, `HK_RESIDUAL.md` §3.3), and `2.4·10⁻⁷` (C5, a
  384-bit quadrature, not certified). The `certificate` fields inside the C4 JSON files still print the earlier midpoint sums
  `1.2·10⁻¹³` / `4.2·10⁻¹¹`; those are estimates, not bounds (§5).
* **Relation to the stored weights.**
  * C4 value and derivative agree with the XP weights to `1.9·10⁻¹⁷` and `1.1·10⁻¹⁶` of the maximum: the same functional.
  * C5 is a different gauge (max relative difference 0.37, amplification 1.91 vs 1.80). The LR image is the same.
* `*_table.txt` tabulates the weights at `ν = k/50`, `0 ≤ ν ≤ 60`, for convenience only (beyond 60 the weights are `< 10⁻²⁵` of their
  maximum).
* **Cutoff.** HK_CLOSE's erf cutoff (`ν_c = 150`, `a = 6`) changes these weights by less than `e⁻³⁰⁰` relative. They already decay
  super-Gaussianly, so they are admissible cut-off weights in HK_CLOSE's sense.

## 4. Result D Abel weights (`d_abel_X_300.txt`, `d_abel_X_1000.txt`)
* **Content.** `X = χ(h₋+h₊)`, `Y = χn`, and `ψ₂(ν_k)`, with `ν_k = (k+½)dν`:
  * `dν = 0.02`, `ν_max = 327.68` (300);
  * `dν = 0.04`, `ν_max = 1310.72` (1000).
* **Equation solved.** `C[aY] − i tanh t S[bY] = C[ψ₂]/cosh t`, `X = aY`, with:
  * `a = (R₋+R₊)/n`, `b = (R₋−R₊)/n`;
  * `n = (|R₋|+|R₊|)·n_q/|n_q|`, `n_q = R₋^{q²}+R₊^{q²}`.
* **Gauge.** The min-norm gauge is `Y − ψ₂` of minimal `L²(dν)` norm (2-dimensional null space). The choice of `n`, including its phase,
  is part of the definition of the gauge.
* **Functional.** `Re Σ_k dν X_k Ψ̂_{EE,ν_k}/(h₋+h₊)(ν_k)`.
* **Independent reproduction.** The weights were re-solved independently (`d_abel.py`): they agree to `≤ 1.5·10⁻⁶` (ν < 2) and
  `≤ 3·10⁻¹¹` (ν > 25) of max ψ₂, and `dk₀` to `10⁻¹⁰`.

## 5. Errata (v4): fields inside files that are out of date

The files are left byte-identical so that the hashes quoted in the Supplement and in `MANIFEST.txt` remain valid.

| file, field | says | correct statement | source |
|---|---|---|---|
| `nf02_psi.json`, `note` | "certified margin 3/20 except region (+-, J>=400, 3/8<u<=1)" | the flat action of NF02 is **≥ N/5 on every flat heavy state** (both channels, all J, all m ≥ M) [certified]; the corner `+−`, J ≥ 400 is closed by a backward-endpoint lemma (its code unreviewed); sharp ratio `0.2000000103` at `+−`, J = 5, m = M | `CERT_CLOSE.md` §1; UT3 §5.1 |
| `q_deshift_continuum_c4v.json`, `c4d.json`, `certificate` | `eps_W = 1.205e-13` / `4.2e-11` | certified upper bounds `2.673·10⁻¹³` / `9.058·10⁻¹¹` | `HK_RESIDUAL.md` §3.3; UT3 N8 |
| §3 of this README (before v4) | ε_W `1.2·10⁻¹³` / `4.2·10⁻¹¹` "certificates" | as above | — |
| `q3_full_delta0.2.json` | stored doubles define the flat target | the flat target is Q3\* (`q3star_factored.json`); the stored doubles' flat action is negative at `+−`, m = M, J ≥ 1990 | `CERT_CLOSE.md` §2.1 |

## 6. Other published material

* `k2/`: Result Q's k = 2 Abel weight `X^Q` and the continuum k = 2 weights of both functionals (own `README.md`, `SHA256SUMS`).
* `../fig/make_fig.py`: Fig. 1 of the Letter, computed from `nf02_psi.json` and `q3star_factored.json` only (exact term-by-term
  integrals in mpmath); it prints the headline constants (`τ₀ = 2.89126·10⁻⁸`, `2L_AdS = 92.8135` at Δ_g = 1000) as a check.
