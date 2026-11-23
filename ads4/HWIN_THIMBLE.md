# H_win by per-wave thimble integration (Result D, Ω^EE fixed: NF02, min-norm gauge)

*20–21 Nov 2026, H_win-thimble agent. Unrefereed, not committed. No existing document or code file edited. New code
`ads4/code/spinfunc/pw_*.py`, caches `spinfunc/cache/pw_*`. ≤ 2 cores (`OMP_NUM_THREADS=1`), < 1 GB per process. Project modules
imported read-only (in-process wrappers only).*

Grades: [P] proved; [D] derived; [NE] numerical evidence; [M] model; [I] open.

**STATUS: final for this session. Nothing is running.**

> **Update (N2_SUBREGGE.md, unreviewed): N2's premise is right (X, not ψ₂, multiplies the window data), but the claimed leverage collapse does not happen. A dir–cv estimator systematic at (1100, 100) is now the dominant uncertainty.**
> * **What the rows are [DV].** Window rows Γ are raw images / (R₋+R₊), and X multiplies them; Re ΣX𝒜[P_β] = Σψ₂P_β holds to 2·10⁻¹². Γ has poles at the real zeros of R₋+R₊ (ν = 1.171, 6.643), which explains HLJ_CLOSE's "local Γ singular at every J". H_win should be stated on the regular aΓ_c with weight |Y|.
> * **T_B with the exact weight:** (1100, 100): 1409 → 1439 (|X|) / 1152 (|Y|). Over all sample states the change is −28 % to +3 %. LARGESPIN's "1.1" and "3.2" were artefacts (ψ₂ + |X−ψ₂| overcounts; a different floor and structure were switched in; signed ring tolerance 31 vs measured 4.5–6.6). The D_h excess is an a_h normalisation, already in de-shift form (cross-checked against k = 3, which shares a_h). The published Γ − κ contains a stripped-κ artefact of ±120–460 c₂.
> * **Corrected H_LJ(1000) leverage at (1100, 100)** (cv estimator, raw-form refit): **3.0–4.2 central; ≥ 4.1 at 2σ with the full ring class.** (364, 64): 20–36. **With the dir estimator: 0.8–1.4.** This dir–cv systematic is unresolved and dominant.
> * **New explicit hypothesis R1 [I] for EE beyond leading Regge:** ring-type window corrections enter in de-shift form. It holds exactly for the a_h mechanism and for the scalar. In stripped form they would cost 27–228× the ring action at J = 150–344.
> * **Subleading Regge:** CHMRS (3.20) is checked against the exact scalar skeleton (9 states). All its terms sit on the state's own ring, with cost 6.5·10⁻⁴ (J = 100) to ≤ 6·10⁻⁵ (J ≥ 150) of the floor, so margins are ×17 to ×10⁴. The graviton window data fit the same ring classes, with non-ring remainder consistent with zero. The **spinning formula is not derived [I]** (obstacle: the EE = Σ_b identity decomposes the kernel, not the spinning block; it needs ⟨TTO⟩ structure ratios). At each finite 1/m² order the corrections are ring-located or index-shifted at the state's β [D-sketch].
> * **Corrections:** the "leverage 1.1" in the LARGESPIN box above is superseded. Quote "3.0–4.2 (cv), ≥ 4.1 at 2σ; 0.8–1.4 with dir". HWIN_THIMBLE: refit the D_h class in de-shift form, after which the (364, 64) remainder is 7–13 c₂.

## 0. Summary (read this first)
1. **Exact per-wave integrands, scalar and graviton** [P structure + NE, §1]. The graviton spin-J contraction is re-derived in closed
   form (a bidegree-(4,4) derivative of `c_J(gg′)^{J/2} f(t)`). It is linear in `f`, so `f = W±` (exact Legendre-Q waves) splits EE
   and C3H exactly. `f = P_J` reproduces the l9/w0 Fourier contraction to `10⁻¹³` and has no `ε²4^J` engine noise at any J.
2. **Per-wave thimbles remove the spin sign problem** [NE, §2–4]. Each wave sits on the n3_hybrid contour of its own saddle; state and
   `Δ_ref = 2·10⁴` share common random numbers. At equal points the errors are 2–11× smaller than with `P_J` on the z2/fw contour
   (scalar `(1100,100)`, ν = 16.3: ×11; EE `(1064,64)`, ν = 16.7: `δΓ ≈ 100` vs `580 c₂`). Raw state/reference ratios agree with the
   independent fw records to `10⁻⁴`. Exceptions:
   * **C3H `+−` fails per wave** in the window: it is a cancellation between the waves, and the graviton waves have Jacobi-Q poles at
     `θ = 0` (§1.3). For this key `P_J` on the same CRN is used instead.
   * **Per-wave bias in the window** [NE]: at `ℓθ₀ ≲ 3` the per-wave scalar misses the exact scalar by `10²–10³ c₂` (≤ 3σ). The
     candidate mechanism is the Q-cut at `θ = 0` swept between the two wave contours [D-sketch]. At `ℓθ₀ ≳ 4` the error falls to
     `≲ 10–30 c₂` (§3.2). The graviton uses the scalar control variate; `|dir − CV|` is carried as a systematic.
3. **The fold limits the window, not the spin phase** [NE, §3]. Below `ν ≈ 10` the contour variance grows for every mode
   (`1 %` at ν = 10, `4–15 %` at 8.3, `> 50 %` at 6–6.7). Widening the proposal makes it worse. **With c = 3 the Bessel window lies
   entirely below ν ≈ 10 for J ≥ 200 (Δ_g = 1000) and J ≥ 64 (Δ_g = 300).** There it is not measurable by this method, and in part
   not even for J = 64, 100 at 1000.
4. **Measured window values** [NE, §4, table §5]:
   * In the sign-free form `|Γ − κ| ≤ T_B`, **H_win holds at (1064, 64)** on ν ∈ [13, 25].
   * It is **violated at (1100, 100)** at the c = 3 window edge: C2 `+−` ≈ −2150 vs T_B = 1409, in the per-wave, `P_J` and fw data
     alike.
   * It is also exceeded in the X-weighted range at (364, 64) (≈ 1100 vs 255) and (400, 100).
   * **Every excess is D_h-shaped.** The `cs` fit gives `c₀ ≈ 3·10³ – 6·10³`. After removing the D_h class, the remainder is
     `≲ 550` (1000, inside T_B) and `≲ 150–250` (300, at the level of T_B ≈ 250).
5. **Grade of H_win** (§6):
   * The sign-free form is **refuted as a sufficient condition at J = 100** (1000) [NE].
   * The refined form, `|Γ − κ − c₀D_h| ≤ T_B` with the D_h class sent to its own, much larger leading-Regge tolerance (HLJ §5.2,
     [D-sketch]), is [NE] on the measured sub-windows: 1000, J = 64, 100, ν ≥ 13; 300, J = 64, X-range ν ∈ [10, 33].
   * Everything else stays [M]/[I]: ν < 10 at every J, J ≥ 150 (1000), J = 34, 48 (300, no basis rows), (400, 100).

## 1. Exact per-wave integrands (`pw_wave.py`, `pw_grav.py`, checks `pw_tests.py` → `cache/pw_tests.txt`)
### 1.1 Scalar [P + NE]
`exp(λF + G + logJac)·Q_J` (l1_crystal at spin J) `= exp(λF + G₀ + logJac)·P_J(t)` point-wise (ratio 1 to `10⁻¹⁰`, J = 8, 20),
`t = Z·W/√(Z²W²)` (HLJ §1.2). Waves `W± = ½[P_J ∓ (2i/π)Q_J]` (Ferrers Q; upward recurrence; vs mpmath `≤ 2·10⁻¹²` up to J = 300;
`W± ≃ e^{±i((J+½)θ−π/4)}/√(2π(J+½)sin θ)` checked). Each wave has only the logarithmic singularity of Q at `θ = 0`.
`√(Z²W²)` is continued (n1_pw replay) from the n3 reference point next to the crystal, so the labels ± are path-consistent.

### 1.2 Graviton [D + NE]
Every ⟨TTO_J⟩ numerator monomial `V₅^{a₅}H₁₅^{b₁₅}H₂₅^{b₂₅}` has `b₁₅ + b₂₅ ≤ 4` (degree count in `z_a, z_b`), so as a polynomial
in the O polarisation `F_a = (v·z₅)^{J−4}·T_a(z₅)` with a quartic `T_a`. With `z^α(v·z)^{J−4} = [(J−4)!/J!]∂_U^α(U·z)^J|_{U=v}` and the
spin-J pairing of two J-th powers `⟨(U·z)^J,(U′·z)^J⟩ = c_J(g(U,U)g(U′,U′))^{J/2} f(t_U)`, `f = P_J`, `c_J = 2^J J!²/(2J)!`
(`g` = inverse lightcone metric on plain coefficients), the contraction is `[(J−4)!/J!]² Σ T_{a,α}T′_{c,β} ∂^α∂^β G_f` — linear in
`f`. `∂^α∂^β G_f` are Taylor coefficients of a bidegree-(4,4) truncated power series in `(δU, δU′)` (f composed through its
t-derivatives to order 8). Spin-2 indices and the D structures (EE: `E⊗E`; C3H: `q⊗o_H`) exactly as w0/l9.
Checks (random points next to the crystal, `ν = 14`, `Δ = 1000 + J`):

| J | EE: pw(P_J)/Fourier − 1 | C3H | (W₊ + W₋)/P_J − 1 |
|---|---|---|---|
| 8 | `2.3·10⁻¹⁴` | `2.6·10⁻¹³` | `4·10⁻⁶ / 6·10⁻⁷` |
| 20 | `2.7·10⁻¹³` | `4.5·10⁻¹³` | `5·10⁻⁹ / 1·10⁻⁸` |
| 40 | `2.0·10⁻⁸` | `1.6·10⁻⁸` | `2·10⁻¹² / 6·10⁻¹³` |

(J = 40: the difference is the Fourier engine's known `ε²4^J` noise, FOLD_WINDOW §1.1 (`3.8·10⁻⁸` there). The `W₊+W₋` loss at small
`Jθ` is the cancellation of the large Q-derivatives.)

### 1.3 A structural limitation [D]
A split of `P_J` into pieces `∼ e^{±iℓθ}` cannot be analytic at `θ = 0` (`θ ∝ √(1−t)`). For the scalar the waves have a log
singularity (integrable). The graviton contraction differentiates `f` up to 8 times, so the graviton waves have poles
`(1−t)^{−k}`, k ≤ 8, at `θ = 0` (Jacobi-Q), which are not integrable on a real 7-cycle that meets `θ = 0`. On the wave contours
`θ = 0` sits at `≈ ℓθ₀/|ℓ∇θ| ≈ ℓθ₀` Gaussian widths (`|ℓ∇θ| ≈ 1`, lj_saddle), i.e. **inside the window `ℓθ₀ ≲ 3–5` the graviton
waves are only usable where their pole contributions stay small** — empirically EE and C3H `++` at `(1064, 64)`, `ν = 16.7`
(`ℓθ₀ ≈ 2`) were fine, C3H `+−` was not (§2). The scalar is unaffected.

## 2. Method (`pw_core.py`, `pw_points.py`)
* Mode 0: `P_J` on the n3 contour of the J-growing scalar flow function (= `F.saddle(ν_f, Δ, 0)`: its J-terms cancel exactly), i.e.
  the z2/fw production contour. Modes ±1: `W±` on `C(ρ₁ = 3.5, ν_f)` built on the saddle of `Φ± = λF + G₀ + logJac + log W±(t)` at
  `ν_f` (homotopy `λ log W±`, λ: 0 → 1, holomorphic-stencil Newton, `|∇Φ| ≤ 3·10⁻⁹` at the saddle).
* `ν_f`: the n3 rule (`max(ν, 15)` above 13, 13 below) or `PW_NUF`. Proposal: n3_hybrid.adapt on the state's `W₊` (moment-matched
  Gaussian), re-used for `W₋`, for `P_J` and for the reference (CRN, 8 Sobol replicas × 2⁹ × 2 antithetic).
* Graviton on the same contour points; the log replay along the stored flow paths uses a contraction-free copy of the integrand with the
  identical log sequence, the contraction is evaluated at the end points only (≈ 6 ms/point).
* Output in the z2/fw record format; `pw_gamma.py` injects the records into fw_gamma (helicity bases from the fw_k2 EE rows, a_h with the
  fw anchors x = 0.1/0.2, κ-image) — identical Γ conventions to FOLD_WINDOW.

## 3. Precision and validation against the exact scalar (`pw_scal.py`, `cache/pw_scal.txt`) [NE]
`Δc₂ = m²Re(d_st − d_ref)ex₂/E₂` (must vanish), jackknife over 8 replicas, 8192 points per state:

| state, ν | ν_f | waves: `Δc₂` | `P_J`: `Δc₂` | per-mode relative error (waves / P_J) |
|---|---|---|---|---|
| (1100,100), 16.3 | 16.3 | `−607 ± 208` | `−175 ± 2330` | 0.4 % / 0.9 % |
| (1100,100), 13.3 | 15 | `+322 ± 128` | — | 0.6 % |
| (1100,100), 13.3 | 13.3 | `+446 ± 492` | — | 0.7–1.1 % |
| (1100,100), 10.0 | 13 | `−2220 ± 2295` | `−2498 ± 2021` | 1.2 % / 1.2 % |
| (1100,100), 6.0 | 13 | `(−2.3 ± 1.4)·10⁶` | `(−2.7 ± 1.0)·10⁶` | 150 % / 94 % |
| (400,100), 8.3 | 13 | `−3362 ± 2592` | `−3381 ± 6272` | 4–7 % / 7 % |
| (400,100), 6.67 | 13 | `(2 ± 5)·10⁴` | `(−4 ± 3)·10⁴` | 44–130 % / 50 % |
| (1064,64), 16.7 (graviton run) | 16.7 | `−343 ± 202` | `+82 ± 1281` | 0.45 % / 0.8 % |

`ex₂/E₂` varies strongly (the factor is `1.9·10⁴` at (1100,100), 13.3 and `3·10⁵` at 16.3), so the same relative precision maps to very
different `c₂` errors. The `10 c₂` target of the task is **not reached** at 8k points: the achieved scalar validation is `|Δc₂| ≲
200–600` at `ν ≥ 13` (1000), consistent with zero at 1–3σ (the 2.5σ/2.9σ nodes at (1100,100) are re-checked with an independent seed,
§3.1).

### 3.1 Repeat with an independent seed (`_s2`, 16k points) and the (1064, 64) graviton-run scalars
* `(1100, 100)`: `ν = 16.3`: `−63 ± 263` (was `−607 ± 208`); `ν = 13.3`: `+218 ± 109` (was `+322 ± 128`; combined `+265 ± 83`, 3σ).
* `(1064, 64)` (waves, same samples as the graviton, `_w`): `ν = 10 / 13.3 / 16.7 / 20.8 / 25`: `−1048 ± 2589 / −994 ± 426 /
  −343 ± 202 / −449 ± 206 / −95 ± 49`; `P_J` on the same CRN (`_wP`): `−2022 ± 1751 / +426 ± 1057 / +82 ± 1281 / +521 ± 2687 / +4 ± 849`.
* **Reading [NE]:** the per-wave scalar is 2–11× more precise than `P_J`, but over the 7 nodes with `ν ≥ 13` `χ² ≈ 27` (8 dof):
  there is a **systematic of the per-wave estimate at the level `|Δc₂| ≈ 100–1000` (relative `≈ 10⁻³`) in the window**
  (`ℓθ₀ ≈ 2–3`), sign state-dependent. A candidate mechanism [D-sketch]: each wave has the logarithmic cut of `Q_J` at `θ = 0`; the
  two waves are moved to different contours, and the part of the cut swept between them contributes `±∫P_J` over a region of small
  but not negligible weight while `θ = 0` lies `≈ ℓθ₀` Gaussian widths away. It should vanish like the weight there as `ℓθ₀` grows;
  queue C tests `ν = 33, 50` at (1064, 64) and 25, 40 at (1100, 100). For the graviton the scalar control variate (CV) removes the part
  of this bias that is common to scalar and graviton (dir − CV ≈ the scalar `Δc₂`, §4); the remainder is carried as a systematic
  `|Γ_dir − Γ_CV|`.

### 3.2 Exact-scalar validation, all nodes (waves; `Δc₂`, ± jackknife; `ℓθ₀ = (2J+1)ν/m`)

| state | ν (ℓθ₀): `Δc₂` |
|---|---|
| (1040, 40) | 13.3 (1.0): `−359 ± 234`; 20.8 (1.6): `−56 ± 101`; 30 (2.3): `−50 ± 19` |
| (1064, 64) | 10 (1.2): `−1048 ± 2589`; 13.3 (1.6): `−994 ± 426`, `−815 ± 378` (2 seeds); 16.7 (2.0): `−343 ± 202`; 20.8 (2.5): `−449 ± 206`; 25 (3.0): `−95 ± 49`; **33.4 (4.0): `+10 ± 9`**; 50 (6.1): `+27 ± 12` |
| (1100, 100) | 10 (1.8): `−2220 ± 2295`; 13.3 (2.4): `+322 ± 128`, `+218 ± 109`; 16.3 (3.0): `−607 ± 208`, `−63 ± 263`; 20.8 (3.8): `−15 ± 36`; 25 (4.6): `+99 ± 53`, `+92 ± 48`; **40 (7.3): `−16 ± 10`** |
| (364, 64) | 10 (3.6): `−182 ± 239`; 11.7 (4.2): `+242 ± 170`; 13.3 (4.8): `+123 ± 49` |
| (400, 100) | 10 (5.2): `+123 ± 160`; 13.3 (6.9): `+63 ± 89`; 8.3 / 6.67: `±3·10³` / `±5·10⁴` (fold) |

**Reading [NE]:**
* The task's target `|Γ_s − κ| ≲ 10 c₂` is **met only outside the window** (`ℓθ₀ ≥ 4`: `+10 ± 9`, `−16 ± 10`).
* Inside the window (`ℓθ₀ ≲ 3`) the per-wave scalar has a bias of `10²–10³ c₂` (state dependent; up to 3σ) on top of statistical
  errors of `50–400` (8–16k points). This is still 2–20× below T_B(c = 3) at Δ_g = 1000.
* At Δ_g = 300 the scalar nodes are all at `ν ≥ 10`, i.e. above the window.

## 4. Graviton `Γ − κ` in the window (fw_gamma conventions, `pw_gamma.py`, `cache/pw_gamma.txt`) [NE]
### 4.1 `(1064, 64)`, Δ_g = 1000 (`m = 1060.5`, ν_B = 24.7 / 41.1 for c = 3 / 5; T_B = 1913 (C2), 14989 (C3H), c = 3)
a_h from the fw anchor `x = 0.2` (and the fw bulk nodes 0.4, 0.6, which enter nothing else); bases fw_k2 EE rows (`sv3 = 5·10⁻¹⁰`).
CV estimator, `++` (`+−` within `≤ 25` of it for C2), `c₂`, ± jackknife:

| ν | 10.0 | 13.3 | 16.7 | 20.8 | 25.0 |
|---|---|---|---|---|---|
| C2, waves | `−1451 ± 2579` | `−209 ± 216` | `−67 ± 26` | `−64 ± 36` | `−86 ± 56` |
| C2, P_J (same CRN) | `−3750 ± 1682` | `−859 ± 539` | `−70 ± 174` | `−3 ± 329` | `−78 ± 398` |
| C2, waves, dir estimator | `−1875 ± 4867` | `−834 ± 433` | `−228 ± 89` | `−305 ± 121` | `−442 ± 194` |
| C3H `++`, waves | `−4184 ± 3689` | `−249 ± 766` | `−335 ± 722` | `+134 ± 208` | `−230 ± 101` |
| C3H `+−`, P_J (waves fail, §1.3) | `+2024 ± 3584` | `+316 ± 940` | `−980 ± 622` | `−367 ± 1263` | `−185 ± 646` |
| FOLD_WINDOW (fw n3, C2 `++`, CV) | `−2653 ± 1619` | `−111 … −344, +112` (± 80–385) on 13.3–33.4 | | | |

* C2 on `ν ∈ [13.3, 25]`: `|Γ − κ| ≤ 210 ± 216` (CV), `≤ 830 ± 430` (dir): **inside `T_B = 1913` by ≥ 2σ including the dir − CV
  systematic**; 3–8× more precise than FOLD_WINDOW at 16.7–25. `ν = 10`: not decisive (`±2600`).
* C3H: every node inside `T_B(C3H) = 14989` (worst `−4184 ± 3689` at `ν = 10`, i.e. `≤ 1.2·10⁴` at 2σ).

### 4.2 `(1100, 100)`, Δ_g = 1000 (`m = 1093.9`, ν_B = 16.3 / 27.2; T_B(C2) = 1435 (`++`) / 1409 (`+−`), T_B(C3H) ≈ 12 000, c = 3;
c = 5: 734 (C2))
Anchors: the fw direct nodes `x = 0.1, 0.2`; bases fw_k2 EE rows. `Γ − κ` (CV, ± jk), C2:

| ν | 10.0 | 13.3 | 16.3 | 20.8 | 25.0 |
|---|---|---|---|---|---|
| C2 `++`, waves | `−1397 ± 3395` | `+1224 ± 287` | `+919 ± 41` | `−334 ± 11` | `−1682 ± 23` |
| C2 `+−`, waves | `−389 ± 3416` | `−902 ± 292` | `−2154 ± 43` | `−1783 ± 12` | `0 ± 25` |
| C2 `++` / `+−`, P_J (same CRN) | `−6254 / −5254 (±3700)` | `+25 / −2061 (±1075)` | `+1155 / −1917 (±517)` | `−399 / −1848 (±110)` | `−1786 / −104 (±290)` |
| C2 `++` / `+−`, fw records only (ν = 13.3, 16.7, 20.8, 25, 33.4) | — | `+1262 / −1053 (±570)` | `+371 / −2123 (±172)` (16.7) | `−1494 / −1972 (±48)` | `−2616 / −154 (±106)` |
| C3H `++` waves / `+−` P_J | `(−7 ± 8)·10⁶` / `+2064 ± 3922` | `−1403 ± 7.7·10⁴` / `+528 ± 3093` | `+1.1·10⁴ ± 1.4·10⁴` / `+1427 ± 1235` | `−941 ± 307` / `+2188 ± 614` | `−2177 ± 107` / `+905 ± 564` |

* **The per-wave raw data agree with the independent fw records** where both exist (state/reference helicity ratios at ν = 20.8, 25:
  `1·10⁻⁴`; at 13.3: `2·10⁻³`, within errors), at 3–30× smaller error (C2).
* **But `Γ − κ` itself is ≈ ±(1–3)·10³ at the window edge and above it, in all data sets** (fw's own records included), and depends on the
  per-row a_h fit by ±(300–1500) (node sets/estimators: the a_h fit on x ≤ 0.25, ν ≥ 15 is ill-conditioned, FOLD_WINDOW §2.6).
  C2 `+−` at `ν = 16.3–16.7` (the c = 3 window edge): `−2154 ± 43` (waves), `−1917 ± 517` (P_J), `−2123 ± 172` (fw) — **above
  T_B(c=3) = 1409 in every variant** [NE]; c = 5 (T_B = 734) is exceeded at 13–25 by several nodes.
* The excess is **D_h-shaped**: the `cs` fit (`c₀D_h + s xD_h′/(2J+1) + r`, fw_gamma) on ν ∈ [12, 34] (pw nodes 13.3, 16.3, 20.8,
  25 + fw nodes 16.7, 33.4; CV) gives C2 `c₀ = +6532 / +4887`, constant `r = +389 ± 117 / −51 ± 140`, rms residual `157 / 217`
  (`++ / +−`); fw records alone: `c₀ ≈ +4470 / +4180`, `r = +250 ± 232 / −91 ± 240`, residual 229 / 157. So the window correction is
  dominantly `c₀D_h` — the leading-Regge-shaped class (an a_h-type normalisation), as HLJ §5.2 [D-sketch] anticipated — and **what
  is left after the D_h class (`|r| + residual ≲ 550`, ≤ 800 at 2σ) is inside T_B(c=3) = 1409**. C3H `+−` (P_J/fw): `r ≈ +1000 ± 380`,
  inside T_B(C3H) ≈ 12 000; C3H `++` per wave is unusable at ν ≤ 16 (§1.3) and the fw/P_J values (`|Γ−κ| ≲ 4000`) are inside 12 000.
* (1064, 64) for comparison (pw nodes, ν ∈ [12, 34]): C2 `c₀ ≈ −120 … −140`, `r = −423 ± 472 / −372 ± 483`, residual ≤ 24.

### 4.3 `(364, 64)`, Δ_g = 300 (`m = 356.7`, ν_B = 8.3 / 13.8; T_B(C2) = 295 / 286, with W_X 255 / 247 (c = 3))
The window `ν < 8.3` is below the fold (§3): not measured. Nodes 10, 11.7, 13.3 (the X-weighted range `[ν_B, 25]`), with the fw nodes
15–33.4 kept (a_h, PW_KEEPFW=1). C2 `Γ − κ` (CV):

| ν | 10.0 | 11.7 | 13.3 | fw: 15 … 33.4 |
|---|---|---|---|---|
| C2 `++`, waves | `−797 ± 218` | `−1132 ± 109` | `−1109 ± 29` | `−529, +173, +759, +976, +566, −70, −847, −14` (± 11–50) |
| C2 `++`, P_J (same CRN) | `−982 ± 299` | `−952 ± 233` | `−905 ± 231` | same |
| C3H `++` waves / `+−` P_J | `−833 ± 803` / `+1839 ± 972` | `−1091 ± 470` / `+932 ± 264` | `−912 ± 138` / `+763 ± 256` | `|·| ≤ 1350` |

* Waves and P_J agree within 1σ; waves 2–8× more precise.
* `|Γ − κ| ≈ 800–1130` on 10–13.3 **exceeds the sign-free T_B(c=3, +X) = 247–255** — again entirely in the D_h class: the `cs` fit on
  all 11 nodes has `c₀ = +2840 / +2805` (C2 `++ / +−`), `r = +0.5 ± 50 / +3.7 ± 53`, rms residual `136 / 145`; P_J: `r = +11 ± 39`,
  residual 105. C3H: `c₀ ≈ ±(3.5–5.4)·10³`, `r = +135 ± 114 / +105 ± 423`, residual 204 / 902 (P_J: `+40 ± 57 / −82 ± 76`, 315 / 243)
  vs T_B(C3H) ≈ 2800. Remainder after the D_h class: C2 `≲ 150–250` — at the level of T_B (≈ 250), not clearly inside.

### 4.4 `(400, 100)`, Δ_g = 300 (`m = 385.6`, ν_B = 5.8 / 9.6; T_B(C2, +X) = 191 / 179 (c = 3))
* The window `ν < 5.8` is below the fold: not measured.
* Per-wave nodes ν = 10, 13.3 (scalar check `+123 ± 160`, `+63 ± 89`), C2 `++ / +−` (CV):
  * waves: `−509 ± 511 / −1626 ± 517` (ν = 10), `−866 ± 62 / +371 ± 71` (13.3);
  * P_J: `+952 ± 1147 / −192 ± 1172`, `−1037 ± 366 / +193 ± 380`.
* These exceed T_B(+X) ≈ 185 by up to ×9.
* **The D_h decomposition is not reliable here.** The fw nodes 16.7–33.4 that fix a_h have references overwritten by the (1100, 100)
  run (FOLD_WINDOW §2.6). The `cs` fit still gives `c₀ ≈ ±(1.3–7)·10³` and `r ≈ −40 ± 70 … +140 ± 140`. [I].

## 5. Table: `|Γ − κ|` vs T_B (c = 3; C2 unless stated; CV; ± jackknife) [NE where measured]

| Δ_g, J (state) | window `ν < ν_B` | measured sub-range | `max|Γ − κ|` there | T_B (C2 / C3H) | after the D_h class: `|r|` + rms resid | verdict |
|---|---|---|---|---|---|---|
| 1000, 64 (1064, 64) | 24.7 | 13.3–25 (+10) | `209 ± 216` (C3H `++` `335 ± 722`, `+−` `980 ± 622`); ν = 10: `1451 ± 2579` | 1913 / 14989 | `≈ 400 + 20` | **inside** on 13–25 [NE]; ν < 13 [I] |
| 1000, 100 (1100, 100) | 16.3 | 13.3–16.3 (+10) | **`2154 ± 43`** (`+−`, 16.3); `1224 ± 287` (13.3); ν = 10: `±3400` | 1409 / 12066 | `≈ 390 + 160 / 50 + 220` | sign-free **exceeded**; refined form inside [NE]; ν < 13 [I] |
| 1000, 150 | 11.3 | — | not measured (window at the fold; no basis rows for Δ = 1150) | 1184 / 11402 | — | [M] |
| 1000, 200 | 8.8 | — | not measured (below the fold) | 1114 / 12581 | — | [M] |
| 1000, 300 | 6.3 | — | not measured (below the fold; no rows) | 1057 / 18594 | — | [M] |
| 300, 34 | 14.4 | — | not computed (no K₂ basis rows for Δ = 334; 10 ≤ ν ≤ 14.4 would be accessible) | 387 / 3052 | — | [I] |
| 300, 48 | 10.6 | — | not computed (no rows; window at the fold) | ≈ 300 / 2850 | — | [I] |
| 300, 64 (364, 64) | 8.3 | none (X-range 10–33.4) | X-range: `1132 ± 109` (11.7) | 255 (+X) / 2806 | `≈ 0–50 + 140` | sign-free exceeded in the X-range; refined ≈ at T_B [NE, weak]; window [I] |
| 300, 100 (400, 100) | 5.8 | none (X-range 10, 13.3) | X-range: `1626 ± 517` (`+−`, 10) | 179 (+X) / 3505 | unreliable (fw refs overwritten) | [I] |
| cross-check (1064, 64) vs fw | | 13.3–25 | pw C2 `−64 … −209`; fw `−111 … −344, +112` | | | consistent; pw 3–8× more precise |

## 6. Grades
* Exact per-wave integrands (scalar W± identity; graviton closed-form contraction, linear in f) — [P] structure, implementation [NE]
  (`10⁻¹³` against the Fourier engine).
* "Per-wave thimbles remove the spin-J sign problem" — [NE] at `ν ≥ 13`: 2–11× smaller errors at equal cost.
* Per-wave window bias `10²–10³ c₂` at `ℓθ₀ ≲ 3` — [NE]; the cut-sweep mechanism is [D-sketch].
* "The window ν < ν_B lies below the fold for J ≥ 200 (1000), J ≥ 64 (300); there no contour of the n3/pw family gives < 10 %
  errors" — [NE] (ν_f = 13 with adapted or widened proposals; the n4 thimbles of EE_DG300 are no better).
* **H_win, sign-free form** (`max_{ν<ν_B}|Γ − κ| ≤ T_B`):
  * [NE] satisfied at (1064, 64) on ν ∈ [13, 25];
  * [NE] **violated** at (1100, 100) on its measured sub-window (2154 ± 43 vs 1409), and in the X-weighted range at (364, 64) and
    (400, 100).
  * So the sign-free T_B condition **cannot hold as stated for J ≳ 100 at 1000 and J ≥ 64 at 300**.
* **H_win, refined form** (D_h-class part handled by the leading-Regge c-class tolerance of HLJ §5.2, remainder `≤ T_B`):
  * [NE] on the measured sub-windows at 1000 (J = 64, 100, ν ≥ 13);
  * [NE, weak] at 300 (J = 64, X-range: remainder ≈ T_B);
  * [I] at (400, 100);
  * [M] for ν < 10 at every J and for J ≥ 150 (1000), J ≥ 34 (300).
  * The refined form needs HLJ §5.2's [D-sketch] (that the c₀D_h class is harmless; c-class sensitivity was computed at 1000 only).
* **Net effect on H_LJ:**
  * 1000: stays [D + NE + M]. The window clause must be restated in the refined form; with that, J = 64, 100 are supported above ν ≈ 13
    [NE].
  * 300: stays [D] (bulk) + [D-sketch / M] (window). The measured X-range values make the D_h-class argument necessary, not
    optional.

## 7. What would close the gaps (value order)
1. A contour that works below the fold, ν ≈ 2–10, at < 10⁻³. This is the real obstruction for J ≥ 150 (1000) and J ≥ 48 (300).
   * Widening the proposal fails; ν_f = 13 with an adapted proposal fails below ν ≈ 9.
   * Untested: complex ν_f (PW_NUF accepts `a+bj`; the saddle path code exists in n3_flow); exact 1D fibre quadrature (§1.3 of HLJ:
     the fibre is a Beta-function integral); the n4 A± pair with the wave factor.
2. Make the D_h decomposition part of the H_win statement (the `cs` class with its own tolerance) and compute the c-class sensitivity
   at Δ_g = 300.
3. K₂ basis rows (order-0 M0 suffices for the basis; a_h needs anchors at x = 0.1, 0.2) for (1150, 150), (1300, 300), (334, 34),
   (348, 48). Anchors are bulk nodes, where the per-wave method works best.
4. Fresh references for (400, 100) (FOLD_WINDOW §2.6) to make its D_h decomposition reliable.
5. Remove the per-wave window bias: put both waves on one contour inside the window (ℓθ₀ ≲ 3) or quantify the cut sweep.

## Hand-over
* Code (`ads4/code/spinfunc/`): `pw_wave.py` (waves; self-test), `pw_grav.py` (exact per-wave spin-J contraction), `pw_core.py`
  (evaluator, wave saddles, block integration on n3_hybrid contours), `pw_points.py` (CRN state/reference nodes; resumable),
  `pw_scal.py` (exact-scalar check), `pw_probe.py` (helicity-ratio precision), `pw_gamma.py` (Γ − κ in fw_gamma conventions),
  `pw_tests.py` (point-wise identities), queues `pw_queueA.sh`, `pw_queueB.sh` (re-run to resume after a restart).
* Caches `cache/pw_pts_J<J>_D<Δ><tag>.pkl` (`_w`: waves + graviton; `_wP`: P_J on the same CRN; `_v`, `_s2`: scalar-only
  independent seeds; `_s`, `_sP`, `_sS`, `_g0/_g1/_tA/_tB`: tests), logs `cache/pw_queue{A,B,C}.log`, `cache/pw_scal.txt`,
  `cache/pw_gamma.txt` (+`.pkl`).
* **Nothing is running.** The queues finished: A = (1064, 64) and (364, 64) nodes; B = (1100, 100) and (400, 100); C = scalar
  validation. Re-running a queue skips finished nodes.
* pw_gamma options: `PW_KEEPFW=1` keeps all fw nodes at other ν (needed at (364, 64) and (400, 100) for a_h); `FW_MINNODES=3` for
  short rows. Ranges: `pw_gamma.py J:Δ 12 34` for the `cs` fits quoted in §4.2.
* Commands: `PW_TAG=_w python3 pw_scal.py 100:1100`; `PW_TAG=_w python3 pw_probe.py 64:1064`;
  `FW_EST=cv PW_TAGS=_w python3 pw_gamma.py 64:1064` (and `PW_TAGS=_wP` for the P_J records).
