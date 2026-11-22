# Result Q: closing H_OO^Q and H_fwd^Q by direct measurement

*21–22 Nov 2026, forward-close agent. Unrefereed, not committed. No existing document or code file was edited; no `bj_*` file was
touched, imported or run. New code: `ads4/code/spinfunc/qf_*.py`, `qf_queue{A,C,D}.sh` (`qf_queueB.sh` was cut short after its first two
jobs; its last job moved to `qf_queueD.sh`). Caches: `spinfunc/cache/qf_*` (git-ignored). Read-only imports: `po_ee_odd` (parity-odd
`⟨TTO⟩` prefactor on the analytic-radial engine), `bg_eeq`, `bg_ee`, `bg_geom`, `bg_ser`, `bg_qan` (projection, `ffn`), `bg_scalar`,
`bg_scalarJ`, `bf_qan` (`mass`, `Rj`, `T`, `weights`), `p1_tto_odd`, `l5_spin`, `l6_spin_action`, `l7_dc2_analysis`, `l8_analysis`,
`resultq/stage6_model`, `stage6_stress`, `stage6_scan`, `common_q`. Existing caches (`bg_eeq_*`, `bh_eeq_*`, `po_eeq_*`,
`bh_helbasis.pkl`, `bh_shadow.pkl`) were read, never written. Compute: ≤ 2 worker processes, `ulimit -v 3 GB`, `OMP/OPENBLAS = 1`,
≤ 1.4 GB per process.*

Grades: **[P]** proved; **[D]** derived; **[DV]** derived and checked numerically; **[NE]** numerical evidence; **[M]** model; **[I]**
open. `c₂` = `m² ×` relative correction. `s_d = 0.04 Σ_ν χ_d(ν) Re Γ_K(ν)` (Stage 7 §2; BELOW_FOLD2 §4). "Engine" = the analytic-radial
`1/Δ` engine of BELOW_FOLD2 on the real section D; its EE beyond-LO equivalence with the crystal cycle is **assumed** (BELOW_FOLD2 §5.4),
here as in BELOW_FOLD2/3 and PARITY_ODD_CLOSE.

@@SUMMARY@@

---

## 1. Which forward channels can still bind: per-channel tolerances (`qf_tol.py` → `cache/qf_tol.txt`) [M]

On Stage 6's HYB state set (logic of `po_ffscan`), per channel `(h, J)`: the K = 0 minimum of FF/N, and the **failure threshold**
`s_fail = min over the channel's states of FF/(a·2u³)` (the `s_d` at which some state of the channel reaches FF = 0; the measured
universal size is `|s_d| ≈ 0.51`). OO rows are Stage 6's `++` rows (the model's identification), so OO J has the `++` J threshold.

| channel | δ = 0.05: K = 0 min @ m/Λ | `s_fail` | δ = 0.2: K = 0 min @ m/Λ | `s_fail` |
|---|---|---|---|---|
| `++`/OO J = 0 | +0.050 @ 1.000 | **+0.58** | +0.200 @ 1.000 | **+1.78** |
| `++`/OO J = 2 | +0.050 @ 1.10 | **+0.83** | +0.200 @ 1.13 | **+2.81** |
| `++`/OO J = 4 | +0.274 @ 1.57 | +11.3 | +0.391 @ 1.63 | +15.0 |
| `++`/OO J = 6 / 8 | +0.377 @ 2.00 / +0.434 @ 2.46 | +28 / +52 | +0.492 @ 2.11 / +0.567 @ 2.59 | +36 / +66 |
| `+−` J = 4 | +0.079 @ 1.013 (300); +0.059 @ 1.004 (1000) | +1.04 / +0.75 | +0.228 / +0.209 @ gap | +2.38 / +2.10 |
| `+−` J = 5 | +0.050 @ 1.20 | +1.78 | +0.200 @ 1.19 | +5.72 |
| **`+−` J = 6** | +0.050 @ **1.475** | **+1.99** | +0.207 @ gap 1.020 (300); +0.202 @ 1.006 (1000) | +9.5 / +9.4 |
| `+−` J = 7 | +0.146 @ 1.83 | +8.9 | +0.400 @ 1.69 | +20.0 |
| `+−` J = 8 | +0.181 @ 1.026 (300); +0.088 @ 1.008 (1000) | +16.6 / +10.4 | +0.372 @ 1.026; +0.290 @ 1.008 | +30 / +27 |
| `+−` J = 9 | +0.225 @ 2.34 | +25 | +0.416 @ 2.22 | +41 |
| `+−` J = 10 … 13 | ≥ +0.243 | ≥ +35 | ≥ +0.354 | ≥ +51 |

(Δ_g = 300 and 1000 agree except where both are given.) **Only five channels have `s_fail` within a factor 4 of the measured size:
`++`/OO J = 0, 2, `+−` J = 4, 5, 6.** `++` J = 0, 2, `+−` J = 4, 5 were measured before (BELOW_FOLD2/3, PARITY_ODD_CLOSE); this note
measures **OO J = 0, 2** (§2) and **`+−` J = 6** (§3), the last two near-binding channels, and then checks the rest (`+−` J = 7, 8, 9, 11;
`++` J = 6, 8) for an order-of-magnitude surprise (`|s_d| ≳ 9`). Every other channel needs `|s_d| > 8.9` (δ = 0.05) / `> 15` (δ = 0.2)
to fail, i.e. 17–30× the measured universal size.

## 2. H_OO^Q: the parity-odd channels at even J = 0, 2, measured (`qf_eeq.py odd`, `qf_qan.py oo`) [NE + M]

### 2.1 Construction
`po_ee_odd.images` (PARITY_ODD_CLOSE §2.1, read-only) builds the EE k = 4 integrand with the conserved **symmetric parity-odd** `⟨TTO⟩`
structure of `p1_tto_odd` on both sides. At even J this is CHL's `X^O_same` (one structure for J = 0 and every even J ≥ 2; checked at all
Δ used: `conserved_odd` returns one symmetric solution, 2 monomials at J = 0, 15 at J = 2). The image is 1 × 1, so the ν-independent shadow
coefficient drops out with the per-state constant `c`. The `S⁰` sum (both `n₀` sheets) is built into the engine's ψ ∈ (−π, π]. Chain:
`qf_qan.py oo` = `po_qan5`'s chain with `γ = γ_{++} = 0.75` (the OO channel's HYB row is `++`: Stage 3 §2.1, Stage 6 §1.1). The `++` image
at **the same Δ** (BELOW_FOLD2's `bg_eeq` caches, or new `qf_eeq even` images at 1.3Λ) goes through the **same** chain (REFCORR = 1 for both),
so the comparison is like for like. Cost: J = 0 1.5 min, J = 2 11 min per state (depth 5, ≤ 0.9 GB); references 11 s / 83 s.

Channel identification [DV]: `Re Γ_exact/Re Σ_j R_j(ν − ij)² = 0.73–0.78` for ν ≥ 4 in OO (the `++` covariant constant 0.75, not the
`+−` 20.75), at every state.

### 2.2 Result

| state | Δ_g | `s_d` H0 | H2 | H3 | H4 | SB | XP | FF/N δ = 0.05 (K = 0 → meas., H0; SB) | δ = 0.2 |
|---|---|---|---|---|---|---|---|---|---|
| **OO J = 0, gap `1.000Λ`** (Δ = 300) | 300 | **+0.509 ± 0.016** | +0.510 | +0.509 | +0.508 | +0.525 | +0.515 | +0.0500 → **+0.0061 ± 0.0013** (+0.0048) | +0.2000 → **+0.1427 ± 0.0017** (+0.1410) |
| **OO J = 0, gap** (Δ = 1000) | 1000 | **+0.510 ± 0.002** | +0.510 | +0.510 | +0.509 | +0.525 | +0.515 | +0.0500 → **+0.0060 ± 0.0001** (+0.0047) | +0.2000 → **+0.1426 ± 0.0002** (+0.1409) |
| OO J = 0, `1.05Λ` (Δ = 314.925) | 300 | +0.510 ± 0.014 | +0.510 | +0.510 | +0.508 | +0.525 | +0.515 | +0.1477 → +0.1124 (+0.1113) | +0.2891 → +0.2436 (+0.2423) |
| OO J = 0, `1.30Λ` (Δ = 389.55) | 300 | +0.510 ± 0.010 | +0.511 | +0.510 | +0.509 | +0.525 | +0.516 | +0.5147 → +0.5022 (+0.5018) | +0.6027 → +0.5876 (+0.5871) |
| **OO J = 2, `1.05Λ`** (Δ = 314.935) | 300 | **+0.511 ± 0.010** | +0.512 | +0.511 | +0.510 | +0.526 | +0.517 | +0.0611 → **+0.0267 ± 0.0007** (+0.0256) | +0.2214 → **+0.1812 ± 0.0008** (+0.1800) |
| **OO J = 2, `1.05Λ`** (Δ = 1049.93) | 1000 | **+0.510 ± 0.001** | +0.510 | +0.510 | +0.509 | +0.525 | +0.515 | +0.0611 → **+0.0267 ± 0.0001** (+0.0257) | +0.2214 → **+0.1813 ± 0.0001** (+0.1801) |
| OO J = 2, `1.30Λ` (Δ = 389.56) | 300 | +0.511 ± 0.007 | +0.512 | +0.511 | +0.510 | +0.526 | +0.517 | +0.1697 → +0.1532 (+0.1527) | +0.2837 → +0.2645 (+0.2639) |

* errors = size of the a⁴ term (the a⁵ term is ≈ 1/15 of it, BELOW_FOLD2); per-order `s_d` (H0) J = 0, Δ_g = 300: `+1.662 −3.949 +2.516
  +0.265 +0.016`; J = 2: `+10.907 −11.414 +0.891 +0.116 +0.010` (large cancelling low orders, then ×10 per order, as for `++`).
* **Sign: unfavourable, exactly as `++`.** The OO channel carries `+U(ν)`, i.e. `σ_OO = σ_{++} = +1`: per harmonic `Γ_K/Σ|R_j|` =
  `−1.03 −4.58 −8.30 −10.41 −9.92 −7.55 −2.11 +1.45 −0.58` at ν = 0.5, 1, 1.5, 2, 2.5, 3, 4, 6, 10 (OO J = 0, Δ_g = 300): equal to the
  `++` profile at the same state to ≤ 0.002 (§2.3), and the same function `U` across J = 0, 2, 1.0–1.3Λ and both Δ_g to ≤ 0.07.
* **m-independence within the channel [NE]:** `s_d` = +0.509 / +0.510 / +0.510 at 1.00 / 1.05 / 1.30Λ (J = 0), +0.511 / +0.511 at 1.05 / 1.30Λ
  (J = 2), spread ≤ 0.001 ≪ truncation.
* Covariant part (bg_covfit rule): `δγ = +0.015` (2 % of `γ_{++}`) at every OO state; the rest (+0.495 H0) is the universal non-covariant `U`.

### 2.3 Is OO = `++` for the k = 4 functional? The helicity lemma beyond order 0 [NE]

Same Δ, same chain, same projection; the two images are **different functions** (ratio `I_OO/I_{++}` = 15.788 (J = 0, Δ = 300),
19.244 (J = 2, Δ = 314.9); reference 15.997 / 19.248 — the basis-normalisation constant `κ(Δ)` of PARITY_ODD_CLOSE §1.1 (F)):

| state | reference ratio `I_OO/I_{++}`: ν-variation (ν < 41) | state ratio: `m²·max|r(ν)/r(0) − 1|` | `max|Γ_K^OO − Γ_K^{++}|/Σ|R_j|` (ν < 41) | `s_d(OO) − s_d(++)` (every gauge) |
|---|---|---|---|---|
| J = 0, Δ = 300 | `1.0·10⁻⁹` | 0.0042 c₂ | `1.3·10⁻³` | `−4.7·10⁻⁶` |
| J = 0, Δ = 314.9 | `1.0·10⁻⁹` | — | `1.2·10⁻³` | `−4.1·10⁻⁶` |
| J = 0, Δ = 389.5 | `1.0·10⁻⁹` | — | `6.6·10⁻⁴` | `−2.3·10⁻⁶` |
| J = 0, Δ = 1000 | `1.0·10⁻⁹` | — | `7.7·10⁻⁵` | `−1.1·10⁻⁶` |
| J = 2, Δ = 314.9 | `1.1·10⁻⁹` | 0.0054 c₂ | `1.1·10⁻³` | `−4.3·10⁻⁶` |
| J = 2, Δ = 389.6 | `1.1·10⁻⁹` | — | `6.2·10⁻⁴` | `−2.4·10⁻⁶` |
| J = 2, Δ = 1049.9 | `1.1·10⁻⁹` | — | `7.1·10⁻⁵` | `−1.0·10⁻⁶` |

* **Order 0:** the reference images (Δ = 2·10⁴) are proportional to 10⁻⁹ — the order-0 lemma (PARITY_ODD_CLOSE §1.1 (C)) holds for the
  k = 4 kernel as predicted (the Lemma-B factor `R` is polarisation-blind, so the saddle argument is unchanged) [DV].
* **Beyond order 0:** `Γ_K^OO − Γ_K^{++}` is ≤ `1.3·10⁻³` c₂ per harmonic at Δ_g = 300 and falls like `Δ^{−2.3…−2.5}` (300 → 389.5 → 1000:
  ×2.0, ×17), i.e. **`Γ_OO − Γ_{++} = O(m⁻²) c₂`**, a relative `O(Δ⁻⁴)` effect — the **same scaling** PARITY_ODD_CLOSE §1.1 (D) found for
  the k = 2, 3 functionals (`X/D = O(Δ⁻⁴)`). Projected: `|s_d(OO) − s_d(++)| ≤ 5·10⁻⁶`, gauge-independent, 3000× below the truncation
  error and 10⁵× below the failure threshold.
* **Verdict: the helicity lemma's consequence `Γ_OO = Γ_{++}` holds for the k = 4 forward functional** — exactly at order 0 [DV], and
  to `O(m⁻²) c₂` (here ≤ 10⁻³ c₂, `|Δs_d| ≤ 5·10⁻⁶`) beyond, at J = 0 and 2, both Δ_g [NE]. Why the beyond-order-0 corrections cancel in
  the OO/`++` difference remains [I] (as for k = 2, 3).

**Consequence for Result Q.** OO J = 0 at the gap is a **co-binding** state of the full set with `++` J = 0: FF/N = +0.1427 ± 0.0017 (δ = 0.2;
SB +0.1410) and +0.0061 ± 0.0013 (δ = 0.05; SB +0.0048), identical to `++` J = 0 to 10⁻⁶. OO J = 2: +0.0267 / +0.1812 at 1.05Λ (channel
minimum with the m-independent `s_d`, §4: +0.0196 / +0.1665 at 1.089 / 1.118Λ, the `++` J = 2 values). **H_OO^Q is no longer an assumption at
J = 0, 2: it is measured** [NE + M; NE given D ≡ crystal for EE beyond LO, M for the HYB rest of the heavy action].

## 3. H_fwd^Q: the remaining forward channels

### 3.1 Fixed order or J/Δ resummation? (`qf_jcheck.py` → `cache/qf_jcheck.txt`) [DV for the scalar; M for the transfer to EE]
The engine expands at fixed J in `a = 1/Δ`; the spin dependence enters through `P_J`-type factors whose expansion parameter is `J²/Δ`.
A J/Δ-resummed EE engine exists only as the other agent's `bj_*` work (even J, parity-even; not touched, not imported). The question
is therefore whether fixed order K = 4 is adequate at J ≤ 9–11 near the gap. On the **scalar skeleton**, where both forms exist
(`bg_scalar` fixed order vs `bg_scalarJ` N = 12, R = 4, the latter validated to ≤ 0.04 c₂ against the exact skeleton up to J = 32),
`m²[log(R_fixed/R_resum) − const]`, max |·| on ν ≤ 4 / 4–10 / 10–18 (ν ≈ 6 excluded: zero of `g + g̃`, ill-conditioned log):

| (Δ, J) | J²/Δ | K = 4 | K = 6 |
|---|---|---|---|
| (304, 4) | 0.05 | 0.11 / 0.21 / 3.7 | 0.004 / 0.37 / 0.18 |
| (306, 6) | 0.12 | 0.30 / 1.4 / 24 | 0.004 / 0.20 / 0.36 |
| (307, 7) | 0.16 | 0.45 / 2.6 / 46 | 0.004 / 0.12 / 0.76 |
| (308, 8) | 0.21 | 0.65 / 4.5 / 79 | 0.004 / 0.04 / 1.7 |
| (309, 9) | 0.26 | 0.91 / 7.2 / 125 | 0.005 / 0.05 / 3.5 |
| (311, 11) | 0.39 | 1.7 / 15 / 263 | 0.009 / 0.26 / 13 |
| (441.9, 7 / 8 / 9) | 0.11–0.18 | ≤ 0.31 / ≤ 2.5 / ≤ 40 | ≤ 0.001 / ≤ 0.09 / ≤ 0.55 |
| (1007–1011, 7 / 9 / 11) | 0.05–0.12 | ≤ 0.05 / ≤ 0.5 / ≤ 7.4 | ≤ 0.001 / ≤ 0.015 / ≤ 0.03 |
| (336, 24) (control) | 1.71 | 27 / 500 / 5·10⁴ | 0.4 / 11 / 2.7·10³ |

(K = 6 column and the control row: ν ≈ 6 not excluded.)

The fixed-order error grows like `(J²/Δ)·x²` and sits at large ν, where the derivative weight `χ_d` has almost nothing (93 % of `Σ|χ_d|`
at ν < 2.5, 7 % at 2.5–10, < 10⁻³ above). Projection proxy (`δΓ ≈ T_st·L_scalar`, i.e. assuming the EE image carries the scalar's
J/Δ truncation error; ν > 18 extrapolated ∝ ν²): `δs_d` (H0 / SB) = +0.006 / −0.047 (306, 6), +0.011 / −0.087 (307, 7), +0.017 / −0.15
(308, 8), +0.026 / −0.22 (309, 9), +0.051 / −0.44 (311, 11); at 1.475Λ ≤ +0.009 / −0.08 (J ≤ 9); at Δ_g = 1000 ≤ +0.002 / −0.016 (J ≤ 11).
SB is the gauge that weights the large-ν tail, hence its larger number. **Verdict:** fixed order K = 4 is adequate for `s_d` at J ≤ 9
(systematic ≤ 0.03 H0 / ≤ 0.22 SB at Δ_g = 300; ≤ 0.016 at Δ_g = 1000) against failure thresholds ≥ 8.9 for J ≥ 7 and 1.99 for J = 6.
"Large J" in the sense of needing the resummation starts around J ≈ 11 at Δ_g = 300; there the channels are measured at **Δ_g = 1000**,
where J²/Δ ≤ 0.12 and the fixed-order error is ≤ 0.02 — `s_d` being Δ_g-independent (BELOW_FOLD2/3, PARITY_ODD_CLOSE, §2 here), the
Δ_g = 1000 value stands for the channel [M]. Caveat from the measurement (§3.3): at J = 7 the Δ_g = 300 and 1000 values differ by
0.074, more than this proxy (0.011) and the a⁴ term (0.031) — the EE spin structure (`P_J` derivatives up to 8th order) is more J-sensitive
than the scalar; the proxy is therefore a lower estimate, and the J ≥ 7 channels are quoted from Δ_g = 1000.

### 3.2 Construction and runs
* **Odd J (`+−`, parity-odd):** `qf_eeq.py odd` = `po_ee_odd.images` (1 × 1), chain `qf_qan.py pmodd` = `po_qan5` (γ = γ_{+−} = 20.75).
* **Even J ≥ 4 (`++`, `+−`, parity-even):** `qf_eeq.py even` = `bg_eeq.images`, all four `(a, c)` entries of the 2 × 2 prefactor matrix,
  × inverse shadow matrix (`l6_spin_action.shadow`), in the order-0 helicity basis `T0(Δ)` at the same Δ (`qf_helbasis.py` = `bh_helbasis`
  with its own cache), chain `qf_qan.py even` = `bh_qan4`. **Validation:** on BELOW_FOLD3's J = 4 images the copied chain reproduces
  `s_d(++) = +0.517 ± 0.017`, `s_d(+−) = −0.516 ± 0.039`, FF/N `+0.118 / +0.278`, mixing `m²κ ≈ +1.4` exactly.
* Depth: K = 4 (depth 5) where the channel is near-binding or the comparison needs it (`+−` J = 6 gap, J = 7 gaps); K = 3 (depth 4) for the
  rest, whose thresholds (`s_fail ≥ 8.9`) dwarf a K = 3 truncation. Costs (depth 5 / depth 4): J = 6 2 × 2: 105 / 39 min; J = 7 odd: 38 / 14
  min; J = 9 odd: — / 19 min; J = 11 odd: — / 22 min; peak RSS 1.6 GB.
* Masses: the gap state of each J at Δ_g = 300 (`Δ = 300 + J + 0.0298`) and/or 1000 (`Δ = 1000 + J + 0.0090`), and `m = 1.475Λ`
  (Δ_g = 300), which is `+−` J = 6's δ = 0.05 binding mass and inside the requested 1.2–1.5Λ window for every channel.

### 3.3 Result [NE + M]

| channel, state | Δ_g | K | `s_d` H0 | H2 | H3 | H4 | SB | XP | FF/N δ = 0.05 (K=0 → meas.; SB) | δ = 0.2 |
|---|---|---|---|---|---|---|---|---|---|---|
| **`+−` J = 6, gap `1.020Λ`** | 300 | 4 | **−0.564 ± 0.012** | −0.565 | −0.564 | −0.563 | −0.579 | −0.570 | +0.0775 → +0.0877 (+0.0880) | +0.2077 → **+0.2200** (+0.2203) |
@@J6_1000@@
| **`+−` J = 6, `1.475Λ`** | 300 | 3 | **−0.545 ± 0.010** | −0.546 | −0.545 | −0.544 | −0.560 | −0.551 | +0.0497 → **+0.0633** (+0.0637) | +0.3454 → +0.3626 (+0.3631) |
| `++` J = 6, gap `1.020Λ` | 300 | 4 | +0.530 ± 0.100 | +0.531 | +0.530 | +0.529 | +0.545 | +0.536 | +0.7811 → +0.7699 | +0.8497 → +0.8380 |
| `++` J = 6, `1.475Λ` | 300 | 3 | +0.565 ± 0.53 | +0.565 | +0.565 | +0.564 | +0.580 | +0.570 | +0.5635 → +0.5532 | +0.7488 → +0.7377 |
| **`+−` J = 7, gap `1.007Λ`** | 1000 | 4 | **−0.517 ± 0.001** | −0.518 | −0.517 | −0.516 | −0.532 | −0.523 | +0.6571 → +0.6627 | +0.6805 → +0.6880 |
| `+−` J = 7, gap `1.023Λ` | 300 | 4 | −0.591 ± 0.031 | −0.591 | −0.591 | −0.589 | −0.605 | −0.596 | +0.6574 → +0.6641 | +0.6846 → +0.6935 |
| `+−` J = 7, `1.475Λ` | 300 | 3 | −0.539 ± 0.28 | −0.539 | −0.539 | −0.538 | −0.554 | −0.544 | +0.3603 → +0.3700 | +0.4687 → +0.4795 |
@@J8@@
| **`+−` J = 9, gap `1.009Λ`** | 1000 | 3 | **−0.507 ± 0.32** | −0.507 | −0.507 | −0.506 | −0.521 | −0.512 | +0.6633 → +0.6676 | +0.7790 → +0.7831 |
| `+−` J = 9, gap `1.030Λ` | 300 | 3 | −0.43 ± 1.67 | −0.43 | −0.43 | −0.43 | −0.44 | −0.43 | +0.6812 → +0.6856 | +0.8069 → +0.8109 |
| `+−` J = 9, `1.475Λ` | 300 | 3 | −0.48 ± 1.01 | −0.48 | −0.48 | −0.48 | −0.49 | −0.48 | +0.7536 → +0.7588 | +0.8278 → +0.8340 |
| **`+−` J = 11, gap `1.011Λ`** | 1000 | 3 | **−0.48 ± 0.67** | −0.48 | −0.48 | −0.48 | −0.49 | −0.48 | +0.3709 → +0.3738 | +0.8407 → +0.8433 |

Errors = size of the last included order (K = 3 rows: the a³ term; the projected series converges geometrically, ratio 9–30 per order,
so these errors are conservative by roughly that factor). Per-order `s_d` (H0): `+−` J = 6 gap `+37.75 −37.37 −0.86 −0.10 +0.012`;
J = 7 (1000) `+59.30 −57.90 −1.85 −0.065 −0.0008`; J = 9 (1000) `+111.7 −102.7 −9.16 −0.32`; J = 11 (1000) `+176.4 −158.2 −18.0 −0.67`:
the leading two orders grow like J² and cancel to ~10 %, then the series converges.

**Reading.**
* **Sign and size: every new `+−` measurement is favourable and of the universal size**, `s_d = −0.48 … −0.59` (all gauges, spread
  ≤ 0.016), i.e. `σ_{+−} = −1` with `|U|`-projection 0.51 ± 0.08. Per harmonic, the clean states reproduce `−U(ν)`: `+−` J = 7 (Δ_g = 1000)
  `+1.03 +4.58 +8.32 +10.46 +9.98 +7.53 +1.75 −1.83 +0.99` at ν = 0.5 … 10 vs `++` J = 0 `−1.03 −4.58 −8.30 −10.41 −9.92 −7.55 −2.11 +1.45
  −0.58`. The parity-odd `+−` exchange at J = 7, 9, 11 thus carries the same correction as at J = 5 (PARITY_ODD_CLOSE) and as the even
  `+−` at J = 4, 6.
* `++` J = 6: `+0.53 ± 0.10` (gap) / `+0.57 ± 0.53` (1.475Λ): the universal `+U`, unfavourable, irrelevant (`s_fail ≥ 28`).
* **Helicity-even part** `(s_d^{++} + s_d^{+−})/2` at J = 6: `−0.017` (gap), `+0.010` (1.475Λ) (J = 4: `+0.0001`): `σ_h`-universality holds to
  ~3 % of `U` at J = 6.
* **Δ_g-dependence at larger J [NE; origin I].** At Δ_g = 1000 the `+−` values are universal (−0.510 / −0.512 / −0.517 / −0.507 / −0.48 at
  J = 4 / 5 / 7 / 9 / 11). At Δ_g = 300 they drift with J: `s_d(300) − s_d(1000)` = −0.006 (J = 4), −0.027 (J = 5), −0.074 ± 0.031 (J = 7),
  @@J6DIFF@@. At J = 7 the difference is 2.4× the a⁴ term and sits in orders k = 1–3 (`−54.98 −4.41 −0.52` at 300 vs `−57.90 −1.85 −0.065`
  at 1000), i.e. it is a J²/Δ-enhanced finite-Δ term that the fixed-order series does compute, not obviously a truncation artefact (the
  scalar proxy of §3.1 gives only +0.011 for (307, 7)). Whether it is physical (`J²/m²` corrections to `U`) or a slow (J²/Δ)ⁿ tail beyond
  the last order cannot be decided from these data [I]. It is **favourable** in sign, and ≤ 0.08 against thresholds ≥ 8.9 for J ≥ 7. The
  Δ_g = 1000 numbers (J²/Δ ≤ 0.12) are the clean measurement of each channel; at J = 9 the Δ_g = 300 rows (J²/Δ = 0.18–0.26, K = 3) are
  only consistency checks (errors ±1–1.7). Per-harmonic a⁴ terms scale as expected (`+−` J = 7: ≤ 7.8 per `|R|` at Δ_g = 300, ≤ 0.77 at
  1000, ratio ≈ (Δ ratio)²).
* **Per-harmonic caution [NE].** At J ≥ 6 and Δ_g = 300 the per-harmonic a⁴ (a³) term is comparable to or larger than the harmonic value at
  ν ≈ 2–4 (`++` J = 6 gap: −10.4 at ν = 3 vs the value −6.5), dominated by a covariant `∝ Σ_j R_j(ν − ij)²` piece (a J²/Δ renormalisation
  of γ_h; covariant fit `δγ = +0.035 / +0.19` for `++` J = 6 at gap / 1.475Λ). Its derivative projection is small (the projected series
  converges), but the per-harmonic `U`-profile at J ≥ 6, Δ_g = 300 should not be read off these fixed-order images.
* **m-dependence within a channel [NE]:** `+−` J = 6: −0.564 (1.02Λ) / −0.545 (1.475Λ); J = 7: −0.591 (1.02Λ) / −0.539 (1.475Λ) (Δ_g = 300,
  the difference within the J/Δ error); J = 9: −0.43 / −0.48 (both ± 1); OO/`++` J = 0: 0.509 / 0.510 / 0.510 (1.0 / 1.05 / 1.3Λ); J = 2:
  0.511 / 0.511 (1.05 / 1.3Λ). No m-dependence beyond the truncation/J-Δ errors anywhere between the gap and 1.475Λ.
* **Helicity mixing at J = 6 [NE + M].** `m²κ ≈ +49` (`T0` from x ≥ 0.15) / `+210` (x ≥ 0.25) — basis-dependent, because the order-0
  Laplace basis at J = 6 near the gap is less sharp (rank-2 residual `7.9·10⁻⁵`, order-0 off-diagonal `2.6·10⁻²`, vs `1.7·10⁻⁶`, `1.5·10⁻³`
  at J = 4). The **diagonal** `s_d` are basis-independent to `10⁻⁴` (both bases). `s_off` (H0 … SB) = `−0.01 … −0.21` / `+0.10 … −0.76`; in
  the 2 × 2 positivity condition this enters as `(a·2u³ s_off/N)² ≤ (FF/N)_{++}(FF/N)_{+−}`: at most `(0.03·0.76)² = 5·10⁻⁴` against
  `0.77 × 0.088 = 0.068` (δ = 0.05, gap) — negligible.

@@REST@@

