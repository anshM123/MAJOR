# Referee report: `ads4/SPINNING_E2E.md` and `ads4/code/e2e/`

Referee: independent, hostile. Date: 16 Nov 2026. Object: commit `df18d80`. I did not edit anything under review.
Scratch scripts and their outputs are in `ads4/code/referee_e2e/`. They run on a copy of `e2e/` with freshly compiled `.so`.

## Summary verdicts

| | Question | Verdict |
|---|---|---|
| A | Analytic derivation (§4, §6.1) | **Correct as far as I could check.** I found no error in the leading formula, in `b^{(Δ34)}`, in (★), in the conserved-structure normalisation, or in the Weyl-reflection argument for `⟨OφOφ⟩`. That includes the flip component, where I tried to break it. Two gaps need closing: the exact mode symmetry is asserted but never proved (m3), and term-by-term integration of the pointwise expansion is assumed (m4). The extension of the `O(1/m)` argument to `⟨TTTT⟩` is **unproven** (M2). |
| B | Numerics | **Sound.** I reproduced the cached tensors from scratch to the last bit and reproduced the table entries. Precision is adequate. I found no tuned constant and no circularity. Mutation tests are detected at `O(1)` and `O(1/m)`, as they should be. |
| C | Honesty of the "confirmed" labels | **Mostly honest, with one material overstatement.** Item 11 (`−2l²P′/P` "confirmed: −8.00") is not supported as stated for `l = 2` (M1). Items 7 and 9 overgeneralise in their summary forms (m1, m2). Items 3, 5 and 14 are acceptable. |
| D | Physical applicability to the AdS4 graviton bound | **Narrow.** It establishes the absolute leading one-pair heavy Regge moment, per `λ²`, at `k = 4, 6` and `ν ≪ m`. It does not establish the `⟨TTTT⟩` formula, its `±` structure, the mean-field or positivity normalisation, the graviton-relevant `k = 2, 3`, or uniformity in `ν/m` (M3). |
| E | Novelty | **Incremental and technical.** I found no prior paper that does this specific computation. The ingredients are standard or close to it (see E). |

No fatal issues.

## (A) Analytic derivation

1. **Leading formula §4.5.** Every `O(1)` ingredient is pinned by the numerics: factors of 2, `4^{l+1+Δφ}`, signs, `|T_l|`, the `M` table, and the powers of `m` and `sinh β`. A multiplicative error would show up as `R ↛ 1`. My mutation "×2" gives `R → 1.97…1.99`. The `(−1)^{J+l}` and `T_1 < 0` conventions in `predict.py` are consistent with the `|T_1|` of the text.
2. **`b^{(Δ34)}`.** It is even in `Δ34`, so `b^{(1−Δφ)} = b^{(Δφ−1)}`, and the label refers to the `b^l` seed `(l+1, Δφ+l)`. Three things constrain it:
   * the word-by-word identity `c_w(Δ,J)/(b_w c_w(3−Δ,−1−J)) = const`, which fixes it modulo period-2 functions of `Δ34`;
   * direct quadrature at 0.1 %, which excludes an `O(1)` period-2 factor;
   * the `O(1/m)` end-to-end discrimination, which I re-ran (see B).

   The extension is consistent with the standard `±Δ34/2` splitting of the light-transform and shadow coefficients (Kravchuk–Simmons-Duffin; Karateev–Kravchuk–Simmons-Duffin, arXiv:1809.05111). It should be cited as such rather than presented as new.
3. **(★).** Valid. In the CHMRS gauge fixing, the residual group consists of Lorentz transformations and `x → λx, y → y/λ`. The latter scales all four frame points by `1/λ`, so it is a *global* dilatation and acts trivially on the polarisations. A Lorentz-covariant tensor kernel therefore passes through verbatim. This argument is only sketched in §3.1; it should be written out in one sentence.
4. **Conserved structure and normalisation.** Uniqueness, the `L_0`, `ℓ_i` map and the OPE-limit link `λ²(−½)^J W_diff/L_0²` are each checked independently. The OPE-limit test uses a separate contraction formula and converges linearly. I accept these as [DVP].
5. **Weyl antisymmetry, no `O(1/m)` (§6.1). I tried to break it.**
   * *Pointwise.* At fixed `(X, Y)` the frame vectors `(ŷ, e_∥, e_⊥)` are unchanged by `(x, y, m) → (−x, −y, −m)`. So `z_1` and `z_3` are expanded in the *same* frame, and the swap `z_1 ↔ z_3` is a genuine swap of mode labels.
   * *Analyticity in `1/m`.* Every factor of `b𝒯` has the form `(1 + O(1/m))^{(linear in m)}` or `|m|^{power}`, so `b𝒯` is genuinely analytic in `1/m` at fixed `(X, Y)`. There is no Stokes problem in continuing `m → −m`.
   * `out_vertex.txt` confirms the antisymmetry numerically: for example `+1.29835` against `−1.29805`.
   * *At the moment level, flip component.* The mode `(l, l)` is invariant under the swap, so its `O(1/m)` term vanishes on its own, without any further input.
   * *At the moment level, preserved component.* The mode `(l, −l)` maps to `(−l, l)`. The argument then needs the *exact* identity `Π̂^{(l,−l)} = Π̂^{(−l,l)}`. This follows from the `(12) ↔ (34)` symmetry of the s-channel partial wave combined with the improper map `x → −x`, which is harmless for a parity-even structure of even total rank.

   So the argument does **not** rely on a polarisation assignment that fails for the flip component. It does, however, rely on (i) that exact mode symmetry, which the text asserts but only checks numerically (m3), and (ii) exchanging the `1/m` expansion with the `(X, Y)` integral (m4).

## (B) Numerics: what I ran

The scripts are in `ads4/code/referee_e2e/`. I ran a representative subset rather than the full `run_all.sh fast`, which takes more than an hour mostly in validation steps. Both `.so` files compiled cleanly with `gcc -O3/-O2 -fopenmp`, plus `-lquadmath` for the quad build.

* **Fresh recomputation without the cache** (`recompute.py`).
  * `l = 1`, `(44.565, 8)`: 103 s.
  * `l = 2`, `(45.861, 20)`: 126 s.
  * Fresh against cached moments: `0` to `10⁻¹⁶`.
  * Table entries reproduced: `l = 1`, `k = 6`, `ν = 0`, `(1,−1)`: `0.947946` (doc §5.2: `0.947946`); `l = 2`, `(2,−2)`: `0.918418` (`out_t4…l2`: `0.918418`).
* **Precision** (`prec.py`). I compared double and quad evaluation of the full `l = 2` stripped tensor at the largest case, `(115.11, 52)`, at four moment-dominant `(u, t)` points. The relative difference is `≤ 5·10⁻⁹`.
* **No circularity.**
  * `predict.py` uses only `b`, `γ_a`, `F^ν_μ`, the analytic `M` table and `|T_l|`. Nothing comes from the blocks and nothing is fitted.
  * The two-mode extraction uses the *predicted* `M`s. But the raw `(l,−l)` mode ratios alone already tend to 1, and the other `(2l+1)² − 2` modes overdetermine the extraction.
  * Moments with `ν ≠ 0` are complex, because `F^ν_{2l}` is complex. Only the *ratios* are real. §5.1 should say so (m6).
* **Mutation tests** (`mutate.py`). I applied each mutation to the prediction and refit the `k = 6` ratios with the doc's ansatz (A) and with a four-parameter ansatz `c_1…c_4`.

| mutation | `c_1` (preserved), `l = 1` | `c_1` (preserved), `l = 2` | `R` at `m ≈ 100` |
|---|---|---|---|
| none | `−0.007 … +0.001` (A); `≤ 0.020` (4-par) | `−0.008 … +0.002` (A); `≤ 0.025` (4-par) | 0.977–0.997 |
| `b^{(0)}` | `−0.57 … −0.65` | `−0.57 … −0.65` | shifted by about 0.6 % |
| `b^{(l+1−Δφ)}`, the natural wrong `Δ34` | `−0.51 … −0.58` | `+1.03 … +1.13` | |
| `×(1 + 1/m)` | `+1.00` | `+1.00` | |
| `×2` | — | — | `→ 1.96–1.99` |
| `×(1 + 3/m²)` | undetected (as expected) | undetected | |

  The shadow-coefficient mutation is clearly detected. The test is genuine, but only at `O(1/m)`: the end-to-end test cannot see `O(1/m²)` errors, which is correct for a leading-order claim.
* **Weak spot.** For the `l = 2` flip component at `β ≈ 0.4`, `c_1 = −0.27` (A) or `−0.31` (4-par) when all five masses are used. The `b^{(0)}` mutation gives `−0.85`. The numerics therefore barely separate the correct formula from a wrong one there. The `−0.20` quoted in the doc drops the `m = 102` point.

## (C) Labels; numbered issues

### Major

**M1. Item 11 and §6.3: the `l = 2` "confirmation" of `−2l²P′/P` is an extrapolation, not a test. The claim "valid for all β tested" is false for `l = 2`.**

The fit basis `{1, cosh β, P′/P, sinh²β P′/P}` on two values of `β` is equivalent to a separate fit `c_β + d_β P′/P` at each `β`. "`a_2`" is then just the linear-in-`sinh²β` extrapolation of `d_β` to `β = 0`. I refit per `β` (`perbeta.py`):

| | `d_β` at `β = 0.4` | `d_β` at `β = 1` | `d_β` at `β = 2` |
|---|---|---|---|
| `l = 1` | `−2.01` | `−2.01` | `−4.1` (errors are a few units here) |
| `l = 2` | `−8.34` | `−10.79` | `−35.8` |

For `l = 1`, fixing `d = −2` leaves residuals of at most 0.1, so this is a genuine confirmation at `β = 0.4` and `1`. For `l = 2`, fixing `d = −8` leaves residuals up to 6.8 at `β = 1`, against a quoted error of `≲ 0.3`.

*Fix:* change the label to "l = 1 confirmed at β = 0.4, 1; l = 2: β-dependent coefficient (−8.3, −10.8), consistent with −8 only as β → 0 under a fitted sinh²β term". Delete "valid for all β tested (not only β → 0)" for `l = 2`, and correct §0 item 6 and STATUS row 7 accordingly.

**M2. §6.1, last sentence: "applies verbatim … hence to `⟨TTTT⟩`" is unproven.**

The argument used two facts:
* the shadow transform maps the *unique* conserved structure to the same analytic family, with a scalar coefficient `1/b`;
* the same coupling appears on both sides.

For `⟨TTO⟩` there are two parity-even conserved structures, plus a parity-odd one. The shadow transform acts as a matrix. Even if that matrix is diagonal, for example in the CHL helicity basis, its eigenvalues `κ_h(Δ,J)` differ. Ratios such as `κ_h/κ_{h′}` that are not Weyl-even would give `O(1/m)` terms weighted by `λ_hλ_{h′}`.

*Fix:* downgrade the statement to a conjecture, or prove that the shadow matrix in the differential basis is `(1/b) × (Weyl-covariant)`.

**M3. Scope for the AdS4 bound.**

STATUS row 6 ("absolute leading formula in terms of OPE coefficients confirmed end-to-end … no fitted constant") is literally true for `⟨OφOφ⟩`. Read as progress on "the four-stress-tensor Regge formula", it overstates, because none of the following is touched:
* (a) the `±` structure `|c^{++}|² ± |c^{+−}|²`, and in particular the opposite-sign `g_1²S_{13}S_{24}` (`B_3`) term, which cannot appear in parity-even `⟨OφOφ⟩`. The `S` factors are checked only by `H²` quadrature;
* (b) the 2×2 map from `λ^{(a)}_{TTO}` to `c^{h_1h_2}`, and Lemma V for both vertices;
* (c) the mean-field (`2P/m^{2k}`) or unitarity normalisation. The standard basis is singular at double twists, and reality and positivity of `λ²` in this basis are not shown;
* (d) `k = 2, 3`, which are the moments the graviton sum rules use. The doc itself says they are inaccessible here. Whether `1 + O(1/m²)` survives term-by-term integration at low `k` is untested (m4);
* (e) uniformity in `ν/m`. A bound functional with support at `ν ∼ m` needs the spinning analogue of the uniform CHMRS (2.60), and §7 shows that the single-moment formula fails at `ν/m ≈ 0.3`;
* (f) the `TT → φφ` channel;
* (g) any `Δφ` other than 1.75.

*Fix:* rewrite STATUS row 6 as "one-spinning-pair leading heavy Regge moment (k = 4, 6, ν ≪ m), absolute per λ²: confirmed; ⟨TTTT⟩ formula: not established".

### Minor

* **m1.** "`|c_1| ≤ 0.008`" depends on the ansatz. With a `1/m³` term included, the preserved values reach 0.025 and the flip values 0.06. The conclusion is unchanged, since these are still small against the mutation value of 0.57, but the range should be quoted. Item 9 and §0 item 5 should state that the fits are inconclusive for the `l = 2` flip at small `β` (see the weak spot in B).
* **m2.** §0 item 4 and §0.1(2) generalise "`m²(R−1)` converging" to `β = 2`, `k = 4` at `m ≈ 100`, and (in item 7) `k = 2`. There:
  * `β = 2` has `m ≤ 70`;
  * `l = 2`, `β = 2`, `ν = 3` flip gives `m²(R−1) ≈ 3, 176, 375, 454`, which is not converged;
  * `k = 4` at the largest `m` is marked `**`;
  * at `k = 2` only `R_pres` for `m ≲ 40` is usable.

  Mark these as "consistent", not [DVE].
* **m3.** Prove the exact moment symmetries `(q_1,q_3) → (q_3,q_1)` and `→ (−q_1,−q_3)` that §6.1 uses. They are currently stated as fact and checked only to about `10⁻⁵`.
* **m4.** State that the `1/m` expansion is exchanged with the `(X,Y)` integral. Give the non-uniform region (`u → 0`, `u^{k−1/2}` behaviour) and the condition on `k` for the exchange to be legitimate.
* **m5.** Reproducibility notes.
  * `data_modes/` and `*.pkl` are git-ignored, so `run_all.sh fast` on a fresh clone still recomputes everything (about 2.5 h). `fast` only skips `t4b`.
  * §12 lists `out_casimir_terms.txt`, `out_structures.txt` and `out_laplace_table.txt`, which are not in the repository.
* **m6.** §5.1 "imaginary parts < 10⁻⁵ of the real parts" refers to the ratios. The moments themselves are complex for `ν ≠ 0`.
* **m7.** Repeat one end-to-end case at `Δφ = 2.5`, for which `gen_phi_*_2p5.py` already exists, to exclude `Δφ`-dependent compensating errors.
* **m8.** Item 14 is fine. "[D]" there means only that the polynomial was identified. Uniformity in `J` is numerical, at `J ≥ l`, three values of `m` and one `k`.

## (D) What this establishes for the `⟨TTTT⟩` Regge formula

**Established**, for one conserved spinning Regge pair (`J` or `T`), a parity-even exchange `(Δ, J)` with the unique conserved `⟨OφO⟩` structure, `Δφ = 1.75`, `k = 4, 6` and `ν = O(1)`:
* the leading heavy value of all tensor Regge moments, with absolute normalisation per `λ²` (standard-basis structure). Both helicity components at the spinning pair are confirmed, with corrections `1 + O(1/m²)` in Weyl-odd variables;
* the helicity-preserved harmonic is `P_{−½+iν}(cosh β)` and the flip harmonic is `F^ν_{2l}`. This is the one-pair factor `R` of the earlier conditional formula.

**Not established** (see M3): the `⟨TTTT⟩` formula itself; its `|c^{++}|² − |c^{+−}|²` term; `2P/m^{2k}`; positivity normalisation; `k = 2, 3`; `ν ∼ m`; the second channel. The "no `O(1/m)`" property for `⟨TTTT⟩` is conjectural (M2). Only the flat-limit `β → 0` value of the helicity `1/m²` term is supported for `l = 2` (M1).

## (E) Novelty

* **CHMRS** (arXiv:2106.10274) covers scalars only. **CLS-D** is arXiv:2311.04271, Chang–Landau–Simmons-Duffin, "Spinning dispersive CFT sum rules and bulk scattering". It treats spinning functionals in the bulk-point regime only and explicitly leaves the Regge regime to future work.
* So the spinning, curved-regime, fixed-`ν`, absolute heavy Regge moment with exact-block verification appears not to be in the literature.
* Closely related prior work that should be cited:
  * **CHL** (Caron-Huot–Li, arXiv:2102.08160): the 3d helicity basis, mean-field OPE data diagonal in it for currents and stress tensors, and Lorentzian inversion for currents. This is directly relevant to the missing MFT normalisation and to M2.
  * Costa–Hansen–Penedones, arXiv:1707.07689: Regge limit of `⟨JJOO⟩` and `⟨TTOO⟩`, with the same spinning Regge pair in impact-parameter form.
  * Karateev–Kravchuk–Simmons-Duffin: weight-shifting operators (arXiv:1706.07813) for differential-basis shadow coefficients, and harmonic analysis (arXiv:1809.05111) for `b` with `Δ34`.

Novelty is incremental and technical. It is useful as a verified building block, not as a standalone result.
