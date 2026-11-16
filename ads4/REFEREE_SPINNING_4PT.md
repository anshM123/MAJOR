# Referee report: `SPINNING_4PT.md` (two spinning Regge pairs, `⟨JJJJ⟩`; `⟨TTTT⟩` partial)

Referee: independent, adversarial. Date: 16 Nov 2026. Object: `ads4/SPINNING_4PT.md`, `ads4/code/fourpt/` (not `fourpt/tttt/`).
Nothing under review was edited. Referee scripts and logs: `ads4/code/referee_4pt/` (a copy of `fourpt/` with its own caches, so the
authors' `out_*.txt` were not overwritten; `ref_analysis.py`, `ref_sflip.py`, `log_*.txt`).

## 1. What was re-run / reproduced

| item | how | result |
|---|---|---|
| `t1_structures_l1` (5 words, rank 5, transversality, 2-dim conserved space, exact word-wise shadow map) | rerun from scratch (caches deleted) | identical; also checked by hand: `v_1(AE)·b^{(0)}/b^{(2)} = −(Δ+J−2)(Δ−J−3)Y = v_1(AE)` at `(3−Δ,−1−J)` |
| `t1b_helicity_l1` (CHL `T^{E,h}` inside `span(v_1,v_2)`, closed forms) | rerun; CHL eq. `VVstcor` re-read from the arXiv source | identical; the CHL matrix in `chl_map.chl_vectors` is transcribed correctly (`n = (Δ−J−2)/2`, `J̃ = −1−J`) |
| `gen_l1` (25 word pairs, ≤ 2084 terms), `t2_validate_l1` | regenerated symbolic code from scratch (612 s), rerun | conservation `≤ 2·10⁻⁸` (control `0.7…1.6`), pair exchange `≤ 1.6·10⁻¹⁴`, `(H,H)` reduction `≤ 4·10⁻¹⁵`: reproduced |
| `t3_vertex_l1`, `t8_mft_l1` | rerun | identical outputs |
| block-tensor modes `t4` for (β≈0.4, `m=101.8`) and (β≈1, `m=39.3`) | regenerated independently into the referee cache | see §1.1 |
| table entries of §5 | recomputed from the regenerated modes | see §1.1 |
| own analysis on the authors' 10 cached cases | `ref_analysis.py`, `ref_sflip.py` | §3 |

### 1.1 Independent regeneration of two cases

Regenerated from scratch in `referee_4pt/fourpt/cache/` (symbolic set-up 12.6 min; `t4`: 124 s and 1458 s under load ≈ 16):
β≈0.4, `m = 12.4` `(14.166, 2)` and β≈1, `m = 39.3` `(45.861, 20)`. All 81 modes × 4 basis waves × ν = 0, 3 agree with the authors'
caches **bit for bit** (`ref_compare.py`), and the recomputed entries equal `out_t5_analysis_l1.txt`, e.g. β≈1, `m=39.3`, `k=6`, ν=0:
`W[v_1,v_1]` oo/ss `0.8282/1.0205`, `W[v_2,v_2]` `Σ_+/Σ_−` `0.9540/0.9614`; ν=1: `W[ĥ_s,ĥ_s]` ss `0.9607`; β≈0.4, `m=12.4`, ν=0:
`W[v_2,v_2]` ss `0.7050`, `W[ĥ_o,ĥ_o]` oo `−0.1237`. (The `m=101.8` case was started and abandoned for machine load; the large-`m` entries were
checked instead through the authors' caches with the referee's own analysis, §2(c).) This shows reproducibility of a deterministic
pipeline, not independence of method; the independent checks are the raw-mode / one-scale comparisons and the mutations below.
CHL closed forms re-checked at new exact points with different random kinematics: `(29/3, 5/2)` (half-integer `J`), `(−132/13, −8)`,
`(−62/7, −16/3)` (reflected): exact (`ref_chl_newpts.py`).

## 2. Claim-by-claim

**(a) Exact `⟨JJJJ⟩` partial waves.** Reproduced (table above). The validations are the right ones (conservation at three points with a
non-conserved control; pair exchange; scalar `(H,H)` reduction; exact mode symmetries). Conservation at `x_4` is not tested (covered by
the pair-exchange identity together with `x_2`; acceptable). Only parity-even exchanges and parity-even `⟨JJO⟩` structures are built; see issue M4.

**(b) Shadow map, helicity structures.** Correct. `m² = Δ²−3Δ−J²−J+2 = A_+²−A_−²` *is* the heavy-mass variable of the one-pair study
(`A_± = m(cosh, sinh)(β/2)`), so the CHL coefficients `f_opp = J(J+1)(m²+2)`, `f_same = −(Δ−1)(Δ−2)(m²−2)` are polynomials in the heavy mass
and are Weyl-even (checked). CHL's helicity-dependent factors were checked against the source: eq. `shadow T` (helicity basis) carries
`(2−Δ)_{|h_1+h_2|}/(Δ−1)_{|h_1+h_2|}`, eq. `hard shadow` (E/O basis) carries `diag((2−Δ)_{|J_1−J_2|}/(Δ−1)_{|J_1−J_2|}, (2−Δ)_{J_1+J_2}/(Δ−1)_{J_1+J_2})`
(= `1` and `(2−Δ)_2/(Δ−1)_2` for currents); the free data `free data` are `diag((−J)_2/(J+1)_2, (Δ−1)_2/(2−Δ)_2)` as quoted. The consistency check of
§2.3 (Euclidean shadow acting covariantly on `ĥ_h` reproduces CHL's ratio) works. But see m1/m2: the statement "no mixing, CHL eigenvalues
are normalisation artefacts" is basis bookkeeping, true by Hilbert 90 for *any* involutive `T`; the content is elsewhere (MFT normalisation).
The quoted formula is from eq. `shadow T`, not `hard shadow` (minor label).

**(c) End-to-end.** Confirmed, and more cleanly than the paper states it: the *raw* modes (no class fit at all) give the same ratios
(e.g. β≈0.4, `m=101.8`, ν=0: `W[ĥ_s,ĥ_s]` mode `(1,−1,−1,1)`: `0.99272`; `W[ĥ_o,ĥ_o]` mode `(1,−1,1,−1)`: `0.98374`), and a
*one-scale* fit (predicted pattern of all 81 modes, a single free scale) has uniform-weight residual `1.9·10⁻³` with scale `0.9932`
(`W[v_2,v_2]`, ν=0). So "no fitted constant" holds. The LSQ residuals `0.009…0.07` quoted in the table are artefacts of the
`1/|pred|` weighting (uniform-weight residual with 8 free classes: `0.9…1.6·10⁻³`); they should not be presented as a goodness-of-fit measure.

**(d) `±` structure / `S` test.** The test is genuine. Referee mutation (`ref_sflip.py`, `ref_analysis.py`): scaling `S → 0.5S, 2S` at
both pairs gives `R(Σ_−) = 3.98, 0.248` (i.e. `1/s²`, detected); `S = 0` raises the residual to `0.59`; `S_R → −S_R` makes the one-scale
residual `0.99` (vs `1.9·10⁻³`). With all 8 classes free the `±S_R` fits have *identical* residuals (same span) and differ only by the
relabelling `oo ↔ ss` — so the sign of the `S`-term is fixed only together with the helicity labels. Two caveats: (i) in `W[v_2,v_2]`
`T_oo/T_ss ≈ (A_−/A_+)⁴ ≈ 1.6·10⁻³` at β≈0.4, so `R(Σ_−) ≡ R(Σ_+) ≡ R_ss` to 4 digits — the agreement of `Σ_−` there is not separate
evidence; the independent information is `R_oo` (and the exact-helicity combos). (ii) Identification: classes are labelled in frames where pair
24 is referred by the inversion tensor, so `h_2 = h_1` is CHL-*opposite* (`J_z = ±2`) — this is correct and consistent with flat space
(all-incoming `c^{++}`: `J_z = 0`, all `J`; `c^{+−}`: `J ≥ 2`). With that dictionary the derived term is
`−g_1²S_LS_R(|c_same|² − |c_opp|²)`, the *opposite* sign to the conjecture `+g_1²S_{13}S_{24}(|c^{++}|²−|c^{+−}|²)` of
`SPINNING_E2E.md` §8 read with physical labels. The two agree only if E2E's `c^{++}` means `h_1 = h_2` *in the inversion-referred frame*,
which is CHL/flat *opposite* helicity. "Up to the sign convention of `S`" is not an explanation: a common sign flip of `S` leaves `S_LS_R`
invariant and `S` is ν-even in the `(±1,∓1)` modes (checked: `S = ∓(a−1)/8`); only a pair-24 helicity relabelling flips it. See M2.

**(e) MFT normalisation.** The relative helicity weight is right: I re-derived `Q_opp, Q_same` by hand from CHL `free data`, the
exact `a_h`, the factor `2^{(J+1)/2}n` of `VVstcor` (`n = α/2` cancels the `1/α` of `a_h`, so the non-Weyl-even `α²` drops out — not stated
in the text, should be) and the heavy `V_h`; both are Weyl-even, `→ m⁴/4`. But `t8` drops `K_c² 2^{J+1} c^E(Δ,J)^s` as "h-independent"
and then substitutes the scalar factor in the *paper's* normalisation (`P_FK(1,1)/b^{(0)}`); the conversion CHL-scalar→paper-scalar was never
shown to equal CHL-`VV`→paper-`VV` (frame factor `K_c`, `2^J`, `(−1)^J`, embedding `(−2)^{…}` conventions). A `J`-dependent mismatch
would not be an "overall constant" and could be non-Weyl-even. See M3. For a positivity-type bound a common *positive constant* is
harmless; an `m`- or `J`-dependent common factor is not (it changes the relative weight of different heavy states and the `1/m` claim).

**(f) No `O(1/m)`.** The [D] argument (reflection + Lemma S + Assumption T, Weyl-even rational change of basis) is correct as stated.
Fits: the quoted `|c_1| ≲ 0.04` (β≈0.4) is ansatz-dependent: with `(c_1,c_2,c_3)` on 4 masses the raw modes give `c_1 = +0.05…+0.19`
(`ref_analysis.py`) — still ≥ 20× below the mutations. Referee's own mutations: CHL-type non-Weyl-even helicity factor
`Δ(Δ−1)/((3−Δ)(2−Δ))` on the right same-helicity vertex: `c_1 = −3.97` (ν=0), `−3.92` (ν=3); `A_−² → J²` in `V_2(opp)`: `c_1 = +10.2`
in `W[v_2,v_2]` oo — detected. **The `c_1 = −2.6` "weak spot" is a class-fit artefact, not a property of the data:** the raw modes of
`W[ĥ_opp,ĥ_same]` (which at leading order contain only the pf/fp classes) give `c_1(A) = −0.009` at ν=0 (`−0.006` at ν=1). The paper's
explanation ("1/m² leakage amplified by ~70") is not shown and would produce a `1/m²`, not a `1/m`, pattern; the right fix is to quote the raw-mode fit.
So the dismissal is justified in conclusion, wrongly argued.

**(g) Stage 4 (Hilbert 90).** Correct: `r_k(Δ,J) r_k(3−Δ,−1−J) = 1` (checked: `u = (Δ+J)/2 → 1−u`, `w = (Δ−J−1)/2 → 1−w`), so
`Φ[c](Δ) = T(Δ̃)c(Δ̃)` is a σ-semilinear involution and Speiser's lemma applies; on a helicity line `κκ̃ = 1 ⇒ κ = f(Δ̃)/f(Δ)`. Caveats:
(i) the result is unavoidable bookkeeping — *every* such `κ` is an "artefact", so the proposition carries no physical information
about `O(1/m)`; all content is in (v): Weyl-evenness of the *physical* (positive) normalisation, which is not supplied by Hilbert 90 and for
`l = 2` is open; (ii) it needs helicity lines to be 1-dimensional in the parity-even conserved space (true for `J ≥ 2l` generic; low `J` and
parity-odd sectors not discussed); (iii) uses CHL's conformal invariance of helicity for conserved operators [Q].

**(h) `⟨TTTT⟩`.** Labels are mostly honest ([NE], "partial", "inconclusive"). The 14-word count with one relation agrees with the 13
parity-even `(2,2,J)` structures in 3d. Two over-statements: (i) the Richardson "dominant component" test is the ratio
`T(+,−,+,−)/T(+,−,−,+)`, i.e. of right-vertex parity partners whose leading values are equal by construction — a weak probe; say so;
(ii) §7.4/§8 "Established: … the `±` algebra … for `l = 2`" — this follows only if the Weyl-covariant families can be chosen helicity-pure
and their MFT normalisation is Weyl-even, both [ND]. Summary item 7 and STATUS should read "ingredients checked at three points; `±`
algebra and no-`O(1/m)` for stress tensors conjectural".

## 3. Numbered issues

**Fatal:** none found.

**Major**
- **M1 (label / framing).** "Shadow matrix = (1/b)×Weyl reflection, no mixing" and "CHL eigenvalues are normalisation artefacts" are
  presented as results; by Hilbert 90 they hold for any basis choice after renormalisation and carry no physics. Fix: state the invariant
  content (word-wise `b^{(Δ_34)}` mechanism verified; `T(Δ̃)T(Δ) = 1`; helicity commutes with `T`) and move the weight to §4.4.
- **M2 (sign / dictionary).** The `±` formula is derived and verified in the paper's frame labels, but its relation to the
  `SPINNING_E2E.md` §8 conjecture and to the flat-space `B_3` term certified in `FLAT_GRAVITON_BOUND.md` is waved through ("up to the sign
  convention of `S`"). Fix: state the result as `R_LR_R(|c_same|²+|c_opp|²) − g_1²S_LS_R(|c_same|²−|c_opp|²)` with `S_R` defined by
  `M(−ν)` in the inversion-referred frame, and derive explicitly (flat limit of `S_LS_R`) which flat-space sign this is; correct the E2E §8 formula if needed.
- **M3 (MFT normalisation).** "[D from Q], MFT normalisation Weyl-even" covers only the *relative* helicity weight. The `h`-independent
  factor `K_c² 2^{J+1} c^E_s` (CHL conventions) vs `P_FK(1,1)/b^{(0)}` (paper conventions) is asserted, not mapped. Fix: compute the
  CHL→paper conversion of the full `⟨VVO⟩` structure (including the `(x·ε_3)^J` prefactor) and of `T_00`, or check one OPE limit; until then
  mark the absolute/`J`-dependent normalisation [I] and the "no `O(1/m)` from normalisation" statement as covering the helicity ratio only.
- **M4 (scope).** Only parity-even exchanged operators. Physical `⟨JJJJ⟩` (and `⟨TTTT⟩`) also exchanges parity-odd operators with two
  parity-odd `⟨JJO⟩` structures (CHL `T^{O,h}`, MFT `c^O`). The `±` structure, the end-to-end match and the `O(1/m)` verdict are untested
  there. Fix: say so in §0 and STATUS; "settled for currents" → "settled for parity-even exchanges".

**Minor**
- m1. §2.3 quotes `(2−Δ)_{|h_1+h_2|}/(Δ−1)_{|h_1+h_2|}` as eq. `hard shadow`; it is eq. `shadow T`.
- m2. Status row 3: "exact at 10 points incl. 4 at reflected labels"; `out_t1b_helicity_l1.txt` shows 6 checked points (2 reflected, 1 skipped);
  the six fitting points are not in the code. Either add the fit script or say "6 points".
- m3. LSQ residuals in the §5 table are weighting artefacts (uniform: `~10⁻³`); report the one-scale (zero-parameter) residual instead.
- m4. `R(Σ_−)` in `W[v_2,v_2]` is numerically identical to `R_ss` (`T_oo/T_ss ≈ 10⁻³`); do not cite it as an independent `S` test.
- m5. `|c_1| ≲ 0.04/0.2` depends on the ansatz (`c_3/m³` included: up to `0.19`); quote the range over ansätze.
- m6. §6.2 weak spot (i): replace the leakage explanation by the raw-mode fit (`c_1 = −0.009`).
- m7. §4.4: state that the `n = α/2` prefactor of `VVstcor` cancels the non-Weyl-even `1/α` of `a_h`.
- m8. TTTT Richardson test: name the ratio used (`T(+−+−)/T(+−−+)`), and weaken summary item 7.

## 4. Verdicts

* **Analytic correctness:** sound (word-wise shadow map, exact CHL identification, heavy vertex, `±` algebra, reflection argument, Hilbert 90).
  Gaps: sign dictionary to flat space (M2), absolute/`J`-dependent MFT normalisation (M3), parity-odd exchanges (M4).
* **Numerical soundness:** good. Validations and two regenerated cases reproduce; raw modes and a one-scale fit confirm the factorised
  formula with no free constant; `S`-term magnitude tested to ~1% by mutation; no `O(1/m)` in raw modes (the −2.6 is a fit artefact).
* **Label honesty:** mostly honest; over-statements M1, m2, m4, m5, §7.4/item 7.
* **Physical applicability:** gives the leading two-pair heavy Regge moments of `⟨JJJJ⟩` for parity-even exchanges at `k = 4, 6`,
  `ν ≪ m`, incl. the opposite-sign term, and removes the `κ_h` objection for currents. It does **not** give: `⟨TTTT⟩` end-to-end, `k = 2, 3`
  (the moments the graviton sum rules use), `ν ~ m`, parity-odd exchanges, the light side, the flat-space sign dictionary (M2), the
  absolute normalisation (M3). The overall *constant* is harmless for a positivity bound; a `J`- or `m`-dependent common factor is not.
* **Novelty:** incremental. The leading heavy-state formula is the CLS-D bulk-point/"scattering-crystal" dictionary (arXiv:2311.04271,
  which already treats `⟨JJJJ⟩`/photons and `⟨TTTT⟩`/gravitons at leading order) specialised to Regge moments; structures, shadow
  coefficients and MFT data are CHL (2102.08160) and KKSD (1809.05111). New: exact end-to-end numerics with both pairs spinning, the
  explicit `(v_1,v_2) ↔` CHL map, and the absence of `O(1/m)` for currents. A search (2023–2026) found no paper doing this.
