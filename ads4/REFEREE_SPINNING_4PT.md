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

---

# Round 2 (16 Nov 2026, after commit `5df911a`)

Scope: M2 (new §4.5), M3 (new §4.6), text of M1, M4 and the minor items. Referee scripts: `ads4/code/referee_4pt/ref_sign_frames.py`,
`ref_sign_moments.py`; copies of `fourpt/tttt/{t1_mft_l1,t2_mft_heavy_l1}.py` run in `referee_4pt/fourpt/tttt/` (fresh caches; nothing tracked was touched).

## M2: sign of the opposite-sign term. Resolved.

1. **Frame bookkeeping, checked with the actual vectors.** I took the positions (`e2e/moments.frame_positions`) and polarisations
   (`moments4.pol_grid4`) exactly as they are fed to the blocks, at random grid points and angles. I evaluated `D_5 = det[z_in, z_out, n̂]`,
   `n̂ = (x_in − x_out)/|·|` (past-pointing timelike at both pairs; in = 1, 2 lie in the past), and `D_1 = z_in·z_out + (z_in·n̂)(z_out·n̂)` with
   no frame conventions. Results: `D_1 = cos(θ_out − θ_in)` at both pairs, `D_5 = −sin(θ_out − θ_in)` at **both** pairs, so
   `D_5^{13}D_5^{24} = +sin·sin` (ratio `+1.000000`, 6 random points). Any normalisation of `ε` multiplies the product by `c² > 0`, so
   this is convention-free. The authors' "two orientation flips cancel (`det I = −1`, `I x̂ = −x̂`)" is therefore correct.
2. **Kernel action, recomputed with my own weights.** I derived the mode weights numerically on a 64-point circle grid instead of copying
   `w_odd = ∓i/2`, and applied them to the exact cached modes. `K_odd/K_even` = `+0.96122 / −0.96127` (same/opp, β≈0.4, `m=101.8`,
   ν=0) and `+0.96126 / −0.96127` (β≈1). These agree with `out_t4_sign_l1.txt` to all digits, with `Im/Re < 10⁻¹⁶`. Controls:
   * `W[v_1,v_1]` gives `+0.0040, +0.0080` at `m = 102, 72`. This is `∝ 1/m²`, as predicted, since `V_1` weighs both helicities equally.
   * `W[v_2,v_2]` gives `0.958` (β≈0.4) and `0.878` (β≈1), equal to `0.961·(A_+⁴−A_−⁴)/(A_+⁴+A_−⁴)` at both β.
3. **Labels and crossing.** "same" is CHL `T^{++}+T^{−−}`. Its defining property is `J_z = h_1 − h_2 = 0` for both momenta incoming, and it
   exists for all `J`. The Fock computation independently confirms that `ĥ_opp` is absent for `J = 0, 1`. This matches the flat channel
   with `D_{++} = P_J = d^J_{00}` and `σ_{++} = +1` in `FLAT_GRAVITON_BOUND.md` (0.1). The identification goes through `J_z` and the
   Wigner function, so it does not depend on whether helicities are counted all-incoming or in/out.
4. **Attempts to break it.**
   * Using `x_out − x_in`, swapping the arguments of `ε`, or flipping its orientation, each applied at both pairs, leaves the product unchanged.
   * Only a rule that differs between the two pairs flips the sign, and that is not a covariant kernel. It is precisely the inversion-referred
     `S_R` of §4.3.
   * The physically invariant content: under one covariant rule, the odd⊗odd kernel acts with opposite signs on the two helicity channels,
     `K_odd ≈ +0.96 K_even` on `J_z = 0` and `≈ −0.96 K_even` on `|J_z| = 2`, with `J_z` as in item 3. The overall sign of an AdS odd kernel
     is fixed when the flat functional is uplifted, so positivity transfers with `σ_h` exactly as in flat space.
   * Still open, and correctly labelled [I] by the authors: the explicit uplift map (which flat sum rule ↔ which AdS kernel combination, with
     positive coefficient), and the light side.
5. **New error (minor, R2-1).** §4.5 states "`K_even ± K_odd ≥ 0` per helicity, as unitarity requires". This is false beyond ν = 0:
   * on the opposite-helicity wave `K_odd/K_even = −1.00477` (β≈0.4, ν=3) and `−1.01059` (β≈1, ν=3), so `K_even + K_odd < 0`;
   * at β≈1, ν=3, `K_even` itself is negative (`P_{−½+3i}(cosh 1) < 0`).

   The ν-harmonic kernels are not pointwise positive, and unitarity says nothing here. Delete the clause, or restrict it to the ν→0
   leading value and drop "as unitarity requires".

## M3: MFT normalisation. Resolved.

* **Fock data reproduced.** `t1_mft_l1` was rerun in the referee copy with fresh caches: the same exact rationals, `G_offdiag = 0` on every
  `n ≥ 1` primary, and the `n = 0` degeneracy. This reproduces their code; it is not an independent Fock implementation.
* **Independent check of the closed forms.** I used the Fitzpatrick–Kaplan generalised-free-field formula (written by hand, not from their
  code), *without* the Bose factor `1+(−1)^J`:
  * `p[ĥ_opp]` and `p[ĥ_same]` both give `P_FK = 3/245` at `(Δ,J) = (6,2)`;
  * `p[ĥ_opp]` gives `P_FK = 6.184·10⁻⁴` at `(8,2)` and `3.527·10⁻³` at `(7,3)`;
  * at `n = 0`, `(4,2)`: `p = 1/900 = p_opp + p_same` (closed forms).

  All agree, so the convention is distinguishable currents. For identical currents the even-`J` data acquire the factor 2; this should be
  stated.
* **No circularity.** CHL enters nowhere: the inputs are `mft/fock.py` Gram forms, the conserved basis and `f_h`, which were identified
  from CHL in round 1 but are now *tested* by the exact diagonality of the Fock Gram matrix, and FK.
* **Absolute constant.** I checked analytically that `t2_mft_heavy_l1` gives `Q_h → 1`, using `P_FK/b^{(0)} = (8/π²)A_+A_−`,
  `2J+1 = 2A_−`, `2Δ−3 = 2A_+`, `a_{2,2} = π`. The rerun reproduces `m²(Q−1) = 66.30/7.03` (β=0.4). So `1/4` per helicity quadruple and
  `C^± = (2P/m^{2k})(|c_same|² ± |c_opp|²)` follow, with Weyl-even factors.
* **Remaining caveat (minor, R2-2).** §4.6 is headed [D], but the closed forms are identified from 25 primaries (`N ≤ 6`) and then used at
  `n ~ 50`, so they are [ID]. The two-pair lightcone block test (`1+O(ε)`, control 12.6) covers only the `n = 0` family, where
  `ĥ_opp = ĥ_same`. The helicity-*split* absolute normalisation is therefore block-tested only through the shared convention chain. Fix:
  label the closed forms [ID] (or derive them), and if feasible test one `n = 1` family at block level (next lightcone order).

## M1, M4, minors

* M1 is fixed (§0 item 2, §0 item 6, §8).
* M4 is fixed in the text (§0 item 5, §6.3).
* m1–m8 are applied as listed in §11; spot-checked m1, m2 and m6.
* **R2-3 (minor):** the `STATUS.md` row 6 is stale. It still says "shadow matrix exactly (1/b)×Weyl reflection; no `O(1/m)` for currents;
  overall two-pair MFT constant not fixed; referee running". It should say: parity-even exchanges only; Hilbert-90 bookkeeping; constant
  `1/4` fixed by Fock data; refereed (two rounds).

## Round-2 verdicts

No fatal and no major issues remain. Minor: R2-1 (false unitarity clause), R2-2 ([ID] label; `n ≥ 1` block test), R2-3 (`STATUS.md`).

* **Analytic correctness:** sound. The sign is resolved and the normalisation is derived in own conventions. Open, and labelled: the uplift
  map and the parity-odd sector.
* **Numerical soundness:** good. Everything checked reproduces, and the sign is confirmed from the raw vectors and modes with independent weights.
* **Label honesty:** good after R2-1 and R2-2.
* **Physical applicability:** unchanged from round 1, except that the `±` sign and the absolute MFT constant are now established for
  parity-even `⟨JJJJ⟩`. Not given: `⟨TTTT⟩`, `k = 2, 3`, `ν ~ m`, parity-odd exchanges, the light side.
* **Novelty:** incremental, as in round 1.
