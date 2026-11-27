# m1_direct: the action of Ω^EE on light blocks by an independent route (referee point M1, decisive check)

27 Nov 2026. Unrefereed. Nothing outside `ads4/code/m1_direct/` was edited; nothing was committed.
Read-only imports: `pcyc/pcyc_dlib.py`, `pcyc/pcyc_dquad.py` (its `Rule` class only; its `run()`, which writes to `pcyc/cache`, is never
called), `rd/rd_core.py`, `spinfunc/l4_tto_even.py` (`monomials`, `conserved_even`; they write nothing), `m1_remainder/m1_weight.py`,
`fullpos/fp_poly.py` (one cross-check of Q_k). Data read: `prl/data/k2/k2_continuum_D_{300,1000}_rep.json`, `prl/data/nf02_psi.json`,
`no_eft/cache/nl1_*.txt`, `pcyc/cache/dq_D20_S_*.npz` (one reproduction check).

Grades: [P] proved, [D] derived, [DV] derived and checked numerically, [NE] numerical evidence, [E] estimate, [I] open.

## 0. Bottom line

**(a) Direct versus ladder.** The direct route is the literal definition of the functional: CLS-D's per-block integral
(Lorentzian shadow representation of the block, causal domain 4>5>3>0>1>2 = the real chain D′ of `PCYC_CLOSE.md`), integrated
per ν by a 7-D product rule, exact in Δ, with no light-ray matrix elements, no CDH polynomials and no rest-frame reduction. It is
compared with the per-ν ladder polynomial of (fp-H), normalisation-free (ratios of states at the same ν).

* **k = 2 piece (Ĉ′_EE), J = 0** [DV]. Agreement to 2.3·10⁻⁵ (Δ = 20, T = 342), 7.5·10⁻⁶ (Δ = 30), 4.6·10⁻⁷ (Δ = 100) on ν ≤ 6, against
  2·10⁻², 8.5·10⁻³, 6.4·10⁻⁴ for the control without the ladder's 1/T term. The fitted coefficient of the ladder's 1/T term is
  κ = 1.00016 / 1.00002 / 1.00013 (ladder: 1). At Δ = 10 (T = 72) the agreement is 4·10⁻³, limited by quadrature (control 10⁻¹).
  O(1/T²) terms of natural size are excluded at Δ = 20 by a factor ≈ 20.
* **k = 2 piece, spin J** [NE]. J = 2: 4.5·10⁻⁴ at Δ = 30 (control 10⁻²), 5.6·10⁻³ at Δ = 16 (control 5.6·10⁻²). J = 8 at Δ = 30
  (twist 22, both ladders b = 0, 4; polynomials of degree 11 and 7 in ν²): 6–9·10⁻⁴ per structure pair, against 10⁻²–2.4·10⁻¹ for
  the leading-Regge control. J = 4: consistent but not discriminating at the grids run.
* **k = 3 piece (Ĉ₃, j′ = 1, o_H q)** [NE, weaker]. J = 0, where the ladder image is exactly ν-independent: the direct ratio to a
  heavy reference is ν-independent to 1.8·10⁻⁴ at Δ = 20 (5.8·10⁻³ → 1.5·10⁻³ → 1.8·10⁻⁴ under grid refinement), i.e. any 1/T
  term is < 1 % of the size of the k = 2 one. Higher spin at k = 3 is not yet discriminating (§4.4).
* **The per-state numbers are reproduced independently** [DV]. Quadrature of the *published* continuum weight against my own
  ladder polynomials (no certified moment balls) gives log₁₀ĥ to all printed digits of `no_eft`: e.g. Δ_g = 1000, twist 2: J = 8
  +1.918·10⁴² (b = 0), J = 20 +2.475·10⁸⁴ (b = 0; the brief's "10⁸⁰" is 10^{84.4}); Δ_g = 300: J = 8 6.92·10³², J = 20 2.14·10⁶².
  The k = 2 and k = 3 parts have opposite signs and are of the same size (ĥ⁽³⁾ ≈ −2.4 ĥ⁽²⁾ at twist 2), so both pieces matter.
* **The defining integral does not converge below twist 5** [D + NE]. On a block of twist τ = Δ − J the D′ integrand behaves like
  B^{(τ−5)/2 − 1} dB at the light cones of the exchanged operator (B → 0), for k = 2 and k = 3 and every J tested (decay rate
  1 − (5+J)/Δ in the pcyc variable t, measured to 4 digits at J = 0, 2, 4, 8). So CLS-D's per-block integral converges absolutely
  iff **τ > 5**. For the twist-2 and twist-4 double traces, which carry the 10⁶⁸¹, the action of Ω^EE is defined only by analytic
  continuation in Δ. The ladder is analytic for τ > 1 and equals the direct value on τ > 5, so it *is* that continuation.
* **Small gap** [DV]. A toy functional of reach ν ≤ 6 (Δ_g ≈ 8) gives smeared direct/ladder differences −4·10⁻³ (Δ = 10, 12-grid),
  −1·10⁻⁵ (20), −3·10⁻⁶ (30), against −8·10⁻², −2·10⁻², −7·10⁻³ for the leading-Regge control. With such a gap nothing is large;
  the 10⁴²–10⁶⁹¹ at Δ_g = 1000 is the same exact polynomial evaluated where the weight reaches, ν ≈ Λ + 100.

**(b) Verdict.** The 10⁶⁸¹ is **not an artefact of continuing (fp-target) outside its domain**. The ladder is the exact per-ν image
of the functional wherever the functional's own integral converges, light states included, and the unique continuation below.
The per-state light actions (10⁴², 10^{84}, …, 10⁶⁹¹) are properties of Ω^EE as defined: a light (r → 1) state's ladder is
undamped up to its end, and the weight reaches ν ≈ Λ + 100. What remains of M1 is the other alternative of `m1_remainder` §5:
the step that writes the light side of the finite-c_T sum rule as the block-wise sum R = Σ_{𝓛′} p_O ω[O] (H_K applied block by
block to the light tower). For τ ≤ 5 that step has no convergent block integral to rest on (above), and `m1_remainder` §5 shows
that the block-wise R is incompatible with the sum rule at any realistic c_T. The present check does not decide whether the
block-wise light side fails, or whether Einstein-like light data violate an input; it removes the continuation as the explanation.

**(c) Honest size of R.** No sound method in hand gives a remainder smaller than the term-wise one. With the block-wise light side,
|c_T²R| ≈ 10^{679}–10^{682} and c_T* ≈ 10^{680}–10^{683} at Δ_g = 1000 (`m1_remainder`, now with the per-state input confirmed).
Without it, R is not controlled at all at finite c_T. Either way Eq. (5) is a statement about c_T → ∞ at fixed Δ_g, and even that
limit uses the block-wise light side (Theorem A), which for τ ≤ 5 is a continuation, not a convergent sum of integrals (§5).

**(d)** Replacement text: §7.

## 1. Ledger

| step | what | file / output | grade |
|---|---|---|---|
| 1 | read m1_remainder, no_eft/NOTE, Supplement (Order of limits, L0, H_K, defs), sec_fullpos, xcheck_ladder, PCYC_CLOSE, RD_CLOSE, EESTEP_CLOSE, REPRO_ENGINE | — | — |
| 2 | how a "direct row" is computed; can it run below the gap (§2) | — | [D] |
| 3 | own ladder polynomials (fp-zy, fp-U, fp-H; k = 3 via the d = 1 CDH recurrence and ψ₃) | `md_ladder.py`, `md_k3.py` | [DV] (Q_k = fp_poly to 10⁻¹⁶) |
| 4 | per-state ĥ by ν-quadrature of the published weight, k = 2 + k = 3; vs no_eft | `md_hhat_quad.py`, `md_hhat_total.py` → `out/hhat_total.txt` | [DV] |
| 5 | direct D′ integrand, k = 2, J = 0 (Φ_EE2 = Q₂Φ_scal, Q₂ = Q_k4/R) and its boundary exponents | `md_dquad.py`, `md_probe.py` | [D + DV] |
| 6 | spin-J integrand, k = 2 and k = 3 (C3H), own contractions; unit tests | `md_spinJ.py`, `md_test_spinJ.py`, `md_probe3.py` | [DV] (J = 0 vs rd_core to 12 digits; contraction identity to 10⁻¹⁵) |
| 7 | B-zero decay rate = 1 − (5+J)/Δ ⇒ convergence iff τ > 5 | `md_decay.py` → `out/decay_*.txt`, `out/log_decay.txt` | [D + NE] |
| 8 | J = 0, k = 2: ladder vs direct, κ fit, grid convergence | `md_compare_J0.py`, `md_kappa.py`, `md_conv.py` → `out/compare_J0.txt` | [DV] |
| 9 | spin J, k = 2 | `md_dquadJ.py`, `md_compare_J.py` → `out/compare_J.txt` | [NE] |
| 10 | k = 3, J = 0, 2 | `MD_KIND=C3 md_dquadJ.py`, `md_compare_k3.py` → `out/compare_k3.txt` | [NE] |
| 11 | toy small-gap functional | `md_smallgap.py` → `out/smallgap.txt` | [DV] |
| 12 | verdict, honest R, replacement text | §5–7 | — |

## 2. What a "direct row" is, and why it cannot be used below the gap

* **The direct rows of Secs. pos–allgap** evaluate, per state, the beyond-leading-Regge correction Γ(ν) of the per-ν image at nodes
  x = ν/m < 1 (glossary: "node", "row", "K₂ row", "direct row"). Each row is a QMC or saddle (Laplace) evaluation of the CLS-D
  integral in variables scaled by a = 1/Δ around the crystal saddle, as a ratio to a much heavier reference state at the same ν.
* **They fail at the light-cone edge.** Direct quadrature fails as x → 1 (Table tab:computed: "x = 0.95 failed"; H_edge), and the
  part of the weight with ν > m is never sampled.
* **For a light state m is small.** Twist 2 has m = 0, so every ν of the weight sits at x = ∞; twist 4, J = 8 has m ≈ 6. A direct
  row of the published kind does not exist for any light state that matters for M1.
* **xcheck_ladder's 116/116 agreement** tests (fp-H) only at heavy states and only on ν < m. It cannot see terms that are
  exponentially small there.
* **What is valid below the gap** is the integral that the rows approximate: CLS-D's per-block integral over the real chain D′,
  exact in Δ. `PCYC_CLOSE` §1 showed that D′ *is* the physical chain, and checked its scalar version against the exact CHMRS Mellin
  skeleton at Δ = 20, 30, 50 to 10⁻⁸. That is the independent route used here, at light Δ and small ν.

## 3. Method

### 3.1 Integrand [D + DV]
* **Scalar part.** Φ_scal is the scalar k = 2, J = 0 CLS-D integrand on D′ in the pcyc coordinates (t_x, r_x, s, t_y, r_y, φ, α),
  where B^{Δ/2} = e^{−t_x−t_y} exactly (`pcyc_dlib.point`, read-only).
* **k = 2 EE (Ĉ′_EE: j′ = 0, f = v′−u′, dual EE).**
  Φ = Φ_scal · u′v′ · prefĥ_J /(m₁₂² m₃₄²) · |m₁₂m₃₄/(m₁₅m₂₅m₃₅m₄₅)|^{J/2}.
  * At J = 0 this is rd_core's Q_k4 divided by R = ((v′−u′)²−1)/(4u′v′). Q_k4 is checked [P + DV] against repro_engine.
  * At spin J the prefactor is ⟨F_a, Ê², G_c⟩: the conserved parity-even ⟨TTO_Δ⟩ (a) and ⟨TTÕ_{3−Δ}⟩ (c) numerators, with
    spin-2 contractions at x₂, x₄ and the spin-J contraction at x₅ = e.
  * Contractions are by null-circle Fourier modes with κ_n = (−1)^{l+n}2^{−l}C(2l, l+n). The identity ⟨T, (w′·)^l⟩ = T(w′) is
    checked to 10⁻¹⁵ at l = 2, 4, 8.
  * The J-dependent m-powers are those of `w0_ee_common._gs`. At J = 0 the spin code reproduces Q₂ up to a constant (12 digits).
* **k = 3 (Ĉ₃: j′ = 1, f = −1, dual o_H q).**
  Φ = Φ_scal · [u′v′/(v′−u′)] · prefĥ₃,J/(m₁₂²m₃₄²) · |…|^{J/2} · sign det[z₀; x₁₀; x₃₀].
  * This is `l9_points.build_spin`, kind 'C3H', re-implemented point-wise. o_H uses the 5×5 embedding ε, and the j′ = 1 structure
    is the j′ = 0 one times that locally constant sign (H0_SPINNING_FUNCTIONAL §21.6).
  * X₂₀, X₄₀, X₂₄ < 0 on all of D′ (2·10⁴ samples), so the √ factors carry only a constant phase. Both pieces are invariant under
    the transverse reflection used to halve the α range (checked).
* **ν-dependence.** A^{iν/2} is handled exactly as in pcyc_dquad: atom powers in the 1-D weights, exp(iν ε/2) at the nodes.

### 3.2 Quadrature
* **Rule.** Gauss–Laguerre in t_x, t_y with weight t^p e^{−λt}, λ = 1 − (5+J)/Δ (the measured decay); Gauss–Jacobi in r_x, r_y, s,
  φ, α with the measured boundary exponents (`md_probe.py`, `md_probe3.py`).
* **Grids.** 8⁴·6·12·8 up to 12⁴·8·20·12, plus a per-dimension scan at Δ = 10. r_y dominates the error at small Δ.
* **Reproduction.** My wrapper reproduces pcyc's cached scalar D′ integral at Δ = 20 exactly.
* **Comparison.** All comparisons are grid-matched ratios to a heavy J = 0 reference (Δ = 300) at the same ν, with one free complex
  constant per state (or the 2-dim span {S₀, S₄} for J ≥ 4). This cancels every ν-only factor (N_com, A_{k,ν}, CLS-D constants) and
  every Δ-only factor (B|a|², sin², shadow coefficient). It therefore tests exactly the ν-profile at fixed state, i.e. the content
  of (fp-H) beyond the state normalisation.

## 4. Results

### 4.1 Per-state numbers, independent of the moment balls [DV] (`out/hhat_total.txt`)
ĥ = ĥ⁽²⁾ (ν-quadrature of the published weight × my S_b) + ĥ⁽³⁾ (my Q_k). In brackets: no_eft (moment balls).

| Δ_g | state | ĥ⁽²⁾ | ĥ⁽³⁾ | ĥ | log₁₀\|ĥ\| [no_eft] |
|---|---|---|---|---|---|
| 1000 | twist 2, J 8, b 0 | −1.31·10⁴² | +3.23·10⁴² | +1.918·10⁴² | 42.283 [+42.283] |
| 1000 | twist 2, J 8, b 4 | −4.56·10²⁵ | −1.51·10²⁶ | −1.969·10²⁶ | 26.294 [−26.294] |
| 1000 | twist 2, J 20, b 0 | −2.01·10⁸⁴ | +4.48·10⁸⁴ | +2.475·10⁸⁴ | 84.394 |
| 1000 | twist 2, J 20, b 4 | −1.31·10⁷¹ | −2.74·10⁷¹ | −4.05·10⁷¹ | 71.607 |
| 1000 | twist 3, J 8, b 0 | | | +1.128·10⁴¹ | 41.052 [+41.052] |
| 1000 | twist 4, J 8 / 20, b 0 | | | +1.25·10⁴⁰ / +2.87·10⁸¹ | 40.098 [+40.098] / 81.458 |
| 1000 | twist 998, J 8 (heavy-side control) | 2.08·10⁵ | −5.35·10⁴ | 1.543·10⁵ | 5.188 [+5.188] |
| 300 | twist 2, J 8 / 20, b 0 | | | +6.92·10³² / +2.14·10⁶² | 32.840 [+32.840] / 62.331 |

### 4.2 k = 2, J = 0: direct versus ladder [DV] (`out/compare_J0.txt`, `out/kappa`, reference Δ = 300, ν ≤ 6)
The ladder image is P = z₂(s)/3 + 2y₂(s)/T exactly (z₂ = (s²−46s+72)/12, y₂ = −(s³−82s²+840s−1008)/144). LR = the same without
the 1/T term. The residual is max_ν |I₁P₂ − cI₂P₁|/rms.

| Δ (T) | grid | ladder | LR control | κ (ladder: 1) |
|---|---|---|---|---|
| 10 (72) | 10·16·10·10·24·16·10 (tuned) | 4.1·10⁻³ | 9.8·10⁻² | 1.028 (12-grid: 1.049; 10-grid: 0.990; quadrature-limited) |
| 12 (110) | 10⁴·6·16·10 | 3.7·10⁻³ | 6.1·10⁻² | 0.986 |
| 16 (210) | 10⁴·6·16·10 | 4.4·10⁻⁴ | 3.3·10⁻² | 0.9993 |
| 20 (342) | 12⁴·8·20·12 | 2.3·10⁻⁵ | 2.0·10⁻² | 1.00016 |
| 30 (812) | 12⁴·8·20·12 | 7.5·10⁻⁶ | 8.5·10⁻³ | 1.00002 |
| 100 (9702) | 12⁴·8·20·12 | 4.6·10⁻⁷ | 6.4·10⁻⁴ | 1.00013 |
| 7, 8 | up to 12 | not resolved: τ = 7, 8 is close to the convergence threshold τ = 5, where the integrand decays like B^{(τ−5)/2−1} | | |

Grid convergence (`md_conv.py`): 1.3·10⁻³ between the 10- and 12-grids at Δ = 30, 100; 1.7·10⁻² at Δ = 10.

### 4.3 k = 2, spin J [NE] (`out/compare_J.txt`; residual per (a, c) entry; ν ≤ 6)
The (a, c) matrix of structure pairs must lie in span{S₀(ν), S₄(ν)} (S₀ alone for J = 2). r = C/T is the ladder's impact
parameter variable; LR = the same ladder with U_{b,k} → r^k and the 1/T terms dropped.

| J | Δ (twist, r) | grid | ladder residual per entry | LR control |
|---|---|---|---|---|
| 2 | 30 (28, 0.008) | 8⁴·6·12·8 | 4.5·10⁻⁴ | 9.9·10⁻³ |
| 2 | 16 (14, 0.029) | 8⁴·6·12·8 | 5.6·10⁻³ | 5.6·10⁻² |
| 2 | 10, 12 | 8⁴·6·12·8 | not converged (0.1–2) | same |
| 8 | 30 (22, 0.089) | 10⁴·6·16·10 | 7.3, 9.3, 7.2, 5.7·10⁻⁴ (3.0–4.2·10⁻⁴ on ν ≤ 4) | 1.4·10⁻², 2.4·10⁻¹, 3.7·10⁻², 1.0·10⁻¹ |
| 8 | 30 | 8⁴·6·12·8 | 2.1, 2.5, 2.1, 1.8·10⁻² | 3.3·10⁻², 0.23, 4.4·10⁻², 0.11 |
| 8 | 20 (12, 0.21) | 8⁴·6·12·8 | 0.07–0.17 | 0.07–1.6 (not discriminating at this grid) |
| 4 | 30, 14 | 8⁴·6·12·8 | 2·10⁻⁴–3·10⁻³; 3·10⁻³–0.1 | similar: the 2-dim span is too flexible at this grid |

The J = 8, Δ = 30 state has twist 22 ≪ Δ_g. Its residual falls by 30× from the 8- to the 10-grid while the control's stays at
10⁻²–10⁻¹, so the J = 8 ladder (U_{b,k} for b = 0, 4, k ≤ 8; S₀ and S₄ of degree 11 and 7 in s) is the direct image to ≲ 10⁻³.

### 4.4 k = 3 [NE] (`out/compare_k3.txt`)
The ladder's k = 3 image is S3_b(s) = σ_b Σ_k D_k² U_{b,k} q_k(s). At J = 0 it is exactly q₀ = 1, so the direct ratio to the J = 0
reference must be ν-independent at every Δ. (The LR control coincides with the ladder at J = 0.) Below, the residual is
max_ν |ratio/fit − 1| on ν ≤ 6.

| Δ (T) | 8⁴·6·12·8 | 10⁴·6·16·10 | 12⁴·8·20·12 |
|---|---|---|---|
| 20 (342) | 5.8·10⁻³ | 1.5·10⁻³ | 1.8·10⁻⁴ |
| 30 (812) | 2.4·10⁻³ | 6.5·10⁻⁴ | — |
| 100 (9702) | 1.8·10⁻⁴ | — | — |

* **The remainder is quadrature error.** It scales like 1/T at fixed grid and has the same ν-shape at Δ = 20, 30, 100. That is the
  signature of quadrature error in the integrand's Δ-dependent part, which integrates to a ν-independent constant. It falls by 30×
  from the 8- to the 12-grid.
* **Size.** At the 12-grid, T × residual = 0.06. A genuine 1/T term of the size of the k = 2 one (c(ν) = 6y₂/z₂ = 2.6 to 8 on
  1 ≤ ν ≤ 8) would show as ≳ 10⁻², 50× larger.
* **Spin.** J = 2 (Δ = 16, 30) was run only on the 8-grid: residual 1.2·10⁻², 2.5·10⁻³, not yet discriminating from its LR control.
* **Grade.** The k = 3 piece is therefore [NE] at J = 0 and untested at higher spin. The giant light values have
  ĥ⁽³⁾ ≈ −2.4 ĥ⁽²⁾ (§4.1), so the k = 3 ladder deserves the same J ≥ 2 test as k = 2.
* **Is it the paper's functional?** Ĉ₃'s O(1/m²) structure depends on the o_H/o_V design choice (H0_SPINNING_FUNCTIONAL §21.6). The
  integrand used here is the project's C3H (o_H), the one the direct rows used.

### 4.5 The convergence threshold τ > 5 [D + NE] (`out/log_decay.txt`)
* **The exponent.** In the pcyc variable t (B = e^{−2t/Δ}), d log|Φ|/dt at large t equals −(1 − (5+J)/Δ) to 3–4 digits: J = 0, 2, 4, 8
  at Δ = 20, and J = 0, 8 at Δ = 40, six base points each, k = 2. The same holds for k = 3 at J = 0, 2.
* **The divisor count.** This is the exponent of RD_CLOSE's divisor table at m₃₅ = 0 or m₄₅ = 0: scalar (Δ−3)/2, plus −2 from the
  ⟨TTO⟩ polarisation structures, plus −J/2 from the spin-J three-point denominators. Hence |Φ| ∝ B^{(τ−5)/2−1}dB, integrable iff τ > 5.
* **Which states it excludes.** It excludes exactly the twist-2 and twist-4 double traces [TT]_{0,J}, [TT]_{1,J}, at γ = 0 and at
  finite c_T alike.
* **Grade.** The divergence statement is for absolute convergence of the per-block integral in CLS-D's own representation. A
  conditionally convergent or resummed definition is not excluded, but none is used anywhere in the paper.

## 5. Interpretation

1. **The ladder is exact on light states where the functional is defined.**
   * At J = 0 its Δ-dependence is exactly 1 + c(ν)/T. Δ = 20, 30 are light by any standard (τ ≪ Δ_g), and there the direct integral
     confirms the coefficient to 2·10⁻⁴ with no 1/T² term.
   * The J = 2, 8 tests confirm the ladder's finite-T structure (U_{b,k} = Π(C−u_j)/(T−u_j), not r^k) at r = C/T ≤ 0.09, where the
     leading-Regge truncation is off by 10⁻²–10⁻¹; the converged J = 8 residual is ≲ 10⁻³.
   * No e^{−cΔ}-type term is missing. At Δ = 10–20 such a term would be O(1), and the residuals are 10⁻²–10⁻⁵ against controls
     10–1000× larger.
2. **The giant numbers are exact polynomials evaluated far out.**
   * The per-ν image of a state of spin J is a polynomial of degree ≈ J+2 in ν² (times state-independent rational factors).
   * The direct integral confirms it on ν ≤ 6–8. Analyticity in ν, together with polynomial identity on an interval, fixes it everywhere.
   * At Δ_g = 1000 the weight extends to ν ≈ Λ + 100, so J = 8 gives Λ^{~14} ≈ 10⁴². For J ≳ Λ/2 the partial sums of an undamped
     (r = 1) ladder sample the back-to-back point and give e^{πΛ/2} (m1_remainder §3).
3. **Twist ≤ 5 is defined only by continuation.**
   * The light towers that dominate R (twist 2, and twist 4 at 10⁻⁶ of it) lie below τ = 5. There CLS-D's per-block integral diverges,
     so "the action of Ω^EE on [TT]_{0,J}" is only the analytic continuation in Δ. That continuation is the ladder.
   * Row L0's "continuation [D-sketch]" can therefore be upgraded: [DV on τ > 5; unique continuation below].
   * The same fact undermines the *bookkeeping* that sums these continued values into the light side. H_K's Tonelli argument
     (R3: Σ_O ∫|K||dDisc G_O| ≤ ∫|K| dDisc𝒢) needs a convergent per-block integral, and in the representation that defines ω[O] there
     is none for τ ≤ 5. Whether the integral of the same block in the cross-ratio representation of R3 converges was not determined
     here [I].
4. **What the check does not resolve.**
   * `m1_remainder` §5: a block-wise R of size 10⁶⁸¹/c_T² is incompatible with Σ_heavy ≥ 0 and bounded heavy actions at realistic c_T.
     With the per-state values now confirmed, this tension is genuine. It must be resolved by the light-side bookkeeping (block-wise
     H_K on the light tower), by a failure of an input in such theories, or by a cancellation among all O(1/c_T²) light families that
     the tree-level twist-2 estimate does not see.
   * Light towers of twist 6+ have convergent per-block integrals and are still large term-wise (about 10⁻¹² of twist 2 at the ladder
     edge, per m1_remainder §2). So the τ ≤ 5 divergence alone does not dissolve the tension [I].

## 6. Honest size of R and the c_T threshold
* **Term-wise (block-wise light side).** R₂ = c_T²R ≈ 10^{679}–10^{682} and c_T* ≈ 10^{680}–10^{683} at Δ_g = 1000 [E, inputs now
  checked]: the per-state actions [DV here], the tree-level weights [E, m1_remainder]. At Δ_g = 300: c_T* ≈ 10^{209}–10^{211}.
* **Without the block-wise light side.** No estimate. The functional's light side is then −Ĉ_θ[light correlator] at O(1/c_T²), a
  one-loop quantity no one has evaluated. The original "power Δ_g^p/c_T²" has no derivation behind it either.
* **No smaller sound number exists.** A smaller threshold would need either (i) the one-loop light side evaluated as a whole (CHMRS-type
  arc/Laurent coefficient of the O(1/c_T²) light correlator), or (ii) a functional of reach Λ′ ≈ (2/π)ln c_T (m1_remainder §6), with
  its own positivity certificate.

## 7. Proposed replacement text (not applied)

**Supplement, "Order of limits", replace from "The size of the light-side remainder is not a power of Δ_g" to "The statements hold
as c_T → ∞.":**

> The size of the light-side remainder is not a power of $\Delta_g$. On a light block the action of $\Omega^{EE}$ is the exact ladder
> formula~\eqref{eq:fp-target}: a direct evaluation of CLS-D's per-block integral over the physical chain, exact in $\Delta$, agrees
> with it on light states ($\Delta=10$--$100$, $J=0,2,8$) to $10^{-2}$--$10^{-7}$, including its $1/T_\Delta$ terms, so the
> continuation is not an approximation. For twist $\le5$, which includes the twist-2 and twist-4 double traces, CLS-D's per-block
> integral does not converge (the integrand behaves as $B^{(\tau-5)/2-1}$ at the light cone of the exchanged operator), and the
> action is defined only by this continuation. Because the weight reaches $\nu\approx\Lambda+100$ and a light ladder is undamped, the
> continued actions grow to $e^{\pi\Lambda/2}10^{4.5}T_\Delta$ (L0). Summed block by block with tree-level Einstein data they give
> $|c_T^2R|\approx10^{679}$--$10^{682}$ at $\Delta_g=1000$, so the inequality above would constrain the exact $t_4$, $\alpha_T$ only
> for $c_T\gtrsim10^{681}$. A remainder of this size is incompatible with the sum rule at smaller $c_T$ (Sec.~\ref{sec:HK}, R3 and
> Theorem FP), so the block-wise light side, which for twist $\le5$ sums continued values rather than convergent integrals, is not
> controlled at any physically realised $c_T$. Eq.~(5) is a statement about $c_T\to\infty$ at fixed $\Delta_g$, with or without a bulk
> effective theory; the angular bound, Eq.~(4), carries no such constant.

**Supplement, Scope, replace "and any bound on the size of the $1/c_T$ remainder below $c_T\approx10^{690}$ (row L0)" with:**

> and any control of the light remainder at finite $c_T$: block by block it is $\sim e^{\pi\Lambda/2}c_T^{-2}$ ($c_T^*\approx10^{681}$ at
> $\Delta_g=1000$, row L0), the per-block actions of the twist-2 and twist-4 double traces exist only as analytic continuations, and a
> remainder of this size is inconsistent with the sum rule (open)

**Table inputs, row L0, grade column, replace "[NE] ($J\le1000$, ball arithmetic with the moment balls of Theorem FP); continuation to
light states and $J\to\infty$ asymptotics [D-sketch]" with:**

> values [NE] ($J\le1000$, moment balls; reproduced by $\nu$-quadrature of the published weight); exact ladder image on light states
> [DV] where CLS-D's per-block integral converges (twist $>5$; $k=2$: $J=0$ to $2\times10^{-5}$ at $\Delta=20$, $J=2,8$ [NE]; $k=3$ [NE]);
> for twist $\le5$ the action is the unique continuation in $\Delta$; the block-wise sum as light side of the finite-$c_T$ sum rule
> [I]; $J\to\infty$ asymptotics [D-sketch]

and in the statement column, after "no sign", add: "; the defining per-block integral converges only for twist $>5$".

**Letter, replace "Bounded term by term, their remainder is of order $10^{686}/c_T^2$, so Eq.~\eqref{eq:cubic} constrains the actual
$\tau$ only for $c_T\gtrsim10^{690}$ unless their sum over spin cancels; it is a statement about the planar limit." with:**

> The functional's action on these states is exact but grows like $e^{\pi\Lambda/2}$, so their remainder is not controlled at any
> realistic $c_T$ (term by term it is $\sim10^{681}/c_T^2$), and Eq.~\eqref{eq:cubic} is a statement about the planar limit; the angular
> bound below has no such constant.

## 8. Files
* `md_ladder.py`: z_m, y_m, U_{b,k}, the per-ν k = 2 ladder polynomial (sympy, exact).
* `md_k3.py`: q_n, Q_k, the k = 3 per-ν polynomial and ĥ⁽³⁾.
* `md_hhat_quad.py`, `md_hhat_total.py`: independent per-state values.
* `md_dquad.py`: direct D′ integral, k = 2, J = 0 (kind E2) and scalar (S).
* `md_spinJ.py`: spin-J integrands, k = 2 (`ratio`) and k = 3 (`ratio3`).
* `md_dquadJ.py`: direct D′ integral, spin J; `MD_KIND=C3` selects k = 3.
* `md_test_spinJ.py`, `md_probe.py`, `md_probe3.py`, `md_decay.py`: tests and boundary exponents.
* `md_compare_J0.py`, `md_kappa.py`, `md_compare_J.py`, `md_compare_k3.py`, `md_conv.py`, `md_smallgap.py`: comparisons.
* `run_*.sh`: the queues used.
* `out/`: all results (`dq*_*.npz` raw integrals, `*.txt` tables, `log_*` run logs).
