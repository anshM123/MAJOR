# Scope: CFT-exact thermal sum rule (CDP shear sum rule) + eta/s, Hall channel, inverse gap statement

Status: 27 Nov 2026, research agent, not refereed. Nothing in `prl/` or elsewhere was edited; nothing committed.
Grades: **[cited]** quoted from the named source (read from its arXiv LaTeX where an equation is given);
**[computed]** by a script in this directory (recorded output `out_*.txt`); **[estimated]** argument or inference, no
independent computation.

| script | what it does | runtime |
|---|---|---|
| `free_fields_sumrule.py` | direct one-loop shear sum rule for a free Dirac fermion and a free real scalar in 2+1D (pair continuum + transport peak) | <5 s |
| `einstein_ads4_sumrule.py` | numerical integral of `rho_T - rho_0` on the AdS4 black brane (Einstein) | ~5 min |
| `hall_uv_sign.py` | Q = C C C~ coupling: Hall correlator h(i Omega) from Omega -> 0 (eta_H) to Omega -> infinity (UV/OPE constant), same conventions as `../scope_thermal/hall/hall_ef.py` | ~2 min |
| `gap_inverse.py` | inversion of tau <= tau_*(Delta_g) in terms of a measured sum-rule deviation | <1 s |

---

## 0. Summary

| item | verdict | grade |
|---|---|---|
| CDP formula, d = 3: `S_E = 3/8 + t4/32` | **correct**; t4 normalisation = HM/Buchel = this paper | [cited] + [computed] |
| Einstein value 3/8 | numerically 0.37499996 | [computed] |
| coefficient 1/32 | free Dirac (t4=-4) gives exactly 1/4, free scalar (t4=+4) exactly 1/2, so the slope is 1/32 | [computed] |
| corrected KSS `N = (4pi eta/s - 1) + (12/5)(S_E - 3/8) = 216 alpha4 + ...` | **yes**, t4 cancels with the stated signs; remainder `-(17/3840)t4^2 + O(alpha_T^2) + O(Lambda^-8)`; bulk-EFT statement | [computed algebra] on [cited]/[computed] inputs |
| inverse gap statement | `Delta_g <= 2.2401 + 5/(9 delta)` (the coordinator's `2.26 + 5/(9 delta)` is valid); not a bulk-EFT statement, but it inherits assumption (ii) | [computed] |
| Hall sum rule, magnitude `|S_H| = |alpha_T|/32` | **confirmed in the bulk at linear order** (S_H = -45 lam_Q, i.e. 7.5|alpha3~| = 240|alpha3~|/32) | [computed] (bulk); CFT derivation [estimated] |
| **Hall sign: does Z cancel?** | **NO.** `S_H = -(5/12)(4 pi eta_H/s)` for the cubic odd coupling, so in `Z = V + (12/5)(S_E - 3/8 + i S_H)` the cubic imaginary parts **add**. The cancelling combination is `Z' = V + (12/5)(S_E - 3/8 - i S_H) = 216(alpha4 + i alpha4~) + ...` | [computed], two independent routes |
| novelty | no prior combination of the shear sum rule with eta/s; no gap-dependent bound on the sum rule; no CFT3 Hall-channel sum rule found | [cited] (searches) |

---

## 1. The CDP sum rule

**Reference** [cited]: S. Dutta Chowdhury, J. R. David, S. Prakash, *Spectral sum rules for conformal field theories
in arbitrary dimensions*, arXiv:1612.00609 (JHEP 07 (2017) 119). Main result eq. (1.5) (label `mainres`), repeated in
t2,t4 form as eq. (4.3) (label `bounds`):

    lim_{eps->0+} (1/pi) Int_{-inf}^{inf} d omega  delta rho(omega)/(omega - i eps)
        = [ d(d-1)/(2(d+1)) + (3-d) t2/(2(d-1)) + (2 + 3d - d^2) t4/((d-1)(d+1)^2) ] P,
    rho = Im G_R,  G_R(t,x) = i theta(t) <[T_xy(t,x), T_xy(0)]>,  delta rho = rho_T - rho_{T=0},  p = 0.

At d = 3: the t2 coefficient vanishes, the constant is 3/4, the t4 coefficient is 2/32 = 1/16; with eps = 2P:

    S_E := (1/(pi eps)) Int_{-inf}^{inf} d omega [rho_T - rho_0]/omega = 3/8 + t4/32 .        [cited, algebra checked]

* **Integration range**: the whole real line (rho is odd, so this equals (2/(pi eps)) Int_0^inf). Writing Int_0^inf with
  1/(pi eps) would halve everything. The coordinator's S_E with `(delta rho_+ + delta rho_-)` over (0, inf) equals this
  (`delta rho_+ + delta rho_- = 2 Im delta g`).
* **Sign convention**: CDP's G_R = +i theta <[ , ]> gives rho/omega -> +eta > 0. Holographic check below uses the same.
* **t4 normalisation** [cited + computed]: CDP take t2, t4 from Buchel-Escobedo-Myers-Paulos-Sinha-Smolkin 0911.4257
  (their eq. (4.1)-(4.2)); their d = 3 collider window is `-4 <= t4 <= 4`; their CS-matter result
  `delta G = -(P/4)(cos theta - 3)` puts the free fermion (theta = 0) at P/2 and the free boson at P, i.e. free scalar
  t4 = +4, free fermion t4 = -4. That is the Hofman-Maldacena normalisation used by this paper (letter.tex l. 77-79,
  `DICTIONARY_T4.md` §1: free scalar +4, Majorana -4, `|t4| <= 4`). Li-Lü-Mai 1901.03349 (source of
  `4 pi eta/s - 1 = -(3/40) t4`) use the same t4 (collider `|t4| <= 4`; ECG sign checked in `../scope_thermal`).
  So **CDP t4 = HM t4 = Buchel t4 = this paper's t4**, no rescaling.
* CDP's own consequence [cited]: with the collider bounds, `P/2 <= delta G_R(0) <= P`, i.e. `1/4 <= S_E <= 1/2` for
  parity-even CFT3s (their eq. `3dcaus`). Bounding the sum rule by bounds on t4 is therefore CDP's idea; the
  gap-dependent version is not in the literature (Sec. 5).

### 1.1 Independent checks of 3/8 and 1/32

* **Einstein, AdS4 black brane** [computed, `einstein_ads4_sumrule.py`]: `S_E = (1/pi) Int_0^inf [|phi_h/phi_bdy|^2 - w^2] dw`
  (r_+ = 1, with eta/s = 1/4pi reproduced at w -> 0). Integral to w = 30 plus an OPE tail fit:
  **S_E = 0.3749999601** (3/8 - 4e-8). The integrand decays as `w^(-3.99)`, matching the OPE: the first thermal
  operator after T in T x T on the brane is the double trace [TT] (Delta = 6), so delta G ~ w^-3 and delta rho/w ~ w^-4.
  This is also the holographic result of Gulotta-Herzog-Kaminski 1010.4806 and CDP §3 (`d/(2(d+1)) eps`), and the
  d = 3 analogue of Romatschke-Son 0903.3946 (2/5 in d = 4) [cited].
* **Free fields** [computed, `free_fields_sumrule.py`, exact to 1e-10]. At p = 0 the spectral function is a pair
  continuum plus a zero-frequency transport peak, with `(1/pi) Int rho_peak/w = W = Sum_k |<k|T_xy|k>|^2 (-dn/dE)`.
  Under any weak interaction W becomes the hydrodynamic peak `rho/w = eta/(1 + w^2 tau^2)`, of weight eta/tau.

  | theory | t4 | pair part / eps | peak W / eps | S_E | CDP 3/8 + t4/32 |
  |---|---|---|---|---|---|
  | free Dirac (= 2 Majorana) | -4 | -1/8 (Pauli blocking) | 3/8 | **1/4** | 1/4 |
  | free real scalar | +4 | +1/8 (Bose enhancement) | 3/8 | **1/2** | 1/2 |

  In d = 3 the parity-even <TTT> is spanned by the free-boson and free-fermion structures, so `S_E = a + b t4` exactly.
  The two free points alone fix **a = 3/8, b = 1/32**, and the Einstein point (t4 = 0, 3/8) is an independent third
  check. The t4 coefficient had not been checked independently before: CDP's conclusions say so, and the follow-up
  1711.01027 (Dutta Chowdhury) checks only the t2 (Gauss-Bonnet) coefficient in d >= 4 [cited].
* The free scalar matches even though phi^2 (Delta = 1) has a nonzero (IR-divergent) thermal vev. This is consistent:
  the scalar's tail coefficient in the T_xy T_xy channel is proportional to `(D - Delta - 2)`, which vanishes at
  Delta = 1 in D = 3 (Lucas-Sierens-Witczak-Krempa 1704.05461, eq. `ABvisc`, coefficient B) [cited].

### 1.2 Conditions, and whether a holographic CFT with a gap meets them

1. **Analyticity / hydrodynamics** [cited]: G_R analytic in the UHP; `G_R(0)|_T = P` from hydrodynamics (CDP eq.
   `hydro`), with contact terms fixed by `<T_xy> = -(1/sqrt g) d ln Z/d g_xy`. The LHS itself is contact-term free,
   since it is `delta G(0) - delta G(i infinity)`. The split of the RHS into `P - J` needs matching conventions, and the
   free-field check above confirms the combination in contact-term-free form. At finite c_T in 2+1D, long-time tails make
   eta(omega) ~ log omega as omega -> 0, which is integrable; the sum rule survives, though eta itself does not exist
   at finite c_T.
2. **Operators** [cited + estimated]: CDP assume that no operator other than T with Delta <= d in T x T gets a thermal
   vev. In d = 3, thermal one-point functions vanish for odd spin, and unitarity puts non-conserved spin-2 at
   Delta > 3 and spin >= 4 at Delta >= 5. So the condition is exactly: **no scalar (or, if P is broken, pseudoscalar)
   O in T x T with Delta_O <= 3 and b_O = <O>_beta != 0.** Refinements:
   * For 1 < Delta < 3 the unsubtracted S_E **diverges**: the tail goes as `w^(3-Delta)` with a nonzero imaginary part.
     One must subtract `B <O>_T (i w)^(3-Delta)` in the integrand. The RHS is then unchanged, still 3/8 + t4/32 (CDP §2
     end; Witczak-Krempa 1501.03495; LSW 1704.05461 §7) [cited].
   * Delta = 1 is harmless (B ∝ D - Delta - 2 = 0) [cited + computed free scalar]. Delta = 3 (marginal) can shift the
     RHS (R-charge chemical-potential cases, CDP §2 and their ref. David:2011hy) [cited].
   * Delta > 3, e.g. single traces between 3 and the gap: the sum rule is unaffected and converges as
     `w^(2-Delta)` [estimated].
3. **Holographic CFT with a gap** [estimated]: at large c_T a single-trace scalar with lambda_TTO != 0 (bulk phi W^2)
   is sourced on the black brane, so b_O ∝ lambda_TTO at the same order, and its contribution ∝ lambda_TTO^2 is not
   1/c_T-suppressed. The CDP condition is therefore equivalent to **no single-trace scalar of Delta <= 3 coupling to TT**.
   * This is implied by the "no light single trace in T x T" assumption of the 2.89e-8 bound.
   * It is **stronger** than the tau_* assumption, which excludes only Delta <= 1.
   * Double traces with Delta <= 3 need Delta_O <= 3/2 and enter at O(1/c_T).
   * Meltzer-Perlmutter 1712.04861 (§ applications) note that phi W^2 thermal vevs (Myers-Sierens-Witczak-Krempa;
     Lucas et al.) are parametrically small without higher-spin fields [cited].
4. **Depends on nothing but t4 and eps**: **yes**, under 2. The RHS is P - J, where J is the T-exchange OPE coefficient
   times <T>_beta. R^4, alpha4, alpha4' and all higher-point bulk couplings do not change <TTT>, so they do not change S_E.
   The sum rule is non-perturbative in 1/c_T: CDP use only the OPE and hydrodynamics. **S_E receives no alpha4 term**
   [cited + estimated].

---

## 2. Corrected-KSS combination (real part)

Inputs: `4 pi eta/s - 1 = -(3/40) t4 + 216 alpha4 - (17/3840) t4^2 + ...` (bulk EFT, `../scope_thermal/README.md`;
Li-Lü-Mai (7.15) for the t4, t4^2 terms) and `S_E - 3/8 = t4/32` (exact). Since (12/5)(1/32) = 3/40:

    N := (4 pi eta/s - 1) + (12/5)(S_E - 3/8) = 216 alpha4 - (17/3840) t4^2 + O(alpha_T^2) + O(Lambda^-8).   [computed algebra]

Signs check out: t4 > 0 lowers eta/s and raises S_E (free scalar side), consistent with ECG (t4 = -1260 mu < 0,
eta/s = 1 + 189 mu/2). With flat-space alpha4 >= 0 (CHLPS; not established in AdS by the paper), the result is
`4 pi eta/s >= 1 - (12/5)(S_E - 3/8) - (17/3840) t4^2 - O(Lambda^-8)`. Caveats:
* **Not CFT-exact.** S_E is exact, but eta/s needs the bulk EFT on the black brane and power counting of the higher-point
  couplings (`../scope_thermal` §4.1, §5). The t4^2 and alpha_T^2 pieces are not cancelled; they are `<= 1.4e-6` at
  tau = tau_*(1000) and negligible under positivity.
* A scalar with phi W^2 coupling and 3 < Delta < Delta_g shifts eta/s at O(lambda^2) **without** shifting S_E. If heavy, it
  acts like a positive alpha4 [estimated]; in general it must be excluded.

## 3. Inverse statement [computed, `gap_inverse.py`]

    delta_E := |S_E - 3/8| = |t4|/32 <= tau/32 <= tau_*(Delta_g)/32
    <=>  (Delta_g - 2) sqrt(1 - 49/(2 Delta_g + 9)^2) <= 5/(9 delta_E)
    =>   Delta_g <= 2.2401 + 5/(9 delta_E)      (sup of the correction is 0.2401, at Delta_g = 8.25)

tau_* is strictly decreasing, so the exact form is `Delta_g <= g_*(32 delta)`. The approximation `2 + 5/(9 delta)` is
within 0.24 everywhere; for example delta = 0.01 gives 57.64 (exact) against 57.56. The collider bound caps delta at 1/8.
At Delta_g = 1000 the allowed deviation is `|S_E - 3/8| <= 5.6e-4` without positivity and `<= 9.0e-10` with it.

**Assumptions** (complete list):
* everything tau_* needs: (i) higher-spin gap; (ii) c_T -> infinity with a tree-level description below the gap,
  uniformly in spin; detector identities as operator identities; convergence of the mixed detector commutator; no
  scalar with Delta <= 1 in T x T;
* the CDP conditions (§1.2): no scalar or pseudoscalar with 1 < Delta <= 3 and nonzero thermal vev in T x T, or else the
  subtracted version; G_R(0) = P (hydrodynamics in the thermal state on R^2);
* the measured S_E must be its c_T = infinity value, because tau_* holds at leading order in 1/c_T although S_E itself
  is exact.

**Is it a pure CFT statement?** It uses no black-brane bulk EFT and no power counting of higher-point couplings: the
thermal side is exact CFT. It does inherit assumption (ii) of the vacuum bound, which is phrased as a tree-level bulk
description of the light data. Fair wording: "CFT statement at leading order in 1/c_T, under the same assumptions as
tau_*, plus the absence of light (Delta <= 3) scalars with thermal vevs in T x T."

---

## 4. Hall channel (coordinator's extra items)

Setup as given: `T1 = (Txx - Tyy)/2`, `T2 = Txy`, `G^R = [[g, -h], [h, g]]`, `G_± = g ± i h`, `rho_± = Im G_±`.

### 4.1 Exact statements about h [computed derivations; standard complex analysis]

* Reality: `h(t)` is real, so `h(-w) = h(w)*`; `rho_+(-w) = -rho_-(w)` (confirmed); `delta rho_+ - delta rho_- = 2 Re delta h`.
* Euclidean antisymmetry: `G_E^[AB](w_n)` is real and odd in w_n. Hence `delta h(0) = 0`, and **no temperature-dependent
  contact term exists in the odd channel at p = 0**: a local term ∝ P would have to be even. The T = 0 part, the
  gravitational Chern-Simons contact term ~ kappa_g w^3, is T-independent and drops out of delta.
* UV: `delta h(i Omega) -> J_H` (real constant, the T-exchange OPE coefficient with the parity-odd <TTT> structure).
* **Coordinator's S_H**: `lim_Omega (Omega/pi eps) Int_0^inf 2 Re delta h/(w^2 + Omega^2) = Re delta h(i Omega -> inf)/eps = J_H/eps`
  (Poisson representation; needs delta h bounded in the UHP, which the OPE gives). **Correct**, but it is an Abel
  average of a function that tends to the constant J_H, not a convergent sum rule. The equivalent **convergent** form
  (Kramers-Kronig on delta h - J_H, which decays like w^-3) is

      S_H = J_H/eps = (1/(pi eps)) Int_{-inf}^{inf} eta_H(w) dw,     eta_H(w) := -Im delta h(w)/w,   eta_H(0) = eta_H.

  This is the parity-odd analogue of CDP: the frequency-integrated Hall viscosity is fixed by alpha_T and eps.
* **IR identity**: `eta_H = (1/pi) Int_0^inf (delta rho_+ - delta rho_-)/w^2 dw` is **correct**, with eta_H defined by
  `h = -i eta_H w + O(w^2)` (contour around the simple pole of h/w^2 at 0). In the conventions of `hall_ef.py`
  (`<T^xy> = ... + eta_H phi_v`), this is the same eta_H [computed: the Omega -> 0 limit below gives 108].
* These two identities do not relate J_H to eta_H. That relation needs a model; here, the bulk (4.3).

### 4.2 (a) |S_H| = |alpha_T|/32

* **Bulk, linear order** [computed, `hall_uv_sign.py`]: for `L ⊃ lam Q` (Q = C C C~, eps_{vrxy} = +r^2), the large-Omega
  limit gives **S_H = -44.99 ± 0.03 per unit lam** (fit `-89.98 + 27.8/Omega^2 + ...` for Pi = 2 S_H, Chebyshev N = 120).
  With `lam = -alpha3~/6`, `|S_H| = 7.5 |alpha3~| = 240|alpha3~|/32 = |alpha_T|/32`, using the paper's
  `tau = 240|ghat3|` (DICTIONARY_T4 §4, `|alpha_T| = 240|alpha3~|`, graded [D] there). The magnitude claim is
  therefore **confirmed for holographic CFTs at linear order**. Machinery check: the Omega -> 0 limit reproduces
  `4 pi eta_H/s = 108 lam` (107.88 at Omega = 1e-3, linear approach) and sigma(r) = 432/r^5 from `hall_ef.py`.
* **CFT derivation** (rotation R_k on the hard legs, Jain-John 2107.00695 eq. (5.1)) [estimated]: not re-derived. It is
  plausible, and the bulk number agrees, but it is not an independent CFT proof. The CDP-normalisation of alpha_T
  (their 1707.03007 flux matrix) relative to the bulk alpha3~ is itself graded [D] in DICTIONARY_T4. A free-field
  check is impossible because parity-odd <TTT> needs interactions. A possible CFT check is large-N CS-matter, which lies
  on the circle `t4 + i alpha_T = 4 e^{i theta}`, prediction `|S_H| = |sin theta|/8`. But there psi-bar psi
  (Delta = 2, nonzero thermal vev when P is broken) needs a subtraction, and N = infinity gives delta-function
  transport peaks.
* **Light (pseudo)scalars**: the same caveat as §1.2 applies to the odd channel (parity-odd <TTO> structures), with a
  subtraction needed for Delta < 3.

### 4.3 (b) The relative sign: **Z as written does not cancel** [computed, two routes]

The ratio `S_H / (4 pi eta_H/s)` is convention-free once both are defined from the same h, as in the coordinator's
setup: orientation and helicity labelling flip both. In the conventions of `hall_ef.py`, where
`V = (4 pi eta/s - 1) - i 4 pi eta_H/s = 18 ghat3 + 216(alpha4 + i alpha4~)`:

| Omega (h(i Omega) for lam = 1, units of eps/2) | 1e-3 | 0.5 | 1 | 2 | 5 | 10 | 20 | 40 | inf (fit) |
|---|---|---|---|---|---|---|---|---|---|
| Pi(oo) | 0.1079 (= 107.88 Omega) | 28.02 | 22.78 | -13.62 | -72.47 | -87.08 | -89.58 | -89.92 | **-89.98** |

So h(i Omega) starts at `+eta_H Omega` and ends at `J_H` of the **opposite sign**:

    S_H = -(5/12) (4 pi eta_H/s)     (cubic odd coupling; -0.41659 numerically vs -5/12 = -0.41667)

* Im Z = `-4 pi eta_H/s + (12/5) S_H = -2 (4 pi eta_H/s)`: the cubic parts **add (double)**.
* The cancelling combination is

      Z' = V + (12/5)(S_E - 3/8 - i S_H) = 216 (alpha4 + i alpha4~) + O(t4^2, alpha_T^2, Lambda^-8).

**Second route, a structural argument** [estimated, agrees with the computation]: the helicity-+ correlator
`G_+ = g + i h` has IR expansion `P + i w (eta - i eta_H)` and UV data `S_+ := (G_+(0) - G_+(i inf))/eps = S_E - i S_H`.
The IR mode, `e^{-iwv}` in EF, and the UV mode, the continued ingoing mode localised near the boundary, are both
continuations of the ingoing solution, so they have the same duality chirality. Writing
`C^3 + i C^2 C~ -> ghat3 C_+^3 + ghat3* C_-^3`, both limits of G_+ depend on a single complex coupling, which makes them
holomorphic with the real-coupling (parity-even) coefficients: `4 pi (eta - i eta_H)/s - 1 = 18 g` and
`S_+ - 3/8 = -7.5 g`. Then `(4 pi eta_+/s - 1) + (12/5)(S_+ - 3/8) = 0` at cubic order. Its complex conjugate is
exactly Z', not Z.

**Disk statement for Z'** [computed algebra + cited]: with `G4 = alpha4 + alpha4'`,
`|Z' - 108 G4| <= 108 G4  <=>  alpha4~^2 <= alpha4 alpha4'`. That is the flat-space positivity condition
`|ghat4| <= g4` as corrected in `../scope_thermal/hall/README.md` (with `ghat4 ∝ alpha4 - alpha4' ± 2i alpha4~`). Note
that `RESULT_Q_STAGE0.md` quotes `alpha4~^2 <= 4 alpha4 alpha4'`; with that weaker form the disk does not follow. So
the disk holds for Z', is flat-space-conditional, and depends on resolving that factor-2.

**Recommendation**: either state only the real-part (corrected-KSS) result, or state the complex one with
`S_E - 3/8 - i S_H` (equivalently, helicity-resolved `eta_+ = eta - i eta_H`, `S_+ = S_E - i S_H`) and cite the bulk
sign computation.

### 4.4 (c) Gap inequality with the Hall channel [computed]

`delta_th = sqrt((S_E - 3/8)^2 + S_H^2) = tau/32` uses only magnitudes, so it is **unaffected by the sign issue**.
`(Delta_g - 2) sqrt(1 - 49/(2Delta_g + 9)^2) <= 5/(9 delta_th)` is correct, and `Delta_g < 2.26 + 5/(9 delta_th)` is
valid; the sharp constant is 2.2401. Additional assumption beyond §3: `|S_H| = |alpha_T|/32` as a CFT statement, which
is verified only holographically at linear order (4.2). Without it, use delta_E (even channel), which is CFT-exact.

### 4.5 (d) Novelty of a CFT3 Hall-channel sum rule [cited: searches]

None found. Related work:
* Bradlyn-Goldstein-Read 1207.7021: Hall viscosity drops out of the viscosity (f-sum-type) sum rule.
* Golkar-Nguyen-Son 1309.2638, and JHEP 01 (2016) 021: odd stress spectral sum rules in quantum Hall states,
  involving the shift and the Hall viscosity; these are gapped, non-relativistic systems with no OPE.
* Gur-Ari-Hartnoll-Mahajan 1605.01122: finite-T Hall conductivity in CS-matter; they use the longitudinal conductivity
  sum rule only.
* Saremi-Son 1103.4851 and Hoyos et al. 1906.04752: holographic eta_H, which needs an axion or a magnetic field.

INSPIRE full-text searches for "Hall viscosity" + "sum rule" (51 hits, all quantum Hall/superfluid) and for "parity-odd" +
"spectral sum rule" found nothing on CFT3 thermal stress tensors.

---

## 5. Novelty of the even-channel ideas [cited: searches]

* All 12 INSPIRE citers of CDP 1612.00609 were checked (1711.01027, 1707.03007, 1704.05461, 1712.04861, 1905.01355,
  1711.11300, 2101.11635, 2611.07394, 2611.24531, reviews). None combines S_E with eta/s, and none uses t4 bounds
  beyond CDP's own collider window `1/4 <= S_E <= 1/2`.
* David-Prakash-Dutta Chowdhury follow-ups (1707.03007 parity-odd collider; 1711.01027 higher-derivative check, d >= 4,
  t2 only; 1812.07774 parity-odd bootstrap): no viscosity, no gap.
* Prokhorov-Teryaev 2604.04222 ("Causality, the KSS bound, and a novel sum rule"): relates eta/s to the sound speed for
  Unruh radiation; unrelated to t4 or the shear sum rule.
* **Verdict**: (1) the t4-cancelling combination N, (2) a gap-dependent bound on the shear sum rule, (3) the inverse
  "sum rule measures the string scale", and (4) any CFT3 Hall-channel sum rule all appear **new**. Each is a few lines
  of algebra on top of CDP + Li-Lü-Mai + the paper's tau bound; the Hall sign (4.3) is the one non-trivial new
  computation.
