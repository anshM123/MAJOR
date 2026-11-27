# Hall viscosity of the AdS4 black brane from parity-odd curvature terms: independent check

Status: 27 Nov 2026, research agent, not refereed. No edits to `prl/`, nothing committed.
Grades: **[computed]** by a script here; **[cited]** from the named source; **[I]** inferred with argument; **[failed]** = claim not confirmed.

| file | content | runtime |
|---|---|---|
| `hall_ef.py [ORIENT]` | explicit-index (no bivector matrices) eps^2-truncated curvature of the time-dependent two-polarisation shear perturbation in ingoing EF; cross terms of P, Q; chi-EOM at O(omega); eta_H for P, Q, IP, theta(r)P, and higher-point odd operators; reciprocity check | ~5 s |
| `out_hall_ef_orient+1.txt`, `out_hall_ef_orient-1.txt` | recorded output for eps_{vrxy} = +r^2 and -r^2 | |
| `quartic_IP.py` | flat 4-graviton contact amplitudes of I^2, P^2, IP (same conventions as `../quartic_flat.py`) | <1 s |

## 1. Coefficients (item 1)

### Setup
The metric is `ds^2 = -f dv^2 + 2 dv dr + r^2 (e^{eps H})_ij dx^i dx^j`, with `H = [[phi, chi], [chi, -phi]]`, `phi, chi` functions of `(v, r)`, and `f = r^2 - 1/r`.
* `e^{eps H}` is expanded to `eps^2`, which is exact at that order (det = 1).
* The Levi-Civita tensor is `eps_{vrxy} = ORIENT * r^2`.
* Every Weyl contraction is an explicit index sum over sparse component dictionaries. This does not use the collaborator's bivector trick.

### Raw quadratic terms [computed]
With `ORIENT = +1`, the terms of the form `phi_v x (r-derivatives of chi)` are:
* `P^(2) ⊃ -(6/r^2) phi_v chi_rr`. This **matches** the claim.
* `Q^(2) ⊃ (54/r^6) phi_v chi_r + (18/r^5) phi_v chi_rr`. This **matches** the claim.

The full cross-term lists, including `phi_r chi_v`, `phi_rr chi_vv` and the rest, are in the output files.

### Extraction
* The O(omega) part of the chi equation of motion on the ingoing solution is `E_chi = E_chi^Einstein + lam sigma(r) phi_v`.
* Cross terms with no v-derivative vanish identically, which the code checks. The reason: a v-independent configuration is invariant under the orientation-reversing t -> -t, so it carries no pseudoscalar.
* Terms with undifferentiated phi, or with only r-derivatives of phi, are absent (checked). So no O(omega) radial profile of phi enters.
* The Einstein part is exactly `d_r(r^2(f chi_r + chi_v)) + d_v(r^2 chi_r)`, with normalisation N = 1 (checked).
* Therefore `d_r Pi_min = lam sigma phi_v` with `Pi_min = -r^2(f chi_r + chi_v)`, and `4 pi eta_H/s = lam ∫_1^∞ sigma dr`. The normalisation is calibrated by eta, using `T^{xy} = Pi_min(∞)/16piG` and `delta T^{xy} = -p chi - eta chi_v + eta_H phi_v`.

### Results (ORIENT = +1; all signs flip with ORIENT = -1)

| operator | sigma(r) | 4 pi eta_H/s per unit coupling | claim | grade |
|---|---|---|---|---|
| P = C.C~ (constant coefficient) | 0 | **0** (topological) | — | [computed] |
| Q = C C C~ | 432/r^5 | **+108** | +108 | **confirmed** [computed] |
| I P | -3024/r^8 | **-432** | -432 | **confirmed** [computed] |
| theta(r) P | -6 theta'' | **6 theta'(1)** (theta decaying at infinity; checked for theta = r^-k, k = 0..7) | 6 theta'(1) | **confirmed** [computed] |
| P^2, I, J, I^2 | no phi-chi cross terms at eps^2 | **0** | 0 | confirmed [computed] |

The collaborator's conserved momentum `Pi_chi = Pi_min + (108 k3/r^4 - 432 k4/r^7) phi_v` is exactly mine: its added term K(r) satisfies `K' = -sigma` and `K(∞) = 0`.

**Reciprocity** [computed]. The phi-channel response to `chi_v` is -108 (Q) and +432 (IP). The response matrix is therefore antisymmetric, as Onsager requires for a Hall coefficient.

**Saremi-Son check** [cited + computed]. Saremi-Son (arXiv:1103.4851, their Kubo section) find `eta_A = -(lam/8piG) r^4 f_SS' theta'/4 |_{r_h}` for `L ⊃ -(lam/4) theta *RR`.
* Their `*RR` equals `-R R~`, which equals `-P` on Einstein spaces. So `theta_ours = lam theta_SS / 4`.
* Their result then becomes `16piG eta_A = -6 theta_ours'(1)`, which matches my |6 theta'(1)|.
* The sign agrees if their convention is `eps^{vrxy} = +1/sqrt(-g)`, i.e. ORIENT = -1. I did not verify their epsilon convention.
* IP equals `theta P` with `theta = k4 I0 = 12 k4/r^6`, so this cross-check covers the IP coefficient directly.

**Mapping to the action.** The action term `-(a3 J + a3~ Q)/6 + (a4 I^2 + a4' P^2 + 2 a4~ I P)/4` gives:

    4 pi eta_H/s = -18 a3~ - 216 a4~      (ORIENT = +1)      CONFIRMED
    4 pi eta/s   = 1 + 18 a3 + 216 a4      (confirmed earlier, ../README.md)

So `(4pi eta/s - 1) - i (4pi eta_H/s) = 18 ghat3 + 216 (a4 + i a4~)`. The cubic part rotates as a single complex number. This is the same duality structure as `t4 + i alpha_T ∝ ghat3`, and it is a structural consistency check.

**Higher-point odd operators.** These are invisible in vacuum TTT and TTTT. Their coefficients [computed] are:
* I Q: +2592
* J P: -648
* I^2 P: -10368
* J Q: +3240

With power counting `c/Lambda^8`, the remainder is about `2.6e-21 c` at `Delta_g = 1000`. The same caveat as for eta applies.

## 2. Signs and orientation (item 2)
* `a3~ = Im ghat3` holds by definition of the action. The action is CHLPS (2.9)-(2.11), with R -> C, which is exact on Einstein spaces up to topological and Ricci terms [cited].
* `|alpha_T| = 240 |a3~|` is the `[D]` entry of `DICTIONARY_T4.md` §4, and its sign is not fixed there. The paper's statement `tau = 240|ghat3|` (letter.tex l. 79, "In the bulk tau = 240|ghat3|") rests on it.
* Given that statement, `||V_cubic|| = 18 sqrt(a3^2 + a3~^2) = 18 |ghat3| = (3/40) tau` is an **identity** (not merely a bound), and `4pi|eta_H|/s = (3/40)|alpha_T|` holds exactly at linear order.
* **Not independently computed here**: the CFT-side normalisation of alpha_T (the CDP flux matrix). Grade [cited/I].
* The sign relation between eta_H and alpha_T is convention-dependent and unfixed. Only magnitudes should be quoted.

## 3. Physics sanity (item 3)
* **P and T.** Q and IP each contain one epsilon tensor. They are therefore odd under P and under T separately, and even under PT. Hall viscosity requires both P and T to be broken. Parity-odd ⟨TTT⟩ (Giombi-Prakash-Yin 1104.4317; Chern-Simons-matter theories) breaks P and T in the same way, so the picture is consistent [I].
* **Is eta_H a genuine transport coefficient?** Jensen, Kaminski, Kovtun, Meyer, Ritz and Yarom (arXiv:1112.4498) treat the Hall viscosity eta~ of a parity-violating normal fluid in 2+1D as a first-order transport coefficient. It is frame invariant (their frame-invariant version equals eta~), dissipationless, and unconstrained by entropy production. It has a k = 0 Kubo formula [cited]. At k = 0 the only local parity-odd contact term, the 3D gravitational Chern-Simons term, contributes at O(omega^3). So the O(omega) Hall coefficient is scheme-free [I]. In their own holographic model eta~ = 0, because tensor modes decouple.
* **Radial-momentum argument.** It is valid at O(omega). The odd terms make Pi_min non-conserved with an explicit source `sigma(r) phi_v`. The correct treatment integrates this source, which is equivalent to adding `K(r) phi_v` with `K(∞) = 0`. That is what both the collaborator and I do [computed].
  * The cross Lagrangian has no undifferentiated-chi terms, so there is no residual ambiguity.
  * The odd terms' contribution to the boundary canonical momentum vanishes, since C ~ r^-3.
  * Regularity of `chi_r` at r = 1 fixes `Pi_min(1) = -chi_v`.
* **Topological?** Constant-coefficient P gives exactly 0 [computed]; this is the Pontryagin density. Q is **not** topological: it gives a nonzero eta_H here, and it produces a nonzero on-shell all-plus three-graviton amplitude (`../../dictionary/three_point_vertices.py` (d)). A total derivative cannot do either. IP is not topological either: it is the theta(r) P mechanism with theta = I0.

## 4. Quartic inequality (item 4): a factor-2 error, harmless numerically
`quartic_IP.py` [computed] finds the following. For the Lagrangian `(1/16piG)(1/4)(a4 I^2 + a4' P^2 + 2 a4~ I P)`, the all-plus contact amplitude is `prop. to a4 - a4' - 2i a4~`, normalised exactly as CHLPS's `(1/2) ghat4 (s^2+t^2+u^2)^2` (ratio c = 1.000000 checked). The MHV amplitude is `prop. to a4 + a4'`, and IP does not contribute to it.

So **`ghat4 = 8piG (a4 - a4' ± 2i a4~)`**. CHLPS's printed `+ i a4~` (2201.06602, after eq. (2.10); the source line was read) disagrees by a factor 2 for its own Lagrangian normalisation.

The same result follows from self-dual algebra:
* write `a = C+·C+`, so that `I = a + a*` and `P = ∓i(a − a*)`;
* in the Lagrangian, the coefficient of `a^2` is `(a4 − a4' ∓ 2i a4~)/4`, and the coefficient of `a a*` is `(a4 + a4')/2`.

Consequences:
* The flat constraint `|ghat4| <= g4` reads `(a4 - a4')^2 + 4 a4~^2 <= (a4 + a4')^2`, i.e. **`a4 a4' >= a4~^2`** (with `a4, a4' >= 0`). It is not the collaborator's `(a4-a4')^2 + a4~^2 <= G4^2`.
* Their constraint is weaker, so their bound `||V_q|| <= 144 sqrt3 G4 = 249 G4` remains **valid but loose**. The algebra giving `a4^2 + a4~^2 <= (4/3)G4^2` from their constraint is itself correct.
* The sharp statement is `a4^2 + a4~^2 <= a4 (a4 + a4') <= G4^2`, so `||V_q|| <= 216 G4` and `4pi|eta_H|/s|_quartic <= 216 |a4~| <= 108 G4`.
* These are flat-space positivity statements, which hold up to loop effects. They are not established in AdS by the paper.

## 5. Numbers (Delta_g = 1000, Lambda = 998.50)
* **Cubic:** `(3/40) x 2.8914e-8 = 2.1686e-9`.
* **Quartic:** `216 G4 = 2.02e-14`, using `G4 = 92.8/Lambda^6 = 9.36e-17`. The collaborator's form gives `249 G4 = 2.34e-14`.
* **Total:** `||V|| <= 2.1686e-9`. The collaborator's 2.176e-9 corresponds to rounding tau up to 2.9e-8. It is conservative; quote **2.17e-9**.
* **Hall alone:** `|eta_H|/s <= [(3/40) tau + 108 G4]/(4 pi) = 1.726e-10`. The collaborator's 1.732e-10 again reflects the rounding. Without positivity, `(3/40) tau_*/(4pi) = 1.06e-4`.
* **Remainder from quintic odd operators:** about `2.6e-21 x` (natural coefficient).

## 6. Novelty (item 5)
I searched the following:
* Saremi-Son 1103.4851 (axion θRR̃).
* Chen, Lee et al. 1110.0793, "A holographic model for Hall viscosity" (axion θRR̃, charged brane).
* 1311.4882 (spontaneous parity breaking, axion).
* 1306.5486 (dynamical Chern-Simons, charged).
* Jensen-Kaminski-Kovtun-Meyer-Ritz-Yarom 1112.4498 (hydrodynamics; holographic eta~ = 0).
* Hoyos, Peña-Benitez, Witkowski 1906.04752. They generate eta_H only with a magnetic field, through curvature-F couplings, and say explicitly that "the Hall viscosity in existing holographic models ... is likely to vanish, although this question has not been studied systematically".
* Web searches 2015-2026 for parity-odd cubic or quartic gravity together with Hall viscosity, and for Hall viscosity together with parity-odd ⟨TTT⟩.

**Not found:**
1. eta_H from axion-free `C C C~` or `C^2 C C~` on a neutral brane.
2. Any relation between eta_H and the parity-odd ⟨TTT⟩ coefficient alpha_T.
3. Any gap-dependent bound on eta_H.

All three appear new. The result `4 pi eta_H/s = ∓(3/40) alpha_T` at linear order is the Hall analogue of Li-Lü-Mai's `-(3/40) t4` for eta/s. That shear analogue is known; the Hall version is not, as far as I can find.

## 7. Verdict (item 6)
* **Correct?** Yes, with one correction.
  * The claimed coefficients 108, -432, 0, 0 and 6 theta'(1) are reproduced by an independent explicit-index computation, with reciprocity checked and consistency with Saremi-Son.
  * The quartic positivity map has a factor-2 slip inherited from CHLPS's printed ghat4. It makes their quartic bound loose, not wrong, and has no effect on the numbers.
  * The `||V_cubic|| = (3/40) tau` identity is exact given the paper's `tau = 240|ghat3|`. Its alpha_T half rests on the [D]-graded dictionary entry.
* **Caveats are identical to the eta/s case** (`../README.md` §5):
  * bulk-EFT power counting for higher-point couplings (I Q, J P, ... are invisible in TTTT but contribute);
  * leading order in 1/c_T;
  * no light single traces;
  * eta_H at c_T -> infinity.
  * A further one: eta_H requires microscopic P and T breaking, and in every known holographic CFT3 with a large gap, alpha_T is zero or exponentially small. ABJM is parity invariant; in ABJ the parity violation is a discrete holonomy effect. So the relation is a prediction with no current test case.
* **PRL strength.** Better than eta/s alone. "One complex coupling controls the energy-flux anisotropy, the shear viscosity and the Hall viscosity, `(4pi eta/s - 1) - i 4pi eta_H/s = 18 ghat3`, and a higher-spin gap bounds all three" is a clean, new, cond-mat-legible statement. The Hall-viscosity–⟨TTT⟩_odd link appears to be new physics, not just a corollary.
  * It is still a corollary of the tau bound plus an EFT assumption, so it belongs in the letter as a highlighted consequence (one paragraph plus one equation), not as a replacement headline.
  * Hostile referee points to pre-empt: (i) power counting for higher-point couplings is assumed, not derived; (ii) no known example has alpha_T != 0 with a large gap; (iii) quote 2.17e-9 and 1.73e-10, not the rounded-up values.
