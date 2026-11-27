# All-orders shear and Hall response: verification of the collaborator's claims

Status: 27 Nov 2026, research agent, not refereed. No edits to `prl/`, nothing committed.

Grades:
* **[computed]**: by a script here.
* **[cited]**: from the named source.
* **[I]**: inferred, with the argument given.
* **[failed]**: claim not confirmed or wrong.

Conventions: `L_AdS = 1`, `r_h = 1`, `f = r^2 - 1/r`, `eps_vrxy = +r^2`, `V = (4pi eta/s - 1) - i 4pi eta_H/s`, so `dV_H = 4pi eta_H/s`. Responses are per unit coupling of the operator added to `R + 6`.

| script | purpose |
|---|---|
| `static_engine.py` | Static diagonal-gauge shear engine. Handles arbitrary N, f, theta-profiles and phi''' jets; generic membrane flux for `sum_k alpha_k (phi^(k))^2`. |
| `shear_theta.py` | Shear response of `theta(r) I` and `theta(r) J` for a **generic** profile theta (local jets at the horizon). |
| `ef_engine.py`, `hall_theta.py` | EF engine (core of `../hall/hall_ef.py`). Hall response of `theta(r) P`, `theta(r) Q` for generic theta; integration-by-parts locality test. |
| `hall_gradient_direct.py` | `grad I . grad P` computed directly in EF (third-order jets), with no hand integration by parts. |
| `gradient_direct.py` | Direct flux for `(grad I)^2`, `(grad J)^2`, `U = (grad C)^2`, `I U`. Also the background check of `U = 1/2 box I + 6 I + 3 J`, and the first-law route (slow). |
| `boxC_wald.py` | `box C` on the brane; Iyer-Wald entropy of U and I U, with no identity used. |
| `remainder.py` | Exact rational evaluation of E(x), H(x); brute-force majorant check; exact gradient-tower sums. |
| `second_order_shear.py`, `wald_J_corrected.py`, `einstein_alpha_generic.py`, `hall_second_order.py`, `ef_engine_lam.py` | O(coupling^2) terms: shear and mixed Hall |

Recorded outputs: `out_*.txt`.

---

## Item 1. Algebraic F(I, J, P, Q): CONFIRMED [computed]

**Why only the r-profiles of F_y matter.** The shear perturbation has `y^(1) = 0` for all four invariants, in both the static gauge and the time-dependent EF gauge (checked). The O(eps^2) action of F is therefore `sum_i F_{y_i}(r) y_i^(2)`. Only the profiles `theta_i(r) = F_{y_i}` on the brane enter.

**Generic-profile results.**
* **Shear** (`shear_theta.py`, generic theta with local horizon jets):
  * `theta I`: `eta1/eta0 = 4 theta(1) - 6 theta'(1)`.
  * `theta J`: `eta1/eta0 = -102 theta(1) + 18 theta'(1)`.
  * Wald: `s1/s0 = 4 theta(1)` and `6 theta(1)`.
  * Hence `dV_E = -6 theta_I'(1) - 108 theta_J(1) + 18 theta_J'(1)`.
* **Hall** (`hall_theta.py`, generic theta, exact integration by parts):
  * `sigma_P = -6 theta''`.
  * `sigma_Q = 18 (r^2 theta'' - 9 r theta' + 24 theta)/r^5 = d/dr[ 18 (r theta' - 6 theta)/r^4 ]`. This is an exact total derivative with **no bulk remainder**.
  * Hence `dV_H = 6 theta_P'(1) + 108 theta_Q(1) - 18 theta_Q'(1)`. This assumes the profiles do not grow at infinity, which holds for F_y built from decaying invariants.

**The collaborator's formulas follow.** Substitute `theta_i' = F_{iI} I' + F_{iJ} J'`, with `I' = -72` and `J' = -108` at the horizon (P' = Q' = 0):
* `dV_E = 432 F_II - 648 F_IJ - 1944 F_JJ - 108 F_J`.
* `dV_H = 108 F_Q + 1296 F_QI + 1944 F_QJ - 432 F_PI - 648 F_PJ`.

No `F_PP`, `F_PQ`, `F_QQ` or `F_P` term can appear.

**Monomial tests.** All agree with the formulas and with my earlier independent computations:

| monomial | channel | value |
|---|---|---|
| I Q | Hall | 2592 |
| J P | Hall | -648 |
| J Q | Hall | 3240 |
| I^2 P | Hall | -10368 |
| I J | shear | -1944 |
| I^3 | shear | 31104 |
| J^2 | shear | -6480 |
| P^2, I, P | both | 0 |

## Item 2. Scalar-gradient extension F(y, X): CONFIRMED, with one caveat [computed + I]

**The reduction argument is correct.** Its premises, all checked:
* `y^(1) = 0`.
* The shear lives in the xy block, so `g^{vv}, g^{vr}, g^{rr}` and `sqrt(-g)` are exactly unchanged (`det e^{eps H} = 1`).
* Background y depends on r only.

Hence `X_ij^(2) = 2 grad y_(i)^(0) . grad y_(j)^(2)` exactly. A bulk integration by parts (no effect on the EOM; boundary terms at infinity decay) gives `F^(2) ≅ sum_i A_i(r) y_i^(2)`. Item 1's generic-profile formulas then apply with `theta_i = A_i`.

**Wald part.** Iyer-Wald gives `E^{abcd} = sum_i A_i dy_i/dR_abcd - sum_i (dF/d grad_e y_i) grad_e(dy_i/dR)`. The second term vanishes at the bifurcation surface, because `grad y = 0` there (`r - r_h ∝ UV` in Kruskal coordinates). Hence `ds = 4 A_I + 6 A_J` [I].

**Direct checks** (no hand integration by parts):

| operator | direct computation | claim / formula |
|---|---|---|
| `(grad I)^2`, shear | flux `e1 = 32832`; Iyer-Wald `w = 4 A_I(1) = 1728` → `dV_E = 31104` | 31104, confirmed |
| `(grad J)^2`, shear | `e1 = -276048`, `w = 6 A_J(1) = 3888` → `-279936` | formula -279936 |
| `grad I . grad P`, Hall | direct EF with third-order jets: `sigma = 7776 (7r^3 - 30)/r^11`, integral **-15552** | -15552, confirmed |

* **The first-law cross-check did not finish.** The general-(N, f) background solve was too slow and was stopped after 40 min. The Wald handling of derivative operators is instead checked two ways:
  * **IU and U, without field redefinitions:** Iyer-Wald plus flux (item 3) agree with the on-shell identity route (U: -324 both ways).
  * **Identity route for the claimed values:** `(grad I)^2 ≅ -2 I U + 12 I^2 + 6 I J` gives `-2(-16200) + 12(864) + 6(-1944) = 31104`, and `grad I . grad P ≅ -2 P U + 12 I P + 6 P J` gives `-2(3240) + 12(-432) + 6(-648) = -15552`.
  * Both equal the A-formula values. The second also equals the direct EF computation.
* **Caveat.** The formulas need the A_i profiles to decay at infinity. They do for any polynomial F in (y, X).

## Item 3. Ten-derivative operators: CONFIRMED [computed]

**Background identities.**
* `U0 = 180 f/r^8 = 180 (r^3 - 1)/r^9`.
* `box I0 = 216/r^6 - 432/r^9`.
* `U0 - (1/2 box I0 + 6 I0 + 3 J0) = 0`.
* `C . box C = -6 I - 3 J` on the brane.
* `(box C)_trtr / C_trtr = -6 - 3/r^3`.

On Einstein spaces the Penrose wave equation makes `C . box C` a combination `a Lambda I + b J`. The independent r-dependences of I0 and J0 on the brane fix `a Lambda = -6` and `b = -3`. So `U = 1/2 box I - C . box C = 1/2 box I + 6 I + 3 J` holds on any Einstein space with `Lambda = -3` [computed + I]. The Ũ identity follows by the chiral (self-dual) decomposition [I].

**IU, shear, with no identity used.**
* Flux on the brane: `e1 = -15768`.
* Iyer-Wald: `E = -2 I box C` at the bifurcation surface (where `U = grad I = 0`), giving `w = 432`.
* **dV_E = -16200. Confirmed.**
* U alone gives `-288 - 36 = -324`, which equals `3 x (-108)` from the identity.

**PU, Hall.** `(PU)^(2) = U0 P^(2)` exactly, so it is `theta P` with `theta = U0`: `6 U0'(1) = 6 x 180 x 3 =` **3240. Confirmed.**

**PŨ.** `P0 = P^(1) = Ũ0 = Ũ^(1) = 0`, so `(PŨ)^(2) = 0` identically. **0. Confirmed.**

**IJ, IQ, PJ, PQ.** -1944 (shear), 2592, -648, 0 (Hall). **Confirmed** (item 1).

**Count of independent structures.**

Ruhdorfer, Serra and Weiler (arXiv:1908.08050, eq. GravityHilbert) give, at mass dimension 10:

    C_L^5 + C_L^3 C_R^2 + C_L^2 C_R^3 + C_R^5 + C_L^4 D^2 + C_L^2 C_R^2 D^2 + C_R^4 D^2

That is 7 real structures. Their explicit basis is `IC, ĨC, IC̃, ĨC̃, FC, FC̃, F̃C̃`, which in our notation is **{J I, Q I, J P, Q P, U I, U P, Ũ P}**. This is exactly the collaborator's set [cited].

The Hilbert series has **no C^3 D^2 term at dimension 8 and no C^3 D^4 term at dimension 10**. So every `grad C grad C C` contraction is redundant (integration by parts, Bianchi, EOM). Every `grad C grad C C C` contraction reduces to the three structures `U I, U P, Ũ P`, plus quintic terms, plus lower-dimension terms. [cited + I]

Caveat: RSW count in flat space (R_mn = 0). With Lambda != 0, each reduction also produces `Lambda x (dimension-8 and dimension-6 operators)` with O(1) coefficients. The count of new structures is unchanged [I]. These shifts move the k = 4 and k = 3 couplings by O(C/Lambda^8). Through t4 (exact cubic data) they are absorbed; through the TTTT-measured quartic coupling (`g~4 = g4 + O(Lambda^-2)`) they appear at the remainder order.

**Coverage at k >= 6.** It is not complete. Two examples:
* `U^2 = (1/2 box I + 6I + 3J)^2` contains `(box I)^2`, which is not of `F(y, grad y . grad y)` form.
* The C^4 D^4 and C^3 D^6 classes, among others, are not covered.

## Item 4. Remainder: algebra CONFIRMED; "all orders" OVERSTATED [computed + I]

**Power counting is consistent.** A monomial `I^m J^n` has `4m + 6n` derivatives, so `k = 2m + 3n`. The hypothesis is `|gamma| <= C x^{k-1}` with `x = 1/Lambda^2`. E(x) and H(x) are generating functions with weight `x^k`. The remainder is `C [E]_{k>=5} / x`, so `[E]_{k>=5}/x^5 < 1945` means `R <= 1945 C/Lambda^8`. That is consistent.

**The generating functions are exact.** A brute-force sum over (m, n) of the absolute values of the four formula terms reproduces E(x) and H(x) **exactly** (relative difference 0, exact rationals). So they are not merely majorants: they equal the sums of the term-wise absolute responses. They do majorize the true `|response|` sums.

**Rational values** at `x = 1/997002`:
* `(E - 108x^3 - 864x^4)/x^5 = 1944.0377 < 1945`. **Confirmed.** It is dominated by IJ (-1944).
* `(H - 108x^3 - 432x^4)/x^5 = 3240.0136 < 3241`. **Confirmed.** IQ gives 2592 and PJ gives 648.

**Gradient tower: exact sums.**
* Monomials with X-degree >= 3 give exactly 0, because X vanishes at the horizon to that order.
* Monomials with X-degree 1 or 2 and k <= 10 (89 of them) give `sum_{k>=6} |dV| x^{k-6}` = 34993 (shear) and 104977 (Hall), both below 129000. **Confirmed.** The terms fall by factors of about 1e-4 to 1e-5 per order beyond that.

**Horizon jets.** `|y'| = (72, 108, 0, 0)`, `|y''| = (504, 1080, 0, 0)`, `X'_ij = f'(1) y_i' y_j' = 3 y_i' y_j'`. **Confirmed.** Only y', y'' and X' enter A and A' at the horizon, because f(1) = 0.

**Totals.** `R_E = 1944.04 (IJ) + 16200 (IU) + 0.13 (gradient k>=6) = 18144.2` against the claimed 18146 (rounded up), and `R_H = 3240.01 + 3240 (PU) + 0.1 = 6480.1` against the claimed 6482. Both **confirmed** as arithmetic. The norm, `sqrt(18146^2 + 6482^2) = 19269`, matches "19300".

**Hostile qualifications.**
1. **Basis dependence.** `(grad I)^2 ≅ -2 I U + 12 I^2 + 6 I J` and `grad I . grad P ≅ -2 P U + 12 I P + 6 P J` (modulo EOM and total derivatives). In the equally natural basis `{(grad I)^2, grad I . grad P, ...}` with the same C, the k = 5 constants become 1944 + 31104 = 33048 (shear) and 3240 + 15552 = 18792 (Hall). The hypothesis "|gamma_O| <= C/Lambda^{2k-2}" must name its basis. The constant is basis dependent by about a factor of 2; the order `C/Lambda^8` is not.
2. **Completeness.** "All orders" holds only for the class `F(I,J,P,Q, grad y . grad y)` plus the k = 5 U-terms. At k >= 6 the general EFT has operators outside this class (`(box I)^2`, `(grad grad C)^2 C^2`, ...). No majorant was given for them. This is uncovered. A real all-orders statement needs a Hilbert-series count per order times a uniform bound on the response per structure.
3. **The hypothesis itself** (gap power counting for couplings with five or more gravitons) is an assumption, not a CFT result. So is the AdS bookkeeping caveat of item 3.
4. **Size.** Natural C is 20-100: the paper's bounds give `|gamma3| Lambda^4 <= 20` and `gamma4 Lambda^6 <= 23`. Then `|R| <= 19300 x 100/9.88e23 = 2e-18`. That is irrelevant next to 2.17e-9, whatever basis caveat applies.

## Item 5. Nonlinear O(ghat3^2) terms [computed]

**Hall** (`hall_second_order.py`). Take kap Q on the O(lam)-corrected J brane, `f1 = -8/r + 24/r^4 - 16/r^7`, `n1 = -6/r^6`:
* `sigma_Q = 432/r^5 + lam 216 (32 r^6 - 462 r^3 + 1035)/r^11`.
* The integral is `108 + 9828 lam`. Dividing by the Wald factor `s_hat = 1 + 6 lam` gives
* **`4 pi eta_H/s = -18 alpha3~ + 255 alpha3 alpha3~ + O(3)`**, with `lam = -alpha3/6` and `kap = -alpha3~/6`.

**Shear** (`second_order_shear.py`, `wald_J_corrected.py`, `einstein_alpha_generic.py`).
* The Einstein kinetic coefficient is exactly `-N r^2 f/2` for **any** static N, f. So the Einstein flux and Einstein Wald entropy are background independent, and the O(lam^2) background is not needed.
* J-part horizon flux on the corrected brane: `102 + 9516 lam` (Einstein = -1).
* J-part Wald entropy: `4 pi (6 + 312 lam)`.
* Result: **`4 pi eta/s = 1 + 18 alpha3 - 255 alpha3^2 = 1 - (3/40) t4 - (17/3840) t4^2`**.
* This reproduces exactly the t4^2 coefficient of Li, Lü, Mai (arXiv:1901.03349 eq. (7.15)), whose `(C_T - a~)` term vanishes for pure C^3. That is an independent check of the second-order method.

**Structural results** [I, argued in the scripts' docstrings].
* `eta` has **no alpha3~^2 term** at tree level. Q is linear in kap; its quadratic action has only phi-chi cross terms; the induced phi is O(kap omega), so its phi_v feeds back only at O(omega^2). Q changes neither the parity-even background nor the Wald entropy.
* `eta_H` has no `alpha3^2` or `alpha3~^2` term (parity).
* `t4 = -240 alpha3` is exact at tree level for pure C^3, because W = 0 on AdS: C_T and L_AdS are unchanged.

**Combined result** (tree level, through second order in the cubic couplings):

    V = (4pi eta/s - 1) - i 4pi eta_H/s = ghat3 (18 - 255 Re ghat3) + 216 (alpha4 + i alpha4~) + O(|ghat3|^3, ...)
      = ghat3 (18 + (17/16) t4) + ...

The cubic sector stays exactly "complex-linear in ghat3" at second order, multiplied by a real factor.

**Explicit nonlinear remainder.**

    |V - 18 ghat3| <= 255 |Re ghat3| |ghat3| <= 255 |ghat3|^2 = (17/3840) tau^2

* At gap 1000 (`tau <= 2.8914e-8`): **3.7e-18**.
* Positivity-free (`tau_* = 0.0178`): **1.4e-6**.
* Mixed cubic-quartic terms, O(alpha3 alpha4) with size about 1e-10 x 1e-16, and O(alpha4^2) terms were **not computed**. Under the bounds they are below 1e-24.

---

## Summary table

| claim | verdict |
|---|---|
| 1. dV_E and dV_H algebraic formulas | **confirmed for generic profiles** [computed]; all monomial tests pass |
| 2. Gradient extension: reduction argument and formula | **confirmed** [computed + I]. Direct `(grad I)^2` = 31104 (flux + Iyer-Wald); direct `grad I . grad P` = -15552 (EF, no hand integration by parts). First-law cross-check not completed; identity-route cross-check passes |
| 3. U identities; IU, PU, PŨ, IJ, IQ, PJ, PQ | **confirmed** [computed]. IU = -16200 with no identity used. RSW count of 7 confirmed from their Hilbert series [cited]. `grad C grad C C` is redundant (no C^3 D^2 term); all `grad C grad C C C` contractions reduce to {UI, UP, ŨP} [cited + I] |
| 4. Generating functions, rational bounds, jets, totals | **arithmetic confirmed** [computed]: 1944.04 < 1945; 3240.01 < 3241; gradient k>=6 sums 34993 and 104977 < 129000; totals 18144 and 6480 (claims 18146 and 6482 are rounded up). **Overstated**: basis dependent (e.g. 33048 and 18792 in the grad-y basis); "all orders" covers only the F(y, grad y . grad y) class beyond k = 5 [failed as an all-operator theorem] |
| 5. Nonlinear terms | **new** [computed]: `4pi eta/s = 1 + 18 alpha3 - 255 alpha3^2`; `4pi eta_H/s = -18 alpha3~ + 255 alpha3 alpha3~`; no `alpha3~^2` in eta; `V = ghat3 (18 - 255 Re ghat3)`; `|V - 18 ghat3| <= (17/3840) tau^2 = 3.7e-18` at gap 1000 |

**What remains uncovered.**
* Operators outside `F(y, grad y . grad y)` at 12 or more derivatives, e.g. `(box I)^2`, `(grad grad C)^2 C^2`, and general C^4 D^4. No uniform majorant was given.
* The power-counting hypothesis itself (couplings with five or more gravitons).
* The AdS bookkeeping of Lambda x (lower-dimension) shifts produced when dimension-10 operators are reduced. Physical t4 absorbs the cubic shifts; the TTTT-measured g4 absorbs the quartic shifts only up to `O(Lambda^-2)`.
* Mixed cubic-quartic nonlinear terms.
* 1/c_T corrections.

None of these can matter numerically next to 2.17e-9 unless a coupling exceeds natural size by about 10^9.
