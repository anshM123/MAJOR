# Scope: can the TTTT coupling bounds be turned into a bound on eta/s (or another thermal observable)?

Status: written 27 Nov 2026 by a research agent. Not refereed. Nothing in `prl/` was edited.
Marks: **[computed]** derived here by a script in this directory. **[cited]** quoted from the source named.
**[I]** inferred, with the argument given.

Scripts (sympy and numpy; run each with `python3 <script>`; recorded outputs are in `out_*.txt`):

| script | what it does | runtime |
|---|---|---|
| `etas_diag.py` | linear-in-coupling `4 pi eta/s - 1` and Wald `s` for `L = R + 6 + lam O` on the AdS4 black brane, for O in {I2, I3, (W W~)^2, W W W~, I2^2, I2 I3, I2^3, I3^2, I2 (W W~)} | ~40 s |
| `thermo_check.py` | O(lam) background, temperature shift and on-shell free energy for I3 and I2^2; checks the first law against the Wald entropies used above | ~15 s |
| `quartic_flat.py` | D=4 identity `W_GW = (1/32)[(W^2)^2 + (W W~)^2]`; four-graviton contact amplitudes, giving the map to CHLPS `g4` | <1 s |
| `numbers.py` | all X(Delta_g) numbers below | <1 s |

---

## 1. Result in one line

At tree level (leading order in 1/c_T), with `L_AdS = 1`, `t4 = -240 Re ghat3`, and CHLPS normalisation of the quartic terms
`(1/16piG)(1/4)[alpha4 (R^2)^2 + alpha4' (R R~)^2 + 2 alpha4~ R^2 R R~]`, `g4/8piG = alpha4 + alpha4'`:

    4 pi eta/s - 1 = -(3/40) t4 + 216 alpha4 + O(t4^2, tau^2, Lambda^-8)
                   = 18 Re ghat3 + 216 alpha4 + ...

The parity-odd couplings (`alpha_T`, the W W W~ term, and `alpha4~`) and `alpha4'` drop out at linear order.

Numbers at `Delta_g = 1000`:
* cubic term `<= (3/40) tau <= 2.17e-9`, from the certified bound `tau <= 2.8914e-8`;
* quartic term `<= 216 g4/(8 pi G) <= 2.0e-14` (needs `0 <= alpha4 <= g4/8piG`, which uses flat-space positivity `alpha4' >= 0`, plus the conditional quartic bound 92.8);
* positivity-free version: `(3/40) tau_* = 4/(3(Delta_g - 2)) [1 - 49/(2Delta_g+9)^2]^{-1/2} = 1.34e-3`.

All of this holds only under the power-counting assumption of Sec. 4: the higher-point couplings must be as small as the gap implies. Vacuum `<TTTT>` data do not fix eta/s.

---

## 2. Coefficients

### 2.1 Table: `L = R + 6 + lam O` on the planar AdS4-Schwarzschild brane (`f = r^2 - 1/r`, horizon at r = 1) [computed, `etas_diag.py`]

At the horizon: `I2 = W_abcd W^abcd = 12`, `I3 = W_ab^cd W_cd^ef W_ef^ab = 12`, `dI2/dr = -72`, `dI3/dr = -108`; `W W~ = W W W~ = 0`.

| O | eta_1/eta_0 | s_1/s_0 (Wald) | 4 pi eta/s - 1 (per unit lam) |
|---|---|---|---|
| I2 = W^2 (Gauss-Bonnet + Ricci terms in 4D) | 4 | 4 | **0** |
| I3 (Weyl^3) | -102 | 6 | **-108** |
| W W W~ (parity odd, alpha_T) | 0 | 0 | **0** |
| I2^2 = (W^2)^2 | 960 | 96 | **864** |
| (W W~)^2 | 0 | 0 | **0** |
| I2 (W W~) (parity odd quartic) | 0 | 0 | **0** |
| I2 I3 (quintic) | -1824 | 120 | **-1944** |
| I3^2 (sextic) | -6336 | 144 | **-6480** |
| I2^3 (sextic) | 32832 | 1728 | **31104** |

**Method.** Only the correction operator's contribution on the Einstein brane is needed at O(lam). With the horizon held fixed, the Einstein term contributes `r_h^2/16piG` to eta (membrane value) and `r_h^2/4G` to s on *any* static background of this form, so the O(lam) background correction cancels in eta/s.

* eta: take the static shear mode `r^2 diag(e^{eps phi}, e^{-eps phi})`. Its radial quadratic Lagrangian reduces, after integrations by parts, to `alpha phi''^2 + beta phi'^2`; the mass term vanishes, as checked. The ingoing solution is `phi = e^{-iw(t+r_*)}(1 + i w g)` with g regular. The O(w) radial flux is then `w Pi[-r_*]` at the horizon, with `Pi = 2 beta chi' - (2 alpha chi'')'`. This step needs `alpha ~ (r-1)^2`, which was checked. The overall normalisation is fixed by Einstein gravity.
* s: Wald entropy, `E^{abcd} = dL/dR_abcd`, with the Riemann tensor treated as independent.

**Cross-checks.**
1. **I3.** `-108 lam`. The CHLPS action term is `-(alpha3/6) I3`, so `lam = -alpha3/6` and the correction is `+18 alpha3`. With `t4 = -240 alpha3` (`DICTIONARY_T4.md`, checked there against ECG and QTG) this is **`-(3/40) t4`**. That equals the D=4 coefficient in Li, Lü, Mai, arXiv:1901.03349, eq. (7.15): `eta/s = (1/4pi)(1 - (3/40) t4 - (45/2)((C_T - a~)/C_T)^2 - (17/3840) t4^2 + (9/320) t4 (C_T-a~)/C_T)` **[cited]**. It also matches the small-mu expansion `1 + 189mu/2` for Einsteinian cubic gravity (Bueno, Cano, Ruipérez, arXiv:1802.00018, eq. after (eta/s); ECG has `alpha3 = (21/4)mu`, so `18 alpha3 = 94.5 mu`) **[cited]**. BCR call the exact ECG result non-analytic in mu. That is a property of the resummed ECG brane; the linear coefficient agrees.
2. **I2^2.** `864 lam`. Garousi and Ghodsi, arXiv:0808.0411 §3.2, compute M-theory on AdS4 x S7 with the quartic Weyl term `gamma W_GW` built from the 4D Weyl tensor and get `eta/s = (1/4pi)(1 + 1728 gamma/L^6)`, where L is the S7 radius, equal to 2 L_AdS4 **[cited]**. Using `W_GW = (1/32)[(W^2)^2 + (W W~)^2]` [computed] and `(W W~)^2 -> 0`, my result gives `864 gamma/32 = 27 gamma/L_AdS^6 = 1728 gamma/(2 L_AdS)^6`, which agrees. Their temperature shift `T = (3 r0/2pi)(1 + 80 gamma)` also matches my `T = T0 (1 + 40 lam)` with `lam = gamma/32` in `L_AdS = 1/2` units (`thermo_check.py`).
3. **First law.** `w - 2 tau1 = k` for I3 and I2^2, where w is the Wald coefficient, tau1 the temperature shift and k the on-shell free energy. So the Wald entropies are consistent with the Euclidean action [computed].
4. **Sign.** Type II strings and M-theory have positive quartic corrections (Buchel, Liu, Starinets; Myers, Paulos, Sinha 0806.2156; Garousi-Ghodsi), and `alpha4 > 0` **raises** eta/s, consistent with the + sign.

### 2.2 Answers to the specific literature questions

* **Weyl^3 (t4):** nonzero at linear order, `-(3/40) t4`. The correction is nonzero because the black brane has Weyl curvature of order one in AdS units at the horizon (`W^2 = 12/L^4`, independent of T). W^3 therefore changes both the shear-mode action and the Wald entropy (-102 against +6). This is not an AdS5-type accident. In D=4 there is one independent parity-even cubic Weyl invariant, because `I2(ECG contraction) = I1/2` (`DICTIONARY_T4.md` §2).
* **W^2 W~ (alpha_T):** exactly zero at linear order [computed]. This follows from parity: the brane and the h_xy shear mode are invariant under the orientation-reversing swap x <-> y, so every pseudoscalar vanishes on them. alpha_T enters eta/s only at O(alpha_T^2) <= O(tau^2) ~ 1e-15. alpha_T would instead enter the **Hall viscosity** eta_H linearly (Sec. 6). That was not computed.
* **R^4:** `+216 alpha4` (CHLPS normalisation), `alpha4'`-independent. For the supersymmetric combination (`alpha4 = alpha4'`): `+108 g4/(8piG L^6)`. Type II normalisation `g4 M^6/8piG = 2 zeta(3)` gives `+108 * 2zeta(3)/(ML)^6`.
* **Myers-Sachdev-Singh** (arXiv:1010.0443): the Weyl coupling gamma there is the `C_abcd F^ab F^cd` term, i.e. `<JJT>`. It controls the conductivity `sigma(0)/sigma(inf) = 1 + 4 gamma`, not `<TTT>`. The paper's bounds say nothing about it.
* **Buchel-Liu-Starinets, AdS5 R^4** (hep-th/0406264): `eta/s = (1/4pi)(1 + 135 gamma)`, `gamma = zeta(3) lambda^{-3/2}/8`. This was revised to `1 + 120 gamma = 1 + 15 zeta(3) lambda^{-3/2}` once the five-form completion of the C^4 term was included (Myers-Paulos-Sinha 0806.2156; Buchel 0805.2683) **[cited, background only; d=4, not used]**.

### 2.3 Normalisation chain (convention-safe statement)

* `ghat3 = alpha3 + i alpha3~` with action `-(1/3!)(alpha3 R^(3) + alpha3~ R~^(3))` (CHLPS 2201.06602 eq. (2.11)).
* `t4 = -240 alpha3`; `tau = 240 |ghat3|`.
* The coupling of I3 is `gamma3 = -alpha3/6`, so the correction is `-108 gamma3 = 18 alpha3 = -(3/40) t4`. The collaborator's `gamma3 = -Re ghat3/6` is confirmed.
* `|4 pi eta/s - 1|_cubic = 18 |Re ghat3| <= 18 |ghat3| = (3/40) tau`.
* Quartic: `(1/16piG)(c1 (W^2)^2 + c2 (WW~)^2)` gives the MHV contact `g4/8piG = 4(c1+c2)` and the all-plus coefficient `prop. to (c1 - c2)` (`quartic_flat.py`, `out_quartic_flat.txt`). Hence `c1 = alpha4/4`, `c2 = alpha4'/4`, which reproduces CHLPS's `g4 = 8piG(alpha4 + alpha4')` independently. The I2^2 coefficient is `gamma4 = c1 = alpha4/4`, so the correction is `864 gamma4 = 216 alpha4`.
* `R` versus `W` in AdS4: on an Einstein space `R^(3) = W^(3) - (6/L^2) W^(2) + const` and `(R^2)^2 = (W^2)^2 + (48/L^4) W^2 + const`. The extra `W^2` pieces are Gauss-Bonnet plus Ricci terms and drop out (row I2 of the table).

---

## 3. Numbers [computed, `numbers.py`, `out_numbers.txt`]

`Lambda = [(Delta_g-2)(Delta_g-1)]^{1/2} = 998.50` at `Delta_g = 1000`.

| input | X = bound on |4 pi eta/s - 1| |
|---|---|
| cubic, certified `tau <= 2.8914e-8` (no light single trace in T x T) | **2.17e-9** |
| quartic, `216 g4/8piG`, `g~4 Lambda^6/8piG <= 92.8` (positivity assumed; no scalar < 0.46 Lambda) | 2.0e-14 (3.6e-14 with light scalars, 164.5) |
| quartic at the type II value `2 zeta(3)/Lambda^6`, supersymmetric combination | 2.6e-16 |
| `(17/3840) t4^2` | < 4e-18 |
| **total (conditional)** | **2.17e-9** |
| positivity-free: `(3/40) tau_*(1000)` | **1.336e-3** |

Positivity-free bound as a function of the gap (exact algebra): `X_*(Delta_g) = (3/40) tau_* = 4/(3(Delta_g - 2)) [1 - 49/(2Delta_g + 9)^2]^{-1/2}`.

| Delta_g | 10 | 30 | 100 | 300 | 1000 | 1e4 | 1e5 | 1e6 |
|---|---|---|---|---|---|---|---|---|
| X_* | 0.172 | 0.0479 | 0.0136 | 4.47e-3 | 1.34e-3 | 1.33e-4 | 1.33e-5 | 1.33e-6 |

Other bounds, valid for Delta_g >= 1000 only:
* analytic all-gap functionals (one further convergence input): `X <= 21 [log((Delta_g+2)/3)]^{1/2}/(Delta_g-2)^2`, which is 5.1e-5 at 1000;
* with the k=3 input: `31.5 (log Delta_g)^{1/2}/(Delta_g-2)^{5/2}`, which is 2.6e-6 at 1000;
* positivity at each gap's own value, using the faded curves: `X = 18 (|ghat3|Lambda^4)_max / Lambda^4`, which gives 2.17e-9 (1000), 1.4e-10 (2000) and 4.0e-12 (5000). Positivity is established only at 1000.

**Reverse reading.** A measured deviation `delta = |4 pi eta/s - 1|` needs `tau >= (40/3) delta`. That requires a higher-spin single trace in T x T below twist `g_*((40/3)delta) ~ 2 + 1.335/delta` (from `g_*(tau) ~ 2 + 17.8/tau`). Examples: delta = 10% needs Delta_g <~ 15; delta = 1% needs Delta_g <~ 136. This uses the positivity-free inputs and the power counting of Sec. 4, and is meaningful only where the quartic and higher terms are small, i.e. Delta_g >> 1.

**Alternative observable: the thermal free energy** [computed, `thermo_check.py`]. At fixed T, with C_T and L_AdS unchanged by any W^{n>=3} term, `f/f_Einstein = 1 + 2 gamma3 + 16 gamma4 + ... = 1 + t4/720 + 4 alpha4 + ...`. This gives `|f/(f_E) - 1| <= tau/720 = 4.0e-11` (conditional) and `2.5e-5` (positivity-free) at Delta_g = 1000. The coefficient is 54 times smaller than for eta/s, and the observable is less famous.

---

## 4. Check of the collaborator's eta/s derivation

| claim | verdict |
|---|---|
| horizon values `I2 = I3 = 12`, `I2' = -72`, `I3' = -108` | **confirmed** [computed] |
| `delta eta = 4F2 - 102F3 - 6 dF2 + 18 dF3` | **confirmed term by term**: my eta_1/eta_0 for I2, I3, I2^2, I2I3, I2^3, I3^2 = 4, -102, 960, -1824, 32832, -6336, all reproduced by this formula |
| `delta s = 4F2 + 6F3` (Wald) | **confirmed** (4, 6, 96, 120, 1728, 144); first law checked independently for I3 and I2^2 |
| `4 pi eta/s - 1 = 432F22 - 648F23 - 1944F33 - 108F3` | **confirmed** (algebra from the two lines above; no F2 term); I3^2 gives -6480 from both |
| per-term table `0, -108, 864, -1944, 31104` | **confirmed** by an independent method. I used a static radial-momentum/flux computation rather than the time-dependent horizon-pole method. Their intermediate `delta I2`, `delta I3` (time-derivative form) were not re-derived |
| agreement with Li et al. 1901.03349 (7.15) | **confirmed** `-(3/40)t4`; also agrees with BCR ECG `189mu/2` |
| `gamma3 = -Re ghat3/6`, cubic part `18|ghat3| = (3/40)tau` | **confirmed** (strictly `18|Re ghat3| <= 18|ghat3|`) |
| `gamma4 = alpha4/4`, `g4/8piG = alpha4 + alpha4'` | **confirmed** independently by the contact-amplitude computation (`quartic_flat.py`); consistent with CHLPS (2.11)/(eq: g4hat Lagrangian) |
| `alpha4, alpha4' >= 0` from flat forward positivity | **correct as a flat-space statement** (CHLPS §3: `|ghat4| <= g4`, so `alpha4, alpha4' >= 0`, up to loop effects `~ G^2/(M^2 p_*^2)`). It is **not** established in AdS by this paper, which bounds `g~4` from above only and conditionally. Numerically irrelevant (Sec. 3) |
| `|4 pi eta/s - 1|_lin <= 2.175e-9 + 2.0e-14` | **confirmed**: 2.1685e-9 + 2.02e-14 |
| `I2 I3` is quintic around AdS, invisible in vacuum TTTT, yet shifts eta/s | **confirmed and important** (Sec. 4.1) |

**Is F(I2,I3) general enough?**
* **6 derivatives.** On-shell, modulo field redefinitions and topological terms, the general cubic is `gamma3 I3 + gamma3~ (W W W~)`. `∇W∇W` and `W box W` reduce on Einstein spaces (Penrose wave equation `box W = Lambda W + W.W`) to `Lambda W^2 + W^3`. The parity-odd term drops out of eta/s [computed]. So F(I2,I3) is complete [I].
* **8 derivatives.** Two parity-even quartic invariants exist in 4D: `(W^2)^2` and `(WW~)^2`. F(I2,I3) omits `(WW~)^2` and the parity-odd `W^2 WW~`, but both contribute 0 [computed]. 8-derivative operators with derivatives (`∇^2 W^3`, `∇W∇W W`) have no new on-shell three-point amplitudes in 4D. They reduce to field redefinitions + `Lambda x` (cubic, absorbed into the physical t4) + quartic Weyl terms, plus remainders with more fields and more derivatives. So F is complete through O(Lambda^-6) [I].
* **10 derivatives.** F(I2,I3) misses quartic derivative operators (`∇^2 W^4`-type, the `g5` of CHLPS, which is in TTTT but not bounded here) and parity-odd quintics. These are O(Lambda^-8) and numerically irrelevant under power counting.

### 4.1 The obstruction: vacuum TTTT does not determine eta/s

True, and general. eta/s probes the bulk Lagrangian at the horizon, where `W^2 L^4 = 12`. That is a point of order one in curvature space, independent of T by conformal invariance. Every `W^n` operator contributes there (the table: I2I3, I2^3, I3^2 all nonzero), but around pure AdS (`W = 0`) an n-curvature operator starts at n-point order. `<TTT>` and `<TTTT>` see only n <= 4. A CFT-only statement about eta/s would therefore need bounds on all higher-point stress-tensor correlators, or a thermal-state argument such as causality of `<TT>` in the thermal state. Neither is available.

**What power counting buys.** If every 2k-derivative coupling is `c_k/Lambda^{2k-2}` in AdS units with `c_k` of order one, the expansion parameter at the horizon is about `(numbers ~ 10-30) x c / Lambda^2 ~ 1e-5 c`. The calibration is that the paper's bounds allow `|gamma3| Lambda^4 <= 20` and `gamma4 Lambda^6 <= 23`, and type II has `g4 M^6/8piG = 2.4`.
* The first ambiguity is O(Lambda^-8): `-1944 gamma5` with `gamma5 = c5/Lambda^8`, i.e. `2.0e-21 c5`, which is `<= 2e-19` for `c5 <= 100`.
* The next is `3e-26 c6` (`I2^3`) and `7e-27 c` (`I3^2`).
* The remainder is therefore about 1e-10 of the cubic term, and the bound `2.2e-9` is unaffected unless some quintic coefficient is about 10^10 times its natural size.

The assumption is natural but cannot be removed: it is an assumption about infinitely many couplings that the paper's sum rules do not see.

---

## 5. Assumptions needed to state |4 pi eta/s - 1| <= X (honest list)

1. **Everything the coupling bound needs**: (i) higher-spin gap; (ii) tree-level bulk EFT below the gap, uniformly in spin; (iii) positivity (certified at 1000 for the cubic functional); no light single trace in T x T for the 2.89e-8 number; the detector-identity inputs for tau_*.
2. **A local bulk EFT valid on the black-brane background, horizon included**, with a single cutoff `M = Lambda/L`. In a CFT on R^2, T is not a scale: the horizon curvature is `W ~ 1/L^2` at every T. So "T << gap" is not the condition. The condition is `L^2 W << (ML)^2`, i.e. `Delta_g >> 1`, plus the power counting of higher-point couplings (Sec. 4.1). This is the new assumption beyond the paper.
3. **Leading order in 1/c_T.** Loop corrections are not controlled. In an EFT with cutoff M they can be UV-enhanced, of order `Delta_g^k/c_T` (e.g. M-theory: the R^4 term is itself a one-loop effect, with a coefficient `~ c_T^{-2/3}`). To reach the 1e-9 level in a real theory you need `c_T` astronomically large (`c_T >~ 1e9 x Delta_g^k`). In 2+1 dimensions, hydrodynamic fluctuations also make the zero-frequency shear viscosity log-divergent at any finite c_T (long-time tails). So "eta" means the c_T -> infinity transport coefficient, taken before omega -> 0. This was not checked quantitatively here; it is a standard 2D-hydrodynamics fact (classic Pomeau-Résibois / Forster-Nelson-Stephen; see Kovtun, Moore, Romatschke, arXiv:1104.1586, for the 3+1 analogue).
4. **No light scalar (or spin-2) single trace coupling to TT.** A light scalar with a `phi W^2` coupling is sourced on the brane (`W^2 != 0`) and shifts eta/s at O(lambda_{TTphi}^2). That is tree level, so not 1/c_T suppressed, and it is not bounded by the cubic functional. Meltzer-Perlmutter, arXiv:1712.04861, already remark qualitatively that sparse-spectrum bounds on `<TTO>` limit such transport effects. The positivity-free tau_* allows light scalars (except Delta <= 1), so the X_* statement additionally needs either none, or their own bound.
5. **Planar black brane is the thermal saddle**, which is automatic on R^2; there is no Hawking-Page transition. The bulk includes no other light fields that condense.
6. **Quartic sign input**: `0 <= alpha4 <= g4/8piG`, which uses flat-space `alpha4' >= 0`. It matters only at the 1e-14 level.

---

## 6. Novelty check

Searched: arXiv sources of CEMZ 1407.5597, AKT 1805.07393, Kundu 2104.11238, AHKT 1610.09378, Meltzer-Perlmutter 1712.04861, CLPS 2201.06602, CLPS-higher-D 2205.01495, CLSD 2311.04271, BMS 0812.2521, plus web searches.

* **Bounds on eta/s from conformal-collider bounds on higher-derivative couplings: old and well known.**
  * Gauss-Bonnet in d=4: eta/s >= (16/25)(1/4pi) (Brigante-Liu-Myers-Shenker-Yaida 0712.0805, 0802.3318; Hofman 0907.1625; Buchel-Myers 0906.2922).
  * Arbitrary d: Buchel-Escobedo-Myers-Paulos-Sinha-Smolkin 0911.4257; Camanho-Edelstein 0911.3160.
  * Quasi-topological gravity: eta/s >= 0.414/4pi (Myers-Paulos-Sinha 1004.2055).
  * **d=3 ECG**: eta/s <= 1.253/4pi from t4 >= -4 (Bueno-Cano-Ruipérez 1802.00018).
  * These are model-specific and O(1). None is gap dependent.
* **The linear map eta/s <-> t4 in D=4 is known**: Li-Lü-Mai 1901.03349 (7.15), `-(3/40) t4` for generic massless cubic gravity. The R^4 piece for the supersymmetric combination in AdS4 is Garousi-Ghodsi 0808.0411. The separation into `(W^2)^2 -> 864` and `(WW~)^2 -> 0`, and the general `F(I2,I3)` formula, I did not find written down. They are easy.
* **Gap dependence**: CEMZ and AKT give parametric suppression (`t4 ~ 1/Delta^4`, AKT §7 d=3 heuristic `|t4| <~ ln Delta/Delta^4`), with no constants and no eta/s statement (`viscos` appears only in their bibliographies). AHKT (1610.09378, intro) explicitly say there is **no universal CFT derivation of eta/s = 1/4pi**. CLPS-higher-D's `|a-c|/c <= 23/Delta^2` combined with Kats-Petrov/BMS `4pi eta/s - 1 ~ -(c-a)/c` would give a d=4 analogue. I found no paper that states it.
* **Dispersive-sum-rule papers 2021-2026** (CHMRS, CLPS, CLSD, Kundu): none mention viscosity (checked sources).

**What would be new**: the first explicit-constant, gap-dependent bound on eta/s for holographic CFTs, in any dimension as far as I can find. It is in d=3/AdS4, with `X = 2.2e-9` at `Delta_g >= 1000` (conditional) and `X_* = 4/(3(Delta_g-2))` at every gap (positivity-free), and a reverse reading (`delta` needs `Delta_g <~ 2 + 1.34/delta`). **What is not new**: the map from couplings to eta/s, which is one line from Li-Lü-Mai plus Garousi-Ghodsi. The new content is entirely the paper's tau bound.

---

## 7. Verdict (hostile)

* **Correct?** Yes, at linear order and under the listed assumptions. The coefficient `-(3/40) t4` was confirmed three ways: my computation, Li-Lü-Mai, and BCR ECG. The R^4 coefficient was confirmed against Garousi-Ghodsi.
* **Non-vacuous?** Yes. The universal collider bound `|t4| <= 4` gives only `|delta| <= 0.3` at linear order, and is meaningless beyond it. Here `X = 2.2e-9` at 1000, and the positivity-free `X_* = 1.3e-3`.
* **A CFT theorem?** **No.** It needs bulk-EFT power counting for the infinitely many higher-point couplings (Sec. 4.1). The statement is "every holographic CFT3 whose bulk EFT has a single cutoff at the gap", not "every CFT3 with a gap". It also holds only at c_T = infinity, and for d=3 eta itself requires the c_T -> infinity limit before omega -> 0.
* **PRL headline?** **Not as the headline; yes as a one-paragraph corollary.** The novelty is the tau bound, and the eta/s translation multiplies it by 3/40. A referee will note that "Einstein-dominated holographic plasmas have eta/s ~ 1/4pi up to 1/Delta^4" is folklore, and that the extra assumption is precisely the kind of power counting the paper claims to replace with theorems. For cond-mat readers the number is striking but physically remote: no laboratory system has `c_T -> infinity` with a gap of 1000. The physically relevant deviations from 1/4pi in real systems are 1/c_T and small-gap effects, which this does not control. As an *illustration* of what "one part in 10^8 of Einstein" means, it is excellent, and it costs one paragraph.
* **Better observable for cond-mat?** The natural one is the **DC conductivity** of a conserved U(1) current, `sigma(0)/sigma(inf) = 1 + 4 gamma` (Myers-Sachdev-Singh; collider `|gamma| <= 1/12`). It is controlled by `<JJT>`, whose deviation from Einstein-Maxwell enters at `1/Delta^2` (four derivatives, `C F F`), much larger than `t4 ~ 1/Delta^4`. An explicit-constant bound on gamma from `<JJJJ>`/`<JJTT>` dispersive sum rules would be a real cond-mat headline, but it is new work, not a corollary. The **Hall viscosity** eta_H of a parity-violating holographic CFT3 would be linear in alpha_T. That makes it the only thermal transport coefficient where the parity-odd half of the bound matters. **Update:** it does. See `hall/README.md`: `4 pi eta_H/s = -18 alpha3~ - 216 alpha4~` (orientation eps_{vrxy} = +r^2), so `4pi|eta_H|/s = (3/40)|alpha_T|` at cubic order. Saremi-Son 1103.4851 need a dynamical axion; Hoyos et al. 1906.04752 need a magnetic field. It is worth one afternoon.

### Draft corollary paragraph (if used)

*Shear viscosity*. The bounds fix the hydrodynamics of the dual plasma. At leading order in `1/c_T`, the Weyl^3 coupling shifts the shear viscosity of the AdS4 black brane by `4pi eta/s - 1 = -(3/40) t4`, the quartic coupling by `+216 alpha4`, and the parity-odd couplings not at all [Li-Lü-Mai; SM]. Hence, if the couplings with five or more gravitons are no larger than the gap implies, every holographic CFT3 with `Delta_g >= 1000` has `|4pi eta/s - 1| <= 2.2 x 10^-9`; without positivity input, `|4pi eta/s - 1| <= (3/40) tau_*(Delta_g) = (4/3)(Delta_g-2)^{-1}[1-49/(2Delta_g+9)^2]^{-1/2}` at every gap (`1.3 x 10^-3` at 1000). Conversely, a plasma whose `eta/s` differs from `1/4pi` by a fraction `delta` must contain a higher-spin state below twist about `2 + 1.3/delta`.
