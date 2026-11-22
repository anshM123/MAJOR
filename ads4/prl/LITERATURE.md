# Literature positioning for the PRL draft (answers referee M13, M14 and minor items 3, 9, 10)

*Written 22 Nov 2026. Read-only pass: no existing file was edited and nothing was committed. The verified bibliography is
`prl/refs_verified.bib`; bib keys below refer to it. Unreviewed.*

**How it was checked.**
* Every paper cited here was looked up on the arXiv API (title, authors, abstract) and on the INSPIRE API (journal, DOI).
* For the eleven papers the positioning hinges on, the full PDF text was read at the passages quoted: 1904.05905, 1904.05892,
  1610.09378, 1805.07393, 2104.11238, 2311.04271 (v1, v2, v3), 2102.08951, 2106.10274, 2103.12728, 1906.07195, 1804.00949,
  2011.05728, 1205.4142, 1711.07343.
* "Abstract only" marks claims taken from an abstract.
* INSPIRE lists 15 papers citing CLS-D (as of 22 Nov 2026). Their titles were checked, and none derives AdS R³ or R⁴ bounds.
* Web searches used paper titles and topics only.

---

## 1. Positioning (citable text, ≈ 700 words)

*What is known.*

**Cubic couplings: parametric causality results.** Causality constraints on the cubic graviton couplings start with CEMZ
[CEMZ].
* Non-Einstein three-graviton structures lead to causality violation in high-energy scattering unless an infinite tower of
  higher-spin states is added.
* In AdS₅/CFT₄ this means |a−c|/c ≲ Δ_gap⁻².

In the CFT, the same suppression follows from Regge-limit causality and unitarity of *mixed* correlators of T with probe
operators:
* [AHKT], d = 4: t₂/c ≲ Δ_gap⁻² and t₄/c ≲ Δ_gap⁻⁴, assuming spin-2 double traces are projected out.
* [CHP]: non-minimal C_TTT ∝ Δ_g⁻² or Δ_g⁻⁴. [MeltzerPerlmutter]: ⟨TTO⟩ suppressed by powers of Δ_gap.
* d = 3, [AKT, AHKTshock]: the holographic null energy condition gives t₄ and α_T ≲ ln Δ_gap/Δ_gap⁴ (their eq. 7.7, written for
  ñ_s−4ñ_f and ñ_odd), with unspecified O(1) constants.

Without a gap:
* The conformal collider bounds [HofmanMaldacena, HLMPR] confine (t₄, α_T) to a disc in d = 3 [CDPrakash]. The parity-odd
  structure exists only in d = 3 [GPY].
* The numerical ⟨TTTT⟩ bootstrap reproduces these bounds and tightens them at moderate gaps [DKKPS].

**Sum rules on ⟨TTTT⟩ itself.** KKSDZ [KKSDZ] proved that ANEC operators commute:
* J₀ ≤ 1 non-perturbatively (§4.2.4, eq. 4.55, generalising [CaronHuotInversion]), with a caveat discussed in §4.2.6.
* J₀ ≤ 2 in the planar limit.

From ⟨T|[ℰ,ℰ]|T⟩ they derived superconvergence sum rules for general d (eq. 5.134).
* These express ((4+d)t₂+4t₄)² through heavy single-trace data.
* Hence t₂ and t₄ vanish in the planar limit if T is the only single trace.
* In §6.1 they argue, non-rigorously and assuming no Regge "outliers", that heavy contributions are suppressed by powers of
  1/Δ_gap.

Belin–Hofman–Mathys [BHM] reach t₂ = t₄ = 0 (a = c) in d = 4 from the same commutator, in the infinite-gap limit. They state
that finite-gap corrections "including numerical factors" are left open. The light-ray OPE [KKSDZlightray] underlies both papers.

So ⟨TTTT⟩ light-ray sum rules that constrain t₄ are *not* new. What has been missing is an explicit-constant, finite-gap
inequality.

**Dispersive sum rules and AdS EFT bounds.**
* CHMRS unified dispersive CFT sum rules [CHMRSdisp].
* In flat space, CHMRS bounded gravitational couplings at small impact parameter [CHMRSswamp], with Caron-Huot–Van Duong
  [CaronHuotVanDuong] for the no-gravity case. This includes a D = 4 infrared-log bound and 0 ≤ g₀M⁶/8πG ≤ 3.000 on the R⁴
  coefficient of ten-dimensional maximal supergravity. Type II string theory sits at 2ζ(3).
* CHMRS uplifted flat bounds to AdS for *scalar* correlators with graviton exchange [CHMRS]. In AdS₄ they obtain a lower bound
  g₂ ≥ −17.64 (8πG/M²) log(0.124 M R_AdS) and two-sided bounds on g₃ at fixed g₂ (eqs. 5.64–5.65).
* Four-dimensional flat graviton bounds with an IR logarithm are in [CHLPS]; D ≥ 5 is in [CHLPSpw].
* Bern–Kosmopoulos–Zhiboedov [BKZ] bound the four-dimensional R³ coefficient by the R⁴ one, |β_R3|² ≤ β⁺_R4/m²_gap (eq. 6.13),
  at leading order in massive exchange.
* Kundu [Kundu] bounds *scalar* contact couplings φ²□ᵏφ² in AdS_D (D ≥ 4), using the Regge limit of a scalar CFT four-point
  function.
  * For even k ≥ 2 the couplings obey positivity, monotonicity and log-convexity.
  * Gravity weakens only the k = 2 condition.
  * He proposes extending the method to stress-tensor Regge correlators, but does not carry it out.
* CLS-D [CLSD] constructed spinning dispersive sum rules for ⟨TTTT⟩ and their flat-space graviton dictionary at the "scattering
  crystal". They list three remaining tasks: the non-shock sum rules, positivity in the Regge regime, and the light-side
  dictionary.

**Supersymmetric determinations.** In ABJ(M) the R⁴ term is *determined* within the CFT, not bounded: by localisation plus Ward
identities in the M-theory regime [CPY] and, as a function of N/k⁵, in the IIA regime, including the ζ(3) term and its genus-one
correction [BCPW]. The 𝒩 = 6 bootstrap with localisation inputs gives rigorous islands across both regimes [BCJP]; see also
[BCPD4R4, CDempseyP]. For the AdS₅ analogue see [CCTZ, CCZ].

*How this work differs.* We take CLS-D's ⟨TTTT⟩ functionals and carry out, *conditionally*, the steps they list, in
AdS₄/CFT₃, for t₄, α_T and one R⁴ coefficient. Compared with what is above:
* **Explicit constants at fixed Δ_g** (e.g. √(t₄²+α_T²) ≤ 2.89·10⁻⁸ at Δ_g = 1000) and closed-form Λ-dependence. [KKSDZ, BHM]
  give vanishing at infinite gap; [AKT] gives scaling only.
* **The d = 3 parity-odd coupling α_T** is bounded together with t₄.
* **A window for a graviton R⁴ coefficient in AdS.** Earlier AdS bounds concern scalar contacts [CHMRS, Kundu], and earlier
  graviton R⁴ bounds are flat [CHMRSswamp, CHLPS, BKZ].

The result is not a theorem:
* Per-state positivity is established only numerically, on finite sets of states.
* It also rests on the analytic inputs of Table I.
* It applies only to theories with a gap in every single trace other than T, which excludes ABJM.
* Item (2) gives no flat-type R³/R⁴ ratio of BKZ form.

---

## 2. Verdict on "no published bound on R³ … from ⟨TTTT⟩ sum rules" (`letter.tex:366-368`)

**False as written. It is true only after qualification.**

1. **⟨TTTT⟩ light-ray sum rules constraining t₄ already exist.**
   * KKSDZ 1904.05905 §5.3, eq. (5.134), derive sum rules from ⟨T(n₄)|ℰ(n₁)ℰ(n₂)|T(n₃)⟩ for general d. These are superconvergence
     sum rules of the same family as CLS-D's.
   * The intro (p. 4–5) concludes "t₂ and t₄ must vanish in the planar limit" when no other single trace is present.
   * §6.1 gives a non-rigorous large-gap estimate, |α_GB(s)|² ≲ 1/(Δ_gap⁴ c_T) for the non-minimal coupling (eq. 6.7), assuming
     no outliers.
   * That is a parametric bound from a ⟨TTTT⟩-type sum rule.
   * Belin–Hofman–Mathys 1904.05892 (eq. 4.24) obtain t₂ = t₄ = 0 from ⟨T|[ℰ₁,ℰ₂]|T⟩ = 0, in d = 4 at infinite gap.
2. **Parametric R³ bounds in d = 3 exist,** but from other correlators.
   * AKT 1805.07393 eq. (7.7): t₄ and α_T ≲ ln Δ_gap/Δ_gap⁴, from the holographic null energy condition, which follows from Regge
     causality of ⟨ψψTT⟩-type correlators [AHKT, AHKTshock].
   * [AHKT, CHP, MeltzerPerlmutter] give the analogues in d = 4 or in general d.
3. **Numerical ⟨TTTT⟩ crossing bounds the d = 3 ⟨TTT⟩ parameter.** DKKPS 1708.05718 do this, with tightened restrictions at
   moderate gaps (abstract). These are not large-gap holographic bounds.
4. **R⁴.**
   * No published *bound* on a graviton R⁴ coefficient in AdS from ⟨TTTT⟩ was found. The AdS bounds are for scalar contacts
     [CHMRS, Kundu]; the graviton R⁴ bounds are in flat space [CHMRSswamp, CHLPS, BKZ].
   * In ABJM the R⁴ term has been *determined* from the stress-tensor-multiplet four-point function, which contains ⟨TTTT⟩ by
     supersymmetry [CPY, BCPW].

**Defensible replacement claim:** "To our knowledge, no explicit-constant, finite-gap bound on t₄, α_T or a graviton R⁴ coefficient
has been derived from ⟨TTTT⟩ sum rules." It must be stated together with items 1–4 above. It also inherits the conditional status
of our result.

---

## 3. ABJM (referee M13)

**What localisation gives.**
* **BCP 1906.07195** (key `BCPW`; three authors), IIA regime.
  * They study the stress-tensor-multiplet four-point function ⟨SSSS⟩ at N → ∞ with fixed μ = N/k⁵.
  * The R⁴ contact Mellin amplitude is fixed by ∂⁴_m F of the mass-deformed S³ partition function (Fermi gas, all orders in 1/N),
    together with Ward identities and bulk locality.
  * In the 't Hooft limit (eq. 1.2) it appears at order ζ(3) λ^{−3/2}/c_T relative to supergravity, with
    λ = N/k − 1/(3k²) − 1/24 (eq. 3.40).
  * The flat-space limit, using L⁸/ℓ_s⁸ = 4π⁴λ² (eq. 4.11), reproduces A_IIA = A_SG(1 + ℓ_s⁶ ζ(3)/32 · stu + g_s² ℓ_s⁶ π²/96 ·
    stu + …) (eqs. 4.14–4.15). At finite μ the genus-one term is also reproduced (eqs. 4.20–4.21).
* **CPY 1804.00949**, M-theory regime, k = 1, 2 (𝒩 = 8).
  * Localisation fixes the c_T^{−5/3} contact term.
  * It matches the 11d R⁴ term ℓ₁₁⁶ stu/(3·2⁷).
* **BCJP 2011.05728** computes no R⁴ number.
  * It gives the 𝒩 = 6 superconformal block decomposition.
  * With localisation inputs (c_T and one protected OPE coefficient), it gives rigorous numerical islands interpolating between
    the M-theory and IIA regimes.
  * The referee's "BCJP in the IIA regime" should therefore read "BCP 1906.07195 in the IIA regime; BCJP for the 𝒩 = 6
    bootstrap".

**Same quantity, same normalisation?**
* **Normalisation: same.**
  * Our g₄/(8πG) = ζ(3)α′³/32, with M² = 4/α′, gives g₄M⁶/8πG = 2ζ(3).
  * This is BCP's stu coefficient ℓ_s⁶ζ(3)/32 (ℓ_s² = α′), and it is exactly the normalisation of CHMRS [CHMRSswamp, eqs.
    3.35, 3.43].
  * The dictionary L²/α′ = π√(2λ) is BCP's eq. (4.11).
  * The c_T convention is the same (free real scalar c_T = 1, BCP eq. 1.3), so 8πG = 256/c_T applies to both.
  * The λ-scaling agrees: (α′/L²)³ = (π√(2λ))⁻³ ∝ λ^{−3/2}.
* **Quantity: the same flat coefficient, not the same CFT object.**
  * BCP fix the coefficient of a degree-4 Mellin amplitude of ⟨SSSS⟩ and compare it with the *ten-dimensional* IIA amplitude.
  * Our g̃₄ is defined by a 4d light-ray functional on ⟨TTTT⟩. Identifying the two requires all of the following:
    * (F)+(P);
    * restricting the 10d amplitude to 4d graviton kinematics;
    * the supersymmetric relation between ⟨SSSS⟩ and ⟨TTTT⟩.
* **Hence 2ζ(3) is not new information.** It is the value the CFT already yields through localisation.
* **What the placement adds is only the ratio to our bound.** That ratio needs Λ ↔ ML, which BCP do not use and which holds only
  at leading order.
* **It cannot test the bound,** because ABJM violates A1 (light KK modes, and light R-singlets in T×T).
* **The 3.3 % / 2.6 % mainly reflects the 4d IR logarithm** in our upper end (17.24 ln(2Λ/3)). The natural flat comparison for a
  10d maximally supersymmetric limit is CHMRS's 0 ≤ g₀M⁶/8πG ≤ 3.000, where 2ζ(3) sits at 80 %.
* **M-theory regime (CPY): not comparable.** There is no string scale. The R⁴ term is suppressed as c_T^{−2/3} relative to
  supergravity, so Λ ↔ ML is undefined.

**c_T versus the flavour coefficient (checked, no confusion).**
* c_T = −(64/π²) ∂²_{m±} log Z (BCP eq. 3.12; ACP eq. 26) holds for masses coupling to currents in the stress-tensor multiplet.
* CDFKS 1205.4142 eq. (4.6) is a different statement: (1/r²)∂²_m Re F = (π²/2)τ_ff, the *flavour*-current two-point coefficient.
* The letter's parenthesis is therefore correct, but it is irrelevant to the placement, and the referee asks for it to be deleted.

**Recommendation.**
* Remove ABJM from the abstract.
* Reduce the paragraph to a sign-and-normalisation check citing [BCPW, CPY] for the CFT determination and [CHMRSswamp] for the
  normalisation (text in §4, item 10).

---

## 4. Sentences in `letter.tex` that must change

| # | file:line | current text (quoted) | proposed replacement |
|---|---|---|---|
| 1 | `letter.tex:35` | "The ABJM value $2\zeta(3)$ lies at $\ABJMHead\%$ of the upper end." | *Delete* (M13). |
| 2 | `letter.tex:45-48` | "For the stress tensor itself this expectation has been made parametric: causality forces the non-Einstein $\langle TTT\rangle$ structures to vanish as the gap grows~\cite{CEMZ,AHKT,MeltzerPerlmutter,CHP}, in three dimensions as $t_4\lesssim \ln\Delta_g/\Delta_g^4$ with an unknown constant~\cite{AKT}." | "For the stress tensor this expectation has been made parametric. Causality of high-energy scattering~\cite{CEMZ} and Regge-limit causality of correlators of $T$ with probe operators~\cite{AHKT,CHP,MeltzerPerlmutter} force the non-Einstein $\langle TTT\rangle$ structures to vanish as powers of $1/\Delta_g$; in three dimensions $|t_4|,|\alpha_T|\lesssim\ln\Delta_g/\Delta_g^4$ with unspecified constants~\cite{AKT,AHKTshock}. On $\langle TTTT\rangle$ itself, commutativity of ANEC operators gives superconvergence sum rules that express these couplings through heavy single-trace data and force them to vanish when $T$ is the only light single trace~\cite{KKSDZ,BHM}, with a heavy contribution argued to be suppressed by powers of $1/\Delta_g$~\cite{KKSDZ}." |
| 3 | `letter.tex:49-51` | "For scalar correlators, Caron-Huot, Mazáč, Rastelli and Simmons-Duffin (CHMRS) turned flat-space dispersive bounds into two-sided $\log\Delta_g$ bounds on AdS$_4$ contact couplings with explicit constants~\cite{CHMRS}." | "Dispersive CFT sum rules~\cite{CHMRSdisp} turn flat-space bounds~\cite{CaronHuotVanDuong,CHMRSswamp} into AdS bounds: for scalar correlators with graviton exchange, Caron-Huot, Mazáč, Rastelli and Simmons-Duffin (CHMRS) obtained in AdS$_4$ a lower bound on the leading contact coupling and two-sided bounds on the next, with explicit constants and a $\log(MR_{\rm AdS})$ from the AdS infrared regulator~\cite{CHMRS}. Kundu derived positivity, monotonicity and log-convexity of scalar contact couplings $\phi^2\Box^k\phi^2$ in AdS from the Regge limit of scalar correlators~\cite{Kundu}." |
| 4 | `letter.tex:51-54` | "… and left ``the remaining steps for uplifting flat space bounds to AdS'' as future work~\cite{CLSD}." | The quoted string **does not occur in CLS-D** (checked in v1, v2, v3). Replace with: "… and listed what remains (sum rules beyond shock amplitudes, positivity in the Regge regime, and the light-side dictionary) before flat-space graviton bounds can be lifted to AdS~\cite{CLSD}." If a quotation is wanted, the verbatim text (v3, §6) is "there are still several remaining tasks". |
| 5 | `letter.tex:54-55` | "Flat four-dimensional graviton bounds with an infrared logarithm are given in Ref.~\cite{CHLPS}." | "Flat four-dimensional graviton bounds with an infrared logarithm are given in Ref.~\cite{CHLPS} (sharp bounds for $D\ge5$ in Ref.~\cite{CHLPSpw}); in four dimensions the $R^3$ coefficient is bounded by the $R^4$ coefficient~\cite{BKZ}, and in ten-dimensional maximal supergravity $0\le g_4M^6/8\pi G\le3.000$~\cite{CHMRSswamp}." |
| 6 | `letter.tex:57-58` | "This Letter reports an explicit-constant, \emph{conditional} execution of that uplift for the cubic ($R^3$) and quartic ($R^4$) graviton couplings in AdS$_4$/CFT$_3$." | "This Letter reports an explicit-constant, \emph{conditional} execution of these steps in AdS$_4$/CFT$_3$ for the cubic couplings ($t_4$ and the parity-odd $\alpha_T$) and one quartic ($R^4$) graviton coupling. Relative to Refs.~\cite{KKSDZ,BHM,AKT}, which give vanishing at infinite gap or parametric scaling, the new elements are explicit constants at fixed $\Delta_g$, the inclusion of $\alpha_T$, and a bound on a graviton $R^4$ coefficient in AdS rather than in flat space~\cite{CHMRSswamp,CHLPS,BKZ} or for scalar contacts~\cite{CHMRS,Kundu}." |
| 7 | `letter.tex:93` | "We write $\langle TTT\rangle$ in the basis $(c_T,t_4,\alpha_T)$~\cite{HofmanMaldacena,Buchel};" | "We write $\langle TTT\rangle$ in the basis $(c_T,t_4,\alpha_T)$~\cite{HofmanMaldacena,Buchel,CDPrakash}, where $\alpha_T$ multiplies the parity-odd structure that exists only in $d=3$~\cite{GPY};" (minor item 3). |
| 8 | `letter.tex:95-97` | "The quartic coupling $g_4$ is the flat $R^4$ coefficient normalised so that the type~II Virasoro--Shapiro amplitude has $g_4/(8\pi G)=2\zeta(3)/M^6$, $M^2=4/\alpha'$." | Add the source: "… normalised as in Ref.~\cite{CHMRSswamp}, so that …". |
| 9 | `letter.tex:105-108` | "… which we take from the Regge bound of Kologlu, Kravchuk, Simmons-Duffin and Zhiboedov (KKSDZ)~\cite{KKSDZ} with its stated loophole (input K1)." | "… which we take from the non-perturbative Regge bound $J_0\le1$ of Kologlu, Kravchuk, Simmons-Duffin and Zhiboedov (KKSDZ)~\cite[Sec.~4.2.4]{KKSDZ}, which generalises the scalar argument of Ref.~\cite{CaronHuotInversion}, with the caveat of their Sec.~4.2.6 (input K1)." (minor item 9). |
| 10 | `letter.tex:338-352` (whole ABJM paragraph; includes `:347-350`) | "In ABJM theory~\cite{ABJM} at strong coupling, … It is a consistency check, not evidence for the theorem." | Footnote-level: "\emph{ABJM.} In type IIA string theory, $g_4M^6/8\pi G=2\zeta(3)$~\cite{CHMRSswamp}. For ABJM theory~\cite{ABJM} this coefficient has been derived within the CFT from localisation and the stress-tensor-multiplet four-point function~\cite{BCPW} (M-theory regime:~\cite{CPY}). With $\Lambda\leftrightarrow ML$ at leading order in $1/(ML)$, it lies inside the window of Eq.~\eqref{eq:quartic}, and $\hat g_3=0$ by ten-dimensional supersymmetry. ABJM violates A1, and its flat limit is ten-dimensional, so this checks only sign and normalisation." Delete the $c_T$/flavour parenthesis (`:347-350`). |
| 11 | `letter.tex:361-363` | "The parametric form of Eq.~\eqref{eq:cubicL} differs from the $\ln\Delta_g/\Delta_g^4$ of Ref.~\cite{AKT}, and the constants are not comparable." | "Eq.~\eqref{eq:cubicL} scales as $\sqrt{\ln\Lambda}/\Lambda^4$; Ref.~\cite{AKT} gives the parametric $\ln\Delta_g/\Delta_g^4$ for $t_4$ and $\alpha_T$ separately, without constants, so no numerical comparison is possible." |
| 12 | `letter.tex:366-370` | "Literature searches (abstract level, plus one extended search) found no published bound on $R^3$ or $R^4$ couplings from $\langle TTTT\rangle$ sum rules, so novelty is asserted only at search level; the relation to gap-scaling arguments for higher-derivative couplings~\cite{Kundu} and to the $\mathcal N=6$ numerical bootstrap~\cite{BCJP} remains to be discussed in detail. % TODO: expand if space allows." | "To our knowledge, no explicit-constant, finite-gap bound on $t_4$, $\alpha_T$ or a graviton $R^4$ coefficient has previously been derived from $\langle TTTT\rangle$ sum rules: the superconvergence sum rules of Refs.~\cite{KKSDZ,BHM} give vanishing in the planar, infinite-gap limit together with a parametric suppression, and the three-dimensional bounds of Ref.~\cite{AKT} are parametric and come from other correlators. Kundu's constraints~\cite{Kundu} concern scalar contact couplings. In ABJ(M) theory the $R^4$ term is determined rather than bounded, by localisation~\cite{CPY,BCPW}, and the $\mathcal N=6$ bootstrap gives rigorous islands for the same theories~\cite{BCJP}; our bounds need no supersymmetry but are conditional." Remove the TODO. |
| 13 | `letter.tex:316` | "… an explicit tree-level $\langle TTTT\rangle$ exists only in momentum space~\cite{Raju,AlbayrakKharel}, …" | Soften the claim of absence, which was not verified exhaustively: "… explicit tree-level results for $\langle TTTT\rangle$ are, to our knowledge, available in momentum space~\cite{Raju,AlbayrakKharel}, …". |
| 14 | `letter.tex:388` | "\bibliography{refs}" | "\bibliography{refs_verified}" (also `supplement.tex:400`), or copy `refs_verified.bib` over `refs.bib` after review. |

**Supplement, same theme:**
* `supplement.tex:396`: "Novelty is asserted only at literature-search level." → "Positioning relative to Refs.~\cite{KKSDZ,BHM,AKT,CHMRS,CLSD} is given in the Letter."
* `supplement.tex:383-385`: the c_T/flavour sentence is correct as written. Keep it only if the letter keeps a pointer to it.

**Outside this task, noted only:** the title (`:14-15`) and M12 ("one quartic coupling"). Row 6 already uses "one quartic
coupling".

---

## 5. Bibliography changes (`refs.bib` → `refs_verified.bib`)

**Totals.**
* **41 entries, all verified on arXiv and INSPIRE; 0 TODO remaining.** That is 22 original entries plus 19 new ones.
* Every original key is unchanged, so the existing `\cite` commands resolve.
* A BibTeX test (`\nocite{*}`, `apsrev4-2`) gives only "missing number" warnings (issue numbers), as before.

**Applied to all 22 original entries:**
* surname-only authors replaced by full names, as printed on arXiv/INSPIRE;
* `@misc` changed to `@article`;
* journal, volume, pages and DOI added;
* `year` set to the publication year.

| key | arXiv | changes beyond the common ones |
|---|---|---|
| HPPS | 0907.0151 | title "Holography from Conformal Field Theory"; JHEP 10 (2009) 079 |
| HofmanMaldacena | 0803.1467 | title "Conformal collider physics: Energy and charge correlations"; JHEP 05 (2008) 012 |
| Buchel | 0911.4257 | title "Holographic GB gravity in arbitrary dimensions"; JHEP 03 (2010) 111 |
| OsbornPetkou | hep-th/9307010 | was entirely TODO. Now: arXiv number; title "Implications of conformal invariance in field theories for general dimensions"; Ann. Phys. 231 (1994) 311–362. arXiv metadata misspells "Petkos"; INSPIRE has "Petkou, A.C." Initials only (given names not verified) |
| CEMZ | 1407.5597 | title "Causality Constraints on Corrections to the Graviton Three-Point Coupling"; JHEP 02 (2016) 020; year 2014 → 2016 |
| AHKT | 1610.09378 | title "Einstein gravity 3-point functions from conformal field theory"; JHEP 12 (2017) 049; year 2016 → 2017. **It is d = 4** (as the repo noted) |
| CHP | 1707.07689 | title "Bounds for OPE coefficients on the Regge trajectory"; JHEP 10 (2017) 197 |
| MeltzerPerlmutter | 1712.04861 | title "Beyond a=c: Gravitational Couplings to Matter and the Stress Tensor OPE"; JHEP 07 (2018) 157; year → 2018 |
| AKT | 1805.07393 | title "A Conformal Collider for Holographic CFTs"; JHEP 10 (2018) 156 |
| KKSDZ | 1904.05905 | title confirmed; JHEP 11 (2020) 096; year 2019 → 2020 |
| CHMRS | 2106.10274 | title confirmed; JHEP 11 (2021) 164 |
| CHLPS | 2201.06602 | title "Causality constraints on corrections to Einstein gravity" added; JHEP 05 (2023) 122 confirmed |
| CLSD | 2311.04271 | title "Spinning dispersive CFT sum rules and bulk scattering" added; JHEP 04 (2025) 016 confirmed |
| CaronHuotLi | 2102.08160 | **title corrected** to "Helicity basis for three-dimensional conformal field theory" (repo had "…CFT"); JHEP 06 (2021) 041 |
| Raju | 1201.6452 | title confirmed; Phys. Rev. D 85, 126008 (2012) |
| AlbayrakKharel | 1908.01835 | title and JHEP 12 (2019) 135 confirmed. The journal prints the title as "…amplitudes. Part II. Gravitons"; the arXiv form is kept |
| ABJM | 0806.1218 | title "N=6 superconformal Chern-Simons-matter theories, M2-branes and their gravity duals"; JHEP 10 (2008) 091 |
| Kundu | 2104.11238 | title "Swampland Conditions for Higher Derivative Couplings from CFT"; JHEP 01 (2022) 176; year → 2022 |
| BCJP | 2011.05728 | title "The 3d 𝒩=6 Bootstrap: From Higher Spins to Strings to Membranes"; JHEP 05 (2021) 083; year → 2021 |
| ACP | 1711.07343 | title "Solving M-theory with the Conformal Bootstrap"; JHEP 06 (2018) 159; year → 2018 |
| BCPW | 1906.07195 | **authors corrected to Binder–Chester–Pufu (three)**; ABJM_THEOREM.md was right. "Binder–Chester–Pufu–Wang" is arXiv:1902.06263 (𝒩 = 4 SYM). Title "AdS₄/CFT₃ from Weak to Strong String Coupling"; JHEP 01 (2020) 034. Key kept; rename to `BCP` at the next letter edit |
| Fmax | 1205.4142 | **authors added**: Closset–Dumitrescu–Festuccia–Komargodski–Seiberg. Title "Contact Terms, Unitarity, and F-Maximization in Three-Dimensional Superconformal Theories"; JHEP 10 (2012) 053. Eq. (4.6) checked: flavour coefficient τ_ff, not c_T |

**New entries (19):**

| key | arXiv | reference | why |
|---|---|---|---|
| BHM | 1904.05892 | Belin, Hofman, Mathys, JHEP 08 (2019) 032 | M14; ⟨T\|[ℰ,ℰ]\|T⟩ ⇒ t₂ = t₄ = 0 (d = 4) |
| CHMRSdisp | 2008.04931 | Caron-Huot, Mazáč, Rastelli, Simmons-Duffin, "Dispersive CFT Sum Rules", JHEP 05 (2021) 243 | M14 |
| CHMRSswamp | 2102.08951 | same authors, "Sharp boundaries for the swampland", JHEP 07 (2021) 110 | M14; source of the g₄ normalisation |
| BKZ | 2103.12728 | Bern, Kosmopoulos, Zhiboedov, J. Phys. A 54 (2021) 344002 | M14 |
| CPY | 1804.00949 | Chester, Pufu, Yin, JHEP 08 (2018) 115 | M13 |
| BCPD4R4 | 1808.10554 | Binder, Chester, Pufu, "Absence of D⁴R⁴ in M-Theory From ABJM", JHEP 04 (2020) 052 | ABJM context |
| CDempseyP | 2412.14094 | Chester, Dempsey, Pufu, JHEP 07 (2025) 096 | ABJM context |
| CaronHuotInversion | 1703.00278 | Caron-Huot, "Analyticity in Spin in Conformal Theories", JHEP 09 (2017) 078 | minor 9 |
| KKSDZlightray | 1905.01311 | Kologlu, Kravchuk, Simmons-Duffin, Zhiboedov, "The light-ray OPE and conformal colliders", JHEP 01 (2021) 128 | task list |
| CDPrakash | 1707.03007 | Chowdhury, David, Prakash, JHEP 11 (2017) 171 | minor 3 (α_T notation) |
| GPY | 1104.4317 | Giombi, Prakash, Yin, JHEP 07 (2013) 105 | minor 3 (parity-odd structure) |
| CGGHJM | 1910.14392 | Chowdhury, Gadde, Gopalka, Halder, Janagal, Minwalla, JHEP 02 (2020) 114 | M12 (classification) |
| CaronHuotVanDuong | 2011.02957 | Caron-Huot, Van Duong, "Extremal Effective Field Theories", JHEP 05 (2021) 280 | flat dispersive bounds |
| CHLPSpw | 2205.01495 | Caron-Huot, Li, Parra-Martinez, Simmons-Duffin, Phys. Rev. D 108 (2023) 026007 | D ≥ 5 graviton bounds |
| DKKPS | 1708.05718 | Dymarsky, Kos, Kravchuk, Poland, Simmons-Duffin, "The 3d Stress-Tensor Bootstrap", JHEP 02 (2018) 164 | ⟨TTTT⟩ bounds on ⟨TTT⟩ in d = 3 |
| HLMPR | 1603.03771 | Hofman, Li, Meltzer, Poland, Rejon-Barrera, JHEP 06 (2016) 111 | collider bounds proof |
| AHKTshock | 1709.03597 | Afkhami-Jeddi, Hartman, Kundu, Tajdini, "Shockwaves from the OPE", JHEP 03 (2019) 201 | basis of AKT's d = 3 bound |
| CCTZ | 2207.01615 | Caron-Huot, Coronado, Trinh, Zahraee, JHEP 02 (2023) 083 | AdS₅ analogue |
| CCZ | 2412.00249 | Caron-Huot, Coronado, Zahraee, JHEP 05 (2025) 220 | AdS₅ analogue |

Kundu (2104.11238), CEMZ, Hofman–Maldacena, AHKT, KKSDZ, CLS-D and Caron-Huot–Li were already in `refs.bib` and are verified
above.

---

## 6. Open items and caveats

* **CCTZ / CCZ** are described from their abstracts only.
* **DKKPS's "tighter restrictions at moderate gaps"** is also from the abstract.
* **The absence of a position-space Einstein ⟨TTTT⟩ in AdS₄ (row 13) was not settled.** One web search found only momentum-space
  results.
* **BCP eq. (4.21)** (the finite-μ R⁴ coefficient) did not extract cleanly from the PDF. It is cited by number only.
* **Search scope for §2.** The §2 verdict rests on:
  * the INSPIRE citation lists of CLS-D (15 citing papers);
  * the papers read above;
  * three web searches.

  A full INSPIRE `refersto` sweep of KKSDZ and AKT (hundreds of citing papers) was not done. A final check of those lists for
  "t₄", "R³" and "superconvergence" in d = 3 before submission is recommended.
* **The referee also cites "Chowdhury et al. 1910.14392"** for the classification of four-graviton contact structures. The
  abstract confirms the topic. The specific statement about D = 4 eight-derivative structures was not checked in the text.
