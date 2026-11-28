# Scoop check: covariant (HRT) entropy inequalities for vacuum multiply connected AdS3

Checked 2026-11-28. Read-only: no repository file was changed apart from this one.

How each claim was checked:
- **[read full text]** means I extracted the PDF text and read or grepped the passages that matter.
- **[abstract]** means I read only the abstract or metadata.
- **[memory]** means the claim comes from my own background knowledge and was not re-checked in this session.

## Bottom line

1. **No scoop found.** Nothing on arXiv up to about 27 Nov 2026 proves or disproves HRT entropy inequalities for
   multiply connected 2+1 bulks. Nothing uses Liouville currents, geodesic currents or intersection numbers on
   the quotient surface for holographic entropy inequalities. The latest statement of the open problem is
   Grimaldi–Headrick, 2610.33905 (24 Oct 2026). It treats the question as open and cites no work on it newer
   than 2407.07165.
2. **There is no earlier version of the Bonahon-current argument.** I found none in physics or in math. The
   closest items, none of which is a scoop:
   - **Czech–Dong.** They used a kinematic-space "intersection number" measure on the universal cover.
   - **Zhang–Chen (1610.07134).** They used a kinematic space quotiented down to static BTZ and multi-boundary
     geometries, with no inequalities.
   - **Martínez-Granado–D. Thurston (2605.04031, May 2026).** Their result is math. It is the main thing a
     referee could use to call Theorem 1 "a short corollary"; see section 2.
3. **Wording corrections needed in our draft** (details in section 1):
   - Grimaldi–Headrick describe the π₁=ℤ case as "**a proof**", not an argument.
   - GGHH's numerics already assume the common-class ½(ℓ_L+ℓ_R) formula, so Theorem H's formula is not new;
     only its rigorous justification is.
   - Our result covers **vacuum** spacetimes only. The open question GGHH state (2407.07165) is for general
     multiply connected bulks. Their own π₁=ℤ argument and Czech–Dong's simply connected proof are not
     restricted to vacuum. We solve the vacuum case of their open question, not the whole question.
4. **PRL calibration.** In this field PRL is rare. I checked about 50 papers on the holographic entropy cone
   (HEC) and HRT inequalities:
   - **2 were in PRL**, both with sharp conceptual punchlines (2025, 2026).
   - **1 was a PRD Letter.**
   - **Every direct predecessor of our result was in JHEP, CQG or PRD**, or is still unpublished.
   - A PRL submission is plausible but uphill; JHEP, CMP or SciPost is the default. Details in section 3.

---

## 1. Papers relevant to Task 1 (most recent first)

Search coverage:
- **arXiv API.** Date-restricted (15 Sep – 28 Nov 2026) and unrestricted abstract searches:
  - "holographic entropy inequalit(y|ies)", "entropy cone", "covariant + entropy cone", HRT + inequalities
  - contraction map, multiboundary/multi-boundary + entropy, "multiply connected" + holographic
  - kinematic space (+BTZ/wormhole/quotient), Crofton + holographic
  - "geodesic current(s)", "Liouville current(s)", "intersection number" + holographic/entropy
  - BTZ/wormhole + entanglement, "extremal surfaces", bit threads, multi-entropy, math.GT/DG + holographic/AdS
- **arXiv listings.** The past-week pages of hep-th, gr-qc, quant-ph and math.GT, with titles grepped.
- **INSPIRE.** Citations of 2407.07165 (19), 1905.03787 (43), 2610.33905 (0), 2602.04888 (1), 2212.11106 (2)
  and 2502.09894 (14). 2206.06946 is not in INSPIRE.
- **Semantic Scholar** citations of the same papers.
- **INSPIRE author feeds since mid-2026** for Headrick, Czech, Hubeny, Bao, Hernández-Cuenca, Grado-White,
  Engelhardt and Rangamani.
- **Two web searches.**

| arXiv id | Date | Authors | Exact claim | Overlap with ours |
|---|---|---|---|---|
| 2611.21726 | 16 Nov 2026 | Ju, Czech | Disproves the claim that every HEI is superbalanced (He–Hubeny–Rangamani 2002.04558). Reports about 2×10⁵ extreme rays with non-extractable mutual information. Says this invalidates several conjectures, including that null reductions of all HEIs pass the majorization test. Nothing on HRT or multiply connected bulks. **[read full text]** | **None** on the claim. It **helps** us: the majorization evidence for HRT (2601.09987, 2601.09989) covers only superbalanced HEIs. Our Theorem 2 covers every contraction-map HEI, including non-superbalanced ones (with LHS−RHS → +∞ at the cut-off). Say so explicitly. Check that our Δ ≥ 0 step (Lemma B, endpoint balance) does not quietly assume superbalance. Formatted for PRL ("5 pages plus appendices"); it discloses use of GPT-6 models. |
| 2611.19035 | 14 Nov 2026 | Anegawa, Iizuka | Multiway maximin bound on the time over which multi-entropy can probe the wormhole interior, for non-rotating BTZ **[abstract]** | None |
| 2611.15854 | 11 Nov 2026 | Caminiti, Pasterski | Swing surfaces and the flat limit of AdS3 RT **[abstract]** | None |
| 2610.33905 | 24 Oct 2026 | Grimaldi, Headrick (Strings 2026 talk) | Quotes in the next two rows **[read full text, §4]** | None. It is the citation for "open". |
| | | | §4: "The validity of higher HEIs for the HRT formula remains an open problem." | |
| | | | Evidence it lists: "A proof that all HEIs hold in 2 + 1 bulk dimensions when the bulk is either simply connected or a two-boundary genus-zero wormhole [Czech–Dong; GGHH]"; numerics on "randomly generated 2+1 dimensional vacuum multiboundary wormhole configurations"; MMI; membrane theory; light-cone perturbations. It names graph models as a possible strategy. It does not mention currents or kinematic space on the quotient. | |
| 2610.07176 | 5 Oct 2026 (v2 14 Oct) | Sahu | Generalised entanglement wedges on a **static** slice of **pure** AdS3. Uses the Crofton formula and kinematic space to derive strong subadditivity of GWs from boundary SSA. States that "quotients of AdS3, such as conical defects and the BTZ black hole" are "not directly described by our construction". **[read full text: scope and discussion]** | None. Shares only the integral-geometry tool. |
| 2610.03929 | 2 Oct 2026 | Headrick, Kasireddy, Rolph | Quantum covariant bit threads, equivalent to the QES prescription **[abstract]** | None |
| 2609.37395 / 2609.37394 | 29 Sep 2026 | Anegawa, Iizuka / Iizuka | Multi-entropy in four-boundary AdS3 wormholes against the Fenchel–Nielsen twist; static RT **[abstract]** | None |
| 2609.22504 | 18 Sep 2026 | Hubeny, Liu | Tripartite form of the two infinite HEI families; RT only **[abstract]** | None |
| 2608.25287 | 26 Aug 2026 | Bao, Ferko | HEIs beyond SSA along RG flows, small light-cone deformations **[abstract]** | None |
| 2602.04888 v2 | Jan 2026 (v2 Apr 2026) | Zhao | Conditional no-short-cut theorem for covariant graph models, assuming exposed regions exist **[read full text: grep]**. No "Crofton", "kinematic", "Liouville", "multiply", "wormhole" or "quotient" in the text. | None. Conditional and general-dimension, as our LITERATURE.md already says. |
| 2601.09987 v2 | Jan 2026 (v2 16 Oct 2026) | Grimaldi, Headrick, Hubeny, Shteyner | Combinatorial framework for HEIs. The null reduction of any superbalanced HEI passes the majorization test. v2 adds "our results apply to contraction-provable HEIs" **[abstract + grep of full text]**. No geometry, currents or wormholes. | None |
| 2407.07165 v2 | Jul 2024; JHEP 01 (2025) 065 | Grado-White, Grimaldi, Headrick, Hubeny | The Czech–Dong universal-cover step is wrong ("an infinite number of surfaces enter on both sides"). Gives an "argument" (their footnote 11: "not … rigorous") for π₁=ℤ with *any* NEC bulk, via finite covers. §3.2 explains why that argument fails for free groups. Numerics on vacuum (3,0) and (1,1) wormholes. **[read full text: §1, §2.1, §3]** | **Partial, on Theorem H only.** Their numerics already compute S(A) as a minimum over homotopy classes of ½(ℓ_L+ℓ_R)-type lengths, with the homology condition imposed combinatorially (eqs. 2.8–2.10, 2.27). The *formula* in Theorem H is their working assumption; what is new in ours is the proof that the homology region sits on an acausal slice. Their footnote 10 also admits self-intersecting phases; we should state why embedded cuts suffice for the common-class functional. |
| 1905.03787 | 2019; JHEP 10 (2019) 177 | Czech, Dong | All contraction-map HEIs for HRT in 2+1 dimensions, using a positive measure from SSA (conditional mutual information) on a discretised kinematic space. Their "intersection numbers" are "topological quantities". §6.1 lifts to the universal cover with a π₁-invariant auxiliary geometry; that step is wrong. Not restricted to vacuum. **[read full text: §1, §6, §7]** | **Partial, on the mechanism**, which our draft already credits. Our Lemma A is the π₁-equivariant crossing count run downstairs on the quotient. Their measure is the SSA density; ours is the Liouville current. That is why ours is limited to vacuum and theirs is not. |
| 1610.07134 | 2016; JHEP 2017 | Zhang, Chen | Kinematic space quotiented for static BTZ and multi-boundary wormholes. Horizon length computed as a kinematic-space geodesic length, with fundamental regions discussed. **[abstract]** | Small. A static precursor of "Crofton on the quotient" (Lemma B). Cite it. |
| 1708.09838 | 2017 | Cresswell, Peet | Kinematic space of conical defects as a quotient **[abstract]** | Small. Same remark. |

No AI-lab or vendor preprint on the covariant cone turned up. **[abstract-level searches]**

## 2. Task 2: does the Bonahon-current argument exist earlier?

I found no paper, in math or physics, stating any of the following:
- that min_c Σ w_j ℓ_{X_j}(c), or min_c i(μ,c) for a geodesic current μ, satisfies contraction-map
  inequalities;
- that HRT lengths in Mess spacetimes equal i(½(L_{X_L}+L_{X_R}), c) and should be used for entropy
  inequalities.

Specific checks:
- **Ju–Czech 2611.21726: none** **[read full text]**. It is purely about combinatorics of the RT cone.
- **Grimaldi–Headrick–Hubeny–Shteyner 2601.09987: none** **[grep of full text]**. No "Crofton", "kinematic",
  "Liouville", "intersection number", "wormhole" or "quotient".
- **Hernández-Cuenca and Bao et al.: none.** Covers 1505.07839, 1903.09148, 2102.07535, 2109.01150 and 2506.18086
  **[abstracts/metadata; 1505.07839 from memory]**. These are graph or combinatorial methods, plus topological
  link models; none uses currents.
- **Mazzoli–Viaggi 2206.06946 and 2212.11106: no entropy content** **[abstract]**. They supply the pleated-surface
  geometry and the "length spectrum dominates pleated surfaces" statement, which our draft already cites.
  Citing papers (INSPIRE and Semantic Scholar): 2610.22329, 2508.14842, 2510.21334, 2411.10352, 2405.01912, all
  in pure math. None is holographic.
- **Geodesic-current literature, 2022–Nov 2026.** All arXiv hits for "geodesic current(s)" / "Liouville
  current(s)" were checked **[titles/abstracts]**, including:
  - Martínez-Granado–Sapir 2610.03359
  - Cremaschi–Martínez-Granado 2610.30420 (intersection form on infinite-type surfaces)
  - Torkaman 2604.05174
  - Sapir's Thurston-metric extension 2210.08130
  - Hensel–Sapir 2109.14768

  None touches minimal cuts or entropy. Nothing in the Thurston/maximal-stretch literature is related.

**The real risk is Martínez-Granado–D. Thurston, "Lengths as intersections on surfaces", 2605.04031 (May 2026)
[read full text: intro, grep].** Their results:
- **Theorem A:** a curve functional is i(μ,·) for a geodesic current μ iff it is additive and satisfies the
  *smoothing* property, i.e. it does not increase under surgery at essential crossings.
- **Corollary B:** every such functional is the length function of a length pseudo-metric on S.

Applied to μ = Σ w_j L_{X_j}, a referee can argue as follows:
- ½(ℓ_L+ℓ_R) is the length spectrum of a pseudo-metric on S, so "HRT for vacuum wormholes equals RT in a
  virtual static metric" would follow from Bao et al.'s static argument.
- Alternatively: Lemma G plus density of weighted multicurves in currents (Bonahon) plus continuity gives
  Theorem 1 in a few lines.

Limits of that objection:
- MGT is stated for **closed** surfaces; ours have funnels and arcs ending on the boundary with cut-offs.
- Their pseudo-metrics can be badly behaved: not the standard topology, possibly discontinuous (their
  caveat citing DRMG25).
- Neither route gives Theorem H, which is the physics input.

So MGT does not scoop us. But it, together with Bonahon 1988, makes Theorem 1 look like a short observation to
anyone who knows currents. Recommended: cite MGT, give the "virtual metric" viewpoint ourselves, and put the
novelty in (a) the Theorem H identification with the acausal-slice homology condition and (b) arcs, funnels,
cut-offs and endpoint balance. Do not sell the integration step as deep. The smoothing property also answers the
self-intersection issue in GGHH's footnote 10: minimisers of i(μ,·) can be taken embedded.

Not checked: Glorieux–Monclair and Bonsante–Seppi on AdS length spectra ½(ℓ_L+ℓ_R) **[memory only]**. The
draft already cites Bonsante–Seppi. I know of nothing there about cuts or inequalities.

## 3. Task 3: where HEC/HRT-inequality papers were published

Publication data from INSPIRE metadata **[abstract/metadata level]**.

| Paper | Venue |
|---|---|
| **PRL** | |
| Czech–Shuai–Wang, "Entropy Inequalities Constrain Holographic Erasure Correction", 2502.12246 | PRL 135, 141603 (2025) |
| Balasubramanian–Kang–Cummings–Murdia–Ross, "Purely GHZ-like Entanglement is Forbidden in Holography", 2509.03621 | PRL 136, 031602 (2026) |
| (Adjacent, quantum not holographic) "Tip of the Quantum Entropy Cone", 2306.00199 | PRL 131 (2023) |
| **Letter-type, not PRL** | |
| Czech–Shuai–Wang–Zhang, infinite families / entanglement-wedge-nesting topology, 2309.15145 | PRD 109, L101903 (2024) |
| Czech–Shuai, average-entropy cone, 2112.00763 | Commun. Phys. 5, 244 (2022) |
| **HRT-specific predecessors (none in PRL)** | |
| Wall, maximin, 1211.3494 | CQG 31 (2014) |
| Callan–He–Headrick, 1204.2309 | JHEP 06 (2012) |
| Czech–Dong, 1905.03787 | JHEP 10 (2019) |
| Caginalp, Vaidya, 1905.00544 | PRD 101 (2020) |
| Rota–Weinberg, 1712.10004 | PRD 97 (2018) |
| GGHH, 2407.07165 | JHEP 01 (2025) |
| Minimax, 2502.09894 | JHEP 05 (2025) |
| Covariant bit threads, 2208.10507 | JHEP 07 (2023) |
| Grimaldi–Headrick–Hubeny, 2508.21823 | SciPost 20 (2026) |
| Czech et al. majorization, 2601.09989 | SciPost 21 (2026) |
| Bao–Mezei late-time, 1811.00019 | unpublished |
| Zhao, 2602.04888 | unpublished |
| GHHS, 2601.09987 | unpublished |
| **Foundational RT-cone papers (none in PRL)** | |
| Bao et al., 1505.07839 | JHEP 2015 |
| Hayden–Headrick–Maloney, MMI, 1107.2940 | PRD 87 (2013) |
| Hernández-Cuenca, 5 regions, 1903.09148 | PRD 100 (2019) |
| Marolf–Rota–Wien, 1705.10736 | JHEP 2017 |
| He–Hubeny–Rangamani, superbalance, 2002.04558 | JHEP 2020 |
| Hernández-Cuenca–Hubeny–Rota, 2204.00075 | JHEP 2022 |
| Czech–Wang, N=7, 2209.10547 | JHEP 2023 |
| Czech–Liu–Yu, 2401.13029 | SciPost 2024 |
| Bao–Furuya–Naskar, 2409.17317 | JHEP 03 (2025) |
| Bao–Furuya–Naskar, 2506.18086 | JHEP 12 (2025) |
| Akers–Hernández-Cuenca–Rath, 2108.07280 | JHEP 2021 |
| He–Hernández-Cuenca–Keeler, 2312.10137 | CMP 405 (2024) |
| Avis–Hernández-Cuenca, 2102.07535 | Discrete Appl. Math. (2023) |
| **Recent, not yet published** | |
| Ju–Czech, 2611.21726 | PRL-formatted; probably aimed at PRL |

Blunt assessment:
- **PRL has taken HEC papers only recently,** and only ones with a crisp physics punchline that needs no
  technical setup: "inequalities constrain erasure correction", "GHZ entanglement is forbidden".
- **Every proof-type paper on the HRT-vs-RT question went to JHEP or CQG,** including Czech–Dong, which is the
  closest precedent and also claimed a full 2+1 resolution.
- **Points in our favour for PRL:**
  - We close an open question stated verbatim by the Headrick–Hubeny group, the likely referees, and repeated in
    an October 2026 Strings review.
  - Combined with Bao–Furuya–Naskar (JHEP 12 (2025) 140: for rational-coefficient HEIs, contraction maps are
    necessary as well as sufficient), "all contraction-map inequalities" can be read as "the full RT cone".
    Cite it as their published claim, not as settled.
- **Points against:**
  - Vacuum only: locally AdS3, so there are no propagating degrees of freedom. A referee will say the HRT problem
    reduces to two hyperbolic surfaces and is "kinematic".
  - The definition depends on the slice: non-smooth acausal slices are allowed, and the smooth, uniformly
    spacelike version is not proved for regions with endpoints. Referees will go after this.
  - The core is mathematics (currents, pleated surfaces), which is hard to compress into four pages.
  - Provenance is AI-assisted with no independent human check. Ju–Czech shows disclosure is accepted; it is not
    a bar, but referees will scrutinise more.
  - MGT 2605.04031 makes Theorem 1 look short (section 2).
- **My call:** JHEP or CMP for the full paper is the safe venue. A PRL letter with a long companion paper is
  defensible only if the headline is stated narrowly ("HRT entropies of all vacuum 2+1 spacetimes, of any
  topology, obey every contraction-map HEI") and the smooth-slice caveat is closed or stated up front. Even then,
  expect a "too specialised / follows from known math" referee report.

## Corrections to make in our documents

1. Correct "argument for π₁ = ℤ" to: GGHH call it an argument and decline to claim rigour, while the 2026
   review lists it as a proof. Never describe it as unproved. Also note that GGHH's π₁=ℤ result and Czech–Dong's
   simply connected proof hold for **non-vacuum** NEC bulks; ours does not.
2. State that GGHH's numerics already assume the Theorem H formula; Theorem H is its rigorous justification.
3. Cite Zhang–Chen 1610.07134 and Cresswell–Peet 1708.09838 (quotient kinematic space), Martínez-Granado–Thurston
   2605.04031 (currents ⇔ smoothing ⇔ pseudo-metrics) and Ju–Czech 2611.21726 (non-superbalanced HEIs: our
   theorem covers them, the majorization evidence does not).
4. Sahu 2610.07176 is the only fall-2026 paper using Crofton/kinematic space for entropy inequalities. It is
   static, pure AdS3, excludes BTZ quotients, and concerns SSA of generalised wedges. No overlap; cite it
   optionally.
