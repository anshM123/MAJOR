# CHANGES v6 → v7 (large-spin results since v6; plain-prose Letter)

*23 Nov 2026. Writing pass. No physics recomputed: every number is quoted from the deliverables below, after checking it there
(J\*(m) recomputed from the formula; Theorem W's 47 states and failing lists read from `code/p4win/cache/states15.txt`; the
O(1/m²) cost fractions from `code/p2con/cache/jstar_full.txt`). Sources: the update boxes of `ads4/UNIFIED_THEOREM.md` through
EXT_FOLD, and `INDEP_CHECK2.md`, `LRINF_TIGHT.md`, `P2_CONTRACT.md`, `P4_WINDOW.md`, `EXT_FOLD.md` (incl. its batch-2 addendum on
J = 31), `REFEREE_REPORT_3.md`, an outside critique (grade tags read like an audit log; positivity paragraph a wall of numbers) and
the authors' instruction to make the Letter read like a PRL. Files edited: `letter.tex`, `supplement.tex`, `macros.tex`,
`fig/make_fig2.py`, `fig/make_fig.py` (two legend labels) and their outputs, `NOTES.md`; new: this file. Nothing outside `ads4/prl/`
touched (`AG_PROOF.md`, `code/agproof/` not touched); nothing committed. L = letter.tex, S = supplement.tex; line numbers are those of
this version.*

## 0. Published PRLs read for voice (authors' instruction 1)

| Letter | what we took |
|---|---|
| Guerrieri, Penedones, Vieira, *Where is String Theory?*, PRL 127, 081601 (2021), arXiv:2102.02847 | opens with one plain sentence about the physics question before any technique; a continuous argument with no section numbers; rigour signalled by verbs ("a rigorous bound", "we estimate", "seems to") rather than labels; captions state the conclusion; numerical detail in appendices |
| Gopakumar, Sinha, Zahed, *Crossing symmetric dispersion relations for Mellin amplitudes*, PRL 126, 211602 (2021), arXiv:2101.09017 | run-in paragraph headings; almost no numbers in the main text (results symbolic, numerics in the supplement); assumptions named in the sentence that uses them ("assuming its nonperturbative existence"); a short forward-looking close |
| Bellazzini, Lewandowski, Serra, *Positivity of amplitudes, weak gravity conjecture, and modified gravity*, PRL 123, 251103 (2019), arXiv:1902.03250 | the main result stated once, as displayed inequalities, then its physical meaning in prose; limitations phrased as explicit conditions ("we content ourselves with assuming …"); candid self-assessment in the conclusion |
| Camanho, Edelstein, Maldacena, Zhiboedov, arXiv:1407.5597 (JHEP; introduction only, as suggested) | physical framing of the problem in short declarative sentences ("This causality violation is occurring when the theory is still weakly coupled"); assumptions listed in words, essentially no numbers |

Applied: run-in headings kept but phrased as the narrative (Introduction → Why it matters → Setting → Mechanism → Functionals →
Result → "Is the functional positive?" → Large spin → A check with ABJM theory → Discussion); status words instead of tags; only the
headline numbers in the text (τ₀ = 2.89·10⁻⁸ i.e. \|ĝ₃\|Λ⁴ ≤ 120, the R⁴ upper end 92.8 / 164.5, the ABJM 2ζ(3) at 2.6 %), plus the
few thresholds a reader needs (0.46Λ, |c₁| ≈ 60/1300, J = 26/42, "a third of the AdS radius"); captions say what the figure shows.

## 1. Tasks

| item | change | where (v7) | source |
|---|---|---|---|
| **1. Plain prose, no tags** | All bracketed grades removed from the Letter (abstract had none; intro, Setting, Mechanism, Functionals, Result, positivity, large spin, ABJM, Discussion, both captions). Each status now in words: "we prove", "certified", "we derive"/"derive and check numerically", "we sketch the argument"/"only sketched", "with error control", "computed inside models or interpolated", "assuming", "a model calculation". One sentence points to the legend (L:82–87). The internal labels A1′/A2 → "(i) a higher-spin gap", "(ii) weak coupling below the gap"; "(Pos)" → "positivity"; H_LJ, W/L′, AG′, [NE±], K₂ absent from the main text and from Fig. 2 (legend in words). New macro `\HLJstatusLetter` (plain status sentence, follows the H_LJ(300) switch). Supplement §S1 gains a paragraph mapping every Letter phrase to its grade and naming (i), (ii), "the two conditional theorems", "the technical assumption", "three magnitudes", "a hypothesis supported by a sketched argument" | L:73–87, 103–117, 150–163, 179–209, 211–252, 255–303, 307–337; macros.tex (`\HLJstatusLetter`); S:65–73 | critique (a); authors' instruction 3 |
| **2. Positivity paragraph** | now "Is the functional positive?", six sentences around Fig. 2: certified (flat limit, both functionals; leading-Regge part for impact parameter ≳ a third of the AdS radius); error control (binding spin-5 state; thirteen further near-gap spins up to J = 31, two only above their lowest mass, J = 29 thin and dependent on the measured window sign; heavier states at small b; Ω^Q_up only at 1.3–3Λ); conditional theorem (large spin, next paragraph); model/interpolated (everything else, including Ω^Q_up's binding states); no negative state, post-hoc disclosure; Δ_g = 300 status. All lists and budgets in Tables S5–S7, S11 | L:255–267 | critique (b); task 2; EXT_FOLD §0 + addendum |
| **3a. Theorem L′** | J\*(m) = ⌈(0.37211 m − 1)/2⌉ (223 at the gap; 250 / 307 / 279 / 372 / 929 at 1.1, 1.3 m_gap, 1.5, 2, 5Λ), Lemma LR∞ [P + certified] with ∫\|ψ₂\| = 10.386029907129, ∫\|ψ₃\| = 12.835622283960, certified r table, β ≥ 64 tail 2.03·10⁻⁷, F_LR > 0 to β = 0.245; tolerances with the certified floor (Λ_tol 181 at J\* → 32 at 10⁶); ε_R ≤ 0.086 via Lemma P4W; grade [D given P_L′(i), AG′, EDGE]; P_L′(i) [DV] at J ≤ 344, [M] above | S:1041–1073, 1191–1245 (Table S10); Letter L:271–291 ("every state … whose impact parameter exceeds about a third of the AdS radius", "leading margin is certified") | LRINF_TIGHT §0–§6; P4_WINDOW §5, §7 |
| **3b. Theorem W** | window form: premises P0, P1′, W, AG′, EDGE; 47 sampled states (21 spins); d/H ≤ 0.49; sampled J\* = 26 (plain edge envelope), 42 with the Ω^U +− factor (Table S9 lists 51/100/46 variants); k = 3 piece included; computed O(1/m²) constants (Λ_B(+−) ≈ 2, \|Λ_B(++)\| ≤ 0.09, Λ_C2 = 1, C1(+−) ≈ 200× harmless, Λ_A ≈ 6–8 normalisation; off-diagonal 0; O(1/m) image 0 ⇒ P1′ [DV] at leading-Regge kernel level; monopole cancels ∓12.8/±12.8 so PW-R1 holds at O(1/m²); cost ≤ 0.41 H at J = 26, 0.05 at J ≥ 64); K_sh premise removed; the flat tower and the 3/16 residual. Not written: "proved for J ≥ 26/42" | S:1074–1099 (§S7.C), 1100–1135 (§S7.D, Lemma P4W and (AG′)), 1136–1190 (§S7.E, Table S9); L:271–291 | P2_CONTRACT §0–§4; P4_WINDOW §0–§8; INDEP_CHECK2 §4.4 |
| **3c. (AG′) as an assumption** | stated precisely in §S7.D with its status ([D] to all orders of the saddle expansion; local rate [NE] at x ≤ 0.86, [M] at x > 0.94; beyond all orders [A]); in the Letter in plain words: "Both rest on one technical assumption … we have not proved it beyond all orders". EDGE and W named in the Letter as "three magnitudes that we estimated only in models". `% AGPROOF` marks at every place that changes if (AG′) is proved (list in §3) | S:1123–1135; L:278–285 | P4_WINDOW §6, §8 |
| **3d. W and L′ overlap; "345–857 uncovered" removed** | Letter no longer mentions the band; Table S11 rows "26 ≤ J ≤ 344" and "223 ≤ J ≤ 344: also Theorem L′"; §S7.G sentence that the v6 band lies inside L′; Not-claimed list (§S11) no longer lists it | S:1257–1262, 1270–1271, 1401–1403 | LRINF_TIGHT §6 |
| **3e. Fig. 2** | `make_fig2.py`: Theorem L′ boundary J + ½ ≥ 0.186055 m (both panels; dotted white line where it overlaps W near the gap); Theorem W from J = 26 (white dashed line at 42); new error-control bars J = 7 (to 1.05Λ), J = 27 (to 1.10Λ), J = 23 (1.03–1.10Λ only), J = 31 (1.05–1.10Λ only, EXT_FOLD addendum); "345 ≤ J < 857: no theorem" bracket removed; star at (1100, 100) removed (premise now computed at all 47 states); hatched "hypothesis only" now m ≥ 1.3Λ, b_A ≤ b < 0.37211Λ; legend and labels in words ("numerics with error control", "conditional theorem", "theorem 1/2", "model / interpolated"). Script prints the L′ twist-edge J\* as a check (223). `make_fig.py`: labels "this work", "flat D = 4 optimum, transplanted", "(2b): light scalars allowed" | fig/make_fig2.py, fig/make_fig.py; captions L:166–177, 292–303 | as above |
| **4. INDEP_CHECK2 fixes** | (i) CHMRS's ε quote cited at their Sec. 1 (Table S2 row A2); (ii) "technical assumption" cited as CHMRS eq. (5.3) (row L7(b)); (iii) "the logic of CLS-D" → "CLS-D restate the same light/heavy split (Sec. 2.2.2) but do not discuss the order of limits" (§S5.D); the Letter says only "Like Refs. [CHMRS, CLS-D], we use the sum rule at finite c_T"; (iv) LR∞'s constant named (∫\|ψ\| certified; v6 table used Σ\|c_n\|/(n+1)); (v) J\* numbers (L′ 223; W 26/42 window form, k = 3 included; INDEP_CHECK2's 64 recorded as the v6-form value); (vi) magic identity: [P] with b symbolic to O(r⁶), every integer b = −4…4 to O(r⁸); (vii) monopole t-independence derived via 1/c partial fractions, D = (15/2)δ₁ | S:246, 263, 584–585, 1009–1030, 1049–1060, 1155–1160; Table S12 (S:1337–1341); S:1423–1430 | INDEP_CHECK2 §1–§4, §6 (D1–D7) |
| **5. v6 fixes kept** | order of limits (Mechanism L:135–139, Functionals L:186–192; S §S5.D); B_eff(Λ) (L:232 "varies slowly with the gap", values in Table S4 and §S4); framing (Λ⁻⁴ scaling, leading order in 1/c_T, no 10⁻⁸ claim at finite c_T, L:88–101); front matter unchanged | as listed | CHANGES_v6 |
| **6. Length** | `scratchpad/tools/wc_prl7.py letter.tex 3 0.9375 1.0385` (wc_prl6 + expands `\HLJstatusLetter`): prose 2715 words, 212 inline formulas (conservative 264), 3 display rows (48), Fig. 1 180, Fig. 2 164: **≈ 3319 APS words (≈ 3371 conservative)**; limit 3750. Abstract unchanged, **563** rendered characters. Text ends ≈ 0.55 down the right column of page 4; Letter 5 pages, Supplement 23 pages | — | — |
| **build** | both compile into `build/` with 0 errors, no undefined references or citations, no overfull boxes; Letter 7 underfull notices; Supplement: hyperref "token not allowed" (math in section titles), cosmetic | — | — |

Other Supplement updates (EXT_FOLD): Table S5 new J = 7 row; J = 23, 27, 31 row rewritten (direct bulk + GK edge values, budgets,
grades); J = 29 row (budget 0.784, L = 1.24 signed; sign-free 0.94 set by the ν = 20 node and the convention below ν = 10; 23
independent streams); Table S6 direct values for J = 7, 23, 27, 31; Table S7 rows H_M (K₂ error 0.6–8 % on seven states, was 4–5 %),
H_edge (J = 31 now GK), H_LJ(1000) (replaced by Theorem W on its samples), H_R′ (range now 5 ≤ b < 0.37211Λ). Reproduction table:
five rows from INDEP_CHECK2. Record corrections "In v7" (§S11).

## 2. Physics wording choices to note

* "Two conditional theorems" in the Letter are Theorems W and L′; Fig. 2 labels them "theorem 1/2". The Letter never says that
  positivity is proved for J ≥ 26 or 42; it says "positivity on every sampled near-gap state of spin J ≥ 26 … Both rest on one technical
  assumption … three magnitudes that we estimated only in models".
* "Less than half of the margin" = P2_CONTRACT's computed O(1/m²) ≤ 0.41 H and P4_WINDOW's window + boundary ≤ 0.49 H.
* "Impact parameter exceeds about a third of the AdS radius" = β ≥ 0.37 (b ≥ 0.37211Λ); F_LR > 0 alone is certified from β = 0.245.
* "Thirteen further near-gap spins up to J = 31" = 0, 2, 4, 6, 7, 8, 12, 16, 21, 23, 27, 29, 31 (23 and 31 only above their lowest mass).

## 3. `% AGPROOF` marks (places that change if (AG′) is proved)

letter.tex: L:7 (header), L:81 (intro, "rest on one technical assumption"), L:285 (large-spin paragraph: the assumption sentences),
L:303 (Fig. 2 caption), L:335 (Discussion, next steps). supplement.tex: S:17 (header), S:73 (§S1 wording map), S:991 (§S7 intro),
S:1123 (statement and status of (AG′)), S:1161 (Theorem W grade), S:1212 (Theorem L′ grade), S:1257, S:1262 (Table S11 rows), S:1271
(not claimed in §S7), S:1403 (§S11 not claimed). If (AG′) is proved, Theorems W and L′ keep the magnitude premises (W, EDGE, and for
L′ P_L′(i) above J = 344) and, for W, P0 below b = 245 and P1′ beyond the leading-Regge kernel.

## 4. Inconsistencies found

1. **LRINF_TIGHT §6** states J\*(1.3Λ) = ⌈(0.37211·1.3Λ − 1)/2⌉ = 241; the formula gives 242 (0.37211·1.3·998.49987 = 483.02).
   Immaterial (not quoted in the paper).
2. **Number of Theorem W states:** v6 said "40 sampled states (even J 26–344, odd 33–301)"; P2_CONTRACT, P4_WINDOW and INDEP_CHECK2 use 47
   (J, m) states on 21 spins (`states15.txt`). v7 says 47.
3. **Theorem W's J\* across deliverables:** INDEP_CHECK2 (64 with k = 3, v6 full-range premise form), LRINF_TIGHT §6 (repeats 64),
   P2_CONTRACT (26/30/38/50 by K₄; "no J\*" with the full-range series), P4_WINDOW (window form, k = 3 included: 26, or 42 with the
   Ω^U +− factor). The window form supersedes the others; the paper quotes 26/42 and records 64 as the v6-form value.
4. **EXT_FOLD verdict vs addendum (J = 31):** verdict table: [NE±] only at 1.075/1.10Λ with L = 1.53/1.91; addendum (batch 2): [NE±]
   from 1.05Λ, L = 1.22/1.55/1.90 at 1.05/1.075/1.10, twist edge +0.862, L = 0.70 (0.85 signed). The brief (before the addendum) listed
   "23/31 at higher masses"; v7 uses the addendum.
5. **Magic-identity residual:** PW_IMAGES (v6 text) 7.9·10⁻⁴⁸ (b ≤ 4), 2.0·10⁻⁴⁷ (b = 7), 1.1·10⁻⁴⁷ (b = 10); INDEP_CHECK2 ≤ 3.7·10⁻⁴⁸ for
   all b incl. 7, 10 — different normalisations (INDEP_CHECK2: relative to max\|LHS\| at the same order). Both quoted with their source.
6. **Theorem L tolerance:** v6 Letter "up to 43 times the scalar one" used the numerical floor; with the certified floor it is 32 at
   J = 10⁶ (LRINF_TIGHT §5). v7 Supplement quotes the certified column.
7. **P_L′ magnitudes:** LRINF_TIGHT says "premise status (P_L′ magnitudes [NE]/[M]) unchanged", while P2_CONTRACT computes the O(1/m²)
   corrections (not scalar-shaped) at J ≤ 344 with cost ≤ 0.005 H for J ≥ 200. v7: P_L′(i) [DV] at the overlap 223 ≤ J ≤ 344, [M] above.
8. **K₂ calibration:** HIGHSPIN's "4–5 % of the last Laplace step" (quoted in v6 Table S7) is 0.6–8 % over seven states (EXT_FOLD).
9. **CHMRS labels (INDEP_CHECK2 D5):** the ε-quote and eq. (5.3) labels did not appear in v6's prl files (only REGGE_TREE had them);
   v7 adds both citations with the corrected locations. The "logic of … CLS-D" sentence did appear (S §S5.D) and is fixed.
10. **(3/16) residual:** P2_CONTRACT lists it as unexplained; P4_WINDOW §3.3 resolves it as a genuine beyond-hybrid O(1/m⁴) term
    ([DV] + exact-fit [NE] −0.12…−0.27 vs −0.1875). v7 follows P4_WINDOW.
11. **UNIFIED_THEOREM box for LRINF_TIGHT** says W "sampled to 344" with J\* = 64 and "Still uncovered: J < 64"; with P4_WINDOW, J < 26.
    The UT boxes were not edited (outside `prl/`).
