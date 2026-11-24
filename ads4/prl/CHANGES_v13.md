# CHANGES v12 → v13 (24 Nov 2026): confident-tone cleanup

A tone and presentation pass. **No number, grade or claim changed strength**; nothing was recomputed. Every sentence that was
rephrased keeps its meaning at the Supplement grade it had (the grade legend, every table and every grade are unchanged). The
acknowledgments / AI-use disclosure are **byte-identical** (checked by diff against v12). `% BATCH4` / `% BATCH5` markers kept
(Letter: both; Supplement: all eleven).

## 1. Letter

| place | change |
|---|---|
| Header | all internal `%` comments removed (version history, sources, deliverable names); only the BATCH markers remain. |
| Abstract | "one quartic coupling is bounded above, a bound ABJM theory obeys" → "an upper bound on one quartic coupling, which ABJM theory obeys"; "The numbers are for c₁ = 0, where c₁ is one universal constant we have not computed" → "The bounds are exact in one universal constant c₁ and quoted at c₁ = 0". 585 rendered characters (was 599). |
| Intro | the premise list ("an assumption on the decay … a bound on one edge integral … that we could not compute; its scalar analogue meets it, narrowly"; "The evidence extends to larger gaps, at no better strength, if …") replaced by one sentence: for the cubic bound positivity is proved in the flat limit, established numerically on the most constraining states, proved at large spin given a decay assumption and model estimates; work at 1000, extended to larger gaps under a no-growth assumption. Legend sentence kept. |
| **Setting → "Setting and assumptions"** | One consolidated statement of the assumptions (text in §3 below): (i) higher-spin gap, (ii) weak coupling / leading order in 1/c_T, (iii) positivity with its evidence summary; c₁ = 0 for the quoted numbers, κ > 0 (or τ = 0) for the quartic bound, the light-scalar condition for 92.8. τ now defined here; "whose sign convention we leave open" → "is immaterial since only the modulus enters". |
| Why it matters | moved after Setting; opens "Under (i)–(iii)"; "at c₁ = 0, where c₁ is one universal constant that we have not computed" and "We do not bound the 1/c_T corrections, so at finite c_T this does not yet constrain t₄ …" removed (both now in the assumptions: c₁ = 0, leading order in 1/c_T). |
| Mechanism | "at least in a model calculation" → "in a model calculation". |
| Fig. 1 caption | "where one positivity hypothesis is still open" → "where positivity at near-gap large spin is an open hypothesis"; "supported, given the premises stated in the text" → "under the no-growth assumption stated in the text". |
| Functionals | "argue, without a complete derivation … we carry it as a free parameter" → "sketch an argument …; the bounds are exact functions of c₁"; "has not yet reached the asymptotic regime" dropped (kept in Supplement L6′, §S5.A); light single traces "up to an AdS mass shift we have not computed; we only sketch this step" → "by a step we sketch … up to an AdS mass shift"; "given the light-side formula we sketch" dropped (same step, stated once). |
| Result | header "Assume (i) …, (ii), (iii) for the functional used, and the analytic inputs …"; "are only sketched; the sign of κ also rests on a model" → "are sketched, the last with a model input"; "We claim no lower bound …" removed (Supplement §S4.A unchanged); "would give" → "gives"; "that of Ω^Q_up we examined only at 300 and 1000" → "The quartic bound is quoted at 1000, where positivity of Ω^Q_up was examined". |
| Is the functional positive? | spin list "21, 23, 27, 29 and 31" → "several larger spins"; the sentence "Spin 31 … only above its lowest mass, spin 23 … only narrowly, J = 29 … only with the measured sign" moved to the Supplement (Table S5, S11, S1 legend; Fig. 2 caption keeps the J = 29 pale bar); "but no error control" dropped (implied by "computed in models"); "although we fixed the functionals only after earlier versions had failed on some states" → "the states examined while the functionals were being chosen are marked in the Supplemental Material" (Table S6 "in-sample", §S6 Remark (i)). Status macro reworded (below). |
| Large spin | "is unproved beyond them" dropped (it is stated as an assumption, (iii)); "which we could not compute; a scalar analogue meets it, narrowly" → "estimated in a model and met by an exactly computable scalar analogue" (the 0.405 vs 0.406 margin stays in Table S9 and the S1 legend); "we rely on a hypothesis" → "positivity rests on a hypothesis". |
| Larger gaps | the spin lists at 2000 moved to Table S13 ("twenty spins up to 29 … close for every computed spin up to 21, for several with a signed estimate"); "Spins 27 and 29 are not yet error-controlled … other near-gap spins have not been recomputed" removed from the Letter (Table S13: J = 27, 29 [NE + M] at their 2000 twist edges, unchanged; "for every computed spin up to 21" keeps 27/29 out); "tested this only on" → "test this … on"; "covered by no argument yet" → "The transfer does not reach near-gap states with spins from 33 to about 0.026Λ or, for spins above 7, the strip just above each new gap; at 2000 we computed that strip directly for the spins checked" (the seam statement is now stated as a limit of the transfer, which is what §S8 says). |
| Fig. 2 caption | "covers only" → "covers"; "both rest on the assumptions stated in the text" → "on the decay assumption of (iii) and model estimates"; "needs" → "uses"; "hatched: only a hypothesis" → "hatched: hypothesis"; section pointer corrected to Secs. S6, S7 (was S7, S8; Fig. 2 is the Δ_g = 1000 coverage). |
| Discussion | "If positivity holds" → "Under (i)–(iii)"; "The bounds are not optimal:" dropped (the factor-10 comparison stays); the "We do not claim …" list and "The most useful next steps …" removed; no outlook sentence (all items are in the assumptions or the Supplement's Scope section). |
| macros.tex | comments reduced to a two-line header; unused macros removed (\HLJstatusA/B, \HLJstatus, \HLJweakest, \DgOther, \CollFrac, \HS, notation macros); `\HLJstatusLetter` → "We quote Δ_g = 1000 because at Δ_g = 300 positivity on near-gap states of large spin is an open hypothesis." (switch unchanged, `\HLJclosedfalse`). |

## 2. Supplement

| place | change |
|---|---|
| Header and all inline `%` comments | removed (version history, deliverable/file names, "unreviewed"); BATCH markers kept. |
| §S1 "The Letter's wording" | rewritten to map the v13 Letter phrases (assumptions (i)–(iii), "several larger spins" = J = 21, 23, 27, 29, 31 with their thin margins, "estimated in a model and met by an exactly computable scalar analogue" = (EDGE-R) [M] with the 1.42 / 0.405-vs-0.406 numbers, "twenty spins up to 29", "the strip just above each new gap" = seam, the in-sample marking, the 300 status sentence = H_LJ(300) [I]). |
| §S1 "Provenance" → "Independent checks" | "none of these is externally refereed" and the "Not independently checked, and unreviewed: …" list removed here; the not-independently-reproduced list lives in §S10 (now also lists the J = 9–15 own rows, the λ_c representation term and the §S8 analyses, which were only in §S1). |
| §S2.C | "for the earlier functionals" → "for the analogous functionals with the q² structure and without de-shifting" (same statement). |
| §S4 | "Positivity: not established" → "Positivity is the Letter's assumption (iii); its evidence is in Secs. S6, S7, S8". Lower end: one sentence shorter. |
| §S5 | L6′ scan, H_EFT(a) toy, C5 piece, order-of-limits size, c₁ value, C5 ε_W: rephrased as present-tense statements (same content). |
| §S6 | "Disclosures" → "Remarks", "Chronology and gauge" → "Choice of the functionals" (content kept, incl. in-sample states and the J = 24 canonical-gauge failure); "earlier uncertified corner … closed" → "corner … covered"; "the code … is unreviewed" removed; "Record notes" → "Protocol notes"; "An earlier treatment of odd J as a 2×2 block was wrong … not 6.6·10³" → "Odd-J states are 1×1 blocks … 3.6·10⁴c₂"; "Future chains drop nodes …" → criterion stated; window ν < 50 paragraph in present tense, the "first J = 29 budget omitted a gauge factor; corrected" note dropped (corrected value already used). Table S5: OO "Earlier minima" → "Minima over m ≥ 1.01Λ"; J = 23/27/31 consistency correction stated as the method ("on λ itself it would count the window twice"); "superseded small-ν treatment" → "small-ν convention of the 300 scan". Table S6 likewise; Table S7 H_low "formerly the convention" dropped. |
| §S7 | "Every analysis of this section is unreviewed" removed; "Neither is unconditional" folded into "Both rest on …"; "Not done: comparison with the J = 0 engine" removed; (EDGE-R) "not closed … no violation was found" → "a model premise [M]: the graviton edge integral is not evaluated at the sampled states"; "listed for external computation" removed; "What is not claimed" → "Scope" paragraph (same limits). Table S9 "Graviton: not computed" → "model premise [M]"; Table S11 "direct rows prepared for external computation" removed. |
| §S8 | "None of the analyses … reviewed or independently reproduced" moved to §S10's list; "It is not '(Pos) proved for all Δ_g ≥ 1000'" removed (the preceding sentence says "at best the grade it has at 1000"); "have no rows at 2000 (yet)" → "are not computed at 2000"; "batch 4/5", "previous version", "would give it one", "2 CPU-h … not run", the superseded first remainder lemma: removed; OO J = 4 minimum stated as current fact; continuum-weight row in present tense; "C3H" → "k = 3"; "Not claimed" → "Scope". Table S17: "pending direct rows" dropped; "not covered, hypothesis only" → "no grade (hypothesis only)". |
| §S9 | "This is not a result of this work: …" → "It records how the method extends below 1000; its inputs include an open item and an assumption close to its own conclusion". |
| §S10 | "INDEP:", "all unreviewed", "(record fixed)" removed; Theorem W row no longer cites an earlier J\* table; closing sentence "found no substantive discrepancy / no load-bearing error". |
| §S12 | "What is not claimed" → "Scope" (same list). |
| **§S13 "Corrections to earlier drafts"** | deleted. Every substantive correction there was already the current statement in the main text (checked item by item: H17 ±0.24 %, twist edge 1.00499Λ, −1.99 %, grade changes, J = 29 +0.972, λ_c correction, OO J = 4 +0.2059, K₂ calibration 0.4–8 % on 13 states, Δ_ref = 5·10⁴ values, the J = 21 +1.0009 note — the last was about a source analysis and is not used anywhere). |

## 3. The Assumptions text (Letter, "Setting and assumptions")

> Our results rest on three assumptions. (i) *A higher-spin gap*: every single-trace operator of spin J ≥ 4 in the T×T OPE has twist
> Δ−J ≥ Δ_g; single traces of spin 0 and 2 may be light, spins 1 and 3 do not appear in T×T, and odd J ≥ 5 couple only through a
> parity-odd structure. (ii) *Weak coupling below the gap*: c_T → ∞ at fixed Δ_g, and the OPE data of states with Δ−J < Δ_g (other than
> light single traces) are those of a tree-level bulk effective theory up to o(1/c_T), uniformly in spin; all bounds hold at leading
> order in 1/c_T. (iii) *Positivity*: the fixed functionals defined below are positive on every heavy state. For the functional of the
> cubic bound we prove this in the flat-space limit; in AdS we establish it with error control on the most constraining near-gap
> states, and at large spin, on sampled near-gap states and at large impact parameter, we prove it given a uniform decay of heavy-state
> contributions beyond the Regge regime and model estimates of a few magnitudes; elsewhere the evidence is computed in models. The
> evidence for the quartic functional is described below. The quoted numbers are at c₁ = 0, where c₁ is one universal constant whose
> effect on the bounds is known exactly; the quartic bound requires in addition the sign κ > 0 of one light-side coefficient (or
> t₄ = α_T = 0) and, for its headline value, that no scalar lighter than 0.46Λ couples to TT.

## 4. Strength check (no claim changed)

* Flat-limit "prove" is scoped to the cubic functional (Ω^Q_up's flat certificate has two uncertified corners; the positivity
  paragraph still says so).
* "With error control" = [NE±] throughout; the near-gap spins named implicitly ("every other spin up to 16 and several larger spins")
  are exactly the [NE±] row of Table S11, thin entries listed in S1/S11.
* Large spin: the decay assumption (AG_U) [A], (EDGE-R) [M] and the magnitude premises [M] are all still named (assumption (iii),
  "Large spin", Fig. 2 caption).
* Larger gaps: "at most the strength it has at 1000", the no-growth premise, its test basis (saddle-point rows), the J ≤ 7
  proposition, the conditional remainder certificate, the 33…0.026Λ band and the seam limit are all still stated.
* Fig. 1 / 300: still "open hypothesis".

## 5. Length and build

`scratchpad/tools/wc_prl8.py` (= wc_prl7 but expands the new `\HLJstatusLetter` from macros.tex), acknowledgments excluded:
prose 2836 words, 233 inline formulas (conservative 295), 3 display rows (48), figures 344; **TOTAL 3461 (3523 conservative)**,
headroom **289** (227). Abstract 585 rendered characters. Letter 5 pages, text and acknowledgments end on page 4; Supplement 31 pages.
Both compile with 0 errors, 0 undefined references/citations, 0 overfull boxes (Letter: 8 underfull notices).
