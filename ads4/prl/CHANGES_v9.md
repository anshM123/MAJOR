# Changes v8 → v9 (24 Nov 2026): cheap fixes after `REFEREE_REPORT_4.md`

No physics was recomputed. Only writing, one figure script, and the Supplement's organisation changed. Nothing committed.

## Letter (`letter.tex`)

| report item | change |
|---|---|
| M1.1 title | "Einstein Gravity from Conformal Field Theory" → **"How Close to Einstein Must Holographic Gravity Be?"** |
| M1.2–5, m3, m2 (abstract) | The abstract now opens with the assumption ("If fixed light-ray functionals … are positive on heavy states"). It leads with \|ĝ₃\|Λ⁴ ≤ 120 in gap units and gives τ ≤ 2.89·10⁻⁸ as "i.e.". "Certified in flat space" now reads "proved in the flat-space limit for the cubic bound and supported numerically in AdS". "Hold up to one universal constant" now reads "The numbers are for c₁ = 0, where c₁ is one universal constant we have not computed". The rhetorical question and the unclear "its" are gone. The headline is stated for a gap in the single traces of T×T (see M5). ABJM is kept, as instructed. 580 rendered characters. |
| M4.1 Introduction | "two conditional theorems that rest on one technical assumption" now names the main premises in plain words: a decay assumption beyond the Regge regime and, near the gap, an edge integral of the graviton that we could not compute, whose scalar analogue meets it narrowly. |
| M6 Introduction | "establish it … for the low-lying heavy states" → "check it numerically, with error control, on the most constraining state and a sample of other low-lying states". "Most low-lying states" was not used, because only 14 of the 30 spins J ≤ 31 are error-controlled. |
| M4.3 Theorem L′ | "on every state of spin below 10⁶ …" → "the leading margin is certified on every state and the remainder was checked on sampled states up to mass 5Λ". |
| M4.2 Theorem W | "every sampled near-gap state of spin J ≥ 26" → "the 47 sampled near-gap states of spin 26 ≤ J ≤ 344, given the edge bound below". "Both rest on one technical assumption" → "Both rest on an assumption". |
| M4.5 | "prove order by order … and derive at all orders" → "prove through eighth order …, check numerically to fortieth order and derive at all orders". This follows S Sec. S7.A: [P] to O(r⁸) for integer b = −4…4 (O(r⁶) with b symbolic), [DV] to O(r⁴⁰), [D] at all orders. |
| M5 | The Result is now stated under "(i) with no light single trace in T×T", where implication (1) is [D]. Light single traces of spin 0 and 2 follow as an extension in one sentence that rests on the sketched step; (2b) moves into that sentence. (2a) keeps only "τ = 0 or κ > 0", because the scalar condition is automatic without light single traces. The Functionals paragraph now says "so, given the light-side formula we sketch, they enter with the favourable sign and the cubic bound also bounds their couplings" instead of "so the cubic bound is unaffected". "Why it matters" and Discussion use the same scope ("gap in the single traces of T×T"). The Discussion's "(for the quartic one, without light scalars …)" is no longer needed. |
| m4 | "up to a small AdS mass shift" → "up to an AdS mass shift we have not computed; we only sketch this step". |
| M6 J = 29 | J = 29 is error-controlled only at 1.04–1.10Λ, sign-free with L ≥ 1.5 (SMALLNU_WINDOW table); its twist edge needs the measured sign (sign-free L = 0.94). The text now reads "Three of these spins are covered only above their lowest mass; at J = 29 the lowest mass is within error control only with the measured sign of one low-frequency contribution". |
| M6 "exactly its flat value" | → "its flat value, to about one part in 10⁶" (ε_W ≈ 1.1·10⁻⁶ [DV] for k = 2). |
| M6 "We show" | → "We derive" (H_K [D]). |
| M6 reproduction sentence | → "Separate code reproduces every constant of the bounds, and a separate check rebuilt the EE-kernel reduction and the counting lemma from the written statements alone." |
| M1 / Why it matters | "up to one universal constant c₁" → "at c₁ = 0, where c₁ is one universal constant that we have not computed". |
| Fig. 2 caption | States that the first theorem covers only the 47 sampled states (dots) and that the hatched band between them is interpolated. Notes that the pale J = 29 segment needs the measured sign. |
| header comment | v9 line added. Internal comments (m12) were not stripped; do that before arXiv. |

## Fig. 2 (`fig/make_fig2.py`, M4.4, m8)

* Theorem W is no longer a solid region. The 47 sampled states (J, m/Λ) from `code/agproof/cache/chi.txt` are dots, drawn over a white
  band hatched in the model colour. The band is labelled "theorem 1 (sampled states)", with a legend entry "between theorem-1 samples".
  Theorem L′ stays solid and is drawn over the band where the two overlap (J ≥ 223).
* J = 29: the error-control bar now runs 1.04–1.10Λ. The twist-edge-to-1.04 segment is pale (model colour), because it holds only with
  the signed window.
* The script still prints J\* = 223 at the twist edge (check passes). The bars for J = 27, 29 and 31 still sit close to the W band (m8,
  partly addressed).

## Supplement (`supplement.tex`)

* **Title**: "Supplemental Material for ``How Close to Einstein Must Holographic Gravity Be?''".
* **M7 (cheap part)**:
  * Version history is removed from the main flow:
    * S1: "The Letter's wording (v7, v8)" → "The Letter's wording", with its map updated to the new Letter phrases.
    * "What moved here from the Letter (v4) … New in v5 … In v6 …" → a short "What is here".
    * The provenance paragraph no longer narrates v6/v7/v8. It lists what was independently checked and what is unreviewed.
    * "What changed since the previous draft" in S3 is deleted; its grade changes are moved to the end.
  * Inline "v7/v8", "since v7", "(v8)", "earlier drafts", "withdrawn" and "superseded" remarks in S2–S8, the tables and S7 are reworded
    neutrally or dropped. Examples: "The v7 premise (AG′) and why it is replaced" → "The pointwise premise (AG′) and why it is not
    used"; the v7-era margin statements; "P_L′(ii) of v6".
  * The "Record fixes" paragraph and the former "record corrections" list are collected in a new final section, **"Corrections to
    earlier drafts"** (`\label{sec:corr}`), headed "for readers of earlier versions only". Its body is unchanged except for the
    additions above. It still contains "In v5 … In v8 …" sentences, which are appropriate in a corrections list.
  * Not done (not cheap): rewriting S5/S7 as derivations, and removing the composite grade notation.
* **M5**: S5.F now separates the two representations.
  * The total light value of an exchange (single trace plus double-trace dressing) is, for the unsubtracted spinning k = 2, 3
    functionals with [TT] at double zeros, the single-trace block's action and a Regge polynomial [D given (E′)].
  * The e^{πν/2} growth concerns a single block *alone* in CHMRS's *subtracted scalar* representation, where the double traces are not
    annihilated and cancel that growth in the total [D-sketch] (DSKETCH_UP §3.3).
  * "Consequences" is now headed "(given (E′))".
  * S1 says that the Letter's Result is stated without light single traces, and that they form an extension resting on (E′).
  * The net-grade line for item (1) says "[D] if no light single trace couples to TT (the Letter's headline)".
* **m9**:
  * S2.B now mentions the ×100 bump of the continuum k = 2 weight at Λ + 30–33 (from the last Gaussians; Re χ is 2·10⁻²² at Λ+28 and
    2.4·10⁻²⁰ at Λ+32 at Δ_g = 1000, figures from the referee's 60-digit evaluation). It is harmless and included in the tail bound.
  * S2.C says that the grids, truncated or not, are tabulations of the continuum weight, not separate functionals: compact ν-support
    fails H_K by the no-go of S5.D. "Both give the same bound" is replaced by "on ν ≤ Λ the tabulation gives the same light
    coefficient".

## Not fixed (not cheap, or out of scope)

* M2, M3: positivity across the spectrum, the graviton EDGE-R, c₁, and a range of gaps. These need computation.
* M7 beyond the cheap part: a derivation-style rewrite of S5/S7, an external review, and the composite grade tags.
* m1/m2 shortening of "Large spin" (it grew slightly). m5 (x ∝ |c₁|/Λ²). m6 (binding-state AdS correction −2·10⁻⁵). m7 (Fig. 1 marker
  code). m9's "certified" qualifier for the leading-Regge floor of the defined (continuum) weight. m10–m13. m14 (AI-tool disclosure:
  the authors must decide; no acknowledgments were added).

## Build and length

* Build: `pdflatex` / `bibtex` ×3 per NOTES.md §1, into the git-ignored `prl/build/`.
  * Letter: 5 pages. Supplement: 25 pages.
  * Both: 0 errors, 0 undefined references, 0 overfull boxes.
* Length (`wc_prl7.py letter.tex 3 0.9375 1.0385`):
  * prose 2850 words and 215 inline formulas (268 conservative);
  * displays 48; figures 344;
  * **TOTAL ≈ 3457 (3510 conservative)**, against a limit of 3750.
* Abstract: 580 rendered characters, against a limit of 600.
