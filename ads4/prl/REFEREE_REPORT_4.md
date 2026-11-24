# Fourth referee report

**Manuscript:** "Einstein Gravity from Conformal Field Theory", v8 (`ads4/prl/letter.tex`, `supplement.tex`, `macros.tex`,
`refs_verified.bib`, `fig/`). Build `prl/build/letter.pdf` (5 pp.) and `supplement.pdf` (25 pp.) of 24 Nov 03:05, newer than every
source file. Switches: `\HLJclosedfalse` (headline Δ_g = 1000), `\HSdraftfalse`.

*Written 24 Nov 2026 as an internal pre-submission review. No existing file was edited and nothing was committed. References are
`file:line` in the current text, relative to `ads4/prl/` unless they say otherwise ("L" = letter.tex, "S" = supplement.tex).*

*Method.* I read the Letter, the Supplement and the figures cold, and wrote §§1–3 and §6 before opening `REFEREE_REPORT{,_2,_3}.md`.
After that I checked report 3 item by item (§4). I then checked the claims against `ads4/UNIFIED_THEOREM_v3.md` and the update boxes of
`ads4/UNIFIED_THEOREM.md`, through NU_TAIL_CHECK, and against the deliverables they cite. Scratch scripts are in
`/tmp/claude-0/-home-user-MAJOR/60eea3a7-07d6-5453-af55-c4346d86f9b5/scratchpad/r4/`. They import no project code, and each ran in under
two minutes on one process.

---

## 1. Summary and recommendation

**What the Letter claims.** Take a CFT₃ with a unique stress tensor and a higher-spin gap Δ_g in T × T (light spin-0 and spin-2 single
traces are allowed), at tree level (c_T → ∞ with tree-level EFT data below the gap). Fixed light-ray functionals on ⟨TTTT⟩ give
explicit bounds, provided those functionals are positive on every heavy state (and on light spin-0 and spin-2 states). At Δ_g = 1000
and leading order in 1/c_T, the bounds are:
* **cubic:** τ = √(t₄²+α_T²) ≤ 2.89·10⁻⁸, i.e. |ĝ₃|Λ⁴ ≤ 120;
* **quartic (2a):** g̃₄Λ⁶/8πG ≤ 92.8, if no scalar single trace lighter than 0.46Λ couples to TT and κ > 0 (or τ = 0);
* **quartic (2b):** g̃₄Λ⁶/8πG ≤ 164.5 without the scalar condition.

Both hold exactly for every value of one universal constant c₁; the numbers above are its c₁ = 0 values. On positivity, the Letter
reports three kinds of support:
* certificates in the flat limit;
* error-controlled numerics on fourteen near-gap spins and a small-impact-parameter strip;
* two "conditional theorems" at large spin.

**Assessment.**
* **The arithmetic and the logic of the implication are sound.** Every headline constant I recomputed from the published weights
  reproduces to the quoted digits (§5). The implication from positivity to the bounds is correct and short.
* **Report 3's main physics objection (M2, the order of limits) is answered properly.** The sum rule is now used at finite c_T, and the
  heavy side only through its sign.
* **The main text is far more readable than v5.** The new coverage figure (Fig. 2) is the right object to show.
* **The v7–v8 work is real progress.** It includes:
  * the computed spinning O(1/m²) corrections;
  * the certified Lemma LR∞, which removes the uncovered 345–857 band;
  * the honest withdrawal of (AG′);
  * the fix of the functional's definition beyond Λ.

The category of the result has not changed, however. The physical content is positivity of a fixed functional on an infinite AdS
spectrum. In AdS that positivity is established at error-control level only on a finite set of near-gap states at one gap value. At
large spin it is reduced to theorems that rest on four kinds of premise:
* a beyond-all-orders assumption (AG_U);
* an edge integral for the graviton (EDGE-R) that was never computed, and whose scalar analogue reaches 0.405 against a tolerance of
  0.406 at the lowest sampled spin;
* several magnitudes taken from models;
* for Theorem W, 47 sampled states.

v8 also *regresses* in the places an editor reads first:
* the new title;
* the abstract's "certified in flat space" and "hold up to one universal constant";
* the Introduction's "rest on one technical assumption".

All three claim more than the Supplement and the internal record support.

**PRL criteria.**
* **Importance and broad interest.** The target is real and current: explicit-constant graviton bounds from ⟨TTTT⟩ are what CLS-D
  list as missing. The EE structure, de-shifting and the k = 4 arc are genuine technical ideas. But the result is conditional and holds
  at one gap value. Its Λ⁻⁴ scaling was known (CEMZ, AKT). Its constant is about 10× weaker in |ĝ₃| than the transplanted flat optimum,
  which CHLPS already remark on. The quartic bound is dominated by the infrared logarithm, and the ABJM check sits at 2.6 % of the bound,
  so it could only fail grossly. For PRL's broad readership the interest is limited.
* **Correctness.** No error found in the implication or the constants. The physical premise (positivity in AdS) is not established, and
  the light-single-trace part of the implication rests on a sketched step (E′) (M5).
* **Claims established at the level stated.** The body is mostly careful. The title, abstract and Introduction are not (M1, M4, M6).
* **Clarity for a general physics reader.** Much improved. The "Mechanism" and "Large spin" paragraphs are still specialist prose (m2).

**Would an editor send it out?** Probably yes, to one or two referees (I would put it at 55–65 %). The new title and abstract are
attractive, the topic is current and the Letter now reads like a PRL. Against it:
* the abstract's own last sentence ("rest on the positivity … supported in AdS") tells an editor the result is conditional;
* the 25-page Supplement is visibly an internal grade ledger with version history;
* both authors are independent researchers.

Some editors will decline on "specialised and conditional".

**What would a strong referee decide?** A referee from the CHMRS/CLS-D community would reject for PRL. They would:
* call the title an overclaim;
* ask for positivity across the spectrum, at least with error control and with an analytic large-spin argument whose premises are
  derived, not assumed;
* ask for c₁ and for the graviton edge integral;
* ask why the light-single-trace case is in the headline when its key step is sketched;
* recommend a long paper.

**Recommendation: reject for PRL in the present form; transfer to PRD (or a JHEP long paper) strongly encouraged.** This is unchanged
from reports 2 and 3. The writing and the large-spin structure have improved, but the result has not changed category. The title and
abstract must be corrected whatever the venue (M1).

---

## 2. Major issues

### M1. The title and abstract claim more than the paper establishes (regression)

**Where:** L:17 (title); L:30-37 (abstract); S:39 (Supplement title).

**Claim.**
* Title: "Einstein Gravity from Conformal Field Theory".
* Abstract: "we show that, at leading order in 1/c_T, its non-Einstein three-point couplings obey √(t₄²+α_T²) ≤ 2.89·10⁻⁸ … Both
  results rest on the positivity of fixed functionals on heavy states, certified in flat space and supported in AdS, and hold up to one
  universal constant."

**Problems.**
1. **The title promises a derivation of Einstein gravity.** The paper delivers conditional upper bounds on one combination of cubic
   couplings and on one of the R⁴ coefficients, at one gap value, assuming an unproven positivity property.
   * The title is also nearly identical to existing titles in exactly this literature:
     * Belin–Hofman–Mathys, "Einstein gravity from ANEC correlators" (cited as Ref. [11]);
     * Afkhami-Jeddi–Hartman–Kundu–Tajdini, "Einstein gravity 3-point functions from conformal field theory" (Ref. [3]).

     Those papers prove parametric statements. A referee will read this title as claiming more than they did, while the paper
     claims less.
   * Report 1 (M12) objected to "Einstein gravity through quartic order" for the same reason. Reports 2 and 3 accepted the narrower
     title "Conditional bounds on R³ and one R⁴ graviton coupling …". The v8 title is a step back past both.
2. **"certified in flat space" is false for the quartic functional.** The sentence attaches it to "fixed functionals" in the plural.
   * For Ω^Q_up the flat certificate fails on two asymptotic regions, where it is [NE] (S:873-875).
   * Ω^Q_up is *proved negative* on flat light scalars below 0.46Λ.
   * Report 3's M4(a) on exactly this point was fixed in v6. The v8 abstract reintroduces it.
3. **"hold up to one universal constant" is ambiguous.** The bounds hold *for every* c₁, in an exact c₁-dependent form. The *number*
   2.89·10⁻⁸ is the c₁ = 0 value and moves by less than 1 % only for |c₁| ≲ 60. As written, a reader can take the number as valid up to
   an unknown O(1) factor, or as valid unconditionally.
4. **The headline number is the scaling-dominated one.** The abstract leads with 2.89·10⁻⁸, which is mostly Λ⁻⁴ (report 3 M4(c); the
   Letter itself says so at L:89-92). The new content, |ĝ₃|Λ⁴ ≤ 120, is not in the abstract.
5. **"we show" precedes the condition, which comes two sentences later.** The quartic sentence also omits both conditions of (2a): the
   scalar condition and κ > 0 (or τ = 0).
6. **"ABJM theory lies inside this bound"** is a consistency check at 2.6 % of the bound, in a regime where the scalar condition holds
   only at sketch level (L:309-311). In an abstract it reads as support. Report 1 (M13) asked for ABJM to be removed from the
   abstract.

**Fix.**
* Retitle. For example:
  * "How close to Einstein gravity? Explicit bounds on graviton couplings from ⟨TTTT⟩ in holographic CFT₃";
  * or "Conditional bounds on the R³ and R⁴ graviton couplings of holographic CFT₃".
* Rewrite the abstract so that it:
  * states the assumption first ("If fixed light-ray functionals are positive on heavy states, …");
  * quotes |ĝ₃|Λ⁴ ≤ 120 (with τ ≤ 2.9·10⁻⁸ as a consequence);
  * says "at c₁ = 0 for a universal constant c₁ not yet computed (the bound moves by <1 % for |c₁| ≲ 60)";
  * restricts "certified in flat space" to the cubic functional;
  * names the scalar condition, or drops the quartic number;
  * moves ABJM to the body.

**Fatal for PRL?** Not alone; it is a few hours of writing. But an overclaiming title is the first thing a strong referee will cite
under "claims not established at the level stated".

### M2. Positivity, the physical content, is still not established in AdS

**Where:** L:72-77, :251-263, :269-288, Fig. 2; S Tables S5–S7, S11, Sec. S7.

**State of the evidence at Δ_g = 1000 (my reading of S:880-893, S:1356-1374, and the EXT_FOLD, AG_PROOF and EDGE_R boxes).**

*With error control (a ×3 safety factor on model terms; an estimate, not a bound):*
* J = 5, the binding state;
* J = 0, 2, 4, 6, 7, 8, 12, 16, 21, 27, 29;
* J = 23 above 1.03Λ and J = 31 above 1.05Λ;
* the strip m ≥ 1.3Λ, b < 5.

Several of these are thin:
* J = 21 and J = 29 have leverage ≈ 1.2 at the twist edge.
* J = 29 reaches 1.24 only with the *signed* window. Sign-free it is 0.94, which fails the Letter's own definition (M6).
* J = 8 is 1.1–1.8 sign-free.
* Most budgets stop at 1.05Λ or 1.10Λ.

*Computed only inside models or interpolated:*
* near-gap J = 9–11, 13–15, 17–20, 22, 24, 25;
* the twist edges of J = 23 and 31;
* all masses 1.1–1.3Λ at J ≤ 48 (H_cont [NE]);
* 5 ≤ b < b_A at m ≥ 1.3Λ;
* *all* binding states of Ω^Q_up, which carry the quartic bound ([NE + M] given H_RS).

*Conditional theorems:*
* Theorem W on 47 sampled near-gap states with 26 ≤ J ≤ 344, under P0, P1′, W, AG_U and EDGE-R (M4);
* Theorem L′ for b ≥ 0.372Λ and J ≤ 10⁶, under P_L′(i), AG_U and EDGE.

*Hypothesis only:* everything between the samples, intermediate impact parameter at m ≥ 1.3Λ, and J > 10⁶.

**Problem.** The result is a bound only if positivity holds on *every* heavy state. Nothing in the paper is a proof of that. The
error-controlled set is a small, non-generic subset, chosen near the expected binding states. The large-spin region is covered by
statements whose premises are assumed or modelled. The Letter is candid about this in its body (L:251-263, :283-288). The claimed
result is still "explicit bounds", while what is established is "explicit values of an implication whose premise is supported on a
finite set". This is the same structural objection as in reports 2 (M1) and 3 (M1, M5).

**Fix.** There is no writing fix. What would change the category:
* error-control every near-gap spin J ≤ 25 over 1.0–1.3Λ, and the Ω^Q_up binding states;
* compute the graviton EDGE-R at the 47 Theorem W states;
* a monotonicity or analytic argument that controls the regions between samples.

**Fatal for PRL?** Yes, together with M3. PRL needs the main claim established. Here the main claim is an implication whose premise is
the hard part.

### M3. Importance and breadth: one gap value, loose constants, one uncomputed constant

**Where:** abstract; L:85-96, :236-246, :317-335; Fig. 1.

**Problems.**
* **One gap value.** The bound exists at Δ_g = 1000 only.
  * Δ_g = 300 has an open hypothesis (L:263).
  * The Λ-forms are values of the implication (L:244-246; S:1403-1404).
  * A constraint at one value of the gap is a data point, not a constraint on a family of theories.
* **Loose constants.**
  * |ĝ₃|Λ⁴ ≤ 120 is about 10× the transplanted flat optimum (L:320-323; S Sec. S5.H), which CHLPS already state for CFT₃ "without a
    derivation".
  * The R⁴ bound is 17.24 ln(2Λ/3) ≈ 112 minus 19.3. It is set by the 4d infrared logarithm, not by the gap physics.
  * The one theory checked sits at 2.6 %.
* **c₁ is uncomputed.**
  * v8 derives that c₁ is universal (S row L3, [D]; DSKETCH_UP, unreviewed). The value is still open; c₁ = 0 rests on F₂ [A].
  * Being universal, it is fixed by one tree-level spinning Witten diagram with a W³ vertex. A referee will ask why the headline is
    advertised before that calculation.

**Fix.** Positivity over a range of gaps, or a proof that (Pos) at Δ_g implies (Pos) at Δ′_g > Δ_g; and compute c₁. Otherwise this is a
long-paper result.

**Fatal for PRL?** Yes, together with M2.

### M4. The large-spin framing in the Introduction and Fig. 2 understates the premises

**Where:** L:75-77, :272-288; Fig. 2 and its caption (L:290-301); S:1224-1294 (Theorem W), S:1296-1328 (Theorem L′), Table S9.

**Claim.** "for states of large spin we reduce it to two conditional theorems that rest on one technical assumption, a uniform decay
rate of heavy-state contributions beyond the Regge regime" (L:75-77).

**Problems.**
1. **Theorem W does not rest on one assumption.** Its premises (S:1229-1236; Table S9 reading, S:1277-1294) are:
   * AG_U beyond all orders [A];
   * EDGE-R for the graviton [M]. No graviton edge integral was computed at any of the 47 states.
   * P0, which is [NE] below b = 245, i.e. at the sampled gap states with J < 139;
   * P1′, which is [D-sketch] beyond the leading-Regge kernel;
   * the window premise W [M];
   * the boundary factor 2 [D-sketch]. A factor-2 error would raise d/H to ≈ 0.9 at J = 300, 344 (S:1290-1291);
   * the +− transfer convention [M].

   The Large-spin paragraph (L:283-287) lists most of these in words, so the Introduction contradicts the Letter's own body.
2. **"J ≥ 26" is set by a scalar analogue that has almost no margin.**
   * At (26, 1.026Λ) the tolerance on the edge integral is 0.406 S_E,t. The exact *scalar's* integral is 0.285 at the tabulated cutoff
     and 0.405 at the worst cutoff of the lemma's grid (S Table S9; I confirmed both in `code/edger/cache/scalar.txt` and
     `code/agproof/cache/chi.txt`).
   * The comparison also assumes that the graviton's +− edge, in the Ω^U-transferred units, behaves like the scalar's in its own
     units. That is a model convention (S:1256-1261).
   * The Letter's "a scalar analogue meets it, narrowly at the lowest spins" is fair. "Positivity on every sampled near-gap state of
     spin J ≥ 26" (L:277) should carry the condition in the same sentence.
3. **Theorem L′ is "every state" in the Letter, but its constants are sampled.**
   * The Letter (L:277-279) says "on every state of spin below 10⁶ whose impact parameter exceeds about a third of the AdS radius".
   * The certified part, the leading-Regge floor, does hold for every state.
   * The remainder ε_R ≤ 0.133 was evaluated at 16 sampled states, J = 223 to 10⁶ and m = m_gap to 5Λ, and is "smooth in J" (S:1309-1316).
     Yet the J\*(m) formula is applied to every m ≥ Λ: at J ≤ 10⁶ that is m up to ≈ 5000Λ, far outside the sampled masses.
   * P_L′(i) is [M] above J = 344, and EDGE is [M].
4. **Fig. 2 draws Theorem W as a solid region.** It fills the whole band 26 ≤ J ≤ 344, m ≤ 1.3Λ, in the same colour as Theorem L′. The
   Supplement says "Theorem W is a statement about the sampled states; the region between them is covered by smooth trends [NE]"
   (S:1293-1294). The caption's "cover the sampled near-gap states" does not undo the visual impression.
5. **"which we prove order by order in the distance from the Regge limit" (L:272-273) is stronger than the grade.** The identity is [P]
   through O(r⁸) (O(r⁶) with b symbolic), [DV] to O(r⁴⁰), and [D] at all orders (S:1066-1072). "Order by order" reads as a proof at
   every order.

**Fix.**
* Introduction: "for states of large spin we reduce it to conditional theorems whose main premises are an assumption on the decay of
  heavy-state contributions beyond the Regge regime and, near the gap, one uncomputed edge integral of the graviton."
* L:277: "on every sampled near-gap state with J ≥ 26, provided the graviton's edge integral is no larger than its scalar analogue".
* L:278: "for every state … with remainders computed at sampled spins and masses".
* L:272: "which we prove through eighth order".
* Fig. 2: draw Theorem W as markers at the 47 sampled states, with the band between them hatched as "smooth trends". Use a lighter
  shade, or a label, for W than for L′, since W alone needs an uncomputed integral.

**Fatal?** No; wording and figure. But it is the paragraph a large-spin expert will read most closely.

### M5. Light single traces: the A1′ headline rests on a sketched step, and the Supplement's account of it is internally inconsistent

**Where:** L:102-104, :199-205, :209-211, :240-242; S Table S2 rows LS1–LS6, S:741-794 (Sec. S5.F), S Table S7 row Pos_light;
DSKETCH_UP §3.2–§3.3.

**Claim.** The Result is stated under assumption (i), which allows light spin-0 and spin-2 single traces. The Letter says that such
states "enter the k = 2, 3 sum rules like heavy states, with a contribution fixed by their Regge behaviour up to a small AdS mass shift;
we sketch the argument". It then says: "We prove that Ω^EE is positive on them at every mass in the flat limit, so the cubic bound is
unaffected" (L:199-203).

**Problems.**
1. **Light exchanges violate the convergence condition at tree level.** A light scalar coupled through φW² has Regge spin 3, and a light
   massive spin 2 has Regge spin 4; both exceed k = 2. Their tree-level light value is defined only through the arc at infinity, i.e.
   CHMRS's light-side formula carried to spinning functionals: input (E′), [D-sketch]. This is exactly the term-wise use of a
   tree-level object above the convergence border that the Letter correctly disclaims for the heavy side (L:184-186).
   * The order-of-limits argument (S:614-626) shows that the light side is O(1/c_T) for the identity, T and double traces at double
     zeros.
   * It does not show this for a light single trace of Regge spin n > k. Its block's per-ν action is, in CHMRS's representation,
     exponentially large in ν (S:786-788; DSKETCH_UP §3.2).
2. **The Supplement says two incompatible things.**
   * S:758: "With the [TT] double zeros it equals the block action."
   * S:786-791: "The exact per-ν action of a single block grows like e^{πν/2} above its mass … the bare block action is exponentially
     far from these numbers".

   DSKETCH_UP §3.3 reconciles them: the growth is a property of CHMRS's *subtracted scalar* kernel, and for the project's unsubtracted
   spinning k = 2, 3 functionals the block action would equal the Regge polynomial *if* (E′) holds. The Supplement does not say this.
   As written, a reader concludes that the light side is either a polynomial or e^{(π/2)Λ} ≈ e^{1570}.
3. **The "so" at L:202 is a non sequitur.** Flat-limit positivity on light states does not make the AdS cubic bound "unaffected". That
   needs (E′) and the size of the AdS shift δ(Δ_L), which is [M] ("natural size ≈ 21").
4. **The headline therefore inherits a sketched step.** The net grade of item (1) with light single traces is [D-sketch] (S:311-314).
   Without them it is [D].

**Fix.** Choose one of:
* **(a)** State the headline under "no light single trace couples to TT" (where implication (1) is derived, and which ABJM-IIA satisfies
  at the orders computed), and present A1′ as an extension that rests on (E′).
* **(b)** Prove (E′) for the spinning unsubtracted functionals.

In either case:
* rewrite S:752-760 and S:786-794 so that the representation dependence is explicit;
* replace "so the cubic bound is unaffected" by "so, given the light-side formula we sketch, they enter with the favourable sign and the
  cubic bound also bounds their couplings".

**Fatal?** No, if (a) is adopted. A CLS-D or CHMRS referee will raise point 1.

### M6. Wording stronger than the Supplement's grades in the positivity sentences

**Where:** L:74-75, :79-80, :155-156, :180, :242-243, :256-259.

**Problems.**
* **L:74-75, "In AdS we establish it numerically, with error control, for the low-lying heavy states."**
  * "Establish" is too strong for an error budget that the Letter itself calls "an estimate, not a bound" (L:80).
  * "The low-lying heavy states" suggests all of them. In fact it is 14 spins at partial mass ranges (M2); J = 9–11, 13–15, 17–20, 22,
    24, 25 and the 1.1–1.3Λ band are model-level.
* **L:256-259, J = 29.** It is listed under "with error control". By the definition at L:79-80 ("every modelled term carries a
  threefold safety factor") it does not qualify: sign-free leverage 0.94 (S:888, Table S5; SMALLNU_WINDOW, EXT_FOLD boxes). The caveat
  that follows ("depends on the measured sign of one low-frequency contribution") is honest, but J = 29 should be counted separately.
  The count then becomes "twelve further spins …, and J = 29 with a signed window".
* **L:155-156, "every heavy state contributes exactly its flat value".** For k = 2 the de-shift identity holds to ε_W ≈ 1.1·10⁻⁶ [DV],
  not exactly. Only C4 has a certified bound.
* **L:180, "We show that the OPE sum may be exchanged".** "Show" is not in the Letter's legend (L:78-82). Use "we derive" (H_K is [D]).
* **L:242-243, "Independent code reproduces every constant, and a third check rebuilt the two most delicate analytic steps".**
  * The reproductions are internal re-implementations, not external work (S:101-103: "none of these is externally refereed").
  * The two steps rebuilt (the EE-kernel reduction and Lemma Y) are not obviously the most delicate. H_K's corner count, the c₁ ASD
    closure, the R⁴ calibration and (E′) are at least as delicate, and none of them was rebuilt (S:1457-1466, "Not reproduced").
  * Write: "Separate code reproduces every constant of the bounds, and a separate check rebuilt the EE-kernel reduction and the
    counting lemma from the written statements alone".

**Fatal?** No.

### M7. The record behind v7–v8 is unreviewed, and the Supplement is an internal ledger

**Where:** S:39-118 (Sec. S1), S:1457-1466, S:1505-1556; the UT update boxes AG_PROOF, DSKETCH_UP, EDGE_R and NU_TAIL_CHECK (all
"unreviewed").

**Problem.**
* Every analysis that changed between v5 and v8 is marked unreviewed in the record. That covers P2_CONTRACT, P4_WINDOW, LRINF_TIGHT
  (its certification excepted), AG_PROOF, EDGE_R, DSKETCH_UP, NU_TAIL_CHECK and REGGE_TREE. INDEP_CHECK2 reproduced only the v5 layer.
* **The Supplement reads as a project log, not as supplemental material:**
  * version history: "What moved here from the Letter (v4)", "New in v5", "In v6 …", "In v7 …", "In v8 …";
  * a "Record corrections" section;
  * grade composites such as "[D-sketch + M]", "[NE + M] given H_RS" and "[P given P0, P1′, W, AG_U, EDGE-R] + [DV]";
  * 25 pages of tables.
* A referee needs the derivations, not the history of their grades.

**Fix.**
* Before submission, strip the version history and the provenance notes.
* Rewrite Sec. S5 (H_K, c₁, the arc, light single traces) and Sec. S7 (Theorems W, L′) as derivations.
* Obtain one genuinely external read of Sec. S7 and of the c₁ closure.

**Fatal?** Not by itself.

---

## 3. Minor issues

**Length (APS counting).**
* **m1. Compliant; the authors' count is confirmed.**
  * *My count.* I detexed from `\maketitle` to `\bibliography`, with figures and displays removed:
    * prose 2534 words, plus 196 inline formulas counted as one word each, giving 2730;
    * captions 205 words (98 + 107);
    * three display rows × 16 = 48;
    * Fig. 1: 243 × 259.2 pt, aspect 0.9375, so 150/0.9375 + 20 = 180;
    * Fig. 2: 243 × 234 pt, aspect 1.0385, so 164.
  * *Totals.* **≈ 3327 APS words** with each inline formula as one word. If long formulas count as about ⅓ word per token, it is
    ≈ 3560. The authors' ≈ 3333 / 3385 is consistent with the first figure. The limit is 3750, so the headroom is ≈ 190–420 words.
  * *Abstract.* ≈ 576–588 rendered characters, against a limit of 600.
  * *Rendered.* 5 pages; the text ends on page 4; 0 errors; no undefined references.
  * *Note.* Most of the headroom M1 and M4 need would come from shortening "Large spin" (410 words, 12 % of the text) to about 200.

**Prose (meant to read as plain human prose with few numbers).**
* **m2. Mostly achieved.** No grade tags remain, and numbers are confined to the Result, Setting and ABJM paragraphs. Remaining problems:
  * *Long sentences.* Mean ≈ 27 words; many exceed 50 (L:72-77, L:251-263, L:317-323). "Is the functional positive?" averages 34
    words per sentence.
  * *Unexplained jargon:* "celestial frequency ν" (L:188), "de-shifting … a generalised Abel equation" (L:154-156), "Regge 'arc'"
    (L:158), "light-cone edge" (L:283-284), and "splits into two waves, analytic in the frequency, that decay at a uniform rate"
    (L:280-281). A general reader cannot use the last one.
  * *Mannered contrasts.* "a number rather than a scaling law" (L:92), "a consistency check rather than a test of positivity"
    (L:95-96) and "a consistency check at the one-digit level, not evidence for positivity" (L:313-314) appear in quick succession.
    One is enough.
  * *Rhetorical question in the abstract.* "but how close must it be?" is answered at one gap, at tree level and conditionally. Drop it,
    or answer it plainly.

**Title and abstract.**
* **m3.** See M1. In addition, "its non-Einstein three-point couplings" has an unclear antecedent ("its" refers to the CFT).

**Wording against the Supplement's grades** (beyond M4–M6).
* **m4. L:200, "up to a small AdS mass shift".** The size of the shift is [M] (S:306, :776-780, "natural" shift ≈ 21). Say "an AdS mass shift we
  have not computed".
* **m5. L:219, "x ∝ |c₁|/Λ²".** Either give x = 19941|c₁|/(Λ²√(16KL_A)) (v6 did), or drop x and keep only the 1 % thresholds (L:236-237).
* **m6. L:146-148, "curvature moves each heavy state's contribution away from its flat value, by amounts that are not small near the
  gap or at large impact parameter".** Good. It would help to say there that the binding state's AdS correction is −2·10⁻⁵ (report 3
  §5 #6), i.e. the binding state is essentially flat.

**Figures.**
* **m7. Fig. 1.**
  * The transplanted flat optimum uses the same marker code (filled = "claimed", open = "one hypothesis open") as the paper's own
    bounds. Use plain markers for the comparison curve.
  * "(2b): light scalars allowed" should add "needs both functionals".
* **m8. Fig. 2** (see M4.4).
  * The lower-panel bars for J = 27, 29 and 31 overlap the Theorem W region and are hard to read.
  * J = 29's bar should be drawn differently (signed window only).

**Supplement details.**
* **m9. The functional's definition beyond Λ (S:150-163, :186-193).**
  * Adopting the continuum k = 2 weights as the definition is right.
  * But S:192-193 says the truncated-grid reading and the continuum weight "both give the same |ĝ₃|² light coefficient and hence the
    same bound". The truncated grid has *compact ν-support*. By the Supplement's own no-go [P] (S:560-562), its leading-Regge image grows
    like e^{β/2}, so H_K fails under J₀ ≤ 1. The truncated grid is therefore a numerical tabulation on ν ≤ Λ + 10, not an admissible
    functional. Say so, and keep "definition" for the continuum weight alone (as the Table S1 caption already does).
  * The continuum weight does not "decay like e^{−(π/2)(ν−Λ)} beyond Λ" (S:153) uniformly. I evaluated the JSON formula in 60-digit
    arithmetic:
    * χ_re falls to 2.1·10⁻²² at Λ + 28, rises to 2.4·10⁻²⁰ at Λ + 32 (the last Gaussian centres sit at Λ + 30.2 with σ = 0.9; |X| ≈
      10⁻¹²), and then decays super-Gaussianly (10⁻⁴¹ at Λ + 40).
    * NU_TAIL_CHECK accounts for this ("last Gaussians ≤ 1.3·10⁻¹²"), and it is harmless.
    * The sentence should mention the bump, since it is 100× the grid floor that motivated the redefinition.
  * The leading-Regge floor is certified for the exact minimal-norm weight. For the defined (continuum) weight it holds up to its
    ε_W = 1.1·10⁻⁶, which is [DV] (S:1131-1132). "Certified" in L:255-256 should be read with that qualifier. This is trivial, but
    it is the only place where "certified" applies to an object other than the one defined.
* **m10. Theorem L′ table (S Table S10).** The mass column stops at 5Λ, while the theorem is stated for every m with J ≤ 10⁶ (M4.3).
  Add a sentence on why ε_R does not grow with m at fixed β.
* **m11. Order of limits (S:614-626).**
  * The argument is fine, but the uniformity in spin of the light double-trace remainder is "[D-sketch] inside A2".
  * Say once, in the Letter's Setting, that A2's "uniformly in spin" is what makes the infinite sum of light double traces o(1/c_T).
  * Cite CHMRS's equivalent statement.
* **m12. Strip internal comments from the source before arXiv.** The source contains notes such as `% Front matter supplied by the
  authors (v5)`, internal file names (L:1-8, 249-250, 265-268) and `% v6: …` markers.

**References.** The list is complete and correctly formatted; JHEP years render as "J. High Energy Phys. 2009 (10), 079".
* **m13. Suggested additions:**
  * Kravchuk–Simmons-Duffin, *Light-ray operators in CFT* (arXiv:1805.00098), for the light transform (L:126-127);
  * Penedones–Silva–Zhiboedov, *Nonperturbative Mellin amplitudes* (arXiv:1912.11100), for Regge boundedness in Mellin space (S
    Sec. S5);
  * Hartman–Jain–Kundu (arXiv:1509.00014), alongside Refs. [3–6] on Regge causality.

  BCJP [33] is cited only in the Discussion; it could replace the second ABJM citation at L:311.

**Editorial.**
* **m14. AI-tool disclosure.** The internal record shows extensive automated-agent work. If that applies to the submission, APS policy
  requires the use of AI tools to be disclosed (for example in an acknowledgment). The Letter has no acknowledgments section.

---

## 4. Report 3: status of each item in the v8 text

| report-3 item | status | evidence (v8) |
|---|---|---|
| **M1.1** near-gap band 345 ≤ J < 857 in no theorem | **addressed** (by new work) | Lemma LR∞ certified to β ≥ 0.37, so Theorem L′ starts at J\* = 223 at the gap; the band lies inside L′ (S:1297-1303, :1380-1381). I confirmed J\* = 223 (m_gap(223) = √(998·1445) = 1200.88 in AdS units; 0.37211 m = 446.86 ⇒ J ≥ 222.93; J = 222 fails). "Heavier states at intermediate impact parameter" are named as uncovered (L:287-288) |
| **M1.2** Theorem W sampled; P0 omitted | **addressed in text, not in the figure** | "every sampled near-gap state" (L:277); "leading margin is certified only at its larger spins" (L:286-287); P0 graded (S:1230-1232). Fig. 2 still shades W as a solid region (M4.4) |
| **M1.3** premise list incomplete in the Letter | **partly** | The Large-spin paragraph names AG_U, EDGE-R and "other magnitudes … in models" (L:279-287). The Introduction says "rest on one technical assumption" (L:75-77), which contradicts it (M4.1). P1′ [D-sketch] and the boundary factor 2 are not mentioned in the Letter |
| **M1.4** Λ_amp from one state | **addressed** (by new work) | O(1/m²) corrections computed at all 47 states, cost ≤ 0.41H (S Sec. S7.C; P2_CONTRACT) |
| **M1.5** Theorem L tolerance for a scalar-shaped vector | **addressed** | No amplitude is quoted in the Letter; the Supplement keeps the shape and light-cone-normalisation caveats (S:1316-1321) |
| **M1.6** uncertainty not reduced in kind | **partly** | O(1/m²) is now computed [DV]. The beyond-all-orders step is the assumption AG_U [A], and EDGE-R for the graviton is uncomputed [M] |
| **M1 fix (Intro wording)** | **not addressed; regressed after v6** | v6 had "on parts of the large-spin region … conditional statements whose premises are measured on one state or modelled". v8: "reduce it to two conditional theorems that rest on one technical assumption" (L:75-77) |
| **M2** order of limits / J₀ for a tree-level result | **addressed** | L:129-133, :180-186; S Table S2 rows A2, R1, H_K; S:614-635 ("Order of limits", "Term-wise use (not made)"). Rests on REGGE_TREE (unreviewed), but the logic is CHMRS's and is correct. The light-single-trace analogue is not covered (M5) |
| **M3** thin [NE±] states | **partly** | The Letter no longer quotes ±0.018; the Supplement says "Monte Carlo error only" (S:888). "Window measured on ν ≥ 10–20" is in the Supplement only. J = 29 is still counted under "with error control", against the Letter's own definition (M6). "Window live again for J ≳ 32" appears only in the Supplement |
| **M4(a)** "proved in the flat limit" | **regressed** | Fixed in v6; the v8 abstract says "certified in flat space" of both functionals (L:35-36). The body is correct (L:74, :253-254) |
| **M4(b)** leading order in 1/c_T | **addressed** | Abstract (L:33), L:85, :93, :210 |
| **M4(c)** "10⁻⁸" mostly Λ⁻⁴ scaling | **addressed in the body, not in the abstract** | L:89-92 is good; the abstract still leads with 2.89·10⁻⁸ and omits |ĝ₃|Λ⁴ ≤ 120 |
| **M4(d)** de Sitter sentence | **addressed** | Reduced to one clause stating that a dS positivity principle is lacking (L:333-335). Acceptable |
| **M5** importance (one gap, loose constants, c₁) | **not addressed** | Unchanged: one gap; 10× in |ĝ₃|; IR-log dominated; c₁ universal but uncomputed. Next steps are listed (L:331-333) |
| **M6** unreviewed layer | **partly** | INDEP_CHECK2 reproduced the v5 layer (identity, monopole, LR∞, W algebra; S Table S12). The v7–v8 layers (P2_CONTRACT, P4_WINDOW, AG_PROOF, EDGE_R, DSKETCH_UP, NU_TAIL_CHECK) are unreviewed (S:1462-1465) |
| **M7** readability; positivity ledger in the main text | **mostly addressed** | The positivity paragraph is 177 words with no leverage values; Fig. 2 coverage map added. "Large spin" is still 410 words of specialist caveats (m1, m2) |
| m1 length | **addressed** | ≈ 3327 APS words (§3 m1) |
| m2 B_eff(Λ) | **addressed** | "B_eff(Λ) ≈ −19.3 varies slowly" (L:229); values in S:342 and Table S4. I confirmed −19.2646 / −19.2197 |
| m3 overlapping spin lists | **addressed** | The Letter has no lists; S Table S11 is disjoint |
| m4 item (ii) mass range | **moot** | Item no longer in the Letter |
| m5 "proved … floating-point constants" | **addressed** (by new work) | LR∞ is now [P + certified]; I recomputed ∫\|ψ₂\| = 10.3860299071292 and ∫\|ψ₃\| = 12.8356222839601 |
| m6 K grade | **addressed** | K = 1 + O(Λ⁻²) in the Letter; grades in S row L2 |
| m7 "So" + [D] for the identity | **partly** | Now "prove order by order … derive at all orders" (L:272-273), stronger than [P] through O(r⁸) (M4.5) |
| m8 Theorem L "any amplitude" | **addressed** | Removed |
| m9 "only at J ≤ 100" | **moot** | Λ_amp is now computed |
| m10 Theorem W J\* constants | **addressed** | J ≥ 26 is tied to EDGE-R in words (L:277, :283-285) |
| m11 c₁ thresholds by item | **addressed** | "about 60 (cubic) or 1300 (quartic)" (L:236-237); I confirmed 59.6 / 1306 |
| m12 x formula | **regressed to "x ∝"** | L:219 (m5 here) |
| m13 composite grades undefined | **addressed** | No tags in the Letter; plain-word legend L:78-82, mapped in S Sec. S1 |
| m14 Fig. 1 | **addressed** | c₁ = 0, K = 1 in the annotation; open diamonds at 300; "nothing claimed" band; implication-only lines |
| m15 two scope points (H_bulk, 1.1–1.3Λ band) | **partly** | Only in the Supplement (S Table S5 large-mass row; Table S7 H_cont); the Letter's "Everywhere else … models or interpolated" covers them implicitly |
| m16 DKKPS | **addressed** | L:107-109 |
| m17 weighted (2a) | **addressed** | L:238-240 |
| m18 ABJM orders | **addressed** | L:309-310 |
| m19 HKT | **addressed** | L:52-53 |
| m20 identity residual | **addressed** | S Sec. S7.A |
| m21 refined quadrature | **addressed** | S Sec. S7.A |
| m22 P0 grade | **addressed** | S:1230-1232 |
| §7.1 settle M2 | **done** | as above |
| §7.2 rewrite main text | **largely done**, but undone at the top | New title and abstract (M1) |
| §7.3 compute spinning O(1/m²) | **done** | P2_CONTRACT; the higher orders are now AG_U [A] + EDGE-R [M] |
| §7.4 close the near-gap band | **partly** | J = 7, 27 added; J = 23, 31 at higher masses; J = 9–11, 13–15, 17–20, 22, 24, 25 and the 1.1–1.3Λ band remain model-level |
| §7.5 compute c₁ | **not done** | Universality derived [D]; value [I] |
| §7.6 independent review | **partly** | v5 layer reproduced; v7–v8 not |
| §7.7 range of gaps | **not done** | — |
| §7.8 re-optimise Ω^Q | **not done** | — |

**Reports 1 and 2.** Report 3's §4 table stands for them, with three changes in v8:
* **R1-M12 (title) has regressed** past its original form (M1).
* **R1-M13 (ABJM in the abstract) has partly regressed.** "ABJM theory lies inside this bound" is back in the abstract. ABJM-IIA is now
  inside A1′, so this is a placement issue, not a scope error.
* **R2-M4 and R3-M4(a) (grade slips in the abstract) have regressed** ("certified in flat space").

Items fixed since report 3 and still fixed: R2-M2 (no lower end advertised; L:243-245), R2-M8 (no corollary; L:245-246), R2-M5
(positioning against CHLPS/CHLPSpw; L:57-67, :320-327).

---

## 5. Spot checks of load-bearing numbers

| # | number (Letter / Supplement) | traced to | my check | verdict |
|---|---|---|---|---|
| 1 | τ₀ = 2.89·10⁻⁸ (1000), 3.13·10⁻⁶ (300); \|ĝ₃\|Λ⁴ ≤ 120 / 104 (L:220, abstract) | `data/nf02_psi.json` | `chk_const.py`, 40-digit mpmath from the exact rationals: K₀ = 1.00000000003; L_A = 896.22432 / 670.91284; √(16L_A) = 119.748 / 103.608; τ₀ = 2.8912626·10⁻⁸ / 3.132055·10⁻⁶ | ✓ |
| 2 | 2L_AdS = 92.81 / 72.04; B_eff ≈ −19.3; (2b) 164.5; c₁² coefficients 5.44·10⁻⁷ / 6.81·10⁻⁵; 1 % thresholds 60 / 1300 (L:225-237) | `data/q3star_factored.json` | 92.81351 / 72.040593; B_eff = −19.2646 / −19.2197; 2L_AdS + 0.08L_A = 164.51146 / 125.71362; 57600 I_Q²/(κΛ⁴) = 5.4428·10⁻⁷ / 6.8145·10⁻⁵ (I_Q = 1.0160366, κ = 0.1099083); cubic threshold 59.6 / 4.61, quartic 1306 / 103 | ✓ |
| 3 | ABJM 2ζ(3) at 2.6 % (L:312) | S Sec. S10 | 2ζ(3)/92.81351 = 2.5903 %; 3.337 % at 300 | ✓ |
| 4 | 0.46Λ scalar threshold (L:204, :222) | S row LS5 | √(0.430430/2.032073) = 0.460237 | ✓ |
| 5 | **Theorem L′: J\* = 223 at the twist edge**; J\*(m) = 279 / 372 / 929 at 1.5 / 2 / 5Λ (S:1297-1303, Table S10) | LRINF_TIGHT; formula | b = (2J+1)Λ/m ≥ 0.37211Λ ⇔ β ≥ 0.37 (2 sinh(0.185) = 0.37211). m_gap(223) = √(998·1445) = 1200.88 gives J ≥ 222.93; J = 222 gives m = 1200.05 and needs J ≥ 222.78 > 222. Hence 223; 279, 372 and 929 also reproduce | ✓ |
| 6 | **Theorem W: sampled J\* = 26 given EDGE-R**; d_χ ≤ 0.489H; tolerance 0.406 S_E,t at (26, 1.026Λ); exact scalar 0.285 at the tabulated cutoff, grid maximum 0.405; J\* = 38 if ≥ 1.2 is required (L:277, :283-285; S Table S9) | `code/agproof/cache/chi.txt`, `code/edger/cache/scalar.txt` | chi.txt: max d_χ/H = 0.489 (J = 344); minimum transferred tolerance 0.406 at (26, 1.026); next 0.629 (30), 0.886 (34), 1.22 (38). scalar.txt: J = 26 own cutoff 0.285, grid (38 cutoffs) min 0.205, median 0.241, max 0.405. So "J\* = 38 at ≥ 1.2" follows. The graviton column is empty at every state, as stated | ✓; **margin 1.0025 at the worst cutoff, and the graviton is uncomputed** (M4.2) |
| 7 | Lemma LR∞ constants: ∫\|ψ₂\| = 10.386029907129, ∫\|ψ₃\| = 12.835622283960, a₁ = 188.896, a₂ = −1023.546 (S:1103-1118) | `nf02_psi.json` | `absint.py` (polynomial roots in (0,1), then a sum of \|∫\| over sign intervals): 10.3860299071292, 12.8356222839601; ψ₂′(0) = 188.896105, ψ₂″(0) = −1023.5457 | ✓ |
| 8 | **c₁ universality: the Barnes map** T ↦ (A²+ν²)/2A, "polynomial top-Regge coefficient ⇒ polynomial in ν² of the same degree" (S:653-657) | DSKETCH_UP §2.3; `code/dsketch/t1b_barnes.py` | Fresh contour integrals at (A, ν) = (5.3, 2.2) and (8.1, 17.0), 25 digits: ⟨T⟩/⟨1⟩ = 3.10660377358491 and 21.8895061728395, equal to (A²+ν²)/(2A) to all digits shown. This is Barnes' first lemma with a = 0, b = j. Also x/\|c₁\| = 1920 I_D/(Λ²√(16KL_A)) with I_D = 10.3859465 gives 1.67026·10⁻⁴ / 2.16006·10⁻³ | ✓ for the scalar CHMRS representation. The spinning Mellin structure it is applied to is [Q] |
| 9 | **Truncated-grid functional:** rows kept/dropped 15425/959 (300) and 25212/7556 (1000); max \|X\| dropped ≤ 5.7·10⁻¹²; continuum agreement on (Λ, Λ+10] ≤ 2.6·10⁻¹⁰; SHA-256 of the four copies (S Table S1; `data/README.md` §7) | `data/truncate_grids.py` | Reran it on scratch copies: identical row counts; dropped max \|X\| 6.76·10⁻¹⁴ / 1.31·10⁻¹⁴ / 5.67·10⁻¹² / 1.22·10⁻¹³; continuum agreement 1.3·10⁻¹⁰ / 8.1·10⁻¹¹ / 2.6·10⁻¹⁰ / 4.0·10⁻¹²; **all four SHA-256 identical to `TRUNC_SHA256SUMS`**; source grid prefix 98d42b28… matches MANIFEST | ✓. Definitional caveats in m9 |
| 10 | Continuum k = 2 weight tail beyond Λ (S:150-155) | `data/k2/k2_continuum_D_1000.json` | `tail.py`, 60 digits, χ_re at Λ+{0, 10, 20, 28, 30, 32, 35, 40, 60}: 2.1·10⁻¹⁸, 1.1·10⁻²¹, 5.5·10⁻²², 2.1·10⁻²², 2.5·10⁻²¹, 2.4·10⁻²⁰, 9.3·10⁻²³, 2.2·10⁻⁴¹, 1.5·10⁻²⁵³ | decays, but with a ×100 bump at Λ+30–33 from the last Gaussians (m9). Harmless, consistent with NU_TAIL_CHECK |
| 11 | 8πG = 256/c_T (L:101) | standard | C_T^E = 3L²/(π³G) in d = 3; C_T^{free scalar} = 3/(32π²); ratio 32L²/(πG) | ✓ |
| 12 | τ = 240\|ĝ₃\| (L:116) | DICTIONARY_T4 | 240 = 2(d−1)d(d+1)(d+2) at d = 3. The t₄ half is [DV]; the α_T half is [D] with an open sign. The collider disk τ ≤ 4 ⇔ \|ĝ₃\| ≤ 1/60 | ✓ |

---

## 6. Physics checks

* **Normalisations.** Both checks pass (§5 #11–12).
  * *8πG = 256/c_T.* Correct with c_T = 1 for a free real scalar.
  * *τ = √(t₄²+α_T²) = 240|ĝ₃|.* Correct for t₄, with an independent literature check [DV]. For α_T it is derived [D], and only the
    modulus enters, so the open sign is harmless.
  * *Collider comparison.* Correct (Chowdhury–David–Prakash disk, τ ≤ 4).
* **Flavour currents versus c_T.** There is no confusion.
  * c_T enters only through 8πG.
  * Spin-1 and spin-3 single traces, R- and flavour currents included, are excluded by the absence of ⟨TTO_J⟩ structures at J = 1, 3
    for identical conserved T in d = 3 (DKKPS, cited at L:108).
  * The Supplement states explicitly that this concerns the ⟨TTJ⟩ structure, not any current two-point coefficient (S:738-740).
  * Odd J ≥ 5 couple only through a parity-odd structure. That is why the binding state (J = 5) is a 1×1 parity-odd block.
* **Order of limits.** Correct.
  * The sum rule is used for the exact correlator at finite c_T, where J₀ ≤ 1 (KKSDZ) and H_K holds with margins ½ and 3/2.
  * (Pos) is a statement about (Δ, J, parity) and holds at every c_T.
  * The light side is O(1/c_T) because [TT] sits at double zeros for k = 2, 3, and through the tree-level double-trace data at k = 4, 5.
  * The bound is the leading coefficient. This is CHMRS's logic (their Sec. 1.1.3), and the "term-wise use (not made)" paragraph is a
    useful clarification.
  * Two loose ends:
    * the uniformity in spin of the light double-trace remainder is only [D-sketch] inside A2 (m11);
    * the argument does not cover light single traces of Regge spin n > k (M5).
* **A1′ versus ABJM.** Sound at the level claimed.
  * In the IIA 't Hooft regime, KK modes on CP³ have J ≤ 2, and the first massive string level (J ≤ 4, twist ≈ ML − 4) sets the spin-4
    gap. So A1′ holds by the spin count.
  * A2 holds at g_s → 0.
  * The scalar condition of (2a) relies on CPY and BCP's statement that only the stress-tensor multiplet is exchanged at tree level
    through R⁴ order. It is a [D-sketch], and the Letter says so (L:308-311). For N = 6, the ⟨TTTT⟩ statement goes beyond what the
    bottom-component four-point function fixes (S:1484-1487); this is also flagged.
  * The identification Λ ↔ ML shifts the placement by −11 % / −3 % relative at 300 / 1000 (S:1494-1497), so "one digit" is right.
* **The light-scalar condition.**
  * Necessary: Ω^Q_up is negative on flat light scalars with m < 0.46024Λ (exact polynomial −0.4304 + 2.0321m²).
  * The weighted form (L:238-240) correctly turns it into a bound when such scalars couple.
  * (2b) removes the condition at the price of positivity of Ω^EE and the ε = 0.040 repair. That repair is [D-sketch], and ε\* = 0.038
    is expected to fail for Δ_L ≲ 6 [M], so 0.040 is necessary.
  * In AdS the threshold applies to μ² = (Δ_L−2)(Δ_L−1) + δ(Δ_L), with δ uncomputed [M]. At Λ² ≈ 10⁶ this is immaterial.
* **c₁.** The structural claims check out, but the value is missing.
  * The exact c₁ forms are correct: I rederived Eq. (S16) and the quartic supremum t₄\* = 57600 I_Q c₁/(κΛ⁶).
  * Universality is now [D]. The Barnes step reproduces (§5 #8). The ASD closure is [D] with an explicit solution [DV], but unreviewed.
  * The value is [I], and c₁ = 0 needs F₂ [A].
  * A referee will insist on the one W³-vertex Witten diagram: c₁ is universal, so one theory fixes it.
* **H_EFT(b).** g̃₄ = g₄ + O(Λ⁻²) needs light double-trace data to be tree-level EFT. This is the standard CHMRS/CLS-D assumption [A],
  stated at L:195-197.
* **κ > 0.** κ is needed only if τ ≠ 0. The flat value is 0.1099, with O(Λ⁻²) corrections, and the sign is [D-sketch + M]. The sign
  matters: with |ĝ₃|Λ⁴ up to 120, κ|ĝ₃|²Λ⁸ can reach ~10³ ≫ 92.8. The Letter states the condition in (2a) (L:223).
* **The large-spin conditional theorems.** The implications are correct and short:
  * Lemma I;
  * Lemma P4W-χ, i.e. Cauchy with an erfc cutoff, using |erfc(a+ib)| ≤ e^{b²}erfc(a);
  * the sharp 2×2 criterion d ≤ f₊₊f₊₋/(f₊₊+f₊₋).

  Assessment of the premises:
  * **AG_U** asks for a split into two analytic waves with the *uniform* rate 0.97Λβ, and is strictly weaker than (AG′). The record's
    finding that (AG′) is pointwise false at the edge, and the switch to AG_U, are to the authors' credit. But AG_U is still an
    assumption beyond all orders. "It holds to all orders of the saddle-point expansion" (L:282) does not control the edge region, where
    the plane-wave representation is shown *not* to supply it at J ≤ 50 (S:1207-1210).
  * **EDGE-R** for the graviton is uncomputed. The direct method fails above ν/Δ ≈ 0.94, and the thimble contour's k = 3 piece is off
    by 200–1000× for an unknown reason (S:1286-1290). J\* = 26 therefore rests on a scalar analogue in transferred units, with margin
    1.42 at the tabulated cutoff and 1.0025 at the worst grid cutoff.
  * **Theorem L′** is solid in its certified part (the floor). Its remainder is computed only at 16 sampled states up to m = 5Λ (M4.3).
* **The functional beyond Λ.**
  * The v8 fix is correct and costless: define the k = 2 pieces by the continuum weights times θ(Λ+100, 10). Every computed number uses
    ν ≤ Λ ≤ m.
  * The continuum weights decay fast enough to beat the exponential growth of heavy images above their mass. NU_TAIL_CHECK bounds the
    tail at ≤ 4.4·10⁻⁷ (D) and 1.6·10⁻⁴ (Q) of a state's action, [DV] with a scalar proxy and [M] for EE.
  * Two cautions (m9):
    * the truncated grid is not itself an admissible functional (compact support);
    * the continuum weight has a small bump at Λ + 30–33.
* **Light single traces.** This is the one structural point a CHMRS-type referee will press (M5). Their light value at k = 2, 3 is an
  arc at infinity whose spinning form, (E′), is sketched. The Supplement mixes two representations without saying so.

---

## 7. What still stands between this and PRL acceptance (ranked by impact ÷ effort)

1. **Correct the title, abstract and Introduction (M1, M4.1, M6).** *Effort: 2–4 hours, no computation.*
   * Retitle.
   * Lead the abstract with the condition and with |ĝ₃|Λ⁴ ≤ 120 at c₁ = 0.
   * Restrict "certified in flat space" to the cubic functional.
   * Fix "one technical assumption", "establish … the low-lying heavy states", J = 29's classification, "prove order by order" and
     "independent".

   This is necessary for any venue. By itself it does not change the PRL verdict.
2. **Restrict the headline to "no light single traces in T×T" and present A1′ as an extension (M5).** *Effort: half a day.* It makes
   implication (1) [D] and removes the most likely technical objection. Proving (E′) instead is weeks to months.
3. **Compute the graviton edge integral EDGE-R.** *Effort: ≈ 90 CPU-h of listed external jobs* (`ext/jobs_batch3*.txt`). The 20
   Theorem W states with m ≥ 1.13Λ need only the direct method. The gap-mass states also need the n3 k = 3 (C3H) bug fixed, which could
   take days to weeks. This turns Theorem W's only [M] premise into data, and makes "J ≥ 26" a statement about the graviton.
4. **Error-control the rest of the near-gap region at Δ_g = 1000 and the Ω^Q_up binding states.** *Effort: ≈ 1–2 CPU-h per spin per
   mass, by HIGHSPIN and SMALLNU, i.e. ≈ 50–100 CPU-h for J = 9–11, 13–15, 17–20, 22, 24, 25 at 3–4 masses up to 1.3Λ, plus J = 29
   sign-free and J = 23, 31 twist edges.* Fig. 2's lower panel would then be uniformly "error control" below J = 26.
5. **Compute c₁.** *Effort: 1–3 weeks; one spinning AdS₄ tree diagram with a W³ vertex.* It removes the "up to one universal constant"
   clause, and one theory suffices.
6. **External review of the v7–v8 layer and a rewrite of the Supplement as derivations (M7).** *Effort: 1–2 weeks.* The priority targets
   are AG_PROOF (Lemmas F and P4W-χ), P2_CONTRACT, DSKETCH_UP (the c₁ ASD closure) and REGGE_TREE.
7. **Positivity over a range of gaps** (for example, monotonicity of the de-shifted corrections in Λ), so that (Pos) at 1000 implies
   (Pos) for all Δ_g ≥ 1000. *Effort: high, months.* This is the single change that most raises breadth (M3).
8. **A proof of AG_U beyond all orders, or an analytic large-spin argument that does not need it.** *Effort: high, possibly out of reach.*
   Without it, large spin stays conditional.
9. **Re-optimise Ω^Q with light-scalar positivity imposed.** *Effort: weeks.* It removes the 0.46Λ condition and the ε-repair.

**Realistic path.** Items 1–2 now (days), then submit a long paper (PRD or JHEP) with items 3–6 as they land. A PRL becomes plausible
only if item 7 succeeds, or if items 3–5 *and* a credible route on item 8 succeed. Positivity would then be established across the
spectrum at error-control level, at more than one gap, with c₁ known.
