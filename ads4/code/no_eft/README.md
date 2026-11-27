# no_eft: removing the tree-level-bulk-EFT clause from the cubic chain

This directory works out how to remove the tree-level-bulk-EFT clause from the cubic chain. It writes nothing outside
`ads4/code/no_eft/` and commits nothing.

## Ledger (final)

| step | status | file |
|---|---|---|
| 0. read letter.tex (assumptions), supplement.tex (Table inputs, HK order of limits, c₁ chain, LS rows, angular route), sec_fullpos.tex | done | — |
| 1. audit of every A2 use, cubic and quartic | done | NOTE.md §2 |
| 2. reduction: finite-c_T inequality, so "leading order" ⇔ c_T R → 0 | done [P] | NOTE.md §1 |
| 3. spin summation: (ii′_E) detector decoupling ⇒ (b); (ii′) ⇒ (ii′_E); (ii′_E) ⇒ 𝓛_g → 0 (τ_\*) | done [D-sketch via Lemma 5; Lemma 6 is D] | NOTE.md §3 |
| 4. numerics: stripped action of Ω^EE on light states | done [NE] | nl1_light_action.py, nl1b_sanity.py; cache/nl1_*.txt |
| 5. can (ii′) be derived from CFT axioms? failing step, literature | done (no; stays [A]) | NOTE.md §5 |
| 6. side check: light J = 0, 2 single traces block-wise | done [NE/D-sketch] | nl1c_light_st.py; NOTE.md §6 |
| 7. proposed edits (exact old/new; OLD strings verified verbatim against the files) | done | PROPOSED_EDITS.md |

## Files

* **NOTE.md** — the self-contained note: audit, lemmas with grades, numerics, what remains assumed.
* **PROPOSED_EDITS.md** — replacement wording, edits E1–E15, NOT applied.
* **nl1_light_action.py** — stripped light action ĥ_b(Δ, J), in ball arithmetic.
  * Run: `python3 nl1_light_action.py 1000` (about 150 s, mostly exact Q_k), or `NL1_DG=300 python3 nl1_light_action.py 300`.
  * Read-only imports from `../fullpos` (fp_common, fp_poly) and its cached moment balls.
* **nl1b_sanity.py** — runs three checks:
  1. heavy states agree with fp_poly.H_ang;
  2. light J = 0, 2;
  3. continuity in twist.
* **nl1c_light_st.py** — the light single-trace J = 0, 2 scan.
* **cache/** — outputs.

## Key numbers [NE]

* **Size of the light action.**
  * Plateau: |ĥ|/T_Δ ≈ 10^{685.4} at Δ_g = 1000 and ≈ 10^{208.7} at Δ_g = 300. This is about e^{πΛ/2}·10^{4–5}.
  * On the light leading-twist double twists (twists 2–4) |ĥ| reaches about 10^{691} by J = 1000.
  * The sign is not definite. Near the gap (twist 998) the values join the heavy ones, |ĥ| ≈ 10².
* **Consequence.** The term-wise light remainder is about 10^{686}/c_T². So "leading order in 1/c_T" is numerically meaningful only
  for c_T ≳ 10^{690}. This holds with or without a bulk EFT.
