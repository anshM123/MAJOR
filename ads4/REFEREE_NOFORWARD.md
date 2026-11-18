# Referee report: `NOFORWARD_FUNCTIONAL.md` and `code/robust/nf_*.py`

Status: independent hostile referee, 18 Nov 2026. Nothing under review was edited. Scratch scripts are in
`ads4/code/referee_noforward/`: `r_contacts.py`, `r_low.py`, `r_scan.py`. Logs are in `referee_noforward/cache/` (git-ignored).

## What was reproduced

| claim | check | result |
|---|---|---|
| `B₂^{(1)}`, `B₃^{(1)}` contact-free | analytic: in (2.26) at `k = 2, 3` the integrand is `(2s−q)·poly` or `poly`, so there are no poles for any polynomial `f`; `r_contacts.py` residues | **correct** for every polynomial `f` in the MHV `f(s,u)`. Contacts in other helicity amplitudes (`g`, `h`) never enter `B^{(1)}`. **Basis-specific:** `B′₂`, `B′₃` (unstripped) do see contacts (`FLAT_UNSTRIPPED.md` item 1) |
| which `B_k` see the isolated data | `r_contacts.py`, own residue code: `P = w^i t⁶` and `P = t^i w⁶` | the order-`i` data enter `B_{2i+4}` and `B_{2i+5}` (and higher), with ratios `+q` vs `−q`. **Not** as stated in §1.2 (see M1) |
| exact rank table | re-ran `nf_contacts.py 7 16 28` | identical: null-space dimensions 64 / 59 / 47 / 35 / 23 / 11 / 3 / 0 |
| NF02 low side | `r_low.py`, independent sympy from `X`, `Y` in the JSON | `k = 88865643959/3233230000`, `A = 3022.3377`, `B = −5327.2096`, `B_AdS = −6552.6620`. `ψ` matches `psi*_exact`. `C₀ < 0`; `m⁴C₂` has no root in `(0,1]` (max −12.3) |
| region C`+−` (`u ∈ [1/4, 3/8]`) | re-ran `nf_cert.py NF02 3/20 Cpm 400 1/4,5/16,3/8` | PASS, log identical to the tracked log |
| M tightest entry | existing log | `4307.8` vs `C_R ≤ 2.267·10⁶`, as stated (thin, but certified) |
| uncertified region | `r_scan.py`, own Gauss–Legendre and Jacobi code, `J ≤ 20000`, 129 `u`-values in `[3/8, 1]` | no violation. Global min `𝔉/N = 0.2000` (`+−`, `J = 5`, `m = M`). In the uncovered region the min is `0.765` (`J = 401`, `u = 1`) [NE] |
| Λ = 100 numbers, ratios 47 / 6.8, 45–60 | arithmetic | correct |

## Verdicts
* **Correctness:** numbers are correct. The obstruction theorem is true (exact rank check plus a corrected induction), but the
  written proof is wrong in detail (M1).
* **Soundness:** the certificates are sound in the stated regions. NF02's heavy positivity is **not** fully certified.
* **Labels:** mostly honest. The headline answer overclaims (M2).
* **Applicability:** this is a flat statement, and the AdS heavy side is untested (stated by the author). The obstruction is
  confined to smeared fixed-`u` `B^{(1)}` sum rules with exact contact cancellation (stated).
* **Novelty:** modest. The forward-free `B₂ + B₃` class is the paper's (3.7)/footnote functional, and the obstruction is
  elementary (the author concedes this). What is new: the necessary conditions `ψ₃(0) = 0` and `(1−p)⁵`, and a positive
  functional certified with a margin.

## Issues
**Fatal:** none.

**M1 (major) — the §1.2 induction is mis-indexed.**
* What the text says:
  * odd `k = 2ℓ+5` uses data of order `≤ ℓ+1`;
  * "the highest odd `k = 2L+5` alone carries the order-`(L+1)` datum";
  * the pair is `(2L+4, 2L+3)`.
* What is true:
  * `[z^ℓ]` uses data of order `≤ ℓ` for both parities;
  * order-2 data are absent from `B₇` and enter `B₈` and `B₉`;
  * the correct pair is `(2ℓ+4, 2ℓ+5)`, with the top `B_K` alone only when `K` is even.
* The repaired argument goes through: the two families give `±` relative signs, which makes the 2×2 system nonsingular, and
  Müntz then applies.

**M2 (major) — the headline overclaims.** The "Answer to the question" says that a robust forward-free functional "exists", and
the commit message says "certified NF02". Neither carries the caveat that `+−`, `J ≥ 400`, `M ≤ m < 1.633M` is [NE] only.

**m1 (minor) — the factor 23 is mis-stated.**
* §2's "≈ 23× the forward-term class (δ = 0: 26.1)" is arithmetically wrong: 860/26.1 ≈ 33.
* §0's "≈ 23" compares against `F_{g3}` (37.84), which is not an optimum.

**m2 (minor) — "every `k ≤ n+4`" is too strong.** It holds only for a generic degree-`n` contact. For example, `su` enters only
`B₄` and `B₅`.

**m3 (minor) — the reproduced footnote functional fails the derived conditions.**
* The paper's (3.7) has `ψ₃(0) = 0.0017 ≠ 0` and only `(1−p)²` decay. It therefore violates both necessary conditions, consistent
  with `EINSTEIN_OBSTRUCTIONS.md` A.4 (it is negative at `m ≈ M`).
* The document "reproduces" this functional without saying so.

**m4 (minor) — the `nf_cert.py` docstring is stale.** It says that all of `1/4 < u ≤ 1` (`+−`) is uncovered and does not list
`Cpm`.

**m5 (minor) — wording on "these families".**
* "Only contact-free forward-free functionals are `∫ψ₂B₂ + ∫ψ₃B₃`" holds in the stripped basis.
* In the `B′` basis the same functionals are combinations of `B′₂ … B′₇`.

**Consistency with M02:** there is no contradiction. The obstruction concerns `L²` weights. The forward terms
`∂_qB₄|₀` and `B₅|₀` see only `g₅` (checked: only the `n = 1` contact survives at `q = 0`).
