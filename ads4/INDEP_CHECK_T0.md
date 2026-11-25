# INDEP_CHECK_T0: independent hostile check of T0_CERT (leading-Regge hybrid action of Ω^EE ≥ c·N_s)

*25 Nov 2026. Independent verifier. No existing file was edited; nothing was committed by me. Scratch code: `ads4/code/t0check/`
(`chk_*.py` are mine; `t0_common.py, t1_curv.py, t2_lrinf.py, t3_flat.py, t0_checks.py, run_all.sh` are byte-identical copies of
`code/t0cert/`, run there so that `t0cert/cache/` was not overwritten). Logs in `t0check/cache/` (its `.gitignore` ignores all but itself).
Compute ≈ 1.5 CPU-h, ≤ 2 processes, `ulimit -v 3000000`, `OMP_NUM_THREADS=1`.
(Note: a periodic "WIP: agent snapshots" commit by another process swept the new `t0check/*.py` files into git; I did not commit.)*

---

## Verdict: **CONFIRMED** (with minor, documentary fixes; no number, constant or [certified] grade changes)

Theorem T0 as stated in T0_CERT §0.1 — hybrid action H ≥ c·N_s, c = 199/1000 for every Δ_g ≥ 300 and c = 1999/10000 for every
Δ_g ≥ 1000, on `++` (even J), `+−` (J ≥ 4), every m ≥ Λ, every b̂ > 0 — **holds**, and the certificate is a genuine proof modulo
the cited inputs it lists (CERT_CLOSE's 𝔉 ≥ N/5, LR∞ + R̂ anchoring, lrt2 M₁ cells, the small-β M₁ bound, Mehler–Dirichlet).

* Reproduction is **bit-identical** to the committed logs (all nine `t{1,2,3}_*.txt` files, modulo timings).
* The algebraic reduction, Lemma 2 (curvature in gap units), Lemma 3, Lemma 4's majorant structure and Lemma 5 were re-derived by hand
  and found correct; Lemma 2 and Lemma 5 were also tested numerically against an **independent** evaluator (₂F₁ conical functions,
  no Mehler–Dirichlet) — remainder |h_eff| ≤ 191 against the certified H_k ≈ 1.2·10⁴.
* Hostile float scan, 3 × 1534 fresh random states: **min H/N_s = 0.19999998 at Δ_g = 1000, 0.20000001 at 5000, 0.19999969 at 300**
  (all at `+−`, J = 4, m = Λ). **Nothing below 0.1999 at Δ_g ≥ 1000; nothing below 0.199 at 300.**

### Issues (severity)

| # | issue | severity | fix |
|---|---|---|---|
| I1 | **OO coverage needs more than lemma (A).** T0 §0.1 says OO is "[certified given (A)]". It also needs (i) the EO entries to cancel in the S⁰ sum (PARITY_ODD_ROUTEA §3; UNCOND_INVENTORY M17 lists the algebraic lemma as **open**), so that OO is a separate 1×1 block; (ii) κ > 0 at every Δ (only per-row values 1.14–2.20 are computed [DV]); (iii) N_s,OO = κ·N_s,++ as the normalisation. UNCOND_INVENTORY §1 explicitly advises *not* to use the Γ_OO = Γ_++ identity in a proof. | minor (scope/wording; the `++`/`+−` sectors are unaffected) | state OO as "[certified given (A), EO cancellation (M17) and κ > 0]" |
| I2 | **Part L at 300 ≤ Δ_g < 1000** uses LR∞ with the R̂ anchor at Λ(300), outside GAP_LARGESPIN's stated range (Δ_g ≥ 1000). I checked that LR∞ (LARGESPIN_L §3, LRINF_TIGHT §1) and the scaling lemma (GAP_LARGESPIN §1.2) are Λ-parametric with no lower bound on Λ, so the use is valid — but T0_CERT should say so. | minor (documentary) | one sentence in §3 Part L |
| I3 | The committed `t0cert/cache/t0_checks.txt` ends with "certified c = 199/1000"; `run_all.sh` as written leaves `T0_C = 1999/10000` exported when `t0_checks.py` runs, so a fresh run prints "1999/10000". The cached file was not produced by the current script (cosmetic; [DV] checks only). | cosmetic | regenerate or set `T0_C` explicitly |
| I4 | H_k is loose by ≈ 60× (actual |h_eff| ≤ 191 vs H₂ ≤ 12124, H₃ ≤ 12273), and T_gap ≈ 29 vs the true worst negative curvature ≈ 8.2 (Λ⁻² units). Harmless (sound, only conservative); it is what keeps Δ_g ≥ 300 from going lower. | none (remark) | — |
| I5 | "Odd J: only `+−`, whose hybrid kernel is the `+−` kernel" is a kinematic/definitional statement (it is how `hyb` is defined), not proved in T0. PARITY_ODD_CLOSE (G)'s [DV]/[NE] concerns a different statement (the engine's antisymmetric block) and is not needed. | none (remark) | — |

No issue affects any certified inequality.

---

## 1. Reproduction

`run_all.sh` (copied verbatim to `t0check/`, run there): **15 min 57 s wall, 15.4 CPU-min** (doc: ≈ 20 CPU-min). All verdicts PASS.
`diff` of every output `.txt` against `t0cert/cache/` (timings stripped): **identical**. Key numbers vs the document:

| run | Part K: cells / min K_g / tail | Part L max ratio; β ≥ 64 | T_gap; T_abs | Part F | doc |
|---|---|---|---|---|---|
| Δ_g ≥ 1000, c = 0.199 | 2840 / [5.0 ± 0.0785] ≥ 4.92 (b̂ = 19.6, (+,−,σ=+,u=1)) / 6.940 | 0.045691 (β = 1.76); 2.03·10⁻⁷ | 28.9879; 2.90751·10⁻⁵ | 584 polys + 672 cells PASS | ✓ |
| Δ_g ≥ 300, c = 0.199 | 2840 / [5 ± 0.415] ≥ 4.585 / 6.492 | 0.19053; 2.27·10⁻⁶ | 29.4355; 3.30357·10⁻⁴ | PASS (shift 3.31·10⁻⁴) | ✓ |
| Δ_g ≥ 1000, c = 0.1999 | 2840 / ≥ 4.92 / 6.928 | 0.045798 | 29.0097; 2.90969·10⁻⁵ | PASS (shift = T_abs from json) | ✓ |

Other doc constants re-derived by hand: S_min = 0.801·2·188.896 − 1.199·50.643 = 241.89 ✓; ε₀(1+c)(2H₂+H₃) = 0.0440 (1000), 0.492 (300) ✓;
584 = 50 (`++`, J = 0…98 even) × 4 + 96 (`+−`, J = 4…99) × 4 ✓; 672 = 84 cells × 4 patterns × 2 kernels ✓.
`t0_checks.py` reproduces (min H/N_s = 0.1999996902 at (300, `+−`, 4, m = Λ)) apart from I3.

## 2. The algebraic reduction

**Lemma 0** [P, re-derived]: with X = B₂ + cv₂, Y = B₃ + cv₃, H − cN_s = X + Y − c·max_{s}(s₂X + s₃Y) = min_{s₂,s₃}[(1−cs₂)X + (1−cs₃)Y],
which is FF_g(c) + (1−cs₂)cv₂ + (1−cs₃)cv₃. ✓ All four sign patterns are checked in every part (`t1` loops (s₂,s₃,σ,u) ∈ {±1}³×{0,1};
`t3` uses `rcert.gens` with δ = c, which enumerates the same four (s₂,s₃)).

**Lemma 1 use**: FF_g = 𝔉 − c(s₂B₂ + s₃B₃) ≥ 𝔉 − cN ≥ (1/5 − c)N ≥ 0 for c ≤ 1/5 ✓ (c = 0.199, 0.1999). It is needed only for b̂ ≥ 16
(where the curvature part is ≥ 0 by K/L); on b̂ ≤ 16 Part F certifies FF_g ≥ T_abs directly.

**Assembly**: b̂ ≤ 16: FF_g ≥ T_abs ≥ −(curvature part)_g for all g ⇒ every pattern ≥ 0. b̂ ≥ 16: FF_g ≥ 0 and curvature part ≥ 0 ✓.

**Object identity** (`chk_consistency.py`): (i) `prl/data/nf02_psi.json` and `robust/functionals_noforward.json` (Part F's input) have
the same exact coefficients and t0_common's ψ = unnormalised/k₀ **exactly** (k₀ = 27.48509817 is the defining exact decimal, so t3's
`K0` scaling of the shift is exact, not a rounding); a = 0 (asserted by `nf_cert.load`), so `rcert`'s u³ forward term vanishes.
(ii) The float hybrid's ψ is Δ_g-independent (identical at 300/1000/5000). (iii) `rcert.poly_u(g, J)/k₀` equals the hybrid's
(1−cs₂)B₂ + (1−cs₃)B₃ from `flat_parts_nf` to ≤ 5·10⁻¹¹ relative over (h, J ≤ 64, u, all patterns); for `+−` both use
P^{(0,8)}_{J−4}(1−2up²) without the (1−up²)⁴ factor (I checked the hypergeometric ratio (n−k+1)(n+k+8)/k² ⇔ ₂F₁(−n, n+9; 1; ·)).

**Matrix claims**: even J, (++, +−) — the hybrid has no off-diagonal (POS_ERRORS §2.1; off-diagonals appear only in the
beyond-leading-Regge `dF_off`), so λ_min of the N_s-normalised diag is min(H₊₊/N₊₊, H₊₋/N₊₋) ≥ c ✓ (cited, consistent with the
code). Odd J: 1×1 `+−` (I5). OO: see I1.

## 3. The curvature lemma (Lemma 2), re-derived

Independent derivation (all steps agree with T0_CERT §2.1):
1. Mehler–Dirichlet with sinh(s/2) = sinh(β/2) sin φ gives P_{−½+iν}(cosh β) = (2/π)∫₀^{π/2}cos(νs)/cosh(s/2)dφ (Jacobian checked);
   with x = b̂ sin φ = 2Λ sinh(s/2): τ = Λs = 2Λ asinh(x/2Λ), w = (1 + x²/4Λ²)^{−½}, and C_k = (2/π)∫[wΦ(τ) − Φ(x)]dφ ✓.
2. Φ = r − a₁/x² (a₁ = ψ′(0); ψ′(1) = 0); the a₁ part is a₁M₁(β)/Λ² with x = Λσ, τ = Λt. The φ → t change of variables gives
   (√2/π)∫K c dt with c = cosh(t/2)/(4sinh²(t/2)) − 1/t² = LRINF_TIGHT §1.1's form ✓.
3. y = x/2Λ ≤ 1 (b̂ ≤ 2Λ). asinh and (1+y²)^{−½} series alternate with decreasing terms (coefficients 1, 1/6, 3/40, … and
   1, 1/2, 3/8, 5/16, …), giving 0 ≤ δτ ≤ 3ε²x⁵/640, 0 ≤ δw ≤ 3ε²x⁴/128, 0 ≤ x − τ ≤ εx³/24 ≤ x/6 ✓. The remainder decomposition
   δw·r + [wδτ + (1−w)εx³/24]r′ + wR₂ with 1 − w ≤ εx²/8, w ≤ 1, |R₂| ≤ ε²x⁶ sup|r″|/1152 reproduces **h̄ exactly** ✓.
4. g₁ = −(x³r)′/24; x³Φ′ = −x²∫(ψ + pψ′)cos by parts, so 3x²Φ + x³Φ′ = x²∫χ cos; two more parts with χ(1) = χ′(1) = 0,
   χ′(0) = a₁ give x²∫χ cos = −a₁ − ∫χ″cos, hence g₁ = (1/24)∫χ″cos(px)dp ✓; then (2/π)∫cos(pb̂ sin φ)dφ = J₀(pb̂) ✓.

**The bound on h is explicit and correct**: h̄ as written; H_k = sup h̄ certified (Lemma 4). I checked the H_sup code: the
crude/asymptotic majorants R_n are valid for all x > 0 (the asymptotic one is the exact finite integration-by-parts sum: odd j ≥ 3
at p = 0, j ≥ 5 at p = 1, derivative factors correct) and decreasing, so evaluating them at the cell's left end and the powers of x
at its right end is sound; the Taylor cells on [0.5, 80] use |r^{(n)}(c)| + h|r^{(n+1)}(c)| + h²/2|r^{(n+2)}(c)| + h³/6·sup|r^{(n+3)}|
with a valid crude sup; the r″ sup on [5x/6, x] covers the Taylor/majorant overlap correctly; on x ≥ 4000 every term is
x^{3−j} (j ≥ 3) or x^{c−a} (a ≥ 6) and so non-increasing ✓; on (0, 10⁻³] every term is increasing ✓.

**Lemma 3** (F_k enclosures): series tail ratio ≤ 1/16 (4/3 safe), derivative tail with 2(M+2)/b ✓; |F″| ≤ ∫|χ″|p² from
|J₀″| = |J₀ − J₂|/2 ≤ 1 ✓; the 1/b² tail via (pJ₀′(bp))′ = −b²pJ₀(bp) with f(1) = f′(1) = 0 ✓.
**Lemma 5** (flat-Regge remainder): re-derived; pf″ + f′ = f₁ + pg, boundary terms g(1)bJ₁(b), g′(1)J₀(b), and
∫_b^∞J₀ = −J₁ + J₀/b − ∫J₀/t² ⇒ |·| ≤ 1 + 2/b ✓. The −f₁/b³ = −a₂/(2b̂³) term matches LR∞'s flat tail (f₁ = ψ″(0)/2) ✓.

**Numerical tests** (`chk_lemma2.py`, `chk_curvpattern.py`, mpmath dps 40, [DV]):

| test | result |
|---|---|
| M₁: Lemma 2's φ-form vs LRINF's t-form, β = 0.05 … 4 | equal to ≤ 10⁻²¹ |
| pointwise \|h(x)\| ≤ h̄(x), k = 2, 3, Λ = 3, 10, 40, 300, x ∈ [0.05, 600] ∩ (0, 2Λ] | **0 violations**, max \|h\|/h̄ = 0.44 (k=2) / 0.68 (k=3) |
| **direct** Λ²C_k (∫ψ[₂F₁(½−iν, ½+iν; 1; −sinh²(β/2)) − J₀]dp, no Mehler–Dirichlet) vs a₁M₁ + F_k/24, at Λ = 40, 298.5, 998.5, b̂ = 1 … 0.95·2Λ | h_eff = Λ²(Λ²C − a₁M₁ − F/24) is Λ-independent (e.g. b̂ = 19.6, k = 3: 190.7 / 180.5 / 180.4) and **\|h_eff\| ≤ 191 ≤ 0.016·H_k** everywhere |
| Lemma 5 at b = 3.5 … 2000 | \|E_fl\|b³ ≤ e(b) everywhere (tightest 1065 vs 3577, k = 3, b = 10) |
| conclusions of K/T/L on actual curvature (all g, σ, u ∈ {0,1}), Λ(300), Λ(1000), Λ(10⁴), b̂ to β = 12 | negative only for b̂ ≤ 12.99; worst negative −8.2/Λ² vs T_abs = 29/Λ²; min Λ²·curv_g on b̂ ≥ 16 is +0.028 (at β = 12, ∝ M₁) |

## 4. Rigour of the certificate

* **Inter-cell bounds — proved.** Part K: second-order Taylor with the proved |F″| bound on 0.1-cells; M₁ by certified cell lower
  bounds (min over all cells meeting (0, β_max], β_max at Λ_min) — I checked the 910 M₁ cells tile [0.1, 64] with **exactly matching
  rational endpoints**, and M1_small(0.1) = 0.041578 < M₁(0.05) = 0.041656 is consistent. Part L: LRINF §1.4 monotone pieces per
  cell (cited [P + certified]). Part F(a): exact; F(b): CERT_CLOSE/ROBUST region-S boxes (mean-value with derivative enclosures,
  Cauchy tail), copied with the shift subtracted from the constant Taylor coefficient only (`cell_shift` is otherwise line-identical
  to `rcert.cell`). H_k: §3 above. Nothing is sampled.
* **Domains — complete.**
  * b̂: [0, 16] (T + F) ∪ [16, 300] (K cells) ∪ [300, 2Λ] (K tail; 300 < 2Λ_min asserted) ∪ β > 2 asinh 1 (L: 564 cells from the
    one containing 1.7627, plus the β ≥ 64 tail) — closed overlaps at 16, 300, 2Λ.
  * J: Part F(a) J < 100 (each J an exact polynomial on u ∈ [0, min(1, (16/(2J+1))²)], i.e. exactly b̂ ≤ 16); F(b) J ≥ 100 via
    u ≤ (b̂/201)², with the (X, u) expansion Π(X − λᵢu) valid for real J (the `++` series terminates at integer J); K/T/L are
    J-independent. The J = 100 handoff is clean.
  * u ∈ [0, 1] (m ∈ [Λ, ∞]): affinity in u at fixed (b̂, Λ, σ, g) with S_g > 0 at both ends ⇒ on the segment ✓; m → ∞ is u = 0.
  * Λ: every Λ-dependence enters via ε ≤ ε₀ = 1/Λ_min², M₁ minimised over β ≤ β_max(Λ_min), T_abs = T_gap/Λ_min², Λ²|E_fl| decreasing
    in Λ, and Λ²N̂(Λ) ≤ Λ_a²N̂(Λ_a). So "every real Δ_g ≥ DGMIN" is covered and the result is uniform as Λ → ∞.
* **Exact arithmetic in the 584 polynomials — yes.** Fractions → `to_fmpz_poly` (positive LCM) → exact fmpz Descartes/bisection
  (`P(0) > 0`, no sign variation on dyadic subintervals); the shift is a rational ≥ k₀·T_abs (rounded up).
  Part F's threshold is read from the json as mid + rad + 10⁻¹⁵ (Arb's printed ball contains the value) ✓.
* **LR∞ anchoring — valid.** GAP_LARGESPIN §1.2: with R̂ = min(R_∞, V/s⁴), every piece of Λ⁴N̂(β; Λ) is non-increasing in Λ
  (Λ⁴T_a ∝ Λ^{3−p}, p ≥ 3; Λ⁴T_b ∝ Λ^{4−p}, p ≥ 4; Λ⁴R̂(Λt) = t⁻⁴(s⁴R̂)(Λt); Λ⁴T_c₂ ∝ Λ^{3−j}, j ≥ 3; the [0, t₁] piece ∝ Λ⁻⁴),
  hence Λ²N̂(Λ) ≤ Λ_a⁴N̂(Λ_a)/Λ² ≤ Λ_a²N̂(Λ_a) for Λ ≥ Λ_a ✓. (Λ_a = lower end of the Λ_min ball ✓.) See I2 for Δ_g < 1000.
* β ≥ 64 tail of Part L: the flat ratio uses M₁ ≥ (√2/π)K₀ℓ/β = 2ℓ/(πβs) and β/s² decreasing ✓.

## 5. Hostile spot checks (`chk_eval.py` + `chk_scan.py`, [DV])

`chk_eval.py` is my own evaluator (flat: three-term recurrence at Gauss nodes; curvature: Lemma-C φ-form with Φ by Gauss (x ≤ 60) or
the exact integration-by-parts sum, valid at every β). It agrees with `ls_core.hyb_fast` to 10⁻¹⁰ on test states.
**Caution found on the way:** evaluating ψ in the *expanded* monomial form (coefficients up to 2.6·10⁶) loses the (1−p)⁵ zero to
round-off, and the `+−` backward peak (~J⁸ at u = 1) turns this into spurious *negative* flat actions at `+−`, J ≳ 50, m = Λ
(a false "H/N_s = −0.17"). ψ must be evaluated in the factored form (1−p)⁵pX/k₀ (as the certificate's exact arithmetic does).
After the fix the evaluator matches the exact Lemma-P polynomials to 10⁻⁹ there.

Fresh random grid per gap (seeds 11/12/13, 1534 states each): every `++` even J ≤ 60 and `+−` J = 4…60 at m = Λ, the twist edge
and 6 random masses (m/Λ ≤ 200); 150 near-forward states (b̂ ∈ [10⁻³, 2]); 500 states J ∈ [60, 3000], b̂ log-uniform to 2·10⁴
(i.e. past 2Λ at 300/1000); 60 huge-J states J ∈ [10⁴, 10⁵]; 120 twist-edge states J ≤ 3000.

| Δ_g | min H/N_s (where) | twist edge min | near-forward min | J ≥ 60 min | huge J min | # < 0.1999 | # < 0.199 |
|---|---|---|---|---|---|---|---|
| 300 | 0.1999996902 (`+−`, J=4, m=Λ) | 0.20427 (`+−`, J=5) | 0.9819 | 0.7849 | 0.99998 | **0** | **0** |
| 1000 | **0.1999999834** (`+−`, J=4, m=Λ) | 0.20090 (`+−`, J=5) | 0.9898 | 0.8010 | 0.9999998 | **0** | 0 |
| 5000 | 0.2000000110 (`+−`, J=4, m=Λ) | 0.20016 (`+−`, J=5) | 0.9942 | 0.7945 | 0.9999997 | **0** | 0 |

(The only sub-0.2 values are `+−`, J = 4, m = Λ: 0.19999969 at 300 and 0.19999998 at 1000.) The eight lowest states at every gap
are `+−`/`++` J = 4–19 at m = Λ, where 𝔉/N = 0.20000001–0.20000006; H/N_s → 𝔉/N as Λ grows, consistent with sharpness. The
c = 0.1999 statement at Δ_g = 1000 has a true margin of only 1.0·10⁻⁴ at its worst state, but it is certified, not sampled.

## 6. Scope: what T0 does NOT cover

1. **Realised vs target weights (M20).** T0 is about the de-shifted (target-weight) hybrid, ψ₂(ν/Λ) exactly. For the realised Ω^EE
   it holds only at Δ_g ∈ {300, 1000}, given ε_W [DV], and only for β ≲ 200 (T0 §4.2, a [certified + DV] corollary; its arithmetic
   checks: 2(1+c)ε̂_W ≈ 2.6·10⁻⁸ ≪ Part-F slack 6·10⁻⁷ / 9·10⁻⁷). M20 stays open.
2. **Beyond leading Regge (δ).** Nothing about δ_h (the Γ/m² corrections, the beyond-leading-Regge off-diagonals) — i.e. nothing about
   (Pos) of the full action.
3. **OO block.** Conditional on helicity lemma (A) [D], **and** on EO cancellation in the S⁰ sum (M17, open) and κ > 0 at all Δ (I1).
4. **Inherited unreviewed code.** Lemma 1 (𝔉 ≥ N/5) rests on CERT_CLOSE, whose `+−`, J ≥ 400, 3/8 ≤ u ≤ 1 corner (`cert/nf_corner.py`,
   Lemma B) is marked unreviewed there; I did not review it. Also cited, not re-derived: LR∞ (LARGESPIN_L §3 / LRINF_TIGHT), the R̂
   scaling lemma, the lrt2 M₁ cells and region-(c) bound, the small-β M₁ bound (GAP_LARGESPIN §3.2), the region-S machinery
   (ROBUST/NOFORWARD), POS_ERRORS' vanishing hybrid off-diagonal, and the definition of the hybrid itself (`ee1000_common.hyb`:
   curvature attached to the "2" of (2 − up²)ψ₂ and to uσψ₃).
5. The statement is about Δ_g ≥ 300 only; below 300 Part L's ratio would exceed 1 (T0 §5).

## Grades

| claim | grade |
|---|---|
| T0, `++`/`+−`, c = 0.199 (Δ_g ≥ 300), c = 0.1999 (Δ_g ≥ 1000) | **[certified]** (confirmed; cited inputs as listed) |
| Lemma 0, Lemma 1 use, Lemma 2, Lemma 3, Lemma 5 | [P] (re-derived) |
| Lemma 4 (H_k) | [certified] (code read; numerically ≈ 60× conservative) |
| OO sector | [certified given (A), M17 EO cancellation, κ > 0] |
| hostile scans, Lemma-2 numerics | [DV] |

## Files
`ads4/code/t0check/`: `chk_consistency.py` (§2), `chk_lemma2.py` (§3), `chk_curvpattern.py` (§3–4), `chk_eval.py` + `chk_scan.py`
(§5); copies of the t0cert scripts for reproduction. Logs: `cache/run_all_repro.out`, `cache/t*_*.txt`, `cache/chk_lemma2.log`,
`cache/chk_curvpattern.log`, `cache/chk_scan_{300,1000,5000}.log`.
