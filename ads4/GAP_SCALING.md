# GAP_SCALING: can (Pos) at Δ_g = 1000 be carried to all Δ_g ≥ 1000 by a scaling argument? (feasibility study)

*24 Nov 2026, gap-scaling agent. Unreviewed, not committed by this agent. Other agents' commits 8df1ce6, 9fc2f83 and 2565789 swept in
intermediate snapshots of `code/gapscale/`; the working tree holds the final versions. No existing file was edited; nothing in `prl/`
was touched.*
* *New code:* `ads4/code/gapscale/` (`run_all.sh`, ≈ 1.6 CPU-h; at most 2 processes; `ulimit -v 3000000`; `OMP_NUM_THREADS=1`).
  * Scripts: `gs_weights.py`, `gs_state.py`, `gs_lcscan.py`, `gs_lrinf.py`, `gs_lrball.py`, `gs_rows.py`, `gs_bind.py`, `gs_edge.py`, `gs_cover.py`.
  * Caches: `cache/`, whose `.gitignore` ignores everything but itself.
* *Read-only imports:* `deshift/d_abel`, `spinfunc/ee1000_common` (and `w2_chain`, `q3_analysis` through it), `lrtight/lrt_common` and
  `lrt3_rcert`, `neargap2/ng_chain` and `ng_bind` (run in-process with their own write guard, output to stdout only), the spinfunc
  row engines (`w1_rows`, `p3_odd_action`, `p2_shadow_odd`, `w0_ee_common`), `agproof/ag_exact`, `p4win/p4w_common`.
  * `gs_rows.py` re-types `neargap2/ng_rows.py`'s odd branch with every write redirected to `gapscale/cache/rows/`.
  * Every script that calls project engines installs a guard that refuses writes outside `gapscale/cache`.

Grades: **[P]** proved; **[certified]** Arb ball arithmetic with proved monotonicity; **[D]** derived; **[D-sketch]**; **[DV]** derived and
checked in floating point; **[NE]** numerical evidence; **[NE±]** with an error budget; **[M]** model; **[I]** open.
Notation as in the Supplement: `Λ = √((Δ_g−2)(Δ_g−1))`, `u = Λ²/m²`, `b = (2J+1)Λ/m`, `sinh(β/2) = (J+½)/m`, `x = ν/m`, `p = ν/Λ`,
`c₂ = m² ×` (relative correction).

---

> **Update (GAP_LARGESPIN.md, unreviewed): the large-spin floor is now certified for EVERY Δ_g ≥ 1000 (no upper limit), and the large-spin theorems' margins grow with the gap.**
> * **Lemma LR∞ for all Δ_g ≥ 1000 [P + certified]:** a scaling lemma makes Λ⁴·(each bound piece) non-increasing in Λ (envelope R̂ = min(R_∞, V/s⁴), V₂ ≤ 684625, V₃ ≤ 789636), plus a 13-step monotone Λ-grid below 2604.28. Thus r < 1 on β ≥ 0.37 and F_LR > 0 on β ≥ 0.245; J\*(m) = ⌈(0.37211m−1)/2⌉, with J\*/Λ → 0.22387 as Δ_g → ∞. This closes the tail above Λ = 10⁶.
>   * Thresholds improve with the gap: β₀ = 0.2375 / 0.1675 / 0.1175 for Δ_g ≥ 2606 / 5000 / 10⁴.
>   * Gap-unit floor: Λ²F_LR ≥ f_U(β) for all Δ_g ≥ 1000 [certified].
> * **ROUTEA curvature window uniform in the gap [D + certified]:** curv > 0 on b̂ ∈ [114.06, 2Λ] for m ≥ Λ (new exact identity ∫₀^{s₁}g₁ = −s₁³r(s₁)/24). **Lemma C′ is monotone in Λ [P]**, so b̂ ≤ 119 carries to every gap. **For m ≥ 1.3Λ the hybrid floor is rigorous on every b̂ ≤ 2Λ at every gap**, which closes the old (119, 245) hole for m in [1.3, 1.5)Λ.
> * **Theorems L′ / W at Δ_g = 1000 / 2000 / 5000 [DV, premises carried with Δ_g-independent AdS-unit constants]:**
>   * Λ²D = −995.3 / −999.1 / −1001.3.
>   * L′ remainder: max d/(F/2) = 0.19 / 0.035 / 0.005; edge tolerance 12.6 / 228 / 9050.
>   * **Theorem W: d_χ/H = 0.489 / 0.246 / 0.100; the tolerated graviton edge integral at the binding state is 0.406 / 3.45 / 57 S_E,t (∝ Λ^3.1).** So the EDGE-R premise gets much easier with the gap.
>   * Caveat: the 5000 C2/vector tolerance entries are not numerically converged (finite-difference β-derivative).
> * **Remaining:** the premises (P_L′ magnitudes, AG_U, EDGE-R, P1′, W) carry over only under the hypothesis of Δ_g-independent AdS-unit constants. The trend between computed gaps is [D-sketch]. The small-J twist-edge seam and β₀ < 0.10 are not addressed.

> **Update (GAP_UNIFORM.md, unreviewed): there is no uniform-in-gap proof for the near-gap corrections, and the literature has none either. The needed inequality is pinned down, and J ≤ 7 is positive at every gap given one explicit "no-growth" premise with large safety factors.**
> * **Literature [Q]:** no explicit uniform-in-Δ_gap bound exists.
>   * 2205.01495: the uplift is one sentence, with an unquantified O(1/Δ_gap⁴) in eq. (30).
>   * CHMRS §5: non-effective ("sufficiently large M", (5.7), (5.11), (5.26)–(5.29)).
>   * Knop–Mazáč AdS₂: their uniformity comes from pointwise kernel positivity against dDisc ((3.25)–(4.9)), which does not transfer to Ω^EE.
>   * CLS-D: leading order only.
>   * van Rees–Zhao 2210.15683: uniformity "not immediate".
> * **Lemma T [P]:** in chart A the node coefficients, the flat parts and N_f are exactly Λ-independent and the Lemma-C bound decreases, so positivity at a larger gap needs only that Γ/m² does not grow. Checked on 366 states at 1000/2000/5000.
> * **Prop. S [D]:** at fixed ξ = ν/Δ the crystal integral is ∫e^{ΔΦ_ξ}g_J with Φ_ξ exactly Δ-independent, so H_mono-A is ordinary fixed-phase Laplace asymptotics.
> * **Lemma R [D-sketch]:** the D′ thimble coefficients are Δ-independent in chart A (removes the finite-Δ Stokes worry). **Theorem E [D-sketch]:** each fixed near-gap state is positive for all Λ ≥ Λ\*(J, μ). Prop. D′ controls only ξ < 0.107.
> * **Coverage:**
>   * **J ∈ {0, 2, 4, 5, 6, 7} positive at every Δ_g ≥ 1000** [P given B_F + windows]. Premise B_F: sup|Γ_Δg| ≤ F_J·(sup|Γ| at 1000), with F_J = 11.8 (J=7), 14.9 (6), 43 (4), 743 (5), 799 (2), 3.4·10⁴ (0). B_F is [NE]: graviton 1000 → 2000 ratio 0.97–1.06; exact scalar 1000 → 10⁴ drifts −2.6 % to +1.2 %.
>   * J = 8 fails with the certified floor (F = 0.85).
>   * J ≥ 12: the sign-free form fails already at 1000.
> * **Main missing step:** explicit Laplace remainder constants C₃(ξ, J). They need interval sup-norm bounds of the spinning amplitudes on a fixed complex polydisc (checked once at Δ₀, valid for all Δ ≥ Δ₀), which needs a pointwise complex evaluator the project lacks.

## 0. Verdict: (b) partially viable

**In one sentence.** Every *certified or proved* layer of the Δ_g = 1000 positivity evidence carries over to every Δ_g ≥ 1000, by
exact invariance or proved monotonicity. The *beyond-leading-Regge* layer, which holds all the [NE±]/[NE + M]/conditional content,
carries over only under a scaling hypothesis. That hypothesis is now tested at up to five gaps, with favourable exponents everywhere
it was tested. It still has to be assumed. One region, the strip between the 1000 and the Δ′_g twist edges, has no Δ_g = 1000 counterpart
at fixed spin. It needs either a second (fixed-β) chart or per-gap computation.

**The precise statement that is within reach** (a "transfer proposition", [D-sketch + NE]; details in §5). For Δ′_g ≥ 1000, write the
per-state action as

  `λ/N_s (J, u; Λ) = Φ_flat(J, u) + C_LR(J, u; Λ) + Γ-cost(J, u; Λ)`

(flat part, leading-Regge curvature, and beyond-leading-Regge correction). Then:

1. **Φ_flat is exactly Λ-independent at fixed (J, u, h)** [P: the weights ψ_k(p) and the flat partial waves depend only on p, y = up², J].
   The NF02 certificate `Φ_flat ≥ N/5` holds for *all* u ≤ 1, so every new state at Δ′_g is covered, including those below the
   1000 twist edge [certified].
2. **The Lemma-C bound on C_LR decreases monotonically in Λ at fixed (J, u)** [P] (§2). The exact C_LR is `g(J, u)/Λ²·(1 + O(Λ⁻²))`
   [DV, 300 → 20000]. So the sign-free floor "certified flat − Lemma C" > 0 on a set of states that grows with Λ.
   * Largest spin with a sign-free floor as m/Λ → 1, worse helicity: J ≤ 13 / 17 / 23 / 35 / 61 at Δ_g = 1000 / 2000 / 5000 / 2·10⁴ /
     10⁵ [certified flat + P] (§2.2).
3. **Lemma LR∞ holds uniformly in the gap** [certified].
   * The certified remainder ratio satisfies sup_{β ≥ 0.37} r ≤ 0.9988 for *every* Λ ∈ [998.5, 10⁶] (Arb Λ-balls, 454 balls).
   * r(β; Λ) falls like Λ⁻²: β₀(r < 1) = 0.37 / 0.2575 / 0.16 at Δ_g = 1000 / 2000 / 5000.
   * So Theorem L′'s floor, and its threshold J\*(m) = ⌈(0.37211 m − 1)/2⌉, hold at every Δ_g ∈ [1000, ≈10⁶], with improving
     constants (§3.1).
4. **Beyond-leading-Regge corrections at fixed (J, u)** cost `∝ Λ⁻ᵖ` with p ≥ 2, *provided* Γ (in c₂ units) at fixed (J, x) has a finite
   limit as m → ∞. Call this H_mono-A [M + NE]. It is the statement that the flat-space (bulk-point) limit is approached at rate 1/m²
   at fixed spin. Tested:
   * **Graviton, binding state J = 5:** Γ(x) on seven K₂ rows Δ = 200 … 2005 stays within |Γ| ≤ 64 c₂ and changes by ≤ 6 c₂ between
     Δ = 1005 and 2005.
   * **Binding state at Δ_g = 2000** (own-mass row computed here; NEARGAP2 budget run unchanged): λ/N_s = +0.20041, sign-free budget
     6.3·10⁻⁴, **L = 320**, against L = 78 at 1000. E_trunc scales as Λ⁻² (×3.92), E_low as Λ⁻⁴ (×15.6) [NE±] (§4).
   * **Exact scalar at J = 26 and 48** (Arb), including the light-cone edge x → 0.97: Γ(x) in c₂ is Λ-independent to ≤ 4 % from 1000 to
     5000. The whole k = 2 cost beyond the hybrid falls like Λ⁻² (J = 26: Λ²S = −0.343 / −0.253 / −0.233 at 1000 / 2000 / 5000)
     [DV for the scalar] (§3.3).
   * DG_SCAN's graviton rows J = 8–24, 300 → 1000: p_emp = 2.83 (range 1.4–3.2) [NE].
5. **Fixed-β chart (near-gap large spin).** In the variables (β, u), the per-state window tolerance grows as Λ^{1.0–1.2} (DG_SCAN), and
   the exact scalar's k = 2 cost relative to the leading-Regge action falls roughly as Λ⁻³. At β ≈ 0.052, m = 1.0257Λ:
   |S/Ω₂| = 5.0·10⁻⁴ / 8.8·10⁻⁵ / 2.3·10⁻⁵ / 7.5·10⁻⁶ at Δ_g = 500 / 1000 / 1500 / 2000 [DV, scalar] (§3.3).
   This is the chart that covers the twist-edge strip at large J. The graviton's edge integral (EDGE-R) is not computed at any gap, so this
   chart inherits that hole unchanged.

**What this does not give.**
* It does not give a *proof* that (Pos) at 1000 implies (Pos) at Δ′_g. H_mono-A/B are statements about beyond-leading-Regge
  corrections, which are numerical at every gap. A scaling argument can at most transfer the 1000 grade ([NE±], [NE + M],
  [P given AG_U, EDGE-R, …]) to larger gaps. It cannot upgrade it.
* The structural premises (AG_U, EDGE-R, P1′, W, the +− transfer convention) are per-state statements independent of Δ_g. They stay
  premises.

**Realistic effort.**
* **A transfer proposition** ("Pos at 1000 ⇒ Pos at all Δ_g ≥ 1000, given H_mono-A/B; certified layers uniform; tested at 2000 and
  5000"): 3–5 days. Items 1–3 are in hand. Most CPU goes to re-running the error-budgeted near-gap spins at 2000 and 5000 (§6).
* **Per-gap repetition of the full 1000 evidence at one more gap:** ≈ 2–3 weeks of agent time and ≈ 300–400 CPU-h, most of it external
  (§6.2). Its graviton-edge step is blocked at 1000 too.

**Recommended framing for the referee's item 7.** "Positivity is established at error-control level at Δ_g = 1000 and, on the
error-budgeted states, at 2000 and 5000. For every Δ_g ≥ 1000 it follows from the 1000 evidence under an explicit, tested scaling
hypothesis. The certified layers hold uniformly." That would answer "one data point" without overclaiming. It is still not "(Pos)
proved for all Δ_g ≥ 1000".

---

## 1. Q1: how the functionals depend on Δ_g

### 1.1 Exact facts [P / D]
* **Flat weights.** ψ_k(p) (NF02; Q3\*) are fixed functions of p = ν/Λ (Supplement §S2: `Ω^EE = ∫dp ψ₂ Λ⁴Ĉ′_{EE,Λp} + …`). The
  rescaling ν → ν/Λ is the *definition* of the family; no correction.
* **Leading-Regge image.** By the de-shift identity, the leading-Regge action of the de-shifted functional is exactly
  `∫ψ_k(p) P_{−1/2+iΛp}(cosh β) dp` for every Λ. So the Λ-dependence of the de-shift weight X enters the heavy side only through
  beyond-leading-Regge corrections, for example the exact monopole cost `D ≈ −990/Λ²` (Supplement §S8.A), and the light side
  (δk(0) = −68.08/Λ², DG_SCAN).
* **The k = 2 de-shift weight X_Λ** is the minimal-norm solution of a Λ-*independent* linear Abel equation (in ν) with source
  ψ₂(ν/Λ). It is therefore a linear functional of the rescaled source. It has a two-scale structure: an outer region p fixed and an
  inner region ν = O(1–50), where the R_∓ kernel varies.

### 1.2 Numerical check (`gs_weights.py` → `cache/weights.txt`) [DV]
The weights were re-solved with the independent solver `deshift/d_abel.Abel` (CGLS, relative residual ≤ 6·10⁻¹⁶, 20 iterations,
≤ 15 s per gap) at 300, 1000, 2000 and 5000.

* **Outer region:** `Λ²[Re X(Λp) − ψ₂(p)] = S(p) + O(Λ⁻²)`, with the first-order symbol `S = ½ψ₂″ − (23/2)(ψ₂/p)′`.
  Residual `sup|Λ²(ReX−ψ₂) − S| / max|S|`:

  | p range | 300 | 1000 | 2000 | 5000 |
  |---|---|---|---|---|
  | p ≥ 0.1 | 0.37 | 0.029 | 0.0071 | 0.0011 |
  | p ≥ 0.2 | 0.037 | 0.0032 | 0.00079 | 0.00013 |

  The residual falls as Λ⁻², so the next term is O(Λ⁻⁴). At p = 0.5: 1651.8 / 1645.1 / 1644.6 / 1644.5 against S = 1644.4.
* **Inner region:** `Λ·X(ν) → Ξ(ν)`, a Λ-independent profile (the Abel solution for the source a₁ν), with O(1/Λ) corrections.
  * At ν = 5: ΛX = 1124 / 1180 / 1192 / 1200; Λ·Im X = 152 / 166 / 169 / 171.
  * Im X is O(1/Λ) and lives at ν ≲ 20.
* **Light-side shift:** `Λ²δk₀ = −68.31 / −68.10 / −68.08 / −68.077`, converging to the symbol value −68.076.
* **Answer to Q1.** Yes: up to explicit corrections, the weights at different gaps are related by ν → ν/Λ.
  * In the outer region the correction is O(Λ⁻²), with the explicit symbol S.
  * In the inner region ν = O(1) (a fixed ν-window whose weight in p-units shrinks like 1/Λ) the profile is Λ-independent.
  * Where the inner region costs anything per state (small-ν windows), the cost scales as Λ⁻⁴ at fixed (J, u) (§4: E_low).
* **Weights at other gaps.**
  * The min-norm grid weights are produced at any Δ_g in seconds (`d_abel.Abel(Dg, dnu, numax)`; `gs_weights.py` saves
    `cache/weights_<Dg>.npz`).
  * The published *continuum* weights exist only for 300 and 1000. `cert/k2_cont.py` (Target(W, Dg), Basis(Λ)) is parametrised by
    Δ_g and could produce them (not run here; CERT_CLOSE's whole compute, including these fits, was ≈ 2.5 CPU-h).
  * Result Q's k = 4, 5 de-shift weights are Λ-independent (functions of ν only).

---

## 2. Q2: how a heavy state's action scales, and whether the curvature bound is monotone

### 2.1 Gap-unit decomposition (`gs_state.py` → `cache/state.txt`)
At fixed (J, u = Λ²/m², h), b = (2J+1)√u is fixed and β = 2 asinh(b/2Λ) → 0. Each column below follows one gap-unit state from
Δ_g = 300 to 20000:

| state (J, m/Λ, h) | flat/N_f (all Λ) | Λ²·C_LR/N_f at 300 / 1000 / 5000 | LC/N_f at 1000 / 2000 / 5000 | sign-free lb at 1000 / 2000 / 5000 |
|---|---|---|---|---|
| 5, 1.0 (limit of twist edges), +− | 0.200000 | 0.7694 / 0.7696 / 0.7696 | 6.7·10⁻⁵ / 1.7·10⁻⁵ / 2.7·10⁻⁶ | 0.19992 / 0.19998 / 0.20000 |
| 8, 1.0, +− | 0.200000 | 19.499 / 19.500 / 19.500 | 3.9·10⁻³ / 9.8·10⁻⁴ / 1.6·10⁻⁴ | 0.1953 / 0.1988 / 0.1998 |
| 12, 1.0, +− | 0.227467 | 70.17 / 70.20 / 70.20 | 0.053 / 0.013 / 0.0021 | 0.165 / 0.211 / 0.225 |
| 16, 1.0, ++ | 1.000000 | 592.96 / 593.44 / 593.49 | 0.78 / 0.19 / 0.031 | 0.127 / 0.676 / 0.940 |
| 21, 1.0, +− | 0.202193 | 184.25 / 184.64 / 184.68 | 0.65 / 0.16 / 0.026 | −0.27 / +0.034 / +0.172 |
| 26, 1.026, ++ | 1.000000 | — / 2383.1 / 2383.7 | 10.6 / 2.7 / 0.42 | −0.83 / −0.45 / +0.40 |
| 32, 1.3, ++ | 1.000000 | 2074 / 2080 / 2080 | 7.5 / 1.9 / 0.30 | −0.76 / −0.30 / +0.54 |

(m/Λ = 1.0 is not twist-allowed at any finite gap. It is the Δ_g → ∞ limit of every fixed-J twist edge, and it is covered by the flat
certificate.)

* **Flat part:** identical at every Λ [P].
* **Curvature:** C_LR = g(J, u)/Λ² to 3–4 digits [DV].
  * The binding state's curvature is +7.8·10⁻⁷ (1000) → +1.9·10⁻⁷ (2000) N_s.
  * Λ²·C_LR loses digits at Δ_g = 2·10⁴ when |C_LR/N| ≲ 10⁻⁹ (mpmath double-precision hypergeometric), e.g. (5, 2Λ).
* **Lemma C bound** `|P_{−1/2+iν}(cosh β) − J₀(νB)| ≤ (1 − sech(β/2)) + min(2, ν(B − β))`, B = 2 sinh(β/2):
  * at fixed (J, u), 1 − sech(β/2) decreases in Λ;
  * ν(B − β) = b p (1 − asinh(s)/s) with s = b/2Λ, and asinh(s)/s is decreasing, so this term decreases in Λ.
  * Hence LC(J, u; Λ) is **pointwise monotone decreasing in Λ** [P].
  * The sign-free floor `lb = (Φ_flat − LC)/(N_f + LC)` is decreasing in LC (since N_f ≥ |Φ_flat|), hence **increasing in Λ** [P].
  * This is the monotonicity the brief asked about. It holds for the leading-Regge layer and nothing else (Lemma C cannot bound Γ;
    SMALLNU_WINDOW).
* **Lemma C′ is not a curvature bound.** It is the light-side analyticity lemma (Supplement §S4.B) and is Δ_g-independent.

### 2.2 Region where the certified-flat-minus-Lemma-C floor is positive (`gs_lcscan.py` → `cache/lcscan.txt`)
Largest J with `N_f/5 − LC > 0` for all allowed J′ ≤ J, worse of the two helicities (certified flat floor). With the quadrature flat
value in place of N_f/5 the numbers are 1.0–1.3× larger (`cache/lcscan.txt`).

| m/Λ | Δ_g = 300 | 1000 | 2000 | 5000 | 2·10⁴ | 10⁵ |
|---|---|---|---|---|---|---|
| 1.0 (limit) | 8 | 13 | 17 | 23 | 35 | 61 |
| 1.05 | 9 | 13 | 17 | 23 | 37 | 65 |
| 1.3 | 10 | 17 | 21 | 29 | 48 | 83 |
| 2.0 | 16 | 27 | 34 | 47 | 76 | 132 |
| 5.0 | 43 | 70 | 89 | 122 | 197 | 339 |

* The region grows like Λ^{1/3}, because LC/N_f ∼ J⁶/Λ² near the gap [certified flat + P monotone; integrals by quadrature, DV].
* Inside it, the hybrid (leading-Regge) floor is rigorous at every Δ_g ≥ the gap where it first holds.
* Outside it, the hybrid floor is a computed number (DG_SCAN's λ₀: min 0.2009 → 0.2004 → 0.2002 at 1000 / 2000 / 5000, decreasing to the
  certified 0.200). It could be certified by interval quadrature of the explicit 1D curvature integral (≈ 1 day; not done).

### 2.3 Answer to Q2
* In gap units the per-state value is `Φ_flat(J, u) + g(J, u)/Λ² + (beyond-Regge)`. The first two terms are smooth, with the flat
  limit as Λ → ∞ [P + DV].
* The certified Lemma-C version decreases monotonically [P].
* "Positivity at 1000 with the curvature bounded ⇒ positivity at all larger gaps" is therefore **true for the leading-Regge layer**, on
  the states where it was established at 1000 [P].
* It is **not** a statement about the full action: everything [NE±] at 1000 is beyond-Regge (§4, §5).

---

## 3. Q3: the large-spin theorems as functions of Λ

### 3.1 Lemma LR∞ / Theorem L′ floor (`gs_lrinf.py`, `gs_lrball.py`) [certified]
* **Per gap.** `lrtight/lrt3_rcert.py` was re-run with Λ patched; Λ-independent inputs (∫|ψ_k|, M₁ cells) were copied. It reproduces
  1000 exactly (β₀ = 0.37 / 0.245):

  | Δ_g | sup_{β≥0.37} r | sup_{β≥1} r | β₀ (r < 1) | β₀ (F_LR > 0) |
  |---|---|---|---|---|
  | 1000 | 0.997 (ball [999.5, 1000.5]) | 0.0162 | 0.37 | 0.245 |
  | 2000 | ≤ 0.276 (sup over β ≥ 0.35) | 0.00387 | 0.2575 | 0.135 |
  | 5000 | ≤ 0.041 (sup over β ≥ 0.35) | 0.00060 | 0.16 | ≤ 0.10 (grid floor) |

  At fixed β, r falls as Λ⁻² (sup_{β≥1}: ×4.19 for Λ×2, ×6.43 for Λ×2.5).
* **Uniformly in the gap.**
  * The same certification was run with Λ an Arb **ball**, on 454 adjacent balls covering [998.4999, 10⁶]:
    `cache/lrball_part1.log`, 998.5–1652 (ratio 1.001 per ball), and `cache/lrball_resume.log`, 1550–10⁶ (ratio ≤ 1.05).
  * Result: **sup_{β≥0.37} r ≤ 0.9988 on every ball**, worst at the first. So Lemma LR∞ with β₀ = 0.37, and Theorem L′'s
    floor with J\*(m) = ⌈(0.37211 m − 1)/2⌉, hold for every Δ_g between 1000 and ≈ 10⁶ [certified].
  * Each ball also certifies its own smaller β₀, e.g. 0.265 on the ball containing Δ_g = 2000 and 0.17 at 5000.
  * Balls wider than ≈ ×1.07 lose all precision (an Arb power of a wide ball contains 0); the script now asserts on NaN.
  * Above 10⁶ a tail argument is needed: every term of `Λ²·|E_k|` is a decreasing function of Λ at fixed β, except the
    `∫|ψ|·t⁴` piece of T_c1, which is small. This is [D-sketch]; not done.
* **Theorem L′'s tolerances** (Λ_tol, K^sh_tol, ε_R) were not recomputed. `p4win`, `agproof` and `lsw` hard-code Δ_g = 1000 (module
  constants). Structurally, in the fixed-β chart the floor ∝ Λ⁻² and the O(1/m²) costs relative to it ∝ 1/m² [D-sketch].

### 3.2 H_LJ window tolerance (DG_SCAN, recomputed there)
min T_B3 (sign-free window tolerance, c₂) = 1201 / 2816 / 8821 at 1000 / 2000 / 5000, i.e. ∝ Λ^{1.2}. At fixed J/Λ it is ∝ Λ [D for the
weights]. Against a Δ_g-independent window correction (≈ 250–450 c₂ measured at 300/1000) the margin grows [M].

### 3.3 Theorem W's corrections and the edge premise, on the exact scalar (`gs_edge.py`) [DV for the scalar]
Exact scalar heavy blocks (Δ_φ = 3, k = 2, Arb Laplace/Barnes; `agproof/ag_exact`, read-only). The correction `G = exact − hybrid`,
`Γ = m²G/E₂ₕ` in c₂ units, and the whole k = 2 cost `S = ∫₀¹ψ₂G dp` relative to the leading-Regge action Ω₂ = ∫ψ₂P dp.

| chart | state (Δ_g: J, m/Λ) | β | Γ(x = 0.95) c₂ | S/Ω₂ | Λ²·S |
|---|---|---|---|---|---|
| A (fixed J, u) | 1000: 26, 1.0257 | 0.0517 | −4.72·10⁴ | −8.8·10⁻⁵ | −0.343 |
| A | 2000: 26, 1.0257 | 0.0259 | −4.91·10⁴ | −1.6·10⁻⁵ | −0.253 |
| A | 5000: 26, 1.0257 | 0.0103 | −4.91·10⁴ | −2.3·10⁻⁶ | −0.233 |
| A | 1000: 48, 1.0470 → 2000: 48, 1.0470 | 0.093 → 0.046 | — | 2.7·10⁻⁵ → 3.6·10⁻⁶ | 0.017 → 0.009 |
| B (fixed β, u) | 500: 13 / 1000: 26 / 1500: 39 / 2000: 52, all at m/Λ = 1.0257 | 0.052 | — | +5.0·10⁻⁴ / −8.8·10⁻⁵ / +2.3·10⁻⁵ / −7.5·10⁻⁶ | — |
| twist edges, fixed J | 26 at 1000 / 2000 / 5000 (x₁ = 0.975 / 0.987 / 0.995) | — | −4.72 / −4.80 / −4.77·10⁴ | −8.8·10⁻⁵ / +3.1·10⁻⁵ / −7.4·10⁻⁶ | — |

* **Chart A.** Γ(x) at fixed (J, x) is Λ-independent in c₂ units *including the light-cone edge* x = 0.90–0.97 (x = 0.97, twist-edge
  states: 1.48 / 1.66 / 1.70·10⁵). So the cost falls as Λ⁻² (Λ²S converges). This is H_mono-A on the exact scalar.
* **Chart B.** At fixed (β, u) the relative cost falls roughly as Λ⁻³, through stronger oscillatory cancellation (b ∝ Λ). Pointwise,
  however, Γ in c₂ grows like J² at fixed β, i.e. the *relative* pointwise correction ~ β² is Λ-independent.
  * Consequence: sign-free per-node budgets (as in the K₂ λ-level error terms) do **not** improve in chart B. Only signed integrals
    do.
  * A first run that cut the integral sharply at p = 0.2 (`cache/edge_run1.txt`, `edge_run2.txt`) showed the opposite trend. That was
    an artefact of the hard cut through an oscillating integrand. The smooth-weight full-range values above supersede it.
* **What this says about EDGE-R (Theorem W's thinnest premise).**
  * The tolerance at (26, 1.026Λ) is set by the +− transfer factor (1 − x₁²)⁻² = 408.
  * In chart B, x₁ is fixed, so the scalar-like edge margin improves with Λ.
  * At the *twist edge* in chart A, x₁ → 1 (x₁ = 1 − J/Δ_g), so the transfer factor grows like (Δ_g/J)². At J = 26 that is 408 → 1550
    → 9400 at 1000 / 2000 / 5000, faster than the scalar cost falls (Λ^{−1.5} at the twist edge). **The near-edge strip at large J
    must therefore be taken in chart B, not chart A.**
  * The graviton's edge integral is not computed at any gap [M]. A scaling argument cannot supply it.

### 3.4 Answer to Q3
* **LR∞ / L′ floor:** improves (β₀ ∝ roughly Λ^{−1/2}: 0.37 → 0.2575 → 0.16) and is now certified uniformly on [1000, 10⁶].
* **Theorem W:** its window tolerances improve (∝ Λ^{1.0–1.2}). Its sampled-state costs transfer in chart A (∝ Λ⁻², scalar-tested).
  Its graviton edge premise is untouched, and at the twist edge its +− transfer convention gets *worse* with Δ_g at fixed J.
* **Theorem L′'s ε_R and tolerances:** not recomputed (hard-coded Δ_g); expected to improve in chart B [D-sketch].
* **Recomputation at 2000/5000:** cheap for LR∞ (done). For Theorem W it needs the p2con O(1/m²) engine and ag_chi/p4win
  parametrised in Δ_g (§6).

---

## 4. Q4: near-gap margins; the binding state at Δ_g = 2000 (`gs_rows.py`, `gs_bind.py` → `cache/bind.txt`)

**Own-mass K₂ row** Δ = 2005, J = 5 (parity-odd +−, 36 nodes × {EE, C3H} × {state, Δ_ref = 2·10⁴}). 142 points, 16.4 CPU-min, engines
unchanged.

**H_mono-A on the graviton at J = 5** (Γ in c₂, K = 2; rows 200, 400, 800, 1600 published, 305/1005 NEARGAP2, 2005 new):

| x | C2: Δ = 200 / 305 / 400 / 800 / 1005 / 1600 / 2005 | C3H: same rows |
|---|---|---|
| 0.3 | 5.6 / 12.9 / 15.2 / 19.0 / 20.3 / 21.3 / 21.8 | −41.0 / −37.9 / −35.4 / −31.5 / −31.1 / −29.9 / −29.7 |
| 0.5 | −13.8 / −5.8 / 4.5 / 19.7 / 19.8 / 24.9 / 25.3 | −11.2 / −16.0 / −24.8 / −37.5 / −36.9 / −40.9 / −40.5 |
| 0.7 | −9.6 / −6.3 / −1.1 / 6.6 / 6.4 / 9.0 / 9.2 | 1.9 / −1.7 / −7.0 / −14.8 / −14.7 / −17.2 / −17.1 |

* Γ converges to an m-independent profile of size ≤ 41 c₂. The approach is O(J²/m) (K₂'s parameter), and 1005 → 2005 differs by ≤ 5.5 c₂.
* So the correction's per-state cost scales as Λ⁻² at fixed (J, u) [NE].

**Budget at Δ_g = 2000** (NEARGAP2's `run_state`, unchanged; gauge factor from the 2000 min-norm weight):

| Δ_g, m/Λ | λ/N_s (own row) | E_trunc | E_low (ν < 50) | E_edge | total | lower bound | L (sign-free / λ-level) |
|---|---|---|---|---|---|---|---|
| 1000, twist edge 1.00499 (reproduced) | +0.20088 | 1.91·10⁻³ | 3.1·10⁻⁴ | 3.7·10⁻⁴ | 2.59·10⁻³ | +0.1983 | 78 / 130 |
| **2000, twist edge 1.00250** | **+0.20041** | 4.87·10⁻⁴ | 2.0·10⁻⁵ | 1.2·10⁻⁴ | 6.3·10⁻⁴ | **+0.1998** | **320 / 565** |
| 2000, 1.00499 (chart-A image of the 1000 edge) | +0.20089 | 4.89·10⁻⁴ | 2.0·10⁻⁵ | 9.2·10⁻⁵ | 6.0·10⁻⁴ | +0.2003 | 335 / 609 |
| 2000, 1.01 / 1.02 / 1.05 | +0.20212 / +0.20561 / +0.22621 | ≈ 5·10⁻⁴ | 2·10⁻⁵ | ≤ 6·10⁻⁵ | ≤ 5.8·10⁻⁴ | — | 352 / 370 / 394 |

* Certified-floor route at 2000: `(N_f/5 − LemmaC)/N_s = +0.19998` (1000: +0.19993).
* The chain built from the *published* rows only, i.e. Γ held from the Δ = 1600 row (an H_mono-A transfer), gives +0.20041, equal to
  the own-row value to < 10⁻⁵.
* **Answer to Q4.**
  * The binding state's margin at 2000 is at least that at 1000 (both → 0.2 + O(Λ⁻¹) from above, with a ×4 smaller budget).
  * The mechanism is general for fixed-J near-gap states: E_trunc ∝ Λ⁻² (×3.92), E_low ∝ Λ⁻⁴ (×15.6), E_edge ∝ Λ^{−1.6} [NE±]. The
    1/Δ_g expansion of the K₂ rows is the reason [D-sketch: leading Laplace coefficients are Δ-independent in c₂ units].
  * The new twist-edge strip (m/Λ ∈ [1.0025, 1.005) for J = 5) is computed directly here. For J ≤ 25 the strip is ≤ 0.025 wide in m/Λ.

---

## 5. The scaling argument stated precisely, and what is in hand

**Charts.**
* **A** (fixed J, u): the Δ′_g state has a Δ_g = 1000 counterpart iff u⁻¹ ≥ 1 + 2J/999.
* **B** (fixed β, u): J′ + ½ = (J + ½)·m′/m. It is quasi-continuous in J, so it is useful only for large J, and it maps twist edges to
  twist edges up to O(1/Λ).

**Hypotheses** (both [M]; tested as listed).
* **H_mono-A:** for every heavy state, Γ(J, x; m) in c₂ units at fixed (J, x) is bounded uniformly in m ≥ m₁₀₀₀(J, u), so that its
  per-state cost scales ≤ (Λ₁₀₀₀/Λ)². Tested on graviton J = 5 (7 rows), DG_SCAN's J = 8–24 (300 → 1000), and exact scalar J = 26, 48
  (1000 → 5000, incl. the edge).
* **H_mono-B:** at fixed (β, u) the signed costs of beyond-hybrid corrections relative to the floor do not grow with Λ, and the window
  correction at fixed (β, ν) is Δ_g-independent (DG_SCAN's H_mono, second clause). Tested on the exact scalar, β ≈ 0.052, 500 → 2000.

**Transfer proposition** [D-sketch given H_mono-A/B; certified parts as marked]. Let Δ′_g ≥ 1000.
* (i) Every state with a chart-A counterpart inherits the 1000 grade of that counterpart, with every beyond-Regge cost multiplied by
  ≤ (Λ₁₀₀₀/Λ′)². The leading-Regge layer is pointwise no worse [P, §2].
* (ii) Every state with b ≥ 0.37211Λ′ has the certified L′ floor [certified, §3.1].
* (iii) Every remaining near-gap state with J ≥ J_B (J_B ≈ 52 at 2000) has a chart-B counterpart in the 1000 near-gap band, and
  inherits its Theorem W statement under H_mono-B. The premises AG_U, EDGE-R, P1′ and W are assumed at Δ′_g as at 1000.
* (iv) The remaining **seam** consists of the small-J twist-edge strips, J ≲ J_B with m/Λ between the Δ′_g and 1000 twist edges.
  * Width ≤ J(1/999 − 1/Δ′_g) in m/Λ; at J = 5 it is [1.0025, 1.005) at 2000.
  * It needs H_cont or direct rows. The binding state's strip is done here.

Coverage bookkeeping without chart B (`gs_cover.py` → `cache/cover.txt`):
* Chart A plus the certified L′ floor leave every spin J ≤ 344 (2000), 424 (3000), 552 (5000) uncovered on [twist edge(Δ′), m₁₀₀₀(J)).
* E.g. J = 300 at 2000: m/Λ ∈ [1.140, 1.265). L′ starts at J\*′ = 303 (twist edge) at 2000, 463 at 5000.
* Chart B is therefore essential, and it rests on H_mono-B, whose decisive part (the graviton +− edge) is the same unknown as at 1000.

**In hand:**
* (i)'s leading-Regge layer and (ii) [P / certified];
* the chart-A test at the binding state [NE±];
* the scalar tests of A and B [DV, scalar];
* DG_SCAN's proxies.

**Not in hand:**
* H_mono-A for the graviton beyond J = 5 at more than two gaps, and at J ≥ 26, where no graviton rows exist at two gaps;
* H_mono-B for the graviton at all;
* the uniformity of H_mono-A near x → 1 at large J²/m (K₂ invalid; exact scalar fine to J = 48);
* the seam for 6 ≤ J ≤ J_B;
* Theorem L′'s ε_R and tolerances at Δ′_g.

---

## 6. Plan and effort

### 6.1 Scaling route (recommended), ≈ 3–5 days of agent time + ≈ 40 CPU-h local + ≈ 60 CPU-h external
1. **Certified layer, uniform in Δ_g** (½ day):
   * write up Lemma-C monotonicity [P] (§2.1);
   * extend the LR∞ ball certification past 10⁶ by a tail lemma, or stop at 10⁶ (Theorem L′ already stops at J = 10⁶);
   * optionally certify the hybrid curvature integral by interval quadrature, to lift J_LC.
2. **H_mono-A on the graviton at 2000 and 5000 for every error-budgeted spin** (J = 0, 2, 4–8, 12, 16, 21, 27, 29; parity-odd OO
   J = 4–10): own-mass K₂ rows at the twist edges and at 1.02/1.05Λ.
   * Odd J: `gs_rows.py` as is, ≈ 17 min/row at J = 5, ∝ J^{0.3–0.5} beyond.
   * Even J: an even-branch copy of ng_rows; ≈ 25–40 min/row.
   * ≈ 16 spins × 2 gaps × 1–2 rows ≈ 30–40 CPU-h local (2 processes: ≈ 1 day wall).
   * Then `gs_bind.py` per state (minutes).
   * Outcome: [NE±] at 2000 and 5000 on the same set as 1000, a measured p per spin, and the small-J seam closed by rows at the
     Δ′_g twist edges.
3. **Direct-node checks at J = 12, 16, 21, 27, 29** at 2000 (the 1000 states whose [NE±] needed direct rows): ext batch with Δ_g = 2000
   n1 rows.
   * ≈ 34 nodes × 5 states × 2.4 min × (J-scaling) ≈ 10–15 CPU-h; the n3 window nodes ≈ 20 CPU-h.
   * **≈ 30–40 CPU-h ≈ 2 h on the 24-core machine** (`ext/run_ext.py`, new jobs file, same `hc_points` engine).
4. **Chart B at 2000** (2 days):
   * parametrise Δ_g in `lsw_common` (Dg = 1000), `ag_exact` (Dg), `ag_flat` (al = 998) and `p4w_common` (L = W.L);
   * rerun P2_CONTRACT's O(1/m²) costs (≈ 4 CPU-h) and ag_chi/P4W-χ (≈ 1.2 + 0.7 CPU-h) at the chart-B images of the 47 Theorem W
     states (J′ ≈ 2J) and at the seam states J ≤ 344 near the 2000 twist edge;
   * exact-scalar EDGE-R at those states (`er_scalar`-type, ≈ 1–2 CPU-h).
   * Outcome: Theorem W at 2000 on sampled states, on the same premises.
5. **Write-up:** the transfer proposition (§5) with H_mono-A/B stated as named hypotheses and the tests tabulated; Fig. 1 with filled
   symbols at 2000/5000.

### 6.2 Per-gap repetition (the alternative), per additional gap
* Light side: free (Λ-forms).
* Weights: seconds (grids); ≲ 2.5 CPU-h (continuum fits).
* Flat certificates: free.
* LR∞: done.
* Near-gap error budgets J ≤ 31: ≈ 40 CPU-h local + ≈ 40 CPU-h external (6.1 steps 2–3).
* Theorem W: ≈ 6 CPU-h after parametrisation (step 4). Plus EDGE-R for the graviton via the ext route: ≈ 90 CPU-h at 1000
  (`ext/jobs_batch3.txt`), ×(1.5–2) at 2000 because node cost ∝ J; **≈ 150–180 CPU-h ≈ 7–8 h on 22 workers**. Still blocked by the n3
  C3H bug (EDGE_R finding 4).
* Theorem L′ tolerances: ≈ 1 CPU-h after parametrisation.
* H_LJ windows / large-mass strip: ≈ 10–20 CPU-h.
* Agent time: ≈ 2–3 weeks; CPU ≈ 300–400 h, of which ≈ 250 h is external.
* It gives a second data point at the same grade as 1000, not a range. With §6.1 it gives two computed gaps plus a tested transfer
  to all gaps.

### 6.3 What would make it a proof
Two conditions, neither of which is available:
* a theorem that the bulk-point expansion of the TTO light-ray functional at fixed (J, x) has remainders bounded uniformly in m with
  explicit constants (H_mono-A as a lemma);
* the same uniformly in β (H_mono-B).

The first is plausible from the Laplace structure of the crystal integral (the K₂ series is in 1/Δ with Δ-independent leading
coefficients in c₂ units) [D-sketch]. The edge x → 1 and large J²/m are where it would be hard. This is the same obstruction as AG_U.

---

## 7. Reproduce
```
cd ads4/code/gapscale && bash run_all.sh      # ~1.6 CPU-h, <= 2 processes; outputs in cache/
```
Individual steps and timings are listed in `run_all.sh`.
* `gs_bind.py` must run after `gs_rows.py` and `gs_weights.py` (it reads `cache/rows/` and `cache/weights_2000.npz`).
* `gs_cover.py` reads the LR-ball logs.
