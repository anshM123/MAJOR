# REVIEW_SLIVER_L5: adversarial referee report on ANALYTIC_SLIVER (L ≥ 5)

*26 Nov 2026. Referee agent. No existing file edited; `prl/` not touched (only read: `prl/data/nf02_psi.json` via `an_common.nf02`);
nothing committed. Review scripts and logs in `ads4/code/review_l5/` (a private copy of `code/analytic/*.py` plus `rv_*.py`; reruns write
only to `code/review_l5/cache/`; the original `rt_loc3.json` cache was not used, so every Loc′ integral was recomputed from scratch).*

## 0. Verdict

**No error that breaks the L ≥ 5 certificate was found.** Every certified run reproduces bit-for-bit. The [5, 5.5] cell survives with
ρ = 0.94683, a 5.3 % margin. The state relations are proved below for every state with L ≥ 4.5, not only spot-checked.

Findings:
* one mis-attribution in the text (E1);
* four rigor gaps (G1–G4). Each is numerically negligible, by at least 10 orders of magnitude against the 5 % margin, but should be
  written out.

| item | verdict |
|---|---|
| Lemma Loc′ (§7.3) | **holds** (all inequality directions and cell monotonicities checked; constants verified exactly; gap G2a at ν > 1300) |
| [5, 5.5] cell | **survives**, ρ = 0.94683 reproduced; no endpoint error; omitted terms ≲ 1e-10·F (G1, G2) |
| state relations | **proved** (§3), valid for all a ≥ 998, b ∈ {0, 4}, d ≤ e^{−4.5} |
| Theorem E′ | holds (proof re-derived step by step) |
| Lemma X | holds numerically at all 126 (b = 0) + 8 (b = 4) test points at L ≈ 5, 5.5 (bound/exact ≥ 13 … 10³⁰); **proof not re-derived** (G4) |
| helicity bound §5.4 | holds (identity checked to 1e-85; bound checked against exact s = 9 values at L = 5) |
| splitting lemma, Theorem R | hold; the bookkeeping omits the θ-cutoff correction (G1) |
| Lemma NB (explicit d) | holds up to a dropped lower-order second-order term (G3, relative ≤ 1e-7) |
| Lemma V | holds (inhomogeneity R(β) = βc(βc + 1 + sd)/d² re-derived exactly) |
| imported constants | used in the right direction: F_cert is a cell minimum; the dressing constant 40 is an upper bound (true constant ~1e-5); sup\|f\| ≤ 21.34 is an upper bound; K₀ ≤ 1.6 (exact 1.0003) |

## 1. Reruns (from source, fresh caches)

| run | cached ρ | rerun ρ | status |
|---|---|---|---|
| `rt_cert.py 5 5.5 0.5 2 1500` | 0.94683 | 0.94683 | reproduced (all lines identical) |
| `rt_cert.py 5.5 6.5 0.5 2 1500` | 0.17942, 0.069807 | same | reproduced (file identical up to timings) |
| `rt_cert.py 6.5 8.0 0.5 2 1500` | 0.041726, 0.031817, 0.027767 | same | reproduced (identical) |
| `rt_cert.py 8 11` + `11 13.5` + `13.5 16` (`0.5 2`) | 16 cells, max 0.046384 | same | reproduced (identical to `rt_cert_final_8_16.txt`) |
| Lemma Loc′ integrals lI₂, lI₃ (6 cells, 5–8) | `rt_loc3.json` | recomputed | identical (e.g. [5, 5.5]: −894.399, −990.651) |
| `rt_rmass2.py` (EPSW panels → unit cells) | `rt_rmass2.json` | regenerated | identical (163 cells, total 3.2958e-9) |

Not rerun:
* the [4.5, 5] cell, which FAILs anyway;
* `rt_rem.py` (TAIL Lemma Loc constants, used only on [8, 16], where E_ρ ≤ 1.5e-26);
* the TAIL+E′ scripts;
* EPSW_CERT and LRINF_TIGHT themselves (imported).

## 2. Lemma Loc′ (§7.3): line-by-line

* **Representation.**
  - The identity ℒ(z^k) = E[B^k]·C(N−k, a)/C(N, a) = U_k was re-derived (B ~ Beta(n+s, a)).
  - The bound |ℒ(R)| ≤ K₀·E_B E_θ^W|R| holds with the normalised θ-measure W/∫W.
  - The four regions partition {B ∈ [0, 1], θ}.
* **Radii.**
  - Upper: |1−ξ| ≤ (1−Bρ) + |θ| ≤ b_{k+1} + 0.51d + θ₁.
  - Lower on tiers: |1−ξ| ≥ 1 − B > b_k.
  - Lower on Th: |1−ξ| ≥ 1 − ρ ≥ 0.49 d_lo.
  - Dressing circles of radius (1−ρ)/2 ≤ 0.255 d_max stay in |ζ| < 1. The code's lower radius 0.49 d_lo − 0.255 d_max is valid.
* **Directions over the L-cell.**
  - Upper radii, weights (1.02d/b)^998, om = 0.51d in frac, and the d³ prefactor all use d_max = e^{−L_a}.
  - Lower radii and the (1−ρ)^{−1} factor in frac use d_lo.
  - In q(θ), ρ is replaced by 1 − 0.51d ≤ ρ, which is the right direction since q increases with ρ.
  - frac is evaluated at a = 998; it is decreasing in a, as checked exactly at a = 998 … 10⁵.
* **c_j bounds.**
  - Numerator at ν_hi, denominator at ν_lo.
  - |B| is decreasing (|B|² = coth(πν)/(4πν)), so B_abs_real(ν_lo) is correct.
  - |νB| is increasing, so nuB_abs(ν_hi) is correct.
  - The tail ratio ((j+½)² + ν²)/((j+1)√((j+1)² + 4ν²)) ≤ (1 + ν/(J+1))² was re-proved, and so was the complex-w version.
* **B-pole (small ν).**
  - B(w) = E(w)/(2iw) by Γ(1+x) = xΓ(x), so Σ_±BΦ = (h(ν) − h(−ν))/(2iν) with h = EΦ.
  - Cauchy with r′ = 0.05 stays inside the box |Re w| ≤ 0.45, |Im w| ≤ 0.05.
  - h is analytic there: Γ(1+2iw) has poles only at Im w ≥ ½, and c_j has poles at Im w ≤ −½.
  - |(1−ξ)^{−iw}| ≤ |1−ξ|^{Im w} e^{π|Re w|/2} is handled by lo^{−y}.
  - Correct.
* **Exact checks** (`rv_state.py`, mpmath):
  - **K₀.** K₀ = ρ^{−n}·₂F₁((a+1)/2, (a+1)/2; 1; ρ²)/C(N, a) = 1.00033 (a = 998) … 1.000003 (a = 10⁵), for L ∈ {4.5, 5, 5.5, 8, 16} and b ∈ {0, 4}.
    The claimed analytic bound is ≈ (π/2)(1 + 3·10⁻⁴) < 1.6.
  - **Lower bound on I_W.** I_W ≥ e^{−0.51}(1−ρ)^{−a}/(π√a) holds with a factor 2.09 to spare.
  - **θ-fraction at θ = 0.02.** Exact −1314.10 vs bound −1301.47 (L = 4.5, a = 998), and the bound holds at all 40 tested states.
  - **P(B < 1−b).** Computed exactly via the binomial identity: log P = −8575 vs the bound −1758 (L = 5, b = 0.04), and the bound holds at
    all 48 states (a up to 2·10⁴, b up to 0.95). The step (a+i)/(N+s+i) ≤ 2a/N needs a ≥ 997.
* **Growth constant.**
  - (A+B)/den at |β| = hypot(1503, 1151) + 2 is 1.9459 (L = 4.5) and 1.9376 (L = 5).
  - The constant is increasing in |β| and d, so the assert at d_max is valid.
  - The near-real (Δ) tail also needs F/den ≤ 1, so that M′ + 1 doubles. Its value is 0.93: this holds, but the hypothesis is implicit.
    State it.
* **Gap G2a.**
  - `rt_loc3.integrals` adds nothing for ν > 1300. It only asserts b(2000) − b(1300) < 8000 and tail(1300) + b(2000) < −30000.
  - There is no statement for ν > 2001, where the bound grows like e^{πν/2}.
  - True (CONTOUR tail: tail(1150) = e^{−9103}, tail(1300) = e^{−45456}, Gaussian beyond), but not written as a proof.

## 3. State relations: proof (L ≥ 4.5, i.e. d ≤ e^{−4.5} = 0.011109)

Exactly, D − C₀ = (N−n)(N+n+s) = a(2N − a + s), so **d = a(2N − a + s)/(N(N+s))** and, with x = a/N,
  a/(Nd) = (N+s)/(2N − a + s),  D d²/a² = (2 − x + s/N)²/(1 + s/N).

1. **Lower bound.** a/(Nd) ≥ ½ ⇔ a + s ≥ 0. This is always true, so a/N ≥ d/2 ≥ 0.49d.
2. **Upper bound.**
   * Step 1: x(2 − x) ≤ x(2 − x + s/N) = d(1 + s/N). Since x < 1 and s/N ≤ 9/1998, this gives x < 1.0046d ≤ 0.01116.
   * Step 2: then N ≥ a/x ≥ 8.9·10⁴, so s/N ≤ 1.01·10⁻⁴.
   * Step 3: so x ≤ d(1 + s/N)/(2 − x) ≤ 1.0001d/1.9888 ≤ **0.5029d ≤ 0.51d**.
3. **T relation.**
   * D d²/a² ≥ (2 − x)²/(1 + s/N) ≥ (2 − 0.5029·0.011109)²/1.0001 = **3.977 ≥ 3.9**.
   * T = D + b(b+1) ≥ D.
   * Also N/a ≥ 179 ≥ 40, which is rt_pp's standing hypothesis "N ≥ 40a".

These are uniform in a ≥ 998 and J, and need nothing beyond d ≤ e^{−4.5}.

An adversarial scan in exact rationals agrees: b ∈ {0, 4}, a ∈ {998, …, 10¹⁰} plus 120 random values, minimal n at L ∈ {4.5, …, 20}.
Worst cases: a/(Nd) ≤ 0.50141 and Dd²/a² ≥ 3.97775. The relations first fail near L ≈ 2.5 (a/(Nd) = 0.5108) and L ≈ 3.0
(Dd²/a² = 3.8998), far below 4.5.

## 4. The [5, 5.5] cell and the turning point

* **Rerun.** ρ = 0.94683 is reproduced exactly; residual 0.947; ψ₂ pieces ≤ 1.4e-5; E_ρ (Loc′) = e^{−894} relative.
* **Cell endpoints.**
  - Every piece is evaluated at d_max = e^{−5}, and each is increasing in d:
    - d³·d^{j−Re α} carries a positive exponent on every contour piece;
    - rec_step's A, B and den, the NB constant g, the dressing d²/(3.9a²), and the helicity d₀ all increase with d.
  - F_cert is the minimum over β ∈ [6.383, 6.884].
  - No piece is taken at the wrong endpoint.
* **E1 (mis-attribution in §7.4, §9).**
  - The text says the unit cell ν ∈ [1031, 1032] contributes 0.94·F at L = 5 (and 16 at L = 4.5). In fact 0.94 is the **whole [898, ∞)
    band**: the log line reads "[898,oo) worst cell 1031: 9.4e-01", which is the band sum together with the arg-max cell.
  - `rv_cell1031.py` gives the per-cell values:
    - [1031, 1032.1] alone: 0.221;
    - [1030, 1031.1]: 0.125;
    - [1029, 1030.1]: 0.104;
    - …;
    - [1020, 1021.1]: 0.014;
    - [1032, 1033.1]: 0.053.
  - Σ over [1020, 1046] = 0.848. The binding mass is spread over ≈ [1000, 1033].
  - So the proposed cure ("finer certified residual masses near ν ≈ 1031") must cover about 30 cells, not one. Soundness is unaffected.
* **Lemma X at the binding frequencies (DV, `rv_exactS.py`).**
  - Exact |Ŝ| comes from the spin recurrence at n ≈ 3–5·10⁵ (L = 5.00, 5.50). It is stable between 120 and 240 bits.
  - b = 0 at τ = 1031: exact |Ŝ(2.5 + 1031i)| = 1817 (a = 998), 66 (a = 1015), 14.7 (a = 1031), 1.7 (a = 1500), against the bound
    9.74e7.
  - b = 0 at τ = 1150: exact 4.6e23 against the bound 7.1e30.
  - b = 4 at τ = 1031: exact 1842 against PP = 1.17e8.
  - The bound holds at every point, with exact/bound ≤ 7.8e-2 (at τ = 948, Theorem E′), ≤ 5e-4 at τ = 1000 and ≤ 2.5e-5 for τ ≥ 1031.
  - So the certified 0.947 overstates the true residual by about 10⁴. That is loose, but not wrong.
* **Omitted terms at L = 5.**
  - G1 (below): ≲ 1e-10·F.
  - G2: e^{−7000} or less.
  - The ψ₃ terms are present.
  - Float accumulation is covered by the 1 + 10⁻⁹ factor.
  - **The cell survives.**

## 5. Gaps (none affects any ρ by more than 1e-10)

* **G1: θ-cutoff correction missing from E_r.**
  - H is built from the realised f_θ (TAIL_AUDIT §2.8: "CONTOUR's f is f_θ itself").
  - EPSW's panel masses (`rt_rmass2.py`, EPSW_BANDS) are for the uncut f; `ew_common` has no θ. FULLPOS3 §4.2 and TAIL R.2 add δE_θ
    explicitly; Route T′ does not.
  - Add d³(2/Λ)·2∫|f_θ − f||K_J|. With ∫|f_θ − f| ≤ 1.66e-23 and sup_{ν≥Λ+20}|f_θ − f| ≤ 3.65e-24 (`tailcert/cache/tc_cutoff.txt`)
    against a K-weight of ≤ 10¹¹·F per unit mass at ν ≤ 1033, the term is ≲ 1e-11·F.
* **G2: far frequencies.**
  - (a) Loc′ above ν = 1300 (see §2).
  - (b) E_r on [1150, 3000] uses the "crude" bound |K_J| ≤ 10³e^{πν/2}ν³d^{−5/2}:
    - The reflection inequality |Ŝ(β)| ≤ Γ(σ)|Γ(β)|⁻¹Ŝ(σ) is correct (U_k ≥ 0, |(β)_k| ≤ (σ)_kΓ(σ)/|Γ(β)|), but only for Re β > 0,
      i.e. j ≤ 2.
    - The j ≥ 3 terms are not controlled.
    - Beyond |β| ≈ 1900 the step constant exceeds 2: it is ≈ 3.3 at ν = 3000.
    - My estimate is e^{39} vs the claimed e^{43} at ν = 3000, L = 5, so the bound is plausible.
  - (c) Nothing is added for ν > 3000.
  - All of this is multiplied by tail(1150) = e^{−9103}, so it is irrelevant numerically. Replace it with a proved e^{cν} bound.
* **G3: Lemma NB second order.**
  - The code keeps (4/15)g²(σ)₅.
  - E[S₂(K)] also has factorial-moment terms ∝ (σ)₄, (σ)₃, … (S₂(k) = Σj²(j+s)², and Σ_{j<k}j⁴ = (k)₅/5 + 1.5(k)₄ + …).
  - These contribute ≈ (4/3)(1.5 + s/2)g²d(σ)₄, which is ≤ 30d/(σ+4) relative to a term of size ~1e-13.
  - Irrelevant, but the displayed formula is not literally proved.
* **G4: Lemma X.**
  - Only the τ-monotonicity of every term in `rt_pp.cross` was audited. Each term uses the worst of t_lo/t_hi (q = |α(t_hi)|/t_lo,
    G, W and y_max at t_hi; δ = c·t_lo^{−2/3} is a free parameter); uniformity in a holds because G_a is decreasing in a.
  - The phase I/II/III estimates themselves were not re-derived. They are supported by the DV in §4 and by `crossing_dv.txt`.
  - This is the one remaining unreviewed analytic input on which the binding cell rests.

## 6. Other lemmas (hypotheses checked where used)

* **Theorem E′.**
  - Checked:
    - Re asinh w ≤ u/√(1−y²);
    - log(d_k/d_{k+1}) ≥ ε_k/2;
    - ε_k = 2/(N−k) − 2/(N+k+1), which is convex and satisfies Δε ≥ ε²/2;
    - the Riemann-sum direction, and the n-independence.
  - Used only for τ ≤ 0.95·998 and σ ∈ [0, σ_max]; the bound is monotone increasing in σ and τ and decreasing in a.
* **Splitting lemma.**
  - κ_J is analytic on the rectangle: S(β) is a terminating ₃F₂, hence entire in β; B is regular for Re ν ≥ ¼; c_j has its poles in
    Im ν < 0.
  - The orientation is correct.
  - The top line Y = 4.25 misses the B poles at 4i and 4.5i by ≥ 0.25.
* **Theorem R.**
  - The main term (FULLPOS3 §4.1(c)) is d³[G^ψ_A(c) + τσ_bG^ψ_D(c)], without the dressing and without r. Route T′ is consistent with this:
    - Mdiff includes the full dressing;
    - E_r bounds |𝒱_bℒ| (|Ŝ| + 1 + dressing);
    - E_ρ uses all of f.
  - Apart from G1, no term is missing.
* **Helicity bound.**
  - d/d₀ = (D+20)/D and U₄ ≥ (1−d₀)⁴(1 − 10⁻⁶) are fine.
  - The 2^i step growth needs (A+B)/den ≤ 2 at |β| ≤ 1155: the value is 1.35, so it holds.

## 7. Files
* `code/review_l5/rv_state.py`: state relations; K₀, I_W, θ-fraction; Beta tail.
* `code/review_l5/rv_exactS.py`: exact |Ŝ| (s = 1 and s = 9) vs PP0/PP.
* `code/review_l5/rv_cell1031.py`: per-cell residual split at L = 5.
* `code/review_l5/rv_rmass2_rerun.py`: regeneration of the residual masses.
* `code/review_l5/logs/*`: all rerun logs.
* `code/review_l5/cache/`: rerun outputs.
