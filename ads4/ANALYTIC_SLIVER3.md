# ANALYTIC_SLIVER3: below L = 3.65 — the remainder is dead, the hump is what is left (Δ_g = 1000)

*27 Nov 2026. Research agent. Running document; unreviewed, not committed. No existing file edited; nothing under `prl/` touched.*
*New code: `ads4/code/analytic3/` (imports `code/analytic/`, `code/analytic2/` read-only; outputs in `cache/`, logs in `logs/`).*
*Grades: **[P]** proof written here; **[certified]** Arb/exact scalar constants, positive bounds accumulated in floats with a final
(1 + 10⁻⁹) (logs: +10⁻⁶); **[DV]** sampled numerics; **[imported]** taken from another document.*
*"AS §x" = ANALYTIC_SLIVER.md, "AS2 §x" = ANALYTIC_SLIVER2.md.*

---

## 0. Verdict

__VERDICT__

---

## 1. What changed against AS2 §8, in one table

| AS2 obstruction | status here |
|---|---|
| 1. Loc″ remainder (log E_ρ = +26 on [3.58, 3.65]) | **removed down to L ≈ 1.7** by Lemma Loc‴ (§2): log E_ρ2 = −2063 on [3.55, 3.65], −1453 on [3, 3.1], −973 on [2.5, 2.6], −437 on [2, 2.1]; it fails only near L ≈ 1.6 (+155 on [1.5, 1.6], J = 1500) |
| 2. contour right edge (0.32·F at 3.65), top-line j-tail assertion (fails ≈ 3.58) | **fixed** (§3): right edge by Lemma F + Lemma B over a-boxes, 1.1·10⁻³·F on [3.5, 3.65]; top line with 160 explicit j and a general growth constant, 1.0·10⁻²·F |
| 3. contour J = 1500 runs out (growth 1.983 vs 1.99) | **dissolved**: Loc‴ needs no REC hypothesis; the near-real tails need (A+B)/den ≤ 2 only up to \|β\| ≈ J + 13, where it is 1.62 (asserted); the oscillatory tails use the actual G |
| 4. termwise e^{νd/2} hump at large twist; state relations near 2.5–3 | **this is now the only obstruction** of the analytic route (§4, §6): it binds at L ≈ 3.4–3.5 in the real-axis residual. State relations are no longer used in the form 0.49d ≤ a/N ≤ 0.51d anywhere in the new code (exact t-ranges); D d² ≥ 3.9a² is used for L ≥ 3.01 and replaced by its exact form below |

---

## 2. Lemma Loc‴: the j-truncation remainder by a product-measure cell decomposition [P + certified arithmetic] (`loc5.py`, `locrun3.py`)

**Why Loc′/Loc″ failed.** Their regions were one-dimensional: the B-tiers used the 998-th moment, P(1−B ≥ b) ≤ (1.02d/b)^{998}, and the θ-tail
region |θ| > θ₁ was bounded with |1−ξ| up to 2 (so |T_J| ≲ Σ|c_j|2^j ≈ e^{1500}) and with the worst phase e^{πν/2}. Both weights are
off by thousands in the exponent: at L = 3.6, a = 998, the exact log P(1−B ≥ 0.069) is −2521 (the 998-th moment Markov bound gives ≈ −430),
and the θ-weight decays like (t/2 sin(θ/2))^{a+1}, i.e. ≈ e^{−4000} at θ ≈ 0.85, where |1−ξ| is still ≤ 1.

**Setting** (TAIL_CERT §4; re-derived in REVIEW_SLIVER_L5 §2). For G analytic in |z| < 1,
ℒ(G) = E_B[(2πi)⁻¹∮_{|z|=ρ}(1−z)^{−a−1}G(Bz)z^{−n−1}dz]/C(N,a), B ~ Beta(n+s, a), ρ = n/N, t := 1 − ρ = a/N; hence
|ℒ(G)| ≤ K₀·E_B E^W_θ|G(Bρe^{iθ})| with E^W the normalised weight W(θ) = |1−ρe^{iθ}|^{−a−1} and K₀ = ρ^{−n}I_W/C(N,a), I_W = (2π)⁻¹∫W.
R_J(ξ;ν) = F_{p,ν}(ξ) − T_J(ξ), T_J = Σ_±B(±ν)Σ_{j<J}c_j(±ν)(1−ξ)^{j−e₀∓iν}, e₀ = p − ½; for |1−ξ| < 1 also R_J = Σ_±B(±ν)Σ_{j≥J}(…).

**Lemma K (K₀).** For every state (a ≥ 998, any n, s, real N): K₀ ≤ (π/2)√(1 + 2/a)e^{1/(6a)} ≤ 1.5727.
*Proof.* sin(θ/2) ≥ θ/π on [0, π] gives |1−ρe^{iθ}|² ≥ t² + (4ρ/π²)θ², so I_W ≤ (1/π)∫₀^∞(t² + 4ρθ²/π²)^{−(a+1)/2}dθ =
t^{−a}(2√ρ)⁻¹·√πΓ(a/2)/(2Γ((a+1)/2)). Robbins' Stirling bounds give C(N,a)ρ^n ≥ t^{−a}e^{−1/(12a)−1/(12n)}/√(2πaρ) (the (N/a)^a(n/N)^n
factors combine exactly to t^{−a}). Gautschi: Γ(x)/Γ(x+½) = x⁻¹Γ(x+1)/Γ(x+½) < x⁻¹√(x+1), x = a/2. Multiply. ∎
(DV, `check_loc5.py`: exact K₀ = 1.00008–1.00034 at L ∈ {3.6, 3, 2}, a ∈ {998, 1500, 4000}, b ∈ {0, 4}.)

**Lemma W (θ-weight).** P_W(|θ| ≥ θ₀) ≤ √(2π(a+1))·t⁻¹·(t/q(θ₀))^{a+1}/erf(X√((a+1)/2)), q(θ)² = t² + 4ρ sin²(θ/2), X = π√ρ/t.
*Proof.* q increases on [0, π], so ∫_{θ₀≤|θ|≤π}W ≤ 2π q(θ₀)^{−a−1}. Lower bound: q² ≤ t² + ρθ² and 1 + x² ≤ e^{x²} give
I_W ≥ (1/π)t^{−a}ρ^{−½}∫₀^X e^{−(a+1)x²/2}dx. ∎ Uniformity: t/q is increasing in t (use t_hi), t⁻¹ at t_lo, and √(a+1)r^{a+1} decreases
in a when r < 0.999 (else the weight 1 is used).

**Lemma C (Chernoff B-weight).** For every integer m, P(1 − B ≥ u) ≤ (a)_m/((N+s)_m u^m) (1−B ~ Beta(a, n+s)). With
(a+i)/(N+s+i) ≤ (1+x)/(K+x), x = i/a, K := 1/t_hi ≤ N/a, the summand f(x) = log((1+x)/((K+x)u)) is increasing, so
Σ_{i<m}f(i/a) ≤ a∫₀^{m/a}f =: aF(m/a), F(X) = (1+X)log(1+X) − (K+X)log(K+X) + K log K − X log u. Taking m = ⌊aX₀⌋ (X₀ = the zero of f):
log P ≤ 998·max(F(X₀ − 1/998), F(X₀)) whenever this is ≤ 0 (F is decreasing on [0, X₀], negative, and aF is then decreasing in a). ∎
(DV: within 1–2 % of the exact log-tail, e.g. −2502 vs exact −2521 at L = 3.6, u = 0.069, a = 998.)

**Lemma D (dressing).** |𝒱_bG(ξ)| = |(1 + τ𝒟_b)G(ξ)| ≤ (1 + 2.6·10⁻⁵)·sup_{|ζ−ξ|≤t/2}|G(ζ)| for |ξ| ≤ ρ.
*Proof.* 𝒟_b = z(1−z)∂²/2 + [(7/2 − 2b)z + ½]∂ + 9/2 − b(b−1) (FULLPOS3 §4.1), Cauchy on the circle of radius r = t/2 (inside |ζ| < 1):
|τ𝒟_bG| ≤ τ(2/r² + 5/r + 7.5)M with τ = 1/T ≤ 1/N², so τ/r² ≤ 4/a², τ/r ≤ 2/a², τ ≤ 1/a²: total ≤ 25.5/a² ≤ 2.6·10⁻⁵. ∎
(TAIL's constant 40 is replaced by 2.6·10⁻⁵; the δ-enlargement of the regions by t/2 stays.)

**Lemma Loc‴.** Partition (u, |θ|) ∈ [0,1] × [0,π] (u = 1 − B) into cells Q = [u₀,u₁) × [θ₀,θ₁). For ζ in the enlarged cell
Q^δ = {ζ : |ζ − (1−u)ρe^{iθ}| ≤ t_hi/2, (u, θ) ∈ Q, t ∈ [t_lo, t_hi]} let lo ≤ |1−ζ| ≤ hi, |arg(1−ζ)| ≤ A and |Im η| ≤ E (η = 2artanh√ζ,
cosh η = (1+ζ)/(1−ζ)) be Arb enclosures. Then for real ν ∈ [ν_lo, ν_hi] (ν_lo ≥ 0.4)
  |𝒱_bℒ(R_J)| ≤ 1.5727·(1 + 2.6·10⁻⁵)·Σ_Q P(1−B ≥ u₀)·P_W(|θ| ≥ θ₀)·min(T_Q, FT_Q),
  T_Q = 2B̄e^{ν_hi A}|c_J|hi^{J−e₀}/(1 − R_t·hi)  (only if R_t·hi < 1),
  FT_Q = (π/2)e^{ν_hi E}lo^{−p} + 2B̄e^{ν_hi A}Σ_{j<J}|c_j|max(lo, hi)^{j−e₀}  (max(lo,hi)^{j−e₀} = lo^{j−e₀} for j < e₀, hi^{j−e₀} else),
with B̄ = |B(ν_lo)| (|B|² = coth(πν)/(4πν) decreasing), |c_j| ≤ Π_{i<j}((i+½)² + ν_hi²)/(√((i+1)² + 4ν_lo²)(i+1)) (Arb), and
R_t ≥ sup_{j≥J}|c_{j+1}/c_j|: R_t = 1 + ν_hi²/(J+1)², or 1 under the tail-ratio lemma of AS2 §5.
*Proof.* B and θ are independent, so the probability of a cell is P_B(u ∈ [u₀,u₁))·P_W(|θ| ∈ [θ₀,θ₁)) ≤ P(1−B ≥ u₀)·P_W(|θ| ≥ θ₀) (Lemmas C, W);
|R(ξ̄)| = |R(ξ)| for real ν (B(−ν) = B̄(ν), c_j(−ν) = c̄_j(ν), F real on the real axis), so |θ| suffices; Lemma D enlarges by t/2.
On a cell: |(1−ζ)^{j−e₀∓iν}| = |1−ζ|^{j−e₀}e^{±ν arg(1−ζ)} ≤ |1−ζ|^{j−e₀}e^{νA} (both branches; |arg(1−ζ)| < π/2 since |ζ| < 1);
|F_{p,ν}(ζ)| = |1−ζ|^{−p}|P_{−½+iν}(cosh η)| ≤ (π/2)e^{ν|Im η|}|1−ζ|^{−p} (FULLPOS3 Lemma B; Re η > 0 and |Im η| < π/2 because
(1+√ζ)/(1−√ζ) maps the disc to the right half-plane; on the boundary Re η = 0 the bound extends by continuity). The tail series converges
for |1−ζ| < 1 with ratio ≤ R_t·hi. K₀ by Lemma K. ∎
*ν ∈ [0, 0.4]* (pole of B at 0): exactly as Loc′ (AS §7.3 item 7): Σ_±B(±ν)Φ(±ν) = (2iν)⁻¹∫_{−ν}^{ν}(EΦ)′, ≤ (1/r′)sup_{Wbox}|EΦ| with
r′ = 0.05, |(1−ζ)^{−iw}| ≤ e^{|Re w|A}max(lo^{−r′}, hi^{r′}), c_j(w) bounded with |a₀+i|² ≤ (i+½+|y|)² + x², |b₀+i| ≥ i+1−2|y|
(the ratio is then < 1 for every j); F itself has no pole and is bounded directly.
*Uniformity.* The geometry depends on the state only through t = a/N ∈ [t_lo, t_hi], computed from the exact relation
t² − (2+e)t + d(1+e) = 0 (e = s/N ≤ 9·0.7/998; t increasing in d and in e). Chernoff uses K = 1/t_hi; Lemma W as stated; the worst a is 998.
Nothing else (no REC, no D d² relation, no 0.51d) enters.
*Cells.* u- and θ-breaks at 0, t_hi/4·1.25^k (≈ 830 cells at L ≈ 3.5); geometry in Arb (sector boxes r·e^{iθ} plus the t_hi/2 square), done once
per L-cell; per ν-cell only |c_j| and B change (vectorised log-sum-exp in floats, +10⁻⁶ in the log).

**Relation to the coordinator's phase-keeping lemma** (message of 27 Nov: Φ(z) = |arg(1−z)| + 2|arg(1+√z)| ≥ |Im β|, and a global
saddle bound (a+1)log(ε/|1−ρe^{iθ}|) + νΦ ≤ (a+1)Λ*(ν/(a+1)) + Λ₂, ≈ e^{464} instead of e^{1616} at ν = 1029). Loc‴ is the cellwise
version of the same idea, with three differences: Im η is enclosed directly (sharper than Φ), the B-direction is weighted by the optimal
Chernoff moment instead of a fixed B-range, and the remainder's own size (|c_J|hi^J on the central cells, the T_J sum on the outer cells) is
optimised jointly per cell, as the coordinator's caveat (b) asks. The coordinator's caveat (a) (the δ-enlargement destroys the phase
near the saddle) does not bite here: on the cells where A ≈ π/2 because of the enlargement, the tail series has hi ≤ 3t and hi^J ≈ e^{−6000}.

**Numbers [certified]** (J = 1500; log of the bound of |𝒱_bℒ(R_J)| per unit ν-cell, L-cell [L, L + 0.15]):
| L | ν = 1 | 300 | 600 | 900 | 1029 | 1149 |
|---|---|---|---|---|---|---|
| 3.5 | −3639 | −3406 | −3005 | −2471 | −2226 | −1999 |
| 3.0 (cell [3, 3.1]) | −3132 | −2899 | −2411 | −1882 | −1654 | −1443 |
Old (Loc″, AS2 §5): +197 at L = 3.5, ν = 1029 (J = 1500); +603 at L = 3.0.

Integrated (`locrun3.py`, FSUP = 21.34 on [0, 1150), the k = 3 piece with p = 4 undressed; same prefactors as AS2's rt2_loc):
| L-cell | log E_ρ2 | log E_ρ3 |
|---|---|---|
| [3.55, 3.65] | −2063 | −2368 |
| [3.5, 3.65] | −2011 | −2321 |
| [3.0, 3.1] | −1453 | −1746 |
| [2.5, 2.6] | −973 | −1289 |
| [2.0, 2.1] | −437 | −834 |
| [1.5, 1.6] | **+155** | −365 |
(F_cert ≈ e^{−13}.) **The remainder is no longer an obstruction anywhere above L ≈ 1.7.** (At L = 1.5 the central cells' e^{νA}|c_J|hi^J
binds; a larger J would push this further, but the j-sums of §4 fail long before.)

DV of the weights (`check_loc5.py`, exact Beta tails by the binomial duality, exact θ-masses by quadrature): every bound lies above the
exact value, by 0.5–2 % of the exponent, at L ∈ {3.6, 3, 2}, a ∈ {998, 1500, 4000}, b ∈ {0, 4}.

---

## 3. Contour pieces (`contour3.py`) [P + certified arithmetic]

* **Bottom, left** (near-real exponents): `rt_cert.piece_bottom/left` unchanged (Lemma NB, Lemma V). Their geometric j-tails assume REC growth
  ≤ 2 and inhomogeneity F/den ≤ 1 up to |β| ≤ J + 4 + Y + 1: asserted at d_hi — values 1.624 and 0.594 at L = 3.55 (in AS2 the assertion was
  made at |β| = hypot(J+3, 1151), a real-axis requirement of Loc′ that Loc‴ no longer has).
* **Top line** ν = t + iY, Re α = ¼ − k₀ (k₀ = 4, Y = 6.25 for ψ₂): base Ŝ(¼ + it), Ŝ(5/4 + it) by `rt_pp.PP` (uniform over a ≥ 998, both
  channels), absolute REC chain (`rt_cert.chain`), **160 explicit j-terms** and the geometric tail with the *actual* growth constant
  G = max(2, (A+B)/den at |β|_max) and Mdiff_j ≤ (1 + τ(A₀ + A₁/d))max(|Ŝ_j|, |Ŝ_{j−1}|) + 1; rb = ratio·d·G < 0.9 asserted.
* **Right edge** ν = Λ + iy (y-cells of ¼): for each a-box of `kernel3.a_ranges` (998 ≤ a ≤ 4500; unit boxes up to Λ + 30, width 2 up to Λ/0.85,
  then 3 %), base powers at the top of the shifted σ-range by min(Theorem E′ / Lemma X, **Lemma F with a σ-box**) and the helicity bound,
  then **Lemma B** (multi-frame) down to Re β_j, same j-sum. Max over a-boxes per y-cell.

| L-cell | ψ₂ bottom | left | top | right | ψ₃ (all) |
|---|---|---|---|---|---|
| [3.5, 3.65] | 1.06e-5 | 1.20e-5 | 1.03e-2 | 1.09e-3 | 3.2e-5 |
| [3.55, 3.65] | __C1__ |
(AS2 at [3.65, 3.75]: right 0.320, top 3.5e-4 — the right edge drops by ×300.)

---

## 4. The real-axis residual: what binds, measured [certified bounds vs DV exact]

`kernel3.py` = AS2's per-cell kernel bound (Theorem E′, Lemma X, Lemma F, helicity, Lemma B, REC tail) with two changes:
* **twist range 998 ≤ a ≤ A* = 4500** — every state with a ≥ a*(L) is covered by the high-twist cell certificate (ANGULAR_ROUTE §73:
  a* = 4500 for β_c ∈ [3.04, 8.24], i.e. L ≈ 2.3–7.3) [**imported**; the analytic certificates of AS, AS2 covered every a];
* a finer a-cover (3 % instead of 15 % beyond the turning region; ×2–6 at ν ≈ 800–900) and, in the turning band ν ∈ [995, 1059], a-boxes of ¼ and
  ν-subcells of ½ (×1.5–2.5 there; Lemma F is very sensitive to the box: 5.7e4 → 2.4e4 → 1.5e4 for boxes 1×1.1, ½×0.55, ¼×0.3 at ν = 1031).

**Where the bound loses** (`probe_loss.py`, `dv_resum.py`; exact Ŝ by the dual-Hahn n-recurrence plus the exact β-recurrence):
| state (L, ν, a) | true K_J | exact termwise Σ\|·\| (hump) | Lemma B with exact M₀, M₁ | certified (box) |
|---|---|---|---|---|
| (3.2, 900, 1061), κ = 0.85 | 150 | 1.1e8 | 4.1e8 | 4.6e9 (a-box 1) / 4.8e10 (old cover) |
| (3.2, 900, 1150), κ = 0.78 | 117 | 2.1e8 | 5.9e8 | 6.9e9 |
| (3.5, 1031, 998), κ = 1.03 | 2.7e5 | 7.4e7 | 1.1e9 | 1.1e11 |
| (3.0, 800, 2000), κ = 0.40 | 32 | 4.3e9 | | |
| (2.5, 800, 1500), κ = 0.53 | 1.3 (×\|B d^{−α}\|) | 2.6e13 (×\|B d^{−α}\|) | | |
So at moderate twist the per-power and chain bounds lose e^{3–5} (Theorem E′ ×10 at κ ≈ 0.85; Lemma F ×80 at the turning point; Lemma B
×2–15), and **the hump itself (true → exact termwise) loses e^{13–30}**. The hump grows like e^{r₊νd/2}; with perfect per-power bounds the
termwise residual would still fail near L ≈ 3.3.

Residual budget estimates (`probe_resid3.py`, sampled ν-cells, interpolated; DV of certified per-cell bounds): [3.5, 3.65]: ≈ 0.7–1.2·F
(turning band ν ∈ [1000, 1035] ≈ 0.5 of it, ν ∈ [800, 1000] ≈ 0.15); [3.2, 3.35]: ≈ 90·F (old cover). Hence the certificate below is run on
the narrower cell [3.55, 3.65].

---

## 5. Certificate

__CERT__

---

## 6. The hump: the one remaining obstruction of the analytic route, and how to remove it (design + DV)

**Mechanism.** K_J = 2Re[B d^{−α}Σ_j c_jd^jŜ(α−j)], |c_j| ≈ (ν/2)^j/j!, phases ≈ i^j. Σ_jc_jd^j = ₂F₁(a₀,a₀;b₀;d) is O(1) exactly
(Euler: |B·₂F₁(a₀,a₀;b₀;y)| ≤ coth(πν)/(2√(1−y)) for real y ∈ [0,1), since |BΓ(b₀)/Γ(a₀)²| = coth(πν)/(2π)), but termwise
Σ|c_j|d^j ≈ e^{νd/2}. Every termwise bound inherits it.

**Frozen-matrix resummation [P-level, exact identity; DV].** With x_j = Ŝ(α−j), V_j = (x_j, x_{j−1}), V_{j+1} = M_jV_j (TAIL E3),
and a frozen M̄ = M_{j*} with eigenvalues λ₁,₂: if x̄ solves the frozen recurrence with x̄₀ = x₀, x̄₋₁ = x₋₁, then exactly
Σ_jc_jd^jx̄_j = C₁H(dλ₁) + C₂H(dλ₂) (H = ₂F₁(a₀,a₀;b₀;·), C₁ + C₂ = x₀, C₁/λ₁ + C₂/λ₂ = x₋₁), and |H(dλ)| ≤ the Euler bound with
e^{ν|arg(1−dλt)|} ≤ e^{νd|Im λ|/(1−d|λ|)} — **no hump**. `dv_resum.py` (exact values):
| state (L, ν, a) | \|S\| = \|Σc_jd^jx_j\| | termwise | \|frozen resummation\| | termwise \|Σc_jd^j(x_j − x̄_j)\| |
|---|---|---|---|---|
| (3.2, 900, 1061), j* = 0 | 2.675 | 4.1e6 | 2.657 | 3.6e5 |
| (3.5, 1031, 998), j* = 7 | 2763 | 1.3e6 | 2469 | 2.6e5 |
| (2.5, 800, 1500), j* = 30 | 1.344 | 2.6e13 | 1.106 | 2.3e12 |
The frozen resummation is right to 1–20 %; **the drift correction is O(ξ)·|S| but termwise it is still humped**, so it must be resummed too.

**Expansion parameter.** The drift of M_j over the hump (j ≲ νd) is relative ~ j/ν in B_j, and the k-th order drift term, resummed through
derivatives H^{(m)}(dλ) (Σ_j c_jd^jλ^{j}j(j−1)⋯ = derivatives of H, again Euler-bounded), is of size ~ ξ^k/k!·|S| with
  **ξ = ν³d²/(4a²) ≈ νd²κ²/4** (0.25–0.38 at the states above; ≤ 274d² for a ≥ 998, ν ≤ 1030: 0.45 at L = 3.2, 2.7 at L = 2.3, 11 at L = 1.6).
This is the same parameter as FULLPOS3's h-expansion (ν³h, h = 1/D ≈ d²/(4a²)): the per-frequency route needs O(ξ + log(hump)/log(1/ξ))
orders of a matrix drift expansion, each order a finite combination of H^{(m)}(dλ_{1,2}) with explicit eigen-projections, and a Cauchy
remainder in the drift parameter ε (M_j(ε) = M̄ + ε(M_j − M̄), |ε| = R) bounded termwise by Lemma B with complex ε. Estimate: order ≈ 8–12
at L ≈ 3, ≈ 20–30 at L ≈ 2.3; uniformity over (ν, d, a)-boxes is free for the resummed terms (Euler bounds hold on boxes) and only the
remainder is box-sensitive. **Not implemented** (new lemma + 2×2 matrix drift calculus; several days). It cannot reach L ≲ 1.7 (ξ ≳ 10,
and Loc‴ itself fails at L ≈ 1.6 with J = 1500).

---

## 7. Which side is cheaper to finish (meet in the middle)

**Analytic (large-L) side.** After this session: remainder solved down to L ≈ 1.7 (§2), contour solved (§3); the single obstruction
is the termwise hump in the real-axis residual (and, later, on the top line). Termwise it ends at L ≈ 3.4–3.5 even with sharper per-power
bounds (§4). Going lower needs the drift-resummation lemma of §6: new mathematics (a 2×2 matrix drift calculus with Euler-bounded
derivatives H^{(m)}(dλ) and a Cauchy remainder), several days, with an expansion parameter ξ ≈ 274d² that reaches 3 near L ≈ 2.2 and
≳ 10 near L ≈ 1.6. **Per-frequency analytic methods cannot reach the small-L end (L ≲ 1.7) at all**, so the small-L machinery has to be
built in any case.

**Small-L side (FULLPOS4 / ANGULAR_ROUTE §§80, 85, 97).** Lemma U at fixed c is already spin-uniform: one degree-1001 polynomial in
y = x/x_* covers every a ≥ 998 (hence every J) at that c, and the restriction a < a*(L) is not needed. DV margins
(P̂₀ − Σ|P̂_ℓ|)/P̂₀ = 0.978 (L = 0.04), 0.972 (0.1), 0.985 (0.5), 0.9976 (2), 0.99999 (7), both channels; the strip constant
K = sup|H_M|/F_cert ≈ 1.3–1.4 for L ≥ 0.5 at η = 0.02, ≈ 5·10³ at L = 0.04. What is missing is entirely *execution*: (1) rigorous jet
enclosures of the integrated G (GL error on a Bernstein ellipse, 6000–8000 bits); (2) the realised EPSW weight instead of NF02
(with ν up to ~1130); (3) an Arb reproduction of the §78/§90 truncation error; (4) the c-direction by complex-strip (Lemma C) cells.
The c-direction cost is very non-uniform in L: at L ≥ 0.5 (K ≈ 1.3) cells of width ≈ η, i.e. ≈ 25 c-points per unit of β; at L ≲ 0.5,
≈ 1000 per unit of β (FULLPOS4 §8.2).

**Comparison.** β_c = 1.47 (L = 0.5), 3.31 (L = 2), 4.4 (L = 3), 5.0 (L = 3.55).
* Small-L side on [0.06, 0.5]: ≈ 1 unit of β × ~10³ points — the dominant cost of the whole programme, unavoidable.
* Small-L side on [0.5, 3.55]: ≈ 3.5 units of β × 25 ≈ 90 c-points, i.e. ≈ 25–50 core-h DV (≈ 3–10× certified): a 10 % marginal cost on
  top of the unavoidable part, with no new mathematics (the margins there are the best of the whole range: 0.985–0.9999).
* Analytic side from 3.55 down to ≈ 2: one new lemma (drift resummation, §6), its rigorous implementation, then ≈ 1–2 core-h per L-cell
  of 0.1 (≈ 15 cells). It saves only the cheapest stretch of the small-L side.

**Recommendation: the small-L side is cheaper to finish; meet at L ≈ 3.55 (here) or 3.65 (AS2).** Concretely: certify Lemma U +
§97 on c-cells covering β ∈ [0.486, 5.1] (L ∈ [0.06, 3.7]), with the overlap [3.55, 3.7] as a consistency check. The analytic route's
value below 3.55 is then only as an independent check; the drift resummation (§6) is the item to do if a second, analytic proof of the
middle range is wanted. A DV run of `fp4_lemU.py` at L ≈ 3.5 (between the tested 2 and 7) was not made here (CPU budget); by
the monotone trend of the margins it should pass with margin ≈ 0.999.

---

## 8. Grades, inputs, files

* **New and unreviewed here:** Lemma Loc‴ (§2) with Lemmas K (K₀ ≤ 1.5727 for every state), W (θ-weight), C (Chernoff B-weight),
  D (dressing factor 1 + 2.6·10⁻⁵); the Arb cell geometry (sector boxes, Im η enclosures); the contour's general-G j-tails and the
  right edge by Lemma F (σ-box) + Lemma B over a-boxes (§3); the a ≤ 4500 restriction and the finer a-cover of the kernel (§4).
* **Taken from AS/AS2 (unreviewed there except where REVIEW_SLIVER_L5 covers it):** Theorem E′, Lemma X, Lemma F, Lemma B, the helicity
  bound, Lemma NB, Lemma V, the splitting lemma, the per-cell kernel assembly (rt2_kernel), the residual bookkeeping incl. G1/G2 (rt2_resid),
  the tail-ratio lemma (AS2 §5, used for R_t).
* **Imported:** a*(L) = 4500 for β_c ∈ [3.04, 8.24] (ANGULAR_ROUTE §73, author's run); LRINF_TIGHT F_cert; FULLPOS3 §4.1(c) and Lemma B;
  EPSW masses (rt_rmass2, rmass3); CONTOUR tail/wtail; TAIL_CERT E1–E5; FSUP = 21.34.
* **Rigor notes.** As AS2: Arb balls → floats by `.upper()`; the Loc‴ log-sum-exp is in double precision over Arb-rounded inputs, with
  +10⁻⁶ added to every log; frames are floating-point choices, rigor does not depend on them.

Files (all new, `ads4/code/analytic3/`):
* lemmas: `loc5.py` (Loc‴), `kernel3.py` (kernel cover, a ≤ 4500, refined turning boxes; solution frames — tried, no gain on boxes);
* assembly: `locrun3.py` (E_ρ), `contour3.py`, `resid3.py`, `cert3.py`;
* checks / probes [DV]: `check_loc5.py` (weights vs exact), `probe_loss.py` (exact vs termwise vs bound), `dv_resum.py` (frozen resummation),
  `proto_loc5.py` (float prototype), `probe_resid.py`, `probe_resid3.py`, `probe_ranges.py`, `probe_solframe.py`, `probe_contour.py`;
* outputs: `cache/loc3_*.json`, `cache/contour3_*.json`, `cache/resid3_*.json`, `cache/cert3_*.txt`; logs in `logs/`.
