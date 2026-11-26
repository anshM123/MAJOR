# FULLPOS3: the open band 0.04 < L < 7 at Δ_g = 1000 — high-twist closure, signed profiles, coverage

*26 Nov 2026. Research agent. Unreviewed, not committed. No existing file was edited; nothing under `prl/` was touched.*
*New code: `ads4/code/fullpos3/` (outputs in `fullpos3/cache/`, whose `.gitignore` ignores everything but itself).*
*Compute: ≤ 2 processes, `ulimit -v 4000000`, `OMP_NUM_THREADS=1`.*

Grades: **[certified]** Arb balls or exact rationals with every remainder proved; **[P]** proved on paper here;
**[DV]** numerical; **[conj]** conjecture. Author/coordinator claims are marked as such and kept separate.

*(Running document; sections marked "in progress" are partial.)*

---

## 0. Verdict

**Not closed.** The infinite open band is reduced to a compact family. What is certified, at Δ_g = 1000, H_ang,b > 0:
1. **Per-J certificates (§1).** Every J ≤ 1000 at every mass, both channels; b = 4 up to J = 1004. Previously J ≤ 766.
2. **High-twist theorem (§4), certified cell by cell.** For every heavy state with β_c in a scanned cell and twist a ≥ a*(cell), both channels:
   * **a* = 3200 for 0.39 ≤ β_c ≤ 1.04 (L ≈ 0.04–0.25), with margins err/F_cert ≤ 0.16;**
   * a* = 4500 for β_c ∈ [7.44, 8.04];
   * a* = 4000 for β_c ∈ [8.24, 8.40] (L ≈ 7), with err/F_cert = 0.91.
   * The other cells (β_c from 1.04 to 7.44) were still being computed at hand-back (§4.6). The author's threshold was a ≥ 2.2·10⁷.
3. **Consequence (low L).** With n ≈ a√c/(1−√c), every state with a < 3200 and L ≲ 0.058 has J ≤ 1000, so items 1 and 2 cover the band
   0.04 < L ≲ 0.058 completely.
4. **What remains (when the scan is complete):** 998 ≤ a < a*(L) (≈ 3200–4500), J ≥ 1001, J ≲ a*·√c/(1−√c) (≤ ~10⁷ at L = 7). This is a compact
   two-parameter family.

**Key new structure (§3, §4.2).** For the polynomial target weight ψ₂ (degree 14, (1−p)⁵ | ψ₂), the Mehler–sin representation with coth and an
exact integration by parts give
|∫ψ₂P_{−½+iΛp}(cosh ξ)dp| ≤ (√2/π)[J_s + C_cor·J_C + e^{Λ|Im ξ|}J_e].
Only the p = 1 endpoint, of relative size ~10⁻⁶, carries the e^{Λ|Im ξ|} growth. The signed finite-mass correction is therefore tiny: exact
ℒ(G) agrees with G(c) to 4·10⁻⁷ at a = 998. This makes the h-expansion's local terms certifiable with plain Cauchy bounds, with no signed
profile evaluation.

**Obstruction below a*.** Per-order absolute local bounds on the two oscillatory components (the ψ₂ endpoint, and the residual weight at
ν ≈ 1000–1130) lose e^{|A|}, while the grouped multiplier is O(1) (§5). The naive grouped Taylor-shift bound fails by 10⁵ at a = 2000 (§5.1).
The integrated radial contour with the full conical function samples |G| ~ 10²⁰⁰ (§5.2). What is needed is an outgoing-wave (Im β ≥ 0)
integral-operator treatment of those two components.

---

## 1. Coverage extension: per-J certificates 767 ≤ J ≤ 1004 (task 3) [certified]

`fp3_polyext.py` continues FULLPOS2's route D. That is FULLPOS_CERT's per-J Bernstein every-mass certificate (`fullpos/fp_poly.py`,
imported read-only), run with the patched K = 1000 moment balls `moments_1000_K1000_an.json`.
* **Earlier runs (FULLPOS2):** odd J 767–849 (b = 4) and even J 768–810 (b = 0, 4), `fullpos2/cache/f2_polyext_*`.
* **This run:** odd J 851–1003 (b = 4) and even J 812–1004 (b = 0, 4), `fullpos3/cache/fp3_polyext_*`.

**Result: every channel with n = J − b ≤ 1000 PASSES.** There are 0 FAIL lines in the four files, and a script check finds no (J, b) missing.
Each polynomial needs 15 cells. The smallest Bernstein lower bound is 7.0·10⁻⁴⁵ (J = 1004, b = 4), with a ball radius ~10⁻⁵¹.

**Certified: H_ang,b > 0 at every mass for every J ≤ 1000 (both channels), and for b = 4 up to J = 1004.** This extends FULLPOS_CERT's
J ≤ 766. Going further needs moments beyond k = 1000.

---

## 2. Coordinator's resummed-ladder test, integrated [DV]

`fp3_resum_scan.py`. R_M(ν) = Σ_k g_k(ν)c^k exp[−Σ_{ℓ≤M}h^ℓL_ℓ(k)] summed without termination (M = 0: the Legendre kernel at c),
against the exact bare 𝒦₃ (b as stated), weight ψ₂, midpoint rule with 2400 points in p, N_s from T0's `hyb_fast`.
Columns: I_M = 2∫|ψ₂|w³|𝒦₃ − R_M|dp / N_s, and the signed difference S_M = 2∫ψ₂w³(R_M − 𝒦₃)dp / N_s.

| state (J, Δ, b) | a | c | L | I₀ | I₁ | I₂ | I₃ | I₄ | I₆ | S₀ |
|---|---|---|---|---|---|---|---|---|---|---|
| 800, 1800, 0 | 998 | 0.198 | 0.22 | 3.6e4 | 5e22 | 7e10 | 9e5 | 8e5 | 2.5e3 | +1.0e-6 |
| 2000, 3000, 0 | 998 | 0.445 | 0.59 | 3.2e4 | 1e22 | 1.5e3 | 1.1e4 | 52 | 3.5 | +2.1e-6 |
| 5000, 6000, 0 | 998 | 0.695 | 1.19 | 3.0e4 | 6e21 | 1.3e3 | 2.4e3 | 37 | 1.8 | +3.1e-6 |
| 800, 3800, 0 | 2998 | 0.044 | 0.045 | 9.1e3 | 70 | 0.97 | 1.9e-2 | 4.7e-4 | 4.2e-7 | −6.4e-8 |
| 2000, 5000, 0 | 2998 | 0.160 | 0.175 | 1.5e4 | 98 | 1.03 | 1.5e-2 | 2.7e-4 | 1.4e-7 | +7.7e-8 |
| 4000, 10000, 0 | 5998 | 0.160 | 0.174 | 4.0e3 | 5.9 | 1.5e-2 | 5.3e-5 | 2.3e-7 | 6.7e-12 | +7.7e-9 |
| 3000, 13000, 4 | 9998 | 0.053 | 0.055 | 9.9e2 | 0.61 | 7.0e-4 | 1.2e-6 | 2.4e-9 | 1.5e-14 | +4.3e-7 |

Findings.
1. **At a ≥ 3000 the truncated-exponent kernel is an excellent per-frequency reference**: I₂ ≈ 1 at a = 2998, I₂ = 0.015 at
   a = 5998, I₁ = 0.6 at a ≈ 10⁴; I_M falls geometrically in M. So the absolute (per-frequency) route works there with M = 3–4.
2. **At a = 998 it fails**: I₁ ~ 10²², and even M = 6 leaves I ~ 2–2500. The failure sits at p ≈ 0.97–0.99 (ν ≈ a, the turning point),
   where the untruncated k > n part (weight e^{−hL₁(k)} ≈ e^{−2a/3}) meets g_k ~ e^{πν/2}; and the ℓ-series in the exponent does not
   converge for k ≳ n.
3. **The signed integral is far better than any absolute one (new, decisive).** Even the plain Legendre kernel at c (M = 0) has
   |S₀| ≤ 3·10⁻⁶ N_s at a = 998, although I₀ ≈ 3·10⁴. I checked this against **exact** moment sums (§3.1): it is real, not a
   quadrature artefact.

---

## 3. The signed structure: ℒ(G) ≈ G(c) to 10⁻⁶ relative [DV + exact check]

### 3.1 Exact check
With the target weight ψ₂ (NF02 polynomial / k₀), the coefficients M_k = 2∫₀¹ψ₂(p)g_k(Λp)dp (g_k = [z^k]F_{3,ν}) are exact rationals
(Λ² = 997002). They are huge: |M_k| = 20.8, 3·10⁶, 2·10¹¹, 2·10⁴⁴, 6·10²⁷⁵, 6·10⁶²⁷ at k = 0, 1, 2, 10, 100, 400. The exact
finite-mass value Σ_{k≤n}U_kM_k (exact rational arithmetic, scratch script) versus G(c) = 2∫ψ₂F_{3,Λp}(c)dp (4096-node Gauss–Legendre):

| (J, Δ), b = 0 | a | c | Σ U_k M_k (exact) | G(c) | relative difference |
|---|---|---|---|---|---|
| 400, 1400 | 998 | 0.0820 | 2.58980118782e-5 | 2.58980003e-5 | 4.5e-7 |
| 300, 1300 | 998 | 0.0536 | 2.94535172712e-5 | 2.94534606e-5 | 1.9e-6 |
| 400, 2400 | 1998 | 0.0279 | 4.57819446639e-5 | 4.57817796e-5 | 3.6e-6 |

So at the twist edge a = 998 the whole finite-mass correction ℒ(G) − G(c) is ~10⁻⁶ of G(c): h·(signed first-order profile) is of
size h·G·O(1), not h·(Λ/d)³·(absolute size) as an absolute majorant gives.

### 3.2 Why [P-level explanation; made quantitative in §4]
G(z) = 2(1−z)^{−3}∫₀¹ψ₂(p)P_{−½+iΛp}(X)dp. By Mehler, P_{−½+iν}(cosh ξ) = (√2/π)∫₀^ξ cos(νt)(cosh ξ − cosh t)^{−½}dt, so the
p-integral enters only through Ψ(t) = ∫₀¹ψ₂(p)cos(Λpt)dp. Since ψ₂ is a polynomial of degree 14 with a fifth-order zero at p = 1,
integration by parts is exact:
∫₀¹ψ₂(p)e^{iωp}dp = Σ_{j=0}^{14}(−1)^j[ψ₂^{(j)}(1)e^{iω} − ψ₂^{(j)}(0)]/(iω)^{j+1}, and ψ₂^{(j)}(1) = 0 for j ≤ 4.
The p = 0 terms do not grow off the real t-axis; the only growth, e^{Λ|Im t|}, carries the factor (Λ|t|)^{−6}. Hence near z = c the
generating function is a smooth function plus an endpoint piece of relative size ~Λ^{−6}e^{Λ·Im ξ}, and its local derivatives are
O(1)^j·j! (smooth part) or Λ^{−6}(Λ√c/d)^j (endpoint part) instead of (Λ√c/d)^j. This is the mechanism behind §3.1, and the
reason the h-expansion's signed low-order terms are tiny.

---

## 2a. Contour identity of ANGULAR_ROUTE §53 [P + DV]
Termwise: with z = uv/((1+u)(1+v)), ∫₀^∞u^{n+s+k−1}(1+u)^{−N−s−k}du = B(n+s+k, a) and (1/2πi)∮v^{k−n−1}(1+v)^{N−k}dv = C(N−k, n−k)
(generalised binomial, real N); multiplied by 𝒞 these give exactly U_k. Dominated convergence (|q| = q₀ < 1 on the v-circle, u/(1+u) < 1)
allows summing over k. Numerically (mpmath, 30 digits, `scratchpad` script): exact vs contour 0.867139188498088 vs 0.867139188498088
at (J,b,Δ,ν) = (3,0,9.3,1.7); 1.30533893267044 vs 1.30533893267044 at (6,4,15.25,0.8). Both have non-integer N.

## 4. High-twist closure (task 0): an h-expansion bound with endpoint-aware sup norms

### 4.1 Exact set-up [P]
Notation of the brief; s = 2b+1, C′ = n(n+s) = C − b(b+1), T′ = N(N+s) = T − b(b+1), h = 1/T′, c = C′/T′ = U₁, τ = 1/T ≤ h.

**(a) Dressing as a differential operator [P].** The brief's dressing term G_C (weight Re[(3+ℓ)F₄ − (3+ℓ/2)F₃ − (ℓ/2)F_{3,ν+i}])
equals ϑG_A, ϑ = z∂_z: by (X²−1)P′_μ = μ(XP_μ − P_{μ−1}) with μ = ℓ, Re[(ℓ/2)P_{ℓ−1}] = −XP/4 − (X²−1)P′/2, and with X² − 1 =
4z/(1−z)² every P- and P′-coefficient matches z∂_zF₃ = 3z(1−z)^{−4}P + 2z(1−z)^{−5}P′. The ν²-weighted term follows from the ODE of
F₃ (§0 of `fp3_common.py`; s_ν = ν²+¼): (ν² + 25/4)F₄/2 = −z(1−z)F₃″/2 − (1−7z)F₃′/2 + 9F₃/2. Hence, per frequency and for any weight,
  **H/w³ = ℒ(G_A + τ𝒟_bG_A + τσ_bG_D),  𝒟_b = z(1−z)∂²/2 + [(7/2 − 2b)z + ½]∂ + 9/2 − b(b−1),**
  G_A = (2/Λ)∫f(ν)F_{3,ν}dν (realised f), G_D = Λ²∫₀¹ψ₃(p)F_{4,Λp}dp,
which is the author's §42 first-moment operator. (Also checked numerically: the operator form and the brief's dressed form agree to all
printed digits at 11 values of c in `fp3_profile_dv.py`; the 12th run, c = 0.632, had an Arb precision failure in the DV quadrature
and is discarded.)

**(b) Ladder factorisation [P].** U_k = c^kE_k, E_k = Π_{j<k}(1 − hu_j/c)/(1 − hu_j), u_j = j(j+s), so E_k ∈ [0,1] and E_k = 0 for
k > n. For k ≤ n, −log E_k(η) = Σ_ℓ η^ℓa_ℓ(k), a_ℓ(k) = ((c^{−ℓ}−1)/ℓ)S_ℓ(k), S_ℓ(k) = Σ_{j<k}u_j^ℓ, a power series in η with
nonnegative coefficients, convergent for |η| < c/u_{k−1}. P_m(k) := [h^m]E_k is a polynomial in k of degree 3m (P₁ = −(d/c)S₁,
S₁(k) = (k)₃/3 + (b+1)(k)₂, which gives the author's first-order profile −(d/X)[(b+1)cG″ + (c²/3)G‴]).
With T_M(k) = Σ_{m≤M}h^mP_m(k) and ρ_M = E_k − T_M,
  **Ĥ := H/w³ = Σ_{m=0}^{M}h^mP_m(ϑ)G(c) + ℛ_M,  ℛ_M = Σ_{k≤n}c^kg_kρ_M(k) − Σ_{k>n}c^kg_kT_M(k),**
G = G_A + τ𝒟_bG_A + τσ_bG_D = Σg_kz^k (all series absolutely convergent since c < 1).

**(c) Main term.** (1−c)³[G^ψ_A(c) + τσ_bG^ψ_D(c)] = 2Ω₂(β_c) + σ_bu_cΩ₃(β_c), cosh β_c = (1+c)/(1−c), with
u_c = Λ²τ/(1−c) ≤ Λ²/m² ≤ 1 on heavy states (b = 0: equality; b = 4: T(1−c) = T(T−C)/(T−20) > m²). So the main term is bounded
below by LRINF_TIGHT's certified **F_cert(β)** (min over u ≤ 1 and σ = ±1, `lrtight/cache/lrt3_rcert.json`).

### 4.2 Endpoint-aware bounds on |G| off the real axis [P]
**Lemma A (Mehler with coth).** For Re ξ > 0, |Im ξ| < π/2 and real ν > 0,
P_{−½+iν}(cosh ξ) = (√2/π)coth(πν)∫₀^∞sin(ν(ξ+u))(2sinh(ξ+u/2)sinh(u/2))^{−½}du
(analytic continuation from real ξ: the integrand is analytic in ξ, sinh(ξ+u/2) has positive real part, |sin(ν(ξ+u))| ≤ e^{ν|Im ξ|}).
For the polynomial weight ψ_k (degree 14, ψ_k(0) = 0, (1−p)⁵ | ψ_k) write coth(πΛp) = 1 + 2/(e^{2πΛp}−1) and integrate the
sine part by parts exactly: ∫₀¹ψ(p)e^{iωp}dp = Σ_{j≤14}(−1)^j[ψ^{(j)}(1)e^{iω} − ψ^{(j)}(0)]/(iω)^{j+1}. In the sine combination the
p = 0 terms with odd j cancel identically, so with t = ξ + u, ω = Λt, |ω| ≥ Λ(Re ξ + u):
|∫₀¹ψ(p)coth(πΛp)sin(Λpt)dp| ≤ A_s(Re t) + C_cor(|Im ξ|) + e^{Λ|Im ξ|}A_e(Re t),
A_s(x) = Σ_{j even ≥ 2}|ψ^{(j)}(0)|/(Λx)^{j+1}, A_e(x) = Σ_{j≥5}|ψ^{(j)}(1)|/(Λx)^{j+1},
C_cor(η) = ∫₀¹|ψ|·2cosh(Λpη)/(e^{2πΛp}−1)dp ≤ Σ_n|ψ_n|Λ^{−n−1}2n!Σ_{m≥1}(2πm − η)^{−n−1}.
With |2sinh(ξ+u/2)sinh(u/2)| ≥ 2sinh(β+u/2)sinh(u/2) (|sinh w| ≥ sinh Re w), β = Re ξ:
  **|∫₀¹ψ_k(p)P_{−½+iΛp}(cosh ξ)dp| ≤ (√2/π)[J_s(β) + C_cor(η)J_C(β) + e^{Λη}J_e(β)],**
J_X(β) = ∫₀^∞A_X(β+u)(2sinh(β+u/2)sinh(u/2))^{−½}du (A_C ≡ 1), all decreasing in β; certified upper bounds on a β-grid
(`fp3_htdata.jtable`: head ∫₀^{u₀} ≤ A(β)·2√u₀/√sinh β, monotone upper sums, exponential tail).
The e^{Λη} growth carries only the p = 1 endpoint coefficients, |ψ^{(5)}(1)| ≈ 10⁶ over (Λβ)⁶: J_e/J_C ~ 10⁻¹¹ at β = 0.6. The
dominant smooth part is the coth correction, ψ′(0)/(12Λ²) ≈ 1.9·10⁻⁵ — i.e. the bound is of the size of the leading-Regge
function itself, not of ∫|ψ||P| (which is 3·10⁴ times larger). This is §3.2 made rigorous.

**Lemma B (crude).** For Re ξ > 0, |Im ξ| < π/2: |P_{−½+iν}(cosh ξ)| ≤ e^{ν|Im ξ|}(|ξ|/Re ξ)P_{−½}(cosh Re ξ) and ≤ (π/2)e^{ν|Im ξ|}
(Mehler–cos along t = sξ; |cosh ξ − cosh sξ| ≥ cosh β − cosh sβ, resp. ≥ (2/π²)(1−s²)|ξ|² from |sinh w| ≥ (2/π)|w| for |Im w| ≤ π/2).
Used for the realised-weight residual r_θ = f_θ − ψ₂(ν/Λ) (with ∫|r_θ|e^{ην} from EPSW_CERT's certified panel sums, a θ-correction and
a Gaussian tail; `fp3_htdata.Rint`), and for the global Taylor-coefficient bound |g_k| ≤ B_k (circle |z| = 1 − 3/(k+3)).

### 4.3 Local terms [P]
For m ≥ 1, P_m(ϑ)G(c) = Σ_jp_{m,j}ϑ^jG(c) and Cauchy in t = log(z/c) gives |ϑ^jG(c)| ≤ j!R^{−j}S(R), S(R) = sup_{|t|=R}|G(ce^t)|, so
|h^mP_m(ϑ)G(c)| ≤ h^m‖P_m‖_R S(R), ‖P‖_R := Σ_j|p_j|j!R^{−j}.
* m ≤ m₁ (= 8): P_m computed exactly as arb polynomials (Faulhaber/Bernoulli for S_ℓ, exponential recursion mP_m = −Σ_i i a_iP_{m−i}).
* m > m₁: with |Q|_ρ = Σ|q_j|ρ^j (submultiplicative): (M1) |F_q|_ρ ≤ ρ^{q+1}B̂((q+1)/ρ)/(q+1), F_q(k) = Σ_{j<k}j^q,
  B̂(y) = Σ|B_r|y^r/r! = 2 + y/2 − (y/2)cot(y/2) (y < 2π); (M2) |P_m|_ρ ≤ η^{−m}exp(Σ_{ℓ≤m}η^ℓ|a_ℓ|_ρ) for every η > 0;
  (M3) ‖P_m‖_R ≤ max(1, (3m)!/(Rρ)^{3m})|P_m|_ρ (j!/x^j is log-convex in j); ρ = 3m/(eR).
* S(R) uses Lemma A for the ψ parts and Lemma B for r_θ; τ𝒟_bG_A on the circle by Cauchy from a larger circle. For a cell of c the
  geometry is computed at the cell centre c₀ on the radius R + δ_c (δ_c ≥ |log(c/c₀)|; maximum principle), with adaptive arcs and
  Arb enclosures of ξ = arccosh((1+z)/(1−z)).

### 4.4 Far remainder [P] (author's §61 lemma; my first version had a cell-monotonicity bug, see below)
ℛ_M = Σ_k g_k(U_k − T_M(k)), with U_k = 𝒰_k(h; c) = Π_{j<k}(c − hu_j)/(1 − hu_j) and T_M(k) its Taylor polynomial in h (= Σ_{m≤M}h^mc^kP_m(k)).
* **Analyticity.** 𝒰_k is analytic in |h| < 1/u_{k−1}. This holds also past termination: for k > n the numerator has the exact zero
  j = n, so 𝒰_k(h_phys) = 0 = U_k.
* **Sup over the cell.** On |h| = ρ (x_j = ρu_j < 1), uniformly in c ∈ [c₋, c₊]: the Möbius map w ↦ (c − w)/(1 − w) takes its extreme
  moduli on |w| = x at w = ±x. Hence sup|𝒰_k| ≤ M_k(ρ) = Π_{j<k}max{(c₊+x_j)/(1+x_j), |c₋−x_j|/(1−x_j), |c₊−x_j|/(1−x_j)}.
* **Two bounds.** Cauchy gives (i) |U_k − T_M(k)| ≤ M_k(ρ)(h/ρ)^{M+1}/(1 − h/ρ) for h ≤ h_max < ρ, and (ii) |U_k − T_M(k)| ≤ c₊^k +
  M_k(ρ)Σ_{m≤M}(h/ρ)^m for any admissible ρ.
* **Coefficient bound.** |g_k| ≤ B_k: Cauchy on |z| = 1 − 3/(k+3) with Lemma B (|P| ≤ (π/2)e^{ν|Im ξ|}); B_k = O(e^{1573}k⁵).
* **Evaluation.** On blocks of k (each factor is quasi-convex in x, so its block maximum is at an end point; B_k is increasing); then a
  geometric tail with ρ = θ/u_k, θ = min(¾, 0.45(1+c₋)), where M_k ≤ q^k with q < 1.
* **Bug in my first version (found by the coordinator, ANGULAR_ROUTE §61).** It used log(1/E_k(ρ; c)) at c = c₊, but that is decreasing in
  c, so it is not an upper bound over a cell; the validity test ρu < c also used c₊. All scan results produced before the fix were
  discarded, and the scan was rerun with the lemma above.
* **Effect of the radius.** Since 𝒰_k is analytic up to 1/u_{k−1} rather than c/u_{k−1}, the "termination barrier" d·n ≳ 2300 of my
  first version disappears. At β ∈ [2, 2.05] the far bound at a₀ = 1500, M = 600 is 2·10⁻⁵⁸ (it was 10²³⁰). At β ∈ [0.39, 0.44] it is
  2·10⁻⁵⁵ at a₀ = 1250, M = 600. At β ∈ [8, 8.2] it still fails at a₀ = 1500 for M ≤ 850, and needs larger M.

### 4.5 The certificate
For channel b, β-cell and a₀: **F_cert(cell) > (1−c)³[ |G^r(c)| + τ|𝒟_bG_A(c)| + Σ_{m=1}^{M}h^m‖P_m‖_RS(R) + |ℛ_M| ]** ⇒ H_ang,b > 0
for all heavy states with β_c in the cell and twist a ≥ a₀ (`fp3_ht.py`).

### 4.6 Scan results (certified; `fp3_ht.py scan`, results `cache/ht_scanA.jsonl`, `cache/ht_scanB.jsonl`) — in progress
Settings:
* β-cells of width 0.05 (β < 1), 0.1 (β < 3) and 0.2 (up to 8.4);
* bisection in a₀ over {998, 1100, …, 2.2·10⁶};
* 81 radii in [10⁻³, 1]·R_max, adaptive stadium arcs (tolerance 2·10⁻⁴ on Im ξ);
* exact P_m for m ≤ 70, majorant beyond, M chosen as the smallest of {100, …, 1600} with far remainder < 10⁻³F_cert.

A PASS at a₀ certifies every a ≥ a₀ in the cell, because all bounds are monotone in h_max. Bisection may miss a smaller passing a₀; that
costs sharpness, not rigour.

| β-cell | L (b = 0) | a*(b = 0) | err/F_cert | a*(b = 4) | err/F_cert |
|---|---|---|---|---|---|
| [0.39, 0.44] | 0.038–0.048 | 3200 | 0.096 | 3200 | 0.117 |
| [0.49, 0.54] | 0.060–0.072 | 3200 | 0.067 | 3200 | 0.075 |
| [0.59, 0.64] | 0.086–0.10 | 3200 | 0.063 | 3200 | 0.067 |
| [0.69, 0.74] | 0.12–0.13 | 3200 | 0.078 | 3200 | 0.081 |
| [0.79, 0.84] | 0.15–0.17 | 3200 | 0.097 | 3200 | 0.097 |
| [0.89, 0.94] | 0.19–0.21 | 3200 | 0.125 | 3200 | 0.123 |
| [0.99, 1.04] | 0.23–0.25 | 3200 | 0.164 | 3200 | 0.159 |
| [7.44, 7.64] | 6.1–6.3 | 4500 | 0.291 | 4500 | 0.291 |
| [7.84, 8.04] | 6.5–6.7 | 4500 | 0.306 | 4500 | 0.306 |
| [8.24, 8.40] | 6.9–7.0 | 4000 | 0.907 | 4000 | 0.907 |

* The bisection is not monotone between neighbouring cells (4500 vs 4000 at the top). Each entry is a certificate for its own a₀; the
  smallest certified value is not necessarily found.
* **The scan was still running at hand-back.** There are two nohup processes over the 60 cells, odd cells ascending and even cells
  descending. Results are appended to `code/fullpos3/cache/ht_scanA.jsonl` and `ht_scanB.jsonl` (one JSON line per cell with the full
  error budget). At ~13 min per cell per process, completion needs about 5 more hours.
* To resume after an interruption, rerun the same two commands: done cells are skipped.
  ```
  SP=<scratch>; python3 fp3_ht.py scan "<cells A>" ht_scanA.jsonl ; python3 fp3_ht.py scan "<cells B reversed>" ht_scanB.jsonl
  ```
  Cells: [0.39, 1.0) step 0.05, [1.0, 3.0) step 0.1, [3.0, 8.4] step 0.2 (the last cell [8.24, 8.4]). A takes cells 0, 2, 4, …; B takes
  cells 1, 3, … in reverse order.


What binds at a₀ < a*:
* the local orders, not the far remainder. Example: β ∈ [0.39, 0.44], a₀ = 2000, where per-order bounds reach 10¹⁵⁸ at m ~ 300;
* the dominant components are the ψ₂ endpoint (J_e·e^{Λη}) and the residual cells at ν ≈ 1000–1130 (EPSW cells plus θ/tail).

Each order gets an e^{A/e} amplification, A ∝ (Λ√c/d)³h.

## 5. Grouped (resummed) multiplier [DV]: why the h-expansion can reach a ≈ 998
For an oscillatory component G_X(ce^t) = e^{iλt}B(t) (carrier λ = Λ√c/(1−c) for the ψ₂ endpoint; ν√c/(1−c) for residual cells at
frequency ν), (𝒬_M(ϑ) − 1)G_X(c) = Σ_i (𝒬_M − 1)^{(i)}(iλ)B^{(i)}(0)/i!, with 𝒬_M(k) = Σ_{m≤M}h^mP_m(k).

The per-order bounds use Σ_m|h^mP_m(iλ)|. The grouped quantity is |𝒬_M(iλ) − 1| → |e^{−A(iλ)} − 1|, A = Σ_ℓh^ℓa_ℓ(iλ)
(scratch script `symbol.py`, 1500–3000 bits):

| c | λ | a₀ | \|A\| | max_m\|h^mP_m(iλ)\| | \|e^{−A} − 1\| |
|---|---|---|---|---|---|
| 0.042 | 210 | 3200 | 4.3 | 16 | 1.66 |
| 0.042 | 222 | 3200 | 5.1 | 32 | 1.12 |
| 0.042 | 210 | 2000 | 10.8 | 9.8·10³ | 1.58 |
| 0.042 | 222 | 2000 | 12.6 | 7.0·10⁴ | 0.064 |
| 0.042 | 210 | 1400 | 20.9 | 9.7·10⁸ | 1.72 |
| 0.042 | 222 | 1400 | 24.5 | 6.2·10¹⁰ | 0.65 |
| 0.042 | 210 | 998 | 38 | 8.1·10¹⁸ | 0.34 |
| 0.042 | 222 | 998 | 44 | 4.6·10²² | 0.18 |
| 0.58 | 1810 | 998 | 71 | 5.6·10³⁴ | 1.67 |

At a₀ = 3200 and 2000 the grouped value is reached by M ≈ 20–50; per-order terms are < 10⁻¹¹ by m = 50. Below 2000 only the converged
limit was computed, not the convergence in M.

The resummed multiplier is O(1) (it is the Darboux multiplier), while the per-order majorant loses up to 10³⁴. A rigorous grouped bound
needs:
* (i) an exact carrier split of the oscillatory components with sup bounds on B free of e^{ΛR}. For the ψ₂ endpoint this is the author's
  §64 A_±; for residual cells, the outgoing-wave split P = P⁺ + P⁻, with |P^±e^{∓iνβ}| bounded uniformly in ν by the Euler-integral
  form of §55;
* (ii) bounds on |(𝒬_M − 1)^{(i)}(iλ)| uniform over (c, h) boxes, via Cauchy in k on E = e^{−A} plus the h-truncation tail.

Estimate: Σ_i|Q^{(i)}(iλ)|R^{−i} ≈ |E|Σ_i(3|A|/(λR))^i. This needs λR ≳ 3|A|, limited by the nonlinear-phase factor
exp[(ΛR²/2)·√(ce^R)(1+ce^R)/(2(1−ce^R)²)] (author's §64b). That suggests a ≈ 1500–2000 at low L is reachable, and a = 998 only
marginally there.

**Turning-point obstruction to a ≈ 998 [P-level remark].** The h-series of the ladder symbol E(k) = U(k)/c^k (Γ-ratio continuation) at complex
index k converges for |k|²h/c ≲ 1, i.e. |k| ≲ √C′ ≈ n. The ψ₂-endpoint carrier is k = iλ, λ = Λ√c/(1−c), and n ≈ a√c/(1−√c), so
λ/n ≈ Λ/(a(1+√c)).
* At a ≈ Λ ≈ 998 the endpoint carrier sits at 0.5–0.83 of the convergence radius: λ = 210 against n ≈ 250 at c = 0.04. This is the
  turning point ν ≈ a of §§47–49, seen from the ladder side.
* Cauchy estimates in k around iλ therefore have radius ≲ n − λ (~40 at c = 0.04, a = 998), too small to control the high k-derivatives.
* So the h-expansion can approach a = 998 only where c is not small, or if the endpoint component (size ~ J_e ~ 10⁻¹¹ relative) is bounded
  non-perturbatively at the turning point, e.g. by an Airy-window bound as in §48 but for the finite-ladder symbol.

### 5.1 Naive grouped carrier bound: tested, fails [DV] (`fp3_grouped_dv.py`)
Setup:
* G_e(ce^t) = Σ_±e^{±iλt}B_±(t), with |B_±| ≤ 2|1−z|^{−3}(√2/π)(1/2)J_e(Re ξ)·e^{Λ|Im φ|} and |Im φ| ≤ (R²/2)max|ξ″|;
* local term bound Σ_±S_B(R)Σ_i|q_i|i!R^{−i}, where q_i are the Taylor coefficients at ±iλ of the grouped polynomial
  𝒬₇₀ − 1 = Σ_{m=1}^{70}h^mP_m(k).

At β = 0.415 (c = 0.042, λ = 213), a₀ = 2000, b = 0, the best bound is 2.5·10⁵ F_cert (R = 0.4); the per-order version gives 2.1·10⁶.

| R | S_B | Σ\|q_i\|i!R^{−i} (grouped) | (per-order) |
|---|---|---|---|
| 0.10 | 3.1e-10 | 9.7e102 | 9.8e102 |
| 0.20 | 2.8e-9 | 1.3e49 | 1.4e49 |
| 0.30 | 1.7e-7 | 2.9e21 | 3.7e21 |
| 0.40 | 7.6e-5 | 3.5e4 | 2.9e5 |
| 0.60 | 3.1e4 | 1.67 | 1.7e5 |

**Why it fails.** The grouped value at the carrier itself is O(1) (q₀ = 𝒬(iλ) − 1), but the high k-derivatives of a degree-D polynomial
approximating e^{−A(k)}, A ≈ αk³, force R ≳ α^{1/3}D^{2/3} (≈ 0.4 here). At that radius the carrier's quadratic phase costs
e^{ΛR²|ξ″|/2} (e⁹–e³³). Grouping in m alone is therefore not enough: Taylor-expanding the amplitude B at t = 0 is the lossy step.
* **Needed instead:** treat the multiplier as an integral operator. For the cubic part, e^{−α∂_t³} acting on e^{iλt}B is an Airy
  convolution of width ~α^{1/3} (cf. the author's §52 Airy form). Then only B on |t| ≲ few·α^{1/3} enters. This is not attempted here.
* **But absolute Airy integrals re-amplify.** ∫Ai(u)e^{pu}du = e^{p³/3} along real or fixed-ray contours samples e^{iλt} at complex shifts
  and again costs e^{O(|A|)} (the author's §52 envelope factor exp(ν³/(6a²))). Only a steepest-descent contour through the saddle of
  Ai(u)e^{iλ(3κ)^{1/3}u}, κ = hd/(3c), avoids it.
  * The saddle shift is t_* ≈ 3κλ² ≈ hdλ²/c (≈ 0.16 at c = 0.042, a = 2000). There the quadratic carrier phase costs only
    e^{Λt_*²|ξ″|/2} ≈ e^{1.4}.
  * For the cubic part this is feasible. The quintic and higher ladder terms are not small at ν ~ a (author §52), so the full-ladder
    symbol is needed, i.e. the author's §53/§59 contour with steepest descent, applied only to the two oscillatory components (ψ₂
    endpoint and residual cells at ν ≳ 900). Their weight relative to the floor is ~10⁻⁵, so a crude O(1)·a^{1/6} bound on their
    finite-ladder transform would suffice.

### 5.2 The §59 radial contour applied to the integrated G: what it needs [analysis, not executed]
* Applied to the full conical G, the radial contour ℒ(G) = −(i𝒞/2π)∫(e^{iBα}/sin α)∫e^{−av}ξ^{s−1}G(ξ²)dv dα samples z = ξ² over much of
  the disk, including the negative half-axis region. For B = 2J+1 ~ 10³ the localisation |α| ≤ 8800/B is void.
  * There |G^ψ| ~ Λ^{−6}e^{Λ|Im arccosh X|} (e.g. ~e^{1190} at z = −½), so ∫|·| ≫ |∫|.
  * The per-frequency non-amplification of §55 rests on the OUTGOING split P = P⁺ + P⁻ with Im β ≥ 0 on each piece's contour.
* **Measured [DV, Arb, 4096-node Gauss–Legendre in p].** The sample points come from the exact map
  ξ = sinh(v/2)/sinh(v/2 + iα), α = x − iε, ε = 10⁻³: for example (v, x) = (2·10⁻³, −1.2·10⁻³) gives z = 0.0865 + 0.162i and Im arccosh X = 0.50. Values:

  | z | G^ψ(z) | Σ\|terms\| |
  |---|---|---|
  | 0.0865 + 0.162i | (4.8 + 2.5i)·10²⁰¹ | 1.4·10²⁰³ |
  | 0.3 + 0.02i | (−1.59 + 0.23i)·10⁹ | 2.2·10¹⁶ |
  | 0.3 | 4.07·10⁻⁵ | 1.6 |

  These points carry O(1) contour weight (e^{−av} with v = 2ε). Hence, with the full conical function, ∫|integrand|/|ℒ(G)| ≳ 10²⁰⁰. The integrated radial
  contour is usable only after the outgoing split.
* Integrated version: G⁺(z) = 2∫ψ₂P⁺_{Λp}dp, the author's §58 𝒫⁺. Its ψ-integral has the same endpoint structure as Lemma A:
  * the p = 0 terms are non-oscillatory, of size ψ′(0)/(Λβ)² ~ the floor;
  * the p = 1 terms are e^{iΛβ}/Λ⁶ and decay for Im β ≥ 0.
  * So on the outgoing contour |G⁺| ≲ few × floor, and ∫|·| ≈ (log B)·few·F_cert. That is the favourable case for a certified 2D quadrature.
* Uniformity over the compact family (a ≤ 4000, J ≤ 9·10⁶) then needs parameter cells in (B, a). The carrier e^{iBα} with |α| ≲ min(π/2, 8800/B)
  gives relative cell widths ~ 1/(B·A₀) ~ 1/8800 in B, so of order 10⁴ cells per decade of B for each a-cell, before using the extracted-carrier
  trick of §57. I did not build this (time).

### 5.3 §54 global contour, literal reading [DV, failed]
A direct mpmath evaluation of (𝒞/π)Re(I_U − iI_T − iI_L) with 𝒲_p as transcribed in ANGULAR_ROUTE §54 (scratch `c54.py`) returns NaN at
(J, b, Δ, ν, p) = (2, 0, 7.5, 0.8, 3), where the exact value is 1.41545100. The likely cause is the diagonal C = A (ξ = 0, β → ∞) and the branch
choices. It needs the author's precise conventions; not pursued further.

## 6. Obstructions and dead ends recorded here
1. **Coordinator's truncated-exponent kernel R_M per frequency (§2).** It fails at a = 998: R₁ ~ 10²² N_s, from the k > n region at ν ≈ a. It works at
   a ≳ 3000.
2. **Far remainder with circular Cauchy contours in h of E_k = U_k/c^k (my first version).** The radius is c/u_{k−1}, which gives a termination barrier
   d·n ≳ 2300, and it had a cell-monotonicity bug (§4.4). Superseded by the 𝒰_k radius-1/u_{k−1} lemma (author §61). At a ≲ 1500 for β ≳ 8 it still
   needs M > 850; the author's vertical-line version (§§62–63) is not implemented here.
3. **Per-order absolute local majorants.** They lose e^{|A|} against the grouped multiplier (§5): ~10⁴ at a = 2000, up to 10³⁴ at a = 998. This sets
   a* ≈ 3200 at low L.
4. **Grouped bound near a ≈ Λ.** The endpoint carrier iλ approaches the convergence radius ≈ n of the h-series of the ladder symbol (the turning point),
   so Cauchy-in-k control degenerates there (§5).
5. **Naive grouped carrier bound (Taylor shift of the grouped polynomial at the carrier).** Fails by 10⁵ at a₀ = 2000 (§5.1).
6. **Per-frequency envelopes and generic Racah norms.** Dead, as stated in the brief and in §51. Not re-tried.

## 7. Reproduce (`ads4/code/fullpos3/`; `ulimit -v 4000000; export OMP_NUM_THREADS=1`; ≤ 2 processes)
```
python3 fp3_polyext.py 851 1003 2; python3 fp3_polyext.py 812 1004 2     # §1 per-J certificates (≈ 2.5 h)
python3 fp3_resum_scan.py 800,1800,0 2400 6                                # §2 (one state per call, ~1 min)
python3 fp3_profile_dv.py 0.0392,0.1,0.5 64 64                             # §3/§4.1 DV profiles (~1 min per c)
python3 fp3_ht.py 0 2 41/20 5000                                           # §4 single certificate (verbose)
python3 fp3_ht.py scan "0.39,0.44;0.44,0.49" out.jsonl                     # §4.6 cell scan (~12 min per cell)
```
| file | content | grade |
|---|---|---|
| `fp3_common.py` | ψ data, Taylor coefficients of F_{p,ν} at c via the ODE (checked by finite differences) | — |
| `fp3_polyext.py` | per-J every-mass certificates | [certified] |
| `fp3_resum_scan.py` | coordinator's R_M test, integrated | [DV] |
| `fp3_profile_dv.py` | h-expansion profiles at c (ψ weights) | [DV] |
| `fp3_htdata.py` | J_s, J_C, J_e tables, C_cor, ∫\|ψ\|e^{κΛp}, residual exponential moments (EPSW cells), F_cert cells | [certified] inputs |
| `fp3_ht.py` | the high-twist certificate (§4): stadium geometry, exact/majorant P_m, far remainder, cell scan | [certified] |
