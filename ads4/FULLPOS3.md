# FULLPOS3: the open band 0.04 < L < 7 at Δ_g = 1000 — high-twist closure, signed profiles, coverage

*26 Nov 2026. Research agent. Unreviewed, not committed. No existing file was edited; nothing under `prl/` was touched.*
*New code: `ads4/code/fullpos3/` (outputs in `fullpos3/cache/`, whose `.gitignore` ignores everything but itself).*
*Compute: ≤ 2 processes, `ulimit -v 4000000`, `OMP_NUM_THREADS=1`.*

Grades: **[certified]** Arb balls or exact rationals with every remainder proved; **[P]** proved on paper here;
**[DV]** numerical; **[conj]** conjecture. Author/coordinator claims are marked as such and kept separate.

*(Running document; sections marked "in progress" are partial.)*

---

## 0. Verdict (running)

* Not closed yet. See §§1–4.

---

## 1. Coverage extension, per-J certificates 767 ≤ J ≤ 1004 (task 3) [certified, in progress]

`fp3_polyext.py` continues FULLPOS2's route D (FULLPOS_CERT's per-J Bernstein certificate, `fullpos/fp_poly.py` imported read-only,
patched K = 1000 moment balls `moments_1000_K1000_an.json`). FULLPOS2 had done odd J 767–849 (b = 4) and even J 768–810 (b = 0, 4),
all PASS (`fullpos2/cache/f2_polyext_*`). This run: odd J 851–1003, then even J 812–1004. Log: `fullpos3/cache/polyext.log`,
results `fullpos3/cache/fp3_polyext_*.txt`. (Status at the end of the document.)

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

## 4. High-twist closure (task 0): an h-expansion bound with endpoint-aware sup norms

### 4.1 Exact set-up [P]
Notation of the brief; s = 2b+1, C′ = n(n+s) = C − b(b+1), T′ = N(N+s) = T − b(b+1), h = 1/T′, c = C′/T′ = U₁, τ = 1/T ≤ h.

**(a) Dressing as a differential operator [P].** The brief's dressing term G_C (weight Re[(3+ℓ)F₄ − (3+ℓ/2)F₃ − (ℓ/2)F_{3,ν+i}])
equals ϑG_A, ϑ = z∂_z: by (X²−1)P′_μ = μ(XP_μ − P_{μ−1}) with μ = ℓ, Re[(ℓ/2)P_{ℓ−1}] = −XP/4 − (X²−1)P′/2, and with X² − 1 =
4z/(1−z)² every P- and P′-coefficient matches z∂_zF₃ = 3z(1−z)^{−4}P + 2z(1−z)^{−5}P′. The ν²-weighted term follows from the ODE of
F₃ (§0 of `fp3_common.py`; s_ν = ν²+¼): (ν² + 25/4)F₄/2 = −z(1−z)F₃″/2 − (1−7z)F₃′/2 + 9F₃/2. Hence, per frequency and for any weight,
  **H/w³ = ℒ(G_A + τ𝒟_bG_A + τσ_bG_D),  𝒟_b = z(1−z)∂²/2 + [(7/2 − 2b)z + ½]∂ + 9/2 − b(b−1),**
  G_A = (2/Λ)∫f(ν)F_{3,ν}dν (realised f), G_D = Λ²∫₀¹ψ₃(p)F_{4,Λp}dp,
which is the author's §42 first-moment operator. (Checked numerically against the dressed profiles: agreement to all printed digits
at 12 values of c, `fp3_profile_dv.py`.)

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

### 4.4 Far remainder [P]
|ℛ_M| ≤ Σ_k c^kB_k·b_k with, for every admissible ρ (ρu_{k−1} < c):
* k ≤ n_min (then hu_{k−1} < c, so E_k(h) is its own Taylor series): |ρ_M(k)| ≤ (h/ρ)^{M+1}/E_k(ρ) (ρ ≥ h), because 1/E_k = exp(A_k)
  has nonnegative Taylor coefficients dominating those of E_k = exp(−A_k);
* all k: |ρ_M(k)| or |T_M(k)| ≤ 1 + Σ_{m≤M}(h/ρ)^m/E_k(ρ).
Evaluated on blocks of k (all factors monotone in k) up to K_T = max(n_min+1, 6(2M+10)/d), then a geometric tail with ρ_k = (d/2)c/u_k.
Every bound is increasing in h and decreasing in n, so a state set {c ∈ cell, a ≥ a₀} is covered by h_max = h(c_lo, a₀),
n_min = n(c_lo, a₀).

### 4.5 The certificate
For channel b, β-cell and a₀: **F_cert(cell) > (1−c)³[ |G^r(c)| + τ|𝒟_bG_A(c)| + Σ_{m=1}^{M}h^m‖P_m‖_RS(R) + |ℛ_M| ]** ⇒ H_ang,b > 0
for all heavy states with β_c in the cell and twist a ≥ a₀ (`fp3_ht.py`).
