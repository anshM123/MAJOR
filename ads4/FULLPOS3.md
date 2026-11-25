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

## 4. High-twist closure (task 0) — in progress
