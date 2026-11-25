# FULLPOS2: the open large-spin region at Δ_g = 1000 — kernel closed forms, route tests, coverage extension

*26 Nov 2026. Research agent. Unreviewed, not committed. No existing file was edited; nothing under `prl/` was touched.*
*New code: `ads4/code/fullpos2/` (outputs in `fullpos2/cache/`, whose `.gitignore` ignores everything but itself).*
*Compute: ≤ 2 processes, `ulimit -v 4000000`, `OMP_NUM_THREADS=1`.*

Grades: **[certified]** Arb balls or exact rationals with every remainder proved; **[P]** proved on paper here (proof written out);
**[DV]** numerical; **[conj]** conjecture. Author/coordinator claims are marked as such and kept separate.

*(This file is updated as work proceeds; sections marked "in progress" are partial.)*

---

## 0. Verdict (running)

* **Not closed.** See §6 for the obstruction list.

---

## 1. Exact kernel identities

### 1.1 Coefficient recurrence for F_{p,ν} [P]
Let H(z) = P_{−½+iν}((1+z)/(1−z)), X = (1+z)/(1−z), s = ν² + ¼, λ(λ+1) = −s.
With X′ = 2/(1−z)² and 1 − X² = −4z/(1−z)², Legendre's equation (1−X²)P_XX − 2XP_X − sP = 0 becomes, after the chain rule,
z(1−z)²H″ + (1−z)²H′ + sH = 0, i.e. (1−z)²θ²H + szH = 0 (θ = z∂_z).
Put H = (1−z)^p G. Then H′ = (1−z)^p[G′ − pG/(1−z)], H″ = (1−z)^p[G″ − 2pG′/(1−z) + p(p−1)G/(1−z)²], and dividing by (1−z)^p:
(z − 2z² + z³)G″ + (1 − (2+2p)z + (1+2p)z²)G′ + (s − p + p²z)G = 0.
Comparing coefficients of z^k (G = Σg_k z^k):
**(k+1)² g_{k+1} = (2k² + 2pk + p − s) g_k − (k+p−1)² g_{k−1},  g₀ = 1, g₋₁ = 0.**
(p = 0 reproduces the d = 0 CDH recurrence; p = 2, 4 give (k+1)²Q⁽⁰⁾_k and D_k²q_k.) Used in `f2_num.py`, `f2_fast.py`, `f2_dress.py`.

### 1.2 The physical k = 2 kernel is the shift-image of a polynomial [DV, 60+ digits]
Φ₂(ν) := w³Σ_{k=0}^{n+1} V_{b,k} g⁽³⁾_k(ν) (V as in ANGULAR_ROUTE §32 with u₋₁ = (b−1)b) satisfies
𝒦_fin(ν) = −((4x+9)(4x+49)/(4608(4x+1)))·Re[R₋(ν)Φ₂(ν+i) + R₊(ν)Φ₂(ν−i)] with 𝒦_fin = `fp_kernel.kfin`.
Checked at (J,Δ,b) ∈ {(12,40,0),(30,1030,4),(100,1100,0)}, ν ∈ {0.7, 5, 23}: differences ≤ 5·10⁻⁶² (`f2_num.py`).
This is the author's §46 weight conversion H⁽²⁾ = ∫g𝒦_fin = (2/Λ)∫fΦ₂ at the kernel level.

### 1.3 Author's dressing formula (§45) [certified: exact rationals]
w⁻³Φ₂ = (1 + (9−b(b−1))/T)𝒦₃ − ((ν²+25/4)/(2T))𝒦₄ + ((7−2b)/T)·Re[(3+ℓ)𝒦₄ − (3+ℓ/2)𝒦₃ − (ℓ/2)𝒦₃(ν+i)], ℓ = −½+iν,
with 𝒦_p = ℒ(F_{p,ν}). Checked **exactly** (Fractions; 𝒦₃(ν+i) as a Gaussian rational) on 40 random states (b ∈ {0,4}, J ≤ 18,
Δ ≤ J+40, rational ν): 0 mismatches; the imaginary part of the right side vanishes identically on every sample (`f2_dress.py`).
So the author's dressing agrees with fp_kernel's 𝒦_fin, via §1.2.

### 1.4 The dual-Hahn closed form (§45) is exact but NOT tame [DV]
The coordinator verified 𝒦_p = U_n((p)_n/n!)Σ_j(−n)_j(½±iν)_j/((p)_j j!²)R_{n−j} exactly. Conditioning (mpmath, 400 digits):

| (J, Δ, b, ν) | 𝒦₃ | max \|term of j-sum\|·prefactor | max \|R_m\| |
|---|---|---|---|
| (60, 1062, 0, 50) | 0.0648 | 17 | 8.6·10¹⁶⁰ |
| (100, 1300, 0, 200) | 0.0717 | 3.0·10¹¹ | 1.0·10²⁴⁴ |
| (200, 1300, 4, 300) | 0.0804 | 9.2·10³⁷ | 1.4·10³⁶⁵ |

R_m (dual Hahn outside the orthogonality range) is astronomically large, and the j-sum cancels by up to 10³⁸ already at J = 200.
**The convolution form does not remove the e^{O(Λ)} cancellation;** it relocates it. The direct k-sum (§1.1 recurrence) is the
cheaper evaluator: O(n) operations at ≈ 0.75ν + 60 digits.

---

## 2. Primary route (per-frequency envelope) — decisive measurement [DV]

Setting: H⁽²⁾ = 2∫₀¹ψ₂(p)Φ₂(Λp)dp (target weight; realised differs by the certified de-shift residual), hybrid
Φ_hyb(ν) = D_b(1−2ν²/m²) + P_{−½+iν}(cosh β) − J₀((2J+1)ν/m), D_b = ₂F₁(−n, n+2b+1; 1; ·) (T0's flat partial wave; it is
ℒ applied with the flat ladder, §34). `f2_scan.py`.

**J = 800, b = 0, b̂ = 450 (Δ−J = 2843, r = 0.048):** N_s = 2.69·10⁻⁵.

| ν | Φ₂ | Φ_hyb | Φ₂ − Φ_hyb |
|---|---|---|---|
| 1 | 0.9390447 | 0.9390443 | 4·10⁻⁷ |
| 20 | −0.0736948 | −0.0736950 | 1.2·10⁻⁷ |
| 50 | −0.150575 | −0.150563 | −1.3·10⁻⁵ |
| 100 | 0.117093 | 0.116861 | 2.3·10⁻⁴ |
| 300 | 0.00371 | 0.00853 | −4.8·10⁻³ |
| 500 | −0.0482 | −0.0117 | −0.037 |
| 776 | −0.0407 | 0.1223 | −0.163 |

* 2∫|ψ₂||Φ₂ − Φ_hyb| = **1.6·10⁴ N_s**; 2∫|ψ₂||Φ_hyb| = 3.4·10⁴ N_s. The signed integrals agree (1.7448·10⁻⁵ both, coarse quadrature).
* The hybrid kernel is pointwise excellent for ν ≲ 20 (10⁻⁶ relative) and **O(1) wrong** for ν ≳ 300: the exact kernel's phase
  drifts from both the flat and the Legendre phase. Leading mechanism [D]: ℒ acts on e^{iνξ} like a cumulant expansion in the
  rapidity ξ, whose third cumulant gives a phase ≈ ν³√r/(3w²T) (≈ 2 rad at ν = 700 here); T0's D_b reproduces its w → 1 part
  ((2J+1)ν³/(6m³)) but not the rest, and higher orders (author §47: the multiplier 𝔑_a has phase iν³/(12a²) + O(ν⁵/a⁴),
  turning point at ν = a) are O(1) at ν ~ Λ in the open region.
* **Consequence.** Any per-frequency absolute envelope E(ν) needs ∫|ψ|E ≤ 0.2N_s, i.e. relative pointwise accuracy ~10⁻⁵ up to
  ν ~ Λ. That requires an essentially exact reference phase. The coordinator's J = 200 measurement (2.0·10⁴ N_s) is the same
  phenomenon. **The primary route as stated (|𝒦_fin − 𝒦_ref| integrated against |g|) is dead** for both the Legendre and the
  hybrid reference; it can only be used for the low band ν ≲ ν₀ (§4).

## 3. The author's Darboux reference (§§46–47) [DV]
R(ν) = 2Re{B(ν)d^{−α}[𝔑_a(α) + ℭ/M]}, α = p−½+iν, d = 1−c, M = J+½+a/2, ℭ = a[(2b−α)𝔑_a(α) + (p−2b)𝔑_a(α−1)],
B(ν) = Γ(iν)4^{−½+iν}/(√πΓ(½+iν)). Compared with the exact bare 𝒦₃ (`f2_darboux.py`, `f2_fast.py`) at the twist edge a = 998, b = 0:

| J | d | ν = 0.5 | ν = 5 | ν = 20 | ν = 50 | ν = 100 | ν = 400 | ν = 700 |
|---|---|---|---|---|---|---|---|---|
| 2000 | 0.555 | 7·10⁻³ | 1.6 | 5.1 | 30 | 30 | — | — |
| 10⁴ | 0.173 | 3·10⁻⁴ | 0.24 | 0.25 | 5.0 | 8.5 | 85 | 57 |
| 10⁵ | 0.0197 | 8·10⁻⁴ | 1.1·10⁻³ | 0.062 | 0.087 | 0.96 | 4.7 | 9.7 |

(entries: |R − 𝒦₃|/|𝒦₃|). At these states the plain Legendre kernel at c is far better (1.6·10⁻⁶ … 3·10⁻² for ν ≤ 50 at J = 2000).
**Reason [P-level remark]:** R keeps only the j = 0 term of the connection expansion
(1−z)^{−p}P_{−½+iν}(X) = Σ_± B(±ν)(1−z)^{−α_±}₂F₁(½∓iν, ½∓iν; 1∓2iν; 1−z). The j-th term is relatively ~ (νd/2)^j/j!, so R is a
valid reference only when νd ≪ 1, i.e. d ≪ 10⁻³ (L ≳ 7) for the ν ~ Λ part of the weight. Moreover each branch separately has a
logarithmic singularity at z = 0 (₂F₁(a,a;2a;y) at y = 1), and ℒ only reads Taylor coefficients at z = 0, so the branch split
cannot be fed to ℒ without the divergent connection sum (the coefficients of ₂F₁(a,a;2a;y) decay only like 1/j while ℒ((1−z)^{−α+j})
grows polynomially (degree n) in j). **The Darboux route is a large-L (tail) tool; it cannot reach the moderate band
0.04 < L ≲ 7, which is where most of the open region lives at J ≤ 10⁶.**

## 4. (in progress)

## 5. Route D: coverage extension (in progress)
`f2_polyext.py` re-runs FULLPOS's per-J every-mass certificate (fp_poly, imported read-only) with the patched K = 1000 moments for
767 ≤ J ≤ 1004 (n ≤ 1000).

## 6. Obstructions (running list)
1. Per-frequency absolute envelopes vs Legendre or the T0 hybrid: lose ~10⁴ (§2).
2. Dual-Hahn convolution: intermediate terms up to 10³⁶⁵ (§1.4).
3. Darboux (j = 0 connection) reference: O(1)–O(10²) pointwise relative error at moderate L (§3).
