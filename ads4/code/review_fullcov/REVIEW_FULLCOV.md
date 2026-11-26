# Referee report: `code/sliver_fullcov/sliver_merge/PROOF.md` (coupled-cap norming certificate)

Scope: the mathematics of PROOF.md §§1–6, and whether `verify_merged.py`, `check_ledgers.py`, the drivers and
`scalar_ref_portable.cpp` implement it faithfully. The ledgers' fast checks and coefficient authentication were
done by the coordinator. The full GMP replay of every shard is running separately and is **not** covered here.
Inherited and not refereed: the angular-action identity, the coefficient inventory, the 1e-14 coherent
allowance and the 2^-140 tail beyond degree 17199 (G2 / §125), and the four direct-triangle cells on [0.258, 0.30].

Every check is in this directory. `rc.py` is an independent fixed-point evaluator of
F_K = (1−c)^3 Σ_{k≤K}(A_k + hT_{k,b})U_k, built from the compact centres; its error is about 4(K+1)2^−W, with W = 256–400.

## Verdicts by section

### §1 Physical cap, coupled rectangle, corner checks: HOLDS
* **Cap.** Set f(t) = √(t(t+s)). Then f′ = (2t+s)/(2f) ≥ 1, because (2t+s)² − 4t(t+s) = s². So a = N − n ≤ f(N) − f(n) = (1−√c)/√h, which gives h ≤ (1−√c)²/998².
* **Secant.** H_*(c) = (1−√c)²/998² is decreasing and convex, since −√c is convex. Hence the secant through upward-rounded endpoint caps dominates it, and physical states lie in y ∈ [0, 1].
  * Checked numerically: 66 random physical states (integer n, a ∈ {998,…,3000}, both b, all 11 new cells) have y ≤ 0.99967.
* **Corner sufficiency.** With h = yH(c), both 1 − p h and 1 + c − 2p h are bilinear in (c, y). The ratio (c − p h)/(1 − p h) is a Möbius map, hence monotone, in each variable separately. So the extremes of |ratio| sit at the corners, and the binding index is j = K−1.
  * This is also what `exact_schedule` and the per-node C++ test (`2Hp ≤ Den+C`) use.
* **Tail usage.** The ordinary truncation tail uses U_k ≤ c^k ≤ u^k and h ≤ Hl (`check_tail`, j = K+1…17199). It is subtracted only in the final physical inequality, after normalisation, so it is used only where valid.
* **§6 group tail.** It uses exact corner maxima of |ratio| (signed numerators, checked ≤ 1) over the rectangle enlarged by η_y. It does not rely on c^k or on termination, so it is valid on the nonphysical enlarged group.
* **Test.** `test_endtoend`: at physical points, F_K agrees with the terminating F_n, n < 17199, far inside the recorded tail.

### §2 c-approximation (Möbius disk, Bernstein ellipse): HOLDS
* At fixed real y, each factor is a real Möbius map of c, because H is affine. Its pole is real and lies outside the closed disk iff the denominator is positive at both real diameter endpoints; this is asserted.
* The image disk is conjugation-symmetric, so its farthest point from 0 is the image of a real endpoint. The prefactors |A + yH(c)T| and |1−c|³ are moduli of affine functions, so their maxima also sit at real endpoints.
* At the real extended endpoints, y ↦ factor is monotone (positive denominator), so maximising over y needs only y ∈ {lo, hi}.
* E_ρ ⊂ disk of half-length (ρ+ρ⁻¹)/2, and |a_n| ≤ 2Mρ^−n gives the tail 2Mρ^−m/(ρ−1). Per-summand choice of ρ is legitimate.
* `group_check` / `exact_check` code exactly this, using pairs (c0 ± R·dc, H(c0 ± R·dc)·y) for y ∈ {lo, hi}, both s, and suffix-free upper coefficients AI, TI.
* **Tests.**
  * Monte-Carlo over complex disk points: max|factor| / endpoint max = 0.99999 ≤ 1.
  * `test_capprox`: true Chebyshev coefficients in c (DCT, N = 1024, cells 5305_6041 and 8936_9000, y in three groups) have tails ≤ about 1e-20 at the certified degrees. The recomputed majorant reproduces the certified E_c (1.8e-14 … 2.0e-12) and dominates at every m tested.

### §2 y-approximation (partial fractions): HOLDS
* **Residue formula.** I re-derived it via v(v+s) = c p_j:
  * c p_j − p_i = (v−i)(v+i+s) and p_j − p_i = (j−i)(j+i+s);
  * the j-factor is (c−1) = (v−j)(v+j+s)/p_j;
  * this gives |ρ_{j,k}| = (2j+s)/p_j · |(−v)_k(v+s)_k| / [j!(k−1−j)!(j+s)_k].
* **Limit and monotonicity.**
  * Limit = (2j+s)|sin πv| Π_{r<s}(j+r)/(v+r) / (π p_j) ≤ (2/j)c^{−(s−1)}, using v ≥ cj. The 1/π slack is dropped.
  * ρ_{j,k+1}/ρ_{j,k} = (p_k − c p_j)/(p_k − p_j) > 1.
  * The log of the infinite product is Σ_r (α^r − β^r)/(r X^r), with X = (ℓ+s/2)², α = p_j + s²/4 and β = c p_j + s²/4. All terms are nonnegative.
  * The integral comparison gives the 4-term bound, which increases in j and decreases in k. Using u in β, l^−8, k = K, j = a, 1/p_j ≤ 1/(a(a+1)) and the suffix max from a+1 are each in the conservative direction.
* **Chebyshev tail.** For 1/(1 − a y), the tail is 2r^{m+1}/(g(1−r)) with g² = 1−a and r = (1−g)/(1+g); it is increasing in a, so using Hl is valid. The T term enters as T_k/p_j, correctly, since h/(1−p_j h) = (1/p_j)[1/(1−p_j h) − 1].
* **`test_residues`.**
  * The closed formula matches exact rational residues to 5e-58.
  * Partial-fraction identity: exact.
  * Monotone in k.
  * In 780 cases, max |ρ_{j,k}| / bound = 0.307 (≈ 1/π).
* **`test_yapprox`** (cell 5305_6041, K = 3193, c ∈ {l, mid, u}, both b). The exact signed pole coefficients are computed with all 3192 poles.
  * The residues are as large as 1e691, and F − Σ Res_j/(1 − p_j H y) is affine to 1e-122, so the representation is confirmed through 813 digits of cancellation.
  * max |Res_j| / (code's per-pole bound) ≤ 2.4e-10.
  * The exact truncation majorant is ≤ Ey(m) for all m tested; at the certified m = 446 it is 6e-62 … 5e-115 against Ey = 6.1e-14.
  * The DCT coefficients agree with the residue prediction.

### §3 Norming: HOLDS
* **Λ_m.** Order the circular distances as t, h−t, h+t, …, mh ± t, with t ≤ h/2 and mh + t ≤ π. Pair k is ≥ (2k−1)π/(2m+1). With |D| ≤ min(2m+1, 1/sin(δ/2)) and sin(δ/2) ≥ δ/π, this gives Λ_m.
  * Numerically the true constant is ≤ 0.59Λ_m (m ≤ 908), and the algebraic constant at the m+1 cosine nodes is ≤ 0.57Λ_m.
* **C_m(η).** This is Markov on [0,1] (2m²); the node errors are given in x units. Actual nodes stay inside [0,1]:
  * c is clamped into [l,u];
  * y′ = floor(…) ≤ ỹ ≤ 1, with ỹ_0 = 1 exactly.
* **§3a.** At each actual c′, the y-nodes are one perturbed set, so the order "y then c" is legitimate. The formula is reproduced.
* **§3b.** Re-derived step by step:
  * |f(c′, y_i) − f(c′, y′)| ≤ 2Ey + γ(M+Ey), with ‖P‖ ≤ M + Ey;
  * c-norming at the ideal y_i, using E_g over the group's ideal range;
  * then exact-node y-norming, giving M ≤ Λ_y(T + γ·max C_g·M + Ey) + Ey, i.e. M ≤ [Λ_y(T+Ey)+Ey]/(1−λ), with λ = Λ_y γ max C_g. M is finite, so this is valid.
  * The rounding dependence of y on c is charged correctly. λ ≤ 3.3e-3 in the archived cells.
* **Short groups.** The tail enters both δ_g and E_g, and the Markov step uses the full-K Ey. Correct.

### §4 Directed scalar arithmetic: HOLDS
* **Schedule invariants.** upper_k ≥ 2^100 · 2^{s_k} Π R_j ≤ 2^100, and the partial-product bound ≤ 2.
* **Per-step error.** In C++ order:
  * the quotient truncation is < 1 unit before the ≤ 2-bit shift, so < 4 after;
  * the A conversion is < 1 (tdiv composes);
  * the 2T conversion times h′/2, plus the shift by 62, is < 2.
  * The total is < 10 units of 2^{s_k−P}, and propagation is ≤ 2^{−s_k}.
* **Coefficient radii.** At most (k+1)² + |k2| + 2 for T, which the stated (K+1)² + K² + 8K + 67 covers. The full K is used even for K_g, which is conservative.
* **Reference error.** The Horner truncation count is ≤ 8 units, ∂R ≤ 2e-3 per coordinate, and both coordinates are floored to 2^−50.
* **Node convention.** R is evaluated at the *actual* y′ = h′/H(c′) (`yv`), not at the node ỹ. This is consistent with §3.
* **`test_scalar_cpp`.** Portable GMP on 18 nodes from three archived shards (P = 180 ref, P = 80 group, P = 80 short group with K_g = 5822):
  * |scalar − independent F_K| ≤ 8e-5 · E_scalar;
  * the printed [lo, hi] encloses the independent F_K − R.
* **Recurrence identity.** The recurrence is the same as the direct method's, S_k = F_k + r_k S_{k+1}, with F_k = A_k + h(B_k + qA_k + σG_k). The +2q sits inside k2, the signs are σ_0 = +1, σ_4 = −1, and the normalisation is (1−c)³/(1+qh)^4.

### §5 Final inequality: HOLDS
* It is (B_R − M)/(1+qHl)^4 − tail − 2^−140 − 1e-14, where B_R is the minimum of the exact Bernstein coefficients of R on [0,1]², and the inequality requires B_R − M > 0.
* `verify_merged` recomputes every term and takes the max of the archived and fresh error bounds.
* Wording issue only: §5 calls the subtracted term "E_short_tail", but the code subtracts the *physical* ordinary tail beyond K. The §6 group tails are inside δ_g and E_g.

### Code faithfulness: HOLDS
* **`verify_merged.py`.** Same formulas as stated in the proof, not looser ones:
  * cap and corner checks;
  * Ey through `residue_degree`, E_c through `group_check`;
  * Λ_m and C_m(η) via `norm_constant`;
  * ε_c = 2e-24 + 2^−61/(u−l) and ε_y = 2e-24 + 2^−61/Hu;
  * γ, λ < 0.1, T and M exactly as in §3b;
  * exact node reconstruction and input matching.
* **Inherited cells.** The three earlier norming cells use Chebyshev–Lobatto nodes and C = 5/2. That rests on T(θ*+t) ≥ ‖T‖cos(mt) (Bernstein–Szegő), with m(π/2n + π√η) < 1.13 and cos 1.13 > 0.4; this holds.

## End-to-end sanity (independent of all archived code)
* **`test_endtoend`.** 16 random and corner points per new cell and channel, biased toward y ≈ 1.
  * max |F_K − R| / M_certified ≤ 0.0114 in all 22 channel-cells. The certified residual over-covers by about 90×.
  * min F_K ≥ B_R − M everywhere sampled.
* **`test_oldcells`.** On the three earlier norming cells, max |F_K − R| ≤ 0.16 × certified residual.
* **Margins.** M ≈ 2e-8 against B_R ≈ 5–10e-6. An error would have to be structural (wrong function, domain or argument), not a looseness in a constant. No structural error was found.

## Overall
No gap or error found in §§1–6 or in the code's implementation of them.
* **Positivity on c ∈ [0.258, 0.9], both b, stands**, conditional on:
  1. the retained upstream inputs: the angular identity, the coefficient table, the 1e-14 coherent allowance and the 2^-140 tail;
  2. the four direct-triangle cells on [0.258, 0.30], which were not refereed here;
  3. the pending full GMP replay of the archived SIMD shards.
* For item 3, I spot-checked only 18 nodes, which agree.
