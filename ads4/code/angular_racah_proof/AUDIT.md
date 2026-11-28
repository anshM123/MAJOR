# Independent audit: uniform heavy feasibility of the angular Racah dual

Audited: `angular_racah_feasibility/angular_racah_feasibility.tex` (SHA-256 matches `SHA256SUMS.txt`), `README.txt`, and `coefficient_ledger.json`. Source conventions were taken from `../angular_opt/{lp_primal,lp_dual,racah_explore}.py` and `../angular_verify/README.md`.

Method: every identity and inequality below was re-derived by hand and then checked with code written from scratch in SymPy and exact `Fraction` arithmetic. The author's verifier was **not** used as input. It was run once, on a copy in the scratchpad, only to confirm that it still passes (44 checks). The audit scripts are `a1.py` to `a11.py` in `/tmp/claude-0/-home-user-MAJOR/60eea3a7-07d6-5453-af55-c4346d86f9b5/scratchpad/audit/`. That scratchpad is temporary. The load-bearing certificate reconstruction is reproduced in the Appendix.

Grades:
- **[V]**: verified exactly (symbolic or integer arithmetic) and checked by hand.
- **[V-hand]**: elementary step, checked by hand.
- **[minor]**: a presentational gap that does not affect correctness.

## Verdict

**The theorem is correct as stated.** I found no mathematical error.
- Every displayed identity reproduces exactly.
- The logical chain closes on the full stated domain: real g ≥ 4, n = ⌊g/2⌋, real Δ ≥ J+g, even J ≥ 0 for b=0, and every integer J ≥ 4 for b=4.
- My own reconstruction of the degree-five certificate gives the same 2 × 169 coefficients as the ledger. All of them are strictly positive.

What remains is scope and presentation, not correctness. The main scope point is item 8: the theorem is a heavy-state theorem only. With n = ⌊g/2⌋ and the scale W_B, the dual is **infeasible on light spin-2 single traces for every g**. Its relative margin there is between about −7% and −44%, which is much worse than the −0.14% quoted in `angular_verify/README.md`; that figure was for n = 53, not 50. A rescaled dual fixes this, and I give an exact proof below.

## Normalization check (definitions and signs)

- `lp_primal.ladder` uses r_{k+1}/r_k = ((k+4)/(k+1))² (C − s(s+1))/(T − s(s+1)) with s = b+k. Since D_{k+1}/D_k = (k+4)/(k+1) for D_k = C(k+3,3), this gives r_k = D_k² U_k^{(b)}.
- `lp_dual` imposes W_4 + Σλ_k r_k^{(4)} ≥ 0 and W_0 − Σλ_k r_k^{(0)} ≥ 0. With λ_k D_k² = −W_B(−n)_k(n)_k/(k!)², these are exactly (main0) and (main4). This is also `racah_explore` with s = 0 and A = −W_4(g+4). **[V]**
- The shifted variables match: C − (b+i)(b+i+1) = X_x − i(i+s) and T − (b+i)(b+i+1) = Y − i(i+s), checked for b = 0, 4 and i ≤ 3. Summing to min(n, x) equals the formal sum, because P_k(X_x) = 0 for k > x. **[V]**
- Every denominator in the finite sum is positive on the domain: Δ − b − i − 2 ≥ g − 1 for i ≤ J−b−1. **[V-hand]**
- Weights: W_4 = (Y + 25/2)/3 and W_B = (g² + 5g − 3/2)/3. **[V-hand]**
- The dual objective is consistent with the independent LP values in `angular_verify`. At g = 10, 20, 40 the dual bound on |c|² (0.0404, 7.29e-4, 1.25e-5) sits above the LP primal optimum (5.04e-3, 1.54e-4, 3.49e-6), as weak duality requires.

## 1. Spin recurrence (Lemma 2.1, Prop. 2.2): correct [V]

- **The P_k identity.**
  - Checked as a polynomial identity for k = 1..8 with symbolic x and s.
  - Checked for **symbolic k** through the three ratio formulas. Each ratio was confirmed from P_k(X_x) = Γ(x+1)Γ(x+s+k)/(Γ(x−k+1)Γ(x+s)) by `gammasimp`, and the combined identity cancels to 0.
  - The x = 0 boundary was checked separately.
- **The coefficient telescoping.**
  - Multiplying by (Y−X)/(2x+s) and using (Y−X)P_{k−1} = [Y−(k−1)(k−1+s)]P_{k−1} − P_k, the coefficient of P_k is a_{k+1}(k+1)²[Y−k(k+s)] − (k² − ν)a_k.
  - This vanishes exactly for the stated a_k. At k = n the top coefficient is (ν − n²)a_n = 0.
  - Because {P_k(X)} is a basis of polynomials in X, this is a complete proof for every x and every real N with Y ≠ i(i+s).
- **Direct check of the recurrence.** I checked it against the finite ₄F₃ for x = 0..8 with symbolic ν (so all n at once), symbolic Y, and s ∈ {1, 9, symbolic}.
- **Initial values.** y_0 = 1, y_1 = 1 − (s+1)ν/Y, and the x = 0 case p_0(y_1 − y_0) + q_0 y_0 = 0. **[V]**
- [minor] For s = 1, p_{−1}^{(1)} = 0·0/0 is indeterminate. The "set to zero" convention is just the statement F_1 = c_0, which follows directly from the finite sum. This is harmless.

## 2. Helicity 0: correct [V]

- **The energy.** I checked:
  - both forms of E_c;
  - the contraction E_c(v,w) = (1+c)v² + r²(u−cv)²/(1−c) when w = (1+r)cv − ru;
  - ∂_c E_c = (u−v)²/(1−c)²;
  - E_{c_0}(1, c_0) = 1 + c_0.
- **The chain.**
  - The recurrence is F_{x+1} = (1+r)c_x F_x − rF_{x−1} with r = x/(x+1) ∈ [0,1).
  - The chain E_{c_x}(F_x, F_{x+1}) ≤ E_{c_x}(F_{x−1}, F_x) ≤ E_{c_{x−1}}(F_{x−1}, F_x) is valid when c_x ≤ c_{x−1} and |c_x| < 1.
  - T − x(x+1) > ν for x ≤ J gives −1 < c_x < 1. The values c_x decrease because T − x(x+1) decreases and stays positive.
  - The final bound is E ≥ (1+c_J)F_J². **[V-hand]**
- **The bound.** (1+c_0)/(1+c_J) − 1 = νC/(T(T−C−ν)) holds exactly. **[V]** I also checked the bound itself on 300 random exact rational points, after confirming there that the recurrence values equal the ₄F₃.
- **Domain inequalities.**
  - T − C ≥ (2g−4)J + (g−1)(g−2) ≥ g² + 5g − 14 ≥ g² ≥ 4n². **[V]**
  - Hence T − C − ν ≥ 3ν, which gives F_J² ≤ 1 + C/(3T).
  - T − B_g = (J−4)(J+2g+1) + 15/2 and 2h − C = (J−4)(J+4g−3) + 4 hold exactly at T = (J+g−1)(J+g−2). Both increase with T. **[V]**
- **The weight comparison.** The step (1 + h/B)² > 1 + C/T is left implicit. It follows from (1+h/B)² ≥ 1 + 2h/B > 1 + C/B > 1 + C/T. **[V-hand]** The margin bound (b0margin), W_0 + W_B F_J > C/9, is also correct.
- **J = 0.** Trivial, since F_0 = 1.
- **J = 2.**
  - The formula for F_2 is exact. **[V]**
  - T ≥ g(g+1) > 4ν gives F_2 > 1 − 3/2 = −1/2.
  - W_0 − W_B/2 ≥ (g² − 3g + 21/2)/6, exactly at T = g(g+1). Its discriminant is 9 − 42 < 0. **[V]**

## 3. Helicity 4, Sonin energy: correct [V]

- **The auxiliary sequence h_j.**
  - h_j + h_{j+1} = 2/p_j was checked by exact integration for j = 0..7. The integrands add to 2t^j(1+t^8).
  - h_j > h_{j+1} > 0 and h_j < H_j were checked. The second uses (1+t)(1−t/2) = 1 + t(1−t)/2 ≥ 1.
  - h_5 = 545393/180180 − 4 log 2 ≈ 0.25435, which is below 4/15 = H_5. **[V]**
- **Exact decrement (Lemma 5.1).**
  - I checked it with p_{j−1}, p_j, p_{j+1}, q_j, q_{j+1}, h_j, y_j, δ_{j+1} all free, h_{j+1} = 2/p_j − h_j, and δ_j, δ_{j+2} eliminated by the recurrence. The difference is identically 0.
  - The square-completion form (soninsquare) is also an exact identity. **[V]**
- **The sign argument.**
  - Since q_{j+1} ≥ q_j > 0 and h_{j+1} > 0, the bracket is at most (h_j − h_{j+1})p_j(1/2 − p_j/q_j).
  - q_j ≤ 2p_j follows from z_j ≤ 1/4 and p_j/(2j+9) − 1/8 = (2j² + 21j − 9)/(8(j+5)(2j+9)) ≥ 0 for j ≥ 1. **[V]**
  - z_j ≤ 1/4 follows from (N−j)(N+j+9) ≥ (2n−2)(2j+2n+7) and (2n−2)(2j+2n+7) − 4n² = (4j+10)n − 4j − 14 ≥ 0 for n ≥ 2. This uses N − j ≥ N − x ≥ g − 2 ≥ 2n − 2. **[V]**
  - q_j increases while its denominator stays positive. All indices used (j+1 ≤ x) lie in the range where this holds. **[V-hand]**
- **Coercivity.**
  - The numerator 2j⁴ + 41j³ + 160j² − 869j − 2502 becomes 2v⁴ + 81v³ + 1075v² + 4806v + 3528 under j = 5 + v. **[V]**
  - So h_j q_j/4 ≤ (2j+9)H_j/16 ≤ 1/3, which gives E_j ≥ (2/3)y_j² for j ≥ 5.
  - The same bound fails for j ≤ 4 (values −4/5, −15/44, −5/39, −47/5460). This is why the tail starts at j = 5 and x = 0..4 are handled separately. That is consistent with the paper.
- **The reduction y_x² ≤ (3/2)E_x ≤ (3/2)E_5.**
  - The decrement is used only for j = 5..x−1.
  - The value y_{x+1} appears only as an algebraic auxiliary; coercivity at j = x eliminates it through the recurrence at x. **[V-hand]**
  - Exact-rational spot checks over 150 random (g, x, N) found E_{j+1} − E_j < 0 throughout and E_j − (2/3)y_j² ≥ 0. A first float run (mpmath, 50 digits) showed spurious violations at n ≈ 80; these were cancellation in the ₄F₃ sum. Exact arithmetic removes them, which is a warning against floating-point checks of this family.

## 4. The E_5 certificate: correct [V]

- **E_5 itself.** E_5 from its **definition**, with the actual y_6 from the finite sum and no recurrence elimination, equals (E5direct). Here p_4 = 65/18 and q_5 = 19ν/(Y−70). **[V]**
- **The polynomials A, P and the rational form.**
  - A = 𝒟 y_5 and P = 𝒟(y_4 − y_5)/(18ν) are polynomials.
  - A equals the factored expression (Aexplicit), and P equals (Pexplicit).
  - E_5(h) = (19A² + h𝓑)/(19𝒟²). **[V]**
- **The polynomials F_H.**
  - I built F_H as 60·19·𝒟²·(4n² + 18n + 25/2)² × [(2/3)((Y+25/2)/(4n² + 18n + 25/2))² − E_5(H/15)], starting from my own E_5.
  - The result is a polynomial identical to the paper's (FH), and (FHequiv) holds. **[V]**
- **Positivity after substitution.**
  - After n = u+1, Y = 4n² + 30n + 36 + t, both H = 0 and H = 4 have exactly the support {(a,b): 0 ≤ b ≤ 12, 0 ≤ a ≤ 24−2b}, with 169 terms each.
  - All coefficients are strictly positive, and the minimum is 760.
  - All 13 row minima in Table 1 and their positions a = 24−2b match.
  - All 338 coefficients agree with `coefficient_ledger.json`. **[V]**
- **Affine in H.** deg_H F_H = 1, so F_0, F_4 > 0 gives F_H > 0 on [0, 4], that is h ∈ [0, 4/15]. **[V]**
- **h_5 ≤ 4/15.** This holds, because 2/6 + 2/14 − 1/7 − 1/15 = 4/15 and h_5 < H_5. **[V]**
- **Domain coverage.**
  - For x ≥ 5: N ≥ x + g − 2 ≥ g + 3 ≥ 2n + 3, so Y = N(N+9) ≥ (2n+3)(2n+12) = 4n² + 30n + 36.
  - E_5 depends on (n, Y) only, so the single inequality in Y covers every physical x ≥ 5 and every real Δ.
  - The denominators are safe: Y ≥ 70 > 52 when n ≥ 1, and Y ≥ 112 when n ≥ 2. **[V-hand]**

## 5. Low spins x = 0..4: correct [V]

- **x = 1, 2, 3.** The closed forms for 1 − y_1, 1 − y_2 and 1 − y_3 (with B_3) are exact.
- **x = 2 domain.** T ≥ (2n+5)(2n+4) gives T − 3n² − 27 ≥ n² + 18n − 7 > 0.
- **x = 3.**
  - B_3(n², T_0) = 4(n⁴ + 297n³ + 4405n² − 1089n − 338) and ∂_T B_3 = 18(3n² + 88n − 11) at T_0, with ∂_T² B_3 = 72. **[V]**
  - The auxiliary claim 4405n² − 1089n − 338 ≥ 2978n² is an equality at n = 1. That is fine, because the claim is non-strict and the n⁴ and n³ terms keep B_3 strictly positive.
- **x = 4.**
  - 1 − y_4 = 26νQ/D_4 holds exactly.
  - The sum-of-squares identity (x4SOS) holds exactly under n = m+2, N = 2n+2+t.
  - The direct expansion has exactly two negative monomials, −28m⁵ and −8m⁵t. The square 2m⁴(m − 2t − 7)² absorbs them, as stated.
  - The domain is covered: N = Δ − 6 ≥ g + 2 ≥ 2n + 2, and Y ≥ 90 makes D_4 > 0. **[V]**
- **Upper bound suffices.** Only y_x ≤ 1 is needed, because W_4(Δ) ≥ W_B on the domain. Saturation occurs at (4, g+4). **[V-hand]**

## 6. Weight conversion: correct [V-hand]

- W_4/W_B = (Y + 25/2)/(g² + 5g − 3/2).
- Since g² + 5g is increasing for g > 0 and g < 2n+2, g² + 5g − 3/2 < 4n² + 18n + 25/2. The inequality is in fact strict.
- Chaining gives y_x² ≤ (3/2)E_5 < ((Y + 25/2)/(4n² + 18n + 25/2))² ≤ (W_4/W_B)², so W_4 − W_B y_x > 0 for x ≥ 5.

## 7. Numerical spot checks (exact rationals): no violations

The two conditions were evaluated in exact `Fraction` arithmetic, using the source normalization directly.

| Set | Points | g | J | Δ − J − g | Min. relative margin |
|---|---|---|---|---|---|
| Random (4 seeds) | 63,000 | rational in [4, 200]; edge values such as 7 − 10⁻³, integer g | 0 to 2000; both ladders, odd J on b=4 | 0, up to 10⁶ in rational steps, and multiples of g | b=0: 0.582 (at J=4, Δ=g+4, g ≈ 21.8); b=4: 0 at (4, g+4) only (saturation) |
| Grid | — | 65 quarter-integers in [4, 20], 26 integers up to 200, plus 199.999, 200.05 and 7 − 10⁻⁶ | up to 2g + 40, plus 5g, 20g and 1001 | 0, 1/3, 1, g/3, g, 4g | b=0: 0.585; b=4 with J ≥ 5: 0.161 (J=8, gap, g = 199.999) |

No point violates either condition.

## 8. Scope: light states are not covered, and the dual as given fails on light spin 2

The theorem is stated, and proved, only on H_0(g) ∪ H_4(g), that is for twist ≥ g. It says nothing about light states, and that is correct for a heavy-state theorem.

Corollary 8.1 does not claim more: its hypotheses put ρ_0 and ρ_4 on H_0 and H_4. **Its premise therefore excludes every light state other than T itself.** That means no light single traces of spin 0 or 2 (which the paper's gap assumption allows) and no light multi-traces. The README's τ ≤ 240/[…] is conditional on this. It should not be read as a bound for a CFT whose only assumption is the spin ≥ 4 gap.

What happens on the light b=0 states:
- **Light spin 0 is automatically fine.** R_{n,0}(0, Δ) = 1, and W_0 = T/3 + 3/2 ≥ 17/12 for every real Δ, since T ≥ −1/4. So the condition holds for any positive scale.
- **Light spin 2 (3 < Δ < g+2) fails with W_B.** Only λ_0, λ_1, λ_2 enter. The relative margin min (W_0 + W_B F_2)/W_0 is:

  | g | 4 | 6 | 10 | 20 | 40 | 100 | 200 | 1000 | g → ∞ |
  |---|---|---|---|---|---|---|---|---|---|
  | margin | −0.385 | −0.443 | −0.414 | −0.293 | −0.195 | −0.124 | −0.099 | −0.078 | 5/9 − 4√2/9 = −0.0730 |

  The worst point is at T ≈ (6 − 3√2)n², that is Δ ≈ 0.66g. `v10_racah_light.py 100 50` independently gives −0.125. The −0.14% in `angular_verify/README.md` was run at n = 53, not at n = ⌊g/2⌋.

**Extending the theorem to light spin 0/2 single traces works.** A rescaled dual fixes it; see section 9, where I prove the collaborator's repair exactly.

**Light multi-traces (A2′) remain completely outside both the theorem and the repair.** They enter every R_k on **both** ladders at twist < g, with r_k up to about D_k². Extending to them needs one of two things:
- (i) a feasibility statement for the Racah dual on the whole unitarity region of both ladders at twist below g. That is not plausible: W_B·|R| can far exceed W at low twist, and `v6` reports enormous dual values there.
- (ii) a separate quantitative bound on the light spectral weight E_light, the c_T → ∞ remainder described in `angular_verify/README.md` (a).

The corollary should state explicitly that it needs this input.

## 9. The collaborator's light single-trace repair (`LIGHT_SINGLE_TRACE_REPAIR.md`): correct, now with an exact proof

**Claim.** Replace W_B by A = 3(√2−1)n², with n = ⌊g/2⌋. The dual λ_k D_k² = −A(−n)_k(n)_k/(k!)² is then feasible on:
- both heavy ladders;
- light scalars;
- light spin-2 states with 3 < Δ < g+2.

It gives τ ≲ 994.4/g³. The note supports the light spin-2 step only numerically (mpmath, g ≤ 10⁴). All parts are confirmed below, and the numerical step is replaced by a proof.

1. **A ≤ W_B for every real g ≥ 4.** [V]
   - Since g ≥ 2n, W_B ≥ (4n² + 10n − 3/2)/3.
   - (4n² + 10n − 3/2)/3 − 3(√2−1)n² = [(26 − 18√2)n² + 20n − 3]/6 > 0, because 26 − 18√2 ≈ 0.544.
2. **Heavy conditions are preserved.** [V-hand]
   - At fixed (J, Δ), W_0 + tR_{n,0} and W_4 − tR_{n,4} are affine in the scale t.
   - Both are positive at t = 0 (W_b > 0), and nonnegative or positive at t = W_B by the theorem. So they hold for every t ∈ [0, W_B], including t = A.
   - This is the note's sign argument, rephrased.
   - The helicity-4 condition is no longer saturated at (4, g+4): it now has strict margin W_B − A > 0.
3. **Light scalars.** R_{n,0}(0, Δ) = 1, and W_0 ≥ 17/12 > 0 for every real Δ. [V]
4. **Light spin-2 window.** [V]
   - For a non-conserved spin-2 state in 3d, unitarity gives Δ > 3. The lightness condition gives Δ < g+2, and Δ ≥ g+2 is heavy and already covered. So the window is T = (Δ−1)(Δ−2) ∈ (2, g(g+1)).
   - Only k ≤ 2 contributes, so with ν = n², R_{n,0}(2, Δ) = F_2 = 1 − 6ν/T + 6ν(ν−1)/(T(T−2)). Its denominators are positive, since T > 0 and T − 2 = Δ(Δ−3) > 0.
   - Exact decomposition:
     3T(T−2)·(W_0 + A F_2) = (T − (6−3√2)ν)²·(T + (3+3√2)ν) + d(T),
     d(T) = (5/2)T² − (18(√2−1)ν + 9)T + 54(√2−1)ν².
     The cubic is A's tangency polynomial: a³ + 9(√2−1)(a² − 6a + 6) = (a − 6 + 3√2)²(a + 3 + 3√2).
   - The discriminant of d is (1512 − 1188√2)ν² + 324(√2−1)ν + 81. Its leading coefficient is ≈ −168.1, and its only positive root is ν ≈ 1.200. So for ν ≥ 4 (n ≥ 2), d(T) > 0 for every real T.
   - Hence W_0 + A F_2 > 0 on the **whole half-line T > 2**, uniformly in g. This holds regardless of the upper end of the window.
   - The relative margin is O(1/n²):
     - Exact critical-point evaluation gives n² × margin → 3.1297 as n → ∞.
     - The minimum over 2×10⁴ random (g, Δ) with g ∈ [4, 2000] is n² × margin = 1.60, at g ≈ 5.4 (n = 2).
     - This matches the note's 1.3e-7 at g = 10⁴ (3.13/5000² = 1.25e-7).
   - So the asymptotic tangency at Δ ≈ √(6−3√2)·n ≈ 0.663g is real, and the scale A cannot be enlarged at leading order.
5. **The 994.4 constant.** [V]
   - L = A n²(n²−1)/400, so τ ≤ 80√3/(√A · n√(n²−1)).
   - With n = g/2 + O(1) this is (640/√(√2−1))·g⁻³(1 + O(1/g)) = 994.41534…/g³. The note's 994.415… is correct.
   - Caveat on reading this number: it is a g → ∞ coefficient. At finite g, τg³ from this bound is 1014.9 at g = 10 and 994.61 at g = 100. When g sits just below an even integer (n = ⌊g/2⌋ is smallest relative to g) it is much larger: 1753 at g = 11.999 and 1055 at g = 101.999. The 960 bound of the original dual has the same floor effect.
6. **Still not covered: light multi-traces.** The repair covers light **single traces of spin 0 and 2** only, which are b=0 and J ≤ 2. Light multi-traces enter both ladders at every k, and they need the separate E_light input described in section 8.

**Grade.** The repair is correct. With items 1–4 it is now an exact statement: the rescaled dual is feasible on heavy ∪ light-scalar ∪ light-spin-2 states for all real g ≥ 4. It costs 960 → 994.4 in the asymptotic coefficient.

## Minor issues (none affects correctness)

1. The terse step in (b0weighted) and the indeterminate p_{−1}^{(1)} convention are discussed in items 2 and 1 above.
2. Prop. 2.2's proof ("the coefficient of each P_k vanishes") is very compressed, but it is correct and complete once the basis argument is spelled out.
3. **Provenance.** The `source_file_hashes` in `verification_summary.txt` (e.g. lp_dual.py `c605…`, 1716 bytes) do not match the files at commit 01b01e6 (`18aa…`, 1715 bytes). They equal the SHA-256 of each file **with one extra trailing newline**. The content is otherwise identical. This is cosmetic, but the hashes as written do not verify against the pinned commit.
4. Remark 7.1 (ρ ≤ 0.50496) is only a g → ∞ necessary condition at fixed x. I checked its algebra: the limits of x = 2 and x = 3, and the roots of 364ρ⁴ − 234ρ² + 36. It is consistent with the finite-g feasibility up to n ≈ 0.52–0.53g that `v7` found, because W_4/W_B > 1 and n²/Y < ρ² at finite g. "Cannot be extended to 0.51" is an asymptotic statement and should be read that way.
5. The abstract's "No … floating-point sign test" is accurate. The only transcendental quantity, h_5, enters only through the rational bound h_5 ≤ 4/15.

## Appendix: core certificate reconstruction (independent of the author's code)

```python
import sympy as S
n,Y,nu,h,H,u,t=S.symbols('n Y nu h H u t')
def yx(x,s=9):                      # finite 4F3, symbolic nu=n^2 and Y=N(N+9)
    tot,c=0,S.Integer(1)
    for k in range(1,x+1):
        i=k-1; c*=((i**2-nu)/S.Integer(k*k))*(x*(x+s)-i*(i+s))/(Y-i*(i+s)); tot+=c
    return S.cancel(1+tot)
y4,y5,y6=yx(4),yx(5),yx(6)
p=lambda j:S.Rational((j+1)*(j+9),2*j+10); q=lambda j:(2*j+9)*nu/(Y-j*(j+9))
E5=y5**2+h*p(4)*p(5)/q(5)*(y5-y4)*(y6-y5)             # Sonin energy from its definition
D=Y*(Y-10)*(Y-22)*(Y-36)*(Y-52); w=4*n**2+18*n+S.Rational(25,2)
F=S.cancel(((S.Rational(2,3)*((Y+S.Rational(25,2))/w)**2-E5.subs(h,H/15)).subs(nu,n**2))*1140*D**2*w**2)
for HH in (0,4):
    P=S.Poly(S.expand(F.subs(H,HH).subs(Y,4*n**2+30*n+36+t).subs(n,u+1)),u,t)
    print(HH,len(P.terms()),min(P.coeffs()))   # -> 169 terms, all >0, min 760
```
