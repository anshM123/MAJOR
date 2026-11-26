# Angular-action sliver: coupled-cap norming certificate

## Scope and upstream hypotheses

This is a computer-assisted sign certificate for the remaining sliver at Δ_g=1000, with b=0 and b=4, a≥998 and the integer physical spins specified in the angular-action identity. It is conditional on the retained angular-action identity, the corrected spectral transport and coefficient inventory, and the upstream coherent-error and ordinary-tail estimates in the full-band reduction package and the repository. It does not rederive the angular identity, prove a different gap parameter, or establish unrelated assumptions of the intended PRL paper.

The authenticated file `coeff_compact60.txt` contains 17,201 pairs of coefficient enclosures at scale 10^60. Every radius is at most two integer units. The original 3,000-digit table has SHA256 `0d09c4ad62d098ad19e8e75fda06e42776e9a5b571e51653768c0b567624754a`; `verify_compact.py` checks enclosure containment for every entry. The compact table SHA256 is `811e15cdd95319309edf64f5cf862bb674f43fbc73f551a002577c4b56b1ee06`. The inherited action-level coherent error is less than 10^-14 on c≤0.9; the inherited ordinary tail beyond degree 17199 is less than 2^-140. The additional tail from shortening each prefix is evaluated explicitly with the compact enclosures.

The earlier independent repository checkpoint is `ads4/ANGULAR_ROUTE.md`, §§115–126, and `ads4/code/gap_prove/ledger.json` on branch `claude/eloquent-cori-nftfjy`. At inspection the latter file had SHA `a4fa33817e148bb246fa35edf2da6947d10448cd`: 14 contiguous positive rows for b=0 to c=.585 and 13 for b=4 to c=.56, both starting at c=.28. A planned check, an isolated passing sample, or a running process is not counted as coverage.

## 1. Rational reduction and physical cap

Let s=2b+1, q=b(b+1), p_j=j(j+s), n=J−b, N=n+a, h=1/[N(N+s)] and c=n(n+s)h. Define

    U_k(c,h)=product_{j<k} (c-p_j h)/(1-p_j h)
    T_{k,b}=((k+1)^2 A_{k+1}+[-k^2+(8-4b)k+9-2b(b-1)+2q] A_k)/2 + sigma_b G_k
    sigma_0=1, sigma_4=-1
    F_K(c,h)=(1-c)^3 sum_{k=0}^K (A_k+h T_{k,b})U_k(c,h).

The normalized physical action is F_K/(1+qh)^4 plus the charged tail and coherent errors. On physical states, U_k is between 0 and c^k up to termination, and zero beyond termination. Thus replacing degree 17199 by K costs at most the explicitly summed absolute coefficient tail using c≤u and h≤Hl; the normalization only decreases the bound.

The derivative of sqrt(t(t+s)) for t>0 is (2t+s)/(2 sqrt(t(t+s)))≥1. Hence

    998≤a=N−n≤sqrt(N(N+s))−sqrt(n(n+s))=(1−sqrt(c))/sqrt(h).

This gives h≤H_*(c)=(1−sqrt(c))^2/998^2. The cap is decreasing and convex. For a cell [l,u], let Hl,Hu be outward-rounded endpoint caps and H(c) their affine secant. Then all physical h are in [0,H(c)]. We work on the coupled rectangle (c,y), y=h/H(c) in [0,1]; the larger constant-h rectangle is not used. Checks at the four corners verify positive denominators and |(c−p_jh)/(1−p_jh)|≤1. The inequalities before division are affine in c at fixed y and affine in y at fixed c; their extreme values occur at corners. Since p_j increases with j, endpoint indices suffice for the pole and ratio check. This includes both b channels and h=0.

## 2. Independent polynomial approximation errors

For fixed real y, each ladder factor is a real Möbius map of c. For a disk centered on the real c midpoint with positive real endpoint denominators, the image disk has real center, so its largest modulus occurs at one of the two real endpoints. A Bernstein ellipse of parameter rho>1 lies in the disk of relative radius (rho+rho^-1)/2. At the two extended real c endpoints, maximizing the modulus in the other real variable y also reduces to endpoints, because its denominator has positive sign. This yields a product majorant M_k(rho) for each summand, including separately maximized prefactors. A degree-m Chebyshev truncation has sup-norm error at most 2 M_k(rho) rho^-m/(rho−1). Different summands can use different rho values; their sum is still degree m. The accepted positive sums use exact rational corner computations and upward-rounded decimal arithmetic. Floating point only selects the degrees/radii; its estimates are not proof data.

The y approximation uses the partial-fraction representation, with possible simple poles h=1/p_j. For a ladder factor product, let v(v+s)=c p_j. Then the magnitude of its pole coefficient in the convention 1/(1-p_jh) is

    |rho_(j,k)| = ((2j+s)/p_j) |(-v)_k(v+s)_k| / [j!(k-1-j)!(j+s)_k],  1≤j<k.

The integer-v case vanishes. Otherwise the infinite-k limit is at most (2/j)c^(-(s−1)): use the gamma-ratio limit, reflection, v≥cj, and the finite product of (j+r)/(v+r) for 1≤r<s. The ratio of the limiting magnitude to the finite one is the product over ell≥k of (p_ell−cp_j)/(p_ell−p_j). Writing p_ell=(ell+s/2)^2−s^2/4, expand its logarithm as nonnegative powers. Retaining powers 1 to 4, replacing c by u, and comparing each decreasing power sum to its integral gives the lower bound

    sum_{r=1}^4 [(p_j+s^2/4)^r−(u p_j+s^2/4)^r]/[r(2r−1)(K+s/2)^(2r−1)].

This increases with j and decreases with K. Therefore division of the limiting magnitude by the exponential of this lower bound gives a valid bound on every finite residue. The coefficient sum uses suffix maxima of |A_k| and |T_k| and at most K−j terms. Grouping 16 j indices gives the conservative bound used in `full_preflight.py`; it uses an upper rational log(2) and a lower square-root enclosure.

For g=sqrt(1−p_j H) and r=(1−g)/(1+g), the Chebyshev expansion of 1/(1−p_jHy) has tail at degree m bounded by 2 r^(m+1)/[g(1−r)]. Use Hl instead of H(c), sum the residue bounds, retain the affine polynomial part exactly and bound (1−c)^3 by one. This proves the uniform y error Ey for every fixed real c. No joint polynomiality is presumed. The c and y approximation errors are rechecked from the archived coefficient data.

## 3. Degree-matched trigonometric norming

For real polynomials P of degree at most m on [0,1], write p(theta)=P((1+cos theta)/2). This is an even trigonometric polynomial of degree m. Take the 2m+1 equally spaced circle points theta_j=2 pi j/(2m+1). Dirichlet interpolation is exact. Its absolute interpolation norm is bounded by

    Lambda_m = 1+2 sum_{k=1}^m 1/(2k−1).

Indeed the nearest node contributes at most one after division by 2m+1, and the others in pairs have circular distances at least (2k−1)pi/(2m+1). The Dirichlet kernel is bounded by 1/|sin(delta/2)|, and sin(delta/2)≥delta/pi for delta in [0,pi]. Evenness means only m+1 distinct cosine points need evaluation.

If these nodes are perturbed by at most eta in [0,1], Markov's inequality gives ||P'||_infinity≤2m^2 ||P||_infinity. Hence the norming constant at the perturbed nodes is

    C_m(eta)=Lambda_m/[1−2m^2 eta Lambda_m],

provided the denominator is positive. Pi is enclosed by exact Machin-series brackets; cosine nodes use rational Taylor polynomials with remainder <10^-40 and rounding to denominator 10^24. The difference from the ideal node is <2·10^-24. The actual c node is rounded to the interior dyadic grid with denominator 2^61, giving eta_c=2·10^-24+2^-61/(u−l). The actual twist is floor(2^61 H(c) y)/2^61; therefore its y-coordinate differs from the ideal by at most eta_y=2·10^-24+2^-61/Hu.

Choose any rational reference polynomial R(x,y) of degree at most two in each coordinate, x=(c−l)/(u−l). The specific R is fitted to nine numerical values for convenience, but no accuracy assumption on that fit is used. Its exact rational degree-two tensor Bernstein coefficients bound it from below. The sum of absolute monomial coefficients is checked ≤10^-3. Subtracting R preserves either one-axis approximation error, because R has degree ≤m.

### 3a. Full grid

Suppose delta bounds the absolute residual F_K−R on all rounded nodes, including all arithmetic errors. For a function approximable by a degree-m polynomial within E, norming gives ||f||≤C_m(delta+E)+E. First apply this in y at each c node, then in c at arbitrary fixed y. Thus

    M=||F_K−R||≤Cc[Cy delta+(Cy+1) Ey]+(Cc+1) Ec.

This is a continuum result, not a sampling argument without interpolation control.

### 3b. Grouped c norming

Partition the m_y+1 distinct ideal y nodes into groups. For each group g, the c approximation degree m_g and error E_g are computed over its enclosed ideal y range. This avoids using the largest c degree at every y. Let C_g=C_{m_g}(eta_c), Lambda_y=Lambda_{m_y}, gamma=2 m_y^2 eta_y, and delta_g bound the residual on that group's rounded scalar nodes. At a fixed c, a degree-m_y approximant to F_K−R has norm at most M+Ey, so Markov's inequality controls the rounding from ideal y to actual y' by

    |f(c,y)−f(c,y')| ≤ 2 Ey+gamma(M+Ey).

For the ideal y node in group g, c-norming now implies

    sup_c |f(c,y)| ≤ C_g[delta_g+2 Ey+gamma(M+Ey)+E_g]+E_g.

Applying exact-node y norming and moving the term proportional to M to the left gives

    T=max_g {C_g[delta_g+2 Ey+gamma Ey+E_g]+E_g}
    lambda=Lambda_y gamma max_g C_g
    M ≤ [Lambda_y(T+Ey)+Ey]/(1−lambda).

The verifier checks lambda<0.1 and that every ideal y node belongs to exactly one group. This resolves the dependence of actual y rounding on c. In particular it would be invalid simply to reverse the two norming orders without charging that dependence.

## 4. Directed scalar arithmetic

For each grid (or y group), take exact upper bounds R_j for the absolute ladder ratios over its real corners. Set v_0=2^100 and s_0=0; at each j round v_j R_j upwards to an integer, multiply by powers of two until it is in (2^99,2^100], and increase the integer shift s_j by the corresponding amount. The check requires each increment ≤2. The inductive bounds are

    2^s_k product_{j<k}|r_j|≤1,
    product_{j=i}^{k-1} 2^(s_(j+1)−s_j)|r_j|≤v_k/v_i≤2.

Represent the backward rational sum at k in scaled integer units 2^(P−s_k), with P=180 in the first two high cells and P=80 in the grouped scans. Truncate coefficient conversions, ratio division, and dressing multiplication toward zero. The ratio quotient is computed before multiplication by 2^increment, so its truncation error is <4 units; the A conversion contributes <1 and the doubled-T conversion plus multiplication contributes <2 units. Ten units per step is conservative. The first bound above means each local unit error contributes at most one unscaled unit at the result. The compact coefficient radii add a separate bound, yielding

    E_scalar ≤ (K+1)[4+Hl((K+1)^2+K^2+8K+67)]/10^60 + 10(K+1)/2^P.

The reference is evaluated by truncated tensor Horner at dyadic coordinate scale 2^50 and integer coefficient scale 10^20. Its derivatives are bounded by 2·10^-3 in either coordinate. Thus reference error is at most 4/(1000·2^50)+8/10^20. Final projected output bounds use floor and ceiling at 10^-20 scale. All of these errors are included in delta or delta_g.

The SIMD code uses base-2^52 signed limbs, not floating-point signed sums. The double-precision quotient estimate in each long-division step is corrected using exact integer remainder tests until 0≤remainder<denominator; failure to reach this invariant throws an error. Thus the returned quotient is exact regardless of estimate accuracy. The denominator and numerator magnitude are ≤2^61. Limb multiplication uses exact 52-bit multiply-add and carry; intermediates fit unsigned 64 bits. The partial-product bound ≤2 proves that a suffix maximum of the scaled A and doubled-T magnitudes plus ten, multiplied by 2(K+1), bounds every exact signed recurrence integer. The code chooses its limb capacity from this bound. A portable GMP implementation with precisely the same truncation order is provided for independent reproduction on machines without AVX512IFMA. Representative subsets are compared as exact integers, not as approximate floats.

## 5. Final sign and coverage test

For each channel, let B_R be the least degree-two tensor Bernstein coefficient of R, and M the certified residual bound. Then the physical normalized action has lower bound

    (B_R−M)/(1+q Hl)^4 − E_short_tail − 2^-140 − 10^-14,

provided B_R−M>0. The scalar recurrence and approximation proofs above apply to the continuous coupled cap, while the short-tail bound is used only on physical states. The sign test is exact rational arithmetic and must be positive for both b=0 and b=4. The independent checker takes the larger of any archived and freshly rounded upper error bounds; tiny upward-rounding differences cannot create a false pass.

The completed input grids, outputs, logs and exact rational ledgers are archived. `verify_merged.py` verifies their counts, input construction, cap and pole checks, analytic error bounds, norm constants, reference bounds, final inequalities and union of intervals. It does not rerun every scalar recurrence; the drivers and portable checker permit that separate reproduction.

## 6. Additional shortened group prefixes

In the final upper-intermediate cells, the scalar grid and the c-direction approximant may use a group-specific degree K_g≤K. This is distinct from the original physical-only tail bound. For the group enlarged by eta_y (to include both ideal and rounded nodes), use exact corner maxima R_j for the absolute rational ladder factors. Then on this entire real group

    |F_K−F_(K_g)| ≤ sum_{k=K_g+1}^K (|A_k|+Hmax |T_k|) product_{j<k} R_j.

The prefactor (1−c)^3≤1 has been discarded. `short_plan.py` bounds these nonnegative terms using upward-rounded products and sums at 50 decimal digits. Every chosen group remainder is <10^-10. The c approximation error is increased by this remainder, and the node residual bound is also increased by it. The full K y approximation is unchanged; consequently the self-consistent rounding argument still applies to F_K itself. `verify_merged.py` recomputes the group tail and uses the larger of its saved and freshly rounded values. No tail estimate relying on positivity or termination is used on the enlarged nonphysical group.

## 7. Executed theorem and final interval union

**FULL_SLIVER_PASS_CONDITIONAL.** The exact merged ledger covers every c in [129/500, 9/10] for both b=0 and b=4, at Δ_g=1000, a≥998 and all physical integer J≥1001. This proves the requested 0.3<L<2.3 in both repository conventions, conditional on the retained angular-action identity, authenticated spectral coefficient enclosures, and retained coherent/ordinary-tail bounds. No remaining interval is untested in the continuum certificate.

The common normalized action lower bound is **5.140000E-6**. Channel-specific minima, rounded downward for display, are 5.166280E-6 (b=0) and 5.148037E-6 (b=4). The 11 new cells contain 5,669,470 directed scalar evaluations across both channels, and all 22 new channel-cell inequalities passed the separate analytic/ledger checker. Seven earlier local cells were rechecked and merged.

The last inspected repository ledger (`ads4/code/gap_prove/ledger.json`, branch `claude/eloquent-cori-nftfjy`, blob SHA `a4fa33817e148bb246fa35edf2da6947d10448cd`) records common coverage from c=.28 to c=.56. The present local calculation extends beyond that endpoint; no changes were made to the GitHub repository.

| c interval (exact decimal endpoints) | b=0 lower bound | b=4 lower bound |
|---|---:|---:|
| [0.2580, 0.2680] | 1.380437E-5 | 1.322035E-5 |
| [0.2680, 0.2800] | 1.355361E-5 | 1.299169E-5 |
| [0.2800, 0.2900] | 1.351826E-5 | 1.297521E-5 |
| [0.2900, 0.3000] | 1.339302E-5 | 1.286763E-5 |
| [0.3000, 0.3500] | 1.050459E-5 | 1.026366E-5 |
| [0.3500, 0.4462] | 8.358087E-6 | 8.362481E-6 |
| [0.4462, 0.5305] | 7.981095E-6 | 7.972298E-6 |
| [0.5305, 0.6041] | 1.043479E-5 | 1.026382E-5 |
| [0.6041, 0.6677] | 9.575659E-6 | 9.449694E-6 |
| [0.6677, 0.7187] | 8.819155E-6 | 8.723942E-6 |
| [0.7187, 0.7596] | 8.152688E-6 | 8.079011E-6 |
| [0.7596, 0.7929] | 7.543729E-6 | 7.485170E-6 |
| [0.7929, 0.8204] | 7.010982E-6 | 6.964830E-6 |
| [0.8204, 0.8434] | 6.529089E-6 | 6.491804E-6 |
| [0.8434, 0.8628] | 6.092176E-6 | 6.062101E-6 |
| [0.8628, 0.8794] | 5.690831E-6 | 5.666348E-6 |
| [0.8794, 0.8936] | 5.338072E-6 | 5.318014E-6 |
| [0.8936, 0.9000] | 5.166280E-6 | 5.148037E-6 |

The table uses downward-rounded lower bounds; the full rational numbers are in `sliver_merge/FULL_COVERAGE_LEDGER.json`. The coverage proof uses exact rational endpoints, not the displayed decimals. Since e^(0.3−20/998²)>1/(1−.258) and e^2.3<10 (checked by rational series bounds), the closed c interval covers the target even if L is defined with the physical w correction log(1+qh).

## Reproducibility and limits

Start with `REPRODUCTION.md`. The mathematical norming and directed-arithmetic proof is in `sliver_merge/PROOF.md`. `python check_ledgers.py` recomputes the analytic bounds and exact sign/coverage ledger, and `python replay_scalar.py --jobs 4` independently reruns every scalar shard using portable GMP. The executed scans use corrected exact integer SIMD division; seven representative node subsets have been compared against GMP as exact integers. The entire slower portable replay has not been rerun here. The original 3,000-digit source table was authenticated and all compact enclosures verified; raw moment regeneration and the angular-action identity remain retained dependencies. See `authenticate_original.py` to check provenance against `coeff_long_17200.txt` in the prior full-band reduction archive. Historical subdirectory READMEs describe their own earlier checkpoints, not the status of this assembled certificate.
