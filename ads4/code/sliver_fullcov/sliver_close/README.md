# Coupled-cap sliver certificate, November 26, 2026

**PARTIAL, CONDITIONAL — NOT A COMPLETE PROOF.** The new signed calculation covers c in [0.30,0.35] in both channels, with normalized physical lower bounds 1.050459132145526e-5 and 1.0263660091837696e-5. In the direct-action convention L=-log(1-c), this extends the previous lower band to L=0.4307829160924542. The 13 further cells through c=.90 have certified approximation errors and pole exclusion only. Their scalar sign scans have not been executed. Coverage for L>=2.3 is retained, not independently recomputed here. The angular-action identity, corrected spectral transport, coherent seed/model error <1e-14, ordinary tail past degree 17199 <2^-140, and previously certified domains are inherited from the prior packages, especially sliver_full_band_reduction_2026-11-26.zip. The certificate uses the direct-action L convention; the alternative physical w convention differs by log(1+qh).

## Formula and coupled cap

Let p_j=j(j+s), s=2b+1, q=b(b+1), h=1/[N(N+s)], c=n(n+s)h. Let sigma_0=1 and sigma_4=-1. Define

    U_k=product_{j<k} (c-h p_j)/(1-h p_j)
    T_k=((k+1)^2 A_(k+1)+(-k^2+(8-4b)k+9-2b(b-1)+2q) A_k)/2 + sigma_b G_k
    F_K=(1-c)^3 sum_{k<=K} (A_k+h T_k) U_k.

The normalized physical action is F_K/(1+qh)^4 plus the inherited and shortened-prefix errors. On physical states U_k is nonnegative, <=c^k, and terminates at n. This justifies the shorter-tail bound using the actual enclosed coefficient table, not an assumption about unphysical continuation. The tail bound charges all coefficients through 17199 plus the inherited tail beyond.

The physical h cap is H_*(c)=(1-sqrt(c))^2/998^2 by the mean-value inequality for sqrt(t(t+s)). H_* is convex. On [l,u], its affine secant through upward rounded endpoint caps Hl,Hu bounds it from above. Substitute h=H(c)y, y in [0,1], with affine H. For all real (c,y) in the cell, each denominator is positive and the absolute ladder ratio is at most one, as checked at the four corners: equivalently 1-p_jh>0 and -1+p_jh<=c-p_jh<=1-p_jh. The first and last conditions are affine in p_j for fixed c,h, so checking j=0,K-1 suffices. The checks apply to both s=1,9.

## Independent one-axis approximation bounds

Fix the other real variable. Each ladder factor is a real fractional-linear function of the variable under consideration. On a disk centered on the real axis containing no pole, its image is a disk with real center; the largest modulus occurs at one of the two real endpoints (constant maps are harmless). The normalized ellipse of radius parameter rho lies in the disk of radius (rho+1/rho)/2. At the two extended real endpoints, maximize again in the other real variable, also a fractional-linear function with positive denominator. This reduces factor bounds to four exact real corners. The prefactors (1-c)^3 and A_k+hT_k are bounded by their absolute separate maxima. Each summand may use a different rho, since its degree-m truncated Chebyshev series is itself a polynomial. Summing these truncations gives error bounded by the sum of 2M_k(rho_k)rho_k^(-m)/(rho_k-1), using the standard Chebyshev coefficient estimate for analytic functions on an ellipse. All arithmetic in the accepted positive sums is directed upwards. Floating-point planning only selects the radii and degrees; it is not a proof input.

For later twist cells, a stronger independent y approximation uses partial fractions from the previous residue certificate. For 1<=j<k, let v(v+s)=c p_j. The residue in the convention 1/(1-p_jh) is

    |rho_jk|=((2j+s)/p_j)|(-v)_k(v+s)_k|/[j!(k-1-j)!(j+s)_k].

Its infinite-k limit is <=(2/j)c^(-(s-1)) by reflection and the finite gamma-product bound; the integer-v case gives zero. The ratio of the limit to the finite residue is product_{ell>=k}(p_ell-cp_j)/(p_ell-p_j). For u>=c, p_ell=(ell+s/2)^2-s^2/4. Expanding the logarithms of the factors into nonnegative terms and comparing each decreasing power sum with its integral gives, for k<=K,

    log(product) >= sum_{r=1}^4 [(p_j+s^2/4)^r-(u p_j+s^2/4)^r]/[r(2r-1)(K+s/2)^(2r-1)].

The bound is increasing in j and decreasing in K. Use this damping in the residue sum, taking the coefficient-specific suffix maxima of |A_k| and |T_k|. Group 16 j indices at a time and use conservative extrema for both s. With g=sqrt(1-p_j Hl), r=(1-g)/(1+g), the exact Chebyshev Poisson expansion of the pole on y in [0,1] has tail at most 2r^(m+1)/(g(1-r)); using Hl gives a uniform upper bound for h=H(c)y. The affine partial-fraction part is retained exactly. Normalization (1-c)^3 can only decrease the bound. Rational log(2) and square-root enclosures and upward-rounded positive sums certify the final bound. This is an approximation in y for each fixed real c, not an assertion of joint polynomiality or positivity.

## Sign test

For a real polynomial p of degree m, the Bernstein-Szego inequality for p(cos(theta)) gives |p(theta_node)|>=||p||_infty cos(m delta) at an angular grid point within delta of an absolute maximizer, provided m delta<pi/2. The driver rigorously encloses pi using Machin's identity, computes cosine nodes by rational Taylor series, accounts for rational rounding by the square-root angular perturbation bound, and verifies m delta<1.13 with cos(1.13)>2/5. Thus C=5/2 is valid as norm inflation in each axis. For a real-valued function approximated by a polynomial within epsilon, its uniform absolute bound is C times its maximum node absolute value plus (C+1)epsilon. For each channel take R as the midpoint of the global node-value bounds, delta as their half-range (including arithmetic error). Apply the y norming bound at each c node, then c at arbitrary fixed y. It follows that

    F_K >= R-C^2 delta-C(C+1)epsilon_y-(C+1)epsilon_c.

Only if this is positive and the physical normalized lower bound after charging the shorter tail, inherited tail and coherent error remains positive is a channel declared PASS_CONDITIONAL. A rational scalar checker uses truncated signed integer Horner operations at scale 10^60, with exact directed 10^-20 output endpoints; real ladder ratio <=1 ensures errors do not amplify. The driver includes a conservative coefficient and operation error bound <10^-45. The compact table is verified against the original 3000-digit corrected table using exact integer enclosure containment (17201 pairs, radii <=2 at 10^-60 scale). The compact table SHA256 is 811e15cdd95319309edf64f5cf862bb674f43fbc73f551a002577c4b56b1ee06. The original long table SHA256 is 0d09c4ad62d098ad19e8e75fda06e42776e9a5b571e51653768c0b567624754a.

## Reproduction

Requires Python 3, NumPy, SciPy, C++17, GMP and a 64-bit unsigned long. From this directory:

    python verify_compact.py /path/to/coeff_long_17200.txt
    g++ -std=c++17 -O3 norm_scalar.cpp -lgmpxx -lgmp -o norm_scalar
    python norm_driver.py cell_30_35.json

The first command checks containment against the prior archive (optional if the compact table has already been authenticated). The executed first cell has 1,035,440 scalar evaluations and both signs pass; its four shard inputs/outputs/logs and ledger are included. The full_preflight.json ledger and the other 13 cell plans specify approximately 18.17 million further scalar evaluations. `python verify_all.py` recomputes the approximation and coverage preflight. `python run_all.py` runs the first cell and all 13 further sign scans, stopping on a failed comparison. This is an executable sufficient proof test, not a guarantee of success. The coarse scout_results.json contains positive samples only and cannot certify a continuum. No processes remain running. A complete proof requires actual positive ledgers for all cells and validation of inherited assumptions; these have not all been obtained here.
