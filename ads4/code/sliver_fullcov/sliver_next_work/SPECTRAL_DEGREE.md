# Finite spectral-degree reduction for the remaining AdS4 angular sliver

**Scope.** Conditional on the upstream raw angular identity, physical state parametrization, and retained repaired spectral mass/pointwise envelopes in the preceding certificates. No sign claim follows solely from the reduction. The exact `full_grid.py` ledger covers c=0.268 to 0.9, which contains all the remaining target L<2.3. A completed signed lower bound for each channel and each rectangle, or a covering refinement, is still required.

## 1. Gaussian data and correlated seed error

Use the exact binary64 repaired inventory (3564 entries, original order). With U=1200, let normalized Gaussian-cutoff moments be m_k = integral of G theta (nu/U)^k over the real line, with the c=0 contribution halved in the final even-moment sum. The integration-by-parts formula in the preceding package is

    m_(k+1) = (center/U)m_k + k(sigma/U)^2 m_(k-1)
               + (sigma^2/U)(e^-_k - e^+_k).

The boundary e sequences have the mass, mean and variance obtained by completing the square (see preceding README). `gaussian_seed.cpp` bounds m_0 on decimal scale 10^800. All 3564 intervals have width <= 10^20 integer units; replace m_0 by the exact midpoint rounded toward zero at that scale. Its difference from the exact m_0 is <= 10^-780. The choice changes the complete moment sequence by the corresponding *ordinary Gaussian* homogeneous moment solution, not independent per-order uncertainties.

`gaussian_fixed.cpp` computes the moments of the substituted data on binary scale S=2^P, with P selected so S>10^(D+20)(6/5)^(2K), here K=6000, D=3000. For k<=12000 the sum of the absolute boundary recurrence coefficients is <0.93 (the mean lies in [-1099,1099]/1200), the main recurrence sum is <0.87, and the coupling coefficient is <0.001. The exact boundary sequences have magnitude <=1 by induction and the main substituted moments have magnitude <=3. Each boundary initial midpoint and rounded mean differs by at most 2/S; if a seed rounds to zero the exact boundary sequence is bounded by 2/S thereafter. Each fixed-point product truncates at most one unit. In the boundary step a rounded mean costs at most four further units (use approximate moment magnitude <=2); 100*.93+4+2<100. In the main step the two boundary errors and three products cost at most 200*.001+3 units; 100*.87+200*.001+3<100. This gives absolute error <=100/S at every moment in each sequence. The weighted shard error is <(100*1220000+3564)/S<10^10/S. Conversion to decimal scale and normalization nu/1000 use exact rational multiplication and directed integer endpoint division. The 12 shards jointly produced 6001 enclosed moments at D=3000, max interval width 24 integer units.

For the seed discrepancy, let delta=1220000*10^-780. The aggregate homogeneous signed Gaussian mass is <=delta, and its absolute exp(B|nu|)-weighted integral is <=2 delta exp(B*1029+B^2*.901^2/2), B=11/7. The rational raw factor satisfies |R| <= sum |r_a|/a, a=1,9,49; the real-frequency raw angular action bound from `raw_seed_budget_full.py` is <10^19 exp(B|nu|), and we allow factor 2 for the additional qh. Holding F_(0,a) fixed adds a pole atom with mass <=delta/a at nu^2=-a/4. The exact inequalities using the imaginary-pole action bounds in `pole_bound_full_results.json` are checked by `verify_new_seeds.py`. The combined new Gaussian seed error in normalized action is <10^-45 (real contribution bound <9.01e-53). The older uncut rational-pole seed and log-seed approximations still have their separately retained correlated budgets; for c<=.56 the combined old/new seed errors and spectral truncation are <10^-18, and for c<=.9 use 10^-14 conservatively. They are **not** asserted as 3000-digit intervals for the true unperturbed moments.

## 2. Spectral truncation with a positive coefficient majorant

Write P(z,lambda)=P_{-1/2+i nu}((1+z)/(1-z)), lambda=nu^2+1/4, and f=(1-P)/lambda (continued at zero). Legendre's ODE gives (zP')'=-lambda P/(1-z)^2, and consequently

    (zf')' = (1-lambda f)/(1-z)^2,
    m^2 f_m = m - lambda sum_{1<=l<m} (m-l) f_l.

In u=-lambda, every coefficient of f in z and u is nonnegative, by induction on m; P=1+u f has the same property. These are finite spectral polynomials at each z degree. For either generating function, deleting terms of degree >M=6000 in nu^2 costs, in the coefficientwise absolute z-series evaluated at r=.91, at most

    ((|nu^2|+1/4)/tau^2)^(M+1) exp(beta*tau),

when |nu^2|+1/4<tau^2, tau=3150, beta=19/5. Indeed all powers lambda^j with j<=M survive unchanged and the sum of the absolute discarded terms of each higher power is bounded by (|nu^2|+1/4)^j. Positivity allows comparison to the series evaluated at u=tau^2. The Laplace integral for the Legendre function at z=r and u=tau^2 yields P<=exp(beta*tau), since the positive Legendre exponent sqrt(tau^2+1/4)-1/2<=tau, and cosh(beta)>(1+r)/(1-r). For f use (P-1)/tau^2<=P. For real nu>1200, the unconditional bound for the entire absolute spectral expansion (and hence the discarded part) is exp(beta*(nu+1)), by taking u=nu^2+1/4. This last bound also applies to f since u>1 in this range.

Let m=M+1, V=1030, t=2m/V. For nonnegative real nu,

    (nu^2+1/4)^m <= (nu+1)^(2m) <= V^(2m) exp(t*(nu+1-V)).

The raw weight factor |R|/Lambda<1; the absolute real Gaussian envelope includes both signs of the center and the new coherent seed variation, bounded conservatively by Gaussian amplitude sum 2*1220001, width <=.901, center <=1029. Its integral in the last exponential is <=2*1220001*2.27*exp(.901^2*t^2/2). This uses the Gaussian MGF and sqrt(2pi)*.901<2.27. The raw log part, including the small K0 seed correction and dropping theta, is bounded by 10(1+nu)^8 exp(-nu^2/2): |y_i|<1, |Rt|/Lambda<1, log(1+nu^2)<=nu^2, and the homogeneous log seed correction has absolute normalized Gaussian multiplier <10^-40. After using (1+nu)^8<=exp(8nu), its weighted integral is <=30, since (t+8)^2/2-t*(V-1)<0. The three rational pole seed corrections and log corrections are represented by atoms at nu^2=-a/4; their total absolute mass is <1 by the earlier cutoff bound 7.1e-24, numerical seed widths <1e-40, the new seed budget, and the finite rational coefficients. For all such atoms |nu^2|+.25<=V^2. The ideal psi3 measure's weighted L1 mass is <6e7, independently checked from its exact Bernstein coefficients (actual upper <1.34e7); its support is |nu|<=sqrt(997002)<999. The far real tail beyond 1200 is bounded by tangents to the two concave quadratic Gaussian exponents in `spectral_cut.py`.

## 3. Transfer of this spectral error through the action

The inherited raw coefficient identities are

    Z_k=12 f_(k+2)-(k+1)f_(k+1)+(k+3)f_(k+3),
    A_k=binom(k+3,3)^2 [Z_k-k Z_(k-1)/(k+3)]/(k+3).

At r=.91 the spectral error E in the l1 norm of f's coefficients gives |delta A_k| <= (k+4)^6 r^(-k-3) E: bound Z by (2k+16)r^(-k-3)E, use the analogous bound for Z_(k-1), and bound the binomial factor by (k+3)^5/36. The dressing |delta B_k| <= (k+5)^8 r^(-k-4) E follows directly from the two coefficient terms (|( -k^2+(8-4b)k+9-2b(b-1))| <= (k+5)^2 for b=0,4). The physical ladder U_k satisfies 0<=U_k<=c^k before termination, zero afterward, on physical states. With H=998^-2 and c<=.9, the normalized action error is bounded by the raw operator constant

    C_raw = (1+20H)*4^6*6!/r^3/(1-.9/r)^7
             + H*5^8*8!/r^4/(1-.9/r)^9.

Here (k+a)^p <= a^p*p!*binom(k+p,p) for a,p positive integer. For the ideal (1-z)^-4 P coefficient series, the constant H/(1-r)^4 suffices. The normalization d^3/(1+qh)^4 is <=1. Both estimates hold equally for any prefix. `spectral_cut.py` checks the product of the total operator bound and the integrated spectral tail by exact rational inequalities: error <2^-1910. This is a separate allowance from the ordinary z-degree truncation.

## 4. Rectangular moment transport

Define I_(m,j) = integrated f_m with a further factor (nu/V)^2j. If T_j are the signed input moments, set T_j=0 for j>M. Define S_(m,j)=sum_(l=1)^m I_(l,j) and W_(m,j)=sum_(l=1)^m S_(l,j). Then

    4 m^2 I_(m,j)=4m T_j - 4V^2 W_(m-1,j+1) - W_(m-1,j).

This computes precisely the spectral-degree-M truncation at j=0 (and the corresponding shifted moments at higher j). Ascending j permits in-place updates, since W_(m-1,j+1) has not yet been changed. For the ideal series, P_m= -V^2 I_(m,1)-I_(m,0)/4, and four cumulative sums implement multiplication by (1-z)^-4. Uniform integer interval errors propagate by positive sums: E_I=ceil((4m E_T+(4V^2+1) E_W)/(4m^2))+1, E_S+=E_I, E_W+=E_S. Ideal P and cumulative sums carry the same directed allowances. `verify_long_transport.py` compares 236 outputs against independent exact rational finite spectral sums, and checks the 1253 old ideal coefficients overlap the corrected new enclosures. Note: an initial ideal transport run had an unsigned integer negation bug; its vector and resulting failed sign tests are invalid and superseded by the corrected code. Raw transport was unaffected. The valid transport uses explicit signed negation of the product, not an unsigned factor.

## 5. Ordinary truncation and coverage

`tight_spectral_tail.py` improves the earlier ordinary coefficient envelope to 2^2337 using the retained residual bands, Gaussian cutoff allowance, ideal mass estimates, and the retained Mehler bounds on |z|=.995. Each band uses its right frequency endpoint in the exponential weight. The low-frequency residual bound includes the ideal psi2 contribution (less than 22 in the normalized mass convention). The exact tail check is the same as the prior direct-prefix proof with this smaller envelope. `full_grid.py` checks 21 contiguous rectangles from c=.268 to .9, K<=16989, tail <2^-110, and every pole gap >.14469 using h <= (1-sqrt(c))^2/998^2 and s<=9. These rectangles are available for the subsequent sign computation, not sign claims. The coefficient file contains all needed A/G entries through index 17200. A physical sign test can also use a shorter K with ordinary tail <10^-7, as the small-c tests do, provided that error is subtracted from the final lower bound. 
