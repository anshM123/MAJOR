# Sliver continuation: coupled norming and direct low-L prefix reduction

This archive is **not** a full-sliver certificate. The inherited angular-action identity, EPSW repaired spectral data and cutoff estimates, LRINF reference floor, low-frequency residual theorem, base-multiplier and localization estimates are inputs. The new low-L theorem certifies the truncation and pole exclusion, not the sign of the rational prefix. The moment-transform test uses the ideal NF02 weights only.

## Reproduce the next analytic slab

Requires Python 3, C++17 and GMP. `boxes_2.6_2.9.json` gives six coupled parameter boxes (two channels apiece). The full signed transfer is enclosed by the sequential one-variable Chebyshev norming checker from the previous archive, with exact rational scalar nodes and directed integer rounding. `results_2.6_2.9.json` and `run_2.6_2.9.log` preserve the actual enclosures. `p3_2.6_2.9_budget.json` contains the residual-only budget, not a full sign claim. The initial p=4 contour at height 11.25 failed; the higher contour results are retained separately.

```
g++ -O3 -std=c++17 scalar.cpp -lgmpxx -lgmp -o scalar
python check_norm.py boxes_2.6_2.9.json --jobs 4 --max-nodes 300000 --output results_2.6_2.9.json
python resid_budget.py results_2.6_2.9.json 0.0000025
for b in 0 4; do
  python ellipse_extend.py results_2.6_2.9.json 10.25 ct_2.6_2.9_3_$b.json 3 $b
  python ellipse_extend.py results_2.6_2.9.json 19.25 ct_2.6_2.9_4_${b}_high.json 4 $b
done
python assemble_26.py
```

The retained floor 2.5e-6 is conservative: the repository's 16 contiguous `lrt3_rcert.json` beta cells on [3.9,4.3] have minimum 2.9539369377261934e-6, and contain the slab's beta range. Input masses are inherited rather than recomputed. The previous package documents the norming and ellipse-extension proofs.

## New low-L direct prefix theorem

Set c=1-d, h=1/[N(N+s)], s=2b+1, q=b(b+1), and tau=h/(1+q h). For a>=998, the exact physical state identity implies h<= (1-sqrt(c))^2/998^2. The physical ladder U_k is zero for k>n, with 0<=U_k<=c^k before termination. Use the rational product for the prefix; every denominator is bounded below on each entire rectangle, so the expression also gives the zero terms at and beyond physical termination (where applicable).

The retained angular-action identity is H_ang=w^3 L(G_A+tau D_b G_A+tau sigma_b G_D), where w=d/(1+qh) and the conical generating functions are those of FULLPOS3. On |z|=.99 take a Cauchy differentiation radius .005. The Mehler estimate on |z|<=.995 is |F_(p,nu)(z)|<=(pi/2)*200^p*exp(pi nu/2). The exact dressing costs a factor <=1+tau(2/.005^2+5/.005+8)<2. For nu<=1100 use the retained mass bounds (2/Lambda) int |f_theta| < 60000001 and tau Lambda^2 int |psi3| < 60000000. These follow from the NF02 coefficient norms <3e7 and EPSW/cutoff, with tau Lambda^2<2. The integrated dressed generating function on the Cauchy circle is bounded by 2^2554: 4*200^4*120000001<2^60, and exp(pi*1100/2)<2^2494, both checked. For nu>1100, the retained repaired Gaussian envelope, with the grouped exact binary64 coefficient ceilings from the preceding package, bounds the exp(pi nu/2)-weighted tail. Use pi/2<11/7. The high-frequency prefix error is less than 12*200^4*gaussian_mass(1100,11/7).

For a cell c in [l,u], choose K so that (u/.99)^(K+1)<=2^-2670. The omitted low-frequency ladder contribution is bounded by 2^2554*(u/.99)^(K+1)/(1-u/.99). Since u<=.56, the geometric denominator factor is <3. Combining with the high-frequency contribution gives error <2^-110, with the physical normalization w^3<=1. The pole-free condition h_cap(l)*(K-1)*(K-1+s)<1 is checked with exact rational arithmetic. `low_prefix.py` selects six cells, verifies the constants and the entire cover, and writes the actual bounds. This provides a finite rational sign problem, not a positivity proof.

The common numerator after clearing the positive factors is

  d^3 sum_(k=0)^K [(1+q h) A_k + h B_k + sigma_b h D_k]
      prod_(j<k)(c-h j(j+s)) prod_(k<=j<K)(1-h j(j+s)),

with denominator (1+q h)^4 prod_(j<K)(1-h j(j+s)). Here A_k=[z^k]G_A, D_k=[z^k]G_D and

  B_k = ((k+1)^2 A_(k+1) + [-k^2+(8-4b)k+9-2b(b-1)] A_k)/2.

Thus a bivariate rational positivity checker with a strict 2^-110 action margin would close this low-L interval. There is no 1000th-order h expansion or 3000-jet computation at each anchor. The exact globally integrated realized coefficients remain to be generated, and the numerator remains to be certified over its rectangles.

## New global moment transform (executed for ideal weights)

For p=3,4 let I_(k,j)=int_0^1 g_(p,k)(Lambda t) psi_(p-1)(t) t^(2j) dt. The conical Taylor recurrence gives exactly

  4(k+1)^2 I_(k+1,j) = (8k^2+8pk+4p-1)I_(k,j)
                           -4 Lambda^2 I_(k,j+1)-4(k-1+p)^2 I_(k-1,j).

I_(-1,j)=0 and I_(0,j) is an exact rational spectral moment. This triangular recurrence requires O(K^2) small-integer operations on fixed-point integers, not expansion/contraction of high-degree rational CDH polynomials. For each row carry a common l-infinity absolute error E_k in fixed-point units: E_(k+1)=ceil(((8k^2+8pk+4p-1+4Lambda^2)E_k+4(k-1+p)^2 E_(k-1))/(4(k+1)^2))+1. The midpoint division is toward zero. General initial moment balls are allowed by `moment_transform.cpp`. The integrated generating coefficients A_k=2 I_(k,0) (p=3, ideal) and D_k=Lambda^2 I_(k,0) (p=4). For the realized A sequence the input moments must instead be those of f_theta on [0,infinity), in the t=nu/Lambda coordinate; use the same prefactor 2.

```
g++ -O3 -std=c++17 moment_transform.cpp -lgmpxx -lgmp -o moment_transform
python test_moment_transform.py 3249 2500
python low_prefix.py
```

The ideal transform was executed for 3,250 entries per channel at 2,500 decimal places; 15 exact-rational independent low-degree checks per p passed, and 10 pointwise ideal-action checks passed. The log's displayed zero error is floating underflow only: all sign decisions used exact integer/rational errors, not zero. The actual ideal moment error for every transported entry is less than 1e-463 before the coefficient prefactors. This is not a continuum ideal-action certificate, nor a repaired-density certificate.

## Suggested exact realized moment acceleration (derived, not executed)

For a Gaussian times the real erf cutoff theta(t)=[erf((t+nu_c)/10)-erf((t-nu_c)/10)]/2, define M_m=int_R t^m exp(-(t-c)^2/(2 sigma^2)) theta(t) dt. The exact generating function is

  sqrt(2 pi) sigma exp(c z+sigma^2 z^2/2)
      * [erf((c+sigma^2 z+nu_c)/sqrt(100+2 sigma^2))
         -erf((c+sigma^2 z-nu_c)/sqrt(100+2 sigma^2))]/2.

Integration by parts gives M_(m+1)=c M_m + m sigma^2 M_(m-1) + sigma^2 int t^m G theta'(t)dt. Since theta' is the difference of two Gaussians divided by 10 sqrt(pi), the final term is exactly the difference of two ordinary Gaussian moments (complete the square, with variance 100 sigma^2/(100+2 sigma^2)). Thus no high-order cutoff quadrature is needed for the polynomial Gaussian component. For the unshifted rational weight in fp_moments, polynomial division leaves only finitely many rational-pole base integrals; the logarithmic component needs its own certified moments. The equivalence to the realized generating coefficients and the rigorous errors of all these inputs must still be implemented/audited before use. In particular one must not simply apply the ideal moment radii to the realized weight.

## Executed additional positivity slab

`ledger_26.json` and `ledger_26.log` record exact rational comparisons, conditional on the retained inputs. All twelve channel boxes (1,057,106 scalar nodes in total) were enclosed. At the higher p=4 contour height 20.25, the formerly failing top-line estimates are below 7.41e-18. The total contour bounds (including p=3) are below 9.80e-10 and 1.487e-9, respectively. The final E/F ratios are <0.091148 and <0.112364 on L in [2.6,2.9]. This joins the preceding [2.9,3.65] ledger. The whole lower region is still open; specifically the direct prefix reduction does not certify its rational numerator. No claims have been added for uncomputed L intervals.
