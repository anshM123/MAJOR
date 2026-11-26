# Physical sliver continuation — six signed certificates and uniform twist compression

**Status: partial, conditional. Not a full-sliver proof.** Both channels are now covered on c in [.268,.30], overlapping the preceding certificate. The remaining physical gap is c in (.30,1-exp(-2.3)), or L in (-log(.70),2.3). No claim is made for the sign of the whole remaining band. The underlying spectral and angular identities and correlated seed allowances are retained from `sliver_full_band_reduction_2026-11-26.zip`.

## Executed sign checks

`new_sign_ledger.json` is the exact-rational verified ledger for the six successful direct checker runs. Each entry charges the physical truncation and coherent seed errors as in the preceding package. The `.in`, `.out`, `.log` and `.json` files are the actual per-channel checker inputs and results. `verify_sign_ledger.py` independently checks the coefficient transfer, pole gap, Taylor remainder, coverage and lower bound algebra. It does not re-execute the C++ recurrence; the commands below do that. The exact coefficient vector is the corrected `coeff_long_17200.txt` from the preceding archive; SHA256:

```
0d09c4ad62d098ad19e8e75fda06e42776e9a5b571e51653768c0b567624754a
```

The vector is not duplicated in this archive. Place it in the same directory after extracting the preceding archive. The retained band masses are copied here. The old 21-cell `full_grid_results.json` is a truncation/pole-gap ledger only.

## A uniform partial-fraction theorem

Write p_j=j(j+s), s=1 or 9, and let K be the retained ordinary truncation degree. Use c=1-exp(-L), h=1/[N(N+s)], and q=b(b+1). For the finite rational prefix put

    U_k(c,h)=prod_{l=0}^{k-1} (c-h p_l)/(1-h p_l),
    S(c,h)=sum_{k=0}^K (A_k+h T_{k,b}) U_k(c,h),
    T_{k,b}=((k+1)^2 A_(k+1)+[-k^2+(8-4b)k+9-2b(b-1)]A_k)/2
             +q A_k+sigma_b G_k,
    sigma_0=1, sigma_4=-1.

The physical normalized prefix is (1-c)^3 S/(1+qh)^4. For each fixed c in (0,1), S is rational in h with distinct simple possible poles h=1/p_j, 1<=j<K, and affine polynomial part. Its residues in the convention 1/(1-p_j h) are

    rho_j(c)=sum_{k=j+1}^K (A_k+T_{k,b}/p_j) rho_(j,k)(c),
    rho_(j,k)(c)=[(1-p_j h)U_k(c,h)]_(h=1/p_j).

The following exact bound is the new ingredient:

    |rho_(j,k)(c)| <= (2/j) c^(-8) exp(-(1-c) p_j/(k+s)).

Proof: let u>=0 solve u(u+s)=c p_j. Then cj<=u<j. From p_u-p_l=(u-l)(u+l+s),

    |rho_(j,k)| = ((2j+s)/p_j) |(-u)_k (u+s)_k| /
                    [j! (k-1-j)! (j+s)_k].

For integer u the residue is zero. Otherwise take k to infinity using the gamma ratio asymptotic; the limit of the magnitude is

    ((2j+s)/p_j) Gamma(j+s) / [Gamma(j+1) |Gamma(-u) Gamma(u+s)|].

Euler reflection and integrality of s bound this by (2/j)c^(-(s-1)): the gamma ratio reduces to a product of (j+l)/(u+l)<=1/c, 1<=l<s, times |sin(pi u)|/pi<=1. For l>=k>j each further factor in the residue magnitude is 1+(1-c)p_j/(p_l-p_j). The log of their infinite product is at least

    sum_{l>=k} (1-c)p_j/p_l >= (1-c)p_j/(k+s),

using log(1+x)>=x/(1+x), p_l-cp_j<=p_l, and comparison of sum 1/(l+s)^2 to its integral. Dividing the limiting magnitude by this product proves the bound. This proof applies to every fixed c and to the complete prefix, not merely a physical sample.

For a retained rectangle c in [l,u], 0<=h<=H and K fixed, the corrected coefficient file has |A_k|+|T_(k,b)|/p_j<=M=10^700, checked with integer endpoint bounds. Consequently

    |rho_j(c)| <= 2 K M l^(-8)/j exp(-(1-u)p_j/(K+9)).

For each j set v=p_j H<1, g=sqrt(1-v), r=(1-g)/(1+g), t=2h/H-1. The exact Poisson-kernel identity gives

    1/(1-p_j h) = (1/g) [1 + 2 sum_{m>=1} r^m T_m(t)].

Here T_m(t) denotes a Chebyshev polynomial, not the angular dressing coefficients. Truncation at degree m gives an absolute error <=2 r^(m+1)/(g(1-r)). Multiply by the residue bound, sum j, and retain the affine polynomial part exactly. The result Q_m(c,h) is polynomial in c and h, degree <=K in c and <=m in h (because all residue and polynomial-part coefficients are polynomials in c). Uniformly on the rectangle,

    |S(c,h)-Q_m(c,h)| <= sum_{j=1}^{K-1} [2 K M l^(-8)/j]
                  exp(-(1-u)p_j/(K+9)) * 2 r_j^(m+1)/(g_j(1-r_j)).

`residue_tail.py` checks this bound by grouping 32 consecutive j's, outward rounding g and r, bounding the exponential by an exact negative power of two using an upper enclosure of log(2), and upward rounding each nonnegative group contribution to units of 1e-25. It checks the coefficient bound against every original integer enclosure and the pole gap and cap for every rectangle. All 21 rectangles pass with m=600; max error <=531/10^25=5.31e-23. This is a uniform error for the unnormalized rational prefix. The positive physical normalization only decreases the error.

## Physical-domain norming consequence (proved reduction, not executed sign scan)

The twist cap is h<=(1-sqrt(c))^2/998^2, and the cap is convex decreasing since its second derivative is positive. On each c interval it is therefore bounded by the affine secant through outward-rounded cap endpoints. One may refine the intervals and use a secant H(c) instead of the whole enclosing constant H (require H(c)<=H). Substitute h=H(c)y, y in [0,1]. Q_600(c,H(c)y) has degrees at most K+600 in c and 600 in y. Multiplying by (1-c)^3 raises the first degree by at most three. This allows tensor polynomial norming on the coupled cap, avoiding the known unphysical negative point in the unrestricted wide box.

For example, take exact cosine grids with N_c>=8(K+603) and N_y>=4800 subdivisions in the two normalized variables, and a constant reference R for (1-c)^3 S. By the Bernstein inequality for trigonometric polynomials, every univariate polynomial p of degree d is bounded by 4/3 times its largest absolute grid value: nearest grid angular distance <=pi/(2N)<1/(4d), and the derivative of p(cos theta) is bounded by d times its supremum. Iterating over both coordinates gives 16/9. If all grid values of (1-c)^3 S-R are bounded by delta, then the continuum lower bound for (1-c)^3 S is

    R - (16/9)(delta+epsilon) - epsilon,

where epsilon=5.31e-23 or the cell-specific smaller bound. Rational enclosures for the cosine nodes and scalar values must be included in delta. Divide a positive lower bound by (1+qH)^4 and subtract the retained physical ordinary-tail and coherent spectral allowances. This is a sufficient criterion, not a necessary one; refine or choose a better reference if it fails. No scan satisfying this criterion was completed. This compression reduces the twist degree, but it does not by itself establish the c-dependent sign or replace the need for signed interval certificates. In particular the numerical scouting and the failed unrestricted rectangular extension are not evidence of continuum positivity.

## Reproduction

Requires Python 3 and GMP C++ (`g++ -O3 direct_sign_triangle.cpp -lgmpxx -lgmp -o direct_sign_triangle`). Use the preceding archive for the corrected coefficient vector and upstream proofs. Execute:

```
python residue_tail.py coeff_long_17200.txt full_grid_results.json
python verify_sign_ledger.py
for stem in next_67_250_7_25_{0,4} next_7_25_29_100_{0,4} next_29_100_3_10_{0,4}; do
  ./direct_sign_triangle < "$stem.in" > "$stem.recheck"
  cmp "$stem.out" "$stem.recheck" || exit 1
done
```

The exact arithmetic checker was executed for all six sign cases and the output ledger was independently verified. The degree-600 error ledger was executed and verified for all 21 rectangles. No ongoing computation remains. The full positivity gap remains open.
