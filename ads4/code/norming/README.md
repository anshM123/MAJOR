# Signed connection transfer: scalar norming, base-multiplier lemma and positivity-cell assembly

## Scope / result

This package implements sequential one-dimensional Chebyshev norming in `(d,h,nu)` using exact rational samples and a C++ GMP integer recurrence. It includes a new elementary spin-uniform bound for the two base multipliers, including the turning point. It includes a contour top-line integration routine and an end-to-end **conditional** positivity-cell assembly on `3.55 <= L <= 3.65`, both `b=0,4`. The new residual calculation is continuous in all state parameters and frequencies; it is not a sampled sign check.

The complete sliver is NOT proved. The assembly imports the already reported contour3 results on this particular L-cell, the LRINF lower floor, EPSW masses, cutoff and Gaussian envelope, Loc-Joint remainder, angular-action identity, and high-twist coverage. The contour3 bounds apply directly to `998 <= a <= 4500`; extending to all `a>=998` uses the prior high-twist coverage (ANGULAR_ROUTE §73). The new base estimates do not require Lemma F or Lemma X. The rest of the L-band needs full continuum coverage, contour bottom/left/right and top-line calculations, and review of retained inputs. The `L1.6` transfer result alone is not a positivity result.

## Reproduction

Requires Python 3 (standard library), C++17 compiler, and GMP development headers/libraries (`libgmp-dev` on Debian/Ubuntu). No network or Python third-party packages required. In this folder:

```bash
g++ -O3 -std=c++17 scalar.cpp -lgmpxx -lgmp -o scalar
python3 check_norm.py cross_sections.json --jobs 4 --output cross_rerun.json
python3 check_norm.py tensor_run1.json --jobs 4 --output tensor1_rerun.json
python3 check_base_identity.py
python3 check_identity.py
python3 check_lowfreq.py
python3 check_norm.py closure5.json --jobs 4 --output closure_rerun.json
python3 assemble.py closure_rerun.json
python3 check_norm.py complex_pilot.json --jobs 4 --output complex_rerun.json
python3 assemble_top.py complex_rerun.json --dlo .040355 --dhi .041169 --lo 895 --hi 905
```

`JOBS=4 ./reproduce.sh` runs these; `LONG=1` adds the expensive L=1.6 test. `python3 verify_inputs.py /path/to/MAJOR` compares inherited input ceilings to the exact JSON source data in the checkout (Decimal/Fraction). `check_norm.py INPUT.json --preflight` gives node counts, degrees, pole exclusion and real-index tail without running scalar samples. `--max-nodes` prevents unintended huge runs. `--output` records complete result objects with exact rational bounds. An `ENCLOSED` row proves only its transfer bound. `PASS` if optional `maxTransfer` requested is an engineering threshold, not a positivity test. `SKIP_NODE_LIMIT` is not coverage. All input endpoints are decimal strings. Timing is measured wall time using 4 concurrent worker processes on the test Xeon Platinum 8573C; worker_seconds sums overlapping worker batch wall times, not CPU seconds.

The historical tensor1 and tensor2 results were computed with the full J0=1500 ellipse majorant. The supplied input files have `ellipseFull:true` to reproduce that setting; the historical log's configuration predates that explicit option. Otherwise, the ellipse majorant only needs to enclose the partial sum through K-1. The real-index tail still uses all 1500 terms. The candidate grid file `tensor_boxes.json` is only a preflight proposal.

### Measured cases

* `cross_sections.json`: d=.44932896411722156 (L approximately .8), h in [0,6.68e-8]. Dressed b=4; at nu=900, 438 nodes, P<7.745,Q<2.509 (2.47s); nu=1031, 635 nodes, P<6131.170,Q<2365.166 (4.19s). Both cover all physical twists. These reproduce earlier cross-section estimates by a slightly different grid.
* `tensor_run1.json`: L approximately 3.2, relative d-width 2%, nu in [895,905], all twists. All four channels; 7296 nodes each, 1.3–2.6s per channel. P<17.48,Q<6.80. This is a continuous three-dimensional box.
* `tensor_run2_wide.json`: L approximately 1.6, relative d-width .2%, nu in [899.5,900.5], all twists; dressed b=4, 209760 nodes, 428.86s wall, P<21.434,Q<7.882. This is a genuine full-twist 3D enclosure but expensive; do not extrapolate the 3.2 runtime to the lower part of the sliver.
* `complex_pilot.json`: same small box as tensor1, horizontal complex frequency with Y=6.25 (p=3) and Y=7.25 (p=4). All four channels, 7296 nodes each, 2.1–3.2s each. `top_pilot.log` integrates only nu_real in [895,905]. E2<1.410e-14, E3<1.241e-18 for this segment and both channels. It is NOT a full top-line integral.
* `closure5.json`: 5 real frequency bands, each both b channels, d in [.02598,.02873], full twist. 137082 scalar evaluations total, K=120. Certified bounds in closure5_results.json. `assembly.log`: total error/floor <=.102850 (b=0), <=.108191 (b=4), floor 1.98e-6. Conditional normalized action lower bounds >1.7763e-6 and >1.7657e-6. Raw exact rational results and per-band budget in assembly_results.json. These new cells overlap the previously claimed L>=3.65 region.

## 1. Transfer recurrence and normalization

Let D=N(N+s), h=1/D, c=n(n+s)/D, d=1-c, s=2b+1, q=b(b+1). The physical cap is

    h <= d^2 / [a^2 (1+sqrt(1-d))^2].

Use the exact full-box cap computed from d_max with a_min=998 and the square root rounded down. Let beta_j=alpha-j, alpha=p-1/2+i nu, X_j=d^(alpha-j) S(alpha-j), S(beta)=L((1-z)^(-beta)), c_j=(1/2-i nu)_j^2/((1-2i nu)_j j!), r_j=c_{j+1}/c_j. Then Z_j=c_j*d^j*X_j, r_-1=1, Z_-1=X_-1/d. The exact weighted recurrence is

    (1-h beta_j(beta_j-s)) Z_(j+1)
      = r_j [d+h beta_j(s-2 beta_j-1)] Z_j
        + r_j r_(j-1) h beta_j(beta_j+1) Z_(j-1).

Two independent basis solutions start with (Z_0,Z_-1)=(1,0),(0,1). The scalar checker sums the basis coefficients. For p=4 the sum is P X_0+Q X_-1. Q_unscaled is divided by d to obtain Q. For p=3 use A0=9/2-b(b-1)-beta(beta+8-4b)/2, A1=beta(beta+9-4b)/2. Sum the terms `(1+qh+h A0) Z_j + h A1 r_(j-1) Z_(j-1)`. The physical dressed result is `(P X0+Q X_-1)/(1+qh)`; for real h>=0 the denominator is at least 1. The source clears denominators exactly, using rational d,h and complex rational frequency. Each division of the Gaussian integer numerator has real and imaginary truncation error at most 1 scaled unit each. For coefficient combination `(A u+B v)/D`, complex l1 error is bounded by `( |A|_1 e_u + |B|_1 e_v ) |D|_1/|D|^2 + 2` in scaled units, rounded upward. `check_identity.py` checks 576 exact weighted recurrence/dressing identities and 48 transfer sums. No floats are used for certificate comparisons.

For p=3 on the real frequency axis the normalized truncated dressed kernel obeys

    d^3 |K_J0(nu)| <= 2 sqrt(d) |B(nu)| (|P| M0 + |Q| M1),

where |B|^2=coth(pi nu)/(4 pi nu). The factor 2 accounts for both outgoing/incoming branches. For p=4 replace sqrt(d) by 1/sqrt(d) and include physical tau*Lambda^2 where needed. On a complex contour use only the outgoing branch, not conjugation. `assemble.py` uses these normalizations.

## 2. Anisotropic sequential norming proof

Let f be one signed transfer coefficient of the partial sum j<K. For one normalized real variable x in [-1,1], freeze the other variables anywhere in their real intervals. On the Bernstein ellipse of parameter rho>1 in x, bound |f|<=M by the directed majorant in `major_axis`. This majorant is uniform in the other real variables. Its proof follows the exact weighted recurrence, bounding each coefficient on the affine complex disk of radius (rho+rho^-1)/2 in x and radius 1 in the other normalized axes. Centered complex products are bounded by the triangle inequality; `D_j` and the denominators of r_j are checked away from zero. Both independent basis solutions have initial magnitudes <=1. The positive recurrence bounds the sum of their absolute magnitudes, including dressing. The majorant for j>=K on the full real box bounds the discarded connection indices K..1499. Only the finite partial sum enters ellipse approximation; no analytic continuation of the asymptotic infinite j-series is assumed.

The Chebyshev truncation at degree m satisfies ||f-p_m||<=epsilon=2 M rho^(-m)/(rho-1). This follows from Cauchy's estimate for the Laurent coefficients after x=(z+1/z)/2, summed geometrically. Choose m by an exact integer inequality. At real Lobatto nodes x_j=cos(j pi/N), a trigonometric polynomial T(theta)=p_m(cos theta) of degree m obeys the sharpened Bernstein inequality for its real part after a constant phase rotation. At a maximizer of |T|, align its phase so the real trigonometric polynomial reaches its maximum M_p. The inequality (T')^2<=m^2(M_p^2-T^2) for the real part yields T(theta)>=M_p cos(m |theta-theta*|) in the relevant neighborhood (integrate arccos(T/M_p)). Hence if every angle is within pi/(3m) of a sample, ||p_m||<=2 max_j |p_m(x_j)|. Consequently ||f||<=2 max_j |f(x_j)|+3epsilon.

The checker uses exact rational approximations to cosine nodes, denominator 10^12. Machin arctan series and the cosine Taylor recurrence at scale 10^50 give node-coordinate error <10^-12. For angles in [0,pi], |cos a-cos b|>=2 sin^2((a-b)/2)>=2(a-b)^2/pi^2. Thus their angular errors are <3e-6. `proof_node_factor` checks m(pi/(2N)+3e-6)<pi/3 using rational pi bounds. No estimate of f at the rounded node via an enormous complex norm is necessary: the norming theorem applies at the actual rational nodes. Endpoints exact. The Machin series uses <100 nonzero terms each; its omitted alternating term is <= one unit and termwise truncations cost <2000 units in pi. Taylor recurrence propagation, including the squared argument rounding, is bounded by 100 times a few thousand units times exp(pi)<24, far below 10^-42; the omitted term after a computed zero is also covered by this bound. Rounding to 12 decimals adds at most 0.5e-12.

Apply the 1D inequality sequentially in d,h,nu. This produces a tensor grid but **does not require a simultaneously large complex polydisk**. Each dimension has its own degree and radius and is extended by f_bound <- 2*f_bound+3epsilon; fixed dimensions need no extension. Q is normed in its unscaled form: sample |Q| is multiplied by d_max to bound Q_unscaled, extend that function, add the real discarded tail, then divide by d_min. Finally add the real tail to P. All square roots, absolute values, epsilon calculations and pole gaps are directed integer operations. Reported decimal approximations are for inspection; exact rational strings are authoritative.

## 3. New state-local eigenframe bound for base multipliers, including turning

For s=1, fixed real N and integer n, let eta(x)=(2x+1)/[(N-x)(N+x+1)] = 1/(N-x)-1/(N+x+1). The exact dual-Hahn recurrence, valid also for noninteger N, is

    S_(k+1) = S_(k-1) + 2 alpha eta(k) S_k, S_0=S_-1=1.

The terminating ladder identity follows by comparing coefficients of (alpha)_j/j!: U_(k+1,j)-U_(k-1,j) = 2 eta(k) j [U_(k,j-1)-U_(k,j)] (with terminated coefficients zero), using alpha*(alpha)_j=(alpha)_(j+1)-j(alpha)_j. `check_base_identity.py` checks 132 exact spin steps and 48 helicity identities, including noninteger N.

Take alpha=sigma+i nu, sigma>0, nu>0; set kappa=sigma/nu, y(x)=nu eta(x), w=(kappa+i)y, r=sqrt(1+w^2) with Re r>0, lambda_+=w+r, lambda_-=w-r=-1/lambda_+. Then |lambda_+|=exp(u)>=1, u=Re asinh(w). Let T(y)=[[lambda_+,lambda_-],[1,1]]. The spin-step matrix is T diag(lambda_+,lambda_-) T^-1. Rather than use endpoint differences or an Airy window, control coordinate changes by the exact continuous matrix ODE. Submultiplicativity and Gronwall give

    ||T(y2)^-1 T(y1)||_inf <= exp(int_y1^y2 ||T^-1 T'||_inf dy).

The derivatives lambda'_±=±(kappa+i) lambda_±/r imply the key identity

    ||T^-1 T'||_inf = sqrt(1+kappa^2) cosh(u)/|1+w^2|
                    = sqrt(1+kappa^2)/2 * (1/|w+i| + 1/|w-i|).

Here cosh(u)=(|w+i|+|w-i|)/2 by the ellipse identity for w=sinh(u+iv). Thus the integral to an endpoint Y>=y(n) is exactly

    A(Y,kappa) = [asinh(((1+kappa^2)Y+1)/kappa)
                 +asinh(((1+kappa^2)Y-1)/kappa)]/2.

This is finite across y=1 for every sigma>0. The initial norm is <=1.01 when |alpha|/[N(N+1)]<.002: |r0-1|<=|w0|^2 and ||T_0^-1(1,1)||_inf <= (1+max|lambda_±|)/(2|r0|). The final row norm is |lambda_+|+|lambda_-|=|w+i|+|w-i|. The eigenvalue growth is normalized by d^sigma exactly using continuous spin: u(y) is nondecreasing on y>=0 (sum of distances is even convex); int_0^n eta(x) dx=-log d; eta'(x)>=eta(x)^2. Therefore

    sigma log d + sum_(k<n) u(y(k))
       <= int_0^n (u(y(x))-kappa*y(x)) dx
       <= nu int_0^Y max(u(y)-kappa*y,0)/y^2 dy.

For y<1, integrating horizontally in w yields u<=kappa*y/sqrt(1-y^2); for y>1, u<=acosh(y)+kappa*y/sqrt(y^2-1), since |1+(x+iy)^2|>=|1-y^2|. For Y>1 use acosh(y)<=sqrt(2(y-1)), drop the negative normalizing term and bound 1/y<=1. The resulting elementary growth bound is

    B(Y)=sigma log(2/(1+sqrt(1-min(1,Y)^2)))
        + 1_(Y>1) [ (2 nu/3)*(Y-1)*sqrt(2(Y-1))
                         + sigma*sqrt(2(Y-1)) ].

So the state-local bound is

    |d^alpha S_n(alpha)| <= 1.01 (|w(Y)+i|+|w(Y)-i|) exp(A+B).

To bound the actual turning point, put t'=a/(N+1/2). Exactly t'(2-t')=d(1-1/(2N+1)^2), y(n)=(nu/a)*2(1-t')/(2-t'). Take d_eff>=d_min*(1-1/(4 a_min^2)), R=sqrt_up(1-d_eff) and Y=(nu_max/a_min)*2R/(1+R). This is materially sharper at low L than the unlocalized bound nu/a. For a frequency/exponent interval, use kappa_min=sigma_min/nu_max, kappa_max=sigma_max/nu_min and qmax=sqrt(1+kappa_max^2). The integrand for A is bounded using denominator kappa_min and numerator qmax, so A <= (qmax/sqrt(1+kappa_min^2))*A(Y,kappa_min); the code safely drops the divisor. Use sigma_max,nu_max in B and kappa_max in the final row norm. Positive-power starting exponents have sigma>=1/4. For small nu, use the rising-factorial comparison `|(sigma+i nu)_k| <= (sigma)_k exp(nu^2*(sigma^-2+sigma^-1)/2)` and U_k<=c^k. Negative/small real exponents are reached by finitely many exact downward beta recurrences from exponents >=1/4. For s=1 the denominator modulus is >=1-h_max*(|Re beta|+1)^2 and the numerator absolute terms are bounded by `1+2h_max|beta|^2/d_min`, `h_max |beta||beta+1|/d_min^2`.

For b=4 use the exact helicity identity (n'=n+4,N'=N+4): U_k^(9)=U_(k+4)^(1)/U_4^(1), hence L^(9)(F)=L^(1)(z^4 F)/U_4. Expand z^4=(1-(1-z))^4. The shifted d0=d/(1+20h) and U_4>= (1-d_max-12 h_max)^4 (each ladder denominator <=1 and numerator >=1-d_max-12h_max). This gives the weighted sum of five s=1 bounds used in base_bounds.py, with factor (d/d0)^Re beta. All transcendental constants in base_bounds.py are bounded via integer/Fraction square roots, atanh log series, and exp series; no special-function sampling is used for the certificate. The helicity bound can be loose for small L; it is an upper bound, not an optimized estimate.

## 4. Conditional complete positivity assembly on [3.55,3.65]

The d-box [.02598,.02873] contains exp(-L) throughout the target cell (checked directed). `closure5` uses five overlapping real-frequency boxes for nu>=159 through 1040.1, both helicities, all twists. For each of 16 sub-bands in band_masses.csv, combine the enclosing transfer result with the sharper state-local M0,M1 evaluated on that sub-band. The retained mass ceilings come from `ads4/code/analytic2/cache/rmass3_40_898.json` and `ads4/code/review_l5/cache/rt_rmass2.json`, which contain earlier reproduced EPSW upper bounds. Each of the 881 source rows (integer starts 159 through 1039) is included once with a factor 1.001 to cover stored decimal rounding. Their full [T0,T1] intervals are inside the corresponding sub-band. The 0–160 result is inherited from the previous checker, rerun here by check_lowfreq.py, bounded by 1.46e-7. The overlaps at 159/160 and 1040/1040.1 are harmless double charges. The `2/Lambda` prefactor is bounded by 2/998. The cutoff L1 error <=1.34e-14 is charged again against the largest finite-band normalized kernel bound (thus no false assertion |f_theta|<=|f|). 

Above 1040 use the inherited Gaussian pointwise envelope `|f_theta(nu)| <= exp(21-(nu-1029)^2/(2*.901^2)) + exp(-nu^2/3)`. On [1040,1050],[1050,1070],[1070,1100],[1100,1150] use the entire positive recurrence majorant (no samples needed), local M0,M1 and Gaussian mass bounds from the tangent-line exponential integral. Above 1150 use the unsplit full conical kernel with bound `d^3 |V F3| <= 200 exp(pi nu/2)` by the Beta-Cauchy representation, K0<1.6, Cauchy radius t/2, d<=2t, and Mehler bound. Integrate the Gaussian tail with w=11/7>pi/2; no underflow is used (exact powers of two bound exponentials of large magnitude). The Loc-Joint integrated connection remainder is retained as <exp(-46)<1.1e-20. The contour pieces are the conservative rounded-up entries of `ads4/code/analytic3/cache/contour3_3.55_3.65_k4.json`, not recomputed here. That source gives F=1.9866191481e-6; take F_lower=1.98e-6. Then assembly checks total error < floor by exact rational comparison. The contour bound applies to a<=4500; use §73's high-twist certificate for a>=4500.

## 5. Complex top line support

`assemble_top.py` accepts certified norming results with horizontal complex frequency and verifies complete segment coverage for each available channel. It subdivides integration into pieces of length <=5, bounds the exact NF02 polynomial on each complex rectangle using a centered Taylor coefficient norm, and uses the base bound for negative real exponents. For nu=x+iY, B(nu)=2^(2i nu-1)/sqrt(pi) * Gamma(i nu)/Gamma(i nu+.5). Shift z=i nu by integer k>Y, bound the finite product |z+j+.5|/|z+j| at the x interval endpoints, and use the positive beta integral with real part sigma=k-Y>0: |Gamma(z+k)/Gamma(z+k+.5)|<=Gamma(sigma)/Gamma(sigma+.5)<=sqrt(sigma+1)/sigma (Gautschi). For the infinite-twist outgoing term, Euler's integral gives

    |B(nu) H(d)| <= coth(pi x)/(2 pi) * int_0^1 [u(1-u)]^(Y-.5) (1-du)^(-Y-.5) du
                  <= 4^(-Y+.5) (1-d)^(-Y-.5)/pi  (x>=1/4,Y>=1/2).

On the top line the finite and infinite terms may be bounded in modulus separately. The physical d^3 factor is d^(3-Re alpha), the split contour prefactors are 4/Lambda for p=3 and 2 tau Lambda for p=4, with tau<=hcap. The complex pilot has Y=6.25,7.25 respectively and integrates only x in [895,905]. `full_top:false` explicitly records the missing remainder of the contour. A complete top-line cover can be assembled by passing a norming results JSON tiled across x=1/4..999 (the latter exceeds Lambda). This does not certify bottom/left/right; those require comparison with the infinite-twist term or the retained contour method. The [3.55,3.65] complete assembly uses the repository's earlier entire contour bound instead.
