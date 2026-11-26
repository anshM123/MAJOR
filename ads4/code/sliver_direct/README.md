# Direct realized-weight sliver certificate

**Scope:** conditional on the retained angular-action identity, physical state parametrization and spectral-envelope certificates in MAJOR. This is not a proof of the entire sliver. The new direct sign certificates cover both channels on `0.05 <= c <= 0.09`, where `c=1-exp(-L)`. The prior upper coverage begins at L=2.3. The intervening band is still open. The attempted analytic slab [2.15,2.3] was interrupted with seven of ten transfer boxes enclosed; no sign assertion is made for it.

## Reproduction

Requirements: Python 3 with SymPy; C++17, GMP/GMPXX, MPFR. The code uses the public MPFR 4.x ABI declarations when its headers are unavailable. On Linux with only the MPFR runtime use `-l:libmpfr.so.6` in place of `-lmpfr`.

```
g++ -O3 -std=c++17 gaussian_fast_adjacent.cpp -lgmpxx -lgmp -lmpfr -o gaussian_fast_adjacent
g++ -O3 -std=c++17 seed_integrals.cpp -lgmpxx -lgmp -lmpfr -o seed_integrals
g++ -O3 -std=c++17 log_constants.cpp -lgmpxx -lgmp -lmpfr -o log_constants
g++ -O3 -std=c++17 raw_moment_transform.cpp -lgmpxx -lgmp -o raw_moment_transform
g++ -O3 -std=c++17 moment_transform.cpp -lgmpxx -lgmp -o moment_transform
g++ -O3 -std=c++17 direct_sign_cheb.cpp -lgmpxx -lgmp -o direct_sign_cheb
python inventory_cut.py weights.json
python pole_bound.py
python raw_seed_budget_actual.py
python pole_bound_full.py
python raw_seed_budget_full.py
python make_quad.py
./seed_integrals < seeds_input.txt > seed_values.txt
./log_constants > log_constants.txt
python run_full_gauss.py 1254 1600 4
python make_coefficients.py
python full_prefix.py
for b in 0 4; do
  python direct_sign_cheb.py .05 .065 $b
  python direct_sign_cheb.py .065 .09 $b
done
python verify_direct_results.py
```

The archive includes the full moment table, coefficient vector and seed values, so the last four sign computations can be rerun without recomputing the spectral integration. The full inventory calculation took 475.25 seconds wall time with four workers, while sharing CPU with the incomplete analytic slab calculation. The two channels in the first direct cell took about 46 seconds each and those in the second about 72 seconds each, using separate workers. Times are not guarantees for other hardware.

## 1. Exact inventory and global moments

`weights.json` is the transfer of the repaired Gaussian center/width/amplitude triples and the three logarithmic coefficients from repository blob `d7e513e1b6aacbaaab6ac3adde098163b86558d0` at `ads4/prl/data/k2/k2_continuum_D_1000_rep.json`. The seven included base64 chunks contain exact little-endian binary64 amplitudes. All 3564 triples were transferred in original order; the FNV1a checksum of their concatenated little-endian binary64 representation is 1930466845, checked at source and destination. In the calculation each value is converted into its exact binary rational. The source inventory is retained input, not a new proof of the physical spectral construction.

Use the raw unshifted weight in `fullpos/fp_moments.py`:

```
g = -(R(nu^2) B + Rt(nu^2) L) theta / Lambda
R = -96/(4nu^2+1) + (5184/5)/(4nu^2+9) - (4704/5)/(4nu^2+49)
Rt = -(9216/5)/(4nu^2+9) + (55296/5)/(4nu^2+49)
L = (1+nu^2) (y0+y1 nu^2+y2 nu^4) log(1+nu^2) exp(-nu^2/2)
```

Here V=1000, Lambda=sqrt(997002), theta(nu)=[erf((nu+nu_c)/10)-erf((nu-nu_c)/10)]/2, nu_c=Lambda+100. For a Gaussian center c and width s put D=100+2s^2. Full-line normalized moments M_k=int G theta (nu/V)^k satisfy

```
M_0=sqrt(2 pi)s[1-(erfc((nu_c-c)/sqrt(D))+erfc((nu_c+c)/sqrt(D)))/2]
M_(k+1)=(c/V)M_k+k(s/V)^2 M_(k-1)+(s^2/V)(E^-_k-E^+_k).
```

The ordinary Gaussian sequences E^± have means `(100c±2s²nu_c)/(DV)`, variance `100s²/(D V²)` and mass `sqrt(2s²/D) exp(-(c∓nu_c)^2/D)`. This is integration by parts and completing the square. For each center c>0 the full-line even moments account for the two symmetric half-line contributions; c=0 is halved. Signed amplitudes are included before summation. Every operation is outward-directed MPFR arithmetic. `run_full_gauss.py` sums shard endpoints exactly. The complete output `S_1254_1600.txt` has 1255 moment intervals at scale 10^1600.

### Adjacent-center erfc acceleration

For x_new=x_old-delta with delta>=0, the exact increment is

```
erfc(x_new)-erfc(x_old)=(2/sqrt(pi))*delta*exp(-x_old²)*int_0^1 exp(2*x_old*delta*t-delta²*t²) dt.
```

For the Taylor coefficients b_0=1, b_-1=0, `(k+1)b_(k+1)=2*x_old*delta*b_k-2*delta²*b_(k-1)`. The Cauchy bound for truncating after degree m on [0,1] is the prefactor times `exp(2*x_old*delta*R+delta²*R²) R^-m/(R-1)`. The code selects a power-of-two radius and degree by directed comparison to 10^(-digits-12). It uses the exact integral of the polynomial and a symmetric interval remainder. For remote nonnegative arguments use Mills' bound `erfc(x)<=exp(-x²)`. This avoids unnecessary high-precision special-function evaluations. The source's Gaussian centers are all nonnegative and below nu_c, so the argument is positive. The far-edge erfc shortcut retains its error interval.

## 2. Correlated seed enclosures

Set `F_(j,a)=int_0^infty B theta (nu/V)^(2j)/(4nu²+a)`. Then `F_(j+1,a)=(S_j-aF_(j,a))/(4V²)`. An error in F_0 is a single atom at nu²=-a/4. To integrate the uncut seeds, the program uses a full-line Gaussian centered at c, normalized variable t, and halves c=0. For c>=10 it expands `1/(4(c+s*t)^2+a)` in t. Partial fractions at the two conjugate poles give the exact geometric Taylor remainder after degree 2k+1. The absolute normalized Gaussian mean of that remainder is at most

```
(1/a) (s²/(c²+a/4))^(k+1) (2k+1)!!.
```

Indeed the final real-variable denominators have modulus >=sqrt(a)/2. Odd powers have zero Gaussian mean. The code stops only when the remainder is below 1e-52, then includes the Gaussian mass `sqrt(2 pi)s`. For c<10, s<=.451, it integrates t in [-16,16] with width-1/4 panels, using the 97 rational cosine nodes and exact interpolatory quadrature weights in `quad.txt`. The weights integrate every polynomial through degree 96 exactly, because they are the integrals of the cardinal Lagrange polynomials; `make_quad.py` constructs these directly with integer arithmetic. Their l1 norm is <3. Taylor truncation at degree 96 on the disk of radius 4 in panel coordinates has error <= M*4^-97/(1-1/4) on the panel. On the complex disk the imaginary part of t is <=.5, so the Gaussian modulus is <2 and the rational denominator modulus is >.79 for a>=1. Thus M<3. The total weighted quadrature discrepancy plus the Gaussian tails at |t|>16 is <1e-53. The code pads by 1e-53 before multiplying by s.

The same quadrature integrates the three log seeds on nu in [0,24]. On a panel disk of radius 4, |Im nu|<=.5, Re(1+nu²)>=.75, and |nu|<25. The principal log is analytic there, and |log(1+nu²)|<9. For the rational log integrands, |4nu²+a|>=a-1 (a=9,49), |y_i|<1, and |exp(-nu²/2)|<2. A common bound M=2e9 is safe. Taylor quadrature error is <= (24/2)*(2+3)*M*4^-97/(1-1/4), and the Gaussian tails at nu>24 with the polynomial envelope are far smaller than 1e-45. The code pads by 1e-45. All six resulting seed widths are below 1e-40.

For the log moments use

```
G_(j+1)=(2j+1)G_j/V²
E_(j+1)=(G_j-E_j)/V²
K_(j+1)=(2j+1)K_j/V²+2E_(j+1)
```

with high-precision G0=sqrt(pi/2), E0=(pi/2)exp(1/2)erfc(1/sqrt(2)), and K0=int log(1+nu²) exp(-nu²/2). Form `L_j=sum_(i=0)^2 y_i(V^(2i+2)K_(j+i+1)+V^(2i)K_(j+i))` and apply the rational-pole recurrence again. The K0 error is a Gaussian moment sequence; the other log seed errors are correlated pole atoms. The chosen seed midpoints are treated as exact in the transport, not as independent uncertain moment intervals. `inventory_cut.py` verifies the full Gaussian cutoff discrepancy (<7.1e-24 in L1 before the denominators). `pole_bound.py`, `raw_seed_budget_actual.py` and `verify_direct_results.py` budget the coherent seed substitutions at action level. For c<=.56, the total Gaussian pole correction including numerical seed errors is <2.4e-20, the log seed correction is <7.35e-26, and dropping theta on the log component costs <1e-400. For c<=.9 the pole correction is <1.8e-15 and log seed correction <1e-20. The sign tests here use the conservative combined action allowance 1e-18, valid only in the low-c domain <=.56, including the prefix truncation error.

`raw_moment_transform.cpp` transports the raw moments by the exact triangular CDH recurrence

```
I_(0,j)=0, I_(1,j)=T_j
4(m+1)^2 I_(m+1,j)=(8m²-1)I_(m,j)-4V² I_(m,j+1)-4(m-1)^2 I_(m-1,j).
```

Then `Z_k=12 I_(k+2,0)-(k+1)I_(k+1,0)+(k+3)I_(k+3,0)` and `A_k=binom(k+3,3)^2 (Z_k-k Z_(k-1)/(k+3))/(k+3)`. The retained raw angular identity gives the dressed coefficients. The ideal G_D is transported separately from exact NF02 psi3 moments. `direct_coefficients.txt` contains 1253 pairs of midpoint/radius enclosures at scale 10^1600. Their rounding-enclosure widths are <1e-375 **for the seed-substituted coefficient vector**; the coherent discrepancy from the true repaired weight is bounded separately at action level, not to 375 digits per coefficient.

## 3. Prefix and direct positivity proof

The retained physical parametrization is `s=2b+1, q=b(b+1), n=J-b, c=n(n+s)h, h=1/[(n+a)(n+a+s)]`, `a>=998, J>=1001`. Hence `h<=c/[(1001-b)(1002+b)]` and `h<=(1-sqrt(c))²/998²`. The exact rational ladder is

```
U_k = product_(j<k) (c-h*j*(j+s))/(1-h*j*(j+s))
```

with exact termination at n. On physical states 0<=U_k<=c^k through n. `full_prefix.py` uses the retained Mehler/weight envelope to prove `|H-H_prefix|<2^-110` on 22 contiguous c-cells tiling [.05,.9]. Coefficients are bounded on R=.9995 by `2^2577 R^-k` for nu<=1100; the weighted Gaussian tail at larger frequency is accounted for separately. K is chosen so `(c_hi/R)^(K+1)<=2^-2700`. Every denominator gap through K-1 is checked. Largest K=17847. This reduction by itself proves no sign.

For either direct cell, choose the rectangular physical spin cap rounded upward to 1e-12 and represent c=c0+dc*u, h=h0+dh*v with u,v in [-1,1], h0=dh. Write `F_k=(1+qh)A_k+h(B_k+sigma G_k)`, where

```
2 B_k=(k+1)^2 A_(k+1)+[-k^2+(8-4b)k+9-2b(b-1)]A_k.
```

The sign polynomial is the nested rational function `S_k=F_k+(c-p_k h)/(1-p_k h) S_(k+1)`, p_k=k(k+s), S_K=F_K. Normalization is `(1-c)^3/(1+q h)^4`. The exact Taylor rectangle of S_0 through orders mc,mh can be evaluated by triangular fixed-point arithmetic in the truncated polynomial ring. After clearing all affine parameter denominators, with common den and c0=C0/den etc., each coefficient V of the ratio times the next nested sum X is computed by

```
(den-p H0)V_ij=(C0-p H0)X_ij+DC X_(i-1,j)-p DH X_(i,j-1)+p DH V_(i,j-1).
```

The l1 error recurrence uses the lower inverse-denominator gap `den-p(H0+DH)` and adds at most one unit per coefficient division. It also includes the input coefficient radii. `direct_sign_cheb.py` bounds the analytic rational sum on a larger complex polydisk by the absolute Horner majorant, checks the denominator gaps exactly, and bounds the omitted Taylor rectangle by

```
M * (Rc^(-mc-1)+Rh^(-mh-1))/[(1-1/Rc)(1-1/Rh)].
```

The choice of radii and degrees is heuristic, but every error inequality is checked in exact rational arithmetic. Convert the **signed** Taylor polynomial first to 80-digit scale, accounting for all rounding. Then convert to tensor Chebyshev using `x^k=2^-k sum_j (2-delta_(j0)) binom(k,(k-j)/2) T_j(x)` over parity-compatible 0<=j<=k. The coefficients of this map are nonnegative and each input column sums to one, so the l1 error does not increase. The C++ code performs the map with integer binomial weights and a common dyadic denominator, counting each truncation. Finally the constant Chebyshev coefficient minus the sum of the absolute nonconstant coefficients, rounding budget and analytic tail is a rigorous lower bound throughout the rectangle. Subtract the physical coherent seed/truncation allowance after normalization.

### Executed result

Four channel-specific cells have passed. `direct_ledger.json` and the exact `.json`, `.out`, `.log`, `.in` files give the rational ledgers. `verify_direct_results.py` independently checks the exact sums of shard endpoints, full inventory byte checksum, numerical seed widths, coefficient radii, norm-one Chebyshev transformation, pole budgets, spin caps, prefix cutoff, contiguous coverage, and the four final inequalities. The low-c checker has not been independently reproduced by the repository coordinator. Earlier failed monomial coefficient bounds are retained under `monomial_fail_*`; they are not counterexamples.

| c cell | b=0 lower | b=4 lower |
|---|---:|---:|
| [.05,.065] | 2.2097043954e-5 | 2.1954815148e-5 |
| [.065,.09] | 1.9052253503e-5 | 1.8852543276e-5 |

These connect with the earlier small-L corner, subject to its retained certificates. The new upper endpoint is L=-log(.91)~.09431068. The prior high-L coverage starts at 2.3. No claim is made for the open interval between them. A larger moment table is necessary to apply the same global direct-prefix procedure all the way to c=.9. The present table is sufficient for more low-c cells, but their sign inequalities have not been run. No process was left running to deliver results later.
