# Sliver continuation: L=2.3–2.6 and correlated raw-weight moment seeds

**Not a full-sliver proof.** The exact angular-action identity, EPSW repaired spectral data and cutoff estimates, LRINF floor, low-frequency residual theorem, base-multiplier and localization bounds are inherited. The low-L direct-prefix reduction from the preceding archive gives a pole-free finite rational function, not a positivity proof. This package does not generate the full repaired moment vector. The repository data are available on the user's checkout; the complete JSON could not be mounted directly in this computation. The Gaussian sample tests are explicitly subset tests.

## Analytic slab

The 16 contiguous LRINF beta cells in [3.6,4.0] have minimum floor 3.4606212402066455e-6. `assemble_23.py` uses the conservative floor 3e-6 and verifies coverage by rational exponential bounds. The physical d interval is [0.0742735,0.1002589]. Signed transfer boxes use h=d^2 x and the physical cap x<=1/[998^2(1+sqrt(1-d))^2]. The norming proof is the same as in the preceding package, with N>=1.4m nodes per nonconstant axis and norm inflation 5/2 instead of 2; the exact test cos(1.13)>2/5 includes rational node rounding. Each axis adds 7/2 times the Chebyshev approximation error. The binary-fixed-point scalar variant rounds each exact complex ratio to l1 error <=2 units, and propagates errors through exact integer multiplications. Its intervals overlap the original decimal implementation on the 1000 archived benchmark inputs. Both implementations are directed; a run can use either. The new Gaussian tail uses the coupled coordinate and smaller bands near 1040. It replaces a failing, too-loose uncoupled tail bound. For p=4 use a normed right edge on Im nu in [0,25], then absolute bounds above 25, and top height 30.25. The p=3 height is 15.25. `ledger_23.json` records only fully assembled channels. Any missing channel or failed row means the slab is not closed.

Reproduction (GMP, Python 3):

```
g++ -O3 -std=c++17 scalar_binary.cpp -lgmpxx -lgmp -o scalar
python check_norm.py boxes_23_p3.json --jobs 4 --max-nodes 300000 --output results_23_p3.json
python check_norm.py p4_right.json --jobs 4 --max-nodes 200000 --output results_23_p4_right.json
python resid_budget.py results_23_p3.json 0.000003
for b in 0 4; do
  python ellipse_extend.py results_23_p3.json 15.25 ct_23_3_$b.json 3 $b
done
python p4_contour_23.py results_23_p4_right.json
python assemble_23.py
```

## Global moment acceleration

Use raw unshifted g, in the notation of the repository `fullpos/fp_moments.py`:

```
g = -(R(nu^2) B + Rt(nu^2) L) theta / Lambda
R = -96/(4nu^2+1) + (5184/5)/(4nu^2+9) - (4704/5)/(4nu^2+49)
Rt = -(9216/5)/(4nu^2+9) + (55296/5)/(4nu^2+49)
L = (1+nu^2) (y0+y1 nu^2+y2 nu^4) log(1+nu^2) exp(-nu^2/2)
```

Set V=1000, Lambda=sqrt(997002). The exact binary64 centers/widths/amplitudes are the data. For one Gaussian centered at c, width s, theta(nu)=(erf((nu+nu_c)/10)-erf((nu-nu_c)/10))/2 with nu_c=Lambda+100, D=100+2s^2. The full-line normalized moments M_k=int G theta (nu/V)^k are given by:

```
M_0=sqrt(2 pi)s [1-(erfc((nu_c-c)/sqrt(D))+erfc((nu_c+c)/sqrt(D)))/2]
M_(k+1)= (c/V)M_k + k s^2/V^2 M_(k-1) + (s^2/V)(E^-_k-E^+_k).
```

The E^± are ordinary Gaussian moments with normalized mean (100c ± 2s²nu_c)/(D V), variance 100s²/(D V²) and mass sqrt(2s²/D) exp(-(c∓nu_c)²/D). This follows from integration by parts and completing the square for G theta'. Symmetric twin centers give the full-line M_(2j) as the half-line moment S_j; c=0 is halved. `gaussian_cut_moments.cpp` implements the recurrence with directed MPFR intervals. For remote cutoff arguments x>=80 it uses 0<=erfc(x)<=exp(-x²), with the output interval widths checked. The original implementation and 36 quadrature diagnostics remain in the archive. `gaussian_high_test.py` checks seven actual repo centers, 3254 even moments through degree 6506 at 2500 decimal places. No full spectrum result is claimed.

```
g++ -O3 -std=c++17 gaussian_cut_moments.cpp -lgmpxx -lgmp -lmpfr -o gaussian_cut_moments_fast
python gaussian_high_test.py
```

The executable reads K, digits, N, then N rows of exact rational center, width, amplitude (each as numerator denominator); outputs j and directed integer lower/upper bounds for S_j scaled by 10^digits. For the full inventory, provide all 3564 rows from `k2_continuum_D_1000_rep.json`, `Fraction(float(...))` on each field as in the retained repository. The computation can be partitioned by centers and the integer endpoints summed exactly. The code has no dependency on python-flint for this step.

### Correlated rational seeds

Let F_(j,a)=int B theta (nu/V)^(2j)/(4nu²+a). Then F_(j+1,a)=(S_j-a F_(j,a))/(4V²). An error eps in F_(0,a) changes every F_j by eps*(-a/(4V²))^j; it is a single evaluation atom at nu²=-a/4. For a=1,9,49, the exact raw CDH recurrence produces polynomial A_k, checked in `pole_bound.py`. Both-channel normalized action norms of the atoms on c<=.56 are bounded by 84,1006,162133. The raw Gaussian seed cutoff difference is <=7.1e-24/a; the inventory test was executed via repository connector with exact binary64 sums, upper numerator 709561660317879566367702162799953560927323381840827401202 / 10^80. `inventory_cut.py` reproduces the test on the checkout JSON. This gives total normalized action error <2.4e-20. The three uncut F0 can thus be computed using the existing Faddeeva formula at moderate precision instead of 2500 digits. Their actual numerical precision error must also be charged via the same atom bounds; e.g. errors <=1e-60 each are negligible. No independent propagation of the seed uncertainty is permitted.

### Logarithmic seeds

With G_j=int (nu/V)^(2j) exp(-nu²/2), E_j=int (nu/V)^(2j) exp(-nu²/2)/(1+nu²), K_j=int (nu/V)^(2j) exp(-nu²/2) log(1+nu²), use

```
G_(j+1)=(2j+1)G_j/V²,
E_(j+1)=(G_j-E_j)/V²,
K_(j+1)=(2j+1)K_j/V²+2 E_(j+1).
```

G0=sqrt(pi/2), E0=(pi/2) exp(1/2) erfc(1/sqrt(2)) must be enclosed at high precision. K0 and the two log rational seed integrals can be computed at certified moderate precision (<=1e-60). Their correlated errors consist of a Gaussian polynomial rational weight (K0 error) and the two corresponding imaginary pole atoms. `raw_seed_budget.py` proves that the log seed errors cost <1e-38 of normalized action on c<=.56; dropping theta on the logarithmic part costs <1e-400. This uses an explicit raw real-frequency kernel bound <1e14 exp(11nu/7) and elementary Gaussian integrals. This is a budget on the physical action, not on independent coefficient errors.

Compute L_j=sum_(i=0)^2 y_i (V^(2i+2) K_(j+i+1)+V^(2i) K_(j+i)). The log rational moments have the same recurrence as F_(j,a), with S_j replaced by L_j. Form the raw moments T_j=int g (nu/V)^(2j) by the displayed partial fractions. Enclose high-precision rounding and the elementary seed G0,E0 errors in the input balls; charge the moderate seed/cutoff errors at action level. The high-precision transport is a triangular recurrence on scaled integers:

```
I_(0,j)=0, I_(1,j)=T_j;
4(m+1)^2 I_(m+1,j)=(8m²-1)I_(m,j)-4V² I_(m,j+1)-4(m-1)^2 I_(m-1,j).
```

`raw_moment_transform.cpp` implements this exactly with directed rounding; `test_raw_transform.py` verifies 33 outputs against 12 signed exact rational atoms. Set Z_k=12 I_(k+2,0)-(k+1) I_(k+1,0)+(k+3) I_(k+3,0), A_k=binom(k+3,3)^2 (Z_k - k Z_(k-1)/(k+3))/(k+3), and B_k=((k+1)^2 A_(k+1)+[-k²+(8-4b)k+9-2b(b-1)]A_k)/2. This is the retained raw angular identity of the repository. After forming the repaired A_k and the ideal p4 coefficients, the direct-prefix theorem reduces the small-L sliver to a bivariate rational positivity test with the action error allowance (including seed errors and truncation). That test has not been executed for the realized coefficients.

All statements remain conditional on the inherited input certificates; no new assertion of positivity is made for any uncomputed channel or interval.

## Final executed slab result

The full `ledger_23.json` has both channels PASS_CONDITIONAL on [2.3,2.6]. All 10 p=3 boxes (1,307,286 scalar evaluations) and two p=4 right-edge boxes (247,104 evaluations) were enclosed. The final exact error/floor ratios are below 0.316041 and 0.568905; the lower normalized actions exceed 2.0518e-6 and 1.2932e-6 respectively. These extend the earlier conditional positivity coverage [2.6,3.65]. The lower region remains open; do not interpret the seven-center moment test as full-spectrum coverage. The original long subset benchmark was interrupted and replaced by the completed accelerated seven-center high-order test. Its partial log contains only the completed 36 numerical comparisons. The full 2500-digit subset output is included for reproducibility. All processes used in the computation have completed or been stopped.
