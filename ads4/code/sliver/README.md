# Sliver positivity: independently executed fixed-spin mass cells

**Scope:** two previously uncovered physical spin channels in `anshM123/MAJOR`, not the complete spin/twist continuum. The proof uses the retained angular-action identity in `ads4/TAIL_CERT.md`, §1 (E5), and the EPSW residual and cutoff bounds. It does not use Lemma C or the complex reference truncation.

## Reproduce

Requires a C++17 compiler and GMP (`gmpxx`, `gmp`). For the main completed wider-cell run:

```sh
g++ -O3 -std=c++17 check_wide.cpp -lgmpxx -lgmp -o check_wide
./check_wide 1002 0
./check_wide 1005 4
```

`check_cell.cpp` reproduces the narrower pilot cells. The source logs record the integer sign tests and decimal approximations to the bounds. The broader cells use `D0 ± 400000`, which contain every `a` in `[2951,3049]` for each specified spin. The near-gap sources use a different scale and 2400-bit precision; `check_gap_low.cpp` is explicitly a LOW-BAND check, not by itself a full repaired-weight sign certificate.

## Mathematical normalization

Let `n=J-b`, `s=2b+1`, `q_b=b(b+1)`, `N=n+a`, `D=N(N+s)`, `C0=n(n+s)`, `T=D+q_b`, `T0=D0+q_b`, `w=1-J(J+1)/T`. For the terminating ladder `U_k=prod_{j<k}(C0-p_j)/(D-p_j)` with `p_j=j(j+s)`, the certified normalized quantity is

    B = (T/T0) H_ang,b / w^3.

Every multiplier in this normalization is positive. In the exact action formula, `G_A` is built from the repaired density and `G_D` from the exact ideal k=3 polynomial. The first moment operator contributes the coefficients

    (D_b F)_k = ((k+1)^2 g_{3,k+1} + (-k^2+(8-4b)k+9-2b(b-1))g_{3,k})/2.

The code multiplies the physical action by the positive factor `T/T0` before forming the polynomial; therefore its Taylor coefficients have no extra varying dressing denominator.

## Certified arithmetic and analytic remainder

On the mass cell write `D=D0+W*t`, `|t|<=1`. Put `E_k=SCAL^k U_k` and `G_{p,k}=SCAL^{-k} g_{p,k}`. `SCAL=15` for the wide and pilot cells, and 3 for the near-gap checks. The code implements the exact recurrences for E in t and G in the Chebyshev coordinate x=2 nu^2/V^2-1. It assembles the signed kernel before taking the coefficient l1 norm. All coefficient operations use GMP integers scaled by `S=2^P`; `div0` and `shr` truncate towards zero. Every such rounding is bounded by one unit of `1/S` (or an explicit multiple when several operations are combined).

The retained E coefficients satisfy `A+W<=D0-p_k`, with `A=SCAL*(C0-p_k)`. Consequently the l1 rounding error of each row is at most `2(n+1)(Q+1)/S`. The G recurrence tracks its l1 rounding error exactly in integer units. The coefficient assembly error follows by the product rule for the l1 norms. The moments of the ideal NF02 weights are evaluated by the exact integration-by-parts recurrence for Chebyshev moments, with a separate error recurrence. The code performs the final strict sign comparison with integers, not the printed decimal values.

For the wide cells, on the larger mass circle `|t|<=R=10`,

    |E_k(t)| <= mu^k,   mu = max(1, SCAL*C0/(D0-R*W)).

The bound follows because `(C0-p)/(D0-R*W-p)` decreases in p when `D0-R*W>C0`. The dressed kernel weights are bounded by `3 mu^k`; this constant is checked by the integer assertions. The Cauchy bound on the discarded mass Taylor coefficient l1 sum is `3 mu^k/((R-1)*R^Q)`. The code sums it against the rigorous spectral coefficient norms; all arithmetic is outward bounded by integer ceilings. The narrow and gap runs have `mu=1`, `R=100`.

For the wide and full-band cells, the repaired-density error is bounded by

    (2/998)*(3e-6)*||K_2|| + 1e-20,

where `||K_2||` is the signed bivariate coefficient l1 upper bound on the mass and frequency rectangle. The `3e-6` constant dominates the EPSW bound `2.681358e-6` and the audited cutoff error. The spectral tail above 1100 is bounded separately: the repaired coefficient inventory has sum of absolute RBF coefficients < 1.22e6, centres <1029 and widths between .449 and .901. The rational factor |Q R_-|<=16 and the cutoff modulus bound give a Gaussian tail bounded by exp(23-(nu-1029)^2/(2*.901^2))+exp(-nu^2/3). For the wide cells use the conical generating function on |z|=.1; for the gap cells use |z|=.5. The Laplace representation of P and Cauchy's coefficient estimate yield a bound <= 1000 exp(pi*nu/2) for the normalized dressed physical kernel, since SCAL^{-1}/|z|<=2/3 and the corresponding weight bounds hold on the real mass cell. The Gaussian exponent at nu=1100 is < -1376 before the prefactor, decreasing thereafter. The resulting normalized tail is far below 1e-20.

## Retained source inputs

- https://github.com/anshM123/MAJOR/blob/claude/eloquent-cori-nftfjy/ads4/TAIL_CERT.md — exact angular action, dressing and generating-function identity.
- https://github.com/anshM123/MAJOR/blob/claude/eloquent-cori-nftfjy/ads4/TAIL_AUDIT.md — audit, including cutoff error <=1.333e-14 in the exp(.02 nu)-weighted L1 norm.
- https://github.com/anshM123/MAJOR/blob/claude/eloquent-cori-nftfjy/ads4/code/epswcert/cache/ew_cert_1000_re_rep.txt — repaired residual mass, including the >Lambda band bound 3.1175e-10.
- https://github.com/anshM123/MAJOR/blob/claude/eloquent-cori-nftfjy/ads4/prl/data/nf02_psi.json — exact NF02 polynomial inputs.

The results do not independently rederive the retained physical action formula or re-run the original EPSW certificates. They do not certify a continuum in spin, and they do not close the large-spin sliver.

## Completed near-gap bandwise certificate

The full-band L1 coefficient test alone fails for the near-gap cells (it takes an absolute frequency norm over an unnecessarily large interval). This is **not** a physical negative action. Three separate kernel norms are used instead:

| J,b | V=1000 | V=1040 | V=1100 |
|---|---:|---:|---:|
| 1002,0 | 15.712 | 5326.947 | 7.308e12 |
| 1005,4 | 131.893 | 80352.761 | 1.099e14 |

For the frequency bands [0,1000], [1000,1040], [1040,1100], use residual masses 3e-6, 4e-10, 1e-22, respectively. The second bound follows from EPSW's >Lambda mass and the audited cutoff correction. For the third, the direct shifted Gaussian bound above and `(1040-1029)^2/(2*.901^2)>74` give an integral bounded by `(.901^2/11)*exp(-51)` plus the negligible log-term tail, hence <1e-22. The mass interval [998,1002] obeys |t|<=1.004; all polynomial l1 and Cauchy remainder bounds grow by a factor <4, as checked with exact rational arithmetic in `verify_gap_bands.py`. This factor is included in the final sign comparison.

Reproduce the near-gap checks:

```sh
for v in gap gap_low gap_mid; do
  g++ -O3 -std=c++17 check_$v.cpp -lgmpxx -lgmp -o check_$v
  ./check_$v 1002 0 > ${v}_b0.log 2>&1
  ./check_$v 1005 4 > ${v}_b4.log 2>&1
  # The gap (V=1100) exit code 2 is expected: its unsplit residual bound fails.
  # The log still contains the rigorous high-band kernel norm.
done
python3 verify_gap_bands.py
```

The bandwise verifier uses exact Fraction arithmetic, padding the printed double conversions by 1e-12 relative to cover printing/conversion error. The ideal lower bounds include a 1e-90 numerical-error allowance, greater than the computed low-band errors. Both combined signs pass:

| J,b | Physical twist interval | Lower bound on normalized action B |
|---|---|---:|
| 1002,0 | [998,1002] | >3.5714e-5 |
| 1005,4 | [998,1002] | >2.9786e-5 |
| 1002,0 | [2951,3049] | >2.4031e-5 |
| 1005,4 | [2951,3049] | >2.2630e-5 |

All four are newly executed continuous mass intervals. The entire sliver is still open: this implementation works at fixed physical spin, and a uniform certificate in the compactified spin coordinate is not implemented. The repository's tail audit certifies L>=20; do not assume all remaining parameters have been covered.
