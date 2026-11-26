# Sliver extension, L in [2.9,3.65]

This is a **conditional certificate**, not a proof of the full sliver. It retains the angular-action identity, normalization, EPSW repaired spectral weights and residual masses, low-frequency and localization remainder certificates, and the LRINF reference lower bound. On [3.5,3.65] the contour bound for a<=4500 is imported from `ads4/code/analytic3/cache/contour3_3.5_3.65_k4.json`; the earlier high-twist certificate supplies larger twists. On [2.9,3.5] all four contour integrals plus the short real extension are newly bounded with the code here, uniformly in a>=998. No physical spin enumeration is used. Neither the uncovered L<2.9 domain nor the earlier large-L certificates are proved anew by this archive.

## Result

The exact rational sign comparisons are in `full_ledger.json` and `full_ledger.log`.

| L cell | b=0 E/F | b=4 E/F | common reference lower bound |
|---|---:|---:|---:|
| [2.9,3.2] | .079228620 | .088400113 | 2.5e-6 |
| [3.2,3.35] | .246056183 | .297553455 | 1.87e-6 |
| [3.35,3.5] | .168162195 | .191352800 | 1.87e-6 |
| [3.5,3.65] | .137500907 | .148859590 | 1.98e-6 |

See `NORMING_BASE.md` for the exact scalar recurrence, base multiplier bound, sequential norming proof, low-frequency bound, normalization and retained assumptions. The prior package was independently rerun and reviewed by the repository coordinator (§109). Those conclusions and inputs are not independently re-proved here.

## Changes in this run

1. **Physically coupled coordinate:** the cap is h <= d^2/[998^2 (1+sqrt(1-d))^2]. In coupled configurations `h0,dh` encode x=h/d^2 rather than h. Use x between zero and the upward rounded maximum cap on the d interval. Rational scalar nodes pass the exact product h=d^2*x to `scalar.cpp`. The old rectangular (d,h) bounds are valid but loose; `resid_combo_2.9_3.2_results.json` uses their first two frequency bands and the coupled bounds for the last three. For example the dressed b=4 [850,950.1] P bound improves from 17584.81 to 25.090 (different but valid norming grids). `major_rect.py` encloses complex coefficients with directed rectangles, checking pole gaps; only one parameter axis is complex at a time. For the ratio bound r_j, the squared modulus envelope has the form (A+t^2)^2/((j+1)^2(C+4t^2)); its derivative has sign equal to that of C+2t^2-2A, so the maximum over a real interval of t^2 is at an endpoint. The disk enclosure is retained as an alternative.
2. **Contour extension:** for the complete finite signed connection sum, define g(w)=P((w+w^(-1))/2) in normalized frequency. On the annulus 1<=|w|<=R, the maximum modulus is logarithmically convex in log|w|. Therefore the real norm U and an outer ellipse bound M yield an inner bound U^(1-lambda) M^lambda with lambda=log(r)/log(R). This applies separately to P and Q, holding d and x fixed real. `ellipse_extend.py` checks the exact containing ellipse for each target rectangle, bounds the full J=1500 sum on the outer ellipse, excludes poles, and rounds log and exp outward. This also encloses tiny decimal endpoint gaps in some historical configurations (the candidate row search tolerance alone is NOT a coverage assertion). If no row applies, use an absolute recurrence majorant with explicit denominator checks. `contour_direct.py` implements the near-real differences, bottom and left edges. It uses the exact beta difference recurrence and inhomogeneity g beta c(beta c+1+sd), g=h/d^2; a negative-binomial factorial-moment telescoping estimate seeds it (see source comments). The coefficient ratio tail bounds the remaining infinite-twist series. On top and right the physical and infinite-twist terms are individually bounded, with the exact polynomial weight and Gamma duplication/Euler integral envelopes. The contours have heights 10.25 (p=3) and 11.25 (p=4). The rational right endpoint 998.5 exceeds Lambda=sqrt(997002); `extra()` explicitly bounds the added real segment. No extrapolation of the real-axis bound off the real axis without a checked ellipse.
3. **Gaussian tail:** use the existing repaired spectral representation and shifted rational factor. In five disjoint groups with positive center ceilings 1025,...,1029, the exact binary64 absolute amplitude sums are < 1220000,.011,.00053,.000024,.00000026 respectively. All 3564 centers and widths were checked against the source inventory (width in (.449,.901)). The conjugates and twins, shifted Gaussian factor <13, shifted theta factor <1.02 and QR_- factor <20 give the common prefactor 1062; the remaining logarithmic term retains the exp(-t^2/3) envelope. Integrate each group by the log-concave tangent bound, using directed powers of two for exponential tails. The integral mass above 1040 is <1.26e-37. The old, coarser Gaussian envelopes are retained in historical budgets, and remain valid. The cutoff correction is independently charged using the retained 1.34e-14 L1 bound.
4. **Exact coverage/sign checks:** `assemble_full.py` checks d slabs contain exp(-L) for the whole interval using directed exponentials and checks the contour intervals agree with the residual intervals. `verify_inputs.py` checks that the source floor rows tile larger beta ranges, where cosh(beta/2)^2=exp(L), and that their lower bounds exceed the chosen floors. It checks the exact weight and mass ceilings. In this session the source data were retrieved from GitHub and the checks performed, rather than using a local MAJOR checkout; the bundled verifier is provided for replay on a checkout. All sign decisions use rational bounds; printed decimal results are illustrative.

## Replay

Requires Python 3, a C++ compiler, and GMP headers/library. In this directory:

```
python assemble_full.py
python check_identity.py
python check_base_identity.py
python check_lowfreq.py
python verify_inputs.py /path/to/MAJOR
```

`assemble_full.py` recombines existing computed enclosures; it does not verify the scalar samples from scratch. For a full recomputation of the new slabs, run `bash reproduce.sh` (can be costly; `JOBS` sets worker count). The optimized majorant implementation may choose different norming degrees or grids from the historical records; fresh bounds need not be bit-identical. The `*_results.json`, `ct_*.json`, `*_budget.json`, and `*.log` are the actual executed historical outputs. A residual-only `PASS_CONDITIONAL` field in intermediate files is not an end-to-end positivity certificate; use `full_ledger.json` only. Historical failed/loose rectangular tests do not count as coverage.

## Remaining work

The uncovered low-L interval is approximately 0.06<=L<2.9 (exact lower endpoint depends on the prior sliver definition). Below .8 the present localization/connection approach has additional unverified hypotheses. The lower-L region is not certified by an overlap or interpolation of the existing results. The existing claims for L>=3.65 and the retained upstream inputs still have their own review obligations. No new large-L proof is included here.
