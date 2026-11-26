# Signed transfer checker — standalone reproduction and width benchmarks

## Scope

This package checks bounds on the **two signed connection-transfer coefficients** for the truncated connection sum of length J0=1500. It is spin-independent. It also checks an independent low-frequency residual estimate. A `PASS` in this package is NOT a certificate of positivity of the full angular action. In particular, the physical base multipliers, the contour assembly, the high-frequency residual integral and the full parameter cover remain to be certified. The localization theorem and the angular-action identity are inherited, not independently rederived here.

All transfer coefficients, parameter endpoints, roundoff bounds, pole exclusions, and Cauchy remainders are computed with BigInt arithmetic. No external libraries or network access are needed. Node.js 22 and Python 3 (standard library only) were used. Timing is wall time per (p,b) case, measured on an Intel Xeon Platinum 8573C container. It is not a benchmark on the user's laptop. The long order-40 case ran concurrently with other package checks.

## Quick reproduction

```
node check_transfer.js pilot.json rerun_pilot.json > rerun_pilot.log
node check_transfer.js combo.json rerun_combo.json > rerun_combo.log
python3 check_lowfreq.py
python3 check_identity.py
```

`pilot.json` reproduces the original six channel-specific box bounds. `combo.json` is the largest combined width test at L approximately 3.2. `lower40.json` is the all-twist test near L=2.5. `complex_top.json` is a top-line transfer test. `reproduce.sh` reruns the whole benchmark set. `benchmarks.csv` summarizes recorded results; the paired JSON and log files preserve the full 80-decimal-place enclosures, input parameters, errors and elapsed time. `majorant_sweep.json` and `combo_preflight.json` deliberately set an extremely stringent remainder tolerance so the checker stops before constructing the large Taylor polynomial; those `FAIL_REMAINDER` entries are *preflight* diagnostics, not unsuccessful positivity tests.

## Input schema / launcher interface

`node check_transfer.js INPUT.json OUTPUT.json` takes an array of boxes, prints JSON lines to stdout, and writes an array to OUTPUT.json if supplied. Use decimal strings for `d0`, `dd`, `h0`, `dh`, `nu0`, `dnu`, `maxRemainder`, and `maxTransfer` (not scientific notation). A box is

- d=d0+dd*t1, h=h0+dh*t2, nu=nu0+dnu*t3, each t_i in [-1,1];
- optional `nuImag`: fixed imaginary part of nu (default zero); optional `nuAxis: "imag"` varies the imaginary rather than the real part of nu, with nu0 the fixed real part. `dnu` is the half-width in the selected frequency direction. Complex frequency boxes have been tested at one top-line location only;
- `m`: total Taylor degree, `K`: number of connection indices computed in signed Taylor arithmetic, `R`: integer Cauchy polydisk radius (>=2);
- `channels`: array of [p,b] pairs with p=3 or 4, b=0 or 4; default all four. For p=3 the physical dressing is included in the numerator; p=4 is undressed;
- `fullTwist: true` makes the checker assert that h in [0,hcap] covers **every a>=998** for every d in the input box;
- `hSplits`: number of exact subintervals in the Taylor polynomial's h direction. This subdivides the final polynomial only, NOT the recurrence or physical spin. Default no subdivision. Values 16 and 32 were useful;
- `maxRemainder` is a fail-fast numerical tolerance for the common Taylor-plus-connection-index remainder. `maxTransfer` is an optional target on each of |P| and |Q|, used to label PASS or FAIL_TRANSFER. Both are adjustable *engineering thresholds*, not physical sign criteria.

An ERROR, FAIL_REMAINDER, or FAIL_TRANSFER is NOT a counterexample. A launcher should subdivide or change m,K,R and rerun. Check every required [p,b] channel and preserve all boxes and bounds. The optional `make_grid.py` generates a proposed geometric d/frequency cover with a rigorously outward-rounded physical hcap. Its boxes are **unverified**. Do not launch the whole grid before measuring actual failures and cost near L=0.8.

## Full-twist coverage

Put c=1-d, s=2b+1, N=n+a, D=N(N+s), h=1/D. The identity

    h/d^2 = 1/[a^2*(1+c+sqrt(4c+s^2*d^2/a^2))]

implies h <= d^2/[998^2*(1+sqrt(1-d))^2]. The right side increases with d, so using d_max gives a single constant hcap for the entire box. The checker rounds the square root down and the resulting quotient up; it requires h0=dh and 2*h0>=hcap for `fullTwist`. Equivalently y=h/hcap is in [0,1]. This avoids computing the exact curved h_max(d) and requires no separate recurrence runs for different twists. The box includes some nonphysical (d,h) combinations, which is safe for an upper bound. The h=0 boundary is included algebraically.

## Mathematical normalization

Let alpha=p-1/2+i*nu, beta_j=alpha-j, and S(beta)=L((1-z)^(-beta)) with the physical terminating ladder. Set X_j=d^(alpha-j) S(alpha-j) for j>=-1. The exact three-term identity is

    (1-h*beta*(beta-s))*X_{j+1}
      = (1+h*beta*(s-2*beta-1)/d)*X_j + h*beta*(beta+1)*X_{j-1}/d^2.

Define c_0=1, r_j=c_{j+1}/c_j=(j+1/2-i*nu)^2/((j+1)*(j+1-2i*nu)); set r_-1=1. **Individual c_j are meromorphic in complex nu**, not entire. The checker explicitly verifies the denominator gaps on the enlarged polydisk.

Use Z_j=c_j*d^j*X_j for j>=0 and Z_-1=X_-1/d. Then

    D_j Z_{j+1}=A_j Z_j+B_j Z_{j-1},
    D_j=1-h*beta_j*(beta_j-s),
    A_j=r_j*(d+h*beta_j*(s-2*beta_j-1)),
    B_j=r_j*r_{j-1}*h*beta_j*(beta_j+1).

The two basis solutions start (Z_0,Z_-1)=(1,0) and (0,1). Thus the sum of their signed coefficients gives P and Q_unscaled, and Q=Q_unscaled/d when multiplying X_-1. For p=4 the output gives |P|,|Q| in

    sum_{j<J0} Z_j=P X_0+Q X_-1.

For p=3 let q=b(b+1), tau=h/(1+qh), A0(beta)=9/2-b(b-1)-beta*(beta+8-4b)/2, A1(beta)=beta*(beta+9-4b)/2. The dressed sum is

    sum_{j<J0} [(1+tau*A0(beta_j))*Z_j + tau*A1(beta_j)*r_{j-1}*Z_{j-1}]
      = (P X_0 + Q X_-1)/(1+q*h).

The checker constructs the numerator, including the factor 1+qh. It bounds |Q_unscaled|/d_min. The denominator is >=1 for physical real h; for any use with complex h, it must be bounded separately.

Thus for real nu>0 and real physical d,h the full normalized **truncated** dressed F3 kernel is bounded by

    2*d^(1/2)*|B(nu)|*(P_bound*M0+Q_bound*M1)/(1+q*h)

where M0>=|X0|, M1>=|X_-1|, B(nu)=Gamma(2i nu)/Gamma(1/2+i nu)^2. The bound for p=4 is 2*d^(-1/2)*|B(nu)|*(P_bound*M0+Q_bound*M1), before the extra tau*Lambda^2 in the physical action. The factor 2 is for the two conjugate branches on the real frequency axis. On complex contour pieces only use the relevant outgoing branch and |d^(-alpha)|, not the real-axis conjugation formula. The localization remainder must still be included. **Do not substitute M0=M1=1:** use appropriate independently valid Theorem E', Lemma F, or other power multiplier bounds at alpha and alpha+1, including helicity, turning points, and complex exponents. For the ideal contour comparison, one also needs the difference from the infinite-twist expression, or the existing bounds on the near-real pieces. A kernel absolute bound alone does not certify the ideal signed integral.

## Proof of arithmetic enclosure

Work in the truncated multivariate polynomial algebra, discarding terms of total degree >m. Its coefficient l1 norm (real and imaginary absolute values) is submultiplicative. All stored polynomial coefficients have scale S=10^80. Sparse products in denominators are exact (asserted at runtime). For A P+B Q over D, let D_0 be the constant coefficient. The inverse norm is bounded by

    |1/D_0|_1 / (1-|1/D_0|_1*||D-D_0||_1),

with every bound rounded outward; a nonpositive gap causes ERROR. Coefficients are solved triangularly in increasing total degree, with toward-zero rounding in each component. The residual in each coefficient is at most two fixed-point units times |D_0|_1. Incoming errors propagate through the exact sparse coefficient norms; multiplication by the inverse norm gives a directed bound on the total coefficient error. This **corrects the earlier prototype's mixed Euclidean/l1 inverse bound**. All baseline cases were rerun with the corrected rounding tracker. `check_identity.py` independently checks the weighted recurrence, dressing and signed transfer normalization by Gaussian-rational arithmetic on 48 finite physical ladders (576 exact step comparisons).

For the analytic remainder, use the unit polydisk in t, enlarged to radius R>1. The routine `major` proves a bound M_R on each signed coefficient function's magnitude (Q_unscaled), by a positive majorant of the exact recurrence including all J0 indices and the dressing. It uses |D_j| >= 1-h_max*|beta_j|*|beta_j-s| and the lower triangle inequality for the r_j denominator, and fails if any denominator is not bounded away from zero. Magnitudes of all complex affine quantities are bounded with directed integer square roots and triangle inequalities. With two basis initial entries bounded by 1, this controls both transfer coefficients. The absolute majorant also bounds the partial sum through K-1. By multivariate Cauchy, the omitted degree >m in that partial sum is bounded by

    M_R * binom(m+3,2) * R^(-m-1) / (1-(m+4)/((m+2)*R)).

The same majorant at R=1 bounds the omitted connection indices K,...,1499 separately. These two errors are added to the total polynomial roundoff error. Q is divided by d_min only after adding all errors. Final coefficient Euclidean magnitudes are rounded upward with integer square roots, so their sum bounds the polynomial on [-1,1]^3.

For `hSplits=N`, each final polynomial is shifted exactly to h-variable t2=(2*l+1-N+u)/N, u in [-1,1], for l=0,...,N-1. For each group of fixed t1,t3 exponents, the h-polynomial is expanded by the binomial theorem using a common integer denominator N^degree. Sum upward-rounded magnitudes of the shifted coefficients; maximize over all N subintervals, then add the unchanged error bounds. This preserves the h-direction phase better than a single monomial coefficient norm. It does not change the number of recurrence evaluations.

## Low-frequency residual theorem and proof

On 0.8<=L<=3.65 the exact state relation d=t*(2-t+s/N)/(1+s/N), t=a/N, implies 0.012<=t<=0.26. Consequently n>=a. Let W(theta)=|1-rho exp(i theta)|^(-a-1), rho=1-t, and B_beta~Beta(n+s,a), v=1-B_beta. The exact coefficient/Beta-Cauchy identity gives |L(F)| <= K0 E_beta E_W |F(B_beta rho exp(i theta))|. Stirling's two-sided bounds and sin(theta/2)>=theta/pi give K0 <= (pi/2)*sqrt(1+2/a)*exp(1/(6a)) <1.6 (n>=a). Here E_W is normalized angular measure. For any z in the unit disk, the retained Mehler bound is |F3(z;nu)| <= (pi/2)|1-z|^(-3) exp(nu |Im beta(z)|), with beta(z)=2 atanh sqrt(z). At real nu, the imaginary phase satisfies tan |Im beta(z)|=sqrt(2(|z|-Re z))/(1-|z|). The differentiation uses radius delta=t/200; tau<=t^2/a^2 yields a factor <1.081 via 1+tau(2/delta^2+5/delta+8). This holds for both channels.

Let r=B_beta*rho, X=2 sqrt(rho) sin(|theta|/2)/t. Then W=t^(-a-1)*(1+X^2)^(-(a+1)/2) and |Im beta(z)|<=atan(t*X/(1-r)); sqrt(B_beta)<=1 was dropped. Use the three regions below:

* Principal region v>=0.75t and r>=1/4: 1-r>=1.555t. On the enlarged Cauchy disk |1-zeta|>=1.55t, and the derivative |beta'(zeta)|=1/(sqrt(|zeta|)*|1-zeta|) bounds the phase increase by 0.0065. Hence the phase is <=atan(X/1.555)+0.0065. d<=2t. In the X coordinate the normalized measure has density proportional to (1+X^2)^(-(a+1)/2)*(1-X^2/Xmax^2)^(-1/2), Xmax=2sqrt(rho)/t. The denominator is at least .999 sqrt(pi/[2(a+1)]), by integrating the Gaussian to X=1. For X<=1/2 the Jacobian factor is <=1.004 and log(1+X^2)>=7X^2/8. Completing the square and enlarging the numerator integration to the whole real line gives central expectation <=2.16 exp(nu^2/4200), using k>=1.55 and a>=998. The tail X>=1/2 contributes <1e-14 for nu<=200: its exponent nu atan(X/k)-(a+1)log(1+X^2)/2 decreases on this interval, and the angular normalization gives prefactor pi sqrt(a)/t exp((a+1)/(2a)). At a=998 the latter times the damping at X=1/2 bounds all a>=998. The checker evaluates the bound with atan(0.5/k)<=0.5/k. Combining gives a coefficient <12.613 in front of exp(.0065nu+nu^2/4200).
* Lower beta region v<0.75t and r>=1/4: by gamma Chernoff for v~Beta(a,n+s), P(v<0.75t)<exp(-36), since m=a/(N+s) and v/m<=.752, and I(v;m)>=-log(v/m)-(1-v/m). Here 1-r>=t, enlarged gap>=.995t, phase increase <=.0101, and the same angular proof with k=1 gives central <=2.16 exp(nu^2/[2*(7/8)*999]) and angular tail <=1. Bounding nu by 200 makes this entire region <.001, uniformly.
* r<1/4 implies v>.662. Since m<=.26, two-sided gamma Chernoff gives P(v>.662)<exp(-1000). On the enlarged disc, |Im beta|<=pi/2 and |1-zeta|>=.995t. Its contribution is negligible, as checked by directed arithmetic.

For the Chernoff formula use independent gamma variables of shapes a,n+s, and optimize the MGF of (1-v)G1-vG2. It gives exponent -a[log(m/v)+(1/m-1)log((1-m)/(1-v))] on either appropriate tail. The code verifies all numerical inequalities used in the three regions using exact rational log series and outward-rounded exponential series. Thus

    d^3 |L((1+tau D_b) F3)| <= 14 exp(0.0065 nu + nu^2/4200), 0<=nu<=200.

The normalized physical action using w^3 differs by a positive factor <=1, so the same bound suffices. `check_lowfreq.py` integrates the retained absolute residual masses using this envelope with the 2/Lambda <=2/998 prefactor. The interval [0,1] uses 1.9e-6, from the EPSW log's [0,.2] and [.2,1] masses. The other rows use conservative ceilings from the source logs, inflated by 1.002 to cover printed rounding. The cutoff L1 error <=1.34e-14 (from the stronger weighted estimate in TAIL_AUDIT) is included. The minimum inherited F_cert across the requested interval exceeds 1.87e-6. The result through nu=160 is E<1.46e-7<0.08 F_cert. The logged ceiling 200.1 in the last row is used only as a monotone upper envelope for frequencies <=200.

Source inputs (branch claude/eloquent-cori-nftfjy in anshM123/MAJOR):
- `ads4/code/analytic3/logs/resid3_3.55_3.65_p0.log`, `..._p1.log`: per-frequency EPSW mass entries; `residual_masses.csv` stores these 199 rows. `extract_masses.py` checks them against a local checkout.
- `ads4/code/epswcert/cache/ew_cert_1000_re_rep.txt`: near-zero and low-frequency residual masses.
- `ads4/TAIL_AUDIT.md`: cutoff correction (weighted L1 upper bound 1.333e-14), and retained K0/Mehler context as in `ANALYTIC_SLIVER3.md` / `FULLPOS3.md`.
- `ads4/code/lrtight/cache/lrt3_rcert.json`: certified floor, inherited.

This package does not recalculate the EPSW integrals, floor certificate, or angular-action identity. The low-frequency result is conditional on those inherited inputs.

## Observed width scaling

The original small boxes pass; new widest combined real-frequency test near L=3.2 covers d=0.0407622 +/-0.000815 (about 4% full width), h from zero to 4.50e-10 (all twists), and nu=900 +/-10. Order 28, K=200, R=8, hSplits=16: dressed b=4 P<5.466, Q<2.062, remainder<8.96e-5. About 27.65 seconds for this channel. At L=2.5 the 2%-wide d, all-twist, 10-wide nu test passes at m=40,R=10,K=260,hSplits=32, but costs about 108 seconds for dressed b=4. The attempted all-twist box near L=0.8 fails the absolute majorant at order 40 even with positive denominator gaps: pole-free does not mean a small Cauchy remainder. Complex-frequency pilot at nu=900+6.25i +/-5 passes with m=24, full twist, narrow d. Wider tests were only executed for p=3,b=4, not all four channels.

If the 4% d and 20-wide nu sizes could be used everywhere, ceil(2.85/log(1.02/.98))*ceil(990/20)=72*50=3600 geometric boxes before channel multiplication, compared with hundreds of millions of narrow 3D boxes. At 27.6 s/channel and four channel computations this is about 110 core-hours, purely an optimistic extrapolation, NOT an established grid runtime or sign certificate. At lower L, tighter adaptive boxes, higher degrees, sharper analytic majorants, or another resummation may be needed. Do not assume the full 0.8–3.65 band passes from these tests.
