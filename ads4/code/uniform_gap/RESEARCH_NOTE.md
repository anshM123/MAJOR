# Uniform-gap continuation of the NF02 / EE angular certificate

## Scope and provenance

This is a partial extension of the fixed-gap result, not a proof of full-spectrum positivity at all gaps. The new family is different from the published endpoint-repaired Gaussian inventory, even at gap 1000. The retained angular-action and light-side identities, normalization, and the other sum-rule hypotheses are those in `anshM123/MAJOR`, branch `claude/eloquent-cori-nftfjy`, particularly `ads4/prl/sec_fullpos.tex`, `ads4/prl/supplement.tex`, `ads4/K_B2A4.md`, `ads4/T0_CERT.md`, `ads4/GAP_LARGESPIN.md`, and `ads4/TAIL_CERT.md`. Exact NF02 coefficients are transcribed from `ads4/prl/data/nf02_psi.json`, blob 7876a9063dc0404c25b81fda538ff26ecde9175d. This work does not rederive the CFT operator identities or establish positivity of the quartic functional.

**Results:** analytic spectral realization for every real gap >=1000; uniform light coefficients; positive exact finite-mass angular action for every physical spin through 64 and for all spins with hat-b <=40; uniform realized hybrid floor (retaining T0 and LRINF inputs). The remaining full-action sign problem is J>=65, hat-b>40. A conditional quantitative rigidity formula is given below. No repo files were changed.

## 1. Exact inverse of the de-shift equation

Let L=Lambda, t=L^(-2), D=997002. Write f(x)=psi2(|x|/L) on |x|<=L and zero outside. Use the rational functions R_-, R_+ in `check_cutoff.py`, with R_+(z)=R_-(-z). Define

P(z)=(z^2+1/4)^2 (z^2+9/4)(z^2+25/4)(z^2+81/4),
F(z)=P(z) cosh(pi*z/2) Gamma(5/4+iz/2) Gamma(5/4-iz/2).

Gamma recurrence and the checked rational identity give F(z+i)/F(z-i)=R_-(z)/R_+(z). Hence B(z)=F(z+i)/R_-(z)=F(z-i)/R_+(z) is even and real on the real axis. The polynomial identity P(z+i)/R_-(z)=-z(2z-7i)^2(2z-3i)^2(2z+3i)(2z+7i)^2/2048 shows that apparent rational poles cancel; the remaining Gamma factors are analytic for |Im z|<3/2.

Put W=f/F on the real line, and

q(z)=(1/4) integral_R W(v) sech(pi*(z-v)/2) dv,    chi(z)=B(z)q(z).

The convolution is analytic in the open strip |Im z|<1, is even, and is real on the real axis. For 0<=y<1, q(x+iy)+q(x-iy) is convolution with (1/2) Re sech(pi*(x-v+iy)/2). Substituting s=sinh(pi*(x-v)/2) turns this measure into the mass-one Cauchy approximate identity (1/pi)*cos(pi*y/2)/(s^2+cos(pi*y/2)^2) ds. W is compactly supported and Lipschitz. The boundary imaginary parts exist as principal values, using subtraction of W(x) at the simple singularity. Thus q(x+i)+q(x-i)=W(x), and

R_-(x-i)chi(x-i)+R_+(x+i)chi(x+i)=f(x).

This identity is exact, not a fit. Boundary values are understood as limits from the strip. Cutoff family used here:

 theta(z)=[erf((z+L+2000)/10)-erf((z-L-2000)/10)]/2.

It is entire, even, and decays super-Gaussian on horizontal strips. The same strip analyticity mechanism as in the original HK argument applies, subject to the retained non-weight assumptions. Uncut chi is only used as a reference for the estimates and finite angular moments, not assumed to be an admissible OPE-interchange weight without cutoff.

## 2. Light coefficient calculation

Let mu[k,j]=integral_0^1 psi_k(p)*p^j dp. For the uncut family the Einstein light coefficient is E=16*LA, LA=integral_0^1 psi2(p)/(p^2+9/(4L^2)) dp. For the cubic squared coefficient use

Q(z)=z^6-(245/4)z^4+(9619/16)z^2-51575/64.

The raw rational k=2 light polynomial in `K_B2A4` equals R_-(z)Q(z+i)+R_+(z)Q(z-i), exactly (symbolic check). Shift the contour using section 1; the Einstein test function has no pole between the two lines and the polynomial test function is entire. The exponential tails justify the shifts for these tests. In the repository's light-side normalization the k=3 polynomial is p^4-(35/2)t*p^2+(569/16)t^2. It follows that

K(t)=-2[mu[2,6]-(245/4)t mu[2,4]+(9619/16)t^2 mu[2,2]-(51575/64)t^3 mu[2,0]]+2[mu[3,4]-(35/2)t mu[3,2]+(569/16)t^2 mu[3,0]].

Exact coefficients in increasing powers of t:

* 888656439590000/888656439561891
* -2176144896481250/67209310555101
* -1537157051698780625/1683770096012004
* 78949596683828125/4716442845972

K(t)>0.99996 on 0<=t<=1/D, with K(1/D)=0.999967523127696... . Set A=psi2'(0), g(p)=(psi2(p)-A*p)/p^2, eps=3/(2L). Then

LA=(A/2)log(1+eps^(-2))+integral g-eps^2 integral g/(p^2+eps^2).

The exact polynomial Bernstein bounds are -511.773<g<=-188.896, integral g=-332.95059691167796... . The error from A log(2L/3)+integral g is bounded by A eps^2/2+(pi/2)eps||g||_infty. Charging the cutoff bound below gives

E_theta <= 3023 log(2L/3)-5307, K_theta >= 0.99995.

If the full angular action is nonnegative on every heavy state at this gap and the retained light-side sum rule (including its treatment of light exchanges and parity) applies, then

sqrt(t4^2+alpha_T^2) <= (240/L^4)*sqrt((3023 log(2L/3)-5307)/0.99995).

This is a conditional all-gap rigidity theorem, not a new unconditional CFT bound. The square-root log is real: A>0 implies LA=A log L+O(1), while K approaches a positive constant. More generally, for a nonzero fixed smooth psi2 with psi2(0)=0, psi2(1)=0 and nonnegative Hankel transform integral psi2(p)J0(bp)dp for all b, the radial function psi2(|x|)/|x| on the unit disk has nonnegative Fourier transform. Gaussian regularization and Fourier inversion imply its absolute value is bounded by its value at zero, so psi2'(0)>0. Thus within this fixed-shape positive-transform class the logarithm cannot simply be removed. This observation does not exclude other functional constructions or gap-dependent shapes.

## 3. Exact angular sign on the infinite gap axis

Use the retained Eq. (fp-target) and generating identity in sec_fullpos.tex. Write m^2=T-C, C=J(J+1), T=(Delta-1)(Delta-2), n=J-b, b=0,4, sigma_b=+1,-1 respectively, w=m^2/T, u=L^2/m^2. Define B_j=C-(b+j)(b+j+1). The physical twist condition implies a=Delta-J-2>=Delta_g-2 and m^2>=a_g(a_g+2J+1), where a_g=Delta_g-2. The ladder acts on z^k by U_k=prod_{j<k} B_j/(m^2+B_j) and vanishes for k>n. Its dressing is the exact operator D_b of sec_fullpos.tex.

The coefficient g[k,p](y) of z^k in (1-z)^(-p) P_{-1/2+i nu}((1+z)/(1-z)), y=nu^2, obeys the three-term ODE recurrence

4(k+1)^2 g[k+1,p]=(8k^2+8pk+4p-1-4y)g[k,p]-4(k-1+p)^2 g[k-1,p], g[0,p]=1, g[-1,p]=0.

Let A_k=2 t^k integral psi2(p)g[k,3](p^2/t)dp and G_k=t^k integral psi3(p)g[k,4](p^2/t)dp (polynomials in t). With gamma=(-k^2+(8-4b)k+9-2b(b-1))/2 the uncut ideal angular action divided by w^3 is exactly

sum_{k=0}^n [u^k prod_{j<k} B_j / prod_{j<k}(1+B_j*t*u)] * [A_k + u/(1+C*t*u)*(t gamma A_k + (k+1)^2 A_{k+1}/2 + sigma_b G_k)].

The formula follows by applying the dressing to the monomials, then shifting the finite polynomial spectral moments using section 1. Clear the positive denominator to get a polynomial in q=D*t and u, with q-degree <=u-degree for each monomial. For J<=32 all 46 physical spin-channel cases pass exact Bernstein/de Casteljau on [0,1]^2. For the next spins the rectangle includes nonphysical masses, and some signs fail there. Use s=998/a_g and zeta=m_gap(J)^2/m^2. Then q=999*s^2/(998+s), u=zeta*(998+s)/(998+(2J+1)s). Clearing the positive remaining factor (998+(2J+1)s)^(n+1) gives an exact polynomial on [0,1]^2. All 48 remaining spin-channel cases through J=64 pass. Total 94; minimum normalized ideal lower bound 2.7185352228409318e-8.

### Unbounded spin and bounded impact

Let v=(J+1/2)^(-2), x=hat-b^2/4, u=x*v, hat-b=(2J+1)L/m, A[j](v)=1-(b+j+1/2)^2*v. For J>=65 and hat-b<=40 use v<=4/131^2, x<=400, q in [0,1]. For a prefix k<=K, clear the positive denominator (1+t*x*(1-v/4))*prod_{j<K}(1+t*x*A[j]). The numerator is

sum_{k<=K} x^k prod_{j<k}A[j] prod_{j=k}^{K-1}(1+t*x*A[j]) * [(1+t*x*(1-v/4+gamma*v))A_k + x*v*((k+1)^2 A_{k+1}/2+sigma_b G_k)].

Here A[j] denotes the ladder factor and A_k the integrated coefficient. All extra terms after physical termination have an exact zero factor; denominator gaps are positive on the entire extended box. Use K=110. Scale q, v and x to unit cube; keep degrees <=14 in q and <=40 in scaled v and keep the entire x polynomial. Discarded powers cannot feed back to lower powers in the backward polynomial recurrence. The positive absolute-coefficient majorant at radii (100,1,1) and (1,16,1), divided by 100^15 and 16^41, bounds the omitted polynomial terms. For the infinite ladder tail, let W_k=X^k times the sum of absolute coefficients of t^k g[k,p](p^2/t) at t=1/D, p^2=1, X=400. For k>=K+1 the recurrence implies W_{k+1}<=alpha W_k+beta W_{k-1} with alpha=X*(3/D+1/(K+2)^2), beta=2X^2/D^2. With r=1/10, alpha/r+beta/r^2<1, hence W_{K+j}<=max(W_K,W_{K+1}/r)*r^j. The weight suprema and quadratic factors of the dressing are bounded by exact Bernstein and geometric-square sums in `impact_tensor.py`. The discarded angular tail is <9.45e-70; polynomial truncation errors are <7.16e-15 and <1.26e-12. Tensor Bernstein coefficients of the retained polynomial prove positivity on the whole cube, one leaf for each channel. After subtracting tails and dividing by the positive denominator upper bound, the normalized ideal lower bounds exceed 0.01519546 in both channels.

## 4. Uniform cutoff estimate in the certified finite-mass regions

Binet's log-Gamma formula (positive real part a) gives |R(a+iy)|<=1/(12a). In particular

pi exp(-2a-1/(6a))(a^2+y^2)^(a-1/2) <= cosh(pi*y)|Gamma(a+iy)|^2 <= 2pi exp(1/(6a))(a^2+y^2)^(a-1/2).

Combining the polynomial factors gives F(x)>=2^-11(1+x)^(23/2), F(x)<=2^13(1+x)^(23/2), |F(x±i)|<=2^14(1+x)^(23/2), |B(x)|<=512(1+x)^(19/2), x>=0. For the shifted Gamma product use Cauchy-Schwarz with a=3/4,7/4; the shifted factor is bounded by 16(1+x)^(3/2). The elementary constants are checked in `check_cutoff.py`. Evenness, the bound ||psi2||<32 and (1+s)^12 exp(-s/2)<=24^12 give |chi(x)R_±(x)|<=2^90 exp(-(x-L)_+), |chi(x)|<=2^90(1+x)^-2 exp(-(x-L)_+). For x<=L+1000 the cutoff defect is <=exp(-10000), and beyond this use the exponential spectral tail. The raw Einstein light rational function is bounded by 2^10, the raw squared-cubic one by 2^30(1+x)^8. For the Einstein term first split x<=L/2 and use the stronger cutoff leakage exp(-(200+L/20)^2), then integrate (1+x)^-2 for x>=L/2. For the cubic term integrate (1+x)^6 and divide by L^7. Both light errors are <2^-500 by the displayed constants.

For angular cutoff transfer write y=(nu+2)/L. For J<=64 the coefficient recurrence multiplied by (k!)^2 gives |(k!)^2 t^k g[k,3]((nu±i)^2)|<=(1+y^2)^k, because the lower-order recurrence majorants sum to <.01. With x<=(64.5)^2 the sum and dressing are bounded by (2+y^2)exp(130(1+y)). For hat-b<=40 use Cauchy at |z|=1600/L^2. The Legendre Laplace integral has base (1+sqrt(z)e^(i theta))(1+sqrt(z)e^(-i theta))/(1-z); the absolute log is <=86/L (elementary log series). Thus t^k|g[k,3]((nu±i)^2)|<=2 exp(86y)*1600^(-k). Since x/1600<=1/4, the geometric sum with the quadratic dressing is bounded by 4. A common envelope for the dressed ladder action is 16(1+y)^2 exp(140(1+y)). This applies to physical finite ladders and is independent of J and L in the stated regions. The cutoff perturbation is bounded by integrating this envelope against (2/L)|chi R_±|(1-theta). Interior <=2^102 exp(-9440); exterior <=2^105 exp(-440), hence total <2^-500. The ideal certified lower bounds survive with room to spare. This estimate is NOT uniform on the missing unbounded-impact region.

## 5. Relative shifted-density defect and realized hybrid floor

For the boundary value write q(x-i)=W(x)/2+i h(x), where h=(1/4) PV integral W(v)csch(pi(x-v)/2)dv. Put G=Fh. The cutoff shifted defect is r=(Re theta(x-i)-1)f(x)-2 Im theta(x-i)G(x). Since W is even, h odd and theta even-real, r(0)=0 exactly.

For j<=3 the polynomial psi2^(j) is divisible by (1-p)^(5-j), with quotient coefficient l1 norm <2^40. On the real axis |(log F)^(j)|<100, j=1,2,3: use the pole-distance bounds for the five quadratic factors, the log-cosh derivatives and the polygamma sum/integral bounds at a=5/4. Thus |F^(j)|/F and F|(1/F)^(j)| are <=2^21. Away from zero, the derivatives of W up to order 3 are bounded by 2^65*(1+L-|v|)/(L F(v)) in its support, and also by 2^65*(1+|v|)/(L F(v)). Its first derivative has a jump of size at most 2^66/(L F(0)), while the regular second derivative is continuous at 0; the high-order endpoint zeros remove boundary jumps at ±L.

Use the difference form h(x)=(1/4) integral_0^infty csch(pi*s/2)[W(x-s)-W(x+s)]ds. For derivatives j<=2 use the regular derivatives W^(j), plus the delta term in the distributional W'' at 0. The ratio F(x)/F(v)<=2^25(1+|x-v|)^12. On s<=1 the mean-value estimate cancels the 1/s kernel, except for the W' jump when x<1, which contributes at most a constant times |log x|. On s>=1 use csch(pi*s/2)<=3 exp(-1.5s), and (1+s)^13 exp(-s/2)<=26^13<2^62. If x>L, support forces s>=x-L, retaining the exponential tail; for the local portion s<=1 one can lose an extra factor 3 to retain it. The W'' delta term is bounded by the same envelope for x>=1 and by a constant times 1/x at small x. Including F derivatives gives the conservative bounds |G^(j)(x)|<=2^181*(1+|L-x|)exp(-(x-L)_+)/L for x>=1, j<=2; for 0<x<=1 use |G|<=2^181/L, |G'|<=2^181(1+|log x|)/L and |G''|<=2^181(1+|log x|+1/x)/L. These follow directly from the split difference integral and the derivative bounds above (the local interval contribution has length <=1, the outer exponential integral <=1). G is continuous at 0 and G(0)=0.

Let C=L+2000. The erf cutoff and its first two derivatives on x-i obey the bounds 8 exp(-(C-x)_+^2/200) for theta-1, theta', theta''. On 0<=x<=1 the imaginary part has the extra factor x by even-real symmetry and the derivative estimate. The product rule therefore bounds |r| and |r''| on x>=1 by 2^200*(1+|L-x|)exp(-(x-L)_+-(C-x)_+^2/200)/L, and on 0<x<=1 by 2^200*(1+|log x|)exp(-(C-1)^2/200)/L. The 1/x singularity cancels; r'(0+) exists with the same endpoint constant without log. Integration below L, between L and L+1000 and above L+1000 gives

L (||r||_1+|r'(0+)|+||r''||_1) <= 2^212 (exp(-5000)+1002 exp(-1000)) < 2^-1000.

For real beta, let K0=(cosh beta-1)^(-1/2), A=||r||_1, B=|r'(0+)|+||r''||_1. Twice integrating the cosine transform by parts gives min(A,B/t^2). Its integral over positive t vanishes since r(0)=0. The Mehler subtraction lemma (TAIL_CERT R.2) gives the relative bound

|integral r(nu)P_{-1/2+i nu}(cosh beta)dnu| <= (sqrt(2)/pi) K0 [B/beta+8B log(2)/beta^2+A eps1(beta)],

eps1=(sinh(beta/2)-beta/2)/(2(cosh beta-cosh(beta/2))). For beta>=12, beta*eps1<=1. The retained LRINF_TIGHT section 1.3 region-c lower bound gives M1>=(sqrt(2)/pi) K0/(10 beta) there: evaluate the subtracted decreasing terms at beta=12, using e^6>400, sqrt(2)<3/2, log 2<1. For beta in [7/4,12], `check_M1.py` independently certifies M1>9e-5 by lower Riemann sums. It uses t=beta(1-v^2), the monotone numerator (t^2 cosh(t/2)-2(cosh t-1))/t^4 (all coefficients after its constant 1/24 are negative), positive increasing denominator sinhc(t/2)^2, and the weight 2 sqrt(beta)/sqrt(sinh(beta(1-v^2/2))*sinhc(beta*v^2/2)). Exponential endpoint values are enclosed by adjacent correctly-rounded Decimal values; other operations are rational/integer. The integral's factor sqrt(2)/pi is bounded below by 2/5.

T0_CERT's c=.199 sign-pattern proof has spare shift >=6e-7 in Part F, curvature margin >=4.58/L^2 in Part K, and >=.81*S_min*M1/L^2 in Part L (S_min>200). In each pattern the realized change is bounded by 2(1+.199)/L times the shifted residual transform. For Parts F,K and the compact portion of Part L, use |P|<=1 and the norm bound; in the tail use the Mehler relative estimate and the lower M1. All losses are far smaller than the retained margins. Hence the realized hybrid floor H_hyb>=.199 N_s (using the realized k=2 hybrid component) follows for every real gap >=1000, every spin and all impact parameters, subject to the T0/Mehler and angular normalization inputs. This is NOT a bound for H_ang-H_hyb, and cannot establish full finite-mass angular positivity in the uncovered region.

## 6. Reproduction and limitations

From this directory, Python 3 with SymPy and NumPy (object integer arithmetic):

```
python check_uniform.py 32
cp results.json results_J32.json
python check_physical.py 48 33
python check_physical.py 64 49
python impact_tensor.py --X 400 --K 110 --mq 14 --mv 40 --Jmin 65 --r 1/10
python check_cutoff.py
python check_M1.py
python check_residual.py
```

The exact JSON outputs of these runs are included. `check_fast.py`/`fast_results.json` contain diagnostic failures in the enlarged nonphysical mass rectangle, not a physical counterexample; the physical transform addresses them. Earlier smaller bounded-impact certificates are included as diagnostics. The main proofs use signed exact polynomials, not floating-point sign samples; decimal outputs are summaries only. The analytic factorization, cutoff and residual estimates have elementary arithmetic checks, not an independent referee or end-to-end recomputation of every analytic inequality. The cutoff is deliberately changed from Lambda+100 to Lambda+2000 to get uniform estimates, so previous full-action certificates for the old family cannot automatically be imported. All-gap full finite-mass positivity still requires a uniform large-spin, unbounded-impact estimate (J>=65, hat-b>40), or a rigorous alternative construction. The original gap-1000 certificate remains separate. No new PRL acceptance or novelty claim is made here.
