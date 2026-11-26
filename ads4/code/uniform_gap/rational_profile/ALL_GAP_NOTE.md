# All-gap positivity by a rational profile, and an exact admissible NF02 realization

## Status and scope

There are two distinct constructions. (A) The original sharp NF02 target has a new exact spectral realization, with no artificial Gaussian cutoff. Its all-gap leading-Regge/hybrid floor follows from the existing ideal T0 theorem. The full finite-mass sign for that target is **not** proved on its entire remaining spin/impact region. New exact fixed-spin certificates extend the earlier J<=64 result to J<=90, and two exact tensor certificates extend the all-spin bounded-impact region from bhat<=40 to bhat<=60. (B) A different, explicitly constructed k=2-only functional has **analytic strictly positive finite-mass angular action on every heavy state at every real gap >=1000 (the strict-positive modification in §6)**, assuming the same angular-action identity as in Theorem FP. It supplies a weaker, but gap-decaying, quantitative Einstein bound, of order sqrt(log gap)/gap^2, under the retained CFT sum-rule inputs. It is not the conjectured sharp sqrt(log gap)/gap^4 bound. The alternative construction is not a numerical extrapolation of the gap-1000 certificate. None of the GitHub files was edited.

Source conventions are `anshM123/MAJOR`, branch `claude/eloquent-cori-nftfjy`: [sec_fullpos.tex](https://github.com/anshM123/MAJOR/blob/claude/eloquent-cori-nftfjy/ads4/prl/sec_fullpos.tex) (fp-target, fp-single, fp-V and fp-LR); [K_B2A4.md](https://github.com/anshM123/MAJOR/blob/claude/eloquent-cori-nftfjy/ads4/K_B2A4.md) (light coefficients); [LEMMAS_CLOSE.md](https://github.com/anshM123/MAJOR/blob/claude/eloquent-cori-nftfjy/ads4/LEMMAS_CLOSE.md), [EESTEP_CLOSE.md](https://github.com/anshM123/MAJOR/blob/claude/eloquent-cori-nftfjy/ads4/EESTEP_CLOSE.md), and [supplement.tex](https://github.com/anshM123/MAJOR/blob/claude/eloquent-cori-nftfjy/ads4/prl/supplement.tex) (K2nu, LRC-a, LRC*, light side, interchange and operator assumptions); [T0_CERT.md](https://github.com/anshM123/MAJOR/blob/claude/eloquent-cori-nftfjy/ads4/T0_CERT.md). The latest angular ledger is ANGULAR_ROUTE.md §§127–129 (blob c2eab93141bed454866e28cb4a949371d9b2eeff): the old fixed-gap sliver was independently reproduced and refereed. The present results use that work's angular formula as an input, not its numerical sign check for the alternative profile.

## 1. A faster-decaying exact inverse of the de-shift equation

Let Lambda^2=(Delta_g-2)(Delta_g-1). Use the repository's R_-(z) and R_+(z)=R_-(-z):

```
R_-(z)=-2(2z-i)(2z+i)^2(2z+3i)(2z+5i)(2z+11i) /
              [z(2z-7i)(2z-3i)(2z+7i)].
```

Let P(z)=(z^2+1/4)^2(z^2+9/4)(z^2+25/4)(z^2+81/4) and

```
Fstar(z) = P(z) cosh(pi*z/2) Gamma(5/4+i*z/2) Gamma(5/4-i*z/2)
Phi(z) = Fstar(z)/cosh(pi*z)
B(z) = Phi(z+i)/R_-(z) = Phi(z-i)/R_+(z).
```

The modification of the preceding uniform-gap note is the division by cosh(pi*z). Its two shifted values at z+i and z-i coincide, because it is 2i-periodic. Gamma recurrence and exact rational algebra give Fstar(z+i)/Fstar(z-i)=R_-(z)/R_+(z); the Gamma/cosh ratio in this computation is (z+3i/2)/(z-3i/2). The identity

```
P(z+i)/R_-(z) = -z(2z-7i)^2(2z-3i)^2(2z+3i)(2z+7i)^2/2048
```

proves the pole assertions below. B is real-even on the real axis. Phi>0 there. B has a double zero at 0; its only possible poles in the closed strip |Im z|<=1 are simple poles at ±i/2. Indeed the shifted Gamma factors have positive real parts in this strip, and the shifted cosh(pi*(z+i)/2) is i*sinh(pi*z/2). By Stirling on bounded horizontal strips, uniformly for |Re z| large,

```
Phi(x) asymptotic to C |x|^(23/2) exp(-pi|x|),
B(x+i y) = O((1+|x|)^(19/2) exp(-pi|x|)), |y|<=1.
```

For an even real target E on the line, set W(x)=E(x)/Phi(x), and define

```
q(z) = (1/4) integral_R W(v) sech(pi*(z-v)/2) dv,
chi(z)=B(z)q(z).
```

We use either E compact and Lipschitz (NF02), or E and E' both O((1+|x|) exp(-pi|x|)) (the rational family below). In the first case W is compact and Lipschitz; in the second W,W'=O((1+|x|)^(-21/2)) (the logarithmic derivative of Phi is bounded at infinity). q is holomorphic in |Im z|<1, even and real on R. For y↑1, q(x+iy)+q(x-iy) is convolution with (1/2) Re sech(pi*(x-v+i y)/2). Under u=sinh(pi*(x-v)/2) this is the mass-one kernel (1/pi)*cos(pi*y/2)/(u^2+cos(pi*y/2)^2) du. Thus its limit is W(x). Individual complex boundary limits exist: subtract W(x) in the odd imaginary singular kernel on |x-v|<1 and use local Lipschitz continuity; outside use uniform exponential decay. Consequently

```
chi(x-i)R_-(x-i)+chi(x+i)R_+(x+i)=E(x)
```

exactly, using B(x-i)=Phi(x)/R_-(x-i) and its reflected identity. For a compact target q(x+iy)=O(exp(-pi|x|/2)) outside its support enlarged by one, uniformly |y|<=1. In the second case q(x+iy)=O((1+|x|)^(-21/2)): on |x-v|<1 the real kernel has bounded mass and the imaginary kernel is bounded after the subtraction; outside, exponential convolution preserves polynomial decay. Hence:

* compact target: chi(x+iy)=O((1+|x|)^(19/2) exp(-3pi|x|/2));
* rational target: chi(x+iy)=O((1+|x|)^(-1) exp(-pi|x|)).

These constants can depend on the gap/profile. No bound uniform in gap is needed for existence and absolute convergence at each fixed gap. In the lower strip the possible simple pole of chi at -i/2 is canceled by the double zero of R_-. Its 1/z pole at 0 is canceled by the double zero of chi. The upper strip uses R_+.

### Position-space convergence, with retained kernel inputs

The real-axis majorant K2nu in LEMMAS_CLOSE §3.1 uses the heavy normalization N_com ~ |nu|^9 exp(-pi|nu|), hence |chi/N_com| is integrable (for rational target O(|nu|^-10)). The leading edge factors in §3.2 are proportional to chi R_-(z)/(E_7(z+i)s(z+i)) in the lower strip, and their reflections in the upper strip, where E_7=2^12*pi*Gamma((13/2-i nu)/2)^2 Gamma((13/2+i nu)/2)^2 and s(w)=32(4w^2+1)(4w^2+9)/((4w^2+25)(4w^2+81)). E_7 has no zeros and no poles in the needed strip and grows as |nu|^11 exp(-pi|nu|) on horizontal lines. The double zero of R_- at -i/2 leaves a simple zero after division by s(z+i), canceling the simple pole of chi. At zero the aforementioned double zero of chi suffices. Thus these edge weights are holomorphic in the requisite open half strips and decay at least as O(|nu|^-10) on fixed horizontal lines for the rational target, faster for compact target. The reflected statement is identical. A shift of depth 3/4 suffices to improve the raw +1/2 Regge edge exponent to -1/4 (no boundary limiting subtlety); the full depth 1 is also available for the moment integrals. With the retained LRC* description of the exact edge as a finite sum of shifted scalar waves with bounded holomorphic hypergeometric multipliers and the retained K1', K2, K3 and LRC-a assumptions, the convergence argument works for these new weights: the shifted contour edge integrals are absolutely convergent, the corner exponent is negative at J_0<=1, and the position-space smeared sum rule can be interchanged with the OPE. This is a check of the **weight-side sufficient hypotheses** in the repository's argument, not an independent proof of K1', K2, K3, LRC* or the physical angular identity. In particular one cannot assume that an arbitrary exponentially decaying spectral weight is admissible; the polynomial rate and the canceled edge pole matter here.

No cutoff is used. All finite polynomial spectral moments can be shifted exactly, by the half-strip holomorphy and exponential decay. For the compact NF02 target E(x)=psi2(|x|/Lambda)1_{|x|<=Lambda}, the ideal shifted density is now exact. This removes the cutoff deformation and coherent spectral error from the new NF02 family. It is a different realization from the old fitted chi; the retained fp-single identity makes its finite-mass action exactly the ideal one.

## 2. Positive rational profiles at every gap and spin

Fix any real Delta_g>=1000 and M=floor(Delta_g-2), so M>=998. For j=1,...,m with m=M+6 take distinct rational rho_j=j/(100 m^2), and put

```
Rprof(z)=prod_j (1-rho_j*z)^(-1),
F(z)=z^3 (1-z)^M Rprof(z).
```

The Taylor coefficients of Rprof are nonnegative; Rprof(0)=1, and 1<=Rprof(z)<=100/99 on 0<=z<=1 because sum rho_j<=1/100. This F is not NF02. It uses only k=2; set the k=3 smearing identically zero.

### 2.1 Exact realization by the new inverse

Let x=(1+z)/(1-z). The desired spherical profile f(x)=(1-z)^3 F(z) is

```
f(x) = 2^(M+3)*(x-1)^3 / prod_j [(1-rho_j)*x+(1+rho_j)].
```

It is a proper rational function with distinct simple poles at -A_j=-(1+rho_j)/(1-rho_j)<-1. Write its exact partial fractions f(x)=sum_j C_j/(x+A_j). The classical Mehler–Fock identity [DLMF 14.20.14](https://dlmf.nist.gov/14.20.E14) is

```
1/(x+A) = pi integral_0^infty nu*tanh(pi*nu)/cosh(pi*nu) * P_{-1/2+i nu}(x)*P_{-1/2+i nu}(A) dnu.
```

It holds for x,A>1 and by continuity at x=1. Define the even real target

```
E(nu)=(Lambda/2)*pi*nu*tanh(pi*nu)/cosh(pi*nu) * sum_j C_j P_{-1/2+i nu}(A_j).
```

Then (2/Lambda) integral E P(x)=f(x). Mehler–Dirichlet at each fixed A_j>1 bounds P and its first frequency derivative uniformly on the real axis (differentiate the cosine under the finite integral). Thus E,E'=O((1+|nu|)exp(-pi|nu|)), and the rational target meets section 1. Though the partial-fraction coefficients may be very large, this is an existence proof with exact identities; it does not rely on evaluating ill-conditioned residues numerically.

The retained generating identity (fp-single/fp-V) gives exactly F(z)=(2/Lambda) integral E(nu)(1-z)^(-3) P_{-1/2+i nu}((1+z)/(1-z)) dnu. Coefficient comparison is valid by dominated holomorphy on a small disk around 0 or by finite polynomial moment shifts. The angular action in channel b is (1-r)^3 L_b[(1+D_b/T)F], with the definitions below, and the overall state prefactor is positive on the heavy range. This is the same angular normalization as Theorem FP, with the k=3 term absent.

### 2.2 Positive basis lemma (proof, all real twists)

Let b=0,4, n=J-b>=0, a=Delta-J-2>=M, N=n+a, s=2b+1, T=(J+a)(J+a+1), r=J(J+1)/T. Define L_b(z^k)=U_k with

```
U_k = (n)_down k * (n+s)_up k / [(N)_down k * (N+s)_up k]  (0<=k<=n),
U_k=0 (k>n).
```

For any integer k>=0 and integer m0 in [0,M+1], L_b[z^k(1-z)^m0]>=0; it is strictly positive if k<=n.

Proof: k>n gives zero. Divide by positive U_k otherwise. The shifted ladder has n'=n-k, N'=N-k, s'=s+2k, same a. If n'=0 the normalized value is 1. For n'>=1, N'>=M+1>=m0. The generalized hypergeometric probabilities

```
p_j = binom(m0,j)*(n')_down j*(a)_down(m0-j)/(N')_down m0, j=0,...,m0
```

are nonnegative since a>=M and m0<=M+1 (even the last factor of (a)_down(m0-j) is nonnegative), and sum to 1 by Vandermonde's falling-factorial identity, valid for real a by polynomial continuation. Their factorial moments are E[(X)_down ell]=(m0)_down ell (n')_down ell/(N')_down ell; this follows by removing ell chosen successes and applying the same Vandermonde identity. Independently take B0~Beta(n'+s',a). Its moments E[B0^ell]=(n'+s')_up ell/(N'+s')_up ell. Thus the normalized shifted-ladder value is

```
sum_{ell=0}^{m0} (-1)^ell binom(m0,ell) (n')_down ell (n'+s')_up ell / [(N')_down ell (N'+s')_up ell]
   = E[(1-B0)^X] > 0.
```

Here X has the above probabilities. Since B0 in (0,1) almost surely, the expectation is strictly positive. This proves the lemma including noninteger masses; it does not assume that U_k is a positive moment sequence (it is not).

### 2.3 Dressing preserves positivity for this profile

The first-moment operator from fp-V is

```
D_b = z(1-z)/2 * d_z^2 + [(7/2-2b)z+1/2]*d_z +9/2-b(b-1).
```

For P_{k,M}=z^k(1-z)^M, the exact product rule gives

```
D_b P_{k,M} = k(k-1)/2 P_{k-1,M+1} + k/2 P_{k-1,M}
              + C_{k,M,b} P_{k,M} + M(M-9+4b)/2 P_{k+1,M-1},
C_{k,M,b}=k(7/2-2b-M)+9/2-b(b-1)-M/2.
```

Only k>=3 is needed. All the off-diagonal coefficients are nonnegative. If k<=n<=J, the potentially negative diagonal satisfies

```
-C_{k,M,b} <= JM+(9/2)J+M/2+15/2 <= Ja+(9/2)(J+a)+15/2 < T/2.
```

Indeed Ja <= (J+a)^2/4 and (9/2)(J+a)+15/2 < (J+a)^2/4 for J+a>=998. If k=n+1 the diagonal ladder value is zero; if k>n+1 all four terms vanish. Apply the positive basis lemma to every term of (1+D_b/T)P_{k,M}. This proves L_b[(1+D_b/T)P_{k,M}]>=0. If n>=2 and k=3, the P_{2,M} term gives a strict positive contribution. Expand Rprof(z)=sum_{j>=0} d_j z^j with d_j>=0,d_0=1; L_b(1+D_b/T) sees only finitely many coefficients, so termwise summation is exact. Multiplication by (1-r)^3>0 proves:

**Theorem RP (conditional on the exact angular-action identity).** For every real Delta_g>=1000, every heavy primary with Delta-J>=Delta_g, and both physical channels, the rational-profile functional has nonnegative action. It is strictly positive for n=J-b>=2, and zero for n=0,1. The angular matrix is PSD by the positive congruence in fp-target. This theorem covers all J, all twists, and all impact parameters; no large-spin or large-gap limit and no numerical sign samples are used. It applies to the new profile, not to NF02.

## 3. Rigidity consequence and rate

The light-side calculation here retains the repository's per-frequency Einstein and W^3-squared light identities, their parity identification, the vanishing of the linear term, light double-trace zeros and the order of limits at fixed gap. Denote the shifted density by E(nu) and the Einstein light numerator by E_light. The repository's rational raw k=2 W^3-squared coefficient is R_-(nu)Q(nu+i)+R_+(nu)Q(nu-i), where

```
Q(nu)=nu^6-(245/4)nu^4+(9619/16)nu^2-51575/64.
```

The first four generating polynomials g_{k,3}(nu^2) satisfy the **exact** identity Q=-100 g_{2,3}-36 g_{3,3}. Therefore, writing F_j=[z^j]F(z), the normalized coefficient in K_B2A4 is

```
K = -(2/Lambda^7) integral_0^infty E Q = (100 F_2+36 F_3)/Lambda^6.
```

Here F_2=0,F_3=1, so K=36/Lambda^6>0 exactly. No k=3 contribution is used. For the Einstein pole D(nu)=1/(nu^2+9/4) (its imaginary poles lie outside the shifts), contour shift gives E_light=16 Lambda integral E D. The resolvent identity

```
integral_1^infty Q_1(x) P_{-1/2+i nu}(x) dx = 1/(nu^2+9/4),
Q_1(x)=(x/2)*log((x+1)/(x-1))-1>0
```

follows by subtracting the two Legendre differential equations in Sturm–Liouville form. At x=1, P(1)=1, (x^2-1)Q_1'(x)→-1 and Q_1(x)(x^2-1)P'(x)→0, while the upper boundary vanishes by their large-x decays. Thus E_light=8 Lambda^2 I with I=integral_1^infty f(x)Q_1(x)dx. All integrals converge absolutely (E has exponential tail and f decays as x^(-M-3)).

At leading order in 1/c_T, the retained sum-rule inequality K tau^2 Lambda^8/57600 <= E_light now gives tau^2 <=12800 I. Substituting x=(1+z)/(1-z) and dropping the negative term of Q_1,

```
I <= (100/99)*integral_0^1 z^3 (1-z)^M (1+z) (-log z) dz
  = (100/99)*[B(4,M+1)(H_{M+4}-H_3)+B(5,M+1)(H_{M+5}-H_4)]
  <= (100/99)*(1+4/(M+5))*6/(M+1)^4 * log((M+4)/3).
```

Here B(5,M+1)/B(4,M+1)=4/(M+5), and the harmonic sum is bounded by the integral of 1/t from 3 to M+4. The rational constant check for M>=998 gives 12800*(100/99)*(1+4/(M+5))*6 < 280^2. Hence

```
tau <= 280 sqrt(log((M+4)/3))/(M+1)^2
    <= 280 sqrt(log((Delta_g+2)/3))/(Delta_g-2)^2,
M=floor(Delta_g-2), every real Delta_g>=1000.
```

**CFT scope:** this bound is conditional on the retained angular/position-space and light-side sum-rule hypotheses, including the treatment of additional light exchanges. If there are no other light single-trace exchanges besides T, no extra light sign condition is needed. If the exact angular identity extends to long light spin-0 and spin-2 primaries (the allowed low-spin exchanges in the manuscript), the new profile actually annihilates spin 0 (the relevant low coefficients vanish) and has positive spin-2 long action: n=2, V_b(F)=9 U_2/(2T), with T>2 for Delta>3, and the prefactor times (1-r)^3 is 2pi G Lambda^4 T/9. The stress tensor short block is instead handled by the retained light formula. The extension to exceptional short/light representations has not been independently derived here; do not silently inherit the manuscript's NF02 light-matter argument for a different profile. The all-gap *heavy* positivity theorem does not need this further light-side issue resolved.

This is a universal decaying bound of order sqrt(log Delta_g)/Delta_g^2. It is weaker in power than the desired sharp NF02 gap-dependent bound of order sqrt(log Delta_g)/Delta_g^4; known parametric causality results are not thereby improved. At a gap of 1000 the older fixed-gap bound is much stronger. It is best viewed as a second explicit, analytic positivity family and a new mechanism for the all-gap program, not a replacement for the original numerical headline. The proof uses no assumption of monotonicity in gap.

## 4. Original NF02 continuation and remaining target

For the exact NF02 realization of §1, all ideal angular certificates are now directly applicable, without the former cutoff error budget. `check_physical.py 90 65` passed all 39 new channel-spin cases, with exact rational Bernstein ledgers, the minimum normalized lower bound >2.07e-11. Together with the earlier 94 cases, this proves all physical channels through J=90, continuously over every mass and every real gap >=1000. The previous all-spin bhat<=40 certificate also applies. The new `impact_tensor.py --X 900 --K 150 --mq 22 --mv 60 --Jmin 91 --r 1/10` run passed both channels (two x-leaves each). At X=900 and J>=91 it verifies bhat<=60 with normalized lower bounds 0.0041216454207107876 and 0.004121645419125774. The polynomial truncation error bounds are 3.521e-17 and 1.822e-12, and the ordinary ladder tail is <2.302e-81. The proof of the tensor recurrence and truncation is unchanged from the preceding uniform-gap note; in particular its positive pole check and x*v<=1 check were executed, and terms beyond physical spin termination have an exact zero. Together with J<=90, this covers every spin at bhat<=60. The sharp NF02 finite-mass sign remains open on J>=91, bhat>60. The hybrid T0 all-gap positivity follows exactly with its retained inputs, no residual. The K polynomial and conditional sharp rigidity formula from the preceding uniform-gap note hold for this exact weight, using the same light identities and with no cutoff allowance. The unresolved full finite-mass NF02 sign on its remaining unbounded-spin/impact domain is still needed for the sharp all-gap theorem.

## 5. Reproduction

`python check_physical.py 90 65` regenerates the new fixed-spin ledger. `python impact_tensor.py --X 900 --K 150 --mq 22 --mv 60 --Jmin 91 --r 1/10` regenerates the new two-channel impact ledger. `python check_rational.py` checks the factorization and edge zero, both raw light identities and Q=-100 g2-36 g3, the exact product-rule identity, 1800 rational hypergeometric distribution/factorial-moment tests (including fractional twists), 128 exact dressed-basis regression tests, the spin inventory and positive ledger endpoints, and the rational 280 constant. These are regression checks of the written analytic proof; the proof of Theorem RP is not an extrapolation of the tests. There has been no independent human or repository referee review of this new family. The new weight is different even at gap 1000, and the full CFT operator identities and sum-rule hypotheses have not been rederived. No claim of PRL acceptance or novelty priority is made.

## 6. Final strengthened form: strict positivity including the lowest spins

The following small modification is the recommended version of Theorem RP. Set epsilon=(M+1)^(-8) and replace F in §2 by

```
F_strict(z) = (z^3 + epsilon*z) (1-z)^M Rprof(z).
```

Its spherical profile still has the same distinct simple poles at x=-A_j; the numerator in §2.1 is replaced by `2^(M+3)*[(x-1)^3+epsilon*(x-1)*(x+1)^2]`. Thus all the weight realization and convergence statements of §1–2.1 are unchanged. The positive basis and dressing proof of §2.2–2.3 holds also for k=1 and k=2: their first off-diagonal factor k(k-1)/2 is nonnegative, and the diagonal bound uses only k<=J. Since the new k=1 term has a strictly positive k/2 contribution to P_{0,M}, the dressed action is strictly positive for **every** n>=0, including n=0,1. No strict lower bound independent of the state is claimed.

The light normalization remains positive with the same final numerical constant. Let `A(z)=(1-z)^M Rprof(z)=1+A1*z+A2*z^2+...`. Then A1=-M+sum rho_j >= -M, and A2=M(M-1)/2-M*sum rho_j+[z^2]Rprof>=0. For F_strict, F_0=0, F_1=epsilon, F_2=epsilon*A1, F_3=1+epsilon*A2. Hence

```
K Lambda^6 =100 F_2+36 F_3 >=36-100 epsilon*M > 35.99.
```

The extra term in I of §3 is bounded above by

```
epsilon*(100/99)*[B(2,M+1)(H_{M+2}-H_1)+B(3,M+1)(H_{M+3}-H_2)]
 <= epsilon*(100/99)*(1+2/(M+3))*2 log((M+4)/3)/(M+1)^2.
```

The last inequality uses log(M+2)<=2 log((M+4)/3), valid for M>=2. Together with §3 this gives

```
I <= C_I * log((M+4)/3)/(M+1)^4,
C_I=(100/99)*[6(1+4/(M+5))+2(1+2/(M+3))/(M+1)^6].
```

Since `tau^2 <= 460800 I/(K Lambda^6)`, the exact rational test at M=998 gives `460800 C_I/35.99 < 280^2`; all bounding factors are non-increasing in M. The same bound `tau <= 280 sqrt(log((M+4)/3))/(M+1)^2` therefore holds in the strict-positive version.

### Low-spin light matter check, where the retained angular identity extends

The strict-positive modification also has a favorable sign on the two allowed types of long light exchange, without requiring a>=M. This uses only their low Taylor coefficients, not extrapolation from heavy states.

* J=0: T L_b[(1+D_b/T)F_strict]=epsilon/2 exactly, since U_0=1 and U_k=0 for k>=1. The physical prefactor times (1-r)^3 is 2pi G Lambda^4 T/9. The apparent T singularity cancels, and this action coefficient is positive even if the scalar's T is negative. Exceptional masses where the representation formula is singular require the physical continuous extension.
* J=2, b=0: for a long unitary spin-2 primary, Delta>3 and T>2. U_1=6/T, U_2=24/[T(T-2)]>0. By fp-single the quantity `T L_b[(1+D_b/T)F_strict]/U_2` equals

```
epsilon*(T-2)*(13*T/48+2) + F_2*(3*T/2+19/2) + (9/2)*F_3
 >= 9/2 + epsilon*((T-2)*(13*T/48+2)-M*(3*T/2+19/2)).
```

The bracket is bounded below by `T^2/4-3*M*T/2-19*M/2-4 >= -9*M^2/4-19*M/2-4 >= -16*(M+1)^2`. Therefore the quantity is > 9/2-16/(M+1)^6 > 0. The short stress tensor at Delta=3 is not included in this long-block argument; its contribution is the separately retained Einstein/W^3 light formula.

This checks the algebraic sign of the long low-spin expressions. It does not independently establish the heavy-to-light extension of the angular formula for all physical representations. The conservative CFT theorem therefore either assumes that extension and the retained light-side identities, or restricts to the case of no additional light single traces. Strict positivity on the heavy spectrum requires neither restriction.
