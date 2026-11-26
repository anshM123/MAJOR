# All-gap mixed-profile positivity: an exact two-channel construction

## Scope

This proves strict positivity of the full **finite-mass angular action** for a new two-channel spectral family at every real gap \(\Delta_g\ge1000\), every physical mass, every physical spin and both \(b=0,4\), given the retained angular-action identity `fp-target` / `fp-V` in the project's `ads4/prl/sec_fullpos.tex`. It improves the earlier rational-profile conditional CFT rigidity rate from \(\Delta_g^{-2}\sqrt{\log\Delta_g}\) to \(\Delta_g^{-5/2}\sqrt{\log\Delta_g}\). It is NOT full-spectrum positivity of the original NF02 functional, and does not prove its sharper \(\Delta_g^{-4}\sqrt{\log\Delta_g}\) bound. The new profile is a deliberate change of functional, even at gap 1000. No repository was edited. This note is not an independent rederivation of the physical operator identities, light-side formulas, or CFT interchange theorem.

Retained sources: [finite-mass angular formula](https://github.com/anshM123/MAJOR/blob/claude/eloquent-cori-nftfjy/ads4/prl/sec_fullpos.tex), [light normalization](https://github.com/anshM123/MAJOR/blob/claude/eloquent-cori-nftfjy/ads4/K_B2A4.md), [k=3 kernel reading](https://github.com/anshM123/MAJOR/blob/claude/eloquent-cori-nftfjy/ads4/S1_S4_FIXES.md), and the preceding conversation's `ALL_GAP_NOTE.md`, §§1–3 and 6 (exact de-shift realization and positive-basis lemma). All sign inequalities below are analytic; `check_final.py` independently checks the rational constants, product rule, direct monomial reduction, and 960 finite regression cases.

## 1. Profiles and the all-gap theorem

Let \(a_g=\Delta_g-2\), \(M=\lfloor a_g\rfloor\ge998\), \(\Lambda^2=(\Delta_g-2)(\Delta_g-1)\), \(\epsilon=(M+1)^{-8}\), and \(A=7M^2/10\).

Choose any \(R(z)=\sum_{j\ge0}r_j z^j\) analytic near the closed unit disk with nonnegative coefficients, \(r_0=1\), and \(R(1)\le100/99\). A spectrally realizable choice is given in §3.

\[
 F=(z^2+\epsilon z)(1-z)^M R(z),\qquad G=A z^2(1-z)^M R(z).
\]

For a heavy state \(a=\Delta-J-2\ge a_g\), let \(n=J-b\ge0\), \(s=2b+1\), \(N=n+a\), \(T=(J+a)(J+a+1)\), \(w=1-J(J+1)/T>0\). The ladder operator acts on monomials by \(\mathcal L_b(z^k)=U_k\), with
\[
 U_k=\frac{(n)_{\downarrow k}(n+s)_{\uparrow k}}{(N)_{\downarrow k}(N+s)_{\uparrow k}}\quad (0\le k\le n),\qquad U_k=0\quad(k>n).
\]
Use the retained first-moment operator
\[
 \mathcal D_b=\frac{z(1-z)}2\partial_z^2+
  \left[(7/2-2b)z+1/2\right]\partial_z+9/2-b(b-1).
\]
The retained angular formula is
\[
 H_b=w^3\mathcal L_b[(1+\mathcal D_b/T)F+\sigma_bG/T],\qquad
 \sigma_0=+1,\quad\sigma_4=-1.
\]

**Theorem:** \(H_b>0\) for every heavy state at every real \(\Delta_g\ge1000\), both channels. In fact for \(n\ge1\),
\[
 H_b\ge\frac{w^3}{25}\mathcal L_b(F)>0.
\]

### Positive-basis lemma (real twists, not just integer mass)

For integer \(k\ge0, 0\le d\le M+1\),
\(\mathcal L_b[z^k(1-z)^d]\ge0\), strictly if \(k\le n\).
For \(k>n\) the value is zero. Otherwise divide by \(U_k>0\) and shift \(n'=n-k,\ N'=N-k,\ s'=s+2k\). If \(n'=0\) the ratio is one. If \(n'\ge1\), define the finite probabilities
\[
 p_j=\binom dj\frac{(n')_{\downarrow j}(a)_{\downarrow(d-j)}}{(N')_{\downarrow d}}.
\]
They are nonnegative because \(a\ge M\), \(d\le M+1\), and \(N'\ge M+1\); they sum to one by Vandermonde's falling-factorial identity (a polynomial identity, valid for real \(a\)). Their factorial moments are \((d)_{\downarrow l}(n')_{\downarrow l}/(N')_{\downarrow l}\). Let \(X\) have these probabilities, independently \(B\sim{\rm Beta}(n'+s',a)\). The normalized ladder value is precisely \(\mathbb E(1-B)^X>0\), by binomial expansion. This proves the lemma without an integer-twist restriction.

### Dressing and helicity domination

For \(P_{k,M}=z^k(1-z)^M\), the exact product rule gives
\[
 \mathcal D_bP_{k,M}=\tfrac{k(k-1)}2 P_{k-1,M+1}
  +\tfrac k2P_{k-1,M}+C_{k,M,b}P_{k,M}
  +\tfrac{M(M-9+4b)}2 P_{k+1,M-1},
\]
where \(C_{k,M,b}=k(7/2-2b-M)+9/2-b(b-1)-M/2\).
For \(k\ge1\), all off-diagonal coefficients are nonnegative and their ladder values are nonnegative by the lemma. If \(k\le n\), write \(x=J+a\ge998\); then
\[
 -C_{k,M,b}\le JM+(9/2)J+M/2+15/2
 \le x^2/4+(9/2)x+15/2 < (1/4+1/100)x^2.
\]
The last inequality follows at 998 by exact rational arithmetic, and improves thereafter. Since \(T=x(x+1)\), the diagonal multiplier is at least \(37/50\). If \(k>n\) its ladder value vanishes and the same lower bound holds trivially. Applying this termwise to the nonnegative coefficients of \(R\) (the finite ladder and first derivative see only finitely many of them),
\[
 \mathcal L_b[(1+\mathcal D_b/T)F]\ge (37/50)\mathcal L_b(F).
\]
Also \(0\le\mathcal L_b(G)\le A\mathcal L_b(F)\) by the positive-basis lemma; \(A/T\le7/10\) since \(T\ge a^2\ge M^2\). Subtraction proves the stated \(1/25\) margin. For \(n=0\), \(F(0)=G(0)=0\) but \(F'(0)=\epsilon>0\), so \(H_b=w^3\epsilon/(2T)>0\) exactly. Thus there are no omitted heavy states in this new profile's positivity theorem.

## 2. Light-side bound (conditional physical implication)

Under the retained exact per-frequency light identities, the cubic coefficient is
\[
 K\Lambda^6=100F_2+36F_3+8G_2,
\]
where \(F_j=[z^j]F\) and \(G_j=[z^j]G\). For \(k=2\), the previous note gives the shifted light polynomial \(Q=-100g_{2,3}-36g_{3,3}\); for \(k=3\), \(4g_{2,4}=\nu^4-(35/2)\nu^2+569/16\) gives the \(8G_2\) contribution in the normalization of `fp-V`. Here \(F=(2/\Lambda)\int E_2F_{3,\nu}\,d\nu\), \(G=\int E_3F_{4,\nu}\,d\nu\), with \(E_3(\nu)=\Lambda\psi_3(\nu/\Lambda)\). All integrals in this statement are on the positive real axis.

Put \(A_0(z)=(1-z)^M R(z)\). Then \(A_1=-M+r_1\ge-M\), \(A_2=\binom M2-Mr_1+r_2\ge0\) since \(\sum_{j\ge1}r_j\le1/99\). Hence
\(F_1=\epsilon,\ F_2=1+\epsilon A_1,\ F_3=A_1+\epsilon A_2,\ G_2=A\).
The exact lower bound is
\[
 K\Lambda^6\ge(28/5)M^2-36M+100(1-\epsilon M)>(11/2)M^2.
\]

The Einstein coefficient is \(E_{\rm light}=8\Lambda^2 I\), where
\[
 I=\int_1^\infty f_2(x)Q_1(x)\,dx,\qquad
 f_2(x)=(1-z)^3F(z),\quad z=(x-1)/(x+1),\quad
 Q_1(x)=(x/2)\log((x+1)/(x-1))-1.
\]
This is the retained Einstein-pole resolvent calculation (previous note §3). The cubic \(k=3\) Einstein value is zero. Changing variables and dropping the negative \(-2(1-z)\) term gives
\[
 I\le\frac{100}{99}\int_0^1(z^2+\epsilon z)(1-z)^M(1+z)(-\log z)\,dz.
\]
Use \(\int_0^1 z^j(1-z)^M(-\log z)dz=B(j+1,M+1)(H_{M+j+1}-H_j)\). Each harmonic difference for \(j=1,2,3\) is at most \(\log\Delta_g\), by \(\Delta_g\ge M+2\) and the integral bound \(\log((M+j+1)/j)\le\log(M+2)\). Consequently
\[
 I\le\frac{100}{99}\frac{\log\Delta_g}{(M+1)^3}
 \left[2(1+3/(M+4))+(1+2/(M+3))/(M+1)^7\right].
\]
At leading order in \(1/c_T\), the retained sum rule gives
\(K\tau^2\Lambda^8/57600\le E_{\rm light}\), hence \(\tau^2\le460800 I/(K\Lambda^6)\). Since \(a_g<M+1\), we may replace \((M+1)^{-3}M^{-2}\) by \((1+1/M)^2 a_g^{-5}\). The resulting squared coefficient is at most
\[
 \frac{460800\cdot2}{11}\frac{100}{99}
 \left[2(1+3/1002)+(1+2/1001)/999^7\right](999/998)^2
 <420^2.
\]
Every bounding factor decreases in \(M\), and the last inequality is checked exactly in `check_final.py`. Therefore, conditional on the physical sum-rule inputs,
\[
 \boxed{\sqrt{t_4^2+\alpha_T^2}\le
 \frac{420\sqrt{\log\Delta_g}}{(\Delta_g-2)^{5/2}}
 \quad(\Delta_g\ge1000).}
\]

**Extra light exchanges:** the new profile also has nonnegative algebraic long spin-0 and spin-2 action if the retained angular identity extends to them. For spin zero \(T\mathcal L[(1+\mathcal D_0/T)F+G/T]=\epsilon/2\), so the full common prefactor cancels the apparent T singularity. For long spin 2, \(b=0,T>2\), \(U_2=24/[T(T-2)]>0\), and the retained monomial action gives
\[
 \frac{T\mathcal L[(1+\mathcal D_0/T)F+G/T]}{U_2}
 = (T-2)(13T/48+2)F_1+(3T/2+19/2)F_2+(9/2)F_3+G_2.
\]
The first two terms are positive, the third is at least \(-9M/2\), and \(G_2=7M^2/10>9M/2\). This is an algebraic check, not a proof of the angular identity for exceptional light/short representations. A conservative CFT statement therefore assumes no additional light single traces, or that they have nonnegative action; the long-representation calculation is useful evidence toward the latter. The heavy positivity theorem itself does not need either light hypothesis.

## 3. Realization and convergence

For an explicit choice with controlled spectral tails, put \(D=M+6\) and
\(\rho_j(t)=(j+t)/(100(D+1)^2)\) for \(1\le j\le D\).
Take any nonnegative \(C^\infty\) probability bump \(\omega(t)\) compactly supported in \(1/4<t<3/4\), and define
\[
 R(z)=\int\omega(t)\prod_{j=1}^D(1-\rho_j(t)z)^{-1}\,dt.
\]
Every product has nonnegative Taylor coefficients and its value at 1 is \(\le100/99\), since \(\sum\rho_j<1/100\). Averaging preserves both properties. The poles are distinct for each \(t\) and uniformly separated on the compact support.

For each \(t\), both \(f_2=(1-z)^3 F\) and \(f_3=(1-z)^4 G\) are proper rational functions of \(x=(1+z)/(1-z)\), with simple poles \(-A_j(t)=-(1+\rho_j)/(1-\rho_j)<-1\). Indeed the number of denominator factors \(D=M+6\) is at least the sum of the numerator power of \(z\) and the power of \(1-z\) for both profiles. Hence the partial fractions have smooth coefficients \(C_{i,j}(t)\). DLMF 14.20.14 ([source](https://dlmf.nist.gov/14.20.E14)) gives
\[
 \frac1{x+A}=\pi\int_0^\infty
 \frac{\nu\tanh(\pi\nu)}{\cosh(\pi\nu)}P_{-1/2+i\nu}(x)P_{-1/2+i\nu}(A)\,d\nu.
\]
It gives \(f_i(x)=\int\widetilde E_i(\nu)P_{-1/2+i\nu}(x)\,d\nu\) with
\(\widetilde E_i=\pi\nu\tanh(\pi\nu)\cosh(\pi\nu)^{-1}
\sum_j\int\omega C_{i,j}P_\nu(A_j(t))dt\).
Use \(E_2=(\Lambda/2)\widetilde E_2,\ E_3=\widetilde E_3\). Since each \(A_j(t)\) is smooth strictly increasing with nonzero derivative, the signed measure under each \(t\)-integral pushes forward to a smooth compactly supported density in \(A>1\). Repeated integration by parts with the self-adjoint Legendre differential operator \(\partial_A((A^2-1)\partial_A)\), whose eigenvalue on \(P_\nu\) is \(-(\nu^2+1/4)\), gives arbitrary inverse powers of \(\nu\), using the Mehler bound on the compact A support. Thus for each fixed gap
\(\widetilde E_i(\nu)=O_{K,M}((1+\nu)^{-K}e^{-\pi\nu})\) for every K, and \(\widetilde E_i=O(\nu^2)\) at zero. Frequency derivatives can be handled by differentiating the same integration-by-parts identity and Mehler integral. In particular the k=2 density meets the exact fast-decaying de-shift construction in the preceding note §1. No approximation of the angular coefficients is used.

For the new noncompact k=3 spectrum, a sufficient physical convergence condition is the retained raw kernel modulus/edge bound and a polynomial-times-\(e^{\pi|\nu|}\) bound for the inverse heavy normalization with no real-axis singularities (or an integrable one at 0). The new rapid polynomial tail meets this sufficient condition. The retained k=3 kernel reading in `S1_S4_FIXES.md` §1 gives useful algebraic support: using the rest-frame monomial transforms in `q1_c3_regge_kernel.py`, factoring out \(K(c,\lambda)=2^{c+2-\lambda}\sqrt\pi\Gamma(1+c/2)\Gamma((c+3)/2-\lambda)\), I obtained the following rational factors, up to overall helicity sign, with \(\delta=1/2+i\nu\):
\[
 C_A=\frac{16i(\delta^2-16\delta+21)}{(\delta-7)(\delta+4)(\delta+6)},
 \qquad C_B=\frac{-2i(\delta-2)(\delta-1)(\delta+1)}{(\delta-5)(\delta+4)}.
\]
The respective reference factors are \(K(7-\delta,1-\delta)\) and \(K(4+\delta,\delta)\). The two rational factors have no zero or pole on the critical line; \(C_AC_B\to32\). For the first numerator, its imaginary part is \(-15\nu\), and its real value at \(\nu=0\) is \(53/4\). The product of the reference Gamma factors has magnitude of order \(\nu^{13}e^{-\pi\nu}\). The omitted odd angular average \(2B(1/2,1-\delta/2)\) is a zero-free gamma ratio of polynomial size. These are calculations **in the retained monomial/kernel conventions**; they do not independently match all the physical normalization conventions or prove the position-space edge bound. Consequently keep the k=3 convergence condition explicit in the physical CFT implication. The angular positivity proof is not contingent on an estimate of those large-frequency tails: its finite spectral moments are exactly those of the profiles, and all converge.

## 4. NF02 obstruction and status

It is tempting to sharpen the profile by taking the positive-basis exponent proportional to \(a_g^2\), to localize at \(z\sim a_g^{-2}\). This cannot be justified by the positive-basis lemma. A physical exact counterexample to its naive quadratic extension: \(a=998,b=4,n=1\) (J=5), \(m=998^2/4=249001\). Then
\[
 \mathcal L_b[(1-z)^m]=1-\frac{10m}{999\cdot1008}=-\frac{741509}{503496}<0.
\]
This is a counterexample to the *extended basis claim*, not to NF02 positivity. Full sharp NF02 coverage remains at the preceding checkpoint: all-gap finite-mass sign for \(J\le90\) or \(\hat b\le60\), exact spectral realization and the retained all-gap hybrid floor. The unbounded region \(J\ge91,\hat b>60\) still requires a uniform signed finite-mass proof. The fixed-gap certificate at 1000 is not transferred to this new family by monotonicity. It remains a separate and much stronger bound at moderate gaps.

Run `python check_final.py` for the exact regression checks. The proof of the infinite domain uses the positive-basis theorem and algebraic inequalities, not extrapolation of these tests.
