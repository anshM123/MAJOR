# Direct sliver continuation — 2026-11-26

## Scope and provenance

This is a delta to the included original `sliver_direct_certificate.zip`, not an unconditional proof of the full MAJOR sliver. Its upstream assumptions are the physical angular-action identity, spectral construction and envelopes, and the coherent seed-error budgets in the original archive and MAJOR `ads4/ANGULAR_ROUTE.md` §§115–116. The direct-action coordinate is c=1-exp(-L). The source coefficient file encloses the seed-substituted vector and has 1253 pairs, indices 0 through 1252. The source spectrum has all 3564 repaired Gaussian triples. No coefficient or weight in the original inventory has been changed. Large-L coverage is retained from the earlier work, not rerun here.

## Exact continuum enclosure

The new direct cells use the original integer GMP tensor-Chebyshev checker. The triangle cells use a new, mathematically equivalent signed Horner evaluator, storing A_k,G_k multiplied by 2^(-k), so the Horner ratio is doubled and the next A coefficient in the dressing is multiplied by two. Coefficient balls are transferred by integer division and an outward unit allowance. The truncation set T={(i,j):a_c i+a_h j<Q} is downward closed. Thus all computed Taylor coefficients within T are unaffected by omitting higher indices. The polydisk radii are R_c=2^a_c and R_h=2^a_h, in the rescaled cell coordinates. The integer majorant v encloses the absolute complex Horner recurrence in units S=10^D; poles are excluded by the checked denominator gaps.

The omitted Taylor coefficients have total absolute sum bounded by (v/S) sum_{m>=Q}(m+1)2^-m = (v/S)(Q+2)2^(1-Q), since the number of pairs of weighted degree m is at most m+1. Per-step directed arithmetic error is bounded by the l1 recurrence in the C++ code; with d=den-pH0, r=pDH, a=2(C0-pH0), dc=2DC, it is E_new=ceil(((|a|+dc+2r)E+terms*d)/(d-r))+Fr. This follows from the nonnegative-coefficient geometric inverse of d-r*y, and at most one integer unit of division error per retained coefficient. The map to tensor Chebyshev has nonnegative column weights summing to one, so does not amplify l1 error. The final constant coefficient less the other absolute coefficients, arithmetic error, and Cauchy tail is a lower bound on the full rectangle; it is then physically normalized and charged the retained coherent seed and prefix errors. The certificate covers all physical spins and twists in each cell, not merely sampled states.

For the final tight cell, the old 2^-110 physical tail requirement is relaxed only locally, to 3e-6, because a larger error can still prove positivity. The exact bound below is 2.590095e-6 (rounded upward). The local sign test subtracts 3e-6+1e-18. It still uses every coefficient in the existing table (K=1251, A_(K+1) required for the dressing). See the final ledger for the actual status. No positivity is inferred from a prefix bound alone.

## Two new prefix bounds

Write A(z) and G(z) for the two integrated spectral generating functions; the retained total absolute spectral masses are bounded by 60000001 and 60000000, respectively. For the repaired spectrum above 1040 the weighted integral is bounded by `(2/998)*gaussian_mass(1040,B)`; the ideal spectrum has no support there. The Gaussian envelope is inherited. Write F_{p,nu}(z)=(1-z)^(-p) P_{-1/2+i nu}((1+z)/(1-z)), p=3 or 4. Use the even-in-beta Mehler-Dirichlet representation

P(cosh beta) = (2/pi) int_0^1 cos(nu beta x) [(1-x^2) sinhc((1+x)beta/2) sinhc((1-x)beta/2)]^(-1/2) dx,

where sinhc(v)=sinh(v)/v. This follows by putting t=beta*x in the real Mehler-Dirichlet integral and using cosh(beta)-cosh(t)=2 sinh((beta+t)/2)sinh((beta-t)/2). Its analytic continuation to |Im beta|<pi is justified by the nonvanishing of sinhc there, with the branch fixed at beta=0. For |Im v|<=B<pi, |sinhc(v)|>=sin(B)/B, by sinh^2(Re v)>=Re(v)^2 and sin^2(Im v)>=(Im v)^2 (sin B/B)^2. Consequently |P|<=(B/sin B)*exp(B nu).

### Tight ordinary-power tail

On |z|<=R=199/200, beta=2 atanh(sqrt(z)), so |Im beta|<pi/2<11/7, and the prefactor is at most pi/2<2. The combined integrated envelope is

M=2(1-R)^(-4)[120000001*2^2358+(2/998)*gaussian_mass(1040,11/7)] < 2^2417,

using 1040*(11/7)<2358 log 2. Cauchy therefore bounds |A_k|+|G_k|<=2^2417 R^(-k). Define B via 2B=(z-z^2)A''+[1+(7-4b)z]A'+[9-2b(b-1)]A. The retained coefficient formula and the bound on its quadratic term show |B_k|<=2(k+4)^2*2^2417 R^(-k). Physical U_k<=c^k, h<=H=998^(-2), q=b(b+1)<=20, and the physical prefactor is <=1. At t=c_hi/R the sum of all omitted normalized terms (k>K) is thus bounded by

2^2417 t^(K+1) [(1+20H)/(1-t) + 2H((K+5)^2/(1-t)+2(K+5)t/(1-t)^2+t(1+t)/(1-t)^3)].

This is <3e-6 at c_hi=.258, K=1251, by exact rational arithmetic. The old coherent seed bounds apply to each finite ordinary-power prefix as well, since their proofs majorize absolute per-k contributions using U_k<=c^k.

### Conformal prefix theorem (exact-weight vector only)

Let w=(1-sqrt(1-z))/(1+sqrt(1-z)), z=4w/(1+w)^2. Then beta=4 atanh(sqrt(w)) and (1-z)^(-p)=((1+w)/(1-w))^(2p). For |w|<=R_c=.99, |Im beta|<=B0=4 atan(sqrt(R_c))<pi, and sin B0=4 sqrt(R_c)(1-R_c)/(1+R_c)^2>.0099. Hence |P|<=320 exp((22/7)nu), since B0/sin B0<320. The combined spectral sup bound is at most

320*199^8 [120000001*2^4716+(2/998)*gaussian_mass(1040,22/7)] < 2^4813,

using 1040*(22/7)<4716 log 2. Thus the transformed coefficient pairs satisfy |Atilde_j|+|Gtilde_j|<=2^4813 R_c^(-j). Their exact conversion, for j>=1, is Atilde_j=sum_(k=1)^j 4^k (-1)^(j-k) binom(j+k-1,j-k) A_k, and likewise for G; the constant coefficients are unchanged.

Define the physical functional T_h by T_h[z^k]=U_k with U_k=0 for k>n. Lagrange inversion of w=z(1+w)^2/4 gives [z^k]w^j=(j/k)4^(-k)binom(2k,k-j) >=0 for k>=j>=1. This implies T_h[w^j]<=w(c)^j. The physical normalized action is T_h applied to (1+qh)A+h(B+sigma G), times (1-c)^3/(1+qh)^4, where sigma=+1 or -1. In the transformed truncation B is obtained by applying the dressing differential operator to the truncated Atilde(w(z)), not by independently truncating a different B series.

For c<=.9, w(c)<=W=13/25, w'(c)<=4 and w''(c)<=26. Nonnegative power-series coefficients of f=w^j and its derivatives give |T_h[B_b f]|<=f''(c)+5 f'(c)+7.5 f(c)<=50(j+1)^2 W^(j-2), j>=1. With h<=998^(-2), the corresponding combined bound for each pair of coefficients (Atilde_j,Gtilde_j) is <=3(j+1)^2 W^(j-2). Thus at t=W/R_c the normalized error of truncating at J is bounded by

3*2^4813*W^(-2)*t^(J+1) [(J+2)^2/(1-t)+2(J+2)t/(1-t)^2+t(1+t)/(1-t)^3].

At J=5333 this is <2^-110, verified by exact rational arithmetic in `prefix_improvements.py`. This replaces the earlier largest ordinary-power prefix K=17847 (17849 pairs including the next coefficient needed by dressing) by 5334 exact coefficient pairs. The bound is a physical-tail result, not a sign certificate. The coherent seed substitutions have not yet been budgeted for the transformed finite approximation. They must be handled separately, or exact-weight coefficient enclosures must be used. No transformed moment table has been generated here.

### Finite-width recurrence for the transformed ladder

The ordinary ladder gives T_h[(z-c)f+h(1-z)D_z(D_z+s)f]=0, where D_z=z*d/dz. This holds coefficientwise even past termination, since U_k=0 beyond n. In w, D_z=w(1+w)/(1-w)*d/dw, and (1-z)D_z(D_z+s)f=w^2 f''+w[1+2w/(1-w^2)+s(1-w)/(1+w)]f'. Put f=w^j(1-w^2)^2 and clear the denominators. With V_j=T_h[w^j], the identity is sum_{i=0}^4 C_i V_(j+i)=0, where

- C0=-c+h*j*(j+s)
- C1=4-2*h*j*(s-1)
- C2=2*c-2*h*j*j-8*h*j-4*h*s-8*h-8
- C3=4+2*h*(j+4)*(s-1)
- C4=-c+h*(j+4)*(j+4-s).

V0=1, V_j=0 for j>n, and 0<=V_j<=w(c)^j. `verify_conformal_recurrence.py` checks the polynomial identity symbolically and on 1980 exact finite-spin recurrence instances. A stable, continuum-certified summation algorithm for this new recurrence remains to be implemented. The recurrence alone is not a positivity argument.

## Reproduction

Extract the included original archive; copy the delta files into `sliver_direct_certificate/`. Keep original SHA256SUMS unchanged (it describes the original archive). The delta's checksum file covers the delta only. From that directory, with Python 3, SymPy, C++17 and GMP:

```
g++ -O3 -std=c++17 direct_sign_cheb.cpp -lgmpxx -lgmp -o direct_sign_cheb
g++ -O3 -std=c++17 direct_sign_triangle.cpp -lgmpxx -lgmp -o direct_sign_triangle
python prefix_improvements.py
python verify_conformal_recurrence.py
for b in 0 4; do
  python direct_sign_cheb.py .09 .12 "$b"
  python direct_sign_cheb.py .12 .15 "$b"
  python direct_sign_triangle.py .15 .18 "$b"
  python direct_sign_triangle.py .18 .223 "$b"
  python direct_sign_tight.py .223 .258 "$b"
done
python verify_extension.py
```

The verification script recomputes transfer of every input coefficient enclosure from the old table, exact majorants and Taylor tails at the chosen radii, real and complex pole gaps, spin caps, continuity of coverage and all final lower-bound arithmetic. It does not replace running the C++ sign calculations, nor independently reproduce the original spectral integrations. The exact `.in` and `.out` files are included for each completed sign computation. Runtime varies by hardware. Read `extension_ledger.json` for actual verified coverage.
