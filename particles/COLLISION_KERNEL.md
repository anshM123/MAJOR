# An exact positive kinematic measure for the two-particle formation kernel

Research result, 28 November 2026. The single-interval input is the explicit
two-massless-particle black-hole formation calculation of Aref'eva,
Khramtsov, and Tikhanovskaya, arXiv:1706.07390v3, equations (3.40), (4.5),
and (4.6). The result below proves an all-time statement about that interval
kernel and its noncrossing-matching functional. Identifying this functional
with physical HRT for every union still requires simultaneous admissibility.

## 1. The exact kernel

The boundary circle has circumference \(2\pi\). The injection points are
\(0,\pi\), the final nonrotating BTZ radius is \(R>0\), and \(t\ge0\) is
the BTZ rest-frame boundary time. This coordinate is not global collision
time: the particles have \(r(t)=R\coth(Rt)\) and approach the final horizon
as \(t\to\infty\). The original global-coordinate collision is not a finite
event on this time axis.

For endpoints whose shorter arc crosses the injection at zero, choose

\[
-\pi<a<0<b<\pi,\qquad b-a\le\pi,
\]

and define

\[
X=e^{Ra},\quad Y=e^{Rb},\quad Z=e^{Rt},\quad
H=\frac{\pi R}{2},\quad E=\coth H.
\]

Thus \(0<X<1<Y\), \(Z\ge1\), and
\(d=\tfrac12\log(Y/X)\le H\). Omit a common endpoint term for the moment.
The direct and crossing lengths are

\[
D_{\mathrm d}=\log\frac{(Y-X)^2}{XY},
\qquad
D_{\mathrm c}=\log\frac{F_1F_2}{XYZ^2},
\tag{1}
\]

where

\[
\begin{aligned}
F_1&=E(X-Z)(Y-Z)+Z(Y-X),\\
F_2&=E(XZ-1)(YZ-1)+Z(Y-X).
\end{aligned}
\tag{2}
\]

Expansion of (2) gives exactly the logarithm argument in the published
crossing formula (3.40); the companion script verifies this identity
symbolically. The published interval kernel is

\[
D_t(a,b)=\min\{D_{\mathrm d},D_{\mathrm c}\}
 +2\log(1/\epsilon).
\tag{3}
\]

If the shorter arc crosses no injection point, use the direct branch. If
it crosses \(\pi\), translate the crossing chart by \(\pi\). The circle
kernel is symmetric under complementation, so for arcs longer than \(\pi\)
use their shorter complement. The endpoint cutoff in (3) is our consistent
dimensionless prescription; no squared-cutoff convention from the source
is imported.

## 2. Positive factors and endpoint monotonicity

Put \(r_i=F_i/[Z(Y-X)]\). Completing the two elementary arithmetic-geometric
mean bounds gives

\[
\begin{aligned}
r_1&=1+E\frac{Z+XY/Z-X-Y}{Y-X},\\
r_2&=1+E\frac{XYZ+1/Z-X-Y}{Y-X},
\end{aligned}
\]

and hence, for either factor,

\[
r_i\ge1-E\frac{\sqrt Y-\sqrt X}{\sqrt Y+\sqrt X}
 =1-E\tanh(d/2)
 \ge\frac1{1+\cosh H}>0.
\tag{4}
\]

There is therefore no hidden logarithm-sign or branch ambiguity. Let
\(h=D_{\mathrm c}-D_{\mathrm d}=\log(r_1r_2)\). Exact differentiation gives

\[
\begin{array}{ll}
\displaystyle \partial_Xr_1=\frac{E(Y-Z)^2}{Z(Y-X)^2},&
\displaystyle \partial_Yr_1=-\frac{E(X-Z)^2}{Z(Y-X)^2},\\[6pt]
\displaystyle \partial_Xr_2=\frac{E(YZ-1)^2}{Z(Y-X)^2},&
\displaystyle \partial_Yr_2=-\frac{E(XZ-1)^2}{Z(Y-X)^2}.
\end{array}
\tag{5}
\]

Consequently

\[
h_a\ge0,\qquad h_b\le0.
\tag{6}
\]

These signs hold throughout the physical shortest-arc patch, including
times after equilibration. They are the ingredient needed to control the
minimum of the two branches; positivity of the separate branches alone
would not suffice.

## 3. Positivity through every branch change

The determinant of each bilinear polynomial is particularly simple:

\[
F_i\partial_X\partial_YF_i
 -(\partial_XF_i)(\partial_YF_i)=Z^2.
\tag{7}
\]

It follows that

\[
\partial_a\partial_b D_{\mathrm c}
 =R^2XY Z^2(F_1^{-2}+F_2^{-2})>0,
\qquad
\partial_a\partial_b D_{\mathrm d}
 =\frac{2R^2XY}{(Y-X)^2}>0.
\tag{8}
\]

For a regular switch curve, the additional mixed measure is formally
\(-\delta(h)h_ah_b\), which is nonnegative by (6). The following argument
also covers tangencies and singular switch loci. Choose smooth concave
functions \(f_\eta\to\min(0,\cdot)\), with
\(0\le f_\eta'\le1\) and \(f_\eta''\le0\). Then

\[
\begin{aligned}
\partial_a\partial_b\bigl(D_{\mathrm d}+f_\eta(h)\bigr)
={}&(1-f_\eta')\,\partial_a\partial_bD_{\mathrm d}
 +f_\eta'\,\partial_a\partial_bD_{\mathrm c}\\
 &+f_\eta''h_ah_b\ \ge0.
\end{aligned}
\tag{9}
\]

Taking the distributional limit proves positivity without a regularity
assumption on the transition locus.

### Injection-point seams

At \(a=0\), \(t=0\), both factors in (4) equal one, so the crossing
and direct values agree. Their tangential \(b\)-derivatives agree as well.
For \(t>0\), the direct branch is selected in a neighborhood of this seam.
One way to see this last fact is to rewrite the crossing argument as

\[
4\left[Q^2+\sinh^2m(E^2-P^2)\right],\quad
P=\sinh d-E\cosh d,\quad Q=E\cosh(Rt)+\cosh m\,P,
\tag{10}
\]

where \(m=R(a+b)/2\), \(|m|\le d\). Here \(P<0\) and

\[
Q\ge \sinh d\,[\cosh d-E\sinh d]\ge0.
\]

Increasing \(t\) strictly increases (10) away from the degenerate boundary
point, while the direct length is fixed. At the seam their initial values
are equal. The same argument applies at \(b=0\) and at the other injection.
Thus patching across these seams creates no negative mixed measure.

### Antipodal seams

At \(b-a=\pi\), changing to the complementary injection chart sends
\((a,b)\mapsto(b-\pi,a+\pi)\). It leaves \(m\) fixed and sends
\(d\mapsto2H-d\). At \(d=H\),

\[
P=-\frac1{\sinh H},\qquad P'=\cosh H-E\sinh H=0.
\]

Equation (10) therefore has \(\partial_dD_{\mathrm c}=0\) there. The
direct branch has \(\partial_dD_{\mathrm d}=2\coth H>0\). After reflection
about \(d=H\), any derivative jump is downward in the separation coordinate.
Since \(\partial_a\partial_b=(R^2/4)(\partial_m^2-\partial_d^2)\), its
mixed measure is nonnegative. This includes the case when the two branches
tie at the antipodal seam.

## 4. The all-time theorem and its graph corollary

**Theorem K.** For every \(R>0\) and \(t\ge0\), the symmetric equal-time
circle interval kernel (3), patched as above, has a nonnegative
distributional kinematic measure away from coincident endpoints:

\[
\boxed{\partial_a\partial_b D_t\ge0.}
\tag{11}
\]

Equivalently, every cyclically ordered quadruple obeys

\[
D_t(a,c)+D_t(b,d)\ge D_t(a,d)+D_t(b,c).
\tag{12}
\]

The proof uses the exact published kernel and elementary algebra. It does
not invoke maximin, the null energy condition, or a guessed Lorentzian
triangle inequality. It is a theorem about the kernel; the geometric
derivation of that kernel remains the stated input.

Fix finitely many distinct boundary endpoints \(x_0,\ldots,x_{m-1}\),
ordered cyclically, and set \(D_{ii}=0\). For a sufficiently small fixed
cutoff, all discrete weights

\[
w_{ij}=\frac12\left[D_{ij}+D_{i+1,j+1}-D_{i,j+1}-D_{i+1,j}\right]
\tag{13}
\]

are nonnegative, with indices interpreted cyclically. Nonadjacent weights
are rectangle integrals of (11); adjacent weights have a positive endpoint
divergence, so sufficiently small cutoff makes the finite set nonnegative.
They reconstruct every interval exactly by telescoping:

\[
D_{ab}=\sum_{i\in[a,b)}\sum_{j\notin[a,b)}w_{ij}.
\tag{14}
\]

Apply the finite disk chord-arrangement theorem from the standing
`strips/FINITE_GRAPH.md` (or the corresponding disk full-cone theorem).
The noncrossing-matching functional

\[
F_t(A)=\min_{\Pi\in\mathrm{NC}(\partial A)}
          \sum_{\{p,q\}\in\Pi}D_t(p,q)
\tag{15}
\]

has one finite nonnegative-capacity graph representing all unions.
Consequently, for every number of parties and every time,

\[
\boxed{\mathbf F_t\in\mathcal C_N^{\mathrm{RT}}.}
\tag{16}
\]

This statement is exact for the defined functional. The arrangement graph
is essential: replacing it by the bare graph of Bell edges between boundary
atoms gives different disconnected-region entropies.

## 5. Precisely what is still needed for physical collision HRT

Every actual HRT cut on a filled-disk homology slice induces a noncrossing
endpoint matching. Its component costs are no smaller than the globally
least single-interval branch costs, provided that the published interval
list exhausts regular admissible component geodesics. This yields a lower
bound by (15). The reverse direction is the missing statement:

> A minimizing matching in (15) admits choices of its direct/crossing
> geodesic representatives that form one embedded acausal multicurve and
> lie on one homology Cauchy slice of the two-particle spacetime.

Endpoint noncrossing by itself does not prove this Lorentzian assertion.
Cone smoothing and finite-group quotient descent do not prove it either.
Once this geometric lemma is supplied, (16) immediately becomes the desired
arbitrary-party black-hole-formation theorem for this specific model.
Until then, (16) and the numerical outputs must be labelled as matching
functional results, not an unconditional all-union HRT theorem.

The independent geometric audit `COLLISION_ADMISSIBILITY.md` supplies further
information about actual branch existence and the obstruction to a direct
vacuum invisible-domain embedding.

## 6. Reproducible numerical audit

Run:

```bash
python code/bilocal_collision.py \
  --facets data/facets_n5_all372.json \
  --output data/collision
```

The computation uses eight distinct ten-endpoint arrangements, three values
\(R=0.5,1,2\), 81 times in \([0,2]\), and 3, 4, or 5 nonempty parties.
For fewer than five parties, the unused parties are empty; evaluation of
all 372 five-party facets therefore also includes the induced lower-party
inequalities. The 372 facet vectors are the complete five-party catalogue
associated with Hernandez-Cuenca, arXiv:1903.09148, with source URL recorded
in the input JSON.

Results:

- 5,832 matching entropy vectors; 2,169,504 facet evaluations.
- Minimum residual: \(-1.4211\times10^{-13}\), consistent with arithmetic
  error at exact saturation; the stated acceptance tolerance is \(10^{-10}\).
- Minimum discrete kinematic weight: \(2.74628\times10^{-4}>0\).
- Maximum reconstruction error in (14): \(2.8422\times10^{-14}\).
- Exact symbolic checks of (1), (5), and (7) all pass.
- Symmetric-interval saturation at \(t=(b-a)/2\) agrees within
  \(1.78\times10^{-15}\).

The samples check the implementation. The all-time cone statement for the
matching functional follows from Sections 2--4, not from time sampling.
The figure labels its multiparty panels accordingly and does not mark a
fictional finite collision time on the BTZ time axis.

## Primary references

- I. Ya. Aref'eva, M. A. Khramtsov, M. D. Tikhanovskaya, *Thermalization after
  holographic bilocal quench*, JHEP 09 (2017) 115,
  https://arxiv.org/abs/1706.07390.
- S. Hernandez-Cuenca, *The Holographic Entropy Cone for Five Regions*,
  Phys. Rev. D 100 (2019) 026004, https://arxiv.org/abs/1903.09148.
- H.-J. Matschull, *Black Hole Creation in 2+1 Dimensions*, Class. Quantum
  Grav. 16 (1999) 1069, https://arxiv.org/abs/gr-qc/9809087.

The algebraic factorization, switch-positivity proof, and all-party matching
corollary are deductions made in this research run; no claim of literature
priority is made without a separate exhaustive novelty review.
