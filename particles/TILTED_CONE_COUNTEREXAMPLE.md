# A static positive-mass cone disproves the naive chiral cone-length formula

28 November 2026. The example was proposed by the cusp agent and its
spacetime admissibility and cutoff normalization checked here. This is a
counterexample to an attempted formula, not to an entropy inequality.

## 1. Geometry and endpoints

Use radius-one conical AdS,

\[
 ds^2=-\cosh^2r\,dt^2+dr^2+\sinh^2r\,d\phi^2,
 \qquad \phi\sim\phi+\theta,\qquad\theta=\frac{3\pi}{2}.
\]

The particle has a positive deficit and zero spin. Take endpoints

\[
 (t_0,\phi_0)=(0,0),\qquad
 (t_1,\phi_1)=\left(\frac{9\pi}{16},\frac{11\pi}{16}\right).
\]

They fit on a spacelike boundary Cauchy graph: the two constant slopes
needed around the circle are \(9/11\) and \(-9/13\), both of magnitude
less than one. Corners can be smoothed while fixing the endpoint values.
For the homology-slice argument below, one may choose the boundary trace
of the constructed spacelike graph; only these endpoint data are needed.

The chiral angular sweeps of the direct marked class are

\[
 \beta_L=\Delta\phi+\Delta t=\frac{5\pi}{4},\qquad
 \beta_R=\Delta\phi-\Delta t=\frac{\pi}{8}.
\]

## 2. The direct spacetime geodesic exists, stays regular, and is acausal

The metric above follows directly from the ambient embedding

\[
 X=(\cosh r\sin t,\sinh r\cos\phi,
       \sinh r\sin\phi,\cosh r\cos t),
 \qquad \langle X,X\rangle=-1,
\]

with signature \(\operatorname{diag}(-,+,+,-)\), followed by the angular
wedge identification. No extrapolation of a cone-length formula is used
to compute the spacetime geodesic. Take its boundary representatives

\[
 p_0=(0,1,0,1),\qquad
 p_1=(\sin\Delta t,\cos\Delta\phi,\sin\Delta\phi,\cos\Delta t).
\]

Both are null boundary representatives. Put

\[
 N=-2\langle p_0,p_1\rangle
 =2(\cos\Delta t-\cos\Delta\phi)>0.
\]

The direct geodesic is

\[
 X(s)=\frac{e^s p_0+e^{-s}p_1}{\sqrt N}.
\]

Writing the physical circumferential radius as \(\rho=\sinh r\), its
minimum radius obeys

\[
 \rho_{\min}^2
 =\frac{1+\cos\Delta\phi}{\cos\Delta t-\cos\Delta\phi}
 =1.2328835904980993>0.
\]

The spatial projection is a positive linear combination of two nonopposite
spatial endpoint vectors. Its continuous angular sweep is exactly
\([0,11\pi/16]\); it never reaches or winds around the particle.

Global acausality survives the cone quotient. In the conformal ultrastatic
description the spatial optical metric is

\[
 h=d\chi^2+\sin^2\chi\,d\phi^2,
 \qquad 0\le\chi\le\frac\pi2,
 \quad \phi\sim\theta.
\]

For any two points on the geodesic, their angular separation is at most
\(11\pi/16<\min(\pi,\theta/2)\). The shortest optical path is therefore
the ordinary spherical-disc path in the direct developed sector. Any
other angular lift has a larger absolute sweep; a route through the apex
is also no shorter. The developed AdS geodesic is acausal, as is also seen
from \(\langle X(s),X(s')\rangle=-\cosh(s-s')<-1\) and its global-time
range of width less than \(\pi\). Hence no shortcut in the cone quotient
causally relates two distinct points of this geodesic.

For an explicit homology-slice existence argument, its compactified
projection into the optical cone is injective. The time function on that
projection is uniformly Lipschitz with some constant \(L<1\): locally this
follows from its uniform spacelike tangent, including the orthogonal
conformal endpoints, and away from the diagonal it follows from strict
acausality and compactness. The McShane extension

\[
 f(x)=\inf_{k\in K}\{t(k)+L\,d_h(x,k)\}
\]

extends this time data to the full compact optical cone with the same
Lipschitz constant. Its graph is acausal and Cauchy and contains the
geodesic. The spatial arc separates the filled disc and gives the required
homology region. The graph can be smoothed away from the prescribed arc;
no particle endpoint or extra homology boundary is introduced.

## 3. It is the shortest physical branch

Use a common radial cutoff \(\rho=\rho_c\) and subtract
\(2\log\rho_c\). The direct renormalized length is

\[
 L_0=\log\left[4\sin\frac{5\pi}{8}\sin\frac\pi{16}\right]
 =-0.32717186862168246.
\]

The other regular angular branch has sweep
\(\Delta\phi-\theta=-13\pi/16\), giving

\[
 L_{-1}=\log\left[4\sin\frac\pi8\sin\frac{11\pi}{16}\right]
 =0.24118665553584978>L_0.
\]

All other continuous lifted angular sweeps have magnitude greater than
\(\pi\), so no regular boundary-to-boundary spacetime geodesic in those
marked classes exists with these endpoint times. Reducing an angular sweep
modulo \(2\pi\) in the developing image would change the lifted endpoint
and therefore the class; it cannot create an extra branch. An elliptic
particle holonomy has no regular closed spacelike geodesic to add.
Consequently the direct curve is the physical HRT minimum in the usual
regular-locus prescription with filled homology.

## 4. The marked chiral cone functional gives a different answer

On a hyperbolic cone of angle \(\theta\), a regular complete geodesic has
continuous angular sweep of magnitude less than \(\pi\). Thus the direct
left sweep \(5\pi/4\) has no regular representative. Its relaxed marked
class infimum runs radially into and back out of the cone point and has
renormalized length \(2\log2\), in the same cutoff convention. The right
sweep is regular and has length \(2\log(2\sin(\pi/16))\).

The half-sum in this class is consequently

\[
 F_0=\log\left[4\sin\frac\pi{16}\right]
 =-0.24799827671149502.
\]

For the other common marked class, the left sweep is \(-\pi/4\) and the
right is \(-11\pi/8\), so

\[
 F_{-1}=\log\left[4\sin\frac\pi8\right]
 =0.42574718219016017.
\]

For embedded arcs in the particle-marked disc, these are the two eligible
simple classes. Consequently their minimum is \(F_0\). Even if one
enlarges the formal list to all integer angular windings, the additional
entries have both relaxed lengths equal to \(2\log2\) and half-sum
\(\log4\), and do not change the minimum. Therefore

\[
 \boxed{\min_m F_m-L_{\mathrm{HRT}}
 =-\log\sin\frac{5\pi}{8}
 =0.07917359191018744>0.}
\]

There is no common marked class with regular representatives in both
chiral cone metrics for these endpoints. Interpreting the chiral formula
using only regular representatives therefore makes it undefined; relaxing
the cone-class minima makes it defined but wrong.

This failure occurs in a static, spinless, positive-mass geometry, on a
perfectly admissible tilted boundary Cauchy slice. A theorem about a
mathematical common-class cone-length functional cannot be identified with
physical HRT by importing the vacuum formula unchanged.

Nor can an endpoint-local change of renormalization convention restore
the branchwise formula. For the same two fixed endpoints the two branch
discrepancies are different:

\[
 F_0-L_0=-\log\sin\frac{5\pi}{8},\qquad
 F_{-1}-L_{-1}=-\log\sin\frac{11\pi}{16}.
\]

An endpoint counterterm would shift both branches by the same amount.
