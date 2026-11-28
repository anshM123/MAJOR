# Independent algebra audit of the proposed two-particle interval kernel

This note audits the root agent's candidate interval formulas only. It does not certify which geometric geodesics are admissible in the collision spacetime, or the simultaneous homology of a multi-interval matching.

Assume \(R>0\), \(H=\pi R/2\), \(E=\coth H\),
\(A=e^{Ra}<1<B=e^{Rb}\), \(T=e^{Rt}\ge1\), and
\(d=\tfrac12\log(B/A)\le H\). Define

\[
F_1=E(A-T)(B-T)+T(B-A),\qquad
F_2=E(AT-1)(BT-1)+T(B-A),
\]
\[
L_\times=\log\frac{F_1F_2}{ABT^2},\qquad
L_0=\log\frac{(B-A)^2}{AB},\qquad
R_i=\frac{F_i}{T(B-A)},\quad h=\log(R_1R_2).
\]

## 1. Denominators stay uniformly positive

For \(s>0\), \(s+AB/s\ge2\sqrt{AB}\). Apply this first to \(s=T\), and then after multiplying by \(AB\), to obtain

\[
R_i\ge1-E\frac{\sqrt B-\sqrt A}{\sqrt B+\sqrt A}
=1-E\tanh(d/2)
\ge1-\coth H\tanh(H/2)
=\frac1{1+\cosh H}>0.
\]

The lower bound is uniform in \(t\ge0\), and it does not require that the unconstrained minimum in T occur at \(T\ge1\).

## 2. Derivative signs and smooth-branch SSA

Direct differentiation gives

\[
\partial_A R_1=\frac{E(B-T)^2}{T(B-A)^2},\qquad
\partial_B R_1=-\frac{E(A-T)^2}{T(B-A)^2},
\]
\[
\partial_A R_2=\frac{E(BT-1)^2}{T(B-A)^2},\qquad
\partial_B R_2=-\frac{E(AT-1)^2}{T(B-A)^2}.
\]

Thus \(h_a\ge0\) and \(h_b\le0\). Each F is bilinear in A,B and satisfies the exact determinant identity

\[
F_iF_{i,AB}-F_{i,A}F_{i,B}=T^2.
\]

Consequently

\[
\partial_a\partial_b\log F_i
=\frac{R^2ABT^2}{F_i^2}>0,\qquad
\partial_a\partial_bL_0=\frac{2R^2AB}{(B-A)^2}>0.
\]

For \(L=\min(L_0,L_\times)=L_0+\min(0,h)\), its distributional mixed derivative is

\[
L_{ab}=\mathbf1_{h>0}(L_0)_{ab}
+\mathbf1_{h<0}(L_\times)_{ab}
-\delta(h)h_a h_b\ge0.
\]

At singular or tangent switch points, the conclusion follows by smooth concave approximations to \(\min(0,h)\): the extra term is \(f''(h)h_a h_b\ge0\), while \(0\le f'\le1\) makes the smooth curvature a positive convex combination of the branch curvatures. This avoids assuming that \(h=0\) is a regular level set.

## 3. The shortest-arc antipodal seam

Write \(m=(a+b)/2\) and
\(P(d)=\sinh d-E\cosh d\). Then

\[
\frac{F_1}{\sqrt{AB}T}=2\{E\cosh[R(m-t)]+P(d)\},
\quad
\frac{F_2}{\sqrt{AB}T}=2\{E\cosh[R(m+t)]+P(d)\}.
\]

At \(d=H\), \(P(H)=-1/\sinh H\) and \(P'(H)=0\), so
\((L_\times)_d=0\), whereas \((L_0)_d=2\coth H>0\).

The complementary shortest-arc chart must be the chart centered on the other collision line, with
\((a',b')=(b-\pi,a+\pi)\). It preserves m and sends \(d\mapsto2H-d\). Thus the normal derivative of the patched minimum has a downward jump or is continuous at d=H. Since
\(\partial_a\partial_b=\tfrac14\partial_m^2-\tfrac{R^2}{4}\partial_d^2\), its seam measure is nonnegative. The argument is invalid if one switches arc branches without the corresponding collision-chart transformation; continuity at the seam must be retained.

## 4. The endpoint-on-collision-line seam

At A=1,

\[
R_1R_2-1
=\frac{E(T-1)^2}{T^2(B-1)^2}
\left[T(B^2-1)-E(B-T)(BT-1)\right].
\]

The bracket divided by T equals

\[
EB(T+T^{-1})+(1-E)B^2-(1+E)
\ge(B-1)\{B+1-E(B-1)\}\ge0,
\]

where the last inequality is exactly \(B\le e^{2H}\). It is strict for T>1. Therefore the direct branch is active in an open neighborhood of the seam at positive time. The B=1 seam follows by endpoint reflection.

At T=1, the seam values match identically as a function of the other endpoint. Their tangential derivatives therefore match, so the mixed derivative has no seam delta (a jump of the normal derivative by itself does not create a mixed delta on a vertical or horizontal interface). Equivalently one may obtain positivity by the \(t\downarrow0\) distributional limit.

## Conclusion

The interval-kernel positivity argument is algebraically sound, including branch switches and chart seams, under the stated geometric candidate formula and shortest-arc chart prescription. The remaining issue is physical identification of this kernel and of the multi-interval matching formula with HRT in the collision spacetime; this calculation does not supply that geometric input.
