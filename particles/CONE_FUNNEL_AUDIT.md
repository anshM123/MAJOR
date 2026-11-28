# Independent audit: funnel replacement and the relaxed common-class cone theorem

28 November 2026. This audits `../cone/SMALL_ANGLE_FUNCTIONAL_THEOREM.md` and the proposed all-angle relaxed extension. It concerns a Riemannian common-class functional, not its identification with spacetime HRT.

## 1. The radial replacement is valid

In a cone chart with angular period theta, the original metric is
\(dr^2+\sinh^2r\,d\phi^2\). For a small a, choose \(H\ge0\) smooth with support in \((a/3,2a/3)\), satisfying

\[
\int e^{-s}H(s)\,ds<1<\int e^sH(s)\,ds.
\]

This is possible because the support is in positive s. Define

\[
f(t)=\sinh t+\int_t^a\sinh(s-t)H(s)\,ds.
\]

Then \(f''=f+H\). Above the support, \(f=\sinh t\); below it,

\[
f(t)=Ae^t+Be^{-t},\qquad
A=\tfrac12\left(1-\int e^{-s}H(s)\,ds\right)>0,
\quad B=\tfrac12\left(\int e^sH(s)\,ds-1\right)>0.
\]

On the support interval, \(\sinh t>0\) and the integral is nonnegative. Thus f is everywhere positive. The metric
\(dt^2+f(t)^2d\phi^2\) is complete at its new end and is smooth with
\(-C\le K=-1-H/f\le-1\), for a finite C depending on the chosen replacement. It is exactly hyperbolic sufficiently far down the new funnel. No curvature bound uniform as \(a\to0\) is needed.

Let \(m=\min f>0\). Choose \(0<c_0<\min(a/3,\operatorname{arsinh}m)\), and a smooth strictly increasing radial diffeomorphism \(r:\mathbb R\to(0,\infty)\) with

\[
r(t)=t\ (t\ge c_0),\qquad
0<r(t)\le c_0\ (t\le c_0),\qquad 0<r'(t)\le1.
\]

Such a map exists: prescribe a positive smooth derivative q, equal to 1 in a small neighborhood below \(c_0\) and above it, with \(q\le1\) and \(\int_{-\infty}^{c_0}q=c_0\). A positive decaying tail and a movable smooth transition supply the desired integral. Integrating q gives r.

For \(t\le c_0\), \(\sinh r(t)\le\sinh c_0<m\le f(t)\). For \(t\ge c_0\), \(r=t\) and \(f(t)\ge\sinh t\). Thus

\[
(r')^2dt^2+\sinh^2r(t)\,d\phi^2
\le dt^2+f(t)^2d\phi^2.
\]

This is a 1-Lipschitz diffeomorphism to the old punctured cone collar, identical near its exterior. It preserves the punctured marking. The angular period imposes no restriction on this construction.

## 2. Below pi, the exact-preservation argument is consistent

For an embedded regular complete geodesic, the developed collar equation is
\(\tanh r\cos(\phi-\phi_0)=\tanh r_0\). If

\[
r_0<\operatorname{arctanh}\{\tanh R\cos(\theta/2)\},
\]

the two points at angular offsets \(\pm\theta/2\) occur before leaving the cone chart, have equal radius, and are identified transversely. This contradicts embeddedness. Thus the stated positive collar exists when \(\theta<\pi\).

Any regular simple geodesic representative therefore survives a sufficiently small replacement unchanged. It is then the unique geodesic in the same class of the smooth negative-curvature surface, proving exact marked-length equality. The old regular-realization theorem is needed for this **attained, exact-preservation statement**. The relaxed extension below avoids that input entirely.

## 3. Precise relaxed functional and all-angle convergence

Let \(\mathcal C(A)\) be the common set of reduced embedded multicurve classes in the punctured marked surface, whose homology trace is A in the filled surface. Particle ends carry no prescribed label. For each old metric j, set

\[
L_j(c)=\inf_{\gamma\in c\;\mathrm{regular,embedded}}\ell_j^{\rm ren}(\gamma),
\qquad
F_{\rm rel}(A)=\inf_{c\in\mathcal C(A)}\sum_jw_j L_j(c).
\]

Representatives may be taken to have vertical geodesic tails in fixed physical endpoint charts, without changing the infimum: replace the tail from a marked ideal point to its first exit from the local chart by the corresponding vertical segment. The remaining finite part is compact and avoids the finite particle set. No old minimizing geodesic or minimizing class is assumed to exist.

Replace every cone at radii tending to zero, obtaining smooth negatively curved metrics \(\widehat h_{j,a}\) on the same punctured marked surface, now with auxiliary funnel ends. Let \(\widehat L_{j,a}(c)\) be their marked geodesic length and \(\widehat F_a\) their common-class infimum with free auxiliary-end colors.

**Lower bound.** The 1-Lipschitz maps preserve the class and the unchanged physical endpoint regulator. Therefore

\[
\widehat L_{j,a}(c)\ge L_j(c),\qquad
\widehat F_a(A)\ge F_{\rm rel}(A).
\]

**Upper bound.** Given \(\delta>0\), choose one common class c with
\(\sum_jw_jL_j(c)<F_{\rm rel}(A)+\delta\). For each j independently, choose a regular embedded representative \(\gamma_j\in c\) with weighted total cost at most another delta above \(\sum_jw_jL_j(c)\). There are finitely many metrics and finitely many particle points. Each truncated representative is compact and disjoint from them, so for all sufficiently small a every replacement misses every \(\gamma_j\). The representatives have exactly their old lengths in the new metrics. Minimization then gives

\[
\widehat F_a(A)\le\sum_jw_j\ell_j(\gamma_j)
<F_{\rm rel}(A)+2\delta.
\]

Consequently

\[
\boxed{\widehat F_a(A)\longrightarrow F_{\rm rel}(A).}
\]

There is no interchange of a limit with an infimum over an uncontrolled set: the proof fixes one delta-optimal class and finitely many representatives before choosing a. For finitely many party unions, the same choice of sufficiently small a works for all of them. Monotonicity in a is unnecessary.

The angle restriction is absent from this argument. For the particle problem it covers every fixed \(0<\theta<2\pi\), including angles at or above pi whose marked infima may hit the cone. In fact the radial construction itself only requires a finite positive angle.

## 4. Uniform endpoint regulator, with the exact local integral

Choose unchanged physical hyperbolic charts
\(U_p=\{z\in\mathbb H^2:|z|>A_p\}\), placing p at infinity, with no other marked endpoint in the chart. The geodesic family

\[
T_p=\{(a,b):a<-A_p,\ b>A_p\}
\]

lies wholly in this chart. A geodesic cut ending at p has terminal vertical line \(x=x_0\) with \(|x_0|<A_p\). Otherwise that entire vertical line stays in the chart and its second ideal endpoint lies in its unmarked local boundary interval, contradicting the prescribed endpoint set. Thus a cutoff \(Y>A_p\) lies in the terminal segment for every class and every replacement parameter.

Put \(u_0=A_p+x_0\) and \(v_0=A_p-x_0\). With the repository's Liouville normalization, the T_p crossing integral below height Y is exactly

\[
I_p(Y,x_0)
=\int_{u_0}^{Y^2/v_0}\!du\int_{v_0}^{Y^2/u}\!\frac{dv}{(u+v)^2}
=\frac12\log\{(Y^2+u_0^2)(Y^2+v_0^2)\}
-\log(2A_pY),
\]

provided \(Y^2>u_0v_0\), which follows from \(Y>A_p\). Hence

\[
I_p(Y,x_0)=\log Y-\log(2A_p)+e_p(Y,x_0),
\qquad 0\le e_p(Y,x_0)\le\frac{2A_p^2}{Y^2}.
\]

The finite constant is independent of both the class and a. After the positive-current subtraction and restoring the exact local horoball cutoff, the Bell coefficient is

\[
t_p(\epsilon)=\sum_jw_j\log\frac1{2A_{j,p}\epsilon},
\]

with the same fixed \(A_{j,p}\) at every replacement radius. Choose epsilon below the fixed finite threshold making these coefficients nonnegative. Every smoothed functional then lies in the same closed fixed-party RT cone at that one regulator. Passing to the limit gives

\[
\boxed{\mathbf F_{\rm rel}^{\epsilon}\in\mathcal C_N^{\rm RT}.}
\]

For endpoint-balanced inequalities, the common Bell terms cancel as usual.

## 5. Marking and peripheral components

The comparison must use isotopy in the **punctured** surface; particle points are filled only for homology. Collapsing winding classes in the old definition would change the function being compared.

Each new funnel corresponds to one puncture and has a free color. A cut component peripheral to a single auxiliary funnel can be deleted and the funnel-side color reversed. This preserves every physical trace and decreases its positive length in every smoothed metric. In the old relaxed functional its peripheral infimum is zero. Removing such components therefore changes neither minimum. No auxiliary end may be silently assigned a fixed purifier label.

All retained simple classes do have regular unique geodesics in the new complete smooth negative-curvature metrics. This is the only class-realization input required in the relaxed proof. There is no reliance on old cone geodesic realization, old minimizer attainment, or joint compactness of old minimizing sequences.

## 6. Scope of the resulting theorem

The audit supports the all-angle **relaxed marked-length functional** theorem. The cone values are finite limits of exact smooth graph-model values. It also supports the exact regular theorem below pi, given the cited old regular-realization input.

It does not prove that a Lorentzian particle spacetime has HRT equal to this functional. In particular, it does not provide acausal homology slices, control collision vertices, or turn a cone-hitting infimum into an attained spacelike extremal curve. These distinctions remain necessary even though the Riemannian functional statement is now substantially stronger.
