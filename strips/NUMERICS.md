# Vaidya–AdS4 strip entropy cone check

## Result

The thin-null-shell AdS4 calculation checks all **372 five-party RT facet rows** on **551 sampled boundary times**, for two strip configurations and their three-/four-party prefixes. This is **1,229,832 evaluated inequalities**. For the smaller numbers of parties, extra labels are set empty, so this uses redundant projections of the complete five-party system.

The smallest computed regulated margin is **−7.11 × 10⁻¹⁴** in area-density units, consistent with floating-point roundoff at saturated facets. No residual falls below the specified numerical tolerance of −2 × 10⁻⁸. This is a sampled numerical consistency check, **not** continuous-time verification or interval-arithmetic certification.

The two-strip mutual information in the disjoint configuration initially increases and then vanishes at

\[
t/z_h=0.2336006115\ldots.
\]

This illustrates a connectivity transition during thermalization while the sampled entropy vectors remain compatible with the static cone.

## Files

- `vaidya_strips.py`: complete computation and figure generation.
- `vaidya_data.npz`: interval areas, turning points, shell ratios, regimes, all subset entropies, all 372 facet margins, exact UV coefficients, and orbit minima.
- `validation.json`: numerical precision checks, multiple-start roots, saturation windows, configurations, minima, and mutual-information transition.
- `configurations.json`: intervals, subset ordering, and facet orbit sizes.
- `vaidya_entropy_cone.png` / `.pdf`: compact physics figure.
- `vaidya_diagnostics.png` / `.pdf`: supplemental figure with facet orbit envelopes and roundoff-scale minimum margins.

From this directory, run:

```sh
python3 -m pip install -r requirements.txt
python3 vaidya_strips.py
```

Dependencies: NumPy, SciPy, Matplotlib. The script resolves the exact facet data at `data/facets_n5_all372.json` relative to its own location. It can also be invoked by a script path from another working directory. The numerical calculation took approximately eight seconds before plotting in this environment.

## Geometry and normalization

\[
ds^2=\frac{1}{z^2}\left[-(1-\Theta(v)z^3)dv^2-2dv\,dz+dx^2+dy^2\right].
\]

The AdS radius and final horizon position are one. Reported entropies use

\[
\mathcal A\equiv A/L_y=4G_N S/L_y.
\]

The true AdS4 strip functional is

\[
\mathcal A=\int dx\,z^{-2}\sqrt{1-fv'^2-2v'z'}.
\]

This is the weighted area functional. It is not the geodesic length in the unrescaled three-dimensional metric.

The primary equations are Liu–Suh, [arXiv:1311.1200](https://arxiv.org/abs/1311.1200), section III.A, equations (3.3)–(3.23). Facet completeness is from Hernández-Cuenca, [arXiv:1903.09148](https://arxiv.org/abs/1903.09148), Table 1. Exact coefficients were obtained from the [author-maintained HEC database](https://github.com/SergioHC95/Holographic-Entropy-Cone).

## Configurations and times

1. Contiguous strips have widths $0.10,0.14,0.18,0.16,0.12$, spanning $0\le x\le0.70$.
2. Disjoint strips have width $0.12$ and gap $0.0741$, spanning $0\le x\le0.8964$.

For each configuration we compute the first three, first four, and all five parties. Contiguous selected intervals are merged before the endpoint matching calculation. Twenty-three different interval widths suffice for all these entropy vectors.

The time set contains a uniform 401-point grid from −0.025 to 0.7723594859, the quench time, every static-admissibility time, and extra points in the narrow competing-branch saturation windows. Deduplication yields 551 times.

## Extremal branches

Write $r=z_c/z_t$, where $z_t$ is the vacuum-side turning point and $z_c$ the shell crossing. For each $r\in[0,1]$, a scalar equation fixes $z_t$ from the strip width. The parameterized family supplies $t(r)$ and the area. The code brackets **all intersections visible in the parameterized scan**, refines roots, and compares their areas with the admissible static solution. It includes the vacuum branch for $t\le0$. Up to three competing candidates occur: two shell-crossing extremal surfaces and one static surface.

The family is **not assumed monotone in boundary time**. Even these short strips possess very narrow folds close to saturation. At the largest width, the static branch becomes admissible at approximately 0.7286410245, and the shell family persists to approximately 0.72948. The smallest-area branch is selected. Refined extremal-curve peaks are included in the bracketing grid.

### Why deep, radially nonmonotone branches are excluded here

The following argument concerns the standard reflection-symmetric shell-crossing strip family, with one regular turning point and the null-shell matching conditions. It is not a replacement for the general symmetry/HRT theorem.

Reflection in $x$ is a separate restriction of this numerical solver from
the transverse $y$-translation theorem. The computation uses the analytic
matching formula and does not independently establish common-slice
existence or symmetry inheritance.

Define

\[
P(z)=(z_t/z)^4-1,\qquad
H(z)=(1-z^3)P(z)+E^2,\qquad
E^2=\frac{z_c^6}{4}P(z_c).
\]

If the outgoing surface initially moves deeper after the shell, it must turn back before reaching the boundary. But

\[
H'(z)=z^{-5}\left[z_t^4(z^3-4)+3z^7\right]>0
\quad\hbox{for }z^3\ge4.
\]

A crossing with $z_c^3\ge4$ cannot produce that return. A monotone outgoing crossing has $z_c^3\le2$ directly from shell refraction. Thus every returning solution in this family has $z_c^3<4$.

The degenerate endpoint $z_c=z_t$ has $E=0$ and is handled by the static
limit. A turning point behind the final horizon cannot connect to the
boundary in this limit, because $H=(1-z^3)P$ is negative on part of the
required path.

For $0<z\le z_c$, $z^3P(z)=z_t^4/z-z^3$ is decreasing. Therefore

\[
\frac{E^2}{z^3P(z)}\le\frac{z_c^3}{4}<1,
\qquad H(z)<P(z).
\]

The width is at least the vacuum integral, with any radial return adding positive width:

\[
\frac w2\ge\int_0^{z_t}\frac{dz}{\sqrt{P(z)}}
=c z_t,\qquad
c=\frac{\sqrt\pi\,\Gamma(3/4)}{\Gamma(1/4)}.
\]

Consequently $z_t\le w/(2c)$. Here the maximum bound is 0.7481595009, and $h(z)\ge0.5812232268$. All these returning solutions lie outside the final horizon. Radial nonmonotonicity would require $z_c^3>2$, contradicting this bound.

At fixed $r$, the width is strictly increasing in $z_t$ in this exterior domain. This follows by differentiating the dimensionless denominator

\[
Q(u)=1-a u^3+\frac{a^2r^2(1-r^4)u^4}{4(1-u^4)},
\quad a=z_t^3,\quad 0\le u\le r,
\]

for which $\partial_aQ<0$. Thus the scalar width solve has a unique root for each $r$; folds in the boundary time still require the explicit minimum-area comparison described above.

## Cutoff and facet evaluation

The numerical check uses a common radial truncation at ε=10⁻⁴. Each interval area is retained as

\[
\mathcal A_\epsilon(w,t)=\frac2\epsilon+a_\epsilon(w,t).
\]

The small-$z$ vacuum and shell integrals below the cutoff are explicitly removed, retaining the finite-cutoff corrections. Candidate surfaces are anchored at the conformal boundary, their areas are truncated at this radial cutoff, and the smallest truncated candidate area is selected. This is not a new extremization with anchor points placed on the cutoff surface. The figures alone use vacuum-subtracted evolution.

The analytic full-cone theorem uses an exact endpoint-counterterm prescription (or compatible exact finite interval data). Its regulator should not be identified with these literal truncated areas at nonzero ε. The numerical calculation is a separate finite-cutoff consistency check; it does not prove equality between the two regulator prescriptions.

For a union with $m_I$ merged connected boundary intervals, every Catalan endpoint matching has the same divergence $2m_I/\epsilon$. A dynamic program finds the least finite matching cost. For an exact integer facet vector $c_I$, the code evaluates

\[
\sum_I c_I\mathcal A_\epsilon(I)
=\frac2\epsilon\underbrace{\sum_Ic_Im_I}_{\text{exact integer}}
+\sum_Ic_Ia_\epsilon(I).
\]

This avoids subtracting large UV-divergent floating-point numbers. The integer UV coefficient is nonnegative in every checked case. It is retained when positive; it is canceled exactly when zero. The cone test therefore uses regulated entropies, not the vacuum-subtracted plot values.

## Validation

- Maximum width-equation residual: approximately $1.2\times10^{-14}$.
- Maximum time-equation residual: approximately $6\times10^{-16}$.
- Gaussian quadrature refinement from 160 to 320 nodes on 48 representative shell points: differences below $9\times10^{-15}$ in width/time and $9\times10^{-16}$ in finite area at the same solution parameters.
- Independent adaptive quadrature of the static-admissibility time agrees within $1.3\times10^{-15}$.
- Forty-eight independent bounded two-variable solves, using four different initial shell ratios at each representative time/width, converge with no failed roots. The largest difference from continuation is below $1.6\times10^{-12}$ in the shell ratio and $6\times10^{-13}$ in the turning point.
- All 372 facet rows were evaluated at every stored time, for all six configurations. Saturated inequalities generate residuals at the $10^{-14}$ level; the diagnostic figure displays these without hiding their signs.

There is no interval-arithmetic proof that every possible fold has been resolved. The full-cone theorem must supply the general analytic result; these data are a reproducible numerical check and physical illustration.

## Suggested figure caption

**Thin-shell thermalization in AdS4.** (a) Vacuum-subtracted strip entropies, in units $4G_NS/L_y$, for several widths with final horizon radius one. (b) Mutual information between neighboring width-0.12 strips separated by 0.0741; the connected phase disappears at $t\simeq0.23360$. At 551 sampled times, all 372 five-party RT facet rows were checked for contiguous and disjoint three-, four-, and five-strip configurations, with no negative margin larger than floating-point roundoff. The largest observed negative residual was $7.11\times10^{-14}$. Competing shell-crossing and static extremal branches were compared, and a common radial cutoff was retained in the cone tests. These numerical samples illustrate, rather than prove, all-time cone containment.
