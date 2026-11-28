# All cone angles: the relaxed common-class functional is in the RT cone

Research theorem, 28 November 2026. This extends the mathematical functional result in `SMALL_ANGLE_FUNCTIONAL_THEOREM.md`. The word **relaxed** is essential: lengths are infima over regular embedded representatives, so a cone-hitting limiting configuration can contribute even when no regular geodesic in that class attains the infimum.

## Statement

Let `h_1,...,h_m` be finitely many marked hyperbolic cone metrics on a common finite-type bordered surface, with finitely many cone points and unchanged hyperbolic collars at physical conformal infinity. Cone angles may be arbitrary positive numbers. Let `w_j>0`.

Use particle-punctured isotopy classes and filled homology, with no curve endpoints on the particles. For a reduced embedded cut class `c`, define

\[
L_j(c)=\inf\{L_j^{\rm ren}(\gamma):
\gamma\text{ is a regular proper embedded representative of }c\}.
\tag{1}
\]

Representatives are geodesic rays in sufficiently small unchanged endpoint charts; this defines the usual renormalized length by subtracting their terminal logarithms. It is equivalent to allowing rectifiable proper representatives with the same ideal asymptotics and then straightening their tails. The finite physical endpoint set and class topology are fixed. No assumption of geodesic attainment is made.

For every physical party union put

\[
F_{\rm rel}^{\rm ren}(A)=
\inf_{c\in\mathcal C(A)}\sum_jw_jL_j(c).
\tag{2}
\]

Then there is a finite graph-cone vector `G_rel` and endpoint constants `b_p(A)=[p in boundary A]` such that

\[
F_{\rm rel}^{\epsilon}(A)
=G_{\rm rel}(A)+\sum_p t_p(\epsilon)b_p(A),
\tag{3}
\]

where `t_p(epsilon)=(sum_j w_j)log(1/epsilon)+O(1)`. In particular, for sufficiently small cutoff, the relaxed functional lies in the full static RT cone, for any finite number of parties.

Here `F_rel^epsilon` means restoring the fixed terminal endpoint logarithms to (2), rather than deleting arbitrary excursions of nongeodesic comparison curves into unrelated endpoint cutoff disks. This is the same regulator convention as the actual geodesic functional whenever regular minimizing geodesics exist.

## Proof

For each particle and metric, perform Lemma 2 of `SMALL_ANGLE_FUNCTIONAL_THEOREM.md` inside a radius `a>0`. There is no small-angle hypothesis in that lemma. Denote the resulting smooth complete negatively curved metric, with an auxiliary hyperbolic funnel at each particle, by `hhat_{j,a}`. Its marking to the original punctured surface is `1`-Lipschitz and is the identity outside the modified collars.

Auxiliary funnel colors are free, and there are no cut endpoints there. Particle-peripheral closed components can be deleted from a minimizing cut; the same deletion is allowed by filled homology downstairs. Thus the common class sets agree after reduction.

Let `Lhat_{j,a}(c)` be the renormalized length of the unique geodesic representative in the smooth replacement metric. The domination map gives

\[
\widehat L_{j,a}(c)\ge L_j(c).
\tag{4}
\]

Conversely, fix `c` and `delta>0`, and choose a regular proper embedded representative `gamma_j` with

\[
L_j^{\rm ren}(\gamma_j)\le L_j(c)+\delta.
\]

Its nonterminal part is compact and disjoint from the finite particle set. Therefore it misses every sufficiently small replacement collar. For those `a`, its full renormalized length is unchanged, and geodesic minimization in the replacement metric gives

\[
\widehat L_{j,a}(c)\le L_j(c)+\delta.
\]

Combining this with (4),

\[
\widehat L_{j,a}(c)\longrightarrow L_j(c)
\qquad(a\downarrow0).
\tag{5}
\]

No monotonicity of the replacement family and no uniformity in `c` is needed.

Put

\[
\widehat F_a(A)=\inf_c\sum_jw_j\widehat L_{j,a}(c).
\]

Equation (4) gives `Fhat_a(A)>=F_rel(A)`. To obtain the opposite limiting inequality, choose one common class `c_delta` whose weighted old infimum is within `delta` of (2). The number of metrics is finite, so (5) applied to this **one** class yields

\[
\limsup_{a\downarrow0}\widehat F_a(A)
\le\sum_jw_jL_j(c_\delta)
\le F_{\rm rel}(A)+\delta.
\]

Let `delta` decrease to zero. Thus

\[
\widehat F_a(A)\longrightarrow F_{\rm rel}(A).
\tag{6}
\]

The vector has finitely many party-union components, so convergence holds as a finite-dimensional vector.

Every smooth replacement satisfies the negative-curvature Crofton and finite-arrangement theorem in Sections 5--7 of the small-angle note. The endpoint subtraction families can be chosen once, independently of `a`, because all replacement metrics agree on the physical hyperbolic collars. Hence

\[
\widehat F_a^\epsilon=G_a+\sum_p t_p(\epsilon)b_p,
\qquad G_a\in\mathcal C_N^{\rm RT},
\]

with the **same** endpoint terms for every `a`. Equation (6) gives a finite vector limit `G_a -> G_rel`. The fixed-party graph entropy cone is closed, so `G_rel` is in that cone. Restoring the nonnegative endpoint edges at sufficiently small cutoff proves (3). `QED`

## Attainment and physical meaning

For cone angles below `pi`, the cone-surface geodesic realization theorem and uniform collar lemma make every relevant reduced simple length attained by a regular geodesic. The preceding relaxed theorem then agrees with the exact finite-graph theorem in `SMALL_ANGLE_FUNCTIONAL_THEOREM.md`; no relaxation remains.

For angles at least `pi`, some punctured simple arc classes have cone-hitting infima. If a proposed HRT prescription minimizes only over regular spacelike extremal curves, its class set or class costs may differ from (1)--(2). Even in a setup where `F_rel<=F_HRT` has separately been established, that comparison alone does not transfer holographic entropy inequalities, because those inequalities contain coefficients of both signs. No such comparison is asserted for an arbitrary physical particle spacetime. In fact, `TILTED_CONE_COUNTEREXAMPLE.md` gives a static spinless positive-cone example with physical HRT strictly smaller than this relaxed functional. A valid Lorentzian identification requires a different argument in that regime.

Likewise, this argument does not supply a fixed-pair common-class description across particle collisions. It proves a mathematical result for marked cone metrics, not an identification of those metrics with all HRT candidates in a collision spacetime.

The local proposition in `COMMON_SIDE_ROUNDING.md` excludes a finite compatible cone-hitting star from being a global minimum whenever the required particle-side resolutions are admissible. It does not establish the compactness and common-star compatibility that would remove the above physical qualification in full generality.

## References and dependence

The explicit funnel replacement, local endpoint subtraction and positive graph construction are proved in the companion note. The analytic scattering/Crofton input is T. Lefeuvre, [arXiv:1805.05155](https://arxiv.org/abs/1805.05155), Lemma 3.1 and Remark 3.1. Closedness and finite graph realization of the fixed-party holographic entropy cone are the standard graph-cone theorem of Bao et al., [arXiv:1505.07839](https://arxiv.org/abs/1505.07839).
