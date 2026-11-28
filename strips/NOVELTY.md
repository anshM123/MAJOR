# Novelty audit: covariant strip entropies in the full RT cone

Checked against primary literature on 2026-11-28. This is a scoped literature audit, not a proof of absence of prior work.

## Main finding

I found no primary source proving the proposed exact statement for arbitrary strip widths and arbitrary time dependence in translationally invariant higher-dimensional bulk geometries. The closest rigorous predecessors are the three-dimensional kinematic-space argument of Czech–Dong, the leading extensive membrane-limit argument of Bao–Mezei, and the conditional graph constructions of Grado-White et al. and Zhao. A rigorous symmetry reduction together with a direct positive graph construction would go beyond the explicit scope of these sources.

The strongest defensible novelty wording is:

> We give an exact higher-dimensional strip-sector realization of covariant entropies by static graph min-cuts, valid at arbitrary strip widths and boundary times under the stated geometric and maximin hypotheses.

This distinguishes the result from a numerical check, a membrane-limit result, and a proof confined to contraction-certifiable inequalities. Whether the work warrants PRL is a separate editorial judgment; the search does not establish that.

## Primary-source comparison

| Source | What is actually established | Relation to the proposed result |
|---|---|---|
| [Czech–Dong, arXiv:1905.03787](https://arxiv.org/abs/1905.03787), JHEP 10 (2019) 177 | In two boundary / three bulk dimensions, extends inequalities with contraction-map proofs to time dependence. The discussion explicitly says the ordering/topological argument does not obviously generalize to higher dimensions. | The full text contains no discussion of transverse translations, higher-dimensional strips, or dimensional reduction. It is the immediate kinematic predecessor, but its explicit conclusion is contraction-based. Its multiply connected extension requires care highlighted by 2407.07165. |
| [Bao–Mezei, arXiv:1811.00019](https://arxiv.org/abs/1811.00019) | Uses the membrane description to transfer static cut-and-paste inequalities to the leading extensive entropy of large regions at late times in black-brane collapse. | It already connects higher-dimensional thermalization to static entropy inequalities, but only in the membrane scaling regime. Section 4 expressly says subleading higher-derivative terms obstruct that proof method. An exact all-width/all-time strip theorem is stronger. |
| [Caginalp, arXiv:1905.00544](https://arxiv.org/abs/1905.00544), PRD 101 (2020) 026010 | Numerical tests of the five-region inequalities in AdS3–Vaidya, including configurations on varying boundary slices, and NEC-violating examples. | It is not an AdS4 theorem, nor an all-party full-cone proof. An AdS4 figure should acknowledge this numerical predecessor. |
| [Grado-White–Grimaldi–Headrick–Hubeny, arXiv:2407.07165v2](https://arxiv.org/html/2407.07165v2), revised 21 May 2025 | Studies the 2+1-dimensional problem; identifies the nontrivial topology issue and proves a cyclic-topology extension; extensive numerical tests. Appendix B proves dimension-independent facts about degenerate HRT surfaces under a Jacobi-eigenvalue genericity hypothesis. | This is directly relevant to symmetry inheritance, but its explicit symmetry corollary is for finite groups. It does not by itself settle continuous transverse translations or remove nongeneric degeneracies. |
| [Grado-White–Grimaldi–Headrick–Hubeny, arXiv:2502.09894v2](https://arxiv.org/abs/2502.09894), JHEP 05 (2025) 104 | Proves minimax properties and a graph construction conditional on a cooperating collection of time-sheets. Gives partial proofs and possible failure scenarios. | It does not prove unconditional full RT/HRT cone equality. Its Theorem 4.1 is explicitly conditional; Conjecture 4.1 is the cooperating assumption. |
| [Zhao, arXiv:2602.04888v2](https://arxiv.org/html/2602.04888v2), *Graph models for covariant holographic entropy I*, revised 27 April 2026 | Conditional graph/no-short-cut theorem assuming nonempty exposed regions for HRT pairs, possibly after timelike-cluster reduction. | Search results still reproduce the earlier universal-resolution abstract. The current theorem is conditional. The abstract, Theorem 1.1, and main construction must be used instead of the stale search snippets. No exact translational-strip specialization was found. |

## Symmetry source: restrictions that must not be omitted

Appendix B of 2407.07165v2 is independent of spacetime dimension. It assumes that for every HRT surface, the principal eigenvalues of the two null diagonal Jacobi operators are not both zero. The authors explain that this excludes continuous HRT families and give a product-metric counterexample when the assumption fails.

Under that hypothesis:

1. Lemma B.1 states that any maximin slice for a boundary region contains all its HRT surfaces.
2. Theorem B.2 reorganizes a **finite** set of HRT surfaces into a nested family with the same union; outermost and innermost entanglement wedges are respectively the union and intersection.
3. Corollary B.3 supplies an invariant HRT surface when a **finite** group of bulk isometries preserves the boundary region.
4. Corollary B.4 says two HRT surfaces for the same region intersect only in coincident connected components in a smooth bulk.

These results support a symmetry argument, but quoting B.3 as an unconditional compact- or continuous-group theorem would be incorrect. A finite-subgroup limiting proof needs compactness and regularity control; a maximin/minimax sandwich can avoid that particular route if proved in the full metric class.

## Regulator warning from newly released literature

[Ju–Czech, arXiv:2611.21726](https://arxiv.org/abs/2611.21726), submitted 16 November 2026, disproves the broad assertion that all holographic entropy inequalities apart from the standard trivial exceptions are superbalanced / free of area-law UV divergences. This does not obstruct a direct graph construction. It does rule out relying on blanket UV-cancellation claims for unknown full-cone facets.

For the strip theorem, retain a common regulator or explicitly prove the positive Bell-lift statement needed by the graph construction. A minimally subtracted finite entropy vector is not automatically an RT-cone point.

## Practical AdS4 thin-shell equations

Primary source: [Liu–Suh, arXiv:1311.1200](https://arxiv.org/abs/1311.1200), PRD 89 (2014) 066012, section III.A, equations (3.3)–(3.23). The formulae below are their $n=d-1=2$ specialization, with AdS radius $L$, transverse coordinate length $L_y$, strip width $w$, and final horizon radius $z_h$.

\[
ds^2=\frac{L^2}{z^2}\left[-f(v,z)dv^2-2dv\,dz+dx^2+dy^2\right],\qquad
f(v,z)=1-\Theta(v)(z/z_h)^3.
\]

For an $x$-parametrized symmetric extremal surface,

\[
\frac{A}{L_y L^2}=\int_{-w/2}^{w/2}\frac{dx}{z^2}
\sqrt{1-fv'^2-2v'z'}.
\]

Let $z_t$ denote the vacuum-side turning point, and $z_c$ the null-shell crossing. Put

\[
h(z)=1-(z/z_h)^3,\qquad g(z)=1-h(z),
\]

\[
E=-\frac{g(z_c)}2\sqrt{(z_t/z_c)^4-1},\qquad
H(z)=h(z)\left[(z_t/z)^4-1\right]+E^2.
\]

For the branch on which $z$ decreases monotonically from the shell to the boundary,

\[
\frac w2=\int_{z_c}^{z_t}\frac{dz}{\sqrt{(z_t/z)^4-1}}
+\int_0^{z_c}\frac{dz}{\sqrt{H(z)}},
\]

\[
t=\int_0^{z_c}\frac{dz}{h(z)}\left(1+\frac{E}{\sqrt{H(z)}}\right),
\]

\[
\frac{A}{2L_yL^2}=\frac1{z_t}\int_{z_c/z_t}^{1}
\frac{du}{u^2\sqrt{1-u^4}}
+z_t^2\int_\epsilon^{z_c}\frac{dz}{z^4\sqrt{H(z)}}.
\]

The apparent horizon pole in the time integral is removable when $E<0$. The shell refraction rules are

\[
v'_+=v'_-,\qquad z'_+=\left(1-\frac12g(z_c)\right)z'_-,
\qquad z'_-=-\sqrt{(z_t/z_c)^4-1}.
\]

The source warns immediately after (3.19) that nonmonotone radial branches require amended quadratures. In Schwarzschild AdS4, competing strip branches and discontinuous saturation can occur. A numerical check must minimize over all relevant branches; a single shooting solution is not sufficient. The fully vacuum and fully black-brane branches must also be included when admissible.

## Complete five-party facet data supplied

[Hernández-Cuenca, arXiv:1903.09148v3](https://arxiv.org/abs/1903.09148), PRD 100 (2019) 026004, Table 1, establishes that the five-party RT cone has 372 facets in eight orbits under $S_6$, including the purifier. Paper-table orbit sizes are $15,20,45,72,10,60,60,90$.

The primary author-maintained [Holographic-Entropy-Cone database](https://github.com/SergioHC95/Holographic-Entropy-Cone) currently supplies exact representatives at `data/n=5/facets.json`. The catalogues are explicitly complete for $n\le5$ and partial for $n=6$.

Downloaded files:

- `facets_n5_upstream.json`: the eight current exact integer representative rows.
- `facets_n5_all372.json`: all distinct $S_6$ images after purification, with explicit subset ordering and orbit index lists.
- `expand_facets.py`: reproducible expansion and orbit-count verification.

Rows use the convention $\sum_I c_I S(I)\ge0$. Subsets are in cardinality-then-lexicographic order $A,B,C,D,E,AB,AC,\ldots,ABCDE$. Upstream row order differs from the paper table; its orbit sizes are $20,45,15,90,72,60,60,10$. For a permutation that maps a subsystem to one containing the purifier, replace it by its complement in all six labels. Exactly 372 rows are obtained after deduplication.

Numerics on three/four strip parties can either use the corresponding smaller catalogues or specialize the five-party inequalities by setting the extra parties empty. A finite grid checks that grid; it is not an all-time proof. A plot should state its grid, branch treatment, cutoff, and minimum residual tolerance.

## Search coverage

Searched both web indices, including stronger current searches, for combinations of “holographic entropy cone”, “strips”, “Vaidya”, “quench”, “translational”, and “dimensional reduction”. Read the primary full texts of 1905.03787, 1811.00019, 1905.00544, 1903.09148, 1211.3494, 1311.1200, and 2502.09894; inspected the latest primary HTML of 2407.07165v2 and 2602.04888v2. No exact prior theorem matching the proposed strip-sector scope emerged. This is evidence supporting a novelty claim, not a guarantee that no unpublished or unindexed antecedent exists.
