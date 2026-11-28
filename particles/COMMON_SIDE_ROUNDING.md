# Common-side rounding at a positive cone

Research proposition, 28 November 2026. This is a local weighted-length statement. It does not by itself prove compactness or global attainment of a minimum over marked curve classes.

## Proposition

Let `h_1,...,h_m` be cone metrics at a commonly marked point `q`, with total angles `0<theta_j<2 pi`, and let `w_j>0`.

Suppose, in each metric, a finite frontier is made of `2r>=2` regular radial geodesic germs ending at `q`. The germs are indexed in the same cyclic order in every metric, and the sectors have the same alternating binary region colors. Assume distinct germs, or a limiting configuration admitting arbitrarily small compatible disjoint perturbations at arbitrarily small extra length.

Then there is a **single common topological resolution** of the frontier, avoiding `q`, whose metricwise representatives have weighted total length strictly less than that of the cone-hitting frontiers. Outside an arbitrarily prescribed cone chart the representatives are unchanged. Their physical boundary trace is unchanged, with `q` filled for homology.

The resolution need not preserve a preassigned punctured isotopy class: it selects one of the particle-side choices. It is an admissible competitor when minimizing over all particle-marked cut classes with filled homology.

## Proof

Fix one metric `h_j0`. The `2r` adjacent sector angles add up to `theta_j0<2 pi`, so at least one adjacent sector has angle `alpha<pi`. Choose its bounding pair of germs.

There are two canonical noncrossing resolutions of an alternating binary star: include a small disk around the point in the selected region and pair the rays around the complementary sectors, or exclude the disk and pair around the selected sectors. Choose the resolution containing the chosen adjacent pair. This fixes a common pairing in every metric, without asking the same sector to be short in those other metrics.

At equal small radius `rho` on the chosen two germs of metric `j0`, their original broken path through the cone has length `2 rho`. The hyperbolic chord within the sector has length

\[
d_\alpha(\rho)=\operatorname{arcosh}
\big(\cosh^2\rho-\sinh^2\rho\cos\alpha\big)
<2\rho.
\]

Its strict saving is

\[
\delta=2\rho-d_\alpha(\rho)>0,
\qquad
d_\alpha(\rho)=2\rho\sin(\alpha/2)+O(\rho^3).
\]

All other pairs in this metric can be rounded in their prescribed sectors at radius `epsilon`, with added length bounded by a constant times `epsilon`. Such a connector follows each radial germ to radius `epsilon`, runs along the radius-`epsilon` circle through the prescribed sector, and returns along the other germ. Its length excess over the cone-hitting path tends to zero. The same construction applies to **all** pairs in every other metric, regardless of whether the chosen sector there has angle above or below `pi`.

The prescribed pairings involve disjoint adjacent sectors, so the local connectors can be made embedded and mutually disjoint. The chosen chord lies in its own sector. The same pairing gives the same particle-marked topological resolution in all metrics.

Choose the remaining rounding radii sufficiently small that the sum of all weighted excesses is less than `w_j0 delta/2`. The total weighted length then decreases by at least `w_j0 delta/2`. Slight smoothing preserves the strict inequality. Including or excluding the small disk changes only the particle's region color, so filled homology and all physical boundary labels are preserved. `QED`

## Exact limitation

This proves that a finite cone-hitting common frontier with compatible local star data cannot be a global minimum when both particle-side resolutions are admissible. It does not prove that a minimizing sequence in different punctured classes converges to such a common finite star, nor that two independently tightened chiral representatives have compatible local pairing data. Those are precisely the additional statements needed before using this proposition in a general all-angle existence theorem or a branched-cover descent.

## Why arbitrary marked pairs need not have joint regular attainment

There is an elementary counterexample to the unrestricted statement that a common-class minimum for two cone metrics must have regular geodesic representatives in both metrics.

Take two copies of a hyperbolic cone of total angle `theta=3 pi/2`, with commonly marked interval endpoints whose angular separations are

\[
\beta_L=5\pi/4,\qquad \beta_R=\pi/8.
\]

The common punctured disk has exactly two embedded anchored arc classes. Their sector sweeps in the two metrics are

| Common class | Left sweep | Right sweep |
|---|---:|---:|
| First side | `5 pi/4` | `pi/8` |
| Other side | `pi/4` | `11 pi/8` |

For an anchored class with sweep `delta`, radial cutoff `r=R` gives

\[
L^{\rm ren}(\delta)=
\begin{cases}
2\log\sin(\delta/2),&0<\delta<\pi,\\
0,&\delta\ge\pi.
\end{cases}
\]

The first line follows from the hyperbolic distance formula after subtracting `2R`. In the second line, the infimum is the broken radial path through the cone: regular representatives can approach it arbitrarily closely, but a regular complete geodesic cannot sweep an angle at least `pi` in the developed cone chart.

Thus the first common class has a cone-hitting left infimum and a regular right geodesic; the second class has the reverse behavior. No class has regular geodesic representatives in both metrics. A positive weighted common-class infimum exists as a number, but it is not jointly regularly attained. Particle-peripheral closed components cannot help because they are removable with filled homology.

This is a counterexample for **arbitrary pairs of marked cone metrics**. It does not establish that these markings arise from a globally admissible particle spacetime with the proposed HRT/chiral-length identification. Its role is to show that an unrestricted mathematical joint-regular-attainment claim cannot be inserted into that identification without additional spacetime input.
