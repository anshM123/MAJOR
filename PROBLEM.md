# PROBLEM

## Chosen target

**Do the entanglement entropies computed by the covariant Hubeny–Rangamani–Takayanagi (HRT) formula obey the
same linear inequalities as those computed by the static Ryu–Takayanagi (RT) formula?**

Precisely, the sub-case attacked here:

> Let `M` be a vacuum (locally `AdS_3`) solution of 2+1-dimensional gravity with negative cosmological
> constant, of arbitrary topology: a rotating BTZ black hole, a rotating multi-boundary wormhole, or a wormhole
> with handles, with independent left- and right-moving moduli (so in general neither static nor
> time-symmetric). Let the parties be unions of intervals on an arbitrary boundary Cauchy slice and of whole
> asymptotic boundaries. Show that the HRT entropies satisfy every holographic entropy inequality that admits
> a contraction-map proof.

## Why this problem

* **Physical content.** The RT inequalities (strong subadditivity, monogamy of mutual information, and the
  infinite hierarchy beyond) are constraints that a quantum state must satisfy to have a classical geometric
  dual. They are proved for time-symmetric slices. Black-hole formation, evaporation and every genuinely
  dynamical process require the covariant formula, for which only strong subadditivity and monogamy are
  proved in general. Whether the "holographic entropy cone" is the same with and without time dependence is a
  basic open question about how spacetime emerges from entanglement (direction 5 of the brief).
* **It is sharply stated and open.** Grado-White, Grimaldi, Headrick and Hubeny (arXiv:2407.07165) found an
  error in the earlier 2+1-dimensional proof when the bulk is multiply connected and write: "The validity of
  the RT inequalities for multiply-connected 2+1 dimensional bulk spacetimes remains an open question." Their
  evidence for vacuum multi-boundary wormholes is numerical. The October 2026 review by Grimaldi and Headrick
  (arXiv:2610.33905) still lists the validity of higher inequalities for HRT as open.
* **It concerns real gravitational spacetimes**, not a toy model: these are exact solutions of Einstein's
  equations containing several rotating black holes behind a common interior.

## What it is not

* It is not a resolution of the black-hole information problem. It does not address evaporation, the Page
  curve, factorisation, or interior reconstruction (directions 1–4, 6 of the brief). After the audit
  (LITERATURE.md) no proposition in those directions was found that is both sharply posed and within reach of
  a rigorous proof here, apart from the secondary target below.
* It is restricted to vacuum solutions in 2+1 dimensions. Matter and higher dimensions are not covered
  (OPEN_GAPS.md).

## Secondary target (recorded, not yet attacked)

Gaussianity of the one-time-point SYK partition function `z` for `q > 2`, which Saad–Shenker–Stanford–Yao
(arXiv:2103.16754, footnote 6) "expect … but have not proven", together with a quantitative form of their
wormhole + half-wormhole decomposition of `z²`. Proposed route: `z` lies in a single Wiener chaos, so the
Nourdin–Peccati fourth-moment bound applies, and the half-wormhole error is a sum of contraction norms that is
bounded by the fourth cumulant. See APPROACHES.md.

## Self-assessed level (brief's scale)

Stated separately for the two parts, because their verification status differs.

* **Mathematical part** (Theorem 1 of PROOFS.md: contraction-map inequalities for the common-class functional
  `min_c ½(ℓ_L + ℓ_R)`; Theorem 3 for whole-boundary parties). Three automated referee passes, no error found.
  This is what the numerical experiments of arXiv:2407.07165 tested and left open as a proof.
* **Gravitational part** (Theorem H: the HRT entropy, defined with acausal bulk slices, equals that
  functional; hence Theorem 2). Proved in PROOFS.md §§6.3–6.7 without any admissibility hypothesis; status of
  the independent audit: see OPEN_GAPS.md G6. Proposition S+ (PROOFS §6.6b; SMOOTH_SLICES.md) now proves
  the identical entropy with smooth, uniformly spacelike completed Cauchy slices, including for interval
  regions. Its new input is a compact extension theorem applied after conformal doubling.

If both parts survive expert verification: **Level 4** (a new theorem settling a significant special case of
an open problem). Not Level 1–3: matter, higher dimensions, and inequalities without contraction maps for
interval regions remain open.
