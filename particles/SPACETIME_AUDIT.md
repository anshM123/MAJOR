# Particle and collision spacetime audit

## Decisive corrections to the proposed blanket extension

1. A spinless massive particle is a timelike cone axis. Its two chiral
   elliptic holonomies have equal angles. Barbot--Bonsante--Schlenker,
   arXiv:1202.5753v2, Proposition 4.20, proves that the associated left and
   right hyperbolic cone metrics have the same cone angle at each particle.
2. With collisions the same paper associates a graph of interaction-free
   slices, a pair of cone metrics to each slice, and compatible double
   surgeries at collisions. It does not give one unchanged pair of marked
   cone surfaces for the entire spacetime. See Introduction 1.6--1.7,
   Proposition 5.9, and Sections 5.3--5.5. Section 5.3 even gives examples
   where a slice is good on one side of a collision and the required
   transverse-field condition fails on the other side.
3. Matschull's explicit black-hole formation solution,
   arXiv:gr-qc/9809087, is a collision of two *massless* particles. Their
   holonomies are parabolic, not ordinary massive elliptic cone holonomies.
4. Holst--Matschull, arXiv:gr-qc/9905030, uses nonzero impact parameter to
   produce a rotating configuration with closed timelike curves, then
   interprets the CTC region as singular/excised. A proof for a filled
   globally hyperbolic cone manifold cannot import that spacetime without
   specifying the causal excision and new homology boundary.

## An exact local obstruction for spinning point quotients

In global AdS3 of radius one,

\[
 ds^2=-(1+r^2)dt^2+\frac{dr^2}{1+r^2}+r^2d\varphi^2.
\]

Let the quotient identify along a helical generator
\(\xi=a\partial_\varphi+b\partial_t\). Its norm is

\[
 \|\xi\|^2=(a^2-b^2)r^2-b^2.
\]

For nonzero spin/time shift \(b\ne0\), it is timelike near the axis. Its
identified orbits are closed timelike curves. If \(a^2>b^2\), the CTC region
includes \(r<|b|/\sqrt{a^2-b^2}\); if \(a^2\le b^2\), the generator does not
become spacelike at large radius. Removing the causally bad region is a
different geometry with an internal boundary, not a filled timelike cone.

The group-theoretic version is equally simple. An AdS point \(x\) is fixed
by \((g_L,g_R)\) iff \(g_L=xg_Rx^{-1}\). A timelike fixed axis thus requires
the elliptic factors to be conjugate, hence to have equal rotation angles.
Unequal chiral angles do not describe the spinless fixed-axis local model.

This does not assert that every physically regular spinning extended source
violates NEC. It shows that arbitrary spinning point identifications cannot
be assumed to be filled globally hyperbolic cone spacetimes. Li,
arXiv:2403.05524, Section 2 and equations (8)--(9), explicitly distinguishes
the CTC region and the causally admissible quotient region.

For a spinless positive-mass worldline, the usual distributional stress
tensor has \(T_{\mu\nu}k^\mu k^\nu=m\int(u\cdot k)^2\delta\,d\tau\ge0\).
For a spinning point source, derivative-of-delta terms and the chosen
causal domain need a separate treatment; positive mass alone does not prove
the proposed combined hypothesis.

## Why the existing Theorem H does not transfer verbatim

The vacuum proof uses a smooth invisible-domain quotient, an injective
developing description, a boundary conjugacy of discrete purely hyperbolic
groups, and a global acausal circle whose chords can be simultaneously
placed on a pleated Cauchy slice. A generic cone holonomy is elliptic of
irrational angle and its cyclic image is nondiscrete. A cone manifold is
constructed by wedge gluing/developing maps, not by the same smooth proper
quotient of the whole invisible domain. Collisions also change the
interaction-free stereographic data by surgery.

The local boundary-to-boundary geodesic identity in the AdS developing chart
remains valid for an individual developed segment. What is unproved is
global admissibility and the identification of its branches with one
common family of filled-surface cut classes through a collision.

The conformal double from Proposition S+ is smooth only away from the
particle axes. Doubling at infinity does not smooth those interior axes.
Bernal--Sanchez's smooth-spacetime theorem therefore cannot be applied
directly to that singular double. The finite-orbifold route in
`FINITE_ORBIFOLD_THEOREM.md` repairs this by performing the entire smoothing
argument on a genuine smooth finite cover and averaging its temporal
function before descent.

## Primary sources inspected

- Barbot--Bonsante--Schlenker, *Collisions of particles in locally AdS
  spacetimes I*, arXiv:1010.3602 (rewritten first part of 0905.1823).
- Barbot--Bonsante--Schlenker, *Collisions of particles in locally AdS
  spacetimes II*, arXiv:1202.5753v2. In particular Proposition 4.20,
  Proposition 5.9, Section 5.3, and Appendix A.
- Matschull, *Black Hole Creation in 2+1 Dimensions*,
  arXiv:gr-qc/9809087.
- Holst--Matschull, *The Anti-de Sitter Gott Universe: A Rotating BTZ
  Wormhole*, arXiv:gr-qc/9905030.
- Li, *Spinning Particle Geometries in AdS3/CFT2*, arXiv:2403.05524.
- Bernal--Sanchez, arXiv:gr-qc/0512095, Theorem 1.2.

No unconditional generic-collision entropy-cone theorem is asserted here.
