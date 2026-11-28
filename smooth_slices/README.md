# Exact checks for Proposition S+

The analytic proof is in `../SMOOTH_SLICES.md` and `../PROOFS.md` §6.6b.
These scripts check local identities, not the global existence theorem.

Requirements: Python 3 and SymPy (verified with SymPy 1.14.0).

```bash
python3 smooth_slices/verify_rotating_strip.py
python3 smooth_slices/verify_boundary.py
```

Each script writes its corresponding `*_verification.json` next to the
script and exits nonzero if an exact identity fails. No numerical samples
or tolerances are used.

- `verify_rotating_strip.py`: 23 identities for the moving frame, determinant,
  holonomy equivariance, induced metric, pairwise separation, observer
  metric, normal, and a surrounding coordinate metric and its conformal
  compactification. Positivity is established analytically using
  `a,b > 0` and `|q| < 1` in the proof note.
- `verify_boundary.py`: 23 identities for the generic rank-one endpoint
  conic, reflected tangent norm, Poincaré endpoint curve, even BTZ conformal
  metric, conserved momenta, and conformal tangent norm.

The globally hyperbolic double, compact ribbon acausality, Cauchy
restriction, peripheral compactness, and homology transport are proved
in the text. Algebra checks do not certify those topological and causal
arguments. Review remains mathematical review, not proof-assistant
formalization.
