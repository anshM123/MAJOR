# Covariant holographic entropy inequalities for vacuum rotating wormholes in 2+1 dimensions

**Authors:** Ansh Mishra and Aryan Senthilkumar

## What this repository contains

A proof that the covariant (Hubeny–Rangamani–Takayanagi) entanglement entropies of vacuum solutions of
2+1-dimensional gravity with negative cosmological constant — rotating BTZ black holes, rotating multi-boundary
wormholes and wormholes with handles, of arbitrary topology — obey the holographic entropy inequalities known
for the static Ryu–Takayanagi formula, together with the supporting mathematics, numerical tests and a draft
manuscript.

| File | Content |
|---|---|
| `MANUSCRIPT.md` | Draft paper: statements, ideas of proofs, status of every claim |
| `PROOFS.md` | Complete proofs |
| `SMOOTH_SLICES.md` | G1 regularity upgrade: smooth uniformly spacelike homology slices through complete HRT cuts, including boundary endpoints |
| `smooth_slices/` | Exact symbolic checks of the rotating strip and conformal endpoint formulas |
| `PROBLEM.md` | The question and why it matters |
| `LITERATURE.md` | Prior work and what was and was not known |
| `APPROACHES.md` | Proof strategies tried, including the ones that failed |
| `DERIVATIONS.md` | Working calculations |
| `EXPERIMENTS.md` | Numerical tests, with commands and results |
| `COUNTEREXAMPLES.md` | Failed conjectures, controls, invalid test set-ups that were caught |
| `OPEN_GAPS.md` | What is not proved, standing assumptions, where an error could hide |
| `code/` | Scripts for all numerical tests |
| `archive_qi/` | Earlier, unrelated exploratory work (not part of the result) |

## Status

The status of every statement is tabulated in `PROOFS.md` §7 and `MANUSCRIPT.md` §9; please read those before
relying on any claim. The work was developed with AI assistance (Claude) and has been checked by automated
adversarial review and numerical falsification tests, not yet by independent human experts, and it is not
formalised in a proof assistant.

Revision 8 adds Proposition S+ in `PROOFS.md` §6.6b. It removes the smooth-slice restriction to
whole-boundary regions by doubling across the timelike conformal boundary, extending the doubled compact
cut, and restricting back with homology preserved. The detailed proof is in `SMOOTH_SLICES.md`; the other
standing assumptions and open questions remain listed in `OPEN_GAPS.md`.

## Reproducing the numerical tests

```bash
pip install numpy scipy
cd code
python3 contraction.py
python3 annulus_test.py 1 3000 25.0
python3 lemmaA_test.py 0 400
python3 separation_check.py 0
python3 spacelike_check.py
python3 fourholed_test.py 3 600 7 fn
python3 -u fourholed_certified.py 3 40
```
