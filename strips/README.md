# Covariant strips and the static entropy cone

Start with `THEOREM.md`. It states the proved theorem, hypotheses, proof, and
scope in one place. The full technical proofs are in `SYMMETRY.md`,
`COMPACT_MAXIMIN.md`, and `FINITE_GRAPH.md`. The numerical example and literature
comparison are in `NUMERICS.md` and `NOVELTY.md`.
`UV_AVERAGING.md` supplies the exact trace-mismatch estimate and its
controlled conformal-boundary limit, completing the asymptotic comparison.

## What is proved

Under the explicit regular maximin, invariant-anchor, and disk-homology
hypotheses, every finite strip entropy vector has an exact positive graph
realization after specified local endpoint completion. All full RT-cone
inequalities hold for that completed vector, for any number of parties and
at arbitrary widths and times. Every endpoint-balanced inequality holds for
the finite renormalized data. The symmetry theorem compares against all
transversely nonsymmetric competitors.

The proof is not an unconditional HRT existence theorem from NEC alone, a
theorem about noninvariant anchors or arbitrary orbit topology, or a claim
that the bare renormalized vector is itself nonnegative. Literal radial
truncation is distinguished from exact local endpoint completion.

## Reproduce the checks

Run these from the repository root:

```sh
python3 strips/code/verify_finite_graph.py --output strips/data/finite_graph_validation.json
python3 strips/code/check_kk_ricci.py
python3 strips/data/expand_facets.py
python3 strips/numerics/vaidya_strips.py
```

The exact arrangement check uses the Python standard library. The Ricci
check uses SymPy. The Vaidya check uses NumPy, SciPy, and Matplotlib, listed
in `numerics/requirements.txt`. The facet expansion uses the standard
library and the author-maintained five-party representative data included
in `data/`.

The numerical script reads its own `numerics/data/facets_n5_all372.json` and
writes its figures and reports beside the script. The archived figure copies
are in `figures/`. The graph command above writes its validation report to
`data/finite_graph_validation.json`.

## Numerical coverage

Two strip arrangements, their 3/4/5-party prefixes, 551 times, and 372
five-party facet rows give 1,229,832 sampled checks. The worst residual is
about $-7.11\times10^{-14}$, with no violation beyond roundoff. The solver
compares competing shell-crossing and static candidates. See `NUMERICS.md`
for the reflection-family restriction, finite radial truncation, branch
coverage, and accuracy checks.

The exact graph test verifies 756 mincut/matching equalities for 3--8
endpoints, including zero weights. The Ricci test checks the general
Kaluza--Klein identity in nontrivial one- and two-fiber examples using exact
arithmetic.

## Provenance

Prepared against repository commit
`acda6f59a240537f521583630ff9e9fac71e145f`, on 28 November 2026. All changes
are additions under `strips/`; existing wormhole and full-cone proofs are
preserved. Primary references and the exact facet-data source are recorded
in the notes. No claim of priority is inferred merely from the absence of a
matching search result.
