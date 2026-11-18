# The parity-odd sector of the route-A′ functional `Ω^{A′}` (binding state `+−`, odd `J`)

*Started 18 Nov 2026 by the parity-odd-sector agent. Unrefereed, not committed. Updated incrementally; the hand-over is the last section.*
Code: `ads4/code/spinfunc/p1_*.py, p2_*.py, …`; caches `spinfunc/cache/p*` (git-ignored). Imports `l1/l4/l5/l6`, `fourpt/tttt/odd3.py`
read-only; nothing else modified.

Marks: **[certified]**, **[D]** derived, **[DV]** derived and checked (exact algebra / round-off), **[NE]** numerical evidence,
**[I]** inferred/assumed, **[M]** model.

## 0. Why this sector matters

For identical stress tensors and odd `J` the only conserved `⟨TT𝒪_{Δ,J}⟩` structure is parity-odd, and it couples only to `+−`
(`PARITY_ODD_HEAVY.md` §2.1). So every odd-`J` `+−` heavy state, including the binding state of Result D (`+−`, `J = 5`, `m ≈ 1.05Λ`;
`F/N_s ≈ 0.225 / 0.205 / 0.202` at `Δ_g = 100 / 300 / 1000`), is a parity-odd exchange. Its route-A′ action was so far *assumed* to be the
stripped flat `+−` density with a positive coefficient (H12′); the `Ĉ₂` crystal corrections of `H0_SPINNING_FUNCTIONAL.md` §22–§23 covered
even `J` only.

## 1. The covariant parity-odd conserved `⟨T T 𝒪_{Δ,J}⟩` (`p1_tto_odd.py`, `p1b_validate_odd3.py`) [DV]

**Construction.** l4/l5 conventions (Lorentzian lightcone components, CLS-D `V`, `H`, `X_ij = (x_i − x_j)²`). `ε(A,B,C,D,E)` := the
determinant of the 5×5 matrix of embedding components `(A⁺, A⁻, u, v, y)`: a real-coefficient polynomial, so the structure is real at real
spacelike Lorentzian configurations and its OPE coefficient is real in a unitary theory (*Lorentzian-real basis*). Basis (CPPR's 3d rule:
odd = `ε_ij` × even):
* `ε(Z₁,Z₂,X₁,X₂,X₅)X₁₂`, `ε(Z₁,Z₅,X₁,X₂,X₅)X₁₅`, `ε(Z₂,Z₅,X₁,X₂,X₅)X₂₅`, each times the even monomials `V₁^{a₁}V₂^{a₂}V₅^{a₅}H₁₂^{b₁₂}H₁₅^{b₁₅}H₂₅^{b₂₅}`
  with the remaining spins, over the even denominator times `(X₁₂X₁₅X₂₅)^{1/2}`;
* conservation of `T₁`, `T₂` by the `d = 3` Todorov operator (exact TPS differentiation, as l4), `1↔2` symmetry, 3d relations quotiented
  numerically. Exact `Δ`-dependent coefficient vector (real to round-off) from the null space.

**Counts** (`cache/p1_tto_odd.txt`; `Δ = 7.3, 23.7, 101.1, 1003.7`; null-space gap `≤ 3·10⁻¹⁴` vs `≥ 0.09`):

| `J` | odd monomials | independent in `d = 3` | conserved, `1↔2` symmetric | conserved, antisymmetric |
|---|---|---|---|---|
| 0 | 2 | 2 | 1 | 0 |
| 1 | 8 | 6 | **0** | 1 |
| 2 | 15 | 9 | 1 | 0 |
| 3 | 19 | 11 | **0** | 1 |
| 4–7, 9 | 21 | 12 | 1 | 1 |

The 12-dimensional `d = 3` odd space at `J ≥ 4` is the one found by `odd3.py` (`PARITY_ODD_HEAVY.md` §4.2). **No symmetric odd structure at
`J = 1, 3`**: there are no `J = 3` (or `J = 1`) states for identical `T`, consistent with the flat-space `+−` threshold `J ≥ 4`. The odd-`J`
states of the near-gap set are therefore `J = 5, 7, 9, …`.

**Validation against `odd3.py` at matched points** (`cache/p1b_validate_odd3.txt`). `odd3.xspace_odd` builds CHL's parity-odd helicity
structures `X^O_same`, `X^O_opp` (arXiv:2102.08160) in position space from their momentum-space definition — independently of any
conservation solve. Our structure was evaluated at odd3's embedding vectors (frame at infinity, 8 random Euclidean configurations, complex
null polarisations; Euclidean vectors mapped to lightcone components by `u = e₁ + ie₂`, `v = −(e₁ − ie₂)`, `y = e₃`):

| `(Δ, J)` | ours / `X^O_same` | ours / `X^O_opp` |
|---|---|---|
| (23.7, 2), (101.1, 4) | constant, spread `≤ 2·10⁻¹³`, `|Im/Re| ≤ 2·10⁻¹⁴` | not proportional (spread `O(10²)`) |
| (23.7, 5), (43.697, 5), (101.1, 7), (43.697, 19), (78.320, 35), (1003.7, 5), (1003.7, 9) | not proportional | constant, spread `≤ 2.5·10⁻¹²`, `|Im/Re| ≤ 4·10⁻¹⁴` |

So: [DV] the conserved symmetric odd structure **is** CHL's `X^O_same` at even `J` and `X^O_opp` at odd `J` (as `KERNELS_D3.md` §5.1
states), including at the two Regge-test points of `PARITY_ODD_HEAVY.md` §4 (`J = 19, 35`); and the proportionality constant is **real**:
the Lorentzian-real `ε`-structure is a real multiple of `X^O = iS` (`S` Euclidean-real). Our `ε` in lightcone components is `2i ×` the
Euclidean one, which is the expected Wick factor. Hence CHL's `X^O` (not `S`) is the real-coefficient (unitarity) basis.
