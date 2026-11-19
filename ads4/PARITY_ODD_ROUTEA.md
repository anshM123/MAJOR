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

## 2. Shadow coefficient of the odd structure (`p2_shadow_odd.py`, `cache/p2_shadow_odd.txt`) [NE, validated]

Same route as `l5_shadow_matrix` (Laplace expansion of the `x₅` shadow integral around `e`, relative `O(1/Δ²)`), fitted **jointly** over
[even…, odd] so that the even–odd block is tested. The odd structure carries `(X₃₄X₃₅X₄₅)^{−1/2}` (LHS) and `(X₃₄X₃₆X₄₆)^{−1/2}` (RHS),
with exactly the `m_ij` sign/branch conventions of `l5_shadow_matrix` (folded into the numerators). The polarisation sample is symmetrised
under the transverse reflection `y → −y`, which maps the configuration to itself, so that the parity check is exact.

| `J` | order-2 fit residual `Δ = 100 / 300 / 1000` | `ν/Δ`-dependence of `β_OO` | polarisation-dependence | even–odd entries | `arg β_OO / π` |
|---|---|---|---|---|---|
| 0 | `3.0·10⁻⁴ / 1.1·10⁻⁵ / 2.9·10⁻⁷` | `≤ 2·10⁻¹³` | `3·10⁻⁴ / 1·10⁻⁵ / 3·10⁻⁷` | `≤ 3·10⁻¹⁶` | `0` (`≤ 6·10⁻⁵`) |
| 2 | `4.8·10⁻⁴ / 1.8·10⁻⁵ / 4.8·10⁻⁷` | `2·10⁻⁴ / 9·10⁻⁶ / 2·10⁻⁷` | `2·10⁻⁴ / 7·10⁻⁶ / 2·10⁻⁷` | `≤ 3·10⁻¹⁶` | `0` |
| 4 | `2.9·10⁻⁴ / 1.2·10⁻⁵ / 3.3·10⁻⁷` | `7·10⁻⁵ / 4·10⁻⁶ / 1·10⁻⁷` | `2·10⁻⁴ / 7·10⁻⁶ / 2·10⁻⁷` | `≤ 4·10⁻¹⁴` | `0` |
| 5 | `1.3·10⁻⁴ / 4.1·10⁻⁶ / 1.0·10⁻⁷` | `2·10⁻⁵ / 6·10⁻⁷ / 2·10⁻⁸` | `9·10⁻⁵ / 3·10⁻⁶ / 8·10⁻⁸` | `0` (1×1) | `±1` |
| 7 | `5.9·10⁻⁴ / 1.8·10⁻⁵ / 4.5·10⁻⁷` | `5·10⁻⁴ / 2·10⁻⁵ / 4·10⁻⁷` | `3·10⁻⁴ / 9·10⁻⁶ / 2·10⁻⁷` | `0` | `±1` |
| 9 | `1.7·10⁻³ / 4.9·10⁻⁵ / 1.2·10⁻⁶` | `4·10⁻⁴ / 2·10⁻⁵ / 4·10⁻⁷` | `2·10⁻³ / 5·10⁻⁵ / 1·10⁻⁶` | `0` | `±1` |

* All three diagnostics scale as `Δ⁻³` (the `K₃` truncation), as for the even matrix (§21.2 of H0).
* `β_EO = β_OE = 0` exactly (parity), and **`β_OO` is real**, with sign `(−1)^J` relative to the even `β` (`+` at even `J`, `−` at odd `J`).
  Reality of `β_OO` is a branch check: a wrong branch of the half-integer powers would give a factor `±i`. (A `(−1)^J` is the expected
  convention sign of the `x₆ = (−e)⁺` two-point contraction; it is compensated in the block by the same sign in the `𝒪`-index contraction,
  see §3.)

## 3. The odd structure in the route-A′ crystal integrand; parity and even–odd interference (`p3_odd_action.py`) [DV, round-off]

**Integrand.** `l6_spin_action.build_spin` / `l9_points.build_spin` (kinds `C2`: `q²/V⁴`, `f = v′−u′`, `j′ = 0`; `C3H`: `o_H q/V⁴`, `f = −1`,
`j′ = 1`, the `o_H` completion adopted in H0 §23.3), with the structure axis extended to [even…, odd] on both sides: the odd numerators
(batched over the `T₂`/`𝒪` polarisation samples; `ε` by cofactors on TPS) times `(m₁₂m₁₅m₂₅)^{−1/2}` (12 side) and `(m₃₄m₃₅m₄₅)^{−1/2}` (34
side, `𝒪̃` at `3−Δ`; the same convention as p2's LHS, so the `𝒪̃` normalisation cancels against `β`). Block `M_ab = Σ_c P_ac (β⁻¹)_bc`.
At `n₀ = −1` the `j′ = 1` celestial sign of `Ĉ₃` is `−1`.

**Checks** (`Δ = 1000`, `x = 0.4`, `J = 4` (3×3) and `J = 5` (1×1), both kinds, all Laplace orders separately):
* **Branch.** `log B` (CLS-D's fast exponent) equals the sum of the individual `log m_ij` to `4·10⁻¹⁶` at both `S⁰` points. So the
  half-integer powers of the odd structure, which are `(Δ+J) → (Δ+J+1)` shifts of the even ones times integer powers, sit on the same
  branch as the `Δ`-dependent powers (whose branch is fixed by the reality of the even per-state constants and of `β`).
* **Parity / interference [DV].** `M(n₀ = −1)/M(n₀ = +1) = +1` on the EE and OO entries and `−1` on every EO entry, to `≤ 2·10⁻¹³`, at
  orders `ε⁰, ε², ε⁴`, for `Ĉ₂` and `Ĉ₃`. The functional is the sum over `S⁰ = {n₀ = ±1}`, so **even–odd interference cancels identically**,
  order by order; EE and OO are doubled (a common factor). This is the crystal version of `PARITY_ODD_HEAVY.md` §2.2/§4.3(a). Hence the
  per-`n₀` mixed entries (which are `≈ i ×` real relative to EE, e.g. `M₀[E₁,O] = 0.03 + 8.2i` vs `M₀[E₁,E₁] = 10.7`) never enter, and the
  `Re` prescription acts on the block-diagonal EE ⊕ OO action. **Phase convention (coordinator's point (1)):** the odd structure used is
  Lorentzian-real (real-coefficient `ε` in lightcone components; §1: a *real* multiple of CHL's `X^O`), so its OPE coefficient is real and
  `Re[Ĉ/N̂]` is taken in a real-coefficient basis; with the EO block absent, "Re" never turns into "Im".

## 4. The normalisation `N` (`p5_N.py`; `cache/p5_*`) [NE + DV ingredients]

**Method (same state, no seed normalisation).** At odd `J` the CHL even structure `X^E_opp` lives in the `1↔2` **antisymmetric** parity-even
conserved space (2-dimensional with `X^E_same`; built like l4 with the opposite symmetry). The crystal action is computed on the structure
set [anti-even (2), odd (1)] (`P3_ANTI=1`; shadow matrix from p2 on the same set; even–odd blocks again exactly zero). `X^E_kind` is
decomposed in our real basis by an exact fit to `chl_xspace.xspace_structure` at 12 matched Euclidean points (residual `≤ 5·10⁻¹⁴`, real
coefficients), `X^O_kind = p1/κ_O` (§1). Then
`N_kind := Ĉ[W(X^O_kind, X^O_kind)] / Ĉ[W(X^E_kind, X^E_kind)]` at the **same** `(ν, Δ, J)`: every common factor (fast exponent, ν-only
factors, scalar reference, `N̂`) cancels. (Integer `Δ` is avoided: CHL's x-space formula has Γ-poles there; rows at `Δ = 1600.37, 3200.37`.)

| state | functional | `N_kind(x)`, `x = 0.2 … 0.7` | across Laplace orders 0 / ≤ 2 |
|---|---|---|---|
| `J = 5`, `Δ = 1600.37` and `3200.37` | `Ĉ₂` | `N_opp = 2.3175632436` at every `x` (`Im ≤ 5·10⁻⁹`) | identical |
| `J = 5`, `Δ = 3200.37` | `Ĉ₃` (`o_H`) | `2.3175632436` | identical |
| `J = 19`, `Δ = 1600.37` | `Ĉ₂` | `2.31756` (`Im 3·10⁻⁶`) | identical |
| `J = 7`, `Δ = 1600.37, 3200.37`; `J = 5`, `Δ = 400.37` | `Ĉ₂` | `2.31756` (`Im ≤ 3·10⁻⁷`) | identical |
| `J = 4`, `Δ = 400.37, 3200.37` | `Ĉ₂` | `N_same = 0.43149` | identical |
| `J = 4`, `Δ = 1600.37` | `Ĉ₂` | `N_same = 0.43149 = 1/2.31756` | identical |
| `J = 19`, `Δ = 43.697` (the §4 Regge-test state; `J²/m ≈ 10`, outside the crystal regime) | `Ĉ₂` | `3.246 / 3.249` (order 0), `3.14 / 3.05` (order 2) | not converged |

**Reading.**
* The CHL-normalised odd and even `+−` blocks are **exactly proportional under both functionals, with a real positive constant**, uniformly in
  `x` and at every Laplace order — the crystal counterpart of `PARITY_ODD_HEAVY.md` §4.3(b) ("`N` is a normalisation of the whole block").
* The constant is **`N∞ = 2.31756…`, independent of `J` (5, 7, 19), of `Δ` (400 / 1600 / 3200: 6–10 digits), of `x`
  and of the functional**, and `N_same = 1/N_opp` exactly. So `N` is a relative normalisation of CHL's position-space odd vs even helicity
  structures (the `opp`/`same` pair is rescaled by `√N∞ = 1.5224` and its inverse) — not dynamics, not `O(1/m)`. I have not identified
  `2.31756` in closed form (it is not one of the obvious `π`, `√2`, `√3`, `log 2` combinations; `mpmath.identify` finds nothing).
* `J = 19` trend (`Ĉ₂`, `x = 0.2, 0.4`): `N = 2.3001 / 2.3162 / 2.31742 / 2.31755` at `Δ = 100.37 / 200.37 / 400.37 / 800.37`; the deviation
  from `2.31756` falls like `Δ^{−3.5}` and its `x`- and order-dependence like the Laplace truncation, i.e. it is the crystal's truncation
  error at `J²/m = O(1)`, not a `Δ`-dependence of `N`. At `Δ = 43.7` (`J²/m ≈ 10`) the crystal fails (3.1–3.25, order-dependent). So the
  crystal finds **no `Δ`-dependence of `N` anywhere it converges**, while the Regge-moment values are `3.77 / 3.30` (`PARITY_ODD_HEAVY.md` §4,
  `J/m ≈ 0.51`). Together with the exact constancy in `k`, `ν`, kernel found there, the most economical reading is that the Regge
  pipeline's `W[O,O]` carries an extra seed-normalisation factor (its own suspicion, §4.3(c)) on top of `N∞`; not verified [I].
  What is established: **`N > 0`** (the odd exchange enters with the same sign as an even `+−` exchange of the same `(Δ, J)`), and `N` is a
  state-*independent* constant in the heavy limit. The expectation "`N = 1`" assumed that CHL's x-space odd and even structures have equal
  helicity moduli; the computed blocks say their moduli differ by `√N∞`. Note that §1's validation fixes the *shape* of `X^O` exactly but
  not its normalisation (the symmetric odd space is one-dimensional at every `J`, and any transverse polynomial in CHL's `A_i, B_i` is
  conserved), so a normalisation slip in the x-space construction of `odd3.py` relative to `chl_xspace.py` cannot be excluded by it; the
  exact reciprocity `N_same N_opp = 1` points to a relative factor between CHL's `same`/`opp` x-space structures rather than to dynamics.
  Not resolved further [I].
* **For route A′ this is irrelevant beyond its sign:** the per-state constant `a₊₋` is normalised to the Regge point of each state (as for
  the even sector, H0 §22–§23), and a positive per-state factor cannot change a sign (the coordinator's note: a constant positive `N` does not
  matter; the computed `N` is constant and positive).
* **Phase convention** (coordinator's point (1)): the physical real-coefficient basis is our Lorentzian-real `ε`-structure = a *real*
  multiple of CHL's `X^O` (§1), so the physical odd block is `λ_O² κ_O² W[X^O, X^O] = λ_O² κ_O² N W[X^E_opp, X^E_opp]` with `λ_O` real: same sign
  as the even `+−` block.

## 5. Leading order: the parity-odd exchange gives the stripped flat `+−` density (`p4_odd_analysis.py` (a)) [NE]

`Ĉ₂` order-0 action of the odd structure at odd `J`, divided by the scalar `k = 2` crystal action at the same `(ν, Δ, J)` and multiplied by
`P_J(1−2x²)` (absolute form, regular through zeros), compared with `c·d^J₄₄(1−2x²)/(1−x²)²` — i.e. `W₊₋ = d^J₄₄/((1−x²)²P_J)` relative to
the scalar, the stripped `+−` shape of §21.3 of H0 (`ρ₊₋/ρ₊₊ = d^J₄₄/((1−x²)⁴P_J)`, `ρ₊₊/scalar = (1−x²)²`):

| `J = 5`, ratio to the `x = 0.5` value | `x = 0.2` | 0.3 | 0.4 | 0.6 | 0.7 |
|---|---|---|---|---|---|
| `Δ = 1600` | 0.96855 | 0.99739 | 0.97097 | 0.97776 | 0.89711 |
| `Δ = 3200` | 0.98425 | 0.99879 | 0.98550 | 0.98892 | 0.94646 |
| Richardson (`2·3200 − 1600`) | 0.99995 | 1.00019 | 1.00003 | 1.00008 | 0.9958 |

* The order-0 deviation is a clean `1/Δ` term (it halves), and its `Δ → ∞` limit is the stripped `+−` shape to `≤ 2·10⁻⁴`
  (`4·10⁻³` at `x = 0.7`, next to a zero of `P₅`). The constant `c` is real (`arg c/π = 0.002` at order 0, `→ 0`).
* Same for `J = 7, 9` (max absolute deviation `2.6·10⁻² → 1.3·10⁻²` and `5.2·10⁻² → 2.6·10⁻²` from `Δ = 1600 → 3200`, i.e. `∝ 1/Δ`, at the
  `x = 0.7` point next to a zero of `P_J`; `arg c/π ≤ 0.003 → 0.001`).
* **`Ĉ₃`** (`o_H`): `(Ĉ₃/Ĉ₂)·(2−y)/x² = 1.998i` at every `x` (`J = 5, 7`, `Δ = 3200`, order 0; spread `≤ 3·10⁻³` except one point next to a
  zero, `1.2·10⁻²`): exactly the `σ₊₋ = −1` image of H0 §23.3's `++` value `−2i`. So the odd exchange has the full route-A′ `+−` LO action
  (`(2−y)(…)d̃^J₄₄` in `Ĉ₂`, `−u(…)d̃^J₄₄` in `Ĉ₃`), the `Ĉ₃` phase being ν-only (absorbed in `N̂₃`).
* So **a parity-odd exchange at odd `J` reproduces the stripped flat `+−` density at leading order** — the claim the binding state rests on.
  (The `1/Δ` part of order 0 is not an `O(1/m)` term of the action: it combines with `K₁`; the full corrections are §6.)

## 6. `O(1/m²)` corrections of the odd sector and the per-state test (`p3_run.py`, `p4_odd_analysis.py` (b), `p6_odd_neargap.py`)

**Method (identical to H0 §22–§23 for the even sector, coordinator's point (2)).**
* Grid: `J = 5, 7, 9` (`J = 3` has no state, §1); `Δ = 100, 200, 400, 800, 1600, 3200`; `x = 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, …, 0.8`;
  each point with its heavy reference `Δ_ref = 2·10⁴` at the **same `ν`** and kind; `Ĉ₂` (`q²`) and `Ĉ₃` (`o_H`), `K₂` (DEG 6), orders kept
  separate; `n₀ = +1` (the OO entry is identical at `n₀ = −1`, §3).
* **Linear (prefactor) form**: the state series is multiplied order by order by `exp(−κ^{ref})` and truncated at `K` (l9's chain; the
  log/cumulant form is not used). This removes every ν-only factor (= `N̂_k(ν)` of §R.3) exactly to the computed order.
* Per-state constant `a₊₋` from `Ĉ₂` on `x ≤ 0.25`, `ν ≥ 15` (common to `Ĉ₃`), the **CHMRS de-shift** of H0 §23.2 applied (scalar `g(ν)`,
  [I] as there), and `Γ_{k,+−}(x; m) = m²[Re G/(a env) − σ D₊₋]` in `c₂` units × LO shape. `Re` acts on the real-coefficient OO entry (§3–§4).
  The literal CLS-D variant (`P4_DESHIFT=0`) is also computed, to flag where the choice matters.
* Near gap: `F = base₊₋ + Σψ₂ w env₂ [Γ₂/m² − (ρ₂ − D ρ₂(0))] + Σψ₃ w u env₃ [Γ₃/m² + (ρ₃ − D ρ₂(0))]` on every odd-`J` state of the
  exact-skeleton set (l9_neargap conventions; base reproduces the binding `+−`, `J = 5` values `0.2256 / 0.2054 / 0.2021`). For odd `J` the
  action is 1×1 (no mixing: one structure, one helicity; EO cancels). `t*` = scale factor on the whole correction at which the odd sector
  first goes negative.

**LO consistency and size of the corrections** (`cache/p4_odd_analysis_C2_C3H.txt`). Asymptotic (`Δ ≥ 400`, stable to `≤ 10` between rows)
graviton-specific absolute corrections on `x ∈ [0.3, 0.7]`: `Ĉ₂`: `|Γ₂| ≤ 30` (`J = 5`), `≤ 160` (`J = 7`), `≤ 145` (`J = 9`); `Ĉ₃`:
`|Γ₃| ≤ 140 / 92 / 191` plus an `m`-independent Regge-point offset of `Ĉ₃` relative to the common `a₊₋` (`Γ₃(x→0) ≈ 250 / 820 / 1510` at
`J = 5 / 7 / 9`; it is the "`k`" pattern of H0 §19.3, here computed). Imaginary parts (dropped by `Re`): `m·|Im| ≤ 0.9` (`Ĉ₂`) and `≤ 0.03`
(`Ĉ₃`) at `Δ ≥ 400` — no `O(1/m)` imaginary mixing exists for odd `J` (1×1). The `Δ = 100` row (`m ≈ 98`) is large: `|Γ| ≤ 1.2·10³ /
3·10³ / 3.3·10³` (`J = 5 / 7 / 9`), and **not converged**: at `K₃` (`J = 5`, `Ĉ₂`) the per-state constant `a₊₋` drops by a factor 0.53 and `Γ₂`
moves to `±4.5·10³` (the small-`ν` problem of H0 §23.4, now also for the odd sector).

**Per-state test** (`cache/p6_odd_neargap_*`; odd-`J` states of the exact-skeleton set; `J = 5, 7, 9` covered, `J ≥ 11` unchanged):

| `Δ_g` | functional | odd sector `n = 1`? | min `F/N_s` (state) | uncorrected | `t*` | max `|ΔF|/N_s` (J = 5 / 7 / 9) |
|---|---|---|---|---|---|---|
| 1000 | CHMRS-type (de-shift; §23 choice) | **yes** | 0.2023 (`J = 5`, `1.01Λ`) | 0.2021 | (−210, > 10⁴) | 0.0005 / 0.0020 / 0.0029 |
| 300 | CHMRS-type | **yes** | 0.2072 (`J = 5`, `1.02Λ`) | 0.2054 | (−57, > 10⁴) | 0.0058 / 0.0176 / 0.0080 |
| 100 | CHMRS-type | yes, **not controlled** | 0.2217 (`J = 5`, `1.05Λ`) | 0.2256 | (−29, +5.6) | 0.0056 / 0.141 / 0.096 |
| 1000 | literal CLS-D (no de-shift) | yes | 0.2023 | 0.2021 | (−66, > 10⁴) | |
| 300 | literal CLS-D | yes | 0.2077 | 0.2054 | (−6.1, > 10⁴) | |
| 100 | literal CLS-D | yes, not controlled | 0.2505 | 0.2256 | (−0.8, > 10⁴) | |
| 1000 / 300 / 100 | `Ĉ₂` only (CHMRS-type) | yes | 0.2021 / 0.2050 / 0.2522 | | (−284, +169) / (−48, +30) / (−7.3, +4.3) | |

* `t*` is the factor on the *whole* computed correction at which an odd-`J` state first turns negative; `> 10⁴` means no positive factor
  breaks it. The computed correction (`t = 1`) keeps the odd sector positive in every variant. **The binding state stays `+−`, `J = 5`,
  and its margin is essentially unchanged** (`ΔF/N_s = +0.0002 / +0.0018 / −0.0039`).
* **Against the §19.3 tolerances** (worst pattern `|δc₂| ≲ 630 / 3.4·10³ / 2.7·10⁴` on `x > 0.3`): the asymptotic odd-sector corrections
  (`≤ 160` for `Ĉ₂`, `≤ 1.5·10³` incl. the `Ĉ₃` offset, which is concentrated at small `x`) are inside by factors `≳ 4 / ≳ 20 / ≳ 170` on `x ≥ 0.3` (`|Γ| ≤ 190`) at
  `Δ_g = 100 / 300 / 1000` — but at `Δ_g = 100` the relevant rows are not converged (above), so that comparison is not a statement about
  the true `m ≈ 100` correction.
* **Where the CHMRS/CLS-D choice matters:** only at `Δ_g ≤ 300` and only in the margin `t*` (lower end −57 vs −6.1 at 300; −29 vs −0.8 at
  100); `n` and the binding state do not depend on it.

**`Δ_g = 100`, `a`-independent check** (`p7_crystal_sign.py`): `a₊₋ > 0` is a common positive factor of `Ĉ₂` and `Ĉ₃`, so the sign of the
crystal's own `ψ`-integrated action at the `Δ = 100` row (`m ≈ 98 ≈ Λ`) does not depend on the unconverged normalisation. `K₂`:
`F_cr/(aN) = +0.191` (`J = 5`; LO alone `+0.200`), `+0.159` (`J = 7`; LO `+0.180`). `K₃` (`J = 5`): `F_cr/(aN) = +0.322` with
`a_{K₃} = 0.533 a_{K₂}`, i.e. the normalisation-free crystal action itself is `F_cr = 1.76·10⁻¹⁰ N` (`K₂`) vs `1.58·10⁻¹⁰ N` (`K₃`): **positive at
both orders and changed by 10 %**, although `a` and the individual `Γ` are not converged. So at the `Δ_g = 100` binding region the sign of
the odd-`J = 5` action is stable under the next Laplace order; its size relative to the model base is not controlled.

## 7. Graded summary

| statement | grade |
|---|---|
| conserved symmetric parity-odd `⟨TT𝒪_{Δ,J}⟩` in `d = 3`, covariant, exact `Δ`-dependent coefficients; counts 1 (`J = 0, 2, ≥ 4`), 0 (`J = 1, 3`) | [DV] (null space; gap `10⁻¹⁴` vs `10⁻¹`) |
| = CHL's `X^O_same` (even `J`) / `X^O_opp` (odd `J`) at matched points, real proportionality constant (Lorentzian-real basis = real multiple of `X^O`) | [DV] (9 states, spread `≤ 2.5·10⁻¹²`) |
| no odd-`J` `+−` state at `J = 3` | [DV] (count) |
| odd shadow coefficient `β_OO` (x₅-saddle, `O(1/Δ²)`): real, `β_EO = 0`, residuals `∝ Δ⁻³` | [NE], validated |
| even–odd interference cancels identically in `Ĉ₂`, `Ĉ₃` (S⁰ sum; EO flips sign, EE/OO equal, every Laplace order) | [DV] (round-off) |
| half-integer powers on the validated `Δ`-branch (`log B` = sum of logs) | [DV] |
| parity-odd odd-`J` exchange gives the stripped flat `+−` density in `Ĉ₂` (`d^J₄₄/(1−x²)²` rel. scalar) and `Ĉ₃` (`Ĉ₃/Ĉ₂ = 2i x²/(2−y)`), `J = 5, 7, 9` | [NE] (`1/Δ` convergence; Richardson `≤ 2·10⁻⁴`) |
| odd/even block ratio for CHL-normalised structures: `N_opp = 2.31756…` exactly constant in `x`, Laplace order, `Δ ≥ 400`, `J = 5, 7, 19`, functional; `N_same = 1/N_opp` | [NE] |
| `N > 0` (odd exchange has the same sign as an even `+−` exchange at the same state) | [NE] |
| identification of `2.31756`; reconciliation with the Regge-moment `3.77 / 3.30` | **[I]** / open (likely a CHL x-space normalisation slip, §4; Regge values at `J/m ≈ 0.5`) |
| odd-sector `O(1/m²)` corrections (linear chain, reference subtraction, CHMRS de-shift), `J = 5, 7, 9`, `Δ ≥ 200` | [NE] |
| **odd sector `n = 1` at `Δ_g = 1000` and `300`** with `Ĉ₂ + Ĉ₃` corrections; binding `+−`, `J = 5` unchanged (0.2023 / 0.2072); `t* ≥ 57` / `≥ 6` (CLS-D) | **[NE + M]** |
| odd sector at `Δ_g = 100` | `n = 1` with the `K₂` numbers, **not controlled** (`Δ = 100` row unconverged at `K₃`); the `a`-independent crystal action at `m ≈ Λ` is positive at `K₂` and `K₃` (`J = 5`, 10 % change) [NE] |
| odd `J ≥ 11` | [ND] (unchanged base) |
| H12′ for route A′ (parity-odd binding state) | shape and sign [NE]; normalisation constant irrelevant (positive, state-independent) |

## 8. Hand-over (state at the end of this round)

* **Running / to finish** (resumable, `p3_run.py`): `K₃` (`L1_DEG=8`) `Δ = 100` row for `J = 7` (`p_queueC.sh`, `cache/p3_run_k3_*.log`).
  `J = 5` `K₃` and the `J = 19` `N` trend are done. When `J = 7` `K₃` is done (`cache/p3_pts_*_J7_k3.pkl`): `python3 p7_crystal_sign.py 7`.
* **Files.** `p1_tto_odd.py` (structures), `p1b_validate_odd3.py`, `p2_shadow_odd.py` (`P3_ANTI=1`: antisymmetric even set at odd `J`),
  `p3_odd_action.py` (+ `p3_run.py`; `L1_DEG=8` → `K₃` caches `*_k3.pkl`), `p4_odd_analysis.py` (LO check, `Γ` grid; `P4_DESHIFT`, `P4_K`,
  `P4_KINDS`), `p5_N.py` (+ `p5_run.sh`), `p6_odd_neargap.py` (`P6_KINDS`, `P6_K3`), `p7_crystal_sign.py`. Caches `cache/p*`.
* **Open.** (1) `Δ_g = 100`: needs a non-perturbative evaluation at `m ≈ Λ` (as for the even sector, H0 §23.4). (2) Closed form of
  `N∞ = 2.31756`, and whether `odd3.xspace_odd` vs `chl_xspace` differ by it. (3) Odd `J ≥ 11`. (4) Combined even + odd `n` at 300/1000
  is `n = 1` given H0 §23.4 (even) and §6 here (odd); a single combined run of `l9_neargap` + `p6` was not made.

## 9. Combined even + odd `n` at `Δ_g = 300 / 1000` (CHMRS-type, `Ĉ₂ + Ĉ₃` corrections, `K₂`) [NE + M]

Even sector: `cache/l9_neargap_C23H_K2_100_300_1000.txt` (H0 §23.4; even `J ≤ 10`, 2×2). Odd sector: §6 here (`J = 5, 7, 9`, 1×1).
The two sectors are block-diagonal (EO cancels exactly, §3), so the combined `n` is the joint condition and the combined `t*` (same scale
factor on all computed corrections) is the intersection of the two intervals.

| `Δ_g` | even sector: min `λ/N_s` (state) | odd sector: min `F/N_s` (state) | **combined `n`** | combined min (binding) | combined `t*` |
|---|---|---|---|---|---|
| 300 | 0.2205 (`J = 4`, `1.02Λ`) | 0.2072 (`+−`, `J = 5`, `1.02Λ`) | **1** | 0.2072 (`+−`, `J = 5`) | (−5.1, +5.0) (set by even `J = 10`) |
| 1000 | 0.2103 (`J = 4`, `1.01Λ`) | 0.2023 (`+−`, `J = 5`, `1.01Λ`) | **1** | 0.2023 (`+−`, `J = 5`) | (−37, +39) (set by the even sector) |

Uncovered: even `J ≥ 12`, odd `J ≥ 11` (unchanged base, all with margin ≥ 0.4 in the base model).

## 10. Hand-over for the next step: odd sector at `Δ_g = 100` with the direct QMC integrator (H0 §24, `n1_*`) — not started

The coordinator asked for the odd `J = 5, 7, 9` `+−` states at `Δ_g = 100` from the direct randomised-QMC evaluation of the exact crystal
integrand (validated by the `Δ_g = 100` agent; works for `ν ≳ 20`, i.e. `x ≥ 0.2`), with `K₂` below `x = 0.2`. Context ran out before this
was started. Recipe (do not edit `n1_*`; copy into `p8_*.py`):
1. **Point-wise loader.** `n1_pw._load(name, swaps)` re-executes a module with `l1_tps` → `n1_pw`. Load `p3_odd_action` with swaps
   `{'l1_tps': PW, 'l1_crystal': PW.Cr, 'l5_spin': PW.S, 'l6_spin_action': PW.A6}`. Needed edits in the copy (`p8_odd_pw.py`):
   remove `assert T.DEG in (6, 8)` and the `T.DEG >= 8` branches (PW has no DEG); `p2_shadow_odd` must stay the *real* module (shadow
   from cache `cache/p2_shadow.pkl`; precompute `P2.shadow_odd(Δ, J)` for the needed `Δ` with the real engine first); `struct_sets` uses
   `P2.even_set` (real) and `P1.conserved_odd` (real, reduced TPS) — fine. `eps_cof`, `vec5X/Z`, `tto_odd_samples` only use
   `T.mul/T.const` and `S.lin_samples/S.bm/S.Vfun/S.Hfun/S.Xfun`, which PW provides (`S.bm` is the point-wise product).
2. **Driver.** Copy `n1_run.point` with `gfun = lambda zl: P8.build(zl, nu, Δ, J, kind)[:3]` (kinds `C2`, `C3H`; `n₀ = +1`) and the
   structure axis = [odd] only at odd `J` (1×1); keep the scalar row (`_both`). Rows: the `Δ` whose `m(Δ, J) = (m/Λ)·Λ(100)` for
   `m/Λ ∈ {1.0x, 1.05, 1.1, 1.2, 1.3}`, `x ∈ {0.2, 0.25, 0.3, 0.4, …, 0.8}`, plus the references `Δ_ref = 2·10⁴` at the same `ν` (or reuse
   the `K₂` references, which are exact enough at `Δ_ref`: check one).
3. **Analysis.** As `n1_analysis.py` (direct values on `x ≥ 0.2`, `K₂` below), but for the 1×1 `+−` element with `a₊₋` from `Ĉ₂` and the
   CHMRS de-shift (`p4_odd_analysis.grid` conventions); then `p6_odd_neargap.py` with the direct grid (add a `P6_DIRECT` switch analogous
   to `P6_K3`). Cost estimate: per point ≈ the even `n1` cost (1×1 instead of 2×2 prefactor), ×2 kinds.
4. **Still running at hand-over:** `p_queueC.sh` (`K₃`, `Δ = 100`, `J = 7`, `Ĉ₂` then `Ĉ₃`; resumable). Afterwards
   `python3 p7_crystal_sign.py 7` gives the `a`-independent `K₂`/`K₃` sign check for `J = 7` (as done for `J = 5`, §6).

## 11. Odd sector at `Δ_g = 100` with the direct QMC integrator (`p8_*.py`): `x ≥ 0.2` now non-perturbative and equal to `K₂`; `n = 1` in the §6/H0-§24.5 protocol, but **not robust** to the treatment of `x < 0.2` (`J = 9`, also `J = 7` at `m ≈ Λ`) — odd sector at 100 undecided

*Written 19 Nov 2026 by the `Δ_g = 100` odd-sector agent, following the §10 recipe. Unrefereed, not committed. Even-sector files (`n1_*`)
imported only, not edited.*

### 11.1 Point-wise odd integrand and driver (`p8_odd_pw.py`) [NE]

* **Loader.** `p3_odd_action.py`'s source is exec'd with the `n1_pw` swaps `{l1_tps → PW, l1_crystal → PW.Cr, l5_spin → PW.S,
  l6_spin_action → PW.A6}` and three asserted textual edits (the `assert T.DEG in (6, 8)`, the `T.DEG ≥ 8` K3 import branch and the
  `L1_DEG` default removed; `EXP` is never used). `p1_tto_odd`, `p2_shadow_odd` are imported first, so they stay the *real* (TPS) modules:
  structure coefficients from `P1.conserved_odd`, shadow `β_OO` from `P2.shadow_odd` (DEG-6 cache). So `P8.build` is p3's integrand,
  line for line, evaluated point-wise. (Done programmatically rather than as a hand copy, so p3 stays the single source.)
* **Driver** = `n1_run.point` with `gfun = P8.build(…, kind, n₀ = +1)[:3]`, structure axis [odd] (1×1), the scalar `Q_J` on the same
  samples (`n1_run._both`), `n1_direct`'s re-centred log contour, wall/winding masks, antithetic scrambled Sobol, `R = 8` replicas.
  Normalisation: full integral `= exp(logN) E` in p3's units (`exp(logC)(M₀ + M₂ + M₄)·β`). State points `2¹⁴·16` (`2¹³·16` on the extra
  rows at `x ≥ 0.3`), references `Δ_ref = 2·10⁴` at the same `ν` with `2¹²·16`. Cost 40–55 s (state) / 10–20 s (reference) per point.

### 11.2 Validation against the `K₂`/`K₃` crystal chain at large `Δ` (`p8_odd_pw.py validate`, `cache/p8_validate.txt`) [NE]

Direct value vs p3's Laplace partial sums for the 1×1 odd element (`2¹³·16` points):

| point (`J`, `Δ`, `x`, kind) | rel. RQMC error | `|dir/K₀−1|` | `|dir/K₁−1|` | `|dir/K₂−1|` | `|M₄/M₀|` (size of the last `K₂` term) |
|---|---|---|---|---|---|
| 5, 1600, 0.5, `Ĉ₂` | `9.4·10⁻⁵` | `3.0·10⁻²` | `3.8·10⁻⁴` | **`8.2·10⁻⁵`** | `3.6·10⁻⁴` |
| 5, 1600, 0.5, `Ĉ₃` | `9.4·10⁻⁵` | `2.9·10⁻²` | `3.8·10⁻⁴` | **`8.0·10⁻⁵`** | `4.2·10⁻⁴` |
| 9, 1600, 0.6, `Ĉ₃` | `9.7·10⁻⁵` | `1.2·10⁻²` | `5.2·10⁻⁴` | **`6.7·10⁻⁵`** | `5.9·10⁻⁴` |
| 5, 400, 0.7, `Ĉ₂` | `9.9·10⁻⁵` | `7.3·10⁻²` | `3.0·10⁻³` | **`9.1·10⁻⁵`** | `3.1·10⁻³` |
| 5, 400, 0.3, `Ĉ₂` (next to the zero of `d⁵₄₄` at `x = 0.316`) | `2.6·10⁻⁴` | `9.1·10⁻²` | `4.5·10⁻²` | `4.7·10⁻³` | `4.2·10⁻²` |
| 7, 400, 0.4, `Ĉ₂` | `1.7·10⁻⁴` | `2.3·10⁻²` | `2.0·10⁻²` | `1.6·10⁻³` | `2.0·10⁻²` |
| 9, 400, 0.4, `Ĉ₃` | `1.1·10⁻⁴` | `1.4·10⁻¹` | `8.8·10⁻³` | `1.0·10⁻³` | `8.4·10⁻³` |

* At `Δ = 1600` (and `Δ = 400`, `x = 0.7`) direct = `K₂` within the RQMC error (`≤ 10⁻⁴`), with the series converging monotonically onto
  it; at the three `Δ = 400` points where `K₂`'s last term is `≥ 10⁻²` the gap is 5–10× smaller than that term, i.e. inside `K₂`'s truncation
  error. **`K₃` at these three points** (`cache/p8_k3check_*`, L1_DEG=8; shadow difference DEG 8/6 `≤ 2·10⁻⁵`) closes the gap by 10×:
  `|dir/K₃−1| = 4.5·10⁻⁴ / 1.6·10⁻⁴ / 7.9·10⁻⁵` vs `|M₆/M₀| = 4.5·10⁻³ / 1.6·10⁻³ / 1.1·10⁻³`. **So the point-wise odd integrand, its normalisation, branch continuation and the
  shadow factor are those of the validated crystal chain** [NE].

### 11.3 Direct values at the `Δ = 100` row (`m = 98.0–98.4 ≈ Λ(100) = 98.50`), `x ≥ 0.2` (`p8_analysis.py`, `cache/p8_analysis_direct_100.txt`) [NE]

`Γ_{k,+−}(x)` (c₂ units × LO shape, CHMRS de-shift, `a₊₋` and `c₀` of the row as in p4; the `K₂` column uses the same `a`):

| `J` | max `|Γ₂^dir − Γ₂^{K₂}|` on `x ≥ 0.2` | max `|Γ₃^dir − Γ₃^{K₂}|` | RQMC error (Γ units) at `x = 0.2` / `x ≥ 0.25` | `|Γ|` there |
|---|---|---|---|---|
| 5 | 25 | 94 (at `x = 0.2`: 785 ± 87 vs 691) | 34, 87 / ≤ 5 | ≤ 10³ |
| 7 | 39 | 42 | 18, 33 / ≤ 12 | ≤ 1.3·10³ |
| 9 | 51 | 64 | 24, 32 / ≤ 5 | ≤ 2·10³ |

* **The `K₂` chain is accurate on `x ≥ 0.2` at `m ≈ Λ` for the odd states**: the gap is `≤ 51` (`Ĉ₂`) / `≤ 94` (`Ĉ₃`, at the RQMC-noisiest
  node, 1.1σ) c₂ units, 7–25× below the §19.3 tolerance scale 630, and of the same size as the even-sector `+−` gap of H0 §24.4 (`≤ 47`).
* Acceptance rule (differs from n1 in one point): a direct node is used if its relative RQMC error is `≤ 5 %` **or** its absolute error is
  `≤ 60` c₂ units. The extra clause admits `J = 7`, `x = 0.2`, `Ĉ₂` (relative error 16 % because the LO shape `d⁷₄₄` nearly vanishes there,
  `D = −0.076`; absolute error 18). All 48 direct nodes of the row are used.

### 11.4 Near-gap odd-sector test at `Δ_g = 100`, main protocol (`P6_DIRECT=direct_100`; `cache/p8_neargap_13.txt`, `p8_neargap_2.txt`) [NE + M]

Main protocol = H0 §24.5's even-sector protocol: the §6 rows (`Δ = 100, 200, …, 3200`, linear in `1/m` between rows), direct values on
`x ≥ 0.2` of the `Δ = 100` row, `K₂` elsewhere; `Ĉ₂ + Ĉ₃` (CHMRS de-shift, `o_H`), `a₊₋` per row as in §6. `p6_odd_neargap.py` gained the switch
`P6_DIRECT` (`1` / any `p8_analysis` grid tag) and a `GRID_OVERRIDE` hook; the default path reproduces §6 exactly (0.2217, `t* = (−29.3, +5.6)`).

| rows | odd `n` | min `F/N_s` (binding) | `J = 5 / 7 / 9` min `F/N_s` (at `m/Λ`) | `t*` |
|---|---|---|---|---|
| §6 (Laplace `K₂`) | 1 | 0.2217 (`J = 5`, 1.05Λ) | 0.2217 (1.05) / 0.5863 (1.09) / 0.9026 (1.35) | (−29.3, +5.64) |
| **direct on `x ≥ 0.2` (`Δ = 100` row), `K₂` below** | **1** | **0.2262 ± 0.0008** (`J = 5`, 1.05Λ) | 0.2262 ± 0.0008 / 0.5890 ± 0.0029 / 0.8974 ± 0.0048 | (−153 ± 110, **+5.84 ± 0.17**) |

* Errors: 30 draws of independent Gaussian perturbations of every direct node (σ = replica error), whole p6 chain rerun. The lower `t*`
  is ill-determined (no state is close to failing for negative `t`); the upper one (`J = 7`, `m ≈ 1.09Λ`, `ΔF/N_s = −0.12`) is stable.
* Replacing `K₂` by the direct integral on `x ≥ 0.2` moves the odd numbers by `≤ 0.005 N_s` and leaves the binding state (`+−`, `J = 5`,
  1.05Λ) and `n = 1` unchanged. **Everything rests on `x < 0.2` (`ν ≲ 20`), where only Laplace values exist (11.6).**

### 11.5 Extra rows `m = 1.05/1.1/1.2/1.3 Λ` and a normalisation phase at non-even-integer `Δ` (`p8_k2rows.py`, `p8_analysis.py`) [NE]

Direct (`x ≥ 0.2`) and `K₂` (`x = 0.05…0.8`) points were also computed on rows at the state masses (`Δ = 105.07 … 129.90`, `J`-dependent;
`J = 5` at 1.05, all `J` at 1.1/1.2/1.3), K2 into separate caches (`cache/p8_k2_*`), so that the §6 caches/numbers are untouched.
* **Phase finding [NE].** On these rows the state/reference ratio carries an extra **overall phase `e^{−iπ(Δ−100)}`** (x-independent; e.g.
  `arg/π` shifts by `+0.92 / −0.11 / +0.17 / +0.33` at `Δ = 105.07 / 110.11 / 119.83 / 129.67`, = `−(Δ−100)` mod 2 to `≤ 0.02`). It is `≡ 1`
  on the even-integer rows used so far (100, 200, …, 2·10⁴), which is why p4/l9/n1 never saw it; with only `c₀ = Re(...)` removed, `Re`
  projects a rotated quantity (`J = 5` at 1.05Λ came out with `a₊₋ < 0`). `p8_analysis` removes it (`P8_PHASE=1`, factor `e^{iπΔ}`, identity on
  the old rows); afterwards `a₊₋` is smooth across rows (`J = 5`: `9.2, 9.9, 11.8, 16.3, 21.9 ·10⁻¹⁰` at 1.0…1.3Λ). Its origin (a `Δ`-dependent
  branch convention of the normalisation) was not traced. **Any future use of non-even-integer `Δ` rows (even sector included) needs it.**
* Direct vs `K₂` on `x ≥ 0.2` of these rows, normalisation-free (`m² max|ΔG/env| / max|G_{K₂}/env|`): `≤ 51` (`J = 5`), `≤ 69` (`J = 7`),
  `≤ 171` (`J = 9`), `Ĉ₂`/`Ĉ₃` alike — the same picture as at `Δ = 100`: **`K₂` is accurate on `x ≥ 0.2` throughout `m = Λ … 1.3Λ`**.

### 11.6 Robustness of `n` to the treatment of `x < 0.2` — it is not robust (`cache/p8_neargap_13.txt`) [NE + M]

| variant (all with direct values on `x ≥ 0.2`) | odd `n` | min `F/N_s` | `J = 5 / 7 / 9` |
|---|---|---|---|
| main (11.4): `K₂` on `x < 0.2`, rows 100, 200, … | 1 | 0.2262 | 0.2262 / 0.5890 / 0.8974 |
| `K₃` on `x < 0.2` of the `Δ = 100` row (`a₊₋` kept from `K₂`; the `K₃` fit gives `a < 0` for `J = 7`) | **0** | −1.74 (`J = 9`, 1.11Λ) | 0.1974 / 0.1881 / **−1.735** |
| `K₂` on `x < 0.2`, plus the extra rows 1.05–1.3Λ (phase-corrected) | **0** | −0.53 (`J = 9`, 1.3Λ) | 0.2505 / 0.7230 / **−0.532** |
| same without direct values (all `K₂`) | 0 | −0.48 (`J = 9`, 1.3Λ) | 0.2462 / 0.7222 / −0.476 |

`a`-free check (no per-state constant, no de-shift): the reference-normalised crystal action `S = (Σψ₂G₂ + Σψ₃G₃)/(|Σψ₂G₂| + |Σψ₃G₃|)` at
each row's own mass (`K₂` | direct on `x ≥ 0.2` + `K₂` | direct + `K₃` below where computed):

| `m/Λ` | `J = 5` | `J = 7` | `J = 9` |
|---|---|---|---|
| ≈ 1.00 (`Δ = 100`) | +0.226 / +0.229 / +0.199 | +0.421 / +0.434 / **+0.002** | **−0.064 / −0.098 / −0.061** |
| 1.05 | +0.264 / +0.267 / +0.222 | — | — |
| 1.10 | +0.377 / +0.379 / (no K3) | +0.729 / +0.724 / +0.634 | +1.000 / +1.000 / **−0.647** |
| 1.20, 1.30 | +1.000 (all) | +1.000 (all) | +1.000 (all; no K3) |

* **Where the sensitivity sits.** On every row the `x < 0.2` nodes differ between `K₂` and `K₃` by `O(10³–10⁴)` c₂ units (e.g. `J = 9`,
  1.1Λ, `x = 0.05`: `Γ₂ = +7737` (`K₂`) vs `−3470` (`K₃`)), and between neighbouring rows by as much (`J = 7`, `x = 0.05`: −951 at 1.0Λ,
  +4368 at 1.1Λ, `K₂`). These nodes have `ν = 5–16`, i.e. they straddle `ν₀ = 5.98` where the CHMRS de-shift is singular (H0 §24.1) and
  lie where neither Laplace order nor the direct method (H0 §24.2) is valid. The `K₂` main result is therefore an interpolation through
  meaningless nodes that happens to be benign on the coarse row set; a finer row set or the next Laplace order breaks `J = 9`.
* **Robust:** the binding state `+−`, `J = 5` (min `F/N_s` 0.197–0.251 in every variant; `S > 0` at every row and order), and everything
  on `x ≥ 0.2`. **Not robust:** `J = 9` (sign flips with the `x < 0.2` treatment), `J = 7` at `m ≈ Λ` (`S → 0` at `K₃`; no `J = 7` state
  lies below 1.07Λ, where `S = +0.63` even with `K₃`).

### 11.7 `p_queueC.sh` result: `K₃`, `J = 7`, `Δ = 100` row (`python3 p7_crystal_sign.py 7`, `cache/p7_crystal_sign_J7_J5.txt`) [NE]

The job finished normally (19 Nov 01:08). `F_cr/(aN)`: `J = 7`: `K₂` +0.1587 (LO only +0.1804), `K₃` +0.0525 but with `a_{K₃} = −9.1·10⁻¹¹`
(vs `a_{K₂} = +1.28·10⁻⁹`), i.e. the normalisation-free crystal action is `F_cr = +2.0·10⁻¹⁰ N` (`K₂`) → `−4.8·10⁻¹² N` (`K₃`): **for `J = 7`
the sign of the crystal's own action at `m ≈ Λ` is not stable under the next Laplace order** (unlike `J = 5`: `+1.76 → +1.58·10⁻¹⁰ N`,
re-confirmed). Split by `x` (absolute units; same p4/p8 grids, scratch computation): the `K₂ → K₃` change is `−8.8` (`Ĉ₂`, `x < 0.2`) `+0.4` (`Ĉ₃`,
`x < 0.2`) `+1.6` (`x ≥ 0.2`) in units of `10⁻¹⁰` against `N_LO = 41`, almost all from the `x = 0.05` node (`ν = 4.9`); the direct values on
`x ≥ 0.2` agree with `K₂` there (shift +0.05). The job is no longer needed; its `K₃` row is used in 11.6.

### 11.8 Combined even + odd at `Δ_g = 100` [NE + M + I]

Even sector: H0 §24.5 (direct on `x ≥ 0.2`, `K₂` below): `n = 1`, min `λ/N_s = 0.1814` (`J = 4`, 1.04Λ), `t* = (−0.8, +1.8)`. Odd sector: 11.4.
Block-diagonal (§3), so:

| protocol | combined `n` | combined min (binding) | combined `t*` |
|---|---|---|---|
| direct on `x ≥ 0.2` at the `Δ = 100` row, `K₂` below, rows 100, 200, … (both sectors) | **1** | **0.1814** (even `J = 4`, 1.04Λ); odd binding 0.2262 (`+−`, `J = 5`, 1.05Λ) | **(−0.8, +1.8)** (set by the even sector; odd `(−153, +5.8)`) |

**Caveat (`x < 0.2`).** In both sectors the `x < 0.2` (`ν ≲ 20`) part is Laplace-only and uncontrolled (H0 §24.5: 25–35 % of `Σ|ψ density|`
for `J ≥ 6`). For the odd sector 11.6 shows explicitly that two equally defensible treatments of that region (`K₃` nodes; finer rows in
`m`) give `n = 0` through `J = 9`. The even sector was not subjected to the same two tests (only `K₃` for `Ĉ₂`, `J = 8, 10`, existed), so its
`n = 1` is at least as exposed. **`n` at `Δ_g = 100` remains [I]; no computed number supports upgrading it, and the odd sector now provides
concrete counter-indications from the uncontrolled region.** `Δ_g = 300, 1000`: unchanged (§9).

### 11.9 Graded summary of this round

| statement | grade |
|---|---|
| point-wise odd integrand (`p8`) = p3's; direct = `K₂` within RQMC error at `Δ = 1600` (`≤ 10⁻⁴`); at `Δ = 400` inside `K₂`'s truncation, `K₃` closes the gap 10× | [NE] |
| direct vs `K₂` on `x ≥ 0.2` at `m = Λ … 1.3Λ`, `J = 5, 7, 9`, both kinds: `≤ 51 / 94` c₂ at `Δ = 100`; `≤ 171` normalisation-free on the extra rows | [NE] |
| odd `n = 1` at `Δ_g = 100` in the H0-§24.5 protocol: min `F/N_s = 0.2262 ± 0.0008` (`+−`, `J = 5`, 1.05Λ), `t* = (−153 ± 110, +5.84 ± 0.17)` | [NE + M] |
| that `n = 1` is robust to the treatment of `x < 0.2` | **refuted as computed**: `K₃` nodes or finer `m`-rows give `n = 0` via `J = 9` (−1.74 / −0.53) |
| binding odd state `+−`, `J = 5`, positive in every variant (0.197–0.251; `a`-free `S > 0` at all rows and orders) | [NE + M] |
| `J = 7` normalisation-free crystal action at `m ≈ Λ`: `K₂` +, `K₃` ≈ 0 (sign unstable); at 1.1Λ positive in all variants | [NE] |
| overall phase `e^{−iπ(Δ−100)}` of the state/reference normalisation at non-even-integer `Δ` (removed by `e^{iπΔ}`) | [NE]; origin not traced |
| combined even + odd `n = 1` at `Δ_g = 100`, min 0.1814 (even `J = 4`), `t* = (−0.8, +1.8)`, protocol of 11.8 | [NE + M] |
| **`n` at `Δ_g = 100`** | **[I]** (unchanged); odd sector undecided, with `J = 9` counter-indications from `ν ≲ 15` |

### 11.10 Hand-over

* **Nothing is running** (all queues finished 19 Nov ≈ 03:45). Scripts: `p8_odd_pw.py` (loader, driver, `validate`), `p8_k2rows.py`
  (Laplace on the extra rows; `L1_DEG=8` → `p8_k3_*`), `p8_analysis.py` (grids; env `P8_DIRECT`, `P8_XROWS`, `P8_LOWK`, `P8_AFIX`, `P8_PHASE`,
  `P8_RELMAX`, `P8_GERRMAX`), `p8_neargap.py NMC PARTS` (variants, RQMC propagation, `a`-free table), `p6_odd_neargap.py` (`P6_DIRECT`),
  queues `p_queueD.sh`, `p_queueE.sh`, `p_queueG_C2.sh`, `p_queueG_C3H.sh`. Caches `cache/p8_pts_<kind>_J<J>.pkl` (direct; key `(Δ, ν)`),
  `cache/p8_k2_*`, `cache/p8_k3_*`, `cache/p8_k3check_*`, `cache/p3_pts_*_J9_k3.pkl` (`x ≤ 0.25`/`0.15`), grids `cache/p4_grid_{direct,
  direct_100, k2x, direct_100_k3low_afix}.pkl`, outputs `cache/p8_validate.txt`, `p8_analysis_*.txt`, `p8_neargap_{13,2}.txt`.
* **Next step (both sectors):** the region `ν ≲ 20` decides `n` at `Δ_g = 100`; it needs the curved-contour (Picard–Lefschetz) evaluation of
  H0 §24.7(1), or the graviton/scalar-ratio route of §24.7(2), now for the odd `J = 7, 9` states as well. Before that, the even sector should
  be put through the two robustness tests of 11.6 (extra `m`-rows with the `e^{iπΔ}` phase fix; `K₃` on `x < 0.2` for both kinds).
