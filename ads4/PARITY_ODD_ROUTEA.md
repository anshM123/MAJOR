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
* The Regge-moment values `3.77 / 3.30` (`PARITY_ODD_HEAVY.md` §4) are at `J/m ≈ 0.51` (the `β = 1` line), far from the fixed-`J` heavy
  limit; the crystal at their `J = 19` state is not converged (above) and the `J = 19` trend in `Δ` is being computed (`cache/p5_J19.log`).
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
`F_cr/(aN) = +0.191` (`J = 5`; LO alone `+0.200`), `+0.159` (`J = 7`; LO `+0.180`). `K₃`: see §8.

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
| odd sector at `Δ_g = 100` | `n = 1` with the `K₂` numbers, **not controlled** (`Δ = 100` row unconverged at `K₃`); `a`-independent crystal sign positive at `K₂` |
| odd `J ≥ 11` | [ND] (unchanged base) |
| H12′ for route A′ (parity-odd binding state) | shape and sign [NE]; normalisation constant irrelevant (positive, state-independent) |

## 8. Hand-over (state at the end of this round)

* **Running / to finish** (resumable, `p3_run.py`): `K₃` (`L1_DEG=8`) `Δ = 100` rows for `J = 5` (`Ĉ₃` running) and `J = 7` (`p_queueC.sh`);
  `N` at `J = 19`, `Δ = 100…800` (`cache/p5_J19.log`). When done: `L1_DEG=6 python3 p7_crystal_sign.py 5 7` (K₂ vs K₃ crystal sign) and
  `P4_K=3 P4_KINDS=C2,C3H python3 p4_odd_analysis.py 5 7`.
* **Files.** `p1_tto_odd.py` (structures), `p1b_validate_odd3.py`, `p2_shadow_odd.py` (`P3_ANTI=1`: antisymmetric even set at odd `J`),
  `p3_odd_action.py` (+ `p3_run.py`; `L1_DEG=8` → `K₃` caches `*_k3.pkl`), `p4_odd_analysis.py` (LO check, `Γ` grid; `P4_DESHIFT`, `P4_K`,
  `P4_KINDS`), `p5_N.py` (+ `p5_run.sh`), `p6_odd_neargap.py` (`P6_KINDS`, `P6_K3`), `p7_crystal_sign.py`. Caches `cache/p*`.
* **Open.** (1) `Δ_g = 100`: needs a non-perturbative evaluation at `m ≈ Λ` (as for the even sector, H0 §23.4). (2) Closed form of
  `N∞ = 2.31756`, and whether `odd3.xspace_odd` vs `chl_xspace` differ by it. (3) Odd `J ≥ 11`. (4) Combined even + odd `n` at 300/1000
  is `n = 1` given H0 §23.4 (even) and §6 here (odd); a single combined run of `l9_neargap` + `p6` was not made.
