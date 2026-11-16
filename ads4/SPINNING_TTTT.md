# `⟨TTTT⟩`: helicity map, mean-field normalisation and the `O(1/m)` question for stress tensors; status of the doubly-spinning blocks

Status: written 16 Nov 2026 by a derivation/numerics agent; **not refereed**; not committed. Code and outputs: `ads4/code/fourpt/tttt/`
(caches in `tttt/cache/`, git-ignored with `fourpt/cache/`). Imports, does not modify, `ads4/code/e2e/`, `ads4/code/mft/`, `ads4/code/fourpt/`.
Builds on `SPINNING_4PT.md` (revised after `REFEREE_SPINNING_4PT.md`), `SPINNING_E2E.md`, `SPINNING_MFT.md`. CHL = Caron-Huot–Li,
arXiv:2102.08160.

Marks: **[D]** derived; **[DV]** derived and verified numerically; **[DVE]** derived and verified end to end on exact blocks; **[NE]**
numerical evidence (exact rationals where stated); **[ID]** closed form identified from exact data and checked at further points; **[Q]**
quoted; **[I]** inferred; **[ND]** not done.

## 0. Summary

1. **`⟨TTO⟩` helicity map, exact at rational points [NE, exact].** Two conserved parity-even structures (plus one `d = 3` relation among the 14
   differential-basis words); exact opposite- and same-helicity lines (momentum-space projection of CHL, calibrated on currents); CHL's
   `T^{E,h}` identified exactly in the conserved basis; the word-wise shadow map sends each helicity line exactly to the same helicity line at
   `(3−Δ,−1−J)` (shadow commutes with helicity for `l = 2`). Two points, `(97/10, 13/3)` and `(53/4, 11/2)`. §3.1–3.2.
2. **MFT normalisation for stress tensors is Weyl-even [D from Q + exact checks].** The shadow factor of CHL's structures is exactly
   `κ_h = 2^{3J+3/2}/A_h` — for `l = 1` *and* `l = 2`, both helicities, both points (12+ digits) — where `A_h` are CHL's MFT helicity factors.
   Hence, with CHL's `l = 2` MFT data [Q], `p^MFT[e_h] = 8^J ×` (scalar factor) `×` (Weyl-even) for the Weyl-covariant helicity families, exactly as
   proved for `l = 1` by exact Fock data. **The MFT (unitarity) normalisation of the graviton-channel helicity couplings adds no `O(1/m)` term;**
   with the reflection argument this gives no `O(1/m)` in the leading heavy `⟨TTTT⟩` Regge moments for parity-even exchanges (conditional on
   CHL's quoted `l = 2` MFT factors and Assumption T; not verified end to end). §3.3.
3. **Currents (goal 4 for `l = 1`) [D, exact Fock data; DV at block level].** Exact two-particle Fock MFT data in our normalisation, helicity-
   diagonal; absolute constant `1/4` per helicity quadruple, `C^± = (2P/m^{2k})(|c_same|² ± |c_opp|²)`, all factors Weyl-even; leading-twist block
   decomposition of `⟨J_1J_3⟩⟨J_2J_4⟩` reproduced with no free constant (`10⁻⁶`). Physical sign of the opposite-sign term `+(|c_same|²−|c_opp|²)` =
   flat `σ_{++} = +1`, `σ_{+−} = −1` (`SPINNING_4PT.md` §4.5–4.6). §2.
4. **`⟨TTTT⟩` blocks and end-to-end moments (goals 1, 3): not done — blocked.** The symbolic word pairs at jet order 8 are too large
   (`ACH × ACH`: 73 095 terms; `AACC × AACC`: not finished after 67 min; 196 pairs), and order-8 jets need extended-precision post-processing. A
   different algorithm (numerical Taylor-mode application of the operators, or a direct recursion) is needed. §4.
5. **Not done:** `l = 2` absolute constant (expected `1/4` per helicity quadruple [I]); a direct `l = 2` Fock computation (would turn CHL's `l = 2`
   MFT factors from [Q] into [D]); closed forms of the `l = 2` helicity lines in `(Δ,J)`; parity-odd exchanges; `k = 2, 3`.

## 1. Tools (all new, `tttt/`)

* **Structures at infinity** (`infframe.py`): `⟨O_1(x)O_2(0)X(∞,ζ)⟩` for every differential-basis word, from the embedding point at infinity
  `P_∞ = (0,1,0)`, `Z_∞ = (0,0,ζ)` (`s_{i5} = −1`, `Z_i·P_∞ = 0`, `Z_∞·P_i = ζ·x_i`), first in the invariants `a·x, b·x, ζ·x, x², a·b, a·ζ, b·ζ`.
  Checked: the finite-point words converge to these values as `x_5 → ∞` (`O(1/L)`, 150-digit test).
* **Conservation at infinity** (`infcons.py` symbolic in `(Δ,J)`; `infcons2.py` exact at rational points by the chain rule in the invariants,
  harmonic projection for spin 2, evaluation on exact 3d configurations with Gaussian-rational null vectors). Reproduces `v_1, v_2` of
  `SPINNING_4PT.md` exactly for `l = 1`. (A bug — one random scale per *component* of the null vectors, which destroys nullness — was found
  and fixed; results before the fix are discarded.)
* **Momentum-space helicity projection** (`ft_helicity.py`): term-wise Fourier transform `FT[(ζ·x)^e|x|^{2t}] = (−i)^e π^{3/2}2^{2t+3+e}
  Γ(t+e+3/2)/Γ(−t) (ζ·p)^e|p|^{−2t−3−2e}` (harmonic polynomial, Bochner), `(a·x) → i a·∂_p`; evaluated at `p = ẑ` with transverse null `a, b`;
  exact version with the Gamma functions reduced to Pochhammer symbols. Calibration on `l = 1`: the lines `T(+,−) = 0` and `T(+,+) = 0`
  are exactly `ĥ_opp` and `ĥ_same` of `SPINNING_4PT.md` §2.3 (exact rationals at two points).
* **CHL structures in momentum space** (`chl_momentum.py`): eq. `Tpm` with the E normalisation of eq. `even odd`, compared with the Fourier
  transform of our structures (frames: `x_CHL = −x`, so `FT_CHL(p) = FT_ours(−p)`).
* **Two-particle Fock space** (`fock2.py`): exact MFT data of `[O O]_{N,J}`, the method of `SPINNING_MFT.md` §2 with two spinning particles
  (`mft/fock.py` one-particle Gram forms).

## 2. Currents: exact MFT data, absolute constant, block-level test (coordinator goal 4 for `l = 1`)

See `SPINNING_4PT.md` §4.6 (written there as part of the referee revision). In short [D from exact Fock data]:
`p^MFT[ĥ_opp] = 8^J P_FK(1,1;Δ,J)/(8Y²(J−1)J(J+1)(J+2))`, `p^MFT[ĥ_same] = 8^J P_FK(1,1;Δ,J)/(8Y²Δ(Δ−1)(Δ−2)(Δ−3))` (all 25 primaries with `n ≥ 1`, `N ≤ 6`,
helicity-diagonal Gram matrix; `n = 0`: one primary carrying the sum); absolute constant `1/4` per helicity quadruple, `C^± =
(2P/m^{2k})(|c_same|² ± |c_opp|²)`, `Q_h = 1 + O(1/m²)` (`out_t2_mft_heavy_l1.txt`); block-level leading-twist test of `⟨J_1J_3⟩⟨J_2J_4⟩` with the exact
doubly-spinning waves: ratio `1 + O(ε)` (`10⁻⁶` at `ε = 10⁻⁵`), no free constant, control `12.6` (`out_t3_lightcone_l1.txt`) [DV].
Physical sign of the opposite-sign term: `+(|c_same|² − |c_opp|²)` for the covariant odd⊗odd kernel (`SPINNING_4PT.md` §4.5, `out_t4_sign_l1.txt`).

## 3. Stress tensors: exact helicity map and the Weyl parity of the MFT normalisation (coordinator goal 2)

**3.1 Conserved space and exact helicity lines [NE, exact] (`t5_helicity_l2.py`, `out_t5_helicity_l2_pointA.txt`, `out_t5_helicity_l2_53_4_11_2.txt`).**
At `(Δ,J) = (97/10, 13/3)` (and `(53/4, 11/2)`, §3.3): conservation at `x_1` and `x_2` in the frame at infinity leaves a 3-dimensional null space
of the 14 words, containing the single `d = 3` relation among the words (found exactly at infinity as well): **two conserved parity-even
`⟨TTO⟩` structures**. The momentum-space projection (calibrated on currents, §1) gives exactly one *opposite*-helicity line (`T(+,−) = 0`) and
one *same*-helicity line (`T(+,+) = 0`) modulo the relation, with exact rational coefficients.

**3.2 Weyl covariance and CHL identification [NE, exact / 30 digits].** For each helicity line `e_h(Δ)`, the word-wise shadow image `T e_h(Δ)`
lies exactly on the same helicity line at `(3−Δ,−1−J)` (the wrong components vanish identically): the shadow map commutes with helicity for
stress tensors too, so `e_h(Δ̃) := T e_h(Δ)` is a Weyl-covariant helicity family. CHL's `T^{E,same}`, `T^{E,opp}` (momentum space, eqs. `Tpm`,
`even odd`) are proportional to `e_same`, `e_opp` with the same ratio `N_h` on both non-zero components — the exact identification of CHL's
`⟨TTO⟩` structures in the differential basis at this point.

**3.3 The Weyl parity of the MFT normalisation [D from Q + exact checks].** With `T^{E,h} = N_h e_h` and the Weyl-covariant `e_h`, the shadow acts on
CHL's structure as `T^{E,h}(Δ) → κ_h T^{E,h}(Δ̃)`, `κ_h = N_h(Δ)/N_h(Δ̃)`. CHL's MFT data are `p[T^{E,h}] = c(1,1)^s(Δ,J)·A_h·(const)` [Q, eq. `free data`],
`A_opp = (−J)_{2l}/(J+1)_{2l}`, `A_same = (Δ−1)_{2l}/(2−Δ)_{2l}`. Then `p[e_h] = N_h² p[T^{E,h}]` and

  `p[e_h](Δ,J)/p[e_h](3−Δ,−1−J) = κ_h² A_h² × [c(1,1)^s(Δ)/c(1,1)^s(Δ̃)]`,

so the normalisation of the Weyl-covariant family is Weyl-even (up to the helicity-independent scalar factor, which is that of the scalar
problem) iff `κ_h A_h` is helicity-independent. Result (`out_t6_kappa_l1.txt`, `out_t5_helicity_l2_*.txt`):

| `l` | `(Δ, J)` | `κ_opp A_opp` | `κ_same A_same` | `2^{3J+3/2}` |
|---|---|---|---|---|
| 1 | (97/10, 13/3) | 23170.4750059 | 23170.4750059 | 23170.4750059 |
| 1 | (53/4, 11/2) | 262144 | 262144 | 262144 |
| 2 | (97/10, 13/3) | 23170.4750059 | 23170.4750059 | 23170.4750059 |
| 2 | (53/4, 11/2) | 262144 | 262144 | 262144 |

(ratio `κ_opp A_opp/κ_same A_same = 1` to `10⁻³⁶` at the `l = 2` point A, to `10⁻¹⁶` at point B). So **`κ_h = 2^{3J+3/2}/A_h` exactly, for currents and for stress tensors**:
the helicity dependence of CHL's shadow factor is exactly the inverse of their MFT helicity factor, and the remainder `2^{3J+3/2} =
2^{3J/2}/2^{3J̃/2}` is the reflection ratio of the convention factor `8^{J/2}` — the same `8^J` that appears in the exact `l = 1` Fock data
(`SPINNING_4PT.md` §4.6), where it is cancelled by the block normalisation `(−1/8)^J`. Hence `p^MFT[e_h] = 8^J × (scalar MFT factor) × (Weyl-even)`
for `l = 2`, as proved for `l = 1` by the exact Fock data. **Consequence: the MFT (unitarity) normalisation of the stress-tensor helicity
couplings introduces no `O(1/m)` term**, and with the reflection argument of `SPINNING_4PT.md` §6.1 (which holds for any `l` once the families
are Weyl-covariant and helicity-pure — now shown) the leading heavy `⟨TTTT⟩` Regge moments with MFT-normalised helicity couplings have
corrections `1 + O(1/m²)` for parity-even exchanges.
Status: [D] given (a) CHL's `l = 2` MFT helicity factors `A_h` [Q] — their `l = 1` analogue is confirmed exactly by our Fock data, the `l = 2` ones
were not recomputed; (b) Assumption T; (c) exactness of the two-point checks (the identity `κ_h A_h = 2^{3J+3/2}` is checked at two points
per `l`, to 12+ digits, not proved symbolically). Not verified end to end (no `⟨TTTT⟩` blocks, §4).

**3.4 What is not done for `l = 2`.** The closed forms of the helicity lines as functions of `(Δ,J)` (only exact values at two points); the heavy
vertex closed forms `V_h` for `l = 2` and hence the absolute two-pair constant (expected `1/4` per helicity quadruple, by the same algebra as for
`l = 1` [I]); a direct `l = 2` Fock computation of `p^MFT` (the helicity lines at double-twist points need the Fourier transform at integer
exponents, where it has contact terms; not attempted); parity-odd exchanges.

## 4. `⟨TTTT⟩` blocks and moments (coordinator goals 1, 3)

**Not done — blocked by cost in the present (symbolic differential-operator) framework [ND].** Measured cost of the symbolic word pairs at
jet order 8 (`feas_l2.py`, `cache/feas_l2.log`; expanded terms of `D^{(12)}_w D^{(34)}_{w′}[K g]`, before code generation):

| word pair | terms | time |
|---|---|---|
| `HH × HH` | 9 | 0 s |
| `AACC × HH` | 22 800 | 45 s |
| `ACH × ACH` | 73 095 | 171 s |
| `AACC × AACC` (4 operators at each vertex) | not finished | > 67 min (still expanding; aborted at the 90-min limit) |

There are 196 pairs (14 × 14, minus nothing: the 3d relation has to be quotiented after evaluation), most of them of the `AACC × AACC` type.
For comparison, the 25 pairs of `⟨JJJJ⟩` have ≤ 2100 terms. Beyond the symbolic size, three further obstacles were identified: (i) the jets of
order 8 in `(P, Q)` need radial derivative sums `Σ(Δ+n)^k P_j^{(l)}` with `k, l ≤ 8` (`e2e/cblock.c` has a binomial table of size 8 and the Python
jet composition `e2e/spinblock.pq_jets` is hard-wired to order 4) — a straightforward extension; (ii) precision: the conserved combination cancels
the leading `(1/r)^{2l}`-enhanced terms of each word (`SPINNING_E2E.md` §2.2), so order-8 jets in double precision are not reliable at `Δ ≳ 30`;
the `__float128` radial sums exist, but the post-processing would have to run in extended precision too; (iii) the number of polarisation
monomials (multi-degree `(2,2,2,2)` in four null polarisations) is in the thousands, so the monomial-coefficient evaluation of `fourpt/spinblock2`
grows accordingly. **Assessment:** an end-to-end `⟨TTTT⟩` test needs a different algorithm — e.g. numerical (Taylor-mode) application of the
`D`-words at each grid point instead of symbolic expansion, or a direct recursion for the conserved doubly-spinning blocks — estimated at days
of development and CPU, not hours. Consequently **goal 3 (end-to-end `⟨TTTT⟩` moments, `±`, `O(1/m)` fits) is not done**; what exists for stress
tensors is the three-point/shadow/normalisation chain of §3 and the heavy-integrand checks of `SPINNING_4PT.md` §7.

## 5. Status

| # | statement | mark |
|---|---|---|
| 1 | Structures at infinity, conservation there (symbolic `l = 1`, exact at rational points `l = 2`), momentum-space helicity projection | [D], checks [DV] (limit `x_5 → ∞`; `l = 1` reproduces `v_1, v_2`, `ĥ_h`) |
| 2 | `⟨JJ⟩` MFT data in the helicity basis, closed forms, absolute constant `1/4`, Weyl-even | [D from exact Fock data] |
| 3 | Two-pair block normalisation `p(−1/8)^J W` (leading twist) | [DV] |
| 4 | Physical sign `+(|c_same|²−|c_opp|²)` of the opposite-sign term (`l = 1`) | [D] + DVE; functional ↔ flat sum rule [I] |
| 5 | `⟨TTO⟩`: 2 conserved structures + 1 relation; exact helicity lines; CHL `T^{E,h}` identified; shadow commutes with helicity | [NE, exact] at two points |
| 6 | `κ_h = 2^{3J+3/2}/A_h` (`l = 1, 2`) | [NE] at two points per `l` (12–36 digits) |
| 7 | MFT normalisation of `l = 2` helicity couplings Weyl-even ⇒ no `O(1/m)` from it | [D from Q] (CHL `l = 2` MFT factors quoted) |
| 8 | No `O(1/m)` in leading heavy `⟨TTTT⟩` moments (parity-even exchanges, MFT-normalised) | [D] conditional on 7 and Assumption T; not end to end |
| 9 | `⟨TTTT⟩` exact blocks, end-to-end moments, `±`, `1/m` fits | [ND] — blocked (§4) |
| 10 | `l = 2` absolute constant; `l = 2` Fock data; parity-odd exchanges; `k = 2, 3` | [ND] |

## 6. Code index (`ads4/code/fourpt/tttt/`)

| file | content | output |
|---|---|---|
| `infframe.py`, `infcons.py`, `infcons2.py` | structures at infinity; conservation (symbolic / exact at rational points) | — |
| `ft_helicity.py`, `chl_momentum.py` | momentum-space helicity projection; CHL `T^{E,h}` | — |
| `fock2.py`, `mft2.py`, `t1_mft_l1.py` | two-particle Fock space; exact MFT data of `[JJ]` in the helicity basis | `out_t1_mft_l1.txt` |
| `t2_mft_heavy_l1.py` | absolute two-pair heavy constant, `l = 1` | `out_t2_mft_heavy_l1.txt` |
| `t3_lightcone_l1.py` | block-level test of the two-pair MFT normalisation | `out_t3_lightcone_l1.txt` |
| `t4_sign_l1.py` | physical sign of the opposite-sign term | `out_t4_sign_l1.txt` |
| `t5_helicity_l2.py` | `⟨TTO⟩` helicity map, Weyl parity of the MFT normalisation | `out_t5_helicity_l2.txt` |
| `t6_kappa_l1.py` | calibration `κ_h A_h` for currents | `out_t6_kappa_l1.txt` |
| `feas_l2.py` | feasibility of the `⟨TTTT⟩` word pairs | `cache/feas_l2.log` |
