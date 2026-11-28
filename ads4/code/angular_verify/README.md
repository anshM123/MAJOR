# Hostile check of the all-k angular rules and the τ ≲ C/Δg³ claim

This directory checks the claim in `code/angular_opt/` independently. That directory has no README; the claim was taken from the task statement and the docstrings. Every script here was written from scratch and imports only `common.py` from this directory. The claimant's code was read only to learn its conventions, and its numbers are compared at the end. Grades: [P] proved here, [D] derived here, [DV] derived and numerically verified, [NE] numerical evidence only, [failed].

Scripts (run from this directory; every one takes under about 3 min):

| script | what it does |
|---|---|
| `common.py` | ladder ratios (supplement.tex, sec:angular item (b)), endpoint weights, admissible (b,J) |
| `v1_freefield.py` | free scalar/fermion test of R_k: signs, normalisation, T only at k=2 |
| `v2_lp.py g N` | first LP (dual with constraint generation, double precision) and the exact 3-state primal |
| `v3_lp_hp.py g N [light] [odd_b0] [contact\|contact3]` | well-conditioned LP: a basis evaluated in 256-bit arithmetic, constraint generation over J ≤ 20000 and twist up to 10⁵g, exact re-check of the certificate |
| `v4_threestate.py` | closed form of the 3-state optimum (sympy) and its asymptotics |
| `v5_primal_diag.py g N [light]` | which states saturate the LP |
| `v6_light_action.py g N` | size of the optimal dual on light (twist < g) states |
| `v7_racah_check.py g n` | exact check, on a grid, of the claimant's closed-form Racah/4F3 dual |
| `v8_ladders.py` | ladder bookkeeping and the "two ladders agree at leading order" statement |
| `v9_light_family.py` | optimum when the J=2 state is a light spin-2 single trace |
| `v10_racah_light.py g n` | Racah dual on light spin-2 single traces |
| `v11_light_tradeoff.py g N` | how small the light remainder E_light must be |

## (a) Derivation of the k-rules

Take the state |T,λ=+2⟩ at rest. At rest T has only the SO(2) charges ±2, because conservation removes T_{0ν} and tracelessness removes the spatial trace. Define E_m = ∫dφ e^{imφ}E(φ)/E. Two facts fix the structure:

- E(φ) commutes with P.
- Rotations act as E(φ) → E(φ+θ), so ⟨O,q|E_m|T,2⟩ ≠ 0 only if q = 2+m.

Commutativity of the E(φ) as a distribution on S¹×S¹, with no contact term on the diagonal, gives ⟨T|[E_m,E_{−m}]|T⟩ = 0. Insert the complete set of states at rest. This is Parseval for the vector E_{∓m}|T⟩: one fixed m, so there is no exchange of a mode sum with a spectral sum. With m = k+2:

- E_{+m}: q = 4+k, the b=4 ladder at index k (needs J ≥ 4+k).
- E_{−m}: q = −k, the b=0 ladder at index k (needs J ≥ k; even J only, see the structure table at supplement.tex l.281–285 and sec_fullpos.tex l.40–46).
- T can appear only at q = ±2, that is only for E_{−4} (k=2), on the b=0 side. Its weight there is |⟨T,−2|E_{−4}|T,2⟩|² = (t4²+α_T²)/16. This uses the 3d collider normalisation ⟨E⟩ = (1/2π)(1 + (t4/4)cos4φ + …).

Hence **R_k: Σ_O|a_{O,4,k}|² − Σ_{O≠T}|a_{O,0,k}|² = δ_{k2}|c|², with |c|² = τ²/16**. This is exactly the claimed rule, with the claimed sign and normalisation. The paper's k=2 rule is this identity with Σ|a_{0,2}|² ≥ 0 dropped. [D]

- **Free fields** (`v1`): E_m acts as multiplication by e^{imφ}(1+(−1)^m)/2. For every k: LHS = RHS = 1 for even k ≠ 2, everything vanishes for odd k, and for k=2 the spin-6 current gives 1 = |c|² (t4 = ±4). The endpoint normalisation 3/2 + 3/2 = 3 matches the same |a|². [DV]
- **Large m.** Commuting positive E(φ) have a joint spectral measure, so the energy-energy correlator (EEC) is a positive measure. Then ‖E_m|T⟩‖² = ∫∫e^{im(φ₁−φ₂)}dμ ≤ ⟨E_0²⟩ = 1 uniformly in m. Every Fourier mode is a legitimate smeared operator, and no condition appears beyond the one already used for m=4. [P given the POVM]
- **What is exact.** R_k is exact as an identity for amplitudes, given commutativity with no contact term. That is the same input as the paper's [E_4,E_{−4}] = 0, and it needs KKSDZ J₀ < 3 at separated points. Writing R_k in terms of a_{O,b,0} uses the ladder identity for every k up to N. The paper's telescoping already uses the full ladder, so this is not a new input. The LP also uses the endpoint rule, so it inherits 𝒟 = 0 (sec:angcross), the first-moment Ward identity and the c₀ contact coefficient. R_k itself has no 𝒟-type defect.
- **Contact-term sensitivity** (`v3 … contact`). If an unknown antisymmetric contact term x(k+2) is added to every R_k (and x′(k+2)³ at g=10), the LP optimum does not change at g=10, 20, 40 [NE]. The cubic variant failed numerically at g=20.

**Light states.** These are the analogue of 𝓛_g and are the main correction to the claim.

1. **Light single traces of spin 0 or 2.** Letter assumption (i) allows them. They have O(1) normalised amplitudes: √c_T × λ_{TTφ} ~ O(1), and A2′ excludes them. They enter only the b=0 side, at k ≤ J ≤ 2, with positive endpoint cost W.
   - For the paper's k=2 rule they are harmless, because that rule needs J ≥ 6.
   - For the all-k LP they must be included, and the claimant's LP omits them. Including them (`v3 … light`, `v9`) raises the LP optimum. A J=2 state below the gap, at D_A ≈ 0.76g, replaces state A. τg³ becomes:

     | g | without light J=2 | with light J=2 |
     |---|---|---|
     | 20 | 397.2 | 411.2 |
     | 40 | 478.4 | 520.4 |
     | 1000 | 582.8 | 684.9 |
     | ∞ | 587.9 | 693.6 |

     The values at g=20 and g=40 are the LP optimum [NE]. The large-g column is the exact family value, so it is a lower limit on the LP.
   - The claimant's Racah dual is infeasible on light spin-2 states at g=100 (margin −0.14%) and g=200 (−1.7%) (`v10`).
   - So the constant moves (588 → ≳ 694) but the scaling does not.
2. **Light multi-traces (A2′).** They enter every R_k on both ladders. Per state, |a_k|²/|a_0|² ≤ C(k+3,3)²·(1+o(1)). At fixed g and c_T → ∞ their contribution vanishes because E_light → 0 [D given A2′].
   - The coefficient is enormous. The optimal dual has σf − λW per unit |a_0|² on twist-2 double twists of about 1.8e3 (g=10), 5.3e6 (g=20) and 8.3e17 (g=40) (`v6`). The Racah dual at g=50 gives 1.4e22, against L = 1.4e6. The paper's k=2 coefficient is ≤ 101.
   - After re-optimising, the LP keeps its value only if E_light ≲ 3e-9 (g=20) or ≲ 3e-20 (g=40). At E_light ≤ 3e-6 (g=40), τg³ is 3087 instead of 478 (`v11`). With no A2′ there is no bound (τ ≤ 10.4 at g=20).
   - The paper's angular bound has "error terms with no large constants" (Letter l.127/138). The all-k bound loses that property. Its light remainder behaves like that of the positivity route, with a constant that appears to grow exponentially in g.
3. **Spins 1 and 3** do not occur in T×T. b=0 at odd J is absent per the paper's structure table, so the `odd_b0` variant (+4% at g=10, none at g=20, 40) is moot.

## (b) Ladders

The LP's ladder is the paper's a_{k+1}/a_k = ((k+4)/(k+1))·√((C−s(s+1))/(T_Δ−s(s+1))), squared, with s = b+k. Its closed form is r_k = C(k+3,3)²(b−J)_k(J+b+1)_k/((b+2−Δ)_k(Δ−1+b)_k). It reproduces the paper's |a_2|²/|a_0|² = 100(C−20)(C−30)/((T−20)(T−30)) exactly (`v8`). [DV]

I did **not** re-derive the ladder identity. It remains the paper's input: [P given light-transform conservation as an operator identity], pipeline-checked for J ≤ 40. The only check here is consistency at twist 1: the ladder has a pole at the last step, which matches free currents having only q = ±J.

"On a heavy state the two ladders agree to leading order in 1/T" is **false as stated**. Both scale as T^{−k}, but at fixed J their coefficients differ by Π_j (J−4−j)(J+5+j)/((J−j)(J+1+j)): at J=6, r₁⁽⁴⁾/r₁⁽⁰⁾ → 0.52 and r₂⁽⁴⁾/r₂⁽⁰⁾ → 0.16. They agree only for J ≫ k. This is harmless, because the LP uses the exact ratios.

## (c) The LP solved independently

The LP is: max |c|² subject to Σ W ρ ≤ 3, R_k for k = 0..N, and ρ ≥ 0, over heavy states with twist ≥ g. It is solved through its minimax dual with constraint generation.

| g | N | LP |c|² | 3-state |c|² | τ g³ | exact re-check |
|---|---|---|---|---|---|
| 10 | 2 / 3 / 4 | ratio 19.8 / 5.39 / 2.58 | 5.0428e-3 | 1264 / 660 / 456 | |
| 10 | ≥ 6 | 5.04283545e-3 | 5.04283545e-3 | 284.05 | 1e-15 |
| 20 | 4 / 8 | ratio 15.6 / 2.05 | | | |
| 20 | ≥ 12 | 1.54049908e-4 | 1.54049908e-4 | 397.17 | 3e-11 |
| 40 | 16 | ratio 2.14 | | 700 | |
| 40 | 32 | 3.49226380e-6 | 3.49226380e-6 | 478.40 | 1e-15 |
| 80 | 48, 64 | [failed numerically] | 6.6644e-8 | (528.70) | violations 6–30% |

- From N ≈ 0.6–0.8g on, the LP optimum equals the 3-state value to 9 digits [NE: the certificate is exact on the binding states and checked in floating point on 2–6×10⁴ states].
- **3-state configuration** (`v4`, exact): A = (b0, J=2, g+2), B = (b4, 4, g+4), C = (b4, 6, g+6).
  - It satisfies every R_{k≠2} exactly and saturates the budget, with ρ_B ∝ 8(g−2)/(3(g+9)) > 0 for g > 2.
  - |c|² = 475200(g−2)/[g(g−1)(g+2)(g+10)(22g³+100g²+249g+411)].
  - Hence |c|²g⁶ → 21600 and τg³ → 4√21600 = 587.88. [P]
  - So 588 is the g → ∞ limit. At finite g the LP lower limit is smaller: τg³ = 284 at g=10 and 582.8 at g=1000, where τ = 5.83e-7. This matches the claimed "≈ 6e-7".
- **Racah dual** (`v7`, exact on a grid with J < 6g plus 5 larger J and 14 twists):
  - Feasible at g=50 for n ≤ 28 (τg³ = 730), at g=100 for n ≤ 53 (834), and at g=200 for n = 104 (877).
  - At larger n it is infeasible at (b4, J=7, gap).
  - This is consistent with the claimed ≈ 941/g³ asymptote and n ≈ 0.53–0.56g. Feasibility for all (J, Δ) is unproven.

## (d) Loopholes

- **Large-m modes.** No problem; see the Bochner bound in (a).
- **Extra assumptions for higher harmonics.** None of a new kind. The quantitative demand on A2′ is very much stronger, as shown in (a) item 2.
- **"No light single trace".** The bound does not need this condition to exist. But light spin-0 and spin-2 single traces must be included as LP columns, and doing so changes the constant to ≳ 694/g³.
- **Unphysical optimum.** The optimum is three operators exactly at the gap, each with a single structure: J=2 with ++ only, and J=4 and J=6 with +− only. The per-state structure counting allows this. That the configuration is unlikely to be crossing-consistent only means the bound may not be sharp. Validity of an upper bound needs the dual only.
- **Gap definition.** The letter defines the gap as twist ≥ Δg for spin ≥ 4 single traces. The LP uses that definition for all heavy states.

## Verdict

- **R_k.** Correct, with the claimed form, signs and normalisation [D, DV in free fields]. They hold under the same commutativity input as the paper's k=2 rule plus the ladder identity, and the paper already uses the ladder identity for every k.
- **C/Δg³ scaling.** The LP built from these rules has optimum ~g⁻³:
  - lower limit 587.9/g³: exact [P];
  - equality with that lower limit at g = 10, 20, 40: [NE];
  - an upper ~g⁻³ dual exists on grids up to g = 200 [NE].

  No all-(J, Δ) dual certificate exists, so **no bound beyond the paper's (160/9)/Δg is proved**.
- **Errors and omissions in the claim:**
  - (i) Light spin-0/2 single traces are omitted. They are allowed by the gap assumption and raise the constant to ≳ 694. The Racah dual as given is slightly infeasible on them for g ≥ 100.
  - (ii) The light-remainder coefficient grows enormously (about 1e18 at g=40). The improved bound holds only as a strict c_T → ∞ statement and loses the "no large constants" feature of the paper's angular bound.
  - (iii) "Ladders agree to leading order in 1/T" is false at fixed J. The LP does not use it.
  - (iv) 588 is an asymptotic constant: the finite-g LP values are 284 to 583 for g = 10 to 1000.
