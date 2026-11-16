# Energy-flux kernels for `T × T → O` in d = 3: `â^O_{m,2}` through OPE data

Status: working document, 16 Nov 2026. Everything below is a computation by one agent, **not refereed**.
Markers: **[FF]** = verified against free fields (Fock space, §2) by the full position-space pipeline;
**[C]** = computed exactly (symbolically in `Δ`) by the pipeline of §3, with the null tests listed there, but with
no independent free-field test of that particular statement; **[S]** = from a source, specialised to d = 3
(equation numbers given) and compared with the pipeline; **[I]** = inferred (pattern in `J`, numerics, or
physical argument; not derived). Code: `ads4/code/kernels/` (map in §9; `run_all.sh` reproduces every number).

Background: `EINSTEIN_PROGRAM.md` §3 (Proposition S, (★)), `EINSTEIN_OBSTRUCTIONS.md` §3.1–3.4, ST1. The present
file does ST1 (a)–(d) at the level of derivation + exact symbolic computation; the hypotheses (E) of
Proposition S (commutativity without contact terms, mode-wise convergence) are **assumed, not examined**.

---

## 0. Summary

1. **Formalism validated on free fields.** Position-space Wightman structure → detector limit → `α`-integral →
   Fourier transform reproduces, for the free scalar and the free Majorana fermion: the Ward identity
   `â^T_{2,2} = 1`, `c = ±1` (`t_4 = ±4`), `|â^{φ²}_{0,2}| = |â^{ψ̄ψ}_{0,2}| = 1`, and
   `|â^{J_s}_{±s,2}| = 1`, `â^{J_s}_{m,2} = 0` for `|m| < s` (`s = 4, 6`), all exactly (§2, §4).
2. **Selection rules and universal `m`-dependence** (all conserved `⟨TTO⟩` structures, `J ≤ 12`, symbolic in
   `Δ`; §5.1). For the source `|T, +2⟩`:
   * `â^O_{m,2} = 0` for `m = 1, 2, 3` (`O ≠ T`);
   * the components `m ≤ 0` are governed by **one complex number** `a^{++}_O := â^O_{0,2}` (even `J` only), the
     components `m ≥ 4` by **one more number** `a^{+−}_O := â^O_{4,2}` (`J ≥ 4`), with
     `|â_{−k,2}|² = |a^{++}_O|² W⁰_k(Δ,J)`, `|â_{4+k,2}|² = |a^{+−}_O|² W⁴_k(Δ,J)`,

     `W⁰_k = C(k+3,3)² (J+1−k)_{2k} / (Δ−1−k)_{2k}`,  `W⁴_k = C(k+3,3)² (J−3−k)_k (J+5)_k / ((Δ−5−k)_k (Δ+3)_k)`.
   * `a^{++}` receives the parity-even "`++`" structure and (with relative phase `i`) the parity-odd
     structure; `a^{+−}` receives the parity-even "`+−`" structure (even `J ≥ 4`) or the single structure of
     odd `J ≥ 5`. This is the CFT₃ image of the 4D helicity couplings `c^{++}_J`, `c^{+−}_J`; the `J`-dependent
     factors of `W⁰`, `W⁴` are exactly the coefficients of `₂F₁(−J, J+1; 1; y)` and `₂F₁(4−J, J+5; 1; y)`.
3. **Kernels in closed form** (§5.2–5.4): scalar (agrees exactly with KKSZ (5.111)/(5.132) and CMT
   (4.10)–(4.11) at d = 3), spin 2, and canonical `++`, `+−`, odd structures for all `J` (closed forms fitted on
   `J ≤ 8`, confirmed out of sample at `J = 9, …, 12` to 30 digits).
4. **Sum rules** (§6). Commutativity is equivalent to the family, `k = 0, 1, 2, …`,

   `(S_k)  δ_{k,2} (t_4² + α_T²)/16 + Σ_{O≠T, J even} |a^{++}_O|² W⁰_k(Δ_O,J_O) = Σ_{O, J ≥ 4+k} |a^{+−}_O|² W⁴_k(Δ_O,J_O) = a_{2+k}`,

   where `a_n` are the Fourier coefficients of the normalised energy–energy correlator in `|T,2⟩`. `(S_2)` is (★);
   `(S_0)` is `Â_0 = Â_4`. The "off-diagonal" relations (O) of `EINSTEIN_OBSTRUCTIONS.md` §3.3 hold identically.
5. **Large `Δ`** (§5.6): `|K_6|²/|K_4|² = W⁴_2 = 100 (J−5)(J−4)(J+5)(J+6)/((Δ−7)(Δ−6)(Δ+3)(Δ+4))` exactly,
   `|K_{−2}|²/|K_0|² = W⁰_2 = 100 (J−1)J(J+1)(J+2)/((Δ−3)(Δ−2)(Δ−1)Δ)`; both `∼ 100 J⁴/Δ⁴` at fixed `J`, both `→ 100 x⁴`
   for `J = xΔ → ∞`. `K_6` and `K_{−2}`, `K_0` belong to different structures (no relation between them).
6. **Double-twist zeros** (§5.7): every squared kernel has a double zero at every `N = ∞` dimension of `[TT]`
   with the corresponding quantum numbers; double-trace contributions to both sides of every `(S_k)` are
   `O(1/c_T)`.
7. **What the sum rules imply** (§7): (i) `|c|² ≤ a_2 · sup_{J_O ≥ 6} W⁴_2(Δ_O, J_O)` — a gap-dependent bound
   *if* `J_O/Δ_O` is bounded for the spin-≥ 6 operators, none otherwise; (ii) an explicit exact four-operator
   solution of all `(S_k)` with positivity has `|c|² ≃ 2700/Δ_gap⁴`, so these constraints can never give a
   rate better than `Δ_gap^{−4}` for `t_4² + α_T²` (the expected physical rate is `Δ_gap^{−8}` up to logs);
   (iii) a linear programme over all operators (spins ≤ 100) finds maxima `≈ 5·10³/Δ_gap⁴`, i.e. numerical
   evidence — **not a proof** — that the full family `(S_k)` plus positivity of the energy correlator does bound
   `|c|²` by a quantity decaying like `Δ_gap^{−4}`. The statement "homogeneous, hence no rate" of
   `EINSTEIN_OBSTRUCTIONS.md` §3.4 is therefore too pessimistic as stated, and the statement that a rate
   `Δ_gap^{−8}` could follow is excluded.

---

## 1. Conventions

* Metric `(−,+,+)`. `q = (E,0,0)`; detector direction `n(θ) = (cos θ, sin θ)`, `z = (1, n)`.
  `e_± = (0,1,±i)`, `e_0 = (1,0,0)`; a null polarisation is `ζ(a,b) = a² e_+/2 + b² e_−/2 + ab e_0`.
* A spin-`J` primary is encoded as `O(x,ζ) = O_{μ_1…μ_J} ζ^{μ_1}⋯ζ^{μ_J}`. Momentum states
  `|O(p),ζ⟩ = ∫d³x e^{ip·x} O(x,ζ)|0⟩`; components `O(p;ζ(a,b)) = Σ_m a^{J+m} b^{J−m} O_m(p)`. With
  `U(R_φ) = e^{−iφ J_z}` one has `U(R_φ) O_m(q)|0⟩ = e^{−imφ} O_m(q)|0⟩` (charge `m`). Helicity states of the stress
  tensor: `|T,±2⟩ ∝ T_{μν} e_±^μ e_±^ν (q)|0⟩`.
* Two-point normalisation: `⟨O(x,ζ)O(0,ζ′)⟩ = (ζ·I(x)·ζ′)^J/(x²)^Δ`, `I_{μν} = η_{μν} − 2x_μx_ν/x²` (unit);
  `⟨T(x,z)T(0,z′)⟩ = C_T (z·I·z′)²/x⁶`. Free real scalar and free Majorana fermion: `C_T = 3/(32π²)`
  (scalar: computed, §4; fermion: input).
* Wightman functions: operators to the left have `t → t − iε`, i.e. `x_{ij}² = −(t_{ij} − iε)² + x⃗_{ij}²` for
  `i` left of `j`.
* Energy flux: `E(n) = lim_{r→∞} r ∫dt n^i T^0{}_i(t, r n)`, so that `∫dθ E(θ) = H` and
  `∫dθ e^{∓iθ}E(θ) = P^x ∓ iP^y`. `U(R_φ)E(θ)U(R_φ)^{−1} = E(θ+φ)`.
* Matrix elements: `⟨O_m(p′)|E(θ)|T_h(p)⟩ = (2π)³ δ³(p−p′) M_m(θ)`, `M_m(θ) = e^{i(h−m)θ} M_m(0)`;
  norms `⟨O_m(p′)|O_m(p)⟩ = (2π)³δ³(p−p′) N_{O,m}`. Normalised amplitudes (rest frame, `E = 1`)

  `â^O_{m,h} := 2π M_m(0) / (N_{O,m} N_{T,h})^{1/2}`,

  so that `â^T_{2,2} = 1` is the Ward identity `∫E = H`. Phases of `â` depend on the phase conventions for the
  states; only `|â|²` and the relative phase of the parity-even and parity-odd contributions are used.
* `Â_m := Σ_O |â^O_{m,2}|²`. Normalised energy–energy correlator in `|T,2⟩`:
  `⟨E(θ_1)E(θ_2)⟩/⟨1⟩ = (E/2π)² 𝓔(θ_{21})`, `𝓔(θ) = Σ_n a_n e^{inθ}`, `a_n = Â_{2−n}`. Commutativity ⟺ `a_n = a_{−n}`
  ⟺ `Â_m = Â_{4−m}`. (In `EINSTEIN_PROGRAM.md` §3 the labels `a_{4}`, `a_{−4}` are interchanged relative to this
  convention, with `U = e^{+iφJ}`; nothing depends on it.)
* OPE data: `⟨0|O(x_4,ζ) T(x_1,z_1) T(x_3,z_3)|0⟩ = Σ_a λ^{(a)}_{TTO} S_a` with `O` unit-normalised, `T` canonical,
  `S_a` the conserved, Bose-symmetric structures of §3.1. Then

  `â^O_{m,2} = C_T^{−1/2} Σ_a λ^{(a)}_{TTO} K^{(a)}_m(Δ,J)`,  `K^{(a)}_m := 2π M^{(a)}_m(0) / (N_{O,m} π²/3)^{1/2}`,

  where `M^{(a)}_m` is computed with the source `T(e_+)` and `π²/3 = 16 N_{J=2,m=2}(Δ=3)` is the norm of
  `T(e_+)|0⟩` per unit `C_T`. `λ` real for real structures (Hermitian `O`).
* `s := (Δ + J)/2` throughout; `(a)_n` is the Pochhammer symbol.

---

## 2. Free fields in Fock space  (`fock_free.py`)

### 2.1 Kinematics common to both theories

`T(q)|0⟩` and every bilinear `O(q)|0⟩` are two-particle states. With `q = (E,0,0)` and massless particles,
`k_1 = (E/2)(1, n(θ))`, `k_2 = (E/2)(1, −n(θ))`: the two particles are back to back, each of energy `E/2`. The
two-body phase space is `∫dθ/(8πE)`, so a state is a wavefunction `ψ(θ)` on the circle, with
`⟨ψ′|ψ⟩ = 𝒩 ∫dθ ψ̄′ψ` (`𝒩` cancels in all normalised quantities), and

`E(θ′)|θ⟩ = (E/2)[δ(θ′−θ) + δ(θ′−θ−π)]|θ⟩`,   `|θ⟩ := |k_1(θ), k_2(θ)⟩`.

For bosons `|θ+π⟩ = |θ⟩`; for fermions `|θ+π⟩ = −|θ⟩` and the spinor `u(k_2)` continued along `θ → θ+π` changes sign,
so that in both cases wavefunctions are `e^{imθ}` with **even** `m`, each once. `U(R_φ)|θ⟩ = |θ+φ⟩`, hence
`ψ = e^{imθ}` has charge `m`. Consequently, with `|m⟩ := ∫dθ e^{imθ}|θ⟩`:

`⟨m′|E(θ)|m⟩ / (⟨m′|m′⟩⟨m|m⟩)^{1/2} = (E/2π) e^{i(m−m′)θ}` for `m − m′` even (0 otherwise),

`⟨m|E(θ_1)E(θ_2)|m⟩/⟨m|m⟩ = (E²/4π)[δ(θ_{12}) + δ(θ_{12} − π)]`.

So **every** normalised amplitude between two-particle states of even charge difference has modulus one, and
the energy–energy correlator is `𝓔(θ) = π[δ(θ) + δ(θ−π)]`: `a_n = (1 + (−1)^n)/2`.

### 2.2 Free real scalar

`φ(x) = ∫d²k/((2π)² 2|k|) (a_k e^{ik·x} + h.c.)`. Improved stress tensor with null `z`:
`T(x,z) = (3/4)(z·∂φ)² − (1/4) φ (z·∂)²φ`. Creation part: `∂ → −ik`, wavefunction
`F_T = −(3/4)(z·k_1)(z·k_2) + (1/8)((z·k_1)² + (z·k_2)²)`. With `z = ζ(a,b)`:
`z·k_1 = (E/4)(A−B)²`, `z·k_2 = −(E/4)(A+B)²`, `A = a e^{iθ/2}`, `B = b e^{−iθ/2}`, hence

`F_T = (E²/16)(a⁴ e^{2iθ} + b⁴ e^{−2iθ})`:  `|T,±2⟩ = (E²/16)|±2⟩`, no other components (conservation, tracelessness).

Currents: `J_s(x,z) = Σ_k (−1)^k C(2s,2k) (z·∂)^{s−k}φ (z·∂)^kφ`, i.e. `Re(√(z·∂_1) + i√(z·∂_2))^{2s}`; wavefunction
`F_s = (±E^s/2)(a^{2s} e^{isθ} + b^{2s} e^{−isθ})`: only `m = ±s` (checked for `s = 0, 2, 4, 6, 8`). Within bilinears,
"only `m = ±s`" is equivalent to "conserved primary": a descendant `(z·∂)^{s−s′}J_{s′}` has wavefunction
`(z·q)^{s−s′}F_{s′} ∝ (ab)^{s−s′}(…)`, which has charges `±s′`.

Results: `⟨T,2|E(θ)|T,2⟩/⟨T,2|T,2⟩ = E/2π`; `⟨T,−2|E(θ)|T,2⟩/norm = (E/2π) e^{4iθ}`, i.e. **`c = +1`, `t_4 = +4`**;
`â^{φ²}_{0,2} = 1`; `â^{J_s}_{±s,2} = (−1)^{s/2}` for `s = 4, 6, 8` (both signs of `m` equal); all other components
vanish. `Â_m = 1` for every even `m`, each produced by exactly one operator:
`m = 2`: `T`; `m = −2`: `T` (this is `|c|² = 1`); `m = 0`: `φ²`; `m = ±s`: `J_s`. (★) reads `1 + 0 = |â^{J_6}_{6,2}|² = 1`.

### 2.3 Free Majorana fermion

`γ^0 = iσ_2, γ^1 = σ_1, γ^2 = σ_3` (real), `ψ(x) = ∫d²k/((2π)²2|k|)(u(k) b_k e^{ik·x} + h.c.)`, `k̸u = 0`,
`u(k) = √(2ω)(sin(θ/2), cos(θ/2))`. `T(x,z) ∝ ψ̄ z̸ (z·∂)ψ`; `γ^0γ^μ` are symmetric, so the antisymmetrised two-particle
wavefunction is `F_T = (1/4)(u_1^T γ^0 z̸ u_2)(z·k_2 − z·k_1) = (iE²/16)(a⁴ e^{2iθ} − b⁴ e^{−2iθ})`.
Currents (even `s`): `J_s = Σ_k (−1)^k C(2s,2k+1) (z·∂)^{s−1−k}ψ̄ z̸ (z·∂)^kψ` (`Im(χλ_1 + iχλ_2)^{2s}` in spinor variables);
`F_s ∝ a^{2s}e^{isθ} − b^{2s}e^{−isθ}` (checked `s = 4, 6`). `ψ̄ψ` (`Δ = 2`, parity odd): `F = u_1^Tγ^0u_2 = −E`, charge 0.

Results: `⟨T,2|E|T,2⟩/norm = E/2π`; **`c = −1`, `t_4 = −4`**; `|â^{ψ̄ψ}_{0,2}| = 1`; `â^{J_s}_{s,2} = −â^{J_s}_{−s,2}`,
modulus 1 (`s = 4, 6`); the relative sign between the `m = +s` and `m = −s` components is opposite to the
scalar. Energy–energy correlator and `a_n` as in §2.1.

**Answer to the question of `EINSTEIN_PROGRAM.md` §4.1:** `|â^{J_s}_{±s,2}|² = 1` is true in both free theories,
for every even `s ≥ 4`, and `|â^{φ²}_{0,2}|² = |â^{ψ̄ψ}_{0,2}|² = 1`.

---

## 3. Formalism  (`structures.py`, `pipeline.py`)

### 3.1 Structures and conservation

Embedding space `X = (X^+, X^-, x^μ) = (1, x², x^μ)`, `Z = (0, 2x·z, z^μ)`, `X·Y = −(X^+Y^- + X^-Y^+)/2 + x·y`. Points:
`4 ↔ O`, `1 ↔` detector `T`, `3 ↔` source `T`. Building blocks (Costa–Penedones–Poland–Rychkov):
`V_1 = V_{1,34}`, `V_3 = V_{3,41}`, `V_4 = V_{4,13}`, `V_{i,jk} = (Z_i·X_j X_i·X_k − Z_i·X_k X_i·X_j)/(X_j·X_k)`;
`H_{ij} = −2[(Z_i·Z_j)(X_i·X_j) − (Z_i·X_j)(Z_j·X_i)]`; `X_{ij} = −2X_i·X_j`;
`ε_{13} = det[Z_1,Z_3,X_1,X_3,X_4]`, `ε_{14} = det[Z_1,Z_4,X_1,X_3,X_4]`, `ε_{34} = det[Z_3,Z_4,X_1,X_3,X_4]`.

* parity even: `V_1^{m_1}V_3^{m_3}V_4^{m_4}H_{13}^{n_{13}}H_{14}^{n_{14}}H_{34}^{n_{34}} / (X_{13}^{5−s} X_{14}^{s} X_{34}^{s})`;
* parity odd: `ε_{ij}X_{ij} × (monomial with the `Z`-degrees at `i`, `j` lowered by one) / (X_{13}^{11/2−s} X_{14}^{s+1/2} X_{34}^{s+1/2})`
  (the factor `X_{ij}` and the half-integer shifts restore homogeneity at the third point; for `J = 0` this is
  CMT (D.3), and the conservation condition `β_1(Δ−3) = β_2(Δ+1)` quoted there is reproduced).

Linearly independent monomials are selected numerically (exact rational points); conservation
`∂_{x_1}·D_{z_1} S = 0` (Todorov operator `D^μ_z = (1/2 + z·∂_z)∂^μ_z − (1/2)z^μ∂_z²`) and symmetry under
`(x_1,z_1) ↔ (x_3,z_3)` are imposed at random rational points, symbolically in `s`. **Counts found [C]** (they
agree with 1708.05718 §2.1 as quoted in `EINSTEIN_OBSTRUCTIONS.md`):

| `J` | independent monomials (even / odd) | conserved + symmetric (even / odd) |
|---|---|---|
| 0 | 3 / 2 | 1 / 1 |
| 2 (generic `Δ`) | 10 / 9 | 1 / 1 |
| 2, `Δ = 3` (`O = T`) | 10 / 9 | 2 / 1 |
| even `J ≥ 4` | 13 / 12 | 2 / 1 |
| 1, 3 | — | 0 / 0 |
| odd `J ≥ 5` | — | 0 / 1 |

(13 = 14 − 1: one Gram-determinant relation in d = 3.) Conserved structures were computed for all `J ≤ 12`.

### 3.2 Detector limit, `α`-integral, Fourier transform

**(F1)** `E(n) = 2∫dα S(X_1 = (0,−α,z), Z_1 = (0,−1,0); …)`, `z = (1,n)`. Derivation: `x_1 = (v/2)z + (u/2)z̄`,
`z̄ = (1,−n)`; at large `v` the flux is `n^iT^0{}_i = T_{μν}z̄^μz̄^ν/4`, and homogeneity gives
`T(x_1,z̄) = (2/v)³(2v)² S((0,−2u,z),(0,−1,0))`; `r = v/2`, `α = 2u`. This is `E = 2L[T](∞,z)` of KKSZ (3.18) **[S]**.
In this limit, with `x := x_4 − x_3`, `P_3 = α + 2z·x_3`, `P_4 = α + 2z·x_4` (verified against the 5×5 determinants and
the embedding expressions at random points):

`V_1 = −z·x/x²`, `V_3 = ((z_3·x)P_3 + (z·z_3)x²)/P_4`, `V_4 = ((ζ·x)P_4 − (z·ζ)x²)/P_3`, `H_{13} = z·z_3`, `H_{14} = z·ζ`,
`H_{34} = (z_3·ζ)x² − 2(z_3·x)(ζ·x)`, `X_{13} = −P_3`, `X_{14} = −P_4`, `X_{34} = x²`,
`ε_{13} = ε(z_3,z,x)`, `ε_{14} = ε(ζ,z,x)`, `ε_{34} = 2(ζ·x)ε(z_3,z,x) + P_3 ε(z_3,ζ,x) + x² ε(z_3,ζ,z)`.

**(F2)** With the Wightman prescription `X_{14} = −(α + 2z·x_4 + iε)`, `X_{13} = −(α + 2z·x_3 − iε)`:

`∫dα X_{14}^{−a} X_{13}^{−b} = −2πi e^{iπa} Γ(a+b−1)/(Γ(a)Γ(b)) · (2z·x + iε)^{1−a−b}`.

(`a + b − 1` is always an integer `≥ 2` here.)

**(F3)** For `p` in the forward cone,

`∫d³x e^{−ip·x} (2z·x + iε)^{−A} (x²)^{−B} = e^{−iπA/2} π^{5/2} 2^{4−2B} / (Γ(B)Γ(A+B−1/2)) · (−p²)^{A+B−3/2} (−4z·p)^{−A}`,

and `x^μ → i∂/∂p_μ` for the polynomial prefactors `(z_3·x)`, `(ζ·x)`, `ε(·,·,x)` (done exactly by expanding
`(−(p+W)²)^γ(−4z·(p+W))^{−A}` in the auxiliary vector `W`). Then `p = (1,0,0)`, `z_3 = e_±`, bra polarisation
`ζ̄(a′,b′) = ζ(b′,a′)`, and the coefficient of `a′^{J+m}b′^{J−m}` is `M_m`.

Each term carries `1/(Γ(a)Γ(b)Γ(B)Γ(A+B−1/2))` with `a ∈ s + ℤ`, `b ∈ 5 − s + ℤ`, `B ∈ s + ℤ` (shifted by 1/2 for odd
structures), so `M_m = G(s) P_m(s)` with `P_m` a polynomial and
`G = −4πi e^{iπs} π^{5/2} 2^{4−2s} / (Γ(s+K_a)Γ(5−s+K_b)Γ(s+K_B)Γ(s−1/2+K_γ))` (integers `K` printed by the code).

**Check of (F1)–(F3) and all normalisations [C]:** for `⟨φ(x_4) T(x_1) φ(x_3)⟩ = λ V_1²/(X_{13}^{5/2}X_{14}^{5/2}X_{34}^{Δ−5/2})`
the chain gives `2πM = −λ (8π/3Δ) ⟨φ(p)|φ(p)⟩`, so `∫E = H` requires `λ = −3Δ/(8π) = −dΔ/((d−1)S_d)`, the known Ward
identity value. (Hand computation; the same formulas are used by the code.)

### 3.3 Norms  **[C]**

`N_{O,m}(Δ,J) = π^{5/2} 2^{4−2Δ−J} C(2J, J+m) (Δ−J−1)_{J−|m|} (Δ−1)_{|m|} / (Γ(Δ+J) Γ(Δ−1/2))`,  `|m| ≤ J`.

Derived symbolically for `J ≤ 5` (`norm_poly`), used as a formula for higher `J`. Positive for `Δ > J + 1`; at
`Δ = J + 1` only `m = ±J` survive. `m`-dependence: `C(2J,J+m)(Δ−J−1)_{J−|m|}(Δ−1)_{|m|}`.

### 3.4 Null tests passed by every conserved structure, `J ≤ 12`, symbolic in `Δ`  **[C]**

* `M_m = 0` for `m = 1, 2, 3` (`O ≠ T`): `∫E = H`, `∫e^{∓iθ}E = P^x ∓ iP^y`. Individual monomials do not satisfy this.
* Source transversality: with a general null source polarisation `ζ(1,d)` only `d⁰` and `d⁴` appear (`J = 0, 2, 4`).
* `M_m(θ) = e^{i(2−m)θ}M_m(0)` at `e^{iθ} = (3+4i)/5` (`J = 4`, both parities).
* `M^{h=−2}_m = ± M^{h=+2}_{−m}` (`+` parity-even structures, `−` parity-odd ones, including odd `J`).
* For `O = T` (`s = 5/2`): `M_{±1} = M_0 = 0`; the parity-odd structure has `M_2 = 0`.

---

## 4. Free-field validation of the pipeline  (`wick_free.py`, `free_checks.py`, `free_canonical.py`)  **[FF]**

Position-space correlators `⟨O T T⟩` are computed by Wick contraction at rational points (scalar:
`G = 1/(4π|x|)`; fermion: Dirac field = two Majorana fields, `S ∝ x̸/|x|³`, normalisation fixed by
`C_T^{Dirac} = 3/(16π²)`), expanded in the monomials of §3.1 (overdetermined exact linear systems, zero
residual), and pushed through §3.2. The two-point functions give `C_T = 3/(32π²)` for the scalar and the unit
normalisations of `φ²`, `ψ̄ψ`, `J_4`, `J_6`.

| theory | `O` | expansion lies in the conserved space | `â_{m,2}` (pipeline) | Fock (§2) |
|---|---|---|---|---|
| scalar | `T` | yes (2-dim at `s = 5/2`) | `â_{2,2} = 1`, `â_{−2,2} = c = 1`, `m = 0, ±1`: 0 | same |
| scalar | `φ²` (`Δ = 1`) | yes | `â_{0,2} = 1` | 1 |
| scalar | `J_4`, `J_6` | yes | `â_{±4,2} = 1`, `â_{±6,2} = −1`; all `\|m\| < s`: 0 | `(−1)^{s/2}` |
| fermion | `T` | yes | `â_{2,2} = 1`, `c = −1` | same |
| fermion | `ψ̄ψ` (`Δ = 2`, odd structure) | yes | `\|â_{0,2}\| = 1`, opposite sign for `h = −2` | modulus 1 |
| fermion | `J_4`, `J_6` | yes | `â_{s,2} = −â_{−s,2}`, modulus 1; `\|m\| < s`: 0 | same |

The scalar Ward identity `â^T_{2,2} = 1` is a genuine check of all normalisations (no free constant); for the
fermion the single constant is fixed by `C_T` and `â^T_{2,2} = 1` then holds. Explicit free data, e.g.
`⟨φ² T T⟩ = (4π)^{−3}(9/8 V_1²V_3² − 9/2 V_1V_3H_{13} + 3/4 H_{13}²)/(X_{13}^{9/2}X_{14}^{1/2}X_{34}^{1/2})` (`φ²` unnormalised).

What free fields do **not** test: non-conserved spinning `O` (`Δ > J + 1`, `J ≥ 2`) and parity-odd structures with
`J ≥ 2`. For these the evidence is §3.4 plus the fact that the closed forms of §5, continued to `Δ = J + 1`,
reproduce the free values (`free_canonical.py`, `J = 4, 6`).

---

## 5. Results

### 5.1 Selection rules and universal `m`-dependence  **[C] for `J ≤ 12`; [I] for `J > 12`**

For every conserved structure (`verify_closed_forms.py`; exact rational identities in `Δ`):

* **(Z)** `P_m = 0` for `m = 1, 2, 3`.
* **(R0)** `P_{−k}/P_0 = C(J,k)(4)_k/(Δ−1−k)_k` for `0 ≤ k ≤ J`, the same for all structures with `P_0 ≠ 0`.
* **(R4)** `P_{4+k}/P_4 = C(J−4,k)(4)_k/(Δ−5−k)_k` for `0 ≤ k ≤ J−4`, the same for all structures with `P_4 ≠ 0`.
* Sectors: `J = 0, 2`: only `m ≤ 0` (both parities). Even `J ≥ 4`: the two parity-even structures have both
  sectors; a unique combination `S^{+−}` has `P_{m≤0} = 0` and a unique combination `S^{++}` has `P_{m≥4} = 0`; the
  parity-odd structure has `P_{m≥4} = 0`. Odd `J ≥ 5`: the single structure has `P_{m≤0} = 0`.

Dividing by the norms of §3.3:

`|â_{−k,2}|²/|â_{0,2}|² = W⁰_k(Δ,J) = C(k+3,3)² (J+1−k)_{2k}/(Δ−1−k)_{2k} = C(k+3,3)² Π_{i=−k}^{k−1} (J+1+i)/(Δ−1+i)`,

`|â_{4+k,2}|²/|â_{4,2}|² = W⁴_k(Δ,J) = C(k+3,3)² (J−3−k)_k (J+5)_k / ((Δ−5−k)_k (Δ+3)_k)`.

`W_0 = 1`; `W⁰_1 = 16 J(J+1)/((Δ−2)(Δ−1))`, `W⁴_1 = 16 (J−4)(J+5)/((Δ−6)(Δ+3))`,
`W⁰_2 = 100 (J−1)J(J+1)(J+2)/((Δ−3)(Δ−2)(Δ−1)Δ)`, `W⁴_2 = 100 (J−5)(J−4)(J+5)(J+6)/((Δ−7)(Δ−6)(Δ+3)(Δ+4))`.
Both are decreasing in `Δ`, increasing in `J`, and positive above the unitarity bound. At `Δ = J + 1`:
`W⁰_k = C(k+3,3)²(J+k)/(J−k)`, `W⁴_k = C(k+3,3)²(J−4)(J+4+k)/((J−4−k)(J+4))`: the top components `k = J`, `k = J−4`
have poles, i.e. `a^{++}`, `a^{+−} → 0` with `|â_{∓J}|²` finite — a conserved current has only `m = ±J`.

Within each sector all components have the same phase; the parity-odd contribution to `a^{++}` is `±i` times
a real multiple of the parity-even one. With `h = −2`: `â_{m,−2}(λ_even, λ_odd) = â_{−m,2}(λ_even, −λ_odd)`.

Remarks. (i) `(J+1−k)_{2k}/k!²` and `(J−3−k)_k(J+5)_k/k!²` are the moduli of the Taylor coefficients of
`d̃^J_{0,0} = ₂F₁(−J,J+1;1;y)` and `d̃^J_{4,4} = ₂F₁(4−J,J+5;1;y)` (2201.06602 (A.3)): the two sectors are the CFT₃ images
of the helicity-conserving and helicity-flipping couplings `c^{++}_{J}`, `c^{+−}_{J}` of two gravitons to a spin-`J`
state, `y = p²/m² ↔ 1/Δ²`, and the weights `[(4)_k]²` come from the detector. (ii) No derivation of (R0), (R4) from
the conformal algebra is given here; Poincaré and dilatation invariance alone allow `2(2J+1)` independent
complex numbers `â_{m,±2}` (see §5.5).

### 5.2 Scalar exchange, `J = 0`  **[FF] + [S]**

`⟨TTO⟩ = λ_e S_e + λ_o S_o`, `S_e = (V_1²V_3² + 4(Δ−4)/(Δ+2) V_1V_3H_{13} + 2(Δ²−6Δ+6)/(Δ(Δ+2)) H_{13}²)/(X_{13}^{5−Δ/2}X_{14}^{Δ/2}X_{34}^{Δ/2})`
(coefficient of `V_1²V_3²` equal to one, as in KKSZ (5.110)), `S_o = ε_{13}X_{13}((Δ+1)V_1V_3 + (Δ−3)H_{13})/(X_{13}^{(11−Δ)/2}X_{14}^{(Δ+1)/2}X_{34}^{(Δ+1)/2})`
(CMT (D.3) with `(β_2−β_1)/4 = −1`). Then `â^O_{0,±2} = C_T^{−1/2}(λ_e K_e ± i λ_o K_o)` up to a common phase, and

`|â^O_{0,2}|² = (λ_e² K_e(Δ)² + λ_o² K_o(Δ)²)/C_T`,

`K_e² = (27/4) π^{9/2} Γ(Δ)Γ(Δ−1/2) / (Γ(2+Δ/2)⁴ Γ(3−Δ/2)² Γ((3+Δ)/2)²)`,

`K_o² = (27/4) π^{9/2} Γ(Δ)Γ(Δ−1/2) / (Γ(2+Δ/2)² Γ((1+Δ)/2)² Γ((3+Δ)/2)² Γ((7−Δ)/2)²)`.

Underlying amplitude: `M_0[S_e] = −3π^{7/2} i e^{iπΔ/2} 2^{−Δ} / (Γ(2+Δ/2)² Γ(3−Δ/2) Γ((3+Δ)/2))`.

* Sources: KKSZ (5.111) at d = 3 has the same modulus (after the factor 4 between `T(e_+)` and the `a⁴`
  component) and (5.132) gives `K_e² = 12 q(Δ)|_{d=3}`, exactly. CMT (4.10)–(4.11): `h^{even}_{3d}{}² = 2K_e²`,
  `h^{odd}_{3d}{}² = 2K_o²`, exactly, with `C^{even}_{TTO}/√C_O = λ_e`, `|C^{odd}_{TTO}|/√C_O = |λ_o|` (the factor 2 is the
  change from linear to helicity polarisations). No `1/(d−3)` pole occurs in the one-point amplitude; the poles
  of KKSZ (5.128), (5.130) belong to the SO(d−1) decomposition of the product, which in d = 3 is replaced by
  the Fourier decomposition of §6.
* Free fields: `K_e(1)² = (1024/27)π²` and `λ_{TTφ²}² = 81/(32768π⁴)` give `|â|² = 1`; `K_o(2)² = (32/3)π²` and the free
  fermion give 1 (§4).
* Zeros: `K_e` at `Δ = 6 + 2n`, `K_o` at `Δ = 7 + 2n` (simple zeros of the amplitude, double zeros of `K²`).
* Positivity: `K_e², K_o² ≥ 0`; both diverge like `(Δ − 1/2)^{−1}` at the scalar unitarity bound.
* Large `Δ`: `K_e² ≃ (27/4)·2^{12} π^{3/2} · 4^{Δ} Δ^{−29/2} sin²(πΔ/2)·(1 + O(1/Δ))` (Stirling; checked numerically at
  `Δ = 40.5, 400.5, 4000.5`: ratio 1.080, 1.0061, 1.0006; not used below).

**The relation `Â_0 = Â_4` (= `(S_0)`):**

`Σ_{scalars} (λ_e²K_e² + λ_o²K_o²)/C_T + Σ_{J ≥ 2 even} |a^{++}_O|² = Σ_{J ≥ 4} |a^{+−}_O|²`.

A scalar (or any `++`) coupling requires `+−` couplings of operators of spin ≥ 4.

### 5.3 Spin 2 and canonical structures for `J ≥ 4`

Normalisation conventions (coefficient vectors are polynomials in `s` without common factor; full tables in
`code/kernels/data/structures_J*.json`, `J ≤ 12`):

* `J = 2`: even structure with entry `V_1²V_3²V_4² = 2(s−5)(s−4)(7s−2)`; odd structure with entry `ε_{14}X_{14}V_1V_3²V_4 = 8s+4`.
* even `J ≥ 4`: `S^{+−}`: entry `H_{14}²H_{34}²V_4^{J−4}` has leading term `8s⁴`; `S^{++}`: entry `V_1V_3V_4^JH_{13}` has leading
  term `8s⁴`; `S^{odd}`: entry `ε_{14}X_{14}V_1V_3V_4^{J−2}H_{34} = J(J−1)(J−2)`.
* odd `J ≥ 5`: entry `ε_{13}X_{13}V_1V_3V_4^{J}` has leading term `8s³`.

Amplitudes per unit coefficient, source `T(e_+)`, `θ = 0` (**[C]**: `J = 2` symbolic; `J ≥ 4`: closed forms fitted
on `J = 4,…,8`, confirmed at `J = 9,…,12`; **[I]** beyond):

`M_0[J=2, even] = 3iπ³ e^{iπs} (Δ−3)(Δ−2) / (Γ(Δ+2) Γ(2+Δ/2) Γ(3−Δ/2))`,
`M_0[J=2, odd] = (3/4)π^{7/2} e^{iπs} 2^{1−Δ} (Δ−3)(Δ−2) / (Γ(2+Δ/2) Γ((Δ+3)/2)² Γ((7−Δ)/2))`,

`M_0[S^{++}] = 3i(−1)^{J/2} e^{iπs} π³ 2^{3−J} (Δ−J−1)_J / (Γ((Δ+J)/2) Γ(Δ+3) Γ(1−(Δ−J)/2))`,
`M_0[S^{odd}] = 3(−1)^{J/2} e^{iπs} π³ 2^{4−J} (Δ−J−1)_J / (Γ((Δ+J+1)/2) Γ(Δ+3) Γ((J+3−Δ)/2))`,
`M_4[S^{+−}] = −3i(−1)^{J/2} e^{iπs} π³ 2^{2−J} (Δ−J−1)_{J−4} / (Γ((Δ+J)/2) Γ(Δ−1) Γ(1−(Δ−J)/2))`,
`M_4[odd J] = −3i(−1)^{(J−5)/2} e^{iπs} π³ 2^{4−J} (Δ−J−1)_{J−4} / (Γ((Δ+J+1)/2) Γ(Δ−1) Γ((J+3−Δ)/2))`.

Squared kernels `K² = 12|M|²/N_{O,m}`:

`|K^{++}_0|² = 27π^{7/2} 2^{4+2Δ−J} (Δ−J−1)_J Γ(Δ+J)Γ(Δ−1/2) / (C(2J,J) Γ((Δ+J)/2)² Γ(Δ+3)² Γ(1−(Δ−J)/2)²)`,
`|K^{odd}_0|² = 27π^{7/2} 2^{6+2Δ−J} (Δ−J−1)_J Γ(Δ+J)Γ(Δ−1/2) / (C(2J,J) Γ((Δ+J+1)/2)² Γ(Δ+3)² Γ((J+3−Δ)/2)²)`,
`|K^{+−}_4|² = 27π^{7/2} 2^{2+2Δ−J} (Δ−J−1)_{J−4} Γ(Δ+J)Γ(Δ−1/2) / (C(2J,J+4)(Δ−1)_4 Γ((Δ+J)/2)² Γ(Δ−1)² Γ(1−(Δ−J)/2)²)`,
`|K^{oddJ}_4|² = 27π^{7/2} 2^{6+2Δ−J} (Δ−J−1)_{J−4} Γ(Δ+J)Γ(Δ−1/2) / (C(2J,J+4)(Δ−1)_4 Γ((Δ+J+1)/2)² Γ(Δ−1)² Γ((J+3−Δ)/2)²)`,
`|K_0[J=2, even]|² = (9/2)π^{7/2} 4^Δ (Δ−3)(Δ−2)Γ(Δ−1/2) / (Γ(Δ+2) Γ(2+Δ/2)² Γ(3−Δ/2)²)`,
`|K_0[J=2, odd]|² = (9/8)π^{9/2} (Δ−3)(Δ−2)Γ(Δ+2)Γ(Δ−1/2) / (Γ(2+Δ/2)² Γ((Δ+3)/2)⁴ Γ((7−Δ)/2)²)`.

Hence, for even `J ≥ 4`, `|a^{++}_O|² = (λ_{++}²|K^{++}_0|² + λ_{odd}²|K^{odd}_0|²)/C_T`, `|a^{+−}_O|² = λ_{+−}²|K^{+−}_4|²/C_T`; for `J = 2`
the same with the two `J = 2` kernels; for odd `J ≥ 5`, `|a^{+−}_O|² = λ²|K^{oddJ}_4|²/C_T`. All other components
follow from §5.1. All squared kernels are non-negative for `Δ ≥ J + 1`. The map
`(λ_{++}, λ_{odd}, λ_{+−}) ↦ (Re a^{++}, Im a^{++}, a^{+−})` (in the common phase convention) is diagonal and invertible away
from the zeros of §5.7, so `(a^{++}_O ∈ ℂ, a^{+−}_O ∈ ℝ)` **are** the normalised OPE data.

Free values (`free_canonical.py`) **[FF]**: for `J_4`: `λ_{+−}/√C_T = 0.0084940…`, `λ_{++}/√C_T = ∓0.0042470…`
(scalar/fermion); for `J_6`: `0.0038575…`, `∓0.0019288…`. The free scalar and the free fermion have the same `λ_{+−}`
and opposite `λ_{++}`, and `λ_{+−}²|K^{+−}_J|²/C_T = λ_{++}²|K^{++}_{−J}|²/C_T = 1` at `Δ = J + 1` (top components
`K_{±J}`, finite at the unitarity bound: `Mtop_pm`, `Mtop_pp` in `closed_forms.py`).

In a parity-invariant theory a parity-even operator of even spin has `(λ_{++}, λ_{+−})`, a parity-odd one only
`λ_{odd}` (it contributes only to `m ≤ 0`), and operators of odd spin only to `m ≥ 4`.

### 5.4 `O = T`  **[FF] for the parity-even part, [C] for the parity-odd part**  (`ttt.py`)

At `Δ = 3`, `J = 2` there are three structures (free boson `B`, free fermion `F`, one parity-odd). The parity-odd
one has `M_2 = 0` and `M_{−2}` purely imaginary relative to the even ones. With `⟨TTT⟩ = n_B B + n_F F + p S^{odd}_T`
(`B`, `F` the free-boson and free-fermion correlators, normalised to equal `C_T`; `C_T` below is the total):

`â^T_{2,2} = 1`,  `c = â^T_{−2,2} = (n_B − n_F)/(n_B + n_F) + i κ p/C_T`,  `κ = −2π²`

for `S^{odd}_T` with coefficients `(6,1,3,1,1;−4,−2)` on `(ε_{13}X_{13}·{V_1V_3V_4², V_1V_4H_{34}, V_3V_4H_{14}, H_{14}H_{34}, V_4²H_{13}};
ε_{14}X_{14}·{V_1V_3²V_4, V_3V_4H_{13}})`. So `t_4 = 4(n_B−n_F)/(n_B+n_F)` (as in 1707.03007 (1.6)) and `α_T = −8π²p/C_T` in
this normalisation; the dictionary between `p` and `p_T` of 1707.03007 was **not** worked out, and the identification
`Im 4c = α_T` rests on (E5) of `EINSTEIN_PROGRAM.md`.

### 5.5 Group theory and counting (items 4(i), 4(ii))

* Rest-frame building blocks for the bra polarisation `ζ̄ = ζ(b′,a′)`, source `e_+`, detector `z(θ)`:
  `ζ̄·z = (a′e^{−iθ/2} − b′e^{iθ/2})²/2`, `ζ̄·e_+ = a′²`, `ζ̄·p = −a′b′`, `z·e_+ = e^{iθ}`, `e_+·p = 0`. Lorentz invariance, scaling
  and homogeneity (degree `J` in `ζ̄`, 2 in `z_3`, `−2` in `z`) allow the parity-even monomials
  `(ζ̄·z)^α(ζ̄·e_+)^β(ζ̄·p)^γ(z·e_+)^{2−β}`, `α+β+γ = J`, `β ≤ 2` (`3J` of them for `J ≥ 2`), and as many parity-odd ones; the
  monomial `(α,β,γ)` has Fourier content `m ∈ [β−α, β+α]`, always with the phase `e^{i(2−m)θ}`. Since there are only
  `2J+1` charges, Poincaré + dilatation invariance gives **no** relation among the `â_{m,2}`: a priori `2J+1`
  independent complex numbers per helicity.
* Conformal invariance and conservation reduce them to: `J = 0, 2`: one complex number (2 structures:
  1⁺ + 1⁻); even `J ≥ 4`: one complex and one real number (3 structures: 2⁺ + 1⁻), carried by `(J+1) + (J−3) = 2J−2`
  non-vanishing components; odd `J ≥ 5`: one real number (1 structure), `J−3` components; `J = 1, 3`: nothing.
  The number of independent amplitudes equals the number of conserved structures in every case.

### 5.6 Large `Δ` and the Regge-like regime (item 4(iii))

All statements are exact consequences of §5.1 (**[C]**, `J ≤ 12`).

* Fixed `J`, `Δ → ∞`: `|K_{4+k}|²/|K_4|² = W⁴_k → C(k+3,3)² (J−3−k)_k(J+5)_k/Δ^{2k}`, `|K_{−k}|²/|K_0|² = W⁰_k →
  C(k+3,3)² (J+1−k)_{2k}/Δ^{2k}`. In particular `|K_6|²` is suppressed by `Δ^{−4}` relative to `|K_4|²`, and
  `|K_{−2}|²` by `Δ^{−4}` relative to `|K_0|²`, with coefficients `100(J−5)(J−4)(J+5)(J+6)` and `100(J−1)J(J+1)(J+2)`.
* `|K_6|²` against `|K_{−2}|²`, `|K_0|²`: no kernel statement exists. `K_6 ≠ 0` only for `S^{+−}` (and odd `J`);
  `K_{−2}`, `K_0 ≠ 0` only for `S^{++}`, `S^{odd}`. The ratio is `(λ_{+−}/λ_{++})²` times
  `|K^{+−}_4|²W⁴_2/(|K^{++}_0|²W⁰_2)`, and in the normalisations of §5.3
  `|K^{+−}_4|²/|K^{++}_0|² = C(2J,J)(Δ−1)_4/(4 C(2J,J+4)(Δ−5)_4) → C(2J,J)/(4C(2J,J+4))`: a constant. No parametric
  hierarchy between the two sectors is implied by the kernels; it would have to come from the OPE data.
* `J = xΔ`, `Δ → ∞`, `0 < x < 1` (fixed `k`): `W⁰_k, W⁴_k → C(k+3,3)² x^{2k}`. No suppression: `W⁴_2 → 100x⁴`, which exceeds
  1 for `x > 0.316`. Near the unitarity bound, `Δ = J + 1 + τ`, `J → ∞`: `W⁴_2 → 100`, while the top component
  dominates, `W⁴_{J−4} ∼ J^{15}/τ`.
* Absolute size: `Γ(1−(Δ−J)/2)^{−2} = sin²(π(Δ−J)/2) Γ((Δ−J)/2)²/π²`; all `K²` grow exponentially in `Δ` at fixed `J`
  (to be compensated by the decay of `λ²`); not analysed further.

### 5.7 Zeros and double-trace operators (item 4(iv))

From §5.2–5.3, as functions of `Δ` at fixed finite couplings in the polynomial bases:

| structure | zeros of the amplitude (simple) | `N = ∞` operators at these dimensions |
|---|---|---|
| `J = 0` even / odd | `Δ = 6+2n` / `Δ = 7+2n` | `T_{μν}∂^{2n}T^{μν}` / `ε^{αβγ}T_{αδ}∂^{2n}∂_βT_γ{}^δ` (CMT §4) |
| `J = 2` even / odd | `Δ = 6+2n` / `Δ = 7+2n` | `[TT]`, twist 4+2n / twist 5+2n |
| even `J ≥ 4`: `S^{++}` and `S^{+−}` | `Δ = J+2+2n` | `[TT]_{n,J}`, twist 2+2n |
| even `J ≥ 4`: `S^{odd}`; odd `J ≥ 5` | `Δ = J+3+2n` | `ε`-type `[TT]`, twist 3+2n |

All components `m` of a given structure share the zero (the `m`-dependence of §5.1 is rational and positive).
Hence `|â^{[TT]}_{m,2}|² = (λ²/C_T)|K_m|² ≃ (λ²/C_T)|∂_ΔK_m|² γ²`. With `λ²_{TT[TT]}/C_T = O(C_T)` and `γ = γ̂/C_T`:

**each double-trace `[TT]_{n,J}` contributes `O(γ̂²/c_T)` to both sides of every `(S_k)`** (first order in `1/c_T`, not
second: KKSZ's "`O(1/N⁴)`" refers to `γ²`, and `λ²/C_T ∼ N²`). The other double-traces (`[TO]`, `[OO′]` with
single-trace `O`) do not sit at zeros but have `λ²/C_T = O(1/c_T)`. Uniformity of these estimates over the infinite
towers in `(n, J)` is **not** shown here. The identification of the right-hand column is by quantum numbers only.

---

## 6. The sum rules through OPE data

Under hypotheses (E), with `a^{++}_O`, `a^{+−}_O` as in §5 and `c = (t_4 + iα_T)/4`:

**(S_k)**, `k = 0, 1, 2, …`:

`a_{2+k} = δ_{k,2}|c|² + Σ_{O ≠ T; J_O even, J_O ≥ k} |a^{++}_O|² W⁰_k(Δ_O,J_O) = Σ_{O; J_O ≥ 4+k} |a^{+−}_O|² W⁴_k(Δ_O,J_O)`.

* Left side: `Â_{−k}`; right side: `Â_{4+k}`; `T` contributes `1` to `Â_2` and `|c|²` to `Â_{−2}`; `Â_1 = Â_3 = 0`.
* `(S_2)` is (★):

  `(t_4²+α_T²)/16 + Σ_{O≠T, J even ≥ 2} |a^{++}_O|² · 100(J−1)J(J+1)(J+2)/((Δ−3)(Δ−2)(Δ−1)Δ) = Σ_{J ≥ 6} |a^{+−}_O|² · 100(J−5)(J−4)(J+5)(J+6)/((Δ−7)(Δ−6)(Δ+3)(Δ+4))`,

  with `|a^{++}_O|² = (λ_{++}²|K^{++}_0|² + λ_{odd}²|K^{odd}_0|²)/C_T`, `|a^{+−}_O|² = λ_{+−}²|K^{+−}_4|²/C_T` (§5.3).
  Refinement of `EINSTEIN_PROGRAM.md` §3: the left side receives only `++`/odd couplings (even spin ≥ 2, any
  parity); the right side only `+−` couplings (spin ≥ 6, including odd spins ≥ 7); the `+−` coupling of a
  spin-4 operator does not enter (★) at all.
* `𝓔(θ) = 1 + 2Σ_{k≥0} a_{2+k} cos((2+k)θ) ≥ 0` (positive measure), so `0 ≤ |a_n| ≤ 1`.
* The off-diagonal relations (O), `B̂_m = B̂_{−m}` with `B̂_m = Σ_O conj(â_{m,−2}) â_{m,2}`, hold **identically**: by
  `â_{m,−2}(λ_e,λ_o) = â_{−m,2}(λ_e,−λ_o)` and the reality properties of §5.1, each operator contributes equally to
  `B̂_m` and `B̂_{−m}`. Commutativity in the two-dimensional helicity space is exactly the family `(S_k)`.
* Free theories: `(S_k)` reads `1 = 1` for even `k` (`J_k` on the left for `k ≥ 4`, `φ²`/`ψ̄ψ` for `k = 0`, `|c|² = 1` for
  `k = 2`; `J_{k+4}` on the right) and `0 = 0` for odd `k`.

---

## 7. Consequences

### 7.1 (a), (b): the planar statement

In a large-`c_T` theory with factorisation, by §5.7: `(S_k)` hold with the sums restricted to **single-trace**
operators, up to `O(1/c_T)`:

`δ_{k,2}|c|² + Σ_{ST, J even} |a^{++}_O|² W⁰_k = Σ_{ST, J ≥ 4+k} |a^{+−}_O|² W⁴_k + O(1/c_T)`.

Hypotheses: (E) (commutativity in the planar limit requires Regge intercept `< 3`, KKSZ (4.79)–(4.80));
uniform convergence of the double-trace towers (not shown); the general-`J` formulas of §5 (**[C]** `J ≤ 12`).
The `O(1/c_T)` term is `Σ_{[TT]}(λ²/C_T)|∂_ΔK|²γ² +` contributions of `[TO]`, `[OO′]`; its coefficient was not
computed.

### 7.2 (d): infinite gap

If no single-trace operator of spin ≥ 6 appears in `T × T`, then at leading order `(S_2)` has vanishing right
side and non-negative terms on the left: `t_4 = α_T = 0` **and** `a^{++}_O = 0` for every single-trace `O` of spin
≥ 2 (even or odd parity, light or heavy). Then `(S_1)` gives `a^{+−}_O = 0` for spin 5, and `(S_0)` reduces to
`Σ_{ST scalars}(λ_e²K_e² + λ_o²K_o²)/C_T = Σ_{ST, J=4} λ_{+−}²|K^{+−}_4|²/C_T`: scalar couplings `⟨TTφ⟩` are allowed only if
balanced by spin-4 single-trace operators (if these are absent too, `λ_{TTφ} = 0` unless `Δ_φ` sits at a zero).
Corrections are `O(1/c_T)`. No assumption on light scalars is needed for `t_4 = α_T = 0`. Exact hypotheses: as in §7.1.

### 7.3 (c): is there a gap-dependent bound?

Setting: all operators of spin > 2 (single-trace, planar limit) have `Δ ≥ Δ_gap`; scalars and non-conserved
spin-2 operators are unrestricted. Inputs: `(S_k)` and `𝓔 ≥ 0`. Three statements of different strength.

**(i) Conditional bound (two lines, rigorous given §5.1).** From `(S_2)`, `(S_0)` and `a_2 ≤ 1`:

`|c|² ≤ Σ_{J≥6}|a^{+−}_O|²W⁴_2 ≤ (sup_{O: J_O≥6} W⁴_2(Δ_O,J_O)) · Σ_{J≥4}|a^{+−}_O|² = a_2 · sup W⁴_2 ≤ sup_{O: J_O ≥ 6} W⁴_2(Δ_O, J_O)`.

If the spin-≥ 6 operators have `J ≤ J_max`: `t_4² + α_T² ≤ 1600 (J_max−5)(J_max−4)(J_max+5)(J_max+6)/((Δ_gap−7)(Δ_gap−6)(Δ_gap+3)(Δ_gap+4))`;
if they satisfy `J ≤ xΔ`: `≤ 1600 x⁴(1 + o(1))`. With a dimension gap only, `sup W⁴_2 ≥ 100` (large spin near the
unitarity bound), and this argument gives nothing: the kernels have **no hierarchy uniform in `J`**. Question
(d) of ST1 is answered in the negative: for an operator with only a `+−` coupling and `J = xΔ → ∞`,
`|â_6|²/Σ_m w_m|â_m|² → 100x⁴/Σ_k w_{4+k}C(k+3,3)²x^{2k}`, a positive constant for any fixed weights for which the sum
converges.

**(ii) Explicit exact solutions with `c ≠ 0` at arbitrary gap (rigorous given §5.1; `family4.py`).** Take
`+−` operators `O_4` (`J = 4`, weight `q_4`) and `O_6` (`J = 6`, weight `q_6`), both at `Δ = Δ_gap`, and `++` operators: a scalar
(weight `p_0`) and a spin-2 operator at `Δ′` (weight `p_2`). All `(S_k)` hold iff

`p_0 + p_2 = q_4 + q_6`,  `p_2 · 96/((Δ′−2)(Δ′−1)) = q_6 · 352/((Δ−6)(Δ+3))`,
`|c|² = q_6 · 26400/((Δ−7)(Δ−6)(Δ+3)(Δ+4)) − p_2 · 2400/((Δ′−3)(Δ′−2)(Δ′−1)Δ′)`,

(`k ≥ 3`: `0 = 0`), and `𝓔 = 1 + 2a_2cos2θ + 2a_3cos3θ + 2a_4cos4θ` with `a_2 = q_4+q_6`, `a_3 = q_6W⁴_1`, `a_4 = q_6W⁴_2`.
Example (`Δ_gap = 20`, `q_4 = 0.29`, `q_6 = 0.16`, `Δ′ = 17`): `p_0 = 0.0127`, `p_2 = 0.4373`, `|c|² = 0.0237`, `min 𝓔 = 0.013 > 0`.
Optimising: `|c|²_max → 5400 a_2/Δ_gap⁴` with `Δ′ = (2/3)^{1/2}Δ_gap`, `a_2 → 1/2`, i.e. `|c|² ≃ 2700/Δ_gap⁴`
(numerical optimum, positivity enforced to `10^{−4}`: `|c|²Δ_gap⁴ = 3995, 3925, 3493, 3119, 2885, 2716` at
`Δ_gap = 15, 20, 30, 50, 100, 1000`). By the cone
property (rescaling `c` and all heavy amplitudes by `t ≤ 1` preserves `(S_k)` and positivity) every smaller `|c|`
is realised too. **Hence `(S_k)` and positivity of `𝓔` can never imply a bound stronger than
`t_4² + α_T² ≲ 4·10⁴/Δ_gap⁴`**; the rate `(log Δ_gap)/Δ_gap⁸` of `EINSTEIN_PROGRAM.md` (P1) is out of reach of route S-B by
itself, whatever the treatment of the Regge regime. Conversely it is **not** true that any `|c| ≤ 1` can be
balanced at arbitrarily large gap by these simple families.

**(iii) Is there an upper bound of this order? Numerical evidence only (`lp_balance.py`).** Linear
programme: maximise `|c|²` over non-negative weights on a grid of `(Δ, J)` (`++`: even `J ≤ J_max`; `+−`: `4 ≤ J ≤ J_max+4`;
`Δ ≥ max(Δ_gap, J+1.001)` for `J > 2`, 23 to 74 values of `Δ` per spin up to `Δ + 4000`), subject to `(S_k)`, `k ≤ J_max`
(exact for a spin-truncated ansatz), and `𝓔(θ) ≥ 0` on 721–1441 points.

| `Δ_gap` | 6 | 7.5 | 10 | 15 | 20 | 30 |
|---|---|---|---|---|---|---|
| `J_max` | 24 | 24 | 80 | 80 | 100 | 100 |
| max `\|c\|²` | 0.924 | 0.854 | 0.475 | 0.109 | 0.0308 | 0.00536 |
| `\|c\|² Δ_gap⁴` | — | — | 4750 | 5540 | 4930 | 4340 |

Stable under `J_max` (10: 0.467/0.473/0.475 at 36/56/80; 15: 0.1084/0.1093 at 40/80; 20: 0.0307/0.0308 at 60/100)
and under refinement of the `Δ` grid (0.469, 0.110, 0.0320 at `Δ_gap` = 10, 15, 20 with 74 values per spin). The
optimum is dominated by the four operators of (ii) plus small higher-spin admixtures that raise `a_2`. A grid
LP gives feasible points, i.e. lower bounds on the true supremum; an upper bound requires a dual functional
valid for all `(Δ, J)`, which was **not** constructed. Heuristic reason to expect one **[I]**: for heavy operators
with `J = xΔ` both sectors have the same profile `C(k+3,3)²x^{2k}`, and a signed measure on `x² ∈ [0,1]` whose moments
vanish except the second does not exist, so the Regge regime cannot produce the `δ_{k,2}` mismatch at leading
order; the mismatch comes from `1/Δ` effects and from low spins, as in (ii).

**Answer to 5(c).** The expectation "no gap-dependent bound because the relations are homogeneous" is not
confirmed in that form: homogeneity is broken by `a_n ≤ a_0 = 1` (positivity of the energy correlator
normalised by the Ward identity). What is established: (★) with `a_2 ≤ 1` alone gives a bound only under an
extra hypothesis on `J/Δ` (i); the full system admits `|c|² ≃ 2700/Δ_gap⁴` (ii), and numerically not more than
`≈ 5·10³/Δ_gap⁴` (iii, unproven). In all cases the result is weaker by `Δ_gap⁴` than the expected `Δ_gap^{−8}`; the missing
factor has to come from the spin-two sum rule with the graviton pole (route S-A), in agreement with the
conclusion of `EINSTEIN_OBSTRUCTIONS.md` §3.4 but for a different reason.

Caveats specific to §7.3: only the state `|T,2⟩` is used for positivity (superpositions of the two helicities,
and of `T` with other operators as in CMT, give further conditions; by the cone property they cannot remove
the family (ii) for small enough rescaling, but they can lower the constants); the solutions of (ii), (iii) are
solutions of necessary conditions, not CFTs.

---

## 8. What was not computed

1. A proof of (R0), (R4), (Z) and of the closed forms of §5.3 for all `J` (verified symbolically for `J ≤ 12`).
   The norm formula §3.3 is proven for `J ≤ 5` only.
2. A free-field or other independent check of the parity-odd structures with `J ≥ 2` and of non-conserved
   spinning exchanges (no free example exists; large-`N` Chern–Simons-matter data were not used).
3. Hypotheses (E): commutativity for momentum eigenstates, absence of contact terms at coincident angles
   (the free theories do have a symmetric `δ(θ_{12})` term), convergence of the mode sums. Not examined.
4. The `O(1/c_T)` terms of §7.1 (anomalous dimensions of `[TT]`, `λ_{TT[TO]}`), and uniformity over the towers.
5. A dual functional proving an upper bound in §7.3 (iii); positivity in superposition states.
6. The dictionary between the canonical structures of §5.3 and other bases in the literature (1708.05718
   App. A, the helicity basis of 2102.08160), and between `p` of §5.4 and `p_T`, `α_T` of 1707.03007. Only the
   scalar structures were matched (KKSZ, CMT).
7. Large-`Δ` asymptotics of the absolute kernels beyond the leading Stirling form quoted in §5.2; OPE-coefficient
   asymptotics.
8. Sources checked for d = 3 results: 2512.09986 (blocks for a *scalar* source, parity-odd included; remarks
   that the commutator vanishes trivially block by block there) and 2608.19322 (positivity bounds, scalar
   source; stress-tensor source listed as future work) — as recorded in `EINSTEIN_OBSTRUCTIONS.md` §3.2; neither
   contains the kernels above. Gillioz et al. (momentum-space Wightman functions) were not used.

---

## 9. Code map (`ads4/code/kernels/`)

| file | content |
|---|---|
| `fock_free.py` | §2: Fock-space wavefunctions, `c`, `â`, energy correlator |
| `structures.py` | §3.1: monomials, independence, conservation and Bose symmetry (sympy, exact) |
| `pipeline.py` | §3.2–3.3: detector limit, `α`-integral, Fourier transform, norms (python-flint, exact in `s`) |
| `wick_free.py`, `free_checks.py` | §4: free correlators in position space, projection, amplitudes |
| `ttt.py` | §5.4 |
| `compute_kernels.py J` | conserved structures and `P_m(s)` for spin `J` (cached) |
| `verify_closed_forms.py` | (R0), (R4), (Z) for `J ≤ 12` |
| `canonical.py`, `closed_forms.py`, `data/structures_J*.json` | §5.2–5.3: canonical structures, closed forms, comparison with KKSZ/CMT |
| `free_canonical.py` | free values of `λ_{++}`, `λ_{+−}`; closed forms at `Δ = J+1` |
| `family4.py`, `lp_balance.py gap Jmax` | §7.3 (ii), (iii) |
| `run_all.sh` | everything, about 10 minutes |
