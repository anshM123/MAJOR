# PARITY_ODD_DZ: the parity-odd, odd-twist `[TT]` double zeros (row L1)

*25 Nov 2026. Research agent, about 4 h. Unreviewed and not committed. No existing file was edited; `prl/`, `code/epswcert/`
and the Supplement were not touched.*

* **Code.** New directory `ads4/code/oddz/` (map in §8). Outputs go to `oddz/cache/`, which ignores everything but its own
  `.gitignore`. `code/kernels/{pipeline,kernels_lib,structures}.py` are imported read-only (`KERNELS_CACHE` redirected to
  `oddz/cache`). `code/angproofs/` is used read-only: its symbolic-J structure vectors (`cache/symbolic_structures.pkl`) and a
  re-run of `ap5_kernels.py` whose output is captured in `oddz/cache/oz1_out.txt`.
* **Compute.** About 1 CPU-h, at most 2 processes, `ulimit -v 3000000`, `OMP_NUM_THREADS=1`.

**Grades.**

| grade | meaning |
|---|---|
| [P] | proved: hand proof, or exact computer algebra with symbolic `Δ` (and symbolic `J` where stated) |
| [P \| X] | proved given the named input X |
| [DV] | derived, and checked exactly at finitely many cases |
| [D] | derived, not machine-checked |
| [Q] / [S] | quoted from the literature / standard |

---

## 0. Verdict

**Proved, at the level L1 needs. The parity-odd, odd-twist part of L1 moves from [D-sketch] to [D], the same grade as its
parity-even part. With it, the cubic implication (1), including the `α_T` bound, is [D] for parity-violating theories too.**

The statement that was sketched is a statement about kernels: the action of the `k = 2, 3` functionals on a parity-odd `⟨TT𝒪⟩`
exchange has a double zero at every parity-odd `[TT]` dimension. It is now a theorem about explicit closed forms. It does not use
CLS-D's double-commutator factor, which the caveat said had not been derived for parity-odd structures.

| statement | grade |
|---|---|
| (Z0) **Which dimensions.** At `c_T = ∞` every parity-odd `[TT]` primary has odd twist, `Δ = J+3+2n` (`J ≥ 4`, and odd `J ≥ 5`) or `Δ = 7+2n` (`J = 0, 2`); every parity-even one has even twist | **[P]** (index counting, §1.4) |
| (Z1) **Extreme components, all `J`, all `n`.** `M_{−J}[S^odd] ∝ cos(πh)Γ(h−½)/(Γ(s+½)Γ(Δ+3))` and `M_{+J}[S^oddJ] ∝ cos(πh)Γ(h−½)/(Γ(s+½)Γ(Δ−1))`, with `h = (Δ−J)/2`. For `h > ½` they vanish **exactly** at `h = 3/2 + n`, i.e. at `Δ = J+3+2n`, and the zeros are **simple**; the derivative is given in closed form | **[P]** (symbolic `J`: ANGULAR_PROOFS §3 closed forms, re-run here, plus the reflection formula; §2) |
| (Z2) **All components `m`, both source helicities `h = ±2`.** Every component has a simple zero at the same points | **[P]** for every `J ≤ 31` (and `J = 40`, `S^odd`) and every `n` (exact symbolic-`Δ` computation per `J`, §3.2, §3.4); **[P \| ladder (L)]** for all `J` (§3.1); (L) is [P given light-transform conservation, D], with an exact pipeline check for `J ≤ 40` (LADDER_ALLJ) |
| (Z3) **First moments `E⁽¹⁾`** (the `k = 2` weight is affine in the light-ray times). Every component has a zero of order ≥ 1 at the same points | **[P]** for every `J ≤ 30` and every `n` (§3.2); **[P \| (L), [z·K, E] = 6iE⁽¹⁾]** for all `J` (§3.1) |
| (Z4) **`J = 0, 2`.** Simple zeros exactly at `Δ = 7+2n`, all components, both helicities, both moments; **no** zero at `Δ = 5` | **[P]** (exact, symbolic `Δ`; §4) |
| (Z5) **Basis regularity.** The canonical odd structures are nonzero at every `Δ₀`, so the GFF OPE coefficients in this basis are finite and the zero is not a basis artefact | **[P]** symbolic `J` for `S^odd`; per `J ≤ 31` for `S^oddJ` (§3.3) |
| (DZ-odd) **The L1 statement.** The contribution of each parity-odd `[TT]_{n,J}` to the `k = 2, 3` light side is `∝ γ²`, i.e. `O(c_T^{−2})` in the functional normalisation, against `O(c_T^{−1})` for the `T` block. Parity-violating admixtures (even couplings of odd-twist `[TT]`, odd couplings of even-twist `[TT]`) are also `O(c_T^{−2})` | **[D]** (§1.2–1.3, §5.2): the zero orders are [P]; the factorisation of the functional through light-ray one-point matrices is C1_FINAL §2.2 [D], and the large-`c_T` counting is standard [D] |
| (LS) **Light-side consistency.** The parity-odd part of the `T` block is not at a zero: `T` has `h = ½`, where `1/Γ(3/2−h) = 1`. Its odd `⟨TTT⟩` amplitude is pure flip, and its first moment has no diagonal element. So `α_T` enters the `k = 2` light side only as `α_T²/16` (flip × flip), with no linear term | **[P]** (exact pipeline at `Δ = 3`, §5.1). This also settles the "`2cos²(πτ/2)` vanishes at `τ = 1`" aside of PARITY_ODD_ROUTEA §12.5 |

**What is unchanged.** The double zeros of the parity-even `[TT]` were [D given Q] through CLS-D's factor. The same method also
re-derives them: the extreme amplitudes carry `1/Γ(1−h) = sin(πh)Γ(h)/π`, with zeros exactly at the even twists (§2(d), §3.2
controls). The remaining links of L1 are the same in both parities and are inherited, not new:
* the factorisation of the functional through light-ray one-point matrices (C1_FINAL §2.2, [D]);
* the uniformity of the light-side remainder over the infinite `[TT]` towers (inside A2, Supplement "Order of limits").

**Relation to PARITY_ODD_ROUTEA §12.5.** That section already argued for the odd-sector factor `2cos²(πτ/2)`, graded
[P given Q] with GFF `τ`-parity [P-sketch]. It rests on CLS-D's multivaluedness bookkeeping applied to odd structures, which
INDEP_CHECK3 did not accept as shown. The route here is independent of it and agrees with it, including the factor:
`|K^odd|² ∝ cos²(πτ/2)`.

---

## 1. What L1 needs

### 1.1 Where the caveat comes from

H0 §15.1, item 1, argues that the double traces drop out of the light side:
* CLS-D's double commutator gives `2 sin²(π(τ−10)/2)` for **parity-even** `⟨TT𝒪⟩` [Q], with a double zero at every even `τ`.
* The polynomial weights `f₂ = v′ − u′` and `f₃ = −1` keep those zeros.
* Hence the even-twist `[TT]` enter `∝ γ²`.

The caveat is the next bullet. In `d = 3`, `[TT]` also contains parity-odd double traces (`ε`-contractions) with **odd** `τ`.
`sin²(π(τ−10)/2) = 1` there. Their double zeros would need the parity-odd structures to shift the phase by `π/2`, giving
`sin²(π(τ±1−10)/2) = cos²(πτ/2)`. CLS-D derive the factor only for parity-even structures, so this was [D-sketch].

INDEP_CHECK3 §1.5 points out that the link is load-bearing. At `O(ĝ₃)` the parity-odd double traces do receive `O(1/c_T)` data.
The same input carries the Einstein term `16πG/(ν²+9/4)` and `c₁ = 0` (L3). So a failure would contaminate the light side at the
order of the `T` block.

### 1.2 Reduction to light-ray one-point amplitudes [D, inherited from C1_FINAL §2.1–2.2]

Put the light-ray base point at spatial infinity, `x₀ = ∞`, as in C1_FINAL. Then:
* `L[T](∞,(1,n)) ∝ E(n)`;
* the `k = 2` weight is affine in the retarded times [P], so the light-ray product is `aE₃E₁ + bE₃E₁⁽¹⁾ + cE₃⁽¹⁾E₁`, with `a, b, c`
  depending only on the external labels `(n₁, n₃, x₂, x₄)`;
* the `k = 3` weight is `f₃ = −1`, so only `E₃E₁` appears.

Since `L[T]|Ω⟩ = 0`, the double commutator equals the Wightman function `⟨T₄ L₃ L₁ T₂⟩`. Its `s`-channel contribution from one
primary `𝒪` is the insertion of the projector onto the `𝒪` family between `L₃` and `L₁`. C1_FINAL §2.2 writes this for `𝒪 = T`;
the argument does not use any property of `T` and holds for every `𝒪`:

    Ω_k[𝒪] = Σ ⟨Ψ₄| A₃ P_𝒪 A₁ |Ψ₂⟩,   A ∈ {E, E⁽¹⁾} (k = 2),  A = E (k = 3),
    P_𝒪 = Σ_m ∫d³p |𝒪(p),m⟩⟨𝒪(p),m| / N_{𝒪,m}(p),

with `Ψ₂ = T(x₂)|Ω⟩`, which lies in the `T` family. The only `Δ_𝒪`-dependent objects are therefore:
* the one-point amplitudes `⟨𝒪(p),m|E(n)|T(p),h⟩ = (2π)³δ³ M_{m,h}(p,n)`;
* the first-moment amplitudes `⟨𝒪|E⁽¹⁾(n)|T⟩`. Since `[P^μ, E⁽¹⁾] = −i n̄^μ E`, these are a `δ′` term `∝ M` plus `δ³·M⁽¹⁾`;
* the norms `N_{𝒪,m}`, which are finite and positive for `Δ > J + 1`.

By Lorentz covariance, everything at timelike `p` and any `n` is a finite linear combination of the rest-frame components
`M_{m,±2}(θ = 0)` and `M⁽¹⁾_{m,±2}(θ = 0)`, `|m| ≤ J`; the dependence on `θ` is the phase `e^{i(h−m)θ}`.

**Consequence.** `Ω_k[𝒪](Δ)` is a finite sum of products `M_A(Δ)·conj(M_B(Δ))` of these components. If every component vanishes
at `Δ₀` to order ≥ 1, then `Ω_k[𝒪]` and `∂_ΔΩ_k[𝒪]` both vanish at `Δ₀`: a double zero. (`E` and `E⁽¹⁾` are Hermitian, so the
`3`-side factor is the complex conjugate of a one-point amplitude of the same kind.)

### 1.3 The precise requirement, and the `1/c_T` counting

**(DZ-odd).** Take a parity-odd conserved structure (`S^odd` for even `J ≥ 4`, `S^oddJ` for odd `J ≥ 5`, the odd structure at
`J = 0, 2`), in a basis that is regular and nonvanishing at `Δ₀` (§3.3). Then every rest-frame component
`M_{m,h}(Δ)` (`|m| ≤ J`, `h = ±2`) and `M⁽¹⁾_{m,h}(Δ)` vanishes to order ≥ 1 at every `Δ₀ ∈ 𝒟_odd(J)`, where

    𝒟_odd(J) = {J+3+2n : n ≥ 0}  (J ≥ 4 even, and odd J ≥ 5),     𝒟_odd(J) = {7+2n : n ≥ 0}  (J = 0, 2).

These sets contain all parity-odd `[TT]` dimensions at `c_T = ∞` (Z0).

**Counting.** There are two equivalent normalisations.

* *Functional normalisation* (unit-normalised `T`, disconnected part `O(1)`).
  * The `T` block enters at `O(1/c_T)`.
  * A `[TT]_{n,J}` has `p = p⁽⁰⁾ + p⁽¹⁾/c_T` and `Δ = Δ₀ + γ̂/c_T`.
  * The `O(1/c_T)` double-trace terms of the correlator are `½p⁽⁰⁾γ̂ ∂_ΔG|_{Δ₀} + p⁽¹⁾G|_{Δ₀}`.
  * Both are annihilated when `Ω_k[𝒪](Δ₀) = ∂_ΔΩ_k[𝒪](Δ₀) = 0`, which is the double zero. What remains is `½p⁽⁰⁾γ̂²Ω″(Δ₀)/c_T²`.
* *State normalisation* (KERNELS_D3 §5.7; canonical `T`, unit `𝒪`).
  * `â = C_T^{−1/2}λK` with `λ_{TT[TT]} = O(C_T)`, so `λ²/C_T = O(C_T)`.
  * With a simple zero, `K(Δ₀ + γ) = ∂_ΔK·γ + O(γ²)` and `γ = O(1/C_T)`.
  * So `|â|² = (λ²/C_T)|∂_ΔK|²γ² = O(1/C_T)`. The `T` contributes `O(1)` (exactly `1` to `Â₂`).
  * Without a zero, `|â|² = O(C_T)` would exceed the bound `|a_n| ≤ 1` of KERNELS_D3 §6. So the zero is forced as well as
    sufficient (see the GFF remark, §5.4).

**The zero order required.** A simple zero of the amplitude suffices: `|K|²` then has a double zero, and the double-trace
contribution is relative order `1/c_T` to the `T` block. A zero of higher order would only help. Note that the double zero must
sit in the **amplitude-squared**: a simple zero of `|K|²` would leave `p⁽⁰⁾γ̂∂_Δ` at `O(1/c_T)`.

### 1.4 (Z0) Twist parity of GFF double traces [P]

At `c_T = ∞`, a `[TT]` primary is a linear combination of `∂^{a}T ∂^{b}T` with all indices contracted by `η` and `ε`. In `d = 3`
these are the only invariant tensors, and `εε` reduces to `η`'s, so parity-odd means exactly one `ε`. With `ℓ` derivatives,
`Δ = 6 + ℓ` and there are `4 + ℓ` indices.

The `ε` cannot carry two free indices, because the free indices are symmetrised. Contracting two `ε` indices on the same `T`, or
on two derivatives of the same `T`, gives zero. With `t` traces:

* **(A) `ε` carries one free index.** `4+ℓ = 2 + (J−1) + 2t`, so `τ = Δ − J = 3 + 2t`.
* **(B) `ε` is fully contracted.** `4+ℓ = 3 + J + 2t`, so `τ = 5 + 2t`.
* **Parity-even (no `ε`).** `4+ℓ = J + 2t`, so `τ = 2 + 2t`.

Since `ℓ ≥ 0`, (A) needs `J ≥ 3 − 2t`. Hence:
* `J = 0`: only (B), with `ℓ = 2t − 1 ≥ 1`, so `Δ ∈ {7, 9, …}`;
* `J = 2`: (A) needs `t ≥ 1` and (B) gives `Δ ≥ 7`, so `Δ ∈ {7, 9, …}`;
* `J ≥ 4`: `Δ ∈ {J+3, J+5, …}`.

The same holds for odd `J`. (Spins `J = 1, 3` have no conserved `⟨TT𝒪⟩` structure, so they are irrelevant.) So the parity-odd
`[TT]` dimensions lie in `𝒟_odd(J)`, and the parity-even ones in `𝒟_even(J) = {J+2+2n}` (`J ≥ 4`) and `{6+2n}` (`J = 0, 2`).
This is an inclusion only. It does not use, and the argument does not need, the existence or multiplicity of primaries at each
point. Degenerate primaries at one `Δ₀` each sit at the zero separately.

---

## 2. Extreme components for symbolic `J` [P] (`oz1_extreme.py` → `cache/oz1_out.txt`)

**(a) Closed forms.** ANGULAR_PROOFS §3.2 computes the extreme detector amplitudes for symbolic `J`. It expands `V₄^{J−k}`
binomially and sums each term by Chu–Vandermonde. Re-run here read-only:

    odd   minusJ: (symbolic-J amplitude)/(claimed closed form) = 1
    oddJ  plusJ:  (symbolic-J amplitude)/(claimed closed form) = 1
    odd   plusJ:  absent sector = 0;   oddJ minusJ: absent sector = 0

So, for all even `J ≥ 4` and all odd `J ≥ 5` respectively (source `T(e₊)`, `θ = 0`, `h = (Δ−J)/2`, `s = (Δ+J)/2`):

    M_{−J}[S^odd]  = 8(−1)^{J/2} e^{iπs} π³ 2^{−J} (J+3)! / (Γ(s+½) Γ(Δ+3) Γ(3/2−h)),
    M_{+J}[S^oddJ] = −8i(−1)^{(J−5)/2} e^{iπs} π³ 2^{−J} (J−1)! / (Γ(s+½) Γ(Δ−1) Γ(3/2−h)).

**(b) Zeros.** By reflection, `1/Γ(3/2−h) = −cos(πh)Γ(h−½)/π` (checked symbolically). Hence

    M_{−J}[S^odd]  = −8(−1)^{J/2} e^{iπs} π² 2^{−J} (J+3)! · cos(πh) Γ(h−½) / (Γ(s+½) Γ(Δ+3)),
    M_{+J}[S^oddJ] =  8i(−1)^{(J−5)/2} e^{iπs} π² 2^{−J} (J−1)! · cos(πh) Γ(h−½) / (Γ(s+½) Γ(Δ−1)).

For `h > ½` (above the unitarity bound), `Γ(h−½)`, `Γ(s+½)`, `Γ(Δ+3)` and `Γ(Δ−1)` are finite and positive. So the zeros are
exactly the zeros of `cos(πh)`, namely `h = 3/2 + n`, i.e. `Δ₀ = J+3+2n`. Each is simple:

    ∂_Δ M_{−J}[S^odd] |_{Δ₀} = −8(−1)^{J/2} e^{iπs₀} π² 2^{−J}(J+3)! · (−(π/2) sin(πh₀)) · n! / ((J+1+n)! (J+5+2n)!),

with `sin(πh₀) = (−1)^{n+1} ≠ 0`. The `S^oddJ` case is the same with `(J+5+2n)! → (J+1+2n)!`. Exact values are asserted in
`oz1_out.txt` for `J ∈ {4,6,10,30}` and `{5,7,11,31}`, `n ∈ {0,1,2,7}` (value 0, derivative nonzero), e.g.
`∂_ΔM_{−4}[S^odd](Δ₀=7) = iπ³/34560`.

**(c) Absolute kernels.** With the ladder, ANGULAR_PROOFS §3.4 gives
`|K^odd_0|² = 3456π^{3/2} cos²(πh) 4^Δ [2^J/C(2J,J)] [Γ(s)/Γ(s+½)] [Γ(h−½)/Γ(h)] Γ(Δ−1)Γ(Δ−½)/Γ(Δ+3)²`. Near `Δ₀`:

    |K^odd_0|² = 864 π^{7/2} (Δ−Δ₀)² 4^{Δ₀} [2^J/C(2J,J)] [Γ(s₀)/Γ(s₀+½)] [n!/Γ(3/2+n)] Γ(Δ₀−1)Γ(Δ₀−½)/Γ(Δ₀+3)² + O((Δ−Δ₀)³),

and the same for `|K^oddJ_4|²`. This is H0 §15.1's expected phase shift made explicit:
`cos²(πh) = cos²(πτ/2) = sin²(π(τ±1−10)/2)`.

**(d) Control.** `S^{++}` and `S^{+−}` carry `1/Γ(1−h) = sin(πh)Γ(h)/π`. Their zeros are at `h = 1+n` (the even-twist `[TT]`), and
`|sin| = 1` at odd twist. The odd structures are nonzero at even twist (`|cos| = 1`). The two zero sets are disjoint and
interlaced; each family sits on the zeros of its own structures only.

---

## 3. All components, both helicities, and the first moments

### 3.1 For all `J`, given the helicity ladder [P | (L)]

KERNELS_D3 (R0), (R4) and ANGULAR_ROUTE §8 give the ladder. LADDER_ALLJ grades it [P] given `D_z·[K, E] = 0` ([D]), and checks it
exactly on the pipeline for `J ≤ 40`, both parities and both helicities. For the source `+2`:

    F_{−k} = F_{−J} · C(J,k)(4)_k (Δ−1−J)_J / ((Δ−1−k)_k (4)_J),   0 ≤ k ≤ J    (S^odd, and J = 0, 2),
    F_{4+k} = F_J · C(J−4,k)(4)_k (Δ−1−J)_{J−4} / ((Δ−5−k)_k (4)_{J−4}),   0 ≤ k ≤ J−4    (S^oddJ);

the other sector vanishes identically. At `Δ₀ = J+3+2n`:
* `(Δ₀−1−k)_k ≥ (2+2n)_k > 0`;
* `(Δ₀−1−J)_J = (2+2n)_J > 0`;
* `(Δ₀−5−k)_k > 0` for `k ≤ J−4`.

So every ratio is finite and positive. All components therefore share the simple zero of the extreme one (Z1); this is the
statement "the `m`-dependence is rational and positive" of KERNELS_D3 §5.7, now with its denominators checked at the zeros.

* **Helicity `−2`.** `F^{(−2)}_m = −F^{(+2)}_{−m}` for odd structures [D by parity covariance; DV J ≤ 40]. So the zeros are the
  same.
* **First moments.** LADDER_ALLJ (A) gives `M^orig = (i/6)(W − 18)H`, with `W` tridiagonal in `m` and polynomial in `Δ`
  ([P | [z·K,E] = 6iE⁽¹⁾ ([D])], DV `J ≤ 8`). So `M⁽¹⁾` vanishes wherever all `H_m` vanish, to at least the same order. The
  midpoint moment differs by the translation piece `(i/2)[(Δ+3) − B₁]H`, which is also linear in `H`.

### 3.2 Exact for every `J ≤ 30` and every `n` [P] (`oz_amps.py`, `oz2_zeros.py` → `cache/oz2_*.txt`)

**Method.** For fixed `J` the pipeline gives, exactly and with `Δ` symbolic, `M_m(s) = G_K(s)P_m(s)`. Here `P_m` is a polynomial
and `G_K` is the common prefactor `∝ 1/(Γ(e₁₄+K_a)Γ(e₁₃+K_b)Γ(e₃₄+K_B)Γ(e₃₄−½+K_γ))`. Hence

    ord_{s₀} M_m = #{Γ-factors of G_K vanishing at s₀} + mult_{s₀}(P_m).

For odd structures `e₁₃ + K_b = 11/2 − s + K_b`, so the `e₁₃` factor vanishes simply at every `s₀` beyond a threshold. Beyond the
largest rational root of `P_m`, the order is exactly 1. So **"simple zero for every `n ≥ 0`" is decided by finitely many exact
checks** per `(J, m)`. The script factorises `P_m` over `ℚ` and computes the order at every `n` up to that point, and the stable
order beyond it.

**Inputs.**
* Structures: the ANGULAR_PROOFS Theorem 1 vectors (conserved and Bose-symmetric for symbolic `J` [P]), specialised to `J`; for
  `J = 0, 2`, the read-only `kernels_lib.conserved` (symbolic `s`).
* Components: all `m`; sources `T(e₊)` and `T(e₋)`.
* Moments: zeroth, midpoint first moment and origin first moment, using the `α`-moment rule of `code/ladder/common.py`.

**Additional exact checks per `J`.**
* the structure's coefficient vector has gcd 1, so there is no common zero at any `Δ₀` (Z5);
* the zeroth moment, source `+2`, equals the closed-form ladder head × ladder factor of §3.1, as an identity in `Δ`
  (`gammasimp` ratio = 1), for **every** `m`;
* the source `−2` amplitudes are the mirror `P^−_m = −P^+_{−m}`;
* control: at the even-twist points, the odd amplitudes do **not** vanish (stable order 0).

**Results** (`ALL CHECKS PASSED: True` in every file):

| sector | `J` | components | zeroth moment | first moments |
|---|---|---|---|---|
| `S^odd` | 4, 6, …, 30, and 40 (§3.4) | `m = −J…0` (`h = +2`), `0…J` (`h = −2`) | simple zero at every `Δ = J+3+2n`, all `n` | order ≥ 1 at every `Δ = J+3+2n`; order 2 at isolated points (`J = 6`: `m = 0`, `n = 0`, origin moment) |
| `S^oddJ` | 5, 7, …, 31 | `m = 4…J` / `−J…−4` (first moments: from `m = 3` / `−3`) | simple zero at every `Δ = J+3+2n` | order ≥ 1; order 2 at isolated points (`J = 17`: `m = ±5`, `n = 1`) |
| odd, `J = 0` | 0 | `m = 0` | simple zero at every `Δ = 7+2n` | simple |
| odd, `J = 2` | 2 | `m = −2…0` / `0…2` | simple zero at every `Δ = 7+2n` | simple |
| controls `S^{++}`, `S^{+−}`, even `J = 0, 2` | 0–12 | all | simple zero at every even-twist `[TT]` point | order ≥ 1 |

The extra order at isolated points only strengthens the double zero. Per-file detail is in `cache/oz2_<name>_J<J>.txt`.

**Grade.** For each `J ≤ 30` (`≤ 31` for odd `J`), and for `J = 40`, this is exact computer algebra in `Δ`, covering every `n`.
So (Z2), (Z3) are **[P] for these `J`**. The only input is the structure count (St) ([S]; [C] for `J ≤ 30`), which identifies
the canonical vector as the unique odd structure. The computation also upgrades the ladder (L) for the odd structures, from
rational-`Δ` checks to identities in `Δ`, at these `J`.

### 3.3 (Z5) Regularity of the basis [P]

The double zero is a statement in a structure basis that is finite and nonvanishing at `Δ₀`. In such a basis the GFF OPE
coefficient is finite. (PARITY_ODD_ROUTEA §12.5 warns that CHL's `X^O` basis carries `tan(πτ/2)`, whose pole eats the zero.)

* **`S^odd`.** Its normalising entry is `ε₁₄X₁₄V₁V₃V₄^{J−2}H₃₄ = J(J−1)(J−2)`. This is independent of `Δ` and nonzero for
  `J ≥ 4`, so the structure never vanishes [P, symbolic `J`].
* **`S^oddJ` and `J = 0, 2`.** For each computed `J`, the coefficient polynomials have gcd 1, hence no common root [P per `J`].
  Individual entries can vanish: for example, the `ε₁₃X₁₃H₁₃V₄^J` entry of `S^oddJ` vanishes at `J = 5`, `Δ = 12`, but the others
  do not.
* The 12 odd monomials of the repo basis are linearly independent (KERNELS_D3 §3.1), so a nonzero vector is a nonzero
  structure.

### 3.4 `J = 40`

The same exact analysis at `J = 40` for `S^odd` (`cache/oz2_odd_J40.txt`; about 11 CPU-min) passes every check:
* simple zeros of all 41 zeroth-moment components at every `Δ = 43+2n`, for both sources;
* order ≥ 1 for both first moments;
* the closed-form × ladder identity, exactly in `Δ`, for all `m`;
* the `h = −2` mirror.

This is the largest `J` checked in this file. Together with LADDER_ALLJ's rational-`Δ` checks to `J = 40`, the ladder used in
§3.1 is tested exactly through `J = 40`.

---

## 4. `J = 0` and `J = 2` [P]

These are single structures with single sectors (`m ≤ 0` for source `+2`), computed exactly with symbolic `Δ` (§3.2).

* **`J = 0`.**
  * `M₀` has simple zeros exactly at `Δ = 7+2n`, for both moments and both sources.
  * Check against KERNELS_D3 §5.2: `12|M₀|²/N₀ ÷ K_o² = 1.0` at `Δ = 37/7, 113/9, 51/4`, with
    `K_o² ∝ 1/Γ((7−Δ)/2)² = cos²(πΔ/2)Γ(Δ/2−5/2)²/π²`.
  * There is **no** zero at `Δ = 5` or `Δ = 3` (`Γ(1)`, `Γ(2)`). That matches (Z0): the lowest parity-odd scalar `[TT]` is
    `ε^{αβγ}T_{αδ}∂_βT_γ{}^δ`, at `Δ = 7`.
* **`J = 2`.**
  * `M_{0,−1,−2}` have simple zeros exactly at `Δ = 7+2n`.
  * The pipeline equals the KERNELS_D3 §5.3 closed form `M₀[J=2,odd] ∝ (Δ−3)(Δ−2)/(Γ(2+Δ/2)Γ((Δ+3)/2)²Γ((7−Δ)/2))` times the
    `J = 2` ladder, exactly in `Δ`, for every `m`.
  * There is **no** zero at `Δ = 5` (`h = 3/2`). That matches (Z0): an odd spin-2 `[TT]` needs at least one derivative, so
    `Δ ≥ 7`.

So in every sector the zero set of the parity-odd amplitude above the unitarity bound is **exactly** the set of possible
parity-odd `[TT]` dimensions: `{J+3+2n}` for `J ≥ 4`, `{7+2n}` for `J = 0, 2`. There are no zeros at other points.

---

## 5. The light side (task 3)

### 5.1 The parity-odd part of the `T` block is not at a zero [P] (`oz3_light.py` → `cache/oz3_out.txt`)

`𝒪 = T` has `J = 2`, `Δ = 3`, so `h = ½` and `τ = 1`. The odd amplitude carries `1/Γ(3/2−h) = 1/Γ(1) = 1` there, so it is
**not** at a zero. The exact pipeline on the unique parity-odd `⟨TTT⟩` structure at `s = 5/2` gives:

| object | result |
|---|---|
| zeroth moment, source `+2` | only `M_{−2} ≠ 0` (`= −2.58386i` per unit coefficient): pure flip `Δh = 4`. `M_{+2} = M_{±1} = M₀ = 0` (Ward identity) |
| first moments, source `+2` | diagonal element `M⁽¹⁾_{+2} = 0` (origin and midpoint) |
| source `−2` | mirror image |
| the generic `J = 2` odd family as `Δ → 3` | `M_{−2} → iπ³/2` (finite, nonzero); `M_{−1}, M₀ ∝ (Δ−3) → 0` (conserved-current selection). The `⟨TTT⟩` odd amplitude is `−1/6` × this limit: a pure normalisation |

**Consequences.**
* `α_T` enters the one-point matrices `Ĥ`, `M̂` only off-diagonally. This was graded [D] in C1_FINAL and INDEP_CHECK3, "no free
  realisation". It is now an exact pipeline statement, [P] at `Δ = 3`.
* So the `O(α_T)` linear light term is flip, and is killed by the diagonal `EE` pairing (L3).
* The quadratic term is flip × flip = preserving. It gives the `α_T²/16` in `|c|² = (t₄²+α_T²)/16` of the `k = 2` light side,
  which is the term behind the new `α_T` bound.

**PARITY_ODD_ROUTEA §12.5 aside.** That section noted that the factor `2cos²(πτ/2)` vanishes at `τ = 1`, so a parity-odd `⟨TTT⟩`
would seem to need "the residue treatment". The amplitude shows it does not. The `cos²` zero of `|K^odd|²` at `h = ½` is cancelled
by the `Γ(h−½)` pole: at the unitarity bound the conserved-current selection moves the weight to the top component `m = −2`,
which is finite and nonzero. The `2cos²(πτ/2)` factor language is valid only for `h > ½`.

### 5.2 Parity-violating admixtures [D]

In a parity-violating theory, parity is broken at `O(1/c_T)`: GFF is parity-invariant, and the odd cubic data `α_T` enter
exchange diagrams with two vertices. In the canonical normalisation of §1.3:

* **Odd-twist `[TT]`.** They acquire parity-even couplings `λ_E = O(1)` (vs `λ_O = O(C_T)`). `K_E ≠ 0` at odd twist, so these
  enter `|â|²` as `λ_E²|K_E|²/C_T = O(1/C_T)`.
* **Even-twist `[TT]`.** They acquire odd couplings `λ_O = O(1)`, which enter as `O(1/C_T)` in the same way.
* **Cross terms.** These are `λ_Eλ_O K_E K_O/C_T = O(1)·O(C_T)·O(1)·O(1/C_T)/C_T = O(1/C_T)`. They are absent anyway from
  `|â_{m,2}|²`, by the relative phase `i` (KERNELS_D3 §5.1), and from the heavy kernels, by parity (PARITY_ODD_HEAVY).
* **Mixing.** Even-twist and odd-twist `[TT]` of the same `J` are non-degenerate at `c_T = ∞`, so there is no degenerate mixing.

All of these are relative `O(1/c_T)` to the `T` block, i.e. `O(c_T^{−2})` in the functional normalisation. This is the same
counting as PARITY_ODD_ROUTEA §12.5 ("`λ_E² = O(1/c_T²)`"). So no parity-odd `[TT]` contribution survives at `O(1/c_T)`. In
particular, contact interactions of either parity change only double-trace data and drop out ("contacts drop" in L1 now holds
in both parities).

### 5.3 Check against H0 §15.1's concern

The concern was that odd-twist `[TT]` sit at `sin²(π(τ−10)/2) = 1`, unless the parity-odd structures shift the phase by `π/2`.
They do. The exact amplitude-squared of every parity-odd structure is `cos²(πτ/2)` × (positive, finite) (§2(c), §3.2). That is
`sin²(π(τ±1−10)/2)`, as H0 anticipated. The factor appears in the amplitude itself, before any multivaluedness argument.

The polynomial weights `f₂`, `f₃` do not enter this argument. They multiply the light-ray product, so they become first moments
(§1.2), and those are covered (Z3).

### 5.4 Why the zeros must be there (GFF remark) [D]

For a generalised free `T`, `E(n)T(p)|0⟩ = [E(n), T(p)]|0⟩` is a c-number times `|0⟩`, so it vanishes at `p ≠ 0`. Therefore
`Σ_a λ^GFF_a M^{(a)}_m(Δ_𝒪) = 0` for every GFF `[TT]` primary `𝒪` and every `m`.

For odd `J`, and for the parity-odd sector at even `J`, there is a single odd structure, so this forces `M^odd_m(Δ_GFF) = 0`
whenever the GFF coefficient is nonzero. Equivalently, without the zeros the sum rules of KERNELS_D3 §6 would read
`O(C_T) ≤ 1`. The computed zero sets coincide exactly with the GFF dimensions (§4), with no zero at `J = 2`, `Δ = 5`, where
there is no GFF operator. This is the expected pattern. It is used here only as a consistency check, because the nonvanishing of
each GFF coefficient is not shown.

---

## 6. What is not shown

* **All `m` for `J > 30`** rests on the helicity ladder (L), i.e. on light-transform conservation `D_z·[K, E] = 0` ([D],
  LADDER_ALLJ). The extreme components and the zero set are [P] for all `J`. The ladder factors at the zeros are proved finite and
  positive for all `J`.
* **First moments for `J > 30`** rest on (L) and `[z·K, E] = 6iE⁽¹⁾` ([D], DV `J ≤ 8`).
* **The structure count (St)** for `J > 30` ([S]; [C] `J ≤ 30`), as everywhere.
* **The factorisation of `Ω_k` through light-ray one-point matrices at `x₀ = ∞`** (C1_FINAL §2.2) is [D], as for the `T` block.
  The parity-even half of L1 has the same dependence through CLS-D's factor.
* **Uniformity** of the `O(c_T^{−2})` estimate over the infinite `[TT]` towers. This is inside A2, as for the parity-even towers
  (Supplement "Order of limits"); nothing new.
* **Multi-trace operators** (`[TTT]`, …) are not at these zeros. Their entry at `O(c_T^{−2})` is the standard large-`c_T` statement
  already used in the Supplement; it is not re-examined here.

---

## 7. Supplement and Letter text

**Supplement, Table S2, row L1, grade column.** Replace
`[D]; the parity-odd, odd-twist $[TT]$ double zeros [D-sketch]` by:

```latex
[D]; double zeros in both parities, from the light-ray amplitudes (not CLS-D's factor): the parity-odd structures give
$\langle\mathcal O|\mathcal E|T\rangle\propto\cos(\pi h)\,\Gamma(h-\tfrac12)$, $h=\tfrac{\Delta-J}2$, with simple zeros exactly at the
parity-odd $[TT]$ dimensions $J+3+2n$ ($J\ge4$) and $7+2n$ ($J=0,2$), so $|K|^2\propto\cos^2(\pi\tau/2)$ [P all $J$, extreme
components]; all $m$, $h=\pm2$ and first moments [P] for $J\le30$, else given the helicity ladder [P given D]; odd-$\tau$ $[TT]$ and
parity-violating admixtures enter at $O(c_T^{-2})$ [D] (\texttt{PARITY\_ODD\_DZ.md})
```

**Letter** (replaces "The cubic implication (1) is derived for parity-invariant theories; for parity-violating ones one step (the
double zeros of parity-odd, odd-twist double traces) is sketched."):

> The cubic implication (1), including the bound on $\alpha_T$, is derived for parity-invariant and parity-violating theories
> alike: the light-ray amplitudes of parity-odd double traces vanish at their dimensions, which we show in closed form.

**Other places that carry the same caveat** (for whoever edits; not touched here):
* Supplement l. 120–122 ("Result and other phrases");
* l. 332 (L1), l. 334 (L3: "given L1's double zeros (parity-odd, odd-twist … [D-sketch])");
* l. 353 ("for parity-violating theories net [D-sketch] through L1");
* l. 390–391 (Net grades);
* l. 675–676 (Sec. c1, chain (ii));
* l. 761 (Table c1, row (ii)).

In each place the parenthetical "(parity-odd, odd-twist … [D-sketch])" can be dropped, or replaced by "[D] in both parities".

---

## 8. Reproduction and code map (`ads4/code/oddz/`)

```
cd ads4/code/oddz;  ulimit -v 3000000; export OMP_NUM_THREADS=1
python3 -B oz1_extreme.py                         # §2: ap5 re-run (read-only), reflection, zeros, derivatives, ladder factors  ~30 s
python3 -B oz2_zeros.py odd 0 2 4 6 ... 30 40     # §3.2, §4: exact zero orders, all m, h = ±2, both moments, closed-form identity
python3 -B oz2_zeros.py oddJ 5 7 ... 31           #   (per J: ~1 min at J = 20, ~3-4 min at J = 28-31, ~11 min at J = 40)
python3 -B oz2_zeros.py pp 0 2 4 ... 12; python3 -B oz2_zeros.py pm 4 6 ... 12     # controls (parity-even, even twist)
python3 -B oz3_light.py                           # §5.1: parity-odd <TTT> at Delta = 3; J = 2 family limit
```

| file | content |
|---|---|
| `common.py` | read-only imports, cache redirect and write guard, ANGULAR_PROOFS symbolic-J vectors at fixed `J`, first-moment rule (as `code/ladder/common.py`) |
| `oz_amps.py` | symbolic-`s` amplitudes `M_m = G_K P_m`, all `m`, sources `±2`, zeroth/mid/origin moments; cached `amps_<name>_J<J>.pkl` |
| `oz2_zeros.py` | exact zero orders at all `n`; gcd/regularity; closed form × ladder identity; `h = −2` mirror; controls |
| `oz1_extreme.py` | symbolic-`J` statements of §2 |
| `oz3_light.py` | §5.1 |
| `cache/oz1_out.txt`, `cache/oz2_<name>_J<J>.txt`, `cache/oz3_out.txt`, `cache/run*.log` | outputs |
