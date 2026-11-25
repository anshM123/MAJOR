# C1_FINAL: `c₁ = 0`. Premise F₂ proved for the `k = 2` weight `f = v′ − u′`; the bounds lose their last caveat

*25 Nov 2026, research agent. Unreviewed, not committed. No existing file was edited. New code: `ads4/code/c1final/`
(`c1f_affine.py`, `c1f_freefields.py`, `c1f_onepoint.py`, `c1f_compose.py`, `c1f_bounds.py`); outputs in `ads4/code/c1final/cache/`
(git-ignored by its own `.gitignore`). Read-only imports: `code/c1/c1_quadratic.py` and, through it, `code/dgscan/bounds.py`.
Literature: CLS-D = Chang–Landau–Simmons-Duffin, arXiv:2311.04271, LaTeX source read as data (eqs. `superconvergence_vdirection`,
`light_transform_definition`, `scalarfunctional_definition`, `cross_ratios`, `subtraction_factor_evenjprime`,
`CFT_building_blocks_VH`). Compute: about 2 CPU-min in total, one process, one thread, `ulimit -v 3000000`.*

**Grades.** [P] proved, including exact algebra; [D] derived; [D-sketch]; [DV] derived and checked numerically or symbolically;
[Q] quoted; [NE]; [M]; [A] assumption; [I] open. A chain has the grade of its weakest link.

---

## 0. Verdict

**`c₁ = 0` exactly. Grade [D].** This is the same grade as the implication "(Pos) ⇒ (1)" that it feeds. Before this note it was
[D-sketch given F₂], with F₂ graded [A] (C1_VALUE §4.2).

* **Route (A) closes.** Premise F₂ is proved for the project's `k = 2` weight `f = v′ − u′`. **F₂: [A] → [D + DV].** The proof is
  the `k = 2` analogue of H7's factorisation F and needs two new facts.
  * **The weight is affine.** On light rays based at `x₀ = ∞`, `f` is affine in the two retarded times and has **no `u₁u₃` term**
    [P, exact]. So the `k = 2` functional contains only the bare light-ray operator `E(n) = ∫du T` and its first moment
    `E⁽¹⁾(n) = ∫du u T`, **at most one first moment per term**.
  * **The first moment has the same helicity pattern as `E`.** On the stress-tensor family,
    `E⁽¹⁾(n) = ½{E(n), −i n̄·D}`, the symmetrised product of `E(n)` with the time-of-arrival operator (`D` is the frame-covariant
    momentum derivative). This holds for every parity-even `⟨TTT⟩` [D + DV]. Hence:
    * the Einstein part of `E⁽¹⁾` is helicity-diagonal;
    * its `t₄` part is purely helicity-flipping (`Δh = ±4`);
    * there is no diagonal `O(t₄)` element.
* **The conclusion.** The EE `(2,4)` dual structure contracts the two external states by the identity [P]. Hence the `O(t₄)`
  `s`-channel `T`-family term of `Ψ_{2,ν}[EE]` vanishes **pointwise**, before any integration [D]. The end-to-end free-field check
  gives zero to `≤ 1.8·10⁻¹¹` of the Einstein integrand. The controls are `O(1)`.
* **The Euclidean-vs-Lorentzian helicity point of C1_VALUE §4.2 is resolved (§4) [D-sketch].** It is not load-bearing.
  * The functional, F₂ and the EE zero are all statements about Lorentzian state helicities.
  * Read in those helicities, the bulk remnant `⟨T⁺T⁺T⁺⟩_{W₊}·Δ/⟨TT⟩` has an `s`-channel `T` cut of the form
    W-vertex (flip) × Einstein vertex (diagonal) = a `(2,4)` flip.
  * Its "preserved ⊗ preserved" label is a label of Euclidean momentum space (SD/ASD), and does not carry over to the cut.
* **Route (B), the direct evaluation (§5).** In the project's `Ω_OPE` prescription the direct evaluation on the explicit
  energy-channel exchange (W₊ vertex × (E − κCS) vertex) is the same pointwise computation. It gives `c₁ = 0`. The numerical floor
  is `≤ 2·10⁻¹¹` (finite differences); there is no other uncertainty.
  * Its value in H7's Lemma-R prescription follows through the identity `Ω_OPE = Lemma R` on tree terms (LIGHTSIDE_DIST [D]).
  * An independent Regge-limit (Witten-diagram) evaluation of the remnant was **not** done.
  * The normalisation (`−Ĉ′_EE|_E = 16πG/(ν²+9/4)`, `8πG = 256/c_T`, `√(t₄²+α_T²) = 240|ĝ₃|`) does not enter a zero.
* **Effect on the headline numbers.** The numbers do not change. The `c₁` caveat and the `c₁` error terms are removed: `x = 0` in
  `R₀(x + √(1+x²))`, and the `c₁²` term of Q is zero.

| `Δ_g` | **Result D**: `√(t₄²+α_T²) ≤` | `t₄` interval | `|ĝ₃|²Λ⁸ ≤` | **Result Q**: `g̃₄Λ⁶/8πG ≤` | removed (per unit `|c₁|`) |
|---|---|---|---|---|---|
| 300 | **3.13206·10⁻⁶** (`K = 1`; 3.13326·10⁻⁶ with `K = k(0)(1+δk)`) | symmetric `±R₀` | 10734.61 | **72.0406** | D: `x = 2.160·10⁻³` (+0.216 %); Q: `+6.81·10⁻⁵c₁²` |
| 1000 | **2.89126·10⁻⁸** (2.89136·10⁻⁸) | symmetric `±R₀` | 14339.59 | **92.8135** | D: `x = 1.670·10⁻⁴` (+0.017 %); Q: `+5.44·10⁻⁷c₁²` |

**Other consequences.**
* **The whole `O(t₄)` light term of `Ĉ′_EE` vanishes for every `ν`.** This is a stronger statement than "the constant `c₁` vanishes".
  * The value of `c₁` therefore no longer rests on `J₀(E×W³) = 2`, on polynomiality, or on the shock / ASD-closure analysis.
    These remain as consistency statements.
* **An `O(α_T)` constant vanishes as well [D].** That is the parity-odd analogue (§2.3(d), §3).
* The `|ĝ₃|²` light term (`W × W` = flip × flip = preserving) is kept by EE, as it must be. This is a control (§6).
* **The zero is specific to EE.** For a dual structure that sees `(2,4)` flips (for example `q²`, C1_COMPUTE §3.2) the same `O(t₄)`
  term would be non-zero. The project's `k = 2` pieces of both Ω^EE and Ω^Q_up use EE.

---

## 1. What has to be shown

`c₁` is the coefficient of the `O(t₄)` light term `16πG c₁t₄` of `−Ĉ′_{EE,ν}` (C1_COMPUTE §1). It is one number for Results D and Q.
The chain used here has three links.

* **(i) Prescription [D, inherited].** H7's Lemma-R value of a tree term equals its term-wise `s`-channel value `Ω_OPE` (LIGHTSIDE_DIST:
  `Ω_k[X] = D_LR[X] + 𝒜_k[X]`, with `D_LR = 0` on microcausal terms).
* **(ii) Only the `T` block [D, inherited].** At `O(1/c_T)` the `k = 2` term-wise light side is the `T` block alone:
  * the identity drops because `L[T]|Ω⟩ = 0`;
  * the `[TT]` double traces have double zeros under the unsubtracted `k = 2` functional (CLS-D §3.5.2 [Q]; EE_FOUNDATIONS §1.4).
    The parity-odd caveat of H0 §15.1 applies.
  * This is the same input that makes the Einstein term the `T`-block pole `16πG/(ν²+9/4)`.
* **(iii) The EE action on the `T` block's `E × t₄` structure vanishes.** This is what F₂ delivers (C1_VALUE §4.2). It is proved
  here (§2).

---

## 2. Proof of F₂

### 2.1 The `k = 2` weight is affine in the light-ray times [P] (`c1f_affine.py` (a))

CLS-D's spinning `k = 2` functional (`J_e = 2`) inserts `f⁺_{2;2} = v′ − u′` into `Ψ⁺ = ∫⟨[O₄, L[O₃](x₀,z₃)][L[O₁](x₀,z₁), O₂]⟩`.
Here `u′ = x₁₂²x₃₄²/(x₁₃²x₂₄²)` and `v′ = x₁₄²x₂₃²/(x₁₃²x₂₄²)`.

**Gauge.** Fix `x₀ = ∞`, i.e. spatial infinity. This is spacelike to `x₂` and `x₄`, as CLS-D's integration region requires. The
light rays are then the generators of future null infinity: `x_{1,3} = (u + r, r n_{1,3})`, `r → ∞`, and `L[T](∞, (1,n)) ∝ E(n)`.

**Result.** With `n̄·x := x⁰ − n·x⃗`, exactly:

    f = 2/((1 − n₁·n₃) x₂₄²) · [ u₁ n̄₃·x₄₂ + u₃ n̄₁·x₂₄ + (n̄₁·x₄)(n̄₃·x₂) − (n̄₁·x₂)(n̄₃·x₄) ].

* The `u₁u₃` coefficient is `u₄u₂ − u₂u₄ = 0` identically. In the sympy check, `∂²f/∂u₁∂u₃ = ∂²f/∂u₁² = ∂²f/∂u₃² = 0`.
* The finite-`r` value converges to this limit as `1/r` (`4.2·10⁻³`, `4.2·10⁻⁶`, `4.2·10⁻⁹` at `r = 10³, 10⁶, 10⁹`).

So at fixed `(n₁, n₃, x₂, x₄)` the light-ray product in `Ψ⁺` is

    f·L₃L₁ ∝ a·E₃E₁ + b·E₃E₁⁽¹⁾ + c·E₃⁽¹⁾E₁,      E⁽¹⁾(n) := ∫du u T(u,n),

with `a, b, c` scalar functions of `(n₁, n₃, x₂, x₄)`. **No term contains two first moments.** The `(1↔3)` image has the same form.
* The `x₂`, `x₄` dependence of the coefficients only changes the scalar kernel multiplying the `(2,4)` structure.
* The celestial and dual structures carry the `ν` dependence. They do not depend on `u₁`, `u₃`.

### 2.2 Factorisation through the `T` family [D]

Since `L[T]|Ω⟩ = 0`, the double commutator equals `⟨O₄ L₃ L₁ O₂⟩`. Its `s`-channel `T` block is the insertion of the projector `P_T`
onto the `T` conformal family between `L₃` and `L₁`. `O₂|Ω⟩` and `⟨Ω|O₄` already lie in the `T` family. Hence, exactly:

    T-block term = ⟨Ψ₄| A₃ P_T A₁ |Ψ₂⟩,    A ∈ {E, E⁽¹⁾},    P_T = Σ_h ∫d³P |P,h⟩⟨P,h|/n(P).

Here `|P,h⟩ := ε_h(P)·T(P)|Ω⟩`, with `h = ±2` the rest-frame helicity of the timelike momentum `P` and `ε_h = e_h e_h` in the
boosted frame `e_h(P) = Λ(P)e_h(rest)`. This is the statement of F₂ except for the helicity pattern of the matrices, which §2.3
supplies. It is exactly H7's F for `k = 3` (`f = −1`, only `E₃E₁`), extended by the first moments.

### 2.3 The one-point matrices of `E` and `E⁽¹⁾` on the `T` family [D + DV]

**(a) Form, from translations [P].** Translations shift the retarded time, so `[P^μ, E⁽¹⁾(n)] = −i n̄^μ E(n)` with `n̄ = (1, n)`.
Therefore, with normalised components:

    (E f)_h = Ĥ f,     (E⁽¹⁾ f)_h = −i Ĥ n̄·∂f − (i/2)(n̄·∂Ĥ) f + M̂ f.

`Ĥ(P,n)` is the HM one-point matrix. `M̂(P,n)` is the only new datum.

**(b) Lemma (helicity pattern) [D + DV].** For every parity-even `⟨TTT⟩`:
* `M̂ = ½{Γ(n̄), Ĥ}`, where `Γ = diag(Γ₊, Γ₋)` is the frame connection along `n̄`, with `Γ₋ = −Γ₊`.
* Equivalently, **`E⁽¹⁾(n)|_T = ½{E(n)|_T, −i n̄·D}`**, with `D = ∂ + iΓ` the frame-covariant derivative.

Consequently, with `Ĥ = Ĥ_E + (t₄/4)Ĥ_W`, where `Ĥ_E` is diagonal and `Ĥ_W` is pure flip (Q8, `d = 3` HM):
* the Einstein part of `E⁽¹⁾` is helicity-diagonal;
* its `t₄` part is pure flip, with `M̂_W = 0`.

**(c) Proof.**
1. **The free fields span the structures.**
   * At fixed `c_T` the parity-even `⟨TTT⟩` of a 3d CFT is a two-dimensional space. It is spanned by the free real boson and the
     free Dirac fermion, which sit at the two ends `t₄ = ±4` of the `d = 3` HM interval (`t₄,B = −t₄,F`; the sign is a convention,
     and §2.3(e) reproduces `|t₄| = 4`).
   * The normalised matrix elements are linear in the structure: `Ĥ = Ĥ_E + (t₄/4)Ĥ_W`, `Ĥ_E = (Ĥ_B+Ĥ_F)/2`,
     `Ĥ_W = (Ĥ_B−Ĥ_F)/2` (with `t₄,B = +4`), and likewise for `M̂`.
   * So it suffices to compute both free theories exactly.
2. **The light-ray operators at `ℐ⁺`.** In a free theory `T|Ω⟩` is a two-particle state with wave function
   `ψ_h(k₁,k₂) = ε_h(k₁+k₂)·t(k₁,k₂)`, where:
   * boson: `t = k₁k₂ + k₂k₁ − η k₁·k₂`;
   * fermion: `t = ū(k₁)γ_{(μ}(k₁−k₂)_{ν)}v(k₂)`.

   From the mode expansion at `ℐ⁺` (stationary phase; the improvement terms are total `u`-derivatives and drop from `∫du` and
   `∫du u`):
   * `E(φ) = (1/8π²)∫dω ω a†a`;
   * `E⁽¹⁾(φ) = (1/8π²)∫dω ω (−i/2)(a†∂_ω a − ∂_ω a† a)`.

   Only the diagonal and first derivative of the symmetric `ω, ω′` kernel enter, so the fermion operator has the same form.
3. **Single phase-space points.** The matrices are then single phase-space points, `k₁ = ωn̄` with `ω = P²/(2P·n̄)`, and
   `M ∝ (−i/2)(ψ*_{h′}∂_ωψ_h − ∂_ωψ*_{h′}ψ_h)`.
4. **The mechanism (exact in both theories).** `∂_ω ψ_h|_{k₂} = (iΓ_h + 1/ω) ψ_h`.
   * The frame derivative of `ε_h` is the SO(2) connection plus a longitudinal piece. The longitudinal piece drops by
     conservation, `P^μ t_{μν} = 0`.
   * `t` is linear in `k₁`. For the fermion the `√ω` of `u(k₁)` and the `(k₁−k₂)` factor add up to `1/ω` as well.
   * There is no admixture of `ψ_{−h}`. The real part cancels in the antisymmetric combination, which leaves `M = ½{Γ, H}`. ∎

**(d) Diagonal elements of any structure, parity-odd included [D].**
* In the rest frame of `P`, rotation covariance makes the diagonal elements of `Ĥ` and `M̂` independent of `φ`.
* Their angular integrals are fixed by the Ward identities:
  * `∫dφ n̄^μ E(n) = P^μ`;
  * `∫dφ E⁽¹⁾(n) = D` (the dilatation flux, `ξ^u = u` at `ℐ⁺`).

  Both depend only on `c_T`.
* So **no `⟨TTT⟩` structure other than the Einstein one has a diagonal element**, in `E` or in `E⁽¹⁾`. This covers `α_T` (the
  parity-odd structure), whose flip elements are harmless (§2.5).

**(e) Checks (`c1f_onepoint.py`; rest frame plus three boosted `P`, five directions each).**
* `|Ĥ_B,diag − Ĥ_F,diag| ≤ 1.9·10⁻¹⁵`.
* `|Ĥ_B,flip + Ĥ_F,flip| ≤ 2.4·10⁻¹⁵`.
* `|Ĥ_flip|/Ĥ_diag = 1` to `2·10⁻¹⁶`. This is `|t₄| = 4`, with the HM phase `e^{−4iφ}` (for example `arg = −1.48` at `φ = 0.37`).
* `|M̂_flip| ≤ 6·10⁻¹⁷` (B and F).
* `|M̂_B,diag − M̂_F,diag| ≤ 1.3·10⁻¹⁰`.
* `M̂ = ½{Γ,Ĥ}` to `2.6·10⁻¹⁶`.
* Ward identities:
  * momentum flux `= P^μ` to `2·10⁻¹⁵`;
  * flip moment `≤ 2·10⁻¹⁶`;
  * `∫dφ M̂_diag ≤ 3·10⁻¹²` (dilatation; scale-invariant frames).
* Sanity: `k̸u = 0` and `P^μt_{μν} = 0` to `4·10⁻¹⁶`.

### 2.4 The EE `(2,4)` contraction is the identity [P] (`c1f_affine.py` (b))

With CLS-D's building blocks, `G := H₂₄ − 2V_{2,04}V_{4,02} → x₂₄²(z₂·z₄)` exactly as `x₀ → ∞` (the ratio is 1). So
`EE ∝ x₂₄⁴(z₂·z₄)²` times `(2,4)`-polarisation-independent factors. The functional therefore pairs the two external states through

    W_{h₄h₂}(P₄,P₂) = ε_{h₄}(P₄) · ε*_{h₂}(P₂),

which has two properties:
* `W = δ_{h₄h₂}` at `P₄ = P₂`;
* `∂W` at coincidence `= iΓ δ`, because the frame connection is diagonal and the longitudinal parts are killed by transversality.

**`W` is diagonal through first order. Only first order is ever needed (§2.1).**

The embedding structures `H`, `V` are invariant under `(X₄, Z₄) → (−X₄, −Z₄)`. So writing the bra at `x₄⁺`, as CLS-D do, introduces
no reflection of the polarisation.

### 2.5 Theorem (F₂) and the zero [D]

**Theorem (F₂).** For the `k = 2` CLS-D functional with `f = v′ − u′`, the `s`-channel `T`-family term factorises as
`Σ_h ⟨Ψ₄|A₃|P,h⟩⟨P,h|A₁|Ψ₂⟩`, with `A ∈ {E, E⁽¹⁾}` and at most one `E⁽¹⁾`. Its one-point matrices are helicity-diagonal for Einstein
and pure flip at `O(t₄)` (and at `O(α_T)`).

**Proof.** §2.1–§2.3.

**Corollary.**
* At `O(t₄)` the `T`-block term is a sum of products (Einstein diagonal) × (flip). Each carries at most one derivative, acting on
  scalars, on the frames (diagonal connection) or on the other matrix (pattern preserved).
* It is therefore a `(2,4)` flip, `h₄ = −h₂`.
* Contracted with `W`, which is diagonal through first order, it vanishes **pointwise** in `(P, n₁, n₃, x₂, x₄)`.
* So it vanishes for every `ν`, every celestial and dual scalar kernel, and after the de-shift.
* The same holds separately for each cross term:
  * W at the `(12T)` vertex and E at the `(T34)` vertex;
  * and the mirror arrangement. ∎

**End-to-end check (`c1f_compose.py`).**
* **Setup.** Free-field `Ĥ`, `M̂` and finite-difference derivatives. Eight random timelike `P`, random `n₁, n₃`, random complex
  Gaussian external profiles, a full Cartesian sum over the `(2,4)` polarisation basis, and the three structures
  `E₃E₁`, `E₃E₁⁽¹⁾`, `E₃⁽¹⁾E₁`.

| quantity (pointwise in `P`) | value |
|---|---|
| `O(t₄)` term `X = I[E₃,W₁] + I[W₃,E₁]` under EE, and each cross term separately, relative to the Einstein term `I[E,E]` | **`≤ 1.75·10⁻¹¹`** (finite-difference floor; `≤ 9·10⁻¹⁶` for `E₃E₁`) |
| control: same `X` with an improper reflection in the `(2,4)` contraction (helicity-flip sensitive) | `2.7·10⁻³ – 1.5`, typically `0.2 – 1.4` |
| control: Einstein term `I[E,E]` under EE | non-zero (`5·10⁻⁴ – 0.6`) |
| control: `|ĝ₃|²` term `I[W,W]` under EE (flip × flip = preserving) | `0.03 – 1.0 × I[E,E]`, non-zero |

---

## 3. `c₁ = 0` and its grade

**The chain.**
* (i) `c₁` (Lemma R) `= Ω_OPE[O(t₄) tree]` [LIGHTSIDE_DIST, D].
* (ii) `= ` the `T` block's `E × t₄` term [DZ: `[TT]` double zeros, inherited].
* (iii) `= 0` under EE [F₂ (§2.5) + identity contraction (§2.4)].

**Hence `c₁ = 0` [D].**

**Inputs of grade [D] or weaker, all inherited except the first two:**
* the `ℐ⁺` mode expansion and the flux identities of §2.3(c)–(d) [D];
* the gauge `x₀ = ∞` and the reading of CLS-D's `x₄⁺` (§2.4) [D];
* LIGHTSIDE_DIST's identity (i) [D];
* (DZ) [D, with the H0 §15.1 parity-odd caveat];
* Q8 [Q + DV];
* convergence of the light-ray integrals on the `T` block. This is the same as for the Einstein term.

Everything specific to this note is [P] or [DV].

**What `c₁ = 0` now does not need.** It is stronger than "`c₁ = 0`":
* C1_COMPUTE §3.3 (the `u`-channel zero);
* DSKETCH §2 (`J₀(E×W³) = 2`, polynomiality, Barnes map);
* the ASD-closure lemma.

Each of these is now a consistency check rather than a link. The `u`-, `t`-channel and contact diagrams have only double traces in
the `s`-channel OPE (term-wise zero by DZ). The `s`-channel exchange's value is its `T`-block cross term, which is zero. So **every
individual `O(ĝ₃)` diagram has zero EE light value**. This includes the energy-channel remnant of C1_VALUE/DSKETCH, `b₊ = 7/300 ≠ 0`
notwithstanding: its Euclidean value is non-zero, but its EE-projected arc vanishes.

**Parity-odd analogue.** An `O(α_T)` constant (C1_COMPUTE §3.5 last bullet) also vanishes [D].
* `α_T` has no diagonal element (§2.3(d));
* the Einstein `E⁽¹⁾` has no flip (§2.3(b)).

The step specific to `α_T` is [D] only: no free field realises the parity-odd structure.

---

## 4. The Euclidean vs Lorentzian helicity point (C1_VALUE §4.2, "unresolved bookkeeping") [D-sketch; not load-bearing]

* **The two statements in tension.** C1_VALUE noted that the bulk remnant sits in Euclidean "preserved ⊗ preserved" (legs 1, 2 SD;
  3, 4 ASD), whereas F₂ makes its `T`-family term a `(2,4)` flip.
* **How they are reconciled.** The functional, the HM matrices and the EE zero are statements about **Lorentzian state
  helicities** (`h` of timelike `P` in its rest frame). Read vertex by vertex in that language:
  * the `W₊` vertex `⟨T₁⁺T₂⁺T⁺⟩` is the `t₄ ± iα_T` part of the one-point matrix, i.e. a **flip** (Q8 [Q], reproduced in §2.3(e));
  * `Δ = ⟨T⁺T⁻T⁻⟩_E − κ⟨…⟩_CS` has Lorentzian content equal to the Einstein vertex. The CS piece is ultralocal and has no Wightman
    or light-ray matrix element. The Einstein vertex is **diagonal**.
  * So the remnant's `s`-channel `T` cut is flip × diagonal, a `(2,4)` **flip**. This agrees with F₂, with `c₁ = 0`, and with the
    per-diagram corollary of §3.
* **Why the Euclidean labels do not carry over.**
  * The SD/ASD labels are eigenvalues of Euclidean momentum-space helicity, defined with `|k|`.
  * For timelike Lorentzian momenta `|k| → ±i√(−p²)`, with opposite branches on in- and out-legs.
  * In the bulk, Lorentzian normalisable modes are superpositions of regular and irregular modes of **opposite** duality (C1V `v2`,
    DSKETCH Lemma M). The remnant `Φ₃₄ = −b₊W₊[reg]` exists only through that reflection.
  * So a Euclidean label pair on `(2,4)` across an `s`-channel cut is not a Lorentzian state-helicity pair.
* **What this changes.** C1_COMPUTE §2's dictionary ("preserved line = one SD + one ASD Weyl tensor") is a statement about bulk Weyl
  types. It should not be read as a statement about state helicities across a cut.
* **What was not done.** This section is [D-sketch]. No explicit continuation of the four-point remnant was carried out. Nothing in
  §2–§3 uses it.

---

## 5. Route (B): what was computed

* **Done [DV] (in the `Ω_OPE` prescription).**
  * The direct action of the `k = 2` EE functional on the explicit energy-channel `T` exchange of C1_VALUE/DSKETCH: the `W₊` vertex
    at `(12T)`, the `(E − κCS)` vertex at `(T34)`, and the mirror arrangement.
  * Exactly, this is the pair of cross terms `I[E₃,W₁]`, `I[W₃,E₁]` of §2.5. These are evaluated with first-principles one-point
    matrices, which are independent of the HM quotation and reproduce it.
  * **Result: `c₁ = 0`.** Numerically the value is zero to the finite-difference floor `≤ 2·10⁻¹¹` relative to the Einstein
    integrand, at every sampled kinematic point.
  * The normalisations listed in the task do not enter a zero. Uncertainty: none beyond the grade of the inherited links in §3.
* **Independent kinematic routes.**
  * rest frame vs three boosted frames;
  * two independent theories (boson and fermion) combined through the structure decomposition;
  * eight random kinematic points × three light-ray structures;
  * Ward identities as normalisation checks.
* **Not done.**
  * A Lemma-R (Regge-limit, CHMRS App. E) evaluation of the energy-channel Witten diagrams themselves, independent of LIGHTSIDE_DIST's
    `Ω_OPE = Lemma R`.
  * A cross-check against JJ21's momentum-space three-point data. The free-field basis replaces it for the light-ray matrices.

  If the identity `Ω_OPE = Lemma R` on tree terms failed, `c₁` in H7's Lemma-R prescription would have to be recomputed as the arc
  of the remnant. It would then be the EE-projected top Regge coefficient of `b₊ × ⟨TTT⟩_{W₊}`. That is the only route by which a
  non-zero `c₁` could re-enter.

---

## 6. Checks (summary)

| check | result | grade |
|---|---|---|
| `f = v′−u′` affine in `(u₁,u₃)`, no `u₁u₃`, closed form | exact; finite-`r` → limit as `1/r` | [P] |
| EE `(2,4)` part at `x₀ → ∞` = `x₂₄²(z₂·z₄)` | exact (ratio 1) | [P] |
| HM calibration (free B/F): diagonal equal, flips opposite, `|flip|/diag = 1`, phase `e^{−4iφ}` | `≤ 2.4·10⁻¹⁵` | [DV] (reproduces Q8) |
| `M̂_flip = 0` (B, F); `M̂_diag,B = M̂_diag,F`; `M̂ = ½{Γ,Ĥ}` | `≤ 6·10⁻¹⁷`; `≤ 1.3·10⁻¹⁰`; `≤ 2.6·10⁻¹⁶` | [DV] (+ exact mechanism, §2.3(c)) |
| momentum-flux and dilatation Ward identities | `≤ 2·10⁻¹⁵`, `≤ 3·10⁻¹²` | [DV] |
| `O(t₄)` `T`-block term under EE, pointwise, each cross term | `≤ 1.75·10⁻¹¹` | [DV] |
| flip-sensitive control / Einstein / `W×W` | `O(1)` / non-zero / non-zero | [DV] |
| bounds at `c₁ = 0` from C1_COMPUTE §4 exact forms | reproduce 3.13206·10⁻⁶, 2.89126·10⁻⁸, 72.04059, 92.81351 | [DV] |

---

## 7. Effect on the headline and on existing documents (no file was edited)

1. **Result D (UNIFIED_THEOREM item 1, RESULT_D_STATEMENT).**
   * Replace `(1 + O(c₁Λ⁻²))` and the `R₀(x+√(1+x²))` form by `R₀`.
   * The bound is `3.13206·10⁻⁶` (300) / `2.89126·10⁻⁸` (1000), with a symmetric `t₄` interval.
   * Drop the "`∓0.22 % / ∓0.017 %` per unit `c₁`" error term (and DG_SCAN's `∓19941/(Λ²√(16L_A))`).
   * The H17 band (`±0.24 % / ±0.021 %`) is unrelated and stays.
2. **Result Q.**
   * Upper end `2L_AdS = 72.0406 / 92.8135` exactly, with no `+6.81·10⁻⁵c₁² / +5.44·10⁻⁷c₁²`.
   * NEARGAP2's "the `c₁²` coefficient has to be recomputed for ε = 0.040" is moot.
3. **Grades.**
   * UT3 / supplement row L3: "`c₁ = 0` given premise F₂ … F₂ [A]" → "**`c₁ = 0` [D]** (F₂ proved: affine `k = 2` weight [P];
     `E⁽¹⁾|_T = ½{E, −i n̄·D}` [D + DV]; EE identity contraction [P])".
   * Referee M5 and RU issue 14 are closed.
4. **Suggested letter sentence.**

   > "The only `O(t₄)` light term vanishes: at `x₀ = ∞` the `k = 2` subtraction factor is affine in the light-ray times, so the
   > stress-tensor exchange enters only through the energy flux and its first moment, both of which are helicity-diagonal for Einstein
   > gravity and helicity-flipping at `O(t₄)`, and the EE dual structure annihilates helicity flips."
5. **C1_VALUE §4.1–4.2, §8 and C1_COMPUTE §3.5–§3.6, §7.**
   * "value [I]" → "`c₁ = 0` [D]".
   * The natural-size estimate `|c₁| ≲ 0.1` [M] is superseded.
   * The "Regge-subcritical" reading of the remnant (C1_VALUE §4.2) is replaced by §4 here.

---

## 8. Where it could still fail

* **The prescription identity (i).**
  * If `Ω_OPE ≠ Lemma R` on the borderline (`n = k = 2`) tree terms of `⟨TTTT⟩`, the Lemma-R `c₁` would be the EE-projected arc of
    the remnant. That was not evaluated (§5).
  * The same identity carries the Einstein pole, so a failure would also affect H7.
* **(DZ) at the borderline.** If the `[TT]` double-trace sum failed to commute with the `1/c_T` expansion at `n = k`, double traces
  could contribute at `O(t₄)`.
  * This is inherited and shared with the Einstein term.
  * C1_VALUE step 4 and EE_FOUNDATIONS §1.4 use the same reading.
* **Parity-odd `⟨TTT⟩`.** For the `α_T` constant only, the diagonal-element argument (§2.3(d)) is [D] without a numerical check.
  It does not affect `c₁` (the `t₄` coefficient).
* **Dual structure.** The zero needs the `(2,4)` contraction to be helicity-preserving through first order. It holds for EE [P]. It
  fails for `q²`-type structures (control).

---

## 9. Reproduce

One core, `ulimit -v 3000000; export OMP_NUM_THREADS=1`. Every script prints `ALL CHECKS PASS: True` where it has checks; all did.

```
cd ads4/code/c1final
python3 -B c1f_affine.py     # §2.1, §2.4: f affine at x0 = inf; EE (2,4) -> identity (sympy, ~1 min) -> cache/c1f_affine.txt
python3 -B c1f_onepoint.py   # §2.3: free B/F one-point matrices of E, E1; HM calibration; Ward; M = {Gamma,H}/2 (~5 s) -> cache/c1f_onepoint.txt
python3 -B c1f_compose.py    # §2.5, §5: O(t4) T-block term pointwise under EE + controls (~12 s)       -> cache/c1f_compose.txt
python3 -B c1f_bounds.py     # §0, §7: headline bounds at c1 = 0 (imports ../c1/c1_quadratic.py read-only) -> cache/c1f_bounds.txt
```

`c1f_freefields.py` is the shared module. It contains the frames, the free-field wave functions, and the closed-form `Ĥ`, `M̂`.

---

## 10. Grade table

| statement | grade |
|---|---|
| `f⁺_{2;2} = v′−u′` at `x₀ = ∞` is affine in the retarded times, no `u₁u₃` term | **[P]** (sympy) |
| `k = 2` light-ray content = `E₃E₁`, `E₃E₁⁽¹⁾`, `E₃⁽¹⁾E₁`; `T`-block = `⟨Ψ₄|A₃P_TA₁|Ψ₂⟩` | [P] + [D] |
| `E⁽¹⁾|_T = ½{E|_T, −i n̄·D}` for every parity-even `⟨TTT⟩` (Einstein diagonal, `t₄` pure flip, `M̂_W = 0`) | **[D + DV]** (free-field basis exact; numerics to `10⁻¹⁶`) |
| no diagonal element for any non-Einstein structure (incl. `α_T`), in `E` and `E⁽¹⁾` | [D] (rotation + Ward) |
| EE `(2,4)` contraction = identity, diagonal through first order | **[P]** |
| **Premise F₂** | **[D + DV]** (was [A]) |
| `O(t₄)` `T`-block term of `Ψ_{2,ν}[EE]` = 0 pointwise, each cross term | [D]; [DV] to `1.8·10⁻¹¹` |
| **`c₁ = 0`** (whole `O(t₄)` light term of `Ĉ′_EE`, all `ν`) | **[D]** (inherits LIGHTSIDE_DIST (i), DZ (ii)) |
| `O(α_T)` constant = 0 | [D] |
| remnant's `s`-channel cut is a Lorentzian `(2,4)` flip; Euclidean SD/ASD labels ≠ state helicities across a cut | [D-sketch] (not load-bearing) |
| headline bounds without `c₁` caveat: D `3.13206·10⁻⁶ / 2.89126·10⁻⁸`, Q `72.0406 / 92.8135` | [D] (same as "(Pos) ⇒ (1), (2)") |
