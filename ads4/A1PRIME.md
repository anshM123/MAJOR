# A1′: replacing the full gap A1 by a higher-spin gap. Regge spins of light exchanges, drop-out per sum rule, and ABJM

*23 Nov 2026, A1′ agent. Unrefereed, not committed by this agent (another agent's commit 646acc6 swept in the `ads4/code/a1p/*.py`
files; their content equals the final versions). No existing file was edited. New code: `ads4/code/a1p/` (`a1p_flat.py`,
`a1p_flat2.py`, `a1p_pieces.py`, `a1p_run.py`); outputs in `ads4/code/a1p/cache/` (git-ignored by its own `.gitignore`). Read-only
imports: `resultq/lsd_core.py` (Setup, kernel K), `resultq/lsd_pos.py` (angle sets). Data read: `prl/data/nf02_psi.json`,
`prl/data/q3_full_delta0.2.json`. Literature read from arXiv LaTeX source (searched by number only): CHLPS 2201.06602 (§5
"Light spin-0 and spin-2 matter fields don't lower the cutoff", App. B), CPY 1804.00949 (§3, exchange content of ⟨SSSS⟩), BCP 1906.07195
(§1–§3). Compute: about 1.2 CPU-h, at most 2 processes, 1 thread, < 0.4 GB each.*

Grades: **[P]** proved; **[D]** derived; **[D-sketch]**; **[DV]** derived and checked; **[NE]** numerical evidence; **[M]** model;
**[Q]** quoted; **[A]** assumption; **[I]** open.

---

## 0. Summary

**The idea works, but not in the form proposed.** Every project functional sees only the helicity-preserving (MHV, `f(s,u)`) sector at
fixed `u`, i.e. the flat sum rules `B_k^{(1)}`. In that sector the Regge spins of light exchanges are not "spin of the exchange":
* a light scalar (either parity) has **n = 3**;
* a light massive spin-2 has **n = 4**;
* spin 1 and spin 3 never occur in `T × T`.

The crossed-channel Regge pole of a massive spin-2 (n = 2) exists only in the double-flip sector, which `EE` annihilates. So:

| sum rule | light J = 0 (n = 3) | light J = 2 (n = 4) | consequence |
|---|---|---|---|
| k = 2, 3 (Result D; the smeared part of Q) | contributes its block action, like a heavy state | same | no drop-out: needs a sign |
| k = 4 (Q) | **drops out** (n < k) | **borderline**: arc fixed by the top Regge coefficient (flat: `12|g₂|²p²/m⁴`) | g₄ value piece untouched; the derivative piece gains a favourable term |
| k = 5 (Q) | **drops out** | **drops out** | — |

**Signs (flat, exact polynomials) [P].**
* **Result D.** NF02 is positive on every light J = 0 and J = 2 state at every mass `0 < m ≤ Λ` (minimum 11.32 and 0.448 per unit
  coupling). Light states therefore enter with the sign of heavy states, and the D bound **holds and is strengthened**:
  `K R²Λ⁸/57600 + Σ_light F_D|g|²(…) ≤ 16L_A`. This is CHLPS's flat mechanism, transported to NF02.
* **Result Q.** The published `Ω^Q_up` is positive on light J = 2 at every mass (≥ 1.97). For light scalars:
  * thanks to the k = 4, 5 drop-out its action is bounded, `F_Q⁰ = −0.4304 + 2.0321 m²/Λ²`;
  * without the drop-out it would be `−u²` at small m;
  * it is still **negative for `m < 0.4602 Λ`**.
* **Repair.** `Ω^Q_up + ε*Ω^EE_θ` with `ε* = 0.038034` is positive on all light states. It holds under a **pure higher-spin gap**, with the
  upper end relaxed:
  * **72.04 → 123.08** (Δ_g = 300) and **92.81 → 160.99** (1000);
  * `c₁²` coefficient 1.26·10⁻⁴ / 1.01·10⁻⁶.

**Engine checks [DV]** (LIGHTSIDE's σ ⊗ T_free family, new Witten pieces, effective spin `n_eff = J₀(σ) − 2`).
* **Drop-out at n < k**, single-trace block included term-wise. Extrapolated symmetrised sums are ≤ 1.3·10⁻⁵ of scale (most ≤ 5·10⁻⁶) for:
  * s + t exchange pairs with degree-6 and degree-8 residues;
  * a degree-6 u-exchange;
  * a degree-6 contact;
  * at k = 3, 4, 5.
* **Borderline at the Result-Q borderline k = 4** (new: the R⁴ analogue, by three routes). The value equals that of the top-Regge term to:
  * 5·10⁻⁷ (contact);
  * 8·10⁻⁸ (crossed exchange);
  * 1.4·10⁻⁶ (direct-channel exchange with its single-trace block).
* **Borderline at k = 2.** LIGHTSIDE's xu4 is reproduced; a new direct-channel case agrees to 1.1·10⁻⁶.

**A1′ statements (§6).**
* **A1′_D = A1′_Q = higher-spin gap:** every single trace in `T × T` with J ≥ 4 has twist ≥ Δ_g. Light J ≤ 2 single traces are allowed at
  any twist.
* **D** needs no coupling condition. **Q** needs none with the ε-repair. The published Q constants need one coupling condition: no
  light scalar below 0.46Λ, or a small summed coupling.
* Net grades: implication [D-sketch] (one new [D-sketch] link: "light value at k = 2, 3 = block action"); flat signs [P]; AdS light-state
  images [D-sketch] at leading order, finite-Δ corrections [I].

**ABJM (§7).**
* In the type-IIA 't Hooft regime the higher-spin gap holds: every J ≥ 4 single trace is a string state with twist ≈ ML.
* Through the orders computed (supergravity + R⁴), **`T × T` contains no light single traces at all** at O(1/c_T) (CPY, BCP [Q]; KK
  orthogonality [Q CHLPS]). The "light R-singlet in `T × T`" of UT §1/§9 and referee M13 is not supported at leading order: the cited
  breathing mode is, if it has Δ = 6 = 2Δ_T as I recall (not re-derived), an even-twist, double-trace-degenerate mode with no exchange
  pole (CPY's (B,+) class).
* So ABJM is inside A1′ and A2 at leading order in 1/λ, and Λ ↔ ML becomes the spin-4 gap of CHLPS/CHMRS.
* What is bounded is the 4d MHV coefficient of the zero-mode gravitons, which equals the 10d IIA R⁴ coefficient at zero internal momentum
  (no KK exchange between zero modes at tree level): `g₄M⁶/8πG₄ = 2ζ(3)`.
* The 3.3 % / 2.6 % placement becomes an **in-scope consistency check**. It is not a sharp test: (Pos) is not established for string
  states, and the bound is far from saturation.

---

## 1. Setting and what "Regge spin" means here

* **The Regge point.** All project functionals are AdS uplifts of the flat fixed-`u` sum rules `B_k^{(1)}` of CHLPS (2201.06602 (2.26)),
  dispersing in `s` at fixed `u = −(p₁+p₃)²`. The `u`-channel is the project's (13) channel (fixed Mellin variable δ₁₃).
* **The sector.** They act only on the MHV function `f(s,u)` of `M(1⁺2⁻3⁻4⁺) = ⟨23⟩⁴[14]⁴f`:
  * `Ω^EE`: flip configurations are annihilated (C1_COMPUTE §3.2: ≤ 4·10⁻¹⁴, ≤ 5·10⁻⁴⁶ for the double flip);
  * Ĉ₄, Ĉ₅ are flip-free (F_NORM §2.4).
* **The Regge spin of a tree term X** is `n_X := 4 + (growth exponent of f_X in s at fixed u)`. Equivalently, it is the growth of the full
  amplitude, since the helicity factor `|⟨23⟩⁴[14]⁴| = t⁴ ~ s⁴`. Checks: Einstein `8πG/(stu)` → 2; R⁴ (`g₄`) → 4; `g₅t` → 5.
* **AdS.** The same count holds for the top Mellin degree at fixed δ₁₃ (Penedones [Q]; spinning [D-sketch]).
* **In LIGHTSIDE's product family** the relevant spin is `n_eff = J₀(σ) − 2` (LIGHTSIDE §2). Every engine check below is stated in
  `n_eff`.
* **The light-side rule (LIGHTSIDE L4) on a tree term:** `Ω_k[X] = D_LR[X] + 𝒜_k[X]`, with:
  * `𝒜_k = 0` if `n_X < k`;
  * at `n_X = k`, a linear functional of the top Regge coefficient.

---

## 2. Classification of light exchanges (goal 1)

### 2.1 Which single traces can appear in `T × T` [Q + D]
* **Spin 1 and spin 3: never.** For identical conserved T in d = 3 there is no `⟨TTO_J⟩` structure at J = 1, 3 (PARITY_ODD_ROUTEA §1,
  quoted in HYPOTHESES_CLOSE §3.1). The flat counterpart is CHLPS §2: "Angular momentum conservation forbids the decay of a state of
  half-integer or odd spin to a pair of gravitons". So R-currents (and flavour currents) never enter `⟨TTTT⟩` at O(1/c_T).
  (This concerns the `⟨TTJ⟩` structure, not any current coefficient; c_T is untouched.)
* **J = 0 and J = 2:** one parity-even structure (`++`) and one parity-odd structure each. **No `+−` structure exists for J < 4**: at even J
  the space is (`++`, `+−`) only from J = 4 on (HYPOTHESES_CLOSE §3.1; flat: `c^{+−}` only for J ≥ 4, CHLPS (2.16), "a particle of spin <
  4 cannot decay to a (+−) two-graviton state", App. B).
* **J ≥ 4 even:** `++`, `+−` and one odd structure. **Odd J ≥ 5:** one structure, `+−`, parity-odd.

### 2.2 Regge spin per channel in the sector the functionals see [D; flat Q for J = 0, 2]
Channel labels are those of `f` at fixed `u`:
* `s = (12)` and `u = (13)` pair opposite helicities, so only `+−` couplings (J ≥ 4) produce poles there;
* `t = (14)` pairs like helicities (`++`).

The vertex structure enters through the on-shell three-point amplitude. For a massive J < 4 state it is unique per helicity sector
(AHH-type `[14]^{4−J}(⟨1X⟩⟨4X⟩)^J`), so different Lagrangian vertices differ only by contact terms. The "Regge spin of an exchange" is
therefore the minimal growth over contact completions, and the bounded coefficients are defined relative to that completion.

| light single trace in `T × T` | bulk cubic coupling (on-shell) | s = (12) | t = (14) | u = (13) | **n** (MHV sector) | flip sectors (invisible) |
|---|---|---|---|---|---|---|
| graviton T, Einstein structure | EH, 2∂ | pole | pole | pole | **2** (s, t, u poles combine; each alone would be 3) | — |
| T with `W³` (t₄, α_T) on both vertices | `W³`, 6∂ | — | pole `su/t` | — | **4** (the `2πG|ĝ₃|²su/t` term) | E×W³: n = 2 [D-sketch, C1] |
| scalar, parity-even | `φW²` (4∂); 2∂ couplings vanish on shell | — | `|g₀|²/(m²−t)` | — | **3** | all-plus `h` |
| scalar, parity-odd | `φWW̃` (4∂) | — | same, `|g₀|²` | — | **3** | `h` |
| spin 1, spin 3 | none (no `⟨TTJ⟩`) | — | — | — | absent | absent |
| massive spin 2 (either parity) | `X·WW`-type, unique `++` amplitude | — | `|g₂|²(1 − 6su/m⁴)/(m²−t)` | none (no `+−` for J < 4) | **4** (CHLPS: "spin-4 at fixed s") | u-channel double flip, n = 2: killed by EE (≤ 5·10⁻⁴⁶), absent from C4/C5 |
| J ≥ 4 even, `++` | — | — | residue a polynomial of degree J/2 in `su` | — | **3 + J/2** | |
| J ≥ 4, `+−` | — | pole | — | `d̃^J_{44}(1+2s/m²)/(m²−u)` | **J** | |
| contacts | `R⁴` (g₄), `D²R⁴` (g₅), … | | | | 4, 5, … | all-plus `ĝ₄`: invisible (F_NORM §2.4) |

**The rule [D].**
* **Crossed (`u`) channel:** n = J. This needs a `+−` coupling, so for massive states J ≥ 4.
* **Direct like-helicity (`t`) channel:** `n = 3 + J/2`. The residue on the pole is `d̃^J_{00} = P_J(1+2u/m²)`. Bose symmetry
  `f(s,u) = f(u,s)` forces it to be written as a polynomial in `su` (degree J/2) and `t`. At fixed `u`, `su ~ s`, so `f ~ s^{J/2−1}`.
  The `4` is the helicity factor.
* **Contact ambiguities.** For J = 0 they are contacts with n ≥ 4; for J = 2, contacts that can "shift the g₄ and g₅ coefficients by
  multiples of |g₂|²" (CHLPS App. B, verbatim).
* **The minimal completion** (CHLPS eq. KK MHV amp) is the one used throughout. With it, **g₄ is the R⁴ coefficient with the light fields
  kept in the amplitude.** Integrating out a light scalar would shift it by `|g₀|²/m²`.

**On the task's example.**
* The `φW²` (χ) exchange has n = 3 (F_NORM §2.3), as above.
* A "two-derivative" scalar coupling to two gravitons has a vanishing on-shell 3-point amplitude: in Einstein frame it is a field
  redefinition, or proportional to `V′(φ₀) = 0`. Its "exchange" is a pure contact, so there is no separate two-derivative case.
* For KK modes, the two-derivative couplings to two zero-mode gravitons vanish by orthogonality (CHLPS §5: "proportional to an overlap
  integral ∫√g χ(y) = 0").

**AdS.**
* The same table holds for the Mellin amplitude's top degree at fixed δ₁₃ [D-sketch: top degree = flat amplitude, Penedones; F_NORM
  A3 checks this for a scalar exchange to 10⁻⁶].
* The helicity selection (no `+−` at J < 4) is an exact CFT statement [Q].

---

## 3. Drop-out per sum rule (goal 2)

| exchange (n) | k = 2 | k = 3 | k = 4 | k = 5 |
|---|---|---|---|---|
| T, Einstein (2) | borderline: H7 pole `16πG/(ν²+9/4)` [D] | Ĉ₃: 0 [D] | drops (Lemma A) [D] | drops [D] |
| T, W³⊗W³ (4) | block action → `−K|ĝ₃|²` term | → `|ĝ₃|²` term | arc: the `4πG|ĝ₃|²p²` slot (κ) | drops |
| light scalar, P-even/odd (3) | **block action** (no kernel poles) | **block action** (n = k, but the unsubtracted k = 3 value is again the block's own) | **drops out** [D; DV toy] | **drops out** [D; DV toy] |
| light massive spin 2 (4) | **block action** | **block action** | **borderline**: arc ∝ top coefficient; flat `−B₄ = 12|g₂|²p²/m⁴` [Q CHLPS]: 0 in the value piece, `+12|g₂|²/m⁴` in `∂_{p²}` | **drops out** [D] |
| light J ≥ 4 (3+J/2, or J) | block action, sign uncontrolled | same | n ≥ 5 > k: **breaks** (contaminates g₄, g₅) unless J = 4 `+−` | breaks for n ≥ 6 |

**Why k = 2, 3 never give a drop-out, and what "n > k breaks" means there [D-sketch].**
* The k = 2, 3 CLS-D functionals are unsubtracted. They put `[TT]` at double zeros and drop contacts (L0 [D]).
* In flat space: the kernels `(2s − p²)` and `1` have no poles, so polynomial parts of `f` contribute nothing. A light pole contributes
  exactly its residue, which is the heavy-state image continued to `m < M` (CHLPS: `B₂|_matter = −|g₀|²(p² − 2m²) − …`; reproduced in
  §4.1).
* In AdS: the light exchange's double-trace part sits at the double zeros. For massive J < 4 there is no preserved⊗preserved `u`-channel
  Regge pole (§2.1). So the light value is the block action `Tr[P_L M_Ω(O_L)]`.
* So at k = 2, 3 a light state is neither zero nor fixed by its Regge data. It acts as one more "heavy" state, and the sum rule needs the
  sign of the functional on it (§4).
* For subtracted k ≥ 4 with n > k, the light value depends on the contact completion through n > k contacts (g₅, g₆, …). That is a true
  breakdown.
* Grade of "light value at k = 2, 3 = block action": [D-sketch]. Flat [P]; AdS uses L0 on the exchange's double traces plus the
  helicity selection.

**At k = 4, 5.**
* Drop-out is LIGHTSIDE's Lemma A applied to the light Witten exchange `W_L`, all channels, single-trace block included: block +
  double traces = 𝒜_k = 0 for n < k. This is the same input F_NORM uses for its J₀ = 3 χ exchange.
* Grade [D]; checked on new pieces [DV] (§5).

---

## 4. Signs: Result D with light J = 0, 2 states, and Result Q (goal 4)

### 4.1 Method and validation (`a1p_flat.py`) [DV]
* Light contributions to the low side are taken from CHLPS (5.x)/(B.1) in the minimal completion (units M = 1 ↔ Λ). With
  `P := m⁴ − 6m²p² + 6p⁴`:
  * `B₂|_m = −|g₀|²(p²−2m²) − (|g₂|²/m⁴)(p²−2m²)P`;
  * `B₃|_m = |g₀|² + (|g₂|²/m⁴)P`;
  * `B₄|_m = −12(|g₂|²/m⁴)p²`;
  * `B₅|_m = 0`.
  (CHLPS print `m²` for `m⁴` in P; a typo, as FLAT_GRAVITON_BOUND §1.3 already noted.)
* A functional is **favourable on light matter iff F[matter] ≥ 0** (CHLPS convention).
* **Validation.** Applied to CHLPS's own `F_{g3}` (ψ₂ = (1−p)³p(65+155p+47p²), ψ₃ = 2(1−p)³(81−1927p+4601p²), −5∂B₄), the code returns
  exactly their printed `2591/210 m² + 1/18` and `2591/210 m⁶ − 239/18 m⁴ − 64573/1540 m² + 330151/4004` (exact rationals).

### 4.2 Result D (NF02) [P, flat]
Per `|g₀|²` and per `|g₂|²/m⁴`, as exact polynomials in `m²` (NF02 normalised to `k(0) = 1`):

| | polynomial | min on `0 < m ≤ Λ` | `m → 0` |
|---|---|---|---|
| `F_D⁰` (J = 0, either parity) | `11.317027 + 20.771893 m²` | 11.317 (m → 0) | `μ₃(0) − μ₂(2)` |
| `F_D²` (J = 2, either parity) | `3 − 6.766265 m² − 6.905463 m⁴ + 20.771893 m⁶` | **0.448298** at `m = 0.675Λ` | `3k(0) = 3` |

* **No real root in `x = m² ∈ [0,1]`** (sympy root isolation on the exact rational polynomials): both are sign-definite [P].
* The `m → 0` limit of the spin-2 action is exactly `3 × k(0)`. A massless-limit `++` spin-2 is the `W³`-type exchange, so it enters with
  the sign of the `−K|ĝ₃|²` term that the functional was normalised to make negative. This is the structural reason D is favourable
  (CHLPS: "The spin-2 contribution with m → 0 is proportional to that from the cubic coupling |ĝ₃|²").
* **Hence (flat, leading order) the D sum rule with light states reads**

      K·R²Λ⁸/57600 + Σ_{light ℓ} F_D^{(J_ℓ)}(m_ℓ/Λ)·w_ℓ  ≤  16 L_A + 16 I_D Λ² c₁ t₄,
      w_ℓ = |g₀(ℓ)|²Λ⁴/(πG)  (J = 0),   |g₂(ℓ)|²Λ⁸/(πG m_ℓ⁴)  (J = 2),

  with every `F_D ≥ 0.448`. **The D bound holds with unchanged constants and is strengthened. It also bounds the summed squared couplings
  of all light J ≤ 2 single traces to TT** (the AdS form of CHLPS's "(g3 with KK)").
* **AdS form [D-sketch].**
  * The statement needed is (Pos_Ω^EE) for light J ≤ 2 states.
  * For `Δ_L ≪ Λ` the functional samples `ν = Λp ≫ Δ_L`, where the exchange looks massless. Its image is the `m → 0` flat value times
    `(1 + O(Δ_L²/Λ²))`.
  * For `Δ_L` comparable to `Δ_g` (still below the gap) it is the heavy-type AdS image near the gap, which the project computes only
    above the gap [I].
* **On the task's premise ("light spin-2 at k = 2 is borderline, n = 2").** That n = 2 is the crossed-channel Regge pole. For a massive
  spin-2 it is pure helicity flip on both lines (no `+−` structure at J = 2), and EE annihilates double flips to ≤ 5·10⁻⁴⁶ (C1 §3.2).
  What EE sees is the direct-channel `++` pole, n = 4 > 2. That contributes its block action, which is positive [P flat]. The sign is
  definite and favourable: coupling² × a positive kernel value, as the task hoped, but for a different reason.

### 4.3 Result Q (`Q3-full-delta0.2`) [P, flat] (`a1p_flat2.py`)
* **The flat light action.** `F_Q = ∫ψ₂B₂ + ∫ψ₃B₃ + a(−∂_{p²}B₄|₀ + B₅|₀) − cB₄|₀` on matter. Equivalently: the k = 2, 3 block actions,
  plus `+12a|g₂|²/m⁴` from the k = 4 spin-2 arc, plus 0 from every drop-out.
* **Check:** `K_Q := −2μ₂(6) + 2μ₃(4) + 4a = 0.879266` and `K_Q/8 = 0.109908 = κ` (file value).

| | polynomial in `m²` | sign on `0 < m ≤ Λ` |
|---|---|---|
| `F_Q²` (J = 2) | `2.637798 + 0.480850 m² − 3.156162 m⁴ + 2.032073 m⁶` (`m → 0`: `3K_Q = 24κ`) | **positive**, min 1.9725 [P] |
| `F_Q⁰` (J = 0) | `−0.430430 + 2.032073 m²` | **negative for `m < 0.460237 Λ`** [P] |

* **The drop-out is what makes this mild.** Without the k = 4, 5 drop-out, the light scalar's action would be the heavy formula continued,
  `FF(u) = 2.032 − 0.430u − u²` (`u = Λ²/m²`). That is negative for every `m < 0.903Λ` and diverges like `−Λ⁴/m⁴`. With the drop-out it is
  bounded (≥ −0.430 per unit coupling) and negative only below 0.46Λ.
* **Spin 2:** favourable, again because its `m → 0` limit is the `|ĝ₃|²` structure with coefficient `+24κ`. So light spin-2 states need
  exactly the sign input κ > 0 that Q already makes (L9).
* **Repair (pure higher-spin gap, no coupling condition).**
  * `F_Q + εF_D ≥ 0` on all light J = 0, 2 states for `ε ≥ ε* = 0.4304302/11.3170270 = 0.0380339`. The worst case is `m → 0`; the combined
    J = 2 action stays ≥ 2.185.
  * Both functionals are positive on heavy states by (Pos). So **`Ω^Q_up + ε*Ω^EE_θ`** gives:

        g̃₄Λ⁶/8πG ≤ 2L_AdS + 2ε*L_A − (κ + ε*/8)R²Λ⁸/57600 + 2(I_Q + ε*I_D)c₁t₄Λ²

  * The supremum over (t₄, α_T) gives:

    | Δ_g | `L_A` (NF02, reproduces `R₀`) | published `2L_AdS` | **repaired upper end** | `c₁²` coefficient |
    |---|---|---|---|---|
    | 300 | 670.913 | 72.041 | **123.076** (×1.708) | 1.2598·10⁻⁴ (was 6.8145·10⁻⁵) |
    | 1000 | 896.224 | 92.814 | **160.988** (×1.735) | 1.0062·10⁻⁶ (was 5.4428·10⁻⁷) |

  * Cost: Q now also needs (Pos_Ω^EE), i.e. item (2) is no longer independent of item (1)'s positivity hypothesis.
* **Published constants with a coupling condition.** With `Ω^Q_up` alone:

      g̃₄Λ⁶/8πG ≤ 2L_AdS + Σ_{light scalars ℓ} max(0, 0.430430 − 2.032073 m_ℓ²/Λ²)·|g₀(ℓ)|²Λ⁴/(8πG) + (κ, c₁ terms as before).

  The correction vanishes if no scalar lighter than 0.4602Λ couples to TT. The D sum rule bounds the correction by the repair cost
  (`2ε*L_A`).
* **Lower end `Ω^Q_low`.** Light J = 0 drops out [D]. The light J = 2 arc is `∝ p²`, which is 0 in the flat value piece. In AdS it is a
  degree-1 polynomial in ν², i.e. an `O(Λ⁻²)·|g₂|²/m⁴` term that joins the existing "−O(Λ⁻²)·(EFT data)" [D-sketch]. `g̃₄ = g₄` (bare,
  minimal completion) + O(Λ⁻²).

---

## 5. Engine checks (goal 3) (`a1p_pieces.py`, `a1p_run.py`) [DV]

**Setting.**
* LIGHTSIDE's term-wise s-channel functional on `T′ = σ ⊗ T_free`: production grid (s₀ = 10⁻⁹, 16-point panels, ratio 1.8, 1248
  nodes); angle set P1 (P2, P3 for the k = 4 borderline, §5.3); δ₁₃ step 0.1, |Im δ₁₃| ≤ 9; residues n ≤ 150.
* "sym" = E(a,b) + (−1)^k E(b,a), summed over the six free Wick classes.
* New σ-pieces are built exactly in sympy (M and ∂M/∂δ₁₂; checked against finite differences to 10⁻¹¹).
* The s-channel single-trace pole (δ₁₂ = 1, i.e. the light block) is included term-wise.

| piece | Mellin amplitude | J₀(σ) | **n_eff** | models |
|---|---|---|---|---|
| X3 = xs3 + xt3 | `Q_s³/(δ₁₂−1) + Q_t³/(δ₁₄−1)` (top `y⁵` cancels in the pair; top `y⁴` coefficient `(19δ₁₃−12)/64`) | 4 | **2** | direct-channel exchange with a high-degree residue, block included |
| X4 = xs4 + xt4 | `Q⁴` residues; top `y⁶` coefficient `(25δ₁₃−16)/256` | 6 | **4** | same, at the k = 4 borderline |
| xu6 | `Q_u³/(δ₁₃−1)` | 6 | **4** | crossed-channel exchange at the k = 4 borderline |
| c6 | `σ₂³` | 6 | **4** | the R⁴ analogue of k = 4 in this family (c4 is the k = 2 one) |
| xu4 | `Q_u²/(δ₁₃−1)` | 4 | 2 | LIGHTSIDE's borderline exchange (reproduction) |
| rX3, rX4, r6 = D⁶/8, r6u = D⁶/(64(δ₁₃−1)), r4u | top-Regge representatives, `D = δ₁₂ − δ₁₄` | | | |

### 5.1 Drop-out, n_eff < k (extrapolated sym; scale = largest single-ordering or single-trace value)

| k | piece (n_eff) | sym N = 50 / 100 / 150 | extrapolated | / scale |
|---|---|---|---|---|
| 3 | X3 (2) | −2.7·10⁻⁹ / −6.0·10⁻¹⁰ / −2.6·10⁻¹⁰ | −2.2·10⁻¹¹ | 1·10⁻⁶ (scale 2.3·10⁻⁵; block residue alone 2.3·10⁻⁵) |
| 4 | X3 (2) | +1.0·10⁻⁸ / +2.1·10⁻⁹ / +7.3·10⁻¹⁰ | +5.4·10⁻¹¹ … +2.3·10⁻¹⁰ | ≤ 4·10⁻⁶ (scale 6.6·10⁻⁵) |
| 5 | X3 (2) | −1.1·10⁻¹¹ / −5.9·10⁻¹³ / +2.6·10⁻¹⁴ | ~10⁻¹⁴ | 1·10⁻⁹ |
| 5 | X4 (4) | −9.4·10⁻⁹ / −2.0·10⁻⁹ / +2.0·10⁻¹⁰ | +8·10⁻¹² … −1.7·10⁻¹⁰ | ≤ 5·10⁻⁶ |
| 5 | c6 (4) | −6.4·10⁻⁷ / −1.4·10⁻⁷ / −3.8·10⁻⁸ | −2.4·10⁻⁹ … −1.3·10⁻⁸ | ≤ 1.3·10⁻⁵ (scale 1.0·10⁻³) |
| 5 | xu6 (4) | −6.7·10⁻⁸ / −1.5·10⁻⁸ / −5.2·10⁻⁹ | −3.7·10⁻¹⁰ … −1.4·10⁻⁹ | ≤ 1.3·10⁻⁵ |

* The extrapolation range brackets two methods: a 3-parameter fit `S∞ + CN^{−q}` on N = 60…150, and `S(100)` + a power tail fitted on
  n = 40…100.
* In every case the single-ordering values and the **single-trace block residue** are O(scale). For example, X3 at k = 4 has
  `st(ab) = −6.6·10⁻⁵`, larger than the total. So the zero is a cancellation between the light block and its own double traces, which is
  the content of the drop-out.
* The same pieces at their borderline give |sym|/scale ≈ 1.4 (§5.2). The drop-out is 5–9 orders below that.

### 5.2 Borderline, n_eff = k: value = value of the top Regge term [(P), DV]

| k | piece vs top-Regge term | sym (N = 150) | rel. diff | extrapolated value | (X − r) partial sums N = 50/100/150 |
|---|---|---|---|---|---|
| 2 | xu4 vs r4u (LIGHTSIDE) | −3.41130060·10⁻³ vs −3.41130088·10⁻³ | 8·10⁻⁸ | −3.4133·10⁻³ | +2.8·10⁻⁹ / +6.5·10⁻¹⁰ / +2.8·10⁻¹⁰ |
| 2 | **X3 vs rX3** (new; block included) | −4.15133034·10⁻⁴ vs −4.15132557·10⁻⁴ | 1.1·10⁻⁶ | −4.1536·10⁻⁴ | −5.7·10⁻⁹ / −1.2·10⁻⁹ / −4.8·10⁻¹⁰ |
| **4** | **c6 vs r6** (contact, R⁴ analogue) | −1.01360095·10⁻¹ vs −1.01360147·10⁻¹ | 5.2·10⁻⁷ | −1.0138·10⁻¹ | +8.5·10⁻⁷ / +1.8·10⁻⁷ / +5.3·10⁻⁸ |
| **4** | **xu6 vs r6u** (crossed exchange) | −1.41993270·10⁻² vs −1.41993282·10⁻² | 8.2·10⁻⁸ | −1.4204·10⁻² | +1.3·10⁻⁸ / +2.8·10⁻⁹ / +1.2·10⁻⁹ |
| **4** | **X4 vs rX4** (direct exchange, block included) | −1.13597403·10⁻³ vs −1.13597244·10⁻³ | 1.4·10⁻⁶ | −1.1356·10⁻³ | −2.0·10⁻⁸ / −4.2·10⁻⁹ / −1.6·10⁻⁹ |

* **Convergence.** Borderline residues decay with p = 2.2–2.3 (cf. LIGHTSIDE's p ≈ 2), so the values are finite.
* **Grid noise.** For degree-6 pieces the residues are clean up to n ≈ 110 and grid-noise-dominated beyond (sign flips at n = 120–150). The
  noise is identical in X and its top-Regge term (the differences keep decaying like N⁻²). So the extrapolated values use n ≤ 100 plus
  a fitted tail (uncertainty ~2·10⁻⁴ relative), while the (P) comparison is unaffected.
* **xu4.** The value here (−3.41130·10⁻³ at N = 150) differs from LIGHTSIDE's −3.41228·10⁻³ by 3·10⁻⁴. That is the δ₁₃ quadrature:
  LIGHTSIDE used h = 0.2, whose trapezoid error near the u-pole they estimate at ~10⁻⁴; here h = 0.1. The (P) agreement is the same.
* **New.** These are the first engine checks of the borderline mechanism **at the Result-Q borderline k = 4**, for all three kinds of tree
  term (contact, crossed exchange, direct exchange with its single-trace block).

### 5.3 Robustness at k = 4 (P2, P3)
Same pieces at P2 (Δφ = 4.1) and P3 (Δφ = 0.3, near the collinear set), k = 4.

| tag | c6 vs r6 | xu6 vs r6u | X4 vs rX4 | value (clean-range extrapolation, c6 / xu6 / X4) |
|---|---|---|---|---|
| P1 | 5.2·10⁻⁷ | 8.2·10⁻⁸ | 1.4·10⁻⁶ | −1.0138·10⁻¹ / −1.4204·10⁻² / −1.1356·10⁻³ |
| P2 | 1.1·10⁻⁶ | 8.2·10⁻⁸ | 4.8·10⁻⁶ | −2.2747·10⁻³ / −2.4563·10⁻⁴ / −3.5454·10⁻⁵ |
| P3 | 4.6·10⁻⁷ (N = 150), 8·10⁻⁷ (N = 100) | 3.7·10⁻⁸ | 9.9·10⁻⁵ at N = 150; **8.8·10⁻⁶ at N = 100** | −1.1154·10² / −5.2123·10¹ / not converged |

* At P3 the degree-8 X4 residues are grid-noise dominated already beyond n ≈ 100: the partial sum jumps from −0.0924 to −0.1029 between
  N = 100 and 150, identically in X4 and rX4.
* So (P) holds there at N = 100, but the absolute P3 value of X4 is not converged on this grid [NE, grid-limited].
* c6 and xu6 converge at P3 (residue decay p = 2.1–2.2).
* The borderline values grow towards the collinear set (P3 ≫ P1 ≫ P2), as LIGHTSIDE found for c4 at k = 2. Their integrability against
  C_ν at a = b is the same open item [I] (EESTEP §4: not separately load-bearing).

### 5.4 What the toy does and does not test
* **It tests** the bookkeeping that Q's drop-out uses: n < k ⇒ zero, with the light single-trace block included; n = k ⇒ fixed by the
  top Regge coefficient. It does so at k = 2–5, for direct and crossed exchanges and contacts.
* **It does not test signs.** The σ-residues are arbitrary microcausal polynomials, so the sign of a toy arc carries no physics. The
  physical signs come from §4 (flat helicity structure; exact polynomials).
* **The product shift `n_eff = J₀(σ) − 2`** (LIGHTSIDE §2) is respected throughout. A σ-piece of Mellin degree 6 is the k = 4 borderline,
  and σ-exchanges of spin 2 (`xs`, `xt`, `xu`; n_eff ≤ 0) drop out at every k ≥ 2 (LIGHTSIDE §3.2, §3.4, reconfirmed by the k = 2–5
  runs there).

---

## 6. A1′: statements and grades (goal 5)

**A1′ (higher-spin gap) [A].** Every single-trace primary in `T × T` with spin `J ≥ 4` has twist `Δ − J ≥ Δ_g`. Single traces of spin 0
and 2 (either parity) may have any dimension. Spins 1 and 3 do not occur in `T × T` (§2.1). This is the CHLPS/CEMZ notion: "States
with spin greater than two … are assumed to have mass above the cutoff". With it, `Λ := √((Δ_g−2)(Δ_g−1))` is the spin-4 gap.

**A1′_D (Result D, item (1)).**
* Statement: A1′, plus (Pos_Ω^EE) for **every** J ≤ 2 single trace in `T × T` regardless of twist.
* Conclusion: item (1) holds with **unchanged constants**, strengthened by `+Σ_light F_D w_ℓ` (§4.2).
* No coupling condition.

| link | grade |
|---|---|
| light value at k = 2, 3 = block action (no drop-out; no u-channel preserved Regge pole for J < 4) | [D-sketch] (flat [P]) |
| flat sign of every light J ≤ 2 action (NF02) | **[P]** (exact polynomials, no root on `[0, Λ]`) |
| AdS light-state action = flat value × (1 + O(Δ_L²/Λ²)) for `Δ_L ≪ Λ` | [D-sketch] |
| AdS action of J ≤ 2 states with twist just below Δ_g | [I] (not computed; the project's near-gap computations start at the gap) |
| net implication "A1′_D + (Pos) ⇒ (1)" | **[D-sketch]**, as before (EESTEP, H17 band etc. unchanged) |

**A1′_Q (Result Q, item (2)).** Two versions.
* **(Q-a) Pure higher-spin gap.** A1′, plus (Pos) for `Ω^Q_up`, `Ω^Q_low` and `Ω^EE_θ` on all states of `T × T`.
  * Upper end `Ω^Q_up + 0.038034 Ω^EE_θ`: **123.08 + 1.26·10⁻⁴c₁² (300), 160.99 + 1.01·10⁻⁶c₁² (1000)**.
  * Lower end unchanged in form.
  * No coupling condition.
* **(Q-b) Published constants.** A1′, plus the coupling condition
  `Σ_{light scalars} max(0, 0.4304 − 2.0321m²/Λ²)|g₀|²Λ⁴/8πG ≤ δ`. Upper end `2L_AdS + δ`; `δ = 0` if no scalar lighter than `0.4602Λ`
  couples to TT. Light spin-2 states need nothing (favourable, given κ > 0).

| link | grade |
|---|---|
| drop-out of light J = 0 at k = 4, 5 and of light J = 2 at k = 5 (Lemma A on the light exchange, block included) | **[D]**; toy [DV] (§5.1) |
| light J = 2 at k = 4: arc ∝ top coefficient | structure [D + DV] (§5.2); flat value `12|g₂|²p²/m⁴` [Q]; AdS normalisation [D-sketch] (same route as F_NORM) |
| k = 2, 3 parts of `Ω^Q_up` on light states = block action | [D-sketch] |
| flat signs: J = 2 favourable, J = 0 negative below 0.4602Λ; ε* = 0.038034 | **[P]** |
| meaning of g₄: bare coefficient, minimal light-exchange completion | [D] (CHLPS App. B convention) |
| net implication | **unchanged from UT §3** ([I] in the a = b link, H_EFT(b) [A] for (2′), κ > 0 if R ≠ 0), plus the [D-sketch] link above |

**Does Result Q hold under a pure higher-spin gap?** Yes, as (Q-a), with constants weaker by ×1.71 / ×1.73. The published 72.04 / 92.81
need the mild condition (Q-b). A dedicated Q functional, re-optimised with light-matter positivity imposed as CHLPS do (their
SDPB constraints on `m_ℓ ∈ (0, M)`), would remove most of the cost [I: not done].

**The KK tower [D].**
* **Finitely many terms.** At fixed Δ_g, a compactification spectrum (discrete, finite multiplicities, Δ_ℓ growing with the level) has
  finitely many single traces below Δ_g. So "the light sum" is a finite sum; the number of modes grows like Δ_g^{d_int}, but finitely.
* **No growth per mode.** Every light mode has J ≤ 2, so n ≤ 4 (bulk KK fields have spin ≤ 2). The per-mode Regge spin does not grow
  with the level.
* **Sums of favourable terms are bounded automatically.** Under D, and Q-a, each term is ≥ 0 and the total is bounded by the light side
  (monotone). Convergence is a consequence, not an extra hypothesis.
* **H_K** covers light blocks like any other OPE term (Rindler positivity per block).
* **Not covered [I].** A spectrum that accumulates below Δ_g would need the per-mode drop-out to be summable (dominated convergence for
  Σ_ℓ 𝒜₄[W_ℓ]). That is not the KK case.

---

## 7. ABJM (goal 6)

**(i) The higher-spin gap holds** (IIA 't Hooft regime, λ ≫ 1, c_T → ∞) [D-sketch].
* The spin ≥ 4 single traces are string states. The lightest is at the first massive level, `M² = 4/α′`, whose highest spin is 4. Its
  twist is ≈ ML − 4 with `ML = 2^{5/4}π^{1/2}λ^{1/4}` at leading order.
* Supergravity KK modes on CP³ have spin ≤ 2.
* Large-spin string states have larger twist (short strings `Δ ~ (ML)√J`, long strings `~ log J`), so the minimum twist over J ≥ 4 is at
  the first level.
* So **Λ = the spin-4 gap = CHLPS/CHMRS's M** at leading order in 1/(ML). This is the M in `g₄M⁶/8πG = 2ζ(3)`, `M² = 4/α′`. Referee M13
  problem 3 (J = 0 gap vs string mass) goes away at leading order.

**(ii) `T × T` has no light single traces at the orders computed [Q + D-sketch].**
* **CPY 1804.00949 §3 (𝒩 = 8)** [Q, verbatim]: "the S × S OPE contains two half-BPS operators: the stress tensor multiplet …, and the
  (B,+) multiplet whose component operators all have even twist. … the only bulk fields that contribute a meromorphic piece in the
  exchange diagrams are those in the stress tensor multiplet". An even-twist exchange must vanish ("would produce third order poles").
* **BCP 1906.07195 (𝒩 = 6, IIA)** [Q]:
  * "The leading tree-level contribution comes from supergravity exchange diagrams. These are corrected by higher derivative contact
    interactions …";
  * `M_SG` "contains an infinite series of poles that correspond to the stress tensor multiplet operators";
  * every other tree term is polynomial (through the R⁴ order they compute).
* **Bulk [Q CHLPS §5 + D-sketch].**
  * KK modes couple to two zero-mode gravitons at two derivatives only through `∫√g χ(y) = 0`.
  * A KK coupling of the `φW²` / `XWW` type is generated only by higher-derivative (α′³R⁴, with internal curvature) terms.
  * Estimate: its contribution relative to the R⁴ term is `~ (α′/L²)³ ~ λ^{−3/2}` [M].
* **The only T-multiplet member in `T × T` is T itself.** The other stress-tensor-multiplet bosons are R-charged or spin 1 (BCP Table:
  the U(1) flavour current is an R-singlet but J = 1).
* **Caveat.** For 𝒩 = 6, `⟨SSSS⟩` fixes only 2 of the 5 superconformal structures (BCP §), so the `⟨TTTT⟩` statement rests on the bulk
  3-point argument, not on `⟨SSSS⟩` alone [D-sketch].
* **The referee's breathing-mode example.** It belongs to the even-twist class if, as I recall for AdS₄×S⁷, it has Δ = 6 = 2Δ_T (I did not
  re-derive this). Then it is degenerate with `T_{μν}T^{μν} = [TT]_{0,0}`, has no exchange pole, and sits at a double zero of every
  functional.
* **Conclusion.** "ABJM has light R-singlet single traces in `T × T`" (UT §1 A1 row, §9; REFEREE M12/M13) is **not supported at
  leading order in 1/λ**. At these orders ABJM-IIA satisfies even the `T × T` form of A1. It satisfies A1′ with room to spare: any light
  coupling generated at higher order has J ≤ 2, and is either favourable (D; Q for J = 2) or covered by (Q-a)/(Q-b) with a
  λ-suppressed δ.

**(iii) What the bound constrains there [D-sketch + Q].**
* The external T is the zero-mode graviton. At tree level no KK state is exchanged between zero modes (internal-momentum conservation /
  orthogonality). So the 4d MHV `f(s,u)` of `⟨TTTT⟩`'s flat limit is the **10d IIA amplitude at zero internal momentum**, with
  `8πG₄ = 8πG₁₀/Vol(CP³)`.
* Hence **`g₄` (4d, our definition) = the 10d R⁴ coefficient restricted to 4d kinematics, `g₄M⁶/8πG₄ = 2ζ(3)`** (CHMRS normalisation;
  BCP's flat limit eqs. 4.14–4.15). It is not a different "4d coefficient".
* The light side of our functionals is the AdS₄ graviton exchange plus contacts. That is exactly ABJM's tree-level light content, so the
  4d IR logarithm `17.24 ln(2Λ/3)` in the upper end is real for this sum rule, though not optimal for a 10d theory. It is why the bound
  is loose: against CHMRS's 10d window, 2ζ(3) sits at 80 %.
* What our functional bounds is `g̃₄ = g₄ + O(Λ⁻²)`, given (F) [D] and H_EFT(b) [A]. H_EFT(b) is natural here: the light EFT is 4d
  supergravity plus 10d-generated contacts, with no light exchanges in `T × T`.

**(iv) Is the 2ζ(3) placement now an in-assumption test? Partly, and only weakly.**
* **What changed.** ABJM (IIA, 't Hooft, at leading order in 1/λ) now satisfies the setting (A1′, A2), and Λ ↔ ML is the spin-4 gap. So
  "2ζ(3) at 3.34 % (Δ_g = 300) / 2.59 % (1000) of the upper end" is a consistency check of a theory inside the hypotheses. It no longer
  merely checks normalisation. The corresponding λ are ≈ 2.6·10⁷ (300) and ≈ 3.2·10⁹ (1000).
* **What did not change.**
  1. (Pos) is assumed, not shown, for ABJM's string states in AdS. Their flat images are positive by 4d unitarity of the restricted
     amplitude and the certified flat positivity of NF02/Q3 on all `m ≥ M` states [D-sketch]. AdS (Pos) is [NE + M] at best, as for any
     theory.
  2. At 3 % of the upper end the check could only have failed grossly. It says nothing about the tightness of the bound.
  3. Λ ↔ ML holds at leading order in 1/(ML). With twist(J = 4) ≈ ML − 4, Λ ≈ ML − 5.5, so `(Λ/ML)⁶` shifts the 3.3 % / 2.6 % by
     about −11 % / −3 % relative. The placement is meaningful to one digit.
  4. A2 needs `g_s → 0` (the 't Hooft limit). BCP's finite-`g_s` genus-one R⁴ term is a loop effect outside A2.
* **Honest status.** An in-scope sanity check [D-sketch for "in scope"], not evidence for (Pos) and not a sharp test. It may be stated in
  the letter as such, and no longer needs the "outside A1" disclaimer at leading order in 1/λ.
* **M-theory regime (CPY): not covered.** There is no string scale, and the spin-4 gap in `T × T` is not at a parametrically separated
  scale.

---

## 8. Letter wording (proposed)

**Assumption paragraph (replaces A1):**
> *Higher-spin gap.* We assume that every single-trace operator of spin J ≥ 4 in T × T has twist at least Δ_g [CEMZ; CHLPS]. Spin 1 and
> 3 do not occur in T × T. Light single traces of spin 0 and 2 are allowed at any dimension. In the helicity-preserving sector that our
> functionals probe, their exchange diagrams have Regge spin 3 (J = 0) and 4 (J = 2); the crossed-channel spin-2 Regge pole of a massive
> state is helicity-flipping and is annihilated by our k = 2 structure. Exchanges of Regge spin below k drop out of the k-th sum rule
> together with their own double traces, so light states do not enter the k = 5 sum rule, and only their top Regge coefficient enters at
> k = 4. In the k = 2, 3 sum rules they enter exactly like heavy states. The Result-D functional is positive on all of them at every mass,
> so the cubic bound holds as stated and in addition bounds the summed squared couplings of light states to two gravitons, as in flat
> space [CHLPS]. The published Result-Q functional is positive on light spin-2 states and on scalars heavier than 0.46Λ. For lighter
> scalars we add 0.038 Ω^EE, which relaxes the upper end to 123 (Δ_g = 300) / 161 (1000). Here g₄ is the R⁴ coefficient with the light
> fields kept in the low-energy amplitude.

**ABJM paragraph (replaces the "outside A1" paragraph):**
> In the type-IIA 't Hooft regime of ABJM the higher-spin gap is the string scale, M² = 4/α′, and T × T contains no light single traces
> at tree level (only the stress-tensor multiplet is exchanged [CPY, BCP]). The flat limit of ⟨TTTT⟩ is the ten-dimensional amplitude at
> zero internal momentum, so g₄M⁶/8πG = 2ζ(3) [CHMRS, BCP]. This lies at 3 % of our upper end. That is a consistency check of a theory
> inside our assumptions at leading order in 1/λ, not a test of the positivity hypothesis; the bound is loose there because of the
> four-dimensional infrared logarithm.

**Table I / hypothesis rows:**
* replace "A1 [A]" by "A1′ higher-spin gap [A]";
* add "light J ≤ 2 states: drop-out at k = 4, 5 [D]; k = 2, 3 action positive (D, all masses; Q, J = 2 and m > 0.46Λ) [P flat, D-sketch AdS]".

---

## 9. What remains

1. **AdS light-state actions at finite Δ_L** for `Ω^EE` and the k = 2, 3 parts of `Ω^Q_up`, for J = 0, 2 and Δ_L from O(1) to Δ_g. The
   flat form is proved; the AdS (Pos_light) is [D-sketch] for Δ_L ≪ Λ and [I] near the gap. The project's exact per-state engine could
   do this for a few (Δ_L, J) rows.
2. **AdS normalisation of the light spin-2 k = 4 arc** (the `12|g₂|²/m⁴` derivative term). Flat [Q]; the AdS route is F_NORM-type
   calibration [D-sketch]. Only the favourable sign is used.
3. **A matter-aware Q functional** (an LP with light-state positivity, as CHLPS do in SDPB), to recover most of the ×1.7.
4. **ABJM.**
   * The bulk argument for `⟨TTTT⟩` beyond the two `⟨SSSS⟩`-visible 𝒩 = 6 structures.
   * The breathing-mode dimension (recalled, not re-derived).
   * The order at which higher-derivative terms first induce `⟨T T O_KK⟩` [M estimate λ^{−3/2} relative].
5. **The k = 2, 3 "block action only" link in AdS** (L0 applied to the light exchange's double traces at n > k) [D-sketch].
6. **The k = 4 borderline robustness** in this note is at three angle sets. LIGHTSIDE's alternative grid was not re-run for the
   degree-6 pieces; the large-n grid noise beyond n ≈ 110 is documented in §5.2.

---

## 10. Reproduce

```bash
cd ads4/code/a1p
(ulimit -v 3000000; OMP_NUM_THREADS=1 python3 a1p_flat.py)        # §4.1-4.3 flat light-matter actions, CHLPS check   -> cache/a1p_flat.txt   (s)
(ulimit -v 3000000; OMP_NUM_THREADS=1 python3 a1p_flat2.py)       # exact root isolation, m0, eps*, repaired constants -> cache/a1p_flat2.txt  (s)
python3 a1p_pieces.py                                             # new Mellin pieces, top coefficients, dM self-check (s)
for k in 2 3 4 5; do (ulimit -v 3000000; OMP_NUM_THREADS=1 python3 a1p_run.py profile $k P1); done   # 5-8 min each -> cache/a1p_prof_k*_P1.pkl
for t in P2 P3; do (ulimit -v 3000000; OMP_NUM_THREADS=1 python3 a1p_run.py profile 4 $t); done      # §5.3
for k in 2 3 4 5; do python3 a1p_run.py report $k P1; done; for t in P2 P3; do python3 a1p_run.py report 4 $t; done   # -> cache/a1p_report_k*_*.txt
```

| file | content |
|---|---|
| `a1p_flat.py` | CHLPS light-matter sum rules; actions of NF02 (D), Q3 (Q up) and Ω^Q_low; CHLPS F_g3 validation; drop-out contrast |
| `a1p_flat2.py` | exact sign certificates (sympy), crossing mass 0.4602Λ, ε*, repaired Q constants, L_A / R₀ reproduction |
| `a1p_pieces.py` | sympy-exact σ-pieces (xs/xt with degree-6, 8 residues, xu6, c6, xu4) and top-Regge representatives |
| `a1p_run.py` | term-wise s-channel functional (LIGHTSIDE engine, read-only) with the single-trace block; reports with tail/3-parameter extrapolation |
| `cache/a1p_report_k*_P*.txt` | the numbers of §5 |
