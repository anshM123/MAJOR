# Session 9: independent five-point contact data under maximal supersymmetry (counting), and its lift to AdS₄ and Mellin space

Status: 27 Nov 2026, research agent. Grades: [P] proof, [D] derivation, [D-sketch], [DV] numerically verified, [Q] quoted from literature, [A] assumption.
Scripts:
* `n8_five_point_count.py` → `out_n8_five_point_count.txt` (D = 10–16);
* `n8_five_point_fast.py` → `out_n8_five_point_fast.txt` (D = 16, 18);
* `n8_four_point_check.py` → `out_n8_four_point_check.txt` (validation).

## 1. What was done, and what was not

* **Not done.** The full 3d 𝒩 = 8 five-point Mellin Ward-identity solution (the 3d analogue of the Rastelli–Zhou / Goncalves AdS₅
  construction) was not attempted. It is a research-paper-sized computation. No 3d 𝒩 = 8 five-point correlator exists in the
  literature (searched).
* **Done: the reduced version (4).** An exact count of independent **local** five-point superamplitudes of 4D 𝒩 = 8 supergravity,
  by mass dimension D (the number of derivatives of the operator).

## 2. The count [DV, exact arithmetic mod p = 2³¹ − 1]

**Method (EFK-style).**
* The five-point MHV superamplitude is δ¹⁶(Q)·f. Here f has little-group weight t_i⁻⁴ on every leg and is totally S₅-symmetric,
  and P = M(1⁻2⁻3⁺4⁺5⁺) = ⟨12⟩⁸f must be a polynomial in spinor brackets.
* For n = 5 the anti-MHV sector is the parity conjugate, so it has the same count.
* Procedure:
  1. enumerate every bracket monomial with weights (+4, +4, −4, −4, −4) and dimension D;
  2. evaluate on random momentum-conserving kinematics over GF(p);
  3. find the function space by rank;
  4. impose f(σ·k) = f(k) for the generators (12) and (12345).

**Validation at four points.** The method reproduces the known superinvariant tower R⁴, D²R⁴, D⁴R⁴, D⁶R⁴, D⁸R⁴, D¹⁰R⁴ with
multiplicities 1, 0, 1, 1, 1, 1. D²R⁴ vanishes on shell; at higher orders the count matches the symmetric polynomials in s, t, u.

**Five points.**

| D | operator | monomials | independent (−−+++) polynomials | S₅-symmetric local superamplitudes |
|---|---|---|---|---|
| 10 | R⁵ (the I₂I₃ ≅ IJ order) | 1 | 1 | **0** |
| 12 | D²R⁵ (the ∇I·∇J order) | 40 | 9 | **0** |
| 14 | D⁴R⁵ | 595 | 36 | **0** |
| 16 | D⁶R⁵ | 4983 | 97 | **0** |
| 18 | D⁸R⁵ | 29397 | 210 | **1** |

This reproduces and sharpens EFK (arXiv:1003.5018): R⁵ and D^{2k}R⁵ for k ≤ 3 have no supersymmetrization. **The first
independent five-point contact is D⁸R⁵.** It is unique, and the count agrees with the weight bound: f needs net ten square
brackets, so D ≥ 18.

The anti-MHV sector gives the parity image. A parity-invariant theory has one real coefficient.

## 3. Lift to gauged 𝒩 = 8 supergravity in AdS₄ [D-sketch]

**Dimension bookkeeping.** In gauged supergravity the gauge coupling g ~ 1/L_AdS carries mass dimension 1, so a superinvariant of
dimension d contains only terms ∂^{d−2j} g^{2j}. A 4D C⁵-type term (10 derivatives, e.g. IJ) can therefore arise only from
superinvariants with d ≥ 10, at order g^{d−10}:

| d | superinvariant | flat-space status | its g^{d−10} C⁵ coefficient |
|---|---|---|---|
| 10 | D²R⁴ type | on-shell trivial (no 4-pt and, by the count, no local 5-pt) | absent at g⁰ |
| 12, 14, 16 | D⁴R⁴, D⁶R⁴, D⁸R⁴ | no independent 5-pt invariant exists | fixed by the **four-point** coefficient through the unique supersymmetric completion [A: uniqueness of the gauged completion] |
| 18 | D⁸R⁵ | the first independent five-point coefficient | appears at g⁸ |

**Natural sizes** (IIA ABJM, Δ_g ~ λ^{1/4}, d-dimensional coefficient ~ (α′/L²)^{(d−2)/2}):
* β₅ = Σ_{k=2,3,4} κ_k·γ_{D^{2k}R⁴} + O(Δ_g⁻¹⁶);
* the leading term is O(Δ_g⁻¹⁰), from D⁴R⁴ (ζ(5)).

**Mellin-space lift [D-sketch].**
* By the flat-space limit (Penedones; Chester–Pufu–Yin for ABJM), a polynomial solution of the five-point Ward identities with top
  Mellin degree p maps to a nonzero local flat superamplitude of matching dimension D(p).
* So for D(p) ≤ 16 there is no independent top-degree contact solution. Lower-degree pieces are AdS-curvature dressings of
  higher-dimension invariants, as in the four-point case, where the 3d 𝒩 = 8 Mellin contact counting of CPY/BCP matches the flat
  superinvariant tower.

**Caveats.**
1. The ABJM flat-space limit is 11D (M-theory) or 10D (IIA), not 4D. A 4D zero count excludes 11D structures that survive
   restriction to 4D kinematics. Evanescent 11D structures (zero for 4D graviton kinematics) cannot produce a 4D C⁵ graviton
   coupling, so they are irrelevant for η/s [D-sketch].
2. The six-point **NMHV** sector, which contains the J² ↔ W₊³W₋³ (−−−+++) structure, was **not counted**. The MHV six-point sector
   needs D ≥ 20 by the same weight bound.

## 4. What Theorem H gives

* For a 3d 𝒩 = 8 SCFT at large c_T with a higher-spin gap, under the lift of §3:
  * the five-point stress-multiplet data relevant to the k = 2 hyperbolic coefficient (β₅ and the gradient Y-type D²R⁵ data) are
    **fixed by four-point (Mellin) data**, up to the single independent D⁸R⁵ coefficient c₁₈;
  * that coefficient enters at O(Δ_g⁻¹⁶) if natural.
* Combining with session 8:

  4πη/s − 1 = 27γ_{R⁴} + (D⁴R⁴, D⁶R⁴, … contributions, all four-point data) + O(c₁₈, six-point NMHV),

  with γ_{R⁴} fixed by localization.
* **What this does *not* give.** Supersymmetry removes the independent five-point data at orders 10 to 16. It does not *bound* c₁₈
  or the six-point NMHV coefficients. An O(1) c₁₈ in AdS units would still shift η/s at O(1).
* So the CFT theorem still needs inputs (a) and (b) of Theorem H for the *independent* superinvariant coefficients. Those now start
  at very high order (D = 18 at five points), where natural power counting gives Δ_g⁻¹⁶. The quintic obstruction of sessions 4–7
  (D = 10) is **removed by supersymmetry**. It is not removed *as an assumption*.

## 5. Remaining CFT-native steps

1. Solve the 3d 𝒩 = 8 five-point Ward identities in Mellin space, and confirm zero independent contacts below the D = 18 degree.
   The flat count predicts this, and the four-point analogue is consistent with CPY/BCP.
2. Count the six-point NMHV sector, both flat and in Mellin space.
3. Bound c₁₈ and the higher independent coefficients. This is a supersymmetric version of the higher-point dispersive problem, now at
   D = 18. Note: with supersymmetry, the five-point *interference* measure is still sign-indefinite, so the Gram argument again
   needs six-point data.
