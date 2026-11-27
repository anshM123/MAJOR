# Session 11: an upper bound on the elastic 3 → 3 contact D⁴R⁶ under 𝒩 = 8? Obstruction analysis

Status: 27 Nov 2026, research agent. Grades: [P] proof, [D] derivation, [D-sketch], [Q] quoted from literature, [A] assumption.
**Result: no bound is obtained.** The obstruction is located precisely below, and the parts that supersymmetry does and does not
remove are separated.

## 1. Route (a), flat space: superconvergence under a complex deformation [D-sketch]

**Ingredients that work.**
* 𝒩 = 8 tree amplitudes fall off as 1/z² under every BCFW shift (Arkani-Hamed–Cachazo–Kaplan 0808.1446; Benincasa–Boucher-
  Véronneau–Cachazo hep-th/0702032) [Q].
* At planar or tree level the six-point amplitude is meromorphic, so there are no multi-particle cuts. The analyticity problem of
  THEOREM_Q5 §3 (Landau singularities) is absent at the order we need.
* A UV completion that is soft at large z gives a superconvergence relation:

  (coefficient of a local term in A(0)) = −Σ_heavy Res + (light residues from lower-point amplitudes).

**Where it fails: positivity of the heavy residues needs forward elasticity.**
* A heavy pole in the shifted three-particle channel has residue A_L(1̂, 2, 3 → X)·A_R(X → 4, 5, 6̂).
* That is |A(123 → X)|² ≥ 0 only if the out-triple is the CPT image of the in-triple, i.e. exact forward kinematics p₄ = −p₁, etc.
* BCFW shifts in exact forward kinematics are degenerate. λ₄ ∝ λ₁, so a [1,4⟩ shift does nothing, and shifting a non-paired leg
  breaks the elasticity.
* No one-parameter null deformation is known that keeps a 3 → 3 configuration exactly elastic and varies s₁₂₃.

**Collinear bundling (the Chandrasekaran–Remmen–Shahbazi-Moghaddam trick)** does not rescue it.
* Grouping (1,2) collinear and 3 separate turns 3 → 3 into a composite 2 → 2. Gravity has no collinear poles, so the limit is
  finite.
* But the W₊³W₋³-type contact is built from ⟨ij⟩ among the three same-helicity legs and [ij] among the other three. It **vanishes**
  when two same-helicity legs are collinear (⟨12⟩ = [12] = 0). Mixed-helicity bundles change the elastic pairing.
* Whether *some* bundling keeps D⁴R⁶ nonzero and elastic needs EFK's explicit polynomial, which was not available here. [Open.]

**The remaining flat-space obstruction.** Exact forward elasticity produces graviton t-channel poles 1/t_{ii′} in every 2 → 2
subprocess with a spectator. These are the leading s²/t Regge terms, and maximal supersymmetry does **not** remove them. They
require CHMRS-type impact-parameter smearing in several momentum transfers at once, with the 4D IR logarithm, possibly squared.
Nothing like this has been constructed.

## 2. Route (a′), CFT reformulation: 3 → 3 as a composite 2 → 2, ⟨T A T A⟩ with A = [TT] [D-sketch]

**Setup.**
* Group the in-state as (two gravitons) + (one graviton), with the pair represented by the double-trace primaries
  A = [TT]_{n,J}. These are actual CFT operators, of dimension 6 + 2n + J + O(1/c_T).
* The elastic 3 → 3 becomes the mixed four-point function ⟨T A T A⟩. The paper's own dispersive four-point machinery applies to it.

**Good features.**
1. **Regge boundedness is available.** Rindler positivity and the chaos-bound argument hold for any wedge operators, composite A
   included.
2. **The heavy s-channel weights are positive.** They are Σ_X |⟨X|T A⟩|². For heavy single-trace X, ⟨X|T[TT]⟩ is the XTTT coupling
   g₃. **This is exactly the six-point diagonal.**
3. **Under 𝒩 = 8 the light side is computable.** The theory below the gap is unique 𝒩 = 8 supergravity: tree six-point, one-loop
   four-point (Alday–Chester–Raj computed AdS₄×S⁷ one-loop), and the R⁴ and D^{2k}R⁴ contacts (four-point data). There is no
   unknown light contact at D < 16 (session 10).

**Where it fails: positivity at subleading order in 1/c_T (Remark 7′ again).**
* In planar normalization the target term Σ|g₃|² enters ⟨T A T A⟩_conn at O(c_T⁻²).
* The same correlator has O(c_T⁻¹) "spectator" terms: graviton exchange and heavy exchange (g₂²) between T and **one constituent**
  of A.
* A positive functional that removes the O(c_T⁻¹) terms must annihilate them. The O(c_T⁻²) remainder then contains Σ|g₃|² (positive)
  **plus** O(c_T⁻²) corrections to the spectator terms: heavy-light loops, 1/c_T mixing of [TT] with heavy double traces, and
  anomalous-dimension corrections to the spectator exchange.
* Those corrections are **not sign-definite** and **not computable from light data**: they involve heavy states in loops. So the
  diagonal cannot be isolated with a sign.
* Supersymmetry makes the purely light corrections computable. It does not control the heavy-in-loop ones.

## 3. Route (b), type II check [not done]

* In type II, a 16-derivative six-graviton contact sits at relative order α′⁷. There the single-valued MZV content is ζ(7) alone:
  ζ(2)ζ(5) and ζ(2)²ζ(3) project to 0.
* So the independent D⁴R⁶ coefficient and the six-point part of D⁸R⁴ (also α′⁷ζ(7)) come with the same transcendental. Separating
  them needs the full α′⁷ six-point closed-string expansion and the EFK basis decomposition.
* That is beyond this session. Sign and size are **not** extracted.
* Expectation [A]: positive, since the UV completion is unitary and the channel is elastic, and of natural size ζ(7)(α′/L²)⁷.

## 4. Verdict (question c)

| obstruction | present? | removed by 𝒩 = 8? |
|---|---|---|
| unknown light contacts below D = 16 | no | yes (EFK, sessions 9–10) |
| light side computable (tree + one-loop supergravity) | — | yes |
| multi-particle analyticity (cuts) | absent at planar tree level | — |
| Regge boundedness | available in CFT (Rindler positivity, composite operators) | — |
| forward elasticity vs one-parameter deformation | **yes** | no |
| graviton t-channel poles in spectator subprocesses (4D IR) | **yes** | no |
| positivity at O(c_T⁻²) contaminated by heavy-in-loop spectator corrections | **yes** | **no** |

**No bound on D⁴R⁶ is obtained [D-sketch for the obstruction].** The decisive obstruction is the last row. It is the same
support-lemma structure that blocked β₅, now for the diagonal at subleading order in 1/c_T.

**Consequence for the 𝒩 = 8 η/s statement.** No explicit p can be justified unconditionally. The statement stays as in session
10: with natural sizes for the independent superinvariant coefficients [A], the deviation is localization-fixed to O(Δ_g⁻¹⁰) and
four-point-fixed to O(Δ_g⁻¹⁴).

**Most promising next step.**
* Build the **leading-order** version of the composite sum rule. Use a functional that annihilates all O(c_T⁻¹) spectator terms
  *and* their O(c_T⁻²) corrections, by choosing it to vanish on the full spectator double-trace family, not only at leading
  dimension. Then test numerically, on the AdS₄ supergravity data, whether its action on the heavy single-trace X block stays
  positive.
* This is a concrete, CFT-native SDP-type problem.
