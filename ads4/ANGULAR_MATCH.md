# ANGULAR_MATCH: order of limits for [E, E⁽¹⁾], the uniform constant, the Ward step, and degenerate multi-traces

*25 Nov 2026. Research agent, about 3.5 h. Unreviewed and not committed. No existing file was edited; BATCH4_FOLD.md and `code/b4fold/`
were not touched.*

* **Code.** New code is in `ads4/code/angmatch/`:
  * `am1_kksz.py`: convergence bookkeeping for weighted light-ray products, and the impact parameter of the defect examples;
  * `am2_constant.py`: W and K(g), symbolic and exact-rational;
  * `am3_ward.py`: the per-harmonic commutator from LADDER_ALLJ §4, and the Ward step.
* **Outputs.** They go to `code/angmatch/cache/`, which ignores everything but its own `.gitignore`.
* **Compute.** Under 10 s in total, one process.
* **Literature.** KKSZ (1904.05905) was read from a local arXiv source copy (scratchpad), at §4.1–4.6, eqs. (4.10)–(4.12), (4.52)–(4.55),
  (4.63)–(4.66) and (4.80)–(4.88). No web access was used.

**Grades.**

| grade | meaning |
|---|---|
| [P] | proved, or exact algebra |
| [D] | derived by hand |
| [D-sketch] | derivation outline; steps not all written |
| [DV] | derived and checked numerically or exactly at finite size |
| [NE] | numerical evidence |
| [Q] | quoted from the literature, checked against the source |
| [A] | assumption |
| [I] | open |

---

## 0. Verdict

| item | status | grade | one line |
|---|---|---|---|
| **1. Order of limits for [E, E⁽¹⁾]** | **partially resolved** | transfer of route N [D]; KKSZ condition [Q]; light-cone edge [D-sketch]; 𝒟(c_T) [I] | The M2 resolution (route N) transfers: the angular inequality is an exact statement at every finite c_T, heavy states enter only through c_T-independent per-state inequalities, and the planar J₀ = 2 never enters. Two things do not come along with it. (i) KKSZ's own counting puts the mixed product's **light-cone edge** at exact marginality unless KKSZ's asymptotic light-cone assumption (4.66) is used. (ii) The **defect 𝒟(c_T)** is the angular analogue of H_K. Nothing in §§0–13 gives 𝒟(c_T) = 0; §§17–18 give it conditionally on the regulated null-plane identity, which works at fixed c_T and so is compatible with route N. |
| **2. Uniform constant** | **resolved** | [P] (W, chain, numbers); [D] (shared normalisation) | W = 100(C−20)(C−30)/((T−20)(T−30)) holds exactly, with a₂ = ladder index k = 2 = helicity m = 6. K(g) = (729/4800)(g−2)²(1−49/(2g+9)²) is a proved lower bound for all J ≥ 6 and twist ≥ g. **K(1000) = 151266.27; bound 0.0178135** (the quoted "0.01782" is a safe round-up of 0.017814). The \|c\|² normalisation is shared exactly. |
| **3. Ward step Ĉ(0) = Ĉ(1)** | **resolved** | [P] given the §12 expression; that expression is now [P] given LADDER_ALLJ §4 | Ĉ(0) = Ĉ(±1) = \|a_{T,0}\|² = 1 exactly, for three reasons: there are no heavy harmonics 0, ±1; the 2→1 coefficient is 1 − \|k\|/2 = 0; and the first-moment operator annihilates both T helicities at (Δ, J) = (3, 2). So c₂ = 0. |
| **4. Degenerate multi-traces (Vv₀ = 0, K(Δ_*)v₁ = 0)** | **resolved for TT-built clusters**, modulo the uniformity already assumed | Vv₀ = 0: [D] (bulk), [D given A2] (CFT); K(Δ_*)v₁ = 0: [D-sketch given A2] | **Vv₀ = 0 holds in the strongest form, ΠVΠ = 0: there are no O(c_T^{−1/2}) splittings at all.** The reason is that every one-particle state is a descendant of the protected T. Both conditions are equivalent to two statements about tree-level ⟨TTTT⟩ (no log²u at O(1/c_T); no s-channel block off the double-trace zeros), and tree-level AdS₄ gravity satisfies both. |

**Note on the record.** ANGULAR_ROUTE §§14–18 were added while this work was in progress. I read §§0–13 as briefed and §§14–18 at the
end; §6 below relates the two.
* §17–18 give a conditional route to 𝒟 = 0 (the regulated null-plane identity). It works at fixed finite c_T, which is exactly what
  route N needs (§1.4).
* §17's "M2 transfer" is the statement checked in §1 here: **confirmed**, with the light-cone caveat added.
* §4.2 here closes §15's "trace-changing directions still to be verified" for a T-only light sector.

**Net effect on the angular bound.** √(t₄²+α_T²) ≤ 4√((3+𝒟)/K(g) + 𝓛_g) = 0.017814 at g = 1000 when 𝒟 = 𝓛_g = 0.
* 𝓛_g → 0 is now **[D given A2]**. A2 is the light-data assumption the Result already makes. It needs the §10 bounded-twist tail lemma (not re-checked here).
* The open physical inputs are now these two:
  1. **𝒟(c_T) = 0 at each finite c_T** (or limsup 𝒟 < ∞ for a weaker bound): [I] on §§0–13; conditional on ANGULAR_ROUTE §18's
     assumption 2;
  2. **separated-direction vanishing of [E, E⁽¹⁾] in momentum eigenstates**, which at the light-cone edge rests on KKSZ (4.66) [A, as in KKSZ].

---

## 1. Order of limits for the mixed commutator

### 1.1 What M2 was, and how it was resolved (REGGE_TREE §0, §4; supplement §HK "Order of limits")

**The worry.** Referee 3's M2: H_K uses the non-perturbative J₀ ≤ 1 for a tree-level statement, while the planar intercept is
J₀^planar = 2 (KKSZ (4.80)). Term by term the corner count then fails at k = 2, 3, 5.

**The resolution ("route N").**
1. The sum rule is applied to the **exact** correlator at finite c_T. There J₀ ≤ 1 (KKSZ (4.55)) and H_K holds with margin ½.
2. (Pos) is a statement about (Δ, J, parity). It holds at every c_T for every state above the gap, multi-traces included.
3. Only the light side is expanded in 1/c_T, under A2.
4. The heavy side is never expanded and need not be analytic in 1/c_T (CHMRS §4.4, the eikonal model).

The resolution is [D given A2], and INDEP_CHECK2 checked it against the sources.

### 1.2 Transfer to the angular inequality [D]

At every finite c_T the angular bound is the following chain.

| step | content | what it uses | c_T-dependence |
|---|---|---|---|
| a | [E(n₃), E⁽¹⁾(n₁)] = 0 for n₃ ≠ n₁, inside ⟨T\|·\|T⟩ | KKSZ (4.87) with J₀ ≤ 1 (§1.3), plus the light-cone edge (§1.3) | exact CFT statement |
| b | endpoint rule Σ_{O≠T,b}[(T−b(b−1))/3 + 3/2]\|a_{O,b,0}\|² = 3 + 𝒟 | §12/§13 algebra (re-derived in §3 below); the "3" is −Ỹ₀\|a_{T,0}\|² with \|a_{T,0}\|² = 1 fixed by the momentum Ward identity | exact; 𝒟 = 𝒟(c_T) |
| c | per state, twist ≥ g, J ≥ 6: \|a_{O,4,2}\|² ≤ [(T−12)/3 + 3/2]\|a_{O,4,0}\|²/K(g) | ladder (§8) and K(g) (§2 below): representation theory of (Δ, J) only | **c_T-independent**, the analogue of (Pos) |
| d | every endpoint term is ≥ 0 | unitarity, T ≥ b(b−1) on ladders | c_T-independent |
| e | fourth-harmonic rule Σ_{J≥6}\|a_{O,4,2}\|² ≥ \|c\|² | [E₄, E₋₄] = 0, the POVM argument of §12; needs only J₀ < 3 | exact |
| f | \|c(c_T)\|² ≤ (3 + 𝒟(c_T))/K(g) + 𝓛_g(c_T) | from b–e | exact at each c_T |
| g | c_T → ∞: \|c_∞\|² ≤ (3 + limsup 𝒟)/K(g) + limsup 𝓛_g | 𝓛_g → 0 from A2 (§4) | the only limit taken |

This is exactly route N.
* The heavy sector enters only through the exact sum rules (b, e) and the per-state (Δ, J) inequalities (c, d). The heavy sum is never expanded in 1/c_T.
* The light side is better behaved than for Ω^EE:
  * The T-family entries (the "3" and \|c\|²) are **exact** at every c_T. They are fixed by the energy Ward identity and the exact ⟨TTT⟩.
  * Only 𝓛_g (light multi-traces) needs A2.
* **The planar intercept J₀^planar = 2, at which the mixed product is exactly marginal, never enters.**

So the order-of-limits question in the narrow sense (c_T → ∞ against Regge) **is resolved by the same argument as M2**.

### 1.3 The KKSZ condition, checked against the source [Q]; and a new light-cone caveat [D-sketch]

**The quoted condition is right.**
* KKSZ §4.6 (titled "Other types of null-integrated operators") defines L^n(y) = ∫dv v^{n+1}T_vv in (4.86). So **E = L^{−1}** and the
  first moment **E⁽¹⁾ = L⁰** ((4.88); L⁰ is a weight-shifting descendant of 𝕃[T], (4.89)–(4.95), consistent with [𝖪_n, E] = −6E⁽¹⁾).
* Eq. (4.87) reads verbatim J₁+J₂ > J₀ + 3 + n + m: "the integral of the Wightman function ⟨Ω|O₄TTO₃|Ω⟩ is absolutely convergent in
  the Regge limit if …".
* With J₁ = J₂ = 2 this gives:
  * EE (n = m = −1): J₀ < 3;
  * EE⁽¹⁾ (n = −1, m = 0): **J₀ < 2**;
  * E⁽¹⁾E⁽¹⁾: J₀ < 1.

  The author's table is correct, and so is the formula as quoted ("+3+n+m" with KKSZ's n, i.e. moment = n + 1).
* KKSZ add: "If J₀ = 1 … we can only prove commutativity … for n + m < 0". Our pair has n + m = −1.

**Margins** (`am1_kksz.py`).
* Under J₀ ≤ 1 ((4.55), from Rindler positivity: the "J₁+J₂ > 2" bound (4.52) becomes 4 − 1 > 2) the Regge margin is **1**.
* At J₀^planar = 2 it is **0**.

**What (4.87) does not cover.** KKSZ: "One should also consider constraints coming from the lightcone limit. A full analysis … is outside
the scope of this work."
* In the α-parametrisation the weight v^{n+1} is |α|^{n+1} at the contour ends α → ±∞. Every condition from those ends therefore shifts.
  The Regge corner reproduces (4.87).
* At the **light-cone edge** α_i → ∞ the condition becomes J₁+J₂ ± (Δ₁−Δ₂) − 2(n_i+1) > 2 − τ.

  | sheet | value of τ | margin for EE⁽¹⁾ (weighted leg) | status |
  |---|---|---|---|
  | first sheet (OPE convergent) | τ = τ₀′ ≥ ½ | **τ₀′ ≥ ½** | rigorous |
  | second sheet, KKSZ's asymptotic light-cone expansion (4.66) | τ = τ₀′ | τ₀′ ≥ ½ | (4.66) is an assumption; KKSZ: "we don't have a general proof" |
  | second sheet, Rindler positivity alone ((4.64)–(4.65)) | τ = 0 | **0: exactly marginal (4 > 4 fails)** | — |
  | double-commutator integral | τ_* ≥ 1 | ≥ 1 | existence only, not commutativity |

* **Consequence.** For EE the Rindler-only light-cone margin is 2, and (4.66) enters only through the first/second-sheet boundary
  loophole. **For EE⁽¹⁾, (4.66) becomes load-bearing at the light-cone edge.** The record graded H_K "no light-cone uniformity assumption
  is used"; that cannot be said here. [D-sketch: exponent shift by my count, `am1_kksz.py`.]
* **Momentum eigenstates.** KKSZ establish commutativity on the dense set created by local operators in the configuration 4 > x, x⁺ > 3
  with 3, 4 spacelike, and state that "this subspace does not include … momentum eigenstates". The other configurations rest on the
  conformal-Regge expectation. The angular rules use ⟨T(p″)|·|T(p)⟩, so they inherit this. So does every event-shape sum rule, including
  [E₄, E₋₄] and KKSZ's own §5 rules. [A, shared]

### 1.4 What does not transfer: the defect 𝒟 is the angular H_K [I]

**What 𝒟 is.**
* In Ω^EE the analytic input is H_K, the interchange of the OPE sum with the functional integral. At finite c_T it is [D] from R1′
  (Rindler-positivity bound on the Regge box, margin ½).
* In the angular route the interchange is between the state sum and the harmonic sum (pairing with Q(θ)). Its failure is **exactly 𝒟**:
  * Σ_n q_n Ĉ(n) = c₀Q(0) = 0;
  * Σ_O Σ_n q_n Ĉ_O(n) = 𝒫 − 3;
  * the difference, 𝒟 = 𝒫 − 3 = lim B_N (§13), is a Fubini defect.

**𝒟 is a large-impact-parameter effect.** In all three abstract defect constructions of ANGULAR_ROUTE (§10: J = L³, Δ − J = L² + 2;
§12: J = L², Δ − J = L + 2; §13: J = 2^ℓ, Δ − J = ⌈J/ℓ⌉ + 2) the carrying states have e^β = (Δ+J)/(Δ−J) → ∞ (≈ 2L, 2L, 2ℓ;
`am1_kksz.py`). So 𝒟 lives in the same corner as H_K. [P for the examples; the general statement is a heuristic.]

**Why J₀ ≤ 1 does not give 𝒟 = 0.**
* KKSZ (4.87) gives absolute convergence at **fixed** separated directions. That yields the vanishing of the commutator, i.e. support at
  coincidence. It does not give uniformity as n₃ → n₁, which is what the high-harmonic tail (and so 𝒟) probes.
* §13's counterexample shows that no smoothness statement about the energy correlators closes this either.

**What route N needs.** Only 𝒟(c_T) = 0 at each finite c_T. No uniformity in c_T is needed, because step (f) then holds with 3 at every
c_T. If instead only sup 𝒟(c_T) < ∞ is available, the bound becomes 4√((3 + limsup 𝒟)/K + 𝓛_g).

**Natural route (not attempted).** An R1′-type bound, uniform on a neighbourhood of n₃ = n₁ and in the boost, for the mixed event shape
⟨T|E(n₃)E⁽¹⁾(n₁)|T⟩ at finite c_T, with J₀ ≤ 1. This is the angular counterpart of supplement §HK's corner count.
* A caution: at tree level the heavy k = 2 moment profile is η⁻² (CHMRS (4.45)). Whether the endpoint weights inherit a comparable profile
  is not known. Route N makes this irrelevant only once 𝒟(c_T) = 0 is shown at finite c_T.

### 1.5 Status of item 1

* **Resolved** [D]:
  * the c_T-against-Regge order of limits (route N transfers);
  * the KKSZ citation and condition ([Q], (4.86)–(4.88)).
* **New caveat** [D-sketch]: the light-cone edge of the mixed product needs KKSZ's (4.66).
* **Still open** [I]: 𝒟(c_T) = 0, the angular H_K.

---

## 2. The uniform constant with the shared normalisation (`am2_constant.py`)

### 2.1 Shared normalisation [D, with the relative factor 1/3 now P]

**Common units.** Both rules are written in the same Wigner-normalised amplitudes a_{O,m} = ⟨O,m|E(n)|T,+2⟩ (Fourier component at
harmonic k = 2 − m), in units where **|a_{T,+2}|² = 1**. This is the momentum Ward identity: the harmonic-0 part of E is P⁰ = 1 in the
rest frame.
* **[E₄, E₋₄] = 0** compares harmonic +4 (m = −2: b = 0 ladder k = 2, and the T flip entry) with harmonic −4 (m = 6: **b = 4 ladder,
  k = 2**). This gives Σ_{J≥6}|a_{O,4,2}|² = |c|² + Σ_{O≠T}|a_{O,0,2}|².
  * The T flip entry is c = (t₄ + iα_T)/4 relative to the diagonal 1: supplement §c1, H = 𝟙 + (t₄/4)(e^{−4iφ}σ₊ + h.c.); KERNELS_D3 §6.
  * Hence **|c|² = (t₄²+α_T²)/16** in the same units.
  * The t₄/4 normalisation is the standard d = 3 one: ⟨E⟩ ∝ 1 + t₄(|ε_{ij}n_in_j|²/ε·ε* − ¼).
* **The endpoint rule** has the same a's and the same unit (the RHS 3 = −Ỹ₀|a_{T,+2}|²).
  * Its relative factor between the translation (𝒟-operator) term and the first-moment term, i.e. the 1/3 in [(T−b(b−1))/3 + 3/2], was
    previously "not independently re-derived".
  * `am3_ward.py` (a) now derives the full §12 per-harmonic form **exactly** from the LADDER_ALLJ §4 identity
    (G = (i/6)WF, 𝒟 = 1 + 2cos φ + sin φ ∂_{θ₂}):

        Σ_O[Ḡ₁F₂ − F̄₂G₁ − (i/2)(F̄₁𝒟F₂ + (𝒟F̄₂)F₁)] = −i Σ_k {|a_k|²[cos kφ + (1+k/2)cos(k+1)φ + (1−k/2)cos(k−1)φ] + (1/3)a_k(Wa)_k cos kφ}

    This is exact Laurent-polynomial algebra for a generic real amplitude vector and a generic real symmetric W; with 1/6 in place of 1/3
    the identity fails.

**So the bound is homogeneous in the normalisation.** Only the ratio |c|²/|a_{T,+2}|² = (t₄²+α_T²)/16 enters.

**ANGULAR_GO (c)3 does not arise here.** Its worry was a common normalisation between the k = 2 form and the k = 3 fourth-harmonic
projection b′, behind the old 47266. The §10/§12 route replaces the k = 3 projection by [E₄, E₋₄], which lives in the same amplitudes as
the endpoint rule.

**Residual grade [D].** The amplitudes are those of LADDER_ALLJ §4, whose identity is [D]. The ∂δ cancellation "family by family" was
not re-derived.

### 2.2 W on the b = 4 ladder [P]

From the §8 closed form F_{4+k} = F₄C(J−4,k)(4)_k/(Δ−5−k)_k and the §4 norms N_m ∝ C(2J,J+m)Γ(Δ+m−1)Γ(Δ−m−1), symbolic in J and Δ:
* |a₁/a₀|² = 16(J−4)(J+5)/((Δ−6)(Δ+3));
* |a₂/a₀|² = 100(J−5)(J−4)(J+5)(J+6)/((Δ−7)(Δ−6)(Δ+3)(Δ+4)) **= 100(C−20)(C−30)/((T−20)(T−30))** (exact; it equals the ANGULAR_ROUTE §6 form and the
  product of the §8 normalised ladder ratios).

Here a_k = a_{m=4+k}. So **a₂ is ladder index k = 2, helicity m = 6, harmonic −4**, and it needs J ≥ 6.

### 2.3 The inequality chain [P]

Let R := [(T−12)/3]/W. Put x = J + ½ and y = Δ − 3/2, so that C = x² − ¼ and T = y² − ¼:

    R = (y² − 49/4)(y² − 81/4)(y² − 121/4) / (300 (x² − 81/4)(x² − 121/4)).

R increases with Δ: d/dT[(T−12)(T−20)(T−30)] has real roots 15.46 and 25.87, both below 30. So the worst case is Δ = J + g, i.e.
y = x + g − 2.

| step | claim | justification |
|---|---|---|
| 1 | (y²−a)/(x²−a) − y²/x² = a(y²−x²)/(x²(x²−a)) ≥ 0 | for y > x > √a, with a = 81/4, 121/4; needs J ≥ 5 |
| 2 | 1 − 49/(4y²) ≥ 1 − 49/(2g+9)² | 2y = 2J + 2g − 3 ≥ 2g + 9 ⇔ J ≥ 6 |
| 3 | min_{v>0}(v+A)⁶/v⁴ = (729/16)A² at v = 2A (second log-derivative 1/(3A²) > 0) | with v = x, A = g − 2 |

Together: **R ≥ (729/4800)(g−2)²(1 − 49/(2g+9)²) = K(g) for all real J ≥ 6 and Δ − J ≥ g > 2.**

Exact-rational scan of R(J, Δ = J+g) over J = 6…12g: R ≥ K everywhere; beyond 12g the step-3 bound increases.

| g | K(g) | min_J R (argmin J) | min incl. +3/2 | (min R)/K − 1 |
|---|---|---|---|---|
| 10 | 9.1537 | 10.4923 (18) | 10.5607 | 1.5e−1 |
| 30 | 117.845 | 119.923 (56) | 119.998 | 1.8e−2 |
| 100 | 1456.97 | 1459.47 (196) | 1459.54 | 1.7e−3 |
| 300 | 13485.33 | 13487.97 (596) | 13488.05 | 2.0e−4 |
| **1000** | **151266.27** | 151268.97 (1996) | 151269.05 | 1.8e−5 |
| 3000 | 1365051.3 | 1365054.0 (5996) | 1365054.0 | 2.0e−6 |

### 2.4 Numbers [P]

* K(1000) = 311491286964/2059225 = **151266.2710**.
* Bound 4√(3/K(1000)) = **0.017813513**. With the exact minimum it is 0.017813354; with the +3/2 as well, 0.017813349. The asymptotic
  value 160/(9g) is 0.0177778. 4√(3/K(g)) = (160/9)/(g−2)·(1 − 49/(2g+9)²)^{−1/2} identically, which is the §18 theorem's form.
* §9's cruder 151250 gives 0.0178145. The "0.01782" quoted in §13 and in the brief is a valid upper rounding of 0.017814.
* Other gaps: 0.1815 (g = 100), 0.05966 (g = 300), 0.005930 (g = 3000).

---

## 3. The Ward step Ĉ(0) = Ĉ(1) (`am3_ward.py`) [P given the §12 expression, which is [P] given LADDER_ALLJ §4]

**Setup.** Take the per-harmonic form of §2.1 and write Ĉ(n) for the exponential Fourier coefficients of the state-summed C(φ). The
form is even in φ by construction: real |a|², common ladder phase, real symmetric W. A signed harmonic k contributes:
* to Ĉ(±k): ½[|a_k|² + a_k(Wa)_k/3];
* to Ĉ(±(k+1)): ½(1+k/2)|a_k|²;
* to Ĉ(±(k−1)): ½(1−k/2)|a_k|²;
* with the obvious change at k = 0.

**The four facts.**
1. **No heavy harmonics 0, ±1** for O ≠ T (momentum Ward identity) [P].
2. **The 2 → 1 coefficient vanishes.** Harmonic k = +2 feeds n = 1 with 1 − k/2 = 0; k = −2 feeds n = −1 with 1 + k/2 = 0. Every
   other heavy harmonic, |k| ≥ 3, feeds only |n| ≥ 2 [P].
3. **The T family has harmonics 0 and ±4 only**, since |h| = 2.
   * Harmonic ±4 (the |c|² entry) feeds |n| ∈ {3, 4, 5}.
   * Harmonic 0 feeds n = 0 with coefficient 1 + d_T/3, and n = ±1 with coefficient 1.
4. **d_T = 0.** At (Δ, J) = (3, 2) the ANGULAR_ROUTE §7 operator annihilates e_{±2}:
   * the diagonal is T + C − 2m² = 2 + 6 − 8 = 0;
   * the unnormalised neighbour coefficients are −(Δ−m−2)(J+m+1) at m = 1 and −(Δ+m−2)(J−m+1) at m = −1, both equal to −(Δ−3)(J+2) = 0;
   * equivalently, the normalised b₁ = b₋₂ = 0.

   So the T family carries no first-moment term at harmonic 0 or at harmonic 4. This also proves the §12 step "the reduced first moment
   vanishes within the T family", previously not re-derived [P, given that the §7 operator applies to the T → T transition; the
   normalisation correction −i(Δ−3)/2 vanishes there].

**Result.** **Ĉ(0) = Ĉ(1) = Ĉ(−1) = |a_{T,0}|² = 1.**

**Checks.**
* On 200 random spectra (T family plus random b = 0 and b = 4 ladders) the identity holds to 0.0 in floating point.
* The q-pairing reproduces −3 + Σ_heavy[(T−b(b−1))/3 + 3/2]|a₀|² to 4·10⁻¹⁴.

**Contact.** C_even is supported at φ = 0, and Ĉ(n) = o(n⁴) (the §13 growth chain, not re-checked). So C_even = c₀δ + c₂δ″ and
Ĉ(n) = (c₀ − c₂n²)/2π. Ĉ(0) = Ĉ(1) then forces **c₂ = 0** (and c₀ = 2π). Q(0) = 0 kills c₀δ.

**Odd parts.** Odd parts (δ′, δ‴) would need sin-terms. The §2.1 identity shows there are none, so no Ĉ(1) = Ĉ(−1) argument is needed for
them. **Item 3 resolved.**

---

## 4. Degenerate multi-trace conditions (§13: Vv₀ = 0 and K(Δ_*)v₁ = 0)

### 4.1 Setting

**Notation.**
* ε = c_T^{−1/2}.
* A cluster is the finite set of exact primaries whose c_T → ∞ dimension is Δ_*, at fixed J. It is finite because each constituent has
  Δ ≥ 3.
* For each primary O_α: a_α = ε⁻¹ Σ_s λ_{α,s}K_s(Δ_α) (KERNELS_D3 §5), with
  * λ_α = λ⁰_α + ελ¹_α + …;
  * Δ_α = Δ_* + εδ_α + O(ε²);
  * δ_α the eigenvalues of V = Π(D − D₀)Π/ε on the degenerate space Π.

**The ε-expansion.**
* The O(ε⁻¹) term is λ⁰_α·K(Δ_*). It vanishes for every vector of GFF multi-particle states, by GFF factorisation: in
  ⟨t t L[t] t⟩_GFF every term contains ⟨0|L[t] = 0 or L[t]|0⟩ = 0.
* The O(1) term is **δ_α λ⁰_α·K′(Δ_*) + λ¹_α·K(Δ_*)**. These are §13's "K′(Δ_*)Vv₀" and "K(Δ_*)v₁".
* If both vanish, |a_α|² = O(ε²) = O(1/c_T).

**Zero structure** (KERNELS_D3 §5.7).
* S^{++} and S^{+−} vanish at twist 2 + 2n; S^{odd} (and odd J) vanish at twist 3 + 2n.
* These are exactly the twists of the parity-even and parity-odd GFF [TT]. So λ⁰ always lies along structures with K_s(Δ_*) = 0.

### 4.2 Vv₀ = 0 holds, in the strong form ΠVΠ = 0, for multi-particle states built from T [D]

**Bulk Hamiltonian.** Radial quantisation gives D = H_AdS = H₀ + εV + ε²V₄ + …, where V is the bulk cubic vertex. All cubic couplings
count here, including the t₄ and α_T ones.
* V changes the particle number by ±1 or ±3. Within an energy shell, ±3 needs a zero-energy triple and is impossible.
* By spectator factorisation, every element ⟨N±1|V|N⟩ inside a degenerate shell reduces to a **resonant 1 → 2 element**
  ⟨k_a k_b|V|k_c⟩ with ω_c = ω_a + ω_b.
* N → N elements vanish (odd number change). So it suffices to show that the resonant 1 → 2 elements vanish.

**Claim.** Every one-particle state is a descendant P^μ…|T⟩ of the single one-particle module, and Δ_T = 3 is exact.
* The exact descendant ψ(ε) = P_ε^n|T_ε⟩ is an exact D-eigenvector with eigenvalue 3 + n, for every ε.
* Expanding D_εψ(ε) = (3+n)ψ(ε) to first order gives **ΠVψ₀ = Π(3+n−D₀)ψ₁ = 0**.
* These ψ₀ span the one-particle subspace at level n. So V maps one-particle states to zero inside any degenerate shell, and by
  hermiticity ⟨k_a k_b|V|k_c⟩ = 0 at resonance.
* Hence **ΠVΠ = 0 on the whole Fock space of T**: no O(c_T^{−1/2}) splitting, δ_α = 0, and Vv₀ = 0 for every cluster.

**Inputs.**
* The perturbative bulk Hamiltonian, with ψ(ε) differentiable at ε = 0. This is the same level as A2.
* Gauge constraints of gravitons are not treated explicitly.
* The argument is [D] at that level.

**It fails when** a light single trace φ ≠ T has its own dimension shifted at O(ε), i.e. an extremal degeneracy Δ_φ ∈ Δ_a + Δ_b + ℕ
with a non-zero coupling. That is excluded here: T is the only light single trace (gap assumption) and Δ_T is protected.

**Classical cross-check** (citation to be checked, not used): this is the CFT statement of the absence of secular (resonant) terms at
second order in AdS gravitational perturbation theory.

### 4.3 Both conditions are tree-level ⟨TTTT⟩ statements [D given A2]

**Setting.** |a_α|² = ε⁻² λ_α^† K(Δ_α)K(Δ_α)^† λ_α is an exact identity at every c_T (conformal kinematics). Summed over a cluster it is
the integral of the PSD matrix K K^† against the cluster's ⟨TTTT⟩ s-channel spectral measure μ.
* Near Δ_*, K_s K_s^* has a **double zero** along the matched structures.
* So Σ_α|a_α|² ≤ ε⁻²[C ∫(Δ−Δ_*)² dμ_matched + ∫ K K^† dμ_unmatched].

| condition | sufficient statement about the O(1/c_T) correlator | why tree-level AdS₄ gravity satisfies it |
|---|---|---|
| (i) **Vv₀ = 0** ⇐ second moment ∫(Δ−Δ_*)²dμ = o(ε²) | an O(ε) splitting with λ⁰ ≠ 0 would put **log²u at O(1/c_T)**, with coefficient ∝ Σ_α(λ⁰_α)²δ_α² ≥ 0 | tree-level Witten diagrams (exchange and contact; D̄-functions) carry single logs only. Independently [D] by §4.2 |
| (ii) **K(Δ_*)v₁ = 0** ⇐ no unmatched weight at O(ε²) | the tree-level s-channel decomposition has **no block at a multi-trace dimension along a structure with K_s(Δ_*) ≠ 0** | at O(1/c_T) the s-channel content is T plus the GFF [TT]_{n,J}, each along its matched structure. Factorisation gives λ¹ = 0 for Wick-ordered multi-particle composites: every O(ε) term is (3-pt connected) × (1-pt or self-contraction) = 0. Triple-trace blocks (e.g. parity-even twist-3 [TTT], where K_e ≠ 0) do not appear at tree level |

**Under A2** the cluster measure equals the tree-level one up to o(ε²) in mass, supported within o(1) of Δ_*. Both bounds then give
**Σ_α|a_α|² = o(1) per cluster**.

Grades:
* (i) is [D] (§4.2) and [D given A2];
* (ii) is [D-sketch given A2]. The step that is a sketch is the structure-resolved statement about the tree-level block decomposition for
  spinning externals with a parity-odd cubic coupling; with parity it is standard.

**Why the 4-point function alone would not have sufficed.** It constrains Σ_α λ_α a_α. That sum vanishes at O(1) automatically
(⟨λ⁰|V|λ⁰⟩ = 0 by number parity), even with O(ε) splittings. The sum rules need Σ|a_α|². That is why the second-moment (log²) statement,
or §4.2, is required.

**Interpretation.** §13 left the "analytic justification (differentiable spectral projections, light transform vs projected expansion)"
open. The measure formulation avoids differentiable projections: it uses only Taylor bounds on K K^† and A2. The remaining analytic
input is A2 itself, read as a statement about the cluster spectral measure. That is its natural reading when degeneracies exist.

### 4.4 From clusters to 𝓛_g → 0 [D given A2 and the §10 tail lemma]

𝓛_g(c_T) = Σ_{twist<g}|a_{O,4,2}|².
* For J ≥ L the §10 bounded-twist tail lemma gives O(g⁶/L⁶), uniformly in c_T. That lemma is not re-checked here.
* For J < L there are finitely many clusters, each o(1) by §4.3.

So limsup_{c_T→∞} 𝓛_g ≤ O(g⁶/L⁶) for every L, i.e. **𝓛_g → 0**. **Item 4 resolved** for TT-built clusters, at the stated grades.

---

## 5. Updated status of the angular bound

| input | before | now |
|---|---|---|
| K(g), W, chain, 0.017814 at g = 1000 | constant verified; chain "not re-checked" | **[P]** |
| \|c\|² normalisation shared by both rules | open (ANGULAR_GO (c)3) | **[D]**; the 1/3 is [P] from LADDER_ALLJ §4 |
| Ward step Ĉ(0) = Ĉ(1); T-family first moment vanishes | not re-checked | **[P]** given LADDER_ALLJ §4 |
| order of limits, c_T against Regge | open | **[D]**, route N transfers |
| KKSZ condition J₀ < 2 for [E, E⁽¹⁾] | citation to be checked | **[Q]**: (4.86)–(4.88), margin 1 under J₀ ≤ 1 |
| light-cone edge of [E, E⁽¹⁾] | not considered | **[A]**: KKSZ (4.66) needed (Rindler-only is marginal) |
| momentum eigenstates | noted | [A], shared with every event-shape sum rule |
| 𝓛_g → 0 (degenerate multi-traces) | "Vv₀ = 0, K(Δ_*)v₁ = 0 open" | **[D given A2 + §10 tail lemma]**; Vv₀ = 0 also [D] independently |
| 𝒟(c_T) = 0 | central open estimate | identified as the angular H_K (large impact parameter); needed only pointwise in c_T. [I] on §§0–13; conditional on §18's assumption 2 |

---

## 6. Relation to ANGULAR_ROUTE §§14–18 (added to the record during this work)

| record item | relation to this file |
|---|---|
| §16 "Order of limits (author)": establish the spectral identity at each finite c_T, then take c_T → ∞; "the matching agent is checking the M2 analogue" | **Checked, §1.2: [D].** Route N transfers step by step. The heavy sector enters only through exact sum rules and (Δ, J) inequalities. |
| §17 "M2 transfer: lim_{c_T→∞}[lim_N 𝒫_N(c_T)], J₀ ≤ 1 < 2" | **Confirmed**, with the KKSZ citation fixed as (4.86)–(4.88). One addition: KKSZ's counting puts the light-cone edge of L⁰ at exact marginality without their (4.66) (§1.3). The record's assumption 2 (regulated commutators exist as local distributions, supported at coincidence) should state explicitly that separated vanishing of [E, E⁽¹⁾] includes this edge. With §17's light-scalar assumption Δ_φ > 1, the first-sheet light-cone margin rises from ½ to τ₀′ = 1. |
| §17–18 null-plane route to 𝒟 = 0 | §1.4 here identifies 𝒟 as the angular H_K and shows route N needs only 𝒟(c_T) = 0 at each finite c_T. The §18 regulator theorem works at fixed c_T, so **if** assumption 2 holds at every finite c_T, the order of limits adds nothing. No uniformity in c_T of the K⁰ bound is required. |
| §15 "No-log mixing … trace-changing directions are still to be verified" | **§4.2 closes the trace-changing directions** for a T-only light sector. Resonant 1 → 2 elements of the cubic vertex vanish because Δ_T is protected, so ΠVΠ = 0 between N- and (N±1)-particle states. With other light single traces, non-extremality must be assumed. |
| §17 "Trace-changing multi-traces … given differentiable fixed-window planar expansion under spectral projection and light transform" | §4.3 replaces the differentiable-projection requirement by A2 read as a statement about the cluster spectral measure. It uses only Taylor bounds on K K^†. |
| §18 theorem, 0.017814 at g = 1000 | Same number. K(g) is [P] here (§2.3), with W re-derived from the §8 closed forms. |

## 7. Reproduce (`ads4/code/angmatch/`; `ulimit -v 3000000; export OMP_NUM_THREADS=1`)

```
python3 -B am1_kksz.py       # < 1 s   -> cache/am1_kksz.txt     (weighted light-ray margins; e^beta of defect examples)
python3 -B am2_constant.py   # ~ 2 s   -> cache/am2_constant.txt (W symbolic; chain; exact scan; bounds)
python3 -B am3_ward.py       # ~ 7 s   -> cache/am3_ward.txt     (LADDER_ALLJ §4 -> §12 form exactly; d_T = 0; C^(0) = C^(1))
```
