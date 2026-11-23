# (F): the R⁴ normalisation of Result Q's light side, by calibration

*23 Nov 2026, (F) agent. Unrefereed, not committed. No existing file was edited. New code: `ads4/code/resultq/fn_scalar.py`,
`fn_spin.py`. Outputs: `resultq/cache/fn_*` (git-ignored). Read-only imports: `exactskel/exactmom.py`, `uniform/skel.py`, `lsd_core`,
`lsd_mellin`, `lsd_pos`. Literature read from arXiv LaTeX sources: CHMRS 2106.10274 (§`sec:lightAction`, §`sec:heavy`), CLS-D 2311.04271
(§2.2, §5.6), CHMRS 2102.08951 (§3.8). Compute: about 25 CPU-min in total, at most 2 processes, 1 thread each, < 1 GB each.*

Grades: **[P]** proved; **[D]** derived; **[D-sketch]**; **[DV]** derived and checked; **[NE]** numerical evidence; **[M]** model;
**[Q]** quoted; **[A]** assumption; **[I]** open.

---

## 0. Summary

**Verdict.** The normalisation is confirmed. No factor changes: **−𝒜₄[R⁴] = 2g₄**, i.e. −Ĉ₄^{ds}|_{R⁴} = 2g₄ in Stage 0's
normalisation. The all-plus R⁴ structure (`ĝ₄`) gets 0.
* **Sign.** The statement must carry the minus sign: with `g̃₄ = −½𝒜₄[light]`, the consistent form is `𝒜₄[R⁴] = −2g₄`.
  * The task's "𝒜₄[R⁴] = 2g₄" drops it.
  * LIGHTSIDE_DIST §4 writes "arc(R⁴) = 2g₄ (−Ĉ₄^{ds} = 2g₄)". The bracket is right; "arc = 2g₄" should read "−arc = 2g₄".
  * The sign agrees with `g₄ ≥ 0`: the heavy value piece is positive, so −light > 0.

**How.** By calibration (route II, §2.3), not through a spinning flat-space-limit computation.
* Apply the same sum rule to an auxiliary tree diagram: the exchange of a heavy scalar χ coupled to `W²` (with either sign of parity).
  * Its Regge spin is 3, below the k = 4 borderline, so Lemma A sets its total term-wise value to zero.
  * Hence "arc of the R⁴ contact it generates" = −(heavy χ block's value).
* That block's value is flat-normalised by construction: it is the same heavy normalisation that every heavy state in Result Q uses.
* It equals `2g₄^χ` by the flat `B₄⁽¹⁾` sum rule applied to the flat χ amplitude.
* The arc is m-independent at fixed coupling, so the identity is exact.
* Varying the phase of the χW² coupling shows that the all-plus coefficient is invisible.

**Grade of (F).** [D-sketch] → **[D]**, given three inputs that are already on Result Q's list (rows 4, 6a and the heavy normalisation):
1. Lemma A in its term-wise form, for a J₀ = 3 exchange;
2. `H_EFT(a)` at leading order in the ν-paired (heavy-normalised) form;
3. the CLS-D heavy dictionary for a J = 0 state at ν = O(1).

The derivation needs **no spinning saddle and no helicity-factor computation**. Route I (the direct one: Penedones plus the spinning light
saddle) remains [I] for spin. It is not needed.

**Checks.**

| check | what | result | grade |
|---|---|---|---|
| A1 (§3.1) | light values of Penedones-normalised contacts under CHMRS's heavy-normalised `C_{k,ν}`, own exact implementation | `−C₄,ν[M₂,₀] = 4.0` (= flat) at ν = 0, 1, 3, 7.5, for Δφ = 3 and 2.3. All other CHMRS closed forms reproduced exactly. `C₄[M₁,₀] = C₄[M₀,₁] = 0` | [DV] |
| A3 (§3.2) | top-degree Mellin coefficient (Penedones) against the coupling defined by the heavy density, on exact AdS data only | `R3 = a_top/(c₄ q/2m⁸)`: 0.932 / 0.9956 / 0.99973 at Δ = 80 / 320 / 1280; fit `r0 = 0.9999990` | [DV] |
| A2 (§3.3) | calibration identity at the functional level: C[J=0 block] + C[EFT expansion to degree R] → 0 | R = 4 (top contact alone): 0.62 / 0.23 / 0.066 / 0.017 at Δ = 20 / 40 / 80 / 160 (∝ 1/m²). R = 9: 7·10⁻³ / 1·10⁻⁵ / 6·10⁻⁹ / 5·10⁻¹⁰ | [DV] |
| B1 (§3.4) | spinning-external position space (LIGHTSIDE_DIST's family): Lemma A for heavy σ-exchanges, single-trace included | symmetrised S_h + DT = 2.6·10⁻⁶ (Δ_h = 13), 1.7·10⁻⁴ (21) of S_h | [DV] |
| B2 (§3.4) | (P) on the EFT polynomials | V[q4] = V[r4]/2 to 3·10⁻⁷; V[q5] = (5/4)V[r5] to 8·10⁻⁸ | [DV] |
| B3 (§3.4) | calibration at **fixed celestial angles** in the same family | **fails**: S_h/(−A₄) = −0.80 / −0.36 / −0.98 (Δ_h = 13 / 21 / 29). The discrepancy comes from Mellin residues at n ~ ⟨δ⟩ | [NE], negative. Not a test of (F) (§3.4) |

**What remains between `g̃₄` and `g₄`.**
* **`H_EFT(b)` [A].** The light double-trace data are those of the tree-level bulk EFT, with no conspiracy in the contact tower. This is
  CHMRS's standing assumption, verbatim (2106.10274 §5.1): *"We assume that the light OPE data agrees with a tree-level AdS effective field
  theory (EFT) with a derivative expansion"*, and *"…double-traces with τ ≤ 2Δφ+2p is well-approximated by a bulk EFT with a truncated
  derivative expansion."* CLS-D §2.2 assume the same: *"One can argue that ω|_light is determined by the Wilson coefficients of the
  low-energy EFT in AdS."* It is therefore standard in exactly the literature Q builds on.
* **`H_EFT(a)` for spinning, J_h ≥ 4 heavy exchanges**, in ν-paired form [D-sketch].
* **O(Λ⁻²) light terms** (`b₅g₅`, curvature constants of higher contacts).
* **`b₀G|ĝ₃|²`**, which is 0 at ĝ₃ = 0.

So: **`g̃₄ = g₄ + O(Λ⁻²)·(EFT data)`, [D given Lemma A, H_EFT, heavy dictionary]**.

**(P), the g₅ cancellation (§2.6), [D-sketch].** The same calibration, applied at order u³ to `Λ²(∂Ĉ₄ − Ĉ₅)^{ds}`, gives:
* the D²R⁴ light value is exactly cancelled **iff** the J ≤ 3 `++` heavy-image correction in that combination vanishes as m → ∞;
* the measured `s_d ≈ +0.51` is J-, Δ- and gauge-independent;
* so exact cancellation and the measured σ_hU correction cannot both hold;
* the consistent reading is RU's "double counting" (the σ_hU term *is* the image of the g₅ non-cancellation).

The δ = 0.2 and δ = 0.05 upper bounds survive either reading (`s_fail > 0` in every channel). The "lower slope 0.49" reading does not
survive, and it is already not claimed.

**Consequences for the quoted bound.**
* **Numbers:** none change. 2L_AdS = 72.04 / 92.81 stands.
* **Statement:** Result Q now bounds the bulk Wilson coefficient: `−O(Λ⁻²) ≤ g₄Λ⁶/8πG ≤ 2L_AdS + …`.
  * Conditions: (F) [D] and `H_EFT` [A, CHMRS-standard].
  * Before, it bounded only `g̃₄` with (F) at [D-sketch].
* **ABJM placement:** unchanged. The ĝ₄ = 0 statement is not needed for it.

---

## 1. Objects and normalisations

| object | definition | source |
|---|---|---|
| flat amplitude | `M(1⁺2⁻3⁻4⁺) = ⟨23⟩⁴[14]⁴ f(s,u)`, `f_low = 8πG/(stu) + … + g₄ + g₅t + …` | CHLPS (2.1), (2.7a) [Q]; Stage 0 |
| flat sum rule | `B₄⁽¹⁾(p²) = ∮ds/(4πi) (2s−p²)/[s(s−p²)] [f(s,−p²) + f(p²−s,−p²)]`: low side `−B₄|₀ = 2g₄`; heavy side `2u²(|c⁺⁺|²+|c⁺⁻|²)` | CHLPS (2.26a), (2.32a) [Q]; one residue [P] |
| type II | `M_susy = 8πG/(stu) + g₀ + …`, `g₀M⁶/8πG = 2ζ(3)`, `M² = 4/α′`. Since `M(1⁺2⁻3⁻4⁺) = ⟨23⟩⁴[14]⁴M_susy`, `f = M_susy` and `g₄ ≡ g₀`. Hence `g₄/8πG = 2ζ(3)/M⁶` | 2102.08951 §3.8 "Maximal supergravity" [Q]; Stage 0 [DV]. My count of the arXiv source puts these lines near (3.31) / (3.39); the task's "(3.43)" could not be matched without compiling |
| `G`, `c_T`, `Λ` | `8πG = 256/c_T`, `Λ = √((Δ_g−2)(Δ_g−1))`, `g₄,AdS = g₄,flat/Λ⁸`, `Λ⁸Ĉ₄(Λp) → B₄(p)` on heavy blocks | UT §1; Stage 12 |
| `Ĉ₄^{ds}` | `∫χ(ν)Ĉ₄(ν)`, where Ĉ₄ = Ω₄ (term-wise s-channel functional, LSD §1.1) divided by the forward heavy normalisation N̂₄(ν). The heavy image of a state is `2u²(|c|²)·P_{−½}(cosh β)` at leading Regge, with `W(0) = ∫χN̂ = 1` | Stage 4 §0, §2 [D + DV] |
| `𝒜₄[X]` | Ĉ₄^{ds} on a tree term X, i.e. its Regge arc (zero for n_X < 4) | LSD §1.2 [D + DV] |
| `g̃₄` | `−½Ĉ₄^{ds}|_light = −½𝒜₄[light content]` (= minus the heavy sum, given H_K) | UT §1; LSD §0 |

**Which normalisations enter (F).**
* `Λ` enters only through the flat-unit bookkeeping (`Λ⁸`), identically on the heavy and the light side.
* `8πG` and `c_T` do not enter (F): both sides are linear in g₄. They enter only when g₄ is quoted as `g₄Λ⁶/8πG`.
* The de-shift enters through `W(0) = 1` (Stage 4: `∫χN̂ = 1.0000000000`; FH: stored sc3 moment `m₀ = 1`).
  * The R⁴ light value is ν-independent (§2.5), so the grid-artefact high moments `m_{n≥4}` (FH §4.4) never act on it.

---

## 2. Derivation

### 2.1 Flat side [P given Q]
* For constant f = g₄: `−B₄⁽¹⁾|_low(0) = Res_{s=0}[(2s)/s² · 2g₄]·(2πi/4πi)·2 = 2g₄`.
* For the flat amplitude of a single heavy state χ, e.g. `f_χ = A/(m² − t)` for a J = 0 `++` state, which only sees the like-helicity
  pair (14)(23):
  * `f_χ` is Regge-bounded (`f_χ ~ 1/s`), so `B₄⁽¹⁾` holds exactly;
  * `B₄|_high[χ] = −B₄|_low[f_χ] = 2f_χ(0,0) =: 2g₄^χ`.

  This is the flat calibration that route II lifts.

### 2.2 Route I (direct): status
1. The arc depends only on the top Regge coefficient (LSD (P)) [D; DV].
2. The top-degree Mellin coefficient of a contact equals Penedones' image of the flat vertex [Q for scalars, Penedones 1011.1485; tested
   here, A3]. Spinning: [D-sketch].
3. The universal Regge kernel, in the heavy normalisation, evaluates the flat `B₄` on the flat amplitude.
   * Scalars: CHMRS's light saddle (`BknuFlatSpaceLimit`) [Q], and exact Mellin algebra (A1) [DV].
   * Spinning: this is CLS-D's stated open problem (*"one has to derive a similar dictionary that relates the low energy EFT contribution…
     we expect … by arguments similar to section 3"*) [I].
   * The helicity factor `⟨23⟩⁴[14]⁴` would have to be evaluated at the Regge saddle in light-ray polarisation. CLS-D's own Regge-limit
     check (§5.6, the `j′ = 4` row of `M⁺_{4g}`) shows that the same coefficient-function combination appears in the Regge projection and
     in the heavy dictionary. That supports step 3, but does not prove it.

**Route I is not needed.**

### 2.3 Route II (calibration) [D given inputs 1–3]

**Auxiliary theory.** Einstein gravity plus one heavy scalar χ of mass m (Δ_χ → ∞), with `ℒ ⊃ c χ (W⁺)² + c̄ χ (W⁻)²` (AdS-covariant,
c ∈ ℂ). Take its single tree-level exchange diagram W_χ (all channels) in ⟨TTTT⟩.

1. **Regge spin [D].** The flat image is:
   * `f_χ ∝ |c|²/(m² − t)` (MHV; only (14)(23) are like-helicity);
   * all-plus `∝ c² Σ_x x⁴/(m² − x)`.

   At fixed u, large s: M ~ s³. In Mellin form, a degree-4 residue over one pole gives the same count. **J₀ = 3 < 4.**
2. **Lemma A, term-wise [D; DV for exchanges incl. single trace].** Ω₄[W_χ] = 0 at every a ≠ b, every ν, hence for every de-shift:

       p_χ Ĉ₄^{ds}[G_χ] + Ĉ₄^{ds}[DT(W_χ)] = 0.

   Evidence: LSD §3.4 (τ = 2 exchanges); B1 here (heavy exchanges, Δ_h = 13, 21); FH and A2 (scalar Mellin).
3. **Heavy dictionary [Q + D, already load-bearing in Q].** For a J = 0 state, β → 0 and `P_{−½}(1) = 1`, and corrections at fixed J are
   O(1/m²) (value piece). So:

       p_χ Ĉ₄^{ds}[G_χ] = 2u²|c^{++}_χ|²(1 + O(m⁻²)) = B₄⁽¹⁾|_high[χ](1 + O(m⁻²)) = 2g₄^χ(1 + O(m⁻²)).

   * This is the CLS-D crystal dictionary (2311.04271 eq. `Blockvspartialwave`, §5.2) with Result Q's heavy normalisation N̂₄.
   * It includes the identification of the block's OPE density with the flat coupling.
   * Scalar analogue: A2 gives `m⁸C_dens → 2` (1.99920 at Δ = 160); A3 gives the density ↔ coupling map [DV].
4. **Tree-level EFT matching [P, Lagrangian algebra].** Integrating out χ gives
   `ℒ_eff = (1/m²)·(c W⁺² + c̄ W⁻²)² + O(∇²/m⁴) + (curvature terms with fewer derivatives)`.
   * Its 8-derivative quartic part is, *as an AdS operator*, the R⁴ contact 𝒞[g₄^χ, ĝ₄^χ].
   * Its flat couplings are those of f_χ's constant term: `g₄^χ ∝ |c|²`, `ĝ₄^χ ∝ c²`.
   * The curvature-induced terms have Mellin degree ≤ 3 (n < 4), so `𝒜₄ = 0` on them [D; LSD].
5. **H_EFT(a) at leading order [D in scalar Mellin; DV scalar (A2); spinning ν-paired form D-sketch].**

       Ĉ₄^{ds}[DT(W_χ)] = 𝒜₄[𝒞[g₄^χ, ĝ₄^χ]] + O(g₄^χ/m²).

   * Mellin reason: the ν-paired functional's double-trace part samples M_χ only at Mellin arguments damped by the γ-factor weight
     (|T| ≪ m²), where the 1/m expansion holds.
   * **B3 shows that this is not true pointwise at fixed celestial angles** (§3.4). The input is the ν-paired form, the one FH and A2 test.

**Conclusion.** `𝒜₄[𝒞[g₄^χ, ĝ₄^χ]] = −2g₄^χ(1 + O(m⁻²))`.
* The left side is the functional on a fixed AdS diagram, linear in (g₄, ĝ₄) and m-independent at fixed couplings.
* Take m → ∞ with `|c|²/m^{…}` scaled to hold g₄^χ fixed. Then **`𝒜₄[R⁴] = −2g₄` exactly**, i.e. **`−Ĉ₄^{ds}|_{R⁴} = 2g₄`**.

### 2.4 The all-plus structure [D]
* With `c → e^{iθ}c`, g₄^χ is fixed while ĝ₄^χ rotates.
* The heavy image is θ-independent: Ĉ₄ is helicity-flip-free, and its flat image is MHV `B₄⁽¹⁾` (Stage 12 §1.3).
* Hence the ĝ₄ coefficient of 𝒜₄ is 0, consistent with `B₄⁽¹⁾ ∌ ĝ₄` (CHLPS). Sums over several χ span all three real R⁴ couplings.

### 2.5 Per ν and through the de-shift [D]
* The calibration can be run per pure moment ν, since the J = 0 heavy image is 2u² for every fixed ν as m → ∞. So the per-ν R⁴ light
  value is the constant −2g₄. This is (P)'s statement for R⁴ (degree 0 in ν²), derived without Mellin input. Scalar analogue: A1.
* Applied directly to the realised `Ĉ₄^{ds}`: the heavy image is `2u²W(0) = 2u²`, so `−Ĉ₄^{ds}|_{R⁴} = 2g₄W(0) = 2g₄`.
  * No per-ν polynomiality is needed.
  * The sc3 grid's high moments are irrelevant.

### 2.6 The g₅ cancellation (P) [D-sketch]
* Run route II at order u³ on `Λ²(∂Ĉ₄ − Ĉ₅)^{ds}`, with χ (J = 0 `++`) and with a J = 2 `++` state (also subcritical).
  * R⁴ contributes nothing: its per-ν value is constant (so ∂ = 0), and n = 4 < 5 in Ĉ₅. The curvature-induced 8-derivative terms
    contribute nothing for the same reason.
  * So at this order the light side is `λ₅g₅^χ`, with `g₅^χ = |c|²u³`.
  * The heavy side is the flat image (0 for J = 0 `++`) plus the curvature correction `δ(J,++)`.
* **Consequence.** Lemma A gives `λ₅ = −δ(0,++)/(|c|²u³) = −δ(2,++)/(|c|²u³)`. Exact (P) (λ₅ = 0) holds iff δ(J ≤ 3, `++`) → 0 at fixed J
  as m → ∞.
* **What is measured.**
  * Q's heavy correction is `FF → FF − a·2u³s_d`, with `s_d = +0.509 / +0.511` (J = 0 / 2) at both Δ_g, in all gauges.
  * It is helicity-odd (`σ_h`) and universal (BF2, FWD).
  * Dimensional analysis: Γ_K is in c₂ units (relative correction × m²), and `Λ²∂_{ν²}(2u²Γ/m²) = 2u³Γ′`. So the correction is
    m-independent at fixed J [D-sketch; not measured at m ≫ Λ].
* **Then.** λ₅ = −2s_d ≈ −1.0 ≠ 0. The light side of `−a(∂Ĉ₄ − Ĉ₅)^{ds}` carries `+2a s_d g₅`. Since `g₅ = ⟨σ_h|c|²u³⟩` by its own
  sum rule, this term exactly cancels the σ_h-universal heavy correction.
* **Equivalently.** The AdS combination acts in the flat limit like `∂B₄ − (1 + 2s_d)B₅`. "σ_hU on the heavy side" and "g₅ not cancelled
  on the light side" are the same term. Counting both is the double count RU suspected (RU issue 8; REVISION §3.3).
* **Effect.**
  * Under the consistent reading the heavy action is the flat-like one (σ_hU removed).
  * Every channel has `s_fail > 0` (Q §3.1): +1.78 / +2.81 / +2.38 / +5.72 / … at δ = 0.2, and +0.58 … at δ = 0.05. So both upper bounds
    survive.
  * The lower slope is the double-counted 0.255, matching CHLPS's 0.26.
* **Decisive cheap-ish test (not run).** The engine's `++` J = 0 derivative-piece correction for a state with m = 3–10 × Λ_functional:
  * calibration + exact (P) predicts → 0;
  * the double-counting reading predicts s_d ≈ 0.51 unchanged.

---

## 3. Checks

### 3.1 A1: light values of contacts (`fn_scalar.py`, `cache/fn_scalar.txt`) [DV]
* **Method.** Own implementation:
  * `B̂_{k,T}[M]` as the exact sum of finite residues (sympy), polynomial in T;
  * `C_{k,ν}` = CHMRS (`CfromBMellin`) prefactor × T-integral done by Barnes' first lemma in the rising-factorial basis;
  * cross-check by direct quadrature on Re T = 1 (agrees to all printed digits);
  * contacts `M_{a,b} = c_{2a+3b}(S²+T²+U²)^a(STU)^b` with CHMRS's Penedones normalisation `c_n`.

| contact, k | Δφ = 3: −C at ν = 0, 1, 3, 7.5 | Δφ = 2.3 | CHMRS |
|---|---|---|---|
| M₁,₀, k = 2 | 2.0 (all ν) | 2.0 | 2 |
| **M₂,₀, k = 4** | **4.0 (all ν)** | **4.0** | **4** (`C2C4Light`) |
| M₃,₀, k = 6 | 8.0 | 8.0 | 2³ |
| M₁,₀, M₀,₁, k = 4 | 0 | 0 | 0 (n < k) |
| M₁,₁, k = 4 | −59.5, −57.5, −41.5, 53.0 | −29.82, … | `2(ν² − 29.75)` / `2(ν² − 14.91)` exact |
| M₂,₀, k = 2 | 40020.5 at ν = 0 | 14799.9144 | `C2C4Light` closed form: 40020.5 / 14799.9144 exact |

* The flat value of `B₄` on `(s²+t²+u²)²` is 4, so the borderline light value equals the flat one **exactly and ν-independently**.
* This is the scalar version of (F), and of (P) for the R⁴-type contact.

### 3.2 A3: Penedones' top coefficient against the heavy-density coupling (exact AdS data, no functional) [DV]
* **Data.** J = 0 exchange with exactmom's density-normalised Mack residues `ρ_n = −½K_n` (Δφ = 3). Under this normalisation
  `m⁸(−1)^J C₄,ν[block]/q → 2`, the flat image of one state.
* **Flat side.** The exchange's flat amplitude, with the coupling read off the heavy density (`λ²/m² = q`), is `λ²Σ_x 1/(m² − x)`. Its
  `(s²+t²+u²)²` coefficient is `g₄^flat = q/(2m⁸)`.
* **AdS side.** The top-degree coefficient is `a_top = −½Σ_n ρ_n/s_n⁵`.
* **Penedones predicts** `a_top = c₄g₄^flat(1 + O(m⁻²))`.

| Δ | 40.5 | 80.5 | 160.5 | 320.5 | 640.5 | 1280.5 | fit `r0 + r1/m² + r2/m⁴` (Δ ≥ 80) |
|---|---|---|---|---|---|---|---|
| R3 | 0.7620 | 0.93234 | 0.98260 | 0.995632 | 0.998909 | 0.9997274 | **r0 = 0.99999901**, r1 = −446 |

* The top-coefficient → flat-amplitude map is **tested**, not assumed.
* Its finite-m corrections are large (r1 ≈ −446; 0.4 % at m = 300). They are a matching property of the generated coefficient, not of
  the contact's light value.

### 3.3 A2: calibration identity at the functional level [DV]
Entries are `(C₄,ν[block] + TP_R)/C₄,ν[block]`, where `TP_R` is the exact light value of the exchange's EFT expansion truncated at degree R.
The R = 4 column is the top contact alone, `−4a_top/c₄` by A1.

| Δ | m⁸C_dens (ν = 0) | R = 4 | 5 | 6 | 7 | 8 | 9 |
|---|---|---|---|---|---|---|---|
| 20.5 | 2.01435 | +0.619 | +0.298 | +0.125 | +0.049 | +0.018 | +6.8·10⁻³ |
| 40.5 | 1.99083 | +0.234 | +0.039 | +5.4·10⁻³ | +7.2·10⁻⁴ | +9.4·10⁻⁵ | +1.2·10⁻⁵ |
| 80.5 | 1.99695 | +0.066 | +3.0·10⁻³ | +1.2·10⁻⁴ | +4.3·10⁻⁶ | +1.6·10⁻⁷ | +6.3·10⁻⁹ |
| 160.5 | 1.99920 | +0.017 | +2.0·10⁻⁴ | +1.9·10⁻⁶ | +1.9·10⁻⁸ | +6.5·10⁻¹⁰ | +4.8·10⁻¹⁰ |

* The same holds at ν = 2 and ν = 5, within ×1.4.
* The leading-order residual falls like 1/m² (×4 per doubling of Δ), and each further order gains ≈ (const/m²).
* This independently reproduces FH §4's mechanism, with the top coefficient and C[contact] computed exactly.

### 3.4 Spinning-external position space: LIGHTSIDE_DIST's family (`fn_spin.py`) [DV + negative NE]
* **Setting.** T′ = σ ⊗ T_free, Result Q's swapped Ψ₂ (borderline k = 2 for σ-degree 4), C × C contour, P1 angles, production grid.
* **Exchange.** A heavy scalar σ-exchange (s + t + u), residues `R_m = ((Δ_h/2−1)_m)²/(m!(Δ_h−½)_m)`.
* **Numerics.** The double-trace (DT) value is evaluated on `X4 = W_h − P₃`, where P₃ is the degree ≤ 3 EFT polynomial (subcritical,
  term-wise value 0). This removes ~⟨δ⟩⁴-fold cancellations.
* **Range.** The grid resolves K(−n) only to n ≈ 170: the residues equal LSD's for n ≤ 150 and are noisy beyond 200. Sums use n ≤ 160
  plus a tail fit.

**B1. Lemma A for heavy exchanges, single trace included [DV].** Symmetrised:

| Δ_h | S_h (heavy block) | DT | S_h + DT |
|---|---|---|---|
| 13 | −1.003556·10⁻⁶ | +1.003554·10⁻⁶ | 2.6·10⁻⁶ of S_h |
| 21 | −6.66725·10⁻⁹ | +6.66614·10⁻⁹ | 1.7·10⁻⁴ (tail-limited) |

Per order the two sides differ by ∓3.10·10⁻⁸ (Δ_h = 13): that is the nonzero D_LR line value, which cancels in the symmetric combination.

**B2. (P) on the EFT polynomials [DV].**
* `V[q4]` with `q4 = δ₁₂⁴+δ₁₄⁴+δ₁₃⁴` equals `V[r4]/2` (top Regge term D⁴/8) to 3·10⁻⁷.
* `V[q5] = (5/4)V[r5]` with `r5 = (2−δ₁₃)D⁴/4` to 8·10⁻⁸.
* So a δ₁₃-dependent top coefficient is also all that matters.

**B3. The calibration at fixed celestial angles fails [NE].**

| Δ_h | ⟨δ⟩ = μ₄/μ₅ | DT / A₄ | S_h / (−A₄) |
|---|---|---|---|
| 13 | −6.0 | −0.80 | −0.80 |
| 21 | −15.7 | −0.36 | −0.36 |
| 29 | −30.3 | −0.98 | −0.98 |
| 41, 57 | −61, −120 | not converged at n ≤ 160 | — |

* **Where the discrepancy sits (`diag`).** `DT[X5] = DT[W_h] − A₄`, in units of |A₄|:
  * at n ≤ ⟨δ⟩/3 it stays at the 1–2 % level: +0.017 (Δ_h = 29, n ≤ 10), +0.020 (41, n ≤ 20), +0.014 (57, n ≤ 40);
  * subtracting the degree-5 term as well leaves −0.043 / −0.009 / −0.001, falling with ⟨δ⟩;
  * then +1.3 … +2.0 builds up over n ≈ 0.5 … 2⟨δ⟩ and does not shrink with Δ_h.
* **Why.**
  * At fixed angles the term-wise sum runs over all double-trace Mellin poles δ₁₂ = −n.
  * The EFT term of degree 4 + j has residues ~ n^{j−2} up to n ~ ⟨δ⟩, so every order contributes O(1)·A₄ there.
  * In CHMRS's ν-paired form, only finitely many kernel poles at bounded Mellin argument enter. That is why A2 converges.
* **Why this is not a test of (F).** In this family:
  * `S_h` is the σ-heavy block ⊗ the whole T_free OPE, not one flat-normalised heavy state;
  * the family is not holographic: σ-double-traces ⊗ twist-1 free bilinears sit at non-double-trace twists of T′ (Δ_{T′} = 5), so they
    carry unsuppressed `sin²` weights;
  * there is no flat dictionary.
* **What it does establish.** `H_EFT(a)` is **not** a pointwise-in-angle statement in the swapped functional's term-wise representation.
  Any spinning proof of route II's step 5 must use the ν-paired, heavy-normalised form. A spinning ν-paired check is [I].

---

## 4. Grades

| claim | before | after |
|---|---|---|
| (F) `−Ĉ₄^{ds}|_{R⁴} = 2g₄` | [D-sketch] | **[D given Lemma A (row 4), H_EFT(a) leading order, ν-paired (row 6a), CLS-D heavy dictionary]** (route II); scalar analogue of every step [DV] (A1–A3) |
| ĝ₄ (all-plus R⁴) gives 0 in Ĉ₄ | Stage 12 table [P given (F),(P)] | **[D]**, route II with a phase rotation |
| (P) for R⁴: per-ν light value constant | [D + DV structure] (LSD) | **[D]** (§2.5); scalar exact [DV] (A1) |
| Penedones top-degree ↔ flat coupling (via heavy density) | assumed | scalar **[DV]**, 10⁻⁶ (A3); spinning [D-sketch], not needed by route II |
| Route I (spinning light saddle, helicity factor) | [D-sketch] | [I]; not needed |
| Lemma A for heavy exchanges incl. single trace (position space) | [DV] for τ = 2 (LSD) | **[DV]**, Δ_h = 13, 21 (B1) |
| H_EFT(a) | [DV] scalar Mellin; [D-sketch] spinning | unchanged, plus: **pointwise (fixed-angle) version false [NE]** in LSD's family (B3); scalar ν-paired leading order [DV] (A2) |
| g₅ cancellation (P) | [D-sketch] | **[D-sketch]: equivalent to the vanishing of the m → ∞ limit of the ++ J ≤ 3 heavy correction; inconsistent with the measured m-independent s_d unless that correction disappears at m ≫ Λ**; δ = 0.2 / 0.05 headlines unaffected |
| `g̃₄ = g₄` | [D-sketch] | **`g̃₄ = g₄ + O(Λ⁻²)` [D given H_EFT (A; CHMRS-standard)]** |

---

## 5. Letter wording (proposed)

> The quartic coupling enters the light side only through the Regge arc of its contact diagram, and its normalisation follows from the
> same sum rule. Apply it to an auxiliary tree-level exchange of a heavy scalar coupled to W². That exchange has Regge spin 3, below the
> k = 4 borderline, so its total vanishes. The arc of the R⁴ contact it generates is therefore minus the heavy block's action, which is
> flat-normalised by construction. This gives −Ĉ₄|_{R⁴} = 2g₄ exactly, and zero for the all-plus structure, without a spinning
> flat-space-limit computation. (The scalar analogue of each step is checked to 10⁻⁶–10⁻¹⁰.) The identification of the CFT quantity g̃₄
> with the bulk Wilson coefficient then rests only on the standard assumption that light double-trace data are those of the tree-level
> bulk EFT [CHMRS], up to O(Λ⁻²).

Hypothesis-table row 5 ("(F)+(P) spinning"):
* **(F) [D]** given rows 4 and 6a;
* **(P) for R⁴ [D]**;
* **(P) for g₅ [D-sketch]**, tied to the σ_hU reading (§2.6).
* The "test" column "a J₀ = 4 spinning position-space correlator" is replaced by "route II (calibration); a ν-paired spinning check of
  H_EFT(a) at leading order".

---

## 6. What remains

1. **Spinning H_EFT(a) at leading order, ν-paired** [D-sketch → DV].
   * It is the only spinning-specific input of route II.
   * B3 shows it cannot be checked pointwise in angle.
   * It needs either the C_ν pairing on an exchange family, or a spinning Mellin (fixed-T-type) representation of Ψ₄.
2. **The m-dependence of the `++` J ≤ 3 derivative-piece correction** at m ≫ Λ_functional. This decides (P) for g₅ against the σ_hU
   reading (§2.6). It is not needed for the δ = 0.2 headline.
3. **H_EFT(b)** [A]: CHMRS-standard, unchanged.
4. **Route I for spin** (CLS-D's open light-side dictionary) [I]: no longer on the critical path.
5. The heavy dictionary at ν = O(1) for spinning blocks (CLS-D derive it at ν, Δ ≫ 1 and *"expect"* continuity, as CHMRS proved for
   scalars) is inherited unchanged. Every heavy image of Q already uses it.

---

## 7. Reproduce

```bash
cd ads4/code/resultq
(ulimit -v 3000000; OMP_NUM_THREADS=1 python3 fn_scalar.py)                 # A1, A3, A2 (~30 s)       -> cache/fn_scalar.txt
for o in ab ba; do (ulimit -v 3000000; OMP_NUM_THREADS=1 python3 fn_spin.py run $o P1 160 &); done; wait     # ~2.5 min each
for o in ab ba; do (ulimit -v 3000000; OMP_NUM_THREADS=1 python3 fn_spin.py runsh $o P1 &); done; wait       # ~3 min each
python3 fn_spin.py report P1      # B3 table                                                     -> cache/fn_spin_P1.txt
python3 fn_spin.py reportsh P1    # B1 (Lemma A for heavy exchanges)                             -> cache/fn_spin_sh_P1.txt
python3 fn_spin.py diag P1        # B2, and where B3's discrepancy accumulates                  -> cache/fn_spin_diag_P1.txt
```

| file | content |
|---|---|
| `fn_scalar.py` | exact `B̂_{k,T}` on polynomials (sympy residues), Barnes T-integral, CHMRS `c_n`; A3 Mack-residue moments (exactmom nodes) vs `q_{Δ,J}` (skel); A2 with exactmom.HeavyC |
| `fn_spin.py` | LIGHTSIDE_DIST machinery (lsd_core.Setup, P1 grid): term-wise residues of r4, r5, q4, q5 and X4 for heavy σ-exchanges (Δ_h = 13, 21, 29, 41, 57); single-trace sums; report / reportsh / diag |
| `cache/fn_spin_P1_{ab,ba}_N600.pkl` | the first run (N = 600; only n ≤ 160 is grid-resolved; no q4/q5). Superseded by the N = 160 files, which reproduce its n ≤ 160 numbers identically |
