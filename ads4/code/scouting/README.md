# Scouting for one PRL-level target (28 Nov 2026)

Written by a research agent. Nothing outside this file was edited, and nothing was committed. One small computation was
run from the scratchpad; its script is reproduced in the Appendix, so it can be re-run.

**Source marks**
* **[R]** full arXiv source read this session.
* **[A]** abstract read this session (INSPIRE or arXiv API).
* **[S]** HTML read through a summarising fetch (partial, paraphrased).
* **[M]** memory, not checked this session.
* **[C]** computed here (Appendix).
* **[I]** my inference.

"Open" means not found in INSPIRE citation lists or in web searches. That is a search-level statement, not a proof that
nobody has done it.

---

## 0. Bottom line (blunt)

1. **No candidate is clearly PRL-plausible within weeks.** The best expected odds of a PRL are about 10–15%. No idea
   here is better than that.
2. **Candidate (f), the Herzog–Stergiou bound `c_S/C_T ≤ 16πζ(3)`.**
   * It would be PRL-level (arguably more) if proven.
   * It is **not provable with our tools**, nor with any mechanism known to me, in weeks or months.
   * The "invisible tail" objection is decisive against every route that uses finitely many vacuum or `T₀` data. That
     covers all of our tools. The paper itself contains an explicit planar example of it (§1.5).
   * The weaker finite-N target (rigorous two-sided Ising bounds) is also out of reach.
   * What *is* tractable is a **test, and possibly a refutation, of the stronger pointwise hyperbolic envelope**, using
     thermal-order CFTs on `S¹×H²`. That is the recommended checkpoint (§4).
3. **Recommendation: ONE target, (f-H).** Look for a violation of the pointwise heat-capacity envelope `Ψ′(x) ≤ Ψ′_s(x)`
   on `S¹×H²` in a CFT with persistent thermal order: the large-N d=3 biconical model.
   * The work is a leading-large-N computation on `H²` with our verified `H²` spectral kernels.
   * Certified interval arithmetic would make a violation rigorous at leading order.
   * There is a hard go/no-go on day 7.
   * **Fallback:** (a′), stress-tensor-source energy correlators in the 3D Ising CFT, using our exact commutator rules.
     It must start at once if (f-H) is killed, because two groups have announced it.
4. **The field moved fast in Aug–Nov 2026, onto our tools.** These papers overlap results already in the Letter or
   Supplement and should be cited there whatever we decide:

| paper | overlap with our work |
|---|---|
| Herzog–Stergiou 2610.33179 [R] | `Ψ′(1)=1` is our `c_H(T₀)=π³C_T/12R²`. `Ψ″(1)=3/4+t₄/240` and the free-scalar boundary-term exception are the heat-capacity analogues of our `d log η_Y/d log T − 1 = −t₄/30` and its edge-term exception. |
| Guo–Li 2611.07394 [R] | Restates the shear sum-rule weight as `(3/4+t₄/16)P`, i.e. our `S_E = 3/8+t₄/32`. Gives the Ising value `t₄ = 3.8893305` (from `n_B, n_F`) and `O(2)` `t₄ = 3.79 ± 0.14`. |
| Ricci–Riembau–Son 2611.31587 [S] | Stress-tensor source in 3D with helicity states `|±2⟩`, Fourier multipoles `c_n, d_n`, the 3D Hofman–Maldacena (HM) analogue `H₄`, and the inequality `|d_n|² ≤ c_{n−2}c_{n+2}`. Says the Ising application is "left for future work". |
| Belin–Borak–Henriksson–Piron–Zhiboedov 2608.26246 [S] | Positivity of energy correlators gives positive quadratic sum rules for `TTO` couplings (spin 2 and spin 4), plus large-`c_T` implications. |
| Meçaj–Moult–(Stoffels)–Walters–Xin 2512.09986, 2608.19322 [S] | EEC-hedron, with bounds on 3D Ising and `O(2)` OPE data from scalar sources. List "spinning sources plus superconvergence" as future work. |
| Hartman–Mathys–Stergiou 2611.23370 [R] | ANEC ⇒ strong 4D a-theorem. Say the 6D light-ray sum rule for `Δa` exists, but they "have not found a way to prove positivity". |

---

## 1. Candidate (f): the free-scalar thermal envelope (evaluated first, as requested)

### 1.1 The paper exists, and what it says [R, full source]

C. P. Herzog and A. Stergiou, *The Free Scalar as a Thermal Envelope for Three-Dimensional CFTs*, arXiv:2610.33179
(23 Oct 2026). It is written in PRL format with a Supplement.

* **Flat conjecture.** `c_S/C_T ≤ 16πζ(3) = 60.4219…`, with `s = c_S T²` and `⟨TT⟩ = C_T 𝓘/x⁶` (one real scalar has
  `C_T = 3/32π²`). Equivalent forms:
  * **Casimir form:** `ε_vac(S¹_L×ℝ, NS)/C_T ≤ 16πζ(3)/3`, using `ε_vac = c_S/3` (Shaghoulian).
  * **Thermal-OPE form:** `a_T ≤ 4ζ(3)Δ` for any scalar's thermal stress-tensor coefficient.
  * **Operator-counting form:**
    `limsup log N(E)/E^{2/3} ≤ (3/2)(4ζ(3) C_T/C_{T,s})^{1/3}`.
* **Stronger pointwise conjecture on `S¹×H²`.** Define `x = T/T₀` with `T₀ = 1/2πR`, `Φ = R²s_H` and
  `Ψ = 12Φ/(π³C_T)`. Then the conjecture is `Ψ′(x) ≤ Ψ′_s(x)` for all `x ≥ 1`, i.e. the C_T-normalised hyperbolic heat
  capacity is bounded by the free scalar's.
  * Integrated, it gives `Ψ(x) − Ψ(1) ≤ Ψ_s(x) − Ψ_s(1)`, whose `x → ∞` limit is the flat conjecture.
  * The triangle relations are `Ψ(1) = 6F/(π⁴C_T)`, `Ψ′(1) = 1` and `Ψ/x² → (3/π⁵)c_S/C_T`.
  * They also observe `Ψ″ ≤ Ψ″_s` in every example, with no principle known to imply it.
* **Combined with Bueno–Casini's `C_T/F` conjecture.** This gives `c_S/F ≤ 8.99492847…`, the free-scalar value.
* **What they prove.** A local theorem near `x = 1`: `Ψ″(1) = 3/4 + t₄/240 ∈ [11/15, 23/30]`, from the HM / collider
  bound via Hung–Myers–Smolkin. The free scalar has `Ψ″_s(1) = 19/15`; the excess `1/2` comes from its `φ²` boundary
  term in the modular Hamiltonian. Hence `Ψ′_s − Ψ′ ≥ ½(x−1) + …`. The theorem needs (i) no dimension-one boundary
  scalar and (ii) right-sided smoothness at `x = 1`.
* **Global analytic checks (full `x > 1`).**
  * free Dirac fermion;
  * Einstein hyperbolic black hole (`Ψ_Ein = y²`);
  * large-N critical `O(N)` (monotone thermal mass);
  * large-N sextic line, including the BMB endpoint (`14πζ(3) = 52.87`).
* **Endpoint checks of the flat bound.**
  * Ising 51.08(13), O(2) 50.84(27), O(3) 50.32(23), from Monte Carlo Casimir amplitudes (Bulgarelli–Caselle–Nada–Panero
    2505.20403) plus bootstrap `C_T`.
  * Planar Chern–Simons-matter (CS-matter) theories, regular-boson range 45.32–60.42.
  * Gross–Neveu `45.32 + 0.08/N_f` and QED₃ `45.32 + 0.90/N_f`.
  * `O(N)` at next-to-leading order: `48.3376 + 1.5557/N`.
  * Planar Einstein-scalar branes: `≤ 4π⁵/27`.
  * Multiscalar fixed points in `4−ε` at `O(ε)`.
* **Loopholes they list.**
  * A first-order transition at `x > 1` would violate the pointwise conjecture.
  * Thermal order raises `c_S` at `O(ε)`.
  * The large-n biconical model *saturates* the flat bound at leading order throughout `3 < d < 4` (Smolkin–Yung
    2608.02720).
  * Their wording on proofs: "A proof … must therefore use constraints specific to conformal field theory"; they suggest
    reflection positivity, modular covariance, or a vacuum variational principle on `S¹×ℝ`.
* **Methods disclosure.** The Supplement says the idea and the calculations were produced by AI agents; the idea itself
  was proposed by a Claude model. Agents "tried, unsuccessfully, to find a proof".

### 1.2 Has it been resolved or challenged?

No, as of 28 Nov 2026 [A]. INSPIRE lists one citing paper:

* **Kusuki–Tang, *Entanglement Bounds on Cardy Growth in Higher Dimensions*, 2610.37449 [A].** It proposes `f ≤ 4π²a`
  in 4D and `f ≤ 2.99830949… F` in 3D.
  * I checked that `2.99830949 = (ζ(3)/2π)/(log2/8 − 3ζ(3)/16π²)` is exactly the free-scalar `f/F` [C].
  * Since `c_S = 3f`, this is **identical** to the Herzog–Stergiou corollary `c_S/F ≤ 8.9949`.
  * They derive it from a stronger Rényi inequality stated in terms of `⟨TTT⟩` and "saturated by every free theory".
    They establish it at weak gauge coupling, at large-N vector models and in holography.
  * So this is an independent parallel conjecture, not a challenge.
* **Buchel, 2610.38226 [A].** Neutral AdS₅ black branes with source-free scalar hair, made "entropically favored" by a
  scalar–Gauss–Bonnet coupling. This is d = 4 and does not cite Herzog–Stergiou. It shows that the "hairy branes lower
  the entropy" intuition behind their holographic check can fail once higher-derivative couplings are allowed. It is
  nowhere near the 60.42 value, so it is not a counterexample.

No counterexample, proof or criticism was found.

### 1.3 Related known results

* **2D.** Modular invariance gives `s = πcT/3`, so `c_S/c` is fixed; the analogue is an equality. In 2D the Casimir
  energy `−πc/6L` is likewise fixed [M]. Kusuki–Tang frame their bound as the higher-d inequality replacing the 2D
  Cardy = entanglement equality [A].
* **`F`, `C_T` and `c_S` are independent in 3D.** Kovtun–Ritz 2008, as cited in Herzog–Stergiou [R].
  * The F-theorem (`F_UV > F_IR`) is proved via strong subadditivity (Casini–Huerta 2012) [M].
  * `C_T` is not monotone under RG flow [M].
  * Bueno–Casini–Lasso–Moreno 2307.05164 conjecture `F(A)/F₀ ≤` free scalar, which implies
    `C_T/F ≤ 3/(4π² log2 − 6ζ(3)) ≈ 0.14887` [A].
* **Values of `c_S/C_T`.**

| theory | `c_S/C_T` | source |
|---|---|---|
| free scalar | 60.4220 | [R] |
| free Dirac fermion | 45.3165 | [R] |
| Einstein holography (`4π⁵/27`) | 45.3363 | [R] |
| critical `O(N)`, `N = ∞` (4/5 of free) | 48.3376 | [R]; reproduced [C] |
| BMB endpoint of the sextic line | 52.869 | [R] |
| **d = 3 biconical model, leading order in 1/N** | **53.1658** | **new, [C], §1.4** |
| 3D Ising | 51.08(13) | [R] |

* **Thermal bootstrap for `O(N)`: Barrat–Marchetto–Miscioscia–Pomoni 2411.00978 [A].**
  * They solve a truncated Kubo–Martin–Schwinger (KMS) condition for `⟨φφ⟩_β` and estimate one-point functions and `f`
    for `N = 1, 2, 3`, in agreement with Monte Carlo.
  * The method is not rigorous: the thermal OPE coefficients `a_O` have no positivity, and the tail is truncated or
    modelled.
  * Follow-ups also have no rigour: the thermal Polyakov bootstrap (Guo–Li–Shen 2609.28627 [A]) and neural-network
    approaches (2605.13183 [A]).
* **The only known global crossing constraint in `d ≥ 3` touching heavy data and thermal one-point functions.**
  Simmons-Duffin–Xu's genus-2 crossing, 2511.07569 [A]. It relates asymptotics of heavy-heavy-heavier OPE coefficients to
  squares of thermal one-point coefficients. Both channels are sums of squares, so it is at least positivity-compatible.
* **Thermal effective action.** In odd `d`, locality gives no relations analogous to the even-`d` ones that fix `a`
  (Advant et al. 2609.03013: "no similar relations at odd d") [A].

### 1.4 New data point: the d = 3 biconical model at leading order [C]

The biconical model is an `O(N)` vector `φ` plus one scalar `χ`. Thermal order in it is the only known mechanism by which
an *interacting* CFT saturates the bound: in `3 < d < 4` the flat direction cancels the thermal masses.

**Computation.** I took the d = 3 leading-order effective potential of Smolkin–Yung 2608.02720 [R]:
* gap solution `m = 2T arcsinh(½e^{2πρ/T})`, with `ρ = φ² + αχ²` and `α = −1`;
* forced sextic coupling `g = 15360`;
* `V_LO = V_free + (N/2)∫_{−∞}^{ρ}m² dρ′ + N g χ⁶/720`.

I evaluated it at their minimum and normalised by `C_T = N·C_{T,s}`, its value at this order.

**Results** (Appendix):
* minimum at `χ̄² = 0.059631 T` and `m = 0.674647 T`, reproducing their `0.0596` and `m² = 0.455`;
* `c_S/c_{S,free} = 0.879909`, i.e. **`c_S/C_T = 53.1658`**;
* check: the symmetric point (`χ = 0`) gives exactly `4/5` (Sachdev's `O(N)` value, `m = 2 log φ`).

**Lemma.** In d = 3 at leading order, for this class of models `V_LO − V_free = ½∫m² + gχ⁶/720 ≥ 0`. So no
`O(N)`-plus-singlet thermal-order model of this type can exceed the free value at leading order [D, for this class].

**Consequences.**
* In d = 3 the RG-generated sextic term removes the flat direction already at leading order. The saturation noted by
  Herzog–Stergiou holds only for `3 < d < 4`.
* The decisive flat-space test is therefore the **next-to-leading `1/N` sign of `c_S/C_T` in the `3 < d < 4` biconical
  model**. That would test the conjecture's `d > 3` extension, which Herzog–Stergiou suggest. It is an NLO large-N
  thermal calculation, which we have not built [I].

### 1.5 Tractability with our tools

* **CHM map plus the modular-variance identity.** Here `S(q) − S(1) = ∫_q^1 v Var_v(K) dv` and
  `c_S/C_T = 4π²R² lim_{q→0} q²∫_q^1 v𝒱(v)dv`.
  * The identity is exact. Its normalisation agrees with Herzog–Stergiou's `Ψ` (`Ψ′(1) = 1` ⇔ `c_H(T₀) = π³C_T/12R²`).
  * But `q → 0` is the high-temperature end. Vacuum data fix only the Taylor coefficients of `Var_v(K)` at `v = 1`,
    namely integrated `n`-point stress-tensor functions.
  * Our own Theorem C and the optimal-recovery results (`code/cft_etas`, `code/math_transfer`) show that even complete
    Euclidean data plus bounded analytic continuation yield only logarithmic control far from the data point.
  * The only positivity, `Var_v(K) ≥ 0` (convexity of `log Tr ρ^q`), gives *lower* bounds on Rényi entropies for
    `q < 1`, never upper ones.
* **`c_H(T₀) = π³C_T/12R²` and `c_H(T₀) = 4πη_Y(T₀)`.**
  * The first is Herzog–Stergiou's `Ψ′(1) = 1` (Perlmutter 2013), so it is not new.
  * The second is ours, but local at `T₀`: it converts the conjecture near `x = 1` into a viscosity statement and nothing
    more.
* **Energy-detector positivity.**
  * It gives `t₄ ≤ 4`, hence the local theorem, which Herzog–Stergiou already have.
  * The next order, `Ψ‴(1)`, involves integrated `⟨TTTT⟩`. Our angular rules constrain `⟨TTTT⟩` only through light-ray
    data, so the reach is still local in `x − 1`.
* **Certified proofs.** These need a positive structure whose extremum is the free scalar. None is known for the
  `x → ∞` asymptotics.
* **The invisible-tail objection holds, and it is concrete.**
  * Adding a positive spectral tail `ε θ(E − E₀) e^{aE^{2/3}}` to the `S²` density of states changes `log Z(β)` at any
    fixed `β` by an amount that can be made arbitrarily small (`∝ ε e^{−βE₀}` times a growth factor). Yet it sets `c_S`
    for `a` above the theory's own coefficient.
  * Herzog–Stergiou exhibit the planar version. In the regular-boson CS-matter theory, the marginal sextic `λ₆` leaves
    planar `⟨TT⟩` and separated-point `⟨TTT⟩` unchanged but moves `c_S/C_T`: 49.03 versus 48.82 at `λ = 1/2`.
  * So **no argument using only `⟨TT⟩`, `⟨TTT⟩`, finitely many low-temperature or `T₀` moments, or the collider disc can
    work**, and that is every tool we have.
* **What global constraint could exclude the tail?** Candidates, none usable now:
  1. Thermal/spatial cycle exchange (`ε_vac = c_S/3`). This is exact, but it only re-expresses the problem as a bound on
     the NS Casimir energy, which is vacuum data of a *different* Hilbert space and not fixed by flat-space OPE data.
     Quantum energy inequalities do not help: on compact space they take the Casimir term as input [M].
  2. Genus-2 crossing (Simmons-Duffin–Xu). It is the only candidate with positivity in both channels. Turning it into an
     upper bound on the `S¹×S²` growth in units of `C_T` is a research programme of months to years, and not ours.
  3. Modular covariance on `S¹×S²`. None is known (Herzog–Stergiou say so).
  4. Positivity of the thermal spectral function (`ω ρ_β(ω,k) ≥ 0`). This is true positivity, but it fixes no upper
     bound on `ε/C_T`. Our shear sum rule is of this type and needs `ε` as input.

### 1.6 Weaker target: rigorous two-sided bounds on `c_S/C_T` for 3D Ising

* **Current state.** 51.08(13), from Monte Carlo Casimir amplitudes plus bootstrap `C_T`. The thermal bootstrap gives
  consistent but non-rigorous estimates.
* **Upper bound.** Rigour fails because thermal OPE coefficients `a_O` are not sign-definite, and the KMS condition
  involves an infinite unknown tail. Positivity exists only for spectral functions and for the spatial-channel
  decomposition. Neither bounds `a_T` (equivalently `c_S/C_T`) from above.
* **Lower bound.** A finite set of known operators says nothing about `β → 0`. Multi-twist families are not rigorously
  established with controlled degeneracy, and a generalised-free-field count overshoots (it gives `E^{3/4}` growth).
* **Verdict.** There is no rigorous two-sided bound within reach. Non-rigorous precision is Monte Carlo's domain. Not
  PRL.

### 1.7 Verdict on (f)

* **Flat conjecture.** PRL-level if proven: yes. Provable in weeks or months with our tools: no (≈0%). No mechanism from
  the CFT axioms is currently known that bounds the high-temperature entropy by `C_T` in 3D.
* **Disproof.** Leading-order large-N thermal-order models in d = 3 are safe (§1.4). Next-to-leading order in `3 < d < 4`
  is the sharpest open test, but it is outside our toolset.
* **Pointwise hyperbolic conjecture.** This is genuinely testable with our tools. Thermal-order CFTs are symmetric at
  `x = 1` (the vacuum) and ordered as `x → ∞` (flat space), so **there must be a symmetry-breaking transition at some
  `x_c > 1` on `S¹×H²`**.
  * If it is first order, latent heat violates `Ψ′ ≤ Ψ′_s` outright. Herzog–Stergiou name this loophole and say no
    example is known.
  * If it is second order, the mean-field heat-capacity jump must be compared with the deficit just below `x_c`.
  * A violation would also be the first Rényi phase transition at `q < 1` in a unitary CFT [I]. As far as I recall, the
    holographic Rényi transitions of Belin–Maloney–Matsuura occur on the low-temperature side, `q > 1` [M].
  * Recommended checkpoint: §4.

---

## 2. The other candidates

### (a) LHC / QCD energy correlators

* **Open?** No. The data-testable content is taken:
  * Ricci–Riembau–Son 2611.31587 bound multipoles of the three-point correlator from LEP one- and two-point data [A/S];
  * Belin et al. 2608.26246 bound contact terms from ALEPH data [S];
  * Riembau–Son 2512.16985 give spinning-source sum rules in QCD [A];
  * energy-conservation sum rules and ANEC commutativity are classic (Koloğlu–Kravchuk–Simmons-Duffin–Zhiboedov
    1904.05905; Korchemsky–Zhiboedov; Belin–Hofman–Mathys 1904.05892) [A/M].
* **Why our rules add nothing here.** They are 3D and CFT-specific. In event data `[E_m, E_{−m}] = 0` is automatic (the
  correlator is symmetric). Their content lies in the decomposition over conformal intermediate states, which does not
  exist in QCD.
* **Result.** None testable. **Tools:** n/a. **Time:** n/a. **Risk:** not novel.

### (a′) Stress-tensor-source energy correlators in 3D Ising and O(2): exact superconvergence rules plus bootstrap `⟨TTT⟩` (better variant)

* **Open?** Yes, but announced by others:
  * Ricci–Riembau–Son: 3D stress-tensor source, helicity states, "left for future work" for Ising;
  * Meçaj et al.: future direction "spinning-source EECs and superconvergence sum rules".
  * Neither uses commutator equalities. Ricci–Riembau–Son use positivity only.
* **Result.** With Ising `t₄ = 3.88933` (`|c|² = t₄²/16 = 0.94543`) [C from `n_B, n_F`], the `k = 2` rule forces
  `Σ_{J≥6} |a_{O,4,2}|² ≥ 0.945`. The Racah linear program then turns this, together with the known spectrum, into
  certified two-sided bounds on `λ_{TTC}`, the spin-6 couplings and the Ising stress-tensor EEC multipoles.
  * Physical message: the Ising graviton is about 97% free-boson-like, which forces near-conserved spin-6 weight.
* **Why non-specialists would care.** Only moderately: energy correlators have collider cachet, but the numbers are
  bootstrap data.
* **Tools.** 2 (helicity ladders, `R_k`, parity-weighted rule) and 3 (exact or interval linear-program certificates).
* **Time.** 2–4 weeks. **Risk.** Scoop (two groups); and the result is incremental (more a JHEP paper than a PRL).

### (b) 2+1D quantum critical points (disc predictions, `t₄` and `α_T` values, QMC)

* **Open?** Mostly no.
  * Guo–Li 2611.07394 [R] already give `t₄` (`O(N)`, GNY, QED₃ at `1/N`; Ising 3.8893 from bootstrap) and spell out the
    shear sum-rule weight `(3/4+t₄/16)P`. They flag Chern–Simons extensions with parity-odd structures as next, and
    Guo–Li–Ren 2611.24531 [A] continue the programme.
  * Our Ising number `S_E = 0.49654` (after the `ε` subtraction) is a one-line corollary.
  * Parity-odd `α_T` for QED₃–Chern–Simons at large `N_f` (models of quantum Hall and semion transitions) was not found,
    but it is the Guo–Li team's stated next step and is not our tool.
  * The CS-matter Kalmeyer–Laughlin transition exists on the fuzzy sphere (2507.19580, PRL) [A], but no fuzzy-sphere
    `⟨TTT⟩` exists.
  * No QMC stress-tensor (shear) data exist: the lattice stress tensor is the obstacle [I].
* **Result.** Placing real critical points on the disc. **Why care:** modest. **Tools:** 1. **Time:** days.
* **Risk.** Already done for the parity-even part; nothing is measurable.

### (c) Hall viscosity in experiments

* **Open?** It is open but disconnected.
  * Graphene Hall viscosity (Berdyugin et al. 2019 [M]) is measured in a magnetic field, in a non-conformal regime.
  * Chiral graviton Raman ("Evidence for chiral graviton modes in fractional quantum Hall liquids"; the title was found
    on INSPIRE, while the authors, Liang et al., Nature 2024, are [M]) measures chiral stress spectral weight in *gapped*
    fractional quantum Hall states. There the Golkar–Nguyen–Son sum rules apply (1309.2638 [A]), not ours.
  * Our `S_H` is a reactive (`Im h`) integral, not a chiral-absorption one.
* **Result.** None measurable. **Time:** n/a. **Risk:** no experimental handle.

### (d) Chaos and modular transport

* **Open?** Nothing sharp found.
  * At `T₀` the dynamics is Rindler, so chaos there is kinematic.
  * The `η/s = 1/4π` analogue on causal horizons is old: Chirco–Eling–Liberati 1005.0475 [A title]. Our
    `c_H(T₀) = 4πη_Y(T₀)` is its `H²` cousin.
  * Butterfly-velocity dependence on `t₄` is holographically known [M].
* **Result.** At most a remark. **Risk:** novelty.

### (e) Certified positivity settling a known conjecture

* **(e1) ANEC ⇒ the 6D a-theorem.**
  * Open: Hartman–Mathys–Stergiou 2611.23370 [R], §1: "Δa in six dimensions can be extracted by a light-ray sum rule, but
    we have not found a way to prove positivity". Dilaton amplitudes are also inconclusive (Elvang et al.).
  * Result if it worked: a proof of the 6D a-theorem. That is far above PRL threshold, and of broad interest.
  * Tools: 3 (finding a dual functional, i.e. adding null constraints until positivity holds) and partly 2.
  * Time: 1–2 weeks to reproduce their construction in 6D, then open-ended.
  * Risk: very high. The obstruction may be real, since `Δa` sits at an odd power `s³` of the dilaton amplitude.
    Hartman–Mathys–Stergiou are surely working on it. Success probability about 5% [I].
* **(e2) Bueno–Casini `C_T/F` bound, or the Kusuki–Tang / Herzog–Stergiou `c_S/F` bound.** These meet the same global
  obstruction as (f), plus `F` is non-local. Not tractable.
* **(e3) Modular-bootstrap gap bounds.** Not our toolset.

---

## 3. Ranking

Scores are my judgement, not measurements.

| rank | target | P(clean result in ≤ 4 weeks) | P(PRL given result) | tool fit | scoop risk |
|---|---|---|---|---|---|
| 1 | **(f-H)** refute or test the pointwise `H²` envelope via thermal order | 0.5 (a decisive answer); 0.2–0.3 (a violation) | 0.3–0.4 (a violation, as Letter or Comment); ≈0 (a confirmation) | good (`H²` kernels, gap equations, interval arithmetic) | low–medium |
| 2 | (a′) 3D Ising stress-tensor EEC with exact rules | 0.7 | 0.1–0.15 | excellent | **high** |
| 3 | (e1) 6D a-theorem positivity | 0.05 | 0.9 | partial | medium |
| 4 | (b) disc placement of real critical points | 0.9 | <0.05 | excellent | already done |
| 5 | (f) prove the flat bound | ≈0 | 1 | none | low |
| 6 | (d), (a)-LHC, (c) | — | — | poor | — |

Expected PRL value is about 0.08–0.1 for each of ranks 1–2 and about 0.05 for rank 3. Rank 1 wins because it decides in
a week, it is what the user wants, and its downside is cheap.

---

## 4. Recommended target and first-week plan: (f-H)

**Statement to test.** In the large-N d = 3 biconical CFT (`α = −1`, `g = 15360`; Smolkin–Yung), is
`Ψ′(x) > Ψ′_s(x)` for some `x > 1` on `S¹×H²` with vacuum Rényi boundary conditions? Equivalently: does the `Z₂`
thermal-order transition at `x_c > 1` carry a heat-capacity spike or latent heat exceeding the free scalar?

**Day 1–2.**
* Reproduce the Herzog–Stergiou `Ψ_s`, `Ψ_Dirac` and `Ψ_{O(N)}` curves, including their `O(N)` gap equation
  `0 = ∫dp[p tanh(πp)/√(p²+m²) coth(π√(p²+m²)/x) − 1]`, with our `H²` kernels (`code/hyp_identity/spectral_kernels.py`).
* Cross-check `Ψ″(1) = 3/4 + t₄/240` against our `η_Y` slope `−t₄/30` and the shared free-scalar exception. This also
  yields the citation paragraph the Letter needs.

**Day 2–5. Leading-N biconical saddle on `H²`.**
* Write the leading-order effective potential with the `H²` tadpole `𝓘(m,x)` in place of the flat one. Fix the
  normalisation by matching the flat limit `m = 2T arcsinh(½e^{2πρ/T})` (Appendix).
* Include the conformal-coupling term `ξR = −1/4R²` for `χ`, and keep `m(1) = 0, χ = 0` at `x = 1`.
* Solve for `(m(x), χ̄(x))`; find `x_c` and the order of the transition; compute `Ψ′(x)` and compare with `Ψ′_s`.

**Day 5–6.**
* If `Ψ′ > Ψ′_s` anywhere, certify the sign at leading order with interval arithmetic (tool 3).
* Estimate the size of the `1/N` corrections; `g` and `α` are fixed only up to `O(1/N)`.

**Day 7: go/no-go.**
* **Go** if there is a certified violation, or a first-order transition. Draft "Rényi phase transitions from thermal order
  violate the free-scalar envelope". Notify nobody until it is checked twice.
* **Kill** if there is no violation, or if the saddle turns out ill-posed (next item). Write a two-paragraph note and
  switch to (a′) on day 8. Its first week:
  1. express Ricci–Riembau–Son's `c_n, d_n, H₄` in our ladder variables;
  2. run the Racah linear program with Ising `t₄ = 3.88933` and the bootstrap spectrum `(ε, C, spin-6, …)`;
  3. certify a lower bound on the spin-≥6 weight.

**Main risk.** On `H²` the boundary carries a finite fraction of the volume. A homogeneous `χ` condensate is
non-normalisable, and the conformal coupling sits exactly at the `H²` marginal-stability point. So the Rényi-boundary-
condition saddle may be inhomogeneous, which makes the large-N gap equation position-dependent and much harder.

**Mitigation.** Day 2 should first settle whether the homogeneous Kubo–Martin–Schwinger state used by Herzog–Stergiou
(for which their `O(N)` checks are homogeneous) admits a homogeneous condensate. If it does not, the conjecture's own
definition of `Ψ` has to be sharpened. That is a publishable remark, but not a PRL.

---

## Appendix: the leading-order biconical computation (re-run with `python3 -I`)

```python
from mpmath import mp, mpf, asinh, exp, quad, zeta, pi, findroot, inf
mp.dps = 30
g = mpf(15360)
m = lambda r: 2*asinh(exp(2*pi*r)/2)            # T = 1, d = 3 LO gap solution (Smolkin-Yung)
def V(chi2):                                      # V_LO/N - V_free/N on phi = 0, rho = -chi2 (alpha = -1)
    return quad(lambda r: m(r)**2/2, [-inf, -chi2]) + g*chi2**3/720
c2 = findroot(lambda c: -m(-c)**2/2 + g*c**2/240, mpf('0.06'))   # d V / d chi2 = 0
fl = -zeta(3)/(2*pi)                              # free real scalar f/T^3
for lab, dv in [('symmetric', V(mpf(0))), ('ordered', V(c2))]:
    k = -(fl + dv)/(-fl); print(lab, k, 16*pi*zeta(3)*k)
```

**Output.**

| quantity | value |
|---|---|
| `χ̄²` | 0.0596309 |
| `m` | 0.674647 |
| symmetric point, `c_S/c_{S,free}` | 0.8 |
| symmetric point, `c_S/C_T` | 48.33758 |
| ordered point, `c_S/c_{S,free}` | 0.8799090 |
| ordered point, `c_S/C_T` | 53.16583 |
| `m(0)` (equals `2 log φ`) | 0.9624237 |

**Inputs assumed.** Smolkin–Yung's `V_LO` normalisation (validated by the 4/5 check above), and `C_T = N·C_{T,s}` at
leading order.
