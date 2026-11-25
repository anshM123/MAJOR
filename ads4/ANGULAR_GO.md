# ANGULAR_GO: does the angular route keep cubic sensitivity? Verdict, first-moment check, bandwidth and gap scaling

*25 Nov 2026, research agent (≈ 3.5 h). Unreviewed, not committed. No existing file was edited. `T0_CERT.md` and `code/t0cert/` were not
touched. New code is in `ads4/code/angular2/` (`ag_common.py`, `ag1_kernels.py`, `ag2_bandwidth.py`, `ag3_firstmoment.py`). Outputs go to
`code/angular2/cache/`, whose `.gitignore` ignores everything but itself. Read-only imports: `code/kb2a4/kb_celestial.py`;
`code/kernels/pipeline.py`, `kernels_lib.py` and `structures.py`, with `KERNELS_CACHE` redirected to `angular2/cache` so that nothing is
written in `code/kernels/`. Compute: about 10 CPU-min, at most 2 processes, `ulimit -v 3000000`.*

Grades: [P] proved or exact algebra; [D] derived; [D-sketch]; [DV] derived and numerically checked; [NE]; [M]; [I] open.
Inputs: ANGULAR_ROUTE §0–§6 (as committed during this session), LOWNU_SCOPE, K_B2A4, C1_FINAL.

---

## 0. Verdict

**Step 1: GO, but conditional and with a ceiling.**

**Cubic sensitivity survives with the correct signs.** The mechanism is the one the authors use in §3–§6, not the one the brief feared.
There is a ceiling, though: a bound obtained this way can scale at best as **τ ∝ 1/Δ_g**, never as Δ_g⁻⁴. The authors' 0.0224 at
Δ_g = 1000 sits at that ceiling.

**Why the plain SP2 form is a NO-GO [D].** Suppose every slot kernel is Δ-free with non-negative weights, which is what "positive harmonic
by harmonic" with SP2 alone would give. Then the certificate never sees Δ, so the bound it gives is the same for every gap.
* Its value is a pure number that must hold in every CFT where the sum rules hold.
* It therefore cannot fall with Δ_g. The best it can reach is a collider-type O(1) bound (§1.3).
* The gap enters only through the **Δ-dependence of the k = 2 first-moment transition**, i.e. T_Δ = (Δ−1)(Δ−2) in ANGULAR_ROUTE §3.
* I confirm that Δ-dependence independently at J = 0 (§4). It is quadratic in Δ, so it can buy at most T_Δ ~ Δ_g² in τ², i.e. τ ∝ 1/Δ_g.

**Why the sensitivity is not killed.** Three facts carry it.
* **The |ĝ₃|² light term occupies a heavy slot [P].** It is the T family's W×W term, which sits in harmonic |m| = 4, channel s = −1, the
  same slot as heavy states with J ≥ 2 and λ = −h.
  * In the k = 3 sum rule it has exactly the heavy kernel k⁽¹⁾(4) ∝ Q₂ (verified in §1.1).
  * The **k = 3 fourth-harmonic projection** (weight Q₂ × CDH density, which isolates |m| = 4 by orthogonality) therefore gives an exact
    sum rule: Σ_heavy (w₊ − w₋)_{|m|=4} = b′|c|².
  * Hence Σ_heavy w_{|m|=4} ≥ b′|c|². This holds whatever the sign of b′, since only |b′| enters [D, given the reduction lemma and the
    flip-entry caveat in §2.3].
* **The Einstein term is at harmonic 0, which no heavy state populates (momentum Ward identity).** So the Einstein term is unconstrained by
  harmonic positivity.
* **The global sign is not free [D].**
  * In an Einstein-like theory (τ ≈ 0) the per-ν sum rule gives Σ_heavy = −(Einstein light value) for any weight.
  * So any correct sign-definite heavy form must lie on the side opposite to the Einstein term.
  * If the computed heavy form came out on the same side, the form would be wrong; it is not a convention to choose.
  * ANGULAR_ROUTE §5–§6's calibration (forward normalisation ∫χ(R₋+R₊) = 15/4 > 0; Einstein −Ω = 139/4 − 48 log 2 > 0; −T_ΔZ₂/6 > 0)
    is consistent with this.

**The 303/5120 is checked [P, given K_B2A4's Θ].**
* With the authors' χ = −ν²(4ν²+9)(4ν²+49)/(32(4ν²+1)cosh πν):
  * the Einstein Lemma-R value ∫χ·[R₋D(ν+i) + R₊D(ν−i)] = 139/4 − 48 log 2 (30 digits);
  * the W×W value ∫χ·E_LR·t₄²Θ = −t₄²·303/81920 = −(303/5120)|c|², where |c|² = t₄²/16.
* It agrees with the h1 reading: ∫χP = −13635/8, and 13635/64 = 240²·303/81920.
* **The sign is right: the W×W term is opposite to the Einstein term. But it is right only because of the AdS (sub-flat) terms of N₁₀.**
  * The flat ν¹⁰ term alone gives +0.0458|c|², the wrong side.
  * With this narrow weight the full value is −0.0592|c|².
* **The term is irrelevant to the bound.** 0.059 is negligible next to the gap-enhanced 47266. The sensitivity lives entirely in the k = 3
  fourth-harmonic projection combined with the T_Δ-enhanced k = 2 heavy form.

**What this gives at Δ_g = 1000.** √(t₄²+α_T²) ≤ 0.0224, conditional (ANGULAR_ROUTE §6). It scales as ≈ 22.4/Δ_g (§3).
* This is about 200× worse than LOWNU's ≈ 1e-4.
* It is about 10⁶× worse than the conditional 2.9e-8.
* It would be **unconditional** in the sense that it needs no beyond-Regge remainder (LF3), **if** the items in §2 close.

---

## 1. Step 1 in detail

### 1.1 Slots and kernels [P, DV] (`ag1_kernels.py`)

**Kernels.**
* The celestial Fourier ratios k(m)/k(2) equal the CDH polynomial Q_{m−2}(ν²) with (d,ℓ) = (1,0), for **both** the j′ = 0 (k = 2-type)
  and the j′ = 1 (k = 3) kernels.
  * Checked for m = 3..6 at ν = 0.7, 1.9, 3.3.
  * Agreement is ≤ 2e-10 for j′ = 0 and ≤ 3e-12 for j′ = 1 (truncation-limited).
  * K_B2A4 had checked only m = 4.
* Parity: k⁽⁰⁾(−m) = k⁽⁰⁾(m) and k⁽¹⁾(−m) = −k⁽¹⁾(m), to 1e-28.

**Slot value.**
* Take a heavy family (Δ, J), helicity λ ∈ [−J, J], external h = ±2, harmonic m = λ − h, with |m| ≥ 2 by the Ward identity.
* The k = 3 dual satisfies d₋₋ = −d₊₊ and its kernel is odd. So the diagonal slot value is A₂(|m|; Δ) + s·A₃(|m|), with s = sgn(h·m).

**Both channel signs occur at every |m| ≥ 2.**
* s = +1 needs J ≥ |m|+2; s = −1 needs J ≥ |m|−2 (table in `cache/ag1_kernels.txt`).
* Hence **k = 3 alone, positive slot by slot, is forced to vanish**. This is the exact form of LOWNU's warning: the k = 3 part is
  channel-odd. Its slot values must be dominated by k = 2 slot by slot.

**The T family.**
* Einstein: λ = h, so m = 0. No heavy slot.
* W×W: λ = −h, so |m| = 4, s = −1. This is the heavy slot (4,−), populated by J ≥ 2.
* The O(t₄) cross terms vanish under EE (c₁ = 0, C1_FINAL).

### 1.2 The fourth-harmonic sum rule [D]

**The projector.** For k = 3 the slot kernels are Δ-free polynomials Q_{|m|−2}(ν²) (K_B2A4 §3: the action factorises exactly into
(ν) × (Δ)). The weight χ₃ = Q₂(ν²)·w_CDH^{(1,0)}(ν) is admissible, since it decays like e^{−πν}·poly. By orthogonality it annihilates every
diagonal slot except |m| = 4.

**The sum rule.** Its light side is only the W×W term. Einstein has no k = 3 part, since k⁽¹⁾(0) = 0. Hence:

    Σ_heavy Σ_{|m|=4} s·|c_{hλ}|² ‖Q₂‖² = −(light W×W)  ⇒  Σ_heavy w_{|m|=4} ≥ |b′|·|c|².

**Caveat [I].** For spinning heavy states the k = 3 dual's **flip entries** (h₄ = −h₂) can contribute 2-D harmonics. K_B2A4 only shows that
they vanish for J = 0, where the flip part is 1↔3 symmetric. Either show that the flip part's Q₂-projection vanishes for all J, or bound
it by the same weights.

### 1.3 Why Δ-free (SP2-type) positivity cannot give a gap-dependent bound [D]

**Setup.** Suppose every heavy slot value is Δ-independent and the weights are non-negative, i.e. ω[O] = Σ_σ w_σ(O) v_σ(χ).

**Consequences.**
* A χ with all v_σ ≥ 0 is positive on **every** non-T family of **every** unitary CFT₃ where the k = 2, 3 sum rules hold. No gap and no
  large N enter.
* So the inequality ω[T] = E(χ) + τ²W(χ) ≤ 0 is universal.
* E and W are Δ_g-free functions of ν, so the bound τ² ≤ −E/W is Δ_g-independent.
* It must also be compatible with every CFT in the class, so at best it is a collider-type O(1) statement.

**Consistency check.** The per-ν sum rule then fixes the heavy content of each slot universally. For example, the heavy (4,±) content is
pinned by the Einstein term's Q-expansion together with τ². That is the same statement in other words.

**What follows.** **The gap must come from somewhere Δ enters the slot kernel.** The only place is the k = 2 first-moment transition. The
flat limit already shows it: the J = 0 heavy k = 2 action ∝ ∫ψ₂(2 − ν²/Δ²), so it is not Δ-free.

### 1.4 Bandwidth: the slot positivity confines the weight to ν = O(1) [D-sketch + NE]

**Alternation.**
* The largest zero of Q_j(ν²) lies at ν ≈ 1.9 j for (1,0), and between 1.4 j and 1.9 j for (0,0), for j ≤ 64.
* Beyond it, sign Q_j = (−1)^j (`ag1`, all j checked).
* At large ν the celestial kernel localises on antipodal detectors (φ₃₁ = π), so harmonic m carries (−1)^m.

**Consequence.**
* Positivity of each smeared slot separately cannot hold for both parities of m when the weight has mass at ν ≫ 1.
* Per-state positivity can, because the heavy λ-distribution resums the alternating series. Slot-wise positivity throws that away.
* So any certificate that needs a single sign of the smeared harmonics for all m forces bandwidth O(1).

**Test (`ag2_bandwidth.py`).** Stretch the authors' density to w_s = 8ν²/(s³cosh(πν/s)).
* z_m^raw keeps one sign for 2 ≤ m ≤ 39 only at **s = 1**.
* Already at s = 1.5 the first few signs alternate. At s ≥ 2, z_m alternates as (−1)^m at every m checked.
* The Einstein value also changes sign (s = 1.5 to 3), and the W×W sign is not stable.

**The authors' weight is at the edge.** Their weight (e-folding 1/π in ν) is essentially the widest that keeps the telescoping, so the
route's bandwidth is ν ≈ 1–3.

---

## 2. Step 2: transition matrices, what is proved, and the all-J gap

### 2.1 The scalar first moment: reproduced [DV, exact symbolic]

`ag3_firstmoment.py` runs `pipeline.amplitude` (read-only) on the J = 0 conserved ⟨T T O_Δ⟩ structures, one in each parity, symbolic in
Δ = 2s.

**Rule used.** The first moment uses the analytic rule I₁^mid(a,b) = (a−b)/(2(a+b−2))·(2z·x)·I₀ [P]:
* per monomial, a = e14 − j, b = e13 − i;
* a − b = 2s−5+i−j and a+b−2 = esum−2−i−j;
* multiplying by (2z·x) means l → l+1 with cf ×2, because the pipeline's pow2(2k−l) already converts (z·x)^l.

**Resolving the ambiguity in the literal rule.** The coordinator's rule (cf ×= (a−b)/(2(a+b−2)), l+1) is **half** the α-moment.
* It is the **retarded-time** moment, because u = α/2.
* Reason: the detector ray crosses x₃ = 0's light cone at α = 0 (X₁₃ = 0) and x₄'s at α = −2z·x (X₁₄ = 0), where u = n̄·x = −z·x.

**Results (u-moment, T helicity +2, n = (1,0)), identical in both parities:**

| quantity | result |
|---|---|
| H (zeroth moment, prefactor G stripped) | even: 3s(s−4)(s−3)(s+1)(2s+3)(2s+5)/256; odd: −3(s+2)(s+3)(2s−7)(2s+1)/128 |
| **M_mid/H (u-moment, midpoint-centred)** | **i(Δ²−7)/6** in both parities. **Authors' claim confirmed.** (α-moment: i(Δ²−7)/3) |
| origin-centred u-moment / H (T at 0, O at x) | i(Δ²−7)/6 − i(Δ+3)/2 = i(Δ²−3Δ−16)/6 |
| authors' normalised-state value | i(Δ²−7)/6 − i(Δ−3)/2 = **i(Δ−1)(Δ−2)/6** |

**The normalisation correction.** The step to −i(Δ−3)/2 is the change from unnormalised to normalised momentum states, ½n̄·∂ log(N_O/N_T)
with N ∝ (p²)^{Δ−3/2}. It vanishes at Δ_O = Δ_T = 3, as it must for C1_FINAL's symmetric T→T form. I confirm its structure and that it
vanishes at Δ = 3. I did **not** independently fix its sign (the n̄·p orientation), so the final i(Δ−1)(Δ−2)/6 is [D-sketch] here.

**Agreement with §3.** It matches the J = 0 entry of the authors' A_{Δ,J} = T_Δ + C_J − 2m² (= T_Δ), with M = (i/6)AH. That supports the
coefficient −6 in [𝖪_n, E] = −6E⁽¹⁾ at J = 0 only.

### 2.2 All-J statements: how to prove them and how long

**(a) Tridiagonal first-moment identity M = (i/6)A_{Δ,J}H, all J [I → D, ≈ 1 day].**
* Route: E⁽¹⁾_{OT} = −(1/6)(𝖪_{O,n}E_{OT} − E_{OT}𝖪_{T,n}).
* 𝖪_n on normalised rest-frame helicity states of a spin-J primary acts with:
  * a scalar part (Δ-dependent, giving T_Δ = (Δ−1)(Δ−2) = Casimir-type);
  * a spin part S_{n̄ν}, which mixes m → m ± 1 with the SO(3)-type coefficients √(C_J − m(m±1)).
* The cross term √(T_Δ − m(m±1)) needs the boost part of 𝖪 acting on the momentum-dependent frame.
* This is a finite representation-theory computation in the Gillioz/KERNELS_D3 §3.3 normalisation, now all-J by ANGULAR_ROUTE §4.
* Plan: derive symbolically, then check with `ag3`-type pipeline runs at J = 1..6 in both parities. This needs `amplitude` with the
  first-moment rule on the conserved spinning structures, about 1–2 CPU-h.

**(b) Ladder ratios a_{k+1}/a_k = ((k+4)/(k+1))√(X_k/Y_k) on b = 0 and b = 4, all J [I, 1–3 days, medium risk].**
* The number of independent conserved ⟨T T O_J⟩ structures is small (2 parity-even + 1 odd for even J ≥ 4). So each ladder is fixed
  by one amplitude, and the ratio is a statement about the m-dependence of the light transform of each structure.
* Two routes:
  * **(i) Closed form.** The pipeline's monomial formula is a finite sum of Γ-ratios. For symbolic J it becomes a terminating
    hypergeometric sum in m. Prove the two-term recursion by creative telescoping (Zeilberger/WZ). First generate data to J ≈ 30
    (pipeline, ≈ 1 CPU-h) to fix X_k, Y_k exactly.
  * **(ii) Ward identity.** Apply the K/boost Ward identity of ⟨O|[𝖪_{n′}, E(n)]|T⟩ for n′ ≠ n. Its rest-frame form gives a first-order
    difference equation in m with the same √(T − m(m±1)) coefficients as A. This is shorter if it works, but I have not checked that it
    closes on one ladder.
* Thursday-feasible: route (i) through J ≈ 30 as [DV], and possibly the WZ proof. A clean all-J [P] by Thursday is about 40 %.

**(c) Physical-matching gap. This is the real remaining risk.**
1. **The heavy action equals the rest-frame quadratic form.** This rests on K_B2A4's reduction lemma (§2 there, [D]). It was tested end to
   end only on the T family (Einstein quartic and flat top coefficient). For heavy families it additionally needs:
   * the origin/midpoint/normalisation bookkeeping of §2.1, with the sign fixed;
   * the EE (2,4) contraction being diagonal through first order. That is [P] and is fine.
2. **The k = 3 flip entries** for spinning heavy O in the fourth-harmonic projection (§1.2 caveat).
3. **Common normalisation of |c_{hλ}|²** between the k = 2 form (z, y, T_Δ) and the k = 3 projection (b′). The constant 47266 combines
   the two. A factor error here moves τ by its square root.
4. **Swap and convergence (H_K)** for χ ~ ν⁴/cosh πν. High harmonics grow at most like e^{πν/2}·poly (CDH asymptotics), so this looks
   fine, but it is not written.
5. **Double traces and contacts at this weight.** DZ holds at leading 1/c_T. ν = i/2 is excluded (its finite-part budget is negative,
   ANGULAR_ROUTE §6).
6. **Parity-odd α_T.** The W×W phase invariance is [D-sketch] (K_B2A4 §0).
7. **Large J at small twist [M-risk].** A ⪰ (√T−√C)² degenerates as Δ − J → 2. The positivity then rests on the actual ladder vector
   rather than the floor (ANGULAR_ROUTE §5: arbitrary-amplitude positivity fails at J = 10⁶). So (b) is load-bearing, not cosmetic.

**Sketch of the theorem.** Assume the conserved-structure ladder identities and the tridiagonal identity hold for all J, and that items
(c)1–6 hold. Then for every heavy family with Δ ≥ Δ_g:
* the k = 2 χ-functional is ≥ κ T_{Δ_g} w_{|m|=4} + (non-negative);
* the k = 3 fourth-harmonic projection forces Σ w_{|m|=4} ≥ b′|c|²;
* so (139/4 − 48 log 2) ≥ (κT_{Δ_g}b′ + 303/5120)|c|².

---

## 3. Step 3: the achievable bound and the bandwidth question

**Bound.**
* τ ≤ 4√((139/4 − 48 log 2)/(47266 + 0.059)) = **0.02237** at Δ_g = 1000. This uses the authors' uniform constant; I did not recompute
  47266.
* Since 47266 ∝ T_{Δ_g} ≈ Δ_g², **τ ≈ 22.4/Δ_g**: about 0.075 at Δ_g = 300 and 0.011 at 2000.

**Answer to the bandwidth question (one paragraph).** The angular route does **not** need LOWNU's e^{πν}-weighted de-shift tail bound.
It never moves a contour for heavy states. The k = 3 slot kernels are polynomials and the k = 2 ones are polynomials times T_Δ-linear
coefficients, so the only analytic requirement is a convergent swap against kernels that grow at most like e^{πν/2}, which a 1/cosh πν
weight satisfies. It has a **harsher bandwidth cap of its own**, however. Its positivity certificate needs a single sign of the smeared
harmonics z_m for all m. The CDH kernels alternate as (−1)^m for ν ≳ 1.9m (antipodal localisation), so that sign survives only for weights
of width ν = O(1): the authors' 1/cosh πν works, and the same weight stretched by 1.5 already fails (`ag2`). The bound therefore cannot
gain the (aΛ)⁸ light-side lever that drives LOWNU. All of its gap dependence comes from T_Δ in the first moment, giving τ ≈ 22/Δ_g ≈ 0.022
at Δ_g = 1000. Recovering Δ_g⁻⁴ would need bandwidth ∝ Δ_g, and that needs per-state cancellations between harmonics, which is the
beyond-Regge remainder problem again in a new basis.

**Plan to Thursday.**
1. **½ day.** Pipeline data for first-moment and zeroth-moment ladders, J = 1..30, both parities (`ag3` extended to spinning conserved
   structures; read-only pipeline). Check the tridiagonal M = (i/6)AH and the b = 0/4 ladder ratios exactly per J. Fix the sign of the
   normalisation shift (§2.1).
2. **1 day.** All-J tridiagonal identity from 𝖪's action (§2.2(a)).
3. **1–2 days, parallel.** WZ proof of the ladder ratios (§2.2(b)(i)).
4. **½ day.** Matching items (c)2–3: the k = 3 flip-entry projection for J = 2, 4 from the pipeline, and the common |c|² normalisation
   behind 47266. Recompute 47266 independently.
5. **Expected result.** A theorem-grade τ ≲ 0.022 at Δ_g = 1000, scaling as 1/Δ_g. It is unconditional on LF3 but much weaker than the
   conditional results. Present it as a complementary unconditional statement, not a replacement.

---

## 4. Reproduce (`ads4/code/angular2/`; `ulimit -v 3000000; export OMP_NUM_THREADS=1`)

```
python3 -B ag1_kernels.py       # ~2 min  -> cache/ag1_kernels.txt   (Fourier ratios = CDH, parity, slots, largest zeros)
python3 -B ag2_bandwidth.py     # ~3 min  -> cache/ag2_bandwidth.txt (z_m sign vs stretch; light side vs stretch)
python3 -B ag3_firstmoment.py   # ~1 min  -> cache/ag3_firstmoment.txt (J=0 first moment, both parities, symbolic Delta)
```

The light-side numbers of §0 (139/4 − 48 log 2, −303/81920·t₄², −13635/8) are reproduced in `ag2` at s = 1 and by the formulas in §0.
