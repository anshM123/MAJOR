# INDEP_CHECK3: independent check of `c₁ = 0` (C1_FINAL) and of the Batch-5 fold at Δ_g = 1000 (BATCH5_FOLD)

*25 Nov 2026, independent checker. Unreviewed, not committed. No existing file was edited. New code: `ads4/code/indep3/`
(`i3_affine.py`, `i3_freefield.py`, `i3_basis.py`, `i3_basis2.py`, `i3_diag.py`, `i3_mkset.py`); outputs in `ads4/code/indep3/cache/`
(git-ignored by its own `.gitignore`). Compute ≈ 0.6 CPU-h, ≤ 2 processes, one thread, `ulimit -v 3000000`.*

**Independence statement.** C1_FINAL's code (`code/c1final/`) was **not read** at any point; the free-field engine, the
phase-space bookkeeping, the frames and the end-to-end contraction were written from scratch (mpmath, 30 digits). For Batch 5, the
basis repair was re-derived and fitted with my own code before I looked at any `b5fold` output; `b5fold/cache/b5_basis.txt` and the
mapped point sets were read **afterwards**, as data, for the comparison. To reproduce the downstream values/budgets/L, I then **ran
BATCH5's own chain** (`b5_state.py`, `b5_window.py`, `b5_odd.py`, `b5_sum.py`) in a sandbox copy under `indep3/cache/sbx/` (sibling
directories symlinked read-only, write guard re-rooted into the sandbox), with the J = 18 even records replaced by **my**
independently mapped ones. That part is a re-run, not an independent re-implementation; I read `b5_guard.py`, `lane.sh` and the
env-var lines of `b5_core.py`/`b5_state.py` only to set the sandbox up. Physics libraries used read-only: `ext/point.py` →
`spinfunc/hc_points.py` (local anchor nodes) and `spinfunc/l5_spin.conserved` (local structure bases).

---

## 0. Table

| claim | independent result | agrees? | grade after check |
|---|---|---|---|
| **1(i)** `f⁺_{2;2} = v′−u′` (CLS-D eq. `subtraction_factor_evenjprime` at `k = 2, J_e = 2`, read from source) is affine in the retarded times at `x₀ = ∞`, closed form as in C1_FINAL §2.1 | sympy from the definitions: leading `r→∞` coefficient minus C1_FINAL's closed form `= 0` exactly; `∂²f/∂u₁² = ∂²f/∂u₃² = ∂²f/∂u₁∂u₃ = 0`; finite-`r` error `2.84·10⁻³, 2.85·10⁻⁶, 2.85·10⁻⁹` at `r = 10³, 10⁶, 10⁹` (exact rationals) | **yes** | [P] |
| **1(ii)** `E⁽¹⁾ = ½{E, −i n̄·D}` on the `T` family for both free fields; Einstein part diagonal, `t₄` part pure flip, `M̂_W = 0` | own free boson + free Dirac engine (two-particle wave functions, `E⁽¹⁾` as `ω(−i/2)(a†∂_ω a − h.c.)`, derivative taken numerically on the wave function — the decomposition is **not** assumed). Ward `∫dφ n̄^μĤ = P^μ` to `3·10⁻³⁰`; `Ĥ_B,diag = Ĥ_F,diag`, `Ĥ_B,flip = −Ĥ_F,flip`, `|flip|/diag = 1` to `3·10⁻³⁰` (HM endpoints `t₄ = ±4`), flip phase `e^{±4iφ}`; `M_flip = 0` exactly, `M_B,diag − M_F,diag ≤ 2·10⁻²⁵`, `M = −½{Γ,Ĥ}` to `3·10⁻²⁵` (minus sign = my convention `ψ_h = ε_h*·t`) | **yes** | [D + DV] |
| 1(ii) "free B, F span parity-even `⟨TTT⟩` at fixed `c_T`" | standard (2 parity-even structures in d = 3); normalised matrices are affine in `t₄`, Einstein = c_T-weighted midpoint | yes | [Q] |
| 1(ii) parity-odd `α_T`: no diagonal element of `Ĥ`, `M̂` (rotation + momentum/dilatation Ward identities) | argument re-derived: rest-frame diagonal elements are φ-independent by rotation covariance (boost frames are rotation covariant), their φ-integrals are fixed by `c_T`; Wigner-rotation terms are diagonal × (zero diagonal). No numerical check possible (no free realisation) | yes | [D] |
| **1(iii)** EE `(2,4)` contraction is the identity, diagonal through first order | sympy: `G(X₀ = I) − x₂₄²(z₂·z₄) = 0`; `G(−X₄,−Z₄) = G(X₄,Z₄)`; finite `x₀` → ratio → 1. Physical reading of the `x₄⁺` bra checked (§1.3): for integer-Δ `T`, `⟨0|T(x⁺, pushed-forward z) = −⟨0|T(x,z)` as bras, so the pairing is the ordinary helicity-diagonal one | **yes** | [P] (+ [D] for the `x₄⁺` reading) |
| **Route B numerical zero** (`≤ 1.75·10⁻¹¹`) with `O(1)` controls | pointwise in `P`, η-contracted over all `(μν)` polarisation pairs, 8 random timelike `P` × 3 structures (`E₃E₁`, `E₃E₁⁽¹⁾`, `E₃⁽¹⁾E₁`), random complex wave packets: **`|X|/|I[E,E]| ≤ 6.9·10⁻²⁴`**, each cross term separately `≤ 3.9·10⁻²⁴`; controls: `W×W` = `0.08 – 2.4 × I[E,E]`, reflected `(2,4)` contraction `0.07 – 3.9` | **yes** (zero reproduced, 13 orders tighter) | [DV] |
| **`c₁ = 0`** (whole `O(t₄)` light term of `Ĉ′_EE`) | the three local links and the zero of the `T`-block cross term are reproduced. The conclusion additionally needs two **inherited** links that I did not (and could not) re-derive: (DZ) only the `T` block survives at `O(1/c_T)` — true for parity-even `[TT]` because `f⁺_{2;2}` is a polynomial (no `1/(u′v′)` shift of the even-twist zeros), but for the **parity-odd `[TT]` (odd twist) H0 §15.1 grades the double zeros [D-sketch]**; and `Ω_OPE = Lemma R` on tree terms (LIGHTSIDE_DIST, [D]). Route B is the same computation as Route A's pointwise statement, so it is not an independent route to `c₁` | **partly** (not refuted) | **[D] given (DZ) and LIGHTSIDE_DIST; the parity-odd-`[TT]` part of (DZ) is [D-sketch]** — the same caveat as the Einstein term |
| **2** even-J basis repair (BATCH5 §1.2), J = 18 | my constrained fit (one O(2) angle + reflection bit per basis, colmax normalisation, re-derived): `θ_G = 3.117435`, refl 0; `θ_H = 4.711468`, refl 1 — **identical** to b5_basis; reference: same up to the trivial `(A,B) → (−A,−B)` gauge. Residual `9.1·10⁻⁸` (state) / `7.9·10⁻⁹` (ref); C3H at the anchor with the EE map `5·10⁻⁹`; **my own unfitted node `x = 0.6`: `2.6·10⁻⁷` (EE), `1.5·10⁻⁷` (C3H)** against MC `3.6·10⁻⁴`. All 108 + 108 mapped records = BATCH5's set `fix` to `≤ 1.0·10⁻¹¹` | **yes** | [DV] |
| 2 basis repair, J = 26 (extra) | same: `≤ 9.8·10⁻¹²` vs `fix`; C3H at anchor `2.7·10⁻⁷` | yes | [DV] |
| **2** J = 18 values/budgets/L (twist edge + 1.02/1.03/1.05) | BATCH5 chain re-run on **my** mapped records: _RERUN18_ | _AGREE18_ | [NE±] |
| **2** J = 23 (twist edge + 4 masses) | re-run of `b5_odd.py 23 main`: _RERUN23_ | _AGREE23_ | [NE±] |
| **2** J = 25 (twist edge + 3 masses) | re-run of `b5_odd.py 25 main`: _RERUN25_ | _AGREE25_ | [NE±] |
| 2 failed-node sets | from raw replicas (`err = √(Σ|R−R̄|²/(7·8))`, re-derived and checked on raw records): J = 18 fails **x = 0.9447 and 0.95** (same on my mapping and on `fix`, and the chain drops both); J = 25 rz flank 0.9453 fails; J = 23 batch-5 nodes none | yes; doc omission (§2.4) | — |
| **2** "25 of 31" | recounted from the table: 0, 2, 4–21, 23, 25, 27, 29, 31 = 25; 33 spins J ≤ 32 minus J = 1, 3 = 31. Stricter: 23 / 18 / 16 — all reproduced | **yes** | — |

**Bottom line.**
* **Claim 1 (`c₁ = 0`): partly reproduced, not refuted.** F₂ (affine weight, helicity pattern of `E⁽¹⁾`, diagonal EE contraction)
  and the pointwise zero are reproduced independently, the zero to `7·10⁻²⁴`. What is not re-derived here, and carries the grade, is
  inherited: (DZ) — whose parity-odd-double-trace part is only [D-sketch] in H0 §15.1 — and `Ω_OPE = Lemma R`. C1_FINAL lists both
  in §8; I recommend that the grade line say "[D] given DZ (parity-odd `[TT]` zeros [D-sketch], H0 §15.1)" rather than an unqualified
  [D], and that "Route B" not be presented as an independent confirmation (it is the numerical form of F₂ in the same prescription).
* **Claim 2 (Batch 5 fold): reproduced** for J = 18 (with my independent basis repair), J = 23, J = 25, and the count. One minor
  documentation slip (§2.4).

---

## 1. Claim 1: methods

### 1.1 Step (i): the weight (`i3_affine.py`)
`x_{1,3} = (u_{1,3} + r, r n_{1,3})`, `u′ = x₁₂²x₃₄²/(x₁₃²x₂₄²)`, `v′ = x₁₄²x₂₃²/(x₁₃²x₂₄²)`, mostly-plus. After eliminating
`|n| = 1`, the `r³` terms of `x₁₄²x₂₃² − x₁₂²x₃₄²` cancel; the ratio of the `r²` coefficients gives
`f(∞) = 2[(u₁ − n̄₁·x₄)(u₃ − n̄₃·x₂) − (u₁ − n̄₁·x₂)(u₃ − n̄₃·x₄)]/((1 − n₁·n₃)x₂₄²)`, whose `u₁u₃` term cancels identically. It equals
C1_FINAL's closed form symbolically. I also confirmed from the CLS-D source that for `J_e = 2, k = 2` the subtraction factor is
`(v′−u′)/(u′v′)^0`, i.e. a polynomial — this is also what makes the parity-even double-trace zeros survive (§1.5).

### 1.2 Step (ii): free-field engine (`i3_freefield.py`)
* `t_B = k₁k₂ + k₂k₁ − η k₁·k₂ − ¼(PP − ηP²)` (improvement fixed by tracelessness), `t_F = ū(k₁)γ_(μ(k₁−k₂)_ν)v(k₂)` with
  `u(k) = √ω (sin φ/2, cos φ/2)`; conservation, tracelessness, symmetry `≤ 4·10⁻³¹`.
* `T`-family basis `ψ_h = ε_h*(P)·t`, `ε_h = e_h e_h`, `e_± = Λ_boost(P)(0,1,±i)/√2`. Norm `n(P) = ∫dφ₁|ψ|²/(−P·n̄₁)`; Lorentz
  invariance `n ∝ m³` checked by brute-force boosted quadrature to `4·10⁻³¹`.
* Derived phase-space bookkeeping: `⟨Ψ′|E(n)|Ψ⟩ = ∫d³P Σ_i [ω/(−P·n̄)] Ψ′*Ψ|_{k_i = ωn̄}` with `ω = P²/(2P·n̄)`, the Jacobian
  `d³P = (1 − k̂₂·n̂)dω d²k₂` giving `−P·n̄`. `E⁽¹⁾` is the same with `(−i/2)(Ψ′*∂_ωΨ − ∂_ωΨ′*Ψ)`, `∂_ω` at fixed spectator momentum,
  4th-order central differences at 30 digits. The improvement term of the scalar enters `T_uu` at `ℐ⁺` as `∂_u²(φ²)`, so it drops
  from **both** `∫du` and `∫du u` (it is a second, not first, `u`-derivative — the note's wording "total u-derivatives" is slightly
  loose but the conclusion holds). Real `ω`-dependent normalisations cancel in the antisymmetric form.
* Mechanism re-derived analytically: in d = 3, `e_+·e_+ = 0` forces `∂e_+ ∈ span(e_+, P)`, so `∂ψ_h = (∓iΓ_h + 1/ω)ψ_h` (the
  longitudinal piece dies by `P^μt_μν = 0`; for the fermion `ε·(k₁−k₂) = 2ε·k₁` on the transverse space, giving `1/(2ω) + 1/(2ω)`),
  hence `M_{h′h} ∝ (Γ_h + Γ_{h′})Ĥ_{h′h}` with `Γ₋ = −Γ₊` — no flip part. Numbers in the table.

### 1.3 Step (iii) and the `x₄⁺` bra
* Embedding conventions of my own (`X = (1, x², x)`, `Z = (0, 2x·z, z)`): `G := H₂₄ − 2V_{2,04}V_{4,02} = x₂₄²(z₂·z₄)` at `X₀ = I`.
* Adversarial point examined: in d = 3 the map `x ↦ x⁺` is orientation-reversing on the cylinder, so one might fear an improper
  reflection of the bra's polarisation (that would turn the pairing into the flip-sensitive control, which is `O(1)`). It does not:
  for `T` (Δ = 3, integer), Wightman functions are single-valued in the embedding variables, `H`, `V` are invariant under
  `(X₄,Z₄) → −(X₄,Z₄)` and the X-degree gives an overall sign, so `⟨0|T(x⁺, pushed-forward z)` = `−⟨0|T(x, z)` as bras. The
  `(z₂·z₄)²` pairing is therefore the ordinary η-contraction of the two states, which is helicity-diagonal at coincident momentum
  (`ε*_{h₄}·ε_{h₂} = δ`) and diagonal at first order (`ε*·∂ε = iΓδ`).

### 1.4 End-to-end zero (Route B), my implementation
Bra/ket are the `(μν)` components of `T` smeared with independent random complex Gaussian packets, `F_h = ε_{h,μν}(P)g(P)`, summed
with weights `η^{μμ}η^{νν}` (exact metric contraction — a real orthonormal basis does not exist in Lorentzian signature).
`E₃E₁⁽¹⁾`: bra replaced by `P_T E₃Ψ₄` (coefficients `Ĥ₃F₄`), then the direct `E⁽¹⁾` density; mirror for `E₃⁽¹⁾E₁`. Structures:
`E = (B+F)/2`, `W = (B−F)/2` (normalised). Control "reflected": the bra components replaced by `R ε R`, `R = diag(1,1,−1)`.
Results in the table (`cache/i3_freefield.txt`, rows in `i3_freefield_rows.json`).

### 1.5 Adversarial points (claim 1)
1. **Is the energy-channel `T` block all of the `O(ĝ₃)` light term?** Only through (DZ). For parity-even `[TT]` (even twist) the
   commutator factor `2 sin²(π(τ−10)/2)` has double zeros and the polynomial `f⁺_{2;2}` (checked in the CLS-D source) does not shift
   them, so their `O(1/c_T)` data drop out term-wise. I first worried that LIGHTSIDE_DIST's borderline toy (σ⊗T_free: a
   `u`-exchange with only "double traces" in the s-channel has a non-zero `k = 2` term-wise value) contradicts this; it does not
   refute it, because that product theory contains odd-twist-shifted operators (`[σσ]_n ⊗` twist-1 free currents) that sit off the
   zeros. But the same mechanism is exactly the **parity-odd `[TT]` caveat** (odd twist, `sin² = 1` unless the parity-odd structures
   shift the phase), which H0 §15.1 grades [D-sketch]. At `O(ĝ₃)` the parity-odd double traces do receive `O(1/c_T)` data, so this
   link is genuinely load-bearing for `c₁` (as it is for the Einstein term).
2. **De-shift.** The de-shifted functional is a linear `ν`/celestial integral of `Ω_k` at fixed `(a,b)` (LIGHTSIDE_DIST §1.1); the
   `O(t₄)` term vanishes pointwise for every `ν`, `(n₁,n₃)` and every scalar kernel (my zero is pointwise in `P` for arbitrary
   packets, i.e. for any `(x₂,x₄)` kernel and its first derivatives), so no constant can be reintroduced. The `k = 4, 5` de-shift
   pieces of `Ω^Q_up` are not part of `c₁`'s definition (C1_COMPUTE §1.1).
3. **Parity-odd `⟨TTT⟩`.** The `O(α_T)` cross term is (Einstein diagonal) × (no diagonal) = flip, killed by the diagonal pairing.
   Relies on the Ward argument only (no numerics possible); [D] is fair.
4. **Affinity matters?** Only first-order derivatives are needed because `f` has no `u₁u₃` term; reproduced.
5. **Not addressed (inherited):** `Ω_OPE = Lemma R` at the borderline. C1_VALUE's non-zero Euclidean `b₊ = 7/300` remnant is
   compatible with `c₁ = 0` only through that identity; no Lemma-R evaluation of the remnant exists. This is the one place a non-zero
   `c₁` could re-enter (as C1_FINAL §5/§8 says).

---

## 2. Claim 2: methods

### 2.1 Basis model (re-derived)
`E[a,c]` is linear in `sol(Δ, J)` (index `a`) and `sol(3−Δ, J)` (index `c`); `conserved` returns `U/colmax|U|`, `U` a real
orthonormal null-space basis fixed only up to `Q ∈ O(2)`. Hence `E_ext = M_Gᵀ E_loc M_H`, `M(Q) = pinv(sol_loc)·(U_locQ)/colmax|U_locQ|`.
**Finding:** a fully general `GL(2)×GL(2)` fit on one anchor node is **not identifiable** (`i3_basis.py`: residual `9·10⁻⁸`, but
the solution family spreads by `O(1)` across starts) — BATCH5's "a fully general 2×2 G, H fit gives the same residual (so the O(2)
form is right)" is therefore weak evidence; the real validation is the unfitted node, which I reproduced independently (`x = 0.6`).
The constrained model has 2 continuous parameters per (J, Δ) and is well determined.

### 2.2 Data
Local anchor nodes computed here with `ext/point.py` (unchanged physics, ≈ 100 s each): J = 18 at `ν = 304.899492` (x = 0.3) and
`609.798983` (x = 0.6), J = 26 at `307.247164`, each with the reference Δ = 2·10⁴ (`cache/loc/`). External records: all
`ext/unpacked/out/N5_J{18,26}_*`.

### 2.3 Downstream (re-run)
_DOWNSTREAM_

### 2.4 Discrepancies (claim 2)
1. **Documentation omission, no effect.** BATCH5 §1.1 lists the failed nodes as "x = 0.95 on every even row, and GK nodes 0.9439
   (J = 19), …". On the mapped records (theirs and mine) J = 18 also fails its GK node **x = 0.9447** (and the chain log shows J = 20
   failing 0.9432). The chain drops them correctly (`b5_state` log: failed `[0.944695, 0.95]`); only the prose list is incomplete.
2. **Weak identifiability argument** (§2.1) — the conclusion stands because the unfitted-node test passes.
3. **Reference-map gauge.** My reference angles differ from b5_basis's by `π` (`(A,B) → (−A,−B)`); `AᵀEB` is unchanged. Harmless.
4. Table columns: for the hs_odd states the "new" column is full λ and the cache's `lam` is λ_c; the offsets (J = 23: 0.022,
   J = 25: 0.030) are consistent between twist edge and masses, so the table is internally consistent.

---

## 3. Reproduce
```
cd ads4/code/indep3        # ulimit -v 3000000; export OMP_NUM_THREADS=1
python3 -B i3_affine.py                          # §1.1, §1.3  (~1 min)        -> cache/i3_affine.txt
python3 -B i3_freefield.py 8                     # §1.2, §1.4  (~10 min)       -> cache/i3_freefield.txt
python3 -B ../../ext/point.py n1 18 1018 304.899492 4242 $PWD/cache/loc/J18_x03   # and 609.798983 -> J18_x06; J=26 307.247164 -> J26_x03
python3 -B i3_basis2.py 18 $PWD/cache/loc/J18_x03 $PWD/cache/loc/J18_x06            # §2.1 (I3_THEIRS selects the comparison set)
python3 -B i3_basis2.py 26 $PWD/cache/loc/J26_x03
python3 -B i3_mkset.py 18 $PWD/cache/loc/J18_x03  # writes set 'mine' into the sandbox copy cache/sbx/ads4/code/b5fold
cd cache/sbx/ads4/code/b5fold && bash lane.sh i3_rerun_A.log B5_SET=mine python3 b5_state.py 18 1.02,1.03,1.05 ';;' python3 b5_window.py 18
bash lane.sh i3_rerun_B.log python3 b5_odd.py 25 main 1.03,1.04,1.05 ';;' python3 b5_odd.py 23 main 1.03,1.05,1.075,1.1
bash lane.sh i3_rerun_C.log python3 b5_sum.py
```
(`i3_basis.py` is the unconstrained-fit diagnostic of §2.1.)
