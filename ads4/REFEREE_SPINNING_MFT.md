# Referee report: `ads4/SPINNING_MFT.md` and `ads4/code/mft/`

Referee: independent, hostile. Date: 16 Nov 2026. Object: commit `12bd065`. Nothing under review was edited.
Scratch outputs: `ads4/code/referee_e2e/mft_exact_fresh.txt`.

## Verdicts

| | Question | Verdict |
|---|---|---|
| Correctness | Fock method, spectrum, closed forms, heavy-limit algebra | **No error found.** |
| Soundness | Exact data and checks | **Sound for the exact data.** Reproduced from an empty cache: 145/145 "closed = 1", in about 2 min. |
| Labels | Honesty of the status marks | **Two overstatements.** The "end-to-end DV" of the MFT normalisation (issue 1). Attributing the M1 resolution to MFT normalisation (issue 3). The odd sign and the `⟨TTTT⟩` statement are honestly labelled [I], but the second has little content (issues 2, 4). |
| Applicability | Use for the AdS4 bound | **One spinning pair only.** Nothing on `⟨TTTT⟩`, two structures per parity, or the identical-particle `[TT]` data. |
| Novelty | Against the literature | **Probably partly known.** See the last section. |

No fatal issues.

## Claim (1): Fock-space MFT coefficients

### Primary count is correct

I checked it against characters. The character of a conserved spin-`l` current is `q^{l+1}(χ_l − qχ_{l−1})/(1−q)³`. So the primaries of `[Oφ]` at level `N` are

`spin l ⊗ Sym^N(vector) ⊖ spin (l−1) ⊗ Sym^{N−1}(vector)`.

For `l = 1, 2` and `N ≤ 3` this gives each `J ∈ [l, N+l]` exactly once, matching §3.1.

Null states:
* For generic `Δφ`, `Δ − J − 1 ≥ Δφ > 0`, so every `X` is long and there are no nulls on the `X` side.
* The null descendants of the conserved `O` are removed by working on the range of the singular kernel matrix. The code asserts that every wavefunction lies in that range, and the pivot-inverse norm equals the pseudo-inverse norm there.

Minor: the code only scans `l ≤ J ≤ N+l`. The statement "none for `J < l` or `J > N+l`" is not tested by the code; cite the character argument instead.

### Norm → OPE coefficient map is consistent with SPINNING_E2E.md

* The structure basis is complete: `{z·x, ζ·x, z·ζ, x²}`, and `ε` appears once.
* `V_1 → −z_1·x` holds exactly in the limit, so the `S_n(∞)` normalisation is fixed.
* The relative factors `(−1)^{l−n}(−½)^n` are pinned by check (iii): the null vector of the divergence matches the `L_n` of `structures.py`.
* `(−½)^J` is pinned by the scalar FK check.
* Hand check: `(N,J) = (0,1)`, `l = 1`. Here `L_1 = 11` and `S_1(∞) = −½ ζ·z`. A unit-normalised `Oφ` then gives `p[𝔏] = 4/121`. The closed form gives `2·(7/11)/38.5 = 4/121`.
* Second hand check: odd `(1,1)` gives `7/11` from both the data and the closed form.

### The 145 matches are not circular

The closed forms were recognised on the data and then confirmed. `Δφ` enters only through `P_FK`, and this is checked at three values of `Δφ`, which is non-trivial. The forms are [NE], as labelled.

### Issue 1 (major, labels): there is no end-to-end test of the MFT normalisation

* Check (3) ("`1.000000000000000` in 72 cases") is an algebraic identity: closed forms × the E2E formula, compared with the simplified form.
* §5 is `R_e2e × Q_E` by construction (the text says "exactly").
* So the new ingredient, the bridge between the Fock `p^MFT` and the block normalisation, is checked only through `S_n(∞)` and the E2E OPE limit. It is never checked through a block decomposition of the MFT correlator.
* Status rows 5/8 and §0.3 ("end-to-end DV") should drop "end-to-end".
* Fix: test the leading-twist `τ = τ_0` family, which is all even, in the lightcone limit of `⟨O_1O_3⟩⟨φ_2φ_4⟩`, using the E2E even blocks.

## Claim (2): Heavy limit

* **Even parity.** The algebra is correct. I verified `a_{l+1,Δφ}/a_{1,Δφ} = l!(½)_l` and `w_l·l!(½)_l = ¼`. The cancellation of `L_0²`, `2^J`, `b` and the "pol" polynomials is exact. [D from NE closed forms] is honest.

### Issue 2 (major): fixing the odd `l = 1` sign by positivity is not legitimate as a derivation in a positivity programme

* The sign decides whether odd operators enter the preserved component as `+|ĉ_O|²` or `−|ĉ_O|²`. That is exactly what a bound would rely on.
* The choice made agrees with flat-space unitarity (preserved = `Σ_h|c^h|²`). The [I] label is honest.
* It can and should be derived, in one of two ways:
  * (a) by the explicit Wick-rotation phase of `ε`;
  * (b) by the identity behind claim 10. The polarisation trace of `⟨O_1O_3⟩` is the character of `I(x_13)` on spin `l`, which is a constant, times `x_13^{−2Δ_O}`. So the spin-traced MFT correlator *is* the scalar MFT `⟨χχ⟩⟨φφ⟩` with `Δ_χ = Δ_O`. Its known density, together with the derived even half, then fixes the odd sign.
* As written, claim 10 for `l = 1` is marked [D], but it depends on the [I] sign. Mark it [I], or derive it via (b).
* The odd shadow coefficient `b^{(2−Δφ)}` and the odd OPE normalisation are "by analogy, not tested". The odd formula therefore has no numerical test at all.
* Odd `l = 2`: [ND], honest.

## Claim (3): The "resolution" of M1 (reproduced)

### The zero-parameter test is genuine

`c_2(2) − c_2(1) + (κ_2−κ_1)·pol = −6P′/P`. Here `pol` is derived and ν-independent, and the `B_3` term cancels between `l = 1` and `l = 2`.

I reproduced it with the full-precision E2E fit outputs:

| `β` | residuals over `ν = 0, 1, 3` |
|---|---|
| 0.4 | `0.039, −0.001, −0.002` |
| 1 | `0.017, −0.008, −0.062` |

By hand at `β = 0.4`, `ν = 3`: `−50.01 + 25.47 + 2.25·25.630 = 33.128`, and `−6P′/P = 33.13`.

Per-`β` `d_β` with skeleton `(2,Δφ)` as reference: `−8.001` and `−8.028`. With skeleton `(3,Δφ)`: `−8.34` and `−10.79`, as in my M1.

### Issue 3 (major, attribution and labels)

MFT normalisation changes `c_2` only by ν-independent terms (`κ_l·pol`, `B_3 coshβ`). It therefore **cannot change any `d_β`**. The "resolution" comes entirely from changing the reference:
* from skeleton `(3,Δφ)` to skeleton `(2,Δφ)`, or
* to the `l = 1` data.

That is a single discrete, post-hoc choice, admitted to be "not understood". It is not a continuous re-fit:
* the difference test has 0 free parameters;
* the per-`β` fits have 2 parameters for 3 points.

But "after MFT normalisation … resolves M1" (§0.5, §6) misattributes the cause.

Fix: state that
* the `l`-dependence `−2(l_2² − l_1²)P′/P` is confirmed with no free parameter;
* the absolute per-`l` form `−2l²P′/P` holds relative to an empirically chosen skeleton `(2,Δφ)` [NE, reference choice unexplained].

## Claim (4): `⟨TTTT⟩` "2P/m^{2k}"

### Issue 4 (minor/major, content)

"Each MFT operator gives exactly `2P/m^{2k}`" is true **by definition** of `C^{pres} := (8/M_pres)C`. "The 2 is the helicity count" is interpretation. The `⟨TTTT⟩` extension ("per-pair ¼ absorbed in the kernel normalisation") is [I], honestly labelled, but it is unfalsifiable as stated.

The real `⟨TTTT⟩` questions are untouched:
* the relative weight and sign of the `±` terms;
* two structures per parity;
* the identical-particle factor;
* positivity.

CHL already gives the MFT `[TT]` data in the helicity basis; use it.

## Minor issues

* The run time is about 2 min, not 15.
* Report the exact-rational spectrum check for `N ≤ 8` (`l = 1`), not only `N ≤ 5`.
* §4.1 check (2) is a Gamma-function identity check, not a block check. Say so.

## Applicability

The work supplies the MFT ("unitarity") normalisation for one spinning Regge pair, even parity, `l = 1, 2`, and odd `l = 1` up to the sign.

It does not supply:
* the doubly spinning `⟨TTTT⟩` density, its 2×2 structure mixing, or its `±` weights;
* odd `l = 2`;
* any test beyond the E2E range (`k = 4, 6`, `ν ≪ m`).

## Novelty

The Fock-space method is a nice, block-free exact route.

The `[Oφ]` MFT coefficients with spinning `O` are, however, given in general `d` by Sleight–Taronna arXiv:1804.09334, which 1807.05941 cites for "`[O_J O]` with `O_J` of arbitrary spin". They are also obtainable from Karateev–Kravchuk–Simmons-Duffin arXiv:1809.05111 (pairing/shadow formalism).

The leading-twist `[Jφ]` and `[Tφ]` coefficients appear in lightcone-bootstrap work (Li–Meltzer–Poland; Hofman–Li–Meltzer–Poland–Rejon-Barrera, arXiv:1603.03771). The structural form "scalar MFT × CHL helicity factor" is CHL's observation transplanted.

The paper must check these papers before claiming the closed forms (especially the conserved, d = 3, parity-odd pieces) as new. Likely new: the d = 3 parity-odd closed forms for conserved `O`, and the exact cancellation with the spinning shadow coefficients in the heavy limit.
