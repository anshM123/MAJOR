# Certified near-gap positivity of the route-A′ D-safe action (scalar-skeleton model), Δ_g = 100

Status: 18 Nov 2026, H13′/certification agent. **Unrefereed. Not committed.**

**What changed.** This replaces `CERTIFIED_NEARGAP.md` (Theorems C, C′), which certified the withdrawn B′-uplift `Ω^{NF,SC}`. The
functional is now route A′, `Ω^{A′} = ∫ψ̃₂Λ⁴Ĉ₂ + ∫ψ₃½Λ⁶Ĉ₃` (`H0_SPINNING_FUNCTIONAL.md` §R, §10, §13).

**Code.**
* New: `ads4/code/exactskel/cert_actionA.py` and `certA_inv.py`.
* Reused read-only: `cert_moments.py`, `cert_quad.py`, `h13_common.py`, `h13a_common.py`, and the `spinfunc/h1`/`h2` pickle and helpers.
* Caches: `exactskel/cache/certA_100_Df{3,2}.pkl`, `.log`, `certA_inv_100.txt` (git-ignored).

**Grades.**
* **[certified]**: Arb (192 bits); every truncation replaced by a proven bound; inequalities checked in Arb or exact rationals.
* **[D]**: derived.
* **[DV]**: derived and numerically validated (floating point).
* **[NE]**: numerical evidence.
* **[M]**: model.
* **[I]**: assumption.

__BODY__
