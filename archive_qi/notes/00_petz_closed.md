# Petz recovery conjecture — closed line (not a new result)

Target considered: I(A;C|B) >= -log F(rho_ABC, Petz_{B->BC}(rho_AB)) (Berta–Seshadreesan–Wilde 2014).

Status found during this session (2026-11-15):
* General-channel form disproved by Bhattacharya, arXiv:2309.02074 (qubit pinching example).
* CMI form listed as SOLVED (negatively) in the QIQCOP Zoo, entry op_87c77263c8bab523,
  three certified three-qubit counterexamples dated Sept 2026 (ratios I/(-log F) ~ 0.991, 0.979, 0.932).

What this session found independently (rediscovery, kept for the record; code in code/petz_*.py, code/mp_cmi.py):
* Two-qubit violation of the general-channel form, verified at 60 digits (gap -0.0127).
* Second-order reduction with A a classical bit: rho = sum_{a=+-} (1/2)|a><a| (x) (sigma + a*eps*delta) gives
  I = (eps^2/2) Q1, -log F = (eps^2/2) Q2 with
  Q1 = g^BKM_sigma(delta) - g^BKM_{sigma_B}(delta_B),  Q2 = (1/2) g^SLD_sigma(delta - Petz(delta_B)).
  For product sigma, Q1 >= 2 Q2 (proved in-session). For nearly rank-2 sigma on two qubits the ratio Q1/Q2
  drops to ~0.688 (numerical infimum of a continuation scan; local optimum only) — lower than the published examples.
  Verified on actual states with an independent 80-digit implementation (e.g. eps=0.01: I=6.3004e-7 < -logF=6.6276e-7).
* Not pursued further: the conjecture is already settled; only the optimal constant remains, which is not a major question.

# Other closed/side lines from the first phase (quantum-information triage; out of scope after re-brief)

* Zhang–Yeung for von Neumann entropies (QIQCOP Zoo, listed unsolved): no violation found by (i) random-start ratio
  maximisation in dims up to (3,3,2,2), ranks 1..16 (all runs saturate at ratio exactly 1), (ii) first-order
  (eps log 1/eps) support-operator test and direct minimisation near saturating one-qubit-per-party stabilizer states,
  (iii) second-order BKM test around rho_A (x) rho_{B1C} (x) rho_{B2D} (max ratio ~0.5). Code: entropy_opt.py,
  zy_search.py, zy_structured.py, zy_second_order.py.
* Triage of three conjectured inequalities (all hold numerically, equality cases found): trace-exponential lower bound
  for matrix-word averages; two-copy p-norm multiplicativity of polarized Werner–Holevo channels (1<p<2, d=3,4);
  quantum Mrs. Gerber lower bound (equality on pure-state side information).
* Transpose-degradable vs degradable channels (open since Brádler–Dutil–Hayden–Muhammad 2010): reformulated as
  TE (states omega_{AEE'} with omega_{AE'} = omega_{AE}^{T_E}) versus symmetric 2-extendibility; SDP support-function
  comparison in dims 2x2, 2x4, 4x2, 3x3 found TE-optimal states always within 1e-7 of the 2-extendible set.
  Suggests the theorem "transpose-degradable => degradable"; NOT proven, not pursued. Code: tdeg_search.py, te_vs_2ext.py.
