#!/usr/bin/env bash
# Reproduce every certified statement of ads4/FLAT_GRAVITON_BOUND.md.
# Prints one PASS/FAIL line per certificate and a final verdict; full logs in logs/.
# Requirements: python3 with python-flint (>=0.6, Arb/FLINT), sympy; numpy/scipy/mpmath only for --sanity.
# Runtime: about 4.5 minutes (cert_P ~60 s, cert_S ~40 s, cert_cg ~20 s).
# NB: cert_jhep PASS means its stated findings were verified, including that the JHEP-printed functional is NEGATIVE.
#   bash run_all.sh            certificates only
#   bash run_all.sh --sanity   also the floating-point checks (not part of the proof)
set -u
cd "$(dirname "$0")"
mkdir -p logs
JP=400; BA=120; BS=1; U1=1/100
overall=PASS
run () {   # name, command...
  local name=$1; shift
  "$@" > "logs/$name.log" 2>&1
  local v
  v=$(grep -E '^VERDICT' "logs/$name.log" | tail -1)
  if [[ "$v" == *PASS ]]; then echo "PASS  $name   ($v)"; else echo "FAIL  $name   (${v:-no verdict})"; overall=FAIL; fi
}
python3 -c "import flint, sympy" || { echo "FAIL: python-flint/sympy missing (pip install python-flint sympy)"; exit 1; }
run cert_P       python3 cert_P.py $JP                 # J < 400, all m >= 1: exact polynomials, Descartes
run cert_C       python3 cert_C.py $JP                 # J >= 400, 1 <= m <= 2: crude bounds
run cert_M       python3 cert_M.py $BA $U1 4000        # m >= 2, bhat >= 120: Mehler expansion + remainder
run cert_S_pp    python3 cert_S.py ++ $BS $BA $JP      # J >= 400, m >= 2, bhat <= 120  (++)
run cert_S_pm    python3 cert_S.py +- $BS $BA $JP      # J >= 400, m >= 2, bhat <= 120  (+-)
run cert_cg      python3 cert_cg.py                    # Clebsch-Gordan identity (exact, symbolic in J)
run cert_low     python3 cert_low.py                   # low-energy side, constants A, B, matter signs, F_g4 g5 term
run cert_ir      python3 cert_ir.py                    # I(lambda), phi-bounds; eps<=1/100 Delta_IR table = ESTIMATE only (see cert_eps)
run cert_eps     python3 cert_eps.py                   # IR: negative set of F^eps localised (eps<=1/20000) + Delta_IR under E
run cert_jhep    python3 cert_jhep.py                  # JHEP-printed (3.3a): low side, heavy action, exact COUNTEREXAMPLE
if [[ "${1:-}" == "--sanity" ]]; then
  python3 sanity_numeric.py > logs/sanity_numeric.log 2>&1 && echo "(sanity checks written to logs/sanity_numeric.log; not part of the proof)"
fi
echo "FINAL VERDICT: $overall"
