#!/usr/bin/env bash
# Reproduces every certified statement of ads4/UNIFORM_POSITIVITY.md.
#   usage:  bash run_all.sh            (certificates + non-rigorous sanity checks, about 1 minute)
#           bash run_all.sh --no-sanity (certificates only, a few seconds)
# Each certificate prints one PASS/FAIL line per box/cell and a final VERDICT line; the full per-box output is
# written to logs/<name>.log and the SUMMARY/VERDICT lines are echoed here.
set -u
cd "$(dirname "$0")"
export PYTHONDONTWRITEBYTECODE=1
mkdir -p logs
overall=0
run_cert () {
  name="$1"; shift
  echo "---- $name ----"
  if python3 "$@" > "logs/$name.log" 2>&1; then :; else overall=1; fi
  grep -E "^(SUMMARY|VERDICT|PASS \(|FAIL)" "logs/$name.log" | head -n 40
  nfail=$(grep -c "^FAIL" "logs/$name.log" || true)
  npass=$(grep -c "^PASS" "logs/$name.log" || true)
  echo "     [$npass PASS lines, $nfail FAIL lines; full log: logs/$name.log]"
  if [ "$nfail" != "0" ]; then overall=1; fi
  if ! grep -q "^VERDICT.*PASS$" "logs/$name.log"; then overall=1; fi
}
run_cert iv_selftest     iv_selftest.py
run_cert C1_S            cert_S.py
run_cert C2_smallM       cert_smallM.py
run_cert C3_kernel       cert_K.py
run_cert C4_expansion    cert_expansion.py
run_cert C5_flat         cert_flat.py
echo "=================================================================="
if [ "$overall" = "0" ]; then
  echo "FINAL VERDICT: ALL CERTIFICATES PASS  (Theorem 1: F_M(beta) > 0 for all M > 0, beta > 0; Theorem 2: fhat > 0 on [0,oo))"
else
  echo "FINAL VERDICT: FAIL (see logs)"
fi
if [ "${1:-}" != "--no-sanity" ]; then
  echo "---- non-rigorous numerical sanity checks (NOT part of the proof) ----"
  python3 sanity_numeric.py 2>&1 | tee logs/sanity_numeric.log
fi
exit $overall
