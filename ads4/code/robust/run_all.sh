#!/usr/bin/env bash
# Reproduces ads4/ROBUST_FUNCTIONAL.md.  Outputs in cache/ (git-ignored).  Machine with 4 cores: about 3-4 hours in total
# (grids ~10 min, LPs ~1 h, certificates ~1.5 h, dense checks ~30 min, AdS scans ~1 h).
#   bash run_all.sh            everything
#   bash run_all.sh cert       certificates only (needs cache/functionals.json, which the LP stage writes)
# Certificates (FLINT/Arb, exact rationals) print 'VERDICT ...: PASS/FAIL'; everything else is floating point ([NE]).
set -u -o pipefail
cd "$(dirname "$0")"
mkdir -p cache/logs
python3 -c "import flint, sympy, scipy, numpy" || { echo "need python-flint, sympy, scipy, numpy"; exit 1; }
stage=${1:-all}
if [[ "$stage" == all ]]; then
  # 1. LPs (floating point): uniform-margin Pareto (no M-proxy), certifiable Pareto (with Lemma-M proxy), envelope-targeted family
  python3 -u pareto.py 4,4 6,6          | tee cache/pareto_a.log
  python3 -u pareto2.py 4 4             | tee cache/pareto2_44.log
  python3 -u pareto2.py 6 6             | tee cache/pareto2_66.log
  python3 -u mtest.py                   | tee cache/mtest.log
  python3 -u target.py 4 4 0.02         | tee cache/target_44.log
  # 2. rationalise the chosen functionals (exact coefficients in cache/functionals.json)
  python3 mkfun.py M01 mtest.json delta=0.1 BA=120.0 Y=2.0
  python3 mkfun.py M02 mtest.json delta=0.2 BA=120.0 Y=2.0
  python3 mkfun.py T1n1.5 target_d4_4_0.02.json C1=1 n=1.5
  python3 mkfun.py T1n1.75 target_d4_4_0.02.json C1=1 n=1.75
  python3 mkfun.py T3n2 target_d4_4_0.02.json C1=3 n=2
fi
# 3. certificates: robust margin delta' = 9/50 (M02), 9/100 (M01); regions P, C, M, S++, S+-; exact low side
for spec in "M02 9/50" "M01 9/100"; do
  set -- $spec
  python3 -u lowcheck.py $1            > cache/logs/$1_lowcheck.log 2>&1
  python3 -u rcert.py $1 0 low         > cache/logs/$1_low.log 2>&1
  python3 -u rcert.py $1 $2 P 400      > cache/logs/$1_P.log 2>&1
  python3 -u rcert.py $1 $2 C 400      > cache/logs/$1_C.log 2>&1
  python3 -u rcert.py $1 $2 M 120 1/100 4000 > cache/logs/$1_M.log 2>&1
  python3 -u rcert.py $1 $2 S++ 120 400 > cache/logs/$1_S++.log 2>&1
  python3 -u rcert.py $1 $2 S+- 120 400 > cache/logs/$1_S+-.log 2>&1
done
grep -h '^VERDICT' cache/logs/*.log
if [[ "$stage" == all ]]; then
  # 4. dense floating-point checks and AdS scans
  python3 -u check.py M02 M01 Fg3 T1n1.5 T1n1.75 T3n2 | tee cache/check.log
  python3 -u ads.py M02 regge 20 100 1000 | tee cache/ads_regge_M02.log
  python3 -u ads.py M01 regge 20 100 1000 | tee cache/ads_regge_M01.log
  python3 -u ads.py Fg3 regge 20 100 1000 | tee cache/ads_regge_Fg3.log
  python3 -u ads.py M02 hybrid 100 1000   | tee cache/ads_hybrid_M02.log
  python3 -u adsT.py                       | tee cache/adsT.log
  python3 -u bound.py                      | tee cache/bound.log
fi
