#!/usr/bin/env bash
# Reproduces ads4/NOFORWARD_FUNCTIONAL.md (single core).  Outputs in cache/nf/ (git-ignored).
#   bash nf_run_all.sh        everything (~2-3 h on one core)
#   bash nf_run_all.sh cert   certificates of the tracked NF02 only (functionals_noforward.json)
set -u -o pipefail
cd "$(dirname "$0")"
mkdir -p cache/nf/logs
if [[ "${1:-all}" == all ]]; then
  python3 -u nf_contacts.py 7 16 32 | tee cache/nf/contacts.log
  for rd in "3 6" "5 6" "5 8" "6 8"; do python3 -u nf_lp.py $rd 0.0 0.05 0.1 0.15 0.2 0.3; done | tee cache/nf/lp.log
  python3 -u nf_lp.py 5 12 0.0 0.1 0.2 | tee -a cache/nf/lp.log
  python3 -u nf_lp.py 11 10 0.0 0.1 0.2 | tee -a cache/nf/lp.log
  python3 -u nf_lp.py 5 8 proxy=120 0.0 0.1 0.15 0.2 | tee -a cache/nf/lp.log
  for rd in "11 10" "12 10" "12 12"; do python3 -u nf_lp.py $rd proxy=120 0.0 0.15 0.2; done | tee -a cache/nf/lp.log   # section 4.1 cost
  python3 nf_mkfun.py NF02 5 8 0.2 _proxy120     # regenerates cache/nf/functionals_nf.json only; the tracked file is not touched
  python3 -u nf_check.py NF02 3000 | tee cache/nf/check.log
fi
python3 -u nf_lowcheck.py NF02               > cache/nf/logs/NF02_lowcheck.log 2>&1
python3 -u nf_cert.py NF02 0 low             > cache/nf/logs/NF02_low.log 2>&1
python3 -u nf_cert.py NF02 3/20 M 120        > cache/nf/logs/NF02_M.log 2>&1
python3 -u nf_cert.py NF02 3/20 Cpp 400      > cache/nf/logs/NF02_Cpp.log 2>&1
python3 -u nf_cert.py NF02 3/20 Cpm 400 1/4,5/16,3/8 > cache/nf/logs/NF02_Cpm.log 2>&1
python3 -u nf_cert.py NF02 3/20 P 400        > cache/nf/logs/NF02_P.log 2>&1
python3 -u nf_cert.py NF02 3/20 S++ 120 400  > cache/nf/logs/NF02_S++.log 2>&1
python3 -u nf_cert.py NF02 3/20 S+- 120 400  > cache/nf/logs/NF02_S+-.log 2>&1
grep -h '^FUNCTIONAL\|^VERDICT' cache/nf/logs/*.log
