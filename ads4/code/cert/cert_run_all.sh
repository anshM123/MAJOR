#!/usr/bin/env bash
# Reproduces ads4/CERT_CLOSE.md (single core, <= 2 processes).  Logs in cache/ (git-ignored).  ~2.5 h in total.
#   bash cert_run_all.sh nf     NF02 flat certificate at margin 1/5 (all regions)            ~45 min (+30 min optional P extension)
#   bash cert_run_all.sh q      Q3* flat certificate pieces at margin 9/50 + float scans      ~45 min
#   bash cert_run_all.sh k2     continuum k = 2 weights, X^Q re-solve, export to prl/data/k2   ~15 min
set -u -o pipefail
cd "$(dirname "$0")"; mkdir -p cache
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1; ulimit -v 3000000
R=../robust
what=${1:-all}
if [[ $what == nf || $what == all ]]; then
  python3 -u nf_exact_min.py                                   > cache/nf_exact_min.log
  python3 -u nf_lemmaB_check.py                                > cache/nf_lemmaB_check.log
  (cd $R && python3 -u nf_cert.py NF02 1/5 P 400)              > cache/NF02_P_1-5.log
  (cd $R && python3 -u nf_cert.py NF02 1/5 M 120)              > cache/NF02_M_1-5.log
  (cd $R && python3 -u nf_cert.py NF02 1/5 Cpp 400)            > cache/NF02_Cpp_1-5.log
  (cd $R && python3 -u nf_cert.py NF02 1/5 Cpm 400 1/4,5/16,3/8) > cache/NF02_Cpm_1-5.log
  (cd $R && python3 -u nf_cert.py NF02 1/5 S++ 120 400)        > cache/NF02_S++_1-5.log
  (cd $R && python3 -u nf_cert.py NF02 1/5 S+- 120 400)        > cache/NF02_S+-_1-5.log
  python3 -u nf_corner.py 1/5 400                              > cache/NF02_corner_1-5_J400.log
  [[ ${NF_PEXT:-0} == 1 ]] && python3 -u nf_Pext.py 1/5 400 800 > cache/NF02_Pext_1-5_400_800.log
fi
if [[ $what == q || $what == all ]]; then
  python3 -u q_cert.py dense 3000                              > cache/Q3_dense.log
  python3 -u q_scan_bhat.py 20000                              > cache/Q3_scan_bhat.log
  python3 -u q_backward.py 20000                               > cache/Q3_backward.log
  for reg in "P 400" "C 400" "S++ 120 400" "S+- 120 400" "M 120 1/100 4000"; do
    set -- $reg; python3 -u q_cert.py 9/50 $reg > cache/Q3s_${1}_9-50.log
  done
fi
if [[ $what == k2 || $what == all ]]; then
  for w in "D 300" "D 1000" "Q 300" "Q 1000"; do set -- $w; python3 -u k2_cont.py $w > cache/k2_${1}_${2}.log; done
  python3 -u k2_xq_check.py 300 1000                           > cache/k2_xq_check.log
  python3 -u k2_export.py && (cd ../../prl/data/k2 && sha256sum -c SHA256SUMS)
fi
grep -h '^VERDICT' cache/*.log
