#!/usr/bin/env bash
# Reproduce FLAT_UNSTRIPPED.md.  Outputs in cache/ (git-ignored).  ~20 s; scan_FU (floating point, [NE]) ~2 min.
#   bash run_all.sh            exact/Arb parts
#   bash run_all.sh --scan     also the floating-point scan of F_U
set -u
cd "$(dirname "$0")"; mkdir -p cache
overall=PASS
for s in sumrules cert_FU_negative; do
  python3 $s.py > cache/$s.log 2>&1
  v=$(grep -E '^VERDICT' cache/$s.log | tail -1)
  if [[ "$v" == *PASS ]]; then echo "PASS  $s   ($v)"; else echo "FAIL  $s   (${v:-no verdict})"; overall=FAIL; fi
done
if [[ "${1:-}" == "--scan" ]]; then
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python3 scan_FU.py > cache/scan_FU.log 2>&1 && echo "(scan_FU [NE] -> cache/scan_FU.log)"
fi
echo "FINAL VERDICT: $overall"
