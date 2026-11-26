#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p rerun
python3 check_identity.py | tee rerun/identity.log
python3 check_lowfreq.py | tee rerun/lowfreq.log
for f in pilot wide full_split wider_d m24 nu20 lower lower40 combo complex_top majorant_sweep combo_preflight; do
  node check_transfer.js "$f.json" "rerun/${f}_results.json" > "rerun/$f.log"
  echo "Ran $f; see output statuses (intentional failures included)"
done
