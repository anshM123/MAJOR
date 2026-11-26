#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
# Recompute in a disposable copy if historical outputs need to be preserved.
JOBS=${JOBS:-3}
g++ -O3 -std=c++17 scalar.cpp -lgmpxx -lgmp -o scalar
# All input boxes have directed continuum majorants; newer checker may choose
# tighter grids. Use the retained historical configs rather than a new float grid.
python check_norm.py resid_3p5_3p65.json --jobs "$JOBS" --max-nodes 150000 --output resid_3p5_3p65_results.json > resid_3p5_3p65.log
python assemble35.py > assembly35.log
for slab in 3.35_3.5 3.2_3.35; do
  python check_norm.py "resid_$slab.json" --jobs "$JOBS" --max-nodes 150000 --output "resid_${slab}_results.json" > "resid_$slab.log"
  python check_norm.py "p4_$slab.json" --jobs "$JOBS" --max-nodes 150000 --output "p4_${slab}_results.json" > "p4_$slab.log"
  python - "$slab" <<'PY'
import json,sys
s=sys.argv[1]
json.dump(json.load(open(f'resid_{s}_results.json'))+json.load(open(f'p4_{s}_results.json')),open(f'combined_{s}.json','w'),indent=2)
PY
  python resid_budget.py "resid_${slab}_results.json" 0.00000187
  for p in 3 4; do for b in 0 4; do
    python ellipse_extend.py "combined_$slab.json" 10.25 "ct_${slab}_${p}_${b}.json" "$p" "$b" > "ct_${slab}_${p}_${b}.log"
  done; done
done
python - <<'PY'
import json
r=json.load(open('resid_2.9_3.2.json'))[:2]
json.dump(r,open('resid_low_2.9_3.2.json','w'),indent=2)
PY
python check_norm.py resid_low_2.9_3.2.json --jobs "$JOBS" --max-nodes 150000 --output resid_low_2.9_3.2_results.json > resid_low_2.9_3.2.log
python check_norm.py resid_coupled_hi_2.9_3.2.json --jobs "$JOBS" --max-nodes 150000 --output resid_coupled_hi_2.9_3.2_results.json > resid_coupled_hi_2.9_3.2.log
python check_norm.py p4_coupled_2.9_3.2.json --jobs "$JOBS" --max-nodes 150000 --output p4_coupled_2.9_3.2_results.json > p4_coupled_2.9_3.2.log
python - <<'PY'
import json
r=json.load(open('resid_low_2.9_3.2_results.json'))+json.load(open('resid_coupled_hi_2.9_3.2_results.json'))
json.dump(r,open('resid_combo_2.9_3.2_results.json','w'),indent=2)
json.dump(r+json.load(open('p4_coupled_2.9_3.2_results.json')),open('combined_2.9_3.2.json','w'),indent=2)
PY
python resid_budget.py resid_combo_2.9_3.2_results.json 0.0000025
for p in 3 4; do for b in 0 4; do
  python ellipse_extend.py combined_2.9_3.2.json 10.25 "ct_2.9_3.2_${p}_${b}.json" "$p" "$b" > "ct_2.9_3.2_${p}_${b}.log"
done; done
python assemble_full.py | tee full_ledger.log
