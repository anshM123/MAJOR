#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
JOBS="${JOBS:-4}"
g++ -O3 -std=c++17 scalar.cpp -lgmpxx -lgmp -o scalar
python3 check_identity.py
python3 check_base_identity.py
python3 check_lowfreq.py
python3 check_norm.py cross_sections.json --jobs "$JOBS" --output rerun_cross.json
python3 check_norm.py tensor_run1.json --jobs "$JOBS" --output rerun_tensor1.json
python3 check_norm.py complex_pilot.json --jobs "$JOBS" --output rerun_complex.json
python3 assemble_top.py rerun_complex.json --dlo .040355 --dhi .041169 --lo 895 --hi 905
python3 check_norm.py closure5.json --jobs "$JOBS" --output rerun_closure.json
python3 assemble.py rerun_closure.json
if [[ "${LONG:-0}" == 1 ]]; then
    python3 check_norm.py tensor_run2_wide.json --jobs "$JOBS" --max-nodes 300000 --output rerun_tensor2.json
fi
