#!/bin/bash
# H_OO (HYPOTHESES_CLOSE2.md): continuation of oo_queueD.sh (resumable; waits for the running oo_points job).
# J = 18 (336, 400), 20 (400), 22 (336): 8 nodes nu ~ 8.3-40; J = 16, 14, 12 (400): 6 nodes nu ~ 8.4-20 (K2 kept at 29.9, 39.8);
# J = 4 (400): validation spot nodes nu ~ 10, 20, 40.
cd "$(dirname "$0")"
while pgrep -f "oo_points.py 18 336" > /dev/null; do sleep 30; done
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OO_M=10 OO_R=8
X336=0.025,0.03,0.035,0.04,0.05,0.0625,0.075,0.1
X400=0.021,0.025,0.0293,0.0335,0.042,0.05,0.075,0.1
X400S=0.021,0.025,0.0293,0.0335,0.042,0.05
python3 oo_points.py 18 336 $X336
python3 oo_points.py 18 400 $X400
python3 oo_points.py 20 400 $X400
python3 oo_points.py 22 336 $X336
python3 oo_points.py 4 400 0.025,0.05,0.1
python3 oo_points.py 16 400 $X400S
python3 oo_points.py 14 400 $X400S
python3 oo_points.py 12 400 $X400S
echo "QUEUE OO-D2 DONE"
