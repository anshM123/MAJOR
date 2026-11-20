#!/bin/bash
# ee1000 ensemble (EE_DG1000_LARGEJ.md §4/§6): the large-J chain in all representation variants, one tag each; then ee1000_tables.py.
# main: shape fit K = 4 (best typical cross-validation, §4) on the nodes 0.05 <= x <= 0.8 (K2 not reliable at nu <~ 50 [300-agent direct; EE_NEARGAP §5] nor at x > 0.8 for
# J >= 40); variants K = 5, 6; x <= 0.9; spline; with the x = 0.025, 0.0375 nodes; 64 p-nodes.
# usage: EE1000_JEVEN=... EE1000_JODD=... ./ee1000_ensemble.sh <prefix>
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1
P=$1
EE1000_MODEL=shape EE1000_K=4 EE1000_XMAX=0.8 EE1000_XMIN=0.05 EE1000_TAG=${P}_main python3 ee1000_chain.py > /dev/null
for K in 5 6; do EE1000_MODEL=shape EE1000_K=$K EE1000_XMAX=0.8 EE1000_XMIN=0.05 EE1000_TAG=${P}_K$K python3 ee1000_chain.py > /dev/null; done
EE1000_MODEL=shape EE1000_K=4 EE1000_XMAX=0.9 EE1000_XMIN=0.05 EE1000_TAG=${P}_x09 python3 ee1000_chain.py > /dev/null
EE1000_MODEL=spline EE1000_XMAX=0.8 EE1000_XMIN=0.05 EE1000_TAG=${P}_spline python3 ee1000_chain.py > /dev/null
EE1000_MODEL=shape EE1000_K=4 EE1000_XMAX=0.8 EE1000_XMIN=0 EE1000_TAG=${P}_withlo python3 ee1000_chain.py > /dev/null
EE1000_MODEL=shape EE1000_K=4 EE1000_XMAX=0.8 EE1000_XMIN=0.05 EE1000_NPN=64 EE1000_TAG=${P}_npn64 python3 ee1000_chain.py > /dev/null
python3 ee1000_tables.py $P main K5 K6 x09 spline withlo npn64
