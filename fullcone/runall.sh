#!/bin/bash
# Batch driver for schottky_cone.py (3 jobs at a time to bound memory).
cd "$(dirname "$0")"
for t in pants torus; do
  for N in 3 4 5; do
    case $N in 3) T=60;; 4) T=40;; 5) T=24;; esac
    python3 -I schottky_cone.py $t $N $T $((N*10+${#t})) > run_${t}_N$N.log 2>&1 &
  done
  wait
done
