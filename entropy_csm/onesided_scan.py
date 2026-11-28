#!/usr/bin/env python3
"""One-sided W = V_+ scan (BBP 1703.08460 sec. 7). Incremental: appends one line per (N_c, N_f) to
scan_monopole_onesided_output.txt and skips points already present.  Usage: python3 onesided_scan.py Nc"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scan_monopole import onesided, Fch

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "scan_monopole_onesided_output.txt")
done = set()
if os.path.exists(OUT):
    for l in open(OUT):
        p = l.split()
        if len(p) > 2 and p[0].isdigit():
            done.add((int(p[0]), int(p[1])))
else:
    with open(OUT, "w") as f:
        f.write("Nc Nf | R_Q eta R[V_+]=2 R[V_-] R[meson] | C_T thr7/8 ratio | F F(Higgs-free) | flags\n")

for Nc in [int(a) for a in sys.argv[1:]]:
    for Nf in range(Nc + 2, 3 * Nc + 5):
        if (Nc, Nf) in done:
            continue
        res = onesided(Nc, Nf)
        if res is None:
            line = f"{Nc} {Nf} | no convergent window"
        else:
            D, eta, RVm, F, CT, edge = res
            dimH = 2 * Nc * Nf - Nc * Nc
            ratio = CT / (4 * dimH)
            flags = []
            if ratio < 7 / 8: flags.append("VIOL")
            if 2 * D < 0.5 or RVm < 0.5: flags.append("UNITARITY")
            if edge: flags.append("F-max-at-edge")
            if F < dimH * Fch(0.5): flags.append("F<F(Higgs)")
            if Nf <= 3 * Nc + 1: flags.append("exists-via-realmass-from-T_M[Nf+1]")
            line = (f"{Nc} {Nf} | {D:.5f} {eta:.5f} 2 {RVm:.4f} {2*D:.4f} | {CT:.4f} {3.5*dimH:.1f} {ratio:.4f} | "
                    f"{F:.5f} {dimH*Fch(0.5):.5f} | {' '.join(flags)}")
        print(line, flush=True)
        with open(OUT, "a") as f:
            f.write(line + "\n")
