"""
nf_nodes (NEARGAP_FILL.md §2): direct vs K2 on the own-mass row Delta = 1000 + J, node by node, in Gamma (c2) (the chain's records:
w2/q3 conventions, as neargap2/ng_val), with the jackknife sigma of the direct value over the 8 RQMC replicas; plus the small-nu
direct values (nu ~ 25 and 50, x = 0.025 / 0.05) that the NU_LOW = 25 variant uses.
Usage: python3 nf_nodes.py J  -> cache/nf_nodes_J<J>.txt
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import nf_core as C
import numpy as np
H = C.H

if __name__ == '__main__':
    J = int(sys.argv[1]); D = 1000.0 + J
    gK = H.grid(J)[D]
    gD = H.grid(J, inject=[D])[D]
    reps = [H.grid(J, inject=[D], jk=j)[D] for j in range(8)]
    keys = ['C2', 'C3H'] if J % 2 else [(kd, h) for kd in ('C2', 'C3H') for h in ('++', '+-')]
    lines = [f"# nf_nodes J={J} row {D:g} (m = {C.mphys(D, J):.3f}): Gamma (c2) K2 vs direct +- jk; K2 nodes {len(gK['xs'])}, direct-row nodes {len(gD['xs'])}"]
    dd, zz, mx = [], [], 0.0
    for i, x in enumerate(gD['xs']):
        s = f"   x={x:.4f} nu={x*gD['m']:7.1f}"
        j0 = int(np.argmin(np.abs(gK['xs'] - x))); on = abs(gK['xs'][j0] - x) < 1e-6
        for k in keys:
            vj = np.array([r[k][int(np.argmin(np.abs(r['xs'] - x)))] for r in reps]); e = float(np.sqrt(7/8*((vj - vj.mean())**2).sum()))
            v = float(gD[k][i]); mx = max(mx, abs(v))
            if on:
                d = v - float(gK[k][j0]); dd.append(abs(d)); zz.append(abs(d)/max(e, 1e-9))
                s += f" | {k}: K2 {gK[k][j0]:+8.1f} dir {v:+8.1f} +- {e:5.1f}"
            else:
                s += f" | {k}: dir {v:+8.1f} +- {e:5.1f}"
        lines.append(s)
    lines.append(f"# direct - K2 on the {len(dd)//len(keys)} common nodes: median |d| {np.median(dd):.0f} c2, max {np.max(dd):.0f}; median |d|/sigma "
                 f"{np.median(zz):.1f}, max {np.max(zz):.1f}; max |Gamma_direct| on the row {mx:.0f} c2")
    print("\n".join(lines))
    open(os.path.join(C.NFC, 'nf_nodes_J%d.txt' % J), 'w').write("\n".join(lines) + "\n")
