"""
ee1000_rawcmp (EE_DG1000_LARGEJ.md §3): K2 convergence test without a full row: at a state point (Delta, nu), the all-orders direct state
matrix M_dir (ee1000_direct, l9 normalisation) vs the K2 Laplace sum M0 + M2 + M4 (and M0, M0 + M2), conserved-structure basis.  The
relative deviation |M_K - M_dir|/max|M_dir| times m^2 is the size of the K-truncation error in Gamma (c2 units) up to the O(1) shape factor
(|D_h| <= 1).  Also the order ratios |M2|/|M0|, |M4|/|M0| (max element).
Usage: python3 ee1000_rawcmp.py "J:Delta,..." [kind]   -> cache/ee1000_rawcmp.txt
"""
import os, sys, glob, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import ee1000_chain as CH
import ee1000_direct as ED
import n1_analysis as N1
for fn in sorted(glob.glob(os.path.join(CH.C, '*shadow_k2*.pkl'))):
    N1.SHK2.update(pickle.load(open(fn, 'rb')))
kind = sys.argv[2] if len(sys.argv) > 2 else 'EE'
lines = []
def pr(s=''):
    print(s, flush=True); lines.append(s)
pr(f"# ee1000_rawcmp ({kind}): state matrices, direct (all orders) vs Laplace K0/K1/K2; deviations relative to max|M_dir|, and x m^2")
for spec in sys.argv[1].split(','):
    J, D = int(spec.split(':')[0]), float(spec.split(':')[1])
    R = CH.ee_even(J) if kind == 'EE' else CH._q2_even(kind, J)
    db = ED.dload(J, kind)
    for (Dk, nu), r in sorted(db.items()):
        if Dk != D:
            continue
        k = [kk for kk in R if kk[0] == D and abs(kk[1] - nu) < 1e-7*nu]
        if not k:
            continue
        st = R[k[0]]
        H, l = N1.Hrep(r, dict(Ti=np.eye(2)), D, J)
        Md = np.exp(l - (st['logC'] - st['logK']))*H.mean(0)
        sd = np.exp(l - (st['logC'] - st['logK']))*H.std(0)/np.sqrt(len(H))
        sc = np.abs(Md).max(); m = st['m']
        dev = [np.abs(Mk - Md).max()/sc for Mk in (st['M0'], st['M0'] + st['M2'], st['M0'] + st['M2'] + st['M4'])]
        pr(f"  J={J} Delta={D:.0f} x={nu/m:.3f} (J^2/m={J*J/m:.2f}): |M_K - M_dir|/max|M_dir| K0 {dev[0]:.2e} K1 {dev[1]:.2e} K2 {dev[2]:.2e}"
           f"  (x m^2: {dev[2]*m*m:.0f}; RQMC {np.abs(sd).max()/sc:.1e} = {np.abs(sd).max()/sc*m*m:.0f});  |M2|/|M0| {np.abs(st['M2']).max()/np.abs(st['M0']).max():.3f}"
           f" |M4|/|M0| {np.abs(st['M4']).max()/np.abs(st['M0']).max():.3f}")
open(os.path.join(CH.C, 'ee1000_rawcmp_%s.txt' % kind), 'w').write("\n".join(lines) + "\n")
