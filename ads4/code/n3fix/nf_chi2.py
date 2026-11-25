"""n3fix/nf_chi2.py: pooled element-wise deviations (patched n3 vs n1, full integrals, units of combined MC error; real and imaginary
parts separately) over all validation nodes, per J, kind, state / reference.  Usage: python3 nf_chi2.py"""
import nf_guard  # noqa: F401
import os, pickle
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
for J, Ds in ((26, 1026.0), (29, 1029.0)):
    for kind in ('EE', 'C3H'):
        ld = lambda d: pickle.load(open(os.path.join(HERE, 'cache', d, 'hc_pts_%s_J%d_D%g.pkl' % (kind, J, Ds)), 'rb'))
        A, B = ld('val_n1'), ld('val_n3')
        for tag, sel in (('state', lambda k: k[0] == Ds), ('ref', lambda k: k[0] == 20000.0)):
            z = []
            for k in A:
                if k in B and sel(k):
                    i1 = np.exp(A[k]['logN'])*np.ravel(A[k]['E']); i3 = np.exp(B[k]['logN'])*np.ravel(B[k]['E'])
                    # per-part errors from the replicas
                    for part in (np.real, np.imag):
                        r1 = part(np.exp(A[k]['logN'])*A[k]['reps'].reshape(A[k]['reps'].shape[0], -1))
                        r3 = part(np.exp(B[k]['logN'])*B[k]['reps'].reshape(B[k]['reps'].shape[0], -1))
                        e = np.sqrt(r1.var(axis=0, ddof=1)/r1.shape[0] + r3.var(axis=0, ddof=1)/r3.shape[0])
                        ok = np.abs(part(i1)) > 3*e                  # only parts that are resolved
                        z += list((part(i3) - part(i1))[ok]/e[ok])
            z = np.array(z)
            print(f"J={J} {kind:4s} {tag:5s}: n = {z.size:3d}, mean z {z.mean():+.2f}, rms z {np.sqrt(np.mean(z**2)):.2f}, max |z| {np.abs(z).max():.2f}, "
                  f"fraction |z| > 2: {np.mean(np.abs(z) > 2):.2f}")
