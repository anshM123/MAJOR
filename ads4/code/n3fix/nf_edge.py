"""n3fix/nf_edge.py: the patched-n3 J = 26 (Delta = 1026) nodes, x = 0.85 ... 0.97 (val_n3 + edge): per kind the dominant element of
the full state integral exp(logN) E_state (normalised by the reference's dominant element at the same nu, which removes the nu-only
normalisation), its relative MC error, the max relative error over elements, the 4 structure-basis entries of E_state/|.|, and
the C3H/EE ratio of the dominant elements (EDGE_R2: C3H ~ -C2 for +-, ~ +C2 for ++ after the helicity map; here in the raw
structure basis).  Usage: python3 nf_edge.py"""
import nf_guard  # noqa: F401
import os, pickle
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
J, Ds = 26, 1026.0
al = Ds - J - 2.0; m = np.sqrt(al*(al + 2*J + 1))
recs = {}
for sub in ('val_n3', 'edge', 'edge2', 'edge3'):
    for kind in ('EE', 'C3H'):
        f = os.path.join(HERE, 'cache', sub, 'hc_pts_%s_J%d_D%g.pkl' % (kind, J, Ds))
        if os.path.exists(f):
            for k, r in pickle.load(open(f, 'rb')).items():
                recs.setdefault((kind, k[0]), {})[k[1]] = r
nus = sorted(recs.get(('EE', Ds), {}))
print(f"# J={J} Delta={Ds:g} m={m:.3f}; patched n3, M = 10, R = 8; ratio = I_state/I_ref (dominant element, complex); errors 1 sigma (replicas)")
print("#   x       kind  ratio (state/ref)              rel.err(ratio)  max rel.err(state)  state entries / |dominant|")
for nu in nus:
    out = []
    dom = {}
    for kind in ('EE', 'C3H'):
        s = recs[(kind, Ds)].get(nu); r = recs[(kind, 20000.0)].get(nu)
        if s is None or r is None:
            continue
        Es, Er = np.ravel(s['E']), np.ravel(r['E'])
        j = int(np.argmax(np.abs(Es)))
        ratio = np.exp(s['logN'] - r['logN'])*Es[j]/Er[j]
        # replica-wise ratio error (state and reference on common random numbers)
        rs = s['reps'].reshape(s['reps'].shape[0], -1)[:, j]; rr = r['reps'].reshape(r['reps'].shape[0], -1)[:, j]
        q = rs/rr; qe = q.std(ddof=1)/np.sqrt(len(q))/abs(q.mean())
        mre = np.max(np.abs(np.ravel(s['err'])))/np.max(np.abs(Es))
        dom[kind] = Es[j]*np.exp(s['logN'])
        out.append(f"  {nu/m:.4f}  {kind:4s}  {complex(ratio):+.4e}   {qe:.2e}        {mre:.2e}          "
                   + ' '.join(f"{complex(v/abs(Es[j])):+.3f}" for v in Es))
    if 'EE' in dom and 'C3H' in dom:
        out.append(f"           C3H/EE (dominant, full integrals) = {complex(dom['C3H']/dom['EE']):+.4f}")
    print('\n'.join(out))
