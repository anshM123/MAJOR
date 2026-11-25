"""n3fix/nf_compare.py: patched-n3 vs n1 at common nodes (state and Delta_ref = 2e4), per kind: full integrals exp(logN) E compared
element-wise in units of the combined MC error, the dominant element's ratio, the state/reference ratio per method (what Gamma uses),
and the scalar on the same samples (exp(logN + lqsc) Es for n1, exp(logNs) Es for n3).  Also the ORIGINAL n3 (EDGE_R) where it exists.
Usage: python3 nf_compare.py J Delta [n1dir n3dir]  -> stdout (and cache/cmp_J<J>.txt via run script redirection)"""
import nf_guard  # noqa: F401
import os, sys, pickle
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
J = int(sys.argv[1]); Ds = float(sys.argv[2])
d1 = sys.argv[3] if len(sys.argv) > 3 else 'val_n1'
d3 = sys.argv[4] if len(sys.argv) > 4 else 'val_n3'
al = Ds - J - 2.0; m = np.sqrt(al*(al + 2*J + 1))


def load(d, kind):
    f = os.path.join(HERE, 'cache', d, 'hc_pts_%s_J%d_D%g.pkl' % (kind, J, Ds))
    return pickle.load(open(f, 'rb')) if os.path.exists(f) else {}


def full(r, scalar=False):
    if scalar:
        ln = r['logNs'] if r.get('logNs') is not None else r['logN'] + r['lqsc']
        return np.exp(ln)*r['Es'], np.exp(ln.real)*abs(r['errs'])
    return np.exp(r['logN'])*np.ravel(r['E']), np.exp(r['logN'].real)*np.abs(np.ravel(r['err']))


print(f"# J={J} Delta={Ds:g}: n1 ({d1}) vs patched n3 ({d3});  dev = |I_n3 - I_n1| / sqrt(err_n1^2 + err_n3^2) per element (max over elements)")
for kind in ('EE', 'C3H'):
    A, B = load(d1, kind), load(d3, kind)
    for key in sorted(k for k in A if k in B and k[0] == Ds):
        nu = key[1]; kr = (20000.0, nu)
        line = f"{kind:4s} x={nu/m:.4f}:"
        rat = {}
        for tag, k in (('state', key), ('ref', kr)):
            if k not in A or k not in B:
                continue
            i1, e1 = full(A[k]); i3, e3 = full(B[k])
            dev = np.abs(i3 - i1)/np.sqrt(e1**2 + e3**2)
            j = int(np.argmax(np.abs(i1)))
            line += (f" {tag}: dev max {dev.max():.2f} | dom n3/n1 {complex(i3[j]/i1[j]):.4f} rel.err n1 {e1[j]/abs(i1[j]):.1e} "
                     f"n3 {e3[j]/abs(i3[j]):.1e};")
            rat[tag] = (i1, i3)
            if kind == 'EE':
                s1, se1 = full(A[k], True); s3, se3 = full(B[k], True)
                line += f" scalar n3/n1 {complex(s3/s1):.4f} (dev {abs(s3 - s1)/np.hypot(se1, se3):.1f});"
        if len(rat) == 2:
            q1 = rat['state'][0]/rat['ref'][0]; q3 = rat['state'][1]/rat['ref'][1]
            j = int(np.argmax(np.abs(rat['state'][0])))
            line += f" state/ref n3 vs n1: {complex(q3[j]/q1[j]):.4f}"
        print(line)
# the original (unpatched) n3 values of EDGE_R, for reference
ER = os.path.join(HERE, '..', 'edger', 'cache')
for sub in ('', 'n3ov'):
    for kind in ('C3H',):
        f = os.path.join(ER, sub, 'hc_pts_%s_J%d_D%g.pkl' % (kind, J, Ds))
        if not os.path.exists(f):
            continue
        db = pickle.load(open(f, 'rb'))
        for key, r in sorted(db.items()):
            if key[0] == Ds and 'n3' in r.get('method', ''):
                i, e = full(r)
                print(f"original n3 (EDGE_R, {sub or 'main'}) {kind} x={key[1]/m:.4f}: |I| max {np.abs(i).max():.3e}, rel err {e.max()/np.abs(i).max():.2f}")
