"""
gc_spin (GAP_C3.md §4.3): certified sup norms of the SPINNING amplitudes (EE k = 2: w0_ee_common even builder; C3H k = 3: l9_points
build_spin) on the Morse polydisc, by running the project builders unchanged on the gc_ball Taylor-model backend, combined with the
certified phase data of the log-factored integrand (gc_cert; Phi_xi is kind-independent) in Lemma L3.

Single-structure spins only (J = 0, 2: na = nc = 1).  The <TTO> coefficient vectors are the float ones of l5_spin.conserved at
Delta = 1002 (frozen): the certificate is for that amplitude; their Delta-dependence is not covered (GAP_C3 §4.3).
Envelope normalisation: h = g / (g(z*)/D_++(xi)), D_++ = P_J(1 - 2 xi^2) (the LO shape of the chain's Gamma).
Usage: GC_DEG=5 python3 gc_spin.py KIND J xi1,xi2 rho1,rho2  -> cache/gc_spin.json
"""
import os, sys, json, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import gc_guard  # noqa: F401
import numpy as np
from scipy.special import eval_legendre
import gc_cert as C
import gc_node as N
import gc_ball as B

DREF = 1002.0


def spin(kind, J, xi, rho):
    t0 = time.time()
    ph, _ = C.cert(0, xi*DREF, DREF, rho)               # phase data (kind- and J-independent part)
    badph = ph['bad']
    R = N.build(kind, J, xi*DREF, DREF, rho)
    eG = B.exp(R['G'])
    g = B.mul(eG, R['prefs'][0, 0])
    tot, _ = B.bound(g)
    Mg = float(np.asarray(tot).ravel()[0].upper())
    c0 = abs(complex(g[0].mid())); c0lo = float(abs(g[0]).lower())
    Dh = abs(eval_legendre(J, 1 - 2*xi*xi))
    env = c0lo/max(Dh, 1e-300)                            # |g(z*)|/|D(xi)| (lower bound -> Mh upper bound)
    out = dict(ph, kind=kind, J=J, Mg_amp=Mg, c0=c0, Dh=Dh, bad=badph + len(B.BAD), qamp=B.WMAX[0],
               Bh=[Mg/env] + [0.0]*B.DEG, remh=0.0, Mg=Mg/env, h0=Dh, t=time.time() - t0)
    return out


if __name__ == '__main__':
    kind, J = sys.argv[1], int(sys.argv[2])
    xis = [float(a) for a in sys.argv[3].split(',')]
    rhos = [float(a) for a in sys.argv[4].split(',')]
    fn = os.path.join(HERE, 'cache', 'gc_spin.json')
    for xi in xis:
        for rho in rhos:
            try:
                o = spin(kind, J, xi, rho)
            except Exception as e:
                print('ERR', kind, J, xi, rho, repr(e)[:200], flush=True)
                continue
            db = json.load(open(fn)) if os.path.exists(fn) else {}
            db['%s_%d_%.4f_%.4f_%d' % (kind, J, xi, rho, B.DEG)] = o
            json.dump(db, open(fn + '.tmp', 'w')); os.replace(fn + '.tmp', fn)
            print(f"{kind} J={J} xi={xi} rho={rho} deg={B.DEG}: bad={o['bad']} qamp={o['qamp']:.3f} Mh_env={o['Mg']:.4g} "
                  f"S3={o['S3']:.3g} t={o['t']:.1f}s", flush=True)
