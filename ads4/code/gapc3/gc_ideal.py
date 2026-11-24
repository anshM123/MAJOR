"""
gc_ideal (GAP_C3.md §5) [NE]: what an IDEAL sup-norm certification could give.  Float estimates (random torus sampling, |v_i| = rho,
plus local ascent) of M3(rho) = sup |Phi - Phi* + v.v/2| and Mh(rho) = sup |g/g_env| on the Morse polydisc, for the scalar (l1_crystal,
Q_J) and the spinning EE / C3H amplitudes (w0 / l9 builders on the n1_pw point backend, even J), fed into the same Lemma L3
(gc_lemma) -> Delta_1^ideal.  Not certified (sampling); used only to show where the sup-norm route can and cannot reach.
Usage: python3 gc_ideal.py KIND J xi rho1,rho2,...   -> cache/gc_ideal.json
"""
import os, sys, json
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import gc_guard  # noqa: F401
os.environ['L1_DEG'] = '6'
SF = os.path.join(HERE, '..', 'spinfunc'); sys.path.insert(0, SF)
import numpy as np
import n1_pw as PW
import gc_node as N

DREF = 1002.0


def evaluator(kind, J, xi):
    nu, D = xi*DREF, DREF
    zs, a, L = N.morse(nu, D)
    LP = L*np.sqrt(D)
    if kind == 'S':
        def f(V):
            zv = PW.points(zs, LP, V)
            lamF, G, Q = PW.Cr.build(zv, nu, D, J, pref=True)
            h = np.exp(G - G[0])*Q
            return lamF/D, h, PW.WINDP.copy()
    else:
        mons, sol, _ = PW.S.conserved(D, J)
        _, solt, _ = PW.S.conserved(3.0 - D, J)
        if kind == 'EE':
            import w0_ee_common as W0
            bld = W0.make_build_even_EE(PW, PW.S, PW.Cr, PW.A6)
            fb = lambda zv: bld(zv, nu, D, J, mons, sol, solt)
        else:
            fb = lambda zv: PW.P9.build_spin(zv, nu, D, J, mons, sol, solt, kind)

        def f(V):
            zv = PW.points(zs, LP, V)
            lamF, Gs, pref = fb(zv)
            h = np.exp(Gs - Gs[0])*pref[0, 0]
            return lamF/D, h, PW.WINDP.copy()
    return f


def probe(kind, J, xi, rho, n=3000, seed=0):
    f = evaluator(kind, J, xi)
    rng = np.random.default_rng(seed)
    V = rho*np.exp(1j*rng.uniform(0, 2*np.pi, (n, 7)))
    Phi, h, w = f(V)
    f3 = Phi[1:] - Phi[0] + 0.5*(V**2).sum(1)
    h0 = h[0]
    # envelope normalisation for the scalar: |Q(z*)| / |P_J(cos th*)| is not available point-wise here: report h/|h0| and h0
    # local ascent on the torus phases for the 20 worst points of each quantity
    best3, besth = np.abs(f3).max(), np.abs(h[1:]).max()
    for tgt in ('f3', 'h'):
        vals = np.abs(f3) if tgt == 'f3' else np.abs(h[1:])
        idx = np.argsort(vals)[-10:]
        for i in idx:
            th = np.angle(V[i]); step = 0.3; cur = vals[i]
            for it in range(40):
                cand = th[None, :] + step*rng.normal(size=(16, 7))
                Vc = rho*np.exp(1j*cand)
                P2, h2, _ = f(Vc)
                v2 = np.abs(P2[1:] - P2[0] + 0.5*(Vc**2).sum(1)) if tgt == 'f3' else np.abs(h2[1:])
                j = int(np.argmax(v2))
                if v2[j] > cur:
                    cur, th = v2[j], cand[j]
                else:
                    step *= 0.7
            if tgt == 'f3':
                best3 = max(best3, cur)
            else:
                besth = max(besth, cur)
    return dict(kind=kind, J=J, xi=xi, rho=rho, M3=float(best3), Mh_abs=float(besth), h0=float(abs(h0)), wind=float(w.max()),
                Mh_rel=float(besth/abs(h0)))


if __name__ == '__main__':
    kind, J, xi = sys.argv[1], int(sys.argv[2]), float(sys.argv[3])
    fn = os.path.join(HERE, 'cache', 'gc_ideal.json')
    for rho in [float(a) for a in sys.argv[4].split(',')]:
        r = probe(kind, J, xi, rho)
        db = json.load(open(fn)) if os.path.exists(fn) else {}
        db['%s_%d_%.3f_%.3f' % (kind, J, xi, rho)] = r
        json.dump(db, open(fn + '.tmp', 'w')); os.replace(fn + '.tmp', fn)
        print(json.dumps(r), flush=True)
