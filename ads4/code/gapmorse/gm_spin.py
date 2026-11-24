"""
gm_spin (GAP_MORSE.md §4.3): certified Taylor-model data of the SPINNING graviton amplitudes (EE k = 2, w0_ee_common even builder;
C3H k = 3, l9_points.build_spin) on the Hermitian ball B_rho, in LOG-FACTORED form:
    g = exp( sum_k p_k Log(f_k/f_k*) - 4 Log(v042/v042*) ) x [num/num*  (EE only: e^{lf} = num/(m13 m24))] x pref(z)
with the same elementary factors f_k as gm_scalar (meas, xsq, FP, dual, cel of l1_crystal; the 3pt exponents of l9/w0:
-(10-J)/2 m12, -J/2 m15, m25, -(7-J)/2 m34, -(3+J)/2 m35, m45), and pref = the project's prefactor (E (x) E or q o_H contracted with the
<TTO> samples; l5_spin / w0 / l9 functions on the gm_ball Taylor-model backend, unchanged).  The <TTO> vectors are the float null
vectors of l5_spin.conserved at Delta = 1002 (frozen, as GAP_C3 §4.3).  Envelope: h = g / (|g(z*)| / |P_J(1 - 2 xi^2)|).
Validated against the project builders run on the same backend (exact Taylor coefficients agree; `check`).
Usage: GC_DEG=8 python3 gm_spin.py KIND J xi rho1,rho2   -> cache/gm_spin.json
"""
import os, sys, json, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import gm_guard  # noqa: F401
import numpy as np
from scipy.special import eval_legendre
import gm_ball as B
import gm_scalar as GS
import gm_cert as C
from flint import acb, arb

DREF = 1002.0
up = C.up


def amp_tm(kind, J, xi, rho, Delta=DREF):
    B.BAD.clear(); B.WMAX[0] = 0.0; B.RHO[0] = rho
    nu = xi*Delta
    zs, a, L = C.N.morse(nu, Delta)
    LP = L*np.sqrt(Delta)
    zc, aa = C.N.saddle_ball(nu, Delta); xib = acb(arb(nu)/arb(Delta))
    zv = B.linsub_ball(zc, LP)
    K = GS.TMB(B)
    snap = {}

    def mark():
        snap['bad'] = len(B.BAD)
    r = GS.build(zv, nu, Delta, 0, K, acb(aa), xi=xib, mark=mark)
    nb0 = len(B.BAD)
    geo, lg, P = r['geo'], r['lg'], dict(r['P'])
    x0, x1, x2, x3, x4, e, z0 = (geo[k] for k in ('x0', 'x1', 'x2', 'x3', 'x4', 'e', 'z0'))
    S = B.S
    for k_ in ('m12', 'm15', 'm25', 'm34', 'm35', 'm45'):
        P.pop(k_, None)
    P.update({'m12': -(10 - J)/2, 'm15': -J/2, 'm25': -J/2, 'm34': -(7 - J)/2, 'm35': -(3 + J)/2, 'm45': -(3 + J)/2})
    if kind == 'EE':
        P['m13'] = P.get('m13', 0.0) - 1.0
        P['m24'] = P.get('m24', 0.0) - 1.0
    Gl = None
    for k, p in P.items():
        if p == 0:
            continue
        f = p*lg(k)
        Gl = f if Gl is None else Gl + f
    v042 = S.Vfun(x0, z0, x4, x2)
    Gl = Gl - 4*K.log1(K.rel(v042))
    grel = B.exp(Gl)
    if kind == 'EE':
        grel = B.mul(grel, K.norm(r['num']))
    mons, sol, _ = S.conserved(Delta, J)
    _, solt, _ = S.conserved(3.0 - Delta, J)
    W2, WJ = S.null_circle(2), S.null_circle(J)
    x10 = [x1[i] - x0[i] for i in range(3)]; x30 = [x3[i] - x0[i] for i in range(3)]
    F = S.tto_samples(x1, x10, x2, W2, e, WJ, J, mons, sol)
    Gt = S.tto_samples(x3, x30, x4, W2, e, WJ, J, mons, solt)
    F = S.fourier(S.fourier(F, 1, 2), 2, J); Gt = S.fourier(S.fourier(Gt, 1, 2), 2, J)
    if kind == 'EE':
        Eb = B.W0.e_factor(B, S, x0, z0, x2, x4, W2)
        pref = B.W0._contract(S, B, S.bm(Eb, Eb), F, Gt, J)
    else:                                            # l9_points.build_spin, kind C3H (lines of the 'C3*' branch, unchanged)
        H20 = S.lin_samples([S.Hfun(x2, S.BASIS[m], x0, z0) for m in range(3)], W2)
        H40 = S.lin_samples([S.Hfun(x4, S.BASIS[m], x0, z0) for m in range(3)], W2)
        H24b = np.array([[S.Hfun(x2, S.BASIS[m], x4, S.BASIS[k]) for k in range(3)] for m in range(3)])
        H24 = np.einsum('im,jk,mkn->ijn', W2, W2, H24b)
        inv = B.recip(B.mul(S.Xfun(x2, x0), S.Xfun(x4, x0)))
        q = S.bm(2*S.bm(H20[:, None], H40[None, :]) + S.bm(S.bm(v042, v042)[None, None], H24), inv)
        D = S.bm(q, B.P9.odd_structure(x0, z0, x2, x4, W2, kind[-1]))
        pref = B.W0._contract(S, B, D, F, Gt, J)
    assert pref.shape[:2] == (1, 1)
    g = B.mul(grel, pref[0, 0])
    return dict(g=g, grel=grel, pref=pref[0, 0], dPhi=r['dPhi'], badPhi=snap['bad'], bad=len(B.BAD), nb0=nb0, q=B.WMAX[0], zv=zv,
                mons=mons, sol=sol, solt=solt)


def check(kind, J, xi, rho=0.004):
    """exact Taylor coefficients of g/g(z*) (centre half) against the project builder on the same backend."""
    R = amp_tm(kind, J, xi, rho)
    zv = R['zv']
    if kind == 'EE':
        bEE = B.W0.make_build_even_EE(B, B.S, B.Cr, B.A6)
        lamF, Gs, prefs = bEE(zv, xi*DREF, DREF, J, R['mons'], R['sol'], R['solt'])
    else:
        lamF, Gs, prefs = B.P9.build_spin(zv, xi*DREF, DREF, J, R['mons'], R['sol'], R['solt'], kind)
    go = B.mul(B.exp(Gs), B.obj(prefs)[0, 0])
    a = R['g'][:B.NM]/R['g'][0]; b = go[:B.NM]/go[0]
    err = max(float(abs(a[i] - b[i]).upper())/max(1.0, float(abs(b[i]).upper())) for i in range(B.NM))
    return err


def cert(kind, J, xi, rho):
    t0 = time.time()
    R = amp_tm(kind, J, xi, rho)
    g = R['g']
    Dh = abs(eval_legendre(J, 1 - 2*xi*xi))
    env = arb(float(abs(g[0]).lower()))/arb(Dh)
    bh = B.bombieri(g[:B.NM][None])[0]
    return dict(kind=kind, J=J, xi=xi, rho=rho, deg=B.DEG, norm=B.NORM, q=R['q'], bad=R['bad'], badPhi=R['badPhi'], Dh=Dh,
                bh=[float((bh[k]/env).upper()) for k in range(B.DEG + 1)], remh=float((abs(g[B.NM])/env).upper()),
                h0=float((abs(g[0])/env).upper()), t=time.time() - t0)


if __name__ == '__main__':
    kind, J = sys.argv[1], int(sys.argv[2])
    if sys.argv[3] == 'check':
        for xi in (0.3, 0.5):
            print(kind, J, xi, 'max rel. coefficient difference vs project builder:', check(kind, J, xi), flush=True)
        sys.exit(0)
    xis = [float(a) for a in sys.argv[3].split(',')]
    rhos = [float(a) for a in sys.argv[4].split(',')]
    fn = os.path.join(HERE, 'cache', 'gm_spin.json')
    for xi in xis:
        for rho in rhos:
            key = '%s_%d_%.4f_%.4f_%d' % (kind, J, xi, rho, B.DEG)
            db = json.load(open(fn)) if os.path.exists(fn) else {}
            if key in db:
                continue
            try:
                o = cert(kind, J, xi, rho)
            except Exception as ex:
                o = dict(kind=kind, J=J, xi=xi, rho=rho, deg=B.DEG, bad=-1, err=repr(ex)[:200])
            db = json.load(open(fn)) if os.path.exists(fn) else {}
            db[key] = o
            json.dump(db, open(fn + '.tmp', 'w')); os.replace(fn + '.tmp', fn)
            print(f"{kind} J={J} xi={xi} rho={rho} deg={B.DEG}: bad={o['bad']} q={o.get('q', -1):.3f} remh={o.get('remh', -1):.3g} "
                  f"bh={[round(x, 2) for x in o.get('bh', [])]} t={o.get('t', 0):.1f}s", flush=True)
