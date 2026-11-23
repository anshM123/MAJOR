"""agproof/ag_exact.py (AG_PROOF.md §3): the EXACT scalar (Delta_phi = 3) k = 2 heavy moment at COMPLEX nu, on vertical lines.

The Arb Laplace/Barnes form of exactskel (h13_laplace.Laplace3H, used read-only through h13_exact, which patches in the FLINT Mack
polynomial) evaluates
    m^4 C_{2,nu} = pre * sum_{m < Mx} g_m(nu) R_m,   g_m(nu) = prod_{l<m} ((A/2 + l)^2 + nu^2/4)/((A + l)(l + 1)) / Gamma(A),
with R_m independent of nu.  Each g_m is a polynomial in nu^2, so the same finite sum continues the moment to complex nu (here with
acb arithmetic).  |g_m(nu)| <= g_m(|nu|), so the truncation rule of Laplace3H (M from numax) is applied with numax = max |nu|; the
truncation is checked by doubling M at one state (flag TRUNC).

Output, per gap state J (Delta = 1000 + J, Delta_phi = 3) and vertical line nu_r + i y (nu_r = p_r Lam, p_r in PR, y in [0, Ymax]):
  ex2(nu), the flat part F(nu) = E_2(x) P_J(1 - 2x^2) (x = nu/m, complex), the curvature part 2(P_{-1/2+i nu}(cosh beta) - J0(bhat p)),
and Arb radii.  The real-axis values are checked against spinfunc/cache/sm_exact_1000.pkl (read-only).
-> cache/exact_lines_<J>.pkl
Usage: python3 ag_exact.py J [J ...]       (env AG_NCONF="N0:nq:Lf" (default 320:960:0.5 = sm_exact's reported configuration),
                                            AG_PREC (bits, default 400), AG_TRUNC=1 doubles M (truncation check))
"""
import os, sys, time, pickle
sys.dont_write_bytecode = True
os.environ.setdefault('OMP_NUM_THREADS', '1')
HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, '..'))
CACHE = os.path.join(HERE, 'cache')
sys.path.insert(0, os.path.join(CODE, 'exactskel'))
import numpy as np
from flint import arb, acb, ctx
import h13_exact as H                     # read-only: patches mack.mack_P_fast -> FLINT/mat Mack polynomial; Laplace3 = Laplace3H

Dg = 1000
LAM = float(np.sqrt((Dg - 2)*(Dg - 1)))
PR = (1.0, 0.95, 0.9, 0.8, 0.6, 0.4)          # vertical lines at nu_r = p_r Lam


def state(J):
    Delta = Dg + J
    al = Delta - J - 2.0; ab = al + 2*J + 1
    m2 = al*ab; m = np.sqrt(m2); eta = (al*al + ab*ab)/(2*al*ab)
    return dict(J=J, Delta=float(Delta), m=m, m2=m2, eta=eta, beta=float(np.arccosh(eta)), mfl=m/LAM)


def build(st, numax, prec, conf, trunc=False):
    N0, nq, Lf = conf
    Ld = max(4.0*N0, 0.5*st['Delta']**2/max(2*3.0 - 1.5, 0.5))*Lf
    M = None
    if trunc:
        M = 2*int(3*numax**2/(2*st['Delta']) + 6*np.sqrt(numax + 1) + 60)
    L3 = H.Laplace3(st['Delta'], st['J'], 3.0, ks=(2,), numax=float(numax), prec=prec, N0=N0, nq=nq, L=Ld, M=M)
    L3.prepare(2)
    return L3


def ex2_complex(L3, nu):
    """m^4 C_{2,nu} at complex nu (acb), = L3.value_fast(2, nu) continued."""
    ctx.prec = L3.prec
    k = 2
    A = 2*L3.Dfa + k - arb(1.5)
    nu = acb(complex(nu).real, complex(nu).imag)
    R = L3.R[k]
    I = acb(0)
    gm = acb(1/A.gamma())
    nu2q = nu*nu/4
    for mm in range(len(R)):
        I += gm*R[mm]
        gm = gm*((A/2 + mm)**2 + nu2q)/((A + mm)*(mm + 1))
    if not hasattr(L3, '_pre_c'):
        from exactmom import log_q
        lq = arb(float(log_q(L3.Delta, L3.J, L3.Df)))
        aDf = 4*L3.Dfa.gamma()**2*(L3.Dfa - arb(0.5)).gamma()**2
        L3._pre_c = arb(2)**(1 - 2*k)*arb.pi()**arb(1.5)*aDf*L3.m2**k*(-lq).exp()
    return ((-1)**L3.J)*L3._pre_c*I


def legendre_c(J, z):
    """P_J(z) for complex z (upward recurrence; |z| = O(1), J <= 400: stable for the growing solution)."""
    p0 = np.ones_like(z); p1 = z.copy()
    if J == 0:
        return p0
    for n in range(1, J):
        p0, p1 = p1, ((2*n + 1)*z*p1 - n*p0)/(n + 1)
    return p1


def flat_c(st, nu):
    x = nu/st['m']
    E2 = (2 - x*x)/(1 - x*x)**2
    return E2*legendre_c(st['J'], 1 - 2*x*x)


def run(J):
    st = state(J)
    conf = tuple(float(v) if i == 2 else int(v) for i, v in enumerate(os.environ.get('AG_NCONF', '320:960:0.5').split(':')))
    prec = int(os.environ.get('AG_PREC', '400'))
    trunc = bool(os.environ.get('AG_TRUNC'))
    beta = st['beta']
    Ymax = 25.0/(0.97*beta)                    # Lemma P4W's rectangle height T (nu units)
    numax = float(np.hypot(LAM, Ymax))*1.02
    t0 = time.time()
    L3 = build(st, numax, prec, conf, trunc)
    tb = time.time() - t0
    lines = {}
    ny = 81
    for pr_ in PR:
        ys = np.linspace(0.0, Ymax, ny)
        nus = pr_*LAM + 1j*ys
        vals = [ex2_complex(L3, v) for v in nus]
        ex = np.array([complex(float(v.real.mid()), float(v.imag.mid())) for v in vals])
        rad = np.array([float(abs(v).rad())/max(float(abs(v).mid()), 1e-300) for v in vals])
        lines[pr_] = dict(ys=ys, nus=nus, ex=ex, rad=rad, flat=flat_c(st, nus))
    # real-axis check against the sm_exact cache
    chk = None
    try:
        db = pickle.load(open(os.path.join(CODE, 'spinfunc', 'cache', 'sm_exact_1000.pkl'), 'rb'))
        key = {40: (40, 0.0, 3.0, 'npf2'), 48: (48, 0.0, 3.0, 'd48'), 64: (64, 0.0, 3.0), 100: (100, 0.0, 3.0)}.get(J)
        if key in db:
            r = db[key]
            sel = np.arange(0, r['npn'], max(1, r['npn']//12))
            mine = np.array([float(ex2_complex(L3, float(r['nus'][i])).real.mid()) for i in sel])
            ref = np.asarray(r['ex'][2])[sel]
            chk = float(np.max(np.abs(mine - ref)/np.maximum(np.abs(ref), 1e-3)))
    except Exception as e:                                   # noqa: BLE001
        chk = 'failed: %r' % (e,)
    rec = dict(st=st, conf=conf, prec=prec, trunc=trunc, M=L3.M, Ymax=Ymax, lines=lines, chk_real=chk,
               t_build=tb, t_total=time.time() - t0)
    fn = os.path.join(CACHE, 'exact_lines_%d%s.pkl' % (J, '_trunc' if trunc else ''))
    pickle.dump(rec, open(fn + '.tmp', 'wb')); os.replace(fn + '.tmp', fn)
    print(f"J={J} m/L={st['mfl']:.4f} beta={beta:.4f} M={L3.M} build {tb:.0f}s total {rec['t_total']:.0f}s; real-axis check vs "
          f"sm_exact: {chk}; max rel rad {max(float(np.max(l['rad'])) for l in lines.values()):.1e}", flush=True)
    return rec


if __name__ == '__main__':
    for a in sys.argv[1:]:
        run(int(a))
