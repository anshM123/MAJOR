"""
n3_exact (doc 25): EXACT scalar comparison object at small nu (5 <= nu <= 20), without the local (Taylor) shift of l2.Xexact.

CLS-D's literal Psi_{2,nu} on a scalar block is N(nu) X(nu), X(nu) = g(nu) C(nu + i) + g(-nu) C(nu - i) (CHMRS magic block, l2).
l2.Xexact shifts C by a third-order Taylor series with the LEADING-ORDER derivatives; at nu = O(5-20) the shift i/m is 5-20 % of x and
the remainder's own shift is not small.  Here C(nu +- i) is evaluated EXACTLY up to an interpolation of the smooth remainder:
   m^4 C_2(nu) = E_2(x) P_J(1 - 2x^2) + 2[P_{-1/2+i nu}(cosh beta) - J_0((2J+1)x)] + E_2(x) rho_2(nu)        (e1 definition)
with every term except rho_2 analytic in closed form (evaluated at complex nu), and rho_2 (exact on the 64 Gauss nodes of e1) continued
to nu +- i by barycentric interpolation in s = p^2 (p = nu/Lam; rho_2 is even in nu) on the NLOC nodes nearest to the target.
interp_check(): the same interpolation on half the nodes, compared at the other half (accuracy statement).
"""
import os, sys, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np, mpmath as mp
from scipy.special import eval_legendre, jv
from scipy.interpolate import BarycentricInterpolator
import l2_scalar_validate as V

NLOC = int(os.environ.get('N3_NLOC', '24'))
_DB = {}


def db(Dg=100, Df=3):
    if (Dg, Df) not in _DB:
        _DB[(Dg, Df)] = pickle.load(open(os.path.join(V.EC, 'e1_%d_Df%d.pkl' % (Dg, Df)), 'rb'))
    return _DB[(Dg, Df)]


def Lam(Dg):
    return np.sqrt((Dg - 2.0)*(Dg - 1.0))


def legP(nu, eta):
    """P_{-1/2 + i nu}(eta), complex nu."""
    a, b = mp.mpf(1)/2 + 1j*mp.mpc(nu), mp.mpf(1)/2 - 1j*mp.mpc(nu)
    return complex(mp.hyp2f1(a, b, 1, (1 - mp.mpf(eta))/2))


def hyb(nu, r):
    """E_2 P_J + 2(P_{-1/2+i nu}(eta) - J_0((2J+1)x)), complex nu."""
    J, m = r['J'], r['m']; x = nu/m
    E2 = (2 - x*x)/(1 - x*x)**2
    PJ = complex(mp.legendre(J, 1 - 2*mp.mpc(x)**2))
    return E2*PJ + 2*(legP(nu, r['eta']) - complex(mp.besselj(0, (2*J + 1)*mp.mpc(x))))


def rho_interp(nu, r, Lm, nodes=None):
    """rho_2 at complex nu by barycentric interpolation in s = (nu/Lam)^2 on the NLOC nodes nearest to Re nu."""
    p = r['nus']/Lm; rho = r['rho'][2]
    idx = np.arange(p.size) if nodes is None else np.asarray(nodes)
    near = idx[np.argsort(np.abs(p[idx] - (nu.real/Lm)))[:NLOC]]
    near = np.sort(near)
    s = p[near]**2
    sc = (s - s.min())/(s.max() - s.min())
    f = BarycentricInterpolator(sc, rho[near])
    st = (((nu/Lm)**2) - s.min())/(s.max() - s.min())
    # BarycentricInterpolator accepts complex evaluation points via real/imag linearity? evaluate the Lagrange form directly
    w = f.wi
    d = st - sc
    if np.any(d == 0):
        return rho[near][np.argmin(np.abs(d))]
    t = w/d
    return np.sum(t*rho[near])/np.sum(t)


def C2(nu, r, Lm, nodes=None):
    x = nu/r['m']; E2 = (2 - x*x)/(1 - x*x)**2
    return (hyb(nu, r) + E2*rho_interp(nu, r, Lm, nodes))/r['m']**4


def X(nu, r, Dg=100, nodes=None):
    """exact X(nu) = g C(nu+i) + g~ C(nu-i) and its pieces (CHMRS-type C at the node for reference)."""
    Lm = Lam(Dg)
    gp, gm = V.gfac(nu)
    Cp, Cm = C2(nu + 1j, r, Lm, nodes), C2(nu - 1j, r, Lm, nodes)
    return gp*Cp + gm*Cm, dict(gp=gp, gm=gm, Cp=Cp, Cm=Cm)


def interp_check(Dg=100, states=((1.0, 0), (1.05, 4), (1.1, 10), (1.2, 6), (1.3, 10)), pr=print):
    """accuracy of the rho_2 continuation: interpolate from the even nodes, evaluate at the odd ones (real nu), and compare
    C(nu +- i) from two disjoint node sets (even/odd) at the target nus (complex: the actual use)."""
    d = db(Dg); Lm = Lam(Dg)
    for key in states:
        if key not in d:
            continue
        r = d[key]; nus = r['nus']; ev = np.arange(0, nus.size, 2); od = np.arange(1, nus.size, 2)
        errs = []; cerr = []
        for i in od:
            if 3.0 <= nus[i] <= 25.0:
                errs.append(abs(rho_interp(nus[i] + 0j, r, Lm, ev) - r['rho'][2][i]))
        for i in range(nus.size):
            if 4.0 <= nus[i] <= 22.0:
                for sgn in (1, -1):
                    a = rho_interp(nus[i] + sgn*1j, r, Lm, ev); b = rho_interp(nus[i] + sgn*1j, r, Lm, od)
                    cerr.append(abs(a - b))
        pr(f"   state m/Lam={key[0]} J={key[1]} (m={r['m']:.2f}): rho_2 leave-half-out at real nodes nu in [3,25]: max |err| "
           f"{max(errs):.2e} (|rho_2| max {np.abs(r['rho'][2][(nus > 3) & (nus < 25)]).max():.2e}); even- vs odd-node continuation to "
           f"nu +- i: max diff {max(cerr):.2e}  [in c2 units x m^2: {max(cerr)*r['m']**2:.3f}]")


if __name__ == '__main__':
    interp_check()
