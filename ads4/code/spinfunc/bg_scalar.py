"""
bg_scalar (BELOW_FOLD2.md §2): analytic-radial 1/Delta expansion of the scalar CLS-D crystal integral (k = 2, Psi_{2,nu} on a scalar
block G_{Delta,J}) on the real section D.

  I(nu; Delta, J) = const(Delta, J) * sum_{k <= K} Delta^{-k} sum_{n <= k/2} c_n(J) M_{k,n}(nu),
  c_n(J) = P_J^{(n)}(1)/n! = (J+n)!/((J-n)! 2^n n!^2)            (spin: P_J(t), t = 1 + delta, delta = O(a^2))
  M_{k,n}(nu) = sum_j (i nu/2)^j int_D [lead_nu] e^{-S} [a^k coefficient of corr_j * delta^n * B-correction]  (closed form, bg_ser)
The M_{k,n} are universal (Delta- and J-independent); state/ref ratios for any (Delta, J) follow.
Usage: python3 bg_scalar.py [K]   -> cache/bg_scalar_K{K}.pkl  (nu grid NUS)
"""
import os, sys, pickle, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from math import factorial, lgamma
import bg_ser as S, bg_geom as GG

C = os.path.join(HERE, 'cache')
NUS = np.round(np.concatenate([np.arange(0.02, 20.0, 0.04), np.arange(20.0, 30.01, 0.25)]), 4)


def cJ(J, n):
    """P_J^{(n)}(1)/n! (real J allowed via Gamma)."""
    if n == 0:
        return 1.0
    return np.prod([(J*(J + 1) - i*(i - 1)) for i in range(1, n + 1)])/(2**n*factorial(n)**2)


def build(K, nus=NUS, verbose=True):
    t = time.time()
    logs, polys, G, _ = GG.scalar_pieces(K + 3)
    A = GG.assemble(logs, polys, K)
    t0, dl = GG.spin_t(G, K)
    Bf = GG.bfac(K)
    out = dict(K=K, nus=np.asarray(nus), E0=A['E0'], E1=A['E1'], lconst=A['lconst'], M={})
    dn = [S.Poly.const(1.0)] + [S.Poly() for _ in range(K)]
    nterms = {}
    for n in range(0, K//2 + 1):
        if n:
            dn = S.lmul(dn, dl.c, K)
        for j, cj in A['corr'].items():
            P = S.lmul(S.lmul(cj, dn, K), Bf.c, K)
            for k in range(K + 1):
                if len(P[k]) == 0:
                    out['M'][(k, n, j)] = np.zeros(len(nus), complex); continue
                out['M'][(k, n, j)] = S.integrate_fast(P[k].E, P[k].C, A['E0'], A['E1'], nus)
                nterms[(k, n, j)] = len(P[k])
        if verbose:
            print(f'  bg_scalar K={K}: n={n} done ({time.time() - t:.0f}s)', flush=True)
    out['nterms'] = nterms
    return out


def load(K=6):
    fn = os.path.join(C, f'bg_scalar_K{K}.pkl')
    if not os.path.exists(fn):
        r = build(K)
        pickle.dump(r, open(fn, 'wb'))
    return pickle.load(open(fn, 'rb'))


def moments(r, Dl, J, Kt=None):
    """per-order contributions m_k(nu) = Delta^{-k} sum_n c_n(J) sum_j (i nu/2)^j M_{k,n,j}, k = 0..Kt."""
    K = r['K'] if Kt is None else Kt
    nus = r['nus']
    out = np.zeros((K + 1, len(nus)), complex)
    for (k, n, j), v in r['M'].items():
        if k <= K:
            out[k] += Dl**(-k)*cJ(J, n)*(0.5j*nus)**j*v
    return out


def ratio(r, Dl, J, Dref=2e4, Kt=None):
    """state/ref (same J) as a function of nu, and the per-order pieces of the state."""
    ms = moments(r, Dl, J, Kt); mr = moments(r, Dref, J, Kt)
    return ms.sum(0)/mr.sum(0), ms, mr


if __name__ == '__main__':
    K = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    fn = os.path.join(C, f'bg_scalar_K{K}.pkl')
    r = build(K)
    pickle.dump(r, open(fn, 'wb'))
    print('terms per (k, n, j):', r['nterms'])
