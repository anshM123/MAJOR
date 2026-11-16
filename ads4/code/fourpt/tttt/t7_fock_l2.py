"""
Direct two-particle Fock computation of the MFT data of [T T]_{N,J} (identical conserved stress tensors, Delta_T = 3) in the helicity basis:
CHL's parity-even T^{E,h} constructed in POSITION space (chl_xspace: term-wise inverse FT, exact Pochhammer ratios, common scalar C(Delta,J)
dropped, double-twist limit Delta -> 6 + N taken as the leading term in eps (C eps^k S_h)), mapped to our frame, checked polynomial and
conserved (fock2 range check), Fock Gram matrix exact.  Outputs per (N, J): orders k_h, Gram matrix (off-diagonal must vanish),
p[S_h] = 2^J/G_hh, and the test of CHL's MFT helicity factors:  p[T_opp]/p[T_same] = A_opp/A_same,
A_opp = (-J)_4/(J+1)_4, A_same = (Delta-1)_4/(2-Delta)_4  (leading terms in eps; Delta = 6 + N - 2 eps).
Also l = 1 (currents) as a control, where the Fock data are known in closed form (t1_mft_l1).
Output: out_t7_fock_l2.txt, cache/fock_l2.pkl
"""
import os, sys, time, pickle
import sympy as sp
from chl_xspace import xspace_structure_limit, to_our_frame, EPS, _lead
from ft_helicity import AX, BX, WX, XX, AB, AW, BW
from infframe import X, Z1, Z2, dot, CACHE
from mft2 import ZETA, gram
out = []
def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)

def to_components(F):
    for term in sp.Add.make_args(F):
        for b, e in term.as_powers_dict().items():
            if b in (WX, XX) and (e < 0 or not sp.sympify(e).is_integer):
                raise AssertionError("non-polynomial term %s" % term)
    zs = ZETA
    csub = {AX: dot(Z1, X), BX: dot(Z2, X), WX: dot(zs, X), XX: dot(X, X), AB: dot(Z1, Z2), AW: dot(Z1, zs), BW: dot(Z2, zs)}
    return sp.expand(F.xreplace(csub))

def A_ratio(l, Dv, Jv):
    D = Dv - 2*EPS
    Ao = sp.rf(-Jv, 2*l)/sp.rf(Jv + 1, 2*l)
    As = sp.expand_func(sp.rf(D - 1, 2*l))/sp.expand_func(sp.rf(2 - D, 2*l))
    return _lead(Ao/As)

def run(l, N, Jv):
    Dv = 2*(l + 1) + N
    t = time.time()
    S = {}; k = {}
    for kind in ('same', 'opp'):
        if kind == 'opp' and Jv < 2*l: continue
        kk, F = xspace_structure_limit(l, Dv, Jv, kind)
        if F == 0: continue
        try:
            S[kind] = to_components(to_our_frame(F)); k[kind] = kk
        except AssertionError:
            pr(f"  l={l} N={N} J={Jv}: CHL {kind} structure not polynomial at this double-twist point (J < 2l) -- skipped")
    if not S: return None
    names = list(S)
    G = gram(l, [S[n] for n in names])
    p = {n: 2**Jv/G[i, i] for i, n in enumerate(names)}
    res = dict(names=names, k=k, G=G, p=p)
    msg = f"  l={l} N={N} J={Jv} n={sp.Rational(2*l + N - Jv, 2)}: orders {k}; "
    if len(names) == 2:
        det = G.det()
        res['det'] = det
        msg += f"G_offdiag = {G[0, 1]}, det = {'0 (degenerate)' if det == 0 else 'nonzero'}; "
        kr, cr = A_ratio(l, Dv, Jv)
        # p[T_opp]/p[T_same] = p[S_opp]/p[S_same] * eps^{2 k_same - 2 k_opp}
        ko = 2*k['same'] - 2*k['opp']; co = sp.nsimplify(p['opp']/p['same'])
        res['test'] = ((ko, co), (kr, cr))
        msg += f"p[T_opp]/p[T_same] = {co} eps^{ko}  vs  A_opp/A_same = {cr} eps^{kr}  -> {'MATCH' if (ko, sp.simplify(co - cr)) == (kr, 0) else 'MISMATCH'}; "
    msg += f"p[S_h] = {p}  ({time.time() - t:.0f} s)"
    pr(msg)
    return res

if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else 'l2'
    allres = {}
    if which == 'l1':
        pr("# control l = 1: CHL x-space T^{E,h} in the two-particle Fock space")
        for N in range(0, 5):
            for Jv in range(N + 2, -1, -2):
                allres[(1, N, Jv)] = run(1, N, Jv)
        pickle.dump(allres, open(os.path.join(CACHE, 'fock_l1.pkl'), 'wb'))
    else:
        pr("# [TT]_{N,J}: Fock data in CHL's x-space helicity basis (fourpt/tttt/t7_fock_l2.py); Delta = 6 + N, n = (4 + N - J)/2")
        Nmax = int(sys.argv[2]) if len(sys.argv) > 2 else 4
        for N in range(0, Nmax + 1):
            for Jv in range(N + 4, -1, -2):
                allres[(2, N, Jv)] = run(2, N, Jv)
                pickle.dump(allres, open(os.path.join(CACHE, 'fock_%s.pkl' % which), 'wb'))
    open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out_t7_fock_%s.txt' % which), 'w').write("\n".join(out) + "\n")
