"""
Heavy vertices of CHL's <TTO> helicity structures (coordinator: step 1 of the third request).
Structures: X_h(Delta,J) (chl_xspace; T^{E,h} = C_2 X_h), exact word coefficients by t8_kappa_exact.decompose at (Delta,J) and at the
reflected labels; the right vertex of the shadow representation is (1/b^(0)) L[T X_h(Delta)] = (1/b^(0)) lambda_h L[X_h(Delta~)]
(T X_h(Delta) = lambda_h X_h(Delta~) exactly).  The stripped integrand at fixed (X, Y) = (m x, m y) (fourpt/t7_vertex_l2.helicity, same
normalisation as fourpt/vertex2 for l = 1: divided by ((-1)^J/(4 b^(0))) 2 4^{2(l+1)} |x|^{2(l+1)} |y|^{2(l+1)} e^{2P.x+2Pbar.y}) gives
   T^{hL hR}(h1,h2,h3,h4) -> V_{hL}(h1,h2) V_{hR}(h3,h4)   (heavy limit).
Points on rays of fixed beta (A_+ = s A_+0, A_- = s A_-0), s = 1, 2, 4, 8 (J non-integer allowed for three-point functions).
Output: out_t12_vertex_l2.txt, cache/t12_*.pkl
"""
import os, sys, time, pickle
import mpmath as mp, sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..')); sys.path.insert(0, os.path.join(HERE, '..', '..', 'e2e'))
from t8_kappa_exact import decompose
from embed2 import shifts, Delta as sD, J as sJ
from t1_structures_l1 import b_ratio
mp.mp.dps = 30
out = []
def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
CACHE = os.path.join(HERE, 'cache')

def coeffs(Dv, Jv, kind):
    """exact word coefficients of X_h(Dv,Jv) and of T X_h(Dv,Jv) = lambda X_h(Dv~,Jv~)"""
    fn = os.path.join(CACHE, 't12_%s_%s_%s.pkl' % (str(Dv).replace('/', '_'), str(Jv).replace('/', '_'), kind))
    if os.path.exists(fn):
        return pickle.load(open(fn, 'rb'))
    ws, s1, p1 = decompose(2, Dv, Jv, kind); s1 = s1.subs({p: 0 for p in p1})
    _, s2, p2 = decompose(2, 3 - Dv, -1 - Jv, kind); s2 = s2.subs({p: 0 for p in p2})
    cL = [s1[i] for i in range(len(ws))]
    cR = [sp.nsimplify(s1[i]/b_ratio(shifts(w)[0] - shifts(w)[1]).subs({sD: Dv, sJ: Jv})) for i, w in enumerate(ws)]   # T X_h(Delta)
    # consistency: cR = lambda s2 (mod the d = 3 relation) -- check proportionality on the largest entry pattern
    res = (ws, cL, cR, [s2[i] for i in range(len(ws))])
    pickle.dump(res, open(fn, 'wb'))
    return res

def to_mp(v):
    v = complex(sp.N(v, 40)) if not v.is_real else None
    return v

if __name__ == "__main__":
    from t7_vertex_l2 import helicity
    rays = {'b0.9': (sp.Rational(195, 4), sp.Rational(41, 2)), 'b0.4': (sp.Rational(2279, 100), sp.Rational(9, 2))}
    which = sys.argv[1] if len(sys.argv) > 1 else 'b0.4'
    scales = [int(x) for x in sys.argv[2:]] or [1, 2, 4, 8]
    Ap0, Am0 = rays[which]
    X = [mp.mpf('1.3'), mp.mpf('0.4'), mp.mpf('-0.5')]; Y = [mp.mpf('0.9'), mp.mpf('-0.3'), mp.mpf('0.2')]
    pr(f"# heavy vertices of CHL's X_h, l = 2 (fourpt/tttt/t12_vertex_l2.py), ray {which}: A_+0 = {Ap0}, A_-0 = {Am0}")
    allres = {}
    for s in scales:
        Dv = sp.Rational(3, 2) + s*Ap0; Jv = s*Am0 - sp.Rational(1, 2)
        t = time.time()
        C = {kind: coeffs(Dv, Jv, kind) for kind in ('same', 'opp')}
        tdec = time.time() - t
        R = {}
        for hL in ('same', 'opp'):
            for hR in ('same', 'opp'):
                cL = [mp.mpc(complex(sp.N(v, 40))) for v in C[hL][1]]
                cR = [mp.mpc(complex(sp.N(v, 40))) for v in C[hR][2]]
                res, leak, m = helicity(cL, cR, mp.mpf(Dv), mp.mpf(Jv), X, Y)
                R[(hL, hR)] = (res, leak)
        allres[s] = (Dv, Jv, m, R)
        pickle.dump(allres, open(os.path.join(CACHE, 't12_res_%s.pkl' % which), 'wb'))
        pr(f"  s={s}: Delta={Dv} J={Jv} m={mp.nstr(m, 8)}  (decompositions {tdec:.0f} s, integrands {time.time() - t - tdec:.0f} s)")
        for (hL, hR), (res, leak) in R.items():
            big = max(abs(v) for v in res.values())
            pr(f"    ({hL},{hR}): max|T_h| = {mp.nstr(big, 10)}, 4^-J m^8 max = {mp.nstr(big*mp.mpf(4)**(-mp.mpf(Jv))*m**8, 10)}; leakage {mp.nstr(leak, 3)}; "
               + " ".join(f"{''.join('+' if x > 0 else '-' for x in h)}:{mp.nstr(mp.re(v/big), 5)}{'' if abs(mp.im(v/big)) < 1e-20 else '%+.2gi' % float(mp.im(v/big))}" for h, v in res.items()))
    open(os.path.join(HERE, 'out_t12_vertex_l2_%s.txt' % which), 'w').write("\n".join(out) + "\n")
