"""
Task A: exact MFT OPE coefficients of [O phi]_{N,J} (O = J: l = 1; O = T: l = 2) from the two-particle Fock space
(fock.py), and comparison with the closed forms

  parity even, structure  L := sum_n L_n S_n  (L_n of structures.py; = 2^J x the differential-basis combination; S_cons = L/L_0):
      p_L(Delta, J) = 2^J P_FK(1, Dphi; Delta, J) / ( c_l (J-l+1)_{2l} (Delta-l-1)_{2l} ),   c_1 = 4, c_2 = 24
  parity odd, l = 1, structure at infinity  eps(z1, zeta, x) (zeta.x)^{J-1} (x^2)^{p/2}:
      p_O(Delta, J) = 2^J P_FK(2, Dphi; Delta, J) J (Delta-1) / ( 2 (J+1) (Delta-2) )
  parity odd, l = 2, structure at infinity (Delta-Dphi) eps (z1.zeta)(zeta.x)^{J-2}(x^2)^{(p+2)/2} - (Delta-Dphi-J-2) eps (z1.x)(zeta.x)^{J-1}(x^2)^{p/2}:
      p_O(Delta, J) = 2^J P_FK(2, Dphi; Delta, J) J (J-1) Delta (Delta-1) / ( 12 (J+1)(J+2)(Delta-2)(Delta-3) )

P_FK(D1, D2; Delta, J) = Fitzpatrick-Kaplan coefficient (lambda^2 = 2^J P_FK for <phi1 phi2 X> in the usual structure),
continued in n = (Delta - D1 - D2 - J)/2.
Usage: python3 mft_exact.py  [writes out_mft_exact.txt via run_all.sh]
"""
import sys, os, pickle, time
import sympy as sp
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here); sys.path.insert(0, os.path.join(here, '..', 'e2e'))
from fock import even_monomials, odd_monomials, conserved_combinations, fock_norm
import structures as S

CACHE = os.path.join(here, 'cache'); os.makedirs(CACHE, exist_ok=True)
D_, J_, Dp_ = sp.symbols('Delta J Dphi')

def Lpoly(l):
    c, se = S.get(l)
    return [sp.factor(sp.simplify(sum(c[w]*se[w][n] for w in c))) for n in range(l+1)]

def PFK(D1, D2, Delta, J):
    n = (Delta - D1 - D2 - J)/2
    h = sp.Rational(3, 2); g = sp.gamma
    rf = lambda a, k: g(a + k)/g(a)
    return sp.simplify(rf(D1-h+1, n)*rf(D2-h+1, n)*rf(D1, J+n)*rf(D2, J+n) /
                       (g(n+1)*sp.factorial(J)*rf(J+h, n)*rf(D1+D2+n-2*h+1, n)*rf(D1+D2+2*n+J-1, J)*rf(D1+D2+n+J-h, n)))

C_L = {1: 4, 2: 24}

def closed_even(l, Delta, J, Dphi):
    return 2**J*PFK(1, Dphi, Delta, J)/(C_L[l]*sp.rf(J-l+1, 2*l)*sp.rf(Delta-l-1, 2*l))

def closed_odd_l1(Delta, J, Dphi):
    return 2**J*PFK(2, Dphi, Delta, J)*J*(Delta-1)/(2*(J+1)*(Delta-2))

def closed_odd_l2(Delta, J, Dphi):
    """structure (Delta-Dphi) eps (z1.zeta)(zeta.x)^{J-2}|x|^{p+2} - (Delta-Dphi-J-2) eps (z1.x)(zeta.x)^{J-1}|x|^p"""
    return 2**J*PFK(2, Dphi, Delta, J)*J*(J-1)*Delta*(Delta-1)/(12*(J+1)*(J+2)*(Delta-2)*(Delta-3))

def compute(l, Dphi, Nmax):
    fn = os.path.join(CACHE, 'fock_l%d_D%s_N%d.pkl' % (l, str(Dphi).replace('/', '_'), Nmax))
    if os.path.exists(fn): return pickle.load(open(fn, 'rb'))
    Lp = Lpoly(l)
    res = {}
    for N in range(Nmax + 1):
        for J in range(l, N + l + 1):
            Delta = l + 1 + Dphi + N
            ev = even_monomials(l, N, J); od = odd_monomials(l, N, J)
            if ev:
                Ln = [L.subs({D_: Delta, J_: J, Dp_: Dphi}) for L in Lp]
                q = sp.expand(sum(Ln[n]*(-1)**(l-n)*sp.Rational(-1, 2)**n*ev[n] for n in ev))
                ce = conserved_combinations(ev, l)
                assert len(ce) == 1
                ks = sorted(ev)
                assert sp.Matrix([[ce[0][k] for k in ks], [Ln[n]*(-1)**(l-n)*sp.Rational(-1, 2)**n for n in ks]]).rank() == 1
                res[('E', N, J)] = sp.Integer(2)**J/fock_norm(q, l, N, Dphi)
            co = conserved_combinations(od, l) if od else []
            if co:
                assert len(co) == 1
                v = co[0]; q = sp.expand(sum(v[k]*od[k] for k in od))
                res[('O', N, J)] = (sp.Integer(2)**J/fock_norm(q, l, N, Dphi), v)
            assert not (ev and co)
    pickle.dump(res, open(fn, 'wb'))
    return res

RUNS = [(1, sp.Rational(7, 4), 8), (1, sp.Rational(5, 2), 5), (1, sp.Rational(13, 10), 5),
        (2, sp.Rational(7, 4), 6), (2, sp.Rational(5, 2), 4), (2, sp.Rational(13, 10), 4)]

if __name__ == "__main__":
    which = sys.argv[1:] and [RUNS[int(i)] for i in sys.argv[1:]] or RUNS
    for (l, Dphi, Nmax) in which:
        t = time.time()
        res = compute(l, Dphi, Nmax)
        print("=== l = %d, Dphi = %s, N <= %d   (%d primaries; %.0f s)" % (l, Dphi, Nmax, len(res), time.time() - t))
        nE = nO = 0; okE = okO = True
        for key in sorted(res, key=lambda k: (k[1], k[2])):
            kind, N, J = key
            Delta = l + 1 + Dphi + N
            if kind == 'E':
                val = res[key]; pred = closed_even(l, Delta, J, Dphi)
                r = sp.nsimplify(sp.simplify(val/pred)); nE += 1; okE &= (r == 1)
                print("  E  N=%d J=%d  tau-tau0=%d  p_L = %s   p_L/closed = %s" % (N, J, N - J + l, val, r))
            else:
                val, v = res[key]
                if l == 1:
                    pred = closed_odd_l1(Delta, J, Dphi); r = sp.nsimplify(sp.simplify(val/pred)); nO += 1; okO &= (r == 1)
                    print("  O  N=%d J=%d  tau-tau0=%d  p_O = %s   p_O/closed = %s" % (N, J, N - J + l, val, r))
                else:
                    v0 = v.get(0, 0); v0pred = -(Delta - Dphi - J - 2)/(Delta - Dphi) if 0 in v else 0
                    pred = closed_odd_l2(Delta, J, Dphi)*(Delta - Dphi)**2
                    r = sp.nsimplify(sp.simplify(val/pred)); nO += 1; okO &= (r == 1) and (sp.nsimplify(v0 - v0pred) == 0)
                    print("  O  N=%d J=%d  tau-tau0=%d  p_O = %s  (n=1 coefficient 1, n=0 coefficient %s; predicted %s)  p_O/closed = %s" % (
                        N, J, N - J + l, val, v0, v0pred, r))
        print("  even: %d primaries, all exact matches: %s" % (nE, okE))
        print("  odd:  %d primaries, all exact matches (incl. structure coefficients for l = 2): %s" % (nO, okO))
