"""md_k3.py -- the k = 3 (Q_k) part of the ladder, own implementation: q_n from the d = 1 CDH recurrence of sec_fullpos.tex
(q_0 = 1, (n+4)^2 q_{n+1} = -(s + 4 - 2(n+2)^2) q_n - n^2 q_{n-1}, s = nu^2 + 1/4), Q_k = int_0^1 psi3(p) q_k(nu^2 = Lam^2 p^2) dp
with psi3 from prl/data/nf02_psi.json (exact rationals), and hhat^(3)_b = sum_k D_k^2 U_{b,k} sigma_b Lam^2 Q_k.
Also the per-nu k = 3 ladder polynomial  S3_b(s) = sigma_b sum_k D_k^2 U_{b,k} q_k(s)  for the direct k = 3 test."""
import os, sys, json
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import sympy as sp
import md_ladder as ML

s = ML.s
DATA = os.path.join(HERE, '..', '..', 'prl', 'data')


def qseq(K):
    q = [sp.Integer(1)]
    qm = sp.Integer(0)
    for n in range(K):
        nxt = sp.expand((-(s + 4 - 2*(n + 2)**2)*q[n] - n*n*(q[n - 1] if n >= 1 else 0))/sp.Integer((n + 4)**2))
        q.append(nxt)
    return q


def psi3():
    d = json.load(open(os.path.join(DATA, 'nf02_psi.json')))
    p3 = d['psi3_exact']; p3 = eval(p3) if isinstance(p3, str) else p3
    co = [sp.Rational(c) for c in p3]
    k0 = sp.Rational(d['k0'])
    return co, k0


def Qk(K, lam):
    co, k0 = psi3()
    q = qseq(K)
    nu2 = sp.Symbol('nu2')
    out = []
    for k in range(K + 1):
        P = sp.Poly(sp.expand(q[k].subs(s, nu2 + sp.Rational(1, 4))), nu2)
        v = 0
        for (i,), c in P.terms():
            v += c*sp.Integer(lam)**i*sum(cn/sp.Integer(2*i + n + 1) for n, cn in enumerate(co))/k0
        out.append(sp.nsimplify(v))
    return out


def S3(b, Delta, J, q):
    Delta = sp.nsimplify(Delta); n = J - b; sig = 1 if b == 0 else -1
    tot = 0
    for k in range(n + 1):
        tot += sp.binomial(k + 3, 3)**2*ML.U(b, k, Delta, J)*q[k]
    return sp.expand(sig*tot)


def hhat3(b, Delta, J, Q, lam):
    Delta = sp.nsimplify(Delta); n = J - b; sig = 1 if b == 0 else -1
    return sum(sp.binomial(k + 3, 3)**2*ML.U(b, k, Delta, J)*sig*lam*Q[k] for k in range(n + 1))


if __name__ == '__main__':
    lam = 998*999
    Q = Qk(22, lam)
    sys.path.insert(0, os.path.join(HERE, '..', 'fullpos'))
    import fp_poly as FP
    Qe = FP.Q_exact(4, lam)
    print('Q_k own vs fp_poly.Q_exact:', [float(Q[k])/float(Qe[k].p)*float(Qe[k].q) for k in range(5)])
