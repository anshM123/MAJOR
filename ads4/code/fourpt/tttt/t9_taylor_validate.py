"""
Validation of the Taylor-mode operator application (taylor4) against the symbolic word pairs of fourpt/embed2.word4:
synthetic seed gt(P,Q) = exp(a P + b Q) (1 + c P Q^2) (exact jets), random numerical values of all 24 dot-product symbols (the operators
treat them as independent).  l = 1: all 25 cached pairs (fourpt/cache/phi4_l1.pkl);  l = 2: a few pairs built symbolically here.
Reports the max relative deviation over z-monomials.  Output: out_t9_taylor_validate.txt
"""
import os, sys, time, pickle
import numpy as np, sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..'))
from taylor4 import Engine4, seed, SN, ZN
from embed2 import word4, shifts, D1, D2, D3, D4, cached
from structures2 import words
out = []
def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
rng = np.random.default_rng(1)
sv = {n: float(rng.uniform(0.5, 2.0)) for n in SN}
zv = {n: float(rng.uniform(-1, 1)) for n in ZN}
a_, b_, c_ = 0.37, -0.61, 0.23
Psym, Qsym = sp.symbols('P Q')
gt = sp.exp(a_*Psym + b_*Qsym)*(1 + c_*Psym*Qsym**2)
lg = {n: np.log(sv[n]) for n in SN}
U = lg['s12'] + lg['s34'] - lg['s13'] - lg['s24']; V = lg['s14'] + lg['s23'] - lg['s13'] - lg['s24']
P0, Q0 = (U + V)/2, (V - U)/2
def Gjets(M):
    return [[float(sp.diff(gt, Psym, i, Qsym, j).subs({Psym: P0, Qsym: Q0})) if i + j <= M else 0.0 for j in range(M + 1)] for i in range(M + 1)]
ZS = [sp.Symbol(n, positive=False) for n in ZN]
def check(l, wL, wR, Phi, M):
    Dx = l + 1
    s1, s2 = shifts(wL); s3, s4 = shifts(wR)
    d = (Dx + s1, Dx + s2, Dx + s3, Dx + s4)
    K0 = (sv['s24']/sv['s14'])**((d[0]-d[1])/2)*(sv['s24']/sv['s23'])**((d[2]-d[3])/2)/(sv['s12']**((d[0]+d[1])/2)*sv['s34']**((d[2]+d[3])/2))
    G = Gjets(M)
    t = time.time()
    E = Engine4(M, sv)
    res = E.apply_words(seed(M, sv, d, G), wL, wR)
    tt = time.time() - t
    # symbolic
    sub = {sp.Symbol(n, positive=True): sv[n] for n in SN}
    sub.update({sp.Symbol('G_%d_%d' % (i, j)): G[i][j] for i in range(M + 1) for j in range(M + 1 - i)})
    P = sp.Poly(Phi.subs(sub), *ZS)
    ref = {m: float(c)*K0 for m, c in P.terms()}
    keys = set(ref) | set(k for k, v in res.items() if abs(v) > 0)
    scale = max(abs(v) for v in ref.values())
    err = max(abs(res.get(k, 0.0) - ref.get(k, 0.0)) for k in keys)/scale
    return err, tt, len(ref)
if __name__ == "__main__":
    pr("# Taylor-mode vs symbolic word pairs (fourpt/tttt/t9_taylor_validate.py)")
    phis = pickle.load(open(os.path.join(HERE, '..', 'cache', 'phi4_l1.pkl'), 'rb'))
    worst = 0
    for (wL, wR), Phi in phis.items():
        err, tt, nm = check(1, wL, wR, Phi, 4)
        worst = max(worst, err)
        pr(f"  l=1 {wL} x {wR}: {nm} z-monomials, max rel. deviation {err:.1e}, Taylor-mode time {tt:.2f} s")
    pr(f"  l=1 worst over 25 pairs: {worst:.1e}")
    for (wL, wR) in (('HH', 'HH'), ('HH', 'ACH'), ('BEH', 'HH'), ('AACC', 'HH'), ('HH', 'ABCE')):
        t = time.time()
        Phi, S, Gs = word4(wL, wR, 8)
        Phi = Phi.subs({D1: 3, D2: 3, D3: 3, D4: 3})
        ts = time.time() - t
        M = sum(1 for c in wL + wR if c != 'H')
        err, tt, nm = check(2, wL, wR, Phi, max(M, 1))
        pr(f"  l=2 {wL} x {wR}: {nm} z-monomials, max rel. deviation {err:.1e}; symbolic {ts:.0f} s, Taylor-mode {tt:.2f} s")
    # order-8 pairs (no symbolic reference possible): the vertex operators commute, D^(12) D^(34) = D^(34) D^(12);
    # applying the right word first exercises the truncation bookkeeping differently
    from taylor4 import tables
    for (wL, wR) in (('AACC', 'AACC'), ('ABCE', 'BBEE'), ('ACH', 'AAEE')):
        Dx = 3; s1, s2 = shifts(wL); s3, s4 = shifts(wR)
        d = (Dx + s1, Dx + s2, Dx + s3, Dx + s4); M = sum(1 for c in wL + wR if c != 'H')
        G = Gjets(M); E = Engine4(M, sv)
        A0 = seed(M, sv, d, G)
        r1 = E.apply_words(A0, wL, wR)
        A = A0
        for kind in reversed(wR): A = E.op(A, kind, 3, 4)
        for kind in reversed(wL): A = E.op(A, kind, 1, 2)
        r2 = {k: v[0] for k, v in A.items() if v.shape[0]}
        keys = set(r1) | set(r2); scale = max(abs(v) for v in r1.values())
        err = max(abs(r1.get(k, 0) - r2.get(k, 0)) for k in keys)/scale
        pr(f"  l=2 {wL} x {wR} (order {M}): {len(r1)} z-monomials; (12)-first vs (34)-first max rel. deviation {err:.1e}")
    open(os.path.join(HERE, 'out_t9_taylor_validate.txt'), 'w').write("\n".join(out) + "\n")
