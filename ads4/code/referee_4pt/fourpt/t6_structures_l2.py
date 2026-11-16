"""
<T T O> (two stress tensors, Delta_1 = Delta_2 = 3), parity even, differential basis of 14 words (structures2.words(2)).
Numerics: (a) at 60 digits (mpmath); (b), (c) in double precision (numpy; 6th-order differences, h = 1e-3, error ~1e-11), words
lambdified from the exact expressions of build_l2.py:
 (a) the d = 3 linear relation among the 14 words (singular values of the word matrix on physical 3d kinematics);
 (b) conservation at x_1 and x_2 (divergence of the traceless spin-2 tensor by 6th-order finite differences, h = 1e-8):
     null space dimension (incl. the relation) -> number of conserved structures;
 (c) Lorentzian shadow map  c_w -> c_w b^(0)/b_w  (b_w = b^{(Delta_34 of the seed)}, Delta_34 in {0, +-2, +-4}): the image of the
     conserved space at (Delta,J) must lie in the conserved space at (3-Delta,-1-J) [test: smallest singular value of the
     enlarged matrix]; control: b^(0) for all words;
Output: out_t6_structures_l2.txt
"""
import os, sys, pickle, functools
import mpmath as mp, sympy as sp
from embed2 import Delta, J, CACHE, shifts
from structures2 import words, logMref
from t1_structures_l1 import b_ratio
mp.mp.dps = 60
out = []
def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)

REL, S, V5 = pickle.load(open(os.path.join(CACHE, 'rel_l2.pkl'), 'rb'))
WS = words(2)
NAMES = ['s12', 's15', 's25', 'z12', 'z15', 'z21', 'z25', 'z51', 'z52', 'z1z2', 'z1z5', 'z2z5']
ARGS = [S[n] for n in NAMES] + [Delta, J]
FQ = {w: sp.lambdify(ARGS, REL[w], modules='mpmath') for w in WS}
import numpy as np
FQn = {w: sp.lambdify(ARGS, REL[w], modules='numpy') for w in WS}
LMn = sp.lambdify(ARGS, logMref(2, S, V5).subs(sp.log(V5), sp.log(sp.Abs(V5))), modules='numpy')
LM = sp.lambdify(ARGS, logMref(2, S, V5).subs(sp.log(V5), sp.log(sp.Abs(V5))), modules='mpmath')
BR = {w: sp.lambdify([Delta, J], b_ratio(shifts(w)[0] - shifts(w)[1]), modules='mpmath') for w in WS}

def md(u, v): return -u[0]*v[0] + u[1]*v[1] + u[2]*v[2]
def sub(u, v): return [u[i] - v[i] for i in range(3)]

def inv_args(x1, x2, x5, n, z1, z2, D, Jv):
    return [-md(sub(x1, x2), sub(x1, x2)), -md(sub(x1, x5), sub(x1, x5)), -md(sub(x2, x5), sub(x2, x5)),
            md(z1, sub(x2, x1)), md(z1, sub(x5, x1)), md(z2, sub(x1, x2)), md(z2, sub(x5, x2)), md(n, sub(x1, x5)), md(n, sub(x2, x5)),
            md(z1, z2), md(z1, n), md(z2, n), D, Jv]

def wordvals(x1, x2, x5, n, z1, z2, D, Jv):
    a = inv_args(x1, x2, x5, n, z1, z2, D, Jv)
    return mp.exp(LM(*a)), [FQ[w](*a) for w in WS]

NULL5 = [[1, mp.cos(t), mp.sin(t)] for t in (2*mp.pi*k/5 for k in range(5))]
def tensor_at(pt, x1, x2, x5, n, zother, D, Jv):
    """traceless symmetric tensor at operator pt (1 or 2) for each word: returns list over words of 3x3 matrices (covariant? -> we
    reconstruct T^{mn} contracted as z_m z_n with the Lorentzian dot of contravariant z: W(z) = T_{mn} z^m z^n)."""
    vals = []
    for z in NULL5:
        if pt == 1:
            M, q = wordvals(x1, x2, x5, n, z, zother, D, Jv)
        else:
            M, q = wordvals(x1, x2, x5, n, zother, z, D, Jv)
        vals.append([M*v for v in q])
    # unknowns (T01, T02, T11, T12, T22), T00 = T11 + T22 (traceless w.r.t. eta: -T00 + T11 + T22 = 0)
    A = mp.matrix(5, 5)
    for k, z in enumerate(NULL5):
        A[k, 0] = 2*z[0]*z[1]; A[k, 1] = 2*z[0]*z[2]; A[k, 2] = z[1]**2 + z[0]**2; A[k, 3] = 2*z[1]*z[2]; A[k, 4] = z[2]**2 + z[0]**2
    Ainv = A**-1
    res = []
    for wi in range(len(WS)):
        sol = Ainv*mp.matrix([vals[k][wi] for k in range(5)])
        T = [[None]*3 for _ in range(3)]
        T[0][1] = T[1][0] = sol[0]; T[0][2] = T[2][0] = sol[1]; T[1][1] = sol[2]; T[1][2] = T[2][1] = sol[3]; T[2][2] = sol[4]; T[0][0] = sol[2] + sol[4]
        res.append(T)
    return res

ETA = [-1, 1, 1]
def divergence_rows(D, Jv, x1, x2, x5, n, others, h=mp.mpf('1e-8')):
    """rows of the linear conservation system (one row per (point, other polarisation, free index)), columns = words"""
    st = [-3, -2, -1, 1, 2, 3]; wt = [-mp.mpf(1)/60, mp.mpf(3)/20, -mp.mpf(3)/4, mp.mpf(3)/4, -mp.mpf(3)/20, mp.mpf(1)/60]
    rows = []
    for pt in (1, 2):
        for zo in others:
            dT = [[[mp.mpf(0)]*3 for _ in range(3)] for _ in WS]      # div_nu per word: sum_mu eta^mu d_mu T_{mu nu}
            for mu in range(3):
                for s_, w_ in zip(st, wt):
                    e = [0, 0, 0]; e[mu] = s_*h
                    xa = [x1[i] + (e[i] if pt == 1 else 0) for i in range(3)]; xb = [x2[i] + (e[i] if pt == 2 else 0) for i in range(3)]
                    Ts = tensor_at(pt, xa, xb, x5, n, zo, D, Jv)
                    for wi in range(len(WS)):
                        for nu in range(3):
                            dT[wi][mu][nu] += w_*Ts[wi][mu][nu]/h
            for nu in range(3):
                rows.append([sum(ETA[mu]*dT[wi][mu][nu] for mu in range(3)) for wi in range(len(WS))])
    return mp.matrix(rows)

def svals(M):
    return sorted(mp.svd_r(M, compute_uv=False), reverse=True) if M.rows >= M.cols else sorted(mp.svd_r(M.T, compute_uv=False), reverse=True)

def nullspace(M, tol):
    U, Sv, V = mp.svd_r(M)
    n = M.cols
    svl = [Sv[i] for i in range(min(M.rows, n))] + [mp.mpf(0)]*(n - min(M.rows, n))
    return [[V[i, j] for j in range(n)] for i in range(n) if svl[i] < tol], svl

X1 = [mp.mpf(1)/3, mp.mpf(2)/7, -mp.mpf(3)/5]; X2 = [-mp.mpf(7)/3, mp.mpf(1)/2, -mp.mpf(1)/5]; X5 = [mp.mpf(9)/4, -mp.mpf(2)/3, mp.mpf(3)/7]; N5 = [1, mp.mpf(3)/5, mp.mpf(4)/5]
OTHERS = [[1, mp.mpf(-5)/13, mp.mpf(12)/13], [1, mp.mpf(8)/17, mp.mpf(-15)/17], [1, mp.mpf(-7)/25, mp.mpf(-24)/25]]

def wordvals_n(x1, x2, x5, n, z1, z2, D, Jv):
    a = inv_args(x1, x2, x5, n, z1, z2, D, Jv)
    return np.exp(LMn(*a)), [FQn[w](*a) for w in WS]

def tensor_at_n(pt, x1, x2, x5, n, zother, D, Jv):
    vals = []
    for z in NULL5f:
        M, q = wordvals_n(x1, x2, x5, n, z, zother, D, Jv) if pt == 1 else wordvals_n(x1, x2, x5, n, zother, z, D, Jv)
        vals.append([M*v for v in q])
    A = np.zeros((5, 5))
    for k, z in enumerate(NULL5f):
        A[k] = [2*z[0]*z[1], 2*z[0]*z[2], z[1]**2 + z[0]**2, 2*z[1]*z[2], z[2]**2 + z[0]**2]
    sol = np.linalg.solve(A, np.array(vals))           # [5, nwords]
    T = np.zeros((len(WS), 3, 3))
    T[:, 0, 1] = T[:, 1, 0] = sol[0]; T[:, 0, 2] = T[:, 2, 0] = sol[1]; T[:, 1, 1] = sol[2]; T[:, 1, 2] = T[:, 2, 1] = sol[3]
    T[:, 2, 2] = sol[4]; T[:, 0, 0] = sol[2] + sol[4]
    return T

NULL5f = [[1.0, np.cos(t), np.sin(t)] for t in (2*np.pi*k/5 for k in range(5))]
def divergence_rows_n(D, Jv, x1, x2, x5, n, others, h=1e-3):
    st = [-3, -2, -1, 1, 2, 3]; wt = [-1/60, 3/20, -3/4, 3/4, -3/20, 1/60]
    rows = []
    for pt in (1, 2):
        for zo in others:
            dT = np.zeros((len(WS), 3, 3))
            for mu in range(3):
                for s_, w_ in zip(st, wt):
                    e = np.zeros(3); e[mu] = s_*h
                    xa = np.array(x1) + (e if pt == 1 else 0); xb = np.array(x2) + (e if pt == 2 else 0)
                    dT[:, mu, :] += w_*tensor_at_n(pt, list(xa), list(xb), x5, n, zo, D, Jv)[:, mu, :]/h
            div = -dT[:, 0, :] + dT[:, 1, :] + dT[:, 2, :]      # [word, nu]
            for nu in range(3):
                rows.append(div[:, nu])
    return np.array(rows)

def cons_space(D, Jv):
    D = float(D)
    M = divergence_rows_n(D, Jv, [float(v) for v in X1], [float(v) for v in X2], [float(v) for v in X5], [float(v) for v in N5],
                          [[float(v) for v in z] for z in OTHERS])
    sc = np.abs(M).max(axis=0)           # column scaling (word coefficients have very different sizes)
    U, sv, Vt = np.linalg.svd(M/sc)
    ns = [Vt[i] for i in range(len(sv)) if sv[i] < 1e-7]
    return ns, list(sv), sc

def cons_space_old(D, Jv):
    M = divergence_rows(D, Jv, X1, X2, X5, N5, OTHERS)
    scale = max(abs(M[i, j]) for i in range(M.rows) for j in range(M.cols))
    ns, sv = nullspace(M/scale, mp.mpf('1e-25'))
    return ns, sv

if __name__ == "__main__":
    pr("# <T T O> parity-even structures, differential basis, 60-digit numerics (fourpt/t6_structures_l2.py)")
    # (a) relation
    import random
    random.seed(5)
    rows = []
    NULLS = [[1, mp.cos(t), mp.sin(t)] for t in (0.3, 1.1, 2.0, 2.9, 4.0, 5.2)]
    for k in range(20):
        x1 = [mp.mpf(random.uniform(-2, 2)) for _ in range(3)]; x2 = [mp.mpf(random.uniform(-2, 2)) for _ in range(3)]
        x5 = [mp.mpf(random.uniform(-2, 2)) for _ in range(3)]
        n, z1, z2 = random.sample(NULLS, 3)
        Mref, q = wordvals(x1, x2, x5, n, z1, z2, mp.mpf('9.7'), 6)
        rows.append(q)
    Mw = mp.matrix(rows)
    for i in range(Mw.rows):
        sc = max(abs(Mw[i, j]) for j in range(Mw.cols))
        for j in range(Mw.cols): Mw[i, j] /= sc
    rel, sv = nullspace(Mw, mp.mpf('1e-30'))
    pr("(a) singular values of the 20 x 14 word matrix on physical 3d kinematics:", [mp.nstr(x, 3) for x in sv])
    rel = rel[0]; piv = max(range(len(rel)), key=lambda i: abs(rel[i])); rel = [x/rel[piv] for x in rel]
    pr("    relation:", {w: mp.nstr(x, 12) for w, x in zip(WS, rel) if abs(x) > 1e-30})
    for (D, Jv) in ((mp.mpf('9.7'), 6), (mp.mpf('13.25'), 5), (mp.mpf('11.5'), 8)):
        Dt, Jt = 3 - D, -1 - Jv
        ns1, sv1, sc1 = cons_space(D, Jv); ns2, sv2, sc2 = cons_space(Dt, Jt)
        ns1 = [v/sc1 for v in ns1]; ns2 = [v/sc2 for v in ns2]           # undo the column scaling
        pr(f"(b) (Delta,J) = ({mp.nstr(D, 5)},{Jv}): conservation singular values (smallest 5): {[f'{x:.1e}' for x in sorted(sv1)[:5]]} -> null space {len(ns1)};"
           f"  at (3-Delta,-1-J): {[f'{x:.1e}' for x in sorted(sv2)[:5]]} -> {len(ns2)}")
        # (c) image of conserved space under T within span(ns2) (which contains the relation)
        def resid(vecs, target):
            sc = np.abs(np.array(target)).max(axis=0)
            A = (np.array(target)/sc).T
            Q, R = np.linalg.qr(A)
            out_ = []
            for v in vecs:
                v = np.array(v, dtype=float)/sc
                pv = Q @ (Q.T @ v)
                out_.append(np.linalg.norm(v - pv)/np.linalg.norm(v))
            return out_
        Tv = [[v[i]/float(BR[w](D, Jv)) for i, w in enumerate(WS)] for v in ns1]
        pr(f"(c)   |T v - proj_(conserved at reflected labels)(T v)|/|T v| for the {len(ns1)} null vectors: {[f'{x:.1e}' for x in resid(Tv, ns2)]};"
           f"  control b^(0) for all words: {[f'{x:.1e}' for x in resid(ns1, ns2)]}")
    open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out_t6_structures_l2.txt'), 'w').write("\n".join(out) + "\n")
