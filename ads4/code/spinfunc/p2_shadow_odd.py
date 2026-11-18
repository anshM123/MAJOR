"""
p2 (parity-odd sector, PARITY_ODD_ROUTEA.md sect. 2): shadow coefficient of the parity-odd conserved <T T O_{Delta,J}> (p1), from the
Laplace expansion of the shadow integral around x5 = e to relative O(1/Delta^2) (l3/l5_shadow_matrix route), jointly with the parity-even
structures, so that the even-odd block is checked to vanish:

   int d^3x5 <O(x6,z6) O(x5,.)> . <T3 T4 O~(x5,.)>_c = sum_b beta_{cb} <T3 T4 O(x6,z6)>_b,   c, b over [even (l5_spin.conserved), odd (p1)].

Odd structures carry the extra factor (X34 X35 X45)^{-1/2} (LHS, m_ij = -x_ij^2 as in l5_shadow_matrix._setup) and (X34 X36 X46)^{-1/2}
(RHS, with the same sign/branch conventions as l5_shadow_matrix.rhs: m34 = -x34^2, m36 = +x36^2, m46 = +x46^2); both are folded into the
numerator TPS / numbers so that one fit covers all structures.
Validation (as l5): fit residual ~ Delta^-3 at order 2; beta independent of the crystal point (xi = nu/Delta) and of the polarisations;
beta_{EO} = beta_{OE} = 0 (parity); beta_OO real.
Usage: python3 p2_shadow_odd.py -> cache/p2_shadow_odd.txt
"""
import os, sys, functools, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import l1_tps as T
import l4_tto_even as L4
from l1_crystal import dot, sub
import l5_spin as S
import l5_shadow_matrix as SM
import p1_tto_odd as P


def _odd_num_tps(xa, za, xb, zb, x5, z5, J, mons, sol):
    """odd numerator (incl. eps.X_ij factor, no denominators) at TPS x5, z5 (3-lists of TPS); xa, xb TPS consts; za, zb numeric."""
    Xa = P.emb_X(L4, xa); Xb = P.emb_X(L4, xb); X5 = P.emb_X(L4, x5)
    Za = P.emb_Z(L4, xa, S.tv(za)); Zb = P.emb_Z(L4, xb, S.tv(zb)); Z5 = P.emb_Z(L4, x5, z5)
    v = P.structures_odd_emb(L4, Xa, Za, Xb, Zb, X5, Z5, 0.0, J, list(mons), numerators=True)
    return sum(sol[i, 0]*v[i] for i in range(len(v)))


def lhs_series(Delta, J, xi, trips, with_even=True):
    x3, x4 = SM.crystal_x3x4(xi)
    zs = np.zeros(T.NV, complex)
    zv = [T.var(i, 0.0) for i in range(T.NV)]
    lamF, G, *_ = SM._setup(zv, Delta, J, x3, x4)
    _, H = T.grad_hess(lamF)
    L = np.linalg.inv(T.csym_cholesky(-H)).T
    for _ in range(2):
        lamF2, G2, *_ = SM._setup(T.linsub(zs, L), Delta, J, x3, x4)
        _, H2 = T.grad_hess(lamF2)
        L = L @ np.linalg.inv(T.csym_cholesky(-H2)).T
    lamF2, G2, x5, x56, X3, X4 = SM._setup(T.linsub(zs, L), Delta, J, x3, x4)
    msq = lambda p, q: -dot(sub(p, q), sub(p, q))
    extra = T.powr(T.mul(T.mul(msq(X3, X4), msq(X3, x5)), msq(X4, x5)), -0.5)       # (m34 m35 m45)^{-1/2}
    monsO, solO, _ = P.conserved_odd(3.0 - Delta, J)
    if with_even and J % 2 == 0:
        monsE, solE, _ = S.conserved(3.0 - Delta, J)
    out = []
    for (z3, z4, z6) in trips:
        a = SM._avec(x56, z6)
        prefs = []
        if with_even and J % 2 == 0:
            pe = SM._tto_at_tps_pol(X3, z3, X4, z4, x5, a, J, monsE, solE)
            prefs += [pe[c] for c in range(pe.shape[0])]
        if solO is not None:
            prefs.append(T.mul(_odd_num_tps(X3, z3, X4, z4, x5, a, J, monsO, solO), extra))
        row = []
        for pf in prefs:
            S0, S2, S4 = T.expand_pref(lamF2, G2, pf)
            row.append((S0, S0 + S2, S0 + S2 + S4))
        out.append(row)
    logI0 = lamF[0] + G[0] + 1.5*np.log(2*np.pi) + np.log(np.linalg.det(L[:3, :3]))
    return logI0, np.array(out)


def rhs(Delta, J, xi, trips, with_even=True):
    x3, x4 = SM.crystal_x3x4(xi)
    x6 = np.array([-1.0, -1.0, 0.0])
    d = lambda p, q: -0.5*(p[0]*q[1] + p[1]*q[0]) + p[2]*q[2]
    m36 = d(x3 - x6, x3 - x6); m46 = d(x4 - x6, x4 - x6); m34 = -d(x3 - x4, x3 - x4)
    den = -((10 - Delta - J)/2)*np.log(m34 + 0j) - ((Delta + J)/2)*(np.log(m36 + 0j) + np.log(m46 + 0j))
    extra = np.exp(-0.5*(np.log(m34 + 0j) + np.log(m36 + 0j) + np.log(m46 + 0j)))
    monsO, solO, _ = P.conserved_odd(Delta, J)
    if with_even and J % 2 == 0:
        monsE, solE, _ = S.conserved(Delta, J)
    r = []
    for (z3, z4, z6) in trips:
        row = []
        if with_even and J % 2 == 0:
            row += list(S.tto_samples(S.tv(x3), S.tv(z3), S.tv(x4), z4[None, :], S.tv(x6), z6[None, :], J, monsE, solE)[:, 0, 0, 0])
        if solO is not None:
            v = np.array(P.structures_odd(L4, x3, z3, x4, z4, x6, z6, 0.0, J, list(monsO), numerators=True))
            row.append(complex(v @ solO[:, 0])*extra)
        r.append(row)
    return den, np.array(r)


def shadow_all(Delta, J, xi=0.5, ntrip=8, seed=7, with_even=True):
    rng = np.random.default_rng(seed)
    trips = [(SM._rnull(rng), SM._rnull(rng), SM._rnull(rng)) for _ in range(ntrip//2)]
    R = np.array([1.0, 1.0, -1.0])                    # transverse reflection y -> -y: maps the configuration to itself; with the
    trips = trips + [(R*a, R*b, R*c) for (a, b, c) in trips]   # reflected triples the even-odd blocks of the fit vanish exactly
    logI0, Sv = lhs_series(Delta, J, xi, trips, with_even)
    den, r = rhs(Delta, J, xi, trips, with_even)
    res = {}
    for o in range(3):
        b, *_ = np.linalg.lstsq(r, Sv[:, :, o], rcond=None)
        resid = np.linalg.norm(r @ b - Sv[:, :, o])/np.linalg.norm(Sv[:, :, o])
        res[o] = (b.T, resid)
    return dict(logK=logI0 - den, b=res[2][0], resid=[res[o][1] for o in range(3)], b_orders=[res[o][0] for o in range(3)])


_DB = os.path.join(HERE, 'cache', 'p2_shadow.pkl')


def shadow_odd(Delta, J):
    """cached: (b_full (c x b over [even..., odd]), logK, resid) at xi = 0.5, 8 triples, seed 7."""
    db = pickle.load(open(_DB, 'rb')) if os.path.exists(_DB) else {}
    key = (float(Delta), int(J), T.DEG)
    if key not in db:
        r = shadow_all(float(Delta), int(J))
        db = pickle.load(open(_DB, 'rb')) if os.path.exists(_DB) else {}
        db[key] = (r['b'], r['logK'], r['resid'])
        tmp = _DB + '.tmp%d' % os.getpid(); pickle.dump(db, open(tmp, 'wb')); os.replace(tmp, _DB)
    return db[key]


if __name__ == '__main__':
    out = []
    def pr(s=""):
        print(s, flush=True); out.append(s)
    pr("# p2: shadow coefficients [even..., odd] from the x5 saddle to O(1/Delta^2); beta = exp(logK) b.  Checks: fit residual per order,")
    pr("#     xi = 0.4 vs 0.6, polarisation seed 7 vs 11, even-odd entries, phase of b_OO (relative to the even b at even J)")
    Js = [int(a) for a in sys.argv[1:]] or [0, 2, 4, 5, 7, 9]
    for J in Js:
        for Dl in (100.0, 300.0, 1000.0):
            A = shadow_all(Dl, J, 0.4, seed=7); B = shadow_all(Dl, J, 0.6, seed=7); Cc = shadow_all(Dl, J, 0.4, seed=11)
            bA = A['b']; bB = np.exp(B['logK'] - A['logK'])*B['b']; bC = np.exp(Cc['logK'] - A['logK'])*Cc['b']
            n = bA.shape[0]; o = n - 1
            dxi = np.abs(bA[o, o] - bB[o, o])/abs(bA[o, o]); dpol = np.abs(bA[o, o] - bC[o, o])/abs(bA[o, o])
            eo = (np.abs(bA[:o, o]).max() + np.abs(bA[o, :o]).max())/np.abs(bA).max() if n > 1 else 0.0
            ph = np.angle(np.exp(A['logK'])*bA[o, o])/np.pi
            phE = np.angle(np.exp(A['logK'])*np.linalg.det(bA[:o, :o]))/np.pi if n > 1 else np.nan
            pr(f"  J={J} Delta={Dl:6.0f}: resid (orders 0/1/2) {A['resid'][0]:.1e}/{A['resid'][1]:.1e}/{A['resid'][2]:.1e}; odd: "
               f"|d xi| {dxi:.1e}, |d pol| {dpol:.1e}; even-odd/max {eo:.1e}; arg(beta_OO)/pi {ph:+.6f}; arg(det beta_EE)/pi {phE:+.6f}")
    open(os.path.join(HERE, 'cache', 'p2_shadow_odd.txt'), 'w').write("\n".join(out) + "\n")
