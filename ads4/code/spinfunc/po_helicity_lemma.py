"""
po_helicity_lemma (PARITY_ODD_CLOSE.md sect. 1.2): the 3-point origin of the OO <-> ++ identity.
LEMMA (tested here): let x1, x2, x5 be any (complex) points and z1 a null polarisation of T1 with V_{1,25}(z1) = 0, i.e. z1 is one of the
two null eigen-directions of the little group SO(2)_C of the three points (helicity eigen-polarisation of T1 in the 3-point frame).  Then
the symmetric conserved parity-odd structure <T1 T2 O>_odd (p1_tto_odd) restricted to z1, as a function of (z2, z5), equals
  +-i (w1 <T1T2O>_E1 + w2 <T1T2O>_E2)   (l5_spin.conserved even basis),
with w a CONSTANT (depends on Delta, J and the basis normalisation only) and the sign fixed by which of the two eigen-directions z1 is.
Tests: (1) random complex configurations, z1 constructed from V_{1,25} = 0 + nullity, both roots -> w and residual; (2) control: generic
null z1 (V1 != 0) -> no relation; (3) the crystal saddle point of the functional (z1 = x1 - x0) has V_{1,25} = 0.
Usage: python3 po_helicity_lemma.py J Delta  -> appends cache/po_helicity_lemma.txt;  python3 po_helicity_lemma.py off J Delta nu
"""
import os, sys
sys.dont_write_bytecode = True
os.environ.setdefault('L1_DEG', '6')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import warnings; warnings.filterwarnings('ignore')
import numpy as np
import po_pwtest as PT

P8 = PT.P8; T = P8.T; S = P8.S; A6 = P8.A6; Cr = P8.Cr
PW = PT.PW


def ldot(p, q):
    return -0.5*(p[0]*q[1] + p[1]*q[0]) + p[2]*q[2]


def eigen_null(x1, x2, x5):
    """the two null z with V_{1,25}(z) = 0: z . (x12/x12^2 - x15/x15^2) = 0, z^2 = 0 (lightcone components)."""
    d2 = x2 - x1; d5 = x5 - x1
    n = d2/ldot(d2, d2) - d5/ldot(d5, d5)                  # V_{1,25} prop. to z . n
    # parametrise null z = (l1^2, l2^2, l1 l2) (u, v, y); z.n = -(l1^2 n_v + l2^2 n_u)/2 + l1 l2 n_y = 0 in t = l1/l2
    a, b, c = -0.5*n[1], n[2], -0.5*n[0]
    out = []
    for t in np.roots([a, b, c]):
        out.append(np.array([t*t, 1.0, t], complex))
    return out


def structures(x1, z1, x2, x5, Delta, J):
    """F[s, n2, n5] samples on the null circles: even (l5_spin.conserved) then odd (p1, times (m12 m15 m25)^(-1/2))."""
    PW.set_batch(1)
    tp = lambda v: [np.array([complex(c)]) for c in v]
    X1, Z1, X2, X5 = tp(x1), tp(z1), tp(x2), tp(x5)
    sets = P8.struct_sets(Delta, J)
    (evD, odD), _ = sets[Delta], sets[3.0 - Delta]
    W2, WJ = S.null_circle(2), S.null_circle(J)
    Fe = S.tto_samples(X1, Z1, X2, W2, X5, WJ, J, *evD)
    msq = A6.msq
    ex = T.powr(T.mul(T.mul(msq(X1, X2), msq(X1, X5)), msq(X2, X5)), -0.5)
    Fo = S.bm(P8.tto_odd_samples(X1, Z1, X2, W2, X5, WJ, J, *odD), ex)
    F = np.concatenate([Fe, Fo], axis=0)[..., 0]
    return F


def fit(F):
    A = np.stack([F[a].reshape(-1) for a in range(F.shape[0] - 1)], axis=1); b = F[-1].reshape(-1)
    w, *_ = np.linalg.lstsq(A, b, rcond=None)
    return w, np.linalg.norm(A @ w - b)/np.linalg.norm(b)


def offsaddle(J, Dl, nu, seed=5):
    """(E) of PARITY_ODD_CLOSE: off the crystal saddle (z = z* + L xi), |V_{1,25}(x10)| != 0 and the F-level relation degrades."""
    zs, L = PT.saddle_L(nu, Dl, 1.0)
    rng = np.random.default_rng(seed)
    lines = [f'# offsaddle J={J} Delta={Dl} nu={nu}: |V1|rel = |z1.n|/(|z1||n|), F-level residual of the constant-w relation']
    for sc in (0.1, 0.4, 1.0):
        xi = sc*rng.standard_normal((4, 7))
        zl = PW.points(zs, L, xi)
        Pg = A6.geom(zl, 1.0)
        for n in range(1, 5):
            x0, x1, x2, e = (np.array([complex(c[n]) for c in Pg[k]]) for k in ('x0', 'x1', 'x2', 'e'))
            z1 = x1 - x0
            nn = (x2 - x1)/ldot(x2 - x1, x2 - x1) - (e - x1)/ldot(e - x1, e - x1)
            V1 = abs(ldot(z1, nn))/np.linalg.norm(z1)/np.linalg.norm(nn)
            w, r = fit(structures(x1, z1, x2, e, Dl, J))
            lines.append(f'  scale {sc} pt {n}: |V1|rel {V1:.2e}  F-level residual {r:.2e}  w {np.round(w, 4)}')
    return lines


if __name__ == '__main__' and sys.argv[1] == 'off':
    L_ = offsaddle(int(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]))
    print('\n'.join(L_))
    open(os.path.join(HERE, 'cache', 'po_helicity_lemma.txt'), 'a').write('\n'.join(L_) + '\n')
elif __name__ == '__main__':
    J, Dl = int(sys.argv[1]), float(sys.argv[2])
    rng = np.random.default_rng(7)
    lines = [f'# po_helicity_lemma J={J} Delta={Dl}']
    for trial in range(4):
        x1, x2, x5 = [rng.standard_normal(3) + 1j*rng.standard_normal(3) for _ in range(3)]
        for k, z1 in enumerate(eigen_null(x1, x2, x5)):
            w, r = fit(structures(x1, z1, x2, x5, Dl, J))
            lines.append(f'  config {trial} root {k}: V1=0 eigen-null z1 -> w = {np.round(w, 6)}, residual {r:.1e}')
        zg = np.array([1.3 + 0.2j, 0.7 - 0.4j, 0.0]); zg[2] = np.sqrt(zg[0]*zg[1])
        w, r = fit(structures(x1, zg, x2, x5, Dl, J))
        lines.append(f'  config {trial} control (generic null z1): w = {np.round(w, 4)}, residual {r:.1e}')
    # the functional's saddle (n0 = +-1): z1 = x1 - x0
    for n0 in (1.0, -1.0):
        zs, L = PT.saddle_L(30.0, Dl, n0)
        zl = PW.points(zs, L, np.zeros((0, 7)))
        Pg = A6.geom(zl, n0)
        x0, x1, x2, e = (np.array([complex(c[0]) for c in Pg[k]]) for k in ('x0', 'x1', 'x2', 'e'))
        z1 = x1 - x0
        n = (x2 - x1)/ldot(x2 - x1, x2 - x1) - (e - x1)/ldot(e - x1, e - x1)
        w, r = fit(structures(x1, z1, x2, e, Dl, J))
        lines.append(f'  crystal saddle n0={n0:+.0f} (nu = 30): z1.z1 = {abs(ldot(z1, z1)):.1e}, |z1.n|/(|z1||n|) = '
                     f'{abs(ldot(z1, n))/np.linalg.norm(z1)/np.linalg.norm(n):.1e} -> w = {np.round(w, 6)}, residual {r:.1e}')
    print('\n'.join(lines))
    open(os.path.join(HERE, 'cache', 'po_helicity_lemma.txt'), 'a').write('\n'.join(lines) + '\n')
