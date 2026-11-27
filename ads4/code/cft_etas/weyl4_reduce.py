"""Reduce the t8t8-type quartic Weyl invariant of Garousi-Ghodsi (arXiv:0808.0411, eq. 3.2),
   W = C^h_{mnk} C_p^{mnq} C_h^{rsp} C_q^{rsk} + (1/2) C^{hk}_{mn} C_{pq}^{mn} C_h^{rsp} C_q^{rsk},
to the 4D basis {I^2, P^2} (I = C.C, P = C.C~), by least squares over random 4D Weyl tensors (Euclidean metric; the I^2
coefficient is signature independent).  Then planar 4 pi eta/s - 1 = 864 * a * gamma from our response table (I^2: 864, P^2: 0)."""
import numpy as np, itertools
rng = np.random.default_rng(1)
n = 4
def random_weyl():
    A = rng.normal(size=(n, n, n, n))
    R = A - A.transpose(1, 0, 2, 3)
    R = R - R.transpose(0, 1, 3, 2)
    R = R + R.transpose(2, 3, 0, 1)
    R = (R + np.einsum('abcd->acdb', R) + np.einsum('abcd->adbc', R))*0 + R  # placeholder
    # impose first Bianchi by projecting: R_abcd -> R - (1/3)(R_abcd + R_acdb + R_adbc)
    B = R + np.einsum('abcd->acdb', R) + np.einsum('abcd->adbc', R)
    R = R - B/3
    g = np.eye(n)
    Ric = np.einsum('abad->bd', R)
    Rs = np.trace(Ric)
    C = (R - 0.5*(np.einsum('ac,bd->abcd', g, Ric) - np.einsum('ad,bc->abcd', g, Ric) - np.einsum('bc,ad->abcd', g, Ric) + np.einsum('bd,ac->abcd', g, Ric))
         + Rs/6*(np.einsum('ac,bd->abcd', g, g) - np.einsum('ad,bc->abcd', g, g)))
    return C
eps = np.zeros((4,)*4)
for p in itertools.permutations(range(4)):
    eps[p] = np.linalg.det(np.eye(4)[list(p)])
rows, vals = [], []
for _ in range(12):
    C = random_weyl()
    assert abs(np.einsum('abad->bd', C)).max() < 1e-12
    Ct = 0.5*np.einsum('abef,efcd->abcd', eps, C)
    I = np.einsum('abcd,abcd', C, C); P = np.einsum('abcd,abcd', C, Ct)
    W = (np.einsum('hmnk,pmnq,hrsp,qrsk', C, C, C, C) + 0.5*np.einsum('hkmn,pqmn,hrsp,qrsk', C, C, C, C))
    rows.append([I**2, P**2]); vals.append(W)
coef, res, *_ = np.linalg.lstsq(np.array(rows), np.array(vals), rcond=None)
print('W = a I^2 + b P^2 with a = %.10f, b = %.10f ; max residual %.2e' % (coef[0], coef[1], np.max(np.abs(np.array(rows) @ coef - vals))))
print('our response table: I^2 -> +864, P^2 -> 0  =>  planar 4 pi eta/s - 1 = 864 a gamma = %.6f gamma  (Garousi-Ghodsi: 1728 gamma)' % (864*coef[0]))
