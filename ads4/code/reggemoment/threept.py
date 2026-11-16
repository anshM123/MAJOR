"""
d = 3 conformal three-point structures <O_1(x1,z1) O_2(x2,z2) O_5(x5,z5)> in physical space (metric (-,+,+)),
embedding-formalism building blocks of Costa-Penedones-Poland-Rychkov (arXiv:1107.3554):

   V_1 := V_{1,25} = [ (z1.x13') ... ]  -- explicitly, with x_ij = x_i - x_j, X_ij = x_ij^2:
   V_{i,jk} = [ (z_i.x_ik) X_ij - (z_i.x_ij) X_ik ] / X_jk
   H_ij     = (z_i.z_j) X_ij - 2 (z_i.x_ij)(z_j.x_ij)
   eps_ij   = det[Z_i, Z_j, X_1, X_2, X_5]  (5x5, embedding X = (1, x^2, x), Z = (0, 2 x.z, z))

   parity even:  V_1^{m1} V_2^{m2} V_5^{m5} H_12^{n12} H_15^{n15} H_25^{n25} / (|X_12|^{e12} |X_15|^{e15} |X_25|^{e25}),
                 m_i = l_i - sum_j n_ij,   e12 = (t1+t2-t5)/2, e15 = (t1+t5-t2)/2, e25 = (t2+t5-t1)/2,  t_i = Delta_i + l_i
   parity odd :  eps_15 (|X_15|/(|X_12||X_25|))^{1/2} x (even monomial with the spins at 1 and 5 lowered by one), same e's.

All separations here are timelike; absolute values |X_ij| are used (CHMRS (3.7)).  The spin l_5 may be non-integer
(Lorentzian shadow operator): then |V_5|^{m5} is used.  Everything is evaluated with mpmath; logarithms are returned
for the scalar factor so that Delta, l ~ 10^3 - 10^5 pose no problem.

Conservation at point 1 (spin 1 or 2, Delta_1 = l_1 + 1) is imposed numerically: divergence by mp.diff.
"""
import mpmath as mp
import itertools

def dot(u, v):
    return -u[0]*v[0] + u[1]*v[1] + u[2]*v[2]
def sub(u, v):
    return [u[i]-v[i] for i in range(3)]
def add(u, v, cu=1, cv=1):
    return [cu*u[i]+cv*v[i] for i in range(3)]

def embX(x):
    return [1, dot(x, x), x[0], x[1], x[2]]
def embZ(x, z):
    return [0, 2*dot(x, z), z[0], z[1], z[2]]

def Vijk(zi, xi, xj, xk):
    xij, xik, xjk = sub(xi, xj), sub(xi, xk), sub(xj, xk)
    return (dot(zi, xik)*dot(xij, xij) - dot(zi, xij)*dot(xik, xik))/dot(xjk, xjk)
def Hij(zi, xi, zj, xj):
    xij = sub(xi, xj)
    return dot(zi, zj)*dot(xij, xij) - 2*dot(zi, xij)*dot(zj, xij)

def blocks(x1, z1, x2, z2, x5, z5):
    B = dict(V1=Vijk(z1, x1, x2, x5), V2=Vijk(z2, x2, x5, x1), V5=Vijk(z5, x5, x1, x2),
             H12=Hij(z1, x1, z2, x2), H15=Hij(z1, x1, z5, x5), H25=Hij(z2, x2, z5, x5),
             X12=dot(sub(x1, x2), sub(x1, x2)), X15=dot(sub(x1, x5), sub(x1, x5)), X25=dot(sub(x2, x5), sub(x2, x5)))
    return B

def eps15(x1, z1, x2, x5, z5):
    return mp.det(mp.matrix([embZ(x1, z1), embZ(x5, z5), embX(x1), embX(x2), embX(x5)]))

def monomials(l1, l2):
    """exponents (n12, n15, n25) of the parity-even monomials (l5 assumed large)."""
    out = []
    for n12 in range(0, min(l1, l2)+1):
        for n15 in range(0, l1-n12+1):
            for n25 in range(0, l2-n12+1):
                out.append((n12, n15, n25))
    return out

def log_scalar(B, dims, spins):
    """log of |V_5|^{l5} / (|X_12|^{e12}|X_15|^{e15}|X_25|^{e25})."""
    (D1, D2, D5), (l1, l2, l5) = dims, spins
    t1, t2, t5 = D1+l1, D2+l2, D5+l5
    return (l5*mp.log(abs(B['V5'])) - (t1+t2-t5)/2*mp.log(abs(B['X12'])) - (t1+t5-t2)/2*mp.log(abs(B['X15']))
            - (t2+t5-t1)/2*mp.log(abs(B['X25'])))

def poly_even(mono, B, spins):
    """monomial divided by |V_5|^{l5} (so it is O(1)); the sign of V_5^{l5} is common to all monomials and dropped."""
    n12, n15, n25 = mono
    l1, l2, l5 = spins
    return (B['V1']**(l1-n12-n15)*B['V2']**(l2-n12-n25)*B['H12']**n12*B['H15']**n15*B['H25']**n25
            / B['V5']**(n15+n25))

def poly_odd(k, B, e15, spins):
    """odd structure number k (spin-1 or spin-2 operator at point 1, scalar at 2), divided by |V_5|^{l5}:
       eps_15 sqrt(|X15|/(|X12||X25|)) V_1^{l1-1-k} H_15^k V_5^{-1-k}."""
    l1, l2, l5 = spins
    return e15*mp.sqrt(abs(B['X15'])/(abs(B['X12'])*abs(B['X25'])))*B['V1']**(l1-1-k)*B['H15']**k/B['V5']**(1+k)

E3 = [[mp.mpf(1), 0, 0], [0, mp.mpf(1), 0], [0, 0, mp.mpf(1)]]
ETA = [-1, 1, 1]

def tensor_at_1(F, l1):
    """symmetric coefficients of the polynomial z1 -> F(z1) of degree l1 (1 or 2): returns s_mu or s_{mu nu}."""
    if l1 == 1:
        return [F(E3[m]) for m in range(3)]
    s = [[0]*3 for _ in range(3)]
    d = [F(E3[m]) for m in range(3)]
    for m in range(3):
        s[m][m] = d[m]
        for n in range(m+1, 3):
            s[m][n] = s[n][m] = (F(add(E3[m], E3[n])) - d[m] - d[n])/2
    return s

def divergence_at_1(Sfun, x1, l1):
    """Sfun(x1, z1) -> value of the full structure.  Returns the conservation defect:
       spin 1: d_mu s^mu;   spin 2: the vector  3 d_mu s^mu_nu - d_nu tr s   (Todorov operator in d = 3)."""
    def comp(x):
        s = tensor_at_1(lambda z: Sfun(x, z), l1)
        return s
    if l1 == 1:
        tot = 0
        for m in range(3):
            f = lambda t, m=m: comp([x1[i] + (t if i == m else 0) for i in range(3)])[m]
            tot += ETA[m]*mp.diff(f, 0)
        return [tot]
    out = []
    for nu in range(3):
        tot = 0
        for m in range(3):
            f = lambda t, m=m: comp([x1[i] + (t if i == m else 0) for i in range(3)])[m][nu]
            tot += 3*ETA[m]*mp.diff(f, 0)
        g = lambda t: (lambda s: sum(ETA[m]*s[m][m] for m in range(3)))(comp([x1[i] + (t if i == nu else 0) for i in range(3)]))
        tot -= mp.diff(g, 0)
        out.append(tot)
    return out

def nullspace(rows, ncols, tol=None):
    """numerical null space via SVD (mpmath). Returns (list of null vectors, singular values)."""
    A = mp.matrix(rows)
    U, S, V = mp.svd_c(A) if any(isinstance(x, mp.mpc) for r in rows for x in r) else mp.svd_r(A)
    sv = [S[i] for i in range(len(S))]
    return sv, V
