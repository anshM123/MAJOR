"""
Free-field (Fock space) computation of energy-flux matrix elements in d = 3.

Conventions (used throughout code/kernels):
  metric (-,+,+);  q = (E,0,0);  detector direction n(theta) = (cos theta, sin theta);
  e_pm = (0,1,+-i), e_0 = (1,0,0);  null polarisation zeta(a,b) = a^2 e_+/2 + b^2 e_-/2 + a b e_0;
  O(q;zeta(a,b)) = sum_m a^{J+m} b^{J-m} O_m(q);  U(R_phi) O_m(q)|0> = e^{-i m phi} O_m(q)|0>  (charge m);
  helicity states of T:  |T,+-2> = T_{mu nu} e_pm^mu e_pm^nu (q)|0>  (coefficients of a^4, b^4).
Two massless particles with total momentum q: k1 = (E/2)(1,n(theta)), k2 = (E/2)(1,-n(theta)).
E(theta')|theta> = (E/2)[delta(theta'-theta) + delta(theta'-theta-pi)]|theta>.

Run:  python3 fock_free.py
"""
import sympy as sp

E, th, a, b = sp.symbols('E theta a b', positive=True)
I = sp.I
k1 = sp.Matrix([E/2, E/2*sp.cos(th), E/2*sp.sin(th)])
k2 = sp.Matrix([E/2, -E/2*sp.cos(th), -E/2*sp.sin(th)])
eta = sp.diag(-1, 1, 1)
dot = lambda u, v: sp.expand((u.T*eta*v)[0])
ep, em, e0 = sp.Matrix([0, 1, I]), sp.Matrix([0, 1, -I]), sp.Matrix([1, 0, 0])
zeta = a**2*ep/2 + b**2*em/2 + a*b*e0
assert sp.simplify(dot(zeta, zeta)) == 0

def fourier_modes(F, J):
    """F(theta; a, b) homogeneous of degree 2J in (a,b). Return {m: coefficient function of theta} for a^{J+m} b^{J-m}."""
    F = sp.expand(sp.expand_trig(sp.expand(F)).rewrite(sp.exp))
    P = sp.Poly(F, a, b)
    out = {}
    for (ea, eb), c in P.terms():
        m = sp.Rational(ea - eb, 2)
        c = sp.simplify(c)
        if c != 0:
            out[m] = c
    return out

# ---------------------------------------------------------------- free real scalar
zk1, zk2 = dot(zeta, k1), dot(zeta, k2)
print("=== free real scalar ===")
# T(x;z) = (3/4)(z.d phi)^2 - (1/4) phi (z.d)^2 phi  (improved stress tensor, xi = 1/8, null z)
# creation part: d -> -i k
F_T = -(sp.Rational(3, 4)*zk1*zk2 - sp.Rational(1, 8)*(zk1**2 + zk2**2))
mT = fourier_modes(F_T, 2)
print("T wavefunction modes {m: coeff}:", mT)
assert set(mT.keys()) == {2, -2}
# spin-s currents: J_s = sum_k (-1)^k C(2s,2k) (z.d)^{s-k} phi (z.d)^k phi
def F_J_scalar(s):
    return sum((-1)**k*sp.binomial(2*s, 2*k)*zk1**(s-k)*zk2**k for k in range(s+1))*(-I)**s
for s in (0, 2, 4, 6, 8):
    ms = fourier_modes(F_J_scalar(s), s)
    print(f"J_{s} wavefunction modes:", ms)
    assert set(ms.keys()) == ({s, -s} if s else {0})
# check J_2 is proportional to T
assert sp.simplify(F_J_scalar(2)/F_T - sp.simplify(F_J_scalar(2)/F_T).subs(th, 0.3)) == 0

def analyse(modesT, label, currents):
    # wavefunction of |T,+2> is modesT[2] = cT2 * e^{2 i theta}; |T,-2> is modesT[-2] = cTm2 e^{-2 i theta}
    cp = sp.simplify(modesT[2]*sp.exp(-2*I*th)); cm = sp.simplify(modesT[-2]*sp.exp(2*I*th))
    assert cp.free_symbols <= {E} and cm.free_symbols <= {E}
    # <psi'|E(theta')|psi> = Nrm * int dtheta conj(psi') psi (E/2)[delta(theta'-theta)+delta(theta'-theta-pi)]
    # <psi'|psi> = Nrm * int dtheta conj(psi') psi   (Nrm: two-body phase space, identical-particle factor; cancels)
    thp = sp.Symbol('thetap', real=True)
    def Emat(psi_out, psi_in):
        f = sp.conjugate(psi_out)*psi_in
        return sp.simplify(E/2*(f.subs(th, thp) + f.subs(th, thp - sp.pi)))
    def norm(psi):
        return sp.simplify(sp.integrate(sp.simplify(sp.conjugate(psi)*psi), (th, 0, 2*sp.pi)))
    psi = {2: cp*sp.exp(2*I*th), -2: cm*sp.exp(-2*I*th)}
    N = {h: norm(psi[h]) for h in psi}
    Mpp = sp.simplify(Emat(psi[2], psi[2])/N[2])
    Mmp = sp.simplify(Emat(psi[-2], psi[2])/sp.sqrt(N[2]*N[-2]))
    print(f"[{label}] <T,2|E(theta)|T,2>/<T,2|T,2>  =", Mpp)
    print(f"[{label}] <T,-2|E(theta)|T,2>/norm    =", sp.simplify(Mmp), "  => c e^{4 i theta} with c =",
          sp.simplify(Mmp/Mpp*sp.exp(-4*I*thp)), " t4 = 4c")
    # normalised amplitudes onto currents:  ahat = 2 pi <O,m|E(0)|T,2> / (E sqrt(N_O N_T))
    for (s, modes) in currents:
        for m, cf in sorted(modes.items()):
            psiO = cf
            NO = norm(psiO)
            amp = sp.simplify(2*sp.pi*Emat(psiO, psi[2]).subs(thp, 0)/(E*sp.sqrt(NO*N[2])))
            print(f"[{label}]   O = spin {s}, m = {m}:  ahat = {amp},  |ahat|^2 = {sp.simplify(sp.Abs(amp)**2)}")
    # energy-energy correlator in |T,2>: product of the two delta-function factors
    print(f"[{label}] EEC: <E(t1)E(t2)>/<1> = (E^2/(4 pi)) [delta(t12) + delta(t12 - pi)];"
          f"  a_n = (1+(-1)^n)/2  (two back-to-back particles of energy E/2)")

analyse(mT, "scalar", [(s, fourier_modes(F_J_scalar(s), s)) for s in (0, 4, 6, 8)])

# ---------------------------------------------------------------- free Majorana fermion
print("=== free Majorana fermion ===")
g0 = sp.Matrix([[0, 1], [-1, 0]]); g1 = sp.Matrix([[0, 1], [1, 0]]); g2 = sp.Matrix([[1, 0], [0, -1]])
gam = [g0, g1, g2]
for mu in range(3):
    for nu in range(3):
        assert gam[mu]*gam[nu] + gam[nu]*gam[mu] == 2*eta[mu, nu]*sp.eye(2)
slash = lambda v: sum(((eta*v)[mu]*gam[mu] for mu in range(3)), sp.zeros(2))
u1 = sp.sqrt(E)*sp.Matrix([sp.sin(th/2), sp.cos(th/2)])
u2 = sp.sqrt(E)*sp.Matrix([sp.sin((th+sp.pi)/2), sp.cos((th+sp.pi)/2)])   # continuous-branch spinor of k2
assert sp.simplify(slash(k1)*u1) == sp.zeros(2, 1) and sp.simplify(slash(k2)*u2) == sp.zeros(2, 1)
bil = sp.simplify((u1.T*g0*slash(zeta)*u2)[0])      # u1^T g0 zslash u2  (g0 gamma^mu symmetric)
# T(x;z) = (i/2) psibar zslash (z.d) psi ; creation part antisymmetrised in the two fermions
F_Tf = sp.Rational(1, 2)*bil*(zk2 - zk1)*sp.Rational(1, 2)
mTf = fourier_modes(F_Tf, 2)
print("T wavefunction modes:", mTf)
assert set(mTf.keys()) == {2, -2}
# currents J_s (even s): bil * P_{s-1}(zk1, zk2), fixed by purity of helicity (= conserved primary), see KERNELS_D3.md
def F_J_fermion(s):
    # J_s = sum_k (-1)^k C(2s,2k+1) (z.d)^{s-1-k} psibar zslash (z.d)^k psi   (even s);
    # in spinor variables this is Im[(chi.lambda_1 + i chi.lambda_2)^{2s}]; the assertion below (only charges +-s)
    # is equivalent to "conserved primary" within fermion bilinears (see KERNELS_D3.md section 2).
    P = sum((-1)**k*sp.binomial(2*s, 2*k+1)*zk1**(s-1-k)*zk2**k for k in range(s))
    return sp.expand(bil*P), [(-1)**k*sp.binomial(2*s, 2*k+1)/sp.Integer(2*s) for k in range(1, s)]
cur_f = []
# parity-odd scalar psibar psi (Delta = 2): wavefunction u1^T g0 u2
F0 = sp.simplify((u1.T*g0*u2)[0])
print("psibar psi wavefunction (charge 0):", F0)
cur_f.append((0, {0: F0}))
for s in (4, 6):
    F, coeffs = F_J_fermion(s)
    ms = fourier_modes(F, s)
    print(f"J_{s} (fermion) coefficients c_k/c_0 =", coeffs, " modes:", ms)
    assert set(ms.keys()) == {s, -s}
    cur_f.append((s, ms))
analyse(mTf, "fermion", cur_f)
