"""t2: the T-block part of the bare k=3 light-ray commutator in d=3  (LIGHT_INPUTS.md Sec. 3).

Input [Q]: KKSDZ arXiv:1904.05905, eq. (energycorrelatorsbulk) and the shock matrix H_{ij,kl}(n) printed below it:
single-T intermediate states give <E(n1) E(n3)>_{eps.T(q)} = (q0/Omega)^2 eps^+ H(n1) H(n3) eps / eps^+ eps.
[D] (momentum conservation of E(n) + completeness of the T family at fixed q): the T-family contribution to
<T|[E(n1),E(n3)]|T> is (q0/Omega)^2 eps^+ [H(n1),H(n3)] eps.
Specialised to d = 3 (transverse R^2; traceless symmetric = helicities +-2). Checks:
 (a) Einstein (t2 = t4 = 0): H = identity  ->  commutator identically zero;
 (b) in d = 3 the t2 structure is not independent of identity + t4 (two parity-even <TTT> structures);
 (c) t4 enters only as a helicity flip ~ w^{+-4}, w = e^{i phi}; [H1,H3] is quadratic in the flip, no linear term,
     and helicity-odd: -2i|flip|^2 sin(4 phi_13) sigma_z, |flip| = |t4|/4 (= CDP off-diagonal; t2 drops out in d = 3).
Exact rational algebra in w (fast).
"""
import sympy as sp

d = 3
t2, t4 = sp.symbols('t2 t4', real=True)
w1, w3 = sp.symbols('w1 w3', nonzero=True)          # w = e^{i phi}

def Hmat(w):
    n = [(w + 1/w)/2, (w - 1/w)/(2*sp.I)]
    D = lambda i, j: 1 if i == j else 0
    P2 = lambda i, j: n[i]*n[j] - sp.Rational(1, d-1)*D(i, j)
    P4 = lambda i, j, k, l: n[i]*n[j]*n[k]*n[l] - (D(i,j)*D(k,l)+D(i,k)*D(j,l)+D(i,l)*D(j,k))/sp.Integer((d+1)*(d-1))
    H = lambda i, j, k, l: (sp.Rational(1, 2)*(D(i,k)*D(j,l)+D(i,l)*D(j,k)) - sp.Rational(1, d-1)*D(i,j)*D(k,l)
            + t2/4*(D(i,k)*P2(j,l)+D(i,l)*P2(j,k)+D(j,l)*P2(i,k)+D(j,k)*P2(i,l))
            + t4*P4(i,j,k,l) - (t2+t4)/(d-1)*(D(i,j)*P2(k,l)+D(k,l)*P2(i,j)))
    e = {+1: [1/sp.sqrt(2), sp.I/sp.sqrt(2)], -1: [1/sp.sqrt(2), -sp.I/sp.sqrt(2)]}
    ebar = {+1: e[-1], -1: e[+1]}
    M = sp.zeros(2, 2)
    for a, ha in enumerate((+1, -1)):
        for b, hb in enumerate((+1, -1)):
            s = sum(ebar[ha][i]*ebar[ha][j]*H(i, j, k, l)*e[hb][k]*e[hb][l]
                    for i in range(2) for j in range(2) for k in range(2) for l in range(2))
            M[a, b] = sp.factor(sp.expand(s))
    return M

if __name__ == "__main__":
    H1, H3 = Hmat(w1), Hmat(w3)
    print("t2: KKSDZ shock matrix in d = 3, helicity basis (+2,-2), w = e^{i phi}:")
    print("   H(w) =", H1)
    print("(a) Einstein t2 = t4 = 0:  H =", H1.subs({t2: 0, t4: 0}))
    print("(b) t2 alone (t4 = 0):     H =", H1.subs({t4: 0}).applyfunc(sp.factor),
          " -> the t2 structure vanishes identically in d = 3 (not an independent <TTT> structure)")
    C = (H1*H3 - H3*H1).applyfunc(lambda z: sp.factor(sp.expand(z)))
    print("(c) [H(w1), H(w3)] =", C)
    flip2 = sp.expand(H1[0, 1]*H1[1, 0])
    print("    flip^2 = H_{+-}H_{-+} =", sp.factor(flip2), " (w-independent)")
    ref = -flip2*(w1**4/w3**4 - w3**4/w1**4)         # = -2i |flip|^2 sin(4 phi13)
    print("    C[0,0] + flip^2 (w1^4/w3^4 - w3^4/w1^4) =", sp.simplify(C[0, 0] - ref), "; C[0,0]+C[1,1] =",
          sp.simplify(C[0, 0] + C[1, 1]), "; off-diagonal:", sp.simplify(C[0, 1]), sp.simplify(C[1, 0]))
    print("    linear term d/dt4 at t2 = t4 = 0:", sp.diff(C, t4).subs({t2: 0, t4: 0}))
    print("Conclusion: the T-block contribution to every bare (k = 3) commutator functional is"
          " -2i|flip|^2 sin(4phi13) sigma_z: zero for Einstein, no linear term, quadratic in the flip amplitude.")
