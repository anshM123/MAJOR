"""Symbolic check: for ds^2 = z^-2(-(1-m(v) z^2) dv^2 - 2 dv dz + dx^2) (Vaidya-AdS3, planar),
compute E_ab = R_ab + 2 g_ab (vanishes in vacuum, Lambda=-1) and its null contraction.
Also verify the 3D identity Rm(tau,k,k,tau) = Ric(k,k) used in Prop. J on a generic metric sample.
usage: python3 -I check_nec_sign.py"""
import sympy as sp
v, z, x = sp.symbols('v z x', real=True)
m = sp.Function('m')(v)
X = [v, z, x]
f = 1 - m*z**2
g = sp.Matrix([[-f, -1, 0], [-1, 0, 0], [0, 0, 1]])/z**2
gi = sp.simplify(g.inv())
n = 3
Gam = [[[sp.simplify(sum(gi[a, d]*(sp.diff(g[d, b], X[c]) + sp.diff(g[d, c], X[b]) - sp.diff(g[b, c], X[d]))
        for d in range(n))/2) for c in range(n)] for b in range(n)] for a in range(n)]
def Riem(a, b, c, d):  # R^a_{bcd}
    r = sp.diff(Gam[a][b][d], X[c]) - sp.diff(Gam[a][b][c], X[d])
    r += sum(Gam[a][c][e]*Gam[e][b][d] - Gam[a][d][e]*Gam[e][b][c] for e in range(n))
    return r
Ric = sp.Matrix(n, n, lambda b, d: sp.simplify(sum(Riem(a, b, a, d) for a in range(n))))
E = sp.simplify(Ric + 2*g)
print("E_ab = R_ab + 2 g_ab =", E)
# outgoing null vector l with g(l,l)=0: l = dv + (f/2)... solve: l=(1, w, 0): g = (-f - 2w)/z^2 = 0 -> w=-f/2
l = sp.Matrix([1, -f/2, 0])
print("E(l,l) =", sp.simplify((l.T*E*l)[0]), "  (NEC <=> this >= 0)")
k = sp.Matrix([0, 1, 0])
print("E(k,k) for k=d_z (ingoing null) =", sp.simplify((k.T*E*k)[0]))
