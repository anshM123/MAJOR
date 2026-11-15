"""Check of the Remark after Proposition P: the planes of the triangles spinning around a closed geodesic tend to a
LIGHTLIKE plane when l_L != l_R, and stay equal to one spacelike plane when l_L = l_R.

Axis: joins E21 and E12; (g_L, g_R) = (diag(e^{a/2},e^{-a/2}), diag(e^{b/2},e^{-b/2})).
Triangle k: span{E12, rho^k p, rho^{k+1} p}.  A basis of its plane that stays non-degenerate as k -> -infinity:
E12, u_k = rho^k p normalised to E21-coefficient 1, and d_k = the normalised difference of u_k and u_{k+1}.
We report the determinant of the Gram matrix of the form N(X,Y)=det(X+Y)-det X-det Y in that basis:
nonzero (negative*positive pattern (1,2)) for a spacelike plane, -> 0 for a lightlike limit."""
import numpy as np
N = lambda x, y: np.linalg.det(x+y) - np.linalg.det(x) - np.linalg.det(y)
E12 = np.array([[0, 1.], [0, 0]])
def rho(x, k, a, b):
    gL = np.diag([np.exp(a*k/2), np.exp(-a*k/2)]); gRi = np.diag([np.exp(-b*k/2), np.exp(b*k/2)])
    return gL @ x @ gRi
for a, b in ((1.0, 1.0), (1.0, 1.6), (1.6, 1.0)):
    p = np.outer([0.7, 1.0], [1.3, 1.0]); out = []
    for k in (-1, -3, -6, -10, -15):
        u = rho(p, k, a, b); u = u/u[1, 0]; v = rho(p, k+1, a, b); v = v/v[1, 0]
        d = u - v; d[0, 1] = 0.0; d = d/np.abs(d).max()          # drop the E12 component (E12 is in the plane anyway)
        B = [E12, u - u[0, 1]*E12, d]
        G = np.array([[N(x, y) for y in B] for x in B]); out.append(np.linalg.det(G))
    print(f"a={a}, b={b}: Gram determinant on the plane of triangle k (k=-1,-3,-6,-10,-15): " + ", ".join("%.3e" % x for x in out))
