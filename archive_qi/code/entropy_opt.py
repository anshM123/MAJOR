"""Adversarial optimisation of linear entropy functionals over multipartite quantum states.

A functional is f(rho) = sum_R c_R S(rho_R)  (natural log), R ranging over subsets of parties.
States: rho = L L^dag / tr(L L^dag), L a D x r complex matrix (rank <= r).  Analytic gradients.
"""
import numpy as np, itertools
from scipy.optimize import minimize

class EntropyFunctional:
    def __init__(self, dims, coeffs):
        """dims: local dimensions; coeffs: dict {tuple_of_party_indices: coefficient}"""
        self.dims = tuple(dims); self.n = len(dims); self.D = int(np.prod(dims))
        self.coeffs = {tuple(sorted(k)): v for k, v in coeffs.items() if v != 0}

    def marg(self, rho, R):
        n = self.n; T = rho.reshape(self.dims + self.dims)
        keep = list(R); drop = [i for i in range(n) if i not in R]
        # einsum partial trace
        a = list(range(n)); b = [i + n for i in range(n)]
        for i in drop: b[i] = a[i]
        out = [a[i] for i in keep] + [b[i] for i in keep]
        d = int(np.prod([self.dims[i] for i in keep]))
        return np.einsum(T, a + b, out).reshape(d, d)

    def embed(self, X, R):
        """X on parties R  ->  X (x) 1 on all parties, as D x D matrix"""
        n = self.n; keep = list(R); drop = [i for i in range(n) if i not in R]
        dk = [self.dims[i] for i in keep]
        T = X.reshape(dk + dk)
        a = list(range(n)); b = [i + n for i in range(n)]
        ops = [T, [a[i] for i in keep] + [b[i] for i in keep]]
        for i in drop: ops += [np.eye(self.dims[i]), [a[i], b[i]]]
        return np.einsum(*ops, a + b).reshape(self.D, self.D)

    def entropies(self, rho):
        out = {}
        for R in self.coeffs:
            w = np.linalg.eigvalsh(self.marg(rho, R)); w = w[w > 1e-300]
            out[R] = float(-(w*np.log(w)).sum())
        return out

    def value(self, rho):
        S = self.entropies(rho); return sum(c*S[R] for R, c in self.coeffs.items())

    def value_grad(self, rho, floor=1e-14):
        """f and Hermitian G with df = tr(G drho) for trace-preserving drho"""
        f = 0.0; G = np.zeros((self.D, self.D), complex)
        for R, c in self.coeffs.items():
            m = self.marg(rho, R); w, V = np.linalg.eigh((m + m.conj().T)/2)
            wc = np.clip(w, floor, None)
            f += c*float(-(w[w > 1e-300]*np.log(w[w > 1e-300])).sum())
            G -= c*self.embed((V*np.log(wc))@V.conj().T, R)
        return f, G

def rho_from_x(x, D, r):
    L = (x[:D*r] + 1j*x[D*r:]).reshape(D, r); M = L@L.conj().T; t = np.trace(M).real
    return M/t, L, t

def make_objective(num, den, D, r, mode):
    """mode 'ratio': minimise num/den (den>0 expected);  mode 'value': minimise num."""
    def obj(x):
        rho, L, t = rho_from_x(x, D, r)
        f1, G1 = num.value_grad(rho)
        if mode == "value": f, G = f1, G1
        else:
            f2, G2 = den.value_grad(rho); f2s = f2 if abs(f2) > 1e-14 else 1e-14
            f = f1/f2s; G = G1/f2s - f1/f2s**2*G2
        Gp = G - np.trace(G@rho).real*np.eye(D)
        gL = 2*(Gp@L)/t
        return f, np.concatenate([gL.real.ravel(), gL.imag.ravel()])
    return obj

def optimise(num, den, D, r, x0, mode="ratio", maxiter=2000):
    obj = make_objective(num, den, D, r, mode)
    res = minimize(obj, x0, jac=True, method="L-BFGS-B", options=dict(maxiter=maxiter, ftol=1e-15, gtol=1e-11))
    return res.fun, res.x

# ---- Zhang-Yeung:  2 I(C:D) <= I(A:B) + I(A:CD) + 3 I(C:D|A) + I(C:D|B), parties (A,B,C,D)=(0,1,2,3)
def I2(X, Y):       return {tuple(X): 1, tuple(Y): 1, tuple(sorted(X+Y)): -1}
def I3(X, Y, Z):    return {tuple(sorted(X+Z)): 1, tuple(sorted(Y+Z)): 1, tuple(Z): -1, tuple(sorted(X+Y+Z)): -1}
def comb(*terms):
    out = {}
    for c, t in terms:
        for k, v in t.items(): out[k] = out.get(k, 0) + c*v
    return out
def zy_parts(perm=(0, 1, 2, 3)):
    A, B, C, D = [[p] for p in perm]
    lhs = comb((2, I2(C, D)))
    rhs = comb((1, I2(A, B)), (1, I2(A, C+D)), (3, I3(C, D, A)), (1, I3(C, D, B)))
    return lhs, rhs
