"""Direct test of the renormalised Crofton identity (FULLCONE.md Prop. 1) on the hyperbolic cylinder
X = H^2/<z -> lam z>, lam = e^ell, marked points on the boundary component (0,inf):
   l^eps(c) = sum_{tau != tau_p, nontrivial} m(tau) i(tau,c) + sum_p e_p(c) log(1/(eps d_p))
for single arcs c = (x_i, lam^m x_k) (winding m times) and for the core geodesic.  The class sum is
truncated at shift M (masses decay like lam^-M; the truncation error is printed)."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np

def linked(a, b, c, d):
    lo, hi = min(a, b), max(a, b); return (lo < c < hi) != (lo < d < hi)

def logcr(a1, a2, b1, b2):  # m((a1,a2)x(b1,b2)) for a1<a2<b1<b2
    return np.log((b1-a1)*(b2-a2)/((b1-a2)*(b2-a1)))

def test(rng, M=40):
    ell = rng.uniform(1.0, 4.0); lam = np.exp(ell); n = int(rng.integers(2, 6))
    x = np.sort(rng.uniform(1, lam, size=n)); s = np.exp(rng.uniform(-1, 1, size=n)); eps = 1e-4
    X = lambda K: x[K % n]*lam**(K//n)            # lifted marked point K (integer index)
    piece = lambda K: (X(K), X(K+1))               # lifted piece K between points K, K+1
    def d_p(i):
        pts = [X(K) for K in range(i-3*n, i+3*n+1) if K != i] + [0.0]   # 0 and inf are limit points; inf -> 0
        q = [-1.0/(z - x[i]) for z in pts] + [0.0]
        return max(q) - min(q)
    # classes: (k, K) with K >= k+2 (piece K to the right of piece k), plus (k, 'inner')
    classes = []
    for k in range(n):
        for K in range(k+2, k + n*M):
            a1, a2 = piece(k); b1, b2 = piece(K)
            classes.append(((k, K), logcr(a1, a2, b1, b2), (np.sqrt(a1*a2), np.sqrt(b1*b2))))
        a1, a2 = piece(k); classes.append(((k, 'in'), np.log(a2/a1), (np.sqrt(a1*a2), -1.0)))
    out = []
    for trial in range(6):
        i = int(rng.integers(0, n)); k = int(rng.integers(0, n)); m = int(rng.integers(0, 3))
        if m == 0 and k <= i: continue
        u, v = x[i], X(k + n*m)
        lifts = [(lam**t*u, lam**t*v) for t in range(-M-5, M+5)]
        def inter(rep):
            a, b = rep
            return sum(linked(a, b, c, d) for c, d in lifts) if b > 0 else sum(linked(a, b, c, d) for c, d in lifts)
        iv = sum(mass*inter(rep) for _, mass, rep in classes)
        ell_eps = np.log((u - v)**2/((eps*s[i])*(eps*s[k]*lam**m)))
        rhs = iv + np.log(1/(eps*s[i]*d_p(i))) + np.log(1/(eps*s[k]*d_p(k)))
        out.append(ell_eps - rhs)
    # core geodesic: length ell, intersects exactly the classes going to the inner boundary, once
    core = sum(mass for (key, mass, rep) in classes if key[1] == 'in') - ell
    tail = max(mass for (key, mass, rep) in classes if key[1] != 'in' and key[1] >= key[0] + n*(M-1))
    return out, core, tail

if __name__ == "__main__":
    rng = np.random.default_rng(11); errs = []; cores = []; tails = []
    for _ in range(40):
        o, c, t = test(rng); errs += o; cores.append(c); tails.append(t)
    print(f"arcs: {len(errs)} cases, max |l^eps - (i(m',c) + sum e_p log(1/(eps d_p)))| = {max(abs(e) for e in errs):.2e}")
    print(f"core: max |sum_tau m(tau) i(tau,core) - ell| = {max(abs(c) for c in cores):.2e};  largest truncated class mass {max(tails):.1e}")
