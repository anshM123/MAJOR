"""Step 4 (proof): extreme-component detector amplitudes M_{-J} (bra polarisation e_+/2) and M_{+J} (bra e_-/2) of the
canonical structures for SYMBOLIC J, by the chain (F1)-(F3) of KERNELS_D3 §3.2, with the V4^{J-k} factor expanded
binomially and the resulting single sums evaluated by Chu-Vandermonde:
    sum_j C(n,j) j^{(r)} / (Gamma(x-j) Gamma(y+j)) = n^{(r)} Gamma(x+y+n-r-1) / (Gamma(x+y-1) Gamma(x-r) Gamma(y+n)).
Source T(e_+) (helicity +2), detector direction theta = 0, rest frame p = (1,0,0).
Output: M_{-J}, M_{+J} / (universal constant U = -4 pi i e^{i pi e14} pi^{5/2}) as closed forms in (s, J)."""
import pickle, itertools, sys, time
from common import *
jj, nn, NN, AA, gg = sp.symbols('j n N A g')

def dot(u, v): return -u[0]*v[0] + u[1]*v[1] + u[2]*v[2]
def det3(a, b, c): return sp.Matrix([a, b, c]).det()
def covec(v): return [-v[0], v[1], v[2]]
def eps_covec(a, b): return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]
def pair(w1, w2): return -w1[0]*w2[0] + w1[1]*w2[1] + w1[2]*w2[2]
def cov_on(w, v): return w[0]*v[0] + w[1]*v[1] + w[2]*v[2]
EP = [0, 1, sp.I]; EM = [0, 1, -sp.I]
P = [1, 0, 0]; Z = [1, 1, 0]; Z3 = EP

class Kin:
    def __init__(self, which):
        self.zeta = [sp.Rational(1, 2)*c for c in (EP if which == 'minusJ' else EM)]
        zt = self.zeta
        self.zz3, self.zze, self.z3ze = dot(Z, Z3), dot(Z, zt), dot(Z3, zt)
        self.E0 = sp.expand(det3(Z3, zt, Z))
        self.forms = {1: covec(Z3), 2: covec(zt), 'e1': eps_covec(Z3, Z), 'e2': eps_covec(zt, Z), 'e3': eps_covec(Z3, zt)}
        one = sp.Integer(1)
        # keys (i,j,k,l,n1,n2,eps): P3^i P4^j (x^2)^k (z.x)^l (z3.x)^n1 (zeta.x)^n2 [eps-form]   (pipeline.Kin.blocks)
        self.B = dict(
            V1={(0, 0, -1, 1, 0, 0, 0): -one},
            V3={(1, -1, 0, 0, 1, 0, 0): one, (0, -1, 1, 0, 0, 0, 0): self.zz3},
            H13={(0, 0, 0, 0, 0, 0, 0): self.zz3},
            H14={(0, 0, 0, 0, 0, 0, 0): self.zze},
            H34={(0, 0, 1, 0, 0, 0, 0): self.z3ze, (0, 0, 0, 0, 1, 1, 0): -2*one},
            e13={(1, 0, 0, 0, 0, 0, 1): -one}, e14={(0, 1, 0, 0, 0, 0, 2): -one},
            e34={(0, 0, 1, 0, 0, 1, 1): 2*one, (1, 0, 1, 0, 0, 0, 3): one, (0, 0, 2, 0, 0, 0, 0): self.E0})

def emul(A, B):
    out = {}
    for ka, ca in A.items():
        for kb, cb in B.items():
            assert not (ka[6] and kb[6])
            k = tuple(x + y for x, y in zip(ka[:6], kb[:6])) + (ka[6] or kb[6],)
            out[k] = sp.expand(out.get(k, 0) + ca*cb)
    return {k: c for k, c in out.items() if c != 0}

def fixed_part(mr, kin):
    m1, m3, k, n13, n14, n34, eps = mr
    R = {(0, 0, 0, 0, 0, 0, 0): sp.Integer(1)}
    for nm, pw in (('V1', m1), ('V3', m3), ('H13', n13), ('H14', n14), ('H34', n34)):
        for _ in range(pw):
            R = emul(R, kin.B[nm])
    if eps:
        R = emul(R, kin.B['e' + eps])
    return R

def Dpoly(kin, n1, eps, n2, A0N, g):
    """Fourier derivative factor, j-polynomial part.  [prod (i w.d_p)^{n}] (-p^2)^g (-4 z.p)^{-A} at p = (1,0,0), with the
    direction zeta to the power N = n2 + j (zeta null, zeta.p = 0), z3 to the power n1, the eps form to the power 0/1:
        = i^{N+ns} prod(n_k!) 4^{-A} sum_q C(g,q)(-1)^q sum_{X^q monomials c tu^r ts^d}
              ff(N, r) (A)_{N-r+R} Lu^{N-r} prod Ls^{rem}/rem!
    Returns (P(j), ns, prod n_k!) where P(j) = sum ... ff(N,r) Gamma(A+N-r+R) Lu^{n2-r} prod Ls^rem/rem!  [the factors
    1/Gamma(A), i^N, 4^{-A}, Lu^j are handled by the caller].  A0N = A + N (independent of j)."""
    small = [(kin.forms[1], n1)] + ([(kin.forms['e%d' % eps], 1)] if eps else [])
    small = [(w, k) for w, k in small if k > 0]
    u = kin.forms[2]
    assert cov_on(u, P) == 0 and sp.expand(pair(u, u)) == 0
    K = len(small)
    tu = sp.Symbol('tu'); ts = sp.symbols('ts0:%d' % max(K, 1))
    W = [(tu, u)] + [(ts[q], small[q][0]) for q in range(K)]
    X = sp.expand(sum(2*t*cov_on(w, P) for t, w in W) + sum(t1*t2*pair(w1, w2) for t1, w1 in W for t2, w2 in W))
    nsmall = [small[q][1] for q in range(K)]
    ns = sum(nsmall)
    Lu = sp.expand(cov_on(u, Z)); Ls = [sp.expand(cov_on(w, Z)) for w, _ in small]
    N = n2 + jj
    total = 0
    for q in range(ns + 1):
        Xq = sp.Poly(sp.expand(X**q), tu, *ts[:K])
        for monom, c in Xq.terms():
            r, d = monom[0], monom[1:]
            if any(d[i] > nsmall[i] for i in range(K)):
                continue
            rem = [nsmall[i] - d[i] for i in range(K)]
            term = sp.binomial(g, q)*(-1)**q*c*sp.ff(N, r)*sp.gamma(A0N - r + sum(rem))*Lu**(n2 - r)
            for i in range(K):
                term *= Ls[i]**rem[i]/sp.factorial(rem[i])
            total += term
    return sp.expand(total), ns, sp.Mul(*[sp.factorial(k) for k in nsmall]), Lu

def chu(x, y, n, r):
    """sum_j C(n,j) j^{(r)} / (Gamma(x-j) Gamma(y+j))"""
    return sp.ff(n, r)*sp.gamma(x + y + n - r - 1)/(sp.gamma(x + y - 1)*sp.gamma(x - r)*sp.gamma(y + n))

def mono_extreme(mr, kin, parity):
    """M/U for one monomial (relative label; V4 power n = J - k), symbolic J and s"""
    from sympy.functions.combinatorial.numbers import stirling
    half = sp.Rational(1, 2) if parity == 'odd' else 0
    e13, e14 = 5 - S + half, S + half
    esum = 6 if parity == 'odd' else 5
    n = J - mr[2]
    zze = kin.zze
    total = 0
    for (iT, jT, kT, lT, n1, n2, eps), cT in fixed_part(mr, kin).items():
        c = esum - 1 - iT - jT
        # j-summand = C(n,j)/(Gamma(e14-jT-j) Gamma(B)) x K0 x Ppoly(j), with
        #   alpha integral : (-zze)^{n-j} (-1)^{iT+n} Gamma(n-j+c) / (Gamma(e14-jT-j) Gamma(e13-iT+n))
        #   Fourier        : 2^{-lT} e^{-i pi A/2} 2^{4-2B} / (Gamma(B) Gamma(g+1)),  A = n-j+c-lT, B = e14-kT-n+j
        #   derivatives    : i^{N+ns} prod(n_k!) 4^{-A} [sum ... Gamma(A+N-r+R)/Gamma(A) ...],  N = n2+j
        # j-dependent powers: (-zze)^{-j} e^{i pi j/2} 4^{-j} i^j 4^{j} Lu^j = (Lu/zze)^j = 1  (Lu = z.zeta = zze)
        # Gamma(n-j+c)/Gamma(A) = rf(n-j+c-lT, lT): polynomial in j.
        g = c - lT + e14 - kT - sp.Rational(3, 2)
        A0N = n + c - lT + n2
        Pj, ns, nfac, Lu = Dpoly(kin, n1, eps, n2, A0N, g)
        assert sp.simplify(Lu - zze) == 0
        K0 = ((-zze)**n*(-1)**(iT + n)/sp.gamma(e13 - iT + n)*sp.Integer(2)**(-lT)
              * sp.exp(-sp.I*sp.pi*(n + c - lT)/2)*sp.Integer(2)**(4 - 2*(e14 - kT - n))/sp.gamma(g + 1)
              * sp.I**(n2 + ns)*nfac*sp.Integer(4)**(-(n + c - lT)))
        poly = sp.Poly(sp.expand(sp.rf(n - jj + c - lT, lT)*Pj), jj)
        x, y = e14 - jT, e14 - kT - n
        for (p_,), cp in poly.terms():
            for r in range(p_ + 1):
                st2 = stirling(p_, r)
                if st2:
                    total += cT*K0*cp*st2*chu(x, y, n, r)
    return total

if __name__ == '__main__':
    pass
def structure_extreme(vec, which, parity):
    kin = Kin(which)
    return sum(c*mono_extreme(mr, kin, parity) for mr, c in vec.items() if c != 0)

def num(expr, sv, Jv):
    import mpmath as mp
    mp.mp.dps = 30
    e = expr.subs({J: Jv, S: sp.Rational(sv)})
    return complex(sp.N(e, 30))

if __name__ == '__main__':
    import mpmath as mp, kernels_lib as kl, canonical
    mp.mp.dps = 30
    structs = pickle.load(open(os.path.join(CACHE, 'symbolic_structures.pkl'), 'rb'))
    names = sys.argv[1:] or ['pp', 'pm', 'odd', 'oddJ']
    res = {}
    for name in names:
        par = 'even' if name in ('pp', 'pm') else 'odd'
        for which in ('minusJ', 'plusJ'):
            t0 = time.time() if 'time' in dir() else 0
            tot = structure_extreme(structs[name], which, par)
            res[(name, which)] = tot
            # numerical validation against the read-only pipeline at fixed J
            Jvals = [4, 6, 14] if name != 'oddJ' else [5, 15]
            for Jv in Jvals:
                try:
                    can = canonical.canonical(Jv)
                except FileNotFoundError:
                    continue
                parity, monos, vecJ, K, Pm = can[name]
                m = -Jv if which == 'minusJ' else Jv
                for sv in (sp.Rational(Jv) + sp.Rational(53, 10), sp.Rational(Jv) + sp.Rational(1007, 100)):
                    pv = sp.N(sp.sympify(Pm.get(m, 0)).subs(kl.S, sv), 30)
                    pipe = kl.G_prefactor(parity, K, mp.mpf(sv.p)/sv.q)*mp.mpc(str(sp.re(pv)), str(sp.im(pv)))
                    half = sp.Rational(1, 2) if par == 'odd' else 0
                    mine = complex(sp.N((-4*sp.pi*sp.I*sp.exp(sp.I*sp.pi*(S + half))*sp.pi**sp.Rational(5, 2)*tot)
                                        .subs({J: Jv, S: sv}), 30))
                    print(f"  {name} {which} J={Jv} s={sv}: pipeline {complex(pipe):.6e}  symbolic-J {mine:.6e}")
    pickle.dump(res, open(os.path.join(CACHE, 'extreme_' + '_'.join(names) + '.pkl'), 'wb'))
