"""
kplus.py: finite-mass outgoing kernel K+_fm(nu) = B(nu) sum_j c_j(nu) S(alpha - j), alpha = p - 1/2 + i nu,
S(beta) = L((1-z)^{-beta}) = sum_{k<=n} U_k (beta)_k / k!  (terminating 3F2), continued by the exact three-term recurrence
[N(N+s) - beta(beta-s)] S(beta-1) = [a(N+n+s) + beta(s-2beta-1)] S(beta) + beta(beta+1) S(beta+1)       (TAIL_CERT E3)
and the infinite-twist analogue S_inf(beta) = d^{-beta}.  [DV evaluator, high precision]
p = 3 for F_{3,nu} (alpha = 5/2 + i nu), p = 4 for F_{4,nu} (alpha = 7/2 + i nu).
"""
from flint import arb, acb, ctx
I = acb(0, 1)

def S_direct(beta, U):
    t = acb(1); s = acb(1)
    for k in range(1, len(U)):
        t = t*(beta + k - 1)/k
        s += U[k]*t
    return s

def Bfac(nu, p=3):
    # B(nu) = Gamma(2 i nu)/Gamma(1/2 + i nu)^2 ; connection: (1-z)^{-p} P = 2Re[B sum_j c_j (1-z)^{j - (p - 1/2) - i nu}]
    return ((2*I*nu).lgamma() - 2*(acb(0.5) + I*nu).lgamma()).exp()

def Kplus(nu, st, U, p=3, tol=1e-60, jmax=200000, inf=False):
    nu = acb(nu)
    n, s, N, a = st['n'], st['s'], st['N'], st['a']
    d = st['d']
    al = acb(p - 0.5) + I*nu
    A = acb(0.5) - I*nu; Bb = acb(1) - 2*I*nu
    if inf:
        Sb = d**(-al); ratio = d
    else:
        Sb = S_direct(al, U); Sb1 = S_direct(al + 1, U)
    DD = N*(N + s); E = arb(a)*(N + n + s)
    cj = acb(1); tot = acb(0); j = 0
    jh = 3*float(abs(nu).upper()) + 200          # past the hump of |c_j| d^j
    mn = None; totmn = None; Kplus.jmin = 0; Kplus.minterm = 0.0
    beta = al
    while True:
        term = cj*Sb
        m = float(abs(term).upper())
        if j > 20 and m < tol*float(abs(tot).lower()):
            break
        if not inf and j > jh:
            if mn is None or m < mn:
                mn = m; totmn = tot; Kplus.jmin = j; Kplus.minterm = m
            elif m > 1e6*mn:
                tot = totmn; break
        tot += term
        cj = cj*(A + j)*(A + j)/((Bb + j)*(j + 1))
        if inf:
            Sb = Sb*d
        else:
            Sm = ((E + beta*(s - 2*beta - 1))*Sb + beta*(beta + 1)*Sb1)/(DD - beta*(beta - s))
            Sb1, Sb = Sb, Sm
        beta = beta - 1
        j += 1
        if j > jmax: raise RuntimeError('jmax')
    return Bfac(nu)*tot, j
