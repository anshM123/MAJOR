"""Referee: recompute A, B_AdS, C1 and the Delta_g=100 number by plain quadrature + fit (no sympy closed forms),
and quantify sensitivity of Omega to O(G) non-polynomial Einstein residues in C-_3 and the forward terms."""
import mpmath as mp
mp.mp.dps = 30
psi2 = lambda q: (1-q)**3*q*(65+155*q+47*q**2)
psi3 = lambda q: 2*(1-q)**3*(81-1927*q+4601*q**2)
k0 = 20 - 2*mp.quad(lambda q: psi2(q)*q**6, [0, 1]) + 2*mp.quad(lambda q: psi3(q)*q**4, [0, 1])
print('k(0) =', k0, ' 330151/12012 =', mp.mpf(330151)/12012)
def L(Lam):
    a = mp.mpf(3)/(2*Lam)
    return mp.quad(lambda q: psi2(q)/(q*q+a*a), [0, a, 10*a, 1])
# fit L(Lam) - 65 log Lam = c0 + c1/Lam + (d log Lam + e)/Lam^2 at large Lam
Ls = [mp.mpf(10)**4, 2*mp.mpf(10)**4, 4*mp.mpf(10)**4, 8*mp.mpf(10)**4]
rows = [[1, 1/x, mp.log(x)/x**2, 1/x**2] for x in Ls]
rhs = [L(x) - 65*mp.log(x) for x in Ls]
sol = mp.lu_solve(mp.matrix(rows), mp.matrix(rhs))
print('fit: c0 =', mp.nstr(sol[0], 12), ' expected 65 log(2/3) - 1171/15 =', mp.nstr(65*mp.log(mp.mpf(2)/3) - mp.mpf(1171)/15, 12))
print('     c1 =', mp.nstr(sol[1], 10), ' expected 30 pi =', mp.nstr(30*mp.pi, 10))
A = 16*65/k0; B = 16*sol[0]/k0; C1 = 16*sol[1]/k0
print('A = %s  B = %s  C1 = %s' % (mp.nstr(A, 8), mp.nstr(B, 8), mp.nstr(C1, 8)))
for Dg in [100, 1000]:
    Lam = mp.sqrt((Dg-2)*(Dg-1))
    exact = 16*L(Lam)/k0
    asym = A*mp.log(Lam) + B + C1/Lam
    print('Dg=%d Lam=%.3f  16 L_AdS/k0 exact = %s, A log+B+C1/L = %s  -> sqrt(t4^2+aT^2) <= %s (exact L) / %s (asym)' %
          (Dg, Lam, mp.nstr(exact, 8), mp.nstr(asym, 8), mp.nstr(240*mp.sqrt(exact)/Lam**4, 5), mp.nstr(240*mp.sqrt(asym)/Lam**4, 5)))
    # linear c1 t4 term: -Omega ⊃ -23690 pi G c1 x cos(phi); shift of x^2 bound for |c1| = 1:
    x0 = mp.sqrt(asym)
    print('    |c1|=1 linear term shifts x^2 by ~ %s  (vs retained C1/Lam = %s)' % (mp.nstr(23690/k0*x0/Lam**2, 4), mp.nstr(C1/Lam, 4)))
    # sensitivity: an O(G) Einstein residue h(nu) = 1/(nu^2+9/4)^3 in -C-_3 (allowed by the flat-limit lemma (F)):
    h = lambda nu: 1/(nu**2+mp.mpf(9)/4)**3
    term3 = mp.mpf(1)/2*Lam**6*mp.quad(lambda q: psi3(q)*h(Lam*q), [0, 1/Lam, 10/Lam, 1])
    lead = 16*L(Lam)*Lam**2
    print('    residue G/(nu^2+9/4)^3 in -C-_3 contributes %s * G to -Omega vs Einstein 16 G Lam^2 L = %s * G (ratio %s)' %
          (mp.nstr(term3, 4), mp.nstr(lead, 4), mp.nstr(term3/lead, 4)))
    h4 = lambda nu: 1/(nu**2+mp.mpf(9)/4)**4
    dh = mp.diff(lambda s: h4(mp.sqrt(s)), 0.0) if False else -4*(mp.mpf(4)/9)**5
    print('    residue G/(nu^2+9/4)^4 in -C+_4: -5 Lam^10 d_{nu^2} -> %s * G (ratio %s)' % (mp.nstr(5*Lam**10*dh, 4), mp.nstr(5*Lam**10*dh/lead, 4)))
