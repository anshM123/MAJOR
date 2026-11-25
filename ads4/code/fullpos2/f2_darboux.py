"""[DV] test of the author's Darboux reference (ANGULAR_ROUTE §46) against the exact bare kernel K_p(nu) = L(F_{p,nu})."""
import sys, os, math
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import mpmath as mp
import f2_num as F

def Kexact(p, nu, J, b, T, tol_digits=30):
    st = F.state(J, b, T); n = st['n']
    s = mp.mpf(nu)**2 + mp.mpf(1)/4
    g0, gm = mp.mpf(1), mp.mpf(0); tot = mp.mpf(0); U = st['U']
    for k in range(n + 1):
        tot += U[k]*g0
        g1 = ((2*k*k + 2*p*k + p - s)*g0 - (k + p - 1)**2*gm)/(k + 1)**2
        gm, g0 = g0, g1
    return tot, st

def Nfrak(a, al):
    return (2*a)**al*mp.beta(al, a + 1)*mp.hyp2f1(a, al, a + al + 1, -1)/mp.gamma(al)

def Bcon(nu):
    nu = mp.mpf(nu)
    return mp.gamma(1j*nu)/(mp.sqrt(mp.pi)*mp.gamma(mp.mpf(1)/2 + 1j*nu))*mp.mpf(4)**(-mp.mpf(1)/2 + 1j*nu)

def Kdarb(p, nu, J, b, Delta, order=1):
    n = J - b; N = Delta - b - 2; a = N - n
    C = J*(J + 1); T = (mp.mpf(Delta) - 1)*(Delta - 2)
    c = (C - b*(b + 1))/(T - b*(b + 1)); d = 1 - c
    M = J + mp.mpf(1)/2 + a/mp.mpf(2)
    al = p - mp.mpf(1)/2 + 1j*mp.mpf(nu)
    br = Nfrak(a, al)
    if order >= 1:
        br += a*((2*b - al)*Nfrak(a, al) + (p - 2*b)*Nfrak(a, al - 1))/M
    return 2*mp.re(Bcon(nu)*d**(-al)*br), d

if __name__ == '__main__':
    p = 3
    for (J, b) in [(int(s.split(',')[0]), int(s.split(',')[1])) for s in sys.argv[1:]]:
        a = 998; Delta = J + a + 2; T = (Delta - 1)*(Delta - 2)
        for nu in (0.5, 5, 20, 50, 100, 200, 400, 700, 1000):
            mp.mp.dps = int(0.75*nu + 60 + 0.01*J)
            ke, st = Kexact(p, nu, J, b, T)
            k0, d = Kdarb(p, nu, J, b, Delta, 0); k1, _ = Kdarb(p, nu, J, b, Delta, 1)
            lr = mp.legenp(-mp.mpf(1)/2 + 1j*nu, 0, (2 - d)/d, type=3).real*d**(-p)
            print(f'J={J} b={b} d={float(d):.4g} nu={nu}: exact={mp.nstr(ke,10)} darb0={mp.nstr(k0,10)} darb1={mp.nstr(k1,10)} '
                  f'LR(c)={mp.nstr(lr,10)}  rel(d1)={mp.nstr(abs(k1-ke)/abs(ke),3)} rel(LR)={mp.nstr(abs(lr-ke)/abs(ke),3)}', flush=True)
