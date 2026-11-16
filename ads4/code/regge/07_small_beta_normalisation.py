r"""Check 7: small-beta behaviour of the kernels at fixed nu (the overlap with the bulk-point regime).
Flat side (arXiv:2201.06602 (A.2)-(A.3)): dtilde^J_{h,h'}(x -> 1) = N^J_{h,h'}, independent of p:
    dtilde_{0,0}, dtilde_{4,4} -> 1;  dtilde_{4,0} -> sqrt((J+4)!/(J-4)!)/4! ~ J^4/4!;  dtilde_{4,-4} -> (J+4)!/((J-4)! 8!) ~ J^8/8!.
AdS side, beta -> 0 at fixed nu (J ~ m beta/2):
    Phi_{0,0}, Phi_{4,4} -> 1;
    P^{-4}_{-1/2+i nu}(cosh beta) -> (beta/2)^4/4!                     [= Phi_{4,0}/(1/2+i nu)_4];
    Phi_{4,-4}/prod_{k<4}(nu^2+(k+1/2)^2) -> (beta/2)^8/8!.
So nu-independence of the beta -> 0 limit (required for agreement with the flat sum rule at all nu << m)
selects: unit normalisation for (0,0) and (4,4); P^{-4} for (4,0); Phi_{4,-4}/prod(nu^2+(k+1/2)^2) for (4,-4)."""
import mpmath as mp
mp.mp.dps = 30
def Phi_int(m, n, nu, b):
    c, sh = mp.cosh(b/2), mp.sinh(b/2)
    def f(th):
        e = mp.exp(1j*th); ep = (e*c - sh)/(c - e*sh)
        return mp.exp(-1j*m*th)*(mp.cosh(b) - mp.sinh(b)*mp.cos(th))**(-0.5 - 1j*nu)*ep**n
    return mp.quad(f, mp.linspace(0, 2*mp.pi, 17))/(2*mp.pi)
def poly(m, nu):
    r = mp.mpf(1)
    for k in range(m): r *= nu**2 + (k + mp.mpf(1)/2)**2
    return r
print("  nu   beta   Phi_44   8! Phi_4,-4/(prod (beta/2)^8)   4! P^-4/(beta/2)^4   [same with nu^8 instead of prod]")
for nu in map(mp.mpf, ('0.3', '1.0', '3.0')):
    for be in map(mp.mpf, ('0.2', '0.05')):
        p44 = mp.cosh(be/2)**8*mp.hyp2f1(4.5+1j*nu, 4.5-1j*nu, 1, -mp.sinh(be/2)**2)
        p4m4 = Phi_int(4, -4, nu, be)
        pm4 = mp.legenp(-0.5+1j*nu, -4, mp.cosh(be), type=3)
        print(f"  {mp.nstr(nu,2):>4} {mp.nstr(be,2):>5}  {mp.nstr(mp.re(p44),8):>10}  {mp.nstr(mp.re(mp.factorial(8)*p4m4/(poly(4,nu)*(be/2)**8)),8):>12}"
              f"  {mp.nstr(mp.re(24*pm4/(be/2)**4),8):>12}   [{mp.nstr(mp.re(mp.factorial(8)*p4m4/(nu**8*(be/2)**8)),6)}]")
