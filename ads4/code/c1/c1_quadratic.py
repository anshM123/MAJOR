r"""
c1_quadratic (C1_COMPUTE.md §4): exact (non-linearised) c1-dependence of the Result D and Result Q bounds, valid for every real c1.

Result D (sum rule with heavy >= 0; light terms: Einstein 16 L_A, E x W3 16 c1 t4 Lam^2 I_D, |g3|^2 term K |g3|^2 Lam^8, K = k(0)(1+dk)):
    16 L_A + 16 I_D Lam^2 c1 t4 - K Lam^8 (t4^2 + aT^2)/57600 >= 0,     I_D = int psi2^NF,  sqrt(t4^2 + aT^2) = 240 |g3|.
  => R := sqrt(t4^2 + aT^2) <= R0 (x + sqrt(1 + x^2)),   R0 = 240 sqrt(16 L_A/K)/Lam^4,   x = 1920 I_D |c1| /(Lam^2 sqrt(K 16 L_A)),
     and t4 in [-R0 (sqrt(1+x^2) - x), R0 (sqrt(1+x^2) + x)] for c1 > 0 (mirror for c1 < 0);  |g3|^2 Lam^8 <= (16 L_A/K)(x + sqrt(1+x^2))^2.
Result Q (upper end U = 2 L_AdS - kappa |g3|^2 Lam^8 + 2 c1 t4 Lam^2 I_Q, kappa > 0, I_Q = int psi2^Q):
    sup over (t4, aT) of U = 2 L_AdS + 57600 I_Q^2 c1^2 /(kappa Lam^4)   (no input from D; attained at t4* = 57600 I_Q c1/(kappa Lam^6)),
    and with D's exact t4-interval imposed: the constrained supremum (endpoint value if |t4*| exceeds D's bound).
Usage: python3 c1_quadratic.py -> cache/c1_quadratic.txt (seconds).  Imports dgscan/bounds.py read-only (exact psi2 polynomials, L integrals).
"""
import os, sys, io
os.environ['OMP_NUM_THREADS'] = '1'
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'dgscan'))
import mpmath as mp
import bounds as B          # read-only: lam_of, L_exact, PSI2_NF, PSI2_Q, INT2_NF, INT2_Q, Qd

mp.mp.dps = 30
OUT = os.path.join(HERE, 'cache', 'c1_quadratic.txt')
buf = io.StringIO()


def pr(s=""):
    print(s, flush=True); buf.write(s + "\n")


DK = {300: mp.mpf('-7.666e-4'), 1000: mp.mpf('-6.830e-5')}       # delta k(0)_EE (RESULT_D_STATEMENT §3)
KAPPA = mp.mpf(B.Qd['kappa'])
C1S = [0, 1, -1, 10, -10, 100, -100, 1000, -1000]


def D_bounds(Dg, c1, K=1):
    Lam = B.lam_of(Dg)
    L16 = 16*B.L_exact(B.PSI2_NF, Lam)
    R0 = 240*mp.sqrt(L16/K)/Lam**4
    x = 1920*B.INT2_NF*abs(c1)/(Lam**2*mp.sqrt(K*L16))
    up = x + mp.sqrt(1 + x*x); dn = mp.sqrt(1 + x*x) - x
    return dict(Lam=Lam, L16=L16, R0=R0, x=x, R=R0*up, t4_same=R0*up, t4_opp=R0*dn, g3=L16/K*up**2, lin=R0*(1 + x))


def Q_bounds(Dg, c1, K=1):
    Lam = B.lam_of(Dg)
    L2 = 2*B.L_exact(B.PSI2_Q, Lam)
    IQ = B.INT2_Q
    eps_free = 57600*IQ**2*mp.mpf(c1)**2/(KAPPA*Lam**4)
    t4s = 57600*IQ*abs(c1)/(KAPPA*Lam**6)
    d = D_bounds(Dg, c1, K)
    T = d['t4_same']
    if t4s <= T:
        eps_con = eps_free
    else:
        eps_con = 2*abs(c1)*IQ*Lam**2*T - KAPPA*T**2*Lam**8/57600
    eps_lin = 2*abs(c1)*IQ*Lam**2*d['R0']          # old linearised |eps_c1| <= 2 I_Q Lam^2 |c1| t4*(c1 = 0), ignoring kappa
    return dict(L2=L2, IQ=IQ, eps_free=eps_free, t4s=t4s, eps_con=eps_con, eps_lin=eps_lin, T=T)


def main():
    pr("# c1_quadratic: exact c1-dependence of the Result D and Q bounds")
    pr(f"# I_D = int psi2^NF = {mp.nstr(B.INT2_NF, 8)};  I_Q = int psi2^Q = {mp.nstr(B.INT2_Q, 8)};  kappa = {mp.nstr(KAPPA, 6)};"
       f"  1920 I_D = {mp.nstr(1920*B.INT2_NF, 7)};  57600 I_Q^2/kappa = {mp.nstr(57600*B.INT2_Q**2/KAPPA, 7)}")
    for Dg in (300, 1000):
        d0 = D_bounds(Dg, 0)
        q0 = Q_bounds(Dg, 0)
        pr(f"\n## Delta_g = {Dg}: Lam = {mp.nstr(d0['Lam'], 8)}, 16 L_A = {mp.nstr(d0['L16'], 8)}, R0 = {mp.nstr(d0['R0'], 6)},"
           f" x per unit c1 = {mp.nstr(D_bounds(Dg, 1)['x'], 6)};  2 L_AdS = {mp.nstr(q0['L2'], 7)};"
           f"  Q coefficient 57600 I_Q^2/(kappa Lam^4) = {mp.nstr(57600*B.INT2_Q**2/(KAPPA*d0['Lam']**4), 6)}")
        pr("| c1 | x | D: sqrt(t4^2+aT^2) <= (exact) | ratio to c1=0 | linearised (1+x) R0 | t4 range (c1>0: [-lo, hi]) |"
           " |g3|^2 Lam^8 <= | Q: eps (no D input) | Q: eps (D-constrained) | old lin. 2 I_Q Lam^2 |c1| R0 | Q upper end 2L+eps (no D) |")
        for c1 in C1S:
            d = D_bounds(Dg, c1); q = Q_bounds(Dg, c1)
            lo, hi = (d['t4_opp'], d['t4_same']) if c1 >= 0 else (d['t4_same'], d['t4_opp'])
            pr(f"| {c1:+d} | {mp.nstr(d['x'], 4)} | {mp.nstr(d['R'], 5)} | {mp.nstr(d['R']/d0['R0'], 6)} | {mp.nstr(d['lin'], 5)} |"
               f" [-{mp.nstr(lo, 4)}, {mp.nstr(hi, 4)}] | {mp.nstr(d['g3'], 6)} | {mp.nstr(q['eps_free'], 4)} | {mp.nstr(q['eps_con'], 4)} |"
               f" {mp.nstr(q['eps_lin'], 4)} | {mp.nstr(q['L2'] + q['eps_free'], 6)} |")
        # with the delta k(0) factor
        K = 1 + DK[Dg]
        pr(f"   with K = k(0)(1+dk) = {mp.nstr(K, 7)}: R0 = {mp.nstr(D_bounds(Dg, 0, K)['R0'], 6)}; c1 = 1000: R = {mp.nstr(D_bounds(Dg, 1000, K)['R'], 6)}")
        # c1 at which the linearisation is off by 1 % and 10 %: x + sqrt(1+x^2) vs 1 + x  -> difference sqrt(1+x^2) - 1 ~ x^2/2
        per = D_bounds(Dg, 1)['x']
        for tol in (mp.mpf('0.01'), mp.mpf('0.1')):
            xs = mp.findroot(lambda X: (X + mp.sqrt(1 + X*X))/(1 + X) - 1 - tol, 0.3)
            pr(f"   linearisation error {mp.nstr(100*tol, 3)} % of the bound at x = {mp.nstr(xs, 5)}, i.e. |c1| = {mp.nstr(xs/per, 5)}")
        pr(f"   x = 1 (bound x (1+sqrt2)) at |c1| = {mp.nstr(1/per, 6)}")
    with open(OUT, 'w') as f:
        f.write(buf.getvalue())


if __name__ == '__main__':
    main()
