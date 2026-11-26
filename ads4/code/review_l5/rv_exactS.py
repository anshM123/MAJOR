"""rv_exactS.py (review): exact |S^(sigma+i tau)| = |d^beta S_n(beta)| by the spin recurrence at fixed N (acb balls), for L ~ 5 states,
compared with rt_pp.PP0 (b=0) / rt_pp.PP (b=4, helicity) bounds.  Usage: python3 rv_exactS.py"""
import sys, math, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, acb, ctx
import rt_pp as PP
ctx.prec = 160

def exact_s1(n, a, sig, tau):
    N = arb(n) + arb(a); al = acb(sig, tau)
    Sm, S = acb(1), acb(1)  # midpoint arithmetic (balls blow up); precision check by rerun at 2x prec
    for k in range(n):
        eps = 2*(2*k + 1)/((N - k)*(N + k + 1))
        Sm, S = S, (Sm + eps*al*S).mid()
    d = 1 - arb(n*(n + 1))/(N*(N + 1))
    return abs(d**al*S), d

def exact_gen(n, a, s, sig, tau):
    """dual-Hahn recurrence in n at fixed N: alpha S_k = B_k S_{k+1} - (B_k + D_k) S_k + D_k S_{k-1}."""
    N = arb(n) + arb(a); al = acb(sig, tau)
    Sm = acb(0); S = acb(1)
    for k in range(n):
        Bk = (N + k + s)*(k + s)*(N - k)/((2*k + s)*(2*k + s + 1))
        Dk = -k*(N + k + s)*(N - k)/((2*k + s - 1)*(2*k + s)) if k > 0 else arb(0)
        Sm, S = S, ((al*S + (Bk + Dk)*S - Dk*Sm)/Bk).mid()
    d = 1 - arb(n*(n + s))/(N*(N + s))
    return abs(d**al*S), d

def direct(n, a, s, al):
    # 3F2(al, -n, n+s; -N, N+s; 1) directly (small n)
    N = arb(n) + arb(a); t = acb(1); tot = acb(1)
    for k in range(n):
        t = t*(al + k)*(-n + k)*(n + s + k)/((-N + k)*(N + s + k)*(k + 1))
        tot += t
    return tot

if __name__ == '__main__':
    # sanity: recurrences vs direct sum
    for (n, a, s) in ((7, 4.5, 1), (9, 6.25, 9), (12, 5, 9)):
        al = acb(1.3, 2.7)
        N = arb(n) + arb(a)
        v1 = direct(n, a, s, al)
        d = 1 - arb(n*(n + s))/(N*(N + s))
        v2, _ = exact_gen(n, a, s, 1.3, 2.7)
        print('check', n, a, s, abs(d**al*v1), v2, flush=True)
    mode = sys.argv[1] if len(sys.argv) > 1 else 'tp'
    rows = []
    if mode == 'tp':
        # L ~ 5 and 5.5 states: N/a ~ 2 e^L / (1 + ...)
        for L in (5.0, 5.5):
            for a in (998, 1000, 1015, 1031, 1060, 1150, 1500):
                N = 2*a/math.exp(-L)        # a/N ~ d/2
                n = int(N - a)
                for (sig, tau) in ((2.5, 1031.0), (2.5, 1031.5), (3.5, 1032.1), (0.25, 1000.0), (1.25, 1040.0), (2.5, 1060), (2.5, 1100), (2.5, 1150), (2.5, 948.0)):
                    v, d = exact_s1(n, a, sig, tau)
                    b = PP.PP0(sig, tau, tau)
                    Lr = -math.log(float(d.mid()))
                    print(f'b=0 L={Lr:.3f} n={n} a={a} beta={sig}+{tau}i: exact {float(v.upper()):.4g}  PP0 bound {b:.4g}  ratio {float(v.upper())/b:.3g}', flush=True)
    else:
        for L in (5.0,):
            for a in (998, 1031):
                N = 2*a/math.exp(-L); n = int(N - a)
                for (sig, tau) in ((2.5, 1031.0), (3.5, 1031.5), (0.25, 900.0), (2.5, 940.0)):
                    v, d = exact_gen(n, a, 9, sig, tau)
                    b = PP.PP(sig, tau, tau, float(d.mid())*1.0001)
                    print(f'b=4 L={-math.log(float(d.mid())):.3f} n={n} a={a} beta={sig}+{tau}i: exact {float(v.upper()):.4g} [rad {float(v.rad()):.1e}]  PP bound {b:.4g}', flush=True)
