"""
fullpos/fp_stripscan.py (FULLPOS_STRIPSCAN.md): UNCERTIFIED scan of the positivity margin across the low-twist bridge strip.
Gap 1000; b in {0, 4}; twist a = Delta - J - 2 in {998, 1500, 3000, 1e4, 3e4, 1e5}; J chosen so that L = -log(1 - r), r = C/T, is
as close as possible to the targets {0.04, 0.1, 0.3, 1, 3, 10, 20, 28} (J even for b = 0).  1 - r = (2aJ + a(a+1))/((J+a)(J+a+1)).
H_ang: exact finite ladder sum with the K = 1000 certified moments (so n = J - b <= 1000); H^(2): its k = 2 part (Q = 0);
H_hyb, N_s: T0's hyb_fast (float).  Beyond n = 1000 the state is skipped and the largest reachable L is reported.
"""
import sys, os, math
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'largespin'))
from flint import arb, fmpq, ctx
import fp_common as FC
import fp_poly as FP


def Lof(J, a):
    return -math.log((2*a*J + a*(a + 1))/((J + a)*(J + a + 1)))


def main():
    Dg = 1000; Kmom = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
    import ls_core as LS
    Z, Y = FP.load_moments(Dg, Kmom)
    ctx.prec = 3000
    lam = (Dg - 2)*(Dg - 1); Lam = math.sqrt(lam)
    Q = FP.Q_exact(Kmom, lam)
    A = FP.Action(Dg, Z, Y, Q)
    A2 = FP.Action(Dg, Z, Y, [fmpq(0)]*(Kmom + 1))
    out = []
    def pr(s):
        print(s, flush=True); out.append(s)
    pr('| b | a | L target | J | L | H_ang | H^(2) | H_hyb | N_s | H_ang/N_s | H^(2)/N_s |')
    pr('|---|---|---|---|---|---|---|---|---|---|---|')
    mins = {}
    for b in (0, 4):
        for a in (998, 1500, 3000, 10**4, 3*10**4, 10**5):
            Jmax = Kmom + b
            Lreach = Lof(Jmax, a)
            for Lt in (0.04, 0.1, 0.3, 1, 3, 10, 20, 28):
                # solve L(J) = Lt (L increasing in J)
                lo, hi = 0.0, 1e16
                if Lof(max(b, 4 if b else 0), a) > Lt:
                    J = max(b, 4 if b == 4 else 0)
                else:
                    for _ in range(200):
                        mid = (lo + hi)/2
                        if Lof(mid, a) < Lt: lo = mid
                        else: hi = mid
                    J = int(round(lo))
                if b == 0 and J % 2: J += 1
                if b == 4 and J < 4: J = 4
                if J - b > Kmom:
                    pr(f'| {b} | {a} | {Lt} | {J} (n = {J-b} > {Kmom}) | — | skipped | | | | | |')
                    continue
                Delta = J + a + 2; T = (Delta - 1)*(Delta - 2); C = J*(J + 1)
                Ha = float(A.H_ang(J, b, T).mid()); H2 = float(A2.H_ang(J, b, T).mid())
                hy = LS.hyb_fast(Dg, J, math.sqrt(T - C)/Lam, '++' if b == 0 else '+-')
                m1 = Ha/hy['Ns']; m2 = H2/hy['Ns']
                mins.setdefault('ang', []).append((m1, b, a, J)); mins.setdefault('h2', []).append((m2, b, a, J))
                pr(f'| {b} | {a} | {Lt} | {J} | {Lof(J, a):.4f} | {Ha:.6e} | {H2:.6e} | {hy["base"]:.6e} | {hy["Ns"]:.6e} | {m1:.5f} | {m2:.5f} |')
            pr(f'|   | {a} | (largest reachable L at n = {Kmom}: {Lreach:.4f}) | | | | | | | | |')
    if mins:
        pr('')
        pr(f'min H_ang/N_s over scanned states = {min(mins["ang"])}')
        pr(f'min H^(2)/N_s over scanned states = {min(mins["h2"])}')
    open(os.path.join(FC.CACHE, 'fp_stripscan.txt'), 'w').write('\n'.join(out) + '\n')


if __name__ == '__main__':
    main()


def edge_points():
    """the largest reachable n = 1000 for a in {998, 1500, 3000} (appended to cache/fp_stripscan_edge.txt)."""
    import ls_core as LS
    Dg = 1000; Kmom = 1000
    Z, Y = FP.load_moments(Dg, Kmom)
    ctx.prec = 3000
    lam = (Dg - 2)*(Dg - 1); Lam = math.sqrt(lam)
    A = FP.Action(Dg, Z, Y, FP.Q_exact(Kmom, lam)); A2 = FP.Action(Dg, Z, Y, [fmpq(0)]*(Kmom + 1))
    out = []
    for b in (0, 4):
        for a in (998, 1500, 3000, 10**4):
            J = Kmom + b
            Delta = J + a + 2; T = (Delta - 1)*(Delta - 2); C = J*(J + 1)
            Ha = float(A.H_ang(J, b, T).mid()); H2 = float(A2.H_ang(J, b, T).mid())
            hy = LS.hyb_fast(Dg, J, math.sqrt(T - C)/Lam, '++' if b == 0 else '+-')
            s = (f'| {b} | {a} | edge of reach | {J} | {Lof(J, a):.4f} | {Ha:.6e} | {H2:.6e} | {hy["base"]:.6e} | {hy["Ns"]:.6e} | '
                 f'{Ha/hy["Ns"]:.5f} | {H2/hy["Ns"]:.5f} |')
            print(s, flush=True); out.append(s)
    open(os.path.join(FC.CACHE, 'fp_stripscan_edge.txt'), 'w').write('\n'.join(out) + '\n')
