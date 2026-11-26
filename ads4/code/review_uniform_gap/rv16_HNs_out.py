"""rv16_HNs_out.py -- normalised margins H/N_s (new family, ideal) in the region the note does NOT cover (J>=65, hat-b>40),
at Delta_g = 1000 (N_s from xcheck_ladder/xl_norm.hybrid, read-only), for planning."""
import sys, os, math
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'xcheck_ladder')))
import mpmath as mp
import xl_norm
from flint import arb, fmpq, ctx
from rv7_sweep import AGtab, S_value
Dg = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
lam2 = (Dg-1)*(Dg-2); ag = Dg-2; t = fmpq(1, lam2)
res = []
for J in [65, 100, 200, 300, 1000]:
    medge2 = ag*(ag+2*J+1)
    bedge = (2*J+1)*math.sqrt(lam2/medge2)
    bl = [x for x in (41, 60, 100, 200, 400) if x < bedge] + [bedge]
    prec = 256 + int(3*bedge)
    tab = AGtab(J+1, [t], prec)
    for bh in bl:
        u2 = fmpq(int(round(bh**2*10**8)), 10**8*(2*J+1)**2) if bh != bedge else fmpq(lam2, medge2)
        u = arb(u2)
        for b in (0, 4):
            if b == 0 and J % 2:
                continue
            S = S_value(J, b, t, u, tab)
            C = J*(J+1); m2 = lam2/float(u2.p)*float(u2.q) if False else lam2*float(u2.q)/float(u2.p)
            w = m2/(m2+C)
            H = w**3*float(S.mid())
            mfl = mp.nstr(mp.sqrt(mp.mpf(int(u2.q))/mp.mpf(int(u2.p))), 25)
            hy = xl_norm.hybrid(Dg, J, mfl, '++' if b == 0 else '+-', npn=max(200, int(bedge/2)), dps=30)
            res.append((J, b, round(bh, 2), H/float(hy['Ns']), float(hy['Hhyb']/hy['Ns'])))
            print(f"J={J} b={b} hat-b={bh:.2f}{' (twist edge)' if bh == bedge else ''}: H/Ns = {res[-1][3]:.6f}   Hhyb/Ns = {res[-1][4]:.6f}", flush=True)
print("MIN:", min(res, key=lambda r: r[3]))
