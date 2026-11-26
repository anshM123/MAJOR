"""rv7_sweep.py -- check (3b): dense numerical sweep (NOT a proof) of the sign of the note's closed-form ideal angular action
S = H/w^3 (validated against Eq. fp-H in rv6), in Arb ball arithmetic (each point's sign is rigorous for that point; the sweep is not).
Regions:  'low'   : all physical J<=64 (b=0 even J, b=4 J>=4), masses from the twist edge (zeta = m_edge^2/m^2 = 1) to m = infinity
          'imp'   : J in {65,100,300,1000}, hat-b <= 40 (physical part)
          'out'   : J in {65,...,300}, 40 < hat-b <= 400 (physical part; region claimed NOT covered)
Gaps: Delta_g in {1000,1500,3000,1e4,1e5,inf}; Lambda^2 = (Dg-1)(Dg-2); twist edge m^2 = a_g(a_g+2J+1), a_g = Dg-2.
Usage: python3 rv7_sweep.py low|imp|out"""
import sys, json, math, time
from flint import arb, fmpq, fmpq_poly, ctx
from fractions import Fraction as Fr
from rv_common import PSI2, PSI3, mom

GAPS = [1000, 1001, 1100, 1500, 2000, 3000, 5000, 10**4, 3*10**4, 10**5, 10**6, None]   # None = infinity (t = 0)
ctx.prec = 64

def gpolys_flint(p, N):
    Y = fmpq_poly([0, 1])
    out = [fmpq_poly([1])]; prev = fmpq_poly([0])
    for k in range(N):
        cur = out[-1]
        nxt = ((8*k*k+8*p*k+4*p-1) - 4*Y)*cur - 4*(k-1+p)**2*prev
        nxt = nxt*fmpq(1, 4*(k+1)**2)
        prev = cur; out.append(nxt)
    return out

class AGtab:
    """A_k(t), G_k(t) for k<=N for a list of t (fmpq or 0), as arb at precision prec."""
    def __init__(self, N, ts, prec):
        ctx.prec = prec
        self.N = N
        g3 = gpolys_flint(3, N+1); g4 = gpolys_flint(4, N+1)
        m2 = [arb(fmpq(x.numerator, x.denominator)) for x in (mom(PSI2, j) for j in range(2*N+4))]
        m3 = [arb(fmpq(x.numerator, x.denominator)) for x in (mom(PSI3, j) for j in range(2*N+4))]
        self.A = {}; self.G = {}
        for t in ts:
            self.A[t] = []; self.G[t] = []
        for k in range(N+2):
            c3 = g3[k].coeffs(); c4 = g4[k].coeffs()
            c3 = [arb(c) for c in c3] + [arb(0)]*(k+1-len(c3)); c4 = [arb(c) for c in c4] + [arb(0)]*(k+1-len(c4))
            for t in ts:
                if t == 0:
                    self.A[t].append(2*c3[k]*m2[2*k]); self.G[t].append(c4[k]*m3[2*k])
                else:
                    ta = arb(t)
                    # sum_j c_{k,k-j} t^j mu[2(k-j)]
                    sA = arb(0); sG = arb(0); tp = arb(1)
                    for j in range(k+1):
                        sA += c3[k-j]*tp*m2[2*(k-j)]; sG += c4[k-j]*tp*m3[2*(k-j)]
                        tp *= ta
                    self.A[t].append(2*sA); self.G[t].append(sG)

def S_value(J, b, t, u, tab, kmax=None):
    """closed form S at t (fmpq or 0), u (arb)"""
    n = J-b; C = J*(J+1); sig = 1 if b == 0 else -1
    A = tab.A[t]; G = tab.G[t]; ta = arb(t) if t != 0 else arb(0)
    fac = u/(1+C*ta*u)
    tot = arb(0); pref = arb(1)
    top = n if kmax is None else min(n, kmax)
    for k in range(top+1):
        gam = arb(fmpq(-k*k+(8-4*b)*k+9-2*b*(b-1), 2))
        tot += pref*(A[k] + fac*(ta*gam*A[k] + arb(fmpq((k+1)**2, 2))*A[k+1] + sig*G[k]))
        Bk = C-(b+k)*(b+k+1)
        pref = pref*u*Bk/(1+Bk*ta*u)
    return tot

def gap_data(Dg):
    if Dg is None:
        return fmpq(0), None
    lam2 = (Dg-1)*(Dg-2)
    return fmpq(1, lam2), lam2

def u_of(J, Dg, zeta):
    """u = L^2/m^2 with m^2 = m_edge^2/zeta; at Dg = inf, u = zeta (u_edge -> 1)."""
    if Dg is None:
        return arb(zeta)
    ag = Dg-2; lam2 = (Dg-1)*(Dg-2)
    medge2 = ag*(ag+2*J+1)
    return arb(fmpq(lam2, medge2))*arb(zeta)

def bhat(J, Dg, u):
    # hat-b = (2J+1) L/m = (2J+1) sqrt(u)
    return (2*J+1)*u.sqrt()

def run_low():
    zetas = [fmpq(1)] + [fmpq(10**6-i, 10**6) for i in (1, 10, 100, 300, 1000, 3000)] + [fmpq(400-i, 400) for i in range(1, 400)] + \
            [fmpq(1, 10**j) for j in (3, 4, 6)] + [fmpq(0)]
    prec = 900
    ts = [gap_data(Dg)[0] for Dg in GAPS]
    t0 = time.time(); tab = AGtab(66, ts, prec); print(f"A,G table built {time.time()-t0:.1f}s", flush=True)
    res = []; glob = None
    for Dg in GAPS:
        t, lam2 = gap_data(Dg)
        for J in range(0, 65):
            for b in (0, 4):
                if (b == 0 and J % 2) or J < b:
                    continue
                vals = []
                for z in zetas:
                    u = u_of(J, Dg, z)
                    S = S_value(J, b, t, u, tab)
                    vals.append((S, z))
                # golden refinement around the smallest grid value (zeta in (0,1])
                Smin, zmin = min(vals, key=lambda v: float(v[0].mid()))
                i = [v[1] for v in vals].index(zmin)
                if 0 < i < len(vals)-1:
                    lo, hi = float(vals[i+1][1]), float(vals[i-1][1])
                    if lo > hi: lo, hi = hi, lo
                    f = lambda zz: float(S_value(J, b, t, u_of(J, Dg, fmpq(int(zz*10**12), 10**12)), tab).mid())
                    gr = (math.sqrt(5)-1)/2; a, c = lo, hi
                    x1, x2 = c-gr*(c-a), a+gr*(c-a); f1, f2 = f(x1), f(x2)
                    for _ in range(40):
                        if f1 < f2: c, x2, f2 = x2, x1, f1; x1 = c-gr*(c-a); f1 = f(x1)
                        else: a, x1, f1 = x1, x2, f2; x2 = a+gr*(c-a); f2 = f(x2)
                    zz = fmpq(int(x1*10**12), 10**12); Sz = S_value(J, b, t, u_of(J, Dg, zz), tab)
                    if float(Sz.mid()) < float(Smin.mid()):
                        Smin, zmin = Sz, zz
                rad_ok = all(float(v[0].rad()) < 1e-30*abs(float(v[0].mid())) for v in vals)
                neg = [float(v[1]) for v in vals if not v[0] > 0]
                row = dict(Dg=Dg if Dg else 'inf', J=J, b=b, Smin=float(Smin.mid()), zeta_at_min=float(zmin),
                           S_edge=float(vals[0][0].mid()), all_positive=(len(neg) == 0 and Smin > 0), nonpos_zetas=neg, tight=rad_ok)
                res.append(row)
                if glob is None or row['Smin'] < glob['Smin']:
                    glob = row
        print(f"Dg={Dg}: done {time.time()-t0:.1f}s; running min {glob}", flush=True)
    json.dump(dict(rows=res, global_min=glob), open('out/rv7_low.json', 'w'), indent=1)
    print("GLOBAL MIN (J<=64, all gaps, all masses scanned):", glob)
    print("any non-positive:", [r for r in res if not r['all_positive']])

def run_imp(Js, bmin, bmax, nb, tag):
    res = []; glob = None; t0 = time.time()
    for J in Js:
        prec = 256 + int(3*bmax)
        kmax = None
        ts = [gap_data(Dg)[0] for Dg in GAPS]
        # terms ~ x^k/(k!)^2: beyond k ~ e*sqrt(x)+60 they are < 2^-prec relative; but keep the full ladder up to Kc
        Kc = min(J, int(2.8*bmax/2) + 120)
        tab = AGtab(Kc+1, ts, prec)
        for Dg in GAPS:
            t, lam2 = gap_data(Dg)
            for b in (0, 4):
                if b == 0 and J % 2:
                    continue
                for i in range(nb+1):
                    bh = bmin + (bmax-bmin)*Fr(i, nb)
                    if bh == 0:
                        continue
                    # u = bh^2/(2J+1)^2 ; physical iff u <= u_edge
                    u = arb(fmpq((bh**2).numerator, (bh**2).denominator))/(2*J+1)**2
                    uedge = u_of(J, Dg, fmpq(1))
                    if not (u <= uedge):
                        continue
                    S = S_value(J, b, t, u, tab, kmax=Kc)
                    # truncation check: compare with Kc-20
                    S2 = S_value(J, b, t, u, tab, kmax=Kc-20)
                    row = dict(Dg=Dg if Dg else 'inf', J=J, b=b, bhat=float(bh), S=float(S.mid()), pos=bool(S > 0),
                               trunc_rel=float(abs((S-S2).mid()))/max(abs(float(S.mid())), 1e-300), rad=float(S.rad()))
                    res.append(row)
                    if glob is None or row['S'] < glob['S']:
                        glob = row
        print(f"J={J} done {time.time()-t0:.1f}s; running min {glob}", flush=True)
    json.dump(dict(rows=res, global_min=glob), open(f'out/rv7_{tag}.json', 'w'), indent=1)
    print("GLOBAL MIN:", glob)
    bad = [r for r in res if not r['pos']]
    print("non-positive points:", len(bad), bad[:20])
    print("max truncation rel:", max(r['trunc_rel'] for r in res), " max ball radius:", max(r['rad'] for r in res))

if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'low':
        run_low()
    elif mode == 'imp':
        run_imp([65, 100, 300, 1000], Fr(0), Fr(40), 160, 'imp')
    elif mode == 'out':
        run_imp([65, 80, 100, 130, 160, 200, 250, 300], Fr(40), Fr(400), 360, 'out')
