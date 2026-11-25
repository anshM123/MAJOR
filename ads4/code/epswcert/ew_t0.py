"""
epswcert/ew_t0.py (EPSW_CERT.md §4): T0's realised-weight corollary at Delta_g = 300, 1000 from the ball-arithmetic certificate of
ew_cert.py, for the PUBLISHED k = 2 weight ('pub') and for the minimally repaired one ('rep', ew_repair.py).

Convention (T0_CERT §4.2): the realised k = 2 weight changes only the leading-Regge image carried by the hybrid's "2"-curvature:
X -> X + Delta, Delta = 2 dI(beta)/Lam (p-units), dI(beta) = int_0^oo r(mu) P_{-1/2+i mu}(cosh beta) dmu = (sqrt2/pi) int_0^beta d(t) K(t) dt,
d = W - Phi, K(t) = (cosh beta - cosh t)^(-1/2), K0 = K(0) = 1/(sqrt2 sinh(beta/2)).  For pattern g, Lam^2 * (realised curvature part)
= Lam^2 curv_g + (1 - c s2) 2 Lam dI.
Bounds on d (ew_cert):  |d(t)| <= eps(t) := min(eps_W, B/t^2).
  (U) |dI(beta)| <= eps_W P_{-1/2}(cosh beta) <= eps_W min(1, beta/(2 sinh(beta/2)))      [int_0^beta K <= (pi/2) beta K0]
  (A) dI(beta) = r(0) K0/sqrt2 + R(beta)   [int_0^oo d = (pi/2) r(0): Fourier inversion, r_e continuous, r_e, d in L^1]
      |R| <= (sqrt2/pi) K0 rhoR(beta),  rhoR = B/beta + 8 B ln(4 cosh(beta/4))/beta^2 + 0.51 (eps_W sinh 2 + 2 B cosh(beta/2)/beta)/(cosh beta - 1)
      (beta >= 64; pieces: K0 int_beta^oo |d|; t in [beta/2, beta] with int_{beta/2}^beta K <= 2 K0 ln(4 cosh(beta/4)); t <= beta/2 with
       K - K0 <= 0.51 K0 (cosh t - 1)/(cosh beta - 1)).
Target floor for beta >= 64 (T0 Part L tail, re-evaluated here with the cited t0cert/gaplarge code at Lam_a = Lam(Dg)):
  Lam^2 curv_g >= (S_g - (1+c) tlA) M1 - (1+c) FL(beta),  M1 >= (sqrt2/pi) K0 L(beta),  beta L(beta) increasing (LRINF_TIGHT §1.5).
Usage: python3 ew_t0.py     -> cache/ew_t0.txt, cache/ew_t0.json
"""
import os, sys, json, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
CODE = os.path.abspath(os.path.join(HERE, '..'))
from fractions import Fraction as Fr
import importlib

LINES = []


def pr(s=''):
    print(s, flush=True); LINES.append(s)


def load_T(Dg, c):
    os.environ['T0_DGMIN'] = str(Dg); os.environ['T0_C'] = str(c)
    sys.path.insert(0, os.path.join(CODE, 't0cert'))
    for m in ('t0_common',):
        if m in sys.modules:
            del sys.modules[m]
    import t0_common as T
    return T


def run(Dg, c, wtags):
    T = load_T(Dg, c)
    G = T.G
    from flint import arb
    A, up, lo = T.A, T.up, T.lo
    PI = arb.pi(); S2 = arb(2).sqrt()
    La = lo(T.LAMMIN)                       # lower end of Lam(Dg): all T0 remainder bounds hold for every Lam >= La
    Lam = T.LAMMIN
    cA = A(c)
    t1tag = f'{Dg}_{str(c).replace("/", "_")}'
    j1 = json.load(open(os.path.join(CODE, 't0cert', 'cache', f't1_curv_{t1tag}.json')))
    t3p = os.path.join(HERE, 'cache', f't3_flat_{t1tag}.json')
    j3 = json.load(open(t3p if os.path.exists(t3p) else os.path.join(CODE, 't0cert', 'cache', f't3_flat_{t1tag}.json')))
    src3 = 'epswcert re-run (copy of t3_flat.py)' if os.path.exists(t3p) else 't0cert'
    assert j3['pass_F'] in (True, 'True') and j1['pass_K'] in (True, 'True')
    tx = j1['Tabs_upper'].strip('[]').split('+/-')
    Tabs = arb(tx[0].strip()) + arb(tx[1].strip()) + arb('1e-15')
    shift = Fr(j3['Tabs']); shiftA = A(shift)
    Kmin = arb(j1['Kmin']) - arb('1e-12'); Ktail = arb(j1['Ktail']) - arb('1e-12')
    # ---------------- Part L cells and tail ingredients (as t2_lrinf)
    TAB = {k: G.build_I(G.W[k], La, hat=True) for k in (2, 3)}
    a12, a13 = A(T.A1[2]), A(T.A1[3])
    Smin = (1 - cA)*2*a12 - (1 + cA)*a13
    FE = {k: T.flat_E_consts(k) for k in (2, 3)}

    def e_of(k, b):
        d = FE[k]
        return A(d['f1'])*(1 + 2/b) + A(d['g1']) + (A(d['dg1']) + d['Iq'])/b

    bc = 2*arb(1).asinh()
    cells = []
    for a, b, m1 in T.CELLS:
        if not A(b) > bc:
            continue
        e2, _ = G.E_up(2, a, b, TAB, La); e3, _ = G.E_up(3, a, b, TAB, La)
        L2N = up(La**2*(2*e2 + e3))
        sa = 2*(A(a)/2).sinh()
        fl = up((2*e_of(2, La*sa) + e_of(3, La*sa))/(La*sa**3))
        K = Smin*m1 - (1 + cA)*L2N - (1 + cA)*fl
        Pbar = A(a)/(2*(A(a)/2).sinh())
        if Pbar > 1:
            Pbar = arb(1)
        cells.append((a, b, K, Pbar))
    tl, ell = G.tail_majorant(TAB, La)
    tlA = tl*A(G.A1MIN)
    pats = []
    for s2 in (1, -1):
        for s3 in (1, -1):
            for sg in (1, -1):
                for u in (0, 1):
                    pats.append((s2, s3, sg, u, (1 - cA*s2)*2*a12 + (1 - cA*s3)*u*sg*a13, 1 - cA*s2))

    def betaL(B):
        ln2 = arb(2).log()
        return 1 - B/(2*(B/2).sinh()) - S2*(B*(-B/2).exp()/A(Fr(47, 10)) + 4/B) - 2*(2*ln2).sqrt()*B/(B - ln2)**2

    def betaFL(B):          # beta * FL(beta) * pi sinh(beta/2)   (FL relative to (sqrt2/pi) K0, times beta); decreasing for beta >= 64
        s = 2*(B/2).sinh()
        return B*(2*e_of(2, La*s) + e_of(3, La*s))/(La*s**3)*PI*(B/2).sinh()

    def betarhoR(B, eps, Bc):
        return Bc + 8*Bc*(4*(B/4).cosh()).log()/B + B*arb('0.51')*(eps*arb(2).sinh() + 2*Bc*(B/2).cosh()/B)/(B.cosh() - 1)

    res = {}
    for wt in wtags:
        ew = json.load(open(os.path.join(HERE, 'cache', f'ew_cert_{Dg}_re' + ('' if wt == 'pub' else '_' + wt) + '.json')))
        eps = arb(ew['epsW'])*(1 + arb('1e-12')); Bc = arb(ew['B'])*(1 + arb('1e-12'))      # floats of upper bounds: widen by 1e-12 rel.
        r0lo = arb(ew['r0_lo']) - arb('1e-18'); r0hi = arb(ew['r0_hi']) + arb('1e-18')
        LamU = arb(Lam.upper())
        pr(f'--- Dg = {Dg}, c = {c}, weight = {wt}: eps_W <= {ew["epsW"]:.4e}, B <= {ew["B"]:.4e}, r(0) in [{ew["r0_lo"]:.6e}, {ew["r0_hi"]:.6e}]')
        # Part F/T (bhat <= 16): slack of the Part-F shift over T_abs must absorb (1+c) 2 eps_W / Lam
        need = (1 + cA)*2*eps/La
        slack = shiftA - Tabs
        okF = bool(slack > need)
        pr(f'  Part F/T (bhat <= 16): shift used ({src3}) - T_abs >= {slack.str(4)}  vs  (1+c) 2 eps_W/Lam <= {float(up(need)):.3e}: {"PASS" if okF else "FAIL"}')
        # Part K (16 <= bhat <= 2 Lam)
        needK = (1 + cA)*2*LamU*eps
        okK = bool(Kmin > needK and Ktail > needK)
        pr(f'  Part K (16 <= bhat <= 2Lam): min K_g >= {Kmin.str(5)} (tail {Ktail.str(5)})  vs  (1+c) 2 Lam eps_W <= {float(up(needK)):.3e}: '
           f'{"PASS" if okK else "FAIL"}')
        # Part L cells (2 asinh 1 < beta <= 64)
        okL = True; worst = None
        for a, b, K, Pbar in cells:
            KR = K - (1 + cA)*2*LamU*eps*Pbar
            if not KR > 0:
                okL = False
            q = KR/(K)
            if worst is None or q.lower() < worst[0]:
                worst = (float(q.lower()), float(a))
        pr(f'  Part L cells (beta in (1.7627, 64], {len(cells)} cells): realised K_g > 0 on all: {okL}; min (realised K)/(target K) >= '
           f'{worst[0]:.6f} (cell at beta = {worst[1]})')
        # beta >= 64
        rminus = arb(0) if r0lo > 0 else -r0lo
        B64 = arb(64)
        base = {}
        for (s2, s3, sg, u, Sg, w2) in pats:
            Sp = Sg - (1 + cA)*tlA
            assert Sp > 0
            base[(s2, s3, sg, u)] = (Sp, w2)
        if not rminus > 0:
            okT = True; mn = None
            for key, (Sp, w2) in base.items():
                v = Sp*betaL(B64) - (1 + cA)*betaFL(B64) - w2*2*LamU*betarhoR(B64, eps, Bc)
                okT &= bool(v > 0)
                mn = v if mn is None or v < mn else mn
            pr(f'  beta >= 64 (r(0) >= 0 certified): min over 16 patterns of beta*G_g at beta = 64 (bounds every beta >= 64 by '
               f'monotonicity) >= {mn.str(5)} > 0: {okT}  -> realised corollary for ALL beta')
            beta1 = None
        else:
            # cells [ba, bb] geometric from 64: lower bound of beta*G_g valid on the cell
            ba = arb(64); okT = True; ncell = 0
            while True:
                bb = ba*arb('1.0005')
                ok = True
                for key, (Sp, w2) in base.items():
                    v = Sp*betaL(ba) - (1 + cA)*betaFL(arb(64)) - w2*2*LamU*((PI/2)*rminus*bb + betarhoR(ba, eps, Bc))
                    if not v > 0:
                        ok = False; break
                if not ok:
                    break
                ba = bb; ncell += 1
            beta1 = float(ba.lower())
            pr(f'  beta >= 64, r(0) < 0: realised K_g > 0 certified on [64, {beta1:.1f}] ({ncell} cells); fails beyond (see next line)')
            # certified failure on ++ states: upper bound of Lam^2 H^R/((sqrt2/pi)K0) at beta2, decreasing in beta except -pi Lam |r(0)|
            res.setdefault('fail', {})
            beta2 = failure(T, Dg, c, La, LamU, a12, a13, tlA, betaFL, betarhoR, eps, Bc, r0hi, beta1)
            res['fail'][wt] = beta2
        res[wt] = dict(F=okF, K=okK, L=okL, tail=okT, beta1=beta1)
    return res


def failure(T, Dg, c, La, LamU, a12, a13, tlA, betaFL, betarhoR, eps, Bc, r0hi, beta1):
    """Certified: for every ++ heavy state (even J, any u in (0,1]) with beta >= beta2, the realised hybrid H^R < 0 (so H^R < c N_s).
    Lam^2 H^R/((sqrt2/pi)K0) <= flat + (2a12 + a13 + tlA)[1/beta + (0.51/24) sinh(beta/2)/(cosh beta - 1)] + FLterm
                                 + 2 Lam [(pi/2) r(0) + rhoR(beta)],    r(0) <= r0hi < 0,
    flat = Lam^2 (T2 + T3)/(bhat^2 - 1) * pi sinh(beta/2), bhat = 2 Lam sinh(beta/2) (Legendre integration by parts, |P_J| <= 1:
    |B_2| + |B_3| <= (T2 + T3)/(4 u J(J+1)) = (T2 + T3)/(bhat^2 - u), T_k = sum of |coefficients| of the u-components of
    p(1 - u p^2) d/dp[(2 - u p^2) psi2/p] (k = 2) and p(1 - u p^2) d/dp[psi3/p] (k = 3)).  Every term except the r(0) one decreases
    in beta for beta >= 64, so the bound at beta2 holds for all beta >= beta2."""
    from flint import arb
    A = T.A
    PI = arb.pi()
    # T2, T3 from exact polynomials (lowest power first)
    def pmul(a, b):
        out = [Fr(0)]*(len(a) + len(b) - 1)
        for i, x in enumerate(a):
            for j, y in enumerate(b):
                out[i + j] += x*y
        return out
    def pder(a):
        return [i*a[i] for i in range(1, len(a))] or [Fr(0)]
    tot = Fr(0)
    for k in (2, 3):
        psi = T.PSI[k]
        assert psi[0] == 0
        pt = psi[1:]                                         # psi/p
        # (2 - u p^2) psi/p for k = 2 ; psi/p for k = 3 : components in powers of u
        comps = [[Fr(2)*x for x in pt], [Fr(0), Fr(0)] + [-x for x in pt]] if k == 2 else [pt]
        Gs = [pder(g) for g in comps]                          # d/dp of each u^j component
        # p (1 - u p^2) * G'  -> u-components
        hs = {}
        for j, g in enumerate(Gs):
            h0 = pmul([Fr(0), Fr(1)], g)
            h1 = pmul([Fr(0), Fr(0), Fr(0), Fr(-1)], g)
            for jj, hh in ((j, h0), (j + 1, h1)):
                cur = hs.get(jj) or []
                if len(cur) < len(hh):
                    cur = cur + [Fr(0)]*(len(hh) - len(cur))
                hs[jj] = [cur[i] + (hh[i] if i < len(hh) else 0) for i in range(len(cur))]
        tk = sum((sum(abs(x) for x in h) for h in hs.values() if h), Fr(0))
        tot += tk
    TT = A(tot)

    def U(B):
        s = (B/2).sinh()
        bh = 2*La*s
        flat = LamU**2*TT/(bh*bh - 1)*PI*s
        curv = (2*a12 + a13 + tlA)*(1/B + arb('0.51')/24*s/(B.cosh() - 1)) + (1 + 0)*betaFL(B)/B
        return flat + curv + 2*LamU*((PI/2)*r0hi + betarhoR(B, eps, Bc)/B)
    # smallest grid beta (step 0.1 %) with U < 0
    B = arb(int(beta1))
    while not U(B) < 0:
        B = B*arb('1.001')
    beta2 = float(B.upper())
    pr(f'  CERTIFIED FAILURE (pub): for every ++ state (even J, m >= Lam) with beta >= {beta2:.1f}: Lam^2 H^R/((sqrt2/pi)K0) <= '
       f'{U(B).str(4)} < 0, i.e. H_realised < 0 <= c N_s   (sum of |coeff| T2 + T3 = {float(tot):.4e})')
    return beta2


def main():
    t0 = time.time()
    pr(__doc__.split('Usage')[0])
    out = {}
    out['300_0.199'] = run(300, Fr(199, 1000), ['pub', 'rep'])
    out['1000_0.199'] = run(1000, Fr(199, 1000), ['pub', 'rep'])
    out['1000_0.1999'] = run(1000, Fr(1999, 10000), ['pub', 'rep'])
    pr(f'\n[{time.time() - t0:.0f}s]')
    json.dump(out, open(os.path.join(HERE, 'cache', 'ew_t0.json'), 'w'), indent=1, default=str)
    open(os.path.join(HERE, 'cache', 'ew_t0.txt'), 'w').write('\n'.join(LINES) + '\n')


if __name__ == '__main__':
    main()
