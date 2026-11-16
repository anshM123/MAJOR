"""
t3: AdS heavy action of M02, M01 and F_g3 near the gap, with the skeleton-calibrated remainder.

Per heavy state (tau >= Dg, J allowed, h in {++, +-}), normal form (ADS4_FUNCTIONAL.md, ROBUST_FUNCTIONAL.md section 5.2):
   F_hyb = F_flat(J, u) + [leading Regge - its flat limit at fixed bhat]   (term by term; robust/ads.py hybrid)
Remainder (transferred from the exact scalar skeleton, Df = 3, relative to the skeleton's own hybrid; remainder.py):
   E2 = int |psi2| (2 - x^2) (|rho_2| + err_2) [controlled nodes],   U2 = int |psi2| (2 - x^2) [uncontrolled nodes]
   E3 = u int |psi3| (|rho_3| + err_3),                              U3 = u int |psi3| [uncontrolled]
   EF = a u^3 (|R4| + e4 + (|R5| + e5)/2)  if the forward series is controlled, else UF = |FW_hyb| (counted as uncontrolled)
   x = p sqrt(u) = nu/m.
Spin-specific 1/m^2 terms (not in the skeleton; SPINNING_E2E.md section 6.3, l = 2 extra piece ~ 17 cosh(beta)/m^2 per spinning pair):
   Sp = Cs cosh(beta)/m^2 N_hyb,  Cs = 0 and 40.
Reported per (functional, Dg, h):
   min F_hyb/N;  min (F_hyb - E - Sp)/N on states with no uncontrolled part;  kappa_crit = min over states of (F_hyb - E - Sp)/U
   (positivity holds on every scanned state if the relative AdS correction on the uncontrolled nodes is below kappa_crit);
   C1_eff(flat) = max (E + EF)/(u N) over controlled states (the envelope constant this remainder corresponds to);
   and the fraction U/N at the gap.
Usage: python3 t3_scan.py Dg [Dg ...]      (writes cache/out_t3_<Dg>.txt)
"""
import sys, os, time, pickle
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', 'robust'))
import numpy as np
from rcommon import ReggeAction, SIGMA, flat_parts, gl01, peval
from check import load
from ads import flatT
import remainder as RM
import skelfast

FUNS = ['M02', 'M01', 'Fg3']
CS = (0.0, 40.0)
NP = 48
AMOD, PM = 20.0, 4.0
pn, pw = gl01(NP)

def mfl_grid():
    return np.unique(np.round(np.concatenate([np.arange(1.0, 1.3, 0.01), np.arange(1.3, 3.0, 0.05), np.geomspace(3.0, 10.0, 14)]), 6))

def jlist(h, mfl, bmax=40.0):
    jmax = int(min(400, (bmax*mfl - 1)/2))
    js = list(range(0, min(jmax, 40) + 1)); j = 41
    while j <= jmax:
        js.append(j); j = int(j*1.06) + 1
    js = [x for x in js if x <= jmax]
    if h == '++':
        js = sorted(set(x + (x % 2) for x in js if x + (x % 2) <= max(jmax, 0)))
    else:
        js = [x for x in js if x >= 4]
    return np.array(js, int)

def skeleton_state(Dg, Lam, mfl, J, RA0):
    """functional-independent remainder data for the heavy state with m = mfl Lam and spin J (None if unphysical)."""
    m = mfl*Lam
    Bq = 2*J + 1.0
    alpha = (-Bq + np.sqrt(Bq*Bq + 4*m*m))/2
    if alpha < Dg - 2 - 1e-9:
        return None
    abar = alpha + 2*J + 1
    m2 = alpha*abar
    eta = (alpha**2 + abar**2)/(2*alpha*abar)
    beta = np.arccosh(eta)
    Delta = alpha + J + 2
    u = 1/mfl**2
    x = pn*np.sqrt(u)
    nus = x*m
    out = {'beta': beta, 'm': m, 'eta': eta}
    Pd = None
    for k in (2, 3):
        rho, err, ctrl = RM.residual_fixedJ(k, x, J, m2)
        need = ~ctrl & (x <= 0.6)
        if need.any():
            if Pd is None:
                Pd = skelfast.Pders_vec(nus, eta)
            r2, e2, c2 = RM.residual_seriesC(k, nus, Delta, J, m2, eta, Pd)
            use = need & c2
            rho = np.where(use, r2, rho); err = np.where(use, e2, err); ctrl = ctrl | use
        out[k] = (rho, np.where(ctrl, err, 0.0), ctrl)
    _, _, W, P0 = RA0.transforms(1.0, beta)
    R4, e4, c4, R5, e5, c5, _, _ = RM.forward_remainders(Delta, J, m2, eta, W, P0)
    out['fw'] = (R4, e4, R5, e5, bool(c4 and c5))
    return out

def run(Dg):
    Lam = np.sqrt((Dg - 2.0)*(Dg - 1.0))
    Fs = {n: load(n) for n in FUNS}
    RAs = {n: ReggeAction(F) for n, F in Fs.items()}
    RA0 = RAs['M02']
    lines = []
    def pr(s):
        print(s, flush=True); lines.append(s)
    pr(f"# t3 Dg = {Dg}, Lam = {Lam:.4f}; states: m/Lam in [1, 10], bhat <= 40; NP = {NP} p-nodes")
    pk = os.path.join(HERE, 'cache', f'skel_t3_{Dg}.pkl')
    cache = pickle.load(open(pk, 'rb')) if os.path.exists(pk) else {}
    BANDS = ((1.0, 1.5), (1.5, 3.0), (3.0, 10.01))
    band = {(n, h, b): dict(EN=0.0, UN=0.0, low=np.inf, mod=np.inf, modloc=None, kap=np.inf) for n in FUNS for h in SIGMA for b in BANDS}
    nearg = {}
    res = {(n, h): dict(minFN=(np.inf, None), minLow={c: (np.inf, None) for c in CS}, kap={c: (np.inf, None) for c in CS},
                        C1=(0.0, None), Umax=(0.0, None), nst=0, nunc=0) for n in FUNS for h in SIGMA}
    t0 = time.time()
    for mfl in mfl_grid():
        u = 1/mfl**2
        m = mfl*Lam
        for h, s in SIGMA.items():
            js = jlist(h, mfl)
            if js.size == 0:
                continue
            flat = {n: flat_parts(F, h, u, int(js.max())) for n, F in Fs.items()}
            for J in js:
                key = (round(mfl, 6), int(J))
                if key not in cache:
                    cache[key] = skeleton_state(Dg, Lam, mfl, int(J), RA0)
                sk = cache[key]
                if sk is None:
                    continue
                beta = sk['beta']
                x = pn*np.sqrt(u)
                rho2, err2, c2 = sk[2]; rho3, err3, c3 = sk[3]
                R4, e4, R5, e5, cF = sk['fw']
                bh = (2*J + 1)/mfl
                for n, F in Fs.items():
                    Js, B2, B3, FW = flat[n]
                    i = int(np.searchsorted(Js, J))
                    I2, I3, W, P0 = RAs[n].transforms(Lam, beta)
                    f2 = flatT(RAs[n].c2, bh); f3 = flatT(RAs[n].c3, bh)
                    dB2 = 2*(I2 - f2); dB3 = s*u*(I3 - f3)
                    dFW = u**3*(2*F.a*m*m*W - 2*F.a*J*(J + 1) - F.a/2) + F.a*s*u**3*(P0 - 1)
                    b2, b3, fw = B2[i] + dB2, B3[i] + dB3, FW[i] + dFW
                    Fh = b2 + b3 + fw
                    N = abs(b2) + abs(b3) + abs(fw)
                    a2 = np.abs(peval(F.psi2, pn))*(2 - x*x)*pw
                    a3 = u*np.abs(peval(F.psi3, pn))*pw
                    E = np.sum(a2*(np.abs(rho2) + err2)*c2) + np.sum(a3*(np.abs(rho3) + err3)*c3)
                    U = np.sum(a2*(~c2)) + np.sum(a3*(~c3))
                    if cF:
                        EF = F.a*u**3*(abs(R4) + e4 + (abs(R5) + e5)/2); UF = 0.0
                    else:
                        EF = 0.0; UF = abs(fw)
                    # pole-model extrapolation on uncontrolled nodes [NE]: |rho| <= A max(1, (1-x^2)^-PM)/m^2, A = 10, PM = 4
                    rmod = AMOD*np.maximum(1.0, (1 - np.minimum(x*x, 0.999))**(-PM))/m**2
                    Um = np.sum(a2*rmod*(~c2)) + np.sum(a3*rmod*(~c3)) + (0.0 if cF else F.a*u**3*AMOD*2/m**2*max(1.0, abs(2*m*m*W)))
                    if mfl < 1.5:
                        key2 = (n, h, round(float(mfl), 3))
                        nearg[key2] = min(nearg.get(key2, (np.inf,))[0:1] + ((Fh - E - EF - Um)/N,)), 
                        nearg[key2] = (min(nearg[key2]), )
                    for b in BANDS:
                        if b[0] <= mfl < b[1]:
                            bb = band[(n, h, b)]
                            bb['EN'] = max(bb['EN'], (E + EF)/N); bb['UN'] = max(bb['UN'], (U + UF)/N)
                            lowm = (Fh - E - EF - Um)/N
                            if lowm < bb['mod']:
                                bb['mod'] = lowm; bb['modloc'] = (round(float(mfl), 3), int(J), round(float(beta), 4))
                            if U + UF == 0:
                                bb['low'] = min(bb['low'], (Fh - E - EF)/N)
                            else:
                                bb['kap'] = min(bb['kap'], (Fh - E - EF)/(U + UF))
                    r = res[(n, h)]
                    r['nst'] += 1
                    loc = (round(float(mfl), 3), int(J), round(float(beta), 4))
                    if Fh/N < r['minFN'][0]:
                        r['minFN'] = (Fh/N, loc)
                    Ut = U + UF
                    if Ut > 0:
                        r['nunc'] += 1
                        if Ut/N > r['Umax'][0]:
                            r['Umax'] = (Ut/N, loc)
                    else:
                        c1 = (E + EF)/(u*N)
                        if c1 > r['C1'][0]:
                            r['C1'] = (c1, loc)
                    for Cs in CS:
                        Sp = Cs*np.cosh(beta)/m**2*N
                        low = Fh - E - EF - Sp
                        if Ut == 0 and low/N < r['minLow'][Cs][0]:
                            r['minLow'][Cs] = (low/N, loc)
                        if Ut > 0:
                            kap = low/Ut
                            if kap < r['kap'][Cs][0]:
                                r['kap'][Cs] = (kap, loc)
                        elif low <= 0 and -np.inf < r['kap'][Cs][0]:
                            r['kap'][Cs] = (-np.inf, loc)
    for n in FUNS:
        for h in SIGMA:
            r = res[(n, h)]
            pr(f"  {n:4s} {h}  states {r['nst']:5d} (with uncontrolled nodes: {r['nunc']:5d})  min F_hyb/N = {r['minFN'][0]:+.4f} at {r['minFN'][1]}")
            for Cs in CS:
                pr(f"        Cs={Cs:4.0f}: min (F_hyb - E - Sp)/N [fully controlled states] = {r['minLow'][Cs][0]:+.4f} at {r['minLow'][Cs][1]};"
                   f"  kappa_crit = {r['kap'][Cs][0]:+.4g} at {r['kap'][Cs][1]}")
            pr(f"        C1_eff(flat-anchored) = max (E+EF)/(u N) = {r['C1'][0]:.4g} at {r['C1'][1]};  max U/N = {r['Umax'][0]:.4f} at {r['Umax'][1]}")
    pr("# per m/Lam band (Cs = 0): max (E+EF)/N, max U/N, min (F-E)/N [controlled], kappa_crit, min (F - E - U_model)/N [pole model A=%g, p=%g] at" % (AMOD, PM))
    for n in FUNS:
        for h in SIGMA:
            for b in BANDS:
                bb = band[(n, h, b)]
                pr(f"  {n:4s} {h} m/Lam in [{b[0]:.1f},{min(b[1],10):.1f}): E/N <= {bb['EN']:.3g}  U/N <= {bb['UN']:.3g}  min(F-E)/N = {bb['low']:+.4f}  kappa_crit = {bb['kap']:.3g}  model: {bb['mod']:+.4f} at {bb['modloc']}")
    pr("# near the gap: min over J of the pole-model margin (F - E - U_model)/N per m/Lam (M02, M01)")
    for n in ('M02', 'M01'):
        for h in SIGMA:
            ms = sorted(set(k[2] for k in nearg if k[0] == n and k[1] == h))
            pr(f"  {n} {h}: " + " ".join(f"{mm:.2f}:{nearg[(n, h, mm)][0]:+.3f}" for mm in ms[:16]))
    pr(f"# time {time.time()-t0:.0f} s, {len(cache)} skeleton states")
    open(os.path.join(HERE, 'cache', f'out_t3_{Dg}.txt'), 'w').write("\n".join(lines) + "\n")
    pickle.dump(cache, open(os.path.join(HERE, 'cache', f'skel_t3_{Dg}.pkl'), 'wb'))

if __name__ == '__main__':
    for Dg in [int(a) for a in sys.argv[1:]] or [100]:
        run(Dg)
