"""
poserr/large_m_map.py (POS_ERRORS.md §2): per-state budget of Omega^EE (Result D) at large mass, on the hybrid states of
spinfunc/ee1000_common.hyb (flat NF02 with exact Wigner functions + exact leading-Regge curvature; read-only import).
For each (Delta_g, m/Lam, J, helicity):
  lam_h   = (B2 + B3 + cv2 + cv3)/Ns          hybrid value (exact functions, quadrature)            [NE]
  flat_h  = (B2 + B3)/N_flat, N_flat = |B2| + |B3| (certified >= 1/5 for every flat state, CERT_CLOSE §1)
  cvb_h   = Mehler-Dirichlet bound on |cv2| + |cv3| (POS_ERRORS Lemma C):  |P_{-1/2+i nu}(cosh b) - J0(2 nu sinh(b/2))|
            <= (1 - sech(b/2)) + min(2, nu (2 sinh(b/2) - b))                                        [P]
  W_h     = sum (|psi2| env(h,x,2) + u |psi3| env(h,x,3)) / (m^2 Ns_h)   : cost in lam_h of a uniform sign-free |Gamma| = 1 c2
            on every nu (split nu < 12 / nu >= 12)
  W_off   = the same for the off-diagonal entry (ee1000_common.dF_off structure, normalised by sqrt(Ns++ Ns+-))
Output: cache/large_m_<Dg>.json and a text summary cache/large_m_<Dg>.txt.
Usage: python3 large_m_map.py Dg [mfl,mfl,...]
"""
import os, sys, json, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'spinfunc'))
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np
import ee1000_common as E
Q = E.Q
CA = os.path.join(HERE, 'cache')
NUSPLIT = 12.0
NUHI = 20.0


def jmax(Dg, mfl):
    """largest J with m = mfl*Lam twist-allowed: m^2 >= (Dg-2)(Dg-1+2J)."""
    Lm = E.lam_of(Dg); m2 = (mfl*Lm)**2
    return int(np.floor((m2/(Dg - 2.0) - (Dg - 1.0))/2))


def jgrid(Dg, mfl, bmax):
    Jm = min(jmax(Dg, mfl), int(np.floor(bmax*mfl/2)))
    if Jm < 0:
        return []
    Js = set(range(0, min(Jm, 60) + 1))
    b = 2.0
    while True:
        J = int(round(b*mfl/2 - 0.5))
        if J > Jm:
            break
        Js.add(J); b *= 1.12
    Js.add(Jm)
    return sorted(J for J in Js if not (J % 2 == 1 and J < 4) and J >= 0)


def npn_of(J, mfl, Lm):
    beta = 2*np.arcsinh((J + .5)/(mfl*Lm))
    ph = max(Lm*beta, (2*J + 1)*np.arcsin(min(1.0, 1/mfl)))
    return int(max(256, min(4096, 64*np.ceil((1.6*ph + 200)/64))))


def mdbound(nu, beta):
    b = 2*np.sinh(beta/2)
    return (1 - 1/np.cosh(beta/2)) + np.minimum(2.0, nu*(b - beta))


def state(Dg, J, mfl):
    Lm = E.lam_of(Dg)
    hs = ['++'] if J < 4 else ['++', '+-']
    npn = npn_of(J, mfl, Lm)
    out = dict(J=J, mfl=mfl, bhat=(2*J + 1)/mfl, npn=npn)
    sts = {}
    for h in hs:
        st = E.hyb(Dg, J, mfl, h, npn); sts[h] = st
        x, m, u = st['x'], st['m'], st['u']
        e2 = Q.A.env(h, x, 2); e3 = Q.A.env(h, x, 3)
        nu = st['pn']*Lm
        beta = 2*np.arcsinh((J + .5)/m)
        dP = mdbound(nu, beta)
        cvb = (np.abs(st['p2'])*2*dP).sum() + (np.abs(st['p3'])*u*dP).sum()
        w2v = np.abs(st['p2'])*e2/m**2/st['Ns']; w3v = u*np.abs(st['p3'])*e3/m**2/st['Ns']; wv = w2v + w3v
        Nf = abs(st['B2']) + abs(st['B3'])
        out[h] = dict(lam=float(st['base']/st['Ns']), Ns=float(st['Ns']), Nf=float(Nf), flat=float((st['B2'] + st['B3'])/Nf),
                      cv=float(st['cv2'] + st['cv3']), cvabs=float(abs(st['cv2']) + abs(st['cv3'])), cvb=float(cvb),
                      W_lo=float(wv[nu < NUSPLIT].sum()), W_hi=float(wv[nu >= NUSPLIT].sum()),
                      W2_lo=float(w2v[nu < NUSPLIT].sum()), W2_hi=float(w2v[nu >= NUSPLIT].sum()),
                      W_mid=float(wv[(nu >= NUSPLIT) & (nu < NUHI)].sum()), W_top=float(wv[nu >= NUHI].sum()),
                      W3_lo=float(w3v[nu < NUSPLIT].sum()), W3_hi=float(w3v[nu >= NUSPLIT].sum()),
                      rig=float((Nf/5 - cvb)/st['Ns']))
    if len(hs) == 2:
        sp, sm = sts['++'], sts['+-']
        x, m, y, u = sp['x'], sp['m'], sp['y'], sp['u']
        eo2 = np.sqrt(Q.A.env('++', x, 2)*Q.A.env('+-', x, 2))
        wo2 = np.abs(sp['p2'])*eo2/m**2/np.sqrt(sp['Ns']*sm['Ns']); wo3 = np.abs(sp['p3'])*u*eo2/(2 - y)/m**2/np.sqrt(sp['Ns']*sm['Ns'])
        wo = wo2 + wo3
        nu = sp['pn']*Lm
        out['W_off_lo'] = float(wo[nu < NUSPLIT].sum()); out['W_off_hi'] = float(wo[nu >= NUSPLIT].sum())
        out['W3_off'] = float(wo3.sum())
        out['W_off_mid'] = float(wo[(nu >= NUSPLIT) & (nu < NUHI)].sum()); out['W_off_top'] = float(wo[nu >= NUHI].sum())
    else:
        out['W_off_lo'] = out['W_off_hi'] = out['W3_off'] = out['W_off_mid'] = out['W_off_top'] = 0.0
    return out


def run(Dg, mfls, bmax=400.0):
    t0 = time.time()
    res = dict(Dg=Dg, Lam=float(E.lam_of(Dg)), states=[])
    for mfl in mfls:
        Js = jgrid(Dg, mfl, bmax)
        for J in Js:
            res['states'].append(state(Dg, J, mfl))
        print('Dg %g mfl %g: %d states, J <= %d, %.0fs' % (Dg, mfl, len(Js), Js[-1] if Js else -1, time.time() - t0), flush=True)
        json.dump(res, open(os.path.join(CA, "large_m_%g%s.json" % (Dg, os.environ.get("LM_TAG", ""))), "w"))
    return res


def run_gap(Dg, Jmax=40):
    """near-gap strip: J = 0..Jmax at m/Lam in {m_gap(J), 1.02, 1.05, 1.1, 1.2} (those >= m_gap)."""
    t0 = time.time()
    res = dict(Dg=Dg, Lam=float(E.lam_of(Dg)), states=[])
    for J in range(0, Jmax + 1):
        if J % 2 == 1 and J < 4:
            continue
        mg = E.mgap(Dg, J)
        for mfl in sorted(set([mg*1.000001] + [v for v in (1.02, 1.05, 1.1, 1.2) if v > mg])):
            res['states'].append(state(Dg, J, mfl))
    json.dump(res, open(os.path.join(CA, 'gap_m_%g.json' % Dg), 'w'))
    print('gap strip Dg %g: %d states, %.0fs' % (Dg, len(res['states']), time.time() - t0))


if __name__ == '__main__' and len(sys.argv) > 2 and sys.argv[2] == 'gap':
    run_gap(float(sys.argv[1]), int(sys.argv[3]) if len(sys.argv) > 3 else 40)
    sys.exit(0)

if __name__ == '__main__':
    Dg = float(sys.argv[1])
    mfls = [float(v) for v in sys.argv[2].split(',')] if len(sys.argv) > 2 else [1.3, 1.5, 2, 3, 5, 10, 20]
    bmax = float(sys.argv[3]) if len(sys.argv) > 3 else 400.0
    run(Dg, mfls, bmax)
