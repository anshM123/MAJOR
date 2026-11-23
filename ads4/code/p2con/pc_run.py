"""
p2con/pc_run.py (P2_CONTRACT.md §1.4-§1.6): the TRUE O(1/m^2) block-side correction of every Theorem-W sample state, through the exact
leading-Regge three-channel EE kernel (pc_engine), with the exact drift-free heavy symbols (pc_sym), and its exact cost.

Per state (lsw/cache/costs_1000.json: J, m, beta, floors; read-only):
  1. symbols at q = e^{beta/2} (pc_gensym cache), drift-free (W1c congruence), seed s2 (CHMRS T, Delta_phi = 3);
  2. raw images of the O(1) and O(eps), O(eps^2) block tensors of every entry on 130 nu nodes (exact; side tables beta-independent);
     checks: O(1) image = N_h x 2pi(R_- P_{l-1} + R_+ P_{l+1}) (published R), helicity-diagonal; O(eps) image = 0;
  3. de-shift reconstruction (PW_IMAGES PW2): raw = R_-^sc F(nu+i) + R_+^sc F(nu-i), F(mu) = a(mu) P_mu + b(mu) P'_mu with
        a, b in span{1, mu^2, 1/(4mu^2 + (2j+1)^2), j = 0..7}   (exact fit; residual and held-out nodes reported);
     normalisation F_0 = c P (so f = F/c, f_0 = P_beta: LARGESPIN_W c2 units);
  4. exact cost (PW4): cost = 2 [ sum psi2 F/c dnu/Lam + D F(i/2) ] / m^2  (and the X-route 2 Re sum X A[F/c] dnu/Lam / m^2 as a check);
  5. class content of f = F/c: large-mu polynomial part  kB eta nu^2 P + kC2 (eta^2-1) nu^2 P' + kA eta P + kC P'  (asymptotic,
     canonical), remainder f_rem = O(1/mu^2) (multiplier tail incl. the monopole Q = F(i/2)), with exact costs of each.
Output: cache/run_states.json (+ log).  Usage: python3 pc_run.py [i0 i1]  (state index range; default all)
"""
import os, sys, json, time, pickle
sys.dont_write_bytecode = True
os.environ.setdefault('OMP_NUM_THREADS', '1')
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', 'n2sub')); sys.path.insert(0, os.path.join(HERE, '..', 'lsw'))
import numpy as np
import sympy as sp
import pc_engine as E, pc_sym as Sy, pc_contract as Ct
from pc_gensym import qrat
import n2_common as N
import lsw_common as W

CACHE = os.path.join(HERE, 'cache')
NODES = np.concatenate([np.linspace(0.05, 6, 40), np.linspace(6.2, 40, 50), np.geomspace(41, 1500, 40)])
HOLD = np.array([0.37, 1.91, 4.43, 9.7, 17.3, 28.1, 63.0, 211.0, 777.0])
JS = list(range(0, 8))
D1000 = json.load(open(os.path.join(HERE, '..', 'pwimg', 'cache', 'pw3_cost.json')))['1000']['D']
cfac = lambda m: 32*(4*m*m + 1)*(4*m*m + 9)/((4*m*m + 25)*(4*m*m + 81))
Rm_, Rp_ = N.Rk()


def get_tables(nus, tag):
    fn = os.path.join(CACHE, 'tables_%s.pkl' % tag)
    if os.path.exists(fn):
        return pickle.load(open(fn, 'rb'))
    t = Ct.tables(nus, E.ALPHAS3, E.ALPHAS3)
    pickle.dump(t, open(fn, 'wb'))
    return t


def basis(m, npoly=2):
    return [m**(2*k) for k in range(npoly)] + [1/(4*m*m + (2*j + 1)**2) for j in JS]


def legP(beta, mu):
    """P_mu, P_mu' (d/deta) for l = -1/2 + i mu (complex mu ok)."""
    eta = np.cosh(beta)
    mu = np.atleast_1d(np.asarray(mu, complex))
    nph = int(600 + 3*np.max(np.abs(mu.real))*beta + 50*np.max(np.abs(mu.imag)))
    P = N.Pnu(beta, mu, nph); Pm = N.Pnu(beta, mu + 1j, nph)
    l = -0.5 + 1j*mu
    return P, l*(eta*P - Pm)/(eta*eta - 1)


def design(nus, beta):
    eta = np.cosh(beta)
    nus = np.asarray(nus, float)
    Rms = Rm_(nus + 0j)/cfac(nus + 1j); Rps = Rp_(nus + 0j)/cfac(nus - 1j)
    l = -0.5 + 1j*nus
    P0, dP0 = legP(beta, nus); Pm = N.Pnu(beta, nus + 1j, int(600 + 3*nus.max()*beta)); Pp = N.Pnu(beta, nus - 1j, int(600 + 3*nus.max()*beta))
    dPm = -l*P0 + eta*dP0; dPp = (l + 1)*P0 + eta*dP0
    cols = []
    for (Xm, Xp) in ((Pm, Pp), (dPm, dPp)):
        bm = basis(nus + 1j); bp = basis(nus - 1j)
        for k in range(len(bm)):
            cols.append(Rms*bm[k]*Xm + Rps*bp[k]*Xp)
    return np.array(cols).T


def fitF(y, A, nus):
    w = 1/(1 + np.asarray(nus)**2)
    Ar = np.vstack([(A*w[:, None]).real, (A*w[:, None]).imag]); yr = np.concatenate([(y*w).real, (y*w).imag])
    sc = np.linalg.norm(Ar, axis=0)
    x, *_ = np.linalg.lstsq(Ar/sc, yr, rcond=None); x = x/sc
    res = np.linalg.norm(Ar@x - yr)/max(np.linalg.norm(yr), 1e-300)
    return x, res


def Feval(x, beta, mu):
    mu = np.atleast_1d(np.asarray(mu, complex))
    nb = len(x)//2
    b = basis(mu); a = sum(x[k]*b[k] for k in range(nb)); bb = sum(x[nb + k]*b[k] for k in range(nb))
    P, dP = legP(beta, mu)
    return a*P + bb*dP, a, bb, P, dP


def F_half(x, beta):
    """F(i/2) = a(i/2) [P_{-1} = 1] + lim b(mu) P'_mu  (the j = 0 pole of b times P'_mu ~ -i eps/(1+eta), mu = i/2 + eps)."""
    eta = np.cosh(beta); nb = len(x)//2
    m = 0.5j
    a = x[0] + x[1]*m*m + sum(x[2 + j]/(4*m*m + (2*j + 1)**2) for j in JS if j != 0)
    r0b = x[nb + 2 + 0]
    a_pole = x[2 + 0]
    return a - r0b/(4*(1 + eta)), a_pole


def costs(x, beta, m):
    nu, Xw, ps, dn, L = W.NU, W.X, W.PS, W.DN, W.L
    F, a, b, P, dP = Feval(x, beta, nu)
    f = F/cfac(nu)
    Q, apole = F_half(x, beta)
    tot = 2*((ps*f).real.sum()*dn/L + D1000*Q.real)/m**2
    # X route: A[f] = (R- f(nu+i) + R+ f(nu-i))/(R- + R+)
    Fm = Feval(x, beta, nu + 1j)[0]/cfac(nu + 1j); Fp = Feval(x, beta, nu - 1j)[0]/cfac(nu - 1j)
    Am = Rm_(nu + 0j); Ap = Rp_(nu + 0j)
    totX = 2*(Xw*(Am*Fm + Ap*Fp)/(Am + Ap)).real.sum()*dn/L/m**2
    # polynomial (ring) part of f and its cost (R1 form, strip-analytic: no defect)
    nb = len(x)//2
    eta = np.cosh(beta)
    kB = x[1]/(32*eta); kC2 = x[nb + 1]/(32*(eta*eta - 1)); kA = (x[0] + 24*x[1])/(32*eta); kC = (x[nb] + 24*x[nb + 1])/32
    poly = (kB*eta*nu**2 + kA*eta)*P + (kC2*(eta*eta - 1)*nu**2 + kC)*dP
    cB = 2*(ps*(kB*eta*nu**2*P)).real.sum()*dn/L/m**2
    cC2 = 2*(ps*(kC2*(eta*eta - 1)*nu**2*dP)).real.sum()*dn/L/m**2
    cA = 2*(ps*(kA*eta*P)).real.sum()*dn/L/m**2
    cC = 2*(ps*(kC*dP)).real.sum()*dn/L/m**2
    cpoly = cB + cC2 + cA + cC
    return dict(total=float(tot), totalX=float(totX), Q=float(Q.real), Qim=float(Q.imag), apole=float(apole), kB=float(kB), kC2=float(kC2),
                kA=float(kA), kC=float(kC), cB=float(cB), cC2=float(cC2), cA=float(cA), cC=float(cC), cpoly=float(cpoly),
                crem=float(tot - cpoly), cmono=float(2*D1000*Q.real/m**2))


def run_state(r, tabN, tabH, log):
    t0 = time.time()
    beta_s = r['beta']; m = r['m']; J = r['J']
    q = qrat(beta_s); beta = float(2*sp.log(q))
    sym = Sy.load_q(q)
    Lt, info = Sy.drift_free(sym)
    S = Sy.numeric_series(Lt, q); seed = Sy.seed_terms(q)
    even = J % 2 == 0
    hs = ('same', 'opp') if even else ('odd',)
    MN = Ct.Mpair(tabN['nus'], beta); MH = Ct.Mpair(tabH['nus'], beta)
    shpN = Ct.lead_shape(tabN['nus'], beta)
    AN = design(tabN['nus'], beta); AH = design(tabH['nus'], beta)
    out = dict(J=J, m=m, mfl=r['mfl'], beta=beta, beta_state=beta_s, q=str(q), fl={h: v['floor'] for h, v in r['fl'].items()})
    out['drift'] = {str(k): complex(sp.N(v)).real for k, v in info.get('even', {}).get('U1', {}).items()}
    Nh = {}
    for h in hs:
        r0 = Ct.raw_image(Sy.block_tensor(S, h, h, 0), tabN, beta, MN)
        rat = r0/shpN
        Nh[h] = float(np.median(rat.real))
        out.setdefault('checks', {})['lead_%s' % h] = float(np.max(np.abs(rat - Nh[h]))/abs(Nh[h]))
    if even:
        r0 = Ct.raw_image(Sy.block_tensor(S, 'same', 'opp', 0), tabN, beta, MN)
        out['checks']['lead_offdiag'] = float(np.max(np.abs(r0/shpN))/Nh['same'])
    pairs = [(h, hp) for i, h in enumerate(hs) for hp in hs[i:]]
    ent = {}
    for (h, hp) in pairs:
        nrm = np.sqrt(Nh[h]*Nh[hp])
        r1 = Ct.raw_image(Sy.block_tensor(S, h, hp, 1), tabN, beta, MN)
        out['checks']['eps1_%s_%s' % (h, hp)] = float(np.max(np.abs(r1/(nrm*Ct.norm_hpm(tabN['nus'])))))
        for piece, sd, kind in (('full', seed, 'both'), ('symbol', None, 'sym'), ('seed', seed, 'seed')):
            if piece != 'full' and h != hp:
                continue
            if kind == 'seed':
                T = {k: Sy.mul_terms(v, seed) for k, v in Sy.block_tensor(S, h, hp, 0).items()}
            else:
                T = Sy.block_tensor(S, h, hp, 2, sd)
            y = Ct.raw_image(T, tabN, beta, MN)/(2*np.pi*nrm)
            yH = Ct.raw_image(T, tabH, beta, MH)/(2*np.pi*nrm)
            x, res = fitF(y, AN, tabN['nus'])
            hold = float(np.max(np.abs(AH@x - yH))/np.max(np.abs(yH)))
            c = costs(x, beta, m)
            c.update(fitres=float(res), hold=hold)
            ent['%s|%s|%s' % (h, hp, piece)] = dict(c, x=[float(v) for v in x])
    out['ent'] = ent
    out['Nh'] = Nh
    out['time'] = time.time() - t0
    log(f"J={J:3d} m/L={r['mfl']:.3f} beta={beta:.4f}  checks: " + ', '.join(f'{k}={v:.1e}' for k, v in out['checks'].items())
        + f"  [{out['time']:.0f}s]")
    for k, v in ent.items():
        log(f"    {k:22s} fit {v['fitres']:.1e} hold {v['hold']:.1e} | total {v['total']:+.3e} (X {v['totalX']:+.3e}) | kB {v['kB']:+.5f} kC2 {v['kC2']:+.5f} "
            f"kA {v['kA']:+.3f} kC {v['kC']:+.3f} Q {v['Q']:+.3f} | cpoly {v['cpoly']:+.2e} crem {v['crem']:+.2e} cmono {v['cmono']:+.1e}")
    return out


def main():
    R = json.load(open(os.path.join(HERE, '..', 'lsw', 'cache', 'costs_1000.json')))
    i0, i1 = (int(sys.argv[1]), int(sys.argv[2])) if len(sys.argv) > 2 else (0, len(R))
    tabN = get_tables(NODES, 'nodes'); tabH = get_tables(HOLD, 'hold')
    logf = open(os.path.join(CACHE, 'run_%d_%d.log' % (i0, i1)), 'w')
    def log(s):
        print(s, flush=True); logf.write(s + '\n'); logf.flush()
    log(f'tables: tail {tabN["tail"]:.1e} / {tabH["tail"]:.1e}')
    res = []
    for r in R[i0:i1]:
        res.append(run_state(r, tabN, tabH, log))
        json.dump(res, open(os.path.join(CACHE, 'run_%d_%d.json' % (i0, i1)), 'w'))


if __name__ == '__main__':
    main()
