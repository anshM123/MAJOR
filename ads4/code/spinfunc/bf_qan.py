"""
bf_qan (BELOW_FOLD.md §4): from the bf_q sums, the exact (all-orders) EE k = 4 forward image of a low-J '++' state on nu < 6, the beyond-
model correction Gamma_K(nu) = Gamma_exact,tot - Gamma_HYB,tot (c2 units, Stage-7 convention: image = sum_j R_j F_j,
F_j = 2[P(nu - ij) + Gamma_j/m^2], Gamma_tot = sum_j R_j Gamma_j), and its projections on the de-shift weights (all gauges).

Normalisation.  total_st(nu) = 2 c I_st(nu)/N(nu), N(nu) := I_ref(nu)/T_ref(nu) (T = sum_j R_j P_{-1/2+i(nu-ij)}(cosh beta), Delta_ref =
2e4), smoothed per replica by a polynomial fit in nu^2 (it is smooth iff the D-integral reproduces the zero 1.0685 of N^_4: checked).
c (one constant per state; the a_h class, invisible to the derivative piece since int chi_d N^_4 = 0) is fixed by Gamma_K(nu_min) = 0.
HYB (Stage 6, '++'): Gamma_HYB,tot = gamma sum_j R_j (nu - ij)^2, gamma = (c4 + 1/2)/2 = 0.75.
Usage: python3 bf_qan.py Delta J [Dg]  -> prints; appends to cache/bf_qan.txt
"""
import os, sys, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
RQ = os.path.join(HERE, '..', 'resultq'); sys.path.insert(0, RQ)
import warnings
warnings.filterwarnings('ignore')
import numpy as np
import mpmath as mp
import bf_q as Q

NUMAX_FIT = 3.6
ZERO = 1.068482


def legP(mu, eta):
    a = mp.mpf(1)/2 + 1j*mp.mpc(mu)
    return complex(mp.hyp2f1(a, 1 - a, 1, (1 - mp.mpf(eta))/2))


def mass(D, J):
    al = D - J - 2.0
    return np.sqrt(al*(al + 2*J + 1)), (al**2 + (al + 2*J + 1)**2)/(2*al*(al + 2*J + 1))


def Rj(nus):
    import stage4_hj as HJ4
    F4 = HJ4.numfuncs('EE3')
    return {j: np.array([complex(F4[j](complex(v))) for v in nus]) for j in F4}


def T(nus, eta, R):
    return sum(R[j]*np.array([legP(v - 1j*j, eta) for v in nus]) for j in R)


def weights():
    """{gauge: (nu, chi_value, chi_deriv)} on dnu = 0.04 midpoints."""
    G = pickle.load(open(os.path.join(RQ, 'cache', 'stage6_gauge.pkl'), 'rb'))
    Cn = pickle.load(open(os.path.join(RQ, 'cache', 'stage6_canon.pkl'), 'rb'))
    out = {}
    nu = G[('C4', 0.04, 'nu')]
    for k in ('H0 min-norm', 'H2 (1+(nu/5)^2)', 'H3 smooth', 'H4 amp-L1'):
        out[k.split()[0]] = (nu, G[('C4v', 0.04)][k]['chi'], G[('C4d', 0.04)][k]['chi'])
    out['SB'] = (Cn['C4d']['nu_04'], Cn['C4v']['chi_04'], Cn['C4d']['chi_04'])
    xv = pickle.load(open(os.path.join(RQ, 'cache', 'stage7_xprec_c4v_d004_N3072.pkl'), 'rb'))
    xd = pickle.load(open(os.path.join(RQ, 'cache', 'stage7_xprec_c4d_d004_N3072.pkl'), 'rb'))
    cv = np.array([float(mp.mpf(x)) for x in xv['weights']]); cd = np.array([float(mp.mpf(x)) for x in xd['weights']])
    out['XP'] = ((np.arange(len(cd)) + 0.5)*0.04, cv, cd)
    return out


def gamma_K(D, J, tag='', deg=8):
    db = pickle.load(open(Q.fn(D, J, tag), 'rb'))
    nus = db['nus']; seeds = sorted(db['reps'])
    m, eta = mass(D, J); mref, etaref = mass(Q.DREF, J)
    R = Rj(nus)
    Tst, Tref = T(nus, eta, R), T(nus, etaref, R)
    NH = sum(R.values())
    HYB = 0.75*sum(R[j]*(nus - 1j*j)**2 for j in R)
    sel = (nus < NUMAX_FIT) & (np.abs(nus - ZERO) > 0.12)
    out = []; Nraw = []
    for s in seeds:
        Ist, Iref = db['reps'][s][(D, 'EE4')], db['reps'][s][(Q.DREF, 'EE4')]
        N = Iref/Tref; Nraw.append(N)
        x = nus[sel]**2
        cr = np.polyfit(x, N[sel].real, deg); ci = np.polyfit(x, N[sel].imag, deg)
        Nf = np.polyval(cr, nus**2) + 1j*np.polyval(ci, nus**2)
        out.append(Ist/Nf)
    tot = np.array(out)                                   # = total_st/(2c)
    c = Tst[0]/tot.mean(0)[0]                             # Gamma_exact(nu_min) = 0 convention, then shift to Gamma_K(nu_min) = 0
    Gex = m*m*(c*tot - Tst[None, :])
    GK = Gex - HYB[None, :]
    GK = GK - (GK.mean(0)[0]/NH[0])*NH[None, :]           # a_h class: Gamma_K(nu_min) = 0 (adds const x N^_4 only)
    return dict(nus=nus, GK=GK, Gex=Gex, m=m, NH=NH, Nraw=np.array(Nraw), sel=sel)


def project(r, W, numax):
    out = {}
    for g, (nu, cv, cd) in W.items():
        idx = [int(round((v - 0.02)/0.04)) for v in r['nus'] if v < numax]
        GK = r['GK'][:, :len(idx)]
        sd = 0.04*np.sum(cd[idx][None, :]*GK.real, axis=1)
        sv = 0.04*np.sum(cv[idx][None, :]*GK.real, axis=1)
        A = 0.04*np.sum(np.abs(cd[idx])*np.abs(r['NH'][:len(idx)]))            # amplification of a uniform |Gamma| on [0, numax)
        Atot = 0.04*np.sum(np.abs(cd)*np.abs(np.real(sum(Rj(nu[:1]).values()))))  # placeholder (not used)
        out[g] = (sd.mean(), sd.std(ddof=1)/np.sqrt(len(sd)), sv.mean(), sv.std(ddof=1)/np.sqrt(len(sv)), A)
    return out


if __name__ == '__main__':
    D, J = float(sys.argv[1]), int(sys.argv[2])
    r = gamma_K(D, J)
    lines = []
    def pr(s=''):
        print(s, flush=True); lines.append(s)
    pr(f'# bf_qan: EE k=4 on D, Delta={D:.4f} J={J} m={r["m"]:.3f}, {r["GK"].shape[0]} replicas')
    Nm = r['Nraw'].mean(0)
    pr('  smoothness of N = I_ref/T_ref through the zero 1.0685 of N^_4 (|N| at nu = 0.94 .. 1.22):  ' +
       ' '.join(f'{abs(Nm[i]):.3e}' for i in range(len(r['nus'])) if 0.9 < r['nus'][i] < 1.25))
    pr('  nu     Gamma_K (Re, Im) +- err      Gamma_exact Re')
    GKm = r['GK'].mean(0); GKe = r['GK'].std(0, ddof=1)/np.sqrt(r['GK'].shape[0])
    for i, v in enumerate(r['nus']):
        if i % 5 == 0:
            pr(f'  {v:5.2f}  {GKm[i].real:+9.2f} {GKm[i].imag:+9.2f}  +- {GKe[i]:7.2f}     {r["Gex"].mean(0)[i].real:+9.2f}')
    W = weights()
    for numax in (2.5, 3.0, 3.5):
        P = project(r, W, numax)
        pr(f'  projections on nu < {numax}: s_d = int chi_d Re Gamma_K (+- QMC);  s_v (a_h-convention dependent)')
        for g, (sd, ed, sv, ev, A) in P.items():
            pr(f'     {g:4s}: s_d = {sd:+8.4f} +- {ed:.4f}   s_v = {sv:+8.3f} +- {ev:.3f}')
    open(os.path.join(HERE, 'cache', 'bf_qan.txt'), 'a').write('\n'.join(lines) + '\n')
    pickle.dump(r, open(os.path.join(HERE, 'cache', 'bf_qan_D%g_J%d.pkl' % (D, J)), 'wb'))
