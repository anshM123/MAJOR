"""
rv_k2eng (RESULT_Q_REVISION.md §1.3): the k = 2 (C'_EE) graviton correction BELOW THE FOLD from the closed-form analytic-radial engine
(bg_eeq EE2 images, already cached by BELOW_FOLD2: state K = 4, reference Delta = 2e4 K = 2; nu < 41 on the dnu = 0.04 grid), projected on
Result Q's smeared weight with the exact Abel (min-norm) weight of Q's psi2.

Chain (the bg_qan chain with the k = 2 kernel; read-only imports bf_qan, bg_eeq, stage4_hj, u3_chi_EE):
  T(nu; eta) = sum_{j = +-1} R_j^{EE1}(nu) P_{-1/2 + i(nu - ij)}(eta)   (Stage 4 exact leading-Regge EE k = 2 kernel; R_j = u2's h_-+)
  Gamma_exact = m^2 (c T_ref I_st/I_ref - T_st),  c: Gamma_exact(nu_min) = 0
  HYB's own beyond-LR term: gamma2 sum_j R_j (nu - ij)^2 (+ the reference's (m/m_ref)^2 copy); gamma2 is FITTED on nu in [2.5, 12]
     (the flat-exact/LR_flat ratio predicts +1/4 without the (2 - y) envelope in the image, -1/4 with it) and reported.
  Gamma_K = Gamma_exact - HYB - const N^ (Gamma_K(nu_min) = 0; N^ = R_-1 + R_+1: the a_h normalisation class)
  correction of the smeared k = 2 piece from the window nu < nu_w:
     dB2_w = (1/Lam) int_0^{nu_w} dnu env(x) Re[ X(nu) Gamma_K(nu)/N^(nu) ] / m^2 = (1/Lam) int env Re[ Y Gamma_K / n ] / m^2
  with X = a Y the min-norm Abel weight of Q's psi2 (u3_chi_EE.SolverAB, chi-normalisation, finite through the real zeros of N^),
  and, for comparison, the same window with the LOCAL weight psi2 and Result D's K2-chain Gamma (rv_smear model) at x = nu/m.
Also printed: Gamma_K / N^ (per harmonic, c2) next to the K2-chain Gamma at the overlapping x (nu = 15 ... 40).
Usage: python3 rv_k2eng.py  -> cache/rv_k2eng.txt (+ cache/rv_k2eng_T_*.pkl, rv_k2eng_Y_<Dg>.pkl)
"""
import os, sys, io, pickle, time
os.environ.setdefault('OMP_NUM_THREADS', '2')
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
SF = os.path.join(HERE, '..', 'spinfunc'); sys.path.insert(0, SF); sys.path.insert(0, os.path.join(HERE, '..', 'h7i14'))
import numpy as np

CACHE = os.path.join(HERE, 'cache')
buf = io.StringIO()


def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); buf.write(s + "\n")


def kernel(nus):
    import stage4_hj as HJ4
    F = HJ4.numfuncs('EE1')
    return {j: np.array([complex(F[j](complex(v))) for v in nus]) for j in F}


def Tcache(nus, eta, R, key):
    import bf_qan as BQ
    fn = os.path.join(CACHE, 'rv_k2eng_T_%s.pkl' % key)
    if os.path.exists(fn):
        d = pickle.load(open(fn, 'rb'))
        if len(d) == len(nus):
            return d
    v = BQ.T(nus, eta, R)
    pickle.dump(v, open(fn, 'wb'))
    return v


def Ysol(Dg, c2, Lam):
    fn = os.path.join(CACHE, 'rv_k2eng_Y_%d.pkl' % Dg)
    if os.path.exists(fn):
        return pickle.load(open(fn, 'rb'))
    import u3_chi_EE as U3
    import v4_chi as C
    dnu, N = C.SETTINGS[Dg]
    Rm, Rp = U3.Rfun('EE')
    S = U3.SolverAB(c2, Lam, dnu, N, Rm, Rp)
    t0 = time.time()
    Y, res, its, _ = S.solve_lsqr()
    lhs, rhs, _ = S.V3Y(Y)
    d = dict(nu=S.nu, Y=Y, n=S.n, a=S.a, phi=S.phi, res=res, its=its, V3=(lhs, rhs), t=time.time() - t0)
    pickle.dump(d, open(fn, 'wb'))
    return d


def main():
    import bf_qan as BQ, bg_eeq as EQ
    import stage6_model as M
    import rv_smear as RS
    F, c2, c3 = M.load('Q3-full-delta0.2')
    psi2 = lambda p: np.polynomial.polynomial.polyval(p, c2)
    models = RS.load_models()
    cases = [(300, 300.0, 0), (300, 314.925, 0), (300, 314.935, 2), (1000, 1000.0, 0), (1000, 1049.93, 2)]
    Ys = {}
    for Dg, D, J in cases:
        Lam = M.lam(Dg)
        if Dg not in Ys:
            Ys[Dg] = Ysol(Dg, c2, Lam)
            y = Ys[Dg]
            pr(f"## Dg={Dg}: Abel min-norm weight of Q's psi2: LSQR residual {y['res']:.1e} ({y['its']} its, {y['t']:.0f}s); "
               f"V3 graviton-pole identity rel {abs(y['V3'][0] - y['V3'][1])/abs(y['V3'][1]):.1e}")
        y = Ys[Dg]
        db = pickle.load(open(EQ.fn(D, J, 4), 'rb'))
        st, rf = db['st'], db['ref']
        n = 1024
        nus = st['nus'][:n]
        m, eta = BQ.mass(D, J); mref, etaref = BQ.mass(EQ.DREF, J)
        R = kernel(nus)
        Tst = Tcache(nus, eta, R, 'D%g_J%d' % (D, J)); Tref = Tcache(nus, etaref, R, 'D%g_J%d' % (EQ.DREF, J))
        NH = sum(R.values())
        S2 = sum(R[j]*(nus - 1j*j)**2 for j in R)
        out = {}
        for Kt in range(0, 5):
            Ist, _ = EQ.value(st, 'EE2', Kt=Kt); Iref, _ = EQ.value(rf, 'EE2')
            tot = Tref*Ist[:n]/Iref[:n]
            cc = Tst[0]/tot[0]
            Gex = m*m*(cc*tot - Tst)
            sel = (nus > 2.5) & (nus < 12) & (np.abs(nus - 6.643) > 0.6)
            # fit gamma2 (and a const N^ term) on the window: Gex ~ g (1 - (m/mref)^2) S2 + k NH
            A = np.array([(1 - (m/mref)**2)*S2[sel], NH[sel]]).T
            coef = np.linalg.lstsq(np.concatenate([A.real, A.imag]), np.concatenate([Gex[sel].real, Gex[sel].imag]), rcond=None)[0]
            out[Kt] = dict(Gex=Gex, g2=coef[0])
        g2 = out[4]['g2']
        pr(f"\n## Dg={Dg} Delta={D} J={J} (m = {m:.3f} = {m/Lam:.4f} Lam): fitted gamma2 = {g2:+.4f} (K=3: {out[3]['g2']:+.4f}); "
           f"flat-exact/LR_flat prediction with the (2-y) envelope: -0.25")
        for gam, lab in ((g2, 'fitted'), (-0.25, 'pred -1/4')):
            res = {}
            for Kt in (3, 4):
                Gex = out[Kt]['Gex']
                HYB = gam*(1 - (m/mref)**2)*S2
                GK = Gex - HYB
                GK = GK - (GK[0]/NH[0])*NH
                res[Kt] = GK
            GK = res[4]
            # Abel-weight window integral (Y/n) and local psi2 window, nu < nu_w
            Yi = np.interp(nus, y['nu'], y['Y'].real) + 1j*np.interp(nus, y['nu'], y['Y'].imag)
            ni = np.interp(nus, y['nu'], y['n'].real) + 1j*np.interp(nus, y['nu'], y['n'].imag)
            x = nus/m
            envx = 2 - x*x
            B2 = None
            for nuw in (10.0, 15.0, 25.0, 41.0):
                s_ = nus < nuw
                dB = 0.04/Lam*np.sum(envx[s_]*(Yi[s_]*GK[s_]/ni[s_]).real)/m**2
                dB3 = 0.04/Lam*np.sum(envx[s_]*(Yi[s_]*res[3][s_]/ni[s_]).real)/m**2
                # Result D K2-chain model on the same window, local psi2 weight
                Gc = models[0][J][0](x[s_], m)[('C2', '++')]
                dBc = 0.04/Lam*np.sum(envx[s_]*psi2(nus[s_]/Lam)*Gc)/m**2
                pr(f"   [{lab}] window nu < {nuw:4.1f}: engine (Abel X^Q) dB2 = {dB:+.3e} (K=3: {dB3:+.3e}); K2-chain model, local psi2: {dBc:+.3e}")
            # per-harmonic comparison at the overlap
            if lab == 'fitted':
                pr("   per harmonic Re Gamma_K/N^ (engine, c2) vs K2-chain Gamma_C2,++ at x = nu/m:   nu: engine [a^4 step] | chain")
                for v in (0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 8.0, 10.0, 12.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0):
                    i = int(np.argmin(np.abs(nus - v)))
                    e4 = (res[4][i]/NH[i]).real; e3 = (res[3][i]/NH[i]).real
                    ch = float(models[0][J][0](np.array([nus[i]/m]), m)[('C2', '++')][0])
                    pr(f"      {nus[i]:6.2f}: {e4:+9.3f} [{e4 - e3:+.3f}] | {ch:+9.3f}")
                pr(f"   max |Re Gamma_K/N^| on nu < 41 excluding |nu - zeros| < 0.3: "
                   f"{np.max(np.abs((GK/NH).real)[(np.abs(nus - 1.1715) > 0.3) & (np.abs(nus - 6.643) > 0.3)]):.2f} c2")
    open(os.path.join(CACHE, 'rv_k2eng.txt'), 'w').write(buf.getvalue())


if __name__ == '__main__':
    main()
