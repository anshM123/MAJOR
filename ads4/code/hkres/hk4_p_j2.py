r"""
hk4_p_j2 (HK_RESIDUAL.md §4): F_NORM §2.6's deciding test for (P) at J = 2 (and J = 0 as the control against EESTEP's direct engine).

EESTEP §5 ran J = 0 with the DIRECT 7-d quadrature (repro_engine/direct_rho.py), which is J = 0 only.  For J = 2 the only engine in the
repo is the analytic-radial one (spinfunc/bg_eeq.py, BELOW_FOLD2 §2.3: exact order by order in a = 1/Delta, conserved <TTO_J> for any
even J), imported read-only; images are computed here (kind EE4 only) and pickled to hkres/cache (bg_eeq's own cache is not touched).
The chain is EESTEP §5(b)'s, unchanged: r_D(nu) = I(D; nu)/I(D; 0.02); r_inf = Richardson of r_6000 and r_20000 in 1/m^2 (same J);
rho = r_D/r_inf, T_ref -> T_inf = sum_j R_j;  Gamma_K = m^2 (c T_inf rho - T_st) - 0.75 sum_j R_j (nu - ij)^2 - const N^_4;
s_d = 0.04 sum_{nu<12} chi_d Re Gamma_K (Stage-6 H0 weight, Delta_g-independent), m^2 = (D-J-2)(D+J-1), cosh beta = 1 + (2J+1)^2/(2m^2).
Truncation: the a-expansion is kept to order K; |s_d(K) - s_d(K-1)| is printed (it multiplies m^2: relative size a^{K-1}).
The C5 piece is NOT run: no engine in the repo evaluates the k = 5 (j' = 1) heavy image at finite m (direct_rho and bg_eeq are k = 4 / 2).
Usage:  python3 hk4_p_j2.py run D J K depth      -> cache/hk4_eeq_D{D}_J{J}_K{K}.pkl
        python3 hk4_p_j2.py report K             -> cache/hk4_p_j2.txt
"""
import os, sys, time, pickle
sys.dont_write_bytecode = True
os.environ.setdefault('OMP_NUM_THREADS', '1')
HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.join(HERE, '..')
sys.path.insert(0, os.path.join(CODE, 'repro_engine'))
sys.path.insert(0, os.path.join(CODE, 'spinfunc'))
import numpy as np

CACHE = os.path.join(HERE, 'cache')
DS = (300.0, 1000.0, 3000.0, 6000.0, 2e4)
OUT = []


def pr(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); OUT.append(s)


def fn(D, J, K):
    return os.path.join(CACHE, 'hk4_eeq_D%g_J%d_K%d.pkl' % (D, J, K))


def run(D, J, K, depth):
    import bg_eeq as EQ
    EQ.C = CACHE                                    # never write into spinfunc/cache
    t0 = time.time()
    rec = EQ.images(D, J, K, depth, kinds=('EE4',))
    pickle.dump(rec, open(fn(D, J, K) + '.tmp', 'wb')); os.replace(fn(D, J, K) + '.tmp', fn(D, J, K))
    print(f'hk4 run D={D:g} J={J} K={K} depth={depth}: {time.time() - t0:.0f}s', flush=True)


def report(K):
    import bg_eeq as EQ
    import sdchain as C
    pr(f'# hk4_p_j2: m-dependence of the ++ derivative-piece correction s_d at fixed functional, analytic-radial images (order K = {K})')
    for J in (0, 2):
        recs = {D: pickle.load(open(fn(D, J, K), 'rb')) for D in DS if os.path.exists(fn(D, J, K))}
        if len(recs) < len(DS):
            pr(f'  J = {J}: missing runs {[D for D in DS if D not in recs]}'); continue
        n = 1024
        nus = recs[DS[0]]['nus'][:n]
        sel = nus < 12.0
        nus = nus[sel]
        R = C.Rj(nus); NH = sum(R.values()); Tinf = sum(R.values())
        res = {}
        for Kt in (K - 1, K):
            rr = {}
            for D, rec in recs.items():
                I, _ = EQ.value(rec, 'EE4', Kt=Kt)
                I = I[:n][sel]
                rr[D] = I/I[0]
            ma2 = C.mass(2e4, J)[0]**2; mb2 = C.mass(6000.0, J)[0]**2
            rinf = (ma2*rr[2e4] - mb2*rr[6000.0])/(ma2 - mb2)
            for D in DS[:4]:
                m, eta = C.mass(D, J)
                Tst = sum(R[j]*C.legP(nus - 1j*j, eta) for j in R)
                tot = Tinf*rr[D]/rinf; c = Tst[0]/tot[0]
                GK = m*m*(c*tot - Tst) - 0.75*sum(R[j]*(nus - 1j*j)**2 for j in R)
                GK = GK - (GK[0]/NH[0])*NH
                full = np.zeros(1024, complex); full[:len(nus)] = GK
                nfull = (np.arange(1024) + 0.5)*0.04
                sd, sv = C.sd(nfull, full, 'H0', 0.0, 12.0)
                win = [C.sd(nfull, full, 'H0', lo, hi)[0] for lo, hi in ((0, 1), (1, 2), (2, 3), (3, 4), (4, 12))]
                res[(Kt, D)] = (sd, sv, win, m)
        pr(f'  J = {J}:   D      m/Lambda_300   s_d^inf (nu<12)   s_v^inf    |s_d(K) - s_d(K-1)|   windows [0,1) [1,2) [2,3) [3,4) [4,12)')
        ys = []
        for D in DS[:4]:
            sd, sv, win, m = res[(K, D)]
            dK = abs(sd - res[(K - 1, D)][0])
            ys.append(sd)
            pr(f'         {D:6.0f}    {m/298.5:7.2f}        {sd:+.5f}        {sv:+.4f}       {dK:.1e}            ' + ' '.join(f'{w:+.4f}' for w in win))
        ys = np.array(ys); Dv = np.array(DS[:4])
        A = np.vstack([np.ones_like(Dv), 1/Dv**2]).T
        coef, *_ = np.linalg.lstsq(A, ys, rcond=None)
        pr(f'         fit s_d^inf = {coef[0]:+.5f} + {coef[1]:.1f}/D^2;  spread over m/Lambda = 1..20: {np.ptp(ys):.1e} '
           f'(relative {np.ptp(ys)/abs(np.mean(ys)):.1e})')
    pr('  EESTEP §5 (direct engine, J = 0): s_d^inf = +0.50814 / +0.50952 / +0.50961 / +0.50945 at D = 300 / 1000 / 3000 / 6000')
    with open(os.path.join(CACHE, 'hk4_p_j2.txt'), 'w') as fh:
        fh.write('\n'.join(OUT) + '\n')


if __name__ == '__main__':
    if sys.argv[1] == 'run':
        run(float(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]))
    else:
        report(int(sys.argv[2]))
