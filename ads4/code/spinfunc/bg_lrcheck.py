"""
bg_lrcheck (BELOW_FOLD2.md §3.3): leading-Regge checks of the analytic-radial images.
 (i)  EE: the reference (Delta_ref = 2e4) image I_ref(nu) vs Stage 4's exact leading-Regge kernel T_ref = sum_j R_j P_{-1/2+i(nu-ij)}(eta_ref)
      (k = 4: EE3; k = 2: EE1): N(nu) = I_ref/T_ref must be smooth and zero-free through the real zeros of N^_4 (1.068482, 4.545932,
      8.169433, 13.094015); printed: |N| around each zero and the max relative curvature jump.
 (ii) scalar, large J (bg_scalarJ, resummed in J/Delta): approach of state/ref to T_st/T_ref (rank-1 scalar kernel sc1) along
      fixed bhat = (2J+1)/m: m^2 [log(R_bg/R_LR)] must stay O(1) (the 1/m^2 corrections) as m grows.
Usage: python3 bg_lrcheck.py Delta J K   (needs cache/bg_eeq_D{Delta}_J{J}_K{K}.pkl)  -> cache/bg_lrcheck.txt
"""
import os, sys, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
RQ = os.path.join(HERE, '..', 'resultq'); sys.path.insert(0, RQ)
import numpy as np
import bf_qan as BQ, bg_eeq as EQ

ZEROS = (1.068482, 4.545932, 8.169433, 13.094015)
lines = []


def pr(s=''):
    print(s, flush=True); lines.append(s)


def kern(nus, eta, key):
    import stage4_hj as HJ4
    F = HJ4.numfuncs(key)
    R = {j: np.array([complex(F[j](complex(v))) for v in nus]) for j in F}
    return BQ.T(nus, eta, R), sum(R.values())


def ee_check(D, J, K):
    db = pickle.load(open(EQ.fn(D, J, K), 'rb'))
    rf = db['ref']
    nus = rf['nus']
    sel = np.where((nus < 16))[0]
    _, eta = BQ.mass(EQ.DREF, J)
    for kind, key in (('EE4', 'EE3'), ('EE2', 'EE1')):
        I, _ = EQ.value(rf, kind)
        T, NH = kern(nus[sel], eta, key)
        N = I[sel]/T
        ph = N/np.abs(N)
        pr(f'## {kind} (Stage-4 kernel {key}), Delta_ref = 2e4, J = {J}: N = I_ref/T_ref on nu < 16')
        pr(f'   |N| range {np.abs(N).min():.4e} .. {np.abs(N).max():.4e};  phase spread of N/|N|: {np.ptp(np.angle(ph*np.conj(ph[0]))):.3e} rad')
        lN = np.log(N)
        d2 = np.abs(np.diff(lN, 2))
        pr(f'   max |second difference of log N| (dnu = 0.04): {d2.max():.3e} at nu = {nus[sel][1:-1][np.argmax(d2)]:.2f}; median {np.median(d2):.3e}')
        for z in ZEROS:
            i = np.argmin(np.abs(nus[sel] - z))
            pr(f'   near zero {z}: |T| = {abs(T[i]):.3e}, |I_ref| = {abs(I[sel][i]):.3e}, |N| = {abs(N[i]):.4e} (neighbours {abs(N[i-1]):.4e}, {abs(N[i+1]):.4e})')
        # sign changes of the real-phased I_ref vs those of T
        rot = I[sel]*np.conj(ph[0])/np.abs(N[0])
        sc = nus[sel][1:][np.sign(rot.real[1:]) != np.sign(rot.real[:-1])]
        scT = nus[sel][1:][np.sign(T.real[1:]) != np.sign(T.real[:-1])]
        pr(f'   sign changes of Re I_ref (phase-fixed): {np.round(sc, 2)};   of Re T_ref: {np.round(scT, 2)}')


if __name__ == '__main__':
    D, J, K = float(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    ee_check(D, J, K)
    open(os.path.join(HERE, 'cache', 'bg_lrcheck.txt'), 'a').write('\n'.join(lines) + '\n')
