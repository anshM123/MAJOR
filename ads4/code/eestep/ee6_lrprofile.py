r"""
ee6_lrprofile (EESTEP_CLOSE.md §1.4, §2): two small facts used in the note.
(a) The leading-Regge EE position kernel L_nu(beta) (r^0 Cauchy coefficient of E_nu, as in ee3) is NOT proportional to P_l(cosh beta):
    the ratio L/P_l depends on nu (the (13) light-ray polarisations make the f-stripped LR kernel a multi-harmonic object).
(b) The nu-independent coefficients C_b(x) of ee3 at fixed rho tend to beta-independent constants (uniformity along the light-cone edge).
Output: cache/ee6_lrprofile.txt.
"""
import os, sys
sys.dont_write_bytecode = True
os.environ.setdefault('OMP_NUM_THREADS', '1')
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np
import ee_blocks as EB
import ee_quad as Q
import ee3_deshift as D3

OUT = []


def pr(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); OUT.append(s)


def L(nu, beta, pa, N=200, MC=32):
    coeffs = [EB.block_coeffs_b(nu, b, N) for b in range(5)]
    Pt = EB.legendre_scaled(nu, beta, N + 1); Ys = [EB.radial_Y(B, m, Pt, beta) for B, m in coeffs]
    R = 0.3*np.exp(-beta); vals = []
    for j in range(MC):
        rc = R*np.exp(2j*np.pi*j/MC); rho = rc*np.exp(beta); rb = rc*np.exp(-beta)
        v = ((1 - rho)*(1 - rb)/((1 + rho)*(1 + rb)))**2
        C = D3.Cb_of(rc, beta, Q.pols(beta, *pa))
        vals.append(sum(C[b]*v**(b/2)*EB.Gh(Ys[b], rho)[0] for b in range(5)))
    return np.mean(vals), Pt[0].real


if __name__ == '__main__':
    pr('# ee6_lrprofile')
    pr('(a) L_nu(beta)/P_l(cosh beta) for the two polarisation sets of ee3 (scalar-like LR would make this nu-independent)')
    for pa in D3.POLS:
        for beta in (1.0, 3.0, 6.0):
            rat = []
            for nu in (0.5, 1.5, 3.0, 6.0):
                l, p = L(nu, beta, pa); rat.append(l.real/p)
            pr(f'    pols {np.round(pa, 2).tolist()}, beta = {beta}: nu = 0.5, 1.5, 3, 6: ' + ' '.join(f'{x:+.4e}' for x in rat))
    pr('(b) max_b |C_b(rho, beta)| at fixed rho vs beta (polarisations boosted with xh; beyond beta ~ 13 boost round-off sets in)')
    for ip, pa in enumerate(D3.POLS):
        for rho in (0.2, 0.5, 0.85):
            row = []
            for beta in (2, 4, 8, 12):
                r = rho*np.exp(-beta)
                row.append(np.max(np.abs(D3.Cb_of(r, beta, Q.pols(beta, *pa)))))
            pr(f'    set {ip}, rho = {rho}: beta = 2, 4, 8, 12: ' + ' '.join(f'{x:.5e}' for x in row))
    with open(os.path.join(HERE, 'cache', 'ee6_lrprofile.txt'), 'w') as fh:
        fh.write('\n'.join(OUT) + '\n')
