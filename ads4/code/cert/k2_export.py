"""Write ads4/prl/data/k2/ (CERT_CLOSE.md §3): X^Q grid weights and the continuum k = 2 weights of D and Q, with README and SHA256SUMS.
Reads: resultq/cache/rv_k2eng_Y_<Dg>.pkl (X^Q grid, read-only), cert/cache/k2_<W>_<Dg>_{re,im}_m2.npz and the k2_cont / k2_xq_check logs.
Usage: python3 k2_export.py"""
import os, sys, json, pickle, hashlib, datetime
import numpy as np
import mpmath as mp
HERE = os.path.dirname(os.path.abspath(__file__)); ADS4 = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, HERE)
import k2_cont as K       # noqa: E402
OUT = os.path.join(ADS4, 'prl', 'data', 'k2')
os.makedirs(OUT, exist_ok=True)


def write_xq(Dg):
    d = pickle.load(open(os.path.join(ADS4, 'code', 'resultq', 'cache', 'rv_k2eng_Y_%d.pkl' % Dg), 'rb'))
    nu, Y, a, phi = d['nu'], d['Y'], d['a'], d['phi']
    X = a * Y
    Lam = K.lam_of(Dg)
    fn = os.path.join(OUT, 'xq_abel_X_%d.txt' % Dg)
    with open(fn, 'w') as f:
        f.write('# Result Q k = 2 Abel weight X^Q = chi (h_- + h_+) of Q3-full-delta0.2 psi2 (min-norm gauge, Y - psi2 of minimal l2 norm on this grid; '
                'u3_chi_EE.SolverAB LSQR, residual %.1e; stored as ads4/code/resultq/cache/rv_k2eng_Y_%d.pkl, RESULT_Q_REVISION §1.3); '
                'Lam = %.10f; grid nu_k = (k + 1/2) %.2f; Y = X/a; psi2 column = psi2^Q(nu/Lam) (no 1/Lam); the smeared k = 2 part of Omega^Q_up is '
                '(1/Lam) Re sum_k dnu X_k Psi^_(EE,nu_k)/(h_- + h_+)(nu_k)\n' % (d['res'], Dg, Lam, nu[1] - nu[0]))
        f.write('# nu_k  Re X  Im X  Re Y  Im Y  psi2(nu_k)\n')
        for i in range(len(nu)):
            f.write('%.3f %.17e %.17e %.17e %.17e %.17e\n' % (nu[i], X[i].real, X[i].imag, Y[i].real, Y[i].imag, phi[i]))
    return fn


def write_cont(W, Dg, info):
    T = K.Target(W, Dg); B = K.Basis(T.Lam)
    xr = np.load(os.path.join(HERE, 'cache', 'k2_%s_%d_re_m2.npz' % (W, Dg)))
    xi = np.load(os.path.join(HERE, 'cache', 'k2_%s_%d_im_m2.npz' % (W, Dg)))
    assert np.allclose(np.array(B.cent), xr['cent']) and np.allclose(np.array(B.cent), xi['cent'])
    # double vs 50-digit evaluation of chi at sample points (rounding check)
    mp.mp.dps = 50
    def chi_mp(nu, x):
        z = mp.mpf(nu); q = z * z / (z * z + mp.mpf(1) / 4) ** 2; s = 0
        for (c, sg), xk in zip(B.cent, x[:len(B.cent)]):
            g = mp.exp(-(z - c) ** 2 / (2 * mp.mpf(sg) ** 2))
            if c > 0: g += mp.exp(-(z + c) ** 2 / (2 * mp.mpf(sg) ** 2))
            s += mp.mpf(float(xk)) * g
        s *= q
        l = (z * z + 1) * mp.log(z * z + 1); e = mp.exp(-z * z / 2)
        for k in range(B.nsing): s += mp.mpf(float(x[len(B.cent) + k])) * z ** (2 * k) * l * e
        return s
    samp = [0.3, 1.0, 5.0, 30.0, 100.0, T.Lam * 0.8, T.Lam - 1]
    cmax = np.abs(K.chi_eval(B, xr['x'], np.linspace(0, T.Lam + 30, 20001))).max()
    rd = max(abs(float(chi_mp(v, xr['x'])) - float(K.chi_eval(B, xr['x'], np.array([v]))[0])) for v in samp) / cmax
    js = dict(
        name='Omega^EE k = 2 continuum de-shift weight (Result D, NF02 psi2/k0)' if W == 'D' else 'Omega^Q_up smeared k = 2 continuum de-shift weight (Q3-full-delta0.2 psi2)',
        Delta_g=Dg, Lam=float(T.Lam),
        target='psi(nu) = psi2^NF02(nu/Lam)/k0, k0 = 27.48509817 (nu <= Lam; 0 beyond)' if W == 'D' else 'psi(nu) = psi2^Q(nu/Lam) (nu <= Lam; 0 beyond)',
        kink_slope_kappa=float(T.kappa),
        identity='W[chi](t) := int_0^oo chi(nu) [R_-(nu) cos((nu+i)t) + R_+(nu) cos((nu-i)t)] dnu = Phi(t) := int_0^oo psi(nu) cos(nu t) dnu; '
                 'R_-+ = prl/data/kernels.json D_k2_EE; X = chi (R_- + R_+) is the weight in the convention of d_abel_X_*.txt',
        form='chi(nu) = chi_re(nu) + i chi_im(nu); each part = Q(nu) sum_k x_k [g_k(nu - c_k) + g_k(nu + c_k)] + sum_{j=0}^{2} y_j nu^(2j) l(nu) exp(-nu^2/2), '
             'Q(nu) = nu^2/(nu^2 + 1/4)^2, g_k(z) = exp(-z^2/(2 sigma_k^2)) (only g(nu) for c_k = 0), l(nu) = (nu^2+1) log(nu^2+1) (principal branch)',
        singularities='poles of order 2 at nu = +-i/2 (double zeros of R_-/R_+ in the crossed strips); branch points at +-i (they carry the |mu| kink of psi); '
                      'double zero at 0 (cancels the 1/nu poles of R_-+); analytic in |Im nu| < 1 otherwise; super-Gaussian decay beyond nu ~ Lam + 30',
        centres_sigma=[[float(c), float(s)] for c, s in B.cent],
        x_re=[repr(float(v)) for v in xr['x'][:len(B.cent)]], y_re=[repr(float(v)) for v in xr['x'][len(B.cent):]],
        x_im=[repr(float(v)) for v in xi['x'][:len(B.cent)]], y_im=[repr(float(v)) for v in xi['x'][len(B.cent):]],
        eps_W_re=float(xr['eps_W']), eps_W_im_null=float(xi['eps_W']), int_abs_psi=float(xr['intpsi']),
        certificate='for all t >= 0: |W[chi_re](t) - Phi(t)| <= eps_W_re, |W[chi_im](t)| <= eps_W_im_null (shifted form, CERT_CLOSE.md §3.2); '
                    'eps_W = int_0^oo |E(mu) - psi(mu)| dmu, E(mu) = 2 Re[chi(mu - i) R_-(mu - i)], by a 4.2e5-point quadrature in double precision '
                    '(grade [DV], not interval-verified); hence |image - target| <= eps_W P_{-1/2}(cosh beta) for all beta',
        precision='coefficients are IEEE doubles (exact as given); double-precision evaluation of chi_re agrees with 50-digit evaluation to %.1e of max|chi_re| at nu = %s'
                  % (rd, ', '.join('%.4g' % v for v in samp)),
        checks=info)
    fn = os.path.join(OUT, 'k2_continuum_%s_%d.json' % (W, Dg))
    json.dump(js, open(fn, 'w'), indent=0)
    # convenience table
    step = 0.02 if Dg < 500 else 0.04
    nu = np.arange(0, T.Lam + 30, step) + step / 2
    cr = K.chi_eval(B, xr['x'], nu); ci = K.chi_eval(B, xi['x'], nu)
    Rm, Rp = K.kern.load_D(); S = Rm.np(nu) + Rp.np(nu)
    X = (cr + 1j * ci) * S
    tf = os.path.join(OUT, 'k2_continuum_%s_%d_table.txt' % (W, Dg))
    with open(tf, 'w') as f:
        f.write('# continuum k = 2 weight %s, Dg = %d, tabulated for convenience only (the JSON formula is the definition); nu_k = (k + 1/2) %.2f, '
                'same nodes as the stored grid weight\n# nu  Re chi  Im chi  Re X  Im X  psi(nu)\n' % (W, Dg, step))
        for i in range(len(nu)):
            f.write('%.3f %.17e %.17e %.17e %.17e %.17e\n' % (nu[i], cr[i], ci[i], X[i].real, X[i].imag, T(nu[i])))
    return fn, tf


def parse_log(W, Dg):
    lines = open(os.path.join(HERE, 'cache', 'k2_%s_%d.log' % (W, Dg))).read().strip().splitlines()
    return [l for l in lines if 'checks' in l]


if __name__ == '__main__':
    files = []
    for Dg in (300, 1000):
        files.append(write_xq(Dg))
    for W in ('D', 'Q'):
        for Dg in (300, 1000):
            files += list(write_cont(W, Dg, parse_log(W, Dg)))
    print('\n'.join(files))
