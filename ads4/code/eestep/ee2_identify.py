r"""
ee2_identify (EESTEP_CLOSE.md §1): the exact EE per-nu kernel as a finite combination of shifted scalar continuous-spin partial waves.

(a) Reduction.  For random polarisation sets: (i) Q(c) = P_even (X20 X40)^4 is a polynomial of degree 8 in c = cos th (Chebyshev tail);
    (ii) the sin-odd part of P integrates to zero against S_nu; (iii) K^EE_nu (direct quadrature with the dressing) = sum_{b=-4}^{4} c_b Psi^(b)_nu
    (direct quadrature with t^b inserted).
(b) Identification of Psi^(b) (prefactor-free, as lrc1):  Q_b(x) = [Psi^(b)_nu(x)/Psi^(b)_nu0(x)] / [Gh_{b'}(x; nu)/Gh_{b'}(x; nu0)] must be
    x-independent for the right b'; tested b' = +b, -b, 0 over 9 points (rho in {0.25, 0.55, 0.85}, beta in {0.6, 1.4, 2.2}).
(c) Prefactor:  Psi^(b)_nu0 / Gh_{b'}(nu0) = kappa r^alpha v^gamma  (v = (1-z)(1-zb)); fitted alpha, gamma.
(d) End to end:  K^EE_nu(x)/K^EE_nu0(x)  vs  [sum_b c_b(x) v^{gamma_b} Gh_{b'}(x; nu)] / [same at nu0]  * N(nu)/N(nu0) (N from the scalar
    b = 0 kernel at one point), over the 9 points and 3 polarisation sets; plus the same with the leading-Regge P_l only (control).
Output: cache/ee2_identify.txt.
"""
import os, sys, time
sys.dont_write_bytecode = True
os.environ.setdefault('OMP_NUM_THREADS', '1')
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np
import ee_blocks as EB
import ee_quad as Q

OUT = []
T0 = time.time()


def pr(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); OUT.append(s)


NUS = np.array([0.4, 0.8, 1.3, 2.7]); NU0 = 0.8; I0 = 1
PTS = [(rho, beta) for rho in (0.25, 0.55, 0.85) for beta in (0.6, 1.4, 2.2)]
rng = np.random.default_rng(11)
POLS = [tuple(rng.uniform(0, 2*np.pi, 4)) for _ in range(3)]
BS = list(range(-4, 5))


def main():
    rows = []
    for (rho, beta) in PTS:
        r = rho*np.exp(-beta)
        parts = Q.kernel_parts(r, beta, nth=192, ds=0.012)
        th = parts['th']
        t, s = Q.t_of_th(th, r, beta)
        Psi = {b: Q.integrate(parts, NUS, t**b) for b in BS}
        G = {bp: np.array([EB.Gh_point(nu, bp, rho, beta, N=200) for nu in NUS]) for bp in BS}
        ee = []
        for ip, pa in enumerate(POLS):
            pol = Q.pols(beta, *pa)
            P = Q.dressing(th, r, beta, pol)
            Podd = 0.5*(P - Q.dressing(-th, r, beta, pol))
            K = Q.integrate(parts, NUS, P)
            Kodd = Q.integrate(parts, NUS, Podd)
            cb, tail, oddsz, q = Q.laurent_cb(r, beta, pol)
            Ksum = sum(cb[b]*Psi[b] for b in BS)
            ee.append(dict(K=K, Kodd=Kodd, Ksum=Ksum, cb=cb, tail=tail, oddsz=oddsz))
        rows.append(dict(rho=rho, beta=beta, r=r, ok=parts['ok'], Psi=Psi, G=G, ee=ee))
        pr(f'  point rho = {rho}, beta = {beta} (r = {r:.3e}, s = {s:.3f}) done [{time.time() - T0:.0f}s]; fixed signs/A > 0: {parts["ok"]}')

    pr('(a) exact reduction to b-shifted insertions (3 random polarisation sets x 9 points)')
    tails = [e['tail'] for rw in rows for e in rw['ee']]
    odd = [np.max(np.abs(e['Kodd'])/np.abs(e['K'])) for rw in rows for e in rw['ee']]
    oddsz = [e['oddsz'] for rw in rows for e in rw['ee']]
    rel = [np.max(np.abs(e['Ksum']/e['K'] - 1)) for rw in rows for e in rw['ee']]
    amp = [max(abs(e['cb'][b]*rw['Psi'][b][I0]) for b in BS)/abs(e['K'][I0]) for rw in rows for e in rw['ee']]
    pr(f'    (i)   Chebyshev tail of Q(c) beyond degree 8, relative: max {max(tails):.1e}')
    pr(f'    (ii)  sin-odd part: size relative to even part up to {max(oddsz):.2f}; its integral / |K^EE|: max {max(odd):.1e}')
    pr(f'    (iii) K^EE (direct) vs sum_b c_b Psi^(b) (direct): max rel diff {max(rel):.1e}  (cancellation in the b-sum: max_b|c_b Psi^(b)|/|K| up to {max(amp):.1e})')

    pr('(b) identification of Psi^(b) = int S_nu t^b: x-spread of Q_b over the 9 points for blocks with b\' = +b, -b, 0')
    gam = {}
    for b in BS:
        line = []
        for bp in sorted(set([b, -b, 0])):
            spreads = []
            for i, nu in enumerate(NUS):
                if i == I0:
                    continue
                qs = np.array([(rw['Psi'][b][i]/rw['Psi'][b][I0])/(rw['G'][bp][i]/rw['G'][bp][I0]) for rw in rows])
                spreads.append(np.max(np.abs(qs/np.mean(qs) - 1)))
            line.append(f"b'={bp:+d}: {max(spreads):.1e}")
        pr(f'    b = {b:+d}: ' + ';  '.join(line))
    pr('(c) prefactor Psi^(b)_nu0 / Gh_{b\'}(nu0) = kappa r^alpha v^gamma, with b\' the identified label')
    for b in BS:
        bp = -b
        lr = np.log([rw['r'] for rw in rows]); lv = np.log([EB.v_cross(rw['rho'], rw['beta']) for rw in rows])
        pf = np.array([rw['Psi'][b][I0]/rw['G'][bp][I0] for rw in rows])
        A = np.vstack([lr, lv, np.ones_like(lr)]).T
        coef, *_ = np.linalg.lstsq(A, np.log(np.abs(pf)), rcond=None)
        res = np.max(np.abs(A @ coef - np.log(np.abs(pf))))
        gam[b] = (coef[0], coef[1], np.exp(coef[2])*np.exp(1j*np.mean(np.angle(pf))))
        pr(f"    b = {b:+d} (b' = {bp:+d}): alpha = {coef[0]:.6f}, gamma = {coef[1]:+.6f}, max residual {res:.1e}, phase spread {np.ptp(np.angle(pf)):.1e}")

    pr('(d) end to end, exact prefactors: K^EE_nu(x) (direct quadrature)  vs  pred = Nk(nu) r^2 sum_b c_b(x) v^{-b/2} Gh_{-b}(x; nu)')
    pr('    (recursion blocks; Nk(nu) = Psi^(0)_nu/(r^2 Gh_0) at the first point only; deviation relative to |K| and to sum_b |c_b Psi^(b)|)')
    rw0 = rows[0]
    Nk = rw0['Psi'][0]/(rw0['r']**2*rw0['G'][0])
    for ip in range(len(POLS)):
        devs = []; devs2 = []; devLR = []; amps = []
        for rw in rows:
            e = rw['ee'][ip]
            v = EB.v_cross(rw['rho'], rw['beta'])
            terms = [e['cb'][b]*rw['r']**2*v**(-b/2)*rw['G'][-b] for b in BS]
            pred = Nk*sum(terms)
            scale = np.abs(Nk)*sum(np.abs(tm) for tm in terms)
            devs.append(np.max(np.abs(e['K']/pred - 1)))
            devs2.append(np.max(np.abs(e['K'] - pred)/scale))
            amps.append(np.max(scale/np.abs(pred)))
            ctl = (e['K']/e['K'][I0])/(rw['G'][0]/rw['G'][0][I0]*(Nk/Nk[I0]))
            devLR.append(np.max(np.abs(ctl - 1)))
        pr(f'    polarisation set {ip} (angles {np.round(POLS[ip], 3).tolist()}): max |K/pred - 1| = {max(devs):.1e}; '
           f'max |K - pred|/sum|c_b Psi_b| = {max(devs2):.1e} (cancellation up to {max(amps):.1e})   [control, scalar nu-profile: {max(devLR):.1e}]')


if __name__ == '__main__':
    pr('# ee2_identify: EE per-nu kernel = finite sum of b-shifted scalar continuous-spin partial waves (EESTEP_CLOSE.md §1)')
    main()
    pr(f'  [{time.time() - T0:.0f}s]')
    with open(os.path.join(HERE, 'cache', 'ee2_identify.txt'), 'w') as fh:
        fh.write('\n'.join(OUT) + '\n')
