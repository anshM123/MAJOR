"""cl_analyse: post-processing of out/box_<tag>.json (no new integrals).
 (1) per-state images vs e^{pi nu/2} nu^{-(A + sL - 1)}  (the t-image mechanism);
 (2) block-wise truncation error |SB_N - CL| vs N (power-law fit) and the J_max needed for |SB_N - CL| <= |CL|;
 (3) spin-data perturbations: c_J -> c_J (1 + eps w_J) for several profiles w_J; the induced change of the light side per unit eps;
 (4) dimension perturbations (first order, per unit g): consistent shift vs inconsistent shifts (from cl_box);
 (5) the smeared functional: the same quantities vs Lam, and their slope in log10 per unit Lam against pi/(2 ln 10) = 0.682.
usage: python3 -B cl_analyse.py [tag] -> out/analyse_<tag>.txt
"""
import os, sys, json
import numpy as np
from math import lgamma
HERE = os.path.dirname(os.path.abspath(__file__))
tag = sys.argv[1] if len(sys.argv) > 1 else 'main'
d = json.load(open(os.path.join(HERE, 'out', 'box_%s.json' % tag)))
P = d['params']; Df, sL, sH = P['Df'], P['sL'], P['sH']; A = 2*Df + 0.5; N = P['N']
L = []
def log(s=''):
    print(s); L.append(s)
Js = np.arange(N + 1)
profiles = {'uniform (all J)': np.ones(N + 1), 'J >= 2 only': (Js >= 2)*1.0, 'J >= 20 only': (Js >= 20)*1.0,
            'alternating (-1)^(J/2)': np.where((Js//2) % 2 == 0, 1.0, -1.0), 'smooth 1/(1+J)': 1/(1 + Js),
            'tail J^-1 beyond 40': np.where(Js > 40, (40/np.maximum(Js, 1))**1.0, 0.0)}
log('# cl_analyse tag=%s  (Df=%s, sL=%s, sH=%s, N=%d)' % (tag, Df, sL, sH, N))
log('\n(1) per-state images:  max_J |I_J(nu)|  vs  E(nu) = e^{pi nu/2} nu^{-(A+sL-1)}  (ratio should be ~ nu-independent at large nu)')
for nu, r in d['nu'].items():
    nu = float(nu); I = np.array(r['I']); big = np.max(np.abs(I))
    E = np.exp(np.pi*nu/2)*nu**(-(A + sL - 1))
    log('  nu=%5.1f  max|I_J| = %.3e  (J=%d)  |CL| = %.3e  ratio max|I|/E = %.3e  max|I|/|CL| = %.2e' % (nu, big, int(np.argmax(np.abs(I))), abs(r['CL']), big/E, big/abs(r['CL'])))
log('\n(2) block-wise truncation error |SB_N - CL| (even N), and power-law fit |SB_N - CL| ~ C N^-p on N in [40, %d]' % N)
for nu, r in d['nu'].items():
    I = np.array(r['I']); S = np.cumsum(I); err = np.abs(S - r['CL'])
    Ns = np.arange(40, N + 1, 2)
    # envelope: running max over windows of 10 to remove oscillation
    env = np.array([err[n - 8:n + 1].max() for n in Ns])
    p, c = np.polyfit(np.log(Ns), np.log(env), 1)
    need = np.exp((np.log(abs(r['CL'])) - c)/p) if env[-1] > abs(r['CL']) else None
    log('  nu=%5.1f  err(N=40)=%.2e err(N=%d)=%.2e  fit p=%.2f  -> N needed for err <= |CL|=%.1e: %s'
        % (float(nu), err[40], N, err[N], -p, abs(r['CL']), ('%.1e' % need) if need else 'reached'))
log('\n(3) spin-data perturbations c_J -> c_J(1 + eps w_J): d(light side)/d eps, compared with |CL|')
for nu, r in d['nu'].items():
    I = np.array(r['I'])
    log('  nu=%5.1f  ' % float(nu) + '  '.join('%s: %+.2e' % (k, np.sum(I*w)) for k, w in profiles.items()) + '   |CL|=%.2e' % abs(r['CL']))
log('\n(3b) large-J data replaced by an asymptotic fit: c_J exact for J <= J0, c_J -> C J^-p beyond (fit on even J in [J0/2, J0]);'
    ' light side with these data minus the exact collective value')
cJ = np.array(d['cJ'])
for J0 in (10, 20, 40):
    ev = np.arange(J0//2 + (J0//2) % 2, J0 + 1, 2)
    sgn = np.sign(cJ[J0])
    pp, cc = np.polyfit(np.log(ev), np.log(np.abs(cJ[ev])), 1)
    capx = cJ.copy()
    big = (Js > J0) & (Js % 2 == 0)
    capx[big] = sgn*np.exp(cc)*Js[big]**pp
    ratio = np.where(cJ != 0, capx/np.where(cJ != 0, cJ, 1), 0.0)
    relerr = np.max(np.abs(ratio[big] - 1))
    out = []
    for nu, r in d['nu'].items():
        I = np.array(r['I'])
        out.append('nu=%s: %+.2e (|CL| %.1e)' % (nu, np.sum(I*(ratio - 1)), abs(r['CL'])))
    log('  J0=%d (fit exponent %.3f, max rel. data error %.1e): ' % (J0, pp, relerr) + '  '.join(out))
log('\n(4) light dimensions moved by g (first order, per unit g), compared with |CL|:')
log('    consistent = toy with sL -> sL + g (residue function and heavy side move with it; sum rule exact);')
log('    pole-only  = pole moved, residue function held fixed;  gamma_J = g c(J) per state, residues fixed (block-wise)')
for nu, r in d['nu'].items():
    dI = np.array(r['dI'])
    log('  nu=%5.1f  consistent %+.2e | pole-only %+.2e | gamma_J=g: %+.2e | g*2/(2+J): %+.2e | g/(1+J)^2: %+.2e | g*(-1)^(J/2): %+.2e | |CL| %.2e'
        % (float(nu), r['dcons'], r['dpole'], dI.sum(), (dI*2/(2 + Js)).sum(), (dI/(1 + Js)**2).sum(),
           (dI*np.where((Js//2) % 2 == 0, 1, -1)).sum(), abs(r['CL'])))
log('\n(5) smeared functional Omega_Lam = int dp psi2(p) C_{Lam p}')
rows = []
for Lam, r in d['smeared'].items():
    I = np.array(r['I']); S = np.cumsum(I); dI = np.array(r['dI'])
    rows.append((float(Lam), abs(r['CL']), np.max(np.abs(I)), np.sum(np.abs(I)), abs(S[-1] - r['CL']), abs(r['dcons']),
                 abs((dI*2/(2 + Js)).sum()), abs(np.sum(I*profiles['J >= 2 only']))))
    log('  Lam=%4.0f  |Omega_light| (collective) = %.3e | max state %.3e | sum|states| %.3e | trunc err N=%d %.2e | '
        'd/dg consistent %.2e | d/dg gamma_J=2g/(2+J) %.2e | d/deps spin (J>=2) %.2e'
        % ((float(Lam), abs(r['CL']), np.max(np.abs(I)), np.sum(np.abs(I)), N, abs(S[-1] - r['CL']), abs(r['dcons']),
            abs((dI*2/(2 + Js)).sum()), abs(np.sum(I*profiles['J >= 2 only'])))))
rows = np.array(rows)
if len(rows) >= 2:
    names = ['collective', 'max state', 'sum|states|', 'trunc err', 'consistent d/dg', 'inconsistent d/dg', 'spin d/deps']
    log('  slopes d log10(.)/d Lam between the last two Lam (pi/(2 ln10) = 0.682):  ' +
        '  '.join('%s %.3f' % (n, (np.log10(rows[-1, i + 1]) - np.log10(rows[-2, i + 1]))/(rows[-1, 0] - rows[-2, 0])) for i, n in enumerate(names)))
open(os.path.join(HERE, 'out', 'analyse_%s.txt' % tag), 'w').write('\n'.join(L) + '\n')
