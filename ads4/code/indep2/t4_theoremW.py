"""Target 4: Theorem W (LARGESPIN_W §4.0).  (i) the 2x2 / Gershgorin algebra (exact + random tests);
(ii) own Lemma-I costs c_{s,h} = [2 Omega_2[S] + sigma_h u Omega_3[S]]/m^2  (de-shift form, R1), u = Lam^2/m^2, Omega_k[S] =
int_0^1 psi_k(p) S(Lam p; beta) dp with exact NF02 psi_k (prl/data/nf02_psi.json), S in {A: eta P, B: eta nu^2 P, C1: P'(eta),
C2: (eta^2-1) nu^2 P'(eta), C3: (eta^2-1) P'(eta), sh+-k: Re P_{nu-+ik}, 4a: eta^2 nu^4 P (/m^4), 4b: eta sinh(beta) nu^4 d_beta P (/m^4)};
(iii) the J* entries from these costs, with the floors f_h taken (as data) from code/lsw/cache/costs_1000.json.
Output: cache/t4_theoremW.txt
"""
import json
from fractions import Fraction
import numpy as np

out = open('cache/t4_theoremW.txt', 'w')
def log(*a):
    s = ' '.join(str(x) for x in a); print(s); out.write(s + '\n'); out.flush()

# ---------------------------------------------------------------- (i) algebra
rng = np.random.default_rng(3)
viol = 0; tight = 0; ntest = 200000
for _ in range(ntest):
    f1, f2 = rng.uniform(0.01, 1, 2)
    H = f1 * f2 / (f1 + f2)
    d = H * rng.uniform(0, 1)
    a = f1 + d * rng.uniform(-1, 1); b = f2 + d * rng.uniform(-1, 1)
    z = d * rng.uniform(0, 1) * np.exp(2j * np.pi * rng.uniform())
    M = np.array([[a, z], [np.conj(z), b]])
    if np.linalg.eigvalsh(M).min() < -1e-14:
        viol += 1
# tightness: a = f1-d, b = f2-d, |z| = d at d slightly above H -> indefinite
for _ in range(1000):
    f1, f2 = rng.uniform(0.01, 1, 2); H = f1 * f2 / (f1 + f2); d = H * 1.001
    M = np.array([[f1 - d, d], [d, f2 - d]])
    tight += np.linalg.eigvalsh(M).min() < 0
log(f'(i) 2x2: {ntest} random Hermitian matrices with |diag - f_h| <= d, |offdiag| <= d, d <= H: {viol} indefinite; '
    f'worst case at d = 1.001 H indefinite in {tight}/1000 (criterion is sharp)')
# Gershgorin 3x3 and the "lowers tolerances by <= 1.5" claim
viol3 = 0
for _ in range(100000):
    f = rng.uniform(0.01, 1, 3); d = f.min() / 3 * rng.uniform(0, 1)
    M = np.diag(f + d * rng.uniform(-1, 1, 3)).astype(complex)
    for (i, j) in ((0, 1), (0, 2), (1, 2)):
        M[i, j] = d * rng.uniform() * np.exp(2j * np.pi * rng.uniform()); M[j, i] = np.conj(M[i, j])
    viol3 += np.linalg.eigvalsh(M).min() < -1e-14
log(f'(i) 3x3 Gershgorin (d <= min f/3): {viol3}/100000 indefinite')
# exact worst-case 3x3 tolerance vs Gershgorin and vs H(2x2)
def exact_tol3(f):
    lo, hi = 0.0, min(f)
    for _ in range(60):
        d = (lo + hi) / 2
        M = np.diag(np.array(f) - d) + d * (np.ones((3, 3)) - np.eye(3))   # worst real case: all off-diag = +d
        ok = np.linalg.eigvalsh(M).min() >= 0
        # check sign patterns (phases) too: worst is a frustrated or unfrustrated pattern
        M2 = np.diag(np.array(f) - d) + d * np.array([[0, 1, 1], [1, 0, -1], [1, -1, 0]])
        ok = ok and np.linalg.eigvalsh(M2).min() >= 0
        lo, hi = (d, hi) if ok else (lo, d)
    return lo
for f in ([1, 1, 1], [1, 1, 10], [1, 10, 10], [1, 3, 1], [0.5, 1, 1]):
    Hf = f[0] * f[1] / (f[0] + f[1])
    log(f'    floors (f++, f+-, fOO) = {f}: H(2x2 of ++,+-) = {Hf:.4f}, Gershgorin min f/3 = {min(f)/3:.4f}, ratio H/(min f/3) = {Hf/(min(f)/3):.3f}; exact 3x3 uniform-d tolerance = {exact_tol3(f):.4f}')

# ---------------------------------------------------------------- (ii) costs
js = json.load(open('../../prl/data/nf02_psi.json'))
k0 = float(js['k0'])
psi = {2: np.array([float(Fraction(c)) for c in js['psi2_exact']]) / k0, 3: np.array([float(Fraction(c)) for c in js['psi3_exact']]) / k0}
Lam = np.sqrt(998.0 * 999.0)

def mehler(mu, beta, nper=4.0):
    mu = np.asarray(mu, complex)
    n = int(max(300, nper * np.max(np.abs(mu.real)) * beta / np.pi + 200))
    x, w = np.polynomial.legendre.leggauss(n)
    U = np.sqrt(beta); u = 0.5 * U * (x + 1); wu = 0.5 * U * w
    s = beta - u * u
    den = np.sqrt(2 * np.sinh((beta + s) / 2) * np.sinh(u * u / 2))
    g = 2 * u / den
    wt = wu * g * np.sqrt(2) / np.pi
    return np.cos(np.outer(mu, s)) @ wt

NP = 4000
px, pw = np.polynomial.legendre.leggauss(NP); px = 0.5 * (px + 1); pw = 0.5 * pw
psiv = {k: np.polyval(psi[k][::-1], px) for k in (2, 3)}

def costs_state(beta, m):
    nu = Lam * px
    eta = np.cosh(beta); sh = np.sinh(beta)
    P = mehler(nu, beta).real
    Pp1 = mehler(nu + 1j, beta)            # P_{-1/2+i(nu+i)} = P_{mu-1}
    Pm1 = mehler(nu - 1j, beta)
    Pp2 = mehler(nu + 2j, beta); Pm2 = mehler(nu - 2j, beta)
    mu = -0.5 + 1j * nu
    e2P1 = (mu * (eta * P - Pp1)).real     # (eta^2-1) P'(eta)
    dP = e2P1 / (eta * eta - 1)            # P'(eta)
    dbP = e2P1 / sh                        # d_beta P
    S = {'A ': eta * P, 'B ': eta * nu**2 * P, 'C1': dP, 'C2': nu**2 * e2P1, 'C3': e2P1,
         'sh+1': Pm1.real, 'sh-1': Pp1.real, 'sh+2': Pm2.real, 'sh-2': Pp2.real,
         '4a': eta**2 * nu**4 * P, '4b': eta * sh * nu**4 * dbP}
    u = Lam**2 / m**2
    res = {}
    for k, Sv in S.items():
        O2 = np.sum(pw * psiv[2] * Sv); O3 = np.sum(pw * psiv[3] * Sv)
        sc = m**4 if k in ('4a', '4b') else m**2
        # [k=2 only (Route W convention), ++ with k=3, +- with k=3]
        res[k] = [2 * O2 / sc, (2 * O2 + u * O3) / sc, (2 * O2 - u * O3) / sc]
    # leading: 2 Omega2 + sigma u Omega3 (for the record)
    res['LR'] = [2 * np.sum(pw * psiv[2] * P) + u * np.sum(pw * psiv[3] * P), 2 * np.sum(pw * psiv[2] * P) - u * np.sum(pw * psiv[3] * P)]
    return res


# stripped (R2) costs: the normalised raw image A[.] is replaced by S itself; cost = 2 Re sum_k X_k S(nu_k) dnu/Lam / m^2 (k = 2 only)
import sympy as _sp
_kj = json.load(open('../../prl/data/kernels.json'))['D_k2_EE']
_nus = _sp.Symbol('nu')
_Rm = _sp.lambdify(_nus, _sp.sympify(_kj['R_minus'], locals={'nu': _nus}), 'numpy')
_Rp = _sp.lambdify(_nus, _sp.sympify(_kj['R_plus'], locals={'nu': _nus}), 'numpy')
_a = np.loadtxt('../../prl/data/d_abel_X_1000.txt')
gnu = _a[:, 0]; gX = _a[:, 1] + 1j * _a[:, 2]; gdnu = gnu[1] - gnu[0]
gchi = gX / (_Rm(gnu + 0j) + _Rp(gnu + 0j))
def mehler_chunk(mu, beta, chunk=2048):
    return np.concatenate([mehler(mu[i:i + chunk], beta) for i in range(0, len(mu), chunk)])
def costs_stripped(beta, m):
    nu = gnu
    eta = np.cosh(beta); sh = np.sinh(beta)
    P = mehler_chunk(nu + 0j, beta).real
    Pp1 = mehler_chunk(nu + 1j, beta); Pm1 = mehler_chunk(nu - 1j, beta)
    Pp2 = mehler_chunk(nu + 2j, beta); Pm2 = mehler_chunk(nu - 2j, beta)
    mu = -0.5 + 1j * nu
    e2P1 = (mu * (eta * P - Pp1)).real
    S = {'A ': eta * P, 'B ': eta * nu**2 * P, 'C1': e2P1 / (eta * eta - 1), 'C2': nu**2 * e2P1, 'C3': e2P1,
         'sh+1': Pm1.real, 'sh-1': Pp1.real, 'sh+2': Pm2.real, 'sh-2': Pp2.real,
         '4a': eta**2 * nu**4 * P, '4b': eta * nu**4 * e2P1}
    res = {}
    for k, Sv in S.items():
        sc = m**4 if k in ('4a', '4b') else m**2
        res[k] = 2 * float(np.real(np.sum(gX * Sv)) * gdnu / Lam) / sc
    return res

states = json.load(open('../lsw/cache/costs_1000.json'))
log('\n(ii) own Lemma-I costs (R1, k = 2 part 2 Omega_2[S]/m^2; O(1/m^4) shapes /m^4) vs Route W cache costs_1000.json "cost"[k][0] (R1 entry)')
mine = []
worst_rel = {}
for st in states:
    c = costs_state(st['beta'], st['m'])
    mine.append(c)
    line = []
    for k in st['cost']:
        ref = st['cost'][k][0]; v = c[k][0]
        ratio = ref / v if v != 0 else float('nan')
        line.append(f"{k.strip()} {v:+.3e}/{ref:+.3e}")
        if k not in ('4a', '4b'):
            worst_rel[k] = max(worst_rel.get(k, 0), abs(v - ref) / abs(ref))
    lr = 2 * np.sum(pw * psiv[2] * mehler(Lam * px, st['beta']).real)
    hk = '++' if '++' in st['fl'] else '+-'
    r4 = st['cost']['4a'][0] / c['4a'][0]
    log(f"  J={st['J']:3d} m/L={st['mfl']:.3f} b={st['beta']:.4f}: 2Om2 mine {lr:.6e} / Ns {st['fl'][hk]['Ns']:.6e}; " + ', '.join(line) +
        f";  4a ratio theirs/mine = {r4:.4e} (m^2 = {st['m']**2:.4e})")
strip = []
wr2 = {}
for st in states:
    cs = costs_stripped(st['beta'], st['m'])
    strip.append(cs)
    for k in st['cost']:
        ref = st['cost'][k][1] / (st['m']**2 if k in ('4a', '4b') else 1)
        wr2[k] = max(wr2.get(k, 0), abs(cs[k] - ref) / abs(ref))
log('  stripped (R2) costs vs cache "cost"[k][1] (4a/4b cache divided by m^2), worst relative difference: ' + ', '.join(f'{k.strip()}: {v:.1e}' for k, v in wr2.items()))
log('  worst relative difference per structure (R1, k = 2): ' + ', '.join(f'{k.strip()}: {v:.1e}' for k, v in worst_rel.items()))
k3 = max(abs(c[k][1] - c[k][0]) / abs(c[k][0]) for c in mine for k in ('A ', 'B ', 'C1', 'C3'))
k3c2 = max(abs(c['C2'][1] - c['C2'][0]) / abs(c['C2'][0]) for c in mine)
log(f'  effect of the omitted k = 3 piece (u Omega_3[S]) on the costs: up to {k3:.2f} relative for A, B, C1, C3; up to {k3c2:.2f} for C2')

# ---------------------------------------------------------------- (iii) J* from own costs
kap = {'A ': 6.46, 'B ': 11 / 6, 'C1': 1 / 6, 'C2': 1 / 6, 'C3': 11.75 / 6}
def jstar(Lams, Ksh, K4, theta, use_mine=True, variant=0, m4fix=1.0, ringR2=False, shR2=False, verbose=False):
    fails = []
    for st, c, cs in zip(states, mine, strip):
        cc = c if use_mine else st['cost']
        hs = (variant,) if (use_mine and variant == 0) else ((1, 2) if use_mine else (0,))
        dh = []
        for h in hs:
            fac4 = m4fix if use_mine else m4fix / st['m']**2
            g = (lambda k: abs(cs[k])) if (use_mine and ringR2) else (lambda k: abs(cc[k][h]))
            gs = (lambda k: abs(cs[k])) if (use_mine and shR2) else (lambda k: abs(cc[k][h]))
            d = sum(Lams[k] * kap[k] * g(k) for k in kap)
            d += Ksh * sum(gs(k) for k in ('sh+1', 'sh-1', 'sh+2', 'sh-2'))
            d += K4 * fac4 * (g('4a') + g('4b'))
            dh.append(d)
        d = max(dh)
        if st['J'] % 2 == 0:
            fpp, fpm = st['fl']['++']['floor'], st['fl']['+-']['floor']
            H = fpp * fpm / (fpp + fpm)
        else:
            H = st['fl']['+-']['floor']
        if d > theta * H:
            fails.append(st['J'])
        if verbose and st['J'] in verbose:
            log(f"      J={st['J']} m/L={st['mfl']:.3f}: d = {d:.3e}, theta H = {theta*H:.3e}, d/(theta H) = {d/(theta*H):.3f}")
    Js = sorted(set(s['J'] for s in states))
    fs = sorted(set(fails))
    js_ = None
    for J0 in Js:
        if all(J < J0 for J in fs):
            js_ = J0; break
    return js_, fs
log('\n(iii) J* with own costs (floors from costs_1000.json); claimed values from LARGESPIN_W §4.3 / P2_SPINNING §7.5 in brackets')
cases = [((1, 10, 1), 1.0, '30'), ((1, 10, 1), 0.5, '38'), ((3, 10, 10), 1.0, '50'), ((3, 10, 10), 0.5, '100'),
         ((10, 30, 10), 1.0, '64'), ((10, 30, 10), 0.5, 'none (fails 300, 344)')]
for (La, Ksh, K4), th, claim in cases:
    Lams = {k: La for k in kap}
    j1, f1 = jstar(Lams, Ksh, K4, th); j2, f2 = jstar(Lams, Ksh, K4, th, use_mine=False)
    j3, f3 = jstar(Lams, Ksh, K4, th, variant=1)
    j4, f4 = jstar(Lams, Ksh, K4, th, shR2=True)
    j5, f5 = jstar(Lams, Ksh, K4, th, shR2=True, ringR2=True)
    log(f'    ring R1 + index-shifted STRIPPED (own costs): J* = {j4}, failing {f4};   all stripped (R2 row): J* = {j5}, failing {f5}')
    log(f'  (Lam_amp, K_sh, K4) = ({La}, {Ksh}, {K4}), theta = {th}: own k=2 costs: J* = {j1}, failing {f1} | own costs incl. k=3 (worse helicity): J* = {j3}, failing {f3} | [claimed {claim}]  (their cached costs, 4a/4b divided by m^2: J* = {j2}, failing {f2})')
log('  margins at the boundary states, (3, 10, 10), theta = 1, own costs with k = 3 included:')
jstar({k: 3 for k in kap}, 10, 10, 1.0, variant=1, verbose=(46, 50, 51, 64))
log('  same, k = 2 only:')
jstar({k: 3 for k in kap}, 10, 10, 1.0, verbose=(46, 50, 51, 64))
# P2_SPINNING §7.5 model constants (Lambda_C1 not stated in §7.5; taken = 1 and = 3.17 as bracket)
for LC1 in (1.0, 3.17):
    Lams = {'A ': 3.17, 'B ': 1.0, 'C1': LC1, 'C2': 1.0, 'C3': 2.03}
    for th, claim in ((1.0, '38'), (0.5, '50')):
        j1, f1 = jstar(Lams, 10, 4, th)
        log(f'  P2 §7.5 model (Lam_A, Lam_B, Lam_C1, Lam_C2, Lam_C3) = (3.17, 1, {LC1}, 1, 2.03), K_sh = 10, K4 = 4, theta = {th}: J* = {j1}, failing {f1} [claimed {claim}]')
# Gershgorin ratio on the actual floors
rat = [ (st['fl']['++']['floor'] * st['fl']['+-']['floor'] / (st['fl']['++']['floor'] + st['fl']['+-']['floor'])) / (min(st['fl']['++']['floor'], st['fl']['+-']['floor']) / 3)
        for st in states if st['J'] % 2 == 0]
fr = [st['fl']['++']['floor'] / st['fl']['+-']['floor'] for st in states if st['J'] % 2 == 0]
log(f'\nGershgorin on the sampled even-J floors: f++/f+- in [{min(fr):.3f}, {max(fr):.3f}]; H/(min(f++,f+-)/3) in [{min(rat):.3f}, {max(rat):.3f}] '
    f'(the claim "lowers tolerances by <= 1.5" needs f_OO >= min(f++, f+-) and f++ ~ f+-)')
