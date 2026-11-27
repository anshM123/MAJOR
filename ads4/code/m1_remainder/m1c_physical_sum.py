r"""
m1_remainder/m1c_physical_sum.py  [NE; exponents computed, O(1) constants A_n, kappa_n set to 1 and flagged]

The ACTUAL term-wise light remainder for tree-level Einstein data, with physical weights, at Delta_g = 1000:

   R(c_T) = sum_{O light} (2 pi G Lam^4/9) sum_b |a_{O,b,0}|^2 hhat_b(O),     8 pi G = 256/c_T,
   R_2   := c_T^2 R = (64 Lam^4/9) sum_O sum_b [c_T |a_{O,b,0}|^2] hhat_b(O)       (c_T-independent at leading order).

 * hhat_b: exact per-state values from m1a (moments of ../fullpos, J <= 1000, every J; even J used).
 * c_T |a|^2: the per-state weights of ../a2prime_einstein/s4_assemble.py (read-only logic, re-implemented here):
       |a_{O,b,0}|^2 = C_T lamhat_X^2 |K_X(Delta,J)|^2,  lamhat_X^2 = A_n 4^{-Delta} 2^J s^{p_X},  gamma = kappa_n/(c_T J),
   K_X the exact closed forms of ../kernels/closed_forms.py (they carry sin^2(pi h) / cos^2(pi h)).  A_n = kappa_n = 1.
   Families (twist, structure, ladder): (2,pp,0) (2,pm,4) (4,pp,0) (4,pm,4) (6,pp,0) (3,odd,0) and the twist-2 N-structure
   on ladder 4 with B_J = 1/J (a2prime S3(c)).
 * Threshold: Lemma 2 of no_eft/NOTE.md, tau^2 <= tau_0^2 (1 + X),  X = pi c_T |R| /(512 n_0 L_A),  16 L_A = 14339.59;
   X <= 1 needs c_T >= pi |R_2| /(512 n_0 L_A).  n_0 (the T-block normalisation) is not tracked: quoted with n_0 = 1.
 * Also: E0 := c_T sum_light sum_b w |a|^2 (unweighted detector share) for the angular route, L_g <= 70 E_light.
Output: out/m1c_physical_sum.txt
"""
import os, sys, json, math
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'kernels'))
import mpmath as mp
import closed_forms as cf

mp.mp.dps = 120
OUT = os.path.join(HERE, 'out')
CT = mp.mpf(10)**40           # only to make gamma tiny; results are multiplied back to the c_T-independent combination
DG = int(os.environ.get('M1_DG', '1000'))
KMAX = 1000 if DG == 1000 else 300
SUF = '' if DG == 1000 else f'_dg{DG}'
LAM = math.sqrt((DG - 2)*(DG - 1))
L16 = {1000: '14339.59', 300: '10734.61'}[DG]


def ctA2(name, h0, J, kappa=1):
    """c_T |a_{O,b,0}|^2 at leading order (c_T-independent)."""
    gam = mp.mpf(kappa)/(CT*J)
    h = mp.mpf(h0) + gam/2
    D = J + 2*h; s = (D + J)/2
    p = mp.mpf(1)/2 if name in ('pp', 'pm') else mp.mpf(3)/2
    lam2 = mp.mpf(4)**(-D)*mp.mpf(2)**J*s**p
    K2 = {'pp': cf.Kpp2, 'pm': cf.Kpm2, 'odd': cf.Kodd2}[name](D, J)
    return CT*CT*lam2*K2, D


def ctA2_N(J, B):
    D = J + 2; s = (D + J)/2
    eps = mp.mpf('1e-15')
    strip = cf.Kpm2(J + 2 + 2*eps, J)/mp.sin(mp.pi*(1 + eps))**2
    return (B**2)*mp.mpf(4)**(-D)*mp.mpf(2)**J*s**mp.mpf(0.5)*strip, D


N0 = {1000: mp.mpf('4.53987e-7'), 300: mp.mpf('6.64037e-6')}


def load(tw, b):
    d = json.load(open(os.path.join(OUT, f'm1a_hhat_{tw}_b{b}{SUF}.json')))
    return {int(J): (v[0], v[1]) for J, v in d.items()}


def main():
    lines = []
    def pr(s):
        print(s, flush=True); lines.append(s)
    FAM = [('pp', 1, 2, 0), ('pm', 1, 2, 4), ('pp', 2, 4, 0), ('pm', 2, 4, 4), ('pp', 3, 6, 0), ('odd', mp.mpf(3)/2, 3, 0)]
    tot_signed = mp.mpf(0); tot_abs = mp.mpf(0); E0 = mp.mpf(0)
    per = {}
    pr('family (struct, h, twist, ladder): sum_{J even<=1000} c_T|a|^2 hhat  [signed, abs], dominant J, and c_T|a|^2 T at J=500, 1000')
    for name, h0, tw, b in FAM + [('N', 1, 2, 4)]:
        H = load(tw, b)
        s_sig = mp.mpf(0); s_abs = mp.mpf(0); best = (0, None)
        for J in range(max(4, b), KMAX + 1, 2):
            if J not in H:
                continue
            sg, lg = H[J]
            if name == 'N':
                w, D = ctA2_N(J, mp.mpf(1)/J)
            else:
                w, D = ctA2(name, h0, J)
            val = w*sg*mp.mpf(10)**lg
            s_sig += val; s_abs += abs(val)
            if abs(val) > best[0]:
                best = (abs(val), J)
            T = (D - 1)*(D - 2)
            E0 += w*((T - b*(b - 1))/3 + mp.mpf(3)/2)
        per[(name, tw, b)] = (s_sig, s_abs)
        tot_signed += s_sig; tot_abs += s_abs
        w5, D5 = (ctA2_N(500, mp.mpf(1)/500) if name == 'N' else ctA2(name, h0, 500))
        w10, D10 = (ctA2_N(1000, mp.mpf(1)/1000) if name == 'N' else ctA2(name, h0, 1000))
        pr(f'  {name:4s} h={float(h0):3.1f} tw={tw} b={b}: signed {mp.nstr(s_sig, 5):>12s}  abs {mp.nstr(s_abs, 5):>12s}  '
           f'max term at J={best[1]}   c_T|a|^2 T: J=500 {mp.nstr(w5*(D5-1)*(D5-2), 3)}, J=1000 {mp.nstr(w10*(D10-1)*(D10-2), 3)}')
    pref = 64*mp.mpf(LAM)**4/9
    R2s = pref*tot_signed; R2a = pref*tot_abs
    pr(f'\nR_2 = c_T^2 R (J <= 1000, even J, A_n = kappa_n = 1):  signed {mp.nstr(R2s, 6)}  (log10 |.| = {float(mp.log10(abs(R2s))):.2f}),'
       f'  term-wise abs {mp.nstr(R2a, 6)} (log10 = {float(mp.log10(R2a)):.2f})')
    # tail J > 1000: per-state |hhat|/T assumed <= 10^{685.7} (plateau sup seen at J <= 1000), weights continued
    PLATEAU = {}
    for tw_ in (2, 3, 4, 6):
        vals = [v[1] - math.log10((J + tw_ - 1)*(J + tw_ - 2)) for b_ in (0, 4) for J, v in load(tw_, b_).items() if J > 0.6*KMAX]
        PLATEAU[tw_] = max(vals)
    pr(f'plateau sup log10|hhat/T| over J > 0.6 Jmax: {PLATEAU}')
    tail = mp.mpf(0)
    for name, h0, tw, b in FAM:
        for J in range(KMAX + 2, 60001, 2) if tw == 2 else range(KMAX + 2, 12001, 2):
            w, D = ctA2(name, h0, J)
            T = (D - 1)*(D - 2)
            tail += w*T*mp.mpf(10)**PLATEAU[tw]
    pr(f'tail J > 1000 (if the plateau |hhat|/T persisted, per family sup values): pref*tail = {mp.nstr(pref*tail, 4)} '
       f'(log10 {float(mp.log10(pref*tail)):.2f})')
    L_A = mp.mpf(L16)/16
    n0 = N0[DG]
    for lab, R2 in (('signed J<=1000', R2s), ('abs J<=1000', R2a), ('abs + tail', R2a + pref*tail)):
        cstar = mp.pi*abs(R2)/(512*L_A)
        pr(f'threshold (X = 1) from {lab:16s}: n_0 = 1: c_T* = {mp.nstr(cstar, 4)} (log10 {float(mp.log10(cstar)):.2f});'
           f'  n_0 = {mp.nstr(n0, 4)} (T-block normalisation, m1_n0): log10 c_T* = {float(mp.log10(cstar/n0)):.2f}')
    pr(f'\nangular route: E0 = c_T sum w|a|^2 over the same families (J <= 1000) = {mp.nstr(E0, 6)};  L_g <= 70 E0/c_T;'
       f'  16 L_g <= tau_*^2 = 3.17e-4 needs c_T >= {mp.nstr(16*70*E0/mp.mpf("3.17e-4"), 4)}')
    open(os.path.join(OUT, f'm1c_physical_sum{SUF}.txt'), 'w').write('\n'.join(lines) + '\n')


if __name__ == '__main__':
    main()
