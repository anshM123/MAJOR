"""PARITY_ODD_DZ §5: the parity-odd part of the light side.  O = T (J = 2, Delta = 3, s = 5/2, h = 1/2, tau = 1):
(1) the parity-odd <TTT> structure (unique at s = 5/2): zeroth-moment components m = -2..2 and first-moment components,
    source helicity +2 and -2.  Expected: M_2 = M_{+-1} = M_0 = 0, M_{-2} != 0 (pure flip, Delta h = 4); first moments: no
    diagonal (m = 2) element.
(2) the generic J = 2 parity-odd family at Delta -> 3: m = -2 component finite and nonzero (h = 1/2 is NOT a zero of the
    odd amplitude, 1/Gamma((7-Delta)/2) = 1 there), m = 0, -1 components -> 0 (conserved-current selection), and the
    family's m = -2 amplitude is proportional to the <TTT> odd one.
Output: cache/oz3_out.txt"""
import mpmath as mp
from common import *
from oz_amps import amps
from oz2_zeros import G_expr
mp.mp.dps = 30
out = []
pr = lambda x: (print(x, flush=True), out.append(str(x)))
sig = sp.Rational(5, 2)
monos, basis = kl.conserved(2, 'odd', sval='5/2')
pr(f'(1) parity-odd <TTT> structures at s = 5/2: {len(basis)}')
res = {}
for src in ('plus', 'minus'):
    kin = pl.Kin(src)
    e0 = [pl.scri_mono(m, kin) for m in monos]
    for mo, ex in (('zero', e0), ('mid', [moment_expr(e, 'odd', 'mid') for e in e0]),
                   ('orig', [moment_expr(e, 'odd', 'orig') for e in e0])):
        K = pl.Kbounds(ex, 'odd')
        polys = [pl.amplitude(e, kin, 'odd', sig=sig, K=K)[0] for e in ex]
        tot = kl.combine(polys, [sp.Rational(c) for c in basis[0]])
        comp = {m: kl.poly_to_sym(v).subs(S, sig) for m, v in pl.tcoeffs(tot, 2).items()}   # moment_expr keeps a symbolic s
        G = complex(sp.N(G_expr('odd', K).subs(S, sig), 25))
        vals = {m: complex(sp.N(comp.get(m, 0), 25))*G for m in range(-3, 4)}
        res[(mo, src)] = vals
        pr(f'    {mo:4s} moment, source {src:5s}: ' + ', '.join(f'm={m:+d}: {mp.nstr(mp.mpc(vals[m]), 8)}' for m in range(-3, 4)))
z = res[('zero', 'plus')]
pr(f'    => zeroth moment, source +2: only m = -2 nonzero (pure flip): {all(abs(z[m]) < 1e-25 for m in (-1, 0, 1, 2, 3)) and abs(z[-2]) > 1e-6}')
pr(f'    => first moments, source +2: diagonal m = +2 element = {mp.nstr(abs(mp.mpc(res[("orig","plus")][2])), 3)} (orig), '
   f'{mp.nstr(abs(mp.mpc(res[("mid","plus")][2])), 3)} (mid)')
# (2) generic J = 2 odd family at Delta -> 3
r = amps('odd', 2)
K, P = r[('zero', 'plus')]
G = G_expr('odd', K)
pr('(2) generic J = 2 parity-odd family, source +2, near Delta = 3 (s = 5/2 + eps/2):')
for Dv in (sp.Rational(3) + sp.Rational(1, 10**6), sp.Rational(3) + sp.Rational(1, 10**9)):
    sv = (Dv + 2)/2
    v = {m: complex(sp.N((G*P.get(m, 0)).subs(S, sv), 25)) for m in (-2, -1, 0)}
    pr(f'    Delta = 3 + {mp.nstr(mp.mpf(Dv.p)/Dv.q - 3, 3)}: ' + ', '.join(f'M_{m:+d} = {mp.nstr(mp.mpc(v[m]), 10)}' for m in (-2, -1, 0)))
lim = sp.limit(sp.gammasimp(G*P[-2]), S, sig)
pr(f'    exact limit Delta -> 3 of M_-2(family) = {sp.nsimplify(lim)} ; <TTT>_odd M_-2 = {mp.nstr(mp.mpc(z[-2]), 15)}')
pr(f'    ratio <TTT>_odd / family limit = {mp.nstr(mp.mpc(z[-2])/complex(sp.N(lim, 25)), 15)}  (a constant normalisation)')
write('oz3_out.txt', out)
