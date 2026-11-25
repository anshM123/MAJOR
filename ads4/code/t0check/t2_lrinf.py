"""
t0cert/t2_lrinf.py (T0_CERT.md §3, Part L): bhat > 2 Lam (beta > 2 asinh 1), every Lam >= Lam(DGMIN).

curv_k = Omega_LR,k(beta) - Omega_flat,k(bhat) = a1_k M1(beta)/Lam^2 + E_k - E_fl,k  with
  * Lemma LR-inf (LRINF_TIGHT / GAP_LARGESPIN U1, [P + certified]): Omega_LR,k = a1 M1/Lam^2 - a2/(2 bhat^3) + E_k, and with the
    envelope R^ anchored at Lam_a = Lam(DGMIN):  Lam^2 (2|E2| + |E3|) <= Lam_a^2 N^(beta; Lam_a) for every Lam >= Lam_a
    (gl_common.build_I(hat=True), E_up per beta-cell; tail_majorant for beta >= 64);
  * T0_CERT Lemma 5 (flat Regge remainder, only |J_n| <= 1 and int_0^inf J0 = 1): Omega_flat,k(b) = -a2/(2 b^3) + E_fl,
    |E_fl| <= e_k(b)/b^3, so Lam^2 |E_fl| <= e_k(bhat_min)/(Lam_a s^3), s = 2 sinh(beta/2), bhat_min = Lam_a s_a.
Sign-pattern curvature: K_g >= S_min M1 - (1+c)[Lam^2(2|E2| + |E3|)] - (1+c)[Lam^2(2|E_fl,2| + |E_fl,3|)] > 0,
   S_min = (1-c) 2 a1_2 - (1+c) a1_3 (worst g, sigma, u <= 1).
Output: cache/t2_lrinf_<DGMIN>_<c>.{txt,json}
"""
import os, sys, json, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fractions import Fraction as Fr
from flint import arb
import t0_common as T
from t0_common import A, up, lo, G

lines = []; pr = T.pr_to(lines)
c = T.C


def main():
    t0 = time.time()
    pr(__doc__.split('Output')[0])
    La = lo(T.LAMMIN)
    TAB = {k: G.build_I(G.W[k], La, hat=True) for k in (2, 3)}
    Smin = A((1 - c)*2*T.A1[2] - (1 + c)*T.A1[3])
    FE = {k: T.flat_E_consts(k) for k in (2, 3)}
    for k in (2, 3):
        pr(f'k={k}: flat-remainder constants |f1| = {float(FE[k]["f1"]):.6f}, |g(1)| = {float(FE[k]["g1"]):.6f}, '
           f'|g\'(1)| = {float(FE[k]["dg1"]):.6f}, int|p g\'\' + g\'| <= {FE[k]["Iq"].str(6)}')

    def e_of(k, b):
        d = FE[k]
        return A(d['f1'])*(1 + 2/b) + A(d['g1']) + (A(d['dg1']) + d['Iq'])/b

    bc = 2*arb(1).asinh()
    ok = True; worst = None; n = 0
    for a, b, m1 in T.CELLS:
        if not A(b) > bc:
            continue
        e2, _ = G.E_up(2, a, b, TAB, La); e3, _ = G.E_up(3, a, b, TAB, La)
        L2N = up(La**2*(2*e2 + e3))
        sa = 2*(A(a)/2).sinh()
        bmin = La*sa
        fl = up((2*e_of(2, bmin) + e_of(3, bmin))/(La*sa**3))
        K = Smin*m1 - A(1 + c)*L2N - A(1 + c)*fl
        n += 1
        if not K > 0:
            ok = False; pr(f'  FAIL beta-cell [{float(a)}, {float(b)}]: K >= {K.str(5)}')
        rel = (A(1 + c)*(L2N + fl))/(Smin*m1)
        if worst is None or rel.upper() > worst[0].upper():
            worst = (up(rel), float(a))
    pr(f'beta-cells in (2 asinh 1, 64]: {n}; all K_g > 0: {ok}; max over cells of (1+c)(LR-inf + flat remainders)/(S_min M1) <= '
       f'{worst[0].str(5)} (cell at beta = {worst[1]})')
    tl, ell = G.tail_majorant(TAB, La)
    B = arb(64); s = 2*(B/2).sinh()
    fr = (2*e_of(2, La*s) + e_of(3, La*s))*arb.pi()*B/(2*ell*La*s**2)
    Kt = Smin - A(1 + c)*tl*A(G.A1MIN) - A(1 + c)*fr
    okt = bool(Kt > 0)
    pr(f'tail beta >= 64: Lam^2(2|E2|+|E3|)/(A1min M1) <= {tl.str(5)}; flat remainder/M1-floor <= {up(fr).str(5)}; '
       f'K_g/M1 >= {Kt.str(6)} > 0: {okt}')
    verdict = ok and okt
    pr(f'\nVERDICT t2_lrinf (DGMIN = {T.DGMIN}, c = {c}): Part L {"PASS" if verdict else "FAIL"}  ({time.time()-t0:.0f}s)')
    tag = f'{T.DGMIN}_{str(c).replace("/", "_")}'
    json.dump(dict(DGMIN=str(T.DGMIN), c=str(c), pass_L=verdict, max_rel=float(worst[0].upper()), tail=float(Kt.lower())),
              open(os.path.join(T.CACHE, f't2_lrinf_{tag}.json'), 'w'), indent=1)
    open(os.path.join(T.CACHE, f't2_lrinf_{tag}.txt'), 'w').write('\n'.join(lines) + '\n')


if __name__ == '__main__':
    main()
