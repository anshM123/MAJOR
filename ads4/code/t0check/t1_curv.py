"""
t0cert/t1_curv.py (T0_CERT.md §3, Parts K and T): certified sign of the curvature part of every sign-pattern action, b-hat <= 2 Lam,
for every Lam >= Lam(DGMIN).

For a sign pattern g = (s2, s3), sigma = +-1, u in [0, 1] (affine in u -> u in {0, 1}), with (T0_CERT Lemma 2/3)
   Lam^2 sum_k (1 - c s_k) w_k curv_k = S_g M1(beta) + (1-c s2) 2 I2(bhat) + (1-c s3) u sigma I3(bhat) + eps * (...),
   S_g = (1-c s2) 2 a1_2 + (1-c s3) u sigma a1_3 > 0,  I_k = F_k/24,  |eps (...)| <= EPS0 [(1-c s2) 2 H2 + (1-c s3) u H3],
Part K: K_g := S_g M1min(beta_max) + ... - EPS0 (...) > 0 on bhat-cells of [B_POS, B_C] and on the tail [B_C, 2 Lam]
        (|F_k(b)| <= CT_k/b^2, M1 >= min over beta <= 2 asinh 1);
Part T: on bhat in [0, B_POS], T_gap := max(0, -min_g (K_g without the M1 term)), the threshold the flat action must beat
        (FF_g(c) >= T_gap/Lam^2 suffices; Lam^2 >= LAM2MIN).
Output: cache/t1_curv_<DGMIN>_<c>.{txt,json}
"""
import os, sys, json, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fractions import Fraction as Fr
from flint import arb
import t0_common as T
from t0_common import A, up, lo

lines = []; pr = T.pr_to(lines)
c = T.C
B_C = Fr(os.environ.get('T1_BC', '300'))
HW = Fr(os.environ.get('T1_H', '1/10'))


def main():
    t0 = time.time()
    pr(__doc__.split('Output')[0])
    assert 2*B_C <= 2*float(T.LAMMIN.lower()) or True
    assert A(B_C) < 2*T.LAMMIN, 'B_C must be <= 2 Lam_min'
    pr(f'DGMIN = {T.DGMIN}, Lam_min^2 = {T.LAM2MIN}, EPS0 = {T.EPS0.str(8)}, c = {c}, B_POS = {T.B_POS}, B_C = {B_C}, cell width {HW}')
    T.init_L()
    H = {}
    for k in (2, 3):
        H[k], tl, wh = T.H_sup(k)
        pr(f'k={k}: a1 = {float(T.A1[k]):.6f}; int|chi\'\'| <= {T.L0[k].str(6)}, int|chi\'\'|p^2 <= {T.L2[k].str(6)}; '
           f'H_k = sup hbar <= {H[k].str(6)} (attained near {wh}; asymptotic tail value {tl.str(5)})')
    a1 = {k: A(T.A1[k]) for k in (2, 3)}
    pats = [(s2, s3, sg, u) for s2 in (1, -1) for s3 in (1, -1) for sg in (1, -1) for u in (0, 1)]

    def Kvals(I2, I3, M):
        out = []
        for s2, s3, sg, u in pats:
            c2 = 1 - c*s2; c3 = (1 - c*s3)*u*sg
            S = A(c2)*2*a1[2] + A(c3)*a1[3]
            assert S > 0
            eps = T.EPS0*(A(c2)*2*H[2] + A(abs(c3))*H[3])
            out.append(((s2, s3, sg, u), S*M + A(c2)*2*I2 + A(c3)*I3 - eps))
        return out

    # ---------------- Part T: [0, B_POS]
    Tgap = arb(0); Targ = None; minK0 = None
    b = Fr(0)
    while b < T.B_POS:
        bb = min(T.B_POS, b + HW)
        I2 = T.F_cell(2, b, bb)/24; I3 = T.F_cell(3, b, bb)/24
        for pat, v in Kvals(I2, I3, arb(0)):
            if (-v).upper() > Tgap.upper():
                Tgap = up(-v); Targ = (float(b), pat)
        b = bb
    Tabs = up(Tgap/A(T.LAM2MIN))
    pr(f'\nPart T (bhat in [0, {T.B_POS}]): T_gap = max(0, -K_g w/o M1) <= {Tgap.str(6)} at {Targ}; '
       f'threshold T_abs = T_gap/Lam_min^2 <= {Tabs.str(6)}')
    # ---------------- Part K: [B_POS, B_C]
    okK = True; worst = None; nc = 0
    b = Fr(T.B_POS)
    rows = []
    while b < B_C:
        bb = min(B_C, b + HW)
        I2 = T.F_cell(2, b, bb)/24; I3 = T.F_cell(3, b, bb)/24
        bmax = 2*(A(bb)/(2*lo(T.LAMMIN))).asinh()
        M = T.M1min_upto(bmax)
        for pat, v in Kvals(I2, I3, M):
            if not v > 0:
                okK = False
                pr(f'  FAIL cell [{float(b)}, {float(bb)}] {pat}: K >= {v.str(5)}')
            if worst is None or v.lower() < worst[0].lower():
                worst = (v, float(b), pat)
        nc += 1
        if int(b) != int(bb) and int(bb) in (20, 30, 50, 100, 200, 300):
            rows.append((float(bb), float(worst[0].lower())))
        b = bb
    pr(f'Part K cells: {nc} cells of width {HW} on [{T.B_POS}, {B_C}], all K_g > 0: {okK}; min lower bound {worst[0].str(6)} at '
       f'bhat = {worst[1]}, pattern (s2, s3, sigma, u) = {worst[2]}')
    # ---------------- tail [B_C, 2 Lam]
    CT = {k: T.tail_const(k) for k in (2, 3)}
    Itail = {k: up(CT[k]/(24*A(B_C)**2)) for k in (2, 3)}
    Mt = T.M1min_upto(2*arb(1).asinh())
    Smin = (1 - c)*2*T.A1[2] - (1 + c)*T.A1[3]
    Kt = A(Smin)*Mt - A(1 + c)*(2*Itail[2] + Itail[3]) - T.EPS0*A(1 + c)*(2*H[2] + H[3])
    okT = bool(Kt > 0)
    pr(f'Part K tail bhat in [{B_C}, 2 Lam]: |I_k| <= CT_k/(24 B_C^2) = {Itail[2].str(4)} / {Itail[3].str(4)} (CT = {CT[2].str(5)} / '
       f'{CT[3].str(5)}); S_min = {float(Smin):.4f}; M1 >= {Mt.str(6)} on beta <= 2 asinh 1; K >= {Kt.str(6)} > 0: {okT}')
    verdict = okK and okT
    pr(f'\nVERDICT t1_curv (DGMIN = {T.DGMIN}, c = {c}): Part K {"PASS" if verdict else "FAIL"}; Part T threshold T_abs <= {Tabs.str(8)}'
       f'  ({time.time()-t0:.0f}s)')
    tag = f'{T.DGMIN}_{str(c).replace("/", "_")}'
    json.dump(dict(DGMIN=str(T.DGMIN), c=str(c), B_POS=str(T.B_POS), B_C=str(B_C), Tgap_upper=Tgap.upper().str(30, radius=False) if hasattr(Tgap.upper(), 'str') else str(Tgap.upper()),
                   Tabs_upper=str(Tabs.upper()), Kmin=float(worst[0].lower()), Kmin_at=worst[1], Ktail=float(Kt.lower()), H2=float(H[2].upper()),
                   H3=float(H[3].upper()), pass_K=verdict),
              open(os.path.join(T.CACHE, f't1_curv_{tag}.json'), 'w'), indent=1)
    open(os.path.join(T.CACHE, f't1_curv_{tag}.txt'), 'w').write('\n'.join(lines) + '\n')


if __name__ == '__main__':
    main()
