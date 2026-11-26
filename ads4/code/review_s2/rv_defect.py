"""rv_defect.py [DV, adversarial]: Lemma F's per-step quantities over a box, tested at random interior/corner points p.
Re-implements tpf3.lemmaF3_log's per-step bounds (same formulas, same frames via tpf3.build) and records, for each step i and column c,
  bnd_delta[i][c] >= sup_B |delta_{c,i}(p)|,  bnd_R[i][c] >= sup_B |R_{c,i}(p)|,  e_i >= sup_B ||E_i(p)||_inf,
then evaluates the exact delta_{c,i}(p), |R_{c,i}(p)| and E_i(p) = Q_{i+1}(p) M_i(p) F_i(p) (mpmath, 50 digits) at sample points and
reports the worst ratio exact/bound (must be <= 1)."""
import sys, math, random
sys.path.insert(0, '/home/user/MAJOR/ads4/code/review_s2/ws/code/analytic2')
import numpy as np
import mpmath as mp
from flint import arb, acb, ctx, acb_mat
import tpf3 as TF
from tpframes import up, arb_iv
mp.mp.dps = 50
NP = 4


def bounds(s_lo, s_hi, t_lo, t_hi, a_lo, a_hi, i_lo, i_hi, K, i_s):
    ctx.prec = 96
    p0 = [(t_lo + t_hi)/2, (a_lo + a_hi)/2, (i_lo + i_hi)/2, (s_lo + s_hi)/2]
    D = [(t_hi - t_lo)/2, (a_hi - a_lo)/2, (i_hi - i_lo)/2, (s_hi - s_lo)/2]
    Da = [arb(x) for x in D]
    box = [arb_iv(t_lo, t_hi), arb_iv(a_lo, a_hi), arb_iv(i_lo, i_hi), arb_iv(s_lo, s_hi)]
    alB = acb(box[3], box[0]); al0 = acb(arb(p0[3]), arb(p0[0]))
    W0 = []; DW0 = []; D2B = []; WB = []
    for i in range(K):
        m0 = arb(p0[1]) + (K - i); mB = box[1] + (K - i)
        ed0 = TF.eps_derivs(m0, arb(p0[2])); edB = TF.eps_derivs(mB, box[2])
        d10, _ = TF.w_derivs(ed0, al0); _, d2B = TF.w_derivs(edB, alB)
        W0.append(ed0['e']*al0/2); DW0.append(d10); D2B.append(d2B); WB.append(edB['e']*alB/2)
    cols = TF.build(p0[3], p0[0], p0[1], p0[2], K, i_s)
    _c = TF._c
    def Rball(cg):
        R0, gg = cg; v = _c(R0)
        for l in range(NP): v = v + _c(gg[l])*acb(arb(0, D[l]))
        return v
    def Rabs(cg):
        R0, gg = cg; R0a = _c(R0); lin = arb(0); tot1 = arb(0)
        for l in range(NP):
            lin += abs((R0a.conjugate()*_c(gg[l])).real)*Da[l]; tot1 += abs(_c(gg[l]))*Da[l]
        return (abs(R0a)**2 + 2*lin + tot1**2).sqrt()
    FB = [acb_mat([[Rball(cols[i][0]), Rball(cols[i][1])], [acb(1), acb(1)]]) for i in range(K + 1)]
    QB = [f.inv() for f in FB]
    bd = []; bR = []; eN = []
    for i in range(K):
        if i + 1 == i_s:
            E = QB[i + 1]*acb_mat([[2*WB[i], acb(1)], [acb(1), acb(0)]])*FB[i]
            eN.append(max(up(abs(E[r, 0]) + abs(E[r, 1])) for r in range(2))); bd.append(None); bR.append(None); continue
        rowd = []; rowR = []; rsum = arb(0); Rmax = 0.0
        for c in range(2):
            (R0, gl), (R1, gl1) = cols[i][c], cols[i + 1][c]
            R0a, R1a = _c(R0), _c(R1)
            bound = abs(2*W0[i] + 1/R0a - R1a)
            for l in range(NP):
                bound += abs(2*DW0[i][l] - _c(gl[l])/R0a**2 - _c(gl1[l]))*Da[l]
            RB = Rball(cols[i][c])
            for l in range(NP):
                for mm in range(NP):
                    if D[l] == 0 or D[mm] == 0: continue
                    bound += abs(2*D2B[i][l][mm] + 2*_c(gl[l])*_c(gl[mm])/RB**3)*Da[l]*Da[mm]/2
            Ra = Rabs(cols[i][c]); rowd.append(up(bound)); rowR.append(up(Ra)); rsum += Ra*bound; Rmax = max(Rmax, up(Ra))
        qrow = max(up(abs(QB[i + 1][0, 0])), up(abs(QB[i + 1][1, 0])))
        eN.append(Rmax + qrow*up(rsum)); bd.append(rowd); bR.append(rowR)
    return p0, cols, bd, bR, eN


def exact_at(p, p0, cols, K, i_s):
    tau, a, io, sig = [mp.mpf(x) for x in p]
    al = mp.mpc(sig, tau)
    dp = [mp.mpf(p[l]) - mp.mpf(p0[l]) for l in range(NP)]
    def R(cg):
        R0, gg = cg
        return mp.mpc(R0) + sum(mp.mpc(gg[l])*dp[l] for l in range(NP))
    def w(i):
        m = a + K - i; B = 2 - (m - 1)*io
        return (2/m - 2*io/B)*al/2
    out = []
    for i in range(K):
        Fi = mp.matrix([[R(cols[i][0]), R(cols[i][1])], [1, 1]])
        Fn = mp.matrix([[R(cols[i + 1][0]), R(cols[i + 1][1])], [1, 1]])
        Mi = mp.matrix([[2*w(i), 1], [1, 0]])
        E = mp.inverse(Fn)*Mi*Fi
        nE = max(abs(E[r, 0]) + abs(E[r, 1]) for r in range(2))
        dl = [abs(2*w(i) + 1/R(cols[i][c]) - R(cols[i + 1][c])) for c in range(2)]
        Rs = [abs(R(cols[i][c])) for c in range(2)]
        out.append((nE, dl, Rs))
    return out


if __name__ == '__main__':
    random.seed(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
    K = 200
    boxes = [(2.5, 2.5, 1031.0, 1032.1, 998.0, 999.0), (3.5, 3.5, 1003.0, 1004.1, 998.0, 999.0), (2.5, 2.5, 1040.0, 1041.1, 1020.0, 1021.0),
             (4.5, 4.5, 1020.0, 1021.1, 1005.0, 1006.0), (2.5, 2.5, 1045.0, 1046.1, 1030.0, 1031.0)]
    for (sl, sh, tl, th, al_, ah) in boxes:
        for (La, Lb) in ((3.65, 3.75), (4.0, 4.5)):
            import rt2_kernel as RK
            il, ih = RK.iota_range(al_, ah, math.exp(-Lb), math.exp(-La))
            for i_s in (K, K - 20):
                try:
                    p0, cols, bd, bR, eN = bounds(sl, sh, tl, th, al_, ah, il, ih, K, i_s)
                except ZeroDivisionError:
                    print('singular frames', (tl, al_, La, i_s)); continue
                worst = dict(E=0.0, d=0.0, R=0.0)
                pts = [[tl, al_, il, sl], [th, ah, ih, sh], [tl, ah, il, sh], [th, al_, ih, sl]] + \
                      [[random.uniform(tl, th), random.uniform(al_, ah), random.uniform(il, ih), random.uniform(sl, sh)] for _ in range(4)]
                for p in pts:
                    ex = exact_at(p, p0, cols, K, i_s)
                    for i in range(K):
                        nE, dl, Rs = ex[i]
                        worst['E'] = max(worst['E'], float(nE/eN[i]))
                        if bd[i] is not None:
                            for c in range(2):
                                worst['d'] = max(worst['d'], float(dl[c]/bd[i][c]) if bd[i][c] > 0 else 0.0)
                                worst['R'] = max(worst['R'], float(Rs[c]/bR[i][c]))
                print(f'box tau[{tl},{th}] a[{al_},{ah}] sig {sl} L[{La},{Lb}] i_s={i_s}: max exact/bound - 1: ||E|| {worst["E"]-1:+.3e}, '
                      f'|delta| (ratio) {worst["d"]:.4f}, |R| {worst["R"]-1:+.3e}   (log sum e_i = {sum(math.log(x) for x in eN):.3f})', flush=True)
