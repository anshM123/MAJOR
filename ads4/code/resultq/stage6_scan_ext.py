"""Result Q, Stage 6 (RESULT_Q_STAGE6.md §1.4): extension of region A (HYB) to LARGE J at small u.

Reason: the flat Stage-0 action of both functionals dips at u ~ 1e-4, bhat ~ 55-85 (delta = 0.05: +0.025 at m = 90.6 Lam, J = 2689; delta = 0.2:
+0.185 at J = 2980), a band the Stage-0 LP never constrained (its grid stopped at J = 2500).  There the scaling function 2 int psi2 J0 is tiny
and the O(u) terms dominate, so the finite-J (flat-exact minus Bessel) corrections that LR drops could be of leading size; HYB is therefore
used for J > 3000 too.  (Outcome, RESULT_Q_STAGE6.md §1.6: at J > 3000 HYB and LR agree to <= 0.002 -- the dip is carried by both.)
Region A': m/Lam = x in geomspace(20, 600, 100), every integer J in (3000, min(100 x, 40000)] (bhat <= 200), twist gap imposed; HYB and LR on the
same states.  Region B': LR (Stage 5 beta_quantities, dense small-beta grid) for J > 3000 with bhat > 200 or m/Lam > 600 (the complement).
Output: cache/stage6_scan_ext_<delta>.txt
"""
import os, sys, io, time, pickle
os.environ.setdefault('OMP_NUM_THREADS', '2')
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import common_q as CQ
import stage6_model as M
import stage5_positivity as S5P
from stage6_scan import BMAX, Best, fmt

CACHE = M.CACHE
buf = io.StringIO()
XE = np.geomspace(20, 600, 100)


def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); buf.write(s + "\n")


def flat_ext(name):
    fn = os.path.join(CACHE, 'stage6_flatext_%s.pkl' % name.split()[0])
    if os.path.exists(fn):
        return pickle.load(open(fn, 'rb'))
    F, c2, c3 = M.load(name)
    out = {'++': [], '+-': []}
    for x in XE:
        u = 1/x**2; Jm = int(min(100*x, 40000))
        for h in ('++', '+-'):
            Js, B2, B3 = CQ.parts(h, u, Jm, [c2], [c3])
            sel = Js > 3000
            if h == '++':
                sel &= Js % 2 == 0
            out[h].append((u, Js[sel], B2[sel, 0], B3[sel, 0]))
        print('  flat ext x = %.1f (J <= %d)' % (x, Jm), flush=True)
    pickle.dump(out, open(fn, 'wb'))
    return out


def regionB(F, name, Lam, Dg, best):
    a, c = F['a'], F['c']
    F_, c2, c3 = M.load(name)
    T = [S5P.CosT(c2), S5P.CosT(np.convolve(c2, [0, 0, 1])), S5P.CosT(c3)]
    betas = np.unique(np.concatenate([np.linspace(1.0, 12.0, 221), np.linspace(12.0, 30.0, 91),
                                      2*np.arcsinh(np.geomspace(20.0, 20000.0, 1200)/(2*Lam))]))
    us = np.concatenate([np.geomspace(1e-9, 0.05, 200), np.linspace(0.0505, 1.0, 100)])
    for beta in betas:
        I20, I22, I3, P0, D, H5 = S5P.beta_quantities(beta, Lam, T)
        bh = 2*Lam*np.sinh(beta/2)
        for h in ('++', '+-'):
            s = M.SIG[h]; c4 = M.C4[h]
            J = 0.5*(bh/np.sqrt(us) - 1)
            ok = (J > 3000) & (us <= M.u_gap(J, Dg)) & ((bh > 200) | (us < 1/600.**2))
            if not np.any(ok):
                continue
            u = us[ok]
            B2 = 2*I20 - u*I22; B3 = s*u*I3
            B4 = 2*u*u*P0; dB4 = -2*u*u*Lam**2*D + c4*u**3*P0; B5 = s*u**3*P0
            FW = -a*(dB4 - B5); FF = B2 + B3 + FW - c*B4
            N = np.abs(B2) + np.abs(B3) + np.abs(FW) + c*np.abs(B4)
            best.add(('LRB\'', 'J>3000 & (bhat>200 or m>600Lam)'), FF/N, dict(h=h, J=J[ok], x=1/np.sqrt(u), b=np.full(len(u), bh)))


if __name__ == '__main__':
    deltas = [float(a) for a in sys.argv[1:]] or [0.05, 0.2]
    for delta in deltas:
        name = M.FUNCS[delta]; F, c2, c3 = M.load(name)
        t0 = time.time()
        fl = flat_ext(name)
        pr(f"\n######## delta = {delta}: '{name}' region A' (J > 3000, bhat <= 200, m/Lam in [20, 600]) [{time.time() - t0:.0f}s]")
        for Dg in (100, 300, 1000):
            Lam = M.lam(Dg); tab = M.Tables(name, Lam, BMAX[Dg])
            best = Best(); prof = []
            for h in ('++', '+-'):
                for (u, Js, B2f, B3f) in fl[h]:
                    sel = Js <= M.jmax_gap(u, Dg)
                    if not np.any(sel):
                        continue
                    J = Js[sel].astype(float); x = 1/np.sqrt(u)
                    FF, N, Nrm, P = M.hyb(F, h, u, J, B2f[sel], B3f[sel], tab, Lam)
                    FFl, Nl, _ = M.lr(F, h, u, J, tab, Lam)
                    info = dict(h=h, J=J, x=np.full(len(J), x), b=P['bh'])
                    best.add(('HYB', "A' all"), FF/N, info); best.add(('LR', "A' all"), FFl/Nl, info)
                    i = int(np.argmin(FF/N)); prof.append((x, h, (FF/N)[i], J[i], P['bh'][i], (FFl/Nl)[i], np.max(np.abs(FF/N - FFl/Nl))))
            regionB(F, name, Lam, Dg, best)
            pr(f"=== Dg = {Dg}:")
            for (mod, lab), t in sorted(best.d.items()):
                pr(f"   {mod:5s} {lab:36s}: {fmt(t)}")
            pr("   profile (per m/Lam: channel, HYB min, J, bhat; LR at that state; max |HYB - LR| at this m/Lam):")
            for (x, h, v, J, b, vl, dmax) in prof:
                if h == '+-' or True:
                    pass
            P = sorted(prof, key=lambda r: r[2])[:6]
            for (x, h, v, J, b, vl, dmax) in sorted(P):
                pr(f"     m/Lam {x:7.2f} {h}: HYB {v:+.4f} (J {J:.0f}, bhat {b:.1f}); LR there {vl:+.4f}; max|HYB-LR| {dmax:.3f}")
        open(os.path.join(CACHE, f'stage6_scan_ext_{delta}.txt'), 'w').write(buf.getvalue())
