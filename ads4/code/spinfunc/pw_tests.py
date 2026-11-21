"""
pw_tests (HWIN_THIMBLE.md): point-wise identities of the per-wave integrands.
  1. scalar: exp(lamF + G + logJac) Q_J (l1_crystal, J) == exp(lamF + G_0 + logJac) P_J(t) (pw_core), up to a J-only sign;
  2. graviton: pw_grav pref (mode 0, f = P_J) == the Fourier-contraction pref (z2-patched EE builder / l9 C3H builder) at random
     points next to the crystal; and pref(+1) + pref(-1) == pref(0).
Usage: python3 pw_tests.py [J ...]  -> stdout, cache/pw_tests.txt
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
os.environ.setdefault('L1_DEG', '6')
import numpy as np
import n1_direct as D
import n1_pw as PW
import l1_crystal as Cr
import pw_core as C
import pw_grav as PG
import pw_wave as WV


def run(Js, pr):
    import z2_points as Z2
    Z2.patch()
    rng = np.random.default_rng(5)
    for J in Js:
        nu, Delta = 14.0, 1000.0 + J
        c, _ = Cr.saddle(nu, Delta)
        Wp = 0.02*(rng.normal(size=(4, 7)) + 1j*rng.normal(size=(4, 7)))
        Wp[0] = 0
        PW.set_batch(4)
        zl, lj = D.zmap([Wp[:, i].copy() for i in range(7)], c, pw=True)
        if J <= 24:
            lam, G, Q = PW.Cr.build(zl, nu, Delta, J, pref=True)
            lam0, G0 = PW.Cr.build(zl, nu, Delta, 0)
            zw, zz = C.zw_geom(zl)
            t = zw/np.sqrt(zz)
            r = np.exp(lam + G - lam0 - G0)*Q/WV.wave_derivs(J, t, 0, 0)[0]
            pr(f"J={J} scalar identity: ratio {np.round(r, 12)}")
        for kind in ('EE', 'C3H'):
            ev = C.WEv(nu, Delta, J, c, np.zeros(7, complex), np.eye(7, dtype=complex), 0, kind, nuf=nu)
            PW.set_batch(4)
            zl, lj = D.zmap([Wp[:, i].copy() for i in range(7)], c, pw=True)
            lamF, Gs, pref = PW.P9.build_spin(zl, nu, Delta, J, ev.mons, ev.sol, ev.solt, kind)
            ref = np.exp(lamF + Gs)[None, None]*pref
            out = {}
            for md in (0, 1, -1):
                ev.mode = md
                PW.set_batch(4)
                Pg, pr_ = PG.grav(ev, Wp.astype(complex))
                out[md] = np.exp(Pg - lj)[None, None]*pr_
            rel = np.abs(out[0]/ref - 1).max()
            rs = np.abs((out[1] + out[-1])/out[0] - 1).max()
            pr(f"J={J} {kind}: max rel |pref_pw(P_J)/pref_Fourier - 1| = {rel:.1e};  |(W+ + W-)/P - 1| = {rs:.1e};  "
               f"|W+|/|P| {np.abs(out[1][0, 0]/out[0][0, 0]).round(3)}")


if __name__ == '__main__':
    lines = []
    def pr(s=""):
        print(s, flush=True); lines.append(s)
    run([int(a) for a in sys.argv[1:]] or [8, 20, 40], pr)
    open(os.path.join(HERE, 'cache', 'pw_tests.txt'), 'a').write("\n".join(lines) + "\n")
