"""
np2_onepart.py -- the regulated null-plane identity in one-particle wavepackets of the d = 3 free massless field
(scalar or Dirac: the one-body operators are identical, np1_kernels.py).

For real profiles F (Re / Im of f_m, optionally smoothly cut off at |y| ~ R) and A = A(F):
    G(t)   = <psi,[A, A_t] psi>/i  computed by exact free evolution  A_t = e^{iHt} A e^{-iHt}  (FFT, per K-slice)
    I_eta  = -int rho_eta'(t) G(t) dt,  rho_eta(t) = e^{-t^2/2eta^2}/(sqrt(2pi) eta)   (hat rho >= 0, even)
and compares lim I_eta with
    W      = (1/2) <psi, A(F'^2) psi>              (the Ward / f-sum prediction)
    2||H^{1/2} A psi||^2 - 2 Re <A psi, A H psi>     (the quadratic-form side)
Also: the spectral form of the author's regulator theorem, the m-scaling, a smooth vs sharp cutoff.
Output: cache/np2_onepart.txt
"""
import numpy as np, os, sys, time
sys.path.insert(0, os.path.dirname(__file__))
from lfcore import Slice, gauss_legendre, profile, fm, I_eta_from_G

out = []
def P(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); out.append(s)

def cov_gauss(p0, sig, y0=0.0):
    def chi_q(K, q):
        h = q**2/(4*K)
        px, py = h - K, q
        return np.exp(-((px - p0[0])**2 + (py - p0[1])**2)/(2*sig**2))/np.sqrt(2*K)*np.exp(-1j*q*y0)
    return chi_q

def build_state(chi_q, nK=48, dy=0.02, L0=120.0, tailtol=1e-26):
    # K-support scan
    Ks = np.geomspace(1e-7, 5, 4000)
    qs = np.linspace(-4, 4, 4001)
    mx = np.array([np.max(np.abs(chi_q(K, qs))**2) for K in Ks])
    keep = Ks[mx > 1e-30*mx.max()]
    Klo, Khi = keep.min(), keep.max()
    Kn, wK = gauss_legendre(Klo, Khi, nK)
    slices = []
    for K, w in zip(Kn, wK):
        L = L0
        while True:
            s = Slice(K, w, L, dy)
            c = s.from_q(chi_q(K, s.q))
            tail = np.sum(np.abs(c[np.abs(s.y) > 0.75*L])**2)/np.sum(np.abs(c)**2)
            if tail < tailtol or L > 5000:
                break
            L *= 2
        s.chi = c
        slices.append(s)
    nrm = sum(s.w*np.real(s.ip(s.chi, s.chi)) for s in slices)
    for s in slices:
        s.chi = s.chi/np.sqrt(nrm)
    return slices, (Klo, Khi)

def quantities(slices, Ffun, etas=(), Gt=None):
    """Ffun(y) -> (F, F') real.  Returns dict."""
    W = Hn = Rc = PK = 0.0
    spec_parts = []
    for s in slices:
        F, Fp = Ffun(s.y)
        v = s.K*F*s.chi
        W += s.w*0.5*s.K*np.real(np.sum(Fp**2*np.abs(s.chi)**2))*s.dy
        Hn += s.w*2*s.hform(v)
        Rc += s.w*2*np.real(s.ip(v, s.K*F*s.H(s.chi)))
        PK += s.w*s.K*np.real(s.ip(s.chi, s.chi))
        spec_parts.append((s.w, s.h, np.abs(np.fft.fft(v))**2*s.dy/s.N))
    res = dict(W=W, Hn=Hn, Rc=Rc, PK=PK)
    def G(ts):
        tot = np.zeros(len(ts))
        for s in slices:
            F, _ = Ffun(s.y)
            a = s.K*F
            ch = max(1, int(2**21//s.N))
            for i0 in range(0, len(ts), ch):
                tt = ts[i0:i0 + ch]
                Uc = s.evolve(s.chi, tt)
                Uac = s.evolve(a*s.chi, tt)
                tot[i0:i0 + ch] += s.w*2*np.imag(np.sum(np.conj(Uac)*a[None, :]*Uc, axis=1))*s.dy
        return tot
    res['G'] = G
    res['I'] = {eta: I_eta_from_G(G, eta) for eta in etas}
    res['Ispec'] = {eta: sum(w*2*np.sum(h*np.exp(-0.5*(eta*h)**2)*pw) for (w, h, pw) in spec_parts) - Rc for eta in etas}
    return res

if __name__ == '__main__':
    t0 = time.time()
    etas = (0.4, 0.2, 0.1, 0.05, 0.025, 0.0125)
    states = {
        'S1 backward packet p0=(-0.6,0.3) sig=0.12': cov_gauss((-0.6, 0.3), 0.12),
        'S2 packet p0=(0.2,0.5) sig=0.10 (reaches small K)': cov_gauss((0.2, 0.5), 0.10),
        'S3 = S1 displaced to y0=1.5': cov_gauss((-0.6, 0.3), 0.12, y0=1.5),
    }
    built = {}
    for name, cq in states.items():
        sl, (Klo, Khi) = build_state(cq)
        built[name] = sl
        P(f"\n##### {name}:  K in [{Klo:.3g}, {Khi:.3g}], {len(sl)} K-nodes, N_y in [{min(s.N for s in sl)}, {max(s.N for s in sl)}]")
        P("  profile        lim-check:  I_eta for eta = " + ", ".join(str(e) for e in etas))
        for m in (0, 1, 2, 4):
            for part in ('re', 'im'):
                if m == 0 and part == 'im':
                    continue
                r = quantities(sl, lambda y: profile(y, m, part), etas)
                Is = [r['I'][e] for e in etas]
                P(f"  {part} f_{m}:  I_eta = " + "  ".join(f"{x:.10f}" for x in Is))
                P(f"           W=(1/2)<A(F'^2)> = {r['W']:.10f}   2||H^1/2 A psi||^2 = {r['Hn']:.10f}   2Re<Apsi,AHpsi> = {r['Rc']:.10f}")
                P(f"           I(eta_min) - W = {Is[-1] - r['W']:.2e}   [2||H^1/2Apsi||^2 - 2Re<..>] - W = {r['Hn'] - r['Rc'] - r['W']:.2e}"
                  f"   (I - W)/eta^2 at last two etas: {(Is[-2]-r['W'])/etas[-2]**2:.4f}, {(Is[-1]-r['W'])/etas[-1]**2:.4f}")
                P(f"           spectral form 2 int lam rhohat(eta lam) dnu_v - 2Re<v,AHpsi>: " +
                  "  ".join(f"{r['Ispec'][e]:.10f}" for e in etas))
        P(f"  <A(1)> = <P_K> = {r['PK']:.10f}")

    # ---------------- m-scaling (complex A_m and its Re / Im parts), state S1 and S3
    P("\n##### m-scaling of ||H^{1/2} A_m psi||^2 (A_m = A(f_m), complex), and of W_re + W_im = 2<A(y^2 + m^2)>")
    for name in ('S1 backward packet p0=(-0.6,0.3) sig=0.12', 'S3 = S1 displaced to y0=1.5'):
        sl = built[name]
        PK = sum(s.w*s.K*np.real(s.ip(s.chi, s.chi)) for s in sl)
        Ay2 = sum(s.w*s.K*np.real(np.sum(s.y**2*np.abs(s.chi)**2))*s.dy for s in sl)
        P(f"  [{name}]  <A(1)> = {PK:.8f},  <A(y^2)> = {Ay2:.8f}")
        P("    m    ||H^1/2 A_m psi||^2   m^2<A(1)>     remainder      W_re+W_im      2<A(y^2+m^2)>   ||H^1/2A_re||^2+||H^1/2A_im||^2")
        rows = []
        for m in (0, 1, 2, 3, 4, 6, 8, 12, 16):
            Hc = sum(s.w*s.hform(s.K*fm(s.y, m)[0]*s.chi) for s in sl)
            rr = quantities(sl, lambda y: profile(y, m, 're'))
            ri = quantities(sl, lambda y: profile(y, m, 'im'))
            rows.append((m, Hc))
            P(f"  {m:3d}   {Hc:16.8f}   {m*m*PK:12.6f}   {Hc - m*m*PK:12.6f}   {rr['W'] + ri['W']:12.8f}   {2*(Ay2 + m*m*PK):12.8f}"
              f"   {(rr['Hn'] + ri['Hn'])/2:12.8f}")
        ms = np.array([r[0] for r in rows], float); Hs = np.array([r[1] for r in rows])
        cfit = np.polyfit(ms, Hs, 2)
        P(f"    quadratic fit ||H^1/2 A_m psi||^2 = {cfit[0]:.8f} m^2 + {cfit[1]:.8f} m + {cfit[2]:.8f}  (exact one-particle law:"
          f" m^2<A(1)> + m<(1+y^2)K j> + const; max fit residual {np.max(np.abs(np.polyval(cfit, ms) - Hs)):.1e})")

    # ---------------- cutoff profiles: smooth (C^inf) vs sharp
    P("\n##### cutoff profiles on S2 (the state reaching small K), m = 2, Re part")
    sl = built['S2 packet p0=(0.2,0.5) sig=0.10 (reaches small K)']
    r_full = quantities(sl, lambda y: profile(y, 2, 're'), etas)
    P(f"  uncut: W = {r_full['W']:.10f}, I(eta=0.0125) = {r_full['I'][0.0125]:.10f}")
    for R in (1.0, 2.0, 4.0, 8.0):
        r = quantities(sl, lambda y: profile(y, 2, 're', R=R), etas)
        P(f"  smooth cutoff R={R:4.1f}: W = {r['W']:.10f}, I_eta = " + "  ".join(f"{r['I'][e]:.8f}" for e in etas) +
          f",  I(eta_min)-W = {r['I'][etas[-1]] - r['W']:.1e}")
    etas2 = (0.1, 0.05, 0.025, 0.0125, 0.00625, 0.003125)
    def sharp(y, R=2.0):
        f, fp = profile(y, 2, 're')
        ind = (np.abs(y) < R).astype(float)
        return f*ind, fp*ind      # F' here omits the delta functions: (F')^2 is not defined for a jump
    r = quantities(sl, sharp, etas2)
    P("  SHARP cutoff R=2 (A psi has a jump, A psi not in Dom(H^1/2)): I_eta = " +
      "  ".join(f"{r['I'][e]:.5f}" for e in etas2))
    Iv = np.array([r['I'][e] for e in etas2])
    P("    successive ratios I(eta/2)/I(eta) = " + "  ".join(f"{x:.4f}" for x in Iv[1:]/Iv[:-1]) +
      "   (eta^{-1/2} growth -> ratio sqrt2 = 1.4142): liminf I_eta = infinity, as the regulator theorem requires")
    P(f"\n[time {time.time() - t0:.1f}s]")
    open(os.path.join(os.path.dirname(__file__), 'cache', 'np2_onepart.txt'), 'w').write('\n'.join(out) + '\n')
