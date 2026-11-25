"""
np4_k0.py -- the author's K^0 estimate (ANGULAR_ROUTE §17-18) in the d = 3 free scalar, in the ORIGINAL (celestial) frame.

    psi_m = T(F_m)|0>,  F^_m(P) = exp(-|P - P0|^2/(2 delta_m^2)),  P0 = (1,0,0),  delta_m = eps/(1+|m|)   (helicity +2)
    E_m   = int dtheta e^{i m theta} E(theta)  -> one-body multiplication by |p| e^{i m theta_p}
    K^0   = -|p| Laplacian_p  on one-particle wavefunctions with measure d^2p/|p|   (P^0 = |p|; (P^0+K^0)/2 >= Delta_phi = 1/2,
            ground state e^{-|p|}; [K^0, P^0] = -(2 p.grad + 1) = -2i D).

Checks, for w = E_m psi_m:
  (i)   ||(K^0)^{1/2} E_m psi_m||^2 / ||psi_m||^2 = O(m^2)  (author's claim) -- value of the constant;
  (ii)  family decomposition: at fixed total P the two-particle states with CM angular momentum j span the conserved-current
        family J_|j| (Delta = |j|+1).  Target fraction ||w_{m+2}||^2/||w||^2 (author: stays >= 1/2 with Doppler error O(eps));
  (iii) coercivity chain <P^0>_w <K^0>_w >= (sum_O Delta_O ||w_O||^2)^2 ... and sum_O Delta_O^2 |a_O|^2 ~ (m+3)^2 (exactly
        saturating the O(m^2) family bound, since free event shapes are delta functions: |a_N|^2 = 1 for every even N).
Measure: d^2p1 d^2p2/(|p1||p2|) = d^3P dtheta*/sqrt(s)  (CM angle theta*).
Output: cache/np4_k0.txt
"""
import numpy as np, os, time
out = []
def P(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); out.append(s)

def run(m, eps=0.15, nG=22, nth=None):
    d = eps/(1 + abs(m))
    nth = nth or max(64, 8*(abs(m) + 4))
    x, wx = np.polynomial.hermite.hermgauss(nG)
    X0, X1, X2 = np.meshgrid(x, x, x, indexing='ij')
    Wg = (wx[:, None, None]*wx[None, :, None]*wx[None, None, :]).ravel()
    P0 = 1 + d*X0.ravel(); Px = d*X1.ravel(); Py = d*X2.ravel()   # weight e^{-|P-P0|^2/d^2} = |F^|^2
    Wg = Wg*d**3
    s2 = P0**2 - Px**2 - Py**2
    ok = (s2 > 0) & (P0 > 0)
    P0, Px, Py, Wg, s2 = P0[ok], Px[ok], Py[ok], Wg[ok], s2[ok]
    rs = np.sqrt(s2)
    th = 2*np.pi*np.arange(nth)/nth
    # boost from rest to P: velocity b = Pvec/P0, gamma = P0/rs
    g = P0/rs; bx, by = Px/P0, Py/P0; b2 = bx**2 + by**2
    def boost(e, kx, ky):   # boost (e,kx,ky) by velocity (bx,by)
        bk = bx[:, None]*kx + by[:, None]*ky
        gg = g[:, None]
        coef = np.where(b2[:, None] > 0, (gg - 1)/np.maximum(b2[:, None], 1e-300), 0.0)
        e2 = gg*(e + bk)
        kx2 = kx + coef*bk*bx[:, None] + gg*bx[:, None]*e
        ky2 = ky + coef*bk*by[:, None] + gg*by[:, None]*e
        return e2, kx2, ky2
    half = (rs/2)[:, None]*np.ones_like(th)[None, :]
    c, s = np.cos(th)[None, :], np.sin(th)[None, :]
    e1, p1x, p1y = boost(half, half*c, half*s)
    e2, p2x, p2y = boost(half, -half*c, -half*s)
    # wavefunction pieces (F^ factor is in the quadrature weight; its gradient appears explicitly)
    z1 = (p1x + 1j*p1y)/2; z2 = (p2x + 1j*p2y)/2
    Pi = -0.75*z1*z2 + 0.125*(z1**2 + z2**2)
    t1 = np.arctan2(p1y, p1x); t2 = np.arctan2(p2y, p2x)
    Em = e1*np.exp(1j*m*t1) + e2*np.exp(1j*m*t2)
    w = Em*Pi
    dP0 = (P0 - 1)[:, None]; dPx = Px[:, None]; dPy = Py[:, None]
    def grad(ei, pix, piy, ti, zo, zi):
        # gradient wrt p_i of  w * Fhat  divided by Fhat
        ux, uy = pix/ei, piy/ei                       # p-hat
        tx, ty = -np.sin(ti), np.cos(ti)              # theta-hat
        gE = np.exp(1j*m*ti)
        gEx, gEy = gE*(ux + 1j*m*tx), gE*(uy + 1j*m*ty)
        gFx = -(dP0*ux + dPx)/d**2; gFy = -(dP0*uy + dPy)/d**2
        gPi = (-0.75*zo + 0.25*zi)
        gPx, gPy = gPi*0.5, gPi*0.5j
        gx = gEx*Pi + Em*gFx*Pi + Em*gPx
        gy = gEy*Pi + Em*gFy*Pi + Em*gPy
        return gx, gy
    g1x, g1y = grad(e1, p1x, p1y, t1, z2, z1)
    g2x, g2y = grad(e2, p2x, p2y, t2, z1, z2)
    meas = (Wg/rs)[:, None]*(2*np.pi/nth)
    norm_psi = np.sum(meas*np.abs(Pi)**2)
    nw = np.sum(meas*np.abs(w)**2)
    P0w = np.sum(meas*(e1 + e2)*np.abs(w)**2)
    K0w = np.sum(meas*(e1*(np.abs(g1x)**2 + np.abs(g1y)**2) + e2*(np.abs(g2x)**2 + np.abs(g2y)**2)))
    # K^0 of psi itself (same formula with E_m -> 1)
    # family decomposition in the CM angle
    wf = np.fft.fft(w, axis=1)/nth                    # coefficient of e^{i j theta*}
    js = np.fft.fftfreq(nth, 1.0/nth).astype(int)
    pj = np.sum((Wg/rs)[:, None]*2*np.pi*np.abs(wf)**2, axis=0)
    frac_target = pj[js == m + 2].sum()/pj.sum()
    D2 = np.sum((np.abs(js) + 1)**2*pj)/pj.sum()
    D1 = np.sum((np.abs(js) + 1)*pj)/pj.sum()
    return dict(d=d, npsi=norm_psi, nw=nw, P0w=P0w, K0w=K0w, frac=frac_target, D2=D2, D1=D1, nth=nth)

if __name__ == '__main__':
    t0 = time.time()
    P("m   delta_m    ||E_m psi||^2/||psi||^2  ||K0^1/2 E_m psi||^2/||psi||^2   /(1+m)^2    <P0><K0>/||w||^4   (sum D||w_O||^2)^2   "
      "sum D^2||w_O||^2   (m+3)^2   target frac")
    rows = []
    for m in (0, 2, 4, 6, 8, 12, 16, 24, 32):
        r = run(m)
        a = r['K0w']/r['npsi']
        coer = r['P0w']*r['K0w']/r['nw']**2
        rows.append((m, a))
        P(f"{m:3d}  {r['d']:.5f}   {r['nw']/r['npsi']:14.8f}   {a:18.6f}   {a/(1+m)**2:12.6f}   {coer:14.4f}   {r['D1']**2:12.4f}   "
          f"{r['D2']:12.4f}   {(m+3)**2:6d}   {r['frac']:.6f}")
    P("eps-dependence at m = 8 (packet width delta = eps/9):")
    for eps in (0.05, 0.1, 0.15, 0.2):
        r = run(8, eps=eps)
        P(f"   eps={eps:4.2f}: ||K0^1/2 E_8 psi||^2/||psi||^2 = {r['K0w']/r['npsi']:.4f},  x eps^2/(1+m)^2 = "
          f"{r['K0w']/r['npsi']*eps**2/81:.5f},  target frac = {r['frac']:.6f}")
    P("odd m: E_m psi ~ 0 at rest (back-to-back pair, e^{im theta} + e^{im(theta+pi)}); e.g. m=3: ||E_3 psi||^2/||psi||^2 = "
      f"{(lambda r: r['nw']/r['npsi'])(run(3)):.3e}")
    P(f"[time {time.time() - t0:.1f}s]")
    open(os.path.join(os.path.dirname(__file__), 'cache', 'np4_k0.txt'), 'w').write('\n'.join(out) + '\n')
