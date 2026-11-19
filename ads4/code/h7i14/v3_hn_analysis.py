"""v3 (H7_I14_LIGHTSIDE.md §18-§19): (HN) for C^_2 on <TTTT> from the computed leading-Regge kernel (v3_regge_kernel).

Per nu (real, and nu +- i, nu +- 2i by analytic continuation of the closed forms):
  * f-stripped image  c0^T(nu) P_l      (fit; residual = test of l-isotypy)        and the scalar c0^sc(nu)
  * eta-kernel image  h-(nu) P_{l-1} + h0(nu) P_l + h+(nu) P_{l+1}    (fit; residual = test of the 3-term structure)
  * both helicity classes (h1, h2) = (+,-) [++] and (+,+) [+-]
Derived quantities:
  a_-+(nu) := 2 h_-+^T / (c0^sc g, c0^sc g~)       (= 1 for the scalar functional; common nu-only factors cancel)
  u_T := (h- - h+)/(h- + h+),  sigma_T := Im u_T     (forward-normalised local shift coefficient; scalar: u_sc = (g - g~)/(g + g~))
  exact (HN) test for the scalar-g de-shift with the scalar normalisation:  a_-(mu - i) == a_+(mu + i)  (then N^_2 = a_-(mu - i))
  light-side (HN) defect of the forward-normalised scalar-g de-shift (first order): eps(nu) = Dsigma 2 nu (nu^2 + 9/4)/((nu^2+5/4)^2 + 4 nu^2)
Output: cache/v3_hn_analysis.txt, cache/v3_uT_table.pkl (u_T, a_-+, h's on a nu grid; used by v3_heavy_impact.py).
Run time: ~10-15 min on one core.
"""
import os, sys, io, time, pickle
os.environ.setdefault('OMP_NUM_THREADS', '2')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import v3_regge_kernel as V

buf = io.StringIO()
def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); buf.write(s + "\n")


def betas_for(nu):
    nr = abs(np.real(nu))
    return np.linspace(0.12, min(3.0, max(0.6, 30/max(nr, 1e-9))), 16)


def coefs(nu, hel=(1, -1)):
    bs = betas_for(nu)
    (hm, h0, hp), res3 = V.fit3(nu, bs, [V.heavy_image(nu, b, *hel) for b in bs])
    cT, resT = V.fit1(nu, bs, [V.heavy_image(nu, b, *hel, with_eta=False) for b in bs])
    cs, resS = V.fit1(nu, bs, [V.heavy_image(nu, b, 1, 1, spin=False, with_eta=False) for b in bs])
    (sm, s0, sp), resS3 = V.fit3(nu, bs, [V.heavy_image(nu, b, 1, 1, spin=False) for b in bs])
    gp, gm = V.gfac(nu)
    return dict(nu=nu, hm=hm, h0=h0, hp=hp, cT=cT, cs=cs, sm=sm, sp=sp, gp=gp, gm=gm,
                am=2*hm/(cs*gp), ap=2*hp/(cs*gm), res=max(res3, resT, resS, resS3))


if __name__ == '__main__':
    t0 = time.time()
    # ---------------------------------------------------------------- 1. helicity classes and h0
    pr("1. Both helicity classes, h0, 3-term residual (spinning Psi_{2,nu}[q^2], a = 7):")
    for nu in (0.05, 0.5, 2.0, 5.0, 12.0):
        A = coefs(nu, (1, -1)); B = coefs(nu, (1, 1))
        pr(f"   nu={nu:5}: max fit resid {max(A['res'], B['res']):.1e}; |h0|/(|h-|+|h+|) = {abs(A['h0'])/(abs(A['hm'])+abs(A['hp'])):.1e}"
           f" (++), {abs(B['h0'])/(abs(B['hm'])+abs(B['hp'])):.1e} (+-); |h-(++)/h-(+-) - 1| = {abs(A['hm']/B['hm'] - 1):.1e},"
           f" |h+(++)/h+(+-) - 1| = {abs(A['hp']/B['hp'] - 1):.1e}")
    # ---------------------------------------------------------------- 2. table on the real axis
    grid = np.concatenate([[0.05, 0.25], np.arange(0.5, 10.01, 0.5), np.arange(11, 20.01, 1.0), [22.5, 25, 30, 35, 40, 50, 60, 80, 100,
                                                                                               140, 200, 300]])
    tab = []
    pr("\n2. Table (++ = +-).  a_-+ = 2h_-+/(c0^sc g, g~) [scalar: 1, 1];  r = c0^T/c0^sc;  u = (h- - h+)/(h- + h+);  sigma = Im u;"
       "  eps = first-order light-side defect of the forward-normalised scalar-g de-shift")
    pr("    nu     resid    a_-                     a_+                     |h-/h+|   Re u_T    sigma_T     sigma_sc    nu^3 Dsigma   eps")
    for nu in grid:
        c = coefs(nu)
        uT = (c['hm'] - c['hp'])/(c['hm'] + c['hp']); us = (c['gp'] - c['gm'])/(c['gp'] + c['gm'])
        ds = uT.imag - us.imag
        eps = ds*2*nu*(nu*nu + 2.25)/((nu*nu + 1.25)**2 + 4*nu*nu)
        c.update(uT=uT, us=us, ds=ds, eps=eps); tab.append(c)
        pr(f"  {nu:6.2f}  {c['res']:.0e}  {c['am'].real:+9.4f}{c['am'].imag:+9.4f}i  {c['ap'].real:+9.4f}{c['ap'].imag:+9.4f}i  "
           f"{abs(c['hm']/c['hp']):7.4f}  {uT.real:+8.4f}  {uT.imag:+10.5f}  {us.imag:+10.5f}  {nu**3*ds:+11.2f}  {eps:+9.2e}")
    pickle.dump(tab, open(os.path.join(HERE, 'cache', 'v3_uT_table.pkl'), 'wb'))
    # ---------------------------------------------------------------- 3. large-nu asymptotics
    pr("\n3. Large-nu fits (nu = 40 ... 300):")
    sel = [c for c in tab if c['nu'] >= 40]
    nus = np.array([c['nu'] for c in sel])
    for lab, y in (("nu*sigma_T", np.array([c['nu']*c['uT'].imag for c in sel])),
                   ("nu*sigma_sc", np.array([c['nu']*c['us'].imag for c in sel])),
                   ("nu^3*Dsigma", np.array([c['nu']**3*c['ds'] for c in sel])),
                   ("nu^2*Re u_T", np.array([c['nu']**2*c['uT'].real for c in sel]))):
        M = np.array([np.ones_like(nus), 1/nus, 1/nus**2]).T
        sol = np.linalg.lstsq(M, y, rcond=None)[0]
        pr(f"   {lab:13s} = {sol[0]:+.4f} {sol[1]:+.3f}/nu {sol[2]:+.1f}/nu^2")
    # ---------------------------------------------------------------- 4. exact (HN) test for the scalar-g de-shift, scalar normalisation
    pr("\n4. Exact (HN) test for C^_2^{cur} = A_g^{-1}[Psi^T/N_sc] (same normalisation as the scalar Psi): need a_-(mu - i) = a_+(mu + i)")
    for mu in (0.7, 2.0, 4.0, 9.0, 20.0):
        cm = coefs(mu - 1j); cp = coefs(mu + 1j)
        pr(f"   mu={mu:5}: a_-(mu-i) = {cm['am']:.6g}, a_+(mu+i) = {cp['ap']:.6g}, ratio - 1 = {abs(cm['am']/cp['ap'] - 1):.2e}"
           f"   [resid {max(cm['res'], cp['res']):.0e}]")
    # ---------------------------------------------------------------- 5. forward graviton/scalar ratio vs the crystal nu-only factor (H0 §22.1)
    pr("\n5. Forward (beta = 0) ratio R(nu) = (h-^T + h+^T)/(h-^sc + h+^sc): fit log R = p log nu + c + a1/nu + a2/nu^2 (+ a3/nu^3),"
       " nu = 40...1000.  H0 §22.1 found the crystal graviton/scalar nu-only factor exp(22i/nu + 82/nu^2 + ...)")
    nus = np.array([40, 60, 80, 100, 140, 200, 300, 450, 700, 1000.])
    vals = []
    for nu in nus:
        bs = np.linspace(0.02, min(3.0, 30/nu), 12)
        (hm, h0, hp), r1 = V.fit3(nu, bs, [V.heavy_image(nu, b, 1, -1) for b in bs])
        (sm, s0, sp_), r2 = V.fit3(nu, bs, [V.heavy_image(nu, b, 1, 1, spin=False) for b in bs])
        vals.append((hm + hp)/(sm + sp_))
    Lg = np.log(np.array(vals)); Lg = Lg.real + 1j*np.unwrap(Lg.imag)
    Mx = np.array([np.log(nus), np.ones_like(nus), 1/nus, 1/nus**2, 1/nus**3]).T
    for cols in (4, 5):
        sol = np.linalg.lstsq(Mx[:, :cols], Lg, rcond=None)[0]
        pr(f"   {cols} terms: p = {sol[0]:.4f}, a1 = {sol[2]:.3f}, a2 = {sol[3]:.2f}" + (f", a3 = {sol[4]:.1f}" if cols == 5 else ""))
    pr(f"\n[{time.time() - t0:.0f}s]")
    open(os.path.join(HERE, 'cache', 'v3_hn_analysis.txt'), 'w').write(buf.getvalue())
