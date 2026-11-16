"""
d = 3 leading-Regge form of the Chang-Landau-Simmons-Duffin functional (arXiv:2311.04271, eq. (3.61), `spinningfunctional_definition')
for O_1 = O_3 conserved of spin l and O_2 = O_4 scalars, expressed through the tensor Regge moments Pihat^{(q1,q3)}.

Light transforms fix the polarisations of O_1, O_3 to the null directions z_1 = x_{1+ 0}, z_3 = x_{30}; with x_1 = -y, x_3 = y these are
yhat + n and yhat - n (n a unit vector orthogonal to y): antipodal points th_1 = phi, th_3 = phi + pi of the celestial circle.
The celestial structure <P~_{d1}(z_1) P~_{d3}(z_3) P_{delta,j'}(z_0)> (delta = 1/2 + i nu) gives the weight |sin(th_0 - phi)|^{-delta}
(times sgn(sin) for the parity-odd structure j' = 1), and the dual structure gives (-z_0.xhat)^{-1/2+i nu} at leading order in the
Regge limit.  Hence, up to a (k, nu)-dependent normalisation that is not tracked,
    Psi_{k,nu,j'=0}  ~  sum_{q1+q3 even} (-1)^{q3} c_{(q1+q3)/2}(delta)/c_0(delta)  Pihat^{(q1,q3)}_{k_r,nu},
    c_n/c_0 = (delta/2)_{|n|}/(1-delta/2)_{|n|}          [Fourier coefficients of |sin psi|^{-delta}],
    Psi_{k,nu,j'=1}  ~  sum_{q1+q3 odd} (-1)^{q3} d_{(q1+q3-1)/2}(delta)/d_0(delta)  Pihat^{(q1,q3)}_{k_r,nu},
    d_n/d_0: Fourier coefficients of sgn(sin psi)|sin psi|^{-delta} e^{-i psi}.
Heavy action:  A_pres(nu) P_{-1/2+i nu}(cosh beta) + A_flip(nu) F^nu_{2l}(beta),  A = sum w M (laplace.M_table).
"""
import sympy as sp, mpmath as mp
from laplace import M_table, nu, a, g1

def weights(l, parity=0):
    delta = sp.Rational(1, 2) + sp.I*nu
    w = {}
    for q1 in range(-l, l+1):
        for q3 in range(-l, l+1):
            if (q1 + q3) % 2 != parity: continue
            if parity == 0:
                n = abs(q1 + q3)//2
                c = sp.rf(delta/2, n)/sp.rf(1 - delta/2, n)
            else:
                # sgn(sin psi)|sin psi|^{-delta} = sum_m d_m e^{i(2m+1) psi}; d_m/d_0 = ((delta+1)/2)_m/((3-delta)/2)_m (m>=0), d_{-m-1} = -d_m
                mm = (q1 + q3 - 1)//2
                sgn = 1
                if mm < 0:
                    mm = -mm - 1; sgn = -1
                c = sgn*sp.rf((delta + 1)/2, mm)/sp.rf((3 - delta)/2, mm)
            w[(q1, q3)] = (-1)**q3*c
    return w

def heavy_response(l, parity=0):
    T = M_table(l); w = weights(l, parity)
    pres = flip = 0
    for (q1, q3), ww in w.items():
        for h1 in (1, -1):
            for h3 in (1, -1):
                R, S, mu = T[(q1, q3, h1, h3)]
                val = ww*(R + g1*S)
                if mu == 0: pres += val
                else: flip += val
    return sp.factor(sp.cancel(sp.together(sp.expand(pres)))), sp.factor(sp.cancel(sp.together(sp.expand(flip))))

if __name__ == "__main__":
    out = []
    def pr(*x):
        s = " ".join(str(v) for v in x); print(s, flush=True); out.append(s)
    # numerical check of the Fourier coefficients
    mp.mp.dps = 20
    nv = mp.mpf('1.3'); dl = mp.mpf(1)/2 + 1j*nv
    c = lambda n: mp.quad(lambda p: mp.sin(p)**(-dl)*mp.exp(-2j*n*p), [0, mp.pi/2, mp.pi])/mp.pi
    pr("# check c_n/c_0 = (delta/2)_n/(1-delta/2)_n at nu = 1.3:", [mp.nstr(c(n)/c(0), 8) for n in (1, 2)], [mp.nstr(mp.rf(dl/2, n)/mp.rf(1-dl/2, n), 8) for n in (1, 2)])
    d = lambda m: mp.quad(lambda p: mp.sin(p)**(-dl)*mp.exp(-1j*(2*m+1)*p), [0, mp.pi/2, mp.pi])/mp.pi
    pr("# check d_m/d_0 at nu = 1.3:", [mp.nstr(d(m)/d(0), 8) for m in (1, 2)], [mp.nstr(mp.rf((dl+1)/2, m)/mp.rf((3-dl)/2, m), 8) for m in (1, 2)])
    for l in (1, 2):
        for parity in (0, 1):
            pres, flip = heavy_response(l, parity)
            pr(f"l = {l}, j' = {parity}:")
            pr("   A_pres(nu; a) =", pres)
            pr("   A_flip(nu; a) =", flip)
            for kk in (2, 4):
                aa = 2*(l+1) + kk - 1
                rr = sp.cancel(sp.together((flip/pres).subs(a, aa))) if pres != 0 else None
                pr(f"   a = {aa}: A_flip/A_pres =", sp.factor(rr) if rr is not None else None)
                if rr is not None and not rr.has(g1):
                    eps = sp.Symbol('eps', positive=True)
                    pr("        large-nu behaviour of A_flip/A_pres (eps = 1/nu):", sp.series(sp.cancel(rr.subs(nu, 1/eps)), eps, 0, 6))
    open('out_clsd_kernel.txt', 'w').write("\n".join(out) + "\n")
