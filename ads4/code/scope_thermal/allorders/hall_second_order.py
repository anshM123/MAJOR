"""Item 5 (Hall): O(lam*kap) term in eta_H/s, kap Q on the O(lam) J-corrected brane.  [computed]"""
import io, contextlib, time
import sympy as sp
t0 = time.time()
with contextlib.redirect_stdout(io.StringIO()):
    import ef_engine_lam as E
r, lam = E.r, E.lamb
print("engine built in %.0f s" % (time.time() - t0))
print("background I, J, P, Q (to O(lam)):", [sp.simplify(E._tr(t[0])) for t in (E.I_, E.J_, E.P_, E.Q_)])
with contextlib.redirect_stdout(io.StringIO()) as buf:
    sig = E.hall(E.cross_part(E._tr(sp.expand(E.r ** 2 * E._tr(E.Nb * 1) * 0 + E.sqrtg[0] * E.Q_[2]))), 'Q', integrate=False)
print(buf.getvalue() if False else "", end="")
sig = sp.expand(sig)
s0, s1 = sp.simplify(sig.coeff(lam, 0)), sp.simplify(sig.coeff(lam, 1))
i0 = sp.integrate(s0, (r, 1, sp.oo)); i1 = sp.integrate(s1, (r, 1, sp.oo))
print("sigma_Q = ", s0, " + lam *", s1)
print("int sigma = ", i0, " + lam *", i1)
# 4 pi eta_H/s = int sigma / s_hat,  s_hat = 1 + 6 lam (Wald, fixed r_h, ../etas_diag.py / second_order.py)
kap, a3, a3t = sp.symbols('kappa alpha3 alpha3t')
res = sp.expand(sp.series(kap * (i0 + lam * i1) / (1 + 6 * lam), lam, 0, 2).removeO())
print("4 pi eta_H/s = ", res)
print("in alphas (lam = -alpha3/6, kap = -alpha3t/6):", sp.expand(res.subs({lam: -a3 / 6, kap: -a3t / 6})))
