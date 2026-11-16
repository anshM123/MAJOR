r"""Check of Lemma M (ads4/POSITIVITY_LEMMA.md): for phi_nu(beta) = P_{(2-d)/2 + i nu}(cosh beta) (spherical function on
H^{d-1}), the ratio r(beta) = -[d phi/d(nu^2)]_{nu=0} / phi_0(beta) (eq. (5.35) of arXiv:2106.10274 up to a constant)
satisfies  r'(beta) = \int_0^beta w / w(beta),  w = phi_0^2 sinh^{d-2}.   d = 3: phi_nu = 2F1(1/2-i nu, 1/2+i nu; 1; (1-cosh b)/2);
d = 4: phi_nu = sin(nu b)/(nu sinh b), r = b^2/6;  d = 5: phi_nu = 2F1(3/2 - i nu, 3/2 + i nu; 2; (1-cosh b)/2)."""
import mpmath as mp
mp.mp.dps = 30
def phi(d, nu, b):
    rho = mp.mpf(d-2)/2
    return mp.re(mp.hyp2f1(rho-1j*nu, rho+1j*nu, mp.mpf(d-1)/2, (1-mp.cosh(b))/2))
def r(d, b, e=mp.mpf('1e-6')):
    return -(phi(d, e, b)-phi(d, 0, b))/(e*e)/phi(d, 0, b)
for d in (3, 4, 5):
    print(f"d = {d}")
    for b in map(mp.mpf, ('0.3', '1', '3', '8')):
        w = lambda s: phi(d, 0, s)**2*mp.sinh(s)**(d-2)
        formula = mp.quad(w, [0, b])/w(b); h = mp.mpf('1e-4'); num = (r(d, b+h)-r(d, b-h))/(2*h)
        extra = f"   r = {mp.nstr(r(d,b),8)}, b^2/6 = {mp.nstr(b*b/6,8)}" if d == 4 else ""
        print(f"   beta={mp.nstr(b,3):>4}:  r'(numerical) = {mp.nstr(num, 10)}   formula = {mp.nstr(formula, 10)}{extra}")
