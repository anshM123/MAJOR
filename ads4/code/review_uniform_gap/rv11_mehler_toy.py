"""rv11_mehler_toy.py -- numerical sanity test (not a proof) of the relative Mehler bound quoted in the note Sec.5:
 |int_0^inf r(nu) P_{-1/2+i nu}(cosh b) dnu| <= (sqrt2/pi) K0 [B/b + 8 B log2/b^2 + A eps1(b)],
 K0 = (cosh b - 1)^(-1/2), A = ||r||_1, B = |r'(0+)| + ||r''||_1, eps1 = (sinh(b/2)-b/2)/(2(cosh b - cosh(b/2))),
on toy residuals with r(0) = 0. LHS by the Mehler-Dirichlet form (sqrt2/pi) int_0^b rhat(t)/sqrt(cosh b - cosh t) dt."""
import mpmath as mp
mp.mp.dps = 40
toys = {
 'nu e^-nu': (lambda n: n*mp.exp(-n), lambda t: (1-t*t)/(1+t*t)**2),
 # r = nu e^{-nu} cos(3 nu): rhat = (1/2) Re[1/(1-i(t+3))^2 + 1/(1-i(t-3))^2]
 'nu e^-nu cos3nu': (lambda n: n*mp.exp(-n)*mp.cos(3*n), lambda t: mp.re(1/(1-1j*(t+3))**2 + 1/(1-1j*(t-3))**2)/2),
 # r = sin(nu) e^{-nu/2}: rhat = Re int sin(nu) e^{-nu/2} cos(nu t) = (1/2)[ (1+t)/((1/4)+(1+t)^2) + (1-t)/((1/4)+(1-t)^2) ]
 'sin(nu) e^-nu/2': (lambda n: mp.sin(n)*mp.exp(-n/2), lambda t: ((1+t)/(mp.mpf(1)/4+(1+t)**2) + (1-t)/(mp.mpf(1)/4+(1-t)**2))/2),
}
def lhs(rhat, b):
    b = mp.mpf(b)
    # t = b(1-v^2)
    f = lambda v: rhat(b*(1-v*v))*2*b*v/mp.sqrt(2*mp.sinh(b*(1-v*v/2))*mp.sinh(b*v*v/2)) if v > 0 else rhat(b)*2*b/mp.sqrt(b*mp.sinh(b))
    return mp.sqrt(2)/mp.pi*mp.quad(f, mp.linspace(0, 1, 21))
for name, (r, rhat) in toys.items():
    A = mp.quad(lambda n: abs(r(n)), mp.linspace(0, 80, 321))
    d1 = mp.diff(r, mp.mpf('1e-30'))
    B = abs(d1) + mp.quad(lambda n: abs(mp.diff(r, n, 2)), mp.linspace(0, 80, 321))
    # check rhat formula and int rhat = 0
    chk = abs(mp.quad(lambda n: r(n)*mp.cos(n*mp.mpf('0.7')), mp.linspace(0, 120, 481)) - rhat(mp.mpf('0.7')))
    print(f"toy {name}: A={mp.nstr(A,6)} B={mp.nstr(B,6)} rhat check {mp.nstr(chk,3)}")
    for b in [12, 15, 20, 30, 45, 60]:
        bm = mp.mpf(b)
        K0 = 1/mp.sqrt(mp.cosh(bm)-1)
        e1 = (mp.sinh(bm/2)-bm/2)/(2*(mp.cosh(bm)-mp.cosh(bm/2)))
        R = mp.sqrt(2)/mp.pi*K0*(B/bm + 8*B*mp.log(2)/bm**2 + A*e1)
        L_ = lhs(rhat, b)
        print(f"   beta={b}: |LHS| = {mp.nstr(abs(L_),6)}  bound = {mp.nstr(R,6)}  ratio = {mp.nstr(abs(L_)/R,4)}")
