# r2_special.py -- numerical checks of the special-function identities used by ALL_GAP_NOTE.
#  (a) Mehler-Fock expansion of 1/(x+A)  (note sec 2.1, "DLMF 14.20.14")
#  (b) int_1^oo f(x) P_{-1/2+i nu}(x) dx with f=1/(x+A) equals pi P(A)/cosh(pi nu) (the transform behind (a))
#  (c) resolvent identity int_1^oo Q_1(x) P_{-1/2+i nu}(x) dx = 1/(nu^2+9/4)
#  (d) Phi>0 on R, Phi asymptotics C|x|^{23/2}e^{-pi|x|}; B real-even; B zeros/poles in the strip; B = Phi(z+i)/R_- = Phi(z-i)/R_+
#  (e) Mehler-Dirichlet bounds |P_{-1/2+i nu}(A)|<=P_{-1/2}(A), |d/dnu P|<= arccosh(A) P_{-1/2}(A)
import mpmath as mp
mp.mp.dps = 30
Pl = lambda nu, x: mp.legenp(-0.5+1j*nu, 0, x, type=3)
def Q1(x): return x/2*mp.log((x+1)/(x-1))-1
print('(a) Mehler-Fock 1/(x+A) = pi int_0^oo nu tanh(pi nu)/cosh(pi nu) P(x) P(A) dnu')
for (x, A) in [(1, 1.5), (2, 1.5), (5, 1.01), (1.3, 3.0), (20, 1.0005)]:
    x = mp.mpf(x); A = mp.mpf(A)
    val = mp.pi*mp.quad(lambda n: n*mp.tanh(mp.pi*n)/mp.cosh(mp.pi*n)*mp.re(Pl(n, x)*Pl(n, A)), [0, 2, 5, 10, 20, 40])
    print('   x=%-6s A=%-7s  rhs=%s  1/(x+A)=%s  rel.err=%.2e' % (x, A, mp.nstr(val, 15), mp.nstr(1/(x+A), 15), abs(val*(x+A)-1)))
print('(b) int_1^oo P(x)/(x+A) dx = pi P(A)/cosh(pi nu)')
for (nu, A) in [(0.3, 1.5), (2.0, 1.01), (5.0, 3.0)]:
    nu = mp.mpf(nu); A = mp.mpf(A)
    lhs = mp.quad(lambda x: mp.re(Pl(nu, x))/(x+A), [1, 2, 10, 100, 1e4, mp.inf])
    rhs = mp.pi*mp.re(Pl(nu, A))/mp.cosh(mp.pi*nu)
    print('   nu=%s A=%s lhs=%s rhs=%s rel=%.2e' % (nu, A, mp.nstr(lhs, 14), mp.nstr(rhs, 14), abs(lhs/rhs-1)))
print('(c) int_1^oo Q_1 P_{-1/2+i nu} dx = 1/(nu^2+9/4)')
for nu in [0, 0.5, 1.7, 4.0, 9.0]:
    nu = mp.mpf(nu)
    # substitute x=(1+z)/(1-z): dx = 2dz/(1-z)^2, Q_1 = (1+z)(-log z)/(2(1-z)) - 1
    f = lambda zz: mp.re(Pl(nu, (1+zz)/(1-zz)))*((1+zz)*(-mp.log(zz))/(2*(1-zz))-1)*2/(1-zz)**2
    lhs = mp.quad(f, mp.linspace(0, 1, 30))
    print('   nu=%-4s lhs=%s 1/(nu^2+9/4)=%s rel=%.2e' % (nu, mp.nstr(lhs, 15), mp.nstr(1/(nu**2+2.25), 15), abs(lhs*(nu**2+2.25)-1)))
print('    Q_1>0 on (1,oo):', all(Q1(mp.mpf(x)) > 0 for x in [1.0001, 1.1, 2, 10, 100, 1e4]))
print('    (x^2-1)Q_1\'(x) -> -1 at x->1:', mp.nstr((mp.mpf('1.000001')**2-1)*mp.diff(Q1, mp.mpf('1.000001')), 8))

# (d) the new inverse
Pp = lambda w: (w**2+0.25)**2*(w**2+2.25)*(w**2+6.25)*(w**2+20.25)
Fstar = lambda w: Pp(w)*mp.cosh(mp.pi*w/2)*mp.gamma(1.25+1j*w/2)*mp.gamma(1.25-1j*w/2)
Phi = lambda w: Fstar(w)/mp.cosh(mp.pi*w)
Rm = lambda w: -2*(2*w-1j)*(2*w+1j)**2*(2*w+3j)*(2*w+5j)*(2*w+11j)/(w*(2*w-7j)*(2*w-3j)*(2*w+7j))
Rp = lambda w: Rm(-w)
B1 = lambda w: Phi(w+1j)/Rm(w)
B2 = lambda w: Phi(w-1j)/Rp(w)
print('(d) new inverse')
mn = min(mp.re(Phi(mp.mpf(x))) for x in mp.linspace(-30, 30, 601))
print('   min Phi on [-30,30] grid:', mp.nstr(mn, 6), ' max |Im Phi|:', mp.nstr(max(abs(mp.im(Phi(mp.mpf(x)))) for x in mp.linspace(-30, 30, 61)), 3))
for x in [20, 40, 80, 160]:
    x = mp.mpf(x)
    print('   Phi(x)/(x^{23/2}e^{-pi x}) at x=%s : %s' % (x, mp.nstr(mp.re(Phi(x))/(x**11.5*mp.exp(-mp.pi*x)), 10)))
worst = 0
for w in [mp.mpc(0.3, 0.2), mp.mpc(-1.7, -0.6), mp.mpc(2.5, 0.95), mp.mpc(7, -0.3), mp.mpc(0.01, 0.49)]:
    worst = max(worst, abs(B1(w)/B2(w)-1))
print('   max |B1/B2-1| at complex points in the strip:', mp.nstr(worst, 3))
print('   B real-even on R:', mp.nstr(max(abs(mp.im(B1(mp.mpf(x))))/abs(B1(mp.mpf(x))) for x in [0.1, 1, 3, 9]), 3),
      mp.nstr(max(abs(B1(mp.mpf(x))-B1(mp.mpf(-x)))/abs(B1(mp.mpf(x))) for x in [0.1, 1, 3, 9]), 3))
print('   B(x)/x^2 as x->0 :', [mp.nstr(mp.re(B1(mp.mpf(t)))/mp.mpf(t)**2, 10) for t in ['1e-3', '1e-5']])
for p in [0.5j, -0.5j]:
    h = mp.mpf('1e-6')
    print('   (z-%s)B(z) near %s :' % (p, p), mp.nstr(h*B1(p+h), 10), mp.nstr(h*B1(p+h*1j)/1j, 10), '(constant => simple pole)')
# no other poles: scan |B| on a grid of the closed strip away from +-i/2 (look for blow-ups)
mx = 0
for xr in mp.linspace(-6, 6, 61):
    for yi in mp.linspace(-1, 1, 41):
        w = mp.mpc(xr, yi)
        if abs(w-0.5j) < 0.05 or abs(w+0.5j) < 0.05 or abs(w) < 1e-9: continue
        mx = max(mx, abs(B1(w)))
print('   max|B| on strip grid away from +-i/2 (|x|<=6):', mp.nstr(mx, 5))
for x in [10, 20, 40]:
    for y in [-1, -0.5, 0, 0.7, 1]:
        w = mp.mpc(x, y)
        r = abs(B1(w))/(x**9.5*mp.exp(-mp.pi*x))
        print('   |B(x+iy)|/(x^{19/2}e^{-pi x}) x=%s y=%s : %s' % (x, y, mp.nstr(r, 6)))
# (e) Mehler-Dirichlet bounds
A = mp.mpf('1.00001')
P0 = mp.re(Pl(0, A))
mxP = max(abs(Pl(n, A)) for n in mp.linspace(0, 200, 201))
mxD = max(abs(mp.diff(lambda t: mp.re(Pl(t, A)), n)) for n in mp.linspace(0, 200, 41))
print('(e) A=1.00001: max|P|/P_{-1/2}(A)=%s, max|dP/dnu|/(arccosh(A)P_{-1/2}(A))=%s' % (mp.nstr(mxP/P0, 6), mp.nstr(mxD/(mp.acosh(A)*P0), 6)))
