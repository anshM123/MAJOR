"""rv1_symbolic.py -- independent symbolic checks of RESEARCH_NOTE Sec.1 (Gamma/rational identities) and Sec.2
(light polynomial identity, Einstein test function).  Written from the note's statements; no code of the package used."""
import sympy as sp
import mpmath as mp
z = sp.symbols('z')
I = sp.I
Rm = -2*(2*z-I)*(2*z+I)**2*(2*z+3*I)*(2*z+5*I)*(2*z+11*I)/(z*(2*z-7*I)*(2*z-3*I)*(2*z+7*I))
Rp = Rm.subs(z, -z)
P = (z**2+sp.Rational(1,4))**2*(z**2+sp.Rational(9,4))*(z**2+sp.Rational(25,4))*(z**2+sp.Rational(81,4))
out = []
def rep(s):
    print(s); out.append(s)

# (a) F(z+i)/F(z-i): cosh ratio = -1 ; Gamma ratio = (3/4 - i z/2)/(3/4 + i z/2)  (Gamma recurrence, derived by hand)
ratio_F = -P.subs(z, z+I)/P.subs(z, z-I)*(sp.Rational(3,4)-I*z/2)/(sp.Rational(3,4)+I*z/2)
d = sp.simplify(sp.together(ratio_F - Rm/Rp))
rep(f"(1a) F(z+i)/F(z-i) - R_-/R_+  (symbolic, after Gamma/cosh recurrences) = {d}")
# numeric check of the same with actual Gamma, cosh at random complex points (no hand recurrence)
mp.mp.dps = 40
def Fn(w):
    Pw = (w**2+mp.mpf(1)/4)**2*(w**2+mp.mpf(9)/4)*(w**2+mp.mpf(25)/4)*(w**2+mp.mpf(81)/4)
    return Pw*mp.cosh(mp.pi*w/2)*mp.gamma(mp.mpf(5)/4+1j*w/2)*mp.gamma(mp.mpf(5)/4-1j*w/2)
Rmn = sp.lambdify(z, Rm, 'mpmath'); Rpn = sp.lambdify(z, Rp, 'mpmath')
mx = 0
for w in [mp.mpc(0.3,0.2), mp.mpc(2.7,-0.4), mp.mpc(11.1,0.9), mp.mpc(-5.3,0.1), mp.mpc(40.2,-0.33)]:
    lhs = Fn(w+1j)/Fn(w-1j); rhs = Rmn(w)/Rpn(w)
    mx = max(mx, abs(lhs/rhs-1))
rep(f"(1a) numeric max |F(z+i)/F(z-i) / (R_-/R_+) - 1| at 5 complex points (mpmath Gamma, 40 dps) = {mp.nstr(mx,3)}")
# polynomial identity
target = -z*(2*z-7*I)**2*(2*z-3*I)**2*(2*z+3*I)*(2*z+7*I)**2/2048
d2 = sp.simplify(sp.together(P.subs(z, z+I)/Rm - target))
rep(f"(1a) P(z+i)/R_-(z) - [-z(2z-7i)^2(2z-3i)^2(2z+3i)(2z+7i)^2/2048] = {d2}")
# B even/real: B(x)=F(x+i)/R_-(x) vs F(x-i)/R_+(x), and B(-x)=B(x), Im B(x)=0
mx1 = mx2 = mx3 = 0
for x in [mp.mpf('0.37'), mp.mpf('1.9'), mp.mpf('7.3'), mp.mpf('55.5'), mp.mpf('300.25')]:
    B1 = Fn(x+1j)/Rmn(x); B2 = Fn(x-1j)/Rpn(x); B3 = Fn(-x+1j)/Rmn(-x)
    mx1 = max(mx1, abs(B1/B2-1)); mx2 = max(mx2, abs(B3/B1-1)); mx3 = max(mx3, abs(mp.im(B1)/mp.re(B1)))
rep(f"(1a) B two forms agree: {mp.nstr(mx1,3)}; even: {mp.nstr(mx2,3)}; |Im B/Re B| on R: {mp.nstr(mx3,3)}")
# R_+(z) = conj(R_-(conj z)) (needed for 'LHS = 2 Re[R_-(x-i)chi(x-i)]')
d3 = sp.simplify(Rp - sp.conjugate(Rm.subs(z, sp.conjugate(z))))
w = sp.symbols('w', real=True)
d3 = sp.simplify((Rp - sp.conjugate(Rm.subs(z, sp.conjugate(z)))).subs(z, w))
rep(f"(1c) R_+(x) - conj(R_-(x)) on real x: {d3}")

# Sec.2: light polynomial.  K_B2A4 Sec.4.5: P_K(nu) = -(4nu^2+1) N10 /(8(4nu^2+9)(4nu^2+49))
N10 = 1024*z**10 - 195328*z**8 + 8229504*z**6 - 90780512*z**4 + 251603732*z**2 - 114141375
PK = -(4*z**2+1)*N10/(8*(4*z**2+9)*(4*z**2+49))
Q = z**6 - sp.Rational(245,4)*z**4 + sp.Rational(9619,16)*z**2 - sp.Rational(51575,64)
d4 = sp.simplify(sp.together(Rm*Q.subs(z, z+I) + Rp*Q.subs(z, z-I) - PK))
rep(f"(2) R_-Q(z+i)+R_+Q(z-i) - P_KB2A4(z) = {d4}")
# also: is the forward factor R_-+R_+ consistent with K_B2A4's (16nu^4-728nu^2+969)?
rep(f"(2) R_- + R_+ = {sp.factor(sp.simplify(Rm+Rp))}")
# Einstein: K_B2A4 Sec.4.4: R_- D(nu+i) + R_+ D(nu-i) = -32(16nu^4-440nu^2+273)/((4nu^2+9)(4nu^2+49)); find D=c/(nu^2+9/4)
c = sp.symbols('c')
D = c/(z**2+sp.Rational(9,4))
expr = sp.together(Rm*D.subs(z, z+I) + Rp*D.subs(z, z-I) + 32*(16*z**4-440*z**2+273)/((4*z**2+9)*(4*z**2+49)))
num = sp.Poly(sp.numer(expr), z)
sol = sp.solve(num.coeffs(), c)
rep(f"(2) Einstein: R_-D(z+i)+R_+D(z-i) = Lemma-R Einstein function iff D = c/(z^2+9/4) with c = {sol}")
open('out/rv1_symbolic.log', 'w').write('\n'.join(out)+'\n')
