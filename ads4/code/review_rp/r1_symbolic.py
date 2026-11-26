# r1_symbolic.py -- independent symbolic checks of the algebraic claims of ALL_GAP_NOTE (Theorem RP).
# Written from the paper (sec_fullpos.tex fp-U, fp-single, fp-V; prl/data/kernels.json R_-,R_+; K_B2A4.md P(nu))
# and the note's own definitions.  Run: python3 r1_symbolic.py
import sympy as sp
z, nu, k, M, b, T = sp.symbols('z nu k M b T')
I = sp.I
out = []
def rep(name, ok, extra=''):
    out.append((name, ok)); print(('PASS ' if ok else 'FAIL ') + name + (' : ' + str(extra) if extra != '' else ''))

# ---------- 1. operator D_b and the product rule on P_{k,M}=z^k(1-z)^M -------------------------------------
Db = lambda f, bb: z*(1-z)/2*sp.diff(f, z, 2) + ((sp.Rational(7, 2)-2*bb)*z+sp.Rational(1, 2))*sp.diff(f, z) + (sp.Rational(9, 2)-bb*(bb-1))*f
P = lambda kk, mm: z**kk*(1-z)**mm
C = k*(sp.Rational(7, 2)-2*b-M)+sp.Rational(9, 2)-b*(b-1)-M/2
claim = k*(k-1)/2*P(k-1, M+1)+k/2*P(k-1, M)+C*P(k, M)+M*(M-9+4*b)/2*P(k+1, M-1)
diff = sp.simplify(sp.expand(sp.powsimp((Db(P(k, M), b)-claim)/(z**(k-1)*(1-z)**(M-1)), force=True)))
rep('product rule D_b P_{k,M} (symbolic k,M,b)', diff == 0, diff)

# ---------- 2. V_b(z^k) from the operator equals V_{b,k} of fp-single --------------------------------------
# D_b z^k = k^2/2 z^{k-1} + [-k(k-1)/2+(7/2-2b)k+9/2-b(b-1)] z^k ; fp-single: U_k[1+((k+3)^2/2-u_{k-1})/T]+k^2/(2T)U_{k-1}
uk1 = (b+k-1)*(b+k)
lhs = -k*(k-1)/2+(sp.Rational(7, 2)-2*b)*k+sp.Rational(9, 2)-b*(b-1)
rep('diagonal coefficient of D_b z^k = (k+3)^2/2-u_{k-1}', sp.expand(lhs-((k+3)**2/2-uk1)) == 0)
rep('off-diagonal coefficient of D_b z^k = k^2/2', sp.expand(Db(z**k, b).subs(z, z) - (k**2/2*z**(k-1)+lhs*z**k)).equals(0) or
    sp.simplify(sp.expand((Db(z**k, b)-(k**2/2*z**(k-1)+lhs*z**k))/z**(k-1))) == 0)

# ---------- 3. R_-, R_+ (kernels.json) and the new inverse ----------------------------------------------------
Rm = lambda w: -2*(2*w-I)*(2*w+I)**2*(2*w+3*I)*(2*w+5*I)*(2*w+11*I)/(w*(2*w-7*I)*(2*w-3*I)*(2*w+7*I))
Rp = lambda w: -2*(2*w-11*I)*(2*w-5*I)*(2*w-3*I)*(2*w-I)**2*(2*w+I)/(w*(2*w-7*I)*(2*w+3*I)*(2*w+7*I))
rep('R_+(z)=R_-(-z) (kernels.json)', sp.simplify(Rp(z)-Rm(-z)) == 0)
Pp = lambda w: (w**2+sp.Rational(1, 4))**2*(w**2+sp.Rational(9, 4))*(w**2+sp.Rational(25, 4))*(w**2+sp.Rational(81, 4))
idt = sp.simplify(Pp(z+I)/Rm(z) - (-z*(2*z-7*I)**2*(2*z-3*I)**2*(2*z+3*I)*(2*z+7*I)**2/2048))
rep('P(z+i)/R_-(z) = -z(2z-7i)^2(2z-3i)^2(2z+3i)(2z+7i)^2/2048', idt == 0, idt)
# Fstar(z+i)/Fstar(z-i) = R_-/R_+ :  cosh(pi(z+i)/2)/cosh(pi(z-i)/2) = -1 ;
# Gamma(5/4+i(z+i)/2)Gamma(5/4-i(z+i)/2)/[Gamma(5/4+i(z-i)/2)Gamma(5/4-i(z-i)/2)]
#   = Gamma(3/4+iz/2)Gamma(7/4-iz/2)/[Gamma(7/4+iz/2)Gamma(3/4-iz/2)] = (3/4-iz/2)/(3/4+iz/2)
gr = (sp.Rational(3, 4)-I*z/2)/(sp.Rational(3, 4)+I*z/2)
ratio = sp.simplify(Pp(z+I)/Pp(z-I)*(-1)*gr - Rm(z)/Rp(z))
rep('Fstar(z+i)/Fstar(z-i) = R_-(z)/R_+(z) (Gamma recurrence + rational algebra)', ratio == 0, ratio)
rep('Gamma/cosh ratio equals (z+3i/2)/(z-3i/2) as note says', sp.simplify(-gr-(z+3*I/2)/(z-3*I/2)) == 0, sp.simplify(-gr))

# ---------- 4. raw light identities -------------------------------------------------------------------------
Q = nu**6-sp.Rational(245, 4)*nu**4+sp.Rational(9619, 16)*nu**2-sp.Rational(51575, 64)
N10 = 1024*nu**10-195328*nu**8+8229504*nu**6-90780512*nu**4+251603732*nu**2-114141375
PK = -(4*nu**2+1)*N10/(8*(4*nu**2+9)*(4*nu**2+49))          # K_B2A4 sec 4.5 : X_theta*(light) = theta chi P / Lambda^6
lhs = sp.simplify(sp.together(Rm(nu)*Q.subs(nu, nu+I)+Rp(nu)*Q.subs(nu, nu-I)))
rep('R_-Q(nu+i)+R_+Q(nu-i) = K_B2A4 P(nu)', sp.simplify(lhs-PK) == 0, sp.factor(lhs-PK))
D = lambda w: 1/(w**2+sp.Rational(9, 4))
ein = sp.simplify(Rm(nu)*D(nu+I)+Rp(nu)*D(nu-I) + 32*(16*nu**4-440*nu**2+273)/((4*nu**2+9)*(4*nu**2+49)))
rep('R_-D(nu+i)+R_+D(nu-i) = -32(16nu^4-440nu^2+273)/((4nu^2+9)(4nu^2+49))', ein == 0, ein)
# uniqueness of Q: degree-6 even Q with R_-Q(nu+i)+R_+Q(nu-i)=P.  Solve generally
cs = sp.symbols('c0:7')
Qg = sum(cs[j]*nu**j for j in range(7))
eqn = sp.numer(sp.together(Rm(nu)*Qg.subs(nu, nu+I)+Rp(nu)*Qg.subs(nu, nu-I)-PK))
sol = sp.solve(sp.Poly(sp.expand(eqn), nu).coeffs(), cs, dict=True)
rep('Q is the unique degree<=6 solution', len(sol) == 1 and sp.expand(Qg.subs(sol[0])-Q) == 0, sol)

# ---------- 5. generating polynomials g_{k,3}: coefficients of (1-z)^{-3} P_{-1/2+i nu}((1+z)/(1-z)) --------
# Pfaff: P_{-1/2+i nu}((1+z)/(1-z)) = (1-z)^{a} 2F1(a,a;1;z), a=1/2-i nu.  Independent of the note's recurrence.
a = sp.Rational(1, 2)-I*nu
Kmax = 6
ser = sp.series((1-z)**(a-3)*sum(sp.rf(a, j)**2/sp.factorial(j)**2*z**j for j in range(Kmax+1)), z, 0, Kmax+1).removeO()
g = [sp.simplify(sp.expand(ser.coeff(z, j))) for j in range(Kmax+1)]
for j in range(4):
    print('   g_%d,3 =' % j, sp.expand(g[j]))
rep('g_{k,3} real even polynomials in nu', all(sp.simplify(sp.im(sp.expand(gj.subs(nu, sp.Symbol('x', real=True))))) == 0 for gj in g[:4]))
rep('Q = -100 g_{2,3} - 36 g_{3,3}', sp.expand(Q+100*g[2]+36*g[3]) == 0, sp.expand(Q+100*g[2]+36*g[3]))
# leading coefficient of g_3 (used in my reading: only nu^6 survives when F0=F1=F2=0)
rep('[nu^6] g_{3,3} = -1/36', sp.expand(g[3]).coeff(nu, 6) == sp.Rational(-1, 36), sp.expand(g[3]).coeff(nu, 6))
# g_1: check heavy-formula/J=0 statement uses F_1 = (2/Lambda) int E (11/4 - nu^2)
rep('g_{1,3} = 11/4 - nu^2', sp.expand(g[1]-(sp.Rational(11, 4)-nu**2)) == 0, g[1])

# ---------- 6. Section-6 algebra (J=0, J=2 light checks) ---------------------------------------------------
eps, F1, F2, F3 = sp.symbols('epsilon F1 F2 F3')
def Vk(kk, bb, J, TT, Uf):
    u = lambda j: (bb+j)*(bb+j+1)
    um1 = (bb-1)*bb
    ukm1 = um1 if kk == 0 else u(kk-1)
    val = Uf(kk)*(1+(sp.Rational((kk+3)**2, 2)-ukm1)/TT)
    if kk >= 1: val += sp.Rational(kk**2, 2)/TT*Uf(kk-1)
    return val
def Ufac(bb, J, TT):
    def Uf(kk):
        if kk > J-bb: return 0
        r = sp.Integer(1)
        for j in range(kk):
            r *= (J*(J+1)-(bb+j)*(bb+j+1))/(TT-(bb+j)*(bb+j+1))
        return r
    return Uf
Fs = {0: 0, 1: eps, 2: F2, 3: F3}
V0 = sum(Vk(kk, 0, 0, T, Ufac(0, 0, T))*Fs[kk] for kk in range(2))
rep('J=0: T V_0(F_strict) = eps/2', sp.simplify(T*V0-eps/2) == 0, sp.simplify(T*V0))
U2 = Ufac(0, 2, T)(2)
V2 = sum(Vk(kk, 0, 2, T, Ufac(0, 2, T))*Fs[kk] for kk in range(4))
target = eps*(T-2)*(13*T/48+2)+F2*(3*T/2+sp.Rational(19, 2))+sp.Rational(9, 2)*F3
rep('J=2,b=0: T V/U_2 = eps(T-2)(13T/48+2)+F2(3T/2+19/2)+9F3/2', sp.simplify(T*V2/U2-target) == 0, sp.simplify(T*V2/U2-target))

print('\nSUMMARY: %d checks, %d failed' % (len(out), sum(1 for _, o in out if not o)))
