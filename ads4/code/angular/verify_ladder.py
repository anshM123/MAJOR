"""ANGULAR_ROUTE §8 check: closed-form ladders F_{-k}, F_{4+k} solve the light-transform-conservation recurrence
m(m-4)(T+C-2m^2)F_m = (m-3)(m-4)(D-m-2)(J+m+1)F_{m+1} + m(m-1)(D+m-2)(J-m+1)F_{m-1}, with F_1=F_2=F_3=0, F_{+-(J+1)}=0;
and the normalised ratio with a = F/sqrt(N_m) (N_m from §4, symmetric in m)."""
import sympy as sp
D = sp.symbols('Delta')
def check(J, Dv=None):
    Dx = D if Dv is None else Dv
    T = (Dx-1)*(Dx-2); C = J*(J+1); F = {}
    for k in range(0, J+1): F[-k] = sp.binomial(J, k)*sp.rf(4, k)/sp.rf(Dx-1-k, k)
    for k in range(0, J-3): F[4+k] = sp.binomial(J-4, k)*sp.rf(4, k)/sp.rf(Dx-5-k, k)
    for m in (1, 2, 3): F[m] = 0
    g = lambda m: F.get(m, 0)
    for m in range(-J, J+1):
        lhs = m*(m-4)*(T + C - 2*m*m)*g(m)
        rhs = (m-3)*(m-4)*(Dx-m-2)*(J+m+1)*g(m+1) + m*(m-1)*(Dx+m-2)*(J-m+1)*g(m-1)
        assert sp.simplify(lhs - rhs) == 0, (J, m)
    N = lambda m: sp.binomial(2*J, J+m)*sp.gamma(Dx+m-1)*sp.gamma(Dx-m-1)
    for b, idx in ((0, lambda k: -k), (4, lambda k: 4+k)):
        kmax = J - b
        for k in range(0, kmax):
            s = b + k
            r = (g(idx(k+1))/sp.sqrt(N(idx(k+1))))/(g(idx(k))/sp.sqrt(N(idx(k))))
            target = sp.Rational(k+4, k+1)*sp.sqrt((C - s*(s+1))/(T - s*(s+1)))
            if Dv is not None: assert abs(sp.N(r/target) - 1) < 1e-20, (J, b, k, sp.N(r/target))
for J in range(4, 11): check(J)                       # symbolic Delta
print("recurrence: closed-form ladders exact for J=4..10, symbolic Delta")
for J in (4, 7, 12, 25, 40):
    for Dv in (sp.Rational(2*J+11, 2), sp.Integer(J+1000)): check(J, Dv)
print("normalised ratio a_{k+1}/a_k = (k+4)/(k+1) sqrt(X_k/Y_k) with a = F/sqrt(N): OK (J up to 40, two Delta each)")
