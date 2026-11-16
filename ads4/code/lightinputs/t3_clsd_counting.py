"""t3: Regge counting of the CLS-D spinning subtraction factors (arXiv:2311.04271, eqs. (87),(88),(106)-(111)).

[Q] At the scattering crystal (s = m^2, u = -nu^2):
  f^+_{k;Je} = -nu^{2(k-2Je+1)} (2m^2-nu^2)/(m^2(m^2-nu^2)) / (m^2(m^2-nu^2))^{(k-2Je)/2},   k even,
  f^-_{k;Je} = -nu^{2(k-2Je)}  nu^2/(m^2(m^2-nu^2))           / (m^2(m^2-nu^2))^{(k-2Je-1)/2}, k odd.
CLS-D read the change of Regge spin of the kernel off the large-m power (their text after (88)): m^{-2n} <-> +n.
[D] Both f^+ and f^- give n = k - 2Je + 1.  CLS-D's text says k - 2Je for f^-; with their own (108) it is k - 2Je + 1.
The bare light-ray commutator converges for J0 < J1 + J3 - 1 = 2Je - 1 (KKSDZ), hence with the subtraction factor:
J0 < 2Je - 1 + (k - 2Je + 1) = k for BOTH parities (CLS-D state k - 1 > J0 for odd j').
Cross-check: Je = 2, k = 3: f^- = -1 (their (111)), so the k = 3 odd graviton functional IS the bare commutator,
which converges for J0 < 3 = k (KKSDZ), not J0 < 2.
"""
import sympy as sp
m, nu = sp.symbols('m nu', positive=True)

def fplus(k, Je):
    return -nu**(2*(k-2*Je+1))*(2*m**2-nu**2)/(m**2*(m**2-nu**2))/(m**2*(m**2-nu**2))**sp.Rational(k-2*Je, 2)
def fminus(k, Je):
    return -nu**(2*(k-2*Je))*nu**2/(m**2*(m**2-nu**2))/(m**2*(m**2-nu**2))**sp.Rational(k-2*Je-1, 2)

def mpower(expr):
    # leading power of m at large m
    t = sp.symbols('t', positive=True)
    e = sp.simplify(expr.subs(m, 1/t))
    lead = sp.limit(sp.log(sp.Abs(e))/sp.log(t), t, 0)
    return -lead   # expr ~ m^{power}

print("t3: Regge-spin increase n of the CLS-D kernels (expr ~ m^{-2n}); convergence needs J0 < (2Je-1) + n")
print(" Je  k  parity   n(from (108))   CLS-D text n   convergence J0 <   unsubtracted?")
for Je in [0, 1, 2]:
    for k in range(2, 8):
        if Je == 0 and k % 2:
            continue
        f = fplus(k, Je) if k % 2 == 0 else fminus(k, Je)
        n = -mpower(f)/2
        text = (k + 1 if Je == 0 else k - 2*Je + 1) if k % 2 == 0 else k - 2*Je
        base = 2*Je - 1 if Je > 0 else -1
        # unsubtracted: no inverse power of (u'v'), i.e. exponent (k-2Je+2)/4 resp (k-2Je+1)/4 ... <= 0
        expo = sp.Rational(k-2*Je+2, 4) if (k % 2 == 0 and (k-2*Je-2) % 4 == 0) else \
               sp.Rational(k-2*Je+4, 4) if k % 2 == 0 else \
               sp.Rational(k-2*Je+1, 4) if (k-2*Je-3) % 4 == 0 else sp.Rational(k-2*Je+3, 4)
        print(f" {Je}   {k}  {'even' if k%2==0 else 'odd '}      {n}               {text}             {base + n}"
              f"            {'yes' if expo <= 0 else 'no'} ((u'v')^-{expo})")
print("Je=2: f^+_{2;2}|crystal =", sp.simplify(fplus(2, 2)), ";  f^-_{3;2}|crystal =", sp.simplify(fminus(3, 2)), " (CLS-D (111))")
