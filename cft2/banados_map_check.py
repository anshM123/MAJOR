"""Symbolic check of the Banados <- Poincare map (Roberts 1204.1982 form).
Poincare:  ds^2 = (du^2 + dw dwb)/u^2        (w, wb independent light-cone coordinates)
Banados:   ds^2 = dz^2/z^2 + L dx^2 + Lb dxb^2 + (1/z^2 + z^2 L Lb) dx dxb
Map:  w  = f  - 2 z^2 f'^2 fb''/(4 f' fb' + z^2 f'' fb'')
      wb = fb - 2 z^2 fb'^2 f''/(4 f' fb' + z^2 f'' fb'')
      u  = 4 z (f' fb')^{3/2}/(4 f' fb' + z^2 f'' fb'')
Claim: pull-back is Banados with L = -(1/2){f,x}, Lb = -(1/2){fb,xb}.  Checked for concrete f, fb."""
import sympy as sp
z, x, xb = sp.symbols('z x xb', positive=True)
def check(f, fb):
    f1, f2 = sp.diff(f, x), sp.diff(f, x, 2); g1, g2 = sp.diff(fb, xb), sp.diff(fb, xb, 2)
    D = 4*f1*g1 + z**2*f2*g2
    w = f - 2*z**2*f1**2*g2/D; wb = fb - 2*z**2*g1**2*f2/D; u = 4*z*(f1*g1)**sp.Rational(3, 2)/D
    X = [z, x, xb]
    J = lambda F: [sp.diff(F, v) for v in X]
    Ju, Jw, Jwb = J(u), J(w), J(wb)
    g = sp.zeros(3)
    for i in range(3):
        for j in range(3):
            g[i, j] = (Ju[i]*Ju[j] + sp.Rational(1, 2)*(Jw[i]*Jwb[j] + Jw[j]*Jwb[i]))/u**2
    schw = lambda F, v: sp.diff(F, v, 3)/sp.diff(F, v) - sp.Rational(3, 2)*(sp.diff(F, v, 2)/sp.diff(F, v))**2
    L, Lb = -schw(f, x)/2, -schw(fb, xb)/2
    target = sp.Matrix([[1/z**2, 0, 0], [0, L, (1/z**2 + z**2*L*Lb)/2], [0, (1/z**2 + z**2*L*Lb)/2, Lb]])
    # numeric spot check at random points (symbolic simplification is slow)
    import random
    worst = 0
    for _ in range(5):
        sub = {z: random.uniform(0.05, 0.4), x: random.uniform(-0.5, 0.5), xb: random.uniform(-0.5, 0.5)}
        diff = (g - target).subs(sub).evalf(30)
        worst = max(worst, max(abs(d) for d in diff))
    return worst
print("f=exp(2x), fb=tan(xb):", check(sp.exp(2*x), sp.tan(xb)))
print("f=x+0.3 sin x, fb=exp(xb)+xb:", check(x + sp.Rational(3, 10)*sp.sin(x), sp.exp(xb) + xb))
