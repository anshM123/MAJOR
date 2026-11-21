import mpmath as mp
mp.mp.dps=20
def P(nu,b): return mp.hyp2f1(mp.mpf(1)/2-1j*nu, mp.mpf(1)/2+1j*nu, 1, -mp.sinh(b/2)**2).real
def M1(b):  # NF02_LARGE_B form
    return mp.quad(lambda nu: 2*nu*P(nu,b)/(mp.e**(2*mp.pi*nu)+1), [0,2,6,15])
def K1(t): return mp.nsum(lambda m: 8*mp.pi*m*t/(4*mp.pi**2*m**2+t**2)**2,[1,mp.inf])
def Q1(b):  # Lemma Q exterior representation
    return mp.sqrt(2)/mp.pi*mp.quad(lambda t: K1(t)/mp.sqrt(mp.cosh(t)-mp.cosh(b)), [b, b+1, b+5, b+30, mp.inf])
def M1int(b): # interior representation of NF02_LARGE_B
    return mp.sqrt(2)/mp.pi*mp.quad(lambda t: (mp.cosh(t/2)/(4*mp.sinh(t/2)**2)-1/t**2)/mp.sqrt(mp.cosh(b)-mp.cosh(t)),[0,b])
for b in [0.3,1.0,2.5,6.0]:
    b=mp.mpf(b); print(float(b), mp.nstr(M1(b),12), mp.nstr(Q1(b),12), mp.nstr(M1int(b),12))
