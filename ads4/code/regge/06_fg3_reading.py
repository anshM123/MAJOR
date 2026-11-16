r"""Check 6: which reading of arXiv:2201.06602 (3.3a) reproduces the positivity claimed there (their Fig. 3)?
Heavy densities (2.33), (2.35a) with stripped partial waves (A.3); units M = 1.
A: (1-p)^3 p only on the B_2 polynomial;  B: (1-p)^3 p on both;  C: (1-p)^3 on both, p only on B_2.
Only C is non-negative on all sampled (m, J) in both channels, so the B_3 weight is 2 (1-p)^3 (81 - 1927 p + 4601 p^2)."""
import mpmath as mp
mp.mp.dps = 20
def F(m, J, ch, reading):
    a = 0 if ch == '++' else 4
    d = lambda p: mp.hyp2f1(a-J, J+a+1, 1, p*p/m**2)
    sgn = 1 if ch == '++' else -1
    psi2 = lambda p: (1-p)**3*p*(65+155*p+47*p*p)
    if reading == 'A': psi3 = lambda p: 2*(81-1927*p+4601*p*p)
    else: psi3 = lambda p: (1-p)**3*p*2*(81-1927*p+4601*p*p)
    I = mp.quad(lambda p: (psi2(p)*(2*m*m-p*p) + sgn*psi3(p))*d(p)/m**6, mp.linspace(0,1,max(4,int(J/m)+4)))
    dp1 = (J-a)*(J+a+1)/mp.mpf(2)
    return I + 5*(4*dp1-1)/m**10
for reading in 'AB':
    print("reading", reading)
    for ch in ('++','+-'):
        neg = []
        for m in (1, 1.05, 1.2, 1.5, 2, 3, 5):
            for J in ([0,2,4,6,8,10,14,20,30,40,60] if ch=='++' else [4,5,6,7,8,10,14,20,30,40,60]):
                v = F(mp.mpf(m), J, ch, reading)
                if v < 0: neg.append((m, J, mp.nstr(v*mp.mpf(m)**10, 4)))
        print("  ", ch, "negative cases:", neg[:12], "count", len(neg))
def F2(m, J, ch, psi3):
    a = 0 if ch == '++' else 4
    d = lambda p: mp.hyp2f1(a-J, J+a+1, 1, p*p/m**2)
    sgn = 1 if ch == '++' else -1
    psi2 = lambda p: (1-p)**3*p*(65+155*p+47*p*p)
    I = mp.quad(lambda p: (psi2(p)*(2*m*m-p*p) + sgn*psi3(p))*d(p)/m**6, mp.linspace(0,1,max(4,int(J/m)+4)))
    dp1 = (J-a)*(J+a+1)/mp.mpf(2)
    return I + 5*(4*dp1-1)/m**10
print("reading C: psi3 = (1-p)^3 * 2(81-1927p+4601p^2)")
psi3C = lambda p: (1-p)**3*2*(81-1927*p+4601*p*p)
for ch in ('++','+-'):
    neg = []; mn = None
    for m in (1, 1.02, 1.05, 1.1, 1.2, 1.35, 1.5, 2, 3, 5, 8, 16):
        for J in ([0,2,4,6,8,10,14,20,30,40,60,100] if ch=='++' else [4,5,6,7,8,9,10,14,20,30,40,60,100]):
            v = F2(mp.mpf(m), J, ch, psi3C)*mp.mpf(m)**10
            if mn is None or v < mn[0]: mn = (v, m, J)
            if v < 0: neg.append((m, J, mp.nstr(v, 4)))
    print("  ", ch, "negative cases:", neg[:12], "count", len(neg), " min m^10 F:", mp.nstr(mn[0], 5), "at", mn[1:])
