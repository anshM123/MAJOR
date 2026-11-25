"""
Audit of ANGULAR_ROUTE §21 item 2: collinear block lower bound and exponent bookkeeping.

(a) k^{(d)}_{2h}(z) := z^h 2F1(h+d, h+d; 2h; z), the SL(2) block of <A B B^dag A^dag> with d = h_B - h_A.
    Claim: k^{(d)}_{2h}(1 - S^-2) >~ 4^h S^{1/2+2d} uniformly for h in [S, 2S].
    We print min/max over h in [S,2S] of k/(4^h S^{1/2+2d}) for S = 8 .. 4096.
(b) GFF calibration (A != B, one Wick contraction): (1-z)^{-2 h_B} z^{-(h_A+h_B)}... decomposed into
    k^{(d)} with the known coefficients c_n = (2h_A)_n (2h_B)_n / (n! (2h_A+2h_B+n-1)_n); check that
    sum_{h in [U,2U]} c 4^h  ~  U^{2(h_A+h_B) - 1/2}  (the exponent 11/2 - M of the claim).
"""
import mpmath as mp
mp.mp.dps = 30

def k(h, d, z):
    return z**h*mp.hyp2f1(h + d, h + d, 2*h, z)

print("(a) block lower bound: ratio k/(4^h S^(1/2+2d)) over h in [S,2S], z = 1 - S^-2")
for d in (-1, -0.5, 0, 0.5, 1, 2):
    row = []
    for S in (8, 32, 128, 512, 2048, 4096):
        z = 1 - mp.mpf(1)/S**2
        vals = [k(mp.mpf(S)*(1 + j/16.), d, z)/(mp.mpf(4)**(S*(1 + j/16.))*mp.mpf(S)**(0.5 + 2*d)) for j in range(17)]
        row.append((S, float(min(vals)), float(max(vals))))
    print(f"  d = {d:5}: " + "  ".join(f"S={S}: [{a:.3g},{b:.3g}]" for S, a, b in row))

print("(b) GFF calibration of the crossed-identity exponent")
for hA, hB in ((2.5, 2.5), (0.5, 0.5), (2.5, 2.0), (2.0, 2.0), (1.0, 1.0), (1.5, 0.5)):
    d = hB - hA
    out = []
    for U in (16, 64, 256, 1024):
        tot = mp.mpf(0)
        for n in range(int(U - hA - hB), int(2*U - hA - hB) + 1):
            if n < 0:
                continue
            c = mp.rf(2*hA, n)*mp.rf(2*hB, n)/(mp.factorial(n)*mp.rf(2*hA + 2*hB + n - 1, n))
            tot += c*mp.mpf(4)**(hA + hB + n)
        out.append((U, tot/mp.mpf(U)**(2*(hA + hB) - 0.5)))
    print(f"  (hA,hB)=({hA},{hB}): sum_[U,2U] c 4^h / U^(2(hA+hB)-1/2) = " + ", ".join(f"U={U}:{float(r):.4g}" for U, r in out))
