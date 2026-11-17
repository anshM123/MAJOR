"""Referee r2: refine T3 of r1 at nu=20 (D1^2 flip/pres, a=9) with finer quadrature."""
import sys
from r1_direct import response, WD
nu, a = 20.0, 9
for Nu, Neta, npan in ((4096, 300, 20),):
    pD = (response(2, a, nu, WD, 1, -1, Nu=Nu, Neta=Neta, npan=npan) + response(2, a, nu, WD, -1, 1, Nu=Nu, Neta=Neta, npan=npan))/2
    fD = (response(2, a, nu, WD, 1, 1, Nu=Nu, Neta=Neta, npan=npan) + response(2, a, nu, WD, -1, -1, Nu=Nu, Neta=Neta, npan=npan))/2
    N1 = (4*nu**2 + 137)*(4*nu**2 + 561)
    s = f"r2 nu=20 Nu={Nu} Neta={Neta}x{npan}: pres {pD:.6e} claimed {N1/(4*(4*nu**2+81)*(4*nu**2+169)):.6e}; |flip/pres| {abs(fD/pD):.5f} claimed 0.69993"
    print(s); open(sys.argv[1], 'w').write(s + "\n")
