"""Target 2 (side): J*(m) table of LARGESPIN_L §3.3: smallest J with beta >= 1, sinh(beta/2) = (J+1/2)/m,
m^2 = (Delta-J-2)(Delta+J-1) (CHMRS (d = 3); Route W's cost cache uses the same), gap: Delta - J = 1000.  Output: cache/t2b_jstar_m.txt"""
import math
Lam = math.sqrt(998 * 999)
s1 = 2 * math.sinh(0.5)
def mgap(J, tau=1000):
    D = J + tau
    return math.sqrt((D - J - 2) * (D + J - 1))
out = open('cache/t2b_jstar_m.txt', 'w')
def log(s): print(s); out.write(s + '\n')
log(f'2 sinh(1/2) = {s1:.7f}; Lam = {Lam:.4f}')
for c in (1.0, 1.1, 1.3):
    J = 0
    while (2 * J + 1) < s1 * c * mgap(J):
        J += 1
    beta = 2 * math.asinh((J + 0.5) / (c * mgap(J)))
    log(f'm = {c} m_gap(J): J* = {J} (beta = {beta:.5f}; J*-1 gives beta = {2*math.asinh((J-0.5)/(c*mgap(J-1))):.5f})')
for c in (1.5, 2.0, 5.0):
    m = c * Lam
    J = math.ceil((s1 * m - 1) / 2)
    log(f'm = {c} Lam: J* = {J}   (with the rounded constant 1.0421: {math.ceil((1.0421*m-1)/2)})')
for (beta0, claim) in ((1.5, (1737, 2034, 2702, 1232, 1642, 4105)), (2.0, (3188, 3776, 5113, 1760, 2347, 5867))):
    s = 2 * math.sinh(beta0 / 2); res = []
    for c in (1.0, 1.1, 1.3):
        J = 0
        while (2 * J + 1) < s * c * mgap(J): J += 1
        res.append(J)
    for c in (1.5, 2.0, 5.0):
        res.append(math.ceil((s * c * Lam - 1) / 2))
    log(f'J at beta >= {beta0}: {res}  claimed {claim}')
