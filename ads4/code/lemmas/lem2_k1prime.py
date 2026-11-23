r"""
lem2_k1prime (LEMMAS_CLOSE.md §2): K1' -- the Rindler/Cauchy-Schwarz Regge bound made uniform over the whole Regge box, with its edge
behaviour, checked on exactly known unitary correlators (scalar externals).

Frame (Caron-Huot 1703.00278 §2.2; KKSDZ 1904.05905 §4.2.3-4.2.4): x4 = (rho, rhob) = -x3, x1 = (1, 1) = -x2 in light-cone coordinates,
z = 4 rho/(1+rho)^2.  Regge sheet: 0 < rho < 1 < rhob (rhob continued around 1 with -+ i0).  Write rho_b := 1/rhob.  The Regge box of
H7 §14 / HK_CLOSE is B_c = {0 < rho <= c, 0 < rho_b <= c}; the large-eta corner (light-cone edges inside the box) is rho_b/rho -> 0
or rho/rho_b -> 0 at fixed max(rho, rho_b).
K1' (statement checked):  sup_{B_c} |G_Regge| <= G_E(c, c),   sup_{B_c} |dDisc G| <= 2 G_E(c, c),
with G_E(rho, rho_b) the Euclidean (first-sheet) correlator, monotone in each variable (positive s-channel coefficients, identical ops).
Edge behaviour: G_E(rho, rho_b) -> 1 (identity) as either variable -> 0.
Also shown: a bound |dDisc G| <= B r^{1-J0} with the TRUE intercept J0 < 1 is NOT uniform towards the edges (2d Ising sigma, J0 = 1/2);
the uniform statement holds at the Rindler value J0 = 1 (r^0), which is what K1' asserts and what H_K uses.
Examples: 2d Ising <sigma sigma sigma sigma> (exact):  G = (1 + sqrt(rho rhob))/((1 - rho^2)(1 - rhob^2))^{1/4};
GFF of dimension D (d-independent):  G = 1 + (z zb)^D + (z zb/((1-z)(1-zb)))^D.
Output: cache/lem2_k1prime.txt (seconds).
"""
import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = []


def pr(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); OUT.append(s)


def ising(rho, rhob_c):
    """sigma 4-pt, stripped, at complex (rho, rhob) with explicit branch: rhob_c complex (rhob -+ i0 handled by caller)."""
    return (1 + np.sqrt(rho*rhob_c))/((1 - rho**2)*(1 - rhob_c**2))**0.25


def ising_regge(rho, rb, sgn):
    """second sheet: rhob = 1/rb > 1 approached from Im rhob = sgn*0 (power (1 - rhob^2)^{-1/4} continued around rhob = 1)."""
    rhob = 1/rb
    ph = np.exp(-sgn*1j*np.pi/4)          # (1 - rhob^2)^{1/4} = |..|^{1/4} e^{sgn i pi/4}
    return (1 + np.sqrt(rho*rhob))/((1 - rho**2)**0.25*np.abs(1 - rhob**2)**0.25)*ph


def gff_E(rho, rb, D):
    z = 4*rho/(1 + rho)**2; zb = 4*rb/(1 + rb)**2
    return 1 + (z*zb)**D + (z*zb/((1 - z)*(1 - zb)))**D


def gff_regge(rho, rb, D, sgn):
    rhob = 1/rb
    z = 4*rho/(1 + rho)**2; zb = 4*rhob/(1 + rhob)**2     # zb analytic through rhob = 1; (1 - zb)^{-D} picks the phase
    one_m_zb = ((1 - rhob)/(1 + rhob))**2                   # = (1 - zb) on the second sheet up to the phase e^{-+ 2 pi i}
    return 1 + (z*zb)**D + (z*zb/((1 - z)*one_m_zb))**D*np.exp(sgn*2j*np.pi*D)


def check(label, GE, GR, c, J0true=None):
    g = np.logspace(-9, np.log10(c), 220)
    R, B = np.meshgrid(g, g, indexing='ij')
    up, dn = GR(R, B, +1), GR(R, B, -1)
    E = GE(R, B)
    dd = E - 0.5*(up + dn)
    bound = GE(np.array(c), np.array(c))
    pr(f'  {label}, c = {c}: G_E(c,c) = {float(np.real(bound)):.6f}')
    pr(f'     sup_box |G_Regge| / G_E(c,c) = {np.max(np.abs(up))/np.real(bound):.4f}  (must be <= 1);  pointwise |G_R| <= G_E(rho, rho_b) '
       f'everywhere: {bool(np.all(np.abs(up) <= np.real(E)*(1 + 1e-12)))}')
    pr(f'     sup_box |dDisc| / (2 G_E(c,c)) = {np.max(np.abs(dd))/(2*np.real(bound)):.4f};  min dDisc/G_E = {np.min(np.real(dd)/np.real(E)):.1e} (>= 0 up to float cancellation: Rindler positivity)')
    for e in (1e-3, 1e-6, 1e-9):
        i = np.argmin(abs(g - e))
        pr(f'     edge rho_b = {g[i]:.0e}: max_rho |G_R| = {np.max(np.abs(up[:, i])):.5f}, G_E -> {np.max(np.real(E[:, i])):.6f};  '
           f'corner rho = rho_b = {g[i]:.0e}: |G_R| = {abs(up[i, i]):.3e}, dDisc = {np.real(dd[i, i]):.6f}')
    if J0true is not None:
        r = np.sqrt(R*B)
        for J0 in (J0true, 1.0):
            q = np.abs(up)/r**(1 - J0)
            pr(f'     sup_box |G_R| / r^(1-J0), J0 = {J0}: {np.max(q):.3e}  (at rho = {R.flat[np.argmax(q)]:.1e}, rho_b = {B.flat[np.argmax(q)]:.1e})')


if __name__ == '__main__':
    pr('# lem2_k1prime: uniform Regge-box bound (K1\'), LEMMAS_CLOSE.md §2')
    for c in (0.5, 0.9):
        check('2d Ising sigma (Delta = 1/8)', lambda R, B: ising(R, B), lambda R, B, s: ising_regge(R, B, s), c, J0true=0.5)
    for D in (0.5, 1.3, 3.0):
        check(f'GFF Delta = {D}', lambda R, B: gff_E(R, B, D), lambda R, B, s: gff_regge(R, B, D, s), 0.9)
    pr('  Ising |G_R| in closed form on the box: (sqrt rho + sqrt rho_b)/((1-rho^2)(1-rho_b^2))^{1/4} <= G_E = (1 + sqrt(rho rho_b))/(...)'
       ' <=> (1 - sqrt rho)(1 - sqrt rho_b) >= 0.  At fixed rho, rho_b -> 0: |G_R| -> sqrt(rho)/(1-rho^2)^{1/4} while r^{1/2} -> 0:'
       ' the true-intercept form r^{1-J0}, J0 = 1/2, is not uniform towards the edge.')
    with open(os.path.join(HERE, 'cache', 'lem2_k1prime.txt'), 'w') as f:
        f.write('\n'.join(OUT) + '\n')
