"""Result Q, Stage 2 (RESULT_Q_STAGE12.md §2.4): leverage of light-side inputs on the g4 bound, before and after the exact light side.

Bound (flat units M = 1; Stage 0): with F = int psi2 B2 + int psi3 B3 - a (dB4|_0 - B5|_0) - c B4|_0, c = 1/2, heavy F >= 0:
   g4/(8 pi G) <= 2L - kappa |g3|^2 - a (d_dB4 - d_B5)/(8 pi G) - d_B4/(16 pi G)
where d_X is any EXTRA term in the light side -X|_low (beyond the flat values 2g4, 4 pi G|g3|^2 + g5, g5).
AdS-unit -> flat-unit dictionary (Lam = M R): G_AdS = G_f/Lam^2, |g3|^2_AdS = |g3_f|^2/Lam^8, g_n,AdS = g_n,f/Lam^{2n-2}... and the
flat-normalised lifts Lam^8 C4, Lam^10 d_{nu^2}C4 (d_{p^2} = Lam^2 d_{nu^2}), Lam^10 C5.
Items:
 (i)   test residue r 16 pi G (nu^2+9/4)^{-4} in -C4(0)             -> shift  r (4/9)^4 Lam^6
 (ii)  test residue r 16 pi G (nu^2+9/4)^{-4} in -d_{nu^2}C4 (i.e. derivative of (i)) -> shift a r 8 (4/9)^5 Lam^8
 (iii) test residue r 16 pi G (nu^2+9/4)^{-5} in -C5(0)             -> shift a r 2 (4/9)^5 Lam^8
 (iv)  test linear-t4 term r 16 pi G t4 in -C4(0), t4 <= 240 x/Lam^4 (x^2 = |g3|^2 Lam^8)  -> shift r 240 x Lam^2
 (v)   exact-light-side residual b0 G|g3|^2 in -C4(0) (b0 an O(1) AdS number, not computed) -> shift b0 x^2/(16 pi Lam^2)
 (vi)  contact contamination b5 g5 in -C4(0) (|g5| <= g4 M^-2)       -> relative shift <= b5/(2 Lam^2)
Output: cache/stage12_leverage.txt.  Instant.
"""
import os, io, json, math

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, 'cache')
buf = io.StringIO()


def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); buf.write(s + "\n")


def main():
    F = json.load(open(os.path.join(CACHE, 'functionals_q.json')))
    q = 4. / 9.
    D_LA = {100: 7508.5, 300: 10734.6, 1000: 14339.6}          # Result D: 16 L_A (|g3|^2 Lam^8 <= 16 L_A / k(0))
    pr("# stage12_leverage: shift of the bound on g4 Lam^6/(8 pi G) per unit of each light-side input (flat units M = 1)")
    for nm in ('Q0-full-robust', 'Q0-min-robust (no B3)'):
        f = F[nm]; a, A, B = f['a'], f['A'], f['B']
        pr(f"## {nm}: a = {a:.4f}, c = 1/2, A = {A:.3f}, B = {B:.3f}")
        pr("   Lam    signal   (i) C4 res   (ii) dC4 res   (iii) C5 res   (iv) t4-lin(x=10)   (v) b0, x^2=110 / D-bound   (vi) b5 rel")
        for Lam in (100, 300, 1000):
            S = A * math.log(2 * Lam / 3) + B
            i1 = q ** 4 * Lam ** 6
            i2 = a * 8 * q ** 5 * Lam ** 8
            i3 = a * 2 * q ** 5 * Lam ** 8
            i4 = 240 * 10 * Lam ** 2
            i5a = 110 / (16 * math.pi * Lam ** 2)
            i5b = D_LA[Lam] / (16 * math.pi * Lam ** 2)
            i6 = 1 / (2 * Lam ** 2)
            pr(f"  {Lam:5d}  {S:7.2f}   {i1:.2e}     {i2:.2e}      {i3:.2e}       {i4:.2e}           {i5a:.1e} / {i5b:.1e}          {i6:.1e}")
    pr("   (i)-(iv): what an UNKNOWN O(G) or O(G t4) residue of unit strength would do (Stage-0 'leverage'); entries are absolute shifts of")
    pr("   the bound on g4 Lam^6/(8 pi G) per unit r; the signal itself is 35-80.  (v), (vi): per unit b0, b5 (natural size O(1)-O(10^2)).")
    pr("   With the exact light side (RESULT_Q_STAGE12 §2): (i)-(iv) are multiplied by r = 0 exactly (Lemma A); only (v), (vi) remain.")
    # exact AdS IR log from the exact graviton pole of the k = 2 piece (Result D: -C'_EE = 16 pi G/(nu^2 + 9/4), nu = Lam p)
    import numpy as np
    x, w = np.polynomial.legendre.leggauss(200)
    pr("## exact AdS graviton-pole integral (replaces Stage 0's [M] translation log(M/m_IR) -> log(2 Lam/3)):")
    pr("   g4 Lam^6/(8 pi G) <= 2 L_AdS(Lam) - kappa x^2 - (b0 terms),  L_AdS = int_0^1 psi2(p) dp/(p^2 + 9/(4 Lam^2))  [P given the pole (D)]")
    pr("   functional                    Lam    2 L_AdS (exact)   A log(2Lam/3) + B   diff     2zeta(3)/(2 L_AdS)")
    for nm in ('Q0-full-robust', 'Q0-min-robust (no B3)'):
        f = F[nm]; c2 = np.array(f['psi2'], float)
        for Lam in (100, 300, 1000):
            # split [0,1] at a few multiples of the pole scale for accuracy
            e = 1.5 / Lam
            edges = [0, e / 10, e, 10 * e, 100 * e, 1.0] if 100 * e < 1 else [0, e / 10, e, 10 * e, 1.0]
            edges = sorted(set(min(1.0, v) for v in edges))
            L = 0.
            for lo, hi in zip(edges[:-1], edges[1:]):
                pp = lo + (hi - lo) * (x + 1) / 2
                L += ((hi - lo) / 2 * w * np.polynomial.polynomial.polyval(pp, c2) / (pp ** 2 + 9 / (4 * Lam ** 2))).sum()
            Sm = f['A'] * math.log(2 * Lam / 3) + f['B']
            pr(f"   {nm:28s} {Lam:5d}   {2 * L:12.4f}      {Sm:12.4f}      {2 * L - Sm:+.4f}   {2 * 1.2020569031595942 / (2 * L):.4f}")
    pr("## (vii) Result D's open linear-t4 term in the k = 2 piece, -C'_EE ⊃ 16 pi G c1 t4 (nu-independent, c1 not computed [I]):")
    pr("   flat units: 16 pi G_f Lam^2 t4 c1 int psi2, t4 = 240 x/Lam^4  => shift of the g4 bound = 480 |c1| x int(psi2)/Lam^2")
    for nm in ('Q0-full-robust', 'Q0-min-robust (no B3)'):
        c2 = np.array(F[nm]['psi2'], float)
        Ipsi = float(np.polynomial.polynomial.polyval(1.0, np.polynomial.polynomial.polyint(c2)))
        row = ", ".join(f"Lam={Lam}: {480 * 10.5 * Ipsi / Lam ** 2:.2e} (x=10.5) / {480 * math.sqrt(D_LA[Lam]) * Ipsi / Lam ** 2:.2e} (x=D-bound)" for Lam in (100, 300, 1000))
        pr(f"   {nm}: int psi2 = {Ipsi:.4f};  per unit c1: {row}")
    open(os.path.join(CACHE, 'stage12_leverage.txt'), 'w').write(buf.getvalue())


if __name__ == '__main__':
    main()
