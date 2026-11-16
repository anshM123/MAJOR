"""t4: (a) leverage of hypothetical O(G) Einstein residues that (F) alone does not exclude (referee M1);
(b) the bound with a reduced smearing scale Lam_eff = Lam/n when heavy positivity is only available for m >= n Lam (referee M2);
(c) size of the linear-t4 term against the retained C1/Lam (referee m1).

(a) Test residues, each with coefficient 16 pi G (same as the Einstein pole), decaying fast enough to be allowed by (F):
    C+_2 : r(nu) = (nu^2+9/4)^-2        enters  Lam^4 int psi2 r
    C-_3 : h(nu) = (nu^2+9/4)^-3        enters  (1/2) Lam^6 int psi3 h
    C+_4 : h(nu) = (nu^2+9/4)^-4        enters  5 Lam^10 |d_{nu^2} h|_0
    C-_5 : h(nu) = (nu^2+9/4)^-5        enters  (5/2) Lam^10 h(0)
    printed as ratio to the leading Einstein term 16 pi G Lam^2 L_AdS(Lam).
"""
import mpmath as mp

psi2 = lambda p: (1 - p) ** 3 * p * (65 + 155 * p + 47 * p ** 2)
psi3 = lambda p: 2 * (1 - p) ** 3 * (81 - 1927 * p + 4601 * p ** 2)
k0 = mp.mpf(330151) / 12012
A = 1040 / k0
B = 16 * (65 * mp.log(mp.mpf(2) / 3) - mp.mpf(1171) / 15) / k0
C1 = 5765760 * mp.pi / 330151


def LAdS(Lam):
    a = mp.mpf(1.5) / Lam
    return mp.quad(lambda p: psi2(p) / (p * p + a * a), [0, a, 10 * a, 1])


def main():
    print('=== (a) leverage of residues not excluded by (F) (coefficient 16 pi G each), ratio to 16 pi G Lam^2 L_AdS ===')
    for Dg in [20, 100, 1000]:
        Lam = mp.sqrt((Dg - 2) * (Dg - 1))
        lead = Lam ** 2 * LAdS(Lam)
        c2 = Lam ** 4 * mp.quad(lambda p: psi2(p) * (Lam ** 2 * p ** 2 + 2.25) ** -2, [0, 1.5 / Lam, 15 / Lam, 1])
        c3 = 0.5 * Lam ** 6 * mp.quad(lambda p: psi3(p) * (Lam ** 2 * p ** 2 + 2.25) ** -3, [0, 1.5 / Lam, 15 / Lam, 1])
        c4 = 5 * Lam ** 10 * 4 * mp.mpf(2.25) ** -5        # |d/d(nu^2) (nu^2+9/4)^-4| at 0 = 4 (9/4)^-5
        c5 = 2.5 * Lam ** 10 * mp.mpf(2.25) ** -5
        print('Delta_g=%5d Lam=%8.2f  C+2 r: %.3e   C-3 h: %.3e   C+4 fwd: %.3e   C-5 fwd: %.3e'
              % (Dg, Lam, c2 / lead, c3 / lead, c4 / lead, c5 / lead))
    print('  => a C-_3 residue must be < ~1e-3/Lam^3 of the pole coefficient, forward residues < ~Lam^-8, to be harmless:')
    print('     only an exact zero (the spinning carry-over (E) of CHMRS (eftformula)) is usable.')

    print('\n=== (b) bound with Lam_eff = Lam/n (heavy positivity assumed only for m >= n Lam_eff) ===')
    print('  sqrt(t4^2+alpha_T^2) <= 240 sqrt(A log Lam_eff + B + C1/Lam_eff) / Lam_eff^4   (nan: right side <= 0)')
    ns = [1, 6.3, 7.4, 8.1, 22]
    print('  %8s ' % 'Delta_g' + ''.join('   n=%-6s' % n for n in ns))
    for Dg in [50, 100, 300, 1000, 10000]:
        Lam = mp.sqrt((Dg - 2) * (Dg - 1))
        row = []
        for n in ns:
            Le = Lam / n
            rhs = A * mp.log(Le) + B + C1 / Le
            row.append('%10.2e' % (240 * mp.sqrt(rhs) / Le ** 4) if rhs > 0 else '%10s' % 'nan')
        print('  %8d ' % Dg + ' '.join(row))
    print('  (n^4 scaling up to the log; Lam_eff >~ 20 is also needed for the O(Lam_eff^-2) light-side terms to be small)')

    print('\n=== (c) linear t4 term: x -> x - 431 c1 cos(phi)/Lam^2 ; x^2 shift vs C1/Lam ===')
    for Dg in [50, 100, 300, 1000]:
        Lam = mp.sqrt((Dg - 2) * (Dg - 1))
        x0 = mp.sqrt(A * mp.log(Lam) + B + C1 / Lam)
        dx = 23690 / (2 * k0) / Lam ** 2
        print('Delta_g=%5d  x0=%.3f  |dx|/|c1| = %.3e  |d(x^2)|/|c1| = %.3f   C1/Lam = %.3f' % (Dg, x0, dx, 2 * x0 * dx, C1 / Lam))


if __name__ == '__main__':
    main()
