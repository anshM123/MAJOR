"""
gapscale/gs_cover.py (GAP_SCALING.md §5): coverage bookkeeping of the two scaling charts at a gap Delta' > 1000, near the gap
(m <= 1.3 Lam'), with exact integer twist edges.
  * chart A (fixed J, fixed u = Lam^2/m^2): the Delta' state (J, u) has a Delta_g = 1000 counterpart iff u^-1 >= 1 + 2J/999
    (i.e. m/Lam >= m_gap(1000, J)/Lam(1000)); its evidence transfers under H_mono-J (corrections ~ Lam^-2).
  * Theorem L' floor at Delta' (certified, gs_lrball): b >= 2 Lam' sinh(beta0'/2), beta0' = the certified threshold of the ball
    containing Lam' (r < 1), i.e. J >= J*'(m) = ceil((2 sinh(beta0'/2) m - 1)/2).
  * SEAM: near-gap states at Delta' with neither (they need per-gap work).
The beta0' values are read from cache/lrball_part1.log and cache/lrball_resume.log (written by gs_lrball.py).
Usage: python3 gs_cover.py -> cache/cover.txt
"""
import os, re, math
HERE = os.path.dirname(os.path.abspath(__file__))
C = os.path.join(HERE, 'cache')


def balls():
    out = []
    for fn in ('lrball_part1.log', 'lrball_resume.log'):
        p = os.path.join(C, fn)
        if not os.path.exists(p):
            continue
        for line in open(p):
            m = re.match(r'\s*([\d.]+)\s+([\d.]+) \|.*\| (\S+) \| (\S+)$', line)
            if m and m.group(3) != 'None':
                a, b = float(m.group(1)), float(m.group(2))
                num, _, den = m.group(3).partition('/')
                out.append((a, b, float(num)/float(den or 1)))
    return out


def beta0_of(Lam, B):
    c = [b0 for a, b, b0 in B if a <= Lam <= b]
    return max(c) if c else None


def main():
    B = balls()
    L1 = math.sqrt(998*999)
    lines = ["# gs_cover: near-gap (m <= 1.3 Lam') coverage at Delta' by chart A (counterpart at 1000) and Theorem L' floor (J >= J*')"]
    for Dg in (1500, 2000, 3000, 5000, 10000):
        Lam = math.sqrt((Dg - 2)*(Dg - 1))
        b0 = beta0_of(Lam, B)
        if b0 is None:
            lines.append(f"Dg={Dg}: (reference gap)"); continue
        s = 2*math.sinh(b0/2)
        seamJ = []; seam_states = 0; area = 0.0
        jcap = int(0.345*(Dg - 1))
        for J in range(0, jcap + 1):
            mg = math.sqrt((Dg - 2)*(Dg - 1 + 2*J))/Lam          # twist edge at Delta'
            mA = math.sqrt(1 + 2*J/999) if Dg > 1000 else mg       # chart-A counterpart needs m/Lam >= mA
            # Theorem L' covers m/Lam <= (2J+1)/(s Lam)   (b >= s Lam  <=>  (2J+1) Lam/m >= s Lam)
            mL = (2*J + 1)/(s*Lam)
            lo, hi = mg, min(mA, 1.3)
            if J >= 26 or Dg == 1000:
                pass
            # uncovered: [mg, min(mA, 1.3)) minus (m <= mL)
            ulo = max(lo, mL)
            if hi > ulo and hi > lo:
                seamJ.append(J); area += (hi - ulo)
        lines.append(f"Dg={Dg}: Lam'={Lam:.1f}, certified beta0'={b0:.4f} (b >= {s:.4f} Lam'), J_cap={jcap}; J*'(twist edge) ~ "
                     f"{next(J for J in range(10**6) if 2*J + 1 >= s*math.sqrt((Dg - 2)*(Dg - 1 + 2*J)))} (at 1.3 Lam': "
                     f"{next(J for J in range(10**6) if 2*J + 1 >= s*1.3*Lam)}); seam spins: "
                     + (f"{len(seamJ)} spins, J = {seamJ[0]}..{seamJ[-1]}, total m/Lam-width {area:.3f}" if seamJ else "none"))
        if seamJ:
            ex = [J for J in (5, 26, 50, 100, 150, 200, 250, 300) if J in seamJ]
            for J in ex:
                mg = math.sqrt((Dg - 2)*(Dg - 1 + 2*J))/Lam; mA = math.sqrt(1 + 2*J/999); mL = (2*J + 1)/(s*Lam)
                lines.append(f"     J={J}: twist edge {mg:.4f}, chart A from {mA:.4f}, L' up to {mL:.4f} -> uncovered m/Lam in "
                             f"[{max(mg, mL):.4f}, {min(mA, 1.3):.4f})")
    open(os.path.join(C, 'cover.txt'), 'w').write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == '__main__':
    main()
