"""Flat certificate for Result Q's upper functional Q3-full-delta0.2 (CERT_CLOSE.md §2): FF >= DELTA N on every flat heavy state.

Flat action (prl/data/q3_full_delta0.2.json, RESULT_Q_STAGE0 / common_q.heavy_action), u = 1/m^2, JJ = J(J+1):
    FF = int psi2 (2 - u p^2) D_h + s u int psi3 D_h + a (2 JJ - c4 + s) u^3 - 2 c u^2,     N = |B2| + |B3| + |FW| + 2 c u^2,
with psi2, psi3 the stored IEEE doubles read as EXACT binary rationals, a, c likewise.  FF - DELTA N >= 0  <=>  for all
(s2, s3) in {+-1}^2 the sign-pattern action of rcert.gens (forward sign sF = sign(2JJ - c4 + s)) minus 2 c (1 + DELTA) u^2 is >= 0.
This is ROBUST_FUNCTIONAL.md's certificate (rcert.py, read-only import) with one extra term, -cq u^2, cq = 2c(1 + DELTA), inserted:
  P   J < JP, all u in [0,1]         Lemma-P polynomial + (-cq) u^2, exact Descartes
  C   J >= JP, u >= 1/4             rcert.region_C with Phi -> Phi - 4 cq   (m^2 <= 4)
  M   u <= 1/4, bhat >= BA          rcert.region_M with eta -> eta + cq/(aphi bhat^2/2)
  S   J >= JP, u <= 1/4, bhat <= BA rcert.region_S (Cauchy cells) with the u^2 Taylor coefficient lowered by cq; small-b bound lowered by cq u^2
Usage: python3 q_cert.py DELTA REGION [args]     REGION in dense, P, C, M, S++, S+-
"""
import sys, os, json, time, hashlib
from fractions import Fraction as Fr
import flint
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'robust'))
import nf_cert            # noqa: E402  (patches rcert to centred interval evaluation; harmless for P/C)
import rcert              # noqa: E402
from rcert import A, fq, gens   # noqa: E402

DATA = os.path.join(HERE, '..', '..', 'prl', 'data', 'q3_full_delta0.2.json')


def load_q():
    d = json.load(open(DATA))
    psi2 = [Fr(float(x)) for x in d['psi2']]; psi3 = [Fr(float(x)) for x in d['psi3']]
    a = Fr(float(d['a'])); c = Fr(float(d['c']))
    canon = 'psi2=' + ','.join(map(str, psi2)) + ';psi3=' + ','.join(map(str, psi3)) + ';a=%s;c=%s' % (a, c)
    print('FUNCTIONAL %s from prl/data/q3_full_delta0.2.json (file sha256[:16] %s); sha256(exact binary coefficients) = %s'
          % (d['name'], hashlib.sha256(open(DATA, 'rb').read()).hexdigest()[:16], hashlib.sha256(canon.encode()).hexdigest()), flush=True)
    return psi2, psi3, a, c


# ------------------------------------------------------------------ dense floating-point scan (not a certificate)
def dense(F, Jcap=3000):
    import numpy as np
    sys.path.insert(0, os.path.join(HERE, '..', 'resultq'))
    from common_q import heavy_action
    from rcommon import jmax_for
    psi2, psi3, a, c = F
    p2 = np.array([float(x) for x in psi2]); p3 = np.array([float(x) for x in psi3])
    us = np.unique(np.concatenate([np.logspace(-6, 0, 1201), np.linspace(0.005, 0.05, 901), np.linspace(0.7, 0.9, 401)]))
    best = {}
    for h in ['++', '+-']:
        for u in us:
            Js, FF, N = heavy_action(p2, p3, float(a), float(c), h, u, jmax_for(u, 400.0, Jcap))
            q = FF / N; j = int(np.argmin(q))
            if h not in best or q[j] < best[h][0]:
                best[h] = (q[j], int(Js[j]), u)
        print('dense %s: min FF/N = %.6f at J = %d, u = %.6g (m = %.4f M)' % (h, best[h][0], best[h][1], best[h][2], best[h][2] ** -0.5), flush=True)
    return best


# ------------------------------------------------------------------ P
def region_P(F, delta, JP):
    psi2, psi3, a, c = F
    cq = 2 * c * (1 + delta)
    ok_all = True; cnt = 0; t0 = time.time(); worst = None
    for h in ['++', '+-']:
        Js = range(0, JP, 2) if h == '++' else range(4, JP)
        for J in Js:
            sF = 1 if (2 * J * (J + 1) - rcert.C4[h] + rcert.SIGMA[h]) >= 0 else -1
            for g in gens(psi2, psi3, a, delta, h, sF):
                co = rcert.poly_u(g, J)
                while len(co) < 3: co.append(Fr(0))
                co[2] -= cq
                P, L = rcert.fbc.to_fmpz_poly(co)
                ok, n = rcert.positive_on_01(P)
                cnt += 1; ok_all &= ok
                if not ok or J < 8 or J % 100 == 0:
                    print('P %s J=%d deg=%d intervals=%d %s' % (g.label, J, P.degree(), n, 'PASS' if ok else 'FAIL'), flush=True)
    print('VERDICT q_cert_P (delta=%s, JP=%d, %d polynomials, %.0fs): %s' % (delta, JP, cnt, time.time() - t0, 'PASS' if ok_all else 'FAIL'))


# ------------------------------------------------------------------ C (copy of rcert.region_C with Phi -> Phi - 4 cq)
def region_C(F, delta, JP):
    from math import comb
    flint.ctx.prec = 200
    psi2, psi3, a, c = F
    cq = 2 * c * (1 + delta)
    ok_all = True
    for h in ['++', '+-']:
        for g in gens(psi2, psi3, a, delta, h, 1):
            G0 = rcert.div_omp3(g.W0); G1 = rcert.div_omp3(g.W1)
            gh = lambda x: A(rcert.ev(G0, x).abs_upper()) + A(rcert.ev(G1, x).abs_upper())
            Phi = 2 * fq(g.aphi) * JP * (JP + 1) - fq(g.cphi) - 4 * fq(cq)
            if h == '++':
                tot = rcert.integral_upper(lambda x: gh(x) * (1 - x) ** 3, 0, 1, 4000)
                ok = bool(Phi > 64 * A(tot)) and bool(4 * fq(g.aphi) * JP > 0)
                print('C %s J>=%d: 2aJJ-c-4cq = %s > 64 int(|W0|+|W1|) <= %s : %s' % (g.label, JP, Phi.str(6), (64 * A(tot)).str(6), ok))
            else:
                c_a = rcert.integral_upper(lambda x: gh(x) / ((1 - x) * (1 + x) ** 4), 0, Fr(1, 2), 2000)
                gam = rcert.sup_on(lambda x: gh(x) / (1 + x) ** 4, Fr(1, 2), 1, 4000)
                g1 = rcert.sup_on(gh, Fr(1, 2), 1, 4000)
                Bn = A(comb(JP + 4, 8)); assert bool(Bn > A(4) ** 4 / A(3) ** 4)
                U = A(c_a) + A(gam) * Bn.log() / 4 + A(g1) / (4 * (A(3) / 2) ** 4)
                ok1 = bool(Phi > 64 * U)
                mono = bool(2 * fq(g.aphi) * (2 * JP + 1) > 16 * A(gam) * 8 / A(JP - 3))
                ok = ok1 and mono
                print('C %s J>=%d: U<=%s, 2aJJ-c-4cq=%s > 64U=%s: %s ; monotone %s' % (g.label, JP, U.str(5), Phi.str(6), (64 * U).str(6), ok1, mono))
            ok_all &= ok
    print('VERDICT q_cert_C (delta=%s, JP=%d): %s' % (delta, JP, 'PASS' if ok_all else 'FAIL'))


# ------------------------------------------------------------------ M (copy of rcert.region_M with eta -> eta + cq/(aphi b^2/2))
def region_M(F, delta, BA, u1, NP):
    flint.ctx.prec = 128
    psi2, psi3, a, c = F
    cq = 2 * c * (1 + delta)
    b = fq(BA); ok_all = True
    neg = lambda x: (-x) if x < 0 else Fr(0)
    for h in ['++', '+-']:
        for g in gens(psi2, psi3, a, delta, h, 1):
            assert g.W0[0] == 0
            alpha = g.W1[0]
            q0 = 2 * (g.W0[2] if len(g.W0) > 2 else 0)
            q1 = 2 * (g.W1[2] if len(g.W1) > 2 else 0) + (16 * g.W0[0] if h == '+-' else 0)
            q2 = 16 * g.W1[0] if h == '+-' else Fr(0)
            gam0, gam1, gam2 = -q0 / 2, -q1 / 2, -q2 / 2
            assert gam0 > 0, 'leading 1/b^3 coefficient not positive'
            CR1 = rcert.CR_bound(g, 0, u1, NP); CRm = rcert.CR_bound(g, u1, Fr(1, 4), NP)
            dl = 8 * fq(u1).sqrt() / b if h == '+-' else A(0)
            bm = b * (1 - dl)
            pi1 = fq(neg(alpha)) + fq(neg(gam1)) / bm ** 2
            aph = fq(g.aphi)
            eta = (fq(max(g.aphi / 2 + g.cphi, 0)) * fq(u1) + fq(neg(gam2)) / bm ** 3 + fq(cq)) / (aph / 2 * b ** 2)
            val = fq(gam0) * b / (1 + dl) ** 3 - CR1 / (1 - dl) ** 4 - pi1 ** 2 / (2 * aph * (1 - eta) * (1 - dl) ** 2)
            c1 = bool(val > 0) and bool(eta < 1)
            bm2 = b - 4 if h == '+-' else b
            pi1m = fq(neg(alpha)) + fq(neg(gam1)) / bm2 ** 2
            etam = (fq(max(g.aphi / 2 + g.cphi, 0)) / 4 + fq(neg(gam2)) / bm2 ** 3 + fq(cq)) / (aph / 2 * b ** 2)
            A2 = aph / 2 * b ** 2 * (1 - etam); A1 = pi1m / bm2; A0 = CRm / bm2 ** 4
            uu = fq(u1)
            c2 = bool(2 * A2 * uu > A1) and bool(A2 * uu ** 2 - A1 * uu - A0 > 0) and bool(etam < 1)
            ok = c1 and c2
            print('M %s: gamma0=%.5g alpha=%.5g gamma1=%.5g aphi=%.4g; C_R(u<=u1)<=%s C_R(u1..1/4)<=%s; b^4 FF >= %s (u<=u1) %s ; u in [u1,1/4] %s'
                  % (g.label, float(gam0), float(alpha), float(gam1), float(g.aphi), CR1.str(5), CRm.str(5), val.str(5), c1, c2), flush=True)
            ok_all &= ok
    print('VERDICT q_cert_M (delta=%s, BA=%s, u1=%s): %s' % (delta, BA, u1, 'PASS' if ok_all else 'FAIL'))


# ------------------------------------------------------------------ S (copies of rcert.cell / small_b / region_S with the -cq u^2 term)
def cell(g, ak, bc, hh, JP, Tmax, cq, kappa=3):
    r = min(Fr(2), bc / 2); assert hh <= r / 10
    u0 = ((bc + hh) / (2 * JP + 1)) ** 2
    rho = fq(kappa * u0); rr = fq(r)
    M = rcert.cauchy_M(g, fq(bc), rr, rho)
    tb = rcert.tail_bound(g, ((bc + r) ** 2) / 4, Tmax)
    KS, IS = rcert.KS, rcert.IS
    gm = [[None] * KS for _ in range(IS)]
    for k in range(KS):
        pk = ak[k]; fi = A(1)
        for i in range(IS):
            err = A((tb / rr ** i).upper())
            gm[i][k] = pk(fq(bc)) / fi + (-err).union(err)
            pk = pk.derivative(); fi = fi * (i + 1)
    a2 = 2 * g.aphi
    gm[0][2] += fq(a2 * bc * bc / 4); gm[1][2] += fq(a2 * 2 * bc / 4); gm[2][2] += fq(a2 / 4)
    gm[0][3] -= fq(g.aphi / 2 + g.cphi)
    gm[0][2] -= fq(cq)                                   # the -2c(1+delta) u^2 term (exact, b-independent)
    q = fq(hh) / rr
    def tail_fn(ub):
        x = fq(ub) / rho
        return M * (q ** IS + x ** KS) / ((1 - q) * (1 - x))
    stats = {'boxes': 0}
    ok = rcert.certify_box(gm, rcert.deriv_t(gm), rcert.deriv_u(gm), tail_fn, -hh, hh, Fr(0), u0, 0, stats)
    return ok, stats['boxes'], M, u0


def small_b(g, BS, JP, Tmax, cq):
    ok, low = rcert.small_b(g, BS, JP, Tmax)
    u = fq((Fr(BS) / (2 * JP + 1)) ** 2)
    low = low - fq(cq) * u ** 2
    return bool(low > 0), low


def region_S(F, delta, h, BS, BA, JP, Tmax=Fr(10 ** 6)):
    flint.ctx.prec = rcert.PREC_S
    psi2, psi3, a, c = F
    cq = 2 * c * (1 + delta)
    ok_all = True; t0 = time.time(); ncell = 0
    for g in gens(psi2, psi3, a, delta, h, 1):
        flint.ctx.prec = 200
        assert rcert.tmax_ok(g, fq(Tmax)), 'Tmax'
        flint.ctx.prec = rcert.PREC_S
        ok, low = small_b(g, BS, JP, fq(Tmax), cq)
        print('S %s small-b (b<=%s): lower bound %s %s' % (g.label, BS, low.str(6), 'PASS' if ok else 'FAIL'), flush=True)
        ok_all &= ok
        ak = rcert.build_ak(g)
        b = Fr(BS); fails = 0
        while b < BA:
            r = min(Fr(2), (b + Fr(1, 10)) / 2); hh = min(r / 10, Fr(1, 10)); bc = b + hh
            while hh > min(Fr(2), bc / 2) / 10:
                hh = hh / 2; bc = b + hh
            if bc + hh > BA:
                hh = (BA - b) / 2; bc = b + hh
            ok, nbox, M, u0 = cell(g, ak, bc, hh, JP, fq(Tmax), cq); ncell += 1
            if not ok:
                fails += 1
                print('S %s cell b in [%s, %s] umax=%.3e M<=%s FAIL' % (g.label, float(b), float(b + 2 * hh), float(u0), M.str(4)), flush=True)
            ok_all &= ok
            b = b + 2 * hh
        print('S %s: cells done, fails=%d (%.0fs)' % (g.label, fails, time.time() - t0), flush=True)
    print('VERDICT q_cert_S %s (delta=%s, b in [%s,%s], %d cells, %.0fs): %s' % (h, delta, BS, BA, ncell, time.time() - t0, 'PASS' if ok_all else 'FAIL'))


if __name__ == '__main__':
    F = load_q()
    if sys.argv[1] == 'dense':
        dense(F, int(sys.argv[2]) if len(sys.argv) > 2 else 3000); sys.exit()
    delta, region = Fr(sys.argv[1]), sys.argv[2]
    if region == 'P': region_P(F, delta, int(sys.argv[3]) if len(sys.argv) > 3 else 400)
    elif region == 'C': region_C(F, delta, int(sys.argv[3]) if len(sys.argv) > 3 else 400)
    elif region == 'M': region_M(F, delta, Fr(sys.argv[3]) if len(sys.argv) > 3 else Fr(120), Fr(sys.argv[4]) if len(sys.argv) > 4 else Fr(1, 100),
                                 int(sys.argv[5]) if len(sys.argv) > 5 else 4000)
    elif region.startswith('S'): region_S(F, delta, region[1:], Fr(1), Fr(sys.argv[3]) if len(sys.argv) > 3 else Fr(120),
                                          int(sys.argv[4]) if len(sys.argv) > 4 else 400)
