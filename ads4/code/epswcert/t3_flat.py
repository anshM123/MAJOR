"""
t0cert/t3_flat.py (T0_CERT.md §3, Part F): on every flat state with bhat = (2J+1) sqrt(u) <= B_POS (u in (0, 1]; ++: even J >= 0;
+-: J >= 4), every sign pattern g = (s2, s3) of the c-shifted NF02 action satisfies
      FF_g(c) := int (W0 + u W1) D_h(1 - 2 u p^2) dp  >=  T_abs        (W0 = 2(1-c s2) psi2, W1 = -(1-c s2) p^2 psi2 + (1-c s3) sigma psi3)
with T_abs the Part-T threshold of t1_curv (read from its json; argv[1] may override).  Unnormalised NF02 coefficients
(robust/functionals_noforward.json, as CERT_CLOSE), so the shift is k0 * T_abs with k0 = 27.48509817 (prl/data/nf02_psi.json).

  (a) J < JP: exact Lemma-P polynomial in u (rcert.poly_u), minus the shift, on u in [0, u_J], u_J = min(1, (B_POS/(2J+1))^2),
      rescaled to [0, 1]; positivity by exact Descartes/bisection (rcert.positive_on_01).
  (b) J >= JP: CERT_CLOSE/ROBUST region S (rcert: (bhat, u) Taylor/Cauchy cells, u <= ((bhat+h)/(2JP+1))^2, and the small-bhat
      series for bhat <= 1), copied with the constant shift subtracted.  nf_cert's centred interval evaluation is applied (import).
Read-only imports: robust/{nf_cert, rcert}, flatbound/common (through rcert).  Output: cache/t3_flat_<DGMIN>_<c>.{txt,json}
"""
import os, sys, json, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'robust'))
from fractions import Fraction as Fr
import flint
import nf_cert                      # read-only: patches rcert to centred interval evaluation (as in CERT_CLOSE's runs)
import rcert as R
from rcert import A, fq, iv

lines = []
def pr(s=''):
    print(s, flush=True); lines.append(s)

DGMIN = Fr(os.environ.get('T0_DGMIN', '1000'))
C = Fr(os.environ.get('T0_C', '199/1000'))
B_POS = Fr(os.environ.get('T0_BPOS', '16'))
JP = int(os.environ.get('T3_JP', '100'))
K0 = Fr('27.48509817')


def region_a(F, shift):
    psi2, psi3, a = F
    ok_all = True; cnt = 0; t0 = time.time()
    for h in ['++', '+-']:
        Js = range(0, JP, 2) if h == '++' else range(4, JP)
        for J in Js:
            uJ = min(Fr(1), (B_POS/(2*J + 1))**2)
            for g in R.gens(psi2, psi3, a, C, h, 1):
                P = R.poly_u(g, J)
                Q = [co*uJ**i for i, co in enumerate(P)]
                Q[0] -= shift
                Pz, L = R.fbc.to_fmpz_poly(Q)
                ok, n = R.positive_on_01(Pz)
                cnt += 1
                if not ok or J < 6 or J % 25 == 0:
                    pr('  (a) %s J=%d u<=%.4g deg=%d intervals=%d %s' % (g.label, J, float(uJ), Pz.degree(), n, 'PASS' if ok else 'FAIL'))
                ok_all &= ok
    pr('(a) J < %d: %d polynomials, %s (%.0fs)' % (JP, cnt, 'PASS' if ok_all else 'FAIL', time.time() - t0))
    return ok_all


def small_b_shift(g, BS, Tmax, shift):
    ok, low = R.small_b(g, BS, JP, Tmax)
    return bool(low - fq(shift) > 0), low - fq(shift)


def cell_shift(g, ak, bc, hh, Tmax, shift, kappa=3):
    """rcert.cell (region S) with the constant `shift` subtracted from the (t, u) Taylor model's constant term."""
    r = min(Fr(2), bc/2); assert hh <= r/10
    u0 = ((bc + hh)/(2*JP + 1))**2
    rho = fq(kappa*u0); rr = fq(r)
    M = R.cauchy_M(g, fq(bc), rr, rho)
    tb = R.tail_bound(g, ((bc + r)**2)/4, Tmax)
    gm = [[None]*R.KS for _ in range(R.IS)]
    for k in range(R.KS):
        pk = ak[k]; fi = A(1)
        for i in range(R.IS):
            err = A((tb/rr**i).upper())
            gm[i][k] = pk(fq(bc))/fi + (-err).union(err)
            pk = pk.derivative(); fi = fi*(i + 1)
    a2 = 2*g.aphi
    gm[0][2] += fq(a2*bc*bc/4); gm[1][2] += fq(a2*2*bc/4); gm[2][2] += fq(a2/4)
    gm[0][3] -= fq(g.aphi/2 + g.cphi)
    gm[0][0] -= fq(shift)
    q = fq(hh)/rr
    def tail_fn(ub):
        x = fq(ub)/rho
        return M*(q**R.IS + x**R.KS)/((1 - q)*(1 - x))
    stats = {'boxes': 0}
    ok = R.certify_box(gm, R.deriv_t(gm), R.deriv_u(gm), tail_fn, -hh, hh, Fr(0), u0, 0, stats)
    return ok, stats['boxes']


def region_b(F, shift, Tmax=Fr(10**6), BS=Fr(1)):
    psi2, psi3, a = F
    ok_all = True; t0 = time.time(); ncell = 0
    for h in ['++', '+-']:
        for g in R.gens(psi2, psi3, a, C, h, 1):
            flint.ctx.prec = 200
            assert R.tmax_ok(g, fq(Tmax))
            flint.ctx.prec = R.PREC_S
            ok, low = small_b_shift(g, BS, fq(Tmax), shift)
            pr('  (b) %s bhat <= %s: FF - shift >= %s %s' % (g.label, BS, low.str(5), 'PASS' if ok else 'FAIL'))
            ok_all &= ok
            ak = R.build_ak(g)
            b = BS; fails = 0
            while b < B_POS:
                r = min(Fr(2), (b + Fr(1, 10))/2); hh = min(r/10, Fr(1, 10)); bc = b + hh
                while hh > min(Fr(2), bc/2)/10:
                    hh = hh/2; bc = b + hh
                if bc + hh > B_POS:
                    hh = (B_POS - b)/2; bc = b + hh
                ok, nbox = cell_shift(g, ak, bc, hh, fq(Tmax), shift); ncell += 1
                if not ok:
                    fails += 1; pr('  (b) %s cell [%s, %s] FAIL' % (g.label, float(b), float(b + 2*hh)))
                ok_all &= ok
                b = b + 2*hh
            pr('  (b) %s: cells on [%s, %s] done, fails = %d (%.0fs)' % (g.label, BS, B_POS, fails, time.time() - t0))
    pr('(b) J >= %d, bhat <= %s: %d cells, %s (%.0fs)' % (JP, B_POS, ncell, 'PASS' if ok_all else 'FAIL', time.time() - t0))
    return ok_all


def main():
    tag = f'{DGMIN}_{str(C).replace("/", "_")}'
    if len(sys.argv) > 1:
        Tabs = Fr(sys.argv[1])
    else:
        d = json.load(open(os.path.join(HERE, 'cache', f't1_curv_{tag}.json')))
        tx = d['Tabs_upper'].strip('[]').split('+/-')         # arb string 'mid +/- rad' -> rational upper bound
        Tabs = Fr(tx[0].strip()) + (Fr(tx[1].strip()) if len(tx) > 1 else 0) + Fr(1, 10**15)
    shift = Fr(K0*Tabs).limit_denominator(10**12)
    if shift < K0*Tabs:
        shift += Fr(1, 10**12)
    pr(__doc__.split('Read-only')[0])
    pr(f'DGMIN = {DGMIN}, c = {C}, B_POS = {B_POS}, JP = {JP}; T_abs <= {float(Tabs):.6e}; shift (unnormalised) = {float(shift):.6e}')
    F = nf_cert.load('NF02')
    okb = region_b(F, shift)
    oka = region_a(F, shift)
    v = oka and okb
    pr(f'\nVERDICT t3_flat (DGMIN = {DGMIN}, c = {C}): Part F {"PASS" if v else "FAIL"}')
    json.dump(dict(DGMIN=str(DGMIN), c=str(C), Tabs=str(Tabs), JP=JP, pass_F=v), open(os.path.join(HERE, 'cache', f't3_flat_{tag}.json'), 'w'))
    open(os.path.join(HERE, 'cache', f't3_flat_{tag}.txt'), 'w').write('\n'.join(lines) + '\n')


if __name__ == '__main__':
    main()
