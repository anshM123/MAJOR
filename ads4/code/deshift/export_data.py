"""deshift/export_data.py -- DESHIFT_CLOSE.md §5: write the functional data for the Supplement to ads4/prl/data/ (plain text / JSON).
Inputs are read-only: robust/functionals_noforward.json, resultq/functionals_q3.json, resultq/functionals_q7_xprec.json,
resultq/functionals_q4_deshift.json, resultq/functionals_q5_deshift.json, spinfunc/cache/u2_hclosed.pkl, spinfunc/cache/u3_chi_EE_<Dg>.pkl,
and this directory's cache/ab_<piece>_<tag>.json (continuum weights).
Usage: python3 export_data.py <tag_c4v> <tag_c4d> <tag_c5>
"""
import os, sys, json, pickle, hashlib
import numpy as np
from flint import arb, acb, ctx
import kern
import ab_check as AC

ctx.prec = 384
DATA = os.path.join(kern.ADS4, 'prl', 'data')
os.makedirs(DATA, exist_ok=True)
manifest = []


def put(name, text):
    fn = os.path.join(DATA, name)
    open(fn, 'w').write(text)
    manifest.append((name, os.path.getsize(fn), hashlib.sha256(text.encode()).hexdigest()[:16]))


def main(tags):
    C = kern.CODE
    # 1. NF02 (Result D, Omega^EE)
    nf = json.load(open(os.path.join(C, 'robust', 'functionals_noforward.json')))['NF02']
    put('nf02_psi.json', json.dumps(dict(
        name='NF02', source='ads4/code/robust/functionals_noforward.json',
        definition='psi_k(p) = (1-p)^5 p X_k(p)/k0 on 0 <= p <= 1 (k = 2, 3), i.e. psi2(p) = sum_n psi2_exact[n] p^n / k0, '
                   'psi3(p) = sum_n psi3_exact[n] p^n / k0 (exact rationals, lowest power first); p = nu/Lam.',
        k0='27.48509817', normalisation='k(0) := -2 int psi2 p^6 dp + 2 int psi3 p^4 dp = 1 (with the 1/k0)',
        r=nf['r'], X=nf['X'], Y=nf['Y'], psi2_exact=nf['psi2_exact'], psi3_exact=nf['psi3_exact'], note=nf['note']), indent=1))
    # 2. Q3-full-delta0.2 (Result Q, flat functional)
    q = json.load(open(os.path.join(C, 'resultq', 'functionals_q3.json')))['Q3-full-delta0.2']
    put('q3_full_delta0.2.json', json.dumps(dict(
        name='Q3-full-delta0.2', source='ads4/code/resultq/functionals_q3.json',
        definition=('flat Stage-0 action on a D = 4 heavy state (u = m_flat^-2-normalised, y = u p^2): '
                    'FF = int_0^1 psi2(p) (2 - y) D_h dp + s u int_0^1 psi3(p) D_h dp + a (2J(J+1) - c4 + s) u^3 - 2 c u^2, '
                    'D_++ = P_J(1 - 2y), D_+- = P^(0,8)_(J-4)(1 - 2y), (s, c4) = (+1, 1) for ++ and (-1, 41) for +-; '
                    'Stage-0 norm N = |B2| + |B3| + |FW| + 2 c u^2; psi2, psi3 polynomials in p (lowest power first, IEEE doubles '
                    'given as repr strings = exact binary values); A = 2 psi2\'(0); upper end 2L = A log(2Lam/3) + B_eff'),
        **{k: (repr(v) if isinstance(v, float) else ([repr(x) for x in v] if isinstance(v, list) else v)) for k, v in q.items()}),
        indent=1))
    # 3. kernels
    kq4 = json.load(open(os.path.join(C, 'resultq', 'functionals_q4_deshift.json')))['EE3']['R_j']
    kq5 = json.load(open(os.path.join(C, 'resultq', 'functionals_q5_deshift.json')))
    Rm, Rp = kern.load_D()
    put('kernels.json', json.dumps(dict(
        note='leading-Regge heavy images (exact rationals in nu, sympy syntax, I = sqrt(-1)); the de-shift identity is W[w](t) = Phi(t)',
        Q_k4_EE3=dict(R_j=kq4, image='sum_{j=+-1,+-3} R_j(nu) P_{-1/2+i(nu-ij)}(cosh beta)',
                      W='W[chi](t) = int_0^oo chi(nu) sum_j R_j(nu) cos((nu - i j) t) dnu'),
        Q_k5_C5_pp=dict(R_j=kq5['R_j_k5_pp'], image='Re[ c(nu)/Nhat5(nu) sum_{j=0,+-2} R_j(nu) P_{-1/2+i(nu-ij)}(cosh beta) ], '
                        'Nhat5 = sum_j R_j; +- class = -(++)', W='W[c](t) = Re int_0^oo (c/Nhat5) sum_j R_j cos((nu - i j) t) dnu',
                        raw_image_nu0=kq5['raw_image_nu0'], R_0_k3_pp=kq5['R_0_k3_pp']),
        D_k2_EE=dict(R_minus=str(Rm.expr), R_plus=str(Rp.expr),
                     abel='C[a Y] - i tanh t S[b Y] = C[psi2]/cosh t, X = a Y, a = (R_- + R_+)/n, b = (R_- - R_+)/n, '
                          'n = (|R_-| + |R_+|) nq/|nq|, nq = R_-^{q2} + R_+^{q2} (below)',
                     q2_minus=str(kern.Rat(eval(pickle.load(open(os.path.join(C, 'spinfunc', 'cache', 'u2_hclosed.pkl'), 'rb'))['q2'][0],
                                                     {**__import__('sympy').__dict__})).expr),
                     q2_plus=str(kern.Rat(eval(pickle.load(open(os.path.join(C, 'spinfunc', 'cache', 'u2_hclosed.pkl'), 'rb'))['q2'][1],
                                                    {**__import__('sympy').__dict__})).expr))), indent=1))
    # 4. stored XP grid weights (Omega^Q_up as fixed on 22 Nov)
    X = json.load(open(os.path.join(C, 'resultq', 'functionals_q7_xprec.json')))
    for key in ('C4_value', 'C4_deriv'):
        e = X[key]
        put(f'q_deshift_xp_{key}.txt', f"# Result Q Stage-7 XP de-shift weight {key} (functionals_q7_xprec.json); grid nu_k = (k + 1/2)/25, "
            f"k = 0..{e['N'] - 1}; functional = (1/25) sum_k chi_k Psi_(4,nu_k) (midpoint rule)\n# k  nu_k  chi_k (64 digits)\n" +
            "".join(f"{k} {(k + 0.5)/25:.2f} {v}\n" for k, v in enumerate(e['chi'])))
    e = X['C5']
    put('q_deshift_xp_C5.txt', "# Result Q Stage-7 XP de-shift weight C5 (functionals_q7_xprec.json); grid nu_k = (k + 1/2)/25; "
        "functional = (1/25) sum_k Re[c_k Psi_(5,nu_k)/Nhat5(nu_k)]\n# k  nu_k  Re c_k  Im c_k\n" +
        "".join(f"{k} {(k + 0.5)/25:.2f} {a} {b}\n" for k, (a, b) in enumerate(zip(e['c_re'], e['c_im']))))
    # 5. continuum (cut-off) weights
    for piece, tag in zip(('c4v', 'c4d', 'c5'), tags):
        P, x, d = AC.load(piece, tag)
        chk = open(os.path.join(kern.CACHE, f'ab_check_{piece}_{tag}.txt')).read()
        meta = dict(piece=piece, source=f'ads4/code/deshift/cache/ab_{piece}_{tag}.json (ab_solve.py), certified by ab_check.py',
                    formula=('C4 (chi, even, real on R): chi(nu) = Q(nu) B(nu), Q = nu^2 (nu^2+1)(nu^2+4) / [(nu^2+9)^p (4nu^2+1)(4nu^2+25)'
                             '(4nu^2+(3+2sqrt2)^2)], p = 1 (value) or 3 (derivative); C5: D(nu) = c(nu)/Nhat5(nu) = Q(nu) B(nu), '
                             'Q = i nu (nu^2+1)/(nu^2+4) * prod_a [i/(nu - i a)], a in {-1.195925515207959757416778302592622585275..., '
                             '1.5, sqrt(73) - 7.5}; c = D Nhat5 on the real line.  B(nu) = x_0 g(nu) + sum_{k>=1} [x_k^re (g(nu-kh) + g(nu+kh))'
                             ' + x_k^im i (g(nu-kh) - g(nu+kh))] with g(z) = exp(-z^2/(2 sigma^2)) (C4: real x_k^re only).'),
                    h=d['h'], sigma=d['sigma'], K=d['K'], labels=d['labels'], coefficients=d['coeffs'],
                    certificate=[l for l in chk.splitlines() if l.startswith(('(1)', '(2)', '(3)', '(5)', '(6)'))])
        put(f'q_deshift_continuum_{piece}.json', json.dumps(meta, indent=0))
        S = AC.Sol(P, x)
        rows = []
        for k in range(0, 3001):
            nu = arb(k)/50
            if P.sym == 'even':
                w = S.weight(nu); rows.append(f"{k/50:.2f} {w.real.str(30, radius=False)}\n")
            else:
                z = acb(nu); Dv = S.u(z); c = Dv*P.Nhat(z) if k > 0 else acb(0)     # D(0) = 0 exactly (Q has a factor nu)
                rows.append(f"{k/50:.2f} {c.real.str(25, radius=False)} {c.imag.str(25, radius=False)} "
                            f"{Dv.real.str(25, radius=False)} {Dv.imag.str(25, radius=False)}\n")
        hdr = ("# continuum de-shift weight " + piece + " tabulated at nu = k/50 (0 <= nu <= 60; |w| < 1e-25 max|w| beyond); exact "
               "definition in q_deshift_continuum_" + piece + ".json\n# nu  " + ("chi(nu)" if P.sym == 'even' else "Re c  Im c  Re D  Im D") + "\n")
        put(f'q_deshift_continuum_{piece}_table.txt', hdr + "".join(rows))
    # 6. Result D Abel weights (stored u3, min-norm gauge) at Dg = 300, 1000
    for Dg in (300, 1000):
        st = pickle.load(open(os.path.join(C, 'spinfunc', 'cache', f'u3_chi_EE_{Dg}.pkl'), 'rb'))
        dnu = st['nu'][1] - st['nu'][0]
        put(f'd_abel_X_{Dg}.txt', f"# Result D k = 2 Abel weight X = chi (h_- + h_+) (u3_chi_EE_{Dg}.pkl; NF02 psi2, EE kernel, min-norm "
            f"gauge = Y - psi2 of minimal l2 norm on this grid); Lam = {st['Lam']:.10f}; grid nu_k = (k + 1/2) {dnu:.2f}; Y = X/a; "
            f"functional = Re sum_k dnu X_k Psi^_(EE,nu_k)/(h_- + h_+)(nu_k)\n# nu_k  Re X  Im X  Re Y  Im Y  psi2(nu_k)\n" +
            "".join(f"{n:.3f} {Xv.real:.17e} {Xv.imag:.17e} {Yv.real:.17e} {Yv.imag:.17e} {ph:.17e}\n"
                    for n, Xv, Yv, ph in zip(st['nu'], st['X'], st['Y'], st['phi'])))
    put('MANIFEST.txt', "".join(f"{n}  {s} bytes  sha256[:16] {h}\n" for n, s, h in manifest))
    print("\n".join(f"{n} {s}" for n, s, h in manifest))


if __name__ == '__main__':
    main(sys.argv[1:4])
