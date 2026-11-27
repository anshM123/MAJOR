import sys, ast
for f in sys.argv[1:]:
    print(f)
    for l in open(f):
        if l.startswith('{'):
            r = ast.literal_eval(l)
            print(r['J0'], r['b'], r['a0'], 'J=', r['J'], 'bh/lam=%.1f ell=%.4f H=%.4g dev=%.3e last=%.1e rad=%.1e' % (r['bhat_over_lam'], r['ell'], r['H'], r['dev'], r['last'], r['Srad']))
