import json, sys, math
for fn in sys.argv[1:]:
    d = json.load(open(fn))
    ok = [r for r in d if math.isfinite(r['S']) and (not r['trunc'] or abs(r['last']) < 1e-25*abs(r['S'])) and r['Srad'] < 1e-10*abs(r['S'])]
    bad = [r for r in d if r not in ok]
    neg = [r for r in ok if r['S'] <= 0]
    print(fn, 'points', len(d), 'valid', len(ok), 'invalid(trunc/prec)', len(bad), 'NONPOS valid', len(neg))
    if ok:
        m = min(ok, key=lambda r: r['S']); print('  min S', m['S'], {k: m[k] for k in ('J','b','a','bhat','ell')})
        okr = [r for r in ok if r['ratio'] is not None and r['r'] < 0.45]
        if okr:
            m2 = min(okr, key=lambda r: r['ratio']); print('  min H/HLR (r<0.5)', m2['ratio'], {k: m2[k] for k in ('J','b','a','bhat','ell')})
            m3 = max(okr, key=lambda r: abs(r['ratio']-1)); print('  max |H/HLR-1| (r<0.5)', m3['ratio'], {k: m3[k] for k in ('J','b','a','bhat','ell')})
        mx = max(ok, key=lambda r: r['ell']); print('  max ell valid', mx['ell'], mx['J'], mx['a'], 'bhat', mx['bhat'])
        mb = max(ok, key=lambda r: r['bhat']); print('  max bhat valid', mb['bhat'], mb['J'], mb['a'])
