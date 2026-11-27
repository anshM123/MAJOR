import json, math, glob
rows = []
for fn in sorted(glob.glob('out/sweep1_*.json')) + sorted(glob.glob('out/sweep2_*.json')):
    d = json.load(open(fn))
    if 'valid' in d[0]:
        ok = [r for r in d if r['valid']]
    else:
        ok = [r for r in d if math.isfinite(r['S']) and (not r['trunc'] or abs(r['last']) < 1e-25*abs(r['S'])) and r['Srad'] < 1e-10*abs(r['S'])]
    neg = [r for r in ok if r['S'] <= 0]
    m = min(ok, key=lambda r: r['S'])
    dev = [abs(r['ratio']-1) for r in ok if r.get('ratio') is not None and r['r'] < 0.45]
    rows.append('| %s | %d / %d | %d | %.3g (J=%d,b=%d,a=%.0f) | %.3g | %.3g | %d | %s |' % (fn[4:], len(ok), len(d), len(neg), m['S'], m['J'], m['b'], m['a'],
        max(r['ell'] for r in ok), max(r['bhat'] for r in ok), max(r['J'] for r in ok), ('%.2e' % max(dev)) if dev else '-'))
hdr = '| file | valid/total | nonpositive | min S=H/w^3 (state) | max ell | max bhat | max J | max abs(H/H_LR-1) (r<0.45) |\n|---|---|---|---|---|---|---|---|'
t = hdr + '\n' + '\n'.join(rows)
open('out/sweep_summary.md', 'w').write(t + '\n'); print(t)
