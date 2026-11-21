import json, os
from common_q import CACHE
for i, l in enumerate(open(os.path.join(CACHE, 'lp_g4.jsonl'))):
    r = json.loads(l)
    tag = '%2d d=%d/%d r=%d delta=%.2f %s%s%s' % (i, r['d2'], r['d3'], r['r'], r['delta'], 'B3' if r['b3'] else 'noB3', ' psi3(0)=0' if r['psi3zero'] else '', ' matter' if r['matter'] else '')
    if r.get('fwd') is False: tag += ' NOFWD'
    if r.get('G3'): tag += ' G3=%g lref=%s' % (r['G3'], r['lref'])
    print(tag, ('A=%.3f B=%.3f kappa=%.4f a=%.4f' % (r['A'], r['B'], r['kappa'], r['a'])) if r.get('status') == 0 else 'status %s' % r['status'])
