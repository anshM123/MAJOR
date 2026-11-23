"""
ext/pack.py -- collect every finished point into one JSON file (results/ext_results.json) that can be committed to git.
Numpy arrays and complex numbers are encoded losslessly; ads4/ext/unpack.py restores the original pickles.
"""
import os, sys, json, pickle, glob
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))


def enc(o):
    if isinstance(o, np.ndarray):
        if np.iscomplexobj(o):
            return {'__nd__': 'c', 'dtype': str(o.dtype), 'shape': list(o.shape), 're': o.real.ravel().tolist(), 'im': o.imag.ravel().tolist()}
        return {'__nd__': 'r', 'dtype': str(o.dtype), 'shape': list(o.shape), 'v': o.ravel().tolist()}
    if isinstance(o, (complex, np.complexfloating)):
        return {'__c__': [float(o.real), float(o.imag)]}
    if isinstance(o, np.generic):
        return o.item()
    if isinstance(o, dict):
        return {'__d__': [[enc(k), enc(v)] for k, v in o.items()]}
    if isinstance(o, tuple):
        return {'__t__': [enc(v) for v in o]}
    if isinstance(o, list):
        return [enc(v) for v in o]
    if isinstance(o, float) and (o != o or o in (float('inf'), float('-inf'))):
        return {'__f__': repr(o)}
    return o


if __name__ == '__main__':
    res = []
    for d in sorted(glob.glob(os.path.join(HERE, 'out', '*', '*'))):
        if not os.path.exists(os.path.join(d, 'done')):
            continue
        for f in sorted(glob.glob(os.path.join(d, '*.pkl'))):
            res.append({'dir': os.path.relpath(d, HERE), 'file': os.path.basename(f), 'data': enc(pickle.load(open(f, 'rb'))),
                        'log': open(os.path.join(d, 'log.txt')).read()[-2000:]})
    os.makedirs(os.path.join(HERE, 'results'), exist_ok=True)
    out = os.path.join(HERE, 'results', 'ext_results.json')
    json.dump(res, open(out, 'w'))
    print(f'packed {len(res)} files -> {out} ({os.path.getsize(out)/1e6:.1f} MB)')
