"""ext/unpack.py -- inverse of pack.py: results/ext_results.json -> pickles under DEST/<dir>/<file> (DEST default: ext/unpacked)."""
import os, sys, json, pickle
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))


def dec(o):
    if isinstance(o, dict):
        if '__nd__' in o:
            a = (np.array(o['re']) + 1j*np.array(o['im'])) if o['__nd__'] == 'c' else np.array(o['v'])
            return a.astype(o['dtype']).reshape(o['shape'])
        if '__c__' in o:
            return complex(*o['__c__'])
        if '__d__' in o:
            return {dec(k): dec(v) for k, v in o['__d__']}
        if '__t__' in o:
            return tuple(dec(v) for v in o['__t__'])
        if '__f__' in o:
            return float(o['__f__'])
    if isinstance(o, list):
        return [dec(v) for v in o]
    return o


if __name__ == '__main__':
    src = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, 'results', 'ext_results.json')
    dest = sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, 'unpacked')
    for r in json.load(open(src)):
        d = os.path.join(dest, r['dir']); os.makedirs(d, exist_ok=True)
        pickle.dump(dec(r['data']), open(os.path.join(d, r['file']), 'wb'))
    print('unpacked to', dest)
