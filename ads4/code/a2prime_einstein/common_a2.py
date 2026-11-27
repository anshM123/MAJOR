"""Shared read-only imports for the A2' Einstein check.  Nothing outside this directory is written."""
import os, sys
sys.dont_write_bytecode = True   # do not write .pyc next to the read-only imports
HERE = os.path.dirname(os.path.abspath(__file__))
ADS4 = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, os.path.join(ADS4, 'kernels'))
sys.path.insert(0, os.path.join(ADS4, 'angproofs'))
OUT = os.path.join(HERE, 'out')
os.makedirs(OUT, exist_ok=True)

class Tee:
    """print to stdout and to out/<name>.txt"""
    def __init__(self, name):
        self.f = open(os.path.join(OUT, name), 'w')
    def __call__(self, *a):
        s = ' '.join(str(x) for x in a)
        print(s); self.f.write(s + '\n'); self.f.flush()
