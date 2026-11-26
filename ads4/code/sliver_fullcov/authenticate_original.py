"""Authenticate the compact table and the saved direct-certificate prefix
against the original corrected coefficient file from the full-band reduction archive.
"""
import sys,hashlib,json,subprocess
from pathlib import Path
root=Path(__file__).resolve().parent
if len(sys.argv)!=2:raise SystemExit('usage: python authenticate_original.py /path/to/coeff_long_17200.txt')
p=Path(sys.argv[1]).resolve();h=hashlib.sha256()
with p.open('rb') as f:
    for x in iter(lambda:f.read(1<<20),b''):h.update(x)
assert h.hexdigest()=='0d09c4ad62d098ad19e8e75fda06e42776e9a5b571e51653768c0b567624754a'
q=root/'physical_sliver_run/coeff_long_17200_prefix.txt'
with p.open('rb') as f,q.open('rb') as g:
    N,D=map(int,f.readline().split());K,E=map(int,g.readline().split());assert N>=K and D==E==3000
    for k in range(K+1):assert f.readline()==g.readline(),k
    assert not g.read()
subprocess.run([sys.executable,str(root/'sliver_merge/verify_compact.py'),str(p)],check=True)
print('PASS original SHA, exact retained prefix, and all 17201 compact coefficient pairs')
