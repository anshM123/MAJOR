#!/usr/bin/env python3
"""Replay every archived signed scalar calculation and compare exact output.
Default is portable GMP; --simd requires AVX512IFMA and the appropriate CPU.
This does not regenerate the upstream 3000-digit spectral coefficient vector.
"""
from pathlib import Path
import argparse,subprocess,concurrent.futures,os
root=Path(__file__).resolve().parent
p=argparse.ArgumentParser();p.add_argument('--jobs',type=int,default=4);p.add_argument('--simd',action='store_true');p.add_argument('--new-only',action='store_true');p.add_argument('--match',default='');a=p.parse_args();assert a.jobs>=1
build=root/'build';build.mkdir(exist_ok=True)
def compile_one(name,src,simd=False):
    exe=build/name
    subprocess.run(['g++','-O3','-std=c++17']+(['-march=native'] if simd else [])+[str(src),'-lgmpxx','-lgmp','-o',str(exe)],check=True)
    return exe
m=root/'sliver_merge';table=m/'coeff_compact60.txt'
new=compile_one('scalar_ref',m/('scalar_ref80.cpp' if a.simd else 'scalar_ref_portable.cpp'),a.simd)
jobs=[]
for src in sorted(m.glob('cell_*.in')):
    if "_cross" in src.name:continue
    out=src.with_suffix('.out');assert out.is_file() and 'COUNT ' in out.read_text();jobs.append((src,out,[str(new),str(table)]))
if not a.new_only:
    for d,s in [('sliver_close','norm_scalar'),('sliver_fast','scalar_scaled')]:
        folder=root/d;exe=compile_one(s,folder/(s+'.cpp'))
        for src in sorted(folder.glob('cell_*.in')):jobs.append((src,src.with_suffix('.out'),[str(exe),str(table)]))
    for d in ('sliver_next_work','physical_sliver_run'):
        folder=root/d;exe=compile_one(d+'_triangle',folder/'direct_sign_triangle.cpp')
        for src in sorted(folder.glob('next_*.in')):jobs.append((src,src.with_suffix('.out'),[str(exe)]))
jobs=[j for j in jobs if a.match in str(j[0])];assert jobs

def run(job):
    src,out,cmd=job
    with src.open('rb') as f:
        got=subprocess.run(cmd,stdin=f,capture_output=True,check=True).stdout
    if got!=out.read_bytes():raise RuntimeError('exact output mismatch: '+str(src))
    print('PASS',src.relative_to(root),flush=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=a.jobs) as ex:list(ex.map(run,jobs))
print('all',len(jobs),'archived scalar shards replayed exactly',flush=True)
