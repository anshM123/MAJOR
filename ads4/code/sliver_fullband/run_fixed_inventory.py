from pathlib import Path
from fractions import Fraction as F
import concurrent.futures,subprocess,json,sys,time,os
sys.set_int_max_str_digits(0)
K=int(sys.argv[1]);D=int(sys.argv[2]);J=int(sys.argv[3]);rows=list(zip(json.load(open('weights.json'))['centres_sigma'],json.load(open('weights.json'))['x_re']));assert len(rows)==3564
st=time.monotonic()
def run(i):
 r=rows[len(rows)*i//J:len(rows)*(i+1)//J]; rr=[' '.join(f'{q.numerator} {q.denominator}'for q in map(F,[*cs,x])) for cs,x in r]; tag=f'fixed_{K}_{D}_{i}'
 with open(tag+'.seeds','w') as out,open(tag+'.seed.log','w') as err:subprocess.run(['./gaussian_seed'],input=str(len(rr))+'\n'+'\n'.join(rr)+'\n',text=True,stdout=out,stderr=err,check=True)
 seed=Path(tag+'.seeds').read_text().splitlines();assert len(seed)==len(rr);assert all(0<=int(a)<=int(b) and int(b)-int(a)<=10**20 for a,b in map(str.split,seed));v=f'{K} {D} {len(rr)}\n'+'\n'.join(' '.join(x) for x in zip(rr,seed))+'\n'
 print('SEEDS_PASS',i,'seconds',time.monotonic()-st,flush=True)
 with open(tag+'.txt.tmp','w') as out,open(tag+'.fixed.log','w') as err:subprocess.run(['./gaussian_fixed'],input=v,text=True,stdout=out,stderr=err,check=True)
 os.replace(tag+'.txt.tmp',tag+'.txt');print('MOMENTS_PASS',i,'seconds',time.monotonic()-st,flush=True);return tag+'.txt'
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:fs=list(pool.map(run,range(J)))
lo=[0]*(K+1);hi=[0]*(K+1)
for f in fs:
 rr=[tuple(map(int,l.split()))for l in open(f)];assert len(rr)==K+1
 for i,(j,l,u)in enumerate(rr):assert j==i and l<=u;lo[i]+=l;hi[i]+=u
with open(f'S_{K}_{D}.txt','w') as f:
 for j in range(K+1):f.write(f'{j} {lo[j]} {hi[j]}\n')
print('PASS_FULL_INVENTORY_SEED_SUBSTITUTED',len(rows),'K',K,'D',D,'maxWidth',max(u-l for l,u in zip(lo,hi)),'seconds',time.monotonic()-st,flush=True)
