import json,subprocess,time,sys,concurrent.futures,os
from fractions import Fraction as F
from pathlib import Path
sys.set_int_max_str_digits(0)
K=int(sys.argv[1]);D=int(sys.argv[2]);J=int(sys.argv[3]) if len(sys.argv)>3 else 4
v=json.load(open('weights.json')); rows=list(zip(v['centres_sigma'],v['x_re']));assert len(rows)==3564
# shard contiguous centers (incremental erf requires sorted centers)
def run(i):
 r=rows[len(rows)*i//J:len(rows)*(i+1)//J];fn=f'g{K}_{D}_{i}.txt';log=fn+'.log'; t=time.monotonic()
 v=f'{K} {D} {len(r)}\n'+''.join(' '.join(f'{q.numerator} {q.denominator}'for q in map(F,[*cs,x]))+'\n'for cs,x in r)
 with open(fn+'.tmp','w') as out,open(log,'w')as err:subprocess.run(['./gaussian_fast_adjacent'],input=v,text=True,stdout=out,stderr=err,check=True)
 os.replace(fn+'.tmp',fn);print('shard',i,'centers',len(r),'elapsed',time.monotonic()-t,flush=True);return fn
st=time.monotonic()
with concurrent.futures.ThreadPoolExecutor(max_workers=J) as pool:fs=list(pool.map(run,range(J)))
lo=[0]*(K+1);hi=[0]*(K+1)
for f in fs:
 rr=[tuple(map(int,l.split()))for l in open(f)];assert len(rr)==K+1
 for i,(j,l,u)in enumerate(rr):assert j==i and l<=u;lo[i]+=l;hi[i]+=u
with open(f'S_{K}_{D}.txt','w') as f:
 for j in range(K+1):f.write(f'{j} {lo[j]} {hi[j]}\n')
print('PASS_FULL_INVENTORY',len(rows),'K',K,'D',D,'maxWidth',max(u-l for l,u in zip(lo,hi)),'seconds',time.monotonic()-st,flush=True)
