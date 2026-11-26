from fractions import Fraction as F
from pathlib import Path
import subprocess,random
# exact signed spectral atoms at rational nu^2, with moments scaled by 1000
rng=random.Random(1325)
K=32;D=200;S=10**D
atoms=[(F(rng.randrange(-50,1000001),rng.randrange(1,9)),F(rng.randrange(-30,31),17)) for j in range(12)]
T=[sum(w*(v/1000000)**j for v,w in atoms) for j in range(K+1)]
lines=f'{K} {D}\n'+''.join(f'{int(x*S)} 1\n' for x in T)
r=subprocess.run([str(Path(__file__).resolve().parent/'raw_moment_transform')],input=lines,text=True,capture_output=True,check=True)
prev=[F(0)]*len(atoms);cur=[F(1)]*len(atoms)
for line in r.stdout.splitlines():
 m,x,e=map(int,line.split());assert abs(sum(w*f for (v,w),f in zip(atoms,cur))*S-x)<=e
 nxt=[((2*m*m-v-F(1,4))*g-(m-1)**2*f)/(m+1)**2 for (v,w),f,g in zip(atoms,prev,cur)];prev,cur=cur,nxt
assert m==K+1
print(f'PASS all {m} raw-transform outputs against independent exact rational recurrence; {len(atoms)} signed atoms')
