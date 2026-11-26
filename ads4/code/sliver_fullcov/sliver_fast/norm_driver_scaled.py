from fractions import Fraction as F
from math import isqrt,ceil
from pathlib import Path
import json,sys,subprocess,concurrent.futures,time,os

def arctan_bracket(q,n=64):
 s=sum((F((-1)**k,(2*k+1)*q**(2*k+1)) for k in range(n)),F(0));e=F(1,(2*n+1)*q**(2*n+1));return s,s+e
lo,hi=arctan_bracket(5);l2,h2=arctan_bracket(239);plo=16*lo-4*h2;phi=16*hi-4*l2
PI0=F('3.1415926535897932384626433832795028841971693993751')
assert max(abs(PI0-plo),abs(PI0-phi))<F(1,10**48) and phi<F(22,7)
NODE_DEN=10**24;DEN=2**61

def cos_nodes(n):
 v=[]
 for j in range(n+1):
  x=PI0*j/n;z=-x*x;term=F(1);s=term
  for k in range(1,27):term*=z/((2*k-1)*(2*k));s+=term
  rem=abs(term*z/(53*54));assert rem<F(1,10**40)
  y=(1+s)/2;u=min(NODE_DEN,max(0,round(y*NODE_DEN)));assert abs(F(u,NODE_DEN)-y)<=F(1,NODE_DEN)
  v.append(u)
 assert v[0]==NODE_DEN and v[-1]==0 and all(v[i]>v[i+1] for i in range(n));return v

def sqrtup(x):return F(isqrt(x.numerator*10**50//x.denominator)+1,10**25)

def ngrid(m,err):
 eta=F(22,7)*sqrtup(err/2);n=ceil(F(7,5)*m)
 while m*(F(22,7)/(2*n)+eta)>=F(113,100):n+=1
 a=F(113,100);assert 1-a*a/2+a**4/24-a**6/720>F(2,5)
 return n

def main(path):
 p=Path(path).resolve();data=json.load(open(p));l,u,Hl,Hu=map(F,(data['l'],data['u'],data['Hl'],data['Hu']));K=data['K'];mc=data['c']['deg'];my=data['y']['deg'];assert l<u and 0<Hu<=Hl
 for c,h in [(l,Hl),(u,Hu)]:
  z=F(isqrt(c.numerator*10**44//c.denominator),10**22);assert z*z<=c and h>=(1-z)**2/998**2
 for s in (1,9):
  for j in range(K):
   pp=j*(j+s)
   for c,h in [(l,F(0)),(u,F(0)),(l,Hl),(u,Hu)]:assert 1-pp*h>0 and abs(c-pp*h)<=1-pp*h
 # exact prefix schedule: choose integer shift_k with 2^shift_k * product(max |r|) <= 1.
 # cumulative product rounded up at fixed 2^100 scale, rescale by factors of 2 while <= 1/2.
 sc=2**100;upper=sc;shifts=[0];sh=0
 for j in range(K):
  mx=max(abs(c-h*j*(j+s))/(1-h*j*(j+s)) for c,h in [(l,F(0)),(u,F(0)),(l,Hl),(u,Hu)] for s in (1,9))
  assert mx<=1
  upper=ceil(upper*mx)
  t=0
  while upper<=sc//2 and t<2:upper*=2;sh+=1;t+=1
  assert sc//2<upper<=sc;shifts.append(sh)
 # rounding errors in coefficient conversion < 1 unit; stage division < 1 unit; T multiply <1; use 10 units.
 err=(K+1)*(F(4)+Hl*((K+1)**2+K*K+8*K+67))/10**60+F(10*(K+1),2**180);assert err<F(1,10**45)
 Ncx=ngrid(mc,F(4,10**24)+F(2,DEN)/(u-l));Ny=ngrid(my,F(4,10**24)+F(2,DEN)/Hu)
 print('grid',Ncx+1,Ny+1,'scalar nodes',2*(Ncx+1)*(Ny+1),'last shift',shifts[-1],flush=True)
 yc=cos_nodes(Ncx);yy=cos_nodes(Ny);rows=[]
 for a in yc:
  c0=l+(u-l)*F(a,NODE_DEN);C=max(ceil(l*DEN),min((u*DEN).numerator//(u*DEN).denominator,(c0*DEN).numerator//(c0*DEN).denominator));c=F(C,DEN);assert l<=c<=u and abs(c-c0)<=F(1,DEN)
  H=Hl+(Hu-Hl)*(c-l)/(u-l);rows.append((C,H.numerator,H.denominator))
 assert len(set(C for C,_,_ in rows))==len(rows)
 wd=p.parent;stem=p.stem;parts=[]
 for worker in range(4):
  r=rows[worker::4];inp=wd/f'{stem}_{worker}.in';out=wd/f'{stem}_{worker}.out';log=wd/f'{stem}_{worker}.log';parts.append((inp,out,log,len(r)))
  inp.write_text(f'{K} {DEN} {len(yy)} {len(r)} 180\n'+' '.join(map(str,shifts))+'\n'+f'{NODE_DEN}\n'+ ' '.join(map(str,yy))+'\n'+''.join(f'{C} {Hn} {Hd}\n' for C,Hn,Hd in r))
 def worker(part):
  inp,out,log,num=part
  with open(inp)as f,open(out,'w')as g,open(log,'w')as h:subprocess.run([str(wd/'scalar_scaled'),os.environ.get('SLIVER_COEFF',str(wd/'coeff_compact60.txt'))],stdin=f,stdout=g,stderr=h,check=True)
  return out
 st=time.monotonic()
 with concurrent.futures.ThreadPoolExecutor(max_workers=4) as exe:
  for res in exe.map(worker,parts):print('done',res.name,flush=True)
 mins={0:10**100,4:10**100};maxs={0:-10**100,4:-10**100};count=0
 for inp,out,log,num in parts:
  for line in out.read_text().splitlines():
   a,*z=line.split()
   if a=='COUNT':assert int(z[0])==num*len(yy)*2;count+=int(z[0]);continue
   b=int(a);v,w=map(int,z);mins[b]=min(mins[b],v);maxs[b]=max(maxs[b],w)
 assert count==2*(Ncx+1)*(Ny+1)
 result=[]
 for b in (0,4):
  a=F(mins[b],10**20)-err;z=F(maxs[b],10**20)+err;R=(a+z)/2;delta=(z-a)/2;E_c=F(data['c']['error']);E_y=F(data['y']['error']);C=F(5,2)
  lower=R-C*C*delta-C*(C+1)*E_y-(C+1)*E_c
  Hphys=lower/(1+b*(b+1)*Hl)**4-F(data['physical_tail'])-F(1,2**140)-F(1,10**14)
  result.append(dict(b=b,min_node=str(a),max_node=str(z),reference=str(R),delta=str(delta),prefix_lower=str(lower),physical_lower=str(Hphys),lower_float=float(Hphys),status='PASS_CONDITIONAL' if Hphys>0 else 'INCONCLUSIVE'))
 print(result,'seconds',time.monotonic()-st,flush=True)
 json.dump(dict(cell=stem,count=count,Nc=Ncx,Ny=Ny,scalar_error=str(err),results=result),open(wd/f'{stem}_ledger.json','w'),indent=2)
if __name__=='__main__':main(sys.argv[1])
