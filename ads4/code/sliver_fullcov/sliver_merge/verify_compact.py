from pathlib import Path
import sys,hashlib
sys.set_int_max_str_digits(0)
assert len(sys.argv)==2, 'supply the original 3000-digit coeff_long_17200.txt'
with open(sys.argv[1]) as f, open(Path(__file__).with_name('coeff_compact60.txt')) as g:
 n,d=map(int,next(f).split());m,e=map(int,next(g).split());assert (n,d,m,e)==(17200,3000,17200,60)
 factor=10**(d-e)
 for k in range(n+1):
  a,r,b,s=map(int,next(f).split());A,R,B,S=map(int,next(g).split());assert min(r,s,R,S)>=0 and max(R,S)<=2
  assert abs(a-A*factor)+r<=R*factor and abs(b-B*factor)+s<=S*factor,k
 assert not f.read().strip() and not g.read().strip()
print('PASS: all 17201 pairs enclosed')
