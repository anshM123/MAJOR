"""lj_critmulti (HLJ_ANALYTIC.md §3.4): the multi-state critical-point search (80 random starts per state/nu/wave) used in §3.4."""
import os, sys; os.environ['L1_DEG']='2'; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lj_crit as C, lj_saddle as S, numpy as np, l1_crystal as Cr
rng=np.random.default_rng(7)
for (Delta,J,nu) in [(1100.,100,10.),(1100.,100,900.),(1300.,300,10.),(1300.,300,300.),(364.,64,10.),(364.,64,200.)]:
  alpha=Delta-J-2.; m=np.sqrt(alpha*(alpha+2*J+1)); beta=np.arccosh(1+(2*J+1)**2/(2*m*m))
  for sgn in (1,-1):
    zc,a=Cr.saddle(nu,Delta)
    zp,okp,_=S.track_J(nu,Delta,J,sgn); Ip=S.info(zp,nu,Delta,J,sgn)
    nok=0; pts=[]
    for sig in [0.05,0.2,0.5,1.0,2.0]:
      for n in range(16):
        z=zc.copy(); z[1:]=zc[1:]+a*sig*(rng.normal(size=6)+1j*rng.normal(size=6)); z[5:]=a*sig*(rng.normal(size=2)+1j*rng.normal(size=2)); z[0]=rng.uniform(-.9,.9)*(1+.3j*rng.normal())
        z,ok=C.newton_free(z,nu,Delta,J,sgn)
        if ok:
          I=S.info(z,nu,Delta,J,sgn)
          if I['gnorm']<1e-6 and np.isfinite(I['Phi']): nok+=1; pts.append((z,I))
    seen=[]; lines=[]
    for z,I in pts:
      if any(np.abs(z-s).max()<1e-6*max(1,np.abs(z).max()) for s in seen): continue
      seen.append(z)
      lines.append('     dRePhi=%+9.3f kt/beta=%+8.4f t=%s x0=%s'%((I['Phi']-Ip['Phi']).real,((I['kappa']-np.log(2))/beta).real,np.round(I['t'],6),np.round(z[0],4)))
    print('Delta=%g J=%d nu=%g wave %+d: main saddle ok=%s kt/beta=%.4f; %d/80 starts converged, %d distinct'%(Delta,J,nu,sgn,okp,((Ip['kappa']-np.log(2))/beta).real,nok,len(seen)))
    print('\n'.join(lines))
