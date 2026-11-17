"""Referee: recompute near-threshold states (Dg=100, M02) from scratch (no skeleton cache) and compare with cached u4 records."""
import pickle, numpy as np, sys, os
OM = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'omegau'); sys.path.insert(0, OM); os.chdir(OM)
from common import load, Transforms, state, flat_unstripped_parts, flat_parts, gl01, KS
from adsaction import action
import remk
from u4_scan import crit
pn,pw=gl01(64); pf,wf=gl01(400)
Dg=int(sys.argv[1]) if len(sys.argv)>1 else 100; Lam=np.sqrt((Dg-2.)*(Dg-1.))
name=sys.argv[2] if len(sys.argv)>2 else 'M02'
F,PH=load(name); TR=Transforms(PH); W0=Transforms({k:np.array([1.0]) for k in KS})
C=pickle.load(open(f'cache/skel_u4_{Dg}_{remk.TAG}.pkl','rb'))
R=pickle.load(open(f'cache/recs_u4_{Dg}_{remk.TAG}.pkl','rb'))
states=[('++',(1.35,0)),('++',(1.35,2)),('++',(1.4,0)),('++',(1.4,2)),('+-',(1.3,4)),('+-',(1.35,4)),('++',(1.01,0)),('++',(1.01,2))]
for h,key in states:
    mfl,J=key; u=1/mfl**2
    st=state(Dg,Lam,mfl,J)
    if st is None: print(key,'unphysical'); continue
    _,W,P0=W0.regge(Lam,st['beta'])
    sk=remk.skeleton_state(st,J,pn*np.sqrt(u),W,P0)
    d=max(np.abs(sk[k][0]-C[key][k][0]).max() for k in KS) if key in C else np.nan
    Jsu,parts,FWu=flat_unstripped_parts(F,PH,h,u,J+2,pf,wf); Jss,B2,B3,FWs=flat_parts(F,h,u,J+2)
    i=int(np.searchsorted(Jsu,J)); i2=int(np.searchsorted(Jss,J))
    r=action(F,PH,TR,h,st,J,mfl,Lam,pn,pw,sk,parts[i],FWu[i],(B2[i2],B3[i2],FWs[i2])); r['m']=st['m']
    old=[x for x in R[(name.replace('F_g3','Fg3'),h)] if abs(x['mfl']-mfl)<1e-6 and x['J']==J]
    o=crit(old[0],0,'abs')[1] if old else np.nan
    # 'trunc': use the truncated series value everywhere (signed, incl. uncontrolled nodes) -- indicative only
    print(f"{h} m/L={mfl} J={J}: rho diff vs cache {d:.1e}; kappa(abs) new {crit(r,0,'abs')[1]:.4f} cached {o:.4f}; Fh/Ns {r['Fh']/r['Ns']:.4f} "
          f"U/Ns {r['U']/r['Ns']:.4f} NU/Ns {r['NU']/r['Ns']:.3f} E/Ns {r['E']/r['Ns']:.4f}; (Fh+Sall)/Ns {(r['Fh']+r['Sall'])/r['Ns']:+.4f}")
