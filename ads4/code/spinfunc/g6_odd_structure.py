"""
g6 (M6, C^_3 ingredient): covariant parity-odd photon (2,4) structure with crystal image O = W2V4 - V2W4 (doc 21.6).
Candidates are the d = 3 (epsilon-dualised) j' = 1 monomials of CLS-D eq. 3pt_convention_spinning: eps(Z2,X2,X4,X0,Z0) H40 (+ 2<->4) and
V042 [eps(Z2,..) V4 - V2 eps(Z4,..)], evaluated at the D = 4 crystal with the g4 polarisation maps.
Usage: python3 g6_odd_structure.py > cache/g6_odd_structure.txt
"""
import os, numpy as np, itertools
HERE='/home/user/MAJOR/ads4/code/spinfunc'
src = open(os.path.join(HERE, 'g4_images.py')).read().split("for x in (0.1")[0]
G = {'__file__': os.path.join(HERE, 'g4_images.py')}; exec(compile(src, 'g4h', 'exec'), G)
H, V, Xij, ldot, cart_to_lc, flat_vectors, EX, EZ = G['H'], G['V'], G['Xij'], G['ldot'], G['cart_to_lc'], G['flat_vectors'], G['EX'], G['EZ']
def vec5(E): return np.array([E[0], E[1], *E[2]], dtype=complex)
def eps(*Es): return np.linalg.det(np.array([vec5(E) for E in Es]))
for x in (0.2, 0.5, 0.8):
    m=1.0; nu, s, (V2,W2,V4,W4), fd, B2, B4 = flat_vectors(x, m)
    a=(m-s)/nu; Xc=np.array([-1j*a,1j*a,0.0]); x4=Xc/(-ldot(Xc,Xc)); x2=-x4
    z0=np.array([1.0,0,0],dtype=complex); Z=[cart_to_lc(0,-1.0,0.0), cart_to_lc(0,0.0,1.0)]
    Oimg=np.array([[fd(B2[i],W2)*fd(B4[j],V4)-fd(B2[i],V2)*fd(B4[j],W4) for j in range(2)] for i in range(2)])
    Qimg=np.array([[fd(B2[i],W2)*fd(B4[j],W4)-fd(B2[i],V2)*fd(B4[j],V4) for j in range(2)] for i in range(2)])
    print(f"x={x}: O image {np.round(Oimg.flatten(),6)}   Q image {np.round(Qimg.flatten(),6)}")
    for sgn in (+1,-1):
        x0=np.array([0,0,sgn*a],dtype=complex)
        X0,X2,X4=EX(x0),EX(x2),EX(x4); Z0=EZ(x0,z0)
        cands={}
        for nm,f in (('e(Z2,Z4,X0,X2,X4)',lambda z2,z4: eps(EZ(x2,z2),EZ(x4,z4),X0,X2,X4)),
                     ('e(Z2,Z4,Z0,X2,X4)',lambda z2,z4: eps(EZ(x2,z2),EZ(x4,z4),Z0,X2,X4)),
                     ('e(Z2,Z4,X0,Z0,X2)',lambda z2,z4: eps(EZ(x2,z2),EZ(x4,z4),X0,Z0,X2)),
                     ('e(Z2,Z4,X0,Z0,X4)',lambda z2,z4: eps(EZ(x2,z2),EZ(x4,z4),X0,Z0,X4)),
                     ('e(Z2,X2,X4,X0,Z0)V4',lambda z2,z4: eps(EZ(x2,z2),X2,X4,X0,Z0)*V(x4,z4,x0,x2)),
                     ('V2 e(Z4,X4,X2,X0,Z0)',lambda z2,z4: V(x2,z2,x0,x4)*eps(EZ(x4,z4),X4,X2,X0,Z0)),
                     ('e(Z2,X2,X4,X0,Z0)H40',lambda z2,z4: eps(EZ(x2,z2),X2,X4,X0,Z0)*H(x4,z4,x0,z0)),
                     ('H20 e(Z4,X4,X2,X0,Z0)',lambda z2,z4: H(x2,z2,x0,z0)*eps(EZ(x4,z4),X4,X2,X0,Z0))):
            img=np.array([[f(Z[i],Z[j]) for j in range(2)] for i in range(2)])
            cands[nm]=img
            r=img.flatten()/np.where(abs(Oimg.flatten())>1e-12,Oimg.flatten(),np.nan)
            print(f"   S0 {sgn:+d} {nm:24s} img {np.round(img.flatten(),6)}  ratio/O {np.round(r,6)}")
print("==== normalised o_H, o_V vs O at S0 -1 (all four entries), and S0 +1 vs reflected O")
for x in (0.1,0.3,0.5,0.7,0.9):
    m=1.0; nu, s, (V2,W2,V4,W4), fd, B2, B4 = flat_vectors(x, m)
    a=(m-s)/nu; Xc=np.array([-1j*a,1j*a,0.0]); x4=Xc/(-ldot(Xc,Xc)); x2=-x4
    z0=np.array([1.0,0,0],dtype=complex); Z=[cart_to_lc(0,-1.0,0.0), cart_to_lc(0,0.0,1.0)]
    Oimg=np.array([[fd(B2[i],W2)*fd(B4[j],V4)-fd(B2[i],V2)*fd(B4[j],W4) for j in range(2)] for i in range(2)])
    for sgn in (-1,+1):
        x0=np.array([0,0,sgn*a],dtype=complex); X0,X2,X4=EX(x0),EX(x2),EX(x4); Z0=EZ(x0,z0)
        v=V(x0,z0,x4,x2)
        oH=lambda z2,z4: (eps(EZ(x2,z2),X2,X4,X0,Z0)*H(x4,z4,x0,z0)/Xij(x4,x0)+H(x2,z2,x0,z0)*eps(EZ(x4,z4),X4,X2,X0,Z0)/Xij(x2,x0))/Xij(x2,x4)
        oV=lambda z2,z4: v*(eps(EZ(x2,z2),X2,X4,X0,Z0)*V(x4,z4,x0,x2)/Xij(x4,x0)-V(x2,z2,x0,x4)*eps(EZ(x4,z4),X4,X2,X0,Z0)/Xij(x2,x0))/Xij(x2,x4)
        Ot = Oimg if sgn==-1 else -Oimg*np.array([[-1,1],[1,-1]])
        for nm,f in (('oH',oH),('oV',oV)):
            img=np.array([[f(Z[i],Z[j]) for j in range(2)] for i in range(2)]); r=(img/Ot).flatten()
            print(f"  x={x} S0 {sgn:+d} {nm}: ratio img/O entries {np.round(r,8)}  (spread {np.ptp(np.abs(r))/np.abs(r).max():.1e}); times 4nu^2(m^2-nu^2)^2: {np.round(r[0]*4*nu**2*(m*m-nu*nu)**2,8)}")
print("==== corrected homogeneity: oH' = [e(Z2,X2,X4,X0,Z0) H40 / sqrt(X24 X20 X40^3) + (2<->4)];  ratio to q-image normalisation")
for x in (0.1,0.3,0.5,0.7,0.9):
    m=1.0; nu, s, (V2,W2,V4,W4), fd, B2, B4 = flat_vectors(x, m)
    a=(m-s)/nu; Xc=np.array([-1j*a,1j*a,0.0]); x4=Xc/(-ldot(Xc,Xc)); x2=-x4
    z0=np.array([1.0,0,0],dtype=complex); Z=[cart_to_lc(0,-1.0,0.0), cart_to_lc(0,0.0,1.0)]
    Oimg=np.array([[fd(B2[i],W2)*fd(B4[j],V4)-fd(B2[i],V2)*fd(B4[j],W4) for j in range(2)] for i in range(2)])
    x0=np.array([0,0,-a],dtype=complex); X0,X2,X4=EX(x0),EX(x2),EX(x4); Z0=EZ(x0,z0)
    n=np.sqrt(Xij(x2,x4)*Xij(x2,x0)*Xij(x4,x0)**3+0j)
    oH=lambda z2,z4: (eps(EZ(x2,z2),X2,X4,X0,Z0)*H(x4,z4,x0,z0)+H(x2,z2,x0,z0)*eps(EZ(x4,z4),X4,X2,X0,Z0))/n
    img=np.array([[oH(Z[i],Z[j]) for j in range(2)] for i in range(2)]); r=(img/Oimg).flatten()
    print(f"  x={x}: ratio {np.round(r[0],8)} spread {np.ptp(np.abs(r))/np.abs(r).max():.1e};  x 4nu^2(m^2-nu^2)^2 = {np.round(r[0]*4*nu**2*(m*m-nu*nu)**2,8)}; X20={Xij(x2,x0):.5f} X24={Xij(x2,x4):.5f} V042={V(x0,z0,x4,x2):.5f}")
