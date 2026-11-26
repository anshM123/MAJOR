#!/usr/bin/env node
'use strict';
const fs=require('fs'),cases=JSON.parse(fs.readFileSync(process.argv[2]||__dirname+'/boxes.json','utf8'));
const S=10n**80n,J=1500,abs=x=>x<0n?-x:x;
const ceil=(x,y)=>{if(x<0n||y<=0n)throw Error('ceil domain');return(x+y-1n)/y},mu=(x,y)=>ceil(x*y,S),du=(x,y)=>ceil(x*S,y);
function sqr(n){if(n<0n)throw Error('negative sqrt');if(n<2n)return n;let x=1n<<BigInt(Math.ceil(n.toString(2).length/2));for(;;){let y=(x+n/x)>>1n;if(y>=x)return x;x=y}}
const hyp=(x,y)=>sqr(x*x+y*y)+1n,hlo=(x,y)=>sqr(x*x+y*y);
function dec(s){let [a,b='']=String(s).split('.');let den=10n**BigInt(b.length),num=BigInt(a+b)*S;if(num%den)throw Error('decimal too long');return num/den}
function fmt(v){let sign=v<0n?'-':'';let a=abs(v).toString().padStart(81,'0');return sign+a.slice(0,-80)+'.'+a.slice(-80)}
function setup(cfg){
 const m=cfg.m,K=cfg.K,R=BigInt(cfg.R),B=m+8;
 if(K<1||K>J||m<1||R<2n)throw Error('bad setup');
 let d0=dec(cfg.d0),dd=dec(cfg.dd),h0=dec(cfg.h0),dh=dec(cfg.dh),n0=dec(cfg.nu0),dn=dec(cfg.dnu),ni=dec(cfg.nuImag||"0"),axis=cfg.nuAxis||"real";
 if(!["real","imag"].includes(axis))throw Error("bad frequency axis");
 if(d0-dd<=0n||d0+dd>=S||h0-dh<0n||n0<0n||(axis==="real"&&n0-dn<0n))throw Error('not real physical box');
 let root=sqr((S-d0-dd)*S),hRequired=ceil((d0+dd)**2n*S,998n**2n*(S+root)**2n);
 const fullTwist=(h0-dh===0n&&h0+dh>=hRequired);if(cfg.fullTwist&&!fullTwist)throw Error('h cap insufficient');
 let ex=[],deg=[],ix=new Map();for(let k=0;k<=m;k++)for(let i=0;i<=k;i++)for(let j=0;j<=k-i;j++){let e=i*B*B+j*B+k-i-j;ix.set(e,ex.length);ex.push(e);deg.push(k)}const N=ex.length,zero=()=>Array(N).fill(0n),shifts=new Map();
 function sh(e,k){if(!shifts.has(e)){let a=new Int32Array(N).fill(-1);for(let i=0;i<N;i++)if(deg[i]+k<=m)a[i]=ix.get(ex[i]+e);shifts.set(e,a)}return shifts.get(e)}
 function sparse(a){let t=new Map;for(let [e,k,r,i] of a){let v=t.get(e)||[e,k,0n,0n];v[2]+=r;v[3]+=i;t.set(e,v)}return [...t.values()].filter(v=>v[2]||v[3]).sort((a,b)=>a[1]-b[1])}
 const con=(r,i=0n)=>[[0,0,r,i]],plus=(...a)=>sparse(a.flat()),times=(a,k)=>a.map(([e,d,r,i])=>[e,d,r*k,i*k]),neg=a=>times(a,-1n);
 function half(a){return a.map(([e,d,r,i])=>{if(r%2n||i%2n)throw Error('non-exact half');return[e,d,r/2n,i/2n]})}
 function prod(a,b){let z=[];for(let [e,d,r,i] of a)for(let [f,g,u,v] of b){let rr=r*u-i*v,ii=r*v+i*u;if(rr%S||ii%S)throw Error('non-exact sparse product');z.push([e+f,d+g,rr/S,ii/S])}return sparse(z)}
 const d=[[0,0,d0,0n],[B*B,1,dd,0n]],h=[[0,0,h0,0n],[B,1,dh,0n]],inu=[[0,0,-ni,n0],[1,1,axis==="imag"?-dn:0n,axis==="real"?dn:0n]],one=con(S);
 function nr(j){let x=plus(con((2n*j+1n)*S/2n),neg(inu));return prod(x,x)}
 function dr(j){return times(plus(con((j+1n)*S),times(inu,-2n)),j+1n)}
 function norm(a){return a.reduce((s,v)=>s+abs(v[2])+abs(v[3]),0n)}
 function invDenBound(a){let v=a.find(x=>x[0]===0),r=v[2],i=v[3];let inv=ceil((abs(r)+abs(i))*S*S,r*r+i*i);let rest=a.filter(x=>x[0]).reduce((s,v)=>s+abs(v[2])+abs(v[3]),0n),gap=S-mu(inv,rest);if(gap<=0n)throw Error('l1 inverse gap failed');return du(inv,gap)}
 function muladd(out,p,a){for(let [e,k,r,i] of a){let shf=sh(e,k);for(let f=0;f<N;f++){let t=shf[f];if(t<0)continue;let u=p.r[f],v=p.i[f];if(!u&&!v)continue;out.r[t]+=r*u-i*v;out.i[t]+=r*v+i*u}}}
 function dense(r=0n){let p={r:zero(),i:zero(),err:0n};p.r[0]=r;return p}
 function solve(A,p,Bb,q,D){let rhs=dense();muladd(rhs,p,A);if(Bb.length)muladd(rhs,q,Bb);let u=dense(),c=D.find(v=>v[0]===0),rr=c[2],ii=c[3],dv=rr*rr+ii*ii;if(!dv)throw Error('zero denominator');let maps=D.filter(x=>x[0]).map(x=>[x,sh(x[0],x[1])]);
 for(let j=0;j<N;j++){let x=rhs.r[j],y=rhs.i[j];let r=(x*rr+y*ii)/dv,i=(y*rr-x*ii)/dv;u.r[j]=r;u.i[j]=i;if(!r&&!i)continue;for(let [[e,k,a,b],map] of maps){let t=map[j];if(t<0)continue;rhs.r[t]-=a*r-b*i;rhs.i[t]-=a*i+b*r}}
 let res=ceil(norm(A)*p.err+norm(Bb)*q.err+2n*BigInt(N)*(abs(rr)+abs(ii)),S);u.err=mu(invDenBound(D),res);return u}
 function addprod(out,p,A){let v=dense();muladd(v,p,A);for(let j=0;j<N;j++){out.r[j]+=v.r[j]/S;out.i[j]+=v.i[j]/S}out.err+=ceil(norm(A)*p.err,S)+2n*BigInt(N)}
 function signedNorm(p){let x=p.err;for(let j=0;j<N;j++)x+=hyp(p.r[j],p.i[j]);return x}
 // Subdivide the already computed polynomial in h. All arithmetic is exact until
 // an upward-rounded Euclidean norm. This does not repeat the recurrence.
 function splitHNorm(p,parts){
   if(!Number.isSafeInteger(parts)||parts<1)throw Error('bad h split count');
   let groups=new Map();for(let l=0;l<N;l++){let i=Math.floor(ex[l]/(B*B)),j=Math.floor(ex[l]/B)%B,k=ex[l]%B,g=i*B+k;if(!groups.has(g))groups.set(g,[]);groups.get(g)[j]=l;}
   let bin=Array.from({length:m+1},()=>Array(m+1).fill(0n));for(let i=0;i<=m;i++){bin[i][0]=bin[i][i]=1n;for(let j=1;j<i;j++)bin[i][j]=bin[i-1][j-1]+bin[i-1][j];}
   let np=BigInt(parts),powers=Array.from({length:m+1},(_,i)=>np**BigInt(i));let mx=0n;
   for(let cell=0;cell<parts;cell++){
     let c=BigInt(2*cell+1-parts),cp=Array.from({length:m+1},(_,i)=>c**BigInt(i));let sum=0n;
     for(let group of groups.values()){
       let degree=group.length-1,den=powers[degree];
       for(let k=0;k<=degree;k++){
         let r=0n,i=0n;for(let j=k;j<=degree;j++){let idx=group[j],w=bin[j][k]*cp[j-k]*powers[degree-j];r+=p.r[idx]*w;i+=p.i[idx]*w;}
         sum+=ceil(hyp(r,i),den);
       }
     }
     if(sum>mx)mx=sum;
   }
   return mx+p.err;
 }
 function major(p,b,R){const s=BigInt(2*b+1),q=BigInt(b*(b+1)),dm=d0+R*dd,hm=h0+R*dh,vr=R*dn;let prev=S,cur=S,total=0n,tail=0n,minDen=S;
 function r(j){if(j<0n)return S;let u=hyp((2n*j+1n)*S/2n+ni,n0)+vr,lo=hlo((j+1n)*S+2n*ni,2n*n0)-2n*vr;if(lo<=0n)throw Error('frequency pole in majorant');return ceil(u*u,(j+1n)*lo)}
 for(let j=0n;j<BigInt(J);j++){let re=(BigInt(2*p-1)*S/2n-j*S-ni),bb=hyp(re,n0)+vr,bc=hyp(re+S,n0)+vr,bs=hyp(re-s*S,n0)+vr;let v=cur;
 if(p===3){let a0=8n*S+ceil(mu(bb,hyp(re+BigInt(8-4*b)*S,n0)+vr),2n),a1=ceil(mu(bb,hyp(re+BigInt(9-4*b)*S,n0)+vr),2n);v=mu(S+q*hm+mu(hm,a0),cur)+mu(mu(hm,a1),mu(r(j-1n),prev))}
 total+=v;if(j>=BigInt(K))tail+=v;if(j===BigInt(J-1))break;
 let den=S-mu(hm,mu(bb,bs));if(den<=0n)throw Error('polydisk denominator not positive');if(den<minDen)minDen=den;
 let aa=du(dm+mu(hm,mu(bb,hyp((s-1n)*S-2n*re,2n*n0)+2n*vr)),den),ab=du(mu(hm,mu(bb,bc)),den);let next=mu(r(j),mu(aa,cur)+mu(mu(ab,r(j-1n)),prev));prev=cur;cur=next;}
 return {total,tail,minDen}}
 function run(p,b){let start=process.hrtime.bigint();const s=BigInt(2*b+1),q=BigInt(b*(b+1));
 let MR=major(p,b,R),M1=major(p,b,1n),bin=BigInt((m+3)*(m+2)/2),num=bin*R*BigInt(m+2),de=R**BigInt(m+1)*(R*BigInt(m+2)-BigInt(m+4));if(de<=0n)throw Error('Cauchy ratio');let rem=ceil(MR.total*num,de)+M1.tail;
 if(cfg.maxRemainder&&rem>dec(cfg.maxRemainder))return {p,b,status:'FAIL_REMAINDER',remainder:fmt(rem),majorant:fmt(MR.total),tail:fmt(M1.tail),minDen:fmt(MR.minDen),seconds:Number(process.hrtime.bigint()-start)/1e9};
 let pa=dense(),p0=dense(S),qa=dense(S),q0=dense(),P=dense(),Q=dense();
 for(let j=0;j<K;j++){let l=BigInt(j),beta=plus(con((BigInt(2*p-1)-2n*l)*S/2n),inu);let rn=nr(l),rd=dr(l),rnp=j?nr(l-1n):one,rdp=j?dr(l-1n):one;
 let den=plus(one,neg(prod(h,prod(beta,plus(beta,con(-s*S))))));let A=prod(prod(plus(d,prod(h,prod(beta,plus(con((s-1n)*S),times(beta,-2n))))),rn),rdp),BB=prod(prod(prod(h,prod(beta,plus(beta,one))),rn),rnp),D=prod(prod(den,rd),rdp);
 if(p===3){let A0=plus(con((9n-2n*BigInt(b*(b-1)))*S/2n),neg(half(prod(beta,plus(beta,con(BigInt(8-4*b)*S))))));let A1=half(prod(beta,plus(beta,con(BigInt(9-4*b)*S))));let ca=plus(one,times(h,q),prod(h,A0)),cb=prod(h,A1);let vp=j?solve(rnp,pa,[],pa,rdp):pa,vq=j?solve(rnp,qa,[],qa,rdp):qa;addprod(P,p0,ca);addprod(P,vp,cb);addprod(Q,q0,ca);addprod(Q,vq,cb)}else{addprod(P,p0,one);addprod(Q,q0,one)}
 if(j+1<K){let pn=solve(A,p0,BB,pa,D),qn=solve(A,q0,BB,qa,D);pa=p0;p0=pn;qa=q0;q0=qn}}
 let rawP=signedNorm(P)+rem,rawQ=du(signedNorm(Q)+rem,d0-dd);let ubP=cfg.hSplits?splitHNorm(P,cfg.hSplits)+rem:rawP,ubQ=cfg.hSplits?du(splitHNorm(Q,cfg.hSplits)+rem,d0-dd):rawQ;let status=cfg.maxTransfer&&(ubP>dec(cfg.maxTransfer)||ubQ>dec(cfg.maxTransfer))?'FAIL_TRANSFER':'PASS';
 return {p,b,status,P_bound:fmt(ubP),Q_bound:fmt(ubQ),P_unsplit:fmt(rawP),Q_unsplit:fmt(rawQ),P_center:[fmt(P.r[0]),fmt(P.i[0])],Q_unscaled_center:[fmt(Q.r[0]),fmt(Q.i[0])],remainder:fmt(rem),rounding:[fmt(P.err),fmt(Q.err)],majorant:fmt(MR.total),unitMajorant:fmt(M1.total),tail:fmt(M1.tail),minDen:fmt(MR.minDen),coeffs:N,seconds:Number(process.hrtime.bigint()-start)/1e9};}
 return {run,fullTwist,hRequired:fmt(hRequired)};
}
let out=[];for(let cfg of cases){let prep=setup(cfg);for(let [p,b] of cfg.channels||[[3,0],[3,4],[4,0],[4,4]]){let res;try{res=prep.run(p,b)}catch(e){res={p,b,status:'ERROR',error:String(e)}}let row={name:cfg.name,box:cfg,fullTwist:prep.fullTwist,hRequired:prep.hRequired,...res};out.push(row);console.log(JSON.stringify(row));}}
if(process.argv[3])fs.writeFileSync(process.argv[3],JSON.stringify(out,null,2)+'\n');
