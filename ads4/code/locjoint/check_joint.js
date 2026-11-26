#!/usr/bin/env node
'use strict';
// Node >=18, no dependencies. All certificate arithmetic is BigInt.
// CSV values are integers scaled by 10^16, except index columns.
const fs=require('fs'),path=require('path');
const out=process.argv[2]||path.join(__dirname,'output');fs.mkdirSync(out,{recursive:true});
const S=10n**16n,a=998n,J=1500n;
const fd=(n,d)=>{if(d<=0n)throw Error('denominator');return n>=0n?n/d:-((-n+d-1n)/d)},cd=(n,d)=>-fd(-n,d),ml=(x,y)=>fd(x*y,S),mh=(x,y)=>cd(x*y,S);
const Q=(x,d=1)=>fd(BigInt(x)*S,BigInt(d)),min=(a,b)=>a<b?a:b,max=(a,b)=>a>b?a:b;
const LN2L=6931471805599453n,LN2U=6931471805599454n,PIU=31415926535897933n;
const check=(ok,why)=>{if(!ok)throw Error(why)};
function lograt(n,d){check(n>0n&&d>0n,'log domain');let k=0n;while(n<d){n*=2n;k--}while(n>=2n*d){d*=2n;k++}let vl=fd((n-d)*S,n+d),vu=cd((n-d)*S,n+d),zl=ml(vl,vl),zu=mh(vu,vu),pl=vl,pu=vu,sl=vl,su=vu;for(let j=3n;j<=39n;j+=2n){pl=ml(pl,zl);pu=mh(pu,zu);sl+=fd(pl,j);su+=cd(pu,j)}return [2n*sl+(k>=0n?k*LN2L:k*LN2U),2n*su+2n+(k>=0n?k*LN2U:k*LN2L)]}
const lnlo=x=>lograt(x,S)[0],lnhi=x=>lograt(x,S)[1];
function isqrt(n){check(n>=0n,'sqrt domain');if(n<2n)return n;let x=1n<<BigInt(Math.ceil(n.toString(2).length/2));while(true){let y=(x+n/x)>>1n;if(y>=x)return x;x=y}}
const sqrtlo=x=>isqrt(x*S),sqrthi=x=>isqrt(x*S)+1n;
function atanSmall(n,d){check(n>=0n&&2n*n<=d,'atan reduced domain');let l=fd(n*S,d),u=cd(n*S,d),zl=ml(l,l),zu=mh(u,u),pl=l,pu=u,al=l,au=u;for(let j=1n;j<=23n;j++){pl=ml(pl,zl);pu=mh(pu,zu);if(j%2n){al-=cd(pu,2n*j+1n);if(j<23n)au-=fd(pl,2n*j+1n)}else{al+=fd(pl,2n*j+1n);au+=cd(pu,2n*j+1n)}}return [al,au]}
function atanI(n,d){check(n>=0n&&n<=d&&d>0n,'atan domain');if(2n*n<=d)return atanSmall(n,d);let [l,u]=atanSmall(d-n,d+n);return [fd(PIU-1n,4n)-u,cd(PIU,4n)-l]}
function atanU(n,d){check(n>=0n&&d>0n,'atanU domain');if(n===0n)return 0n;if(n>d){let [l]=atanI(d,n);return cd(PIU,2n)-l}return atanI(n,d)[1]}
function dec(x){let neg=x<0n;if(neg)x=-x;return (neg?'-':'')+(x/S).toString()+'.'+(x%S).toString().padStart(16,'0')}
function coef(V){let lc=0n,lmax=0n;for(let j=0n;j<J;j++){let x=j+1n,A=(2n*j+1n)**2n;let l0=lograt(A,4n*x*x)[1],l1=lograt(A+4n*V*V,4n*x)[1]-fd(lograt(x*x+4n*V*V,1n)[0],2n);lc+=max(l0,l1);lmax=max(lmax,lc)}return lmax}
const xx=J+1n,Y=1150n**2n;check(xx*xx*(32n*Y-24n)>(4n*Y+1n)**2n,'tail ratio endpoint');
const Vs=[999n,1040n,1150n],coefLogs=Vs.map(coef),Cs=coefLogs.map(x=>cd(x,S));check(Cs[0]<=547n&&Cs[1]<=570n&&Cs[2]<=630n,'coefficient budgets');
const logN=lograt(100n,1n)[1];
const runs=[
{id:'A',begin:Q(12,1000),end:Q(12,100),tNum:28n,tDen:25n,gNum:29n,gDen:25n,lowerBeta:false,localPhase:false,tailCut:Q(99,100),tailExtra:15n,vs:[1150n],cs:[630n],targets:[-90n]},
{id:'B',begin:Q(12,100),end:Q(20,100),tNum:41n,tDen:40n,gNum:53n,gDen:50n,lowerBeta:true,localPhase:true,tailCut:Q(999,1000),tailExtra:20n,vs:[1150n],cs:[630n],targets:[-100n]},
{id:'C',begin:Q(20,100),end:Q(26,100),tNum:203n,tDen:200n,gNum:26n,gDen:25n,lowerBeta:true,localPhase:true,tailCut:Q(999,1000),tailExtra:20n,vs:Vs,cs:[547n,570n,630n],targets:[-177n,-129n,-2n]}];
const report={scale:S.toString(),J:Number(J),a_min:Number(a),coefficient_bounds:Vs.map((v,i)=>({V:Number(v),log_upper_scaled:coefLogs[i].toString(),log_upper:dec(coefLogs[i]),integer_ceiling:Cs[i].toString()})),runs:[]};
function cuts(th,cfg){let q=[0n],v=cd(th,4n);while(v<S){q.push(v);v=cd(cfg.gNum*v,cfg.gDen)}q.push(S);return q}
function doSlab(cfg,tl,th,slab,fdout){
 const u=cuts(th,cfg),x=u,bs=[],ws=[],worst=cfg.vs.map(()=>-(10n**9n)*S),locations=cfg.vs.map(()=>null);
 const mu=fd(tl*S,S+cd(9n*th,998n)),kmU=cd((S-mu)*S,mu),kmL=fd((S-th)*S,th),del=cd(tl,100n);
 for(let i=0;i<u.length-1;i++){const z=u[i],z1=u[i+1];let H=0n;if(z>th)H=lograt(th,z)[0]+ml(kmL,lograt(S-th,S-z)[0]);else if(cfg.lowerBeta&&z1<mu)H=lograt(mu,z1)[0]+ml(kmU,lograt(S-mu,S-z1)[0]);bs.push(min(0n,-a*H))}
 for(let j=0;j<x.length-1;j++){const z=x[j],k=lograt(th*th*S+4n*(S-th)*z*z,th*th*S)[0];ws.push(k<cd(S,a)?0n:min(0n,logN-lnlo(tl)+Q(502,1000)-fd((a+1n)*k,2n)))}
 let rows=[];
 for(let i=0;i<u.length-1;i++)for(let j=0;j<x.length-1;j++){
 const u0=u[i],u1=u[i+1],x0=x[j],x1=x[j+1],r0=ml(S-u1,S-th),r1=mh(S-u0,S-tl);
 const lo=min(S,sqrtlo(ml(S-r1,S-r1)+4n*ml(r0,ml(x0,x0))));const q2=r=>mh(S-r,S-r)+4n*mh(r,mh(x1,x1));const hi=sqrthi(max(q2(r0),q2(r1)));
 const A=atanU(2n*mh(r1,mh(x1,sqrthi(S-ml(x0,x0)))),S-r1+2n*ml(r0,ml(x0,x0)));
 let E=cd(PIU,2n);if(cfg.localPhase){const coslo=S-2n*mh(x1,x1),relo=min(ml(r0,coslo),ml(r1,coslo))-del,R=r1+del;E=atanU(sqrthi(2n*(R-relo)),S-R)}
 const ll=-4n*lnlo(lo),lhi=lnhi(max(S,hi)),ltail=hi<=cfg.tailCut?lnhi(hi):null;
 const row=[slab,i,j,tl,th,u0,u1,x0,x1,r0,r1,lo,hi,A,E,bs[i],ws[j]];
 for(let h=0;h<cfg.vs.length;h++){const V=cfg.vs[h],C=cfg.cs[h];const ft=S+max((cfg.localPhase?2n:4n)*S+(cfg.localPhase?V*E:cd(V*PIU,2n)),(C+16n)*S+V*A+J*lhi)+ll;const tail=ltail!==null?(C+cfg.tailExtra)*S+V*A+J*ltail+ll:null;const v=bs[i]+ws[j]+(tail===null?ft:min(ft,tail));if(v>worst[h]){worst[h]=v;locations[h]=[i,j]}row.push(ft,tail===null?'NA':tail,v)}rows.push(row.join(','))}
 fs.writeSync(fdout,rows.join('\n')+'\n');const count=(u.length-1)**2;return {slab,t_lo:tl.toString(),t_hi:th.toString(),t_lo_decimal:dec(tl),t_hi_decimal:dec(th),beta_cells:u.length-1,angular_cells:x.length-1,count,max_cell_scaled:worst.map(String),max_cell_decimal:worst.map(dec),max_locations:locations};
}
let total=0;
for(const cfg of runs){const beginTime=process.hrtime.bigint(),csv=path.join(out,'cells_'+cfg.id+'.csv'),fdout=fs.openSync(csv,'w');fs.writeSync(fdout,['slab','beta_index','angle_index','t_lo','t_hi','v_lo','v_hi','x_lo','x_hi','r_lo','r_hi','q_lo','q_hi','arg_upper','conical_phase_upper','log_beta_upper','log_angle_upper',...cfg.vs.flatMap(v=>['log_FT_'+v,'log_tail_'+v,'log_weighted_'+v])].join(',')+'\n');let tl=cfg.begin,slab=0,count=0,slabs=[],maxCount=0,worst=cfg.vs.map(()=>-(10n**9n)*S);while(tl<cfg.end){const th=min(cfg.end,cd(cfg.tNum*tl,cfg.tDen));const r=doSlab(cfg,tl,th,slab,fdout);slabs.push(r);count+=r.count;maxCount=Math.max(maxCount,r.count);worst=worst.map((w,i)=>max(w,BigInt(r.max_cell_scaled[i])));console.log(JSON.stringify({run:cfg.id,...r}));tl=th;slab++}fs.closeSync(fdout);const totals=worst.map(w=>w+lograt(BigInt(maxCount),1n)[1]+lograt(8n,5n)[1]);const pass=worst.every((w,i)=>w<cfg.targets[i]*S);check(pass,'joint cell target failed: '+cfg.id);const summary={run:cfg.id,t_slabs:slab,joint_cells:count,max_cells_per_slab:maxCount,V:cfg.vs.map(Number),max_cell_log:worst.map(dec),max_cell_log_scaled:worst.map(String),log_kernel_upper:totals.map(dec),log_kernel_upper_scaled:totals.map(String),pass,elapsed_seconds:Number(process.hrtime.bigint()-beginTime)/1e9,slabs};report.runs.push(summary);total+=count;console.log('SUMMARY '+JSON.stringify({...summary,slabs:undefined}))}
report.total_joint_cells=total;check(total===190031,'unexpected cell count');check(BigInt(report.runs[0].log_kernel_upper_scaled[0])<-85n*S,'A kernel');check(BigInt(report.runs[1].log_kernel_upper_scaled[0])<-95n*S,'B kernel');for(const [i,b] of [-168n,-120n,7n].entries())check(BigInt(report.runs[2].log_kernel_upper_scaled[i])<b*S,'C kernel');report.geometry_and_kernel_pass=true;fs.writeFileSync(path.join(out,'summary.json'),JSON.stringify(report,null,2)+'\n');console.log('PASS: '+total+' joint cells; geometry and kernel budgets passed. See README for analytic assumptions and spectral inputs.');
