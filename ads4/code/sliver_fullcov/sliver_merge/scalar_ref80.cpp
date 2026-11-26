#include <immintrin.h>
#include <gmpxx.h>
#include <vector>
#include <array>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <algorithm>
#include <cstdint>
using Z=mpz_class;using V=__m512i;using Dv=__m512d;using Mask=__mmask8;
constexpr uint64_t B=(1ULL<<52),BM=B-1;constexpr int LM=64;
const V zero=_mm512_setzero_si512(),one=_mm512_set1_epi64(1),mask52=_mm512_set1_epi64(BM),sign52=_mm512_set1_epi64(B/2);
inline V spl(uint64_t x){return _mm512_set1_epi64(x);} inline V low(V x){return _mm512_and_si512(x,mask52);} inline V add(V x,V y){return _mm512_add_epi64(x,y);} inline V sub(V x,V y){return _mm512_sub_epi64(x,y);} inline V hi(V x){return _mm512_srli_epi64(x,52);} inline V mul(V x,V y){return _mm512_mullo_epi64(x,y);} inline V bit(Mask x){return _mm512_maskz_set1_epi64(x,1);} inline Mask ult(V x,V y){return _mm512_cmp_epu64_mask(x,y,_MM_CMPINT_LT);} inline Mask uge(V x,V y){return _mm512_cmp_epu64_mask(x,y,_MM_CMPINT_GE);} inline Mask eq(V x,V y){return _mm512_cmpeq_epi64_mask(x,y);} inline Mask neg(V x){return _mm512_cmp_epi64_mask(x,zero,_MM_CMPINT_LT);} 
struct Divisor{V d,dl,dh; Dv inv;Divisor(V x):d(x),dl(low(x)),dh(hi(x)),inv(_mm512_div_pd(_mm512_set1_pd(1.0),_mm512_cvtepu64_pd(x))){};
 inline V divstep(V &r,V a) const {
  Dv xx=_mm512_fmadd_pd(_mm512_cvtepu64_pd(r),_mm512_set1_pd(double(B)),_mm512_cvtepu64_pd(a));
  V q=_mm512_min_epu64(_mm512_cvttpd_epu64(_mm512_mul_pd(xx,inv)),mask52);
  V pl=low(_mm512_madd52lo_epu64(zero,q,dl));V ph=_mm512_madd52hi_epu64(mul(q,dh),q,dl);
  V rl=low(sub(a,pl)),rh=sub(sub(r,ph),bit(ult(a,pl)));
  Mask m=neg(rh);int it=0;
  while(m){if(++it>8)throw std::runtime_error("division negative correction");q=_mm512_mask_sub_epi64(q,m,q,one);V z=_mm512_mask_add_epi64(rl,m,rl,dl);rh=_mm512_mask_add_epi64(rh,m,rh,add(dh,hi(z)));rl=low(z);m=neg(rh);}
  m=ult(dh,rh)|(eq(dh,rh)&uge(rl,dl));it=0;
  while(m){if(++it>8)throw std::runtime_error("division positive correction");q=_mm512_mask_add_epi64(q,m,q,one);V br=bit(ult(rl,dl));rh=_mm512_mask_sub_epi64(rh,m,rh,add(dh,br));rl=low(_mm512_mask_sub_epi64(rl,m,rl,dl));m=ult(dh,rh)|(eq(dh,rh)&uge(rl,dl));}
  r=add(_mm512_slli_epi64(rh,52),rl);return q;
 }
};
struct Coef{int L,nt;std::array<uint64_t,LM> a,t;bool ts;};
struct Table{std::vector<Coef>v[2];};
int bits(const Z&v){return mpz_sizeinbase(v.get_mpz_t(),2);}
Table tables(const std::vector<Z>&A,const std::vector<Z>&G,int K,int P,const std::vector<int>&sh){Table R; Z S,sc=Z(1)<<P;mpz_ui_pow_ui(S.get_mpz_t(),10,60);std::vector<Z> AA(K+1),T(K+1);for(int k=0;k<=K;k++){AA[k]=A[k]*sc;AA[k]/=S;mpz_tdiv_q_2exp(AA[k].get_mpz_t(),AA[k].get_mpz_t(),sh[k]);}
 for(int i=0;i<2;i++){int b=4*i,q=b*(b+1),sg=i?-1:1;R.v[i].resize(K+1);for(int k=0;k<=K;k++){long k1=(long)(k+1)*(k+1),k2=-(long)k*k+(8-4*b)*(long)k+9-2*b*(b-1)+2*q;T[k]=k1*A[k+1]+k2*A[k]+2*sg*G[k];T[k]*=sc;T[k]/=S;mpz_tdiv_q_2exp(T[k].get_mpz_t(),T[k].get_mpz_t(),sh[k]);}
  Z mx=0;for(int k=K;k>=0;k--){Z v=abs(AA[k])+abs(T[k])+10;if(v>mx)mx=v;Z bd=2*(K+1)*mx;int L=(bits(bd)+1+51)/52+1;if(L>LM-3)throw std::runtime_error("limb length");Coef &x=R.v[i][k];x.L=L;size_t cnt;Z aa=AA[k];if(aa<0)aa+=Z(1)<<(52*L);mpz_export(x.a.data(),&cnt,-1,8,0,12,aa.get_mpz_t());if(cnt>(size_t)L)throw std::runtime_error("a export");x.ts=T[k]<0;Z tt=abs(T[k]);mpz_export(x.t.data(),&cnt,-1,8,0,12,tt.get_mpz_t());x.nt=cnt;if(cnt>(size_t)L)throw std::runtime_error("t export");}
 }return R;}
void eval8(unsigned long C,const uint64_t*hv,int K,int s,const std::vector<int>&sh,const std::vector<Coef>&tab,Z*res){
 V H=_mm512_loadu_si512(hv),CC=spl(C),v[LM],prod[LM+2],Q[LM+2],TP[LM+2];int Lv=1;v[0]=zero;
 for(int k=K;k>=0;k--){const Coef &x=tab[k];int L=x.L;
  // Magnitude of the previous signed value and multiplication by the unscaled ratio numerator.
  Mask sv=uge(v[Lv-1],sign52),sr=0;V d=spl(1ULL<<61),rn=zero;
  if(k<K){V z=mul(H,spl((uint64_t)k*(k+s)));d=sub(d,z);sr=ult(CC,z);rn=_mm512_mask_sub_epi64(sub(CC,z),sr,z,CC);}
  Mask sq=sv^sr;V signv=_mm512_maskz_set1_epi64(sv,BM),carry=zero,cy=bit(sv);V rnl=low(rn),rnh=hi(rn);
  for(int j=0;j<Lv;j++){V a=add(_mm512_xor_si512(v[j],signv),cy);cy=hi(a);a=low(a);V w=_mm512_madd52lo_epu64(carry,a,rnl);carry=add(_mm512_madd52hi_epu64(mul(a,rnh),a,rnl),hi(w));prod[j]=low(w);}
  Divisor dd(d);V rem=carry;for(int j=Lv-1;j>=0;j--)Q[j]=dd.divstep(rem,prod[j]);
  int delta=k<K?sh[k+1]-sh[k]:0;
  // Multiply the dressing coefficient magnitude by H and divide by 2^62.
  carry=zero;for(int j=0;j<x.nt;j++){V a=spl(x.t[j]);V w=_mm512_madd52lo_epu64(carry,a,H);carry=add(_mm512_madd52hi_epu64(zero,a,H),hi(w));TP[j]=low(w);}TP[x.nt]=carry;TP[x.nt+1]=zero;
  // Signed scaled quotient + signed dressing + signed A, all in base 2^52 two's complement.
  V qsign=_mm512_maskz_set1_epi64(sq,BM),tsign=x.ts?mask52:zero;V qcy=bit(sq),tcy=x.ts?one:zero,ac=zero,shiftcarry=zero;
  V sl=spl(delta),srsh=spl(52-delta);
  for(int j=0;j<L;j++){
    V old=j<Lv?Q[j]:zero;V qq=low(_mm512_or_si512(_mm512_sllv_epi64(old,sl),shiftcarry));shiftcarry=_mm512_srlv_epi64(old,srsh);
    V qt=add(_mm512_xor_si512(qq,qsign),qcy);qcy=hi(qt);qt=low(qt);
    V tt=j+1<=x.nt?_mm512_srli_epi64(TP[j+1],10):zero; if(j+2<=x.nt)tt=_mm512_or_si512(tt,low(_mm512_slli_epi64(TP[j+2],42)));
    tt=add(_mm512_xor_si512(tt,tsign),tcy);tcy=hi(tt);tt=low(tt);
    V w=add(add(add(qt,tt),spl(x.a[j])),ac);v[j]=low(w);ac=hi(w);
  } Lv=L;
 }
 alignas(64) uint64_t out[LM][8];for(int j=0;j<Lv;j++)_mm512_store_si512(out[j],v[j]);for(int z=0;z<8;z++){uint64_t limbs[LM];for(int j=0;j<Lv;j++)limbs[j]=out[j][z];mpz_import(res[z].get_mpz_t(),Lv,-1,8,0,12,limbs);if(limbs[Lv-1]>=B/2)res[z]-=Z(1)<<(52*Lv);}
}
int main(int argc,char**argv){try{std::ifstream f(argv[1]);int N,D;f>>N>>D;if(D!=60)throw std::runtime_error("scale");std::vector<Z>A(N+1),G(N+1);Z r,g;for(int k=0;k<=N;k++){f>>A[k]>>r>>G[k]>>g;if(r<0||g<0||r>2||g>2)throw std::runtime_error("radius");}
 int K;unsigned long Den;int Ny,Nc,P;std::cin>>K>>Den>>Ny>>Nc>>P;if(K>=N||K<0||Den!=(1UL<<61)||Ny<1||Nc<1||(P<64||P>180))throw std::runtime_error("settings");std::vector<int>sh(K+1);for(int &x:sh)std::cin>>x;if(sh[0]!=0)throw std::runtime_error("shift0");for(int k=0;k<K;k++)if(sh[k+1]<sh[k]||sh[k+1]-sh[k]>2)throw std::runtime_error("shift");auto T=tables(A,G,K,P,sh);Z Yden;std::cin>>Yden;std::vector<Z>Y(Ny);for(auto&y:Y){std::cin>>y;if(y<0||y>Yden)throw std::runtime_error("node");}
 // Reference polynomials in x=(c-l)/(u-l), y=h/H(c), with integer coefficients at scale 1e20.
 Z ln,ld,un,ud;std::cin>>ln>>ld>>un>>ud; if(ld<=0||ud<=0||ln*ud>=un*ld)throw std::runtime_error("reference endpoints");
 Z ref[2][3][3];for(int i=0;i<2;i++)for(int j=0;j<3;j++)for(int k=0;k<3;k++)std::cin>>ref[i][j][k];
 Z sc=Z(1)<<P,unit;mpz_ui_pow_ui(unit.get_mpz_t(),10,20);Z de=Z(Den)*Den*Den*sc;Z lows[2],highs[2];bool init[2]={false,false};unsigned long long total=0;Z Hbig,lo,hi,ztemp,hn,hd;
 for(int ic=0;ic<Nc;ic++){unsigned long C;std::cin>>C>>hn>>hd;if(C>=Den||C==0||hn<=0||hd<=0)throw std::runtime_error("C/H");Z scaleC=Z(Den-C)*(Den-C)*(Den-C)*unit;Z hnum=hn*Den,hden=hd*Yden;Z xnum=(Z(C)*ld-ln*Den)*ud,xden=Z(Den)*(un*ld-ln*ud);
 if(xnum<0||xnum>xden)throw std::runtime_error("x outside interval");
 Z Rscale=Z(1)<<50, X=xnum*Rscale/xden, XR[2][3];
 for(int i=0;i<2;i++)for(int k=0;k<3;k++){XR[i][k]=ref[i][2][k];XR[i][k]=XR[i][k]*X/Rscale+ref[i][1][k];XR[i][k]=XR[i][k]*X/Rscale+ref[i][0][k];}
 std::vector<uint64_t>H(Ny);for(int j=0;j<Ny;j++){Hbig=hnum*Y[j]/hden;if(!Hbig.fits_ulong_p())throw std::runtime_error("h overflow");H[j]=Hbig.get_ui();unsigned __int128 z=(unsigned __int128)H[j]*(K-1)*(K+8);if(z>=Den||2*z>(unsigned __int128)Den+C)throw std::runtime_error("ratio or pole");}
  for(int j=0;j<Ny;j+=8){alignas(64)uint64_t hv[8];for(int z=0;z<8;z++)hv[z]=H[std::min(j+z,Ny-1)];for(int i=0;i<2;i++){Z vals[8];eval8(C,hv,K,1+8*i,sh,T.v[i],vals);for(int z=0;z<8&&j+z<Ny;z++){ztemp=vals[z]*scaleC;mpz_fdiv_q(lo.get_mpz_t(),ztemp.get_mpz_t(),de.get_mpz_t());mpz_cdiv_q(hi.get_mpz_t(),ztemp.get_mpz_t(),de.get_mpz_t());Z yv=Z(H[j+z])*hd*Rscale/(Z(Den)*hn);if(yv<0||yv>Rscale)throw std::runtime_error("y outside interval"); Z rr=XR[i][2];rr=rr*yv/Rscale+XR[i][1];rr=rr*yv/Rscale+XR[i][0];lo-=rr;hi-=rr;
 if(!init[i]||lo<lows[i])lows[i]=lo;if(!init[i]||hi>highs[i])highs[i]=hi;init[i]=true;total++;if(argc>2)std::cerr<<ic<<" "<<j+z<<" "<<4*i<<" "<<vals[z]<<"\n";}}}
  if((ic+1)%20==0)std::cerr<<"rows "<<ic+1<<"/"<<Nc<<" nodes "<<total<<"\n";
 }for(int i=0;i<2;i++)std::cout<<4*i<<" "<<lows[i]<<" "<<highs[i]<<"\n";std::cout<<"COUNT "<<total<<"\n";
 }catch(const std::exception&e){std::cerr<<e.what()<<"\n";return 1;}return 0;}
