#include <gmpxx.h>
#include <vector>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <algorithm>
using Z=mpz_class;
int main(int argc,char**argv){try{
 std::ifstream f(argv[1]);int N,D;f>>N>>D;if(D!=60)throw std::runtime_error("scale");std::vector<Z>A(N+1),G(N+1),AA,T[2];T[0].resize(N);T[1].resize(N);Z r,g,S;mpz_ui_pow_ui(S.get_mpz_t(),10,60);
 for(int k=0;k<=N;k++){f>>A[k]>>r>>G[k]>>g;if(r<0||g<0||r>2||g>2)throw std::runtime_error("radius");}
 int K;unsigned long Den;int Ny,Nc,P;std::cin>>K>>Den>>Ny>>Nc>>P; if(K>=N||K<0||Den!=(1UL<<61)||Ny<1||Nc<1||P!=180)throw std::runtime_error("settings");
 std::vector<int> shift(K+1);for(int &s:shift)std::cin>>s;if(shift[0]!=0)throw std::runtime_error("shift0");for(int k=0;k<K;k++)if(shift[k+1]<shift[k]||shift[k+1]-shift[k]>2)throw std::runtime_error("shift");
 Z sc=Z(1)<<P; AA.resize(K+1);for(int k=0;k<=K;k++){AA[k]=A[k]*sc;AA[k]/=S;mpz_tdiv_q_2exp(AA[k].get_mpz_t(),AA[k].get_mpz_t(),shift[k]);}
 for(int i=0;i<2;i++){int b=4*i,q=b*(b+1),sg=i?-1:1;for(int k=0;k<=K;k++){long k1=(long)(k+1)*(k+1),k2=-(long)k*k+(8-4*b)*(long)k+9-2*b*(b-1)+2*q;T[i][k]=k1*A[k+1]+k2*A[k]+2*sg*G[k];T[i][k]*=sc;T[i][k]/=S;mpz_tdiv_q_2exp(T[i][k].get_mpz_t(),T[i][k].get_mpz_t(),shift[k]);}}
 Z Yden;std::cin>>Yden;std::vector<Z>Y(Ny);for(auto&y:Y){std::cin>>y;if(y<0||y>Yden)throw std::runtime_error("node");}
 Z unit;mpz_ui_pow_ui(unit.get_mpz_t(),10,20);Z den3=Z(Den)*Den*Den;Z de=den3*sc;Z lows[2]={0,0},highs[2]={0,0};bool init[2]={false,false};Z v,w,hn,hd,cp,scaleC,lo,hi,ztemp;unsigned long long total=0;
 for(int ic=0;ic<Nc;ic++){
  unsigned long C;std::cin>>C>>hn>>hd;if(C>=Den||C==0||hn<=0||hd<=0)throw std::runtime_error("C/H");cp=Z(Den-C)*(Den-C)*(Den-C);scaleC=cp*unit;Z hnum=hn*Den;Z hden=hd*Yden;
  for(auto&y:Y){Z Hbig=hnum*y/hden;if(!Hbig.fits_ulong_p())throw std::runtime_error("h overflow");unsigned long H=Hbig.get_ui();
   for(int i=0;i<2;i++){int s=1+8*i;v=0;
    for(long k=K;k>=0;k--){unsigned long p=k*(k+s);if(k<K){if(H>0&&p>(Den-1)/H)throw std::runtime_error("pole");unsigned long z=H*p,rd=Den-z;bool negative=z>C;unsigned long rn=negative?z-C:C-z;if(rn>rd)throw std::runtime_error("node ratio exceeds one");rn<<=(shift[k+1]-shift[k]);mpz_mul_ui(v.get_mpz_t(),v.get_mpz_t(),rn);mpz_tdiv_q_ui(v.get_mpz_t(),v.get_mpz_t(),rd);if(negative)mpz_neg(v.get_mpz_t(),v.get_mpz_t());}
     mpz_mul_ui(w.get_mpz_t(),T[i][k].get_mpz_t(),H);mpz_tdiv_q_2exp(w.get_mpz_t(),w.get_mpz_t(),62);v+=w;v+=AA[k];}
    ztemp=v*scaleC;mpz_fdiv_q(lo.get_mpz_t(),ztemp.get_mpz_t(),de.get_mpz_t());mpz_cdiv_q(hi.get_mpz_t(),ztemp.get_mpz_t(),de.get_mpz_t());
    if(!init[i]||lo<lows[i])lows[i]=lo;if(!init[i]||hi>highs[i])highs[i]=hi;init[i]=true;total++;
   }
  }
  if((ic+1)%20==0)std::cerr<<"rows "<<ic+1<<"/"<<Nc<<" nodes "<<total<<"\n";
 }
 for(int i=0;i<2;i++)std::cout<<4*i<<" "<<lows[i]<<" "<<highs[i]<<"\n";std::cout<<"COUNT "<<total<<"\n";
 }catch(const std::exception&e){std::cerr<<e.what()<<"\n";return 1;}return 0;}
