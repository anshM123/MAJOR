// Portable GMP recurrence for independently reproducing the SIMD output.
// Every division truncates toward zero except the explicitly directed final bounds.
#include <gmpxx.h>
#include <vector>
#include <array>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <algorithm>
using Z=mpz_class;
Z trunc2(Z a,unsigned long k){mpz_tdiv_q_2exp(a.get_mpz_t(),a.get_mpz_t(),k);return a;}
int main(int argc,char**argv){try{
 if(argc<2)throw std::runtime_error("coefficient path");std::ifstream f(argv[1]);int N,D;f>>N>>D;if(D!=60)throw std::runtime_error("scale");std::vector<Z>A(N+1),G(N+1);Z r,g,S,unit;mpz_ui_pow_ui(S.get_mpz_t(),10,60);mpz_ui_pow_ui(unit.get_mpz_t(),10,20);
 for(int k=0;k<=N;k++){f>>A[k]>>r>>G[k]>>g;if(!f||r<0||r>2||g<0||g>2)throw std::runtime_error("table");}
 int K,Ny,Nc,P;unsigned long Den;std::cin>>K>>Den>>Ny>>Nc>>P;if(K>=N||K<0||Den!=(1UL<<61)||Ny<1||Nc<1||P<64||P>180)throw std::runtime_error("settings");std::vector<int>sh(K+1);for(int &v:sh)std::cin>>v;if(sh[0]!=0)throw std::runtime_error("shift0");for(int k=0;k<K;k++)if(sh[k+1]<sh[k]||sh[k+1]-sh[k]>2)throw std::runtime_error("shift");
 Z sc=Z(1)<<P;std::vector<Z>AA(K+1);std::array<std::vector<Z>,2>T;T[0].resize(K+1);T[1].resize(K+1);
 for(int k=0;k<=K;k++){AA[k]=trunc2(A[k]*sc/S,sh[k]);for(int i=0;i<2;i++){int b=4*i,q=b*(b+1),sg=i?-1:1;long k1=(long)(k+1)*(k+1),k2=-(long)k*k+(8-4*b)*(long)k+9-2*b*(b-1)+2*q;T[i][k]=trunc2((k1*A[k+1]+k2*A[k]+2*sg*G[k])*sc/S,sh[k]);}}
 Z Yden;std::cin>>Yden;std::vector<Z>Y(Ny);for(auto&y:Y){std::cin>>y;if(y<0||y>Yden)throw std::runtime_error("node");}
 Z ln,ld,un,ud;std::cin>>ln>>ld>>un>>ud;if(ld<=0||ud<=0||ln*ud>=un*ld)throw std::runtime_error("reference endpoints");Z ref[2][3][3];for(int i=0;i<2;i++)for(int j=0;j<3;j++)for(int k=0;k<3;k++)std::cin>>ref[i][j][k];
 Z de=Z(Den)*Den*Den*sc,lows[2],highs[2];bool init[2]={false,false};unsigned long long total=0;
 for(int ic=0;ic<Nc;ic++){unsigned long C;Z hn,hd;std::cin>>C>>hn>>hd;if(C==0||C>=Den||hn<=0||hd<=0)throw std::runtime_error("row");Z scaleC=Z(Den-C)*(Den-C)*(Den-C)*unit;Z hnum=hn*Den,hden=hd*Yden;Z xnum=(Z(C)*ld-ln*Den)*ud,xden=Z(Den)*(un*ld-ln*ud);if(xnum<0||xnum>xden)throw std::runtime_error("x");Z Rscale=Z(1)<<50,X=xnum*Rscale/xden,XR[2][3];
 for(int i=0;i<2;i++)for(int k=0;k<3;k++){XR[i][k]=ref[i][2][k];XR[i][k]=XR[i][k]*X/Rscale+ref[i][1][k];XR[i][k]=XR[i][k]*X/Rscale+ref[i][0][k];}
 for(int j=0;j<Ny;j++){Z H=hnum*Y[j]/hden;Z pmax=Z(K-1)*(K+8);if(H*pmax>=Den||2*H*pmax>Z(Den)+C)throw std::runtime_error("ratio or pole");Z yv=H*hd*Rscale/(Z(Den)*hn);if(yv<0||yv>Rscale)throw std::runtime_error("y");
 for(int i=0;i<2;i++){Z v=0;int s=1+8*i;
 for(int k=K;k>=0;k--){if(k<K){Z z=H*k*(k+s),d=Z(Den)-z,rn=Z(C)-z;if(d<=0||abs(rn)>d)throw std::runtime_error("ratio");v*=rn;v/=d;v<<=(sh[k+1]-sh[k]);}v+=AA[k]+trunc2(T[i][k]*H,62);}
 Z n=v*scaleC,lo,hi;mpz_fdiv_q(lo.get_mpz_t(),n.get_mpz_t(),de.get_mpz_t());mpz_cdiv_q(hi.get_mpz_t(),n.get_mpz_t(),de.get_mpz_t());Z rr=XR[i][2];rr=rr*yv/Rscale+XR[i][1];rr=rr*yv/Rscale+XR[i][0];lo-=rr;hi-=rr;if(!init[i]||lo<lows[i])lows[i]=lo;if(!init[i]||hi>highs[i])highs[i]=hi;init[i]=true;total++;if(argc>2)std::cerr<<ic<<" "<<j<<" "<<4*i<<" "<<v<<"\n";
 }}if((ic+1)%20==0)std::cerr<<"rows "<<ic+1<<"/"<<Nc<<" nodes "<<total<<"\n";}
 if(!std::cin)throw std::runtime_error("input read");for(int i=0;i<2;i++)std::cout<<4*i<<" "<<lows[i]<<" "<<highs[i]<<"\n";std::cout<<"COUNT "<<total<<"\n";
 }catch(const std::exception&e){std::cerr<<e.what()<<"\n";return 1;}return 0;}
