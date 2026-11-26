// Spectral-moment-truncated positive summation transport, integer directed error.
// Input: kind M K V2, then M+1 midpoint/radius pairs. Output: m midpoint radius.
#include <gmpxx.h>
#include <iostream>
#include <vector>
using Z=mpz_class;
Z up(const Z&a,unsigned long d){Z b;mpz_cdiv_q_ui(b.get_mpz_t(),a.get_mpz_t(),d);return b;}
int main(){int kind,M,K;unsigned long V2;if(!(std::cin>>kind>>M>>K>>V2)||M<1||K<1||K>30000||(kind!=0&&kind!=4))return 1;std::vector<Z>T(M+2),S(M+2),W(M+2);Z eT=0,r;for(int j=0;j<=M;++j){if(!(std::cin>>T[j]>>r)||r<0)return 1;if(r>eT)eT=r;}
Z eS=0,eW=0,z,q,first,second,ee;std::vector<Z> g(4),ge(4);if(kind==4){for(int k=0;k<4;++k){g[k]=T[0];ge[k]=eT;}std::cout<<0<<' '<<g[3]<<' '<<ge[3]<<'\n';}
for(unsigned long m=1;m<=(unsigned long)K;++m){unsigned long den=4*m*m;Z e=up(4*m*eT+(4*V2+1)*eW,den)+1;
 for(int j=0;j<=M;++j){mpz_mul_ui(z.get_mpz_t(),T[j].get_mpz_t(),4*m);mpz_submul_ui(z.get_mpz_t(),W[j+1].get_mpz_t(),4*V2);z-=W[j];mpz_tdiv_q_ui(q.get_mpz_t(),z.get_mpz_t(),den);S[j]+=q;W[j]+=S[j];if(j==0)first=q;if(j==1)second=q;}
 eS+=e;eW+=eS;
 if(kind==0)std::cout<<m<<' '<<first<<' '<<e<<'\n';else{Z p=-(4*V2*second+first)/4;Z ep=up((4*V2+1)*e,4)+1;g[0]+=p;ge[0]+=ep;for(int j=1;j<4;++j){g[j]+=g[j-1];ge[j]+=ge[j-1];}std::cout<<m<<' '<<g[3]<<' '<<ge[3]<<'\n';}
 if(m%250==0||m==(unsigned long)K){std::cout.flush();std::cerr<<"m "<<m<<"/"<<K<<'\n';}
}return 0;}
