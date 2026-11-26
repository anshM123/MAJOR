#include <gmpxx.h>
#include <vector>
#include <iostream>
#include <stdexcept>
#include <chrono>
using I=mpz_class; I ab(const I&x){return x>=0?x:-x;} I up(const I&n,const I&d){I q;mpz_cdiv_q(q.get_mpz_t(),n.get_mpz_t(),d.get_mpz_t());return q;} I pow10(unsigned n){I q;mpz_ui_pow_ui(q.get_mpz_t(),10,n);return q;}
int main(){try{int K,mc,mh,b,D;I den,C0,DC,H0,DH; if(!(std::cin>>K>>mc>>mh>>b>>D>>den>>C0>>DC>>H0>>DH))throw std::runtime_error("header");if(K<0||mc<1||mh<1||D<1||(b!=0&&b!=4)||DC<=0||H0<DH||DH<=0||den<=0)throw std::runtime_error("settings");std::vector<I>A(K+2),Ar(K+2),G(K+2),Gr(K+2);for(int k=0;k<=K+1;k++){if(!(std::cin>>A[k]>>Ar[k]>>G[k]>>Gr[k])||Ar[k]<0||Gr[k]<0)throw std::runtime_error("coeffs");}
int q=b*(b+1),s=2*b+1,sg=b==0?1:-1;long M=(long)(mc+1)*(mh+1);std::vector<I>X(M),V(M);I E=0,maxE=0;
for(int k=K;k>=0;k--){long k1=(long)(k+1)*(k+1),k2=-(long)k*k+(8-4*b)*(long)k+9-2*b*(b-1);I T=(k1*A[k+1]+k2*A[k])/2+q*A[k]+sg*G[k],Tr=up(k1*Ar[k+1]+ab(I(k2))*Ar[k],I(2))+1+q*Ar[k]+Gr[k];I F0=A[k]+H0*T/den,F1=DH*T/den,Fr=Ar[k]+up((H0+DH)*Tr,den)+2; // sum coefficient radius
if(k==K){X[0]=F0;X[1]=F1;E=Fr;continue;}
I p=(long)k*(k+s),a=C0-p*H0,r=p*DH,d=den-p*H0;if(d<=r)throw std::runtime_error("real pole");
for(int i=0;i<=mc;i++){long o=(long)i*(mh+1);for(int j=0;j<=mh;j++){long l=o+j;I z=a*X[l];if(i)z+=DC*X[l-(mh+1)];if(j)z+=r*(V[l-1]-X[l-1]);V[l]=z/d;}}
V[0]+=F0;V[1]+=F1;E=up((ab(a)+DC+r)*E+M*d,d-r)+Fr;X.swap(V);if(k%200==0)std::cerr<<"k "<<k<<'\n';}
I rem=0;for(long i=1;i<M;i++)rem+=ab(X[i]);std::cout<<"CENTER "<<X[0]<<"\nL1 "<<rem<<"\nERROR "<<E<<"\nLOWER "<<X[0]-rem-E<<"\n";return 0;}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
